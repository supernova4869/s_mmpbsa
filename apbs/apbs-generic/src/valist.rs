// APBS valist.rs - Atom list management
// Port of src/generic/valist.h / valist.c

use std::io::{BufRead, BufReader};
use std::fs::File;

use crate::error::{ApbsError, ApbsResult};
use crate::vatom::Vatom;
use crate::vparam::Vparam;
use quick_xml::events::Event;
use quick_xml::Reader;

/// Container for a list of atom objects
#[derive(Debug, Clone)]
pub struct Valist {
    /// Number of atoms
    pub number: usize,
    /// Geometric center of molecule
    pub center: [f64; 3],
    /// Minimum coordinates across all atoms
    pub mincrd: [f64; 3],
    /// Maximum coordinates across all atoms
    pub maxcrd: [f64; 3],
    /// Maximum atomic radius
    pub maxrad: f64,
    /// Net charge of molecule
    pub charge: f64,
    /// Atom list
    pub atoms: Vec<Vatom>,
}

impl Default for Valist {
    fn default() -> Self {
        Self {
            number: 0,
            center: [0.0; 3],
            mincrd: [f64::MAX; 3],
            maxcrd: [f64::MIN; 3],
            maxrad: 0.0,
            charge: 0.0,
            atoms: Vec::new(),
        }
    }
}

impl Valist {
    pub fn new() -> Self {
        Self::default()
    }

    /// Get number of atoms
    #[inline]
    pub fn number_atoms(&self) -> usize {
        self.number
    }

    /// Get atom by index
    #[inline]
    pub fn get_atom(&self, i: usize) -> &Vatom {
        &self.atoms[i]
    }

    /// Get mutable atom by index
    #[inline]
    pub fn get_atom_mut(&mut self, i: usize) -> &mut Vatom {
        &mut self.atoms[i]
    }

    /// Get center X coordinate
    #[inline]
    pub fn center_x(&self) -> f64 {
        self.center[0]
    }

    /// Get center Y coordinate
    #[inline]
    pub fn center_y(&self) -> f64 {
        self.center[1]
    }

    /// Get center Z coordinate
    #[inline]
    pub fn center_z(&self) -> f64 {
        self.center[2]
    }

    /// Get the atoms slice
    pub fn atoms(&self) -> &[Vatom] {
        &self.atoms
    }

    /// Get mutable atoms slice
    pub fn atoms_mut(&mut self) -> &mut [Vatom] {
        &mut self.atoms
    }

    /// Compute statistics (center, bounds, maxrad, charge)
    pub fn get_statistics(&mut self) {
        if self.atoms.is_empty() {
            return;
        }

        self.mincrd = [f64::MAX; 3];
        self.maxcrd = [f64::MIN; 3];
        self.maxrad = 0.0;
        self.charge = 0.0;

        for atom in &self.atoms {
            for d in 0..3 {
                if atom.position[d] < self.mincrd[d] {
                    self.mincrd[d] = atom.position[d];
                }
                if atom.position[d] > self.maxcrd[d] {
                    self.maxcrd[d] = atom.position[d];
                }
            }
            if atom.radius > self.maxrad {
                self.maxrad = atom.radius;
            }
            self.charge += atom.charge;
        }

        for d in 0..3 {
            self.center[d] = (self.mincrd[d] + self.maxcrd[d]) / 2.0;
        }

        self.number = self.atoms.len();
    }

    /// Read PQR file (whitespace-delimited ATOM/HETATM lines)
    pub fn read_pqr(&mut self, params: Option<&Vparam>, filename: &str) -> ApbsResult<()> {
        let file = File::open(filename).map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let reader = BufReader::new(file);

        self.atoms.clear();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result.map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

            let trimmed = line.trim();
            if trimmed.starts_with("ATOM") || trimmed.starts_with("HETATM") {
                let mut atom = Self::parse_pqr_line(trimmed, line_num + 1, params)?;
                atom.set_atom_id(self.atoms.len() as i32);
                self.atoms.push(atom);
            }
        }

        self.get_statistics();
        Ok(())
    }

    /// Parse a single PQR line
    fn parse_pqr_line(line: &str, line_num: usize, params: Option<&Vparam>) -> ApbsResult<Vatom> {
        let fields: Vec<&str> = line.split_whitespace().collect();

        // PQR has both chain-less and chain-containing variants:
        // ATOM serial atomName resName resSeq x y z charge radius
        // ATOM serial atomName resName chain resSeq x y z charge radius
        // Parse numeric payload from the tail so both layouts are handled.
        if fields.len() < 9 {
            return Err(ApbsError::Parse {
                line: line_num,
                message: format!("Expected at least 9 fields, got {}", fields.len()),
            });
        }

        let atom_name = fields[2];
        let res_name = fields[3];
        let tail = fields.len();

        let x: f64 = fields[tail - 5].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid x coordinate: {}", fields[tail - 5]),
        })?;
        let y: f64 = fields[tail - 4].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid y coordinate: {}", fields[tail - 4]),
        })?;
        let z: f64 = fields[tail - 3].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid z coordinate: {}", fields[tail - 3]),
        })?;
        let charge: f64 = fields[tail - 2].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid charge: {}", fields[tail - 2]),
        })?;
        let radius: f64 = fields[tail - 1].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid radius: {}", fields[tail - 1]),
        })?;

        let mut atom = Vatom::new();
        atom.set_position([x, y, z]);
        atom.set_atom_name(atom_name);
        atom.set_res_name(res_name);

        // If parameters are provided, look up and override charge/radius/epsilon
        if let Some(p) = params {
            if let Some(adata) = p.get_atom_data(res_name, atom_name) {
                atom.set_charge(adata.charge);
                atom.set_radius(adata.radius);
                atom.set_epsilon(adata.epsilon);
            } else {
                atom.set_charge(charge);
                atom.set_radius(radius);
            }
        } else {
            atom.set_charge(charge);
            atom.set_radius(radius);
        }

        Ok(atom)
    }

    /// Read PDB file (whitespace-delimited ATOM/HETATM lines)
    pub fn read_pdb(&mut self, params: Option<&Vparam>, filename: &str) -> ApbsResult<()> {
        let file = File::open(filename).map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let reader = BufReader::new(file);

        self.atoms.clear();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result.map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

            let trimmed = line.trim();
            if trimmed.starts_with("ATOM") || trimmed.starts_with("HETATM") {
                let mut atom = Self::parse_pdb_line(trimmed, line_num + 1, params)?;
                atom.set_atom_id(self.atoms.len() as i32);
                self.atoms.push(atom);
            }
        }

        self.get_statistics();
        Ok(())
    }

    /// Parse a single PDB line
    fn parse_pdb_line(line: &str, line_num: usize, params: Option<&Vparam>) -> ApbsResult<Vatom> {
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() < 5 {
            return Err(ApbsError::Parse {
                line: line_num,
                message: format!("Expected at least 5 fields for PDB, got {}", fields.len()),
            });
        }

        let atom_name = fields[2];
        let res_name = fields[3];

        // Find x, y, z - they follow the residue sequence number
        // In PDB format: ATOM serial atomName resName resSeq x y z ...
        // We need fields[4] = resSeq, fields[5] = x, fields[6] = y, fields[7] = z
        if fields.len() < 8 {
            return Err(ApbsError::Parse {
                line: line_num,
                message: format!("Expected at least 8 fields for PDB coordinates, got {}", fields.len()),
            });
        }

        let x: f64 = fields[5].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid x coordinate: {}", fields[5]),
        })?;
        let y: f64 = fields[6].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid y coordinate: {}", fields[6]),
        })?;
        let z: f64 = fields[7].parse().map_err(|_| ApbsError::Parse {
            line: line_num,
            message: format!("Invalid z coordinate: {}", fields[7]),
        })?;

        let mut atom = Vatom::new();
        atom.set_position([x, y, z]);
        atom.set_atom_name(atom_name);
        atom.set_res_name(res_name);

        // Look up charge/radius from parameter database
        if let Some(p) = params {
            if let Some(adata) = p.get_atom_data(res_name, atom_name) {
                atom.set_charge(adata.charge);
                atom.set_radius(adata.radius);
                atom.set_epsilon(adata.epsilon);
            }
        }

        Ok(atom)
    }

    /// Read XML file
    ///
    /// Supported:
    /// `<atom .../>` with attributes:
    /// - name/atom_name, res_name/residue
    /// - x,y,z
    /// - charge, radius, epsilon (optional if params can fill)
    pub fn read_xml(&mut self, params: Option<&Vparam>, filename: &str) -> ApbsResult<()> {
        let file = File::open(filename).map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let mut reader = Reader::from_reader(BufReader::new(file));
        reader.trim_text(true);
        let mut buf = Vec::new();

        self.atoms.clear();

        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(e)) | Ok(Event::Empty(e)) => {
                    let tag = e.name().as_ref().to_ascii_lowercase();
                    if tag.as_slice() == b"atom" {
                        if let Some(mut atom) = Self::parse_xml_atom(&e, params)? {
                            atom.set_atom_id(self.atoms.len() as i32);
                            self.atoms.push(atom);
                        }
                    }
                }
                Ok(Event::Eof) => break,
                Err(e) => {
                    return Err(ApbsError::Parse {
                        line: 0,
                        message: format!("XML parse error in {}: {}", filename, e),
                    })
                }
                _ => {}
            }
            buf.clear();
        }

        self.get_statistics();
        Ok(())
    }
}

impl Valist {
    fn parse_xml_atom(
        e: &quick_xml::events::BytesStart<'_>,
        params: Option<&Vparam>,
    ) -> ApbsResult<Option<Vatom>> {
        let mut atom_name = String::new();
        let mut res_name = String::new();
        let mut x = None;
        let mut y = None;
        let mut z = None;
        let mut charge = None;
        let mut radius = None;
        let mut epsilon = None;

        for a in e.attributes().flatten() {
            let key = String::from_utf8_lossy(a.key.as_ref()).to_ascii_lowercase();
            let val = String::from_utf8_lossy(a.value.as_ref()).to_string();
            match key.as_str() {
                "name" | "atom_name" | "atom" => atom_name = val,
                "res_name" | "residue" | "res" => res_name = val,
                "x" => x = val.parse::<f64>().ok(),
                "y" => y = val.parse::<f64>().ok(),
                "z" => z = val.parse::<f64>().ok(),
                "charge" => charge = val.parse::<f64>().ok(),
                "radius" => radius = val.parse::<f64>().ok(),
                "epsilon" | "eps" => epsilon = val.parse::<f64>().ok(),
                _ => {}
            }
        }

        if atom_name.is_empty() || res_name.is_empty() || x.is_none() || y.is_none() || z.is_none() {
            return Ok(None);
        }

        let mut atom = Vatom::new();
        atom.set_atom_name(&atom_name);
        atom.set_res_name(&res_name);
        atom.set_position([x.unwrap_or(0.0), y.unwrap_or(0.0), z.unwrap_or(0.0)]);

        if let Some(p) = params {
            if let Some(adata) = p.get_atom_data(&res_name, &atom_name) {
                atom.set_charge(adata.charge);
                atom.set_radius(adata.radius);
                atom.set_epsilon(adata.epsilon);
            } else {
                atom.set_charge(charge.unwrap_or(0.0));
                atom.set_radius(radius.unwrap_or(0.0));
                if let Some(eps) = epsilon {
                    atom.set_epsilon(eps);
                }
            }
        } else {
            atom.set_charge(charge.unwrap_or(0.0));
            atom.set_radius(radius.unwrap_or(0.0));
            if let Some(eps) = epsilon {
                atom.set_epsilon(eps);
            }
        }

        Ok(Some(atom))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_valist_pqr() {
        let tmpfile = std::env::temp_dir().join("test_apbs.pqr");
        let mut f = File::create(&tmpfile).unwrap();
        writeln!(f, "ATOM      1  CA  ALA     1       1.000   2.000   3.000 -0.1800  1.9000").unwrap();
        writeln!(f, "ATOM      2  CB  ALA     1       2.000   3.000   4.000  0.0300  1.9000").unwrap();
        writeln!(f, "END").unwrap();
        drop(f);

        let mut valist = Valist::new();
        valist.read_pqr(None, tmpfile.to_str().unwrap()).unwrap();

        assert_eq!(valist.number_atoms(), 2);
        assert!((valist.charge - (-0.15)).abs() < 1e-10);

        let a0 = valist.get_atom(0);
        assert_eq!(a0.atom_name(), "CA");
        assert!((a0.position()[0] - 1.0).abs() < 1e-10);
        assert!((a0.charge - (-0.18)).abs() < 1e-10);

        std::fs::remove_file(tmpfile).ok();
    }

    #[test]
    fn test_valist_pqr_with_chain_id() {
        let tmpfile = std::env::temp_dir().join("test_apbs_chain.pqr");
        let mut f = File::create(&tmpfile).unwrap();
        writeln!(f, "ATOM      1  N   SER A   1      58.907   62.163   46.267     0.184900     1.550000").unwrap();
        writeln!(f, "ATOM   4681  N   LIG A 296      32.021   38.906   40.672    -0.901800     1.550000").unwrap();
        drop(f);

        let mut valist = Valist::new();
        valist.read_pqr(None, tmpfile.to_str().unwrap()).unwrap();

        assert_eq!(valist.number_atoms(), 2);
        let a0 = valist.get_atom(0);
        assert_eq!(a0.atom_name(), "N");
        assert_eq!(a0.res_name(), "SER");
        assert!((a0.position()[0] - 58.907).abs() < 1e-10);
        assert!((a0.position()[1] - 62.163).abs() < 1e-10);
        assert!((a0.position()[2] - 46.267).abs() < 1e-10);
        assert!((a0.charge - 0.1849).abs() < 1e-10);
        assert!((a0.radius - 1.55).abs() < 1e-10);

        let a1 = valist.get_atom(1);
        assert!((a1.position()[0] - 32.021).abs() < 1e-10);
        assert!((a1.charge - (-0.9018)).abs() < 1e-10);
        assert!((a1.radius - 1.55).abs() < 1e-10);

        std::fs::remove_file(tmpfile).ok();
    }

    #[test]
    fn test_valist_xml() {
        let tmpfile = std::env::temp_dir().join("test_apbs.xml");
        let mut f = File::create(&tmpfile).unwrap();
        writeln!(
            f,
            "<molecule><atom name=\"CA\" res_name=\"ALA\" x=\"1.0\" y=\"2.0\" z=\"3.0\" charge=\"-0.18\" radius=\"1.9\"/></molecule>"
        )
        .unwrap();
        drop(f);

        let mut valist = Valist::new();
        valist.read_xml(None, tmpfile.to_str().unwrap()).unwrap();
        assert_eq!(valist.number_atoms(), 1);
        let a0 = valist.get_atom(0);
        assert_eq!(a0.atom_name(), "CA");
        assert!((a0.charge - (-0.18)).abs() < 1e-12);
        std::fs::remove_file(tmpfile).ok();
    }
}
