// APBS vparam.rs - Parameter database for charge/radius assignment
// Port of src/generic/vparam.h / vparam.c

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::error::{ApbsError, ApbsResult};
use crate::vstring::strcasecmp_eq;
use quick_xml::events::Event;
use quick_xml::Reader;

/// Atom data entry in the parameter database
#[derive(Debug, Clone)]
pub struct VparamAtomData {
    pub atom_name: String,
    pub res_name: String,
    pub charge: f64,
    pub radius: f64,
    pub epsilon: f64,
}

impl Default for VparamAtomData {
    fn default() -> Self {
        Self {
            atom_name: String::new(),
            res_name: String::new(),
            charge: 0.0,
            radius: 0.0,
            epsilon: 0.0,
        }
    }
}

/// Residue data entry containing atom data
#[derive(Debug, Clone)]
pub struct VparamResData {
    pub name: String,
    pub atom_data: Vec<VparamAtomData>,
}

/// Parameter database for charge/radius assignment
#[derive(Debug)]
pub struct Vparam {
    /// Residue data indexed by residue name (case-insensitive lookup)
    pub res_data: HashMap<String, VparamResData>,
}

impl Default for Vparam {
    fn default() -> Self {
        Self {
            res_data: HashMap::new(),
        }
    }
}

impl Vparam {
    pub fn new() -> Self {
        Self::default()
    }

    /// Look up residue data by name (case-insensitive)
    pub fn get_res_data(&self, res_name: &str) -> Option<&VparamResData> {
        self.res_data
            .iter()
            .find(|(k, _)| strcasecmp_eq(k, res_name))
            .map(|(_, v)| v)
    }

    /// Look up atom data by residue name and atom name (case-insensitive)
    pub fn get_atom_data(&self, res_name: &str, atom_name: &str) -> Option<&VparamAtomData> {
        let res = self.get_res_data(res_name)?;
        res.atom_data
            .iter()
            .find(|a| strcasecmp_eq(&a.atom_name, atom_name))
    }

    /// Read flat-format parameter file
    /// Format: RESIDUE ATOM CHARGE RADIUS EPSILON
    pub fn read_flat_file(&mut self, filename: &str) -> ApbsResult<()> {
        let file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let reader = BufReader::new(file);

        let mut atoms: Vec<VparamAtomData> = Vec::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result
                .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with('%') {
                continue;
            }

            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() < 5 {
                return Err(ApbsError::Parse {
                    line: line_num + 1,
                    message: format!("Expected 5 fields, got {}", fields.len()),
                });
            }

            let charge: f64 = fields[2].parse().map_err(|_| ApbsError::Parse {
                line: line_num + 1,
                message: format!("Invalid charge: {}", fields[2]),
            })?;
            let radius: f64 = fields[3].parse().map_err(|_| ApbsError::Parse {
                line: line_num + 1,
                message: format!("Invalid radius: {}", fields[3]),
            })?;
            let epsilon: f64 = fields[4].parse().map_err(|_| ApbsError::Parse {
                line: line_num + 1,
                message: format!("Invalid epsilon: {}", fields[4]),
            })?;

            atoms.push(VparamAtomData {
                res_name: fields[0].to_string(),
                atom_name: fields[1].to_string(),
                charge,
                radius,
                epsilon,
            });
        }

        // Group atoms into residues
        for atom in atoms {
            let res_name_upper = atom.res_name.to_uppercase();
            let entry = self.res_data
                .entry(res_name_upper.clone())
                .or_insert_with(|| VparamResData {
                    name: res_name_upper,
                    atom_data: Vec::new(),
                });
            entry.atom_data.push(atom);
        }

        Ok(())
    }

    /// Read XML-format parameter file
    ///
    /// Supported shapes:
    /// 1) `<atom res_name="ALA" atom_name="CA" charge="..." radius="..." epsilon="..."/>`
    /// 2) `<residue name="ALA"><atom name="CA" charge="..." radius="..." epsilon="..."/></residue>`
    pub fn read_xml_file(&mut self, filename: &str) -> ApbsResult<()> {
        let file = File::open(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;
        let mut reader = Reader::from_reader(BufReader::new(file));
        reader.trim_text(true);

        let mut buf = Vec::new();
        let mut cur_res: Option<String> = None;

        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(e)) => {
                    let name = e.name().as_ref().to_ascii_lowercase();
                    if name.as_slice() == b"residue" {
                        for a in e.attributes().flatten() {
                            let k = a.key.as_ref().to_ascii_lowercase();
                            if k.as_slice() == b"name" || k.as_slice() == b"res_name" {
                                cur_res = Some(String::from_utf8_lossy(a.value.as_ref()).to_string());
                            }
                        }
                    } else if name.as_slice() == b"atom" {
                        if let Some(atom) = Self::parse_atom_xml(&e, cur_res.as_deref())? {
                            self.insert_atom(atom);
                        }
                    }
                }
                Ok(Event::Empty(e)) => {
                    let name = e.name().as_ref().to_ascii_lowercase();
                    if name.as_slice() == b"residue" {
                        for a in e.attributes().flatten() {
                            let k = a.key.as_ref().to_ascii_lowercase();
                            if k.as_slice() == b"name" || k.as_slice() == b"res_name" {
                                cur_res = Some(String::from_utf8_lossy(a.value.as_ref()).to_string());
                            }
                        }
                    } else if name.as_slice() == b"atom" {
                        if let Some(atom) = Self::parse_atom_xml(&e, cur_res.as_deref())? {
                            self.insert_atom(atom);
                        }
                    }
                }
                Ok(Event::End(e)) => {
                    let name = e.name().as_ref().to_ascii_lowercase();
                    if name.as_slice() == b"residue" {
                        cur_res = None;
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
        Ok(())
    }

    /// Memory usage in bytes
    pub fn mem_chk(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.res_data.iter().map(|(k, v)| {
                k.len() + v.name.len() + v.atom_data.len() * std::mem::size_of::<VparamAtomData>()
            }).sum::<usize>()
    }
}

impl Vparam {
    fn insert_atom(&mut self, atom: VparamAtomData) {
        let res_name_upper = atom.res_name.to_uppercase();
        let entry = self
            .res_data
            .entry(res_name_upper.clone())
            .or_insert_with(|| VparamResData {
                name: res_name_upper,
                atom_data: Vec::new(),
            });
        entry.atom_data.push(atom);
    }

    fn parse_atom_xml(
        e: &quick_xml::events::BytesStart<'_>,
        inherited_res: Option<&str>,
    ) -> ApbsResult<Option<VparamAtomData>> {
        let mut res_name = inherited_res.unwrap_or("").to_string();
        let mut atom_name = String::new();
        let mut charge = None;
        let mut radius = None;
        let mut epsilon = None;

        for a in e.attributes().flatten() {
            let key = String::from_utf8_lossy(a.key.as_ref()).to_ascii_lowercase();
            let val = String::from_utf8_lossy(a.value.as_ref()).to_string();
            match key.as_str() {
                "res_name" | "residue" | "res" => res_name = val,
                "atom_name" | "name" | "atom" => atom_name = val,
                "charge" => charge = val.parse::<f64>().ok(),
                "radius" => radius = val.parse::<f64>().ok(),
                "epsilon" | "eps" => epsilon = val.parse::<f64>().ok(),
                _ => {}
            }
        }

        if res_name.is_empty() || atom_name.is_empty() || charge.is_none() || radius.is_none() {
            return Ok(None);
        }

        Ok(Some(VparamAtomData {
            atom_name,
            res_name,
            charge: charge.unwrap_or(0.0),
            radius: radius.unwrap_or(0.0),
            epsilon: epsilon.unwrap_or(0.0),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_param_flat_file() {
        let tmpfile = std::env::temp_dir().join("test_apbs.param");
        let mut f = File::create(&tmpfile).unwrap();
        writeln!(f, "ALA CA -0.1800 1.9000 0.1100").unwrap();
        writeln!(f, "ALA CB  0.0300 1.9000 0.1100").unwrap();
        writeln!(f, "GLY CA -0.1800 1.9000 0.1100").unwrap();
        drop(f);

        let mut params = Vparam::new();
        params.read_flat_file(tmpfile.to_str().unwrap()).unwrap();

        // Look up ALA CA
        let atom = params.get_atom_data("ALA", "CA").unwrap();
        assert!((atom.charge - (-0.18)).abs() < 1e-10);
        assert!((atom.radius - 1.9).abs() < 1e-10);

        // Case insensitive
        let atom2 = params.get_atom_data("ala", "cb").unwrap();
        assert!((atom2.charge - 0.03).abs() < 1e-10);

        std::fs::remove_file(tmpfile).ok();
    }

    #[test]
    fn test_param_xml_file() {
        let tmpfile = std::env::temp_dir().join("test_apbs.param.xml");
        let mut f = File::create(&tmpfile).unwrap();
        writeln!(
            f,
            "<params><residue name=\"ALA\"><atom name=\"CA\" charge=\"-0.18\" radius=\"1.9\" epsilon=\"0.11\"/></residue></params>"
        )
        .unwrap();
        drop(f);

        let mut params = Vparam::new();
        params.read_xml_file(tmpfile.to_str().unwrap()).unwrap();
        let atom = params.get_atom_data("ALA", "CA").unwrap();
        assert!((atom.charge - (-0.18)).abs() < 1e-12);
        assert!((atom.radius - 1.9).abs() < 1e-12);
        std::fs::remove_file(tmpfile).ok();
    }
}
