use std::fs::{self, File};
use std::io::Write;
use std::fmt::Formatter;
use std::fmt;

use indicatif::ProgressBar;

use crate::mmpbsa::set_style;

pub struct PDB {
    pub models: Vec<PDBModel>
}

impl fmt::Display for PDB {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "PDB with {} model(s)", self.models.len())
    }
}

impl PDB {
    pub fn new(models: &Vec<PDBModel>) -> PDB {
        PDB { models: models.to_vec() }
    }

    pub fn from(fname: &str) -> PDB {
        let f = fs::read_to_string(fname).unwrap();
        let ms: Vec<&str> = f.split("MODEL").collect();
        let mut models: Vec<PDBModel> = vec![];
        for m in ms {
            let mut t: Vec<&str> = m.split("\n").collect();
            t.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
            if t.len() > 0 {
                models.push(PDBModel::from(m));
            }
        }
        PDB { models }
    }

    pub fn to_pdb(&self, out_file_path: &str) {
        let mut pdb_file = File::create(out_file_path).unwrap();
        writeln!(pdb_file, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        let pb = ProgressBar::new(self.models.len() as u64);
        set_style(&pb);
        for model in &self.models {
            writeln!(pdb_file, "MODEL      {}", model.modelid).unwrap();
            for atom in model.atoms.iter() {
                writeln!(pdb_file, "{}", atom).unwrap();
            }
            writeln!(pdb_file, "ENDMDL").unwrap();
            pb.inc(1);
        }
        writeln!(pdb_file, "END").unwrap();
        pb.finish();
    }
}

#[derive(Clone)]
pub struct PDBModel {
    pub modelid: i32,
    pub atoms: Vec<PDBAtom>
}

impl PDBModel {
    pub fn from(model: &str) -> PDBModel {
        let mut f: Vec<&str> = model.split("\n").collect();
        let modelid = f[0].trim().parse().unwrap_or(1);
        f.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
        let mut atoms: Vec<PDBAtom> = vec![];
        for line in f {
            atoms.push(PDBAtom::from(line));
        }
        PDBModel {
            modelid,
            atoms
        }
    }

    pub fn push_atoms(&mut self, new_atoms: &Vec<PDBAtom>) {
        self.atoms.extend(new_atoms.clone());
    }
}

impl fmt::Display for PDBModel {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "PDB MODEL {}: contains {} atoms", self.modelid, self.atoms.len())
    }
}

#[derive(Clone)]
pub struct PDBAtom {
    typ: String,
    atid: i32,
    pub atname: String,
    pub resname: String,
    chainname: String,
    pub resid: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    occupy: f64,
    bf: f64,
    pub element: String,
    pub charge: String,
}

impl fmt::Display for PDBAtom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.typ.eq("ATOM") || self.typ.eq("HETATM") {
            write!(f, "{:6}{:5} {:4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:2}",
                    self.typ, self.atid, self.atname, self.resname, self.chainname, self.resid, 
                    self.x, self.y, self.z, self.occupy, self.bf, self.element, self.charge)
        }
        else {
            write!(f, "TER")
        }
    }
}

impl PDBAtom {
    pub fn from(line: &str) -> PDBAtom {
        // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
        // ATOM     69  OE2 GLU A   6      29.520 -25.258   3.929  1.00 19.99           O1-
        let typ = line[0..6].trim();
        let atid: i32 = line[9..11].trim().parse().unwrap();
        let atname = line[12..16].trim();
        let resname = line[17..20].trim();
        let chainname = line[21..22].trim();
        let resid: i32 = line[22..26].trim().parse().unwrap();
        let x: f64 = line[30..38].trim().parse().unwrap();
        let y: f64 = line[38..46].trim().parse().unwrap();
        let z: f64 = line[46..54].trim().parse().unwrap();
        let occupy: f64 = line[55..60].trim().parse().unwrap();
        let bf: f64 = line[61..66].trim().parse().unwrap();
        let element = line[70..78].trim();
        let charge = line.get(78..).unwrap_or("  ");
        let charge = if charge.ne("  ") {
            charge.trim()
        } else {
            charge
        };
        return PDBAtom {
            typ: typ.to_string(),
            atid,
            atname: atname.to_string(),
            resname: resname.to_string(),
            chainname: chainname.to_string(),
            resid,
            x,
            y,
            z,
            occupy,
            bf,
            element: element.to_string(),
            charge: charge.to_string(),
        }
    }
}