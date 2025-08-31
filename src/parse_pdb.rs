use std::fs::{self, File};
use std::path::Path;
use std::io::{BufWriter, Write};
use std::fmt::{self, Formatter};

use ndarray::Array2;
use regex::Regex;

#[derive(Clone)]
pub struct PDB {
    pub models: Vec<PDBModel>
}

impl fmt::Display for PDB {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "PDB with {} model(s)", self.models.len())
    }
}

#[allow(dead_code)]
impl PDB {
    pub fn new(models: &Vec<PDBModel>) -> PDB {
        PDB { models: models.to_vec() }
    }

    pub fn from<P: AsRef<Path>>(fname: P) -> PDB {
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

    pub fn simplify_atname(&mut self) {
        self.models.iter_mut().for_each(|m| m.atoms.iter_mut().for_each(|a| {
            // a.atname = a.atname.trim().trim_matches(char::is_numeric).to_string()
            a.atname = a.element.to_string()
        }));
    }

    pub fn to_pdb<P: AsRef<Path>>(&self, out_file_path: P) {
        let pdb_file = File::create(out_file_path).unwrap();
        let mut writer = BufWriter::new(pdb_file);
        writeln!(writer, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        for model in &self.models {
            write!(writer, "{}", model).unwrap();
        }
    }
}

#[derive(Clone)]
pub struct PDBModel {
    pub modelid: i32,
    pub atoms: Vec<PDBAtom>
}

#[allow(dead_code)]
impl PDBModel {
    pub fn from(model: &str) -> PDBModel {
        let mut f: Vec<&str> = model.split("\n").collect();
        let modelid = f[0].trim().parse().unwrap_or(1);
        f.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
        let mut atoms: Vec<PDBAtom> = vec![];
        let re_atname = Regex::new(r"[1-9][A-Z][ABGDEZH0-9][1-9']").unwrap();
        for line in f {
            atoms.push(PDBAtom::from(line, &re_atname));
        }
        PDBModel {
            modelid,
            atoms
        }
    }

    pub fn get_elements(&self) -> Vec<String> {
        self.atoms.iter().map(|a| a.element.to_string()).collect()
    }

    pub fn get_coordinates(&self) -> Array2<f64> {
        let coord: Vec<[f64; 3]> = self.atoms.iter().map(|a| [a.x, a.y, a.z]).collect();
        let coord: Vec<f64> = coord.into_iter().flatten().collect();
        Array2::from_shape_vec((self.atoms.len(), 3), coord).unwrap()
    }

    pub fn insert_atoms(&mut self, pos: usize, new_atom: &PDBAtom) {
        self.atoms.insert(pos, new_atom.clone());
    }

    pub fn push_atoms(&mut self, new_atoms: &Vec<PDBAtom>) {
        self.atoms.extend(new_atoms.clone());
    }

    pub fn to_pdb(&self, out_file_path: &str) {
        let mut pdb_file = File::create(out_file_path).unwrap();
        writeln!(pdb_file, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        for atom in self.atoms.iter() {
            writeln!(pdb_file, "{}", atom).unwrap();
        }
        writeln!(pdb_file, "END").unwrap();
    }
}

impl fmt::Display for PDBModel {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "MODEL      {}", self.modelid).unwrap();
        for atom in self.atoms.iter() {
            writeln!(f, "{}", atom).unwrap();
        }
        writeln!(f, "ENDMDL")
    }
}

#[derive(Clone)]
pub struct PDBAtom {
    typ: String,
    atid: i32,
    pub atname: String,
    pub altloc: String,
    pub resname: String,
    pub chainname: String,
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
            write!(f, "{:6}{:5} {:4}{:1}{:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:2}",
                    self.typ, self.atid, if self.atname.len() == 4 {
                        self.atname.to_string()
                    } else {
                        " ".to_string() + &self.atname
                    }, self.altloc, self.resname, self.chainname, self.resid, 
                    self.x, self.y, self.z, self.occupy, self.bf, self.element, self.charge)
        }
        else {
            write!(f, "TER")
        }
    }
}

impl PDBAtom {
    fn atname_for_pdb(re_atname: &Regex, atname: &str) -> String {
        if atname.len() == 4 {
            if !re_atname.is_match(&atname) {
                format!("{}{}", &atname[3..], &atname[..3])
            } else {
                atname.to_string()
            }
        } else {
            format!(" {:3}", atname.to_string())
        }
    }

    pub fn from(line: &str, re_atname_pdb: &Regex) -> PDBAtom {
        // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
        // ATOM     69  OE2 GLU A   6      29.520 -25.258   3.929  1.00 19.99           O1-
        let typ = line[0..6].trim();
        let atid: i32 = line[9..11].trim().parse().unwrap();
        let atname = PDBAtom::atname_for_pdb(&re_atname_pdb, line[12..16].trim());
        let altloc = line[16..17].trim();
        let resname = line[17..20].trim();
        let chainname = line[21..22].trim();
        let resid: i32 = line[22..26].trim().parse().unwrap();
        let x: f64 = line[30..38].trim().parse().unwrap();
        let y: f64 = line[38..46].trim().parse().unwrap();
        let z: f64 = line[46..54].trim().parse().unwrap();
        let occupy: f64 = line[55..60].trim().parse().unwrap();
        let bf: f64 = line[61..66].trim().parse().unwrap();
        let element = line[70..78].trim();
        let element = if element.is_empty() {
            atname[1..2].to_string()
        } else {
            element.to_string()
        };
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
            altloc: altloc.to_string(),
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