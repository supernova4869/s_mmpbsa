use std::fs::{self, File};
use std::path::Path;
use std::io::Write;
use std::fmt::Formatter;
use std::fmt;

pub struct PDBQT {
    pub models: Vec<PdbqtModel>
}

impl fmt::Display for PDBQT {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "PDBQT with {} model(s)",
        self.models.len()
    )
}
}

#[allow(dead_code)]
impl PDBQT {
    pub fn new(models: &Vec<PdbqtModel>) -> PDBQT {
        PDBQT { models: models.to_vec() }
    }

    pub fn from<P: AsRef<Path>>(fname: P) -> PDBQT {
        let f = fs::read_to_string(fname).unwrap();
        let ms: Vec<&str> = f.split("MODEL").collect();
        let mut models: Vec<PdbqtModel> = vec![];
        for m in ms {
            let mut t: Vec<&str> = m.split("\n").collect();
            t.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
            if t.len() > 0 {
                models.push(PdbqtModel::from(m));
            }
        }
        PDBQT { models }
    }

    pub fn to_pdbqt(&self, out_file_path: &str) {
        let mut pdbqt_file = File::create(out_file_path).unwrap();
        writeln!(pdbqt_file, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        for model in &self.models {
            write!(pdbqt_file, "{}", model).unwrap();
        }
        writeln!(pdbqt_file, "END").unwrap();
    }

    pub fn to_pdb(&self, out_file_path: &str) {
        let mut pdbqt_file = File::create(out_file_path).unwrap();
        writeln!(pdbqt_file, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        for model in &self.models {
            write!(pdbqt_file, "{:-}", model).unwrap();
        }
        writeln!(pdbqt_file, "END").unwrap();
    }
}

#[derive(Clone)]
pub struct PdbqtModel {
    pub modelid: i32,
    pub atoms: Vec<PdbqtAtom>
}

#[allow(dead_code)]
impl PdbqtModel {
    pub fn from(model: &str) -> PdbqtModel {
        let mut f: Vec<&str> = model.split("\n").collect();
        let modelid = f[0].trim().parse().unwrap_or(1);
        f.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
        let mut atoms: Vec<PdbqtAtom> = vec![];
        for line in f {
            atoms.push(PdbqtAtom::from(line));
        }
        PdbqtModel {
            modelid,
            atoms
        }
    }

    pub fn insert_atoms(&mut self, pos: usize, new_atom: &PdbqtAtom) {
        self.atoms.insert(pos, new_atom.clone());
    }

    pub fn to_pdbqt(&self, out_file_path: &str) {
        let mut pdbqt_file = File::create(out_file_path).unwrap();
        writeln!(pdbqt_file, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        writeln!(pdbqt_file, "{}", self).unwrap();
        writeln!(pdbqt_file, "END").unwrap();
    }

    pub fn to_pdb(&self, out_file_path: &str) {
        let mut pdb_file = File::create(out_file_path).unwrap();
        writeln!(pdb_file, "REMARK   Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
        writeln!(pdb_file, "{:-}", self).unwrap();
        writeln!(pdb_file, "END").unwrap();
    }
}

impl fmt::Display for PdbqtModel {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "MODEL      {}", self.modelid).unwrap();
        for atom in self.atoms.iter() {
            if f.sign_minus() {
                writeln!(f, "{:-}", atom).unwrap();
            } else {
                writeln!(f, "{}", atom).unwrap();
            }
        }
        write!(f, "ENDMDL")
    }
}

#[derive(Clone)]
pub struct PdbqtAtom {
    typ: String,
    atid: i32,
    pub atname: String,
    pub resname: String,
    pub chainname: String,
    pub altloc: String,
    pub resid: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    occupy: f64,
    bf: f64,
    pub charge: f64,
    pub attype: String
}

impl fmt::Display for PdbqtAtom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if self.typ.eq("ATOM") || self.typ.eq("HETATM") {
            if f.sign_minus() {
                // Write pdb
                write!(f, "{:6} {:4} {:4}{:1}{:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}           {:2}",
                    self.typ, self.atid, if self.atname.len() == 4 {
                        self.atname.to_string()
                    } else {
                        " ".to_string() + &self.atname
                    }, self.altloc, self.resname, self.chainname, self.resid, 
                    self.x, self.y, self.z, self.occupy, self.bf, if self.atname.len() == 4 {
                        self.atname[1..2].to_string()
                    } else {
                        self.atname[0..1].to_string()
                    }
                )
            } else {
                // Write pdbqt
                write!(f, "{:6} {:4} {:4}{:1}{:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}    {:6.3} {:2}",
                    self.typ, self.atid, if self.atname.len() == 4 {
                        self.atname.to_string()
                    } else {
                        " ".to_string() + &self.atname
                    }, self.altloc, self.resname, self.chainname, self.resid, 
                    self.x, self.y, self.z, self.occupy, self.bf, self.charge, self.attype
                )
            }
        } else {
            write!(f, "TER")
        }
    }
}

impl PdbqtAtom {
    pub fn from(line: &str) -> PdbqtAtom {
        // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
        // ATOM      1  N   ALA A   2      26.338 -25.338  11.581  1.00 42.62     0.614 N 
        let typ = line[0..6].trim();
        let atid: i32 = line[9..11].trim().parse().unwrap();
        let atname = line[12..16].trim();
        let altloc = line[16..17].trim();
        let resname = line[17..20].trim();
        let chainname = line[21..22].trim();
        let resid: i32 = line[22..26].trim().parse().unwrap();
        let x: f64 = line[30..38].trim().parse().unwrap();
        let y: f64 = line[38..46].trim().parse().unwrap();
        let z: f64 = line[46..54].trim().parse().unwrap();
        let occupy: f64 = line[55..60].trim().parse().unwrap();
        let bf: f64 = line[61..66].trim().parse().unwrap();
        let charge: f64 = line[70..76].trim().parse().unwrap();
        let attype = line[77..79].trim();
        return PdbqtAtom {
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
            charge,
            attype: attype.to_string()
        }
    }
}