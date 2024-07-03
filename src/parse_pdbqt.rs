use std::fs;
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

impl PDBQT {
    pub fn new(fname: &str) -> PDBQT {
        let f = fs::read_to_string(fname).unwrap();
        let ms: Vec<&str> = f.split("MODEL").collect();
        let mut models: Vec<PdbqtModel> = vec![];
        for m in ms {
            let mut t: Vec<&str> = m.split("\n").collect();
            t.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
            if t.len() > 0 {
                models.push(PdbqtModel::new(m));
            }
        }
        PDBQT { models }
    }
}

pub struct PdbqtModel {
    pub atoms: Vec<PdbqtAtom>
}

impl PdbqtModel {
    pub fn new(model: &str) -> PdbqtModel {
        let mut f: Vec<&str> = model.split("\n").collect();
        // let modelid = f[0].trim().parse().unwrap_or(1);
        f.retain(|&l| l.starts_with("ATOM") || l.starts_with("HETATM"));
        let mut atoms: Vec<PdbqtAtom> = vec![];
        for line in f {
            atoms.push(PdbqtAtom::new(line));
        }
        PdbqtModel {
            atoms
        }
    }
}

#[derive(Clone)]
pub struct PdbqtAtom {
    pub atname: String,
    pub resname: String,
    pub resid: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub charge: f64,
    pub attype: String
}

impl PdbqtAtom {
    pub fn new(line: &str) -> PdbqtAtom {
        // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
        // ATOM      1  N   ALA A   2      26.338 -25.338  11.581  1.00 42.62     0.614 N 
        // let typ = line[0..6].trim();
        // let atid: i32 = line[9..11].trim().parse().unwrap();
        let atname = line[12..16].trim();
        let resname = line[17..20].trim();
        // let chainname = line[21..22].trim();
        let resid: i32 = line[22..26].trim().parse().unwrap();
        let x: f64 = line[30..38].trim().parse().unwrap();
        let y: f64 = line[38..46].trim().parse().unwrap();
        let z: f64 = line[46..54].trim().parse().unwrap();
        // let occupy: f64 = line[55..60].trim().parse().unwrap();
        // let bf: f64 = line[61..66].trim().parse().unwrap();
        let charge: f64 = line[70..76].trim().parse().unwrap();
        let attype = line[77..79].trim();
        return PdbqtAtom{
            atname: atname.to_string(),
            resname: resname.to_string(),
            resid,
            x,
            y,
            z,
            charge,
            attype: attype.to_string()
        }
    }
}