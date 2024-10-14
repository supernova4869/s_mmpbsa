use std::{fs, io::Write};
use std::path::Path;
use std::fmt::{self, Debug, Display};

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MOL2 {
    pub mol: Molecule,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

impl MOL2 {
    pub fn new(mol: Molecule, atoms: Vec<Atom>, bonds: Vec<Bond>) -> MOL2 {
        MOL2 { mol, atoms, bonds }
    }
    pub fn from(file: &str) -> MOL2 {
        // 读取文件
        let mol2_content = fs::read_to_string(file).unwrap();
        let mut mol2_content: Vec<&str> = mol2_content.split("\n").collect();
        mol2_content.iter_mut().for_each(|s| *s = s.trim());

        // Molecule定位
        let mol_ln = mol2_content.iter().enumerate().find_map(|(index, &s)| {
            if s.eq("@<TRIPOS>MOLECULE") {
                Some(index)
            } else {
                None
            }
        }).unwrap();

        // Molecule字段
        let sys_name = Path::new(file).file_stem().unwrap().to_str().unwrap();
        let num: Vec<i32> = mol2_content[mol_ln + 2].trim().split_whitespace().map(|s| s.parse().unwrap()).collect();
        let at_num = Some(num[0]);
        let bond_num = Some(num[1]);
        let sub_struct_num = Some(num[2]);
        let prop_num = Some(num[3]);
        let set_num = Some(num[4]);
        let sys_type = Some(mol2_content[mol_ln + 3].trim().to_string());
        let at_charge = Some(mol2_content[mol_ln + 4].trim().to_string());
        let mol = Molecule::new(sys_name, at_num, bond_num, sub_struct_num, prop_num, set_num, sys_type, at_charge);

        // Atom字段
        let atom_ln = mol2_content.iter().enumerate().find_map(|(index, &s)| {
            if s.eq("@<TRIPOS>ATOM") {
                Some(index)
            } else {
                None
            }
        }).unwrap();
        let mut atoms: Vec<Atom> = Vec::new();
        for &at_line in &mol2_content[atom_ln + 1 .. atom_ln + 1 + at_num.unwrap() as usize] {
            let mut atom = Atom::from(at_line);
            atom.sub_struct_name = sys_name.to_string();
            atoms.push(atom);
        }
        
        // Bond字段
        let bond_ln = mol2_content.iter().enumerate().find_map(|(index, &s)| {
            if s.eq("@<TRIPOS>BOND") {
                Some(index)
            } else {
                None
            }
        }).unwrap();
        let mut bonds: Vec<Bond> = Vec::new();
        for &bond_line in &mol2_content[bond_ln + 1 .. bond_ln + 1 + bond_num.unwrap() as usize] {
            bonds.push(Bond::from(bond_line));
        }

        // total mol2
        MOL2::new(mol, atoms, bonds)
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Molecule {
    pub sys_name: String,
    at_num: i32,
    bond_num: i32,
    sub_struct_num: i32,
    prop_num: i32,
    set_num: i32,
    sys_type: String,
    at_charge: String,
}

impl Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "@<TRIPOS>MOLECULE\n{}\n{:5}{:6}     1 0 0\n{}\n{}\n\n", 
            self.sys_name, self.at_num, self.bond_num, self.sys_type, self.at_charge)
    }
}

impl Molecule {
    pub fn new(sys_name: &str, at_num: Option<i32>, bond_num: Option<i32>, sub_struct_num: Option<i32>, prop_num: Option<i32>, 
               set_num: Option<i32>, sys_type: Option<String>, at_charge: Option<String>) -> Molecule {
        let sys_name = sys_name.to_string();
        let at_num = at_num.unwrap_or(0);
        let bond_num = bond_num.unwrap_or(0);
        let sub_struct_num = sub_struct_num.unwrap_or(0);
        let prop_num = prop_num.unwrap_or(0);
        let set_num = set_num.unwrap_or(0);
        let sys_type = sys_type.unwrap_or("SMALL".to_string());
        let at_charge = at_charge.unwrap_or("USER_CHARGES".to_string());
        Molecule {
            sys_name, at_num, bond_num, sub_struct_num, prop_num, set_num, sys_type,at_charge
        }
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Atom {
    pub atom_id: usize,
    pub atom_name: String,
    x: f64,
    y: f64,
    z: f64,
    at: String,
    sub_struct_id: i32,
    sub_struct_name: String,
    atom_charge: Option<f64>,
    pub element: String
}

impl Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:7} {:10}{:12.4}{:12.4}{:12.4} {:7}{:3} {:9}{:8.4}\n", 
            self.atom_id, self.atom_name, self.x, self.y, self.z, self.at, 
            self.sub_struct_id, self.sub_struct_name, self.atom_charge.unwrap_or(0.0))
    }
}

impl Atom {
    fn from(line: &str) -> Atom {
        let line: Vec<&str> = line.trim().split_whitespace().collect();
        let atom_id: usize = line[0].parse().unwrap();
        let atom_name: String = line[1].to_string();
        let x: f64 = line[2].parse().unwrap();
        let y: f64 = line[3].parse().unwrap();
        let z: f64 = line[4].parse().unwrap();
        let at: String = line[5].to_string();
        // 修改残基名为MOL, 残基编号为1
        let sub_struct_id = 1;
        let sub_struct_name = "MOL".to_string();
        let atom_charge: Option<f64> = Some(line[8].parse().unwrap());
        let element: Vec<&str> = at.split(".").collect();
        let element = element[0].to_string();
        Atom {
            atom_id, atom_name, x, y, z, at, sub_struct_id, sub_struct_name, atom_charge, element
        }
    }
}

#[derive(Debug)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct Bond {
    bond_id: usize,
    pub a1: usize,
    pub a2: usize,
    bt: String,
}

impl Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:6}{:5}{:5} {}\n", 
            self.bond_id, self.a1, self.a2, self.bt)
    }
}

impl Bond {
    fn from(line: &str) -> Bond {
        let line: Vec<&str> = line.trim().split_whitespace().collect();
        let bond_id: usize = line[0].parse().unwrap();
        let a1: usize = line[1].parse().unwrap();
        let a2: usize = line[2].parse().unwrap();
        let bt: String = line[3].to_string();
        Bond {
            bond_id, a1, a2, bt
        }
    }
}

impl Display for MOL2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = format!("; Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)\n\n");
        out.push_str(format!("{}", self.mol).as_str());
        out.push_str("@<TRIPOS>ATOM\n");
        for a in &self.atoms {
            out.push_str(format!("{}", a).as_str());
        }
        out.push_str("@<TRIPOS>BOND\n");
        for b in &self.bonds {
            out.push_str(format!("{}", b).as_str());
        }
        write!(f, "{}", out)
    }
}

#[allow(dead_code)]
impl MOL2 {
    pub fn output(&self, outfile: &str) {
        let mut file = fs::File::create(outfile).unwrap();
        file.write_fmt(format_args!("{}", self)).unwrap();
        println!("Written to {}", outfile);
    }

    pub fn to_chg(&self, outfile: &str) {
        let mut file = fs::File::create(outfile).unwrap();
        for atom in &self.atoms {
            // H    -7.181892   -0.715999   -1.479709   0.0934874878
            if let Some(c) = atom.atom_charge {
                writeln!(file, " {:4}{:12.6}{:12.6}{:12.6}{:15.10}", atom.atom_name, atom.x, atom.y, atom.z, c).unwrap();
            } else {
                writeln!(file, " {:4}{:12.6}{:12.6}{:12.6}{:15.10}", atom.atom_name, atom.x, atom.y, atom.z, 0.0).unwrap();
            }
        }
        println!("Written to {}", outfile);
    }
}
