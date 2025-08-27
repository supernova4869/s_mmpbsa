use std::fmt::Formatter;
use std::path::Path;
use std::fmt;
use ndarray::Array2;
use std::io::{BufRead, BufReader, BufWriter, Write};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fs::File;
use once_cell::sync::Lazy;

use crate::settings::Settings;

// 预编译所有正则表达式
static RE_NAME: Lazy<Regex> = Lazy::new(|| Regex::new("name\\s*=\\s*\"(.*)\"").unwrap());
static RE_ATOMS_NUM: Lazy<Regex> = Lazy::new(|| Regex::new(r"#atoms\s*=\s*(\d+)").unwrap());
static RE_MOLBLOCK: Lazy<Regex> = Lazy::new(|| Regex::new(r"#molblock\s*=\s*(\d+)").unwrap());
static RE_MOLTYPE: Lazy<Regex> = Lazy::new(|| Regex::new("moltype\\s*=\\s*\\d+\\s*\"(.*)\"").unwrap());
static RE_MOLECULES_NUM: Lazy<Regex> = Lazy::new(|| Regex::new(r"#molecules\s*=\s*(\d+)").unwrap());
static RE_ATNR: Lazy<Regex> = Lazy::new(|| Regex::new(r"atnr\s*=\s*(\d+)").unwrap());
static RE_FUNCTYPE: Lazy<Regex> = Lazy::new(|| 
    Regex::new(r"functype\[(\d*)]=LJ_SR,\s*c6\s*=\s*(\S*),\s*c12\s*=\s*(\S*)").unwrap());
static RE_MOLTYPE_ID: Lazy<Regex> = Lazy::new(|| Regex::new(r"moltype \((\d+)\)").unwrap());
static RE_ATOM_PARAMS: Lazy<Regex> = Lazy::new(|| 
    Regex::new(r".*type=\s*(\d+).*q=\s*([^,]+),.*resind=\s*(\d+).*").unwrap());
static RE_ATOM_NAME: Lazy<Regex> = Lazy::new(|| Regex::new("name=\"(.*)\"").unwrap());
static RE_TYPE_NAME: Lazy<Regex> = Lazy::new(|| Regex::new("name=\"(.*)\",").unwrap());
static RE_RESIDUE: Lazy<Regex> = Lazy::new(|| 
    Regex::new("residue\\[(\\d+)]=\\{name=\"(.+)\",.*nr=([\\d\\-]+).*").unwrap());
static RE_ANGLES_NUM: Lazy<Regex> = Lazy::new(|| Regex::new(r"nr\s*:\s*(\d+)").unwrap());
static RE_ANGLES: Lazy<Regex> = Lazy::new(|| Regex::new(r"\(ANGLES\)\s+(\d+)\s+(\d+)\s+(\d+)").unwrap());
static RE_COORD: Lazy<Regex> = Lazy::new(|| Regex::new(r"\{\s*(.*),\s*(.*),\s*(.*)\}").unwrap());
static RE_ATOM_COUNT: Lazy<Regex> = Lazy::new(|| Regex::new(r"atom \((\d+)\):").unwrap());
static RE_RESIDUE_COUNT: Lazy<Regex> = Lazy::new(|| Regex::new(r"residue \((\d+)\)").unwrap());

pub struct TPR {
    pub name: String,
    pub n_atoms: usize,
    pub molecule_types_num: usize,
    pub molecule_types: Vec<MolType>,
    pub atom_types_num: usize,
    pub lj_sr_params: Vec<LJType>,
    pub molecules: Vec<Molecule>,
    pub temp: f64,
    pub coordinates: Array2<f64>,
}

impl fmt::Display for TPR {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}, with\n{} atoms,\n{} type(s) of molecules,\n{} atom types,\n{} LJ types,\n{}x{} coordinate",
               self.name, self.n_atoms, self.molecule_types_num,
               self.atom_types_num, self.lj_sr_params.len(), 
               self.coordinates.shape()[0], self.coordinates.shape()[1]
        )
    }
}

impl TPR {
    pub fn from(mdp: &str, settings: &Settings) -> TPR {
        let mut name = String::new();
        let mut atoms_num = 0;
        let mut molecule_types_num = 0;
        let mut atom_types_num = 0;
        let mut molecule_types: Vec<MolType> = Vec::new();

        // 使用BufReader高效读取文件
        let file = File::open(mdp).unwrap();
        let mut reader = BufReader::new(file);
        let mut buf = String::with_capacity(256);

        let mut fun_type: Vec<LJType> = Vec::new();
        let mut sigma: Vec<f64> = Vec::new();
        let mut epsilon: Vec<f64> = Vec::new();
        let mut radius: Vec<f64> = Vec::new();

        // 预分配向量容量（稍后根据实际大小调整）
        let mut atom_resids: Vec<usize> = Vec::new();
        let mut atom_types: Vec<usize> = Vec::new();
        let mut atom_radii: Vec<f64> = Vec::new();
        let mut atom_charges: Vec<f64> = Vec::new();
        let mut atom_names: Vec<String> = Vec::new();
        let mut type_names: Vec<String> = Vec::new();
        let mut coordinates: Vec<f64> = Vec::new();
        let mut molecules: Vec<Molecule> = Vec::new();

        // 模拟时间参数
        let mut temperature = 0.0;

        println!("Loading dump file: {}\n", mdp);
        
        // 读取文件行
        while read_line(&mut reader, &mut buf) > 0 {
            let line = buf.trim();
            
            // 使用模式匹配提高可读性和性能
            if line.starts_with("ref-t:") {
                process_ref_t(line, &mut temperature);
            } else if line.starts_with("topology:") {
                process_topology(&mut reader, &mut buf, &mut name, &mut atoms_num, 
                            &mut molecule_types_num, &mut molecule_types);
            } else if line.starts_with("ffparams:") {
                process_ffparams(&mut reader, &mut buf, &mut atom_types_num, 
                            &mut fun_type, &mut sigma, &mut epsilon, 
                            &mut radius, settings);
            } else if line.starts_with("moltype (") {
                process_moltype(&mut reader, &mut buf, &mut molecules, 
                            &mut atom_resids, &mut atom_types, &mut atom_radii,
                            &mut atom_charges, &mut atom_names, &mut type_names,
                            &radius);
            } else if line.starts_with("x (") {
                process_coordinates(&mut reader, &mut buf, atoms_num, &mut coordinates);
            }
        }

        println!("Backup force field radius...");
        let ff_dat = Path::new(mdp).parent().unwrap().join("ff_radius.dat");
        if ff_dat.is_file() {
            std::fs::remove_file(&ff_dat).unwrap();
        }
        let ff_file = File::create(ff_dat).unwrap();
        let mut writer = BufWriter::new(ff_file);
        for r in &atom_radii {
            writeln!(writer, "{:.2}", r).unwrap();
        }
        writer.flush().unwrap();

        println!("System molecular composition:");
        for mol in &molecules {
            println!("Molecule {}: {}", mol.molecule_type_id, mol);
        }

        TPR {
            name,
            n_atoms: atoms_num,
            molecule_types_num,
            molecule_types,
            atom_types_num,
            lj_sr_params: fun_type,
            molecules,
            temp: temperature,
            coordinates: Array2::from_shape_vec((atoms_num, 3), coordinates).unwrap()
        }
    }
}

pub struct MolType {
    pub id: usize,
    pub name: String,
    pub molecules_num: i64,
}

impl MolType {
    fn new(id: usize, name: String, molecules_num: i64) -> MolType {
        MolType {
            id,
            name,
            molecules_num,
        }
    }
}

impl fmt::Display for MolType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Molecule type {}: {}, #{} in the system", self.id, self.name, self.molecules_num)
    }
}

pub struct LJType {
    pub c6: f64,
    pub c12: f64,
}

impl LJType {
    fn new(c6: f64, c12: f64) -> LJType {
        LJType {
            c6,
            c12,
        }
    }
}

pub struct Molecule {
    pub molecule_type_id: usize,
    pub molecule_name: String,
    pub atoms_num: usize,
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
}

impl Molecule {
    fn new(molecule_type_id: usize, molecule_name: String, atoms_num: usize,
           atoms: &Vec<Atom>, residues: &Vec<Residue>) -> Molecule {
        Molecule {
            molecule_type_id,
            molecule_name,
            atoms_num,
            atoms: atoms.to_vec(),
            residues: residues.to_vec(),
        }
    }
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Molecule {}, with {} atoms, {} residues",
               self.molecule_name, self.atoms_num, self.residues.len())
    }
}

#[derive(Clone)]
pub struct Atom {
    pub id: usize,
    pub at_type: String,
    pub type_id: usize,
    pub charge: f64,
    pub resind: usize,
    pub name: String,
    pub radius: f64,
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Atom {}: {} with type {}, charge {}, radius {}, in residue {}",
               self.id, self.name, self.type_id, self.charge, self.radius, self.resind)
    }
}

impl Atom {
    fn new(id: usize, at_type: &str, type_id: usize, charge: f64, residue_index: usize, name: String,
           radius: f64) -> Atom {
        Atom {
            id,
            at_type: at_type.to_string(),
            type_id,
            charge,
            resind: residue_index,
            name,
            radius,
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Residue {
    pub id: usize,
    pub name: String,
    pub nr: i32,
}

impl fmt::Display for Residue {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Residue {}: {} with residue index {}", self.id, self.name, self.nr)
    }
}

impl Residue {
    pub fn new(id: usize, name: String, nr: i32) -> Residue {
        Residue {
            id,
            name,
            nr,
        }
    }
}

// 优化的读取行函数
fn read_line<R: BufRead>(reader: &mut R, buf: &mut String) -> usize {
    buf.clear();
    match reader.read_line(buf) {
        Ok(0) => 0,
        Ok(n) => n,
        Err(_) => 0,
    }
}

// 分解的处理函数

fn process_ref_t(line: &str, temp: &mut f64) {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() > 1 {
        *temp = parts[1].parse().unwrap_or(0.0);
    }
}

fn process_topology<R: BufRead>(
    reader: &mut R, 
    buf: &mut String, 
    name: &mut String, 
    atoms_num: &mut usize, 
    molecule_types_num: &mut usize, 
    molecule_types: &mut Vec<MolType>
) {
    // name
    read_line(reader, buf);
    if let Some(caps) = RE_NAME.captures(&buf) {
        *name = caps[1].trim().to_string().replace(" ", "_");
    }
    println!("System name: {}", name);

    // atom num
    read_line(reader, buf);
    if let Some(caps) = RE_ATOMS_NUM.captures(&buf) {
        *atoms_num = caps[1].trim().parse().unwrap();
    }
    println!("Total atoms number: {}", atoms_num);

    // molecule types num
    read_line(reader, buf);
    if let Some(caps) = RE_MOLBLOCK.captures(&buf) {
        *molecule_types_num = caps[1].trim().parse().unwrap();
    }

    println!("System molecular types:");
    for mt_id in 0..*molecule_types_num {
        loop {
            read_line(reader, buf);
            if buf.trim().starts_with("molblock (") {
                read_line(reader, buf);
                if let Some(caps) = RE_MOLTYPE.captures(&buf) {
                    let mol_name = caps[1].to_string();
                    read_line(reader, buf);
                    if let Some(caps) = RE_MOLECULES_NUM.captures(&buf) {
                        let molecules_num: i64 = caps[1].parse().unwrap();
                        let moltype = MolType::new(mt_id, mol_name, molecules_num);
                        println!("{}", moltype);
                        molecule_types.push(moltype);
                        break;
                    }
                }
            }
        }
    }
}

fn process_ffparams<R: BufRead>(
    reader: &mut R, 
    buf: &mut String, 
    atom_types_num: &mut usize, 
    fun_type: &mut Vec<LJType>, 
    sigma: &mut Vec<f64>, 
    epsilon: &mut Vec<f64>, 
    radius: &mut Vec<f64>, 
    settings: &Settings
) {
    read_line(reader, buf);
    if let Some(caps) = RE_ATNR.captures(&buf) {
        *atom_types_num = caps[1].trim().parse().unwrap();
    }
    println!("Total atom types: {}.", atom_types_num);

    read_line(reader, buf);
    
    // 预分配向量容量
    fun_type.reserve(*atom_types_num * *atom_types_num);
    sigma.reserve(*atom_types_num);
    epsilon.reserve(*atom_types_num);
    radius.reserve(*atom_types_num);

    const INV_SIX: f64 = 1.0 / 6.0;

    for i in 0..*atom_types_num {
        for j in 0..*atom_types_num {
            read_line(reader, buf);
            if let Some(caps) = RE_FUNCTYPE.captures(&buf) {
                let c6: f64 = caps[2].parse().unwrap();
                let c12: f64 = caps[3].parse().unwrap();
                fun_type.push(LJType::new(c6, c12));
                
                if j == i {
                    if c6 != 0.0 && c12 != 0.0 {
                        let ratio = c12 / c6;
                        sigma.push(10.0 * ratio.powf(INV_SIX));
                        epsilon.push(c6 * c6 / (4.0 * c12));
                        radius.push(sigma.last().unwrap() / 2.0);
                    } else {
                        sigma.push(0.0);
                        epsilon.push(0.0);
                        radius.push(settings.radius_ff_default);
                    }
                }
            }
        }
    }
    println!("Total LJ function types: {}", fun_type.len());
}

fn process_moltype<R: BufRead>(
    reader: &mut R, 
    buf: &mut String, 
    molecules: &mut Vec<Molecule>, 
    atom_resids: &mut Vec<usize>, 
    atom_types: &mut Vec<usize>, 
    atom_radii: &mut Vec<f64>, 
    atom_charges: &mut Vec<f64>, 
    atom_names: &mut Vec<String>, 
    type_names: &mut Vec<String>, 
    radius: &[f64]
) {
    let mut atoms: Vec<Atom> = Vec::new();
    let mut residues: Vec<Residue> = Vec::new();
    let offset: usize = molecules.iter().map(|p| p.atoms_num).sum();

    if let Some(caps) = RE_MOLTYPE_ID.captures(&buf) {
        let molecule_type_id: usize = caps[1].parse().unwrap();
        println!("Reading molecule {} information...", molecule_type_id);
        
        read_line(reader, buf);
        let molecule_name = if let Some(caps) = RE_NAME.captures(&buf) {
            caps[1].to_string()
        } else {
            String::new()
        };

        read_line(reader, buf);
        read_line(reader, buf);
        
        let atoms_num: usize = if let Some(caps) = RE_ATOM_COUNT.captures(&buf) {
            caps[1].parse().unwrap()
        } else {
            0
        };

        // 预分配容量
        atom_resids.reserve(atom_resids.len() + atoms_num);
        atom_types.reserve(atom_types.len() + atoms_num);
        atom_radii.reserve(atom_radii.len() + atoms_num);
        atom_charges.reserve(atom_charges.len() + atoms_num);
        atom_names.reserve(atom_names.len() + atoms_num);
        type_names.reserve(type_names.len() + atoms_num);

        // 原子参数
        for _ in 0..atoms_num {
            read_line(reader, buf);
            if let Some(caps) = RE_ATOM_PARAMS.captures(&buf) {
                let atom_type_id: usize = caps[1].parse().unwrap();
                let atom_charge: f64 = caps[2].parse().unwrap();
                let residue_index: usize = caps[3].parse().unwrap();
                
                atom_resids.push(residue_index);
                atom_types.push(atom_type_id);
                atom_radii.push(radius[atom_type_id]);
                atom_charges.push(atom_charge);
            }
        }

        // 原子名称
        read_line(reader, buf);
        for _ in 0..atoms_num {
            read_line(reader, buf);
            if let Some(caps) = RE_ATOM_NAME.captures(&buf) {
                atom_names.push(caps[1].to_string());
            }
        }

        // 原子类型名称
        read_line(reader, buf);
        for _ in 0..atoms_num {
            read_line(reader, buf);
            if let Some(caps) = RE_TYPE_NAME.captures(&buf) {
                type_names.push(caps[1].to_string());
            }
        }

        // 残基
        read_line(reader, buf);
        let res_num: i32 = if let Some(caps) = RE_RESIDUE_COUNT.captures(&buf) {
            caps[1].parse().unwrap()
        } else {
            0
        };

        for _ in 0..res_num {
            read_line(reader, buf);
            if let Some(caps) = RE_RESIDUE.captures(&buf) {
                let id: usize = caps[1].parse().unwrap();
                let name = caps[2].to_string();
                let nr: i32 = caps[3].parse().unwrap();
                residues.push(Residue::new(id, name, nr));
            }
        }

        // 角度处理
        loop {
            read_line(reader, buf);
            if buf.trim().starts_with("Angle:") {
                read_line(reader, buf);
                if let Some(caps) = RE_ANGLES_NUM.captures(&buf) {
                    let angles: i32 = caps[1].parse().unwrap();
                    if angles > 0 {
                        read_line(reader, buf);
                        for _ in 0..angles / 4 {
                            read_line(reader, buf);
                            if let Some(caps) = RE_ANGLES.captures(&buf) {
                                let i: usize = caps[1].parse().unwrap();
                                let j: usize = caps[2].parse().unwrap();
                                let k: usize = caps[3].trim().parse().unwrap();
                                
                                if atom_names[offset + i].starts_with(['H', 'h']) {
                                    let base_name = &atom_names[offset + j];
                                    atom_names[offset + i] = format!("H{}", base_name);
                                }
                                if atom_names[offset + k].starts_with(['H', 'h']) {
                                    let base_name = &atom_names[offset + j];
                                    atom_names[offset + k] = format!("H{}", base_name);
                                }
                            } else {
                                break;
                            }
                        }
                        break;
                    } else { break; }
                }
            }
        }

        // 创建原子
        for id in 0..atoms_num {
            let global_id = id + offset;
            atoms.push(Atom::new(
                global_id,
                &type_names[global_id],
                atom_types[global_id],
                atom_charges[global_id],
                atom_resids[global_id],
                atom_names[global_id].clone(),
                atom_radii[global_id]
            ));
        }

        molecules.push(Molecule::new(molecule_type_id, molecule_name, atoms_num, &atoms, &residues));
    }
}

fn process_coordinates<R: BufRead>(
    reader: &mut R, 
    buf: &mut String, 
    atoms_num: usize, 
    coordinates: &mut Vec<f64>
) {
    println!("Reading coordinate information...");
    
    // 预分配容量
    coordinates.reserve(coordinates.len() + atoms_num * 3);
    
    for _ in 0..atoms_num {
        read_line(reader, buf);
        if let Some(caps) = RE_COORD.captures(&buf) {
            let x: f64 = caps[1].parse().unwrap();
            let y: f64 = caps[2].parse().unwrap();
            let z: f64 = caps[3].parse().unwrap();
            coordinates.push(x * 10.0);
            coordinates.push(y * 10.0);
            coordinates.push(z * 10.0);
        }
    }
}