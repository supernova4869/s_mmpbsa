use std::fmt::Formatter;
use std::path::Path;
use std::{fmt, fs};
use std::io::Write;
use std::io::BufReader;
use ndarray::Array2;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufRead;

use crate::settings::Settings;

pub struct TPR {
    pub name: String,
    pub n_atoms: usize,
    pub molecule_types_num: usize,
    pub molecule_types: Vec<MolType>,
    pub atom_types_num: usize,
    pub lj_sr_params: Vec<LJType>,
    pub molecules: Vec<Molecule>,
    pub dt: f64,
    pub nsteps: u64,
    pub nstxout: u32,
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
    pub fn new(name: &str, n_atoms: usize, molecule_types_num: usize, molecule_types: Vec<MolType>, 
        atom_types_num: usize, lj_sr_params: Vec<LJType>, molecules: Vec<Molecule>, dt: f64, nsteps: u64,
        nstxout: u32, temp: f64, coordinates: Array2<f64>) -> TPR {
            TPR {
                name: name.to_string(),
                n_atoms,
                molecule_types_num,
                molecule_types,
                atom_types_num,
                lj_sr_params,
                molecules,
                dt,
                nsteps,
                nstxout,
                temp,
                coordinates,
            }
    }

    pub fn from(mdp: &str, settings: &Settings) -> TPR {
        let mut name = String::new();
        let mut atoms_num = 0;
        let mut molecule_types_num = 0;
        let mut atom_types_num = 0;
        let mut molecule_types: Vec<MolType> = vec![];

        let file = File::open(mdp).unwrap();
        let mut reader = BufReader::new(file);
        let mut buf = String::from("");

        let mut fun_type: Vec<LJType> = vec![];
        let mut sigma: Vec<f64> = vec![];
        let mut epsilon: Vec<f64> = vec![];
        let mut radius: Vec<f64> = vec![];

        let mut atom_resids: Vec<usize> = vec![];   // residue ids of each atom
        let mut atom_types: Vec<usize> = vec![];    // atom type
        let mut atom_radii: Vec<f64> = vec![];      // atom radius
        let mut atom_charges: Vec<f64> = vec![];    // atom charge
        let mut atom_names: Vec<String> = vec![];   // atom name
        let mut type_names: Vec<String> = vec![];   // atom type name
        
        let mut coordinates: Vec<f64> = vec![];   // atom coordinates

        let mut molecules: Vec<Molecule> = vec![];

        // simulation time parameters
        let mut dt = 0.0;
        let mut nsteps = 0;
        let mut nstxout = 0;
        let mut temp = 0.0;

        println!("Loading dump file: {}\n", mdp);
        loop {
            let bytes = read_line(&mut reader, &mut buf);
            if bytes == 0 {
                break;
            }

            // MD steps
            if buf.starts_with("inputrec:") {
                loop {
                    read_line(&mut reader, &mut buf);
                    if buf.trim().starts_with("dt") {
                        let re = Regex::new(r"dt\s+=\s*(.*)").unwrap();
                        dt = re.captures(&buf).unwrap().get(1).unwrap().as_str().trim().parse().unwrap();
                        read_line(&mut reader, &mut buf);
                        let re = Regex::new(r"nsteps\s+=\s*(.*)").unwrap();
                        nsteps = re.captures(&buf).unwrap().get(1).unwrap().as_str().trim().parse().unwrap();
                    } else if buf.trim().starts_with("nstxout-compressed") {
                        let re = Regex::new(r"nstxout-compressed\s+=\s*(.*)").unwrap();
                        nstxout = re.captures(&buf).unwrap().get(1).unwrap().as_str().trim().parse().unwrap();
                        break
                    }
                }
            }

            if buf.trim().starts_with("ref-t:") {
                let ref_t: Vec<&str> = buf.split(" ").filter_map(|p| match p.trim().len() {
                    0 => None,
                    _ => Some(p)
                }).collect();
                temp = ref_t[1].trim().parse().expect("Failed to get temperature.");
            }

            // molecules define
            if buf.starts_with("topology:") {
                // name
                read_line(&mut reader, &mut buf);       // name="Protein in water"
                let re = Regex::new("name\\s*=\\s*\"(.*)\"").unwrap();
                name = re.captures(&buf).unwrap().get(1).unwrap().as_str().trim().to_string();
                name = name.replace(" ", "_");
                println!("System name: {}", name);

                // atom num
                read_line(&mut reader, &mut buf);       // #atoms = 3218
                let re = Regex::new(r"#atoms\s*=\s*(\d+)").unwrap();
                atoms_num = re.captures(&buf).unwrap().get(1).unwrap().as_str().trim().parse().unwrap();
                println!("Total atoms number: {}", atoms_num);

                // molecule types num
                read_line(&mut reader, &mut buf);
                let re = Regex::new(r"#molblock\s*=\s*(\d+)").unwrap();
                molecule_types_num = re.captures(&buf).unwrap().get(1).unwrap().as_str().trim().parse().unwrap();

                println!("System molecular types:");
                for mt_id in 0..molecule_types_num {
                    loop {
                        read_line(&mut reader, &mut buf);
                        if buf.trim().starts_with("molblock (") {
                            read_line(&mut reader, &mut buf);
                            let re = Regex::new("moltype\\s*=\\s*\\d+\\s*\"(.*)\"").unwrap();
                            let name = re.captures(&buf).unwrap().get(1).unwrap().as_str().to_string();
                            read_line(&mut reader, &mut buf);
                            let re = Regex::new(r"#molecules\s*=\s*(\d+)").unwrap();
                            let molecules_num: i64 = re.captures(&buf).unwrap().get(1).unwrap().as_str().parse().unwrap();
                            let moltype = MolType::new(mt_id, name, molecules_num);
                            println!("{}", moltype);
                            molecule_types.push(moltype);
                            break;
                        }
                    }
                }
            }

            // force field parameters (atom radius here)
            if buf.trim().starts_with("ffparams:") {
                read_line(&mut reader, &mut buf);
                let re = Regex::new(r"atnr\s*=\s*(\d+)").unwrap();
                atom_types_num = re.captures(&buf).unwrap().get(1).unwrap().as_str().parse().unwrap();
                println!("Total atom types: {}.", atom_types_num);

                read_line(&mut reader, &mut buf);
                // functype[0]=LJ_SR, c6= 2.07413384e-03, c12= 1.51207360e-06
                let re = Regex::new(r"functype\[(\d*)]=LJ_SR,\s*c6\s*=\s*(\S*),\s*c12\s*=\s*(\S*)").unwrap();
                for i in 0..atom_types_num {
                    for j in 0..atom_types_num {
                        read_line(&mut reader, &mut buf);
                        let m = re.captures(&buf).unwrap();
                        let c6: f64 = m.get(2).unwrap().as_str().parse().unwrap();
                        let c12: f64 = m.get(3).unwrap().as_str().parse().unwrap();
                        fun_type.push(LJType::new(c6, c12));
                        // calculate σ, ε, radius for each atom
                        if j == i {
                            if c6 != 0.0 && c12 != 0.0 {
                                sigma.push(10.0 * (c12 / c6).powf(1.0 / 6.0)); // nm to A
                                epsilon.push(c6.powi(2) / (4.0 * c12));
                                radius.push(sigma[i] / 2.0); // sigma is diameter
                            } else {
                                sigma.push(0.0);
                                epsilon.push(0.0);
                                radius.push(settings.radius_ff_default);
                            }
                        }
                    }
                }
                println!("Total LJ function types: {}", fun_type.len());
            }

            if buf.trim().starts_with("moltype (") {
                let mut atoms: Vec<Atom> = vec![];
                let mut residues: Vec<Residue> = vec![];
                // currently exist atoms, to avoid duplicate process
                let offset: usize = molecules.iter().map(|p| p.atoms_num).sum();

                let re = Regex::new(r"moltype \((\d+)\)").unwrap();
                let molecule_type_id: usize = re.captures(&buf).unwrap().get(1).unwrap().as_str().parse().unwrap();
                println!("Reading molecule {} information...", molecule_type_id);
                read_line(&mut reader, &mut buf);
                let re = Regex::new("name\\s*=\\s*\"(.*)\"").unwrap();
                let molecule_name = re.captures(&buf).unwrap().get(1).unwrap().as_str().to_string();
                read_line(&mut reader, &mut buf);
                read_line(&mut reader, &mut buf);
                let re = Regex::new(r"atom \((\d+)\):").unwrap();
                let atoms_num: usize = re.captures(&buf).unwrap().get(1).unwrap().as_str().parse().unwrap();

                // atom parameters
                // atom[     0]={type=  0, typeB=  0, ptype=    Atom, m= 1.60000e+01,
                // q=-4.91104e-01, mB= 1.60000e+01, qB=-4.91104e-01, resind=    0, atomnumber= -1}
                let re = Regex::new(r".*type=\s*(\d+).*q=\s*([^,]+),.*resind=\s*(\d+).*").unwrap();
                for _ in 0..atoms_num {
                    read_line(&mut reader, &mut buf);
                    let c = re.captures(&buf).unwrap();
                    let atom_type_id: usize = c.get(1).unwrap().as_str().parse().unwrap();
                    let atom_charge: f64 = c.get(2).unwrap().as_str().parse().unwrap();
                    let residue_index: usize = c.get(3).unwrap().as_str().parse().unwrap();
                    atom_resids.push(residue_index);
                    atom_types.push(atom_type_id);
                    atom_radii.push(radius[atom_type_id]);
                    atom_charges.push(atom_charge);
                }

                // atom names
                read_line(&mut reader, &mut buf);
                // atom[0]={name="O1"}
                let re = Regex::new("name=\"(.*)\"").unwrap();
                for _ in 0..atoms_num {
                    read_line(&mut reader, &mut buf);
                    let name = re.captures(&buf).unwrap().get(1).unwrap().as_str();
                    atom_names.push(name.to_string());
                }

                // atom types
                read_line(&mut reader, &mut buf);
                // type[0]={name="N3",nameB="N3"}
                let re = Regex::new("name=\"(.*)\",").unwrap();
                for _ in 0..atoms_num {
                    read_line(&mut reader, &mut buf);
                    let name = re.captures(&buf).unwrap().get(1).unwrap().as_str();
                    type_names.push(name.to_string());
                }

                // residues
                read_line(&mut reader, &mut buf);
                let re = Regex::new(r"\s*residue \((\d+)\)").unwrap();
                let res_num: i32 = re.captures(&buf).unwrap().get(1).unwrap().as_str().parse().unwrap();
                let re = Regex::new("residue\\[(\\d+)]=\\{name=\"(.+)\",.*nr=([\\d\\-]+).*").unwrap();
                for _ in 0..res_num {
                    read_line(&mut reader, &mut buf);
                    let m = re.captures(&buf).unwrap();
                    let id: usize = m.get(1).unwrap().as_str().parse().unwrap();
                    let name = m.get(2).unwrap().as_str().to_string();
                    let nr: i32 = m.get(3).unwrap().as_str().parse().unwrap();
                    residues.push(Residue::new(id, name, nr));
                }

                // angles
                loop {
                    read_line(&mut reader, &mut buf);
                    if buf.trim().starts_with("Angle:") {
                        read_line(&mut reader, &mut buf);
                        let re = Regex::new(r"nr\s*:\s*(\d+)").unwrap();
                        let angles: i32 = re.captures(&buf).unwrap().get(1).unwrap().as_str().parse().unwrap();
                        if angles > 0 {
                            read_line(&mut reader, &mut buf);
                            let re = Regex::new(r"\(ANGLES\)\s+(\d+)\s+(\d+)\s+(\d+)").unwrap();
                            // assign H types by connection atoms from angle information
                            for _ in 0..angles / 4 {
                                read_line(&mut reader, &mut buf);
                                if re.is_match(&buf) {
                                    let c = re.captures(&buf).unwrap();
                                    let i: usize = c.get(1).unwrap().as_str().parse().unwrap();
                                    let j: usize = c.get(2).unwrap().as_str().parse().unwrap();
                                    let k: usize = c.get(3).unwrap().as_str().trim().parse().unwrap();
                                    if atom_names[offset + i].starts_with(['H', 'h']) {
                                        atom_names[offset + i] = format!("H{}", atom_names[offset + j]);
                                    }
                                    if atom_names[offset + k].starts_with(['H', 'h']) {
                                        atom_names[offset + k] = format!("H{}", atom_names[offset + j]);
                                    }
                                } else {
                                    break;
                                }
                            }
                            break;
                        } else { break; }
                    }
                }

                for id in 0..atoms_num {
                    let id = id + offset;
                    atoms.push(Atom::new(id,
                                        &type_names[id],
                                        atom_types[id],
                                        atom_charges[id],
                                        atom_resids[id],
                                        atom_names[id].to_string(),
                                        atom_radii[id]));
                }

                molecules.push(Molecule::new(molecule_type_id, molecule_name, atoms_num,
                                             &atoms, &residues));
            }

            // coordinates
            // x (3218x3):
            if buf.trim().starts_with("x (") {
                println!("Reading coordinate information...");
                let re = Regex::new(r"\{\s*(.*),\s*(.*),\s*(.*)\}").unwrap();
                for _ in 0..atoms_num {
                    read_line(&mut reader, &mut buf);
                    // x[    0]={ 1.41430e+00,  1.38595e+00,  2.14591e-01}
                    let caps = re.captures(&buf).unwrap();
                    let x: f64 = caps.get(1).unwrap().as_str().parse().unwrap();
                    let y: f64 = caps.get(2).unwrap().as_str().parse().unwrap();
                    let z: f64 = caps.get(3).unwrap().as_str().parse().unwrap();
                    coordinates.push(x * 10.0);
                    coordinates.push(y * 10.0);
                    coordinates.push(z * 10.0);
                }
            }
        }

        println!("Backup force field radius...");
        let ff_dat = Path::new(mdp).parent().unwrap().join("ff_radius.dat");
        if ff_dat.is_file() {
            fs::remove_file(&ff_dat).unwrap();
        }
        let mut ff_dat = File::create(ff_dat).unwrap();
        for r in atom_radii {
            writeln!(ff_dat, "{:.2}", r).unwrap();
        }

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
            dt,
            nsteps,
            nstxout,
            temp,
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

fn read_line(reader: &mut BufReader<File>, buf: &mut String) -> usize {
    buf.clear();
    reader.read_line(buf).unwrap()
}
