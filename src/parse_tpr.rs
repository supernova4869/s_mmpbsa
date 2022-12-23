use std::fmt::Formatter;
use std::fmt;
use std::io::BufReader;
use regex::Regex;
use std::fs::File;
use std::io::BufRead;

pub struct TPR {
    pub name: String,
    pub atoms_num: usize,
    pub molecule_types_num: usize,
    pub molecule_types: Vec<MolType>,
    pub atom_types_num: usize,
    pub lj_sr_params: Vec<LJType>,
    pub molecules: Vec<Molecule>,
    pub dt: f64,
    pub nsteps: u64,
}

impl fmt::Display for TPR {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}, with\n{} atoms,\n{} type(s) of molecules,\n{} atom types,\n{} LJ types",
               self.name, self.atoms_num, self.molecule_types_num,
               self.atom_types_num, self.lj_sr_params.len()
        )
    }
}

impl TPR {
    pub fn new(mdp: &str) -> TPR {
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
        let mut atom_sigmas: Vec<f64> = vec![];     // atom sigma
        let mut atom_epsilons: Vec<f64> = vec![];   // atom epsilon
        let mut atom_charges: Vec<f64> = vec![];    // atom charge
        let mut atom_names: Vec<String> = vec![];   // atom name

        let mut molecules: Vec<Molecule> = vec![];

        // simulation time parameters
        let mut dt = 0.0;
        let mut nsteps = 0;

        println!("Loading dumped tpr file: {}\n", mdp);
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
                        break;
                    }
                }
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

            // force field parameters
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
                        let func_id: usize = m.get(1).unwrap().as_str().parse().unwrap();
                        let c6: f64 = m.get(2).unwrap().as_str().parse().unwrap();
                        let c12: f64 = m.get(3).unwrap().as_str().parse().unwrap();
                        fun_type.push(LJType::new(func_id, c6, c12));
                        // calculate σ, ε, radius for each atom
                        if j == i {
                            if c6 != 0.0 && c12 != 0.0 {
                                sigma.push(10.0 * (c12 / c6).powf(1.0 / 6.0)); // nm to A
                                epsilon.push(c6.powi(2) / (4.0 * c12));
                                radius.push(sigma[i] / 2.0); // sigma is diameter
                            } else {
                                sigma.push(0.0); // nm to A
                                epsilon.push(0.0);
                                radius.push(0.0); // sigma is diameter
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
                    atom_sigmas.push(sigma[atom_type_id]);
                    atom_epsilons.push(epsilon[atom_type_id]);
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

                loop {
                    read_line(&mut reader, &mut buf);
                    if buf.trim().starts_with("residue (") {
                        // residues
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
                    }
                    if buf.trim().starts_with("excls:") {
                        break;
                    }
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
                                         atom_types[id],
                                         atom_charges[id],
                                         atom_resids[id],
                                         atom_names[id].to_string(),
                                         atom_sigmas[id],
                                         atom_epsilons[id],
                                         atom_radii[id]));
                }

                molecules.push(Molecule::new(molecule_type_id, molecule_name, atoms_num,
                                             &atoms, &residues));
            }
        }
        println!("System molecular composition:");
        for mol in &molecules {
            println!("Molecule {}: {}", mol.molecule_type_id, mol);
        }
        TPR {
            name,
            atoms_num,
            molecule_types_num,
            molecule_types,
            atom_types_num,
            lj_sr_params: fun_type,
            molecules,
            dt,
            nsteps,
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
    pub func_id: usize,
    pub c6: f64,
    pub c12: f64,
}

impl LJType {
    fn new(func_id: usize, c6: f64, c12: f64) -> LJType {
        LJType {
            func_id,
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
    pub type_id: usize,
    pub charge: f64,
    pub residue_index: usize,
    pub name: String,
    pub sigma: f64,
    pub epsilon: f64,
    pub radius: f64,
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "Atom {}: {} with type {}, charge {}, radius {}, in residue {}",
               self.id, self.name, self.type_id, self.charge, self.radius, self.residue_index)
    }
}

impl Atom {
    fn new(id: usize, type_id: usize, charge: f64, residue_index: usize, name: String,
           sigma: f64, epsilon: f64, radius: f64) -> Atom {
        Atom {
            id,
            type_id,
            charge,
            residue_index,
            name,
            sigma,
            epsilon,
            radius,
        }
    }
}

#[derive(Clone)]
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
    fn new(id: usize, name: String, nr: i32) -> Residue {
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
