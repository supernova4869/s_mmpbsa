use std::collections::HashMap;
use std::fmt::Debug;
use std::str::FromStr;
use std::fs;
use std::io::Write;
use indicatif::ProgressBar;
use ndarray::{Array1, Array2};
use regex::Regex;
use crate::index_parser::Index;

pub fn fetch_from_tpr(mdp: &String, ndx: &Index, receptor_grp: usize, ligand_grp: usize,
                      rad_type: i32, rad_lj0: f64)
                      -> (Array2<usize>, Array2<usize>, Array2<f64>, Array2<f64>, Array2<f64>, Array2<f64>) {
    // read mdp file
    // let mut qrv = fs::File::create(qrv).expect("Create qrv file failed");
    // qrv.write_all("Protein Ligand\n".as_bytes()).expect("Writing qrv file failed");
    let ndx_rec = &ndx.groups[receptor_grp].indexes;
    let ndx_lig = &ndx.groups[ligand_grp].indexes;
    let mdp = fs::read_to_string(mdp).unwrap();
    let mdp: Vec<&str> = mdp.split("\n").collect();

    // get MD parameters
    let re = Regex::new(r"\s*ffparams:").unwrap();
    let locator = get_md_locators_first(&mdp, &re).unwrap();
    // number of atom types
    let re = Regex::new(r"\s*atnr=(\d+)").unwrap();
    let atnr = re.captures(mdp[locator + 1]).unwrap();
    let atnr = atnr.get(1).unwrap().as_str();
    // qrv.write_all(format!("{}\n", atnr).as_bytes()).expect("Writing qrv file failed");
    let atnr: usize = atnr.parse().unwrap();
    // LJ parameters
    let locator = locator + 3;
    let mut sigma: Array1<f64> = Array1::zeros(atnr);
    let mut epsilon: Array1<f64> = Array1::zeros(atnr);
    let mut rad: Array1<f64> = Array1::ones(atnr) * rad_lj0;

    // println!("Generating qrv file...");
    println!("Reading atom L-J parameters..");
    let pb = ProgressBar::new(atnr as u64);
    for i in 0..atnr {
        // qrv.write_all(format!("{:6}", i).as_bytes()).expect("Writing qrv file failed");
        // get c6 and c12 parameters for each atom
        for j in 0..atnr {
            let re = Regex::new(r".*c6\s*=\s*(.*),.*c12\s*=\s*(.*)").unwrap();
            let m = re.captures(&mdp[locator + i * atnr + j]).unwrap();
            let c6 = m.get(1).unwrap().as_str();
            let c12 = m.get(2).unwrap().as_str().trim();
            // qrv.write_all(format!(" {} {}", c6, c12).as_bytes()).expect("Writing qrv file failed");
            let c6: f64 = c6.parse().unwrap();
            let c12: f64 = c12.parse().unwrap();
            // calculate σ, ε, radius for each atom
            if j == i && c6 != 0.0 && c12 != 0.0 {
                sigma[i] = 10.0 * (c12 / c6).powf(1.0 / 6.0); // 转换单位为A
                epsilon[i] = c6.powi(2) / (4.0 * c12);
                rad[i] = 0.5 * sigma[i]; // sigma为直径
            }
        }
        // qrv.write_all("\n".as_bytes()).expect("Writing qrv file failed");
        pb.inc(1);
    }

    // number of groups
    let re = Regex::new(r"\s*#molblock\s*=\s*(.+?)\s*").unwrap();
    let grp_num: usize = get_md_params_first(&mdp, &re).parse().unwrap();
    // number of molecules
    let re = Regex::new(r"\s*#molecules\s*=\s*(.+?)\s*").unwrap();
    let mol_num: Vec<usize> = get_md_params_all(&mdp, &re);
    // number of atoms
    let re = Regex::new(r"atom \((.+)\):").unwrap();
    let max_atm_num: usize = *get_md_params_all(&mdp, &re).iter().max().unwrap();
    let mut grp_names: Vec<String> = vec![];                // name of each group
    let mut grp_atom_nums: Vec<usize> = vec![];             // atom number of each group

    // locator of molecule types
    let re = Regex::new(r"\s*moltype.+\(").unwrap();
    let locators = get_md_locators_all(&mdp, &re).unwrap();

    // initialize atom information, each molecule per line
    let mut res_ids = Array2::<usize>::zeros((mol_num.len(), max_atm_num));      // residue ids of each atom
    // atom type
    let mut c_atoms = Array2::<usize>::zeros((mol_num.len(), max_atm_num));     // atom type
    // atom radius
    let mut r_atoms = Array2::<f64>::zeros((mol_num.len(), max_atm_num));        // atom radius
    // atom sigma
    let mut s_atoms = Array2::<f64>::zeros((mol_num.len(), max_atm_num));        // atom sigma
    // atom epsilon
    let mut e_atoms = Array2::<f64>::zeros((mol_num.len(), max_atm_num));        // atom epsilon
    // atom charge
    let mut q_atoms = Array2::<f64>::zeros((mol_num.len(), max_atm_num));        // atom charge
    // atom name
    let mut t_atoms = Array2::<String>::default((mol_num.len(), max_atm_num)); // atom name

    // get atom parameters
    for locator in locators {
        let re = Regex::new(r"\s*moltype.+\((\d+)\)").unwrap();
        let mol_id = re.captures(&mdp[locator]).unwrap().get(1).unwrap();
        let mol_id: usize = mol_id.as_str().parse().unwrap();
        let re = Regex::new(r"name=(.*)").unwrap();
        let n = re.captures(&mdp[locator + 1]).unwrap().get(1).unwrap().as_str().trim();
        grp_names.push(n[1..(n.len() - 1)].to_string());
        let re = Regex::new(r"\((.*)\)").unwrap();
        let num = re.captures(&mdp[locator + 3]).unwrap().get(1).unwrap();
        let num: usize = num.as_str().parse().unwrap();
        grp_atom_nums.push(num);
        let locator = locator + 4;

        println!("Reading the {}/{} system atoms information.", mol_id + 1, mol_num.len());
        println!("Reading atom property parameters...");
        let pb = ProgressBar::new(grp_atom_nums[mol_id] as u64);   // progress bar
        for i in 0..grp_atom_nums[mol_id] {
            let re = Regex::new(r".*type=\s*(\d+).*q=\s*([^,]+),.*resind=\s*(\d+).*").unwrap();
            let c = re.captures(&mdp[locator + i]).unwrap();
            let at_type = c.get(1).unwrap();
            let at_type: usize = at_type.as_str().parse().unwrap();
            let res_id = c.get(3).unwrap();
            let res_id: usize = res_id.as_str().parse().unwrap();
            res_ids[[mol_id, i]] = res_id;
            c_atoms[[mol_id, i]] = at_type;
            r_atoms[[mol_id, i]] = rad[at_type];
            s_atoms[[mol_id, i]] = sigma[at_type];
            e_atoms[[mol_id, i]] = epsilon[at_type];
            let q = c.get(2).unwrap();
            let q: f64 = q.as_str().parse().unwrap();
            q_atoms[[mol_id, i]] = q;
            pb.inc(1);
        }
        let locator = locator + grp_atom_nums[mol_id] + 1;     // 不加1是"atom (3218):"行

        // get atom names
        println!("Reading atom names...");
        let pb = ProgressBar::new(grp_atom_nums[mol_id] as u64);
        for i in 0..grp_atom_nums[mol_id] {
            let re = Regex::new("name=\"(.*)\"").unwrap();
            let name = re.captures(&mdp[locator + i]).unwrap();
            let name = name.get(1).unwrap().as_str();
            t_atoms[[mol_id, i]] = name.to_string();
            pb.inc(1);
        }
    }

    // get residues information
    let re = Regex::new(r"\s*residue \((\d+)\)").unwrap();
    let locators = get_md_locators_all(&mdp, &re).unwrap();
    let mut resnums: Vec<usize> = vec![];
    for i in 0..locators.len() {
        let res_num = re.captures(&mdp[locators[i]]).unwrap().get(1).unwrap();
        let res_num: usize = res_num.as_str().trim().parse().unwrap();
        resnums.push(res_num);
    }
    let max_res_num: usize = *resnums.iter().max().unwrap() as usize;
    let mut res_names = Array2::<String>::default((mol_num.len(), max_res_num));

    for (idx, locator) in locators.into_iter().enumerate() {
        let re = Regex::new(".*name=\"(.+)\",.*nr=(\\d+).*").unwrap();
        println!("Reading residues information...");
        let pb = ProgressBar::new(resnums[idx] as u64);
        for i in 0..resnums[idx] {
            let m = re.captures(&mdp[locator + 1 + i]).unwrap();
            let name = m.get(1).unwrap().as_str();
            let nr = m.get(2).unwrap().as_str();
            let nr: i32 = nr.parse().unwrap();
            res_names[[idx, i]] = format!("{:05}{}", nr, name);
            pb.inc(1);
        }
    }

    // assign H types by connection atoms from angle information
    let re = Regex::new(r"^ +Angle:").unwrap();
    let locators = get_md_locators_all(&mdp, &re).unwrap();
    println!("Reading angles...");
    for (mol_id, locator) in locators.into_iter().enumerate() {
        let angles_num: Vec<&str> = (&mdp[locator + 1]).trim().split(" ").collect();
        let angles_num: usize = angles_num[1].parse().unwrap();
        let angles_num = angles_num / 4;
        if angles_num > 0 {
            let re = Regex::new(r"\d+ type=\d+ \(ANGLES\)\s+(\d+)\s+(\d+)\s+(\d+)").unwrap();
            let pb = ProgressBar::new(angles_num as u64);
            for l_num in locator + 3..locator + 3 + angles_num {
                if re.is_match(&mdp[l_num]) {
                    let paras = re.captures(&mdp[l_num]).unwrap();
                    let i = paras.get(1).unwrap();
                    let i: usize = i.as_str().parse().unwrap();
                    let j = paras.get(2).unwrap();
                    let j: usize = j.as_str().parse().unwrap();
                    let k = paras.get(3).unwrap();
                    let k: usize = k.as_str().trim().parse().unwrap();
                    if t_atoms[[mol_id, i]].starts_with(['H', 'h']) {
                        t_atoms[[mol_id, i]] = format!("H{}", t_atoms[[mol_id, j]]);
                    }
                    if t_atoms[[mol_id, k]].starts_with(['H', 'h']) {
                        t_atoms[[mol_id, k]] = format!("H{}", t_atoms[[mol_id, j]]);
                    }
                } else { break; }
                pb.inc(1);
            }
        }
    }

    // fix atom radius
    let mut atom_id = 0;
    for i in 0..grp_num {
        for n in 0..mol_num[i] {
            let pb = ProgressBar::new(grp_atom_nums[i] as u64);
            for j in 0..grp_atom_nums[i] {
                if ndx_rec.contains(&atom_id) || ndx_lig.contains(&atom_id) {
                    match rad_type {
                        1 => {
                            let re = Regex::new(r"([a-zA-Z]+)\d*").unwrap();
                            let res = re.captures(t_atoms[[i, j]].as_str()).unwrap();
                            let res = res.get(1).unwrap().as_str();
                            r_atoms[[n, j]] = get_radi(res);
                        }
                        _ => ()
                    }
                }
                atom_id += 1;
                pb.inc(1);
            }
        }
    }

    println!("Finished reading MD parameters.");
    return (res_ids, c_atoms, r_atoms, s_atoms, e_atoms, q_atoms);
}

fn get_md_locators_first(strings: &Vec<&str>, re: &Regex) -> Result<usize, usize> {
    for (idx, l) in strings.into_iter().enumerate() {
        if re.is_match(l) {
            return Ok(idx);
        }
    }
    return Err(0);
}

fn get_md_locators_all(strings: &Vec<&str>, re: &Regex) -> Result<Vec<usize>, Vec<usize>> {
    let mut locators: Vec<usize> = vec![];
    for (idx, l) in strings.into_iter().enumerate() {
        if re.is_match(l) {
            locators.push(idx);
        }
    }
    if let 0 = locators.len() {
        return Err(locators);
    } else {
        return Ok(locators);
    }
}

#[warn(dead_code)]
fn get_md_params_first(strings: &Vec<&str>, re: &Regex) -> String {
    for line in strings {
        if re.is_match(line) {
            let value = re.captures(line).unwrap();
            let value = value.get(1).unwrap().as_str();
            return value.to_string();
        }
    }
    return String::new();
}

fn get_md_params_all<T: Debug + FromStr>(strings: &Vec<&str>, re: &Regex) -> Vec<T> where <T as FromStr>::Err: Debug {
    let mut values: Vec<T> = vec![];
    for line in strings {
        if re.is_match(line) {
            let value = re.captures(line).unwrap();
            let value = value.get(1).unwrap().as_str();
            let value: T = value.trim().parse().unwrap();
            values.push(value);
        }
    }
    return values;
}

fn get_radi(at_type: &str) -> f64 { // mBondi from AMBER20/parmed/tools/changeradii.py
    let rad_bondi: HashMap<&str, f64> = vec![("C", 1.7),
                                             ("H", 1.2),
                                             ("N", 1.55),
                                             ("HC", 1.3),
                                             ("O", 1.5),
                                             ("HN", 1.3),
                                             ("F", 1.5),
                                             ("HP", 1.3),
                                             ("SI", 2.1),
                                             ("HO", 0.8),
                                             ("P", 1.85),
                                             ("HS", 0.8),
                                             ("S", 1.8),
                                             ("CL", 1.7),
                                             ("BR", 1.85),
                                             ("I", 1.98)].into_iter().collect();
    let at_type = at_type.to_uppercase();
    let mut radius = 1.5;
    if at_type.len() >= 2 {
        let r = rad_bondi.get(&at_type[0..2]);
        if let Some(m) = r {
            radius = *m;
        } else {
            let r = rad_bondi.get(&at_type[0..1]);
            if let Some(m) = r {
                radius = *m;
            }
        }
    } else {
        let r = rad_bondi.get(&at_type[0..1]);
        if let Some(m) = r {
            radius = *m;
        }
    }
    return radius;
}