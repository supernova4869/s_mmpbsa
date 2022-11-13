use std::collections::HashMap;
use std::fmt::Debug;
use std::str::FromStr;
use std::fs;
use std::io::Write;
use std::path::Path;
use indicatif::ProgressBar;
use ndarray::Array2;
use regex::Regex;
use crate::index_parser::Index;
use crate::mmpbsa::gen_file_sha256;
use crate::Parameters;

pub fn gen_qrv(mdp: &String, tpr:&String, ndx: &Index, wd: &Path,
               receptor_grp: usize, ligand_grp: usize,
               qrv: &Path, settings: &Parameters) {
    // read mdp file
    let qrv = qrv.to_str().unwrap();
    let rad_type = settings.rad_type;
    let rad_lj0 = settings.rad_lj0;

    let mut qrv_content = fs::File::create(qrv).expect("Create parameter qrv file failed");
    qrv_content.write_all("Protein Ligand\n".as_bytes()).expect("Writing parameter qrv file failed");
    let ndx_rec = &ndx.groups[receptor_grp].indexes;
    let ndx_lig = &ndx.groups[ligand_grp].indexes;
    let mdp_content = fs::read_to_string(mdp).unwrap();
    if mdp_content.len() == 0 {
        println!("Error with _mdout.mdp: file empty");
    }
    let mdp_content: Vec<&str> = mdp_content.split("\n").collect();

    // get MD parameters
    let re = Regex::new(r"\s*ffparams:").unwrap();
    let locator = get_md_locators_first(&mdp_content, &re);
    // number of atom types
    let re = Regex::new(r"\s*atnr=(\d+)").unwrap();
    let atnr = re.captures(mdp_content[locator + 1]).expect("Parse mdp file error.");
    let atnr = atnr.get(1).unwrap().as_str();
    qrv_content.write_all(format!("{}\n", atnr).as_bytes()).expect("Writing parameter qrv file failed");
    let atnr: usize = atnr.parse().unwrap();
    // LJ parameters
    let locator = locator + 3;
    let mut sigma: Vec<f64> = vec![0.0; atnr];
    let mut epsilon: Vec<f64> = vec![0.0; atnr];
    let mut rad: Vec<f64> = vec![rad_lj0; atnr];

    println!("Generating qrv file...");
    println!("Writing atom L-J parameters..");
    let pb = ProgressBar::new(atnr as u64);
    for i in 0..atnr {
        qrv_content.write_all(format!("{:6}", i).as_bytes()).expect("Writing qrv file failed");
        // get c6 and c12 parameters for each atom
        for j in 0..atnr {
            let re = Regex::new(r".*c6\s*=\s*(.*),.*c12\s*=\s*(.*)").unwrap();
            let m = re.captures(&mdp_content[locator + i * atnr + j]).unwrap();
            let c6 = m.get(1).unwrap().as_str();
            let c12 = m.get(2).unwrap().as_str().trim();
            qrv_content.write_all(format!(" {} {}", c6, c12).as_bytes()).expect("Writing qrv file failed");
            let c6: f64 = c6.parse().unwrap();
            let c12: f64 = c12.parse().unwrap();
            // calculate σ, ε, radius for each atom
            if j == i  && c6 != 0.0 && c12 != 0.0 {
                sigma[i] = 10.0 * (c12 / c6).powf(1.0 / 6.0); // 转换单位为A
                epsilon[i] = c6.powi(2) / (4.0 * c12);
                rad[i] = 0.5 * sigma[i]; // sigma为直径
            }
        }
        qrv_content.write_all("\n".as_bytes()).expect("Writing qrv file failed");
        pb.inc(1);
    }

    // number of groups
    let re = Regex::new(r"\s*#molblock\s*=\s*(.+?)\s*").unwrap();
    let grp_num: usize = get_md_params_first(&mdp_content, &re).parse().unwrap();
    // number of molecules
    let re = Regex::new(r"\s*#molecules\s*=\s*(.+?)\s*").unwrap();
    let mol_num: Vec<usize> = get_md_params_all(&mdp_content, &re);
    // number of atoms
    let re = Regex::new(r"atom \((.+)\):").unwrap();
    let max_atm_num: usize = *get_md_params_all(&mdp_content, &re).iter().max().unwrap();
    let mut sys_names: Vec<String> = vec![];                // name of each system
    let mut sys_atom_nums: Vec<usize> = vec![];             // atom number of each system

    // locator of molecule types
    let re = Regex::new(r"\s*moltype.+\(").unwrap();
    let locators = get_md_locators_all(&mdp_content, &re);

    // initialize atom information, each molecule per line
    // 考虑使用泛型值来优化
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
        let mol_id = re.captures(&mdp_content[locator]).unwrap().get(1).unwrap();
        let mol_id: usize = mol_id.as_str().parse().unwrap();
        let re = Regex::new(r"name=(.*)").unwrap();
        let n = re.captures(&mdp_content[locator + 1]).unwrap().get(1).unwrap().as_str().trim();
        sys_names.push(n[1..(n.len() - 1)].to_string());
        let re = Regex::new(r"\((.*)\)").unwrap();
        let num = re.captures(&mdp_content[locator + 3]).unwrap().get(1).unwrap();
        let num: usize = num.as_str().parse().unwrap();
        sys_atom_nums.push(num);
        let locator = locator + 4;

        println!("Reading the {}/{} system atoms information.", mol_id + 1, mol_num.len());
        println!("Reading atom property parameters...");
        let pb = ProgressBar::new(sys_atom_nums[mol_id] as u64);   // progress bar
        for i in 0..sys_atom_nums[mol_id] {
            let re = Regex::new(r".*type=\s*(\d+).*q=\s*([^,]+),.*resind=\s*(\d+).*").unwrap();
            let c = re.captures(&mdp_content[locator + i]).unwrap();
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
        pb.reset();
        let locator = locator + sys_atom_nums[mol_id] + 1;     // 不加1是"atom (3218):"行

        // get atom names
        println!("Reading atom names...");
        let pb = ProgressBar::new(sys_atom_nums[mol_id] as u64);
        for i in 0..sys_atom_nums[mol_id] {
            let re = Regex::new("name=\"(.*)\"").unwrap();
            let name = re.captures(&mdp_content[locator + i]).unwrap();
            let name = name.get(1).unwrap().as_str();
            t_atoms[[mol_id, i]] = name.to_string();
            pb.inc(1);
        }
        pb.reset();
    }

    // get residues information
    let re = Regex::new(r"\s*residue \((\d+)\)").unwrap();
    let locators = get_md_locators_all(&mdp_content, &re);
    let mut resnums: Vec<usize> = vec![];
    for i in 0..locators.len() {
        let res_num = re.captures(&mdp_content[locators[i]]).unwrap().get(1).unwrap();
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
            let m = re.captures(&mdp_content[locator + 1 + i]).unwrap();
            let name = m.get(1).unwrap().as_str();
            let nr = m.get(2).unwrap().as_str();
            let nr: i32 = nr.parse().unwrap();
            res_names[[idx, i]] = format!("{:05}{}", nr, name);
            pb.inc(1);
        }
        pb.reset();
    }

    // assign H types by connection atoms from angle information
    let re = Regex::new(r"^ +Angle:").unwrap();
    let locators = get_md_locators_all(&mdp_content, &re);
    println!("Reading angles...");
    for (mol_id, locator) in locators.into_iter().enumerate() {
        let angles_num: Vec<&str> = (&mdp_content[locator + 1]).trim().split(" ").collect();
        let angles_num: usize = angles_num[1].parse().unwrap();
        let angles_num = angles_num / 4;
        if angles_num > 0 {
            let re = Regex::new(r"\d+ type=\d+ \(ANGLES\)\s+(\d+)\s+(\d+)\s+(\d+)").unwrap();
            let pb = ProgressBar::new(angles_num as u64);
            for l_num in locator + 3..locator + 3 + angles_num {
                if re.is_match(&mdp_content[l_num]) {
                    let paras = re.captures(&mdp_content[l_num]).unwrap();
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
            pb.reset();
        }
    }

    // output to qrv file
    let mut atom_id_total = 0;
    let mut atom_id_feature = 0;
    for i in 0..grp_num {
        for n in 0..mol_num[i] {
            println!("Writing atoms...");
            let pb = ProgressBar::new(sys_atom_nums[i] as u64);
            for j in 0..sys_atom_nums[i] {
                if ndx_rec.contains(&atom_id_total) || ndx_lig.contains(&atom_id_total) {
                    atom_id_feature += 1;
                    let mut radi: f64;
                    match rad_type {
                        0 => radi = r_atoms[[n, j]],
                        1 => radi = {
                            let re = Regex::new(r"([a-zA-Z]+)\d*").unwrap();
                            let res = re.captures(t_atoms[[i, j]].as_str()).unwrap();
                            let res = res.get(1).unwrap().as_str();
                            get_radi(res)
                        },
                        _ => {
                            println!("Error: radType should only be 0 or 1. Check settings.");
                            return;
                        }
                    }
                    qrv_content.write_all(format!("{:6} {:9.5} {:9.6} {:6} {:9.6} {:9.6} {:6} \"{}\"-1.{} {} {:-6}  ",
                                                  atom_id_feature, q_atoms[[i, j]], radi, c_atoms[[i, j]], s_atoms[[i, j]],
                                                  e_atoms[[i, j]], atom_id_total + 1, sys_names[i], j + 1,
                                                  res_names[[i, res_ids[[i, j]]]], t_atoms[[i, j]]).as_bytes())
                        .expect("Writing qrv file failed.");
                    if ndx_rec.contains(&atom_id_total) {
                        qrv_content.write_all("Rec\n".as_bytes()).expect("Writing qrv file failed.");
                    } else if ndx_lig.contains(&atom_id_total) {
                        qrv_content.write_all("Lig\n".as_bytes()).expect("Writing qrv file failed.");
                    }
                }
                atom_id_total += 1;
                pb.inc(1);
            }
            pb.reset();
        }
    }

    // generate md5 for tpr file
    if Path::new(mdp).is_file() {
        println!("Generating mdp sha...");
        let mut mdp_sha = fs::File::create(wd.join(".mdp.sha")).unwrap();
        mdp_sha.write_all(gen_file_sha256(mdp).as_bytes()).expect("Failed to write sha");
    }
    if Path::new(qrv).is_file() {
        println!("Generating qrv sha...");
        let mut qrv_sha = fs::File::create(wd.join(".qrv.sha")).unwrap();
        qrv_sha.write_all(gen_file_sha256(qrv).as_bytes()).expect("Failed to write sha");
    }

    println!("Finished generating qrv file.");
}

fn get_md_locators_first(strings: &Vec<&str>, re: &Regex) -> usize {
    for (idx, l) in strings.into_iter().enumerate() {
        if re.is_match(l) {
            return idx;
        }
    }
    return 0;
}

fn get_md_locators_all(strings: &Vec<&str>, re: &Regex) -> Vec<usize> {
    let mut locators: Vec<usize> = vec![];
    for (idx, l) in strings.into_iter().enumerate() {
        if re.is_match(l) {
            locators.push(idx);
        }
    }
    return locators;
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