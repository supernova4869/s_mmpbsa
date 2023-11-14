use std::cmp::Ordering;
use std::io::stdin;
use std::path::Path;
use crate::index_parser::Index;
use crate::{get_input_selection, settings::Settings};
use crate::apbs_param::{PBASet, PBESet};
use std::io::Write;
use std::fs::{File, self};
use crate::atom_property::AtomProperty;
use crate::parse_tpr::TPR;
use crate::mmpbsa::{self, get_residues};
use crate::analyzation::{self, get_infile};

pub fn set_para_mmpbsa(trj: &String, tpr: &mut TPR, ndx: &Index, wd: &Path,
                       receptor_grp: usize, ligand_grp: Option<usize>,
                       bt: f64, et: f64, dt: f64, settings: &mut Settings) {
    // atom indexes
    println!("Preparing atom indexes...");
    let ndx_lig = match ligand_grp {
        Some(ligand_grp) => Some(&ndx.groups[ligand_grp].indexes),
        None => None
    };
    let ndx_rec = &ndx.groups[receptor_grp].indexes;
    let ndx_com = match ndx_lig {
        Some(ndx_lig) => {
            match ndx_lig[0] > ndx_rec[0] {
                true => {
                    let mut ndx_com = ndx_rec.to_vec();
                    ndx_com.extend(ndx_lig);
                    ndx_com
                }
                false => {
                    let mut ndx_com = ndx_lig.to_vec();
                    ndx_com.extend(ndx_rec);
                    ndx_com
                }
            }
        }
        None => ndx_rec.to_vec()
    };

    // atom properties
    println!("Parsing atom properties...");
    let mut aps = AtomProperty::new(tpr, &ndx_com);
    
    let (ndx_com_norm, ndx_rec_norm, ndx_lig_norm) = normalize_index(ndx_rec, ndx_lig);
    
    // kinds of radius types
    let radius_types = vec!["ff", "amber", "Bondi", "mBondi", "mBondi2"];
    let mut pbe_set = PBESet::new(tpr.temp);
    let mut pba_set = PBASet::new(tpr.temp);
    loop {
        println!("\n                 ************ MM/PB-SA Parameters ************");
        println!("-10 Return");
        println!(" -3 Output PBSA parameters");
        println!(" -2 Output structural parameters");
        println!(" -1 Output ff parameters");
        println!("  0 Start MM/PB-SA calculation");
        println!("  1 Toggle whether to use Debye-Huckel shielding method, current: {}", settings.use_dh);
        println!("  2 Toggle whether to use entropy contribution, current: {}", settings.use_ts);
        println!("  3 Select atom radius type, current: {}", radius_types[settings.rad_type]);
        println!("  4 Input atom distance cutoff for MM calculation (A), current: {}", settings.r_cutoff);
        println!("  5 Input coarse grid expand factor (cfac), current: {}", settings.cfac);
        println!("  6 Input fine grid expand amount (fadd), current: {} A", settings.fadd);
        println!("  7 Input fine mesh spacing (df), current: {} A", settings.df);
        println!("  8 Prepare PB parameters for APBS");
        println!("  9 Prepare SA parameters for APBS");
        let i = get_input_selection();
        match i {
            -10 => return,
            -1 => {
                let mut paras = File::create(wd.join("paras_ff.txt")).unwrap();
                paras.write_all(format!("Atom types num: {}\n", tpr.atom_types_num).as_bytes()).unwrap();
                paras.write_all("c6:\n".as_bytes()).unwrap();
                for i in 0..aps.c6.shape()[0] {
                    for j in 0..aps.c6.shape()[1] {
                        paras.write_all(format!("{:13.6E} ", aps.c6[[i, j]]).as_bytes()).unwrap();
                    }
                    paras.write_all("\n".as_bytes()).unwrap();
                }
                paras.write_all("c12:\n".as_bytes()).unwrap();
                for i in 0..aps.c12.shape()[0] {
                    for j in 0..aps.c12.shape()[1] {
                        paras.write_all(format!("{:13.6E} ", aps.c12[[i, j]]).as_bytes()).unwrap();
                    }
                    paras.write_all("\n".as_bytes()).unwrap();
                }
                println!("Forcefield parameters have been written to paras_ff.txt");
            }
            -2 => {
                let mut paras = File::create(wd.join("paras_structure.txt")).unwrap();
                paras.write_all(format!("Receptor group: {}\n", 
                    ndx.groups[receptor_grp as usize].name).as_bytes()).unwrap();
                match ligand_grp {
                    Some(ligand_grp) => {
                        paras.write_all(format!("Ligand group: {}\n", 
                            ndx.groups[ligand_grp as usize].name).as_bytes()).unwrap();
                    }
                    None => {
                        paras.write_all("Ligand group: None\n".as_bytes()).unwrap();
                    }
                }
                paras.write_all(format!("Atom radius type: {}\n", radius_types[settings.rad_type]).as_bytes()).unwrap();
                paras.write_all(format!("Atoms:\n     id   name   type        sigma      epsilon   charge   radius   resnum  resname\n").as_bytes()).unwrap();
                for idx in 0..ndx_com_norm.len() {
                    paras.write_all(format!("{:7}{:>7}{:7}{:13.6E}{:13.6E}{:9.2}{:9.2}{:9}{:>9}\n", 
                    aps.atm_index[idx], aps.atm_name[idx], aps.atm_typeindex[idx], aps.atm_sigma[idx], 
                    aps.atm_epsilon[idx], aps.atm_charge[idx], aps.atm_radius[idx], aps.atm_resnum[idx] + 1, 
                    aps.atm_resname[idx]).as_bytes()).unwrap();
                }
                println!("Structural parameters have been written to paras_structure.txt");
            }
            -3 => {
                let mut paras = File::create(wd.join("paras_pbsa.txt")).unwrap();
                paras.write_all(format!("Use Debye-Huckel shielding method: {}\n", settings.use_dh).as_bytes()).unwrap();
                paras.write_all(format!("Use entropy contribution: {}\n", settings.use_ts).as_bytes()).unwrap();
                paras.write_all(format!("Atom radius type: {}\n", radius_types[settings.rad_type]).as_bytes()).unwrap();
                paras.write_all(format!("Atom distance cutoff for MM calculation (A): {}\n", settings.r_cutoff).as_bytes()).unwrap();
                paras.write_all(format!("Coarse grid expand factor (cfac): {}\n", settings.cfac).as_bytes()).unwrap();
                paras.write_all(format!("Fine grid expand amount (fadd): {} A\n", settings.fadd).as_bytes()).unwrap();
                paras.write_all(format!("Fine mesh spacing (df): {} A\n\n", settings.df).as_bytes()).unwrap();
                paras.write_all(format!("PB settings:\n{}\n\n", pbe_set).as_bytes()).unwrap();
                paras.write_all(format!("SA settings:\n{}\n", pba_set).as_bytes()).unwrap();
                println!("PBSA parameters have been written to paras_pbsa.txt");
            }
            0 => {
                println!("Applying {} radius...", radius_types[settings.rad_type]);
                aps.apply_radius(settings.rad_type, tpr, &ndx_com_norm, &radius_types);

                // Temp directory for PBSA
                let mut sys_name = String::from("_system");
                println!("Input system name (default: {}):", sys_name);
                let mut input = String::new();
                stdin().read_line(&mut input).expect("Error input");
                if input.trim().len() != 0 {
                    sys_name = input.trim().to_string();
                }
                let temp_dir = wd.join(&sys_name);
                if let Some(_) = settings.apbs.as_ref() {
                    println!("Temporary files will be placed at {}/", temp_dir.display());
                    if !temp_dir.is_dir() {
                        fs::create_dir(&temp_dir).expect(format!("Failed to create temp directory: {}.", &sys_name).as_str());
                    } else {
                        println!("Directory {}/ not empty. Clear? [Y/n]", temp_dir.display());
                        let mut input = String::from("");
                        stdin().read_line(&mut input).expect("Get input error");
                        if input.trim().len() == 0 || input.trim() == "Y" || input.trim() == "y" {
                            fs::remove_dir_all(&temp_dir).expect("Remove dir failed");
                            fs::create_dir(&temp_dir).expect(format!("Failed to create temp directory: {}.", &sys_name).as_str());
                        }
                    }
                } else {
                    println!("Note: APBS not found. Will not calculate solvation energy.");
                };
                println!("Collecting residues list...");
                let residues = get_residues(tpr, &ndx_com);
                let results = mmpbsa::fun_mmpbsa_calculations(trj, &temp_dir, &sys_name, &aps,
                                                              &ndx_com_norm, &ndx_rec_norm, &ndx_lig_norm, residues,
                                                              bt, et, dt, &pbe_set, &pba_set, settings);
                analyzation::analyze_controller(&results, pbe_set.temp, &sys_name, wd);
            }
            1 => {
                settings.use_dh = !settings.use_dh;
            }
            2 => {
                settings.use_ts = !settings.use_ts;
            }
            3 => {
                println!("Input atom radius type (default mBondi), Supported:{}", {
                    let mut s = String::new();
                    for (k, v) in radius_types.iter().enumerate() {
                        s.push_str(format!("\n{}):\t{}", k, v).as_str());
                    }
                    s
                });
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.rad_type = 3;
                } else {
                    let s = s.trim().parse().expect("Input not valid number.");
                    if s == 0 {
                        settings.rad_type = 0;
                    } else if s < radius_types.len() {
                        settings.rad_type = s;
                    } else {
                        println!("Radius type {} not supported. Will use mBondi instead.", radius_types[s]);
                        settings.rad_type = 3;
                    }
                }
            }
            4 => {
                println!("Input cutoff value (A), default 0 (inf):");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.r_cutoff = f64::INFINITY;
                } else {
                    settings.r_cutoff = s.trim().parse().expect("Input not valid number.");
                    if settings.r_cutoff == 0.0 {
                        settings.r_cutoff = f64::INFINITY;
                    }
                }
            }
            5 => {
                println!("Input coarse grid expand factor, default 3:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.cfac = 3.0;
                } else {
                    settings.cfac = s.trim().parse().expect("Input not valid number.");
                }
            }
            6 => {
                println!("Input fine grid expand amount (A), default 10:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.fadd = 10.0;
                } else {
                    settings.fadd = s.trim().parse().expect("Input not valid number.");
                }
            }
            7 => {
                println!("Input fine mesh spacing (A), default 0.5:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.df = 0.5;
                } else {
                    settings.df = s.trim().parse().expect("Input not valid number.");
                }
            }
            8 => {
                let pb_fpath = wd.join("PB_settings.txt");
                pbe_set.save_params(&pb_fpath);
                println!("PB parameters have been wrote to {0}.\n\
                    Edit it and input its path to reload (default: {0}).", &pb_fpath.to_str().unwrap());
                let infile = get_infile(&pb_fpath.to_str().unwrap());
                pbe_set = PBESet::load_params(infile);
            }
            9 => {
                let sa_fpath = wd.join("SA_settings.txt");
                pba_set.save_params(&sa_fpath);
                println!("SA parameters have been wrote to {0}.\n\
                    Edit it and input its path to reload (default: {0}).", &sa_fpath.to_str().unwrap());
                let infile = get_infile(&sa_fpath.to_str().unwrap());
                pba_set = PBASet::load_params(infile);
            }
            _ => println!("Invalid input")
        }
    }
}

// convert rec and lig to begin at 0 and continous
pub fn normalize_index(ndx_rec: &Vec<usize>, ndx_lig: Option<&Vec<usize>>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let offset = match ndx_lig {
        Some(ndx_lig) => ndx_lig[0].min(ndx_rec[0]),
        None => ndx_rec[0]
    };
    let mut ndx_rec: Vec<usize> = ndx_rec.iter().map(|p| p - offset).collect();
    let mut ndx_lig = match ndx_lig {
        Some(ndx_lig) => ndx_lig.iter().map(|p| p - offset).collect(),
        None => ndx_rec.clone()
    };
    let ndx_com = match ndx_lig[0].cmp(&ndx_rec[0]) {
        Ordering::Greater => {
            ndx_lig = ndx_lig.iter().map(|p| p - ndx_lig[0] + ndx_rec.len()).collect();
            let mut ndx_com = ndx_rec.to_vec();
            ndx_com.extend(&ndx_lig);
            ndx_com
        }
        Ordering::Less => {
            ndx_rec = ndx_rec.iter().map(|p| p - ndx_rec[0] + ndx_lig.len()).collect();
            let mut ndx_com = ndx_lig.to_vec();
            ndx_com.extend(&ndx_rec);
            ndx_com
        }
        Ordering::Equal => Vec::from_iter(0..ndx_rec.len())
    };
    (ndx_com, ndx_rec, ndx_lig)
}
