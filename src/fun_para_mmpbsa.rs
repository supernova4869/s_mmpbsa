use std::io::stdin;
use std::path::Path;
use crate::utils::{self, get_input, get_input_selection, get_residue_range_ca};
use crate::index_parser::Index;
use crate::settings::Settings;
use crate::apbs_param::{PBASet, PBESet};
use std::io::Write;
use std::fs::{File, self};
use crate::atom_property::AtomProperties;
use crate::parse_tpr::{Residue, TPR};
use crate::{mmpbsa, parse_xvg};
use crate::analyzation;

pub fn set_para_mmpbsa(tpr: &mut TPR, ndx: &Index, wd: &Path, aps: &mut AtomProperties,
                       ndx_rec: &Vec<usize>, ndx_lig: &Vec<usize>,
                       receptor_grp: usize, ligand_grp: Option<usize>,
                       residues: &Vec<Residue>, settings: &mut Settings) {
    // kinds of radius types
    let radius_types = vec!["ff", "amber", "Bondi", "mBondi", "mBondi2"];
    let mut pbe_set = PBESet::new(tpr.temp);
    let mut pba_set = PBASet::new(tpr.temp);
    let mut ala_list: Vec<i32> = vec![];
    loop {
        println!("\n                 ************ MM/PB-SA Parameters ************");
        println!("-10 Return");
        println!(" -3 Output PBSA parameters");
        println!(" -2 Output LJ parameters");
        println!(" -1 Output structural parameters");
        println!("  0 Start MM/PB-SA calculation");
        println!("  1 Choose electrostatic screening method, current: {}", match settings.elec_screen {
            1 => "Ding's method",
            2 => "Supernova's method",
            _ => "None"
        });
        println!("  2 Select residues list for alanine scanning, current: {:?}", ala_list);
        println!("  3 Select atom radius type, current: {}", radius_types[settings.radius_type]);
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
                let mut paras = File::create(wd.join("paras_atom_properties.txt")).unwrap();
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
                paras.write_all(format!("Atom radius type: {}\n", radius_types[settings.radius_type]).as_bytes()).unwrap();
                paras.write_all(format!("Atoms:\n     id   name   type   charge   radius   resnum  resname\n").as_bytes()).unwrap();
                for ap in &aps.atom_props {
                    paras.write_all(format!("{:7}{:>7}{:7}{:9.2}{:9.2}{:9}{:>9}\n", 
                        ap.id, ap.name, ap.type_id, ap.charge, ap.radius, ap.resid + 1, ap.resname).as_bytes()).unwrap();
                }
                println!("Structural parameters have been written to paras_atom_properties.txt");
            }
            -2 => {
                let mut paras = File::create(wd.join("paras_LJ.txt")).unwrap();
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
                println!("Forcefield parameters have been written to paras_LJ.txt");
            }
            -3 => {
                let mut paras = File::create(wd.join("paras_pbsa.txt")).unwrap();
                paras.write_all(format!("Electrostatic screening method: {}\n", settings.elec_screen).as_bytes()).unwrap();
                paras.write_all(format!("Atom radius type: {}\n", radius_types[settings.radius_type]).as_bytes()).unwrap();
                paras.write_all(format!("Atom distance cutoff for MM calculation (A): {}\n", settings.r_cutoff).as_bytes()).unwrap();
                paras.write_all(format!("Coarse grid expand factor (cfac): {}\n", settings.cfac).as_bytes()).unwrap();
                paras.write_all(format!("Fine grid expand amount (fadd): {} A\n", settings.fadd).as_bytes()).unwrap();
                paras.write_all(format!("Fine mesh spacing (df): {} A\n\n", settings.df).as_bytes()).unwrap();
                paras.write_all(format!("PB settings:\n{}\n\n", pbe_set).as_bytes()).unwrap();
                paras.write_all(format!("SA settings:\n{}\n", pba_set).as_bytes()).unwrap();
                println!("PBSA parameters have been written to paras_pbsa.txt");
            }
            0 => {
                println!("Applying {} radius...", radius_types[settings.radius_type]);
                aps.apply_radius(settings.radius_type, &tpr.get_at_list(), &radius_types, wd);

                // Temp directory for PBSA
                let mut sys_name = String::from("_system");
                println!("Input system name (default: {}):", sys_name);
                let mut input = String::new();
                stdin().read_line(&mut input).expect("Error input");
                if input.trim().len() != 0 {
                    sys_name = input.trim().to_string();
                }
                let temp_dir = wd.join(&sys_name);
                if let Some(_) = settings.apbs_path.as_ref() {
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
                    println!("Note: Since APBS not found, solvation energy will not be calculated.");
                };
                
                // run MM/PB-SA calculations
                println!("Extracting atoms coordination...");
                let (time_list, coordinates) = parse_xvg::read_coord_xvg(wd.join("_MMPBSA_coord.xvg").to_str().unwrap());
                
                let (result_wt, result_as) = mmpbsa::fun_mmpbsa_calculations(&time_list, &coordinates, &temp_dir, &sys_name, &aps,
                                                                &ndx_rec, &ndx_lig, &ala_list, &residues, wd,
                                                                &pbe_set, &pba_set, settings);
                analyzation::analyze_controller(&result_wt, &result_as, pbe_set.temp, &sys_name, wd, settings);
            }
            1 => {
                println!("Input the electrostatic screening method:");
                println!("0: no screening\n1: Ding's method\n2: Supernova's method");
                settings.elec_screen = get_input(1);
            }
            2 => {
                println!("Select the residues for alanine scanning:");
                println!(" 1 Select the residues within the first layer (0-4 A)");
                println!(" 2 Select the residues within the second layer (4-6 A)");
                println!(" 3 Select the residues within the third layer (6-8 A)");
                println!(" 4 Select the residues within specific distance");
                println!(" 5 Directly input the resudues list");
                let i: i32 = get_input_selection();
                let receptor_res: Vec<Residue> = residues.iter().filter_map(|r| if r.id != aps.atom_props[ndx_lig[0]].resid {
                    Some(r.clone())
                } else {
                    None
                }).collect();
                let atom_res = &aps.atom_props.iter().map(|a| a.resid).collect();
                let atom_names = &aps.atom_props.iter().map(|a| a.name.to_string()).collect();
                match i {
                    1 => {
                        let rs = get_residue_range_ca(&tpr.coordinates, ndx_lig, 4.0, 
                            &atom_res, &atom_names, &receptor_res);
                        ala_list = rs.iter().filter_map(|&i| Some(residues[i].nr)).collect();
                    },
                    2 => {
                        let rs = get_residue_range_ca(&tpr.coordinates, ndx_lig, 6.0, 
                            &atom_res, &atom_names, &receptor_res);
                        let inner_rs = get_residue_range_ca(&tpr.coordinates, ndx_lig, 4.0, 
                            &atom_res, &atom_names, &receptor_res);
                        ala_list = rs.iter().filter_map(|&i| if !inner_rs.contains(&i) {
                            Some(residues[i].nr)
                        } else {
                            None
                        } ).collect();
                    },
                    3 => {
                        let rs = get_residue_range_ca(&tpr.coordinates, ndx_lig, 8.0, 
                            &atom_res, &atom_names, &receptor_res);
                        let inner_rs = get_residue_range_ca(&tpr.coordinates, ndx_lig, 6.0, 
                            &atom_res, &atom_names, &receptor_res);
                        ala_list = rs.iter().filter_map(|&i| if !inner_rs.contains(&i) {
                            Some(residues[i].nr)
                        } else {
                            None
                        } ).collect();
                    },
                    4 => {
                        println!("Input the cut-off distance you want to expand from ligand, default: 4 A");
                        let cutoff = get_input(4.0);
                        let rs = get_residue_range_ca(&tpr.coordinates, ndx_lig, cutoff, 
                            &atom_res, &atom_names, &receptor_res);
                        ala_list = rs.iter().filter_map(|&i| Some(residues[i].nr)).collect();
                    },
                    5 => {
                        println!("Input the residues list for alanine scanning:");
                        let rs = get_input("".to_string());
                        ala_list = utils::range2list(rs.as_str());
                    },
                    _ => {}
                }
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
                    settings.radius_type = 3;
                } else {
                    let s = s.trim().parse().expect("Input not valid number.");
                    if s == 0 {
                        settings.radius_type = 0;
                    } else if s < radius_types.len() {
                        settings.radius_type = s;
                    } else {
                        println!("Radius type {} not supported. Will use mBondi instead.", radius_types[s]);
                        settings.radius_type = 3;
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
                let pb_fpath = wd.join("PB_settings.yaml");
                pbe_set.save_params(&pb_fpath);
                println!("PB parameters have been wrote to {0}.\n\
                    Edit it and input its path to reload (default: {0}).", &pb_fpath.to_str().unwrap());
                let pb_fpath = get_input(pb_fpath.to_str().unwrap().to_string());
                pbe_set = PBESet::load_params(pb_fpath);
            }
            9 => {
                let sa_fpath = wd.join("SA_settings.yaml");
                pba_set.save_params(&sa_fpath);
                println!("SA parameters have been wrote to {0}.\n\
                    Edit it and input its path to reload (default: {0}).", &sa_fpath.to_str().unwrap());
                let sa_fpath = get_input(sa_fpath.to_str().unwrap().to_string());
                pba_set = PBASet::load_params(sa_fpath);
            }
            _ => println!("Invalid input")
        }
    }
}
