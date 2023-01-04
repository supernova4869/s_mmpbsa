use std::io::stdin;
use std::path::Path;
use crate::{get_input_selection, index_parser, parameters::Parameters};
use crate::{mmpbsa, analyzation};
use crate::apbs_param::{PBASet, PBESet};
use crate::atom_radius::{Radius, RADIUS_TABLE};
use crate::parse_tpr::TPR;

pub fn set_para_mmpbsa(trj: &String, tpr: &mut TPR, ndx: &String, wd: &Path,
                       complex_grp: usize,
                       receptor_grp: usize,
                       ligand_grp: usize,
                       bt: f64, et: f64, dt: f64,
                       atom_radius: &Radius,
                       settings: &mut Parameters) {
    // save a copy of default force field atom type
    let atom_radius_ff = atom_radius.clone();
    tpr.apply_radius(settings.rad_type, &atom_radius.radii);
    let mut pbe_set = &PBESet::new(tpr.temp);
    let mut pba_set = &PBASet::new();
    loop {
        println!("\n                 ************ MM/PB-SA Parameters ************");
        println!("-10 Return");
        println!("  0 Start MM/PB-SA calculation");
        println!("  1 Toggle whether to use Debye-Huckel shielding method, current: {}", settings.use_dh);
        println!("  2 Toggle whether to use entropy contribution, current: {}", settings.use_ts);
        println!("  3 Select atom radius type, current: {}", RADIUS_TABLE[&settings.rad_type]);
        println!("  4 Input coarse grid expand factor (cfac), current: {}", settings.cfac);
        println!("  5 Input fine grid expand amount (fadd), current: {} A", settings.fadd);
        println!("  6 Input atom distance cutoff for MM calculation (A), current: {}", settings.r_cutoff);
        println!("  7 Input fine mesh spacing (df), current: {} A", settings.df);
        println!("  8 Prepare PB parameters for APBS");
        println!("  9 Prepare SA parameters for APBS");
        let i = get_input_selection();
        match i {
            -10 => return,
            0 => {
                let mut sys_name = String::from("_system");
                println!("Input system name (default: {}):", sys_name);
                let mut input = String::new();
                stdin().read_line(&mut input).expect("Error input");
                if input.trim().len() != 0 {
                    sys_name = input.trim().to_string();
                }
                // 定义results形式, 其中应包含所需的全部数据
                let ndx = index_parser::Index::new(ndx);
                let results = mmpbsa::fun_mmpbsa_calculations(trj, tpr, &ndx, wd, &sys_name,
                                                              complex_grp as usize,
                                                              receptor_grp as usize,
                                                              ligand_grp as usize,
                                                              bt, et, dt,
                                                              pbe_set, pba_set, settings);
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
                    for (k, v) in RADIUS_TABLE.clone().into_iter() {
                        s.push_str(format!("\n{}):\t{}", k, v).as_str());
                    }
                    s
                });
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.rad_type = 1;
                } else {
                    let s = s.trim().parse().unwrap();
                    if s == 0 {
                        settings.rad_type = 0;
                    } else if RADIUS_TABLE.contains_key(&s) {
                        settings.rad_type = s;
                    } else {
                        println!("Radius type {} not supported. Will use mBondi instead.", s);
                        settings.rad_type = 1;
                    }
                }
                tpr.apply_radius(settings.rad_type, &atom_radius_ff.radii);
            }
            4 => {
                println!("Input coarse grid expand factor, default 3:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.cfac = 3.0;
                } else {
                    settings.cfac = s.trim().parse().unwrap();
                }
            }
            5 => {
                println!("Input fine grid expand amount (A), default 10:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.fadd = 10.0;
                } else {
                    settings.fadd = s.trim().parse().unwrap();
                }
            }
            6 => {
                println!("Input cutoff value (A), default 0 (inf):");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.r_cutoff = f64::INFINITY;
                } else {
                    settings.r_cutoff = s.trim().parse().unwrap();
                    if settings.r_cutoff == 0.0 {
                        settings.r_cutoff = f64::INFINITY;
                    }
                }
            }
            7 => {
                println!("Input fine mesh spacing (A), default 0.5:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.df = 0.5;
                } else {
                    settings.df = s.trim().parse().unwrap();
                }
            }
            8 => {
                pbe_set.save_params(wd.join("PB_settings.txt"));
                println!("PB parameters have been wrote to PB_settings.txt. \
                    Edit it and press return to reload it.");
                stdin().read_line(&mut String::new()).unwrap();
                // pbe_set = &pbe_set.load_params(wd.join("PBESet.txt"));
            }
            9 => {
                pba_set.save_params(wd.join("SA_settings.txt"));
                println!("SA parameters have been wrote to SA_settings.txt. \
                    Edit it and press return to reload it.");
                stdin().read_line(&mut String::new()).unwrap();
            }
            _ => println!("Invalid input")
        }
    }
}
