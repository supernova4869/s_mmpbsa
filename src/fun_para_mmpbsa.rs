use std::io::stdin;
use std::path::Path;
use crate::{get_input_value, index_parser, Parameters};
use crate::{mmpbsa, analyzation};
use crate::apbs_param::{PBASet, PBESet};

enum AtomRadius {
    ForceField,
    MBondi,
}


pub fn set_para_mmpbsa(trj: &String, mdp: &String, ndx: &String, wd: &Path,
                                                             complex_grp: usize,
                                                             receptor_grp: usize,
                                                             ligand_grp: usize,
                                                             bt: f64, et: f64, dt: f64,
                                                             settings: &mut Parameters) {
    let atom_rad_type = AtomRadius::MBondi;
    let pbe_set = PBESet::new();
    let pba_set = PBASet::new();
    loop {
        println!("\n                 ************ MM/PB-SA Parameters ************");
        println!("-10 Return");
        println!("  0 Start MM/PB-SA calculation");
        println!("  1 Toggle whether to use Debye-Huckel shielding method, current: {}", settings.use_dh);
        println!("  2 Toggle whether to use entropy contribution, current: {}", settings.use_ts);
        println!("  3 Select atom radius type, current: {}", match atom_rad_type {
            AtomRadius::ForceField => "from force field",
            AtomRadius::MBondi => "mBondi"
        });
        println!("  4 Input atom radius for LJ parameters, current: not support");
        println!("  5 Input coarse grid expand factor, current: {}", settings.cfac);
        println!("  6 Input fine grid expand amount, current: {} A", settings.fadd);
        println!("  7 Input fine mesh spacing, current: {} A", settings.df);
        println!("  8 Prepare PB parameters for APBS");
        println!("  9 Prepare SA parameters for APBS");
        let i = get_input_value();
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
                let results = mmpbsa::do_mmpbsa_calculations(trj, mdp, &ndx, wd, &sys_name,
                                                             complex_grp as usize,
                                                             receptor_grp as usize,
                                                             ligand_grp as usize,
                                                             bt, et, dt,
                                                             &pbe_set, &pba_set,
                                                             &settings);
                // analyzation::analyze_controller(wd, &sys_name, results);
            }
            1 => {
                settings.use_dh = !settings.use_dh;
            }
            2 => {
                settings.use_ts = !settings.use_ts;
            }
            // 在这里加半径选择, 之前生成qrv的半径逻辑删掉
            5 => {
                println!("Input coarse grid expand factor");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.cfac = 3.0;
                } else {
                    settings.cfac = s.trim().parse().unwrap();
                }
            }
            6 => {
                println!("Input fine grid expand amount (A)");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.fadd = 10.0;
                } else {
                    settings.fadd = s.trim().parse().unwrap();
                }
            }
            7 => {
                println!("Input fine mesh spacing (A)");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.df = 0.5;
                } else {
                    settings.df = s.trim().parse().unwrap();
                }
            }
            8 => {
                println!("Current PB settings:");
                println!("{}", pbe_set);
            }
            9 => {
                println!("Current SA settings:");
                println!("{}", pba_set);
            }
            _ => println!("Invalid input")
        }
    }
}