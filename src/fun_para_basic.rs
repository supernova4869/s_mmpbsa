use std::io::stdin;
use std::path::Path;
use std::thread;
use colored::*;
use crate::parameters::Config;
use crate::settings::Settings;
use crate::utils::{append_new_name, get_input, get_input_selection, make_ndx};
use crate::{confirm_file_validity, convert_cur_dir, set_program};
use crate::fun_para_system;
use crate::parse_tpr::TPR;

fn list_basic_programs(settings: &mut Settings) {
    println!(" -5 Set number of parallel kernels (only used by PBSA calculations), current: {}", settings.nkernels);
    println!(" -4 Set delphi path, current: {}", match &settings.delphi_path {
        Some(s) => s.to_string(),
        None => String::from("Not set")
    });
    println!(" -3 Set apbs path, current: {}", match &settings.apbs_path {
        Some(s) => s.to_string(),
        None => String::from("Not set")
    });
    println!(" -2 Set PBSA kernel, current: {}", match &settings.pbsa_kernel {
        Some(s) => s.to_string(),
        None => String::from("Not set")
    });
    println!(" -1 Toggle whether debug mode, current: {}", settings.debug_mode);
}

fn set_basic_programs(opt: i32, settings: &mut Settings) {
    match opt {
        -1 => settings.debug_mode = !settings.debug_mode,
        -2 => {
            println!("Input PBSA kernel (if empty, means not to do PBSA calculation):");
            let s: String = get_input_selection().unwrap();
            if s.eq("apbs") || s.eq("delphi") {
                settings.pbsa_kernel = Some(s)
            } else {
                settings.pbsa_kernel = None
            }
        }
        -3 => {
            println!("Input APBS path:");
            let s: String = get_input_selection().unwrap();
            match set_program(&Some(s), "apbs", settings) {
                Some(s) => settings.apbs_path = Some(s),
                None => settings.apbs_path = None
            }
        }
        -4 => {
            println!("Input Delphi path:");
            let s: String = get_input_selection().unwrap();
            match set_program(&Some(s), "delphi", settings) {
                Some(s) => settings.delphi_path = Some(s),
                None => settings.delphi_path = None
            }
        }
        -5 => {
            let num_cores = thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1);
            
            // 使用一半的核心
            let threads_to_use = (num_cores / 2).max(1);
            
            println!("Input number of parallel kernels (default: system cores num, currently {}):", threads_to_use);
            settings.nkernels = get_input(threads_to_use).max(1);
        }
        _ => {}
    }
}

pub fn set_para_basic_tpr(tpr_path: &String, trj_path: &Option<String>, ndx_path: &Option<String>, 
                            config: &Option<Config>, wd: &Path, settings: &mut Settings) {
    let mut trj = trj_path.clone().unwrap_or(String::new());
    let mut ndx = ndx_path.clone().unwrap_or(String::new());
    let mut tpr = TPR::from(&tpr_path);
    println!("\nFinished loading tpr file: {}", tpr);

    loop {
        println!("\n                 ************ MM/PB-SA Files ************");
        println!("-10 Exit program");
        list_basic_programs(settings);
        if !trj.is_empty() && !ndx.is_empty() {
            println!("{}", "  0 Go to next step (complete)".green().bold());
        } else {
            println!("{}", "  0 Go to next step (incomplete)".red().bold());
        }
        println!("  1 Assign trajectory file (xtc, trr, pdb, gro), current: {}", match trj.len() {
            0 => "undefined",
            _ => trj.as_str()
        });
        println!("  2 Assign index file (ndx), current: {}", match ndx.len() {
            0 => "undefined",
            _ => ndx.as_str()
        });
        let i = get_input_selection();
        match i {
            Ok(-1) => settings.debug_mode = !settings.debug_mode,
            Ok(0) => {
                if trj.is_empty() {
                    println!("Please assign trajectory file, including xtc, trr, gro or pdb.");
                } else if ndx.is_empty() {
                    println!("Please assign index file.");
                } else {
                    // go to next step
                    fun_para_system::set_para_trj(&trj, &mut tpr, &ndx, config, &wd, &tpr_path, settings);
                }
            }
            Ok(1) => {
                println!("Input trajectory file path, default: ?md.xtc (\"?\" means the same directory as tpr):");
                trj.clear();
                stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                if trj.trim().is_empty() {
                    trj = "?md.xtc".to_string();
                }
                trj = convert_cur_dir(&trj, tpr_path);
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr", "gro", "pdb", "pdbqt"], tpr_path);
            }
            Ok(2) => {
                println!("Input index file path, default: ?index.ndx (\"?\" means the same directory as tpr):");
                println!("Note: if no index file prepared, the default index.ndx will be generated according to tpr.");
                ndx.clear();
                stdin().read_line(&mut ndx).expect("Failed while reading index file");
                if ndx.trim().is_empty() {
                    ndx = "?index.ndx".to_string();
                }
                ndx = convert_cur_dir(&ndx, tpr_path);
                if !Path::new(&ndx).is_file() {
                    let tpr_path = append_new_name(tpr_path, ".tpr", "");
                    make_ndx(&vec!["q"], wd, settings, &tpr_path, "", &ndx);
                }
                ndx = confirm_file_validity(&mut ndx, vec!["ndx", "pdbqt"], tpr_path);
            }
            Ok(-10) => break,
            Ok(other) => {
                set_basic_programs(other, settings);
            },
            Err(_) => {}
        };
    }
}
