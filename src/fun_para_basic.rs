use std::io::stdin;
use std::path::Path;
use std::process::exit;
use std::{env, thread};
use colored::*;
use crate::parameters::Config;
use crate::settings::Settings;
use crate::utils::{get_input, get_input_selection, make_ndx};
use crate::{confirm_file_validity, convert_ask_dir};
use crate::fun_para_system;
use crate::parse_tpr::TPR;

fn list_basic_programs(settings: &mut Settings) {
    println!(" -4 Set number of parallel kernels, current: {}", settings.nkernels);
    println!(" -3 Toggle whether do PBSA calculations, current: {}", settings.calc_pbsa);
    println!(" -2 Toggle whether do MM calculations, current: {}", settings.calc_mm);
    println!(" -1 Toggle whether debug mode, current: {}", settings.debug_mode);
}

fn set_basic_programs(opt: i32, settings: &mut Settings) {
    match opt {
        -1 => settings.debug_mode = !settings.debug_mode,
        -2 => settings.calc_mm = !settings.calc_mm,
        -3 => settings.calc_pbsa = !settings.calc_pbsa,
        -4 => {
            let num_cores = thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1);
            
            // 使用一半的核心
            let threads_to_use = (num_cores / 2).max(1);
            
            println!("Input number of parallel kernels (default: {}, currently {}):", 
                threads_to_use, settings.nkernels);
            settings.nkernels = get_input(threads_to_use).max(1);
        }
        _ => {}
    }
}

pub fn set_para_basic_tpr(tpr_dump_path: &String, trj_path: &Option<String>, 
                            tpr_path: &String, ndx_path: &Option<String>, 
                            config: &Option<Config>, settings: &mut Settings) {
    let mut trj = trj_path.clone().unwrap_or("".to_string());
    let mut ndx = ndx_path.clone().unwrap_or("".to_string());
    let mut tpr = TPR::from(&tpr_dump_path);
    println!("\nFinished loading tpr file: {}", tpr);
    if config.is_some() {
        settings.calc_mm = config.as_ref().unwrap().program_set.calc_mm;
        settings.calc_pbsa = config.as_ref().unwrap().program_set.calc_pbsa;
        let trj = config.as_ref().unwrap().program_set.trj.clone();
        let ndx = config.as_ref().unwrap().program_set.ndx.clone();
        let ndx = if ndx.is_empty() {
            "index.ndx".to_string()
        } else {
            ndx
        };
        if !Path::new(&trj).is_file() {
            println!("Not valid trajectory file in config: {}. Check again.", trj);
            exit(0);
        }
        if !Path::new(&ndx).is_file() {
            println!("{} not found. Generating default index.ndx.", ndx);
            make_ndx(&vec!["q"], &env::current_dir().unwrap(), settings, &tpr_path, "", &ndx);
        }
        fun_para_system::set_para_trj(&trj, &mut tpr, &ndx, config, &tpr_path, settings);
    }

    loop {
        println!("\n                 ************ MM-PBSA Files ************");
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
                    fun_para_system::set_para_trj(&trj, &mut tpr, &ndx, config, &tpr_path, settings);
                }
            }
            Ok(1) => {
                println!("Input trajectory file path, default: ?md.xtc (\"?\" means the same directory as tpr):");
                trj.clear();
                stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                trj = trj.trim().to_string();
                if trj.is_empty() {
                    trj = "?md.xtc".to_string();
                }
                trj = convert_ask_dir(&trj, tpr_path);
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr", "gro", "pdb", "pdbqt"], tpr_dump_path);
            }
            Ok(2) => {
                println!("Input index file path, default: ?index.ndx (\"?\" means the same directory as tpr):");
                println!("Note: if no index file prepared, the default index.ndx will be generated according to tpr.");
                ndx.clear();
                stdin().read_line(&mut ndx).expect("Failed while reading index file");
                ndx = ndx.trim().to_string();
                if ndx.is_empty() {
                    ndx = "?index.ndx".to_string();
                }
                ndx = convert_ask_dir(&ndx, tpr_path);
                if !Path::new(&ndx).is_file() {
                    make_ndx(&vec!["q"], &env::current_dir().unwrap(), settings, &tpr_path, "", &ndx);
                }
                ndx = confirm_file_validity(&mut ndx, vec!["ndx", "pdbqt"], tpr_dump_path);
            }
            Ok(-10) => break,
            Ok(other) => {
                set_basic_programs(other, settings);
            },
            Err(_) => {}
        };
    }
}
