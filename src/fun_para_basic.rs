use std::io::stdin;
use std::path::Path;
use crate::parse_pdbqt::PDBQT;
use crate::settings::Settings;
use crate::utils::get_input_selection;
use crate::{check_apbs, check_delphi, confirm_file_validity, convert_cur_dir};
use crate::fun_para_system::{set_para_trj, set_para_dock};
use crate::parse_tpr::TPR;

pub fn set_para_basic(infile: &String, wd: &Path, settings: &mut Settings) {
    if infile.ends_with("tpr") || infile.ends_with("dump") {
        let mut trj = String::new();
        let mut ndx = String::new();
        let mut tpr = TPR::new(&infile, &settings);
        println!("\nFinished loading tpr.");
    
        loop {
            println!("\n                 ************ MM/PB-SA Files ************");
            println!("-10 Exit program");
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
            println!("  0 Go to next step");
            println!("  1 Assign trajectory file (xtc or trr), current: {}", match trj.len() {
                0 => "undefined",
                _ => trj.as_str()
            });
            println!("  2 Assign index file (ndx), current: {}", match ndx.len() {
                0 => "undefined",
                _ => ndx.as_str()
            });
            let i = get_input_selection();
            match i {
                -1 => settings.debug_mode = !settings.debug_mode,
                -2 => {
                    println!("Input PBSA kernel (if empty, means not to do PBSA calculation):");
                    let s: String = get_input_selection();
                    if s.eq("apbs") || s.eq("delphi") {
                        settings.pbsa_kernel = Some(s)
                    } else {
                        settings.pbsa_kernel = None
                    }
                }
                -3 => {
                    println!("Input APBS path:");
                    let s: String = get_input_selection();
                    match check_apbs(&Some(s)) {
                        Some(s) => settings.apbs_path = Some(s),
                        None => settings.apbs_path = None
                    }
                }
                -4 => {
                    println!("Input Delphi path:");
                    let s: String = get_input_selection();
                    match check_delphi(&Some(s)) {
                        Some(s) => settings.delphi_path = Some(s),
                        None => settings.delphi_path = None
                    }
                }
                0 => {
                    if trj.len() == 0 {
                        println!("Trajectory file not assigned.");
                    } else if ndx.len() == 0 {
                        // 可能要改, 以后不需要index也能算
                        println!("Index file not assigned.");
                    } else {
                        // go to next step
                        set_para_trj(&trj, &mut tpr, &ndx, &wd, &infile, settings);
                    }
                }
                1 => {
                    println!("Input trajectory file path, default: ?md.xtc (\"?\" means the same directory as tpr):");
                    trj.clear();
                    stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                    if trj.trim().is_empty() {
                        trj = "?md.xtc".to_string();
                    }
                    trj = convert_cur_dir(&trj, &settings);
                    trj = confirm_file_validity(&mut trj, vec!["xtc", "trr"], &settings);
                }
                2 => {
                    println!("Input index file path, default: ?index.ndx (\"?\" means the same directory as tpr):");
                    ndx.clear();
                    stdin().read_line(&mut ndx).expect("Failed while reading index file");
                    if ndx.trim().is_empty() {
                        ndx = "?index.ndx".to_string();
                    }
                    ndx = convert_cur_dir(&ndx, &settings);
                    ndx = confirm_file_validity(&mut ndx, vec!["ndx"], &settings);
                }
                -10 => break,
                _ => println!("Error input.")
            };
        }
    } else { // pdbqt
        // FUCK!!! I will combine it with above sooner or later
        let receptor = PDBQT::new(&infile);
        println!("\nFinished loading receptor pdbqt.");
        println!("{}", receptor);
        let mut ligand = String::new();

        loop {
            println!("\n                 ************ MM/PB-SA Files ************");
            println!("-10 Exit program");
            println!(" -2 Toggle whether debug mode, current: {}", settings.debug_mode);
            println!(" -1 Set apbs path, current: {}", match &settings.apbs_path {
                Some(s) => s.to_string(),
                None => String::from("Not set")
            });
            println!("  0 Go to next step");
            println!("  1 Current receptor file: {}", infile);
            println!("  2 Assign ligand file (pdbqt), current: {}", match ligand.len() {
                0 => "undefined",
                _ => ligand.as_str()
            });
            let i = get_input_selection();
            match i {
                -2 => settings.debug_mode = !settings.debug_mode,
                -1 => {
                    println!("Input APBS path (if empty, means not to do PBSA calculation):");
                    let s: String = get_input_selection();
                    match check_apbs(&Some(s)) {
                        Some(s) => settings.apbs_path = Some(s),
                        None => settings.apbs_path = None
                    }
                }
                0 => {
                    if ligand.len() == 0 {
                        println!("Ligand file not assigned.");
                    } else {
                        // go to next step
                        let ligand = PDBQT::new(&ligand);
                        set_para_dock(&receptor, &ligand, &wd, settings);
                    }
                }
                2 => {
                    println!("Input ligand file path, default: ?DSDP_out.pdbqt (\"?\" means the same directory as receptor):");
                    ligand.clear();
                    stdin().read_line(&mut ligand).expect("Failed while reading ligand file");
                    if ligand.trim().is_empty() {
                        ligand = "?DSDP_out.pdbqt".to_string();
                    }
                    ligand = convert_cur_dir(&ligand, &settings);
                    ligand = confirm_file_validity(&mut ligand, vec!["pdbqt"], &settings);
                }
                -10 => break,
                _ => println!("Error input.")
            };
        }
    }
}