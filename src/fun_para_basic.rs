use std::io::stdin;
use std::path::Path;
use crate::settings::Settings;
use crate::utils::{append_new_name, get_input_selection, make_ndx};
use crate::{check_apbs, check_delphi, confirm_file_validity, convert_cur_dir};
use crate::fun_para_system::set_para_trj;
use crate::parse_tpr::TPR;

pub fn set_para_basic(tpr_path: &String, wd: &Path, settings: &mut Settings) {
    let mut trj = String::new();
    let mut ndx = String::new();
    let mut tpr = TPR::new(&tpr_path, &settings);
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
                    println!("Index file not assigned.");
                } else {
                    // go to next step
                    set_para_trj(&trj, &mut tpr, &ndx, &wd, &tpr_path, settings);
                }
            }
            1 => {
                println!("Input trajectory file path, default: ?md.xtc (\"?\" means the same directory as tpr):");
                trj.clear();
                stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                if trj.trim().is_empty() {
                    trj = "?md.xtc".to_string();
                }
                trj = convert_cur_dir(&trj, tpr_path);
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr"], tpr_path);
            }
            2 => {
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
                    make_ndx("q", wd, settings, &tpr_path, &ndx);
                }
                ndx = confirm_file_validity(&mut ndx, vec!["ndx"], tpr_path);
            }
            -10 => break,
            _ => println!("Error input.")
        };
    }
}