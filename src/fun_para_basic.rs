use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};
use std::io::stdin;
use std::path::Path;
use crate::settings::Settings;
use crate::utils::{append_new_name, get_input_selection, make_ndx};
use crate::{confirm_file_validity, convert_cur_dir, set_program};
use crate::fun_para_system::{set_para_trj, set_para_trj_pdbqt};
use crate::parse_tpr::TPR;

pub fn set_para_basic_tpr(tpr_path: &String, wd: &Path, settings: &mut Settings) {
    let mut trj = String::new();
    let mut ndx = String::new();
    let mut tpr = TPR::from(&tpr_path, &settings);
    println!("\nFinished loading input file.");

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
                match set_program(&Some(s), "apbs") {
                    Some(s) => settings.apbs_path = Some(s),
                    None => settings.apbs_path = None
                }
            }
            -4 => {
                println!("Input Delphi path:");
                let s: String = get_input_selection();
                match set_program(&Some(s), "delphi") {
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
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr", "pdbqt"], tpr_path);
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
                    make_ndx(&vec!["q"], wd, settings, &tpr_path, "", &ndx);
                }
                ndx = confirm_file_validity(&mut ndx, vec!["ndx", "pdbqt"], tpr_path);
            }
            -10 => break,
            _ => println!("Error input.")
        };
    }
}

pub fn set_para_basic_pdbqt(init_receptor_path: &String, wd: &Path, settings: &mut Settings) {
    let mut receptor_path = String::from(init_receptor_path);
    let mut ligand_path = String::new();
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
        println!("  1 Assign docking receptor file, current: {}", match receptor_path.len() {
            0 => "undefined",
            _ => receptor_path.as_str()
        });
        println!("  2 Assign docking ligand file, current: {}", match ligand_path.len() {
            0 => "undefined",
            _ => ligand_path.as_str()
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
                match set_program(&Some(s), "apbs") {
                    Some(s) => settings.apbs_path = Some(s),
                    None => settings.apbs_path = None
                }
            }
            -4 => {
                println!("Input Delphi path:");
                let s: String = get_input_selection();
                match set_program(&Some(s), "delphi") {
                    Some(s) => settings.delphi_path = Some(s),
                    None => settings.delphi_path = None
                }
            }
            0 => {
                if ligand_path.len() == 0 {
                    println!("Ligand file not assigned.");
                } else if receptor_path.len() == 0 {
                    println!("Receptor file not assigned.");
                } else {
                    // prepare pdbqt files
                    let complex_path = pdbqt2pdb(&receptor_path, &ligand_path, wd, settings);
                    // prepare index file
                    make_ndx(&vec!["q"], wd, settings, &complex_path, "", "MMPBSA_index.ndx");
                    // go to next step
                    set_para_trj_pdbqt(&complex_path, &wd, settings);
                }
            }
            1 => {
                println!("Input docking receptor file path, default: ?receptor.pdbqt (\"?\" means the same directory as initial input):");
                receptor_path.clear();
                stdin().read_line(&mut receptor_path).expect("Failed while reading receptor file");
                if receptor_path.trim().is_empty() {
                    receptor_path = "?receptor.pdbqt".to_string();
                }
                receptor_path = convert_cur_dir(&receptor_path, &init_receptor_path);
                receptor_path = confirm_file_validity(&mut receptor_path, vec!["pdbqt"], &init_receptor_path);
            }
            2 => {
                println!("Input docking ligand file path, default: ?ligand.pdbqt (\"?\" means the same directory as initial input):");
                ligand_path.clear();
                stdin().read_line(&mut ligand_path).expect("Failed while reading ligand file");
                if ligand_path.trim().is_empty() {
                    ligand_path = "?DSDP_out.pdbqt".to_string();
                }
                ligand_path = convert_cur_dir(&ligand_path, &init_receptor_path);
                ligand_path = confirm_file_validity(&mut ligand_path, vec!["pdbqt"], &init_receptor_path);
            }
            -10 => break,
            _ => println!("Error input.")
        };
    }
}

fn pdbqt2pdb(receptor_path: &String, ligand_path: &String, wd: &Path, settings: &Settings) -> String {
    let receptor_file_path = Path::new(receptor_path);
    let receptor_file_stem = receptor_file_path.file_stem().unwrap().to_str().unwrap();
    let ligand_file_path = Path::new(ligand_path);
    let ligand_file_stem = ligand_file_path.file_stem().unwrap().to_str().unwrap();
    let out_file_stem = format!("{}_{}", receptor_file_stem, ligand_file_stem);
    let out_file_name = append_new_name(&out_file_stem, ".pdb", "MMPBSA_");
    let out_file_path = wd.join(out_file_name);
    let out_file_path = out_file_path.to_str().unwrap();
    let mut pml_file = fs::File::create(wd.join(String::from(&out_file_stem) + ".pml")).unwrap();
    writeln!(pml_file, "cmd.load(r\"{}\", \"Protein\")", receptor_file_path.to_str().unwrap()).unwrap();
    writeln!(pml_file, "cmd.load(r\"{}\", \"Ligand\")", ligand_path).unwrap();
    writeln!(pml_file, "cmd.h_add(\"all\")").unwrap();
    writeln!(pml_file, "cmd.save(r\"{}\", selection=\"(all)\", state=0)", out_file_path).unwrap();
    writeln!(pml_file, "quit").unwrap();
    println!("\nLoading docking results files with PyMOL...");
    let result = Command::new(settings.pymol_path.as_ref().unwrap())
        .args(vec!["-cq", wd.join(String::from(out_file_stem).to_string() + ".pml").as_os_str().to_str().unwrap()])
        .stdout(Stdio::null())
        .spawn();
    match result {
        Ok(mut child) => {
            child.wait().ok();
        }
        Err(_) => {
            eprintln!("The configured PyMOL '{}' not found.", settings.pymol_path.as_ref().unwrap());
        }
    }
    println!("Finished loading docking results files.");
    return out_file_path.to_string()
}