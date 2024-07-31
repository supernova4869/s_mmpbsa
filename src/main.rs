mod index_parser;
mod mmpbsa;
mod parse_tpr;
mod parse_pdbqt;
mod analyzation;
mod fun_para_basic;
mod fun_para_system;
mod fun_para_mmpbsa;
mod atom_radius;
mod apbs_param;
mod prepare_apbs;
mod settings;
mod atom_property;
mod coefficients;
mod utils;

use std::fs;
use std::env;
use std::fs::File;
use std::io::{stdin, Write};
use std::path::Path;
use std::process::Command;
use regex::Regex;
use settings::{Settings, get_base_settings, get_settings_in_use};

fn main() {
    welcome();
    // initialize parameters
    let mut settings = match get_settings_in_use() {
        Some(settings_file) => {
            Settings::from(&settings_file)
        }
        None => {
            println!("Note: settings.ini not found, will use 1 kernel.");
            Settings::new()
        }
    };
    let programs = check_basic_programs(&settings);
    settings.gmx_path = programs.0;
    settings.apbs_path = programs.1;
    settings.delphi_path = programs.2;

    let args: Vec<String> = env::args().collect();
    let mut infile: String = String::new();
    match args.len() {
        1 => {
            println!("Input path of tpr, dump or pdbqt file, e.g. D:/Shinichi/Shiho.tpr or D:/Conan/Ai.pdbqt");
            println!("Hint: input \"o\" to simply load last-opened file");
            loop {
                stdin().read_line(&mut infile).expect("Failed to get input file.");
                if infile.trim() == "o" {
                    infile = settings.last_opened.to_string();
                    if infile.len() == 0 {
                        println!("Last-opened tpr or mdp not found.");
                    }
                }
                if !Path::new(infile.trim()).is_file() {
                    println!("Not file: {}, input again:", infile.trim());
                    infile.clear();
                } else {
                    break;
                }
            }
        }
        2 => {
            infile = args[1].to_string()
        }
        _ => {}
    }

    infile = confirm_file_validity(&infile, vec!["tpr", "dump", "pdbqt"], &settings);

    settings.last_opened = fs::canonicalize(Path::new(&infile))
        .expect("Cannot convert to absolute path.").display().to_string();
    change_settings_last_opened(&infile);

    // dump tpr and do nothing with pdbqt
    if infile.ends_with("tpr") || infile.ends_with("dump") {
        infile = get_dump(&infile, &settings);
    }

    match settings.debug_mode {
        true => println!("Debug mode open."),
        false => println!("Debug mode closed."),
    }

    // go to next step
    fun_para_basic::set_para_basic(&infile, &Path::new(&infile).parent().unwrap(), &mut settings);
}

fn get_dump(infile: &String, settings: &Settings) -> String {
    // get dump file or dumpped tpr
    let tpr_dump_path = fs::canonicalize(Path::new(&infile)).expect("Cannot get absolute tpr path.");
    let tpr_dump_name = tpr_dump_path.file_stem().unwrap().to_str().unwrap();
    let tpr_dir = tpr_dump_path.parent().expect("Failed to get tpr parent path");
    let dump_path = tpr_dir.join(tpr_dump_name.to_string() + ".dump");
    println!("Currently working at path: {}", Path::new(&tpr_dir).display());

    // It names tpr but exactly dump file _(:qゝ∠)_
    let tpr_name = match infile.ends_with(".tpr") {
        true => {
            println!("Found tpr file: {}", infile);
            let gmx = settings.gmx_path.as_ref().unwrap();
            let dump_to = dump_path.to_str().unwrap().to_string();
            dump_tpr(&infile, &dump_to, gmx);
            dump_to
        }
        false => {
            println!("Found dump file: {}", infile);
            infile.to_string()
        }
    };
    tpr_name
}

fn welcome() {
    println!("\
        ========================================================================\n\
        | s_mmpbsa: Supernova's tool of calculating binding free energy using  |\n\
        | molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method |\n\
        ========================================================================\n\
        Website: https://github.com/supernova4869/s_mmpbsa\n\
        Developed by Jiaxing Zhang (zhangjiaxing7137@tju.edu.cn), Tianjin University.\n\
        Version 0.4, first release: 2022-Oct-17, current release: 2024-Jul-25\n");
    println!("Usage 1: run `s_mmpbsa` and follow the prompts.\n\
        Usage 2: run `s_mmpbsa Miyano_Shiho.tpr` to load tpr file.\n\
        Usage 3: run `s_mmpbsa Miyano_Shiho.dump` to load dumped tpr file.\n\
        Usage 4: run `s_mmpbsa Haibara_Ai.pdbqt` to load receptor pdbqt file.\n");
}

// 把ext_list改成enum
pub fn confirm_file_validity(file_name: &String, ext_list: Vec<&str>, settings: &Settings) -> String {
    let mut f_name = String::from(file_name);
    loop {
        f_name = convert_cur_dir(&f_name, &settings).trim().to_string();
        if !Path::new(&f_name).is_file() {
            println!("Not valid file: {}. Input file path again.", f_name);
            f_name.clear();
            stdin().read_line(&mut f_name).expect("Failed to read file name.");
            continue;
        }
        // check extension
        let file_ext = Path::new(&f_name).extension()
            .expect("Input file has no extension.")
            .to_str().expect("Extension not valid Unicode.");
        for i in 0..ext_list.len() {
            if file_ext != ext_list[i] {
                continue;
            } else {
                return f_name.trim().to_string();
            }
        }
        println!("Not valid {:?} file, currently {}. Input file path again.", ext_list, file_ext);
        f_name.clear();
        stdin().read_line(&mut f_name).expect("Failed to read file name.");
    }
}

fn get_built_in_gmx() -> Option<String> {
    if cfg!(windows) {
        Some(env::current_exe().expect("Cannot get current s_mmpbsa program path.")
            .parent().expect("Cannot get current s_mmpbsa program directory.")
            .join("programs").join("gmx")
            .join("win").join("gmx.exe").to_str()
            .expect("The built-in gromacs not found.").to_string())
    } else {
        println!("Built-in gromacs not supported on Linux.");
        None
    }
}

fn get_built_in_apbs() -> Option<String> {
    if cfg!(windows) {
        Some(env::current_exe().expect("Cannot get current s_mmpbsa program path.")
            .parent()
            .expect("Cannot get current s_mmpbsa program directory.")
            .join("programs").join("apbs")
            .join("win").join("apbs.exe").to_str()
            .expect("The built-in apbs not found.").to_string())
    } else if cfg!(unix) {
        Some(env::current_exe().expect("Cannot get current s_mmpbsa program path.")
            .parent()
            .expect("Cannot get current s_mmpbsa program directory.")
            .join("programs").join("apbs")
            .join("linux").join("apbs").to_str()
            .expect("The built-in apbs not found.").to_string())
    } else {
        println!("Built-in apbs not found for current system.");
        None
    }
}

fn get_built_in_delphi() -> Option<String> {
    if cfg!(windows) {
        Some(env::current_exe().expect("Cannot get current s_mmpbsa program path.")
            .parent()
            .expect("Cannot get current s_mmpbsa program directory.")
            .join("programs").join("delphi")
            .join("win").join("delphi.exe").to_str()
            .expect("The built-in delphi not found.").to_string())
    } else if cfg!(unix) {
        Some(env::current_exe().expect("Cannot get current s_mmpbsa program path.")
            .parent()
            .expect("Cannot get current s_mmpbsa program directory.")
            .join("programs").join("delphi")
            .join("linux").join("delphi").to_str()
            .expect("The built-in delphi not found.").to_string())
    } else {
        println!("Built-in delphi not found for current system.");
        None
    }
}

fn check_basic_programs(settings: &Settings) -> (Option<String>, Option<String>, Option<String>) {
    (check_gromacs(&settings.gmx_path), check_apbs(&settings.apbs_path), check_delphi(&settings.delphi_path))
}

fn check_gromacs(gmx: &Option<String>) -> Option<String> {
    match gmx {
        Some(gmx) => {
            let gmx = match gmx.as_str() {
                "built-in" => {
                    get_built_in_gmx()
                }
                _ => Some(gmx.to_string())
            };
            match gmx {
                Some(gmx) => {
                    match check_program_validity(gmx.as_str()) {
                        Ok(p) => {
                            println!("Using Gromacs: {}", gmx);
                            Some(p)
                        }
                        Err(_) => {
                            println!("Note: {} not valid. Will try built-in gmx.", gmx);
                            match get_built_in_gmx() {
                                Some(gmx) => {
                                    match check_program_validity(gmx.as_str()) {
                                        Ok(p) => {
                                            println!("Using Gromacs: {}", gmx);
                                            Some(p)
                                        }
                                        Err(_) => {
                                            println!("Warning: no valid Gromacs program in use.");
                                            None
                                        }
                                    }
                                }
                                _ => {
                                    println!("Warning: no valid Gromacs program in use.");
                                    None
                                }
                            }
                        }
                    }
                }
                None => None
            }
        }
        _ => {
            println!("Warning: no valid Gromacs program in use.");
            None
        }
    }
}

fn check_apbs(apbs: &Option<String>) -> Option<String> {
    match apbs {
        Some(apbs) => {
            let apbs = match apbs.as_str() {
                "built-in" => {
                    get_built_in_apbs()
                }
                _ => Some(apbs.to_string())
            };
            match apbs {
                Some(apbs) => {
                    match check_program_validity(apbs.as_str()) {
                        Ok(p) => {
                            println!("Using APBS: {}", apbs);
                            match fs::remove_file(Path::new("io.mc")) {
                                _ => Some(p)
                            }
                        }
                        Err(_) => {
                            println!("Warning: no valid APBS program in use.");
                            None
                        }
                    }
                }
                None => None
            }
        }
        None => None
    }
}

fn check_delphi(delphi: &Option<String>) -> Option<String> {
    match delphi {
        Some(delphi) => {
            let delphi = match delphi.as_str() {
                "built-in" => {
                    get_built_in_delphi()
                }
                _ => Some(delphi.to_string())
            };
            match delphi {
                Some(delphi) => {
                    match check_program_validity(delphi.as_str()) {
                        Ok(p) => {
                            println!("Using delphi: {}", delphi);
                            Some(p)
                        }
                        Err(_) => {
                            println!("Warning: no valid Delphi program in use.");
                            None
                        }
                    }
                }
                None => None
            }
        }
        None => None
    }
}

fn check_program_validity(program: &str) -> Result<String, ()> {
    let output = Command::new(program).arg("--version").output();
    match output {
        Ok(output) => {
            match output.status.code() {
                Some(0) => Ok(program.to_string()),
                Some(13) => Ok(program.to_string()),    // Fuck APBS cannot return 0 without input
                Some(1) => Ok(program.to_string()),    // Currently do not know delphi's test command
                _ => {
                    Err(())
                }
            }
        }
        Err(_) => Err(())
    }
}

pub fn convert_cur_dir(p: &String, settings: &Settings) -> String {
    if p.starts_with('?') {
        let last_opened = &settings.last_opened;
        if last_opened.len() != 0 {
            Path::new(last_opened).parent()
                .expect("Cannot get path of last-opened file.")
                .join(p[1..].to_string()).to_str()
                .expect("Path of last-opened file not valid unicode.").to_string()
        } else {
            p.to_string()
        }
    } else {
        p.to_string()
    }
}

fn change_settings_last_opened(tpr_mdp: &String) {
    if get_base_settings().is_file() {
        let base_settings = fs::read_to_string(get_base_settings())
            .expect("Cannot read settings.ini.");
        let re = Regex::new("last_opened.*\".*\"").unwrap();
        let last_opened = fs::canonicalize(Path::new(&tpr_mdp))
            .expect("Cannot convert to absolute path.").display().to_string();
        let settings = re.replace(base_settings.as_str(),
            format!("last_opened = \"{}\"", &last_opened));
        let mut settings_file = File::create(get_base_settings())
            .expect("Cannot edit settings.ini.");
        settings_file.write_all(settings.as_bytes()).expect("Cannot write to settings.ini.");
    }
}

fn dump_tpr(tpr: &String, dump_to: &String, gmx: &str) {
    let tpr_dump = Command::new(gmx).arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(dump_to).expect("Cannot create md.dump.");
    outfile.write(tpr_dump.as_bytes()).expect("Cannot write md.dump.");
    println!("Dumped tpr file to {}", dump_to);
}