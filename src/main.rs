mod index_parser;
mod mmpbsa;
mod parse_tpr;
mod analyzation;
mod fun_para_basic;
mod fun_para_trj;
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
use crate::parse_tpr::TPR;
use settings::{Settings, get_base_settings};
use crate::settings::init_settings;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut tpr_dump = String::new();        // may be dump file
    let mut trj = String::from("");
    let mut ndx = String::from("");

    welcome();
    // initialize parameters
    let mut settings = init_settings();
    let programs = check_basic_programs(settings.gmx, settings.apbs);
    settings.gmx = programs.0;
    settings.apbs = programs.1;

    match args.len() {
        1 => {
            println!("Input path of .tpr or .dump file, e.g. D:/Study/ZhangYang.tpr or D:/Study/ZhangYang.dump");
            println!("Hint: input \"o\" to simply load last-opened .tpr or .dump file");
            loop {
                stdin().read_line(&mut tpr_dump).expect("Failed to read tpr or dumped file.");
                if tpr_dump.trim() == "o" {
                    tpr_dump = settings.last_opened.to_string();
                    if tpr_dump.len() == 0 {
                        println!("Last-opened tpr or mdp not found.");
                    }
                }
                if !Path::new(tpr_dump.trim()).is_file() {
                    println!("Not file: {}, input again:", tpr_dump.trim());
                    tpr_dump.clear();
                } else {
                    break;
                }
            }
        }
        2 => tpr_dump = args[1].to_string(),
        _ => {
            for i in (1..args.len()).step_by(2) {
                match args[i].as_str() {
                    "-f" => { trj = args[i + 1].to_string() }
                    "-s" => { tpr_dump = args[i + 1].to_string() }
                    "-n" => { ndx = args[i + 1].to_string() }
                    _ => {
                        println!("Omitted invalid option: {}", args[i])
                    }
                }
            }
        }
    }
    tpr_dump = confirm_file_validity(&mut tpr_dump, vec!["tpr", "dump"], &settings);

    settings.last_opened = fs::canonicalize(Path::new(&tpr_dump))
        .expect("Cannot convert to absolute path.").display().to_string();
    change_settings_last_opened(&tpr_dump);

    // get mdp or dump tpr
    let tpr_dump_path = fs::canonicalize(Path::new(&tpr_dump)).expect("Cannot get absolute tpr path.");
    let tpr_dump_name = tpr_dump_path.file_stem().unwrap().to_str().unwrap();
    let tpr_dir = tpr_dump_path.parent().expect("Failed to get tpr parent path");
    let dump_path = tpr_dir.join(tpr_dump_name.to_string() + ".dump");
    println!("Currently working at path: {}", Path::new(&tpr_dir).display());

    // It names tpr but exactly dump file _(:qゝ∠)_
    let tpr_name = match tpr_dump.ends_with(".tpr") {
        true => {
            println!("Found tpr file: {}", tpr_dump);
            let gmx = settings.gmx.as_ref()
                .expect("Gromacs not configured correctly.").as_str();
            let dump_to = dump_path.to_str().unwrap().to_string();
            dump_tpr(&tpr_dump, &dump_to, gmx);
            dump_to
        }
        false => {
            println!("Found dump file: {}", tpr_dump);
            tpr_dump.to_string()
        }
    };

    let mut tpr = TPR::new(&tpr_name, &settings);
    println!("\nFinished loading tpr.");

    // go to next step
    fun_para_basic::set_para_basic(&trj, &mut tpr, &ndx, &tpr_dir, tpr_name.as_str(), &mut settings);
}

fn welcome() {
    println!("\
        ========================================================================\n\
        | s_mmpbsa: Supernova's tool of calculating binding free energy using  |\n\
        | molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method |\n\
        ========================================================================\n\
        Website: https://github.com/supernova4869/s_mmpbsa\n\
        Developed by Jiaxing Zhang (zhangjiaxing7137@tju.edu.cn), Tian Jin University.\n\
        Version 0.2, first release: 2022-Oct-17, current version: 2024-Apr-5\n");
    println!("Usage 1: run `s_mmpbsa` and follow the prompts.\n\
        Usage 2: run `s_mmpbsa WangBingBing.tpr` to directly load tpr file.\n\
        Usage 3: run `s_mmpbsa WangBingBing.dump` to directly load dumped tpr file.\n\
        Usage 4: run `s_mmpbsa -f md.xtc -s md.tpr -n index.ndx` to assign all files.\n\
        Usage 5: run `s_mmpbsa -f md.xtc -s md.dump -n index.ndx` to assign all files.\n");
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

fn check_basic_programs(gmx: Option<String>, apbs: Option<String>) -> (Option<String>, Option<String>) {
    (check_gromacs(gmx), check_apbs(apbs))
}

fn check_gromacs(gmx: Option<String>) -> Option<String> {
    match gmx {
        Some(gmx) => {
            let gmx = match gmx.as_str() {
                "built-in" => {
                    get_built_in_gmx()
                }
                _ => Some(gmx)
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

fn check_apbs(apbs: Option<String>) -> Option<String> {
    match apbs {
        Some(apbs) => {
            let apbs = match apbs.as_str() {
                "built-in" => {
                    get_built_in_apbs()
                }
                _ => Some(apbs)
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

fn check_program_validity(program: &str) -> Result<String, ()> {
    let output = Command::new(program).arg("--version").output();
    match output {
        Ok(output) => {
            match output.status.code() {
                Some(0) => Ok(program.to_string()),
                Some(13) => Ok(program.to_string()),    // Fuck APBS cannot return 0 without input
                _ => Err(())
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

fn dump_tpr(tpr: &String, dump_to: &String, gmx: &str) {
    let tpr_dump = Command::new(gmx).arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(dump_to).expect("Cannot create md.dump.");
    outfile.write(tpr_dump.as_bytes()).expect("Cannot write md.dump.");
    println!("Dumped tpr file to {}", dump_to);
}