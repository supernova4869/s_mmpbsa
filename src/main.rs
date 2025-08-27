mod index_parser;
mod mmpbsa;
mod parse_tpr;
mod parse_xvg;
mod parse_pdb;
mod parse_gro;
mod parse_pdbqt;
mod parse_mol2;
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

use std::{fs, io};
use std::env;
use std::fs::File;
use std::io::{stdin, Write};
use std::path::Path;
use std::process::{exit, Command};
use analyzation::SMResult;
use regex::Regex;
use settings::{Settings, get_base_settings, get_settings_in_use};
use utils::get_input;

fn main() {
    let version = 0.8;
    welcome(&version.to_string(), "2025-Aug-27");
    let mut settings = env_check();
    match settings.debug_mode {
        true => println!("Debug mode on.\n"),
        false => println!("Debug mode off.\n"),
    }

    let args: Vec<String> = env::args().collect();
    let mut input: String = String::new();
    let mut ligand = String::new();
    match args.len() {
        1 => {
            println!("Input path of tpr file, e.g. D:/md.tpr");
            println!("Or, input path of docking receptor file, e.g. D:/receptor.pdbqt");
            println!("Hint: input \"o\" to simply load last-opened file.");
            println!("Hint: input \"a\" to start analyzation mode.");
            stdin().read_line(&mut input).expect("Failed to get input file.");
            if input.trim().eq("o") {
                input = settings.last_opened.to_string();
                if input.len() == 0 {
                    println!("Last-opened tpr not found.");
                }
            } else if input.trim().eq("a") {
                input = env::current_dir().unwrap().to_str().unwrap().to_string();
                println!("Input path of working dir with .sm results (default: current dir):");
                let temp = get_input("".to_string());
                if temp.len() != 0 {
                    input = temp;
                }
            } else {
                input = input.trim().to_string();
            }
        }
        2 => {
            input = args[1].to_string();
        }
        3 => {
            input = args[1].to_string();
            ligand = args[2].to_string();
        }
        _ => {}
    }

    if Path::new(&input).is_file() {
        let in_file = confirm_file_validity(&input, vec!["tpr", "pdbqt"], &input);
        change_settings_last_opened(&mut settings, &in_file);
        if in_file.ends_with("tpr") {
            let in_file = get_dump(&in_file, &settings);
            fun_para_basic::set_para_basic_tpr(&in_file, &Path::new(&in_file).parent().unwrap(), &mut settings);
        } else { // pdbqt
            let wd = fs::canonicalize(Path::new(&in_file)).unwrap();
            fun_para_basic::set_para_basic_pdbqt(&in_file, &ligand, &wd.parent().unwrap(), &mut settings);
        };
    } else if Path::new(&input).is_dir() {
        let wd = Path::new(&input);
        let sm_list: Vec<String> = fs::read_dir(wd).unwrap().into_iter().filter_map(|f| {
            let f = f.unwrap().path();
            if let Some(ext) = f.extension() {
                if ext.to_str().unwrap().eq("sm") {
                    Some(f.to_str().unwrap().to_string())
                } else {
                    None
                }
            } else {
                None
            }
        }).collect();
        if !sm_list.is_empty() {
            println!("Please input MD temperature (default: 298.15):");
            let temperature = get_input(298.15);
            println!("Please input system name (default: system):");
            let sys_name = get_input("system".to_string());
            println!("Loading MM/PB-SA results...");
            let result_wt = sm_list.iter().find(|f| {
                let f_name = Path::new(f).file_name().unwrap().to_str().unwrap();
                f_name.starts_with(&format!("_MMPBSA_{}", sys_name)) && f_name.ends_with("_WT.sm")
            }).ok_or_else(|| {
                println!("The required _MMPBSA_{}_WT.sm file not found. Please check.", sys_name);
                exit(0);
            });
            let result_wt = SMResult::from(result_wt.unwrap());
            let result_as: Vec<SMResult> = sm_list.iter().filter(|&f| {
                let f_name = Path::new(f).file_name().unwrap().to_str().unwrap();
                f_name.starts_with(&format!("_MMPBSA_{}", sys_name)) && !f_name.ends_with("_WT.sm")
            }).map(|f| SMResult::from(f)).collect();
            analyzation::analyze_controller(&result_wt, &result_as, temperature, &sys_name, wd, &settings);
        } else {
            println!("There is no MM/PB-SA results at {}. Please run MM/PB-SA calculations first.", &input);
        }
    } else if input.eq("--version") {
        utils::show_famous_quotes();
    } else {
        println!("Input {} not file or directory. Please check.", Path::new(&input).to_str().unwrap());
    }
}

fn welcome(version: &str, today: &str) {
    println!("\
        ========================================================================\n\
        | s_mmpbsa: Supernova's tool of calculating binding free energy using  |\n\
        | molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method |\n\
        ========================================================================\n\
        Website: https://github.com/supernova4869/s_mmpbsa\n\
        Developed by Supernova (zhangjiaxing7137@tju.edu.cn), Tianjin University.\n\
        Version {}, first release: 2022-Oct-17, current release: {}\n", version, today);
    println!("Usage 1: run `s_mmpbsa` and follow the prompts.\n\
        Usage 2: run `s_mmpbsa Haibara_Ai.tpr` to load MD tpr file.\n\
        Usage 3: run `s_mmpbsa Miyano_Shiho.pdbqt Kudo_Shinichi.pdbqt` (receptor first) to load docking results.\n");
}

pub fn confirm_file_validity(file_name: &String, ext_list: Vec<&str>, tpr_path: &str) -> String {
    let mut f_name = String::from(file_name);
    loop {
        f_name = convert_cur_dir(&f_name, &tpr_path).trim().to_string();
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

fn get_built_in_gmx() -> String {
    env::current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent().expect("Cannot get current s_mmpbsa program directory.")
        .join("programs").join("gmx")
        .join(if cfg!(windows) {"win"} else {"linux"}).join("gmx")
        .display().to_string()
}

fn get_built_in_apbs() -> String {
    env::current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent()
        .expect("Cannot get current s_mmpbsa program directory.")
        .join("programs").join("apbs")
        .join(if cfg!(windows) {"win"} else {"linux"}).join("apbs")
        .display().to_string()
}

fn get_built_in_delphi() -> String {
    env::current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent()
        .expect("Cannot get current s_mmpbsa program directory.")
        .join("programs").join("delphi")
        .join(if cfg!(windows) {"win"} else {"linux"}).join("delphi")
        .display().to_string()
}

fn get_built_in_antechamber() -> String {
    env::current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent()
        .expect("Cannot get current s_mmpbsa program directory.")
        .join("programs").join("amber")
        .join(if cfg!(windows) {"win"} else {"linux"})
        .join("bin").join("antechamber")
        .display().to_string()
}

fn get_built_in_sobtop() -> String {
    env::current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent()
        .expect("Cannot get current s_mmpbsa program directory.")
        .join("programs").join("sobtop")
        .join("sobtop")
        .display().to_string()
}

fn get_built_in_obabel() -> String {
    env::current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent()
        .expect("Cannot get current s_mmpbsa program directory.")
        .join("programs").join("openbabel")
        .join("obabel")
        .display().to_string()
}

fn set_program(p: &Option<String>, name: &str, settings: &Settings) -> Option<String> {
    if let Some(p) = p {
        let p = if p.eq("built-in") {
            match name {
                "gromacs" => get_built_in_gmx(),
                "apbs" => get_built_in_apbs(),
                "delphi" => get_built_in_delphi(),
                "antechamber" => get_built_in_antechamber(),
                "sobtop" => get_built_in_sobtop(),
                "obabel" => get_built_in_obabel(),
                _ => String::from("")
            }
        } else {
            utils::get_program_path(p).unwrap()
        };
        if !p.is_empty() {
            if settings.debug_mode {
                println!("Checking {} validity...", name);
            }
            match check_program_validity(p.as_str()) {
                Ok(p) => {
                    println!("Using {}: {}", name, p);
                    if name.eq("apbs") {
                        if Path::new("io.mc").is_file() {
                            fs::remove_file("io.mc").ok();
                        }
                    }
                    Some(p)
                }
                Err(_) => {
                    println!("Warning: no valid {} program in use.", name);
                    None
                }
            }
        } else {
            None
        }
    }
    else {
        None
    }
}

fn check_program_validity(program: &str) -> Result<String, ()> {
    let version_arg = match program {
        "obabel" => "-V",
        _ => "--version"
    };
    let output = Command::new(program).arg(version_arg).output();
    match output {
        Ok(output) => {
            // println!("{}", output.status.code().unwrap());
            match output.status.code() {
                Some(0) => Ok(program.to_string()),
                Some(13) => Ok(program.to_string()),    // APBS
                Some(127) => Ok(program.to_string()),    // APBS
                Some(1) => Ok(program.to_string()),    // delphi
                Some(24) => Ok(program.to_string()),    // sobtop
                Some(69) => Ok(program.to_string()),    // ?
                _ => {
                    Err(())
                }
            }
        }
        Err(_) => Err(())
    }
}

pub fn convert_cur_dir(p: &String, tpr_path: &str) -> String {
    if p.starts_with('?') {
        Path::new(tpr_path).parent()
            .expect("Cannot get path of tpr file.")
            .join(p[1..].to_string()).to_str()
            .expect("Path of tpr file not valid unicode.").to_string()
    } else {
        p.to_string()
    }
}

fn change_settings_last_opened(settings: &mut Settings, tpr: &String) {
    settings.last_opened = fs::canonicalize(Path::new(&tpr))
        .expect("Cannot convert to absolute path.").display().to_string();
    if get_base_settings().is_file() {
        let base_settings = fs::read_to_string(get_base_settings())
            .expect("Cannot read settings.ini.");
        let re = Regex::new("last_opened.*\".*\"").unwrap();
        let settings = re.replace(base_settings.as_str(),
            format!("last_opened = \"{}\"", &settings.last_opened));
        let mut settings_file = File::create(get_base_settings())
            .expect("Cannot edit settings.ini.");
        settings_file.write_all(settings.as_bytes()).expect("Cannot write to settings.ini.");
    }
}

fn get_dump(tpr_path: &String, settings: &Settings) -> String {
    // get dumpped tpr
    let tpr_dump_path = fs::canonicalize(Path::new(&tpr_path)).expect("Cannot get absolute tpr path.");
    let tpr_dump_name = tpr_dump_path.file_stem().unwrap().to_str().unwrap();
    let tpr_dir = tpr_dump_path.parent().expect("Failed to get tpr parent path");
    let dump_path = tpr_dir.join(tpr_dump_name.to_string() + ".dump");
    println!("Currently working at path: {}", Path::new(&tpr_dir).display());
    let gmx = settings.gmx_path.as_ref().unwrap();
    let dump_to = dump_path.to_str().unwrap().to_string();
    dump_tpr(&tpr_path, &dump_to, gmx);
    dump_to
}

fn dump_tpr(tpr: &String, dump_to: &String, gmx: &str) {
    let tpr_dump = Command::new(gmx).arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(dump_to).expect("Cannot create md.dump.");
    outfile.write(tpr_dump.as_bytes()).expect("Cannot write md.dump.");
    println!("Dumped tpr file to {}", dump_to);
}

fn env_check() -> Settings {
    // initialize parameters
    let mut settings = match get_settings_in_use() {
        Some(settings_file) => {
            println!("Note: found settings.ini at {}.", settings_file.display());
            Settings::from(&settings_file)
        },
        None => {
            println!("Note: no settings.ini found.");
            Settings::new()
        }
    };

    // Set global parallel kernels
    rayon::ThreadPoolBuilder::new()
        .num_threads(settings.nkernels)
        .build_global()
        .unwrap();

    // check necessary dat path
    if !Path::new(env::current_exe().unwrap().parent().unwrap().join("dat/").as_path()).is_dir() {
        println!("Error: the dat/ folder with atom radius not found, please check and retry.");
        io::stdin().read_line(&mut String::new()).unwrap();
        std::process::exit(0);
    }
    settings.gmx_path = set_program(&settings.gmx_path, "gromacs", &settings);
    settings.apbs_path = set_program(&settings.apbs_path, "apbs", &settings);
    settings.delphi_path = set_program(&settings.delphi_path, "delphi", &settings);
    settings.antechamber_path = set_program(&settings.antechamber_path, "antechamber", &settings);
    settings.sobtop_path = set_program(&settings.sobtop_path, "sobtop", &settings);
    settings.obabel_path = set_program(&settings.obabel_path, "obabel", &settings);
    settings
}