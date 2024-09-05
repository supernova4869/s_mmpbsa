mod index_parser;
mod mmpbsa;
mod parse_tpr;
mod parse_xvg;
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
    welcome("2024-Sep-3");
    let mut settings = env_check();
    match settings.debug_mode {
        true => println!("Debug mode on.\n"),
        false => println!("Debug mode off.\n"),
    }

    let args: Vec<String> = env::args().collect();
    let mut input: String = String::new();
    match args.len() {
        1 => {
            println!("Input path of tpr file, e.g. D:/Conan/Ai.tpr");
            println!("Hint: input \"o\" to simply load last-opened file; input \"a\" to start analyzation.");
            stdin().read_line(&mut input).expect("Failed to get input file.");
            if input.trim().eq("o") {
                input = settings.last_opened.to_string();
                if input.len() == 0 {
                    println!("Last-opened tpr not found.");
                }
            } else if input.trim().eq("a") {
                input = Path::new(&settings.last_opened).parent().unwrap().to_str().unwrap().to_string();
                println!("Input path of working dir with .sm results (default: {}):", input);
                let temp = get_input("".to_string());
                if temp.len() != 0 {
                    input = temp;
                }
            } else {
                input = input.trim().to_string();
            }
        }
        2 => {
            input = args[1].to_string()
        }
        _ => {}
    }

    if Path::new(&input).is_file() {
        let tpr_file = confirm_file_validity(&input, vec!["tpr"], &input);
        change_settings_last_opened(&mut settings, &tpr_file);
        let tpr_file = get_dump(&tpr_file, &settings);
        fun_para_basic::set_para_basic(&tpr_file, &Path::new(&tpr_file).parent().unwrap(), &mut settings);
    } else if Path::new(&input).is_dir() {
        let wd = Path::new(&input);
        let sm_list: Vec<String> = fs::read_dir(wd).unwrap().into_iter().filter_map(|f| {
            let f = f.unwrap().path();
            if f.extension().unwrap().to_str().unwrap().eq("sm") {
                Some(f.to_str().unwrap().to_string())
            } else {
                None
            }
        }).collect();
        if !sm_list.is_empty() {
            println!("Please input MD temperature (default: 298.15):");
            let temperature = get_input(298.15);
            println!("Please input system name (default: _system):");
            let sys_name = get_input("_system".to_string());
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
    } else {
        println!("Input {} not file or directory. Please check.", Path::new(&input).to_str().unwrap());
    }
}

fn welcome(today: &str) {
    println!("\
        ========================================================================\n\
        | s_mmpbsa: Supernova's tool of calculating binding free energy using  |\n\
        | molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method |\n\
        ========================================================================\n\
        Website: https://github.com/supernova4869/s_mmpbsa\n\
        Developed by Supernova (zhangjiaxing7137@tju.edu.cn), Tianjin University.\n\
        Version 0.4, first release: 2022-Oct-17, current release: {}\n", today);
    println!("Usage 1: run `s_mmpbsa` and follow the prompts.\n\
        Usage 2: run `s_mmpbsa Miyano_Shiho.tpr` to load tpr file.\n");
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

fn set_program(p: &Option<String>, name: &str) -> Option<String> {
    if let Some(p) = p {
        let p = if p.eq("built-in") {
            match name {
                "gromacs" => get_built_in_gmx(),
                "apbs" => get_built_in_apbs(),
                "delphi" => get_built_in_delphi(),
                _ => None
            }
        } else {
            Some(p.to_string())
        };
        if let Some(p) = p {
            match check_program_validity(p.as_str()) {
                Ok(p) => {
                    println!("Using {}: {}", name, p);
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
            Settings::from(&settings_file)
        }
        None => {
            println!("Note: settings.ini not found, will use 1 kernel.");
            Settings::new()
        }
    };
    // check necessary dat path
    if !Path::new(env::current_exe().unwrap().parent().unwrap().join("dat/").as_path()).is_dir() {
        println!("Error: the dat/ folder with atom radius not found, please check and retry.");
        io::stdin().read_line(&mut String::new()).unwrap();
        std::process::exit(0);
    }
    settings.gmx_path = set_program(&settings.gmx_path, "gromacs");
    settings.apbs_path = set_program(&settings.apbs_path, "apbs");
    fs::remove_file("io.mc").ok();
    settings.delphi_path = set_program(&settings.delphi_path, "delphi");
    settings
}