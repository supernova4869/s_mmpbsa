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
mod parameters;
mod atom_property;

use std::fs;
use std::env;
use std::fs::File;
use std::io::{stdin, Write};
use std::path::Path;
use std::process::Command;
use std::str::FromStr;
use regex::Regex;
use chrono::prelude::Local;
use crate::parse_tpr::TPR;
use parameters::Parameters;
use crate::parameters::init_settings;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut tpr = String::new();        // may also be dump file
    let mut trj = String::from("");
    let mut ndx = String::from("");

    welcome();
    // initialize parameters
    let mut settings = init_settings();
    let programs = check_basic_programs(settings.gmx.as_str(), settings.apbs.as_str());
    settings.gmx = programs.0;
    settings.apbs = programs.1;

    match args.len() {
        1 => {
            println!("Input path of .tpr or .dump file, e.g. D:/Study/ZhangYang.tpr or D:/Study/ZhangYang.dump");
            println!("Hint: input \"o\" to simply load last-opened .tpr or .dump file");
            loop {
                stdin().read_line(&mut tpr).expect("Failed to read tpr or dumped file.");
                if tpr.trim() == "o" {
                    tpr = settings.last_opened.to_string();
                    if tpr.len() == 0 {
                        println!("Last-opened tpr or mdp not found.");
                    }
                }
                if !Path::new(tpr.trim()).is_file() {
                    println!("Not file: {}, input again:", tpr.trim());
                    tpr.clear();
                } else {
                    break;
                }
            }
        }
        2 => tpr = args[1].to_string(),
        _ => {
            for i in (1..args.len()).step_by(2) {
                match args[i].as_str() {
                    "-f" => { trj = args[i + 1].to_string() }
                    "-s" => { tpr = args[i + 1].to_string() }
                    "-n" => { ndx = args[i + 1].to_string() }
                    _ => {
                        println!("Omitted invalid option: {}", args[i])
                    }
                }
            }
        }
    }
    tpr = confirm_file_validity(&mut tpr, vec!["tpr", "dump"], &settings);

    change_settings_last_opened(&tpr);
    settings.last_opened = tpr.to_string();

    // working directory (path of tpr location)
    let wd = Path::new(&tpr).parent().unwrap();
    println!("Currently working at path: {}", fs::canonicalize(Path::new(&tpr)).unwrap().as_path().parent().unwrap().display());
    // get mdp or dump tpr
    let tpr = match tpr.ends_with(".tpr") {
        true => {
            println!("Found tpr file: {}", tpr);
            let p = tpr[..tpr.len() - 4].to_string() + ".dump";
            dump_tpr(&tpr, &p, settings.gmx.as_str());
            p
        }
        false => {
            println!("Found dump file: {}", tpr);
            tpr.to_string()
        }
    };

    let mut tpr = TPR::new(tpr.as_str(), &settings);
    println!("\nFinished reading tpr.");

    let mut atm_radius: Vec<f64> = Vec::new();
    for mol in &tpr.molecules {
        for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
            for atom in &mol.atoms {
                atm_radius.push(atom.radius);
            }
        }
    }
    let atom_radius = atom_radius::Radius::new(0, atm_radius);

    // go to next step
    fun_para_basic::set_para_basic(&trj, &mut tpr, &ndx, wd, &atom_radius, &mut settings);
}

fn welcome() {
    println!("\
        ========================================================================\n\
        | super_mmpbsa: Supernova's tool of calculating binding free energy by |\n\
        | molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method |\n\
        ========================================================================\n\
        Website: https://github.com/supernovaZhangJiaXing/super_mmpbsa\n\
        Developed by Jiaxing Zhang (zhangjiaxing7137@tju.edu.cn), Tian Jin University.\n\
        Version 0.1, first release: 2022-Oct-17\n\
        Current time: {}\n\n\
        Usage 1: run `super_mmpbsa` and follow the prompts.\n\
        Usage 2: run `super_mmpbsa WangBingBing.tpr` to directly load tpr file.\n\
        Usage 3: run `super_mmpbsa WangBingBing.dump` to directly load dumped tpr file.\n\
        Usage 4: run `super_mmpbsa -f md.xtc -s md.tpr -n index.ndx` to assign all files.\n\
        Usage 5: run `super_mmpbsa -f md.xtc -s md.dump -n index.ndx` to assign all files.\n",
             Local::now().format("%Y-%m-%d %H:%M:%S").to_string());
}

pub fn get_input_selection<T: FromStr>() -> T {
    loop {
        let mut input = String::from("");
        stdin().read_line(&mut input).expect("Error input.");
        match input.trim().parse() {
            Ok(num) => return num,
            Err(_) => {
                println!("Error input, input again.");
                continue;
            }
        };
    }
}

pub fn get_input_value(default: &String) -> String {
    let mut input = String::from("");
    stdin().read_line(&mut input).expect("Error input.");
    if input.is_empty() {
        default.to_string()
    } else {
        input.trim().to_string()
    }
}

// 把ext_list改成enum
pub fn confirm_file_validity(file_name: &String, ext_list: Vec<&str>, settings: &Parameters) -> String {
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
        let file_ext = Path::new(&f_name).extension().unwrap().to_str().unwrap();
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
    if cfg!(windows) {
        env::current_exe().unwrap().parent().unwrap().join("programs").join("gmx")
            .join("win").join("gmx.exe").to_str().unwrap().to_string()
    } else {
        println!("Currently not supported.");
        String::new()
    }
}

fn check_basic_programs(gmx: &str, apbs: &str) -> (String, String) {
    let mut gmx_path: String = String::new();
    let mut apbs_path: String = String::new();
    let mut gmx = gmx.to_string();
    if gmx == "built-in" {
        gmx = get_built_in_gmx();
    }
    match check_program_validity(gmx.as_str()) {
        Ok(p) => {
            gmx_path = p;
            println!("Using Gromacs: {}", gmx);
        }
        Err(_) => {
            println!("Note: Gromacs not configured correctly: {}. Now trying default gmx.", gmx);
            match check_program_validity("gmx") {
                Ok(p) => {
                    gmx_path = p;
                    println!("Using Gromacs: gmx");
                }
                Err(_) => {
                    println!("Note: default gmx invalid. Now trying built-in gmx of super_mmpbsa.");
                    match check_program_validity(get_built_in_gmx().as_str()) {
                        Ok(p) => {
                            gmx_path = p;
                            println!("Using Gromacs: {}", get_built_in_gmx().as_str());
                        }
                        Err(_) => {
                            println!("Warning: no valid Gromacs program in use.");
                        }
                    }
                }
            }
        }
    }

    // apbs
    match check_program_validity(apbs) {
        Ok(p) => {
            apbs_path = p;
            println!("Using APBS: {}", apbs);
            fs::remove_file(Path::new("io.mc")).unwrap();
        }
        Err(_) => {
            println!("Note: APBS not configured correctly. Now trying default apbs.");
            match check_program_validity("apbs") {
                Ok(p) => {
                    apbs_path = p;
                    println!("Using APBS: apbs");
                    fs::remove_file(Path::new("io.mc")).unwrap();
                }
                Err(_) => {
                    println!("Warning: no valid APBS program in use.");
                }
            }
        }
    }

    return (gmx_path, apbs_path);
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

pub fn convert_cur_dir(p: &String, settings: &Parameters) -> String {
    if p.starts_with('?') {
        let last_opened = &settings.last_opened;
        if last_opened.len() != 0 {
            Path::new(last_opened).parent().unwrap().join(p[1..].to_string()).to_str().unwrap().to_string()
        } else {
            p.to_string()
        }
    } else {
        p.to_string()
    }
}

fn change_settings_last_opened(tpr_mdp: &String) {
    // change settings.ini last opened file
    let settings_file = env::current_exe().unwrap().parent().unwrap().join("settings.ini");
    if settings_file.is_file() {
        let settings = fs::read_to_string(&settings_file).unwrap();
        let re = Regex::new("last_opened.*\".*\"").unwrap();
        let last_opened = fs::canonicalize(Path::new(&tpr_mdp)).unwrap().display().to_string();
        let settings = re.replace(settings.as_str(), format!("last_opened = \"{}\"",
                                                             &last_opened));
        let mut settings_file = File::create(settings_file).unwrap();
        settings_file.write_all(settings.as_bytes()).unwrap();
    }
}

fn dump_tpr(tpr: &String, dump_to: &String, gmx: &str) {
    let tpr_dump = Command::new(gmx).arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(dump_to).unwrap();
    outfile.write(tpr_dump.as_bytes()).unwrap();
    println!("Dumped tpr file to {}", fs::canonicalize(dump_to).unwrap().display());
}