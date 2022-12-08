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

use std::fs;
use std::env;
use std::fs::File;
use std::io::{stdin, Write};
use std::path::Path;
use std::process::Command;
use std::str::FromStr;
use regex::Regex;
use toml;
use toml::Value;
use chrono::prelude::Local;
use crate::parse_tpr::TPR;

pub struct Parameters {
    rad_type: i32,
    rad_lj0: f64,
    mesh_type: i32,
    grid_type: i32,
    use_dh: bool,
    use_ts: bool,
    cfac: f64,
    fadd: f64,
    df: f64,
    nkernels: i32,
    preserve: bool,
    gmx: String,
    apbs: String,
    last_opened: String,
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut tpr_mdp = String::new();
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
            println!("Input path of .tpr or dumped .mdp file, e.g. D:/Study/ZhangYang.tpr or D:/Study/ZhangYang_dump.mdp");
            println!("Hint: input \"o\" to simply load last-opened .tpr or dumped .mdp file");
            loop {
                stdin().read_line(&mut tpr_mdp).expect("Failed to read tpr or dumped file.");
                if tpr_mdp.trim() == "o" {
                    tpr_mdp = settings.last_opened.clone();
                    if tpr_mdp.len() == 0 {
                        println!("Last-opened tpr or mdp not found.");
                    }
                }
                if !Path::new(tpr_mdp.trim()).is_file() {
                    println!("Not file: {}, input again:", tpr_mdp.trim());
                    tpr_mdp.clear();
                } else {
                    break;
                }
            }
        }
        2 => tpr_mdp = args[1].to_string(),
        _ => {
            for i in 1..args.len() {
                match args[i].as_str() {
                    "-f" => { trj = args[i + 1].to_string() }
                    "-s" => { tpr_mdp = args[i + 1].to_string() }
                    "-n" => { ndx = args[i + 1].to_string() }
                    _ => {
                        if i % 2 == 1 {
                            println!("Omitted invalid option: {}", args[i])
                        }
                    }
                }
            }
        }
    }
    tpr_mdp = confirm_file_validity(&mut tpr_mdp, vec!["tpr", "mdp"], &settings);

    change_settings_last_opened(&tpr_mdp);
    settings.last_opened = tpr_mdp.to_string();

    // working directory (path of tpr location)
    let wd = Path::new(&tpr_mdp).parent().unwrap();
    println!("Currently working at path: {}", fs::canonicalize(Path::new(&tpr_mdp)).unwrap().as_path().parent().unwrap().display());
    // get mdp or dump tpr to mdp
    let mut mdp_path = String::from(&tpr_mdp).clone();
    if tpr_mdp.ends_with(".tpr") {
        mdp_path = tpr_mdp[..&tpr_mdp.len() - 4].to_string() + "_dumped.mdp";
        dump_tpr(&tpr_mdp, &mdp_path, settings.gmx.as_str());
    }
    let tpr = TPR::new(mdp_path.as_str());
    println!("\nFinished reading tpr file:\n{}.", tpr);

    // go to next step
    // 修改接口, 开始下面流程
    // fun_para_basic::set_para_basic(&mut trj, &mdp_path, &mut ndx, wd, &mut settings);
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
        Usage 3: run `super_mmpbsa WangBingBing_dumped.mdp` to directly load dump file.\n\
        Usage 4: run `super_mmpbsa -f md.xtc -s md.tpr -n index.ndx` to assign all needed files.\n\
        Usage 5: run `super_mmpbsa -f md.xtc -s md_dumped.mdp -n index.ndx` to assign all needed files.\n",
             Local::now().format("%Y-%m-%d %H:%M:%S").to_string());
}

pub fn get_input_value<T: FromStr>() -> T {
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

fn init_settings() -> Parameters {
    let mut params = Parameters {
        rad_type: 1,
        rad_lj0: 1.2,
        mesh_type: 0,
        grid_type: 1,
        use_dh: true,
        use_ts: true,
        cfac: 3.0,
        fadd: 10.0,
        df: 0.5,
        nkernels: 4,
        preserve: true,
        gmx: String::from("gmx"),
        apbs: String::from("apbs"),
        last_opened: String::new(),
    };
    if Path::new("settings.ini").is_file() {
        // Find settings locally
        let settings = fs::read_to_string("settings.ini").unwrap();
        let settings = Regex::new(r"\\").unwrap().replace_all(settings.as_str(), "/").to_string();
        let settings: Value = toml::from_str(settings.as_str()).expect("Error with settings.ini grammar");
        read_user_settings(&mut params, &settings);
        println!("Note: found settings.ini in the current path. Will use {} kernels.", params.nkernels);
    } else {
        let super_mmpbsa_path = env::current_exe().unwrap()
            .parent().unwrap().join("settings.ini");
        if super_mmpbsa_path.is_file() {
            let settings = fs::read_to_string(super_mmpbsa_path).unwrap();
            let settings = Regex::new(r"\\").unwrap()
                .replace_all(settings.as_str(), "/").to_string();
            let settings: Value = toml::from_str(settings.as_str()).unwrap();
            read_user_settings(&mut params, &settings);
            println!("Note: found settings.ini in super_mmpbsa directory. Will use {} kernels.", params.nkernels);
        } else {
            println!("Note: settings.ini not found. Will use 1 kernel.");
        }
    }
    println!("(Currently multi-threading not yet utilized)");
    return params;
}

fn read_user_settings(params: &mut Parameters, settings: &Value) {
    params.rad_type = parse_param(settings, "radType", params.rad_type);
    params.rad_lj0 = parse_param(settings, "radLJ0", params.rad_lj0);
    params.mesh_type = parse_param(settings, "meshType", params.mesh_type);
    params.grid_type = parse_param(settings, "gridType", params.grid_type);
    params.cfac = parse_param(settings, "cfac", params.cfac);
    params.fadd = parse_param(settings, "fadd", params.fadd);
    params.df = parse_param(settings, "df", params.df);
    params.nkernels = parse_param(settings, "nkernels", params.nkernels);
    // String type cannot move
    params.preserve = match settings.get("preserve").unwrap().to_string().as_str() {
        "y" => true,
        "Y" => true,
        _ => false
    };
    let mut gmx = settings.get("gmx").unwrap().to_string();
    gmx = gmx[1..gmx.len() - 1].to_string();
    if gmx.len() == 0 {
        gmx = "gmx".to_string();
    }
    params.gmx = gmx;
    let mut apbs = settings.get("apbs").unwrap().to_string();
    apbs = apbs[1..apbs.len() - 1].to_string();
    if apbs.len() == 0 {
        apbs = "apbs".to_string();
    }
    params.apbs = apbs;
    let mut last_opened = settings.get("last_opened").unwrap().to_string();
    last_opened = last_opened[1..last_opened.len() - 1].to_string();
    if last_opened.len() == 0 {
        last_opened = String::new();
    }
    params.last_opened = last_opened;
}

fn parse_param<T: FromStr>(settings: &Value, key: &str, default: T) -> T {
    match settings.get(key).unwrap().to_string().parse::<T>() {
        Ok(v) => v,
        Err(_) => default
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