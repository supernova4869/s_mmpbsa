mod index_parser;
mod mmpbsa;
mod parse_tpr;
mod analyzation;

use std::fs;
use std::env;
use std::fmt::{Display, format};
use std::fs::File;
use std::io::{stdin, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::rc::Rc;
use std::str::FromStr;
use regex::Regex;
use toml;
use toml::Value;
use xdrfile::{Frame, XTCTrajectory};
use chrono::prelude::Local;

pub struct Parameters {
    rad_type: i32,
    rad_lj0: f64,
    mesh_type: i32,
    grid_type: i32,
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
    let mut use_dh = true;
    let mut use_ts = true;

    // start workflow
    welcome();
    // initialize parameters
    let mut settings = init_settings();
    let programs = check_basic_programs(settings.gmx.as_str(), settings.apbs.as_str());
    settings.gmx = programs.0;
    settings.apbs = programs.1;
    if settings.apbs.len() == 0 {
        println!("APBS invalid. Press any key to exit.");
        let mut temp = String::new();
        stdin().read_line(&mut temp).unwrap();
        return;
    }

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

    // working directory (path of tpr location)
    let wd = Path::new(&tpr_mdp).parent().unwrap();
    println!("Currently working at path: {}", fs::canonicalize(Path::new(&tpr_mdp)).unwrap().as_path().parent().unwrap().display());
    // get mdp or dump tpr to mdp
    let mut mdp_path = String::from(&tpr_mdp).clone();
    if tpr_mdp.ends_with(".tpr") {
        mdp_path = tpr_mdp[0..tpr_mdp.len() - 4].to_string() + "_dumped.mdp";
        dump_tpr(&tpr_mdp, &mdp_path, settings.gmx.as_str());
    }
    loop {
        println!("\n                 ************ SuperMMPBSA functions ************");
        println!("-2 Toggle whether to use entropy contribution, current: {}", use_ts);
        println!("-1 Toggle whether to use Debye-Huckel shielding method, current: {}", use_dh);
        println!(" 0 Ready for MM-PBSA calculations");
        println!(" 1 Assign trajectory file (xtc or trr), current: {}", match trj.len() {
            0 => "undefined",
            _ => trj.as_str()
        });
        println!(" 2 Assign index file (ndx), current: {}", match ndx.len() {
            0 => "undefined",
            _ => ndx.as_str()
        });
        println!(" 3 Exit program");
        let i = get_input_value();
        match i {
            -2 => { use_ts = !use_ts; }
            -1 => { use_dh = !use_dh; }
            0 => {
                if trj.len() == 0 {
                    println!("Trajectory file not assigned.");
                } else if ndx.len() == 0 {
                    // 可能要改, 以后不需要index也能算
                    println!("Index file not assigned.");
                } else {
                    mmpbsa_calculation(&trj, &mdp_path, &ndx, &wd, use_dh, use_ts, &settings);
                }
            }
            1 => {
                println!("Input trajectory file path (if in the same directory with tpr, then simply input (e.g.) `?md.xtc`:");
                stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                trj = convert_cur_dir(&trj, &settings);
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr"], &settings);
            }
            2 => {
                println!("Input index file path (if in the same directory with tpr, then simply input (e.g.) `?index.ndx`::");
                stdin().read_line(&mut ndx).expect("Failed while reading index file");
                ndx = convert_cur_dir(&ndx, &settings);
                ndx = confirm_file_validity(&mut ndx, vec!["ndx"], &settings);
            }
            3 => break,
            _ => println!("Error input.")
        };
    }
}

fn welcome() {
    println!("SuperMMPBSA: Supernova's tool of calculating binding free energy using\n\
        molecular mechanics Poisson-Boltzmann surface area (MM-PBSA) method.\n\
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

fn dump_tpr(tpr: &String, dump_to: &String, gmx: &str) {
    let tpr_dump = Command::new(gmx).arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(dump_to).unwrap();
    outfile.write(tpr_dump.as_bytes()).unwrap();
    println!("Finished loading tpr file, md parameters dumped to {}", fs::canonicalize(dump_to).unwrap().display());
}

fn mmpbsa_calculation(trj: &String, mdp: &String, ndx: &String,
                      wd: &Path, use_dh: bool, use_ts: bool, settings: &Parameters) {
    let mut complex_grp: i32 = -1;
    let mut receptor_grp: i32 = -1;
    let mut ligand_grp: i32 = -1;
    let xtc = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let frames: Vec<Rc<Frame>> = xtc.into_iter().map(|p| p.unwrap()).collect();
    let mut bt: f64 = frames[0].time as f64;
    let mut et: f64 = frames[frames.len() - 1].time as f64;
    let mut dt: f64 = (frames[1].time - frames[0].time) as f64;
    let ndx = index_parser::Index::new(ndx);
    loop {
        println!("\n                 ************ MM-PBSA calculation ************");
        println!("-10 Return");
        println!("  0 Do MM-PBSA calculations now!");
        println!("  1 Select complex group, current:            {}", match complex_grp {
            -1 => String::from("undefined"),
            _ => format!("{}): {}, {} atoms",
                         complex_grp,
                         ndx.groups[complex_grp as usize].name,
                         ndx.groups[complex_grp as usize].indexes.len())
        });
        println!("  2 Select receptor groups, current:          {}", match receptor_grp {
            -1 => String::from("undefined"),
            _ => format!("{}): {}, {} atoms",
                         receptor_grp,
                         ndx.groups[receptor_grp as usize].name,
                         ndx.groups[receptor_grp as usize].indexes.len())
        });
        println!("  3 Select ligand groups, current:            {}", match ligand_grp {
            -1 => String::from("undefined"),
            _ => format!("{}): {}, {} atoms",
                         ligand_grp,
                         ndx.groups[ligand_grp as usize].name,
                         ndx.groups[ligand_grp as usize].indexes.len())
        });
        println!("  4 Set start time of analysis, current:      {} ns", bt / 1000.0);
        println!("  5 Set end time of analysis, current:        {} ns", et / 1000.0);
        println!("  6 Set time interval of analysis, current:   {} ps", dt);
        println!("  7 Prepare PB parameters: polar");
        println!("  8 Prepare SA parameters: non-polar");
        let i = get_input_value();
        match i {
            -10 => return,
            0 => {
                let mut sys_name = String::from("_system");
                println!("Input system name (default: {}):", sys_name);
                let mut input = String::new();
                stdin().read_line(&mut input).expect("Error input");
                if input.trim().len() != 0 {
                    sys_name = input.trim().to_string();
                }
                // 定义results形式, 其中应包含所需的全部数据
                let results = mmpbsa::do_mmpbsa_calculations(&trj, mdp, &ndx, wd, &sys_name,
                                                             use_dh, use_ts,
                                                             complex_grp as usize,
                                                             receptor_grp as usize,
                                                             ligand_grp as usize,
                                                             bt, et, dt,
                                                             &settings);
                analyzation::analyze_controller(&sys_name, results);
            }
            1 => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input complex group num:");
                complex_grp = get_input_value();
            }
            2 => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input receptor group num:");
                receptor_grp = get_input_value();
            }
            3 => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input ligand group num:");
                ligand_grp = get_input_value();
            }
            4 => {
                println!("Input start time (ns), should be divisible of {} ns:", dt / 1000.0);
                let mut new_bt = get_input_value::<f64>() * 1000.0;
                while (new_bt - bt) % dt != 0.0 || new_bt > frames[frames.len() - 1].time as f64 || new_bt < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_bt / 1000.0);
                    println!("Input start time (ns) again, should be divisible of {} ns:", dt / 1000.0);
                    new_bt = get_input_value::<f64>() * 1000.0;
                }
                bt = new_bt;
            }
            5 => {
                println!("Input end time (ns), should be divisible of {} ns:", dt / 1000.0);
                let mut new_et = get_input_value::<f64>() * 1000.0;
                while (new_et - et) % dt != 0.0 || new_et > frames[frames.len() - 1].time as f64 || new_et < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_et / 1000.0);
                    println!("Input end time (ns) again, should be divisible of {} ns:", dt / 1000.0);
                    new_et = get_input_value::<f64>() * 1000.0;
                }
                et = new_et;
            }
            6 => {
                println!("Input interval time (ns), should be divisible of {} ns:", dt / 1000.0);
                let mut new_dt = get_input_value::<f64>() * 1000.0;
                while new_dt % dt != 0.0 {
                    println!("The input {} ns is not a valid time step.", new_dt / 1000.0);
                    println!("Input interval time (ns) again, should be divisible of {} ns:", dt / 1000.0);
                    new_dt = get_input_value::<f64>() * 1000.0;
                }
                dt = new_dt;
            }
            _ => println!("Invalid input")
        }
    }
}

fn get_input_value<T: FromStr>() -> T {
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
fn confirm_file_validity(file_name: &String, ext_list: Vec<&str>, settings: &Parameters) -> String {
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
            .join("gmx.exe").to_str().unwrap().to_string()
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
        // gmx = get_built_in_gmx(&env::current_exe().unwrap());
        gmx = get_built_in_gmx();
    }
    match check_program_validity(gmx.as_str(), "GROMACS version:") {
        Ok(p) => {
            gmx_path = p;
            println!("Using Gromacs: {}", gmx);
        }
        Err(_) => {
            println!("Warning: Gromacs not configured correctly: {}. Now trying default gmx.", gmx);
            match check_program_validity("gmx", "GROMACS version") {
                Ok(p) => {
                    gmx_path = p;
                    println!("Using Gromacs: gmx");
                }
                Err(_) => {
                    println!("Warning: default gmx invalid. Now trying built-in gmx of super_mmpbsa.");
                    match check_program_validity(
                        get_built_in_gmx().as_str(),
                        "GROMACS version",
                    ) {
                        Ok(p) => {
                            gmx_path = p;
                            println!("Using Gromacs: {}", get_built_in_gmx().as_str());
                        }
                        Err(_) => {
                            println!("Error: no valid Gromacs program in use.");
                        }
                    }
                }
            }
        }
    }

    // apbs
    match check_program_validity(apbs, "Version") {
        Ok(p) => {
            apbs_path = p;
            println!("Using APBS: {}", apbs);
            fs::remove_file(Path::new("io.mc")).unwrap();
        }
        Err(_) => {
            println!("Warning: APBS not configured correctly. Now trying default apbs.");
            match check_program_validity("apbs", "Version") {
                Ok(p) => {
                    apbs_path = p;
                    println!("Using APBS: apbs");
                    fs::remove_file(Path::new("io.mc")).unwrap();
                }
                Err(_) => {
                    println!("Error: no valid APBS program in use.");
                }
            }
        }
    }

    return (gmx_path, apbs_path);
}

fn check_program_validity(program: &str, target_output: &str) -> Result<String, String> {
    let test_program = Command::new(program).arg("--version").output();
    match test_program {
        Ok(test_program) => {
            let test_program = String::from_utf8(test_program.stdout)
                .expect(format!("Getting {} output failed.", program).as_str());
            if test_program.find(target_output) == None {
                Err(program.to_string())
            } else {
                Ok(program.to_string())
            }
        }
        Err(_) => {
            Err(program.to_string())
        }
    }
}

fn init_settings() -> Parameters {
    let params: Parameters;
    if Path::new("settings.ini").is_file() {
        // Find settings locally
        let settings = fs::read_to_string("settings.ini").unwrap();
        let settings = Regex::new(r"\\").unwrap().replace_all(settings.as_str(), "/").to_string();
        let settings: Value = toml::from_str(settings.as_str()).expect("Error with settings.ini grammar");
        params = read_settings(&settings);
        println!("Note: found settings.ini in the current path. Will use {} kernels.", params.nkernels);
    } else {
        let super_mmpbsa_path = env::current_exe().unwrap()
            .parent().unwrap().join("settings.ini");
        if super_mmpbsa_path.is_file() {
            let settings = fs::read_to_string(super_mmpbsa_path).unwrap();
            let settings = Regex::new(r"\\").unwrap()
                .replace_all(settings.as_str(), "/").to_string();
            let settings: Value = toml::from_str(settings.as_str()).unwrap();
            params = read_settings(&settings);
            println!("Note: found settings.ini in super_mmpbsa directory. Will use {} kernels.", params.nkernels);
        } else {
            params = Parameters {
                rad_type: 1,
                rad_lj0: 1.2,
                mesh_type: 0,
                grid_type: 1,
                cfac: 3.0,
                fadd: 10.0,
                df: 0.5,
                nkernels: 4,
                preserve: true,
                gmx: String::from("gmx"),
                apbs: String::from("apbs"),
                last_opened: String::new(),
            };
            println!("Note: settings.ini not found. Will use 1 kernel.");
        }
    }
    println!("(Currently multi-threading not yet utilized)");
    return params;
}

fn read_settings(settings: &Value) -> Parameters {
    let rad_type = settings.get("radType").unwrap().to_string().parse().unwrap();
    let rad_lj0 = settings.get("radLJ0").unwrap().to_string().parse().unwrap();
    let mesh_type = settings.get("meshType").unwrap().to_string().parse().unwrap();
    let grid_type = settings.get("gridType").unwrap().to_string().parse().unwrap();
    let cfac = settings.get("cfac").unwrap().to_string().parse().unwrap();
    let fadd = settings.get("fadd").unwrap().to_string().parse().unwrap();
    let df = settings.get("df").unwrap().to_string().parse().unwrap();
    let nkernels = settings.get("nkernels").unwrap().to_string().parse().unwrap();
    let preserve = match settings.get("preserve").unwrap().as_str().unwrap() {
        "y" => true,
        "Y" => true,
        _ => false
    };
    let mut gmx = settings.get("gmx").unwrap().to_string();
    gmx = gmx[1..gmx.len() - 1].to_string();
    if gmx.len() == 0 {
        gmx = "gmx".to_string();
    }
    let mut apbs = settings.get("apbs").unwrap().to_string();
    apbs = apbs[1..apbs.len() - 1].to_string();
    if apbs.len() == 0 {
        apbs = "apbs".to_string();
    }
    let mut last_opened = settings.get("last_opened").unwrap().to_string();
    last_opened = last_opened[1..last_opened.len() - 1].to_string();
    if last_opened.len() == 0 {
        last_opened = String::new();
    }
    return Parameters {
        rad_type,
        rad_lj0,
        mesh_type,
        grid_type,
        cfac,
        fadd,
        df,
        nkernels,
        preserve,
        gmx,
        apbs,
        last_opened,
    };
}

fn convert_cur_dir(p: &String, settings: &Parameters) -> String {
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
    let settings = env::current_exe().unwrap().parent().unwrap().join("settings.ini");
    let settings = fs::read_to_string(&settings).unwrap();
    let re = Regex::new("last_opened.*\".*\"").unwrap();
    let last_opened = fs::canonicalize(Path::new(&tpr_mdp)).unwrap().display().to_string();
    let settings = re.replace(settings.as_str(), format!("last_opened = \"{}\"",
                                                         &last_opened));
    let mut settings_file = File::create(env::current_exe().unwrap().parent().unwrap()
        .join("settings.ini")).unwrap();
    settings_file.write_all(settings.as_bytes()).unwrap();
}