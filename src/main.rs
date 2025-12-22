mod parse_ndx;
mod mmpbsa;
mod parse_tpr;
mod parse_pdb;
mod parse_gro;
mod read_xtc;
mod analyzation;
mod fun_para_basic;
mod fun_para_system;
mod fun_para_mmpbsa;
mod atom_radius;
mod parameters;
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
use colored::Colorize;
use regex::Regex;
use settings::{Settings, get_base_settings, get_settings_in_use};
use utils::get_input;
use clap::Parser;

use crate::parameters::Config;

#[derive(Parser)]
#[command(name = "s_mmpbsa")]
#[command(disable_version_flag = true)]
struct Cli {
    /// input xtc file path
    #[arg(short = 'f', long, value_name = "md.xtc", default_value = None)]
    xtc: Option<String>,

    /// input tpr file path
    #[arg(short = 's', long, value_name = "md.tpr", default_value = None)]
    tpr: Option<String>,

    /// input ndx file path
    #[arg(short = 'n', long, value_name = "index.ndx", default_value = None)]
    ndx: Option<String>,
    
    /// enter analyzation mode
    #[arg(short, long)]
    analyze: bool,
    
    /// generate template config file
    #[arg(short = 'p', long, value_name = "template")]
    template: bool,
    
    /// assign config file path
    #[arg(short, long, value_name = "config.yaml")]
    config: Option<String>,
    
    /// show version info
    #[arg(short = 'V', long)]
    version: bool,
}

fn main() {
    let cli = Cli::parse();
    let compile_date = "2025-Dec-21";
    welcome(&env!("CARGO_PKG_VERSION"), compile_date);
    
    // Show version info
    if cli.version {
        utils::show_famous_quotes();
        exit(0);
    }

    // Build template
    if cli.template {
        // 在当前路径生成一个config模板
        let config = Config::new();
        config.save(&Path::new(&env::current_dir().unwrap()).join("config.yaml"));
        println!("Template config has been written to {}", "config.yaml".cyan().bold());
        println!("You can edit and use it by `{}`", "s_mmpbsa -c config.yaml".red().bold());
        exit(0);
    }

    let mut settings = env_check();
    match settings.debug_mode {
        true => println!("Debug mode on.\n"),
        false => println!("Debug mode off.\n"),
    }

    // Analyzation mode
    if cli.analyze {
        println!("Input path of working dir with .sm results (default: current dir):");
        let input = get_input("".to_string());
        let wd = if input.is_empty() {
            env::current_dir().unwrap()
        } else {
            Path::new(&input).to_path_buf()
        };
        if !wd.is_dir() {
            println!("Input {} not directory. Please check.", Path::new(&input).to_str().unwrap());
            get_input("".to_string());
            exit(0);
        }
        let sm_list: Vec<String> = fs::read_dir(&wd).unwrap().into_iter().filter_map(|f| {
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
            analyzation::analyze_controller(&result_wt, &result_as, temperature, &sys_name, &wd, &settings);
        } else {
            println!("There is no MM/PB-SA results at {}. Please run MM/PB-SA calculations first.", &input);
        }
    };

    let tpr = if cli.tpr.is_some() {
        cli.tpr.unwrap().to_string()
    } else {
        println!("Input path of tpr file, e.g. D:/md.tpr");
        println!("Hint: input \"o\" to simply load last-opened file.");
        let mut input = String::new();
        stdin().read_line(&mut input).expect("Failed to get input file.");
        if input.trim().eq("o") {
            input = settings.last_opened.to_string();
            if input.is_empty() {
                println!("Last-opened tpr not found.");
            }
        }
        input
    };

    let trj = cli.xtc;
    let ndx = cli.ndx;
    
    // Config file
    let config = if cli.config.is_some() {
        println!("Loading config file: {}", cli.config.as_deref().unwrap().cyan().bold());
        let mut config = Config::load(&cli.config.as_deref().unwrap());
        loop {
            match config {
                Ok(_) => break,
                Err(e) => {
                    println!("Error format with config file, details:\n{}", e.to_string().red().bold());
                    println!("Edit {} again and press ENTER to reload", &cli.config.as_deref().unwrap());
                    get_input("".to_string());
                    config = Config::load(&cli.config.as_deref().unwrap());
                }
            };
        }
        config.ok()
    } else {
        None
    };
    if Path::new(&tpr).is_file() {
        let tpr = confirm_file_validity(&tpr, vec!["tpr"], &tpr);
        change_settings_last_opened(&mut settings, &tpr);
        let tpr = get_dump(&tpr, &settings);
        fun_para_basic::set_para_basic_tpr(&tpr, &trj, &ndx, &config, &Path::new(&tpr).parent().unwrap(), &mut settings);
    } else {
        println!("Input {} not file. Please check.", Path::new(&tpr).to_str().unwrap());
        println!("Press ENTER to exit.");
        get_input("".to_string());
        exit(0);
    }
}

fn welcome(version: &str, today: &str) {
    println!("\
        ========================================================================\n\
        | s_mmpbsa: Supernova's tool of calculating binding free energy using  |\n\
        | molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method |\n\
        ========================================================================\n\
        Website: https://github.com/supernova4869/s_mmpbsa\n\
        Latest documentation: https://s-mmpbsa.readthedocs.io/en/latest/\n\
        Developed by Supernova (zhangjiaxing7137@tju.edu.cn), Tianjin University.\n\
        Version {}, first release: 2022-Oct-17, current release: {}\n", version, today);
    println!("Usage 1: run `s_mmpbsa` and follow the prompts.\n\
        Usage 2: run `s_mmpbsa -s Haibara_Ai.tpr` to load MD tpr file.\n");
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

fn set_program(p: &Option<String>, name: &str, settings: &Settings) -> Option<String> {
    if let Some(p) = p {
        let p = if p.eq("built-in") {
            match name {
                "gromacs" => get_built_in_gmx(),
                "apbs" => get_built_in_apbs(),
                "delphi" => get_built_in_delphi(),
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
    let output = Command::new(program)
        .arg("--version")
        .output()
        .map_err(|_| ())?;
    
    // 定义可接受的退出码
    let valid_codes = [0, 1, 13, 127];
    
    match output.status.code() {
        Some(code) if valid_codes.contains(&code) => Ok(program.to_string()),
        _ => Err(())
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
    let cur_path = env::current_exe().unwrap();
    let cur_path = cur_path.parent().unwrap();
    if !Path::new(cur_path.join("dat/").as_path()).is_dir() {
        println!("Error: the dat/ folder with atom radius not found, please check and retry.");
        io::stdin().read_line(&mut String::new()).unwrap();
        std::process::exit(0);
    }
    settings.gmx_path = set_program(&settings.gmx_path, "gromacs", &settings);
    settings.apbs_path = set_program(&settings.apbs_path, "apbs", &settings);
    settings.delphi_path = set_program(&settings.delphi_path, "delphi", &settings);
    settings
}