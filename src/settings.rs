use std::{env, fs};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use regex::Regex;
use toml::Value;

pub struct Settings {
    pub rad_type: usize,
    pub rad_ff_default: f64,
    pub use_dh: bool,
    pub use_ts: bool,
    pub cfac: f64,
    pub fadd: f64,
    pub r_cutoff: f64,
    pub df: f64,
    pub nkernels: i32,
    pub preserve: bool,
    pub gmx: Option<String>,
    pub apbs: Option<String>,
    pub last_opened: String,
}

pub fn init_settings() -> Settings {
    let mut params = Settings {
        rad_type: 3,
        rad_ff_default: 1.5,
        use_dh: true,
        use_ts: true,
        cfac: 3.0,
        fadd: 10.0,
        r_cutoff: 0.0,
        df: 0.5,
        nkernels: 1,
        preserve: true,
        gmx: None,
        apbs: None,
        last_opened: String::new(),
    };
    
    match find_settings_in_use() {
        Some(settings) => {
            let settings = fs::read_to_string(settings).unwrap();
            let settings = Regex::new(r"\\").unwrap().replace_all(settings.as_str(), "/").to_string();
            let settings: Value = toml::from_str(settings.as_str()).expect("Error with settings.ini's grammar.");
            read_user_settings(&mut params, &settings);
            println!("Note: found settings.ini, will use {} kernels.", params.nkernels);
        }
        None => {
            println!("Note: settings.ini not found, will use 1 kernel.");
        }
    }

    return params;
}

pub fn find_settings_in_use() -> Option<PathBuf> {
    if Path::new("settings.ini").is_file() {
        Some(Path::new("settings.ini").to_path_buf())
    } else {
        let settings_file = env::current_exe()
            .expect("Cannot get current super_mmpbsa program path.")
            .parent().expect("Cannot get current super_mmpbsa program directory.")
            .join("settings.ini");
        if settings_file.is_file() {
            Some(settings_file)
        } else {
            None
        }
    }
}

fn read_user_settings(params: &mut Settings, settings: &Value) {
    params.rad_type = parse_param(settings, "radType", params.rad_type);
    params.rad_ff_default = parse_param(settings, "radDef", params.rad_ff_default);
    params.cfac = parse_param(settings, "cfac", params.cfac);
    params.fadd = parse_param(settings, "fadd", params.fadd);
    params.r_cutoff = parse_param(settings, "r_cutoff", params.r_cutoff);
    if params.r_cutoff == 0.0 {
        params.r_cutoff = f64::INFINITY;
    }
    params.df = parse_param(settings, "df", params.df);
    params.nkernels = parse_param(settings, "nkernels", params.nkernels);
    // String type cannot move
    params.preserve = match settings.get("preserve").unwrap().to_string().as_str() {
        "y" => true,
        "Y" => true,
        _ => false
    };
    let gmx = settings.get("gmx").unwrap().to_string();
    params.gmx = Some(gmx[1..gmx.len() - 1].to_string());
    let apbs = settings.get("apbs").unwrap().to_string();
    params.apbs = Some(apbs[1..apbs.len() - 1].to_string());
    let last_opened = settings.get("last_opened").unwrap().to_string();
    params.last_opened = last_opened[1..last_opened.len() - 1].to_string();
}

fn parse_param<T: FromStr>(settings: &Value, key: &str, default: T) -> T {
    match settings.get(key) {
        Some(v) => match v.to_string().parse::<T>() {
            Ok(v) => v,
            Err(_) => default
        }
        None => default
    }
}