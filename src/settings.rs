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
    pub if_alanine_scanning: bool,
}

pub fn init_settings() -> Settings {
    let mut settings = Settings {
        rad_type: 3,
        rad_ff_default: 1.5,
        use_dh: true,
        use_ts: true,
        cfac: 3.0,
        fadd: 10.0,
        r_cutoff: 0.0,
        df: 0.5,
        nkernels: 1,
        preserve: false,
        gmx: None,
        apbs: None,
        last_opened: String::new(),
        if_alanine_scanning: false
    };
    
    match find_settings_in_use() {
        Some(settings_file) => {
            let settings_file = fs::read_to_string(settings_file).unwrap();
            let settings_file = Regex::new(r"\\").unwrap().replace_all(settings_file.as_str(), "/").to_string();
            let setting_values: Value = toml::from_str(settings_file.as_str()).expect("Error with settings.ini's grammar.");
            read_user_settings(&mut settings, &setting_values);
            println!("Note: found settings.ini, will use {} kernels.", settings.nkernels);
        }
        None => {
            println!("Note: settings.ini not found, will use 1 kernel.");
        }
    }

    return settings;
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

fn read_user_settings(settings: &mut Settings, setting_values: &Value) {
    settings.rad_type = parse_param(setting_values, "radType", settings.rad_type);
    settings.rad_ff_default = parse_param(setting_values, "radDef", settings.rad_ff_default);
    settings.cfac = parse_param(setting_values, "cfac", settings.cfac);
    settings.fadd = parse_param(setting_values, "fadd", settings.fadd);
    settings.r_cutoff = parse_param(setting_values, "r_cutoff", settings.r_cutoff);
    if settings.r_cutoff == 0.0 {
        settings.r_cutoff = f64::INFINITY;
    }
    settings.df = parse_param(setting_values, "df", settings.df);
    settings.nkernels = parse_param(setting_values, "nkernels", settings.nkernels);
    settings.preserve = match setting_values.get("preserve").unwrap().to_string()[1..2].to_string().as_str() {
        "y" => true,
        "Y" => true,
        _ => false
    };
    let gmx = setting_values.get("gmx").unwrap().to_string();
    settings.gmx = Some(gmx[1..gmx.len() - 1].to_string());
    let apbs = setting_values.get("apbs").unwrap().to_string();
    settings.apbs = Some(apbs[1..apbs.len() - 1].to_string());
    let last_opened = setting_values.get("last_opened").unwrap().to_string();
    settings.last_opened = last_opened[1..last_opened.len() - 1].to_string();
    settings.if_alanine_scanning = match setting_values.get("alanine_scanning").unwrap().to_string()[1..2].to_string().as_str() {
        "y" => true,
        "Y" => true,
        _ => false
    };
}

fn parse_param<T: FromStr>(setting_values: &Value, key: &str, default: T) -> T {
    match setting_values.get(key) {
        Some(v) => match v.to_string().parse::<T>() {
            Ok(v) => v,
            Err(_) => default
        }
        None => default
    }
}

pub fn get_base_settings() -> PathBuf {
    return env::current_exe().unwrap().parent().unwrap().join("settings.ini");
}
