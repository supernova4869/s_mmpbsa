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
    pub fix_pbc: bool,
    pub preserve: bool,
    pub gmx: Option<String>,
    pub apbs: Option<String>,
    pub last_opened: String,
    pub if_alanine_scanning: bool,
    pub debug_mode: bool,
}

impl Settings {
    pub fn new() -> Settings {
        Settings {
            rad_type: 3,
            rad_ff_default: 1.5,
            use_dh: true,
            use_ts: true,
            cfac: 3.0,
            fadd: 10.0,
            r_cutoff: 0.0,
            df: 0.5,
            nkernels: 1,
            fix_pbc: true,
            preserve: false,
            gmx: Some("gmx".to_string()),
            apbs: None,
            last_opened: String::new(),
            if_alanine_scanning: false,
            debug_mode: false,
        }
    }

    pub fn from(settings_file: &PathBuf) -> Settings {
        let settings_file = fs::read_to_string(settings_file).unwrap();
        let settings_file = Regex::new(r"\\").unwrap().replace_all(settings_file.as_str(), "/").to_string();
        let setting_values: Value = toml::from_str(settings_file.as_str()).expect("Error with settings.ini's grammar.");
        let default_settings = Settings::new();
        
        // Read settings
        let rad_type = parse_param(&setting_values, "radType", default_settings.rad_type);
        let rad_ff_default = parse_param(&setting_values, "radDef", default_settings.rad_ff_default);
        let cfac = parse_param(&setting_values, "cfac", default_settings.cfac);
        let fadd = parse_param(&setting_values, "fadd", default_settings.fadd);
        let r_cutoff = parse_param(&setting_values, "r_cutoff", default_settings.r_cutoff);
        let r_cutoff = if r_cutoff == 0.0 {
            f64::INFINITY
        } else {
            r_cutoff
        };
        let df = parse_param(&setting_values, "df", default_settings.df);
        let nkernels = parse_param(&setting_values, "nkernels", default_settings.nkernels);
        let fix_pbc = parse_param(&setting_values, "fix_pbc", "\"y\"".to_string());
        let fix_pbc = match fix_pbc[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let preserve = parse_param(&setting_values, "preserve", "\"y\"".to_string());
        let preserve = match preserve[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let gmx = parse_param(&setting_values, "gmx", "gmx".to_string());
        let gmx = Some(gmx[1..gmx.len() - 1].to_string());
        let apbs = parse_param(&setting_values, "apbs", "".to_string());
        let apbs = Some(apbs[1..apbs.len() - 1].to_string());
        let last_opened = parse_param(&setting_values, "last_opened", "\"\"".to_string());
        let last_opened = last_opened[1..last_opened.len() - 1].to_string();
        let if_alanine_scanning = parse_param(&setting_values, "alanine_scanning", "\"y\"".to_string());
        let if_alanine_scanning = match if_alanine_scanning[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let debug_mode = parse_param(&setting_values, "debug_mode", "\"y\"".to_string());
        let debug_mode = match debug_mode[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };

        println!("Note: found settings.ini, will use {} kernels.", nkernels);

        Settings {
            rad_type,
            rad_ff_default,
            use_dh: true,
            use_ts: true,
            cfac,
            fadd,
            r_cutoff,
            df,
            nkernels,
            fix_pbc,
            preserve,
            gmx,
            apbs,
            last_opened,
            if_alanine_scanning,
            debug_mode,
        }
    }
}

pub fn get_settings_in_use() -> Option<PathBuf> {
    if Path::new("settings.ini").is_file() {
        Some(Path::new("settings.ini").to_path_buf())
    } else {
        let settings_file = get_base_settings();
        if settings_file.is_file() {
            Some(settings_file)
        } else {
            None
        }
    }
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
    env::current_exe()
            .expect("Cannot get current s_mmpbsa program path.")
            .parent().expect("Cannot get current s_mmpbsa program directory.")
            .join("settings.ini")
}
