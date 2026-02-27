use std::{env, fs};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use regex::Regex;
use toml::Value;

#[allow(dead_code)]
pub struct Settings {
    pub radius_type: usize,
    pub r_cutoff: f64,
    pub fix_pbc: bool,
    pub elec_screen: usize,
    pub calc_mm: bool,
    pub calc_pbsa: bool,
    pub gmx_path: Option<String>,
    pub apbs_path: Option<String>,
    pub cfac: i32,
    pub fadd: f64,
    pub df: f64,
    pub chg_m: usize,
    pub pymol_path: Option<String>,
    pub nkernels: usize,
    pub debug_mode: bool,
    pub last_opened: String,
}

impl Settings {
    pub fn new() -> Settings {
        Settings {
            elec_screen: 1,
            radius_type: 3,
            r_cutoff: 0.0,
            fix_pbc: true,
            gmx_path: Some("gmx".to_string()),
            apbs_path: None,
            calc_mm: true,
            calc_pbsa: true,
            cfac: 3,
            fadd: 10.0,
            df: 0.5,
            chg_m: 0,
            pymol_path: None,
            nkernels: 1,
            debug_mode: false,
            last_opened: String::new(),
        }
    }

    pub fn from(settings_file: &PathBuf) -> Settings {
        let settings_file = fs::read_to_string(settings_file).unwrap();
        let settings_file = Regex::new(r"\\").unwrap().replace_all(settings_file.as_str(), "/").to_string();
        let setting_values: Value = toml::from_str(settings_file.as_str()).expect("Error with settings.ini's grammar.");
        let default_settings = Settings::new();
        
        // Read settings
        let elec_screen = parse_param(&setting_values, "screen_method", default_settings.elec_screen);
        let radius_type = parse_param(&setting_values, "radius_type", default_settings.radius_type);
        let r_cutoff = parse_param(&setting_values, "r_cutoff", default_settings.r_cutoff);
        let r_cutoff = if r_cutoff == 0.0 {
            f64::INFINITY
        } else {
            r_cutoff
        };
        let fix_pbc = parse_param(&setting_values, "fix_pbc", "\"y\"".to_string());
        let fix_pbc = match fix_pbc[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let gmx_path = parse_param(&setting_values, "gmx_path", "\"built-in\"".to_string());
        let gmx_path = Some(gmx_path[1..gmx_path.len() - 1].to_string());
        let apbs_path = parse_param(&setting_values, "apbs_path", "".to_string());
        let apbs_path = Some(apbs_path.trim_start_matches('\"').trim_end_matches('\"').to_string());
        let calc_mm = parse_param(&setting_values, "calc_mm", "\"y\"".to_string());
        let calc_mm = match calc_mm[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let calc_pbsa = parse_param(&setting_values, "calc_pbsa", "\"y\"".to_string());
        let calc_pbsa = match calc_pbsa[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let cfac = parse_param(&setting_values, "cfac", default_settings.cfac);
        let fadd = parse_param(&setting_values, "fadd", default_settings.fadd);
        let df = parse_param(&setting_values, "df", default_settings.df);
        let chg_m = parse_param(&setting_values, "chg_m", 0);
        let pymol_path = parse_param(&setting_values, "pymol_path", "".to_string());
        let pymol_path = Some(pymol_path.trim_start_matches('\"').trim_end_matches('\"').to_string());
        let nkernels = parse_param(&setting_values, "n_kernels", default_settings.nkernels);
        let debug_mode = parse_param(&setting_values, "debug_mode", "\"y\"".to_string());
        let debug_mode = match debug_mode[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let last_opened = parse_param(&setting_values, "last_opened", "\"\"".to_string());
        let last_opened = last_opened[1..last_opened.len() - 1].to_string();

        Settings {
            elec_screen,
            radius_type,
            r_cutoff,
            fix_pbc,
            gmx_path,
            apbs_path,
            calc_mm,
            calc_pbsa,
            cfac,
            fadd,
            df,
            chg_m,
            pymol_path,
            nkernels,
            debug_mode,
            last_opened,
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
        Some(v) => v.to_string().parse::<T>().unwrap_or(default),
        None => default
    }
}

pub fn get_base_settings() -> PathBuf {
    env::current_exe()
            .expect("Cannot get current s_mmpbsa program path.")
            .parent().expect("Cannot get current s_mmpbsa program directory.")
            .join("settings.ini")
}
