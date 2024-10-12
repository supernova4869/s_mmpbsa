use std::{env, fs};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use regex::Regex;
use toml::Value;

pub struct Settings {
    pub radius_type: usize,
    pub radius_ff_default: f64,
    pub r_cutoff: f64,
    pub fix_pbc: bool,
    pub elec_screen: usize,
    pub gmx_path: Option<String>,
    pub cfac: f64,
    pub fadd: f64,
    pub df: f64,
    pub pbsa_kernel: Option<String>,
    pub apbs_path: Option<String>,
    pub delphi_path: Option<String>,
    pub pymol_path: Option<String>,
    pub gaussian_dir: Option<String>,
    pub sobtop_dir: Option<String>,
    pub multiwfn_dir: Option<String>,
    pub nkernels: i32,
    pub debug_mode: bool,
    pub last_opened: String,
}

impl Settings {
    pub fn new() -> Settings {
        Settings {
            elec_screen: 1,
            radius_type: 3,
            radius_ff_default: 1.5,
            r_cutoff: 0.0,
            fix_pbc: true,
            gmx_path: Some("gmx".to_string()),
            cfac: 3.0,
            fadd: 10.0,
            df: 0.5,
            pbsa_kernel: None,
            apbs_path: None,
            delphi_path: None,
            pymol_path: None,
            gaussian_dir: None,
            sobtop_dir: None,
            multiwfn_dir: None,
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
        let radius_ff_default = parse_param(&setting_values, "radius_default", default_settings.radius_ff_default);
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
        let gmx_path = parse_param(&setting_values, "gmx_path", "gmx".to_string());
        let gmx_path = Some(gmx_path[1..gmx_path.len() - 1].to_string());
        let cfac = parse_param(&setting_values, "cfac", default_settings.cfac);
        let fadd = parse_param(&setting_values, "fadd", default_settings.fadd);
        let df = parse_param(&setting_values, "df", default_settings.df);
        let pbsa_kernel = parse_param(&setting_values, "pbsa_kernel", "".to_string());
        let pbsa_kernel = Some(pbsa_kernel[1..pbsa_kernel.len() - 1].to_string());
        let apbs_path = parse_param(&setting_values, "apbs_path", "".to_string());
        let apbs_path = Some(apbs_path[1..apbs_path.len() - 1].to_string());
        let delphi_path = parse_param(&setting_values, "delphi_path", "".to_string());
        let delphi_path = Some(delphi_path[1..delphi_path.len() - 1].to_string());
        let pymol_path = parse_param(&setting_values, "pymol_path", "".to_string());
        let pymol_path = Some(pymol_path[1..pymol_path.len() - 1].to_string());
        let gaussian_dir = parse_param(&setting_values, "gaussian_dir", "".to_string());
        let gaussian_dir = Some(gaussian_dir[1..gaussian_dir.len() - 1].to_string());
        let sobtop_dir = parse_param(&setting_values, "sobtop_dir", "".to_string());
        let sobtop_dir = Some(sobtop_dir[1..sobtop_dir.len() - 1].to_string());
        let multiwfn_dir = parse_param(&setting_values, "multiwfn_dir", "".to_string());
        let multiwfn_dir = Some(multiwfn_dir[1..multiwfn_dir.len() - 1].to_string());
        let nkernels = parse_param(&setting_values, "n_kernels", default_settings.nkernels);
        let debug_mode = parse_param(&setting_values, "debug_mode", "\"y\"".to_string());
        let debug_mode = match debug_mode[1..2].to_string().as_str() {
            "y" => true,
            "Y" => true,
            _ => false
        };
        let last_opened = parse_param(&setting_values, "last_opened", "\"\"".to_string());
        let last_opened = last_opened[1..last_opened.len() - 1].to_string();

        println!("Note: found settings.ini, will use {} kernels.", nkernels);

        Settings {
            elec_screen,
            radius_type,
            radius_ff_default,
            r_cutoff,
            fix_pbc,
            gmx_path,
            cfac,
            fadd,
            df,
            pbsa_kernel,
            apbs_path,
            delphi_path,
            pymol_path,
            gaussian_dir,
            sobtop_dir,
            multiwfn_dir,
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
