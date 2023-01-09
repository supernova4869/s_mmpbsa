use std::{env, fs};
use std::path::Path;
use std::str::FromStr;
use regex::Regex;
use toml::Value;

pub struct Parameters {
    pub rad_type: usize,
    pub rad_default: f64,
    // pub mesh_type: i32,
    // pub grid_type: i32,
    pub use_dh: bool,
    pub use_ts: bool,
    pub cfac: f64,
    pub fadd: f64,
    pub r_cutoff: f64,
    pub df: f64,
    pub nkernels: i32,
    pub preserve: bool,
    pub gmx: String,
    pub apbs: String,
    pub last_opened: String,
}

pub fn init_settings() -> Parameters {
    let mut params = Parameters {
        rad_type: 1,
        rad_default: 1.2,
        // mesh_type: 0,
        // grid_type: 1,
        use_dh: true,
        use_ts: true,
        cfac: 3.0,
        fadd: 10.0,
        r_cutoff: 0.0,
        df: 0.5,
        nkernels: 4,
        preserve: true,
        gmx: String::new(),
        apbs: String::new(),
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
    params.rad_default = parse_param(settings, "radDef", params.rad_default);
    // params.mesh_type = parse_param(settings, "meshType", params.mesh_type);
    // params.grid_type = parse_param(settings, "gridType", params.grid_type);
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
    match settings.get(key) {
        Some(v) => match v.to_string().parse::<T>() {
            Ok(v) => v,
            Err(_) => default
        }
        None => default
    }
}