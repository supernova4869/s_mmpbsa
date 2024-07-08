use std::collections::HashMap;
use std::env::current_exe;
use std::fs;
use crate::atom_property::AtomProperty;
use crate::parse_tpr::TPR;
use indicatif::{ProgressBar, ProgressStyle};

impl AtomProperty {
    // ff_radius would not be used
    pub fn apply_radius(&mut self, radius_type: usize, at_list: &Vec<String>, total_at_num: usize, radius_types: &Vec<&str>) {
        let pb = ProgressBar::new(total_at_num as u64);
        pb.set_style(ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:50.cyan/cyan} {percent}% {msg}").unwrap()
            .progress_chars("=>-"));
        let rad_type = radius_types[radius_type];
        let radii_table = get_radii_map(rad_type);
        for (i, r) in &mut self.atom_radius.iter_mut().enumerate() {
            *r = get_radii(&radii_table, &at_list[i]);
            pb.inc(1);
        }
        pb.finish();
    }
}

// get atom radius from dat
pub fn get_radii(radii_table: &HashMap<String, f64>, at_type: &str) -> f64 {
    if at_type.len() >= 2 {
        match radii_table.get(&at_type[0..2]) {
            Some(&m) => m,
            _ => {
                match radii_table.get(&at_type[0..1]) {
                    Some(&m) => m,
                    _ => radii_table["*"]
                }
            }
        }
    } else {
        match radii_table.get(at_type) {
            Some(&m) => m,
            _ => radii_table["*"]
        }
    }
}

impl TPR {
    pub fn get_at_list(&self) -> Vec<String> {
        let mut atom_radius: Vec<String> = vec![];
        for mol in &self.molecules {
            for _ in 0..self.molecule_types[mol.molecule_type_id].molecules_num {
                for atom in &mol.atoms {
                    let at = atom.name.to_uppercase();
                    atom_radius.push(at);
                }
            }
        }
        atom_radius
    }
}

pub fn get_radii_map(rad_type: &str) -> HashMap<String, f64> {
    let mut radii_table: HashMap<String, f64> = HashMap::new();
    let radii_file = current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent().expect("Cannot get current s_mmpbsa program directory.")
        .join("dat").join(format!("{}.dat", &rad_type))
        .to_str().expect("The atom radius data files (dat/) not found.").to_string();
    let radii_file_content = fs::read_to_string(radii_file)
        .expect(format!("Error reading atom radius data file: {}.dat", rad_type).as_str());
    for l in radii_file_content.split("\n").filter(|p| !p.trim().starts_with("//") && !p.trim().is_empty()) {
        let k_v: Vec<&str> = l.split(":").collect();
        radii_table.insert(k_v[0].to_string(), k_v[1].trim().parse::<f64>().unwrap());
    }
    radii_table
}

#[derive(Clone)]
pub struct AD4param {
    pub sigma: f64,
    pub eps: f64,
    pub sigma_hb: f64,
    pub eps_hb: f64
}

pub fn get_ad4_map() -> HashMap<String, AD4param> {
    let mut ad4_map: HashMap<String, AD4param> = HashMap::new();
    let radii_file = current_exe().expect("Cannot get current s_mmpbsa program path.")
        .parent().expect("Cannot get current s_mmpbsa program directory.")
        .join("dat").join("AD4_parameters.dat")
        .to_str().expect("The atom radius data files (dat/) not found.").to_string();
    let radii_file_content = fs::read_to_string(radii_file)
        .expect("Error reading atom radius data file: AD4_parameters.dat");
    for l in radii_file_content.split("\n").filter(|p| !p.trim().starts_with("//") && !p.trim().is_empty()) {
        let k_v: Vec<&str> = l.trim().split(" ").filter(|s| !s.is_empty()).collect();
        let ad4 = AD4param {
            sigma: k_v[1].parse::<f64>().unwrap() / 10.0, // in nm
            eps: k_v[2].parse::<f64>().unwrap() * 4.18,    // in kJ/mol
            sigma_hb: k_v[3].parse::<f64>().unwrap() / 10.0, // in nm
            eps_hb: k_v[4].parse::<f64>().unwrap() * 4.18,    // in kJ/mol
        };
        ad4_map.insert(k_v[0].to_string(), ad4);
    }
    ad4_map
}

// get param from AD4
pub fn get_ad4_param(ad4_map: &HashMap<String, AD4param>, at_type: &str) -> AD4param {
    ad4_map.get(at_type).unwrap().clone()
}