use std::collections::HashMap;
use std::env::current_exe;
use std::fs;
use crate::atom_property::AtomProperty;
use crate::parse_tpr::TPR;
use indicatif::{ProgressBar, ProgressStyle};

impl AtomProperty {
    // ff_radius would not be used if radius_type not 0
    pub fn apply_radius(&mut self, radius_type: usize, tpr: &TPR, total_at_num: usize, radius_types: &Vec<&str>) {
        let mut idx = 0;
        
        let pb = ProgressBar::new(total_at_num as u64);
        pb.set_style(ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:50.cyan/cyan} {percent}% {msg}").unwrap()
            .progress_chars("=>-"));
        match radius_type {
            0 => {
                for mol in &tpr.molecules {
                    if idx >= total_at_num {
                        break;
                    }
                    for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                        if idx >= total_at_num {
                            break;
                        }
                        else {
                            for atom in &mol.atoms {
                                if idx >= total_at_num {
                                    break;
                                } else {
                                    self.atm_radius[idx] = atom.radius;
                                    idx += 1;
                                    pb.inc(1);
                                    pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
                                }
                            }
                        }
                    }
                };
            }
            _ => {
                let mut radii_table: HashMap<&str, f64> = HashMap::new();
                let rad_type = radius_types[radius_type];
                let radii_file = current_exe().expect("Cannot get current s_mmpbsa program path.")
                    .parent().expect("Cannot get current s_mmpbsa program directory.")
                    .join("dat").join(format!("{}.dat", &rad_type))
                    .to_str().expect("The atom radius data files (dat/) not found.").to_string();
                let radii_file = fs::read_to_string(radii_file)
                    .expect(format!("The atom radius data file not found: {}.dat", rad_type).as_str());
                for l in radii_file.split("\n").filter(|p| !p.trim().starts_with("//")) {
                    let k_v: Vec<&str> = l.split(":").collect();
                    radii_table.insert(k_v[0], k_v[1].trim().parse::<f64>().unwrap());
                }

                for mol in &tpr.molecules {
                    if idx >= total_at_num {
                        break;
                    }
                    for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                        if idx >= total_at_num {
                            break;
                        }
                        for atom in &mol.atoms {
                            if idx >= total_at_num {
                                break;
                            }
                            self.atm_radius[idx] = get_radii(&radii_table, &atom.name.as_str().to_uppercase());
                            idx += 1;
                            pb.inc(1);
                            pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
                        }
                    }
                }
            }
        }
        pb.finish();
    }
}

// get atom radius from dat
pub fn get_radii(radii_table: &HashMap<&str, f64>, at_type: &str) -> f64 {
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
