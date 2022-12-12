use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::parse_tpr::{Atom, TPR};

lazy_static! {
    pub static ref RADIUS_TABLE: HashMap<usize, String> = HashMap::from ([
        (0, "ff".to_string()),
        (1, "mBondi".to_string()),
    ]);
}

impl Atom {
    pub fn apply_radius(&mut self, atom_radius_type: usize, ff: f64) {
        match atom_radius_type {
            0 => {
                self.radius = ff;
            }
            1 => {
                self.radius = get_mbondi(self.name.as_str());
            }
            _ => {}
        }
    }
}

// store all atom radius data
pub struct Radius {
    pub radius_type: usize,
    pub radius_name: String,
    pub radii: Vec<f64>,
}

impl Radius {
    pub fn new(radius_type: usize, radii: Vec<f64>) -> Radius {
        Radius { radius_type, radius_name: RADIUS_TABLE[&radius_type].to_string(), radii }
    }
}

impl TPR {
    // ff_radius would not be used if atom_radius_type not 0
    pub fn apply_radius(&mut self, atom_radius_type: usize, ff_radius: &Vec<f64>) {
        let mut idx = 0;
        for mol in &mut self.molecules {
            for _ in 0..self.molecule_types[mol.molecule_type_id].molecules_num {
                for mut atom in &mut mol.atoms {
                    atom.apply_radius(atom_radius_type, ff_radius[idx]);
                    idx += 1;
                }
            }
        }
    }
}

// get atom radius, returns 1.5 if not specified
pub fn get_mbondi(at_type: &str) -> f64 { // mBondi from AMBER20/parmed/tools/changeradii.py
    let rad_bondi: HashMap<&str, f64> = HashMap::from([
        ("C", 1.7), ("H", 1.2), ("N", 1.55), ("HC", 1.3),
        ("O", 1.5), ("HN", 1.3), ("F", 1.5), ("HP", 1.3),
        ("SI", 2.1), ("HO", 0.8), ("P", 1.85), ("HS", 0.8),
        ("S", 1.8), ("CL", 1.7), ("BR", 1.85), ("I", 1.98),
    ]);
    let at_type = at_type.to_uppercase();
    let mut radius = 1.5;
    if at_type.len() >= 2 {
        let r = rad_bondi.get(&at_type[0..2]);
        if let Some(&m) = r {
            radius = m;
        } else {
            let r = rad_bondi.get(&at_type[0..1]);
            if let Some(&m) = r {
                radius = m;
            }
        }
    } else {
        let r = rad_bondi.get(&at_type[0..1]);
        if let Some(&m) = r {
            radius = m;
        }
    }
    return radius;
}