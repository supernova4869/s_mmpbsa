use std::collections::HashMap;

use ndarray::Array2;
use crate::{atom_radius::{get_radii, get_radii_map}, parse_tpr::TPR};

#[derive(Clone)]
pub struct AtomProperties {
    pub c6: Array2<f64>,
    pub c12: Array2<f64>,
    pub at_map: HashMap<String, usize>,
    pub radius_type: String,
    pub atom_props: Vec<AtomProperty>
}

#[derive(Clone)]
pub struct AtomProperty {
    pub charge: f64,
    pub radius: f64,
    pub type_id: usize,
    pub id: usize,
    pub name: String,
    pub resname: String,
    pub resid: usize,
}

impl AtomProperty {
    pub fn change_atom(&mut self, new_type_id: Option<&usize>, new_name: &str, radius_type: &str) {
        if let Some(&new_type_id) = new_type_id {
            self.type_id = new_type_id;
            self.name = new_name.to_string();
            let radii_table = get_radii_map(radius_type);
            self.radius = get_radii(&radii_table, new_name);
        }
    }
}

impl AtomProperties {
    pub fn from_tpr(tpr: &TPR, ndx_com: &Vec<usize>) -> AtomProperties {
        // c6 and c12
        let mut c6: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        let mut c12: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        for i in 0..tpr.atom_types_num {
            for j in 0..tpr.atom_types_num {
                c6[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c6;
                c12[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c12;
            }
        }

        let mut atom_props: Vec<AtomProperty> = vec![];

        let mut cur_atom_id = 0;
        let mut resid_offset = 0;      // residues number that has been overpast

        // HashMap to store the first occurrence index of each string
        let mut at_map: HashMap<String, usize> = HashMap::new();
        let mut index = 0;
        for mol in &tpr.molecules {
            for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                for atom in &mol.atoms {
                    if ndx_com.contains(&cur_atom_id) {
                        atom_props.push(AtomProperty {
                            charge: atom.charge,
                            radius: atom.radius,
                            type_id: atom.type_id,
                            id: atom.id,
                            name: atom.name.to_string(),
                            resname: mol.residues[atom.resind].name.to_string(),
                            resid: atom.resind + resid_offset,
                        });
                        if !at_map.contains_key(&atom.at_type) {
                            // If the string is not in the map, insert it with the current index
                            at_map.insert(atom.at_type.to_string(), index);
                            index += 1;
                        }
                    }
                    cur_atom_id += 1;
                }
                resid_offset += mol.residues.len();
            }
        }

        // normalize resid
        let first_resid = atom_props[0].resid;
        for ap in atom_props.iter_mut() {
            ap.resid -= first_resid;
        }

        AtomProperties {
            c6,
            c12,
            at_map,
            radius_type: "ff".to_string(),
            atom_props
        }
    }
}