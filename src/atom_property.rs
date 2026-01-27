use std::collections::{BTreeSet, HashMap};
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
    pub fn from_tpr(tpr: &TPR, ndx_com: &BTreeSet<usize>) -> AtomProperties {
        // Initialize c6 and c12 arrays
        let atom_types_num = tpr.atom_types_num;
        let mut c6 = Array2::zeros((atom_types_num, atom_types_num));
        let mut c12 = Array2::zeros((atom_types_num, atom_types_num));
        
        // Use zip and iterators to fill arrays
        for (idx, param) in tpr.lj_sr_params.iter().enumerate() {
            let i = idx / atom_types_num;
            let j = idx % atom_types_num;
            c6[[i, j]] = param.c6;
            c12[[i, j]] = param.c12;
        }

        // Pre-allocate atom_props with estimated capacity
        let mut atom_props = Vec::with_capacity(ndx_com.len());
        let mut at_map = HashMap::with_capacity(atom_types_num);
        let mut index = 0;
        let mut cur_atom_id = 0;
        let mut resid_offset = 0;

        for mol in &tpr.molecule_blocks {
            let mol_type = &tpr.molecule_types.iter().find_map(|mt| if mt.molecule_name.eq(&mol.name) {
                Some(mt)
            } else {
                None
            }).unwrap();
            for _ in 0..mol.molecules_num {
                for atom in &mol_type.atoms {
                    if ndx_com.contains(&cur_atom_id) {
                        // Get residue once to avoid multiple lookups
                        let residue = &mol_type.residues[atom.resind];
                        
                        atom_props.push(AtomProperty {
                            charge: atom.charge,
                            radius: atom.radius,
                            type_id: atom.type_id,
                            id: atom.id,
                            name: atom.name.to_string(),
                            resname: residue.name.to_string(),
                            resid: atom.resind + resid_offset,
                        });

                        // Use entry API for more efficient HashMap insertion
                        at_map.entry(atom.at_type.clone())
                            .or_insert_with(|| {
                                let old_index = index;
                                index += 1;
                                old_index
                            });
                    }
                    cur_atom_id += 1;
                }
                resid_offset += mol_type.residues.len();
            }
        }

        // Normalize resid more efficiently
        if let Some(first_resid) = atom_props.first().map(|ap| ap.resid) {
            for ap in &mut atom_props {
                ap.resid -= first_resid;
            }
        }

        AtomProperties {
            c6,
            c12,
            at_map,
            radius_type: "ff".to_string(),  // TODO: Can this be a &'static str?
            atom_props,
        }
    }
}
