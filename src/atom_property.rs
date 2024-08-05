use std::collections::HashMap;

use ndarray::{Array1, Array2};
use crate::{atom_radius::{get_ad4_map, get_ad4_param, get_radii, get_radii_map}, parse_pdbqt::PDBQT, parse_tpr::TPR};
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Clone)]
pub struct AtomProperties {
    pub c6: Array2<f64>,
    pub c12: Array2<f64>,
    pub c10: Array2<f64>,
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
        let c10: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        for i in 0..tpr.atom_types_num {
            for j in 0..tpr.atom_types_num {
                c6[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c6;
                c12[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c12;
            }
        }

        let mut atom_props: Vec<AtomProperty> = vec![];

        let mut cur_atom_id = 0;
        let mut resid_offset = 0;      // residues number that has been overpast

        let pb = ProgressBar::new(tpr.atom_types_num as u64);
        pb.set_style(ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:50.cyan/cyan} {percent}% {msg}").unwrap()
            .progress_chars("=>-"));

        let mut at_list = vec![];
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
                        at_list.push(atom.at_type.to_string());
                    }
                    cur_atom_id += 1;
                    pb.inc(1);
                    pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
                }
                resid_offset += mol.residues.len();
            }
        }

        // normalize resid
        let first_resid = atom_props[0].resid;
        for ap in atom_props.iter_mut() {
            ap.resid -= first_resid;
        }

        pb.finish();
        
        // HashMap to store the first occurrence index of each string
        let mut at_map: HashMap<String, usize> = HashMap::new();
        let mut ordered_atom_types = Vec::new();
        
        let mut atom_type_id: Array1<usize> = Array1::zeros(ndx_com.len());
        let mut index = 0;
        
        for (i, s) in at_list.iter().enumerate() {
            if !at_map.contains_key(s) {
                // If the string is not in the map, insert it with the current index
                at_map.insert(s.to_string(), index);
                ordered_atom_types.push(s);
                index += 1;
            }
            atom_type_id[i] = at_map[s];
        }

        AtomProperties {
            c6,
            c12,
            c10,
            at_map,
            radius_type: "ff".to_string(),
            atom_props
        }
    }

    pub fn from_pdbqt(receptor: &PDBQT, ligand: &PDBQT) -> AtomProperties {
        let mut atoms = receptor.models[0].atoms.to_vec();
        atoms.extend(ligand.models[0].atoms.to_vec());
        let atom_num = atoms.len();

        // 将每个原子的原子类型进行简并编号
        let at_list: Vec<String> = atoms.iter().map(|a| a.attype.to_owned()).collect();
        
        // HashMap to store the first occurrence index of each string
        let mut at_map: HashMap<String, usize> = HashMap::new();
        let mut ordered_atom_types = Vec::new();
        
        let mut atom_type_id: Array1<usize> = Array1::zeros(atom_num);
        let mut index = 0;
        for (i, s) in at_list.iter().enumerate() {
            if !at_map.contains_key(s) {
                // If the string is not in the map, insert it with the current index
                at_map.insert(s.to_string(), index);
                ordered_atom_types.push(s);
                index += 1;
            }
            atom_type_id[i] = at_map[s];
        }

        // AD4 parameters for all atoms
        let ad4_map = get_ad4_map();

        // c6 and c12 parameters for all atoms
        let mut c6: Array2<f64> = Array2::zeros((at_map.len(), at_map.len()));
        let mut c12: Array2<f64> = Array2::zeros((at_map.len(), at_map.len()));
        let mut c10: Array2<f64> = Array2::zeros((at_map.len(), at_map.len()));
        for i in 0..at_map.len() {
            for j in 0..at_map.len() {
                let sigmai = get_ad4_param(&ad4_map, ordered_atom_types[i]).sigma;
                let epsi = get_ad4_param(&ad4_map, ordered_atom_types[i]).eps;
                let sigmaj = get_ad4_param(&ad4_map, ordered_atom_types[j]).sigma;
                let epsj = get_ad4_param(&ad4_map, ordered_atom_types[j]).eps;
                let sigmaij = (sigmai + sigmaj) / 2.0;
                let epsij = f64::sqrt(epsi * epsj) * 0.1560;
                c12[[i, j]] = epsij * sigmaij.powi(12);
                c6[[i, j]] = 2.0 * epsij * sigmaij.powi(6);

                // for H-bond
                let hb_h_ids: Vec<usize> = at_map.iter().filter_map(|(k, &v)| {
                    if k.eq("HD") || k.eq("HS") {
                        Some(v)
                    } else {
                        None
                    }
                }).collect();
                let hb_a_ids: Vec<usize> = at_map.iter().filter_map(|(k, &v)| {
                    if k.eq("NA") || k.eq("NS") || k.eq("OA") || k.eq("OS") || k.eq("SA") {
                        Some(v)
                    } else {
                        None
                    }
                }).collect();

                // From autodock:
                // Note that the Rij_hb value is non-zero for heteroatoms only, and zero for H atoms;
                // to obtain the length of an H-bond, look up Rij_hb for the heteroatom only; 
                // this is combined with the Rii value for H in the receptor, in AutoGrid.
                // For example, the Rij_hb for OA-HD H-bonds will be (1.9 + 1.0) Angstrom, 
                // and the weighted epsij_hb will be 5.0 kcal/mol * FE_coeff_hbond.
                if hb_h_ids.contains(&i) && hb_a_ids.contains(&j) {
                    let sigma_heavy_hb = get_ad4_param(&ad4_map, ordered_atom_types[j]).sigma_hb;
                    let eps_heavy_hb = get_ad4_param(&ad4_map, ordered_atom_types[j]).eps_hb;
                    let sigmaij_hb = sigmai / 2.0 + sigma_heavy_hb;
                    let epsij_hb = epsi * f64::sqrt(eps_heavy_hb) * 0.0974;
                    c12[[i, j]] = 5.0 * epsij_hb * sigmaij_hb.powi(12);
                    c10[[i, j]] = 6.0 * epsij_hb * sigmaij_hb.powi(10);
                } else if hb_h_ids.contains(&j) && hb_a_ids.contains(&i) {
                    let sigma_heavy_hb = get_ad4_param(&ad4_map, ordered_atom_types[i]).sigma_hb;
                    let eps_heavy_hb = get_ad4_param(&ad4_map, ordered_atom_types[i]).eps_hb;
                    let sigmaij_hb = sigmaj / 2.0 + sigma_heavy_hb;
                    let epsij_hb = epsj * f64::sqrt(eps_heavy_hb) * 0.0974;
                    c12[[i, j]] = 5.0 * epsij_hb * sigmaij_hb.powi(12);
                    c10[[i, j]] = 6.0 * epsij_hb * sigmaij_hb.powi(10);
                }
            }
        }

        let mut atom_props : Vec<AtomProperty> = vec![];

        // 每个原子对应的残基编号
        let resid_offset_receptor = receptor.models[0].atoms[0].resid;
        let resid_offset_ligand = ligand.models[0].atoms[0].resid;
        let mut atom_resid: Vec<usize> = receptor.models[0].atoms.iter().map(|at| (at.resid - resid_offset_receptor) as usize).collect();
        let resid_receptor_last = atom_resid.last().unwrap().to_owned();
        let resid_ligand = ligand.models[0].atoms.iter()
            .map(|at| (at.resid - resid_offset_ligand) as usize + resid_receptor_last + 1);
        atom_resid.extend(resid_ligand);
        let atom_resid: Array1<usize> = Array1::from_vec(atom_resid);

        for (i, atom) in atoms.iter().enumerate() {
            atom_props.push(AtomProperty {
                charge: atom.charge,
                radius: get_ad4_param(&ad4_map, &atom.attype).sigma * 5.0,
                type_id: atom_type_id[i],
                id: i,
                name: atom.atname.to_string(),
                resname: atom.resname.to_string(),
                resid: atom_resid[i],
            });
        }

        AtomProperties {
            c6,
            c12,
            c10,
            at_map,
            radius_type: "AD4_parameters".to_string(),
            atom_props
        }
    }
}