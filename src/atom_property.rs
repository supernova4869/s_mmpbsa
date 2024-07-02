use std::{collections::{HashMap, HashSet}, ops::Deref};

use ndarray::{Array1, Array2};
use crate::{atom_radius::{get_ad4_map, get_ad4_param, get_radii, get_radii_map, AD4param}, parse_pdbqt::PDBQT, parse_tpr::TPR};
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Clone)]
pub struct AtomProperty {
    pub c6: Array2<f64>,
    pub c12: Array2<f64>,
    pub c10: Array2<f64>,
    pub atom_charge: Array1::<f64>,
    pub atom_radius: Array1::<f64>,
    pub atom_type_id: Array1<usize>,
    pub atom_id: Array1<usize>,
    pub atom_name: Array1<String>,
    pub atom_resname: Array1<String>,
    pub atom_resid: Array1<usize>,
}

impl AtomProperty {
    pub fn from_tpr(tpr: &TPR, ndx_com: &Vec<usize>) -> AtomProperty {
        // c6 and c12
        let mut c6: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        let mut c12: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        let mut c10: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        for i in 0..tpr.atom_types_num {
            for j in 0..tpr.atom_types_num {
                c6[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c6;
                c12[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c12;
            }
        }

        let mut atom_charge: Array1<f64> = Array1::zeros(ndx_com.len());
        let mut atom_radius: Array1<f64> = Array1::zeros(ndx_com.len());
        let mut atom_type_id: Array1<usize> = Array1::zeros(ndx_com.len());
        let mut atom_id: Array1<usize> = Array1::zeros(ndx_com.len());
        let mut atom_name: Array1<String> = Array1::default(ndx_com.len());
        let mut atom_resname: Array1<String> = Array1::default(ndx_com.len());
        let mut atom_resid: Array1<usize> = Array1::zeros(ndx_com.len());

        let mut index_total = 0;
        let mut index = 0;
        let mut resid_offset = 0;      // residues number that has been overpast

        let mut size: u64 = 0;
        for mol in &tpr.molecules {
            for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                size += mol.atoms.len() as u64;
            }
        }
        let pb = ProgressBar::new(size);
        pb.set_style(ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:50.cyan/cyan} {percent}% {msg}").unwrap()
            .progress_chars("=>-"));

        for mol in &tpr.molecules {
            for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                for atom in &mol.atoms {
                    if ndx_com.contains(&index_total) {
                        atom_charge[index] = atom.charge;
                        atom_radius[index] = atom.radius;
                        atom_type_id[index] = atom.type_id;
                        atom_id[index] = atom.id;
                        atom_name[index] = atom.name.to_string();
                        atom_resname[index] = mol.residues[atom.resind].name.to_string();
                        atom_resid[index] = atom.resind + resid_offset;
                        index += 1;
                    }
                    index_total += 1;
                    pb.inc(1);
                    pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
                }
                resid_offset += mol.residues.len();
            }
        }

        // normalize resid
        let atom_resid = &atom_resid - atom_resid[0];

        pb.finish();

        AtomProperty {
            c6,
            c12,
            c10,
            atom_charge,
            atom_radius,
            atom_type_id,
            atom_id,
            atom_name,
            atom_resname,
            atom_resid
        }
    }

    pub fn from_pdbqt(receptor: &PDBQT, ligand: &PDBQT) -> AtomProperty {
        let mut atoms = receptor.models[0].atoms.to_vec();
        atoms.extend(ligand.models[0].atoms.to_vec());
        let atom_num = atoms.len();

        // 将每个原子的原子类型进行简并编号
        let at_list: Vec<String> = atoms.iter().map(|a| a.attype.to_owned()).collect();
        
        // HashMap to store the first occurrence index of each string
        let mut at_map = HashMap::new();
        let mut ordered_atom_types = Vec::new();
        
        let mut atom_type_id: Array1<usize> = Array1::zeros(atom_num);
        let mut index = 0;
        for (i, s) in at_list.iter().enumerate() {
            if !at_map.contains_key(s) {
                // If the string is not in the map, insert it with the current index
                at_map.insert(s, index);
                ordered_atom_types.push(s);
                index += 1;
            }
            atom_type_id[i] = at_map[s];
        }
        println!("{}", atom_type_id);

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
                let hb_h_ids: Vec<usize> = at_map.iter().filter_map(|(&k, &v)| {
                    if k.eq("HD") || k.eq("HS") {
                        Some(v)
                    } else {
                        None
                    }
                }).collect();
                let hb_a_ids: Vec<usize> = at_map.iter().filter_map(|(&k, &v)| {
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
                    let epsij_hb = f64::sqrt(epsi * eps_heavy_hb) * 0.0974 * 0.1560;
                    c12[[i, j]] = 5.0 * epsij_hb * sigmaij_hb.powi(12);
                    c10[[i, j]] = 6.0 * epsij_hb * sigmaij_hb.powi(10);
                } else if hb_h_ids.contains(&j) && hb_a_ids.contains(&i) {
                    let sigma_heavy_hb = get_ad4_param(&ad4_map, ordered_atom_types[i]).sigma_hb;
                    let eps_heavy_hb = get_ad4_param(&ad4_map, ordered_atom_types[i]).eps_hb;
                    let sigmaij_hb = sigmaj / 2.0 + sigma_heavy_hb;
                    let epsij_hb = f64::sqrt(epsj * eps_heavy_hb) * 0.0974 * 0.1560;
                    c12[[i, j]] = 5.0 * epsij_hb * sigmaij_hb.powi(12);
                    c10[[i, j]] = 6.0 * epsij_hb * sigmaij_hb.powi(10);
                }
            }
        }

        let atom_charge: Array1<f64> = Array1::from_iter(atoms.iter().map(|at| at.charge));
        let atom_radius: Array1<f64> = Array1::from_iter(atoms.iter().map(|at| get_ad4_param(&ad4_map, &at.attype).sigma * 5.0));
        let atom_id: Array1<usize> = Array1::from_iter(0..atoms.len());
        let atom_name: Array1<String> = Array1::from_iter(atoms.iter().map(|at| at.atname.to_string()));
        let atom_resname: Array1<String> = Array1::from_iter(atoms.iter().map(|at| at.resname.to_string()));
        // 每个原子对应的残基编号
        let resid_offset_receptor = receptor.models[0].atoms[0].resid;
        let resid_offset_ligand = ligand.models[0].atoms[0].resid;
        let mut atom_resid: Vec<usize> = receptor.models[0].atoms.iter().map(|at| (at.resid - resid_offset_receptor) as usize).collect();
        let resid_receptor_last = atom_resid.last().unwrap().to_owned();
        let resid_ligand = ligand.models[0].atoms.iter()
            .map(|at| (at.resid - resid_offset_ligand) as usize + resid_receptor_last + 1);
        atom_resid.extend(resid_ligand);
        let atom_resid: Array1<usize> = Array1::from_vec(atom_resid);

        AtomProperty {
            c6,
            c12,
            c10,
            atom_charge,
            atom_radius,
            atom_type_id,
            atom_id,
            atom_name,
            atom_resname,
            atom_resid
        }
    }
}