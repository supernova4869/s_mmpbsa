use ndarray::{Array1, Array2};
use crate::parse_tpr::TPR;

pub struct AtomProperty {
    pub c6: Array2<f64>,
    pub c12: Array2<f64>,
    pub atm_charge: Array1::<f64>,
    pub atm_radius: Array1::<f64>,
    pub atm_typeindex: Array1<usize>,
    pub atm_sigma: Array1::<f64>,
    pub atm_epsilon: Array1::<f64>,
    pub atm_index: Array1<usize>,
    pub atm_name: Array1<String>,
    pub atm_resname: Array1<String>,
    pub atm_resnum: Array1<usize>,
}

impl AtomProperty {
    pub fn new(tpr: &TPR, ndx_com: &Vec<usize>) -> AtomProperty {
        // c6 and c12
        let mut c6: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        let mut c12: Array2<f64> = Array2::zeros((tpr.atom_types_num, tpr.atom_types_num));
        for i in 0..tpr.atom_types_num {
            for j in 0..tpr.atom_types_num {
                c6[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c6;
                c12[[i, j]] = tpr.lj_sr_params[i * tpr.atom_types_num + j].c12;
            }
        }

        let mut atm_charge: Array1::<f64> = Array1::zeros(ndx_com.len());
        let mut atm_radius: Array1::<f64> = Array1::zeros(ndx_com.len());
        let mut atm_typeindex: Array1<usize> = Array1::zeros(ndx_com.len());
        let mut atm_sigma: Array1::<f64> = Array1::zeros(ndx_com.len());
        let mut atm_epsilon: Array1::<f64> = Array1::zeros(ndx_com.len());
        let mut atm_index: Array1<usize> = Array1::zeros(ndx_com.len());
        let mut atm_name: Array1<String> = Array1::default(ndx_com.len());
        let mut atm_resname: Array1<String> = Array1::default(ndx_com.len());
        let mut atm_resnum: Array1<usize> = Array1::zeros(ndx_com.len());

        let mut idx = 0;
        let mut idx_com = 0;
        for mol in &tpr.molecules {
            for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                for atom in &mol.atoms {
                    if ndx_com.contains(&idx) {
                        atm_charge[idx_com] = atom.charge;
                        atm_radius[idx_com] = atom.radius;
                        atm_typeindex[idx_com] = atom.type_id;
                        atm_sigma[idx_com] = atom.sigma;
                        atm_epsilon[idx_com] = atom.epsilon;
                        atm_index[idx_com] = atom.id;
                        atm_name[idx_com] = atom.name.to_string();
                        atm_resname[idx_com] = mol.residues[atom.residue_index].name.to_string();
                        atm_resnum[idx_com] = atom.residue_index;
                        idx_com += 1;
                    }
                    idx += 1;
                }
            }
        }
        AtomProperty {
            c6,
            c12,
            atm_charge,
            atm_radius,
            atm_typeindex,
            atm_sigma,
            atm_epsilon,
            atm_index,
            atm_name,
            atm_resname,
            atm_resnum
        }
    }
}