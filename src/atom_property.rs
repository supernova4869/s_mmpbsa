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

        let mut idx_total = 0;
        let mut idx = 0;
        let mut resind_offset = 0;      // residues number that has been overpast
        for mol in &tpr.molecules {
            for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
                for atom in &mol.atoms {
                    if ndx_com.contains(&idx_total) {
                        atm_charge[idx] = atom.charge;
                        atm_radius[idx] = atom.radius;
                        atm_typeindex[idx] = atom.type_id;
                        atm_sigma[idx] = atom.sigma;
                        atm_epsilon[idx] = atom.epsilon;
                        atm_index[idx] = atom.id;
                        atm_name[idx] = atom.name.to_string();
                        atm_resname[idx] = mol.residues[atom.residue_index].name.to_string();
                        atm_resnum[idx] = atom.residue_index + resind_offset;
                        idx += 1;
                    }
                    idx_total += 1;
                }
                resind_offset += mol.residues.len();
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