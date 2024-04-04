use ndarray::{Array1, Array2};
use crate::parse_tpr::TPR;
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Clone)]
pub struct AtomProperty {
    pub c6: Array2<f64>,
    pub c12: Array2<f64>,
    pub atm_charge: Array1::<f64>,
    pub atm_radius: Array1::<f64>,
    pub atm_typeindex: Array1<usize>,
    pub atm_index: Array1<usize>,
    pub atm_name: Array1<String>,
    pub atm_resname: Array1<String>,
    pub atm_resid: Array1<usize>,
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
        let mut atm_index: Array1<usize> = Array1::zeros(ndx_com.len());
        let mut atm_name: Array1<String> = Array1::default(ndx_com.len());
        let mut atm_resname: Array1<String> = Array1::default(ndx_com.len());
        let mut atm_resid: Array1<usize> = Array1::zeros(ndx_com.len());

        let mut idx_total = 0;
        let mut idx = 0;
        let mut resind_offset = 0;      // residues number that has been overpast

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
                    if ndx_com.contains(&idx_total) {
                        atm_charge[idx] = atom.charge;
                        atm_radius[idx] = atom.radius;
                        atm_typeindex[idx] = atom.type_id;
                        atm_index[idx] = atom.id;
                        atm_name[idx] = atom.name.to_string();
                        atm_resname[idx] = mol.residues[atom.resind].name.to_string();
                        atm_resid[idx] = atom.resind + resind_offset;
                        idx += 1;
                    }
                    idx_total += 1;
                    pb.inc(1);
                    pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
                }
                resind_offset += mol.residues.len();
            }
        }

        pb.finish();

        AtomProperty {
            c6,
            c12,
            atm_charge,
            atm_radius,
            atm_typeindex,
            atm_index,
            atm_name,
            atm_resname,
            atm_resid
        }
    }
}