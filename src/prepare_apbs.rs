use std::collections::BTreeSet;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;
use std::path::Path;
use ndarray::{Array1, ArrayView2};
use apbs_generic::valist::Valist;
use crate::apbs_runner::MemMols;
use crate::parameters::*;
use crate::atom_property::AtomProperties;
use crate::settings::Settings;

pub fn prepare_pqr(cur_frm: usize, times: &Vec<f64>,
                   temp_dir: &Path, sys_name: &str, coord: &ArrayView2<f64>,
                   ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>,
                   aps: &AtomProperties) {
    let f_name = format!("{}_{}ns", sys_name, times[cur_frm]);
    let mut pqr_com = if ndx_lig.is_some() {
        Some(File::create(&temp_dir.join(format!("{}_com.pqr", f_name))).expect("Error: Failed to write pqr file"))
    } else {
        None
    };
    let mut pqr_lig = if ndx_lig.is_some() {
        Some(File::create(&temp_dir.join(format!("{}_lig.pqr", f_name))).expect("Error: Failed to write pqr file"))
    } else {
        None
    };
    let mut pqr_rec = File::create(&temp_dir.join(format!("{}_rec.pqr", f_name))).expect("Error: Failed to write pqr file");
    
    // loop atoms and write pqr information (from pqr)
    for at_id in 0..aps.atom_props.len() {
        let index = aps.atom_props[at_id].id;
        let at_name = &aps.atom_props[at_id].name;
        let resname = &aps.atom_props[at_id].resname;
        let resnum = aps.atom_props[at_id].resid;
        let x = coord[[at_id, 0]];
        let y = coord[[at_id, 1]];
        let z = coord[[at_id, 2]];
        let q = aps.atom_props[at_id].charge;
        let r = aps.atom_props[at_id].radius;
        let atom_line = format!("ATOM  {:5} {:-4} {:3} X {:3}    {:8.3} {:8.3} {:8.3} {:12.6} {:12.6}\n",
                                index, at_name, resname, resnum, x, y, z, q, r);

        // write qrv files
        // if has ligand
        if let Some(pqr_com) = &mut pqr_com {
            pqr_com.write_all(atom_line.as_bytes()).unwrap();
        }
        if let Some(pqr_lig) = &mut pqr_lig {
            if let Some(ndx_lig) = ndx_lig {
                if ndx_lig.contains(&at_id) {
                    pqr_lig.write_all(atom_line.as_bytes()).unwrap();
                }
            }
        }
        if ndx_rec.contains(&at_id) {
            pqr_rec.write_all(atom_line.as_bytes()).unwrap();
        }
    }
}

/// Builds the APBS input text for one frame in memory.
///
/// The returned text references molecules as `mol 1` (complex), `mol 2`
/// (receptor) and `mol 3` (ligand), matching [`build_molecules`].
pub fn build_apbs_input_text(ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>, coord: &ArrayView2<f64>,
                  atom_radius: &Array1<f64>, pbe_set: &PBESet, pba_set: &PBASet,
                  f_name: &str, settings: &Settings) -> String {
    let mut input_apbs = String::from("read\n");
    if ndx_lig.is_some() {
        input_apbs.push_str(&format!("  mol pqr {}_com.pqr\n", f_name));
        input_apbs.push_str(&format!("  mol pqr {}_rec.pqr\n", f_name));
        input_apbs.push_str(&format!("  mol pqr {}_lig.pqr\n", f_name));
    } else {
        input_apbs.push_str(&format!("  mol pqr {}_rec.pqr\n", f_name));
    }
    input_apbs.push_str("end\n\n");

    let (rec_box, lig_box, com_box) =
        gen_mesh_edges(ndx_rec, ndx_lig, coord, atom_radius);

    let mut pbe_set_vacuum = PBESet::from(pbe_set);
    pbe_set_vacuum.sdie = 1.0;

    let mut pba_set_modify = PBASet::from(pba_set);
    pba_set_modify.gamma = 1.0;
    pba_set_modify.bconc = 0.0;

    if let Some(com_box) = com_box {
        input_apbs.push_str(&prepare_apbs_content(format!("{}_com", f_name).as_str(), 1,
            com_box, settings, pbe_set, &pbe_set_vacuum, &pba_set_modify));
        input_apbs.push_str(&prepare_apbs_content(format!("{}_rec", f_name).as_str(), 2,
            rec_box, settings, pbe_set, &pbe_set_vacuum, &pba_set_modify));
        input_apbs.push_str(&prepare_apbs_content(format!("{}_lig", f_name).as_str(), 3,
            lig_box.unwrap(), settings, pbe_set, &pbe_set_vacuum, &pba_set_modify));
    } else {
        input_apbs.push_str(&prepare_apbs_content(format!("{}_rec", f_name).as_str(), 1,
            rec_box, settings, pbe_set, &pbe_set_vacuum, &pba_set_modify));
    }
    input_apbs
}

/// Builds the complex/receptor/ligand molecule lists that the APBS input
/// references as `mol 1/2/3`, straight from the in-memory coordinates.
///
/// Coordinates are in Angstrom, exactly what a PQR file would have carried,
/// so the solver sees identical input without the disk round-trip.
pub fn build_molecules(aps: &AtomProperties, coord: &ArrayView2<f64>,
                       ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>) -> MemMols {
    fn valist_from(aps: &AtomProperties, coord: &ArrayView2<f64>,
                   indices: impl Iterator<Item = usize>) -> Valist {
        let mut alist = Valist::new();
        for (i, at_id) in indices.enumerate() {
            let ap = &aps.atom_props[at_id];
            let mut atom = apbs_generic::vatom::Vatom::new();
            atom.set_position([coord[[at_id, 0]], coord[[at_id, 1]], coord[[at_id, 2]]]);
            atom.set_charge(ap.charge);
            atom.set_radius(ap.radius);
            atom.set_atom_name(&ap.name);
            atom.set_res_name(&ap.resname);
            atom.set_atom_id(i as i32);
            alist.atoms.push(atom);
        }
        alist.get_statistics();
        alist
    }

    let mut mols = MemMols::new();
    match ndx_lig {
        Some(ndx_lig) => {
            mols.insert("1".to_string(), Arc::new(valist_from(aps, coord, 0..aps.atom_props.len())));
            mols.insert("2".to_string(), Arc::new(valist_from(aps, coord, ndx_rec.iter().copied())));
            mols.insert("3".to_string(), Arc::new(valist_from(aps, coord, ndx_lig.iter().copied())));
        }
        None => {
            mols.insert("1".to_string(), Arc::new(valist_from(aps, coord, ndx_rec.iter().copied())));
        }
    }
    mols
}

fn get_bounds(ndx: &BTreeSet<usize>, coord: &ArrayView2<f64>, atom_radius: &Array1<f64>) -> [f64; 6] {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    
    for &p in ndx {
        let x = coord[[p, 0]];
        let y = coord[[p, 1]];
        let z = coord[[p, 2]];
        let r = atom_radius[p];
        
        // 左边界
        min_x = min_x.min(x - r);
        min_y = min_y.min(y - r);
        min_z = min_z.min(z - r);
        
        // 右边界
        max_x = max_x.max(x + r);
        max_y = max_y.max(y + r);
        max_z = max_z.max(z + r);
    }
    
    [min_x, min_y, min_z, max_x, max_y, max_z]
}

pub fn gen_mesh_edges(
    ndx_rec: &BTreeSet<usize>, 
    ndx_lig: &Option<BTreeSet<usize>>, 
    coord: &ArrayView2<f64>,
    atom_radius: &Array1<f64>
) -> ([f64; 6], Option<[f64; 6]>, Option<[f64; 6]>) {
    
    // 一次遍历计算受体的边界
    let rec_box = get_bounds(ndx_rec, coord, atom_radius);
    
    let (lig_box, com_box) = if let Some(ndx_lig) = ndx_lig {
        // 一次遍历计算配体的边界
        let lig_box = get_bounds(ndx_lig, coord, atom_radius);
        
        // 计算组合边界
        let com_box = [
            rec_box[0].min(lig_box[0]),  // min_x
            rec_box[1].min(lig_box[1]),  // min_y
            rec_box[2].min(lig_box[2]),  // min_z
            rec_box[3].max(lig_box[3]),  // max_x
            rec_box[4].max(lig_box[4]),  // max_y
            rec_box[5].max(lig_box[5]),  // max_z
        ];
        
        (Some(lig_box), Some(com_box))
    } else {
        (None, None)
    };
    
    (rec_box, lig_box, com_box)
}

pub fn prepare_apbs_content(file: &str, mol_index: i32, box_: [f64;6],
                settings: &Settings, pbe_set: &PBESet, pbe_set0: &PBESet, pba_set: &PBASet) -> String {
    let cfac = settings.cfac;
    let fadd = settings.fadd;
    let df = settings.df;

    let [min_x, min_y, min_z, max_x, max_y, max_z] = box_;

    let x_len = (max_x - min_x).max(0.1);
    let x_center = (max_x + min_x) / 2.0;
    let y_len = (max_y - min_y).max(0.1);
    let y_center = (max_y + min_y) / 2.0;
    let z_len = (max_z - min_z).max(0.1);
    let z_center = (max_z + min_z) / 2.0;

    let c_x = x_len * cfac as f64;
    let c_y = y_len * cfac as f64;
    let c_z = z_len * cfac as f64;
    let f_x = (x_len + fadd).min(c_x);
    let f_y = (y_len + fadd).min(c_y);
    let f_z = (z_len + fadd).min(c_z);

    let n_x = ((((f_x / df).round() - 1.0) / 32.0).round() * 32.0 + 1.0).max(33.0) as i32;
    let n_y = ((((f_y / df).round() - 1.0) / 32.0).round() * 32.0 + 1.0).max(33.0) as i32;
    let n_z = ((((f_z / df).round() - 1.0) / 32.0).round() * 32.0 + 1.0).max(33.0) as i32;

    let mg_set = "mg-auto";

    let xyz_set = format!("  {mg_set}\n  mol    {mol_index:7}\
        \n  dime   {n_x:7}  {n_y:7}  {n_z:7}\
        \n  cglen  {c_x:7.3}  {c_y:7.3}  {c_z:7.3}\
        \n  fglen  {f_x:7.3}  {f_y:7.3}  {f_z:7.3}\
        \n  fgcent {x_center:7.3}  {y_center:7.3}  {z_center:7.3}\
        \n  cgcent {x_center:7.3}  {y_center:7.3}  {z_center:7.3}\n");

    return format!("\nELEC name {}_SOL\n\
    {}\n\
    {}\n\
    end\n\n\
    ELEC name {}_VAC\n\
    {}\n\
    {}\n\
    end\n\n\
    APOLAR name {}_SAS\n  \
    mol    {:7}\n{}\n\
    end\n\n\
    print elecEnergy {}_SOL - {}_VAC end\n\
    print apolEnergy {}_SAS end\n\n", file, xyz_set, pbe_set.to_string(), file,
                   xyz_set, pbe_set0.to_string(), file, mol_index,
                   pba_set.to_string(), file, file, file);
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::apbs_runner::run_apbs_in_process_text;
    use ndarray::Array2;

    fn tiny_aps(charges: [f64; 3]) -> AtomProperties {
        let names = ["CA", "CB", "OW"];
        let resnames = ["ALA", "ALA", "SOL"];
        AtomProperties {
            c6: Array2::zeros((0, 0)),
            c12: Array2::zeros((0, 0)),
            at_map: std::collections::HashMap::new(),
            radius_type: "mBondi".to_string(),
            atom_props: names
                .iter()
                .enumerate()
                .map(|(i, &name)| crate::atom_property::AtomProperty {
                    charge: charges[i],
                    radius: 1.4,
                    type_id: 0,
                    id: i,
                    name: name.to_string(),
                    at_type: "C".to_string(),
                    resname: resnames[i].to_string(),
                    resid: i,
                })
                .collect(),
        }
    }

    /// The in-memory path (input text + in-memory molecules) must run the real
    /// solver end to end and yield per-atom results for every calculation of
    /// the complex/receptor/ligand decomposition.
    #[test]
    fn in_memory_input_runs_full_pbsa_decomposition() {
        let aps = tiny_aps([0.4, -0.2, -0.2]);
        let coord = Array2::from_shape_vec((3, 3), vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 9.0, 0.0, 0.0])
            .unwrap();
        let ndx_rec: BTreeSet<usize> = [0usize, 1].into_iter().collect();
        let ndx_lig: Option<BTreeSet<usize>> = Some([2usize].into_iter().collect());

        let mols = build_molecules(&aps, &coord.view(), &ndx_rec, &ndx_lig);
        assert_eq!(mols["1"].number_atoms(), 3, "complex = every atom");
        assert_eq!(mols["2"].number_atoms(), 2, "receptor subset");
        assert_eq!(mols["3"].number_atoms(), 1, "ligand subset");
        // The in-memory atoms carry the same coordinates and charges that the
        // PQR files used to hold.
        let a0 = mols["2"].get_atom(0);
        assert_eq!(a0.position, [0.0, 0.0, 0.0]);
        assert_eq!(a0.charge, 0.4);
        assert_eq!(a0.radius, 1.4);

        let pbe_set = PBESet::new(298.15);
        let pba_set = PBASet::new(298.15);
        let settings = Settings::new();
        let radius: Array1<f64> = Array1::from_iter(aps.atom_props.iter().map(|a| a.radius));
        let text = build_apbs_input_text(&ndx_rec, &ndx_lig, &coord.view(), &radius,
            &pbe_set, &pba_set, "t_0ns", &settings);
        assert!(text.contains("mol pqr t_0ns_com.pqr"));
        assert!(text.contains("ELEC name t_0ns_com_SOL"));
        assert!(text.contains("APOLAR name t_0ns_com_SAS"));

        let run = run_apbs_in_process_text(&text, &mols)
            .expect("in-memory PB/SA run must succeed");
        let names: Vec<&str> = run.calcs.iter().map(|c| c.name.as_str()).collect();
        for want in [
            "t_0ns_com_SOL", "t_0ns_com_VAC",
            "t_0ns_rec_SOL", "t_0ns_rec_VAC",
            "t_0ns_lig_SOL", "t_0ns_lig_VAC",
            "t_0ns_com_SAS", "t_0ns_rec_SAS", "t_0ns_lig_SAS",
        ] {
            assert!(names.contains(&want), "missing {want} in {names:?}");
        }
        for calc in &run.calcs {
            let expect = if calc.name.contains("_com_") {
                3
            } else if calc.name.contains("_rec_") {
                2
            } else {
                1
            };
            assert_eq!(calc.per_atom.len(), expect, "{}", calc.name);
            assert!(calc.per_atom.iter().all(|v| v.is_finite()), "{}", calc.name);
        }
    }

    /// Without a ligand only the receptor molecule and its SOL/VAC/SAS set
    /// exist.
    #[test]
    fn in_memory_input_without_ligand_has_receptor_only() {
        let aps = tiny_aps([0.4, -0.2, -0.2]);
        let coord = Array2::from_shape_vec((3, 3), vec![0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 9.0, 0.0, 0.0])
            .unwrap();
        let ndx_rec: BTreeSet<usize> = (0..3).collect();

        let mols = build_molecules(&aps, &coord.view(), &ndx_rec, &None);
        assert_eq!(mols.len(), 1);
        assert_eq!(mols["1"].number_atoms(), 3);

        let pbe_set = PBESet::new(298.15);
        let pba_set = PBASet::new(298.15);
        let settings = Settings::new();
        let radius: Array1<f64> = Array1::from_iter(aps.atom_props.iter().map(|a| a.radius));
        let text = build_apbs_input_text(&ndx_rec, &None, &coord.view(), &radius,
            &pbe_set, &pba_set, "t_0ns", &settings);
        assert!(text.contains("mol pqr t_0ns_rec.pqr"));
        assert!(!text.contains("_com"));
        assert!(!text.contains("_lig"));

        let run = run_apbs_in_process_text(&text, &mols)
            .expect("in-memory PB/SA run must succeed");
        let names: Vec<&str> = run.calcs.iter().map(|c| c.name.as_str()).collect();
        for want in ["t_0ns_rec_SOL", "t_0ns_rec_VAC", "t_0ns_rec_SAS"] {
            assert!(names.contains(&want), "missing {want} in {names:?}");
        }
        assert_eq!(names.len(), 3, "no complex/ligand calcs expected: {names:?}");
    }
}
