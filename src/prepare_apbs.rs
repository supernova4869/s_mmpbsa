use std::collections::BTreeSet;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use ndarray::{Array1, ArrayView2};
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

pub fn write_apbs_input(ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>, coord: &ArrayView2<f64>,
                  atm_radius: &Array1<f64>, pbe_set: &PBESet, pba_set: &PBASet,
                  temp_dir: &PathBuf, f_name: &String, settings: &Settings) {
    let mut input_apbs = File::create(temp_dir.join(format!("{}.apbs", f_name))).unwrap();
    writeln!(input_apbs, "read").expect("Failed writing apbs file.");
    if ndx_lig.is_some() {
        writeln!(input_apbs, "  mol pqr {}_com.pqr", f_name).expect("Failed writing apbs file.");
        writeln!(input_apbs, "  mol pqr {}_rec.pqr", f_name).expect("Failed writing apbs file.");
        writeln!(input_apbs, "  mol pqr {}_lig.pqr", f_name).expect("Failed writing apbs file.");
    } else {
        writeln!(input_apbs, "  mol pqr {}_rec.pqr", f_name).expect("Failed writing apbs file.");
    }
    writeln!(input_apbs, "end\n").expect("Failed writing apbs file.");
    
    let (rec_box, lig_box, com_box) =
        gen_mesh_edges(ndx_rec, ndx_lig, coord, atm_radius);

    let mut pbe_set0 = PBESet::from(pbe_set);
    pbe_set0.sdie = 1.0;

    if let Some(com_box) = com_box {
        input_apbs.write_all(dim_apbs(format!("{}_com", f_name).as_str(), 1,
                                    com_box, settings, pbe_set, &pbe_set0, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
        input_apbs.write_all(dim_apbs(format!("{}_rec", f_name).as_str(), 2,
                                    rec_box, settings, pbe_set, &pbe_set0, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
        input_apbs.write_all(dim_apbs(format!("{}_lig", f_name).as_str(), 3,
                                    lig_box.unwrap(), settings, pbe_set, &pbe_set0, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
    } else {
        input_apbs.write_all(dim_apbs(format!("{}_rec", f_name).as_str(), 1,
                                    rec_box, settings, pbe_set, &pbe_set0, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
    }
}

fn get_lb(ndx: &BTreeSet<usize>, axis: usize, coord: &ArrayView2<f64>, atm_radius: &Array1<f64>) -> f64 {
    ndx.iter().map(|&p| coord[[p, axis]] - atm_radius[p])
        .fold(f64::INFINITY, |prev, curr| prev.min(curr))
}

fn get_ub(ndx: &BTreeSet<usize>, axis: usize, coord: &ArrayView2<f64>, atm_radius: &Array1<f64>) -> f64 {
    ndx.iter().map(|&p| coord[[p, axis]] + atm_radius[p])
        .fold(f64::NEG_INFINITY, |prev, curr| prev.max(curr))
}

pub fn gen_mesh_edges(ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>, coord: &ArrayView2<f64>,
                       atm_radius: &Array1<f64>) -> ([f64; 6], Option<[f64; 6]>, Option<[f64; 6]>) {
    let min_x_rec = get_lb(ndx_rec, 0, coord, atm_radius);
    let min_y_rec = get_lb(ndx_rec, 1, coord, atm_radius);
    let min_z_rec = get_lb(ndx_rec, 2, coord, atm_radius);
    let max_x_rec = get_ub(ndx_rec, 0, coord, atm_radius);
    let max_y_rec = get_ub(ndx_rec, 1, coord, atm_radius);
    let max_z_rec = get_ub(ndx_rec, 2, coord, atm_radius);

    let rec_box = [min_x_rec, min_y_rec, min_z_rec, max_x_rec, max_y_rec, max_z_rec];

    let (lig_box, com_box) = if let Some(ndx_lig) = ndx_lig {
        let min_x_lig = get_lb(ndx_lig, 0, coord, atm_radius);
        let min_y_lig = get_lb(ndx_lig, 1, coord, atm_radius);
        let min_z_lig = get_lb(ndx_lig, 2, coord, atm_radius);
        let max_x_lig = get_ub(ndx_lig, 0, coord, atm_radius);
        let max_y_lig = get_ub(ndx_lig, 1, coord, atm_radius);
        let max_z_lig = get_ub(ndx_lig, 2, coord, atm_radius);
        let min_x_com = min_x_rec.min(min_x_lig);
        let min_y_com = min_y_rec.min(min_y_lig);
        let min_z_com = min_z_rec.min(min_z_lig);
        let max_x_com = max_x_rec.max(max_x_lig);
        let max_y_com = max_y_rec.max(max_y_lig);
        let max_z_com = max_z_rec.max(max_z_lig);

        (Some([min_x_lig, min_y_lig, min_z_lig, max_x_lig, max_y_lig, max_z_lig]), 
        Some([min_x_com, min_y_com, min_z_com, max_x_com, max_y_com, max_z_com]))
    } else {
        (None, None)
    };

    return (rec_box, lig_box, com_box);
}

pub fn dim_apbs(file: &str, mol_index: i32, box_: [f64;6],
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

    // 格点数为32的倍数, apbs的特殊要求
    let t = 32.0;
    let n_x = ((f_x / df / t).round() * t) as i32 + 1;
    let n_y = ((f_y / df / t).round() * t) as i32 + 1;
    let n_z = ((f_z / df / t).round() * t) as i32 + 1;

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