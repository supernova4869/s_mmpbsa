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
                  atom_radius: &Array1<f64>, pbe_set: &PBESet, pba_set: &PBASet,
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
        gen_mesh_edges(ndx_rec, ndx_lig, coord, atom_radius);

    let mut pbe_set_vacuum = PBESet::from(pbe_set);
    pbe_set_vacuum.sdie = 1.0;

    if let Some(com_box) = com_box {
        input_apbs.write_all(dim_apbs(format!("{}_com", f_name).as_str(), 1,
                                    com_box, settings, pbe_set, &pbe_set_vacuum, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
        input_apbs.write_all(dim_apbs(format!("{}_rec", f_name).as_str(), 2,
                                    rec_box, settings, pbe_set, &pbe_set_vacuum, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
        input_apbs.write_all(dim_apbs(format!("{}_lig", f_name).as_str(), 3,
                                    lig_box.unwrap(), settings, pbe_set, &pbe_set_vacuum, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
    } else {
        input_apbs.write_all(dim_apbs(format!("{}_rec", f_name).as_str(), 1,
                                    rec_box, settings, pbe_set, &pbe_set_vacuum, pba_set).as_bytes()).
            expect("Failed writing apbs file.");
    }
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

    let n_x = ((((f_x / df).round() - 1.0) / 32.0).round() * 32.0 + 1.0) as i32;
    let n_y = ((((f_y / df).round() - 1.0) / 32.0).round() * 32.0 + 1.0) as i32;
    let n_z = ((((f_z / df).round() - 1.0) / 32.0).round() * 32.0 + 1.0) as i32;

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