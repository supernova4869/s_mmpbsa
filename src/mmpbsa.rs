use crate::index_parser::Index;
use xdrfile::*;

pub fn do_mmpbsa_calculations(trj: &String, tpr: &String, ndx: &Index, use_dh: bool, use_ts: bool,
                          complex_grp: usize, receptor_grp: usize, ligand_grp: usize) {
    let complex = &ndx.groups[complex_grp].indexes;
    let receptor = &ndx.groups[receptor_grp].indexes;
    let ligand = &ndx.groups[ligand_grp].indexes;
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let atoms = get_atoms_trj(trj, ligand, 9501.0, 10000.0, 500.0);
    for i in 0..atoms.len() {
        println!("frame {}: {:?}\n", i, atoms[i]);
    }
}

fn get_atoms_trj(trj: XTCTrajectory, at_ids: &Vec<i32>, begin: f64, end: f64, dt: f64) -> Vec<Vec<[f32; 3]>> {
    let mut coord_matrix:Vec<Vec<[f32; 3]>> = vec![];
    let mut start_t = 0.0;
    for result in trj.into_iter() {
        let frame = result.unwrap();
        let ts: f64 = frame.time as f64;
        if ts > begin && start_t == 0.0 {
            start_t = ts;
        }
        let mut c: Vec<[f32; 3]> = vec![[0.0; 3]; at_ids.len()];
        if ts >= begin && ts <= end && (ts - start_t) % dt == 0.0 {
            for i in 0..at_ids.len() {
                c[i] = frame.coords[at_ids[i] as usize];
            }
            coord_matrix.push(c);
        }
    }
    return coord_matrix;
}