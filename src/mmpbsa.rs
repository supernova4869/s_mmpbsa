use crate::index_parser::Index;
use xdrfile::*;

pub fn do_mmpbsa_calculations(trj: &String, tpr: &String, ndx: &Index, use_dh: bool, use_ts: bool,
                          complex_grp: usize, receptor_grp: usize, ligand_grp: usize) {
    let complex = &ndx.groups[complex_grp].indexes;
    let receptor = &ndx.groups[receptor_grp].indexes;
    let ligand = &ndx.groups[ligand_grp].indexes;
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let atoms = get_atoms_trj(trj, ligand, 0.0, 0.0);
    for i in 0..atoms.len() {
        println!("frame {}: {:?}", i, atoms[i]);
    }
}

fn get_atoms_trj(trj: XTCTrajectory, at_ids: &Vec<i32>, begin: f64, end: f64) -> Vec<Vec<[f32; 3]>> {
    let mut coords:Vec<Vec<[f32; 3]>> = vec![];
    for (idx, result) in trj.into_iter().enumerate() {
        // 增加时间限制
        let mut c: Vec<[f32; 3]> = vec![[0.0; 3]; at_ids.len()];
        let frame = result.unwrap();
        for i in 0..at_ids.len() {
            c[i] = frame.coords[at_ids[i] as usize];
        }
        coords.push(c);
    }
    return coords;
}