use crate::index_parser::Index;
use xdrfile::*;

pub fn do_mmpbsa_calculations(trj: &String, tpr: &String, ndx: &Index, use_dh: bool, use_ts: bool,
                          complex_grp: usize, receptor_grp: usize, ligand_grp: usize) {
    let complex = &ndx.groups[complex_grp].indexes;
    let receptor = &ndx.groups[receptor_grp].indexes;
    let ligand = &ndx.groups[ligand_grp].indexes;
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    // atoms: 三维数组, 原子坐标x帧号
    let atoms = get_atoms_trj(trj, ligand, 9501.0, 10000.0, 500.0);
    for i in 0..atoms.len() {
        println!("frame {}: {:?}\n", i, atoms[i]);
    }
    // 0. 定义apbs所需的参数
    // 1. 预处理轨迹: 复合物完整化, 团簇化, 居中叠合, 然后生成pdb文件
    // 2. 获取每个原子的电荷, 半径, LJ参数, 然后生成qrv文件
    // 3. Mpdb>pqr, 输出apbs, 计算MM, APBS
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