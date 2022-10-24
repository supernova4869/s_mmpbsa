use crate::index_parser::Index;
use xdrfile::*;
use crate::Parameters;
use ndarray::{ArrayBase, Dim, Ix, OwnedRepr};

pub fn do_mmpbsa_calculations(trj: &String, tpr: &String, ndx: &Index, use_dh: bool, use_ts: bool,
                              complex_grp: usize, receptor_grp: usize, ligand_grp: usize,
                              bt: f64, et: f64, dt: f64, params: &Parameters) {
    let complex = &ndx.groups[complex_grp].indexes;
    let receptor = &ndx.groups[receptor_grp].indexes;
    let ligand = &ndx.groups[ligand_grp].indexes;
    let gmx = &params.gmx;
    let apbs = &params.apbs;

    // 0. MMPBSA parameters
    let rad_type = params.rad_type;
    let rad_lj0 = params.rad_lj0;
    let mesh_type = params.mesh_type;
    let grid_type = params.grid_type;
    let cfac = params.cfac;
    let fadd = params.fadd;
    let df = params.df;

    // 1. 预处理轨迹: 复合物完整化, 团簇化, 居中叠合, 然后生成pdb文件
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let (coordinates, boxes) =
        get_atoms_trj(trj);   // frames x atoms(3x1)
    // pbc whole
    // println!("{:?}", coordinates);
    println!("{:?}", boxes);
    // ...
    // 2. 获取每个原子的电荷, 半径, LJ参数, 然后生成qrv文件
    // 3. Mpdb>pqr, 输出apbs, 计算MM, APBS
}

fn get_atoms_trj(trj: XTCTrajectory) -> (ArrayBase<OwnedRepr<[f32; 3]>, Dim<[Ix; 2]>>,
                                         ArrayBase<OwnedRepr<[f32; 3]>, Dim<[Ix; 2]>>) {
    let mut coord_matrix: Vec<[f32; 3]> = vec![];
    let mut box_size: Vec<[f32; 3]> = vec![];
    let mut num_frames = 0;
    let mut num_atoms = 0;
    for result in trj.into_iter() {
        let frame = result.unwrap();
        let atoms = frame.coords.to_vec();
        for a in atoms {
            coord_matrix.push(a);
        }
        for i in 0..3 {
            box_size.push(frame.box_vector[i]);
        }
        num_frames += 1;
        num_atoms = frame.num_atoms();
    }
    let coord_matrix = ndarray::Array::from_shape_vec(
        (num_frames, num_atoms), coord_matrix).unwrap();
    let box_size = ndarray::Array::from_shape_vec(
        (num_frames, 3), box_size).unwrap();
    return (coord_matrix, box_size);
}