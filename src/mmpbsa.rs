use std::fs;
use std::path::Path;
use crate::index_parser::Index;
use xdrfile::*;
use crate::Parameters;
use ndarray::{array, Array, ArrayBase, Dim, Ix, NdIndex, OwnedRepr, Shape};
use fs::File;
use std::cmp::max;
use std::io::Write;
use regex::{Match, Regex};
use std::str::FromStr;
use std::fmt::Debug;

pub fn do_mmpbsa_calculations(trj: &String, mdp: &String, ndx: &Index, wd: &Path, use_dh: bool, use_ts: bool,
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
    // pbc whole 先不写, 先默认按照已经消除了周期性来做后续处理, 之后再看周期性的事
    // 2. 获取每个原子的电荷, 半径, LJ参数, 然后生成qrv文件
    let qrv_path = "_Protein_in_water.qrv";
    let qrv_path = wd.join(qrv_path);
    let qrv_path = qrv_path.to_str().unwrap();
    gen_qrv(mdp, ndx, wd, receptor_grp, ligand_grp, qrv_path, rad_type, rad_lj0);
    // PB部分的参数从哪来?是否从前一菜单转移?
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
    let coord_matrix = Array::from_shape_vec(
        (num_frames, num_atoms), coord_matrix).unwrap();
    let box_size = Array::from_shape_vec(
        (num_frames, 3), box_size).unwrap();
    return (coord_matrix, box_size);
}

fn gen_qrv(mdp: &String, ndx: &Index, wd: &Path, receptor_grp: usize, ligand_grp: usize, qrv: &str, rad_type: i32, rad_lj0: f64) {
    // read mdp file
    let mut qrv = File::create(qrv).expect("Create qrv file failed");
    qrv.write_all("Protein Ligand\n".as_bytes()).expect("Writing qrv file failed");
    let ndx_rec = &ndx.groups[receptor_grp].indexes;
    let ndx_lig = &ndx.groups[ligand_grp].indexes;
    let mdp = fs::read_to_string(mdp).unwrap();
    let mdp: Vec<&str> = mdp.split("\n").collect();

    // get MD parameters
    let re = Regex::new(r"\s*ffparams:").unwrap();
    let locator = get_md_locators_all(&mdp, &re)[0];
    // number of atom types
    let re = Regex::new(r"\s*atnr=(\d+)").unwrap();
    let atnr = re.captures(mdp[locator + 1]).unwrap();
    let atnr = atnr.get(1).unwrap().as_str();
    qrv.write_all(format!("{}\n", atnr).as_bytes()).expect("Writing qrv file failed");
    let atnr: usize = atnr.parse().unwrap();
    // LJ parameters
    let mut locator = locator + 3;
    let mut sigma: Vec<f64> = vec![0.0; atnr];
    let mut epsilon: Vec<f64> = vec![0.0; atnr];
    let mut rad: Vec<f64> = vec![rad_lj0; atnr];

    // println!("Generating qrv file...");
    // for i in 0..atnr {
    //     qrv.write_all(format!("{:6}", i).as_bytes()).expect("Writing qrv file failed");
    //     // 获取每个原子的C6和C12
    //     for j in 0..atnr {
    //         let re = Regex::new(r".*c6\s*=\s*(.*),.*").unwrap();
    //         let re = Regex::new(r".*c12\s*=\s*(.*?)\s*$").unwrap();
    //         let c6 = re.captures(mdp[locator]).unwrap().get(1).unwrap().as_str();
    //         let c12 = re.captures(mdp[locator]).unwrap().get(1).unwrap().as_str();
    //         qrv.write_all(format!(" {} {}", c6, c12).as_bytes()).expect("Writing qrv file failed");
    //         let c6: f64 = c6.parse().unwrap();
    //         let c12: f64 = c12.parse().unwrap();
    //         // 计算每种原子类型的σ, ε, 半径
    //         if j == i {
    //             sigma[i] = 0.0;
    //             epsilon[i] = 0.0;
    //             rad[i] = rad_lj0;
    //             if c6 * c12 != 0.0 {
    //                 sigma[i] = 10.0 * (c12 / c6).powf(1.0 / 6.0); // 转换单位为A
    //                 epsilon[i] = c6.powi(2) / (4.0 * c12);
    //                 rad[i] = 0.5 * sigma[i]; // sigma为直径
    //             }
    //         }
    //     }
    //     locator += 1;
    //     qrv.write_all("\n".as_bytes()).expect("Writing qrv file failed");
    // }
    // println!("Finished generating qrv file.");

    // number of systems(?)
    let re = Regex::new(r"\s*#molblock\s*=\s*(.+?)\s*").unwrap();
    let mol_num: i32 = get_md_params_first(&mdp, &re).parse().unwrap();
    // number of molecule types
    let re = Regex::new(r"\s*moltype\s*=\s*(.+?)\s*").unwrap();
    let Imol: Vec<i32> = get_md_params_all(&mdp, &re);
    // number of molecules
    let re = Regex::new(r"\s*#molecules\s*=\s*(.+?)\s*").unwrap();
    let Nmol: Vec<i32> = get_md_params_all(&mdp, &re);
    // number of atoms
    let re = Regex::new(r"atom \((.+)\):").unwrap();
    let Natm: Vec<i32> = get_md_params_all(&mdp, &re);
    let maxNatm = Natm.iter().max().unwrap();
    // Name = Vector{String}()
    // Natm = Vector{Int32}()

    // locator of molecule types
    let re = Regex::new(r"\s*moltype.+\(").unwrap();
    let locators = get_md_locators_all(&mdp, &re);


    // 先看看ndarray再决定下面怎么写
    let arr: ArrayBase<OwnedRepr<[f32; 3]>, Dim<[Ix; 2]>> = Array::zeros((3, 100));
    println!("{:?}", arr);
    // let resID= Array::zeros();

    // Catm = zeros(Int32, length(Nmol), maxNatm)
    // Ratm = zeros(Float64, length(Nmol), maxNatm)
    // Satm = zeros(Float64, length(Nmol), maxNatm)
    // Eatm = zeros(Float64, length(Nmol), maxNatm)
    // Qatm = zeros(Float64, length(Nmol), maxNatm)
    // Tatm = Array{String}(undef, length(Nmol), maxNatm)

    // get atom parameters
}

fn get_md_locators_first(strings: &Vec<&str>, re: &Regex) -> usize {
    for (idx, l) in strings.into_iter().enumerate() {
        if re.is_match(l) {
            return idx;
        }
    }
    return 0;
}

fn get_md_locators_all(strings: &Vec<&str>, re: &Regex) -> Vec<usize> {
    let mut locators: Vec<usize> = vec![];
    for (idx, l) in strings.into_iter().enumerate() {
        if re.is_match(l) {
            locators.push(idx);
        }
    }
    return locators;
}

fn get_md_params_first(strings: &Vec<&str>, re: &Regex) -> String {
    for line in strings {
        if re.is_match(line) {
            let value = re.captures(line).unwrap();
            let value = value.get(1).unwrap().as_str();
            return value.to_string();
        }
    }
    return String::new();
}

fn get_md_params_all<T: Debug + FromStr>(strings: &Vec<&str>, re: &Regex) -> Vec<T> where <T as FromStr>::Err: Debug {
    let mut values: Vec<T> = vec![];
    for line in strings {
        if re.is_match(line) {
            let value = re.captures(line).unwrap();
            let value = value.get(1).unwrap().as_str();
            let value: T = value.trim().parse().unwrap();
            values.push(value);
        }
    }
    return values;
}