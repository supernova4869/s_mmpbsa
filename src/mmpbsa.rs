use std::fs;
use std::path::Path;
use crate::index_parser::Index;
use xdrfile::*;
use crate::Parameters;
use ndarray::{Array, Array1, Array2, Array3, ArrayBase, Dim, Ix, NdIndex, OwnedRepr, Shape, s};
use std::cmp::{max, min};
use std::collections::HashSet;
use std::f32::consts::E;
use std::fs::File;
use std::io::{stdin, Write};
use std::ops::{Add, Range};
use std::process::{Command, Termination};
use std::rc::Rc;
use indicatif::ProgressBar;
use regex::Regex;
use crate::gen_qrv::gen_qrv;

pub fn do_mmpbsa_calculations(trj: &String, mdp: &String, ndx: &Index, wd: &Path, use_dh: bool, use_ts: bool,
                              complex_grp: usize, receptor_grp: usize, ligand_grp: usize,
                              bt: f64, et: f64, dt: f64, settings: &Parameters) {
    // get charge, radius, LJ parameters of each atoms and generate qrv files
    let mut pid = String::from("_Protein_in_water");
    println!("Input system name (default: {}):", pid);
    let mut input = String::new();
    stdin().read_line(&mut input).expect("Error input");
    if input.trim().len() != 0 {
        pid = input.trim().to_string();
    }
    let qrv_path = String::from(pid.as_str()) + ".qrv";
    let qrv_path = wd.join(qrv_path);
    let qrv_path = qrv_path.to_str().unwrap();
    // gen_qrv(mdp, ndx, receptor_grp, ligand_grp, qrv_path, rad_type, rad_lj0);
    // Mpdb>pqr, 输出apbs, 计算MM, APBS
    do_mmpbsa(trj, mdp, ndx, qrv_path, wd, pid,
              complex_grp, receptor_grp, ligand_grp,
              bt, et, dt,
              settings, use_dh, use_ts);
}

fn get_atoms_trj(frames: &Vec<Rc<Frame>>) -> (Array3<f64>, Array3<f64>) {
    let num_frames = frames.len();
    let num_atoms = frames[0].num_atoms();
    let mut coord_matrix = Array3::<f64>::zeros((num_frames, num_atoms, 3));
    let mut box_size = Array3::<f64>::zeros((num_frames, 3, 3));
    for (idx, frame) in frames.into_iter().enumerate() {
        let atoms = frame.coords.to_vec();
        for (i, a) in atoms.into_iter().enumerate() {
            for j in 0..3 {
                coord_matrix[[idx, i, j]] = a[j] as f64;
            }
        }
        for (i, b) in frame.box_vector.into_iter().enumerate() {
            for j in 0..3 {
                box_size[[idx, i, j]] = b[j] as f64;
            }
        }
    }
    return (coord_matrix, box_size);
}

fn do_mmpbsa(trj: &String, mdp: &String, ndx: &Index, qrv: &str, wd: &Path, pid: String,
             complex_grp: usize, receptor_grp: usize, ligand_grp: usize,
             bt: f64, et: f64, dt: f64,
             settings: &Parameters, use_dh: bool, use_ts: bool) {
    // Running settings
    let gmx = &settings.gmx;
    let apbs = &settings.apbs;
    let rad_type = settings.rad_type;
    let rad_lj0 = settings.rad_lj0;
    let mesh_type = settings.mesh_type;
    let grid_type = settings.grid_type;
    let cfac = settings.cfac;
    let fadd = settings.fadd;
    let df = settings.df;

    // PB部分的参数从哪来?是否从前一菜单转移?
    let PBEset0 = "\
    \n  temp  298.15      # 温度\
    \n  pdie  2           # 溶质介电常数\
    \n  sdie  1       # 溶剂介电常数, 真空1, 水78.54\
    \n  \
    \n  npbe              # PB方程求解方法, lpbe(线性), npbe(非线性), smbpe(大小修正)\
    \n  bcfl  mdh         # 粗略格点PB方程的边界条件, zero, sdh/mdh(single/multiple Debye-Huckel), focus, map\
    \n  srfm  smol        # 构建介质和离子边界的模型, mol(分子表面), smol(平滑分子表面), spl2/4(三次样条/7阶多项式)\
    \n  chgm  spl4        # 电荷映射到格点的方法, spl0/2/4, 三线性插值, 立方/四次B样条离散\
    \n  swin  0.3         # 立方样条的窗口值, 仅用于 srfm=spl2/4\
    \n  \
    \n  srad  1.4         # 溶剂探测半径\
    \n  sdens 10          # 表面密度, 每A^2的格点数, (srad=0)或(srfm=spl2/4)时不使用\
    \n  \
    \n  ion charge  1 conc 0.15 radius 0.95  # 阳离子的电荷, 浓度, 半径\
    \n  ion charge -1 conc 0.15 radius 1.81  # 阴离子\
    \n  \
    \n  calcforce  no\
    \n  calcenergy comps";
    let PBEset = "
    \n  temp  298.15      # 温度\
    \n  pdie  2           # 溶质介电常数\
    \n  sdie  78.54       # 溶剂介电常数, 真空1, 水78.54\
    \n  \
    \n  npbe              # PB方程求解方法, lpbe(线性), npbe(非线性), smbpe(大小修正)\
    \n  bcfl  mdh         # 粗略格点PB方程的边界条件, zero, sdh/mdh(single/multiple Debye-Huckel), focus, map\
    \n  srfm  smol        # 构建介质和离子边界的模型, mol(分子表面), smol(平滑分子表面), spl2/4(三次样条/7阶多项式)\
    \n  chgm  spl4        # 电荷映射到格点的方法, spl0/2/4, 三线性插值, 立方/四次B样条离散\
    \n  swin  0.3         # 立方样条的窗口值, 仅用于 srfm=spl2/4\
    \n  \
    \n  srad  1.4         # 溶剂探测半径\
    \n  sdens 10          # 表面密度, 每A^2的格点数, (srad=0)或(srfm=spl2/4)时不使用\
    \n  \
    \n  ion charge  1 conc 0.15 radius 0.95  # 阳离子的电荷, 浓度, 半径\
    \n  ion charge -1 conc 0.15 radius 1.81  # 阴离子\
    \n  \
    \n  calcforce  no\
    \n  calcenergy comps";
    let PBAset = "
    \n  temp  298.15 # 温度\
    \n  srfm  sacc   # 构建溶剂相关表面或体积的模型\
    \n  swin  0.3    # 立方样条窗口(A), 用于定义样条表面\
    \n  \
    \n  # SASA\
    \n  srad  1.4    # 探测半径(A)\
    \n  gamma 1      # 表面张力(kJ/mol-A^2)\
    \n  \
    \n  #gamma const 0.0226778 3.84928\
    \n  #gamma const 0.027     0\
    \n  #gamma const 0.0301248 0         # AMBER-PB4 .0072*cal2J 表面张力, 常数\
    \n  \
    \n  press  0     # 压力(kJ/mol-A^3)\
    \n  bconc  0     # 溶剂本体密度(A^3)\
    \n  sdens 10\
    \n  dpos  0.2\
    \n  grid  0.1 0.1 0.1\
    \n  \
    \n  # SAV\
    \n  #srad  1.29      # SAV探测半径(A)\
    \n  #press 0.234304  # 压力(kJ/mol-A^3)\
    \n  \
    \n  # WCA\
    \n  #srad   1.25           # 探测半径(A)\
    \n  #sdens  200            # 表面的格点密度(1/A)\
    \n  #dpos   0.05           # 表面积导数的计算步长\
    \n  #bconc  0.033428       # 溶剂本体密度(A^3)\
    \n  #grid   0.45 0.45 0.45 # 算体积分时的格点间距(A)\
    \n  \
    \n  calcforce no\
    \n  calcenergy total";
    let qrv = wd.join(pid.to_string() + ".qrv");
    let wd = wd.join("temp");
    println!("Temporary files will be placed at {}", wd.display());
    if !wd.is_dir() {
        fs::create_dir(&wd).expect("Failed to create directory.");
    }
    let err = wd.join(pid.to_string() + ".err");
    let pdb = wd.join(pid.to_string() + ".pdb");

    // run MM-PBSA calculatons
    println!("Running MM-PBSA calculatons...");
    let qrv = fs::read_to_string(qrv).unwrap();
    let qrv: Vec<&str> = qrv.split("\n").collect();
    // c6 and c12
    let Atyp = qrv[1];
    let Atyp: usize = Atyp.trim().parse().unwrap();
    let mut C6 = Array2::<f64>::zeros((Atyp, Atyp));
    let mut C12 = Array2::<f64>::zeros((Atyp, Atyp));
    for i in 0..Atyp {
        for j in 0..Atyp {
            let paralj: Vec<&str> = qrv[i + 2].trim().split(" ").collect();
            let c6: f64 = paralj[2 * j + 1].parse().unwrap();
            C6[[i, j]] = c6;
            let c12: f64 = paralj[2 * j + 2].parse().unwrap();
            C12[[i, j]] = c12;
        }
    }

    let ndxCom = &ndx.groups[complex_grp].indexes;
    let ndxPro = &ndx.groups[receptor_grp].indexes;
    let ndxLig = &ndx.groups[ligand_grp].indexes;
    let mut Qatm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut Ratm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut Catm = Array1::<usize>::zeros(qrv.len() - 3 - Atyp);
    let mut Satm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut Eatm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut atm_index = Array1::<i32>::zeros(qrv.len() - 3 - Atyp);
    let mut atm_name = Array1::<String>::default(qrv.len() - 3 - Atyp);
    let mut atm_resname = Array1::<String>::default(qrv.len() - 3 - Atyp);
    let mut atm_resnum = Array1::<usize>::zeros(qrv.len() - 3 - Atyp);
    // resPro has blanks in the left part
    let mut resPro = Array1::<String>::default(qrv.len() - 3 - Atyp);
    // resLig has blanks in the right part
    let mut resLig = Array1::<String>::default(qrv.len() - 3 - Atyp);
    let mut idx = 0;
    for line in &qrv[Atyp + 2..qrv.len() - 1] { // there's a blank line at the end
        let line: Vec<&str> = line
            .split(" ")
            .filter_map(|p| match p.len() {
                0 => None,
                _ => Some(p)
            })
            .collect();
        Qatm[idx] = line[1].parse().unwrap();
        Ratm[idx] = line[2].parse().unwrap();
        Catm[idx] = line[3].parse().unwrap();
        Satm[idx] = line[4].parse().unwrap();
        Eatm[idx] = line[5].parse().unwrap();
        atm_index[idx] = line[0].parse().unwrap();
        atm_name[idx] = line[line.len() - 2].to_string();
        let res = line[line.len() - 3];
        atm_resnum[idx] = res[0..5].parse().unwrap();
        atm_resname[idx] = res[5..].to_string();

        if line[line.len() - 1] == "Rec" {
            resPro[idx] = "R~".to_string() + line[line.len() - 3];
        } else {
            resLig[idx] = "L~".to_string() + line[line.len() - 3];
        }
        idx += 1;
    }

    // fix resnum as start at 0 and no leap
    let diff = &atm_resnum.slice(s![1..atm_resnum.len() - 1]) - &atm_resnum.slice(s![0..atm_resnum.len() - 2]);
    let mut boundaries: Vec<usize> = vec![0];
    for (idx, res_num) in diff.into_iter().enumerate() {
        if res_num > 0 {
            boundaries.push(idx + 1);
        }
    }
    boundaries.push(atm_resnum.len());
    for b in 1..boundaries.len() {
        for i in boundaries[b - 1]..boundaries[b] {
            atm_resnum[i] = b - 1;
        }
    }

    // atom number of receptor and ligand
    let Npro = ndxPro.len();
    let Nlig = ndxLig.len();
    let Ncom = Npro + Nlig;
    let mut temp = 0.0;
    let mut pdie = 0.0;
    let mut sdie = 0.0;
    let mut Nion: usize = 0;
    let mut Qion = Array1::<f64>::zeros(10);
    let mut Cion = Array1::<f64>::zeros(10);

    // polar parameters
    for i in PBEset.split("\n") {
        if i.trim().len() != 0 {
            let paras: Vec<&str> = i
                .split(" ")
                .filter_map(|p| match p.len() {
                    0 => None,
                    _ => Some(p)
                })
                .collect();
            match paras[0] {
                "temp" => temp = paras[1].parse().unwrap(),
                "pdie" => pdie = paras[1].parse().unwrap(),
                "sdie" => sdie = paras[1].parse().unwrap(),
                "ion" => {
                    Nion += 1;
                    Qion[Nion] = paras[2].parse().unwrap();
                    Cion[Nion] = paras[4].parse().unwrap();
                }
                _ => ()
            }
        }
    }

    let mut gamma = 0.0;
    let mut _const = 0.0;
    // apolar parameters
    for i in PBAset.split("\n") {
        if i.trim().len() != 0 {
            let paras: Vec<&str> = i
                .split("\n")
                .filter_map(|p| match p.len() {
                    0 => None,
                    _ => Some(p)
                })
                .collect();
            if paras[0].ends_with("gamma") && paras[1] == "const" {
                gamma = paras[2].parse().unwrap();
                _const = paras[3].parse().unwrap();
            }
        }
    }
    println!("Finished setting parameters.");

    // 1. 预处理轨迹: 复合物完整化, 团簇化, 居中叠合, 然后生成pdb文件
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let frames: Vec<Rc<Frame>> = trj.into_iter().map(|p| p.unwrap()).collect();
    let total_frames = frames.len();
    // pbc whole 先不写, 先默认按照已经消除了周期性来做后续处理, 之后再看周期性的事
    let (coordinates, boxes) = get_atoms_trj(&frames);   // frames x atoms(3x1)

    // border of the whole molecule
    let minX = coordinates.slice(s![.., .., 0]).iter().
        fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    let maxX = coordinates.slice(s![.., .., 0]).iter().
        fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    let minY = coordinates.slice(s![.., .., 1]).iter().
        fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    let maxY = coordinates.slice(s![.., .., 1]).iter().
        fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    let minZ = coordinates.slice(s![.., .., 2]).iter().
        fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    let maxZ = coordinates.slice(s![.., .., 2]).iter().
        fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));

    let mut minXpro = Array1::<f64>::zeros(total_frames);
    let mut minYpro = Array1::<f64>::zeros(total_frames);
    let mut minZpro = Array1::<f64>::zeros(total_frames);
    let mut maxXpro = Array1::<f64>::zeros(total_frames);
    let mut maxYpro = Array1::<f64>::zeros(total_frames);
    let mut maxZpro = Array1::<f64>::zeros(total_frames);

    let mut minXlig = Array1::<f64>::zeros(total_frames);
    let mut minYlig = Array1::<f64>::zeros(total_frames);
    let mut minZlig = Array1::<f64>::zeros(total_frames);
    let mut maxXlig = Array1::<f64>::zeros(total_frames);
    let mut maxYlig = Array1::<f64>::zeros(total_frames);
    let mut maxZlig = Array1::<f64>::zeros(total_frames);

    let mut minXcom = Array1::<f64>::zeros(total_frames);
    let mut minYcom = Array1::<f64>::zeros(total_frames);
    let mut minZcom = Array1::<f64>::zeros(total_frames);
    let mut maxXcom = Array1::<f64>::zeros(total_frames);
    let mut maxYcom = Array1::<f64>::zeros(total_frames);
    let mut maxZcom = Array1::<f64>::zeros(total_frames);

    let bf = (bt / frames[1].time as f64) as usize;
    let ef = (et / frames[1].time as f64) as usize;
    let dframe = (dt / frames[1].time as f64) as usize;
    let total_frames = (ef - bf) / dframe + 1;
    let ef = ef + 1;        // range lefts the last frame

    // println!("Preparing qrv files...");
    // let pb = ProgressBar::new(total_frames as u64);
    // for cur_frm in (bf..ef).step_by(dframe) {
    //     // process each frame
    //
    //     let f_name = format!("{}_{}ns", pid, frames[cur_frm].time / 1000.0);
    //     let pqr_com = wd.join(format!("{}_com.pqr", f_name));
    //     let mut pqr_com = File::create(pqr_com).unwrap();
    //     let pqr_pro = wd.join(format!("{}_pro.pqr", f_name));
    //     let mut pqr_pro = File::create(pqr_pro).unwrap();
    //     let pqr_lig = wd.join(format!("{}_lig.pqr", f_name));
    //     let mut pqr_lig = File::create(pqr_lig).unwrap();
    //
    //     // get min and max atom coordinates
    //     let coordinates = coordinates.slice(s![cur_frm, .., ..]);
    //
    //     let atoms_rec_x_lb: Vec<f64> = ndxPro.iter().map(|&p| coordinates[[p, 0]] - Ratm[p]).collect();
    //     let atoms_rec_y_lb: Vec<f64> = ndxPro.iter().map(|&p| coordinates[[p, 1]] - Ratm[p]).collect();
    //     let atoms_rec_z_lb: Vec<f64> = ndxPro.iter().map(|&p| coordinates[[p, 2]] - Ratm[p]).collect();
    //     let atoms_rec_x_ub: Vec<f64> = ndxPro.iter().map(|&p| coordinates[[p, 0]] + Ratm[p]).collect();
    //     let atoms_rec_y_ub: Vec<f64> = ndxPro.iter().map(|&p| coordinates[[p, 1]] + Ratm[p]).collect();
    //     let atoms_rec_z_ub: Vec<f64> = ndxPro.iter().map(|&p| coordinates[[p, 2]] + Ratm[p]).collect();
    //     minXpro[cur_frm] = atoms_rec_x_lb.iter().
    //         fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    //     minYpro[cur_frm] = atoms_rec_y_lb.iter().
    //         fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    //     minZpro[cur_frm] = atoms_rec_z_lb.iter().
    //         fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    //     maxXpro[cur_frm] = atoms_rec_x_ub.iter().
    //         fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    //     maxYpro[cur_frm] = atoms_rec_y_ub.iter().
    //         fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    //     maxZpro[cur_frm] = atoms_rec_z_ub.iter().
    //         fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    //
    //     let atoms_lig_x_lb: Vec<f64> = ndxLig.iter().map(|&p| coordinates[[p, 0]] - Ratm[p]).collect();
    //     let atoms_lig_y_lb: Vec<f64> = ndxLig.iter().map(|&p| coordinates[[p, 1]] - Ratm[p]).collect();
    //     let atoms_lig_z_lb: Vec<f64> = ndxLig.iter().map(|&p| coordinates[[p, 2]] - Ratm[p]).collect();
    //     let atoms_lig_x_ub: Vec<f64> = ndxLig.iter().map(|&p| coordinates[[p, 0]] + Ratm[p]).collect();
    //     let atoms_lig_y_ub: Vec<f64> = ndxLig.iter().map(|&p| coordinates[[p, 1]] + Ratm[p]).collect();
    //     let atoms_lig_z_ub: Vec<f64> = ndxLig.iter().map(|&p| coordinates[[p, 2]] + Ratm[p]).collect();
    //     minXlig[cur_frm] = atoms_lig_x_lb.iter().
    //         fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    //     minYlig[cur_frm] = atoms_lig_y_lb.iter().
    //         fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    //     minZlig[cur_frm] = atoms_lig_z_lb.iter().
    //         fold(f64::INFINITY, |prev, curr| prev.min(*curr));
    //     maxXlig[cur_frm] = atoms_lig_x_ub.iter().
    //         fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    //     maxYlig[cur_frm] = atoms_lig_y_ub.iter().
    //         fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    //     maxZlig[cur_frm] = atoms_lig_z_ub.iter().
    //         fold(f64::NEG_INFINITY, |prev, curr| prev.max(*curr));
    //
    //     minXcom[cur_frm] = minXpro[cur_frm].min(minXpro[cur_frm]);
    //     minYcom[cur_frm] = minYpro[cur_frm].min(minYpro[cur_frm]);
    //     minZcom[cur_frm] = minZpro[cur_frm].min(minZpro[cur_frm]);
    //     maxXcom[cur_frm] = maxXpro[cur_frm].max(maxXpro[cur_frm]);
    //     maxYcom[cur_frm] = maxYpro[cur_frm].max(maxYpro[cur_frm]);
    //     maxZcom[cur_frm] = maxZpro[cur_frm].max(maxZpro[cur_frm]);
    //
    //     // loop atoms and write pqr information (from pqr)
    //     for at_id in ndxCom {
    //         let at_id = *at_id;
    //         let index = atm_index[at_id];
    //         let at_name = &atm_name[at_id];
    //         let resname = &atm_resname[at_id];
    //         let chain = "X";
    //         let resnum = atm_resnum[at_id];
    //         let coord = coordinates.slice(s![at_id, ..]);
    //         let x = coord[0];
    //         let y = coord[1];
    //         let z = coord[2];
    //         let q = Qatm[at_id];
    //         let r = Ratm[at_id];
    //         let atom_line = format!("ATOM  {:5} {:-4} {:3} X{:4}    {:8.3} {:8.3} {:8.3} \
    //         {:12.6} {:12.6}\n",
    //                                 index, at_name, resname, resnum, x, y, z, q, r);
    //
    //         // write qrv files
    //         pqr_com.write_all(atom_line.as_bytes()).unwrap();
    //         if ndxPro.contains(&at_id) {
    //             pqr_pro.write_all(atom_line.as_bytes()).unwrap();
    //         }
    //         if ndxLig.contains(&at_id) {
    //             pqr_lig.write_all(atom_line.as_bytes()).unwrap();
    //         }
    //     }
    //
    //     pb.inc(1);
    // }

    // calculate MM and PBSA
    let Iion = Cion.dot(&(&Qion * &Qion));
    let eps0 = 8.854187812800001e-12;
    let kb = 1.380649e-23;
    let Na = 6.02214076e+23;
    let qe = 1.602176634e-19;
    let RT2kJ = 8.314462618 * temp / 1e3;
    let kap = 1e-9 / f64::sqrt(eps0 * kb * temp * sdie / (Iion * qe * qe * Na * 1e3));

    let kJcou = 1389.35457520287;
    let Rcut = f64::INFINITY;

    let Nres = atm_resnum[atm_resnum.len() - 1] + 1;

    let dE = Array1::<f64>::zeros(Nres);
    let dGres = Array2::<f64>::zeros((Nres, 2));
    let dHres = Array2::<f64>::zeros((Nres, 2));
    let MMres = Array2::<f64>::zeros((Nres, 2));
    let COUres = Array2::<f64>::zeros((Nres, 2));
    let VDWres = Array1::<f64>::zeros(Nres);
    let dPBres = Array1::<f64>::zeros(Nres);
    let dSAres = Array1::<f64>::zeros(Nres);

    let vdw = Array1::<f64>::zeros(total_frames);
    let pb = Array1::<f64>::zeros(total_frames);
    let sa = Array1::<f64>::zeros(total_frames);
    let cou = Array2::<f64>::zeros((total_frames, 2));
    let mm = Array2::<f64>::zeros((total_frames, 2));
    let dh = Array2::<f64>::zeros((total_frames, 2));

    let PBres = Array1::<f64>::zeros(Nres);
    let SAres = Array1::<f64>::zeros(Nres);

    let Ipro = ndxPro[0];
    let Ilig = ndxLig[0];

    println!("Start MM/PB-SA calculations...");
    let pb = ProgressBar::new(total_frames as u64);
    pb.inc(0);
    for cur_frm in (bf..ef).step_by(dframe) {
        // // MM
        // let coord = coordinates.slice(s![cur_frm, .., ..]);
        // let mut dEcou = Array1::<f64>::zeros(Nres);
        // let mut dEvdw = Array1::<f64>::zeros(Nres);
        // // traverse receptor/ligand atoms to store parameters
        // for i in 0..Npro {
        //     let ii = i + Ipro;
        //     let qi = Qatm[ii];
        //     let ci = Catm[ii];
        //     let xi = coord[[ii, 0]];
        //     let yi = coord[[ii, 1]];
        //     let zi = coord[[ii, 2]];
        //     for j in 0..Nlig {
        //         let jj = j + Ilig;
        //         let qj = Qatm[jj];
        //         let cj = Catm[jj];
        //         let xj = coord[[jj, 0]];
        //         let yj = coord[[jj, 1]];
        //         let zj = coord[[jj, 2]];
        //         let r = f64::sqrt((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2));
        //         if r < Rcut {
        //             let mut Ecou = qi * qj / r / 10.0;
        //             if use_dh {
        //                 Ecou = Ecou * f64::exp(-kap * r);
        //             }
        //             let Evdw = C12[[ci, cj]] / r.powi(12) - C6[[ci, cj]] / r.powi(6);
        //             dEcou[atm_resnum[ii]] += Ecou;
        //             dEcou[atm_resnum[jj]] += Ecou;
        //             dEvdw[atm_resnum[ii]] += Evdw;
        //             dEvdw[atm_resnum[jj]] += Evdw;
        //         }
        //     }
        // }
        // for i in 0..Nres {
        //     dEcou[i] *= kJcou / (2.0 * pdie);
        //     dEvdw[i] /= 2.0;
        // }
        // let Evdw = dEvdw.sum();
        // let Ecou = dEcou.sum();
        //
        // // APBS
        let f_name = format!("{}_{}ns", pid, frames[cur_frm].time / 1000.0);
        // let mut input_apbs = File::create(wd.join(format!("{}.apbs", f_name))).unwrap();
        // input_apbs.write_all("read\n".as_bytes()).expect("Failed to write apbs input file.");
        // input_apbs.write_all(format!("  mol pqr {0}_com.pqr\n  \
        // mol pqr {0}_pro.pqr\n  \
        // mol pqr {0}_lig.pqr\n\
        // end\n\n", f_name).as_bytes())
        //     .expect("Failed to write apbs input file.");
        //
        // if mesh_type == 0 {
        //     // GMXPBSA
        //     input_apbs.write_all(dimAPBS(format!("{f_name}_com").as_str(), 1,
        //                                  minX, maxX, minY,
        //                                  maxY, minZ, maxZ,
        //                                  cfac, fadd, df, grid_type,
        //                                  PBEset, PBEset0, PBAset).as_bytes()).
        //         expect("Failed writing apbs file.");
        //     input_apbs.write_all(dimAPBS(format!("{f_name}_pro").as_str(), 2,
        //                                  minX, maxX, minY,
        //                                  maxY, minZ, maxZ,
        //                                  cfac, fadd, df, grid_type,
        //                                  PBEset, PBEset0, PBAset).as_bytes()).
        //         expect("Failed writing apbs file.");
        //     input_apbs.write_all(dimAPBS(format!("{f_name}_lig").as_str(), 3,
        //                                  minX, maxX, minY,
        //                                  maxY, minZ, maxZ,
        //                                  cfac, fadd, df, grid_type,
        //                                  PBEset, PBEset0, PBAset).as_bytes()).
        //         expect("Failed writing apbs file.");
        // } else if mesh_type == 1 {
        //     // g_mmpbsa
        //     input_apbs.write_all(dimAPBS(format!("{f_name}_com").as_str(), 1,
        //                                  minXcom[cur_frm], maxXcom[cur_frm], minYcom[cur_frm],
        //                                  maxYcom[cur_frm], minZcom[cur_frm], maxZcom[cur_frm],
        //                                  cfac, fadd, df, grid_type,
        //                                  PBEset, PBEset0, PBAset).as_bytes()).
        //         expect("Failed writing apbs file.");
        //     input_apbs.write_all(dimAPBS(format!("{f_name}_pro").as_str(), 2,
        //                                  minXcom[cur_frm], maxXcom[cur_frm], minYcom[cur_frm],
        //                                  maxYcom[cur_frm], minZcom[cur_frm], maxZcom[cur_frm],
        //                                  cfac, fadd, df, grid_type,
        //                                  PBEset, PBEset0, PBAset).as_bytes()).
        //         expect("Failed writing apbs file.");
        //     input_apbs.write_all(dimAPBS(format!("{f_name}_lig").as_str(), 3,
        //                                  minXcom[cur_frm], maxXcom[cur_frm], minYcom[cur_frm],
        //                                  maxYcom[cur_frm], minZcom[cur_frm], maxZcom[cur_frm],
        //                                  cfac, fadd, df, grid_type,
        //                                  PBEset, PBEset0, PBAset).as_bytes()).
        //         expect("Failed writing apbs file.");
        // }

        // let mut apbs_out = File::create(wd.join(format!("{}.out", f_name))).
        //     expect("Failed to create apbs out file");
        // let apbs_result = Command::new(apbs).
        //     arg(wd.join(format!("{}.apbs", f_name))).
        //     current_dir(&wd).output().expect("running apbs failed.");
        // let apbs_output = String::from_utf8(apbs_result.stdout).
        //     expect("Failed to get apbs output.");
        // apbs_out.write_all(apbs_output.as_bytes()).expect("Failed to write apbs output");

        // parse output
        // let apbs_out = fs::read_to_string(apbs_out).unwrap();
        let apbs_out = fs::read_to_string(wd.join(format!("{}.out", f_name))).unwrap();
        let Esol = Array2::<f64>::zeros((3, Ncom));
        let Evac = Array2::<f64>::zeros((3, Ncom));
        let Esas = Array2::<f64>::zeros((3, Ncom));
        let s: Vec<&str> = apbs_out
            .split("\n")
            .filter_map(|p|
                if p.trim().starts_with("CALCULATION ") ||
                    p.trim().starts_with("Per-atom energies:") ||
                    p.trim().starts_with("Atom") ||
                    p.trim().starts_with("Solvent Accessible Surface Area") ||
                    p.trim().starts_with("SASA") {
                    Some(p.trim())
                } else { None }
            )
            .collect();
        println!("{}", s.len());

        pb.inc(1);
    }
    pb.finish();
    println!("MM-PBSA calculation finished.");
}

fn dimAPBS(file: &str, Imol: i32, minX: f64, maxX: f64, minY: f64, maxY: f64, minZ: f64, maxZ: f64,
           cfac: f64, fadd: f64, df: f64, grid_type: i32, PBEset: &str, PBEset0: &str, PBAset: &str) -> String {
    let lenX = (maxX - minX).max(0.1);
    let cntX = (maxX + minX) / 2.0;
    let lenY = (maxY - minY).max(0.1);
    let cntY = (maxY + minY) / 2.0;
    let lenZ = (maxZ - minZ).max(0.1);
    let cntZ = (maxZ + minZ) / 2.0;
    let mut cX = lenX * cfac;
    let mut fX = cX.min(lenX + fadd);
    let mut cY = lenY * cfac;
    let mut fY = cY.min(lenY + fadd);
    let mut cZ = lenZ * cfac;
    let mut fZ = cZ.min(lenZ + fadd);

    let levN = 4;    // split level
    let t: i32 = 2_i32.pow(levN + 1);
    let mut nX: i32 = max((fX / df) as i32, 33);
    let mut nY: i32 = max((fY / df) as i32, 33);
    let mut nZ: i32 = max((fZ / df) as i32, 33);

    if grid_type == 0 { // GMXPBSA method
        let fpre = 1;
        let cfac = 1.7;
        fX = lenX + 2.0 * fadd;
        cX = fX * cfac;
        nX = t * ((fX / (t as f64 * df)).round() as i32 + 1 + fpre) + 1;
        fY = lenY + 2.0 * fadd;
        cY = fY * cfac;
        nY = t * ((fY / (t as f64 * df)).round() as i32 + 1 + fpre) + 1;
        fZ = lenZ + 2.0 * fadd;
        cZ = fZ * cfac;
        nZ = t * ((fZ / (t as f64 * df)).round() as i32 + 1 + fpre) + 1;
    }
    let MGset = "mg-auto";
    let mem = 200 * nX * nY * nZ / 1024 / 1024;       // MB

//		npX=nX; npY=nY; npZ=nZ
//		gmem=4000
//		ofrac=0.1
//		if(mem>=gmem) {
//			while(mem>gmem) {
//				maxN=max(npX, max(npY, npZ))
//					 if(maxN==npX) npX = t*((npX-1)/t-1)+1
//				else if(maxN==npY) npY = t*((npY-1)/t-1)+1
//				else if(maxN==npZ) npZ = t*((npZ-1)/t-1)+1
//				mem = 200*npX*npY*npZ/1024./1024
//			}

//			t=nX/npX; if(t>1) npX = int(t*(1+2*ofrac) + 1.0);
//			t=nY/npY; if(t>1) npY = int(t*(1+2*ofrac) + 1.0);
//			t=nZ/npZ; if(t>1) npZ = int(t*(1+2*ofrac) + 1.0);
//			MGset="mg-para\n  ofrac "ofrac"\n  pdime "npX" "npY" "npZ
//		}
    let XYZset = format!("  {MGset}\n  mol {Imol}\
        \n  dime   {nX:.0}  {nY:.0}  {nZ:.0}        # 格点数目, 所需内存: {mem} MB\
        \n  cglen  {cX:.3}  {cY:.3}  {cZ:.3}        # 粗略格点长度\
        \n  fglen  {fX:.3}  {fY:.3}  {fZ:.3}        # 细密格点长度\
        \n  fgcent {cntX:.3}  {cntY:.3}  {cntZ:.3}  # 细密格点中心\
        \n  cgcent {cntX:.3}  {cntY:.3}  {cntZ:.3}  # 粗略格点中心");

    return format!("\nELEC name {file}\n\
    {XYZset} \n\
    {PBEset} \n\
    end\n\n\
    ELEC name {file}_VAC\n\
    {XYZset}\n\
    {PBEset0} \n\
    end\n\n\
    APOLAR name {file}_SAS\n  \
    mol {Imol}\n{PBAset}\n\
    end\n\n\
    print elecEnergy {file} - {file}_VAC end\n\
    print apolEnergy {file}_SAS end\n\n");
}