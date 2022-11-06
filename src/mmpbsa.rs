use std::fs;
use std::path::Path;
use crate::index_parser::Index;
use xdrfile::*;
use crate::Parameters;
use ndarray::{Array, Array1, Array2, Array3, ArrayBase, Dim, Ix, NdIndex, OwnedRepr, Shape};
use std::cmp::max;
use std::fs::File;
use std::io::stdin;
use std::ops::{Add, Range};
use std::rc::Rc;
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
    let rad_type = &settings.rad_type;
    let rad_lj0 = &settings.rad_lj0;
    let mesh_type = &settings.mesh_type;
    let grid_type = &settings.grid_type;
    let cfac = &settings.cfac;
    let fadd = &settings.fadd;
    let df = &settings.df;

    // PB部分的参数从哪来?是否从前一菜单转移?
    let PBEset = "
temp  298.15      # 温度
pdie  2           # 溶质介电常数
sdie  78.54       # 溶剂介电常数, 真空1, 水78.54

npbe              # PB方程求解方法, lpbe(线性), npbe(非线性), smbpe(大小修正)
bcfl  mdh         # 粗略格点PB方程的边界条件, zero, sdh/mdh(single/multiple Debye-Huckel), focus, map
srfm  smol        # 构建介质和离子边界的模型, mol(分子表面), smol(平滑分子表面), spl2/4(三次样条/7阶多项式)
chgm  spl4        # 电荷映射到格点的方法, spl0/2/4, 三线性插值, 立方/四次B样条离散
swin  0.3         # 立方样条的窗口值, 仅用于 srfm=spl2/4

srad  1.4         # 溶剂探测半径
sdens 10          # 表面密度, 每A^2的格点数, (srad=0)或(srfm=spl2/4)时不使用

ion charge  1 conc 0.15 radius 0.95  # 阳离子的电荷, 浓度, 半径
ion charge -1 conc 0.15 radius 1.81  # 阴离子

calcforce  no
calcenergy comps";
    let PBAset = "
temp  298.15 # 温度
srfm  sacc   # 构建溶剂相关表面或体积的模型
swin  0.3    # 立方样条窗口(A), 用于定义样条表面

# SASA
srad  1.4    # 探测半径(A)
gamma 1      # 表面张力(kJ/mol-A^2)

#gamma const 0.0226778 3.84928
#gamma const 0.027     0
#gamma const 0.0301248 0         # AMBER-PB4 .0072*cal2J 表面张力, 常数

press  0     # 压力(kJ/mol-A^3)
bconc  0     # 溶剂本体密度(A^3)
sdens 10
dpos  0.2
grid  0.1 0.1 0.1

# SAV
#srad  1.29      # SAV探测半径(A)
#press 0.234304  # 压力(kJ/mol-A^3)

# WCA
#srad   1.25           # 探测半径(A)
#sdens  200            # 表面的格点密度(1/A)
#dpos   0.05           # 表面积导数的计算步长
#bconc  0.033428       # 溶剂本体密度(A^3)
#grid   0.45 0.45 0.45 # 算体积分时的格点间距(A)

calcforce no
calcenergy total";
    let qrv = wd.join(pid.to_string() + ".qrv");
    let wd = wd.join("temp");
    println!("Temporary files will be placed at {}", wd.display());
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
            let c6: f64 = paralj[2 * j].parse().unwrap();
            C6[[i, j]] = c6;
            let c12: f64 = paralj[2 * j + 1].parse().unwrap();
            C12[[i, j]] = c12;
        }
    }

    let ndxCom = &ndx.groups[complex_grp].indexes;
    let ndxPro = &ndx.groups[receptor_grp].indexes;
    let ndxLig = &ndx.groups[ligand_grp].indexes;
    let mut Qatm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut Ratm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut Catm = Array1::<i32>::zeros(qrv.len() - 3 - Atyp);
    let mut Satm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    let mut Eatm = Array1::<f64>::zeros(qrv.len() - 3 - Atyp);
    // 考虑将这两个合并成一个
    let mut resPro = Array1::<String>::default(qrv.len() - 3 - Atyp);
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

        if line[line.len() - 1] == "Rec" {
            resPro[idx] = "R~".to_string() + line[line.len() - 3];
        } else {
            resLig[idx] = "L~".to_string() + line[line.len() - 3];
        }
        idx += 1;
    }

    // atom number of receptor and ligand
    let Npro = resPro.len();
    let Nlig = resLig.len();
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

    let Iion = &Cion.dot(&(&Qion * &Qion));
    let eps0 = 8.854187812800001e-12;
    let kb = 1.380649e-23;
    let Na = 6.02214076e+23;
    let qe = 1.602176634e-19;
    let RT2kJ = 8.314462618 * temp / 1e3;
    let kap = 1e-10 / f64::sqrt(eps0 * kb * temp * sdie / (Iion * qe * qe * Na * 1e3));

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
    // pbc whole 先不写, 先默认按照已经消除了周期性来做后续处理, 之后再看周期性的事
    let (coordinates, boxes) = get_atoms_trj(&frames);   // frames x atoms(3x1)

    let mut time = bt;
    let mut Nfrm: usize = 0;
    // 有空用map/filter优化一下
    for f in frames.into_iter() {
        if f.time as f64 == bt {
            break
        }
        Nfrm += 1;
    }

    while time <= et {
        // process each frame
        let f_name = format!("{}~{}ns", pid, time / 1000.0);
        println!("{}", f_name);
        let minXpro = 0.0;
        let minXlig = 0.0;
        let minYpro = 0.0;
        let minYlig = 0.0;
        let minZpro = 0.0;
        let minZlig = 0.0;

        let maxXpro = 0.0;
        let maxXlig = 0.0;
        let maxYpro = 0.0;
        let maxYlig = 0.0;
        let maxZpro = 0.0;
        let maxZlig = 0.0;

        let minXcom = 0.0;
        let minYcom = 0.0;
        let minZcom = 0.0;
        let maxXcom = 0.0;
        let maxYcom = 0.0;
        let maxZcom = 0.0;

        let pqr_com = wd.join(format!("{}_com.pqr", f_name));
        let pqr_pro = wd.join(format!("{}_pro.pqr", f_name));
        let pqr_lig = wd.join(format!("{}_lig.pqr", f_name));

        // read atom coordinates
        let mut coordinates_rec = Array2::<f64>::zeros((ndxPro.len(), 3));
        let mut coordinates_lig = Array2::<f64>::zeros((ndxLig.len(), 3));
        for (idx, at_id) in ndxPro.into_iter().enumerate() {
            for j in 0..3 {
                coordinates_rec[[idx, j]] = coordinates[[Nfrm, *at_id, j]];
            }
        }
        for (idx, at_id) in ndxLig.into_iter().enumerate() {
            for j in 0..3 {
                coordinates_lig[[idx, j]] = coordinates[[Nfrm, *at_id, j]];
            }
        }
        println!("{}", coordinates_rec);
        println!("{}", coordinates_lig);
        // get min and max atom coordinates
        // loop atoms and write pqr information (from pqr)
        // write qrv files
        // calculate MM
        // calculate PBSA by apbs and extract informations
        // write output

        Nfrm += 1;
        time += dt;
    }
}