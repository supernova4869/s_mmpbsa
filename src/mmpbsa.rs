use std::cmp::Ordering;
use std::fs;
use std::path::Path;
use crate::index_parser::Index;
use xdrfile::*;
use crate::parameters::Parameters;
use ndarray::{Array1, Array2, Array3, s};
use std::fs::File;
use std::io::{stdin, Write};
use std::process::Command;
use std::rc::Rc;
use indicatif::{ProgressBar, ProgressStyle};
use crate::parse_tpr::TPR;
use crate::apbs_param::{PBASet, PBESet};
use crate::atom_property::AtomProperty;
use crate::prepare_apbs::{prepare_apbs_inputs, write_apbs};

pub fn do_mmpbsa_calculations(trj: &String, tpr: &TPR, ndx: &Index, wd: &Path, sys_name: &String,
                              complex_grp: usize, receptor_grp: usize, ligand_grp: usize,
                              bt: f64, et: f64, dt: f64,
                              pbe_set: &PBESet, pba_set: &PBASet, settings: &Parameters)
                              -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    // pdb>pqr, output apbs, calculate MM, calculate APBS

    // Running settings
    let apbs = &settings.apbs;
    let mesh_type = settings.mesh_type;
    let use_dh = settings.use_dh;
    let use_ts = settings.use_ts;

    let temp_dir = wd.join(sys_name);
    println!("Temporary files will be placed at {}/", temp_dir.display());
    if !temp_dir.is_dir() {
        fs::create_dir(&temp_dir).expect(format!("Failed to create temp directory: {}.", sys_name).as_str());
    } else {
        println!("Directory {}/ not empty. Clear? [Y/n]", temp_dir.display());
        let mut input = String::from("");
        stdin().read_line(&mut input).expect("Get input error");
        if input.trim().len() == 0 || input.trim() == "Y" || input.trim() == "y" {
            fs::remove_dir_all(&temp_dir).expect("Remove dir failed");
            fs::create_dir(&temp_dir).expect(format!("Failed to create temp directory: {}.", sys_name).as_str());
        }
    }

    // run MM/PB-SA calculations
    println!("Running MM/PB-SA calculations...");
    println!("Preparing parameters...");

    // c6 and c12
    let atom_types_num = tpr.atom_types_num;
    let mut c6: Array2<f64> = Array2::zeros((atom_types_num, atom_types_num));
    let mut c12: Array2<f64> = Array2::zeros((atom_types_num, atom_types_num));
    for i in 0..atom_types_num {
        for j in 0..atom_types_num {
            c6[[i, j]] = tpr.lj_sr_params[i * atom_types_num + j].c6;
            c12[[i, j]] = tpr.lj_sr_params[i * atom_types_num + j].c12;
        }
    }

    // atom indexes
    let ndx_com = &ndx.groups[complex_grp].indexes;
    let ndx_rec = &ndx.groups[receptor_grp].indexes;
    let ndx_lig = &ndx.groups[ligand_grp].indexes;

    // atom number of receptor and ligand
    let atom_num_rec = ndx_rec.len();
    let atom_num_lig = ndx_lig.len();
    let atom_num_com = atom_num_rec + atom_num_lig;

    // atom properties
    let aps = AtomProperty::new(tpr, ndx_com);

    // normalize residue indexes
    let mut ndx_lig: Vec<usize> = ndx_lig.iter().map(|p| p - ndx_com[0]).collect();
    let mut ndx_rec: Vec<usize> = ndx_rec.iter().map(|p| p - ndx_com[0]).collect();
    let ndx_com: Vec<usize> = ndx_com.iter().map(|p| p - ndx_com[0]).collect();
    if ndx_lig[0] > ndx_rec[0] {
        ndx_lig = ndx_lig.iter().map(|p| p - ndx_lig[0] + ndx_rec.len()).collect();
    } else {
        ndx_rec = ndx_rec.iter().map(|p| p - ndx_rec[0] + ndx_lig.len()).collect();
    }

    // 1. pre-treat trajectory: fix pbc
    println!("Reading trajectory file...");
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let frames: Vec<Rc<Frame>> = trj.into_iter().map(|p| p.unwrap()).collect();
    // pbc whole 先不写, 先默认按照已经消除了周期性来做后续处理, 之后再看周期性的事
    println!("Extracting atoms coordination...");
    let (coordinates, boxes) = get_atoms_trj(&frames);   // frames x atoms(3x1)

    let (bf, ef, dframe, total_frames) = get_frames_range(&frames, bt, et, dt);

    prepare_apbs_inputs(&frames, bf, ef, dframe, total_frames, &temp_dir,
                        sys_name, &coordinates, &ndx_com, &ndx_rec, &ndx_lig, &aps);

    // calculate MM and PBSA
    let eps0 = 8.854187812800001e-12;
    let kb = 1.380649e-23;
    let NA = 6.02214076e+23;
    let qe = 1.602176634e-19;
    let RT2kJ = 8.314462618 * pbe_set.temp / 1e3;
    let ion_strength: f64 = pbe_set.ions.iter()
        .map(|ion| ion.charge * ion.charge * ion.conc).sum();
    let kap = 1e-9 / f64::sqrt(eps0 * kb * pbe_set.temp * pbe_set.sdie / (ion_strength * qe * qe * NA * 1e3));

    let kJcou = 1389.35457520287;
    let Rcut = f64::INFINITY;

    // default gamma for apbs calculation is 1
    let gamma = 0.0301248;      // Here is the surface extension constant from AMBER-PB4
    let _const = 0.0;

    let total_res_num = tpr.molecules.iter().map(|mol|
        mol.residues.len() * tpr.molecule_types[mol.molecule_type_id].molecules_num as usize).sum();

    let mut dE: Array1<f64> = Array1::zeros(total_res_num);
    let mut dGres: Array1<f64> = Array1::zeros(total_res_num);
    let mut dHres: Array1<f64> = Array1::zeros(total_res_num);
    let mut MMres: Array1<f64> = Array1::zeros(total_res_num);
    let mut COUres: Array1<f64> = Array1::zeros(total_res_num);
    let mut VDWres: Array1<f64> = Array1::zeros(total_res_num);
    let mut dPBres: Array1<f64> = Array1::zeros(total_res_num);
    let mut dSAres: Array1<f64> = Array1::zeros(total_res_num);

    let mut vdw: Array1<f64> = Array1::zeros(total_frames);
    let mut pb: Array1<f64> = Array1::zeros(total_frames);
    let mut sa: Array1<f64> = Array1::zeros(total_frames);
    let mut cou: Array1<f64> = Array1::zeros(total_frames);
    let mut mm: Array1<f64> = Array1::zeros(total_frames);
    let mut dh: Array1<f64> = Array1::zeros(total_frames);

    let mut pb_res: Array1<f64> = Array1::zeros(total_res_num);
    let mut sa_res: Array1<f64> = Array1::zeros(total_res_num);

    println!("Start MM/PB-SA calculations...");
    let pgb = ProgressBar::new(total_frames as u64);
    set_style(&pgb);
    pgb.inc(0);
    let mut idx = 0;
    for cur_frm in (bf..ef + 1).step_by(dframe) {
        // MM
        let coord = coordinates.slice(s![cur_frm, .., ..]);
        let mut de_cou: Array1<f64> = Array1::zeros(total_res_num);
        let mut de_vdw: Array1<f64> = Array1::zeros(total_res_num);
        for &i in &ndx_rec {
            let qi = aps.atm_charge[i];
            let ci = aps.atm_typeindex[i];
            let xi = coord[[i, 0]];
            let yi = coord[[i, 1]];
            let zi = coord[[i, 2]];
            for &j in &ndx_lig {
                let qj = aps.atm_charge[j];
                let cj = aps.atm_typeindex[j];
                let xj = coord[[j, 0]];
                let yj = coord[[j, 1]];
                let zj = coord[[j, 2]];
                let r = f64::sqrt((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2));
                if r < Rcut {
                    let mut e_cou = qi * qj / r / 10.0;
                    if use_dh {
                        e_cou = e_cou * f64::exp(-kap * r);
                    }
                    let e_vdw = c12[[ci, cj]] / r.powi(12) - c6[[ci, cj]] / r.powi(6);
                    de_cou[aps.atm_resnum[i]] += e_cou;
                    de_cou[aps.atm_resnum[j]] += e_cou;
                    de_vdw[aps.atm_resnum[i]] += e_vdw;
                    de_vdw[aps.atm_resnum[j]] += e_vdw;
                }
            }
        }
        for i in 0..total_res_num {
            de_cou[i] *= kJcou / (2.0 * pbe_set.pdie);
            de_vdw[i] /= 2.0;
        }
        let e_vdw = de_vdw.sum();
        let e_cou = de_cou.sum();

        // APBS
        let f_name = format!("{}_{}ns", sys_name, frames[cur_frm].time / 1000.0);
        write_apbs(&ndx_rec, &ndx_lig, &coord, &aps.atm_radius,
                   pbe_set, &pba_set, &temp_dir, &f_name, settings);
        // invoke apbs program to do apbs calculations
        if !apbs.is_empty() {
            let mut apbs_out = File::create(temp_dir.join(format!("{}.out", f_name))).
                expect("Failed to create apbs out file");
            let apbs_result = Command::new(apbs).
                arg(format!("{}.apbs", f_name)).
                current_dir(&temp_dir).output().expect("running apbs failed.");
            let apbs_output = String::from_utf8(apbs_result.stdout).
                expect("Failed to get apbs output.");
            apbs_out.write_all(apbs_output.as_bytes()).expect("Failed to write apbs output");
        } else {
            println!("Warning: APBS not found. Will not calculate solvation energy.");
        }

        // parse output
        let apbs_info = fs::read_to_string(temp_dir.join(format!("{}.out", f_name))).unwrap();
        let mut Esol: Array2<f64> = Array2::zeros((3, atom_num_com));
        let mut Evac: Array2<f64> = Array2::zeros((3, atom_num_com));
        let mut Esas: Array2<f64> = Array2::zeros((3, atom_num_com));
        let apbs_info = apbs_info
            .split("\n")
            .filter_map(|p|
                if p.trim().starts_with("CALCULATION ") ||
                    p.trim().starts_with("Atom") ||
                    p.trim().starts_with("SASA") {
                    Some(p.trim())
                } else { None }
            )
            .collect::<Vec<&str>>()
            .join("\n");

        // extract apbs results
        let apbs_info: Vec<&str> = apbs_info        // list of apbs calculation results
            .split("CALCULATION ")
            .filter_map(|p| match p.trim().len() {
                0 => None,
                _ => Some(p.trim())
            })
            .collect();
        for info in apbs_info {
            let info: Vec<&str> = info
                .split("\n")
                .collect();

            let sys_idx;    // com: 0, rec: 1, lig: 2
            let n: f64;
            if info[0].contains(format!("{}_com", f_name).as_str()) {
                sys_idx = 0;
                n = atom_num_com as f64;
            } else if info[0].contains(format!("{}_rec", f_name).as_str()) {
                sys_idx = 1;
                n = atom_num_rec as f64;
            } else {
                sys_idx = 2;
                n = atom_num_lig as f64;
            }

            if info[0].contains("_VAC") {
                for (idx, v) in info[1..].into_iter().enumerate() {
                    let v: Vec<&str> = v
                        .split(" ")
                        .filter_map(|p| match p.trim().len() {
                            0 => None,
                            _ => Some(p)
                        }).collect();
                    let v: f64 = v[v.len() - 2].parse().unwrap();
                    Evac[[sys_idx, idx]] = v;
                }
            } else if info[0].contains("_SAS") {
                for (idx, v) in info[1..].into_iter().enumerate() {
                    let v: Vec<&str> = v
                        .split(" ")
                        .filter_map(|p| match p.trim().len() {
                            0 => None,
                            _ => Some(p)
                        }).collect();
                    let v: f64 = v[v.len() - 1].parse().unwrap();
                    Esas[[sys_idx, idx]] = gamma * v + _const / n;
                }
            } else {
                for (idx, v) in info[1..].into_iter().enumerate() {
                    let v: Vec<&str> = v
                        .split(" ")
                        .filter_map(|p| match p.trim().len() {
                            0 => None,
                            _ => Some(p)
                        }).collect();
                    let v: f64 = v[v.len() - 2].parse().unwrap();
                    Esol[[sys_idx, idx]] = v;
                }
            }
        }

        Esol = Esol - Evac;
        let pb_com: f64 = Esol.slice(s![0, ..]).iter().sum::<f64>();
        let sa_com: f64 = Esas.slice(s![0, ..]).iter().sum::<f64>();
        let pb_rec: f64 = Esol.slice(s![1, ..]).iter().sum::<f64>();
        let sa_rec: f64 = Esas.slice(s![1, ..]).iter().sum::<f64>();
        let pb_lig: f64 = Esol.slice(s![2, ..]).iter().sum::<f64>();
        let sa_lig: f64 = Esas.slice(s![2, ..]).iter().sum::<f64>();

        vdw[idx] = e_vdw;
        cou[idx] = e_cou;
        mm[idx] = e_cou + e_vdw;
        pb[idx] = pb_com - pb_rec - pb_lig;
        sa[idx] = sa_com - sa_rec - sa_lig;
        dh[idx] = mm[idx] + pb[idx] + sa[idx];

        // residue decomposition
        for &i in &ndx_rec {
            dPBres[aps.atm_resnum[i]] += Esol[[0, i]] - Esol[[1, i]];
            dSAres[aps.atm_resnum[i]] += Esas[[0, i]] - Esas[[1, i]];
        }
        for &i in &ndx_lig {
            dPBres[aps.atm_resnum[i]] += Esol[[0, i]] - Esol[[2, i]];
            dSAres[aps.atm_resnum[i]] += Esas[[0, i]] - Esas[[2, i]];
        }

        COUres += &de_cou;
        VDWres += &de_vdw;
        pb_res += &dPBres;
        sa_res += &dSAres;

        idx += 1;
        pgb.inc(1);
    }
    pgb.finish();

    // residue time average
    COUres /= total_frames as f64;
    VDWres /= total_frames as f64;
    pb_res /= total_frames as f64;
    sa_res /= total_frames as f64;
    MMres = COUres + VDWres;
    dHres = MMres + pb_res + sa_res;

    // totally time average and ts
    let dH = dh.iter().sum::<f64>() / dh.len() as f64;
    let MM = mm.iter().sum::<f64>() / mm.len() as f64;
    let COU = cou.iter().sum::<f64>() / cou.len() as f64;
    let VDW = vdw.iter().sum::<f64>() / vdw.len() as f64;
    let PB = pb.iter().sum::<f64>() / pb.len() as f64;
    let SA = sa.iter().sum::<f64>() / sa.len() as f64;

    let TdS = mm.iter()
        .map(|&p| f64::exp((p - MM) / RT2kJ))
        .sum::<f64>() / mm.len() as f64;
    let TdS = -RT2kJ * TdS.ln();
    let dG = dH - TdS;
    let Ki = f64::exp(dG / RT2kJ);

    println!("MM-PBSA calculation finished.");
    return (dH, MM, PB, SA, COU, VDW, TdS, dG, Ki);
}

fn get_atoms_trj(frames: &Vec<Rc<Frame>>) -> (Array3<f64>, Array3<f64>) {
    let num_frames = frames.len();
    let num_atoms = frames[0].num_atoms();
    let mut coord_matrix: Array3<f64> = Array3::zeros((num_frames, num_atoms, 3));
    let mut box_size: Array3<f64> = Array3::zeros((num_frames, 3, 3));
    let pb = ProgressBar::new(frames.len() as u64);
    set_style(&pb);
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
        pb.inc(1);
    }
    pb.finish();
    return (coord_matrix, box_size);
}

pub fn set_style(pb: &ProgressBar) {
    pb.set_style(ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/ctan} {pos}/{len} {msg}").unwrap()
        .progress_chars("=>-"));
}

fn get_frames_range(frames: &Vec<Rc<Frame>>, bt: f64, et: f64, dt: f64) -> (usize, usize, usize, usize) {
    // decide frame step according to time step
    let time_step = (frames[1].time - frames[0].time) as f64;
    let bf = ((bt - frames[0].time as f64) / time_step) as usize;
    let ef = ((et - frames[0].time as f64) / time_step) as usize;
    let dframe = match dt.partial_cmp(&time_step).unwrap() {
        Ordering::Less => 1,            // converted trajectory with big time step (e.g., 1 ns)
        _ => (dt / time_step) as usize  // time step partially less than target dt (e.g., initial trajectory)
    };
    let total_frames = (ef - bf) / dframe + 1;
    (bf, ef, dframe, total_frames)
}
