use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use xdrfile::*;
use crate::settings::Settings;
use ndarray::parallel::prelude::*;
use ndarray::{ArrayBase, OwnedRepr, ViewRepr, Dim, Array1, Array2, Array3, s};
use std::process::Command;
use std::rc::Rc;
use std::env;
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{Local, Duration};
use crate::coefficients::Coefficients;
use crate::analyzation::Results;
use crate::parse_tpr::{Residue, TPR};
use crate::apbs_param::{PBASet, PBESet};
use crate::atom_property::AtomProperty;
use crate::prepare_apbs::{prepare_pqr, write_apbs_input};

pub fn fun_mmpbsa_calculations(trj: &String, temp_dir: &PathBuf,
                               sys_name: &String, aps: &AtomProperty,
                               ndx_com: &Vec<usize>, ndx_rec: &Vec<usize>, ndx_lig: &Vec<usize>, 
                               residues: &Vec<Residue>, bt: f64, et: f64, dt: f64,
                               pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings)
                               -> Results {
    // run MM/PB-SA calculations
    println!("Running MM/PB-SA calculations of {}...", sys_name);
    println!("Preparing parameters...");
    
    println!("Reading trajectory file...");
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let frames: Vec<Rc<Frame>> = trj.into_iter().map(|p| p.unwrap()).collect();

    println!("Extracting atoms coordination...");
    let (coordinates, _) = get_atoms_trj(&frames);   // frames x atoms(3x1)

    let (bf, ef, dframe, total_frames) = get_frames_range(&frames, bt, et, dt);

    if let Some(_) = settings.apbs.as_ref() {
        println!("Preparing pqr files...");
        prepare_pqr(&frames, bf, ef, dframe, total_frames, &temp_dir,
                    sys_name, &coordinates, ndx_com, &ndx_rec, ndx_lig, aps);
    }

    // calculate MM and PBSA
    calculate_mmpbsa(&frames, &coordinates, bf, ef, dframe, 
        total_frames, aps, &temp_dir, &ndx_com, &ndx_rec, &ndx_lig, residues,
        sys_name, pbe_set, pba_set, settings)
}

fn get_atoms_trj(frames: &Vec<Rc<Frame>>) -> (Array3<f64>, Array3<f64>) {
    let num_frames = frames.len();
    let num_atoms = frames[0].num_atoms();
    let mut coord_matrix: Array3<f64> = Array3::zeros((num_frames, num_atoms, 3));
    let mut box_size: Array3<f64> = Array3::zeros((num_frames, 3, 3));
    let pb = ProgressBar::new(frames.len() as u64);
    set_style(&pb);
    for (idx, frame) in frames.into_iter().enumerate() {
        for (i, a) in (&frame.coords).into_iter().enumerate() {
            for j in 0..3 {
                coord_matrix[[idx, i, j]] = a[j] as f64 * 10.0;
            }
        }
        for (i, b) in (&frame.box_vector).into_iter().enumerate() {
            for j in 0..3 {
                box_size[[idx, i, j]] = b[j] as f64 * 10.0;
            }
        }
        pb.inc(1);
    }
    pb.finish();
    return (coord_matrix, box_size);
}

pub fn set_style(pb: &ProgressBar) {
    pb.set_style(ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:50.cyan/cyan} {pos}/{len} {msg}").unwrap()
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

fn calculate_mmpbsa(frames: &Vec<Rc<Frame>>, coordinates: &Array3<f64>, bf: usize, ef: usize, 
                    dframe: usize, total_frames: usize, aps: &AtomProperty, temp_dir: &PathBuf,
                    ndx_com_norm: &Vec<usize>, ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>,
                    residues: &Vec<Residue>, sys_name: &String, 
                    pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> Results {
    let mut elec_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    let mut vdw_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    let mut pb_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    let mut sa_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    
    // parameters for elec calculation
    let coeff = Coefficients::new(pbe_set);

    // Time list of trajectory
    let times: Array1<f64> = (bf..=ef).step_by(dframe)
        .map(|p| frames[p].time as f64).collect();

    // start calculation
    env::set_var("OMP_NUM_THREADS", settings.nkernels.to_string());
    let t_start = Local::now();
    
    println!("Calculating MM/PB-SA binding energy...");

    let pgb = ProgressBar::new(total_frames as u64);
    set_style(&pgb);
    pgb.inc(0);
    let mut idx = 0;
    for cur_frm in (bf..=ef).step_by(dframe) {
        // MM
        let coord = coordinates.slice(s![cur_frm, .., ..]);
        if ndx_lig_norm[0] != ndx_rec_norm[0] {
            let (res_elec, res_vdw) = 
                calc_mm(&ndx_rec_norm, &ndx_lig_norm, aps, &coord, residues, &coeff, &settings);
            elec_res.row_mut(idx).assign(&res_elec);
            vdw_res.row_mut(idx).assign(&res_vdw);
        }

        // PBSA
        calc_pbsa(idx, &coord, frames, ndx_rec_norm, ndx_lig_norm, ndx_com_norm,
            &mut pb_res, &mut sa_res, cur_frm, sys_name, temp_dir, aps, pbe_set, pba_set, settings);

        pgb.inc(1);
        pgb.set_message(format!("at {} ns, ΔH={:.2} kJ/mol, eta. {} s", 
                                        times[idx] / 1000.0,
                                        vdw_res.row(idx).sum() + elec_res.row(idx).sum() + pb_res.row(idx).sum() + sa_res.row(idx).sum(),
                                        pgb.eta().as_secs()));

        idx += 1;
    }
    pgb.finish();

    // end calculation
    let t_end = Local::now();
    let t_spend = Duration::from(t_end - t_start).num_milliseconds();
    println!("MM/PB-SA calculation of {} finished. Total time cost: {} s", sys_name, t_spend as f64 / 1000.0);
    env::remove_var("OMP_NUM_THREADS");
    
    // whether remove temp directory
    if !settings.debug_mode {
        if settings.apbs.is_some() {
            fs::remove_dir_all(&temp_dir).expect("Remove dir failed");
        }
    }

    Results::new(
        aps,
        residues,
        ndx_rec_norm,
        ndx_lig_norm,
        &times,
        coordinates.clone(),
        &elec_res,
        &vdw_res,
        &pb_res,
        &sa_res,
    )
}

pub fn get_residues(tpr: &TPR, ndx_com: &Vec<usize>) -> Vec<Residue> {
    let mut residues: Vec<Residue> = vec![];
    let mut idx = 0;
    let mut resind_offset = 0;
    
    let pb = ProgressBar::new(tpr.n_atoms as u64);
    pb.set_style(ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:50.cyan/cyan} {percent}% {msg}").unwrap()
        .progress_chars("=>-"));
    for mol in &tpr.molecules {
        for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
            for atom in &mol.atoms {
                idx += 1;
                if ndx_com.contains(&idx) && residues.len() <= atom.resind + resind_offset {
                    residues.push(mol.residues[atom.resind].to_owned());
                }
                pb.inc(1);
                pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
            }
            resind_offset += mol.residues.len();
        }
    }
    pb.finish();
    residues
}

fn calc_mm(ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>, aps: &AtomProperty, coord: &ArrayBase<ViewRepr<&f64>, Dim<[usize; 2]>>, 
            residues: &Vec<Residue>, coeff: &Coefficients, settings: &Settings) -> (Array1<f64>, Array1<f64>) {
    let kj_elec = coeff.kj_elec;
    let kap = coeff.kap;
    let pdie = coeff.pdie;
    let mut de_elec: Array1<f64> = Array1::zeros(residues.len());
    let mut de_vdw: Array1<f64> = Array1::zeros(residues.len());

    for &i in ndx_rec_norm {
        let qi = aps.atm_charge[i];
        let ci = aps.atm_typeindex[i];
        let xi = coord[[i, 0]];
        let yi = coord[[i, 1]];
        let zi = coord[[i, 2]];
        for &j in ndx_lig_norm {
            if ndx_lig_norm[0] == ndx_rec_norm[0] && j <= i {
                continue;
            }
            let qj = aps.atm_charge[j];
            let cj = aps.atm_typeindex[j];
            let xj = coord[[j, 0]];
            let yj = coord[[j, 1]];
            let zj = coord[[j, 2]];
            let r = f64::sqrt((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2));
            if r < settings.r_cutoff {
                let e_elec = match settings.use_dh {
                    false => qi * qj / r,
                    _ => qi * qj / r * f64::exp(-kap * r)   // doi: 10.1088/0256-307X/38/1/018701
                };
                let r = r / 10.0;
                let e_vdw = (aps.c12[[ci, cj]] / r.powi(6) - aps.c6[[ci, cj]]) / r.powi(6);
                de_elec[aps.atm_resid[i]] += e_elec;
                de_elec[aps.atm_resid[j]] += e_elec;
                de_vdw[aps.atm_resid[i]] += e_vdw;
                de_vdw[aps.atm_resid[j]] += e_vdw;
            }
        }
    }

    de_elec.par_iter_mut().for_each(|p| *p *= kj_elec / (2.0 * pdie));
    de_vdw.par_iter_mut().for_each(|p| *p /= 2.0);

    return (de_elec, de_vdw)
}

fn calc_pbsa(idx: usize, coord: &ArrayBase<ViewRepr<&f64>, Dim<[usize; 2]>>, frames: &Vec<Rc<Frame>>, 
            ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>, ndx_com_norm: &Vec<usize>,
            pb_res: &mut ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, sa_res: &mut ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>,
            cur_frm: usize, sys_name: &String, temp_dir: &PathBuf, 
            aps: &AtomProperty, pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) {
    // From AMBER-PB4, the surface extension constant γ=0.0072 kcal/(mol·Å2)=0.030125 kJ/(mol·Å^2)
    // but the default gamma parameter for apbs calculation is set to 1, in order to directly obtain the surface area
    // then the SA energy term is calculated by s_mmpbsa
    let gamma = 0.030125;
    let bias = 0.0;
    let f_name = format!("{}_{}ns", sys_name, frames[cur_frm].time / 1000.0);
    if let Some(apbs) = &settings.apbs {
        write_apbs_input(ndx_rec_norm, ndx_lig_norm, coord, &aps.atm_radius,
                pbe_set, pba_set, temp_dir, &f_name, settings);
        // invoke apbs program to do apbs calculations
        let apbs_result = Command::new(apbs).arg(format!("{}.apbs", f_name)).current_dir(temp_dir).output().expect("running apbs failed.");
        let apbs_err = String::from_utf8(apbs_result.stderr).expect("Failed to parse apbs output.");
        let apbs_result = String::from_utf8(apbs_result.stdout).expect("Failed to parse apbs output.");
        if settings.debug_mode {
            let mut outfile = File::create(temp_dir.join(format!("{}.out", f_name))).expect("Failed to create output file.");
            outfile.write_all(apbs_result.as_bytes()).expect("Failed to write apbs output.");
            let mut errfile = File::create(temp_dir.join(format!("{}.err", f_name))).expect("Failed to create err file.");
            errfile.write_all(apbs_err.as_bytes()).expect("Failed to write apbs output.");
        }
        // let apbs_result = fs::read_to_string(temp_dir.join(format!("{}.out", f_name))).expect("Failed to parse apbs output.");

        // preserve CALCULATION, Atom and SASA lines
        let apbs_result: Vec<&str> = apbs_result.split("\n").filter_map(|p|
            if p.trim().starts_with("CALCULATION") || p.trim().starts_with("Atom") || p.trim().starts_with("SASA") {
                Some(p.trim())
            } else {
                None
            }
        ).collect();

        // extract apbs results
        let indexes: Vec<usize> = apbs_result.iter().enumerate().filter_map(|(i, &p)| match p.starts_with("CAL") {
            true => Some(i),
            false => None
        }).collect();

        let mut com_pb_sol: Vec<f64> = vec![];
        let mut com_pb_vac: Vec<f64> = vec![];
        let mut rec_pb_sol: Vec<f64> = vec![];
        let mut rec_pb_vac: Vec<f64> = vec![];
        let mut lig_pb_sol: Vec<f64> = vec![];
        let mut lig_pb_vac: Vec<f64> = vec![];
        let mut com_sa: Vec<f64> = vec![];
        let mut rec_sa: Vec<f64> = vec![];
        let mut lig_sa: Vec<f64> = vec![];

        let mut skip_pb = true;     // the first time PB calculation should be wasted
        for (i, &idx) in indexes.iter().enumerate() {
            let st = idx + 1;
            let ed = match indexes.get(i + 1) {
                Some(&idx) => idx,
                None => apbs_result.len()
            };
            if apbs_result[idx].contains(&"_com_SOL") {
                if !skip_pb {
                    apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut com_pb_sol);
                }
                skip_pb = !skip_pb;
            } else if apbs_result[idx].contains(&"_com_VAC") {
                if !skip_pb {
                    apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut com_pb_vac);
                }
                skip_pb = !skip_pb;
            } else if apbs_result[idx].contains(&"_rec_SOL") {
                if !skip_pb {
                    apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut rec_pb_sol);
                }
                skip_pb = !skip_pb;
            } else if apbs_result[idx].contains(&"_rec_VAC") {
                if !skip_pb {
                    apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut rec_pb_vac);
                }
                skip_pb = !skip_pb;
            } else if apbs_result[idx].contains(&"_lig_SOL") {
                if !skip_pb {
                    apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut lig_pb_sol);
                }
                skip_pb = !skip_pb;
            } else if apbs_result[idx].contains(&"_lig_VAC") {
                if !skip_pb {
                    apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut lig_pb_vac);
                }
                skip_pb = !skip_pb;
            } else if apbs_result[idx].contains(&"_com_SAS") {
                apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut com_sa);
            } else if apbs_result[idx].contains(&"_rec_SAS") {
                apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut rec_sa);
            } else if apbs_result[idx].contains(&"_lig_SAS") {
                apbs_result[st..ed].par_iter().map(|&p| parse_apbs_line(p)).collect_into_vec(&mut lig_sa);
            }
        }

        let com_pb: Array1<f64> = Array1::from_vec(com_pb_sol) - Array1::from_vec(com_pb_vac);
        let com_sa: Array1<f64> = Array1::from_vec(com_sa.par_iter().map(|i| gamma * *i + bias / com_sa.len() as f64).collect());
        let rec_pb: Array1<f64> = Array1::from_vec(rec_pb_sol) - Array1::from_vec(rec_pb_vac);
        let rec_sa: Array1<f64> = Array1::from_vec(rec_sa.par_iter().map(|i| gamma * *i + bias / rec_sa.len() as f64).collect());
        let lig_pb: Array1<f64> = Array1::from_vec(lig_pb_sol) - Array1::from_vec(lig_pb_vac);
        let lig_sa: Array1<f64> = Array1::from_vec(lig_sa.par_iter().map(|i| gamma * *i + bias / lig_sa.len() as f64).collect());

        // residue decomposition
        let offset_rec = match ndx_lig_norm[0].cmp(&ndx_rec_norm[0]) {
            Ordering::Less => 0,
            Ordering::Greater => ndx_rec_norm.len(),
            Ordering::Equal => 0
        };
        let offset_lig = match ndx_rec_norm[0].cmp(&ndx_lig_norm[0]) {
            Ordering::Less => 0,
            Ordering::Greater => ndx_lig_norm.len(),
            Ordering::Equal => 0
        };

        if ndx_rec_norm[0] == ndx_lig_norm[0] {
            // if no ligand, pb_com = pb_lig = 0, so real energy is inversed rec_pbsa
            for &i in ndx_com_norm {
                pb_res[[idx, aps.atm_resid[i]]] += rec_pb[i - offset_lig];
                sa_res[[idx, aps.atm_resid[i]]] += rec_sa[i - offset_lig];
            }
        } else {
            for &i in ndx_com_norm {
                if ndx_rec_norm.contains(&i) {
                    pb_res[[idx, aps.atm_resid[i]]] += com_pb[i] - rec_pb[i - offset_lig];
                    sa_res[[idx, aps.atm_resid[i]]] += com_sa[i] - rec_sa[i - offset_lig];
                } else {
                    pb_res[[idx, aps.atm_resid[i]]] += com_pb[i] - lig_pb[i - offset_rec];
                    sa_res[[idx, aps.atm_resid[i]]] += com_sa[i] - lig_sa[i - offset_rec];
                }
            }
        }
    }
}

fn parse_apbs_line(line: &str) -> f64 {
    line.split(":")
        .skip(1)
        .next().expect("Cannot get information from apbs")
        .trim_start()
        .split(" ")
        .next().expect("Cannot get information from apbs")
        .parse().expect("Cannot parse value from apbs")
}