use std::cmp::Ordering;
use std::fs;
use std::path::PathBuf;
use xdrfile::*;
use crate::settings::Settings;
use ndarray::{ArrayBase, OwnedRepr, ViewRepr, Dim, Array1, Array2, Array3, s};
use std::fs::File;
use std::io::Write;
use std::process::Command;
use std::rc::Rc;
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{Local, Duration};
use crate::analyzation::Results;
use crate::parse_tpr::TPR;
use crate::apbs_param::{PBASet, PBESet};
use crate::atom_property::AtomProperty;
use crate::prepare_apbs::{prepare_pqr, write_apbs_input};

pub fn fun_mmpbsa_calculations(trj: &String, tpr: &TPR, temp_dir: &PathBuf,
                               sys_name: &String, aps: &AtomProperty,
                               ndx_com: &Vec<usize>, ndx_rec: &Vec<usize>, ndx_lig: Option<&Vec<usize>>,
                               bt: f64, et: f64, dt: f64,
                               pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings)
                               -> Results {
    // run MM/PB-SA calculations
    println!("Running MM/PB-SA calculations...");
    println!("Preparing parameters...");

    // pre-treat trajectory: fix pbc
    println!("Reading trajectory file...");
    let trj = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let frames: Vec<Rc<Frame>> = trj.into_iter().map(|p| p.unwrap()).collect();
    // pbc whole 先不写, 先默认按照已经消除了周期性来做后续处理, 之后再看周期性的事
    println!("Extracting atoms coordination...");
    let (coordinates, _) = get_atoms_trj(&frames);   // frames x atoms(3x1)

    let (bf, ef, dframe, total_frames) = get_frames_range(&frames, bt, et, dt);

    if let Some(_) = settings.apbs.as_ref() {
        println!("Preparing pqr files...");
        prepare_pqr(&frames, bf, ef, dframe, total_frames, &temp_dir,
                    sys_name, &coordinates, ndx_com, ndx_rec, ndx_lig, aps);
    }

    // must appear here because the tpr and xtc need origin indexes
    let (ndx_com_norm, ndx_rec_norm, ndx_lig_norm) = 
        normalize_index(ndx_com, ndx_rec, ndx_lig);

    // calculate MM and PBSA
    println!("Start MM/PB-SA calculations...");
    let results = calculate_mmpbsa(tpr, &frames, &coordinates, bf, ef, dframe, 
        total_frames, aps, &temp_dir, &ndx_com_norm, &ndx_rec_norm, &ndx_lig_norm, &ndx_com,
        sys_name, pbe_set, pba_set, settings);
    println!("MM/PB-SA calculation finished.");
    results
}

// convert rec and lig to begin at 0 and continous
fn normalize_index(ndx_com: &Vec<usize>, ndx_rec: &Vec<usize>, ndx_lig: Option<&Vec<usize>>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut ndx_lig = match ndx_lig {
        Some(ndx_lig) => ndx_lig.iter().map(|p| p - ndx_com[0]).collect(),
        None => Vec::from_iter(0..ndx_rec.len())
    };
    let mut ndx_rec: Vec<usize> = ndx_rec.iter().map(|p| p - ndx_com[0]).collect();
    let ndx_com = match ndx_lig[0].cmp(&ndx_rec[0]) {
        Ordering::Greater => {
            ndx_lig = ndx_lig.iter().map(|p| p - ndx_lig[0] + ndx_rec.len()).collect();
            let mut ndx_com = ndx_rec.to_vec();
            ndx_com.extend(&ndx_lig);
            ndx_com
        }
        Ordering::Less => {
            ndx_rec = ndx_rec.iter().map(|p| p - ndx_rec[0] + ndx_lig.len()).collect();
            let mut ndx_com = ndx_lig.to_vec();
            ndx_com.extend(&ndx_rec);
            ndx_com
        }
        Ordering::Equal => Vec::from_iter(0..ndx_rec.len())
    };
    (ndx_com, ndx_rec, ndx_lig)
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

fn calculate_mmpbsa(tpr: &TPR, frames: &Vec<Rc<Frame>>, coordinates: &Array3<f64>,
                    bf: usize, ef: usize, dframe: usize,
                    total_frames: usize, aps: &AtomProperty, temp_dir: &PathBuf,
                    ndx_com_norm: &Vec<usize>, ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>,
                    ndx_com: &Vec<usize>,
                    sys_name: &String, pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> Results {
    let residues: Array1<(i32, String)> = get_residues(tpr, ndx_com);
    let mut elec_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    let mut vdw_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    let mut pb_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));
    let mut sa_res: Array2<f64> = Array2::zeros((total_frames, residues.len()));

    let pgb = ProgressBar::new(total_frames as u64);
    set_style(&pgb);
    pgb.inc(0);
    let mut idx = 0;
    for cur_frm in (bf..=ef).step_by(dframe) {
        // MM
        let coord = coordinates.slice(s![cur_frm, .., ..]);
        if ndx_lig_norm[0] != ndx_rec_norm[0] {
            calc_mm(idx, &ndx_rec_norm, &ndx_lig_norm, &aps, &coord, 
                &residues, &mut elec_res, &mut vdw_res,
                &pbe_set, &settings);
        }

        // PBSA
        calc_pbsa(idx, &coord, frames, ndx_rec_norm, ndx_lig_norm, ndx_com_norm,
            &mut pb_res, &mut sa_res, cur_frm, sys_name, temp_dir, aps, pbe_set, pba_set, settings);

        pgb.inc(1);
        pgb.set_message(format!("Will finish at {}", 
            Local::now().checked_add_signed(Duration::seconds(pgb.eta().as_secs() as i64))
            .expect("Failed to get finish time").format("%Y-%m-%d %H:%M:%S").to_string()));
        idx += 1;
    }
    pgb.finish();

    // Time list of trajectory
    let times: Array1<f64> = (bf..=ef).step_by(dframe)
        .map(|p| frames[p].time as f64).collect();

    Results::new(
        tpr,
        times,
        ndx_com,
        elec_res,
        vdw_res,
        pb_res,
        sa_res
    )
}

pub fn get_residues(tpr: &TPR, ndx_com: &Vec<usize>) -> Array1<(i32, String)> {
    let mut residues: Vec<(i32, String)> = vec![];
    let mut idx = 0;
    let mut resind_offset = 0;
    for mol in &tpr.molecules {
        for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
            for atom in &mol.atoms {
                idx += 1;
                if ndx_com.contains(&idx) && residues.len() <= atom.residue_index + resind_offset {
                    let res = &mol.residues[atom.residue_index];
                    residues.push((res.nr, res.name.to_string()));
                }
            }
            resind_offset += mol.residues.len();
        }
    }
    Array1::from_vec(residues)
}

fn calc_mm(idx: usize, ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>, aps: &AtomProperty, coord: &ArrayBase<ViewRepr<&f64>, Dim<[usize; 2]>>, 
            residues: &Array1<(i32, String)>, elec_res: &mut Array2<f64>, vdw_res: &mut Array2<f64>,
            pbe_set: &PBESet, settings: &Settings) {
    let eps0 = 8.854187812800001e-12;
    let kb = 1.380649e-23;
    let na = 6.02214076e+23;
    let qe = 1.602176634e-19;
    let ion_strength: f64 = pbe_set.ions.iter()
        .map(|ion| ion.charge * ion.charge * ion.conc).sum();
    let kap = 1e-10 / f64::sqrt(eps0 * kb * pbe_set.temp * pbe_set.sdie / (ion_strength * qe * qe * na * 1e3));

    let kj_elec = 1389.35457520287;

    let mut de_elec: Array1<f64> = Array1::zeros(residues.len());
    let mut de_vdw: Array1<f64> = Array1::zeros(residues.len());

    for &i in ndx_rec_norm {
        let qi = aps.atm_charge[i];
        let ci = aps.atm_typeindex[i];
        let xi = coord[[i, 0]];
        let yi = coord[[i, 1]];
        let zi = coord[[i, 2]];
        for &j in ndx_lig_norm {
            // if ndx_lig_norm[0] == ndx_rec_norm[0] && j <= i {
            //     continue;
            // }
            let qj = aps.atm_charge[j];
            let cj = aps.atm_typeindex[j];
            let xj = coord[[j, 0]];
            let yj = coord[[j, 1]];
            let zj = coord[[j, 2]];
            let r = f64::sqrt((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2));
            if r < settings.r_cutoff {
                let e_elec = match settings.use_dh {
                    false => qi * qj / r,
                    _ => qi * qj / r * f64::exp(-kap * r)
                };
                let e_vdw = aps.c12[[ci, cj]] / (r / 10.0).powi(12) - aps.c6[[ci, cj]] / (r / 10.0).powi(6);
                de_elec[aps.atm_resnum[i]] += e_elec;
                de_elec[aps.atm_resnum[j]] += e_elec;
                de_vdw[aps.atm_resnum[i]] += e_vdw;
                de_vdw[aps.atm_resnum[j]] += e_vdw;
            }
        }
    }
    for i in 0..residues.len() {
        de_elec[i] *= kj_elec / (2.0 * pbe_set.pdie);
        de_vdw[i] /= 2.0;
        elec_res[[idx, i]] += de_elec[i];
        vdw_res[[idx, i]] += de_vdw[i];
    }
}

fn calc_pbsa(idx: usize, coord: &ArrayBase<ViewRepr<&f64>, Dim<[usize; 2]>>, frames: &Vec<Rc<Frame>>, 
            ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>, ndx_com_norm: &Vec<usize>,
            pb_res: &mut ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, sa_res: &mut ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>,
            cur_frm: usize, sys_name: &String, temp_dir: &PathBuf, 
            aps: &AtomProperty, pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) {

    // From AMBER-PB4, the surface extension constant γ=0.0072 kcal/(mol·Å2)=0.030125 kJ/(mol·Å2)
    // but the default gamma parameter for apbs calculation is set to 1, in order to directly obtain the surface area
    // then the SA energy term is calculated by super_mmpbsa
    let gamma = 0.030125;
    let bias = 0.0;
    let f_name = format!("{}_{}ns", sys_name, frames[cur_frm].time / 1000.0);
    if let Some(apbs) = &settings.apbs {
        write_apbs_input(ndx_rec_norm, ndx_lig_norm, coord, &aps.atm_radius,
                pbe_set, pba_set, temp_dir, &f_name, settings);
        // invoke apbs program to do apbs calculations
        let mut apbs_out = File::create(temp_dir.join(format!("{}.out", f_name))).
            expect("Failed to create apbs out file");
        let apbs_result = Command::new(apbs).
            arg(format!("{}.apbs", f_name)).
            current_dir(temp_dir).output().expect("running apbs failed.");
        let apbs_output = String::from_utf8(apbs_result.stdout).
            expect("Failed to get apbs output.");
        apbs_out.write_all(apbs_output.as_bytes()).expect("Failed to write apbs output");

        // parse output
        let apbs_results = fs::read_to_string(temp_dir.join(format!("{}.out", f_name))).unwrap();

        let mut e_com: Array2<f64> = Array2::zeros((3, ndx_rec_norm.len() + ndx_lig_norm.len()));
        let mut e_rec: Array2<f64> = Array2::zeros((3, ndx_rec_norm.len()));
        let mut e_lig: Array2<f64> = Array2::zeros((3, ndx_lig_norm.len()));

        // preserve CALCULATION, Atom and SASA lines
        let apbs_results = apbs_results
            .split("\n")
            .filter_map(|p|
                if p.trim().starts_with("CALCULATION ") ||
                    p.trim().starts_with("Atom") ||
                    p.trim().starts_with("SASA") {
                    Some(p.trim())
                } else { None }
            )
            .collect::<Vec<&str>>();

        // extract apbs results
        let mut cur_sys = 0;    // com=0, rec=1, lig=2
        let mut cur_item = 0;   // total=0, vac=1, sas=2
        let mut cur_atom = 0;   // current atom index
        for line in apbs_results {
            if line.starts_with("CALCULATION") {
                if line.contains(format!("{}_com", f_name).as_str()) {
                    cur_sys = 0;
                } else if line.contains(format!("{}_rec", f_name).as_str()) {
                    cur_sys = 1;
                } else if line.contains(format!("{}_lig", f_name).as_str()) {
                    cur_sys = 2;
                }
                if line.contains("_VAC") {
                    cur_item = 1;
                } else if line.contains("_SAS") {
                    cur_item = 2;
                } else {
                    cur_item = 0;
                }
                cur_atom = 0;
            } else {
                let v: Vec<&str> = line
                    .split(":")
                    .collect();
                let v: Vec<&str> = v[1]
                    .split(" ")
                    .filter_map(|p| match p.trim().len() {
                        0 => None,
                        _ => Some(p)
                    }).collect();
                match cur_sys {
                    0 => e_com[[cur_item, cur_atom]] = v[0].parse().unwrap(),
                    1 => e_rec[[cur_item, cur_atom]] = v[0].parse().unwrap(),
                    2 => e_lig[[cur_item, cur_atom]] = v[0].parse().unwrap(),
                    _ => ()
                }
                cur_atom += 1;
            }
        }

        let f = |e_arr: &mut Array2<f64>, n_cols: usize|
            for mut col in e_arr.columns_mut() {
                col[0] -= col[1];
                col[2] = gamma * col[2] + bias / n_cols as f64;
            };
        f(&mut e_com, ndx_com_norm.len());
        f(&mut e_rec, ndx_rec_norm.len());
        f(&mut e_lig, ndx_lig_norm.len());

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
        for &i in ndx_com_norm {
            if ndx_rec_norm.contains(&i) {
                pb_res[[idx, aps.atm_resnum[i]]] += e_com[[0, i]] - e_rec[[0, i - offset_lig]];
                sa_res[[idx, aps.atm_resnum[i]]] += e_com[[2, i]] - e_rec[[2, i - offset_lig]];
            } else {
                pb_res[[idx, aps.atm_resnum[i]]] += e_com[[0, i]] - e_lig[[0, i - offset_rec]];
                sa_res[[idx, aps.atm_resnum[i]]] += e_com[[2, i]] - e_lig[[2, i - offset_rec]];
            }
        }

        // if no ligand, pb_com = pb_lig = 0, so pb_rec should be inversed to be real energy
        if ndx_lig_norm[0] == ndx_rec_norm[0] {
            for j in 0..pb_res.shape()[1] {
                pb_res[[idx, j]] = -pb_res[[idx, j]];
                sa_res[[idx, j]] = -sa_res[[idx, j]];
            }
        }
    }
}