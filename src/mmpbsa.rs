use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::env;
use ndarray::parallel::prelude::*;
use ndarray::{s, Array1, Array2, Array3, ArrayView2, Axis};
use indicatif::ProgressBar;
use chrono::{Local, Duration};
use crate::fun_para_system::normalize_index;
use crate::settings::Settings;
use crate::utils::{self, is_amino, set_style};
use crate::coefficients::{self, Coefficients};
use crate::analysis::{SMResult, SMResults};
use crate::parse_tpr::Residue;
use crate::parameters::{PBASet, PBESet};
use crate::atom_property::{AtomProperties, AtomProperty};
use crate::prepare_apbs::{build_apbs_input_text, build_molecules, prepare_pqr};
use crate::apbs_runner::{self, SolverCalcKind, SolverRun};

pub fn fun_mmpbsa_calculations(time_list: &Vec<f64>, time_list_ie: &Vec<f64>, coordinates_ie: &Array3<f64>, 
                               temp_dir: &PathBuf, sys_name: &str, aps: &AtomProperties,
                               ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>,
                               ala_list: &Vec<i32>, residues: &Vec<Residue>, temperature: f64,
                               pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings)
                               -> SMResults {
    println!("Running MM-PBSA calculations of {}...", sys_name);

    // calculate MM and PBSA
    println!("Calculating binding energy for {}...", sys_name);
    let result_wt = calculate_mmpbsa(time_list, time_list_ie, coordinates_ie, aps, &temp_dir, 
        &ndx_rec, &ndx_lig, residues, temperature,
        sys_name, "", pbe_set, pba_set, settings);

    let mut result_ala_scan: Vec<SMResult> = vec![result_wt];
    if ala_list.len() > 0 {
        // main chain atoms number
        let as_res: Vec<&Residue> = residues.iter().filter(|&r| ala_list.contains(&r.nr) 
            && is_amino(&r.name)).collect();     // gly not contain CB, ala no need to mutate
        let exclude_list = ["N", "CA", "C", "O", "CB", "H", "HA", "HB1", "HB2"];
        for asr in as_res {
            // let (new_coordinates, new_aps, new_ndx_rec, new_ndx_lig) = 
            //     ala_mutate(aps, asr, &exclude_list, coordinates, ndx_rec, ndx_lig);
            let (new_coordinates_ie, new_aps, new_ndx_rec, new_ndx_lig) = 
                ala_mutate(aps, asr, &exclude_list, coordinates_ie, ndx_rec, ndx_lig);

            // After alanine mutation
            let mutation = match utils::resname_3to1(&asr.name) {
                Some(mutation) => mutation,
                None => asr.name.to_string()
            };

            let mut new_residues = residues.clone();
            new_residues[asr.id].name = "ALA".to_string();

            let mutation = format!("-{}{}A", mutation, asr.nr);
            let sys_name = format!("{}{}", sys_name, mutation);
            println!("Calculating binding energy for {}...", sys_name);
            let result_as = calculate_mmpbsa(time_list, time_list_ie, &new_coordinates_ie,
                &new_aps, &temp_dir, &new_ndx_rec, &new_ndx_lig, &new_residues, temperature,
                &sys_name, &mutation, pbe_set, pba_set, settings);
            result_ala_scan.push(result_as);
        }
    };

    println!("");
    utils::show_famous_quotes();

    println!("Writing calculation results...");
    let sm_results = SMResults::new(result_ala_scan);
    sm_results.to_bin(Path::new(&format!("MMPBSA_{}.sm", sys_name)));

    sm_results
}

fn ala_mutate(aps: &AtomProperties, asr: &Residue, exclude_list: &[&str], coordinates: &Array3<f64>, 
                ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>)
                -> (Array3<f64>, AtomProperties, BTreeSet<usize>, Option<BTreeSet<usize>>) {
    let mut new_aps = aps.clone();
    let as_atoms: Vec<AtomProperty> = aps.atom_props.iter().filter_map(|a| if a.resid == asr.id {
        Some(a.clone())
    } else {
        None
    }).collect();
    let mut sc_out: Vec<&AtomProperty> = as_atoms.iter().filter(|&a| !exclude_list.contains(&a.name.as_str())).collect();
    let xgs: Vec<AtomProperty> = as_atoms.iter().filter_map(|a| {
        if a.name.eq("CG") || a.name.eq("CG1") || a.name.eq("CG2") || a.name.eq("SG") || a.name.eq("OG") || a.name.eq("OG1") {
            Some(a.clone())
        } else {
            None
        }}).collect();
    // 通过CB定位新的HB
    let cb: Vec<AtomProperty> = as_atoms.iter().filter_map(|a| {
        if a.name.eq("CB") {
            Some(a.clone())
        } else {
            None
        }}).collect();
    // 脯氨酸需要把CD改成H（CD因此保留，不进入删除列表）
    let pro_cd: Option<&AtomProperty> = if asr.name.eq("PRO") {
        sc_out.retain(|&a| a.name.ne("CD"));
        Some(as_atoms.iter().find(|&a| a.name == "CD").unwrap())
    } else {
        None
    };
    // 突变残基中重命名的原子（CG/SG/OG* 与脯氨酸CD）在删除列表之外保留
    for xg in xgs.iter() {
        new_aps.atom_props[xg.id].change_atom(aps.at_map.get("H"), "HB3", "HC", &aps.radius_type);
    }
    if let Some(cd) = pro_cd {
        new_aps.atom_props[cd.id].change_atom(aps.at_map.get("H"), "H", "H", &aps.radius_type);
    }

    // delete other atoms in the scanned residue
    let del_list: Vec<usize> = sc_out.iter().map(|a| a.id).collect();
    let xg_list: Vec<usize> = xgs.iter().map(|a| a.id).collect();
    new_aps.atom_props.retain(|a| !del_list.contains(&a.id) || xg_list.contains(&a.id));
    let retain_id: Vec<usize> = new_aps.atom_props.iter().map(|a| a.id).collect();
    // 每次删除原子后重新排序剩余原子id
    for (i, ap) in new_aps.atom_props.iter_mut().enumerate() {
        ap.id = i;
    };

    // 只拷贝保留的原子：先 select 出子轨迹，再在子集上重建氢原子坐标，
    // 避免对整条轨迹做一次全量克隆
    let mut new_coordinates = coordinates.select(Axis(1), &retain_id);
    let new_index: HashMap<usize, usize> = retain_id.iter()
        .enumerate()
        .map(|(i, &old_id)| (old_id, i))
        .collect();

    // 获取新的HB坐标
    for xg in xgs.iter() {
        let cb_id = new_index[&cb[0].id];
        let xg_id = new_index[&xg.id];
        for layer in 0..new_coordinates.shape()[0] {
            let cb_coords: Array1<f64> = new_coordinates.slice(s![layer, cb_id, ..]).to_owned();
            let hg_coords: Array1<f64> = new_coordinates.slice(s![layer, xg_id, ..]).to_owned();
            let new_hg_coord: Array1<f64> = transform_coordinate(&cb_coords, &hg_coords, 1.09);
            new_coordinates[[layer, xg_id, 0]] = new_hg_coord[0];
            new_coordinates[[layer, xg_id, 1]] = new_hg_coord[1];
            new_coordinates[[layer, xg_id, 2]] = new_hg_coord[2];
        }
    }
    if let Some(cd) = pro_cd {
        let n = as_atoms.iter().find(|&a| a.name == "N").unwrap();
        let n_id = new_index[&n.id];
        let cd_id = new_index[&cd.id];
        // 获取新的HN坐标
        for layer in 0..new_coordinates.shape()[0] {
            let n_coords: Array1<f64> = new_coordinates.slice(s![layer, n_id, ..]).to_owned();
            let hn_coords: Array1<f64> = new_coordinates.slice(s![layer, cd_id, ..]).to_owned();
            let new_hn_coord = transform_coordinate(&n_coords, &hn_coords, 1.07);
            new_coordinates[[layer, cd_id, 0]] = new_hn_coord[0];
            new_coordinates[[layer, cd_id, 1]] = new_hn_coord[1];
            new_coordinates[[layer, cd_id, 2]] = new_hn_coord[2];
        }
    }

    let mut new_ndx_rec = ndx_rec.clone();
    new_ndx_rec.retain(|&x| !del_list.contains(&x) || xg_list.contains(&x));
    let (new_ndx_rec, new_ndx_lig) = normalize_index(&new_ndx_rec, ndx_lig);

    (new_coordinates, new_aps, new_ndx_rec, new_ndx_lig)
}

fn calculate_mmpbsa(time_list: &Vec<f64>, time_list_ie: &Vec<f64>, coordinates_ie: &Array3<f64>, 
                    aps: &AtomProperties, temp_dir: &PathBuf,
                    ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>,
                    residues: &Vec<Residue>, temperature: f64, sys_name: &str, mutation: &str,
                    pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> SMResult {
    let mut elec_atom: Array2<f64> = Array2::zeros((time_list.len(), aps.atom_props.len()));
    let mut vdw_atom: Array2<f64> = Array2::zeros((time_list.len(), aps.atom_props.len()));
    let mut pb_atom: Array2<f64> = Array2::zeros((time_list.len(), aps.atom_props.len()));
    let mut sa_atom: Array2<f64> = Array2::zeros((time_list.len(), aps.atom_props.len()));
    
    // parameters for elec calculation
    let coeff = Coefficients::new(pbe_set);

    // Time list of trajectory
    let times: Vec<f64> = time_list.iter().map(|t| t / 1000.0).collect();
    let times_ie: Vec<f64> = time_list_ie.iter().map(|t| t / 1000.0).collect();

    // extract coordinates for PBSA from initial
    let coordinates: Array3<f64> = {
        // 预先计算需要保留的索引
        let indices: Vec<usize> = times_ie.iter()
            .enumerate()
            .filter_map(|(i, t)| times.contains(t).then_some(i))
            .collect();
        
        if indices.is_empty() {
            // 如果没有有效帧，返回一个空的 Array3
            Array3::zeros((0, coordinates_ie.shape()[1], coordinates_ie.shape()[2]))
        } else {
            // 直接选择需要的帧，避免克隆
            coordinates_ie.select(Axis(0), &indices)
        }
    };

    // setting up environment
    env::set_var("OMP_NUM_THREADS", settings.nkernels.to_string());

    // Decide whether to use parallel
    let total_iterations = if let Some(ndx_lig) = ndx_lig {
        ndx_rec.len() * ndx_lig.len()
    } else {
        ndx_rec.len()
    };
    const PARALLEL_THRESHOLD: usize = 100000;
    let use_parallel = total_iterations > PARALLEL_THRESHOLD;
    if use_parallel {
        println!("Since there are too many atom pairs (> {}), will use parallel computation for MM.", PARALLEL_THRESHOLD);
    }

    let t_start = Local::now();

    println!("Start MM-PBSA calculation...");
    let pgb = ProgressBar::new(time_list.len() as u64);
    set_style(&pgb);
    pgb.inc(0);
    pgb.set_message(format!("at {} ns...", times[0]));

    // Frames whose PB/SA solve failed; they count as zero, which would
    // silently distort the average binding energy if nobody was told.
    let mut failed_pb_frames: Vec<f64> = Vec::new();
    for cur_frm in 0..time_list.len() {
        // MM
        if settings.calc_mm {
            let coord = coordinates.slice(s![cur_frm, .., ..]);
            if let Some(ndx_lig) = ndx_lig {
                let (de_elec, de_vdw) =
                    calc_mm(&ndx_rec, &ndx_lig, aps, &coord, &coeff, &settings, use_parallel);
                elec_atom.row_mut(cur_frm).assign(&de_elec);
                vdw_atom.row_mut(cur_frm).assign(&de_vdw);
            }
        }

        // PBSA
        if settings.calc_pbsa {
            let coord = coordinates.slice(s![cur_frm, .., ..]);
            match calc_pbsa(&coord, &times, ndx_rec, ndx_lig, cur_frm, sys_name, temp_dir,
                            aps, pbe_set, pba_set, settings) {
                Some((de_pb, de_sa)) => {
                    pb_atom.row_mut(cur_frm).assign(&de_pb);
                    sa_atom.row_mut(cur_frm).assign(&de_sa);
                }
                None => failed_pb_frames.push(times[cur_frm]),
            }
        }

        pgb.inc(1);
        pgb.set_message(format!("at {} ns, ΔH={:.2} kJ/mol, eta. {} s",
                                        times[cur_frm],
                                        vdw_atom.row(cur_frm).sum() + elec_atom.row(cur_frm).sum() +
                                        pb_atom.row(cur_frm).sum() + sa_atom.row(cur_frm).sum(),
                                        pgb.eta().as_secs()));
    }
    pgb.finish();

    if !failed_pb_frames.is_empty() {
        report_failed_pb_frames(&failed_pb_frames, time_list.len());
        // Every frame failed: the whole result would be zero-filled noise, so
        // there is no point in writing it out.
        if failed_pb_frames.len() == time_list.len() {
            println!("All frames failed: the PB/SA terms of the whole result would be zero. Aborting.");
            std::process::exit(1);
        }
        if settings.exit_on_error {
            println!("exit_on_error is enabled: aborting because of the PB/SA failures above.");
            std::process::exit(1);
        }
    }

    let mm_ie: Option<Array1<f64>> = if settings.inter_entropy {
        println!("Start IE calculation...");
        if settings.calc_mm {
            if let Some(ndx_lig) = ndx_lig {
                let pgb = ProgressBar::new(coordinates_ie.shape()[0] as u64);
                set_style(&pgb);
                pgb.inc(0);
                let calc_ie_per_frame = |frame: ArrayView2<f64>| {
                    let (de_elec, de_vdw) = calc_mm(&ndx_rec, &ndx_lig, aps, &frame, &coeff, &settings, false);
                    pgb.inc(1);
                    pgb.set_message(format!("eta. {} s", pgb.eta().as_secs()));
                    de_elec.sum() + de_vdw.sum()
                };
                let atoms_ie: Vec::<f64> = coordinates_ie.axis_iter(Axis(0))
                    .into_par_iter()
                    .enumerate()
                    .map(|(i, frame)| 
                        if let Some(frame_id) = times.iter().position(|&x| x == times_ie[i]) {
                            vdw_atom.row(frame_id).sum() + elec_atom.row(frame_id).sum()
                        } else {
                            calc_ie_per_frame(frame)
                        }).collect();
                pgb.finish();
                Some(Array1::from_vec(atoms_ie))
            } else {
                println!("Since no ligand, will not calculate IE.");
                None
            }
        } else {
            println!("Since MM calculation is not enabled, will not calculate IE.");
            None
        }
    } else {
        None
    };    

    // end calculation
    let t_end = Local::now();
    let t_spend = Duration::from(t_end - t_start).num_milliseconds() as f64 / 1000.0;
    println!("MM-PBSA calculation of {} finished. Total time cost: {:.2} s", sys_name, t_spend);
    fs::write(format!(".MMPBSA_time_cost_{}.txt", sys_name), format!("{:.2}", t_spend))
        .expect("Failed to write time cost file.");
    env::remove_var("OMP_NUM_THREADS");

    let atom_res = &aps.atom_props.iter().map(|a| a.resid).collect();
    let atom_names = &aps.atom_props.iter().map(|a| a.name.to_string()).collect();
    SMResult::new(
        atom_names,
        atom_res,
        residues,
        ndx_lig,
        &times,
        &times_ie,
        &coordinates,
        mutation,
        temperature,
        &elec_atom,
        &vdw_atom,
        &pb_atom,
        &sa_atom,
        &mm_ie,
    )
}

fn calc_mm(ndx_rec: &BTreeSet<usize>, ndx_lig: &BTreeSet<usize>, aps: &AtomProperties, coord: &ArrayView2<f64>, 
            coeff: &Coefficients, settings: &Settings, use_parallel: bool) -> (Array1<f64>, Array1<f64>) {
    let n_atoms = aps.atom_props.len();
    let mut de_elec: Array1<f64> = Array1::zeros(n_atoms);
    let mut de_vdw: Array1<f64> = Array1::zeros(n_atoms);
    
    // 预计算常量
    let r_cutoff_sq = settings.r_cutoff.powi(2);
    let a2nm = 1.0 / 10.0;
    
    // 预提取数据
    let rec_props: Vec<_> = ndx_rec.iter().map(|&i| {
        (i, aps.atom_props[i].charge, aps.atom_props[i].type_id, 
        coord[[i, 0]], coord[[i, 1]], coord[[i, 2]])
    }).collect();
    
    let lig_props: Vec<_> = ndx_lig.iter().map(|&j| {
        (j, aps.atom_props[j].charge, aps.atom_props[j].type_id, 
        coord[[j, 0]], coord[[j, 1]], coord[[j, 2]])
    }).collect();
    
    if use_parallel {
        // 并行版本
        let (par_elec, par_vdw): (Vec<f64>, Vec<f64>) = rec_props.par_iter()
            .map(|&(i, qi, ci, xi, yi, zi)| {
                let mut local_elec = vec![0.0; n_atoms];
                let mut local_vdw = vec![0.0; n_atoms];
                
                for &(j, qj, cj, xj, yj, zj) in &lig_props {
                    let dx = xi - xj;
                    let dy = yi - yj;
                    let dz = zi - zj;
                    let r_sq = dx * dx + dy * dy + dz * dz;
                    
                    if r_sq <= r_cutoff_sq {
                        let r = r_sq.sqrt() * a2nm;
                        let r_inv = 1.0 / r;
                        
                        let e_elec = qi * qj * r_inv * coefficients::screening_method(r, coeff, settings.elec_screen);
                        let r6_inv = r_inv.powi(6);
                        let e_vdw = (aps.c12[[ci, cj]] * r6_inv - aps.c6[[ci, cj]]) * r6_inv;
                        
                        local_elec[i] += e_elec;
                        local_elec[j] += e_elec;
                        local_vdw[i] += e_vdw;
                        local_vdw[j] += e_vdw;
                    }
                }
                
                (local_elec, local_vdw)
            })
            .reduce(|| (vec![0.0; n_atoms], vec![0.0; n_atoms]),
                    |(mut acc_elec, mut acc_vdw), (elec, vdw)| {
                        for i in 0..n_atoms {
                            acc_elec[i] += elec[i];
                            acc_vdw[i] += vdw[i];
                        }
                        (acc_elec, acc_vdw)
                    });
        
        de_elec = Array1::from_vec(par_elec);
        de_vdw = Array1::from_vec(par_vdw);
    } else {
        // 串行版本
        for &(i, qi, ci, xi, yi, zi) in &rec_props {
            for &(j, qj, cj, xj, yj, zj) in &lig_props {
                let dx = xi - xj;
                let dy = yi - yj;
                let dz = zi - zj;
                let r_sq = dx * dx + dy * dy + dz * dz;
                
                if r_sq <= r_cutoff_sq {
                    let r = r_sq.sqrt() * a2nm;
                    let r_inv = 1.0 / r;
                    
                    let e_elec = qi * qj * r_inv * coefficients::screening_method(r, coeff, settings.elec_screen);
                    let r6_inv = r_inv.powi(6);
                    let e_vdw = (aps.c12[[ci, cj]] * r6_inv - aps.c6[[ci, cj]]) * r6_inv;
                    
                    de_elec[i] += e_elec;
                    de_elec[j] += e_elec;
                    de_vdw[i] += e_vdw;
                    de_vdw[j] += e_vdw;
                }
            }
        }
    }
    
    // 应用缩放因子
    let elec_scale = coeff.f / coeff.pdie / 2.0;
    let vdw_scale = 0.5;
    
    de_elec = de_elec * elec_scale;
    de_vdw = de_vdw * vdw_scale;
    
    (de_elec, de_vdw)
}

fn calc_pbsa(coord: &ArrayView2<f64>, times: &Vec<f64>,
            ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>, cur_frm: usize, sys_name: &str, temp_dir: &PathBuf,
            aps: &AtomProperties, pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> Option<(Array1<f64>, Array1<f64>)> {
    // the default gamma parameter for apbs calculation is set to 1, in order to directly obtain the surface area
    // then the SA energy term is subsequently calculated
    let f_name = format!("{}_{}ns", sys_name, times[cur_frm]);
    if settings.calc_pbsa {
        // Assemble the APBS input and the molecule lists in memory; the PQR
        // files are only written (and the input text saved next to them) when
        // debug mode keeps intermediate files for inspection.
        let atom_radius: Array1<f64> = Array1::from_iter(aps.atom_props.iter().map(|a| a.radius));
        let input_text = build_apbs_input_text(ndx_rec, ndx_lig, coord,
                &atom_radius, pbe_set, pba_set, &f_name, settings);
        let mem_mols = build_molecules(aps, coord, ndx_rec, ndx_lig);
        if settings.debug_mode {
            prepare_pqr(cur_frm, &times, &temp_dir, sys_name, coord, &ndx_rec, ndx_lig, aps);
            let apbs_input = temp_dir.join(format!("{}.apbs", f_name));
            fs::write(&apbs_input, &input_text).expect("Failed to write apbs input file.");
        }
        // solve the PB/SA calculations in-process with the built-in APBS
        // solver (apbs-rs), no external process is spawned.
        let solver_run = match apbs_runner::run_apbs_in_process_text(&input_text, &mem_mols) {
            Ok(run) => run,
            Err(e) => {
                println!("PBSA calculation failed at {} ns: {}", times[cur_frm], e);
                return None;
            }
        };
        if settings.debug_mode {
            let mut outfile = File::create(temp_dir.join(format!("{}.out", f_name)))
                .expect("Failed to create output file.");
            outfile.write_all(solver_run.log.as_bytes()).expect("Failed to write apbs output.");
        };

        let apbs_results = ApbsResults::from_solver_run(&solver_run);
        let com_pb_sol = apbs_results.com_pb_sol;
        let com_pb_vac = apbs_results.com_pb_vac;
        let rec_pb_sol = apbs_results.rec_pb_sol;
        let rec_pb_vac = apbs_results.rec_pb_vac;
        let lig_pb_sol = apbs_results.lig_pb_sol;
        let lig_pb_vac = apbs_results.lig_pb_vac;
        let com_sa = apbs_results.com_sa;
        let rec_sa = apbs_results.rec_sa;
        let lig_sa = apbs_results.lig_sa;

        let com_pb: Array1<f64> = Array1::from_vec(com_pb_sol) - Array1::from_vec(com_pb_vac);
        let com_sa: Array1<f64> = Array1::from_vec(com_sa.par_iter()
            .map(|&sa| pba_set.gamma * sa + pba_set.bconc).collect());
        let mut rec_pb: Array1<f64> = Array1::from_vec(rec_pb_sol) - Array1::from_vec(rec_pb_vac);
        let mut rec_sa: Array1<f64> = Array1::from_vec(rec_sa.par_iter()
            .map(|&sa| pba_set.gamma * sa + pba_set.bconc).collect());
        let mut lig_pb: Array1<f64> = Array1::from_vec(lig_pb_sol) - Array1::from_vec(lig_pb_vac);
        let mut lig_sa: Array1<f64> = Array1::from_vec(lig_sa.par_iter()
            .map(|&sa| pba_set.gamma * sa + pba_set.bconc).collect());

        let (pb, sa) = if let Some(ndx_lig) = ndx_lig {
            if ndx_rec.iter().min() < ndx_lig.iter().min() {
                rec_pb.append(Axis(0), lig_pb.view()).unwrap();
                rec_sa.append(Axis(0), lig_sa.view()).unwrap();
                (com_pb - rec_pb, com_sa - rec_sa)
            } else {
                lig_pb.append(Axis(0), rec_pb.view()).unwrap();
                lig_sa.append(Axis(0), rec_sa.view()).unwrap();
                (com_pb - lig_pb, com_sa - lig_sa)
            }
        } else {
            (rec_pb, rec_sa)
        };
        if pb.len() != aps.atom_props.len() || sa.len() != aps.atom_props.len() {
            println!("Warning: The number of PB/SA values does not match the number of atoms. \
                This may indicate an issue with the in-process APBS solver results.");
            return None;
        }
        Some((pb, sa))
    } else {
        Some((Array1::zeros(aps.atom_props.len()), Array1::zeros(aps.atom_props.len())))
    }
}

/// Prints a prominent report of the frames whose PB/SA solve failed.
///
/// `exit_on_error` already decided to continue when this runs: the frames
/// are counted as zero, which is honest only if it is impossible to miss.
fn report_failed_pb_frames(failed_ns: &[f64], total: usize) {
    eprintln!();
    eprintln!("========================================================================");
    eprintln!("WARNING: PB/SA failed for {} of {} frames.", failed_ns.len(), total);
    eprintln!("Those frames are counted with zero PB/SA energy, so the PB/SA terms");
    eprintln!("and the average binding energy are NOT reliable.");
    let shown = failed_ns.iter().take(10)
        .map(|t| format!("{t}"))
        .collect::<Vec<_>>()
        .join(", ");
    if failed_ns.len() > 10 {
        eprintln!("Failed frames (ns): {} ... ({} more)", shown, failed_ns.len() - 10);
    } else {
        eprintln!("Failed frames (ns): {}", shown);
    }
    eprintln!("Set `exit_on_error = \"y\"` in settings.ini to stop instead.");
    eprintln!("========================================================================");
    eprintln!();
}

#[derive(Default)]
struct ApbsResults {
    com_pb_sol: Vec<f64>,
    com_pb_vac: Vec<f64>,
    rec_pb_sol: Vec<f64>,
    rec_pb_vac: Vec<f64>,
    lig_pb_sol: Vec<f64>,
    lig_pb_vac: Vec<f64>,
    com_sa: Vec<f64>,
    rec_sa: Vec<f64>,
    lig_sa: Vec<f64>,
}

impl ApbsResults {
    /// Bucket in-memory solver results into the per-complex PB/SA vectors.
    fn from_solver_run(run: &SolverRun) -> ApbsResults {
        let mut results = ApbsResults::default();
        for calc in &run.calcs {
            match calc.kind {
                SolverCalcKind::Elec => {
                    if calc.name.ends_with("com_SOL") {
                        results.com_pb_sol.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("com_VAC") {
                        results.com_pb_vac.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("rec_SOL") {
                        results.rec_pb_sol.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("rec_VAC") {
                        results.rec_pb_vac.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("lig_SOL") {
                        results.lig_pb_sol.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("lig_VAC") {
                        results.lig_pb_vac.extend_from_slice(&calc.per_atom);
                    }
                }
                SolverCalcKind::Apolar => {
                    if calc.name.ends_with("com_SAS") {
                        results.com_sa.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("rec_SAS") {
                        results.rec_sa.extend_from_slice(&calc.per_atom);
                    } else if calc.name.ends_with("lig_SAS") {
                        results.lig_sa.extend_from_slice(&calc.per_atom);
                    }
                }
            }
        }
        results
    }
}

fn transform_coordinate(base: &Array1<f64>, origin: &Array1<f64>, target_length: f64) -> Array1<f64> {
    let v_ch: Array1<f64> = origin - base;
    let cur_len = v_ch.iter().map(|d| d.powi(2)).sum::<f64>().sqrt();
    let lambda = target_length / cur_len;
    let new_v_ch: Array1<f64> = Array1::from_iter(v_ch.iter().map(|r| r * lambda));
    return base + new_v_ch
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::apbs_runner::{SolverCalcResult, SolverRun};

    #[test]
    fn bucket_solver_results_by_calc_name() {
        let run = SolverRun {
            calcs: vec![
                SolverCalcResult {
                    name: "sys_com_SOL".to_string(),
                    kind: SolverCalcKind::Elec,
                    per_atom: vec![11.0, 22.0],
                },
                SolverCalcResult {
                    name: "sys_com_VAC".to_string(),
                    kind: SolverCalcKind::Elec,
                    per_atom: vec![9.9, 8.8],
                },
                SolverCalcResult {
                    name: "sys_rec_SOL".to_string(),
                    kind: SolverCalcKind::Elec,
                    per_atom: vec![3.3],
                },
                SolverCalcResult {
                    name: "sys_rec_VAC".to_string(),
                    kind: SolverCalcKind::Elec,
                    per_atom: vec![2.2],
                },
                SolverCalcResult {
                    name: "sys_lig_SOL".to_string(),
                    kind: SolverCalcKind::Elec,
                    per_atom: vec![50.0],
                },
                SolverCalcResult {
                    name: "sys_lig_VAC".to_string(),
                    kind: SolverCalcKind::Elec,
                    per_atom: vec![7.7],
                },
                SolverCalcResult {
                    name: "sys_com_SAS".to_string(),
                    kind: SolverCalcKind::Apolar,
                    per_atom: vec![44.0, 45.0],
                },
                SolverCalcResult {
                    name: "sys_rec_SAS".to_string(),
                    kind: SolverCalcKind::Apolar,
                    per_atom: vec![6.6],
                },
                SolverCalcResult {
                    name: "sys_lig_SAS".to_string(),
                    kind: SolverCalcKind::Apolar,
                    per_atom: vec![7.7],
                },
            ],
            log: String::new(),
        };
        let r = ApbsResults::from_solver_run(&run);
        assert_eq!(r.com_pb_sol, vec![11.0, 22.0]);
        assert_eq!(r.com_pb_vac, vec![9.9, 8.8]);
        assert_eq!(r.rec_pb_sol, vec![3.3]);
        assert_eq!(r.rec_pb_vac, vec![2.2]);
        assert_eq!(r.lig_pb_sol, vec![50.0]);
        assert_eq!(r.lig_pb_vac, vec![7.7]);
        assert_eq!(r.com_sa, vec![44.0, 45.0]);
        assert_eq!(r.rec_sa, vec![6.6]);
        assert_eq!(r.lig_sa, vec![7.7]);
    }

    fn aps_from_atoms(atoms: &[(&str, &str)], resid: usize) -> AtomProperties {
        AtomProperties {
            c6: Array2::zeros((0, 0)),
            c12: Array2::zeros((0, 0)),
            at_map: HashMap::new(),
            radius_type: "mBondi".to_string(),
            atom_props: atoms
                .iter()
                .enumerate()
                .map(|(i, &(name, resname))| AtomProperty {
                    charge: 0.0,
                    radius: 1.4,
                    type_id: 0,
                    id: i,
                    name: name.to_string(),
                    at_type: "C".to_string(),
                    resname: resname.to_string(),
                    resid,
                })
                .collect(),
        }
    }

    /// `ala_mutate` keeps the backbone and CB, renames the side-chain gamma
    /// atom to HB3 with its coordinate rebuilt along CB -> Xg at 1.09 A, and
    /// deletes the rest of the side chain — on a per-atom subset, leaving
    /// every other frame coordinate untouched.
    #[test]
    fn ala_mutate_leu_rebuilds_hb_and_drops_rest_of_sidechain() {
        let aps = aps_from_atoms(
            &[("N", "LEU"), ("CA", "LEU"), ("C", "LEU"), ("O", "LEU"),
              ("CB", "LEU"), ("CG", "LEU"), ("CD1", "LEU")],
            0,
        );
        // Two frames; CB -> CG is 1.54 A along +z in both.
        let coordinates = Array3::from_shape_vec(
            (2, 7, 3),
            vec![
                0.0, 0.0, 0.0,  0.0, 0.0, 0.5,  0.0, 0.0, 1.0,  0.0, 0.0, 1.5,
                0.0, 0.0, 2.0,  0.0, 0.0, 3.54, 1.0, 0.0, 3.54,
                1.0, 2.0, 3.0,  1.0, 2.0, 3.5,  1.0, 2.0, 4.0,  1.0, 2.0, 4.5,
                1.0, 2.0, 5.0,  1.0, 2.0, 6.54, 2.0, 2.0, 6.54,
            ],
        )
        .unwrap();
        let exclude_list = ["N", "CA", "C", "O", "CB", "H", "HA", "HB1", "HB2"];
        let asr = Residue { id: 0, name: "LEU".to_string(), nr: 1 };
        let ndx_rec: BTreeSet<usize> = (0..7).collect();

        let (new_coordinates, new_aps, new_ndx_rec, new_ndx_lig) =
            ala_mutate(&aps, &asr, &exclude_list, &coordinates, &ndx_rec, &None);

        // CD1 is gone; CG survives as the new HB3.
        assert_eq!(new_coordinates.dim(), (2, 6, 3));
        let names: Vec<&str> = new_aps.atom_props.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, ["N", "CA", "C", "O", "CB", "CG"]);
        assert_eq!(new_aps.atom_props.iter().map(|a| a.id).collect::<Vec<_>>(), (0..6).collect::<Vec<_>>());
        assert_eq!(new_ndx_rec, (0..6).collect::<BTreeSet<usize>>());
        assert_eq!(new_ndx_lig, None);

        // Backbone and CB coordinates are untouched and in the expected slots.
        for layer in 0..2 {
            for atom in 0..5 {
                for d in 0..3 {
                    assert_eq!(
                        new_coordinates[[layer, atom, d]],
                        coordinates[[layer, atom, d]],
                        "frame {layer} atom {atom} coord {d}"
                    );
                }
            }
        }
        // The new HB sits 1.09 A from CB along the old CB -> CG direction.
        for layer in 0..coordinates.shape()[0] {
            let cb = new_coordinates.slice(s![layer, 4, ..]);
            let hb = new_coordinates.slice(s![layer, 5, ..]);
            let dist: f64 = (0..3usize).map(|d| (hb[d] - cb[d]).powi(2)).sum::<f64>().sqrt();
            assert!((dist - 1.09).abs() < 1e-12, "frame {layer}: |HB-CB| = {dist}");
            let cg_old = coordinates.slice(s![layer, 5, ..]);
            let dot: f64 = (0..3usize)
                .map(|d| (hb[d] - cb[d]) * (cg_old[d] - cb[d]))
                .sum();
            assert!(
                dot > 0.0,
                "frame {layer}: HB must lie on the CB -> CG ray"
            );
        }
    }

    /// Proline keeps both ring atoms: CG becomes HB3 (along CB -> CG) and CD
    /// becomes H (along N -> CD).
    #[test]
    fn ala_mutate_pro_rebuilds_both_ring_hydrogens() {
        let aps = aps_from_atoms(
            &[("N", "PRO"), ("CA", "PRO"), ("C", "PRO"), ("O", "PRO"),
              ("CB", "PRO"), ("CG", "PRO"), ("CD", "PRO")],
            0,
        );
        let coordinates = Array3::from_shape_vec(
            (1, 7, 3),
            vec![
                0.0, 0.0, 0.0,  0.0, 0.0, 0.5,  0.0, 0.0, 1.0,  0.0, 0.0, 1.5,
                1.0, 0.0, 0.0,  1.0, 0.0, 1.54, 0.0, 0.0, 1.5,
            ],
        )
        .unwrap();
        let exclude_list = ["N", "CA", "C", "O", "CB", "H", "HA", "HB1", "HB2"];
        let asr = Residue { id: 0, name: "PRO".to_string(), nr: 1 };
        let ndx_rec: BTreeSet<usize> = (0..7).collect();

        let (new_coordinates, _new_aps, _, _) =
            ala_mutate(&aps, &asr, &exclude_list, &coordinates, &ndx_rec, &None);

        // Nothing is deleted: both ring atoms turn into hydrogens.
        assert_eq!(new_coordinates.dim(), (1, 7, 3));
        // HB3: 1.09 A from CB along the old CB -> CG direction ([1, 0, 1.09]).
        for d in 0..3 {
            let want = [1.0, 0.0, 1.09][d];
            assert!((new_coordinates[[0, 5, d]] - want).abs() < 1e-12);
        }
        // H: 1.07 A from N along the old N -> CD direction ([0, 0, 1.07]).
        for d in 0..3 {
            let want = [0.0, 0.0, 1.07][d];
            assert!((new_coordinates[[0, 6, d]] - want).abs() < 1e-12);
        }
    }
}
