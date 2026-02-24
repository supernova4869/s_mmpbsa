use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use crate::fun_para_system::normalize_index;
use crate::settings::Settings;
use crate::utils::{self, is_amino};
use ndarray::parallel::prelude::*;
use ndarray::{s, Array1, Array2, Array3, ArrayView2, Axis};
use std::process::Command;
use std::env;
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{Local, Duration};
use crate::coefficients::{self, Coefficients};
use crate::analyzation::SMResult;
use crate::parse_tpr::Residue;
use crate::parameters::{PBASet, PBESet};
use crate::atom_property::{AtomProperties, AtomProperty};
use crate::prepare_apbs::{prepare_pqr, write_apbs_input};

pub fn fun_mmpbsa_calculations(time_list: &Vec<f64>, time_list_ie: &Vec<f64>, coordinates_ie: &Array3<f64>, 
                               temp_dir: &PathBuf, sys_name: &str, aps: &AtomProperties,
                               ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>,
                               ala_list: &Vec<i32>, residues: &Vec<Residue>, wd: &Path, temperature: f64,
                               pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings)
                               -> (SMResult, Vec<SMResult>) {
    println!("Running MM-PBSA calculations of {}...", sys_name);
    if ala_list.len() > 0 {
        let as_res: Vec<String> = residues.iter().filter_map(
            |r| if ala_list.contains(&r.nr) && is_amino(&r.name) {
                match utils::resname_3to1(&r.name) {
                    Some(mutation) => Some(format!("{}{}A", mutation, r.nr)),
                    None => Some(format!("{}{}A", r.name.to_string(), r.nr))
                }
            } else {
                None
            }).collect();
        println!("Mutations for alanine scanning: {}", as_res.join(", "));
    }

    // calculate MM and PBSA
    println!("Calculating binding energy for {}...", sys_name);
    let result_wt = calculate_mmpbsa(time_list, time_list_ie, coordinates_ie, aps, &temp_dir, 
        &ndx_rec, &ndx_lig, residues, temperature,
        sys_name, "WT", pbe_set, pba_set, settings);
    result_wt.to_bin(&wd.join(format!("_MMPBSA_{}_{}.sm", sys_name, "WT").as_str()));

    let mut result_ala_scan: Vec<SMResult> = vec![];
    if ala_list.len() > 0 {
        // main chain atoms number
        let as_res: Vec<&Residue> = residues.iter().filter(|&r| ala_list.contains(&r.nr) 
            && is_amino(&r.name)).collect();     // gly not contain CB, ala no need to mutate
        let exclude_list = ["N", "CA", "C", "O", "CB", "HN", "HCA", "HCB"];
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

            let mutation = format!("{}{}A", mutation, asr.nr);
            let sys_name = format!("{}_{}", sys_name, mutation);
            println!("Calculating binding energy for {}...", sys_name);
            let result_as = calculate_mmpbsa(time_list, time_list_ie, &new_coordinates_ie,
                &new_aps, &temp_dir, &new_ndx_rec, &new_ndx_lig, &new_residues, temperature,
                &sys_name, &mutation, pbe_set, pba_set, settings);
            result_as.to_bin(&wd.join(format!("_MMPBSA_{}.sm", sys_name).as_str()));
            result_ala_scan.push(result_as);
        }
    };

    // whether remove temp directory
    if !settings.debug_mode {
        if settings.apbs_path.is_some() {
            fs::remove_dir_all(&temp_dir).expect("Remove dir failed");
        }
        fs::remove_file(wd.join("_MMPBSA_coord_ie.xvg")).ok();
        fs::remove_file(wd.join("_MMPBSA_coord_sol.xvg")).ok();
    }

    println!("");
    utils::show_famous_quotes();

    (result_wt, result_ala_scan)
}

fn ala_mutate(aps: &AtomProperties, asr: &Residue, exclude_list: &[&str], coordinates: &Array3<f64>, 
                ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>)
                -> (Array3<f64>, AtomProperties, BTreeSet<usize>, Option<BTreeSet<usize>>) {
    let mut new_coordinates = coordinates.clone();
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
    for xg in xgs.iter() {
        new_aps.atom_props[xg.id].change_atom(aps.at_map.get("HC"), "HC", &aps.radius_type);
        // 获取新的HB坐标
        for layer in 0..new_coordinates.shape()[0] {
            let cb_coords: Array1<f64> = new_coordinates.slice(s![layer, cb[0].id, ..]).to_owned();
            let hg_coords: Array1<f64> = new_coordinates.slice(s![layer, xg.id, ..]).to_owned();
            let new_hg_coord: Array1<f64> = transform_coordinate(&cb_coords, &hg_coords, 1.09);
            new_coordinates[[layer, xg.id, 0]] = new_hg_coord[0];
            new_coordinates[[layer, xg.id, 1]] = new_hg_coord[1];
            new_coordinates[[layer, xg.id, 2]] = new_hg_coord[2];
        }
    }
    // 脯氨酸需要把CD改成H
    // 通过N定位新的HN
    if asr.name.eq("PRO") {
        let n: Vec<AtomProperty> = as_atoms.iter().filter_map(|a| {
            if a.name.eq("N") {
                Some(a.clone())
            } else {
                None
            }}).collect();
        sc_out.retain(|&a| a.name.ne("CD"));
        let cd = as_atoms.iter().find(|&a| a.name == "CD").unwrap();
        new_aps.atom_props[cd.id].change_atom(aps.at_map.get("H"), "HN", &aps.radius_type);
        // 获取新的HN坐标
        for layer in 0..new_coordinates.shape()[0] {
            let n_coords: Array1<f64> = new_coordinates.slice(s![layer, n[0].id, ..]).to_owned();
            let hn_coords: Array1<f64> = new_coordinates.slice(s![layer, cd.id, ..]).to_owned();
            let new_hn_coord = transform_coordinate(&n_coords, &hn_coords, 1.07);
            new_coordinates[[layer, cd.id, 0]] = new_hn_coord[0];
            new_coordinates[[layer, cd.id, 1]] = new_hn_coord[1];
            new_coordinates[[layer, cd.id, 2]] = new_hn_coord[2];
        }
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
    let mut new_ndx_rec = ndx_rec.clone();
    new_ndx_rec.retain(|&x| !del_list.contains(&x) || xg_list.contains(&x));
    let (new_ndx_rec, new_ndx_lig) = normalize_index(&new_ndx_rec, ndx_lig);

    return (new_coordinates.select(Axis(1), &retain_id).clone(), new_aps, new_ndx_rec, new_ndx_lig)
}

pub fn set_style(pb: &ProgressBar) {
    pb.set_style(ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:50.cyan/cyan} {pos}/{len} {msg}").unwrap()
        .progress_chars("=>-"));
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
    if settings.pbsa_kernel.is_some() && settings.apbs_path.is_some() {
        let apbs_path = settings.apbs_path.as_ref().unwrap();
        let apbs_dir = Path::new(apbs_path).parent().unwrap();
        
        // 处理动态库路径
        if cfg!(target_os = "windows") {
            let current_path = env::var("PATH").unwrap_or_default();
            let new_path = format!("{};{}", apbs_dir.display(), current_path);
            env::set_var("PATH", new_path);
        } else {
            env::set_var("LD_LIBRARY_PATH", 
                format!("{}:{}", Path::new(settings.apbs_path.as_ref().unwrap()).parent().unwrap().to_str().unwrap(), 
                env::var("LD_LIBRARY_PATH").unwrap()));
        }
    }

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

    for cur_frm in 0..time_list.len() {
        // MM
        let coord = coordinates.slice(s![cur_frm, .., ..]);
        if let Some(ndx_lig) = ndx_lig {
            let (de_elec, de_vdw) = 
                calc_mm(&ndx_rec, &ndx_lig, aps, &coord, &coeff, &settings, use_parallel);
            elec_atom.row_mut(cur_frm).assign(&de_elec);
            vdw_atom.row_mut(cur_frm).assign(&de_vdw);
        }

        // PBSA
        let coord = coordinates.slice(s![cur_frm, .., ..]);
        if settings.pbsa_kernel.is_some() {
            let (de_pb, de_sa) = 
                calc_pbsa(&coord, &times, ndx_rec, ndx_lig, cur_frm, sys_name, temp_dir, aps, pbe_set, pba_set, settings);
            pb_atom.row_mut(cur_frm).assign(&de_pb);
            sa_atom.row_mut(cur_frm).assign(&de_sa);
        }

        pgb.inc(1);
        pgb.set_message(format!("at {} ns, ΔPBSA={:.2} kJ/mol, eta. {} s", 
                                        times[cur_frm],
                                        vdw_atom.row(cur_frm).sum() + elec_atom.row(cur_frm).sum() + 
                                        pb_atom.row(cur_frm).sum() + sa_atom.row(cur_frm).sum(),
                                        pgb.eta().as_secs()));
    }
    pgb.finish();

    let mm_ie: Array1<f64> = if let Some(ndx_lig) = ndx_lig {
        println!("Start IE calculation...");
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
        Array1::from_vec(atoms_ie)
    } else {
        Array1::zeros(coordinates_ie.shape()[0])
    };

    // end calculation
    let t_end = Local::now();
    let t_spend = Duration::from(t_end - t_start).num_milliseconds() as f64 / 1000.0;
    println!("MM-PBSA calculation of {} finished. Total time cost: {:.2} s", sys_name, t_spend);
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
    let scale_factor = 1.0 / 10.0;
    
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
                        let r = r_sq.sqrt() * scale_factor;
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
                    let r = r_sq.sqrt() * scale_factor;
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
            aps: &AtomProperties, pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> (Array1<f64>, Array1<f64>) {
    prepare_pqr(cur_frm, &times, &temp_dir, sys_name, coord, &ndx_rec, ndx_lig, aps);
                
    // From AMBER-PB4, the surface extension constant γ=0.0072 kcal/(mol·Å2)=0.030125 kJ/(mol·Å^2)
    // but the default gamma parameter for apbs calculation is set to 1, in order to directly obtain the surface area
    // then the SA energy term is calculated by s_mmpbsa
    let gamma = 0.030125;
    let bias = 0.0;
    let f_name = format!("{}_{}ns", sys_name, times[cur_frm]);
    if let Some(pbsa_kernel) = &settings.pbsa_kernel {
        if pbsa_kernel.eq("apbs") {
            if settings.apbs_path.as_ref().is_none() {
                return (Array1::zeros(aps.atom_props.len()), Array1::zeros(aps.atom_props.len()))
            }
            let apbs = settings.apbs_path.as_ref().unwrap();
            write_apbs_input(ndx_rec, ndx_lig, coord, &Array1::from_iter(aps.atom_props.iter().map(|a| a.radius)),
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
            // let apbs_result = fs::read_to_string(temp_dir.join(format!("{}.out", f_name))).expect("Failed to get apbs output file.");
    
            // preserve CALCULATION, Atom and SASA lines
            let apbs_result: Vec<&str> = apbs_result.split("\n").filter_map(|p|
                if p.trim().starts_with("CALCULATION") || p.trim().starts_with("Atom") || p.trim().starts_with("SASA") {
                    Some(p.trim())
                } else {
                    None
                }
            ).collect();
    
            // extract apbs results
            let indexes: Vec<usize> = apbs_result
                .iter()
                .enumerate()
                .filter_map(|(i, p)| p.starts_with("CAL").then_some(i))
                .collect();

            let mut com_pb_sol: Vec<f64> = vec![];
            let mut com_pb_vac: Vec<f64> = vec![];
            let mut rec_pb_sol: Vec<f64> = vec![];
            let mut rec_pb_vac: Vec<f64> = vec![];
            let mut lig_pb_sol: Vec<f64> = vec![];
            let mut lig_pb_vac: Vec<f64> = vec![];
            let mut com_sa: Vec<f64> = vec![];
            let mut rec_sa: Vec<f64> = vec![];
            let mut lig_sa: Vec<f64> = vec![];
    
            let mut skip_pb = true;     // the first time PB calculation should be dropped
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
            let com_sa: Array1<f64> = Array1::from_vec(com_sa.par_iter().map(|&sa| gamma * sa + bias).collect());
            let mut rec_pb: Array1<f64> = Array1::from_vec(rec_pb_sol) - Array1::from_vec(rec_pb_vac);
            let mut rec_sa: Array1<f64> = Array1::from_vec(rec_sa.par_iter().map(|&sa| gamma * sa + bias).collect());
            let mut lig_pb: Array1<f64> = Array1::from_vec(lig_pb_sol) - Array1::from_vec(lig_pb_vac);
            let mut lig_sa: Array1<f64> = Array1::from_vec(lig_sa.par_iter().map(|&sa| gamma * sa + bias).collect());
    
            if let Some(ndx_lig) = ndx_lig {
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
            }
        }
        else {
            println!("Currently Delphi kernel not available.");
            (Array1::zeros(aps.atom_props.len()), Array1::zeros(aps.atom_props.len()))
        }
    } else {
        (Array1::zeros(aps.atom_props.len()), Array1::zeros(aps.atom_props.len()))
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

fn transform_coordinate(base: &Array1<f64>, origin: &Array1<f64>, target_length: f64) -> Array1<f64> {
    let v_ch: Array1<f64> = origin - base;
    let cur_len = v_ch.iter().map(|d| d.powi(2)).sum::<f64>().sqrt();
    let lambda = target_length / cur_len;
    let new_v_ch: Array1<f64> = Array1::from_iter(v_ch.iter().map(|r| r * lambda));
    return base + new_v_ch
}