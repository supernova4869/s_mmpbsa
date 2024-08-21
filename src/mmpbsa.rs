use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use crate::fun_para_system::normalize_index;
// use xdrfile::*;
use crate::settings::Settings;
use crate::utils::resname_3to1;
use ndarray::parallel::prelude::*;
use ndarray::{s, Array1, Array2, Array3, ArrayView2, Axis};
use std::process::Command;
// use std::rc::Rc;
use std::env;
use indicatif::{ProgressBar, ProgressStyle};
use chrono::{Local, Duration};
use crate::coefficients::Coefficients;
use crate::analyzation::Results;
use crate::parse_tpr::Residue;
use crate::apbs_param::{PBASet, PBESet};
use crate::atom_property::{AtomProperties, AtomProperty};
use crate::prepare_apbs::{prepare_pqr, write_apbs_input};

pub fn fun_mmpbsa_calculations(time_list: &Vec<f64>, coordinates: &Array3<f64>, temp_dir: &PathBuf,
                               sys_name: &String, aps: &AtomProperties,
                               ndx_rec: &Vec<usize>, ndx_lig: &Vec<usize>, 
                               ala_list: &Vec<i32>, residues: &Vec<Residue>, wd: &Path,
                               bf: usize, ef: usize, dframe: usize, total_frames: usize,
                               pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings)
                               -> (Results, Vec<Results>) {
    println!("Running MM/PB-SA calculations of {}...", sys_name);
    if ala_list.len() > 0 {
        let as_res: Vec<String> = residues.iter().filter_map(
            |r| if ala_list.contains(&r.nr) && r.name.ne("GLY") && r.name.ne("ALA") {
                match resname_3to1(&r.name) {
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
    let result_wt = calculate_mmpbsa(&time_list, &coordinates, bf, ef, dframe, 
        total_frames, aps, &temp_dir, &ndx_rec, &ndx_lig, residues,
        sys_name, "WT", pbe_set, pba_set, settings);

    let mut result_ala_scan: Vec<Results> = vec![];
    if ala_list.len() > 0 {
        // main chain atoms number
        let as_res: Vec<&Residue> = residues.iter().filter(|&r| ala_list.contains(&r.nr) 
            && r.name.ne("GLY") && r.name.ne("ALA")).collect();     // gly not contain CB, ala no need to mutate
        let mut new_coordinates = coordinates.clone();
        let exclude_list = ["N", "CA", "C", "O", "CB", "HN", "HCA", "HCB"];
        for asr in as_res {
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
            let new_coordinates: Array3<f64> = new_coordinates.select(Axis(1), &retain_id);
            let mut new_ndx_rec = ndx_rec.clone();
            new_ndx_rec.retain(|&x| !del_list.contains(&x) || xg_list.contains(&x));
            let (new_ndx_rec, new_ndx_lig) = normalize_index(&new_ndx_rec, Some(ndx_lig));

            // After alanine mutation
            let mutation = match resname_3to1(&asr.name) {
                Some(mutation) => mutation,
                None => asr.name.to_string()
            };

            let mut new_residues = residues.clone();
            new_residues[asr.id].name = "ALA".to_string();

            let mutation = format!("{}{}A", mutation, asr.nr);
            let sys_name = format!("{}-{}", sys_name, mutation);
            println!("Calculating binding energy for {}...", sys_name);
            let result_as = calculate_mmpbsa(&time_list, &new_coordinates,
                bf, ef, dframe, total_frames, &new_aps, &temp_dir, 
                &new_ndx_rec, &new_ndx_lig, &new_residues,
                &sys_name, &mutation, pbe_set, pba_set, settings);
            result_ala_scan.push(result_as);
        }
    };

    // whether remove temp directory
    if !settings.debug_mode {
        if settings.apbs_path.is_some() {
            fs::remove_dir_all(&temp_dir).expect("Remove dir failed");
        }
        fs::remove_file(wd.join("_MMPBSA_coord.xvg")).unwrap();
    }

    (result_wt, result_ala_scan)
}

pub fn set_style(pb: &ProgressBar) {
    pb.set_style(ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:50.cyan/cyan} {pos}/{len} {msg}").unwrap()
        .progress_chars("=>-"));
}

fn calculate_mmpbsa(time_list: &Vec<f64>, coordinates: &Array3<f64>, bf: usize, ef: usize, 
                    dframe: usize, total_frames: usize, aps: &AtomProperties, temp_dir: &PathBuf,
                    ndx_rec: &Vec<usize>, ndx_lig: &Vec<usize>,
                    residues: &Vec<Residue>, sys_name: &String, mutation: &str,
                    pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> Results {
    let mut elec_atom: Array2<f64> = Array2::zeros((total_frames, aps.atom_props.len()));
    let mut vdw_atom: Array2<f64> = Array2::zeros((total_frames, aps.atom_props.len()));
    let mut pb_atom: Array2<f64> = Array2::zeros((total_frames, aps.atom_props.len()));
    let mut sa_atom: Array2<f64> = Array2::zeros((total_frames, aps.atom_props.len()));
    
    // parameters for elec calculation
    let coeff = Coefficients::new(pbe_set);

    // Time list of trajectory
    let times: Array1<f64> = Array1::from_iter((bf..=ef).into_iter().step_by(dframe).map(|f| time_list[f] / 1000.0));

    // start calculation
    env::set_var("OMP_NUM_THREADS", settings.nkernels.to_string());
    let t_start = Local::now();
    
    let pgb = ProgressBar::new(total_frames as u64);
    set_style(&pgb);
    pgb.inc(0);
    let mut frame_id = 0;
    pgb.set_message(format!("at {} ns...", times[frame_id]));
    for cur_frm in (bf..=ef).step_by(dframe) {
        // MM
        let coord = coordinates.slice(s![cur_frm, .., ..]);
        if ndx_lig[0] != ndx_rec[0] {
            let (de_elec, de_vdw) = 
                calc_mm(&ndx_rec, &ndx_lig, aps, &coord, &coeff, &settings);
            elec_atom.row_mut(frame_id).assign(&de_elec);
            vdw_atom.row_mut(frame_id).assign(&de_vdw);
        }

        // PBSA
        if settings.pbsa_kernel.is_some() {
            let (de_pb, de_sa) = 
                calc_pbsa(&coord, time_list, ndx_rec, ndx_lig, cur_frm, sys_name, temp_dir, aps, pbe_set, pba_set, settings);
            pb_atom.row_mut(frame_id).assign(&de_pb);
            sa_atom.row_mut(frame_id).assign(&de_sa);
        }

        pgb.inc(1);
        pgb.set_message(format!("at {} ns, ΔH={:.2} kJ/mol, eta. {} s", 
                                        times[frame_id],
                                        vdw_atom.row(frame_id).sum() + elec_atom.row(frame_id).sum() + 
                                        pb_atom.row(frame_id).sum() + sa_atom.row(frame_id).sum(),
                                        pgb.eta().as_secs()));

        frame_id += 1;
    }
    pgb.finish();

    // end calculation
    let t_end = Local::now();
    let t_spend = Duration::from(t_end - t_start).num_milliseconds();
    println!("MM/PB-SA calculation of {} finished. Total time cost: {} s", sys_name, t_spend as f64 / 1000.0);
    env::remove_var("OMP_NUM_THREADS");

    let fr: Vec<usize> = (bf..=ef).step_by(dframe).collect();
    Results::new(
        aps,
        residues,
        ndx_lig,
        &times,
        coordinates.select(Axis(0), &fr),
        mutation,
        &elec_atom,
        &vdw_atom,
        &pb_atom,
        &sa_atom,
    )
}

fn calc_mm(ndx_rec: &Vec<usize>, ndx_lig: &Vec<usize>, aps: &AtomProperties, coord: &ArrayView2<f64>, 
            coeff: &Coefficients, settings: &Settings) -> (Array1<f64>, Array1<f64>) {
    let mut de_elec: Array1<f64> = Array1::zeros(aps.atom_props.len());
    let mut de_vdw: Array1<f64> = Array1::zeros(aps.atom_props.len());

    for &i in ndx_rec {
        let qi = aps.atom_props[i].charge;
        let ci = aps.atom_props[i].type_id;
        let xi = coord[[i, 0]];
        let yi = coord[[i, 1]];
        let zi = coord[[i, 2]];
        for &j in ndx_lig {
            if ndx_lig[0] == ndx_rec[0] && j <= i {
                continue;
            }
            let qj = aps.atom_props[j].charge;
            let cj = aps.atom_props[j].type_id;
            let xj = coord[[j, 0]];
            let yj = coord[[j, 1]];
            let zj = coord[[j, 2]];
            let r = ((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2)).sqrt();
            if r <= settings.r_cutoff {
                let e_elec = match settings.use_dh {
                    false => qi * qj / r,
                    true => qi * qj / r * (-coeff.kap * r).exp()   // doi: 10.1088/0256-307X/38/1/018701
                }; // use A for elec
                // use nm for vdW
                let r = r / 10.0;
                let e_vdw = (aps.c12[[ci, cj]] / r.powi(6) - aps.c6[[ci, cj]]) / r.powi(6);
                de_elec[i] += e_elec;
                de_elec[j] += e_elec;
                de_vdw[i] += e_vdw;
                de_vdw[j] += e_vdw;
            }
        }
    }

    de_elec = de_elec * coeff.kj_elec / coeff.pdie / 2.0;
    de_vdw = de_vdw / 2.0;

    return (de_elec, de_vdw)
}

fn calc_pbsa(coord: &ArrayView2<f64>, time_list: &Vec<f64>, 
            ndx_rec_norm: &Vec<usize>, ndx_lig_norm: &Vec<usize>, cur_frm: usize, sys_name: &String, temp_dir: &PathBuf, 
            aps: &AtomProperties, pbe_set: &PBESet, pba_set: &PBASet, settings: &Settings) -> (Array1<f64>, Array1<f64>) {
    prepare_pqr(cur_frm, &time_list, &temp_dir, sys_name, coord, &ndx_rec_norm, ndx_lig_norm, aps);

    // From AMBER-PB4, the surface extension constant γ=0.0072 kcal/(mol·Å2)=0.030125 kJ/(mol·Å^2)
    // but the default gamma parameter for apbs calculation is set to 1, in order to directly obtain the surface area
    // then the SA energy term is calculated by s_mmpbsa
    let gamma = 0.030125;
    let bias = 0.0;
    let f_name = format!("{}_{}ns", sys_name, time_list[cur_frm]);
    if let Some(pbsa_kernel) = &settings.pbsa_kernel {
        if pbsa_kernel.eq("apbs") {
            let apbs = settings.apbs_path.as_ref().unwrap();
            write_apbs_input(ndx_rec_norm, ndx_lig_norm, coord, &Array1::from_iter(aps.atom_props.iter().map(|a| a.radius)),
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
            let mut rec_pb: Array1<f64> = Array1::from_vec(rec_pb_sol) - Array1::from_vec(rec_pb_vac);
            let mut rec_sa: Array1<f64> = Array1::from_vec(rec_sa.par_iter().map(|i| gamma * *i + bias / rec_sa.len() as f64).collect());
            let mut lig_pb: Array1<f64> = Array1::from_vec(lig_pb_sol) - Array1::from_vec(lig_pb_vac);
            let mut lig_sa: Array1<f64> = Array1::from_vec(lig_sa.par_iter().map(|i| gamma * *i + bias / lig_sa.len() as f64).collect());
    
            if ndx_rec_norm[0] < ndx_lig_norm[0] {
                rec_pb.append(Axis(0), lig_pb.view()).unwrap();
                rec_sa.append(Axis(0), lig_sa.view()).unwrap();
                return (com_pb - rec_pb, com_sa - rec_sa)
            } else if ndx_rec_norm[0] > ndx_lig_norm[0] {
                lig_pb.append(Axis(0), rec_pb.view()).unwrap();
                lig_sa.append(Axis(0), rec_sa.view()).unwrap();
                return (com_pb - lig_pb, com_sa - lig_sa)
            } else {
                return (rec_pb, rec_sa)
            }
        }
        else {
            println!("Currently Delphi kernel not available.");
            return (Array1::zeros(aps.atom_props.len()), Array1::zeros(aps.atom_props.len()))
        }
    } else {
        return (Array1::zeros(aps.atom_props.len()), Array1::zeros(aps.atom_props.len()))
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