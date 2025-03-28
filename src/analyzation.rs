use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::process::{exit, Command, Stdio};
use ndarray::{s, Array1, Array2, Array3, Axis};
use serde::{Deserialize, Serialize};
use plotpy::{Barplot, Curve, Plot};
use crate::parse_tpr::Residue;
use crate::settings::Settings;
use crate::utils::{get_input, get_input_selection, get_residue_range_ca, range2list};

#[derive(Clone, Serialize, Deserialize)]
pub struct SMResult {
    pub mutation: String,
    pub atom_names: Vec<String>,
    pub atom_res: Vec<usize>,
    pub residues: Vec<Residue>,
    pub ndx_lig: Vec<usize>,
    pub times: Vec<f64>,
    pub times_ie: Vec<f64>,
    pub coord: Array3<f64>,
    pub dh: Array1<f64>,
    pub mm: Array1<f64>,
    pub pb: Array1<f64>,
    pub sa: Array1<f64>,
    pub elec: Array1<f64>,
    pub vdw: Array1<f64>,
    pub mm_ie: Array1<f64>,
    pub dh_atom: Array2<f64>,
    pub mm_atom: Array2<f64>,
    pub pb_atom: Array2<f64>,
    pub sa_atom: Array2<f64>,
    pub elec_atom: Array2<f64>,
    pub vdw_atom: Array2<f64>,
}

impl SMResult {
    pub fn new(atom_names: &Vec<String>, atom_res: &Vec<usize>, 
               residues: &Vec<Residue>, ndx_lig: &Vec<usize>, 
               times: &Vec<f64>, times_ie: &Vec<f64>, coord: &Array3<f64>, mutation: &str,
               elec_atom: &Array2<f64>, vdw_atom: &Array2<f64>, 
               pb_atom: &Array2<f64>, sa_atom: &Array2<f64>,
               mm_ie: &Array1<f64>) -> SMResult {
        let mut dh: Array1<f64> = Array1::zeros(times.len());
        let mut mm: Array1<f64> = Array1::zeros(times.len());
        let mut pb: Array1<f64> = Array1::zeros(times.len());
        let mut sa: Array1<f64> = Array1::zeros(times.len());
        let mut elec: Array1<f64> = Array1::zeros(times.len());
        let mut vdw: Array1<f64> = Array1::zeros(times.len());
        for t in 0..times.len() {
            elec[t] = elec_atom.row(t).sum();
            vdw[t] = vdw_atom.row(t).sum();
            mm[t] = elec[t] + vdw[t];
            pb[t] = pb_atom.row(t).sum();
            sa[t] = sa_atom.row(t).sum();
            dh[t] = mm[t] + pb[t] + sa[t];
        }

        let mm_atom: Array2<f64> = elec_atom + vdw_atom;
        let dh_atom: Array2<f64> = &mm_atom + pb_atom + sa_atom;

        SMResult {
            mutation: mutation.to_string(),
            atom_names: atom_names.to_vec(),
            atom_res: atom_res.to_vec(),
            residues: residues.to_owned(),
            ndx_lig: ndx_lig.to_owned(),
            times: times.to_owned(),
            times_ie: times_ie.to_owned(),
            coord: coord.to_owned(),
            dh,
            mm,
            pb,
            sa,
            elec,
            vdw,
            mm_ie: mm_ie.to_owned(),
            dh_atom,
            mm_atom,
            pb_atom: pb_atom.to_owned(),
            sa_atom: sa_atom.to_owned(),
            elec_atom: elec_atom.to_owned(),
            vdw_atom: vdw_atom.to_owned(),
        }
    }

    pub fn to_bin(&self, target: &Path) {
        println!("Saving results to {}", target.to_str().unwrap());
        let mut result_as_serialize = std::fs::File::create(target).unwrap();
        serde_pickle::to_writer(&mut result_as_serialize, self, serde_pickle::SerOptions::new()).unwrap();
    }

    pub fn from(result_serialize: &str) -> SMResult {
        let result_deserialize = std::fs::File::open(result_serialize).unwrap();
        serde_pickle::from_reader(&result_deserialize, serde_pickle::DeOptions::new()).unwrap()
    }
}

pub fn analyze_controller(result_wt: &SMResult, result_as: &Vec<SMResult>, temperature: f64, sys_name: &String, wd: &Path, settings: &Settings) {
    let mut results = result_as.clone();
    results.insert(0, result_wt.clone());
    loop {
        println!("\nTime range: {} - {} ns, step = {} ns", result_wt.times[0], result_wt.times.last().unwrap(), if result_wt.times.len() > 1 {
            result_wt.times[1] - result_wt.times[0]
        } else {
            0.0
        });
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!("-1 Write residue-wised binding energy at specific time to pdb file");
        println!(" 0 Exit program");
        println!(" 1 View binding energy summary");
        println!(" 2 Output binding energy by trajectory");
        println!(" 3 Output binding energy by residue at specific time");
        println!(" 4 Output ligand binding energy by atom at specific time");
        // println!("10 Output residue-wised binding energy by time as default names");
        let sel_fun = get_input_selection();
        match sel_fun {
            Ok(-1) => {
                println!("Input the time point (in ns, e.g. 40) to output (default: average):");
                let (tmin, tmax) = get_time_range();
                let ts_ids = get_time_index(&result_wt.times, tmin, tmax);
                if ts_ids.is_empty() {
                    println!("Not valid time.");
                    continue;
                }
                println!("Writing pdb and pml file(s)...");
                for result in &results {
                    let def_name = format!("MMPBSA_binding_energy_{}.pdb", sys_name);
                    write_pdb_with_bf(result, &def_name, &ts_ids, wd, &(0..result.atom_res.len()).collect(), true);
                    let pml_name = format!("MMPBSA_binding_energy_{}.pml", sys_name);
                    let png_name = format!("MMPBSA_binding_energy_{}", sys_name);
                    write_pml(&pml_name, &def_name, &png_name, wd, settings);
                }
                println!("Finished writing pdb file(s) with binding energy information.");
                println!("Finished drawing figures with pml file(s) by PyMOL.");
            },
            Ok(0) => exit(0),
            Ok(1) => {
                println!("Input the time period (in ns, e.g. 0-40) to output (default: average):");
                let (tmin, tmax) = get_time_range();
                for result in &results {
                    analyze_summary(result, temperature, wd, &format!("{}-{}", sys_name, result.mutation), tmin, tmax)
                }
            },
            Ok(2) => {
                for result in &results {
                    analyze_traj(result, wd, &format!("{}-{}", sys_name, result.mutation))
                }
            },
            Ok(3) => {
                println!("Input the time period (in ns, e.g. 0-40) to output (default: average):");
                let (tmin, tmax) = get_time_range();
                let ts_ids = get_time_index(&result_wt.times, tmin, tmax);
                if ts_ids.is_empty() {
                    println!("Not valid time.");
                    continue;
                }
                let (range_des, target_res) = select_res_by_range(result_wt);
                println!("Writing energy file(s)...");
                for result in &results {
                    analyze_res(result, wd, &format!("{}-{}", sys_name, result.mutation), &ts_ids, &range_des, &target_res);
                }
                println!("Finished writing residue-wised binding energy file(s).");
            },
            Ok(4) => {
                for result in &results {
                    analyze_atom(result, wd, &format!("{}-{}", sys_name, result.mutation))
                }
                println!("Finished writing atom-wised binding energy pdb file(s) for ligand.");
            },
            // Ok(10) => {
            //     for result in &results {
            //         output_all_details(result, wd, &format!("{}-{}", sys_name, result.mutation))
            //     }
            // },
            Ok(_) => {},
            Err(_) => {}
        }
    }
}

fn get_time_range() -> (f64, f64) {
    println!("Note: should be time point or period. Time period should be splitted by \"-\", e.g., 3-5");
    let ts = get_input("".to_string());
    if !ts.trim().is_empty() {
        let tm: Vec<&str> = ts.split("-").collect();
        let tmin: f64 = tm.first().unwrap().parse().unwrap();
        let tmax: f64 = tm.last().unwrap().parse().unwrap();
        (tmin, tmax)
    } else {
        (0.0, f64::INFINITY)
    }
}

fn get_time_index(time_range: &Vec<f64>, tmin: f64, tmax: f64) -> Vec<usize> {
    time_range.iter().enumerate().filter_map(|(i, &t)| if t >= tmin && t <= tmax {
        Some(i)
    } else {
        None
    }).collect()
}

fn write_pml(pml_name: &String, def_name: &String, png_name: &String, wd: &Path, settings: &Settings) {
    let mut pml_file = fs::File::create(wd.join(pml_name)).unwrap();
    writeln!(pml_file, "cmd.load(\"{}\", \"complex\")", def_name).unwrap();
    writeln!(pml_file, "select protein, polymer.protein").unwrap();
    writeln!(pml_file, "preset.b_factor_putty(\"protein\", _self=cmd)").unwrap();
    writeln!(pml_file, "# cmd.spectrum(\"b\", minimum=-100, maximum=100)").unwrap();
    writeln!(pml_file, "select ligand, not polymer.protein and not solvent and not name \"NA\" and not name \"CL\"").unwrap();
    writeln!(pml_file, "cmd.spectrum(\"b\", selection=(\"ligand\"), quiet=0)").unwrap();
    writeln!(pml_file, "zoom ligand, 2.5").unwrap();
    writeln!(pml_file, "cmd.disable(\"ligand\")").unwrap();
    writeln!(pml_file, "ray 1920, 1080, async=1").unwrap();
    writeln!(pml_file, "png {}, 1920, 1080, 300, 1, 1", png_name).unwrap();
    writeln!(pml_file, "quit").unwrap();
    let result = Command::new(settings.pymol_path.as_ref().unwrap())
        .args(vec!["-cq", wd.join(pml_name).as_os_str().to_str().unwrap()])
        .stdout(Stdio::null())
        .spawn();
    match result {
        Ok(mut child) => {
            child.wait().ok();
        }
        Err(_) => {
            eprintln!("The configured PyMOL '{}' not found.", settings.pymol_path.as_ref().unwrap());
        }
    }
}

fn write_pdb_with_bf(result: &SMResult, def_name: &String, ts_ids: &Vec<usize>, wd: &Path, atom_range: &Vec<usize>, reverse: bool) {
    let mut f = fs::File::create(wd.join(def_name)).unwrap();
    let coord = &result.coord;
    writeln!(f, "REMARK  Generated by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
    writeln!(f, "REMARK  B-factor column filled with INVERSED receptor-ligand interaction energy (kJ/mol)").unwrap();
    for (id, &res_id) in result.atom_res.iter().enumerate() {
        if atom_range.contains(&id) {
            let ts_id = ts_ids.last().unwrap();
            write_atom_line(&result, id, &result.atom_names[id], res_id, ts_ids, 
                coord[[*ts_id, id, 0]], coord[[*ts_id, id, 1]], coord[[*ts_id, id, 2]], &mut f, reverse);
        }
    }
    writeln!(f, "END").unwrap();
}

fn write_atom_line(result: &SMResult, id: usize, name: &String, res_id: usize, ts_ids: &Vec<usize>, x: f64, y: f64, z: f64, f: &mut File, reverse: bool) {
    let reverse = (1 - 2 * (reverse as i32)) as f64;
    let dh_avg = result.dh_atom.select(Axis(0), &ts_ids).mean_axis(Axis(0)).unwrap();
    writeln!(f, "ATOM  {:5} {:<4} {:<3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           {:<2}", 
                id + 1, name, result.residues[res_id].name, result.residues[res_id].nr, x, y, z, 
                reverse * dh_avg[id], name.get(0..1).unwrap()).unwrap();
}

fn analyze_summary(results: &SMResult, temperature: f64, wd: &Path, sys_name: &String, tmin: f64, tmax: f64) {
    let ts_ids = get_time_index(&results.times, tmin, tmax);
    let ts_ie_ids = get_time_index(&results.times_ie, tmin, tmax);
    if ts_ids.is_empty() {
        println!("Not valid time.");
        return;
    }

    let beta_kj = 1000.0 / 8.314462618 / temperature;
    let dh_avg = results.dh.select(Axis(0), &ts_ids).mean().unwrap();
    let dh_err = results.dh.select(Axis(0), &ts_ids).std(0.0);
    let mm_avg = results.mm.select(Axis(0), &ts_ids).mean().unwrap();
    let mm_err = results.mm.select(Axis(0), &ts_ids).std(0.0);
    let elec_avg = results.elec.select(Axis(0), &ts_ids).mean().unwrap();
    let elec_err = results.elec.select(Axis(0), &ts_ids).std(0.0);
    let vdw_avg = results.vdw.select(Axis(0), &ts_ids).mean().unwrap();
    let vdw_err = results.vdw.select(Axis(0), &ts_ids).std(0.0);
    let pb_avg = results.pb.select(Axis(0), &ts_ids).mean().unwrap();
    let pb_err = results.pb.select(Axis(0), &ts_ids).std(0.0);
    let sa_avg = results.sa.select(Axis(0), &ts_ids).mean().unwrap();
    let sa_err = results.sa.select(Axis(0), &ts_ids).std(0.0);

    // Interactive Entropy
    let mm_ie_avg = results.mm_ie.select(Axis(0), &ts_ie_ids).mean().unwrap();
    let mm_sum: f64 = results.mm_ie.select(Axis(0), &ts_ie_ids).iter().map(|&mm| ((mm - mm_ie_avg) * beta_kj).exp()).sum();
    let tds = -(mm_sum / ts_ie_ids.len() as f64).ln() / beta_kj;
    let dg = dh_avg - tds;
    let ki = f64::exp(dg * beta_kj) * 1e9;    // nM

    println!("\nEnergy terms summary ({}-{} ns):", results.times[ts_ids[0]], results.times[*ts_ids.last().unwrap()]);
    println!("ΔH: {:.3} ± {:.3} kJ/mol", dh_avg, dh_err);
    println!("ΔMM: {:.3} ± {:.3} kJ/mol", mm_avg, mm_err);
    println!("ΔPB: {:.3} ± {:.3} kJ/mol", pb_avg, pb_err);
    println!("ΔSA: {:.3} ± {:.3} kJ/mol", sa_avg, sa_err);
    println!();
    println!("Δelec: {:.3} ± {:.3} kJ/mol", elec_avg, elec_err);
    println!("Δvdw: {:.3} ± {:.3} kJ/mol", vdw_avg, vdw_err);
    println!();
    println!("-TΔS: {:.3} kJ/mol", -tds);
    println!("ΔG: {:.3} kJ/mol", dg);
    println!("Ki: {:.9e} nM", ki);

    let def_name = format!("MMPBSA_{}.csv", sys_name);
    println!("Writing binding energy terms...");
    let mut energy_sum = fs::File::create(wd.join(&def_name)).unwrap();
    write!(energy_sum, "Energy Term,value,std.P,info ({}-{} ns)\n", results.times[ts_ids[0]], results.times[*ts_ids.last().unwrap()]).unwrap();
    write!(energy_sum, "ΔH,{:.3},{:.3},ΔH=ΔMM+ΔPB+ΔSA (kJ/mol)\n", dh_avg, dh_err).unwrap();
    write!(energy_sum, "ΔMM,{:.3},{:.3},ΔMM=Δelec+ΔvdW (kJ/mol)\n", mm_avg, mm_err).unwrap();
    write!(energy_sum, "ΔPB,{:.3},{:.3},(kJ/mol)\n", pb_avg, pb_err).unwrap();
    write!(energy_sum, "ΔSA,{:.3},{:.3},(kJ/mol)\n", sa_avg, sa_err).unwrap();
    write!(energy_sum, "\n").unwrap();
    write!(energy_sum, "Δelec,{:.3},{:.3},(kJ/mol)\n", elec_avg, elec_err).unwrap();
    write!(energy_sum, "ΔvdW,{:.3},{:.3},(kJ/mol)\n", vdw_avg, vdw_err).unwrap();
    write!(energy_sum, "\n").unwrap();
    write!(energy_sum, "-TΔS,{:.3},,(kJ/mol)\n", -tds).unwrap();
    write!(energy_sum, "ΔG,{:.3},,ΔG=ΔH-TΔS (kJ/mol)\n", dg).unwrap();
    write!(energy_sum, "Ki,{:.9e},,Ki=exp(ΔG/RT) (nM)\n", ki).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_traj(results: &SMResult, wd: &Path, sys_name: &String) {
    println!("Writing binding energy terms...");
    let def_name = format!("MMPBSA_{}_traj.csv", sys_name);
    let mut energy_sum = fs::File::create(wd.join(&def_name)).unwrap();
    write!(energy_sum, "Time (ns),ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW\n").unwrap();
    for i in 0..results.times.len() {
        write!(energy_sum, "{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n",
                            results.times[i], results.dh[i],
                            results.mm[i], results.pb[i], results.sa[i],
                            results.elec[i], results.vdw[i]).unwrap();
    }
    println!("Binding energy terms writen to {}", &def_name);

    // ΔH curve
    println!("Plotting binding energy figures...");
    let mut curve = Curve::new();
    curve.set_line_width(2.0);
    curve.draw(&results.times, &results.dh.to_vec());
    let mut plot = Plot::new();
    if cfg!(windows) {
        plot.set_python_exe("python");
    }
    let def_name = format!("MMPBSA_{}_ΔH_traj.png", sys_name);
    match plot.add(&curve)
        .grid_and_labels("Time (ns)", "Binding Energy (kJ/mol)")
        .set_label_x_fontsize(18.0)
        .set_label_y_fontsize(18.0)
        .set_ticks_x_fontsize(14.0)
        .set_ticks_y_fontsize(14.0)
        .save(&wd.join(&def_name)).ok() {
            Some(_) => println!("Figure drawn to {}", &def_name),
            None => println!("Not drawn due to the matplotlib error.")
        };
    println!("Binding energy terms writen to {}", &def_name);
}

fn select_res_by_range(results: &SMResult) -> (String, Vec<usize>) {
    println!("Determine the residue range to output:");
    println!(" 1 Ligand and receptor residues by: CA within 4 A");
    println!(" 2 Ligand and receptor residues by: CA within 6 A");
    println!(" 3 Ligand and receptor residues by: CA within 8 A");
    println!(" 4 Ligand and receptor residues by: CA within a specified distance");
    println!(" 5 Self-defined residue range");
    // 残基范围确定
    let i = get_input_selection().unwrap();
    let mut range_des = String::from("4A");
    let target_res = match i {
        1 => {
            get_residue_range_from_results(results, 4.0)
        },
        2 => {
            range_des = String::from("6A");
            get_residue_range_from_results(results, 6.0)
        },
        3 => {
            range_des = String::from("8A");
            get_residue_range_from_results(results, 8.0)
        },
        4 => {
            println!("Input the cut-off distance you want to expand from ligand, default: 4");
            let cutoff = get_input(4.0);
            range_des = format!("{:.1}A", cutoff);
            get_residue_range_from_results(results, cutoff)
        },
        5 => {
            let mut res_range = String::new();
            loop {
                println!("Input the residue range you want to output (e.g., 1-3, 5), default: all");
                println!("Input \"?\" to view the residues list");
                res_range = get_input(res_range);
                if res_range.eq("?") {
                    results.residues.iter().enumerate().for_each(|(i, r)| {
                        print!("{}{}, ", r.nr, r.name);
                        if (i + 1) % 5 == 0 {
                            println!();
                        }
                    });
                    println!();
                    res_range.clear();
                } else {
                    break;
                }
            }
            range_des = res_range.to_string();
            let res_range: Vec<i32> = match res_range.len() {
                0 => {
                    range_des = "all".to_string();
                    results.residues.iter().map(|r| r.nr).collect()
                },
                _ => range2list(&res_range)
            };
            results.atom_res
                .iter()
                .filter(|&&i| res_range.contains(&(results.residues[i].nr)))    // 用户筛选用nr
                .map(|&i| i)     // 索引用nr
                .collect()
        },
        _ => vec![],
    };
    (range_des, target_res)
}

fn analyze_res(results: &SMResult, wd: &Path, sys_name: &String, ts_ids: &Vec<usize>, range_des: &String, target_res: &Vec<usize>) {
    let def_name = format!("MMPBSA_{}_res_{}.csv", sys_name, range_des);
    let (tar_res_nr, tar_res_name, 
        tar_res_energy, tar_res_energy_err) = get_target_res_data(results, ts_ids, target_res);
    write_res_csv(&tar_res_nr, &tar_res_name, &tar_res_energy, &tar_res_energy_err, wd, &def_name);
    plot_res_csv(&tar_res_nr, &tar_res_name, &tar_res_energy, &tar_res_energy_err, wd, &format!("MMPBSA_{}_res_{}.png", sys_name, range_des));
}

fn analyze_atom(results: &SMResult, wd: &Path, sys_name: &String) {
    let def_name = format!("MMPBSA_{}_ligand.pdb", sys_name);
    write_pdb_with_bf(results, &def_name, &vec![0], wd, &results.ndx_lig, false);
}

fn get_target_res_data(results: &SMResult, ts_ids: &Vec<usize>, target_res: &Vec<usize>) -> (Vec<i32>, Vec<String>, [Vec<f64>; 6], [Vec<f64>; 6]) {
    let res_nr: Vec<i32> = results.residues.iter().enumerate().filter_map(|(id, res)| if target_res.contains(&id) {
        Some(res.nr)
    } else {
        None
    }).collect();
    let res_name: Vec<String> = results.residues.iter().enumerate().filter_map(|(id, res)| if target_res.contains(&id) {
        Some(res.name.to_string())
    } else {
        None
    }).collect();
    let mut dh_res = vec![];
    let mut dh_res_err = vec![];
    let mut mm_res = vec![];
    let mut mm_res_err = vec![];
    let mut pb_res = vec![];
    let mut pb_res_err = vec![];
    let mut sa_res = vec![];
    let mut sa_res_err = vec![];
    let mut elec_res = vec![];
    let mut elec_res_err = vec![];
    let mut vdw_res = vec![];
    let mut vdw_res_err = vec![];
    // 该残基的所有原子所在索引
    let get_cur_res_atom_ids = |target_resid| results
        .atom_res
        .iter()
        .enumerate()
        .filter_map(|(i, &a)| if a == target_resid {
        Some(i)
    } else {
        None
    }).collect::<Vec<usize>>();
    for (id, _res) in results.residues.iter().enumerate() {
        if !target_res.contains(&id) {
            continue;
        }
        let atom_ids = get_cur_res_atom_ids(id);
        let atom_energy = |arr: &Array2<f64>| {
            arr.select(Axis(1), &atom_ids).select(Axis(0), &ts_ids)
        };
        let dh_resi: Array1<f64> = atom_energy(&results.dh_atom).sum_axis(Axis(1));
        dh_res.push(dh_resi.mean().unwrap());
        dh_res_err.push(dh_resi.std(0.0));
        let mm_resi: Array1<f64> = atom_energy(&results.mm_atom).sum_axis(Axis(1));
        mm_res.push(mm_resi.mean().unwrap());
        mm_res_err.push(mm_resi.std(0.0));
        let pb_resi: Array1<f64> = atom_energy(&results.pb_atom).sum_axis(Axis(1));
        pb_res.push(pb_resi.mean().unwrap());
        pb_res_err.push(pb_resi.std(0.0));
        let sa_resi: Array1<f64> = atom_energy(&results.sa_atom).sum_axis(Axis(1));
        sa_res.push(sa_resi.mean().unwrap());
        sa_res_err.push(sa_resi.std(0.0));
        let elec_resi: Array1<f64> = atom_energy(&results.elec_atom).sum_axis(Axis(1));
        elec_res.push(elec_resi.mean().unwrap());
        elec_res_err.push(elec_resi.std(0.0));
        let vdw_resi: Array1<f64> = atom_energy(&results.vdw_atom).sum_axis(Axis(1));
        vdw_res.push(vdw_resi.mean().unwrap());
        vdw_res_err.push(vdw_resi.std(0.0));
    }
    (res_nr, res_name, [dh_res, mm_res, pb_res, sa_res, elec_res, vdw_res], 
        [dh_res_err, mm_res_err, pb_res_err, sa_res_err, elec_res_err, vdw_res_err])
}

fn write_res_csv(tar_res_nr: &Vec<i32>, tar_res_name: &Vec<String>, 
                 tar_res_energy: &[Vec<f64>; 6], tar_res_energy_err: &[Vec<f64>; 6], 
                 wd: &Path, def_name: &String) {
    let mut res_energy_file = fs::File::create(wd.join(def_name)).unwrap();
    res_energy_file.write_all("id,name,ΔH,ΔH_std.P,ΔMM,ΔMM_std.P,ΔPB,ΔPB_std.P,ΔSA,ΔSA_std.P,Δelec,Δelec_std.P,ΔvdW,ΔvdW_std.P\n".as_bytes()).unwrap();
    for tar_res_id in 0..tar_res_energy[0].len() {
        write!(res_energy_file, "{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n", 
            tar_res_nr[tar_res_id], tar_res_name[tar_res_id],
            tar_res_energy[0][tar_res_id], tar_res_energy_err[0][tar_res_id],
            tar_res_energy[1][tar_res_id], tar_res_energy_err[1][tar_res_id],
            tar_res_energy[2][tar_res_id], tar_res_energy_err[2][tar_res_id],
            tar_res_energy[3][tar_res_id], tar_res_energy_err[3][tar_res_id],
            tar_res_energy[4][tar_res_id], tar_res_energy_err[4][tar_res_id],
            tar_res_energy[5][tar_res_id], tar_res_energy_err[5][tar_res_id])
            .expect("Error while writing residue-wised energy file");
    }
}

fn plot_res_csv(tar_res_nr: &Vec<i32>, tar_res_name: &Vec<String>, 
                tar_res_energy: &[Vec<f64>; 6], tar_res_energy_err: &[Vec<f64>; 6], 
                wd: &Path, def_name: &String) {
    println!("Plotting residue-wised binding energy figures...");
    let mut bar = Barplot::new();
    let mut plot = Plot::new();
    if cfg!(windows) {
        plot.set_python_exe("python");
    }
    let ytick_labels: Vec<String> = tar_res_nr.iter().enumerate().map(|(i, r)| format!("{}{}", tar_res_name[i], r)).collect();
    let ytick_labels: Vec<&str> = ytick_labels.iter().map(|s| s.as_str()).collect();
    bar.set_errors(&tar_res_energy_err[0])
        .set_horizontal(true)
        .draw_with_str(&ytick_labels, &tar_res_energy[0]); // requires string, as expected
    match plot.add(&bar)
            .set_figure_size_inches(6.4, tar_res_nr.len() as f64 * 0.48)
            .set_ymin(-0.5)
            .set_ymax(tar_res_nr.len() as f64 - 0.5)
            .grid_and_labels("Binding Energy (kJ/mol)", "Residue")
            .set_label_x_fontsize(18.0)
            .set_label_y_fontsize(18.0)
            .set_ticks_x_fontsize(14.0)
            .set_ticks_y_fontsize(14.0)
            .save(&wd.join(&def_name)).ok() {
        Some(_) => println!("Figure drawn to {}", &def_name),
        None => println!("Not drawn due to the matplotlib error.")
    }
}

fn get_residue_range_from_results(results: &SMResult, cutoff: f64) -> Vec<usize> {
    let last_frame = results.times.len() - 1;
    get_residue_range_ca(&results.coord.slice(s![last_frame, .., ..]).to_owned(), 
        &results.ndx_lig, cutoff, &results.atom_res, &results.atom_names, &results.residues)
}


// fn analyze_dh_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
//     println!("Writing binding energy terms...");
//     let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
//     energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//     for res in &results.residues {
//         energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
//     }
//     for i in 0..results.times.len() {
//         energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
//         for dh in &results.dh_atom.row(i) {
//             energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
//         }
//     }
//     energy_res.write_all("\n".as_bytes()).unwrap();
//     println!("Binding energy terms have been writen to {}", &def_name);
// }

// fn analyze_mm_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
//     println!("Writing binding energy terms...");
//     let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
//     energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//     for res in &results.residues {
//         energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
//     }
//     for i in 0..results.times.len() {
//         energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
//         for mm in &results.mm_atom.row(i) {
//             energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
//         }
//     }
//     energy_res.write_all("\n".as_bytes()).unwrap();
//     println!("Binding energy terms have been writen to {}", &def_name);
// }

// fn analyze_pb_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
//     println!("Writing binding energy terms...");
//     let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
//     energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//     for res in &results.residues {
//         energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
//     }
//     for i in 0..results.times.len() {
//         energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
//         for pb in &results.pb_atom.row(i) {
//             energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
//         }
//     }
//     energy_res.write_all("\n".as_bytes()).unwrap();
//     println!("Binding energy terms have been writen to {}", &def_name);
// }

// fn analyze_sa_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
//     println!("Writing binding energy terms...");
//     let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
//     energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//     for res in &results.residues {
//         energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
//     }
//     for i in 0..results.times.len() {
//         energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
//         for sa in &results.sa_atom.row(i) {
//             energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
//         }
//     }
//     energy_res.write_all("\n".as_bytes()).unwrap();
//     println!("Binding energy terms have been writen to {}", &def_name);
// }

// fn analyze_elec_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
//     println!("Writing binding energy terms...");
//     let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
//     energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//     for res in &results.residues {
//         energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
//     }
//     for i in 0..results.times.len() {
//         energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
//         for elec in &results.elec_atom.row(i) {
//             energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
//         }
//     }
//     energy_res.write_all("\n".as_bytes()).unwrap();
//     println!("Binding energy terms have been writen to {}", &def_name);
// }

// fn analyze_vdw_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
//     println!("Writing binding energy terms...");
//     let mut energy_res = fs::File::create(wd.join(def_name)).unwrap();
//     energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//     for res in &results.residues {
//         energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
//     }
//     for i in 0..results.times.len() {
//         energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
//         for vdw in &results.vdw_atom.row(i) {
//             energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
//         }
//     }
//     energy_res.write_all("\n".as_bytes()).unwrap();
//     println!("Binding energy terms have been writen to {}", def_name);
// }

// pub fn output_all_details(results: &SMResult, wd: &Path, sys_name: &String) {
//     analyze_dh_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔH.csv", sys_name));
//     analyze_mm_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔMM.csv", sys_name));
//     analyze_pb_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔPB.csv", sys_name));
//     analyze_sa_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔSA.csv", sys_name));
//     analyze_elec_res_traj(results, wd, &format!("MMPBSA_{}_res_Δelec.csv", sys_name));
//     analyze_vdw_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔvdW.csv", sys_name));
// }
