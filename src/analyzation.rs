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
    pub coord: Array3<f64>,
    pub dh: Array1<f64>,
    pub mm: Array1<f64>,
    pub pb: Array1<f64>,
    pub sa: Array1<f64>,
    pub elec: Array1<f64>,
    pub vdw: Array1<f64>,
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
               times: &Vec<f64>, coord: &Array3<f64>, mutation: &str,
               elec_atom: &Array2<f64>, vdw_atom: &Array2<f64>, 
               pb_atom: &Array2<f64>, sa_atom: &Array2<f64>) -> SMResult {
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
            coord: coord.to_owned(),
            dh,
            mm,
            pb,
            sa,
            elec,
            vdw,
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
    println!("\nTime range: {} - {} ns, step = {} ns", result_wt.times[0], result_wt.times.last().unwrap(), if result_wt.times.len() > 1 {
        result_wt.times[1] - result_wt.times[0]
    } else {
        0.0
    });
    let mut results = result_as.clone();
    results.insert(0, result_wt.clone());
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!("-1 Write residue-wised binding energy at specific time to pdb file");
        println!(" 0 Exit program");
        println!(" 1 View binding energy summary");
        println!(" 2 Output binding energy by trajectory");
        println!(" 3 Output binding energy by residue at specific time");
        println!(" 4 Output ligand binding energy by atom at specific time");
        println!("10 Output residue-wised binding energy by time as default names");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            -1 => {
                println!("Input the time point (in ns) to output (default: average):");
                let ts_id = get_time_points(result_wt);
                println!("Writing pdb and pml file(s)...");
                for result in &results {
                    if let Some(ts_id) = ts_id {
                        let def_name = format!("MMPBSA_binding_energy_{}_{}ns.pdb", sys_name, result.times[ts_id]);
                        write_pdb_with_bf(result, &def_name, ts_id, wd, &(0..result.atom_names.len()).collect(), true, true);
                        let pml_name = format!("MMPBSA_binding_energy_{}_{}ns.pml", sys_name, result.times[ts_id]);
                        let png_name = format!("MMPBSA_binding_energy_{}_{}ns", sys_name, result.times[ts_id]);
                        write_pml(&pml_name, &def_name, &png_name, wd, settings);
                    } else {
                        let def_name = format!("MMPBSA_binding_energy_{}_avg.pdb", sys_name);
                        write_pdb_with_bf(result, &def_name, 0, wd, &(0..result.atom_names.len()).collect(), false, true);
                        let pml_name = format!("MMPBSA_binding_energy_{}_avg.pml", sys_name);
                        let png_name = format!("MMPBSA_binding_energy_{}_avg", sys_name);
                        write_pml(&pml_name, &def_name, &png_name, wd, settings);
                    }
                }
                println!("Finished writing pdb file(s) with binding energy information.");
                println!("The pml file(s) could be loaded by PyMOL.");
            },
            0 => exit(0),
            1 => {
                for result in &results {
                    analyze_summary(result, temperature, wd, &format!("{}-{}", sys_name, result.mutation))
                }
            },
            2 => {
                for result in &results {
                    analyze_traj(result, wd, &format!("{}-{}", sys_name, result.mutation))
                }
            },
            3 => {
                println!("Input the time point (in ns) to output (default: average):");
                let ts_id = get_time_points(result_wt);
                let (range_des, target_res) = select_res_by_range(result_wt);
                println!("Writing energy file(s)...");
                for result in &results {
                    analyze_res(result, wd, &format!("{}-{}", sys_name, result.mutation), ts_id, &range_des, &target_res);
                }
                println!("Finished writing residue-wised binding energy file(s).");
            },
            4 => {
                for result in &results {
                    analyze_atom(result, wd, &format!("{}-{}", sys_name, result.mutation))
                }
                println!("Finished writing atom-wised binding energy pdb file(s) for ligand.");
            },
            10 => {
                for result in &results {
                    output_all_details(result, wd, &format!("{}-{}", sys_name, result.mutation))
                }
            },
            _ => println!("Invalid input")
        }
    }
}

fn get_time_points(result: &SMResult) -> Option<usize> {
    let ts = get_input(-1.0);
    if ts != -1.0 {
        get_time_index(ts, result)
    } else {
        None
    }
}

fn get_time_index(ts: f64, results: &SMResult) -> Option<usize> {
    if ts == results.times[0] {
        Some(0)
    } else {
        if results.times.len() > 1 {
            Some(((ts - results.times[0]) / (results.times[1] - results.times[0])) as usize)
        } else {
            None
        }
    }
}

fn write_pml(pml_name: &String, def_name: &String, png_name: &String, wd: &Path, settings: &Settings) {
    let mut pml_file = fs::File::create(wd.join(pml_name)).unwrap();
    writeln!(pml_file, "cmd.load(\"{}\", \"complex\")", def_name).unwrap();
    writeln!(pml_file, "select protein, polymer.protein").unwrap();
    writeln!(pml_file, "preset.b_factor_putty(\"protein\", _self=cmd)").unwrap();
    writeln!(pml_file, "select ligand, not polymer.protein and not solvent and not name \"NA\" and not name \"CL\"").unwrap();
    writeln!(pml_file, "cmd.spectrum(\"b\", selection=(\"ligand\"), quiet=0)").unwrap();
    writeln!(pml_file, "zoom ligand, 2.5").unwrap();
    writeln!(pml_file, "cmd.disable(\"ligand\")").unwrap();
    writeln!(pml_file, "ray 1920, 1080, async=1").unwrap();
    writeln!(pml_file, "png {}", png_name).unwrap();
    writeln!(pml_file, "quit").unwrap();
    let result = Command::new(settings.pymol_path.as_ref().unwrap())
        .args(wd.join(pml_name).as_os_str().to_str())
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

fn write_pdb_with_bf(result: &SMResult, def_name: &String, ts_id: usize, wd: &Path, atom_range: &Vec<usize>, by_frame: bool, reverse: bool) {
    let mut f = fs::File::create(wd.join(def_name)).unwrap();
    let coord = &result.coord;
    writeln!(f, "REMARK  Generated by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)").unwrap();
    writeln!(f, "REMARK  B-factor column filled with INVERSED receptor-ligand interaction energy (kJ/mol)").unwrap();
    for (id, &res_id) in result.atom_res.iter().enumerate() {
        if atom_range.contains(&id) {
            write_atom_line(&result, id, &result.atom_names[id], res_id, ts_id, 
                coord[[ts_id, id, 0]], coord[[ts_id, id, 1]], coord[[ts_id, id, 2]], &mut f, by_frame, reverse);
        }
    }
    writeln!(f, "END").unwrap();
}

fn write_atom_line(result: &SMResult, id: usize, name: &String, res_id: usize, ts_id: usize, x: f64, y: f64, z: f64, f: &mut File, by_frame: bool, reverse: bool) {
    let reverse = (1 - 2 * (reverse as i32)) as f64;
    if by_frame {
        writeln!(f, "ATOM  {:5} {:<4} {:<3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           {:<2}", 
                    id + 1, name, result.residues[res_id].name, result.residues[res_id].nr, x, y, z, 
                    reverse * result.dh_atom[[ts_id, id]], name.get(0..1).unwrap()).unwrap();
    } else {
        let dh_avg = result.dh_atom.mean_axis(Axis(0)).unwrap();
        writeln!(f, "ATOM  {:5} {:<4} {:<3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           {:<2}", 
                    id + 1, name, result.residues[res_id].name, result.residues[res_id].nr, x, y, z, 
                    reverse * dh_avg[id], name.get(0..1).unwrap()).unwrap();
    }
}

fn analyze_summary(results: &SMResult, temperature: f64, wd: &Path, sys_name: &String) {
    let beta_kj = 1000.0 / 8.314462618 / temperature;
    let dh_avg = results.dh.mean().unwrap();
    let mm_avg = results.mm.mean().unwrap();
    let elec_avg = results.elec.mean().unwrap();
    let vdw_avg = results.vdw.mean().unwrap();
    let pb_avg = results.pb.mean().unwrap();
    let sa_avg = results.sa.mean().unwrap();

    let mm_sum: f64 = results.mm.iter().map(|&mm| f64::exp((mm - mm_avg) * beta_kj)).sum();
    let tds = -(mm_sum / results.mm.len() as f64).ln() / beta_kj;
    let dg = dh_avg - tds;
    let ki = f64::exp(dg * beta_kj) * 1e9;    // nM

    println!("\nEnergy terms summary:");
    println!("ΔH: {:.3} kJ/mol", dh_avg);
    println!("ΔMM: {:.3} kJ/mol", mm_avg);
    println!("ΔPB: {:.3} kJ/mol", pb_avg);
    println!("ΔSA: {:.3} kJ/mol", sa_avg);
    println!();
    println!("Δelec: {:.3} kJ/mol", elec_avg);
    println!("Δvdw: {:.3} kJ/mol", vdw_avg);
    println!();
    println!("TΔS: {:.3} kJ/mol", tds);
    println!("ΔG: {:.3} kJ/mol", dg);
    println!("Ki: {:.9e} nM", ki);

    let def_name = format!("MMPBSA_{}.csv", sys_name);
    println!("Writing binding energy terms...");
    let mut energy_sum = fs::File::create(wd.join(&def_name)).unwrap();
    write!(energy_sum, "Energy Term,value,info\n").unwrap();
    write!(energy_sum, "ΔH,{:.3},ΔH=ΔMM+ΔPB+ΔSA (kJ/mol)\n", dh_avg).unwrap();
    write!(energy_sum, "ΔMM,{:.3},ΔMM=Δelec+ΔvdW (kJ/mol)\n", mm_avg).unwrap();
    write!(energy_sum, "ΔPB,{:.3},(kJ/mol)\n", pb_avg).unwrap();
    write!(energy_sum, "ΔSA,{:.3},(kJ/mol)\n", sa_avg).unwrap();
    write!(energy_sum, "\n").unwrap();
    write!(energy_sum, "Δelec,{:.3},(kJ/mol)\n", elec_avg).unwrap();
    write!(energy_sum, "ΔvdW,{:.3},(kJ/mol)\n", vdw_avg).unwrap();
    write!(energy_sum, "\n").unwrap();
    write!(energy_sum, "TΔS,{:.3},(kJ/mol)\n", tds).unwrap();
    write!(energy_sum, "ΔG,{:.3},ΔG=ΔH-TΔS (kJ/mol)\n", dg).unwrap();
    write!(energy_sum, "Ki,{:.9e},Ki=exp(ΔG/RT) (nM)\n", ki).unwrap();
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
    plot.add(&curve)
        .grid_and_labels("Time (ns)", "Binding Energy (kJ/mol)")
        .set_label_x_fontsize(18.0)
        .set_label_y_fontsize(18.0)
        .set_ticks_x_fontsize(14.0)
        .set_ticks_y_fontsize(14.0)
        .save(&wd.join(&def_name)).unwrap();
    println!("Binding energy terms writen to {}", &def_name);
}

fn select_res_by_range(result_wt: &SMResult) -> (String, Vec<usize>) {
    println!("Determine the residue range to output:");
    println!(" 1 Ligand and receptor residues by: CA within 4 A");
    println!(" 2 Ligand and receptor residues by: CA within 6 A");
    println!(" 3 Ligand and receptor residues by: CA within 8 A");
    println!(" 4 Ligand and receptor residues by: CA within a specified distance");
    println!(" 5 Self-defined residue range");
    // 残基范围确定
    let i: i32 = get_input_selection();
    let mut range_des = String::from("4A");
    let target_res = match i {
        1 => {
            get_residue_range_from_results(result_wt, 4.0)
        },
        2 => {
            range_des = String::from("6A");
            get_residue_range_from_results(result_wt, 6.0)
        },
        3 => {
            range_des = String::from("8A");
            get_residue_range_from_results(result_wt, 8.0)
        },
        4 => {
            println!("Input the cut-off distance you want to expand from ligand, default: 4");
            let cutoff = get_input(4.0);
            range_des = format!("{:.1}A", cutoff);
            get_residue_range_from_results(result_wt, cutoff)
        },
        5 => {
            println!("Input the residue range you want to output (e.g., 1-3, 5), default: all");
            let res_range = get_input(String::new());
            range_des = res_range.to_string();
            let res_range: Vec<i32> = match res_range.len() {
                0 => {
                    range_des = "all".to_string();
                    result_wt.residues.iter().map(|r| r.nr).collect()
                },
                _ => range2list(&res_range)
            };
            result_wt.atom_res
                .iter()
                .filter(|&&i| res_range.contains(&(result_wt.residues[i].nr)))    // 用户筛选用nr
                .map(|&i| result_wt.residues[i].id)     // 索引用id
                .collect()
        },
        _ => {
            println!("Invalid selection");
            vec![]
        }
    };
    (range_des, target_res)
}

fn analyze_res(results: &SMResult, wd: &Path, sys_name: &String, ts_id: Option<usize>, range_des: &String, target_res: &Vec<usize>) {
    if let Some(ts_id) = ts_id {
        let def_name = format!("MMPBSA_{}_res_{}_{}ns.csv", sys_name, range_des, results.times[ts_id]);
        let (tar_res_nr, tar_res_name, tar_res_energy) = get_target_res_data(results, ts_id, target_res);
        write_res_csv(&tar_res_nr, &tar_res_name, &tar_res_energy, wd, &def_name);
        plot_res_csv(&tar_res_nr, &tar_res_name, &tar_res_energy, wd, &format!("MMPBSA_{}_res_{}_{}ns.png", sys_name, range_des, results.times[ts_id]));
    } else {
        let def_name = format!("MMPBSA_{}_res_{}_avg.csv", sys_name, range_des);
        let (tar_res_nr, tar_res_name, tar_res_energy) = get_target_res_avg_data(results, target_res);
        write_res_csv(&tar_res_nr, &tar_res_name, &tar_res_energy, wd, &def_name);
        plot_res_csv(&tar_res_nr, &tar_res_name, &tar_res_energy, wd, &format!("MMPBSA_{}_res_{}_avg.png", sys_name, range_des));
    }
}

fn analyze_atom(results: &SMResult, wd: &Path, sys_name: &String) {
    let def_name = format!("MMPBSA_{}_ligand.pdb", sys_name);
    write_pdb_with_bf(results, &def_name, 0, wd, &results.ndx_lig, false, false);
}

fn get_target_res_data(results: &SMResult, ts_id: usize, target_res: &Vec<usize>) -> (Vec<i32>, Vec<String>, [Vec<f64>; 6]) {
    let res_nr: Vec<i32> = results.residues.iter().filter_map(|res| if target_res.contains(&res.id) {
        Some(res.nr)
    } else {
        None
    }).collect();
    let res_name: Vec<String> = results.residues.iter().filter_map(|res| if target_res.contains(&res.id) {
        Some(res.name.to_string())
    } else {
        None
    }).collect();
    let mut dh_res = vec![];
    let mut mm_res = vec![];
    let mut pb_res = vec![];
    let mut sa_res = vec![];
    let mut elec_res = vec![];
    let mut vdw_res = vec![];
    for res in results.residues.iter() {
        if !target_res.contains(&res.id) {
            continue;
        }
        dh_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.dh_atom[[ts_id, a]])
        } else {
            None
        } ).sum::<f64>());
        mm_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.mm_atom[[ts_id, a]])
        } else {
            None
        } ).sum::<f64>());
        pb_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.pb_atom[[ts_id, a]])
        } else {
            None
        } ).sum::<f64>());
        sa_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.sa_atom[[ts_id, a]])
        } else {
            None
        } ).sum::<f64>());
        elec_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.elec_atom[[ts_id, a]])
        } else {
            None
        } ).sum::<f64>());
        vdw_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.vdw_atom[[ts_id, a]])
        } else {
            None
        } ).sum::<f64>());
    }
    (res_nr, res_name, [dh_res, mm_res, pb_res, sa_res, elec_res, vdw_res])
}

fn get_target_res_avg_data(results: &SMResult, target_res: &Vec<usize>) -> (Vec<i32>, Vec<String>, [Vec<f64>; 6]) {
    let res_nr: Vec<i32> = results.residues.iter().filter_map(|res| if target_res.contains(&res.id) {
        Some(res.nr)
    } else {
        None
    }).collect();
    let res_name: Vec<String> = results.residues.iter().filter_map(|res| if target_res.contains(&res.id) {
        Some(res.name.to_string())
    } else {
        None
    }).collect();
    let mut dh_res = vec![];
    let mut mm_res = vec![];
    let mut pb_res = vec![];
    let mut sa_res = vec![];
    let mut elec_res = vec![];
    let mut vdw_res = vec![];
    for res in results.residues.iter() {
        if !target_res.contains(&res.id) {
            continue;
        }
        dh_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.dh_atom.column(a).sum())
        } else {
            None
        } ).sum::<f64>() / results.times.len() as f64);
        mm_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.mm_atom.column(a).sum())
        } else {
            None
        } ).sum::<f64>() / results.times.len() as f64);
        pb_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.pb_atom.column(a).sum())
        } else {
            None
        } ).sum::<f64>() / results.times.len() as f64);
        sa_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.sa_atom.column(a).sum())
        } else {
            None
        } ).sum::<f64>() / results.times.len() as f64);
        elec_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.elec_atom.column(a).sum())
        } else {
            None
        } ).sum::<f64>() / results.times.len() as f64);
        vdw_res.push(results.atom_res.iter().filter_map(|&a| if a == res.id {
            Some(results.vdw_atom.column(a).sum())
        } else {
            None
        } ).sum::<f64>() / results.times.len() as f64);
    }
    (res_nr, res_name, [dh_res, mm_res, pb_res, sa_res, elec_res, vdw_res])
}

fn write_res_csv(tar_res_nr: &Vec<i32>, tar_res_name: &Vec<String>, tar_res_energy: &[Vec<f64>; 6], wd: &Path, def_name: &String) {
    let mut res_energy_file = fs::File::create(wd.join(def_name)).unwrap();
    res_energy_file.write_all("id,name,ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW\n".as_bytes()).unwrap();
    for tar_res_id in 0..tar_res_energy[0].len() {
        write!(res_energy_file, "{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n", 
            tar_res_nr[tar_res_id], tar_res_name[tar_res_id],
            tar_res_energy[0][tar_res_id],
            tar_res_energy[1][tar_res_id],
            tar_res_energy[2][tar_res_id],
            tar_res_energy[3][tar_res_id],
            tar_res_energy[4][tar_res_id],
            tar_res_energy[5][tar_res_id])
            .expect("Error while writing residue-wised energy file");
    }
}

fn plot_res_csv(tar_res_nr: &Vec<i32>, tar_res_name: &Vec<String>, tar_res_energy: &[Vec<f64>; 6], wd: &Path, def_name: &String) {
    println!("Plotting residue-wised binding energy figures...");
    let mut bar = Barplot::new();
    bar.draw(&(0..tar_res_nr.len()).map(|a| a as f64).collect(), &tar_res_energy[0]);
    let mut plot = Plot::new();
    if cfg!(windows) {
        plot.set_python_exe("python");
    }
    let xticks: Vec<usize> = (0..tar_res_nr.len()).collect();
    let xtick_labels: Vec<String> = tar_res_nr.iter().enumerate().map(|(i, r)| format!("{}{}", tar_res_name[i], r)).collect();
    match plot.add(&bar)
            .set_figure_size_inches(tar_res_nr.len() as f64 * 0.64, 4.8)
            .set_ticks_x_labels(&xticks, &xtick_labels)
            .set_rotation_ticks_x(45.0)
            .grid_and_labels("Residue", "Binding Energy (kJ/mol)")
            .set_label_x_fontsize(18.0)
            .set_label_y_fontsize(18.0)
            .set_ticks_x_fontsize(14.0)
            .set_ticks_y_fontsize(14.0)
            .save(&wd.join(&def_name)).ok() {
        Some(_) => println!("Residue-wised binding energy terms writen to {}", &def_name),
        None => println!("Not drawn due to the matplotlib error.")
    }
}

fn get_residue_range_from_results(results: &SMResult, cutoff: f64) -> Vec<usize> {
    let last_frame = results.times.len() - 1;
    get_residue_range_ca(&results.coord.slice(s![last_frame, .., ..]).to_owned(), 
        &results.ndx_lig, cutoff, &results.atom_res, &results.atom_names, &results.residues)
}


fn analyze_dh_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for dh in &results.dh_atom.row(i) {
            energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_mm_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for mm in &results.mm_atom.row(i) {
            energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_pb_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for pb in &results.pb_atom.row(i) {
            energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_sa_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for sa in &results.sa_atom.row(i) {
            energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_elec_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for elec in &results.elec_atom.row(i) {
            energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_vdw_res_traj(results: &SMResult, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for vdw in &results.vdw_atom.row(i) {
            energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", def_name);
}

pub fn output_all_details(results: &SMResult, wd: &Path, sys_name: &String) {
    analyze_dh_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔH.csv", sys_name));
    analyze_mm_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔMM.csv", sys_name));
    analyze_pb_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔPB.csv", sys_name));
    analyze_sa_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔSA.csv", sys_name));
    analyze_elec_res_traj(results, wd, &format!("MMPBSA_{}_res_Δelec.csv", sys_name));
    analyze_vdw_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔvdW.csv", sys_name));
}
