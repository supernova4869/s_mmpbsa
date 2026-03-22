use std::collections::BTreeSet;
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::process::{exit, Command, Stdio};
use ndarray::{Array1, Array2, Array3, Axis};
use serde::{Deserialize, Serialize};
use plotters::prelude::*;
use crate::parse_tpr::Residue;
use crate::settings::Settings;
use crate::utils::{get_input, get_input_selection, range2list};

#[derive(Clone, Serialize, Deserialize)]
pub struct SMResults {
    pub sm_results: Vec<SMResult>
}

impl SMResults {
    pub fn new(sm_results: Vec<SMResult>) -> SMResults {
        SMResults { sm_results }
    }

    pub fn to_bin(&self, target: &Path) {
        println!("Saving results to {}", target.to_str().unwrap());
        let mut result_as_serialize = std::fs::File::create(target).unwrap();
        serde_pickle::to_writer(&mut result_as_serialize, self, serde_pickle::SerOptions::new()).unwrap();
    }

    pub fn from(result_serialize: &str) -> Result<SMResults, serde_pickle::Error> {
        let result_deserialize = std::fs::File::open(result_serialize).unwrap();
        serde_pickle::from_reader(&result_deserialize, serde_pickle::DeOptions::new())
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SMResult {
    pub temperature: f64,
    pub mutation: String,
    pub atom_names: Vec<String>,
    pub atom_res: Vec<usize>,
    pub residues: Vec<Residue>,
    pub ndx_lig: Option<BTreeSet<usize>>,
    pub times: Vec<f64>,
    pub times_ie: Vec<f64>,
    pub coordinates: Array3<f64>,
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
               residues: &Vec<Residue>, ndx_lig: &Option<BTreeSet<usize>>, 
               times: &Vec<f64>, times_ie: &Vec<f64>, coordinates: &Array3<f64>, 
               mutation: &str, temperature: f64, 
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
            temperature,
            mutation: mutation.to_string(),
            atom_names: atom_names.to_vec(),
            atom_res: atom_res.to_vec(),
            residues: residues.to_owned(),
            ndx_lig: ndx_lig.to_owned(),
            times: times.to_owned(),
            times_ie: times_ie.to_owned(),
            coordinates: coordinates.to_owned(),
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
}

pub fn analyze_controller(sm_results: &SMResults, sys_name: &str, settings: &Settings) {
    loop {
        println!("\nTime range: {} - {} ns, step = {} ns", 
            sm_results.sm_results[0].times[0], 
            sm_results.sm_results[0].times.last().unwrap(), 
            if sm_results.sm_results[0].times.len() > 1 {
                sm_results.sm_results[0].times[1] - sm_results.sm_results[0].times[0]
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
                let ts_ids = get_time_index(&sm_results.sm_results[0].times, tmin, tmax);
                if ts_ids.is_empty() {
                    println!("Not valid time.");
                    continue;
                }
                println!("Writing pdb and pml file(s)...");
                for sm_result in &sm_results.sm_results {
                    let def_name = format!("MMPBSA_binding_energy_{}.pdb", sys_name);
                    write_pdb_with_bf(sm_result, &def_name, &ts_ids, &env::current_dir().unwrap(), 
                            &(0..sm_result.atom_res.len()).collect(), true);
                    let pml_name = format!("MMPBSA_binding_energy_{}.pml", sys_name);
                    let png_name = format!("MMPBSA_binding_energy_{}", sys_name);
                    write_pml(&pml_name, &def_name, &png_name, &env::current_dir().unwrap(), settings);
                }
                println!("Finished writing pdb file(s) with binding energy information.");
                println!("Finished drawing figures with pml file(s) by PyMOL.");
            },
            Ok(0) => exit(0),
            Ok(1) => {
                println!("Input the time period (in ns, e.g. 0-40) to output (default: average):");
                let (tmin, tmax) = get_time_range();
                for result in &sm_results.sm_results {
                    analyze_summary(result, &env::current_dir().unwrap(), 
                            &format!("{}-{}", sys_name, result.mutation), tmin, tmax)
                }
            },
            Ok(2) => {
                for result in &sm_results.sm_results {
                    analyze_traj(result, &env::current_dir().unwrap(), &format!("{}-{}", sys_name, result.mutation))
                }
            },
            Ok(3) => {
                println!("Input the time period (in ns, e.g. 0-40) to output (default: average):");
                let (tmin, tmax) = get_time_range();
                let ts_ids = get_time_index(&sm_results.sm_results[0].times, tmin, tmax);
                if ts_ids.is_empty() {
                    println!("Not valid time.");
                    continue;
                }
                let (range_des, target_res) = select_res_by_range(&sm_results.sm_results[0]);
                println!("Writing energy file(s)...");
                for result in &sm_results.sm_results {
                    analyze_res(result, &env::current_dir().unwrap(), 
                        &format!("{}-{}", sys_name, result.mutation), &ts_ids, &range_des, &target_res);
                }
                println!("Finished writing residue-wised binding energy file(s).");
            },
            Ok(4) => {
                for result in &sm_results.sm_results {
                    analyze_atom(result, &env::current_dir().unwrap(), &format!("{}-{}", sys_name, result.mutation))
                }
                if sm_results.sm_results[0].ndx_lig.is_some() {
                    println!("Finished writing atom-wised binding energy pdb file(s) for ligand.");
                } else {
                    println!("System does not contain ligands.");
                }
            },
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

fn write_pdb_with_bf(result: &SMResult, def_name: &String, ts_ids: &Vec<usize>, wd: &Path, atom_range: &BTreeSet<usize>, reverse: bool) {
    let mut f = fs::File::create(wd.join(def_name)).unwrap();
    let coord = &result.coordinates;
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

fn analyze_summary(results: &SMResult, wd: &Path, sys_name: &String, tmin: f64, tmax: f64) {
    let ts_ids = get_time_index(&results.times, tmin, tmax);
    let ts_ie_ids = get_time_index(&results.times_ie, tmin, tmax);
    if ts_ids.is_empty() {
        println!("Not valid time.");
        return;
    }

    let beta_kj = 1000.0 / 8.314462618 / results.temperature;
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

    let def_name = format!("MMPBSA_{}({}-{}ns).csv", sys_name, results.times[ts_ids[0]], results.times[*ts_ids.last().unwrap()]);
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

    let def_name = format!("MMPBSA_{}_ΔH_traj.png", sys_name);
    let output_path = wd.join(&def_name);

    // 创建绘图区域
    let root = BitMapBackend::new(output_path.to_str().unwrap(), (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    // 找到数据的范围
    let x_min = results.times.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let x_max = results.times.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let y_min = results.dh.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let y_max = results.dh.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    // 添加一些边距
    let x_range = x_max - x_min;
    let y_range = y_max - y_min;
    let x_min = x_min - x_range * 0.05;
    let x_max = x_max + x_range * 0.05;
    let y_min = y_min - y_range * 0.05;
    let y_max = y_max + y_range * 0.05;

    // 创建图表
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Binding Energy - {}", sys_name), ("sans-serif", 24).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(90)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();

    // 配置网格和标签
    chart
        .configure_mesh()
        .x_desc("Time (ns)")
        .y_desc("Binding Energy (kJ/mol)")
        .axis_desc_style(("sans-serif", 22).into_font())
        .x_label_style(("sans-serif", 20).into_font())
        .y_label_style(("sans-serif", 20).into_font())
        .light_line_style(&BLACK.mix(0.0))
        .y_label_formatter(&|v| {
            if *v < 0.0 {
                // 将负号替换为规范的 minus 符号 (U+2212)
                format!("−{:.1}", v.abs())
            } else {
                format!("{:.1}", v)
            }
        }).draw()
        .unwrap();

    // 绘制曲线
    chart
        .draw_series(LineSeries::new(
            results.times.iter().zip(results.dh.iter()).map(|(&x, &y)| (x, y)),
            &BLUE.mix(0.75),
        ))
        .unwrap()
        .label("Binding Energy")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE.mix(0.75)));
    
    // 绘制点（可选，如果需要显示数据点）
    chart
        .draw_series(
            results.times.iter().zip(results.dh.iter()).map(|(&x, &y)| {
                Circle::new((x, y), 3, BLUE.filled())
            })
        )
        .unwrap();

    // 添加图例
    chart
        .configure_series_labels()
        .label_font(("sans-serif", 20).into_font())
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .position(SeriesLabelPosition::LowerRight)
        .draw()
        .unwrap();

    // 保存并完成
    root.present().unwrap();
    println!("Binding energy figure drawn to {}", def_name);
}

fn select_res_by_range(results: &SMResult) -> (String, Vec<usize>) {
    let res_range = loop {
        println!("Input the residue range you want to output (e.g., 1-3, 5), default: all");
        println!("Input \"?\" to view the residues list");
        
        let input = get_input(String::new());
        if input == "?" {
            results.residues.iter().enumerate().for_each(|(i, r)| {
                print!("{:6}{:<3},", i + 1, r.name);
                if (i + 1) % 5 == 0 {
                    println!();
                }
            });
            println!();
        } else {
            break input;
        }
    };

    // target_res 要从 0 开始数
    let (range_des, target_res): (String, Vec<usize>) = if res_range.trim().is_empty() {
        (
            "all".to_string(),
            (0..results.residues.len()).collect()
        )
    } else {
        (
            res_range.to_string(),
            range2list(&res_range).iter().map(|&i| i as usize - 1).collect()
        )
    };

    (range_des, target_res)
}

fn analyze_res(results: &SMResult, wd: &Path, sys_name: &String, ts_ids: &Vec<usize>, range_des: &String, target_res: &Vec<usize>) {
    let def_name = format!("MMPBSA_{}_res_{}.csv", sys_name, range_des);
    let (tar_res_nr, tar_res_name, tar_res_energy) = get_target_res_data(results, ts_ids, target_res);
    write_res_csv(results, &tar_res_nr, &tar_res_name, ts_ids, &tar_res_energy, wd, &def_name);
}

fn analyze_atom(results: &SMResult, wd: &Path, sys_name: &String) {
    let def_name = format!("MMPBSA_{}_ligand.pdb", sys_name);
    if let Some(ndx_lig) = &results.ndx_lig {
        write_pdb_with_bf(results, &def_name, &vec![0], wd, ndx_lig, false);
    }
}

fn get_target_res_data(results: &SMResult, ts_ids: &Vec<usize>, target_res: &Vec<usize>) -> (Vec<i32>, Vec<String>, Array2<f64>) {
    let (tar_res_nr, tar_res_name): (Vec<i32>, Vec<String>) = results.residues
        .iter()
        .enumerate()
        .filter(|(id, _)| target_res.contains(id))
        .map(|(_, res)| (res.nr, res.name.clone()))
        .unzip();
    let mut dh_res = Array2::zeros((target_res.len(), ts_ids.len()));
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
    target_res.iter().enumerate().for_each(|(i, &res_id)| {
        let atom_ids = get_cur_res_atom_ids(res_id);
        let dh_resi: Array1<f64> = results.dh_atom.select(Axis(1), &atom_ids).select(Axis(0), &ts_ids).sum_axis(Axis(1));
        dh_res.row_mut(i).assign(&dh_resi);
    });
    (tar_res_nr, tar_res_name, dh_res)
}

fn write_res_csv(results: &SMResult, tar_res_nr: &Vec<i32>, tar_res_name: &Vec<String>, 
                 ts_ids: &Vec<usize>, tar_res_energy: &Array2<f64>, wd: &Path, def_name: &String) {
    let mut res_energy_file = fs::File::create(wd.join(def_name)).unwrap();
    let columns = format!("resid,resname,{}\n",
        &ts_ids.iter().map(|&i| results.times[i].to_string()).collect::<Vec<String>>().join(","));
    write!(res_energy_file, "{}", columns).unwrap();
    for (i, res_energy) in tar_res_energy.rows().into_iter().enumerate() {
        write!(res_energy_file, "{},{},{}\n", 
            tar_res_nr[i], tar_res_name[i],
            res_energy.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(",")
        ).expect("Error while writing residue-wised energy file");
    }
}
