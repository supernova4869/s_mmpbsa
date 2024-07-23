use std::collections::HashSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::process::exit;
use ndarray::{Array1, Array2, Array3};
use crate::atom_property::AtomProperties;
use crate::parse_tpr::Residue;
use crate::settings::Settings;
use crate::utils::{get_input, get_input_selection, range2list, get_outfile};

pub struct Results {
    pub aps: AtomProperties,
    pub residues: Vec<Residue>,
    pub ndx_lig: Vec<usize>,
    pub times: Array1<f64>,
    pub coord: Array3<f64>,
    pub mutation: String,
    pub dh: Array1<f64>,
    pub mm: Array1<f64>,
    pub pb: Array1<f64>,
    pub sa: Array1<f64>,
    pub elec: Array1<f64>,
    pub vdw: Array1<f64>,
    pub dh_res: Array2<f64>,
    pub mm_res: Array2<f64>,
    pub pb_res: Array2<f64>,
    pub sa_res: Array2<f64>,
    pub elec_res: Array2<f64>,
    pub vdw_res: Array2<f64>,
}

impl Results {
    pub fn new(aps: &AtomProperties, residues: &Vec<Residue>,
               ndx_lig: &Vec<usize>, times: &Array1<f64>, 
               coord: Array3<f64>, mutation: &str,
               elec_res: &Array2<f64>, vdw_res: &Array2<f64>, 
               pb_res: &Array2<f64>, sa_res: &Array2<f64>) -> Results {
        let mut dh: Array1<f64> = Array1::zeros(times.len());
        let mut mm: Array1<f64> = Array1::zeros(times.len());
        let mut pb: Array1<f64> = Array1::zeros(times.len());
        let mut sa: Array1<f64> = Array1::zeros(times.len());
        let mut elec: Array1<f64> = Array1::zeros(times.len());
        let mut vdw: Array1<f64> = Array1::zeros(times.len());
        for idx in 0..times.len() {
            elec[idx] = elec_res.row(idx).sum();
            vdw[idx] = vdw_res.row(idx).sum();
            mm[idx] = elec[idx] + vdw[idx];
            pb[idx] = pb_res.row(idx).iter().sum();
            sa[idx] = sa_res.row(idx).iter().sum();
            dh[idx] = mm[idx] + pb[idx] + sa[idx];
        }

        let mm_res: Array2<f64> = elec_res + vdw_res;
        let dh_res: Array2<f64> = &mm_res + pb_res + sa_res;

        Results {
            aps: aps.to_owned(),
            residues: residues.to_owned(),
            ndx_lig: ndx_lig.to_owned(),
            times: times.to_owned(),
            coord: coord.to_owned(),
            mutation: mutation.to_string(),
            dh,
            mm,
            pb,
            sa,
            elec,
            vdw,
            dh_res,
            mm_res,
            pb_res: pb_res.to_owned(),
            sa_res: sa_res.to_owned(),
            elec_res: elec_res.to_owned(),
            vdw_res: vdw_res.to_owned(),
        }
    }

    // totally time average and ts
    pub fn summary(&self, temperature: f64, settings: &Settings) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let rt2kj = 8.314462618 * temperature / 1e3;

        let dh_avg = self.dh.iter().sum::<f64>() / self.dh.len() as f64;
        let mm_avg = self.mm.iter().sum::<f64>() / self.mm.len() as f64;
        let elec_avg = self.elec.iter().sum::<f64>() / self.elec.len() as f64;
        let vdw_avg = self.vdw.iter().sum::<f64>() / self.vdw.len() as f64;
        let pb_avg = self.pb.iter().sum::<f64>() / self.pb.len() as f64;
        let sa_avg = self.sa.iter().sum::<f64>() / self.sa.len() as f64;

        let tds = match settings.use_ts {
            true => {
                -rt2kj * (self.mm.iter().map(|&p| f64::exp((p - mm_avg) / rt2kj)).sum::<f64>() / self.mm.len() as f64).ln()
            }
            false => 0.0
        };
        let dg = dh_avg - tds;
        let ki = f64::exp(dg / rt2kj) * 1e9;    // nM
        return (dh_avg, mm_avg, pb_avg, sa_avg, elec_avg, vdw_avg, tds, dg, ki);
    }
}

pub fn analyze_controller(result_wt: &Results, result_as: &Vec<Results>, temperature: f64, sys_name: &String, wd: &Path, settings: &Settings) {
    if result_as.is_empty() {
        loop {
            println!("\n                 ************ MM-PBSA analyzation ************");
            println!("-1 Write residue-wised bind energy at specific time to pdb file");
            println!(" 0 Exit program");
            println!(" 1 View binding energy terms summary");
            println!(" 2 Output binding energy terms by trajectory");
            println!(" 3 Output binding energy terms by residue at specific time");
            println!(" 4 Output residue-wised binding energy terms by time as default names");
            let sel_fun: i32 = get_input_selection();
            match sel_fun {
                -1 => {
                    let ts_ids = get_time_points(result_wt);
                    write_energy_to_bf(result_wt, &ts_ids, wd, sys_name);
                    println!("Finished writing pdb file(s) with binding energy information.");
                },
                0 => exit(0),
                1 => analyze_summary(result_wt, temperature, wd, sys_name, settings),
                2 => analyze_traj(result_wt, wd, sys_name),
                3 => analyze_res(result_wt, wd, sys_name),
                4 => output_all_details(result_wt, wd, sys_name),
                _ => println!("Invalid input")
            }
        }
    } else {
        let result_diffs: Vec<Results> = result_as.iter().map(|r| Results::new(
            &r.aps, &r.residues, &r.ndx_lig, &r.times, r.coord.to_owned(), &r.mutation,
            &(&r.elec_res - &result_wt.elec_res), 
            &(&r.vdw_res - &result_wt.vdw_res), 
            &(&r.pb_res - &result_wt.pb_res), 
            &(&r.sa_res - &result_wt.sa_res))).collect();
        loop {
            println!("\n                 ************ MM-PBSA analyzation (Alanine scanning) ************");
            println!("-1 Write ALL residue-wised bind energy DIFF at specific time to pdb file");
            println!(" 0 Exit program");
            println!(" 1 View ALL binding energy terms DIFF summary");
            println!(" 2 Output ALL binding energy terms DIFF by trajectory");
            println!(" 3 Output ALL binding energy terms DIFF by residue at specific time");
            println!(" 4 Output ALL residue-wised binding energy terms DIFF by time as default names");
            let sel_fun: i32 = get_input_selection();
            match sel_fun {
                -1 => {
                    let ts_ids = get_time_points(result_wt);
                    for r_df in &result_diffs {
                        write_energy_to_bf(r_df, &ts_ids, wd, &format!("{}-{}", sys_name, r_df.mutation))
                    }
                    println!("Finished writing pdb file(s) with binding energy information.");
                },
                0 => exit(0),
                1 => {
                    for r_df in &result_diffs {
                        analyze_summary(r_df, temperature, wd, &format!("{}-{}", sys_name, r_df.mutation), settings)
                    }
                },
                2 => {
                    for r_df in &result_diffs {
                        analyze_traj(r_df, wd, &format!("{}-{}", sys_name, r_df.mutation))
                    }
                },
                3 => {
                    for r_df in &result_diffs {
                        analyze_res(r_df, wd, &format!("{}-{}", sys_name, r_df.mutation))
                    }
                },
                4 => {
                    for r_df in &result_diffs {
                        output_all_details(r_df, wd, &format!("{}-{}", sys_name, r_df.mutation))
                    }
                },
                _ => println!("Invalid input")
            }
        }
    }
}

fn get_time_points(result: &Results) -> Vec<usize> {
    println!("Input the time point (in ns) to write pdb (default: all):");
    let ts = get_input(-1.0);
    println!("Writing pdb file(s)...");
    if ts != -1.0 {
        if let Some(ts_id) = get_time_index(ts, result) {
            vec![ts_id]
        } else {
            println!("Error input: {} ns", ts);
            vec![]
        }
    } else {
        Vec::from_iter(0..result.times.len())
    }
}

fn get_time_index(ts: f64, results: &Results) -> Option<usize> {
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

fn write_energy_to_bf(results: &Results, ts_ids: &Vec<usize>, wd: &Path, sys_name: &String) {
    for ts_id in ts_ids {
        write_bf_pdb(results, sys_name, *ts_id, wd);
    }
}

fn write_bf_pdb(results: &Results, sys_name: &String, ts_id: usize, wd: &Path) {
    let mut f = fs::File::create(wd.join(&format!("MMPBSA_binding_energy_{}_{}ns.pdb", sys_name, results.times[ts_id]))).unwrap();
    let coord = &results.coord;
    writeln!(f, "REMARK  The B-factor column is filled with the INVERSED residue-wised binding energy (ΔH), in kcal/mol").unwrap();
    for atom in &results.aps.atom_props {
        let res_id = atom.resid;
        let atom_name = atom.name.as_str();
        write_atom_line(res_id, atom.id, atom_name, &results, 
            coord[[ts_id, atom.id, 0]], coord[[ts_id, atom.id, 1]], coord[[ts_id, atom.id, 2]], &mut f);
    }
}

fn write_atom_line(res_id: usize, atom_id: usize, atom_name: &str, results: &Results, x: f64, y: f64, z: f64, f: &mut File) {
    writeln!(f, "ATOM  {:5} {:<4} {:<3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           {:<2}", 
                 atom_id + 1, atom_name, results.residues[res_id].name, results.residues[res_id].nr, x, y, z, 
                 -results.dh_res[[results.dh_res.shape()[0] - 1, res_id]] / 4.18, atom_name.get(0..1).unwrap()).unwrap();
}

fn analyze_summary(results: &Results, temperature: f64, wd: &Path, sys_name: &String, settings: &Settings) {
    let (dh_avg, mm_avg, pb_avg, sa_avg, elec_avg,
        vdw_avg, tds, dg, ki) = results.summary(temperature, settings);
    println!("Energy terms summary:");
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
    println!("Ki: {:.3} nM", ki);

    let def_name = get_outfile(&format!("MMPBSA_{}.csv", sys_name));
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
    write!(energy_sum, "Ki,{:.3e},Ki=exp(ΔG/RT) (nM)\n", ki).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_traj(results: &Results, wd: &Path, sys_name: &String) {
    println!("Writing binding energy terms...");
    let def_name = &get_outfile(&format!("MMPBSA_{}_traj.csv", sys_name));
    let mut energy_sum = fs::File::create(wd.join(&def_name)).unwrap();
    write!(energy_sum, "Time (ns),ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW,(kJ/mol)\n").unwrap();
    for i in 0..results.times.len() {
        write!(energy_sum, "{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n",
                            results.times[i], results.dh[i],
                            results.mm[i], results.pb[i], results.sa[i],
                            results.elec[i], results.vdw[i]).unwrap();
    }
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_res(results: &Results, wd: &Path, sys_name: &String) {
    println!("Determine the residue range to output:");
    println!(" 1 Ligand and receptor residues within 3 A");
    println!(" 2 Ligand and receptor residues within 5 A");
    println!(" 3 Ligand and receptor residues within a specified distance");
    println!(" 4 Self-defined residue range");
    // 残基范围确定
    let i: i32 = get_input_selection();
    let mut range_des = String::from("3A");
    let target_res = match i {
        1 => {
            get_residue_range(results, 3.0)
        },
        2 => {
            range_des = String::from("5A");
            get_residue_range(results, 5.0)
        },
        3 => {
            println!("Input the cut-off distance you want to expand, default: 3");
            let cutoff = get_input(3.0);
            range_des = format!("{:.1}A", cutoff);
            get_residue_range(results, cutoff)
        },
        4 => {
            println!("Input the residue range you want to output (e.g., 1-3, 5), default: all");
            let res_range = get_input(String::new());
            range_des = res_range.to_string();
            let res_range: Vec<i32> = match res_range.len() {
                0 => {
                    range_des = "all".to_string();
                    results.residues.iter().map(|r| r.nr).collect()
                },
                _ => range2list(&res_range)
            };
            results.aps.atom_props
                .iter()
                .filter(|&i| res_range.contains(&(results.residues[i.resid].nr)))    // 用户筛选用nr
                .map(|i| results.residues[i.resid].id)     // 索引用id
                .collect()
        },
        _ => {
            println!("Invalid selection");
            return
        }
    };
    
    println!("Input the time point (in ns) to output (default: average):");
    let ts = get_input(-1.0);
    if ts * 1000.0 > *results.times.last().unwrap() || (ts < 0.0 && ts != -1.0) {
        println!("Error input: {} ns", ts);
        return;
    }
    println!("Writing energy file(s)...");
    if ts != -1.0 {
        if let Some(ts_id) = get_time_index(ts, results) {
            let def_name = format!("MMPBSA_{}_res_{}_{}ns.csv", sys_name, range_des, results.times[ts_id]);
            write_res_csv(results, ts_id, wd, &target_res, &def_name);
        } else {
            println!("Error input: {} ns", ts);
            return;
        }
    } else {
        let def_name = format!("MMPBSA_{}_res_{}_avg.csv", sys_name, range_des);
        write_res_avg_csv(results, wd, &target_res, &def_name);
    }

    println!("Finished writing residue-wised binding energy file(s).");
}

fn write_res_csv(results: &Results, ts_id: usize, wd: &Path, target_res: &HashSet<usize>, def_name: &String) {
    let mut energy_res = fs::File::create(wd.join(def_name)).unwrap();
    energy_res.write_all("id,name,ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW\n".as_bytes()).unwrap();
    for (i, res) in results.residues.iter().enumerate() {
        if !target_res.contains(&res.id) {
            continue;
        }
        write!(energy_res, "{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n", 
            res.nr, res.name,
            results.dh_res[[ts_id, i]],
            results.mm_res[[ts_id, i]],
            results.pb_res[[ts_id, i]],
            results.sa_res[[ts_id, i]],
            results.elec_res[[ts_id, i]],
            results.vdw_res[[ts_id, i]])
            .expect("Error while writing residue-wised energy file");
    }
}

fn write_res_avg_csv(results: &Results, wd: &Path, target_res: &HashSet<usize>, def_name: &String) {
    let mut energy_res = fs::File::create(wd.join(def_name)).unwrap();
    energy_res.write_all("id,name,ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW\n".as_bytes()).unwrap();
    for (i, res) in results.residues.iter().enumerate() {
        if !target_res.contains(&res.id) {
            continue;
        }
        write!(energy_res, "{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n", 
            res.nr, res.name, 
            results.dh_res.column(i).mean().unwrap(),
            results.mm_res.column(i).mean().unwrap(),
            results.pb_res.column(i).mean().unwrap(),
            results.sa_res.column(i).mean().unwrap(),
            results.elec_res.column(i).mean().unwrap(),
            results.vdw_res.column(i).mean().unwrap())
            .expect("Error while writing residue-wised energy file");
    }
}

fn get_residue_range(results: &Results, cutoff: f64) -> HashSet<usize> {
    let mut res_range: HashSet<usize> = HashSet::new();
    let total_frames = results.times.len() - 1;
    let ligand_x: Vec<f64> = results.ndx_lig.iter().map(|&a| results.coord[[total_frames, a, 0]]).collect();
    let ligand_y: Vec<f64> = results.ndx_lig.iter().map(|&a| results.coord[[total_frames, a, 1]]).collect();
    let ligand_z: Vec<f64> = results.ndx_lig.iter().map(|&a| results.coord[[total_frames, a, 2]]).collect();
    for res in &results.residues {
        let atoms_id: Vec<usize> = results.aps.atom_props.iter()
            .filter_map(|a| if a.resid == res.id {
                Some(a.id)
            } else {
                None
            })
            .collect();
        for i in atoms_id {
            let x = results.coord[[total_frames, i, 0]];
            let y = results.coord[[total_frames, i, 1]];
            let z = results.coord[[total_frames, i, 2]];
            for j in 0..results.ndx_lig.len() {
                if (x - ligand_x[j]).powi(2) + (y - ligand_y[j]).powi(2) + (z - ligand_z[j]).powi(2) < cutoff.powi(2) {
                    res_range.insert(results.aps.atom_props[i].resid);
                }
            }
        }
    }
    res_range
}

fn analyze_dh_res_traj(results: &Results, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for dh in &results.dh_res.row(i) {
            energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_mm_res_traj(results: &Results, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for mm in &results.mm_res.row(i) {
            energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_pb_res_traj(results: &Results, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for pb in &results.pb_res.row(i) {
            energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_sa_res_traj(results: &Results, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for sa in &results.sa_res.row(i) {
            energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_elec_res_traj(results: &Results, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for elec in &results.elec_res.row(i) {
            energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn analyze_vdw_res_traj(results: &Results, wd: &Path, def_name: &String) {
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(def_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for res in &results.residues {
        energy_res.write_all(format!(",{}#{}", res.nr, res.name).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i]).as_bytes()).unwrap();
        for vdw in &results.vdw_res.row(i) {
            energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", def_name);
}

pub fn output_all_details(results: &Results, wd: &Path, sys_name: &String) {
    analyze_dh_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔH.csv", sys_name));
    analyze_mm_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔMM.csv", sys_name));
    analyze_pb_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔPB.csv", sys_name));
    analyze_sa_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔSA.csv", sys_name));
    analyze_elec_res_traj(results, wd, &format!("MMPBSA_{}_res_Δelec.csv", sys_name));
    analyze_vdw_res_traj(results, wd, &format!("MMPBSA_{}_res_ΔvdW.csv", sys_name));
}
