use std::collections::HashSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::process::exit;
use ndarray::{Array1, Array2, Array3};
use crate::atom_property::AtomProperty;
use crate::parse_tpr::Residue;
use crate::settings::Settings;
use crate::utils::{get_input, get_input_selection, range2list, get_outfile};

pub struct Results {
    pub aps: AtomProperty,
    pub residues: Vec<Residue>,
    pub ndx_rec: Vec<usize>,
    pub ndx_lig: Vec<usize>,
    pub times: Array1<f64>,
    pub coord: Array3<f64>,
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
    pub fn new(aps: &AtomProperty, residues: &Vec<Residue>,
               ndx_rec: &Vec<usize>, ndx_lig: &Vec<usize>,
               times: &Array1<f64>, coord: Array3<f64>, 
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
            ndx_rec: ndx_rec.to_owned(),
            ndx_lig: ndx_lig.to_owned(),
            times: times.to_owned(),
            coord: coord.to_owned(),
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
    fn summary(&self, temperature: f64, settings: &Settings) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
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

pub fn analyze_controller(results: &Results, temperature: f64, sys_name: &String, wd: &Path, total_at_num: usize, settings: &Settings) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!("-1 Write residue-wised bind energy to pdb file");
        println!(" 0 Exit program");
        println!(" 1 View binding energy terms summary");
        println!(" 2 Output binding energy terms by trajectory");
        println!(" 3 Output average binding energy terms by residue");
        println!(" 4 Output residue-wised binding energy terms by time as default names");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            -1 => write_energy_to_bf(results, wd, sys_name, total_at_num),
            0 => exit(0),
            1 => analyze_summary(results, temperature, wd, sys_name, settings),
            2 => analyze_traj(results, wd, sys_name),
            3 => analyze_res(results, wd, sys_name),
            4 => output_all_details(results, wd, sys_name),
            _ => println!("Invalid input")
        }
    }
}

fn write_energy_to_bf(results: &Results, wd: &Path, sys_name: &String, total_at_num: usize) {
    println!("Input the time point (in ns) to write pdb (default: all):");
    let ts = get_input(-1.0);
    println!("Writing pdb file(s)...");
    if ts != -1.0 {
        let ts_id = get_time_index(ts, results);
        write_bf_pdb(results, sys_name, ts_id, wd, total_at_num);
    } else {
        for ts_id in 0..results.times.len() {
            write_bf_pdb(results, sys_name, ts_id, wd, total_at_num);
        }
    }
    println!("Finished writing pdb file(s) with binding energy information.");
}

fn write_bf_pdb(results: &Results, sys_name: &String, ts_id: usize, wd: &Path, total_at_num: usize) {
    let mut f = fs::File::create(wd.join(&format!("binding_energy_{}_{}ns.pdb", sys_name, results.times[ts_id] / 1000.0))).unwrap();
    let coord = &results.coord;
    writeln!(f, "REMARK  The B-factor column is filled with the INVERSED residue-wised binding energy (ΔH), in kcal/mol").unwrap();
    for atom_id in 0..total_at_num {
        let res_id = results.aps.atm_resid[atom_id];
        let atom_name = results.aps.atm_name[atom_id].as_str();
        write_atom_line(res_id, atom_id, atom_name, &results, 
            coord[[ts_id, atom_id, 0]], coord[[ts_id, atom_id, 1]], coord[[ts_id, atom_id, 2]], &mut f);
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
                            results.times[i] / 1000.0, results.dh[i],
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
    let target_residues = match i {
        1 => {
            get_residue_range(results, 3.0)
        },
        2 => {
            get_residue_range(results, 5.0)
        },
        3 => {
            println!("Input the cut-off distance you want to expand, default: 3");
            let cutoff: f64 = get_input(3.0);
            get_residue_range(results, cutoff)
        },
        4 => {
            println!("Input the residue range you want to output (e.g., 1-3, 5), default: all");
            let res_range = get_input(String::new());
            let res_range: Vec<i32> = match res_range.len() {
                0 => results.residues.iter().map(|r| r.nr).collect(),
                _ => range2list(&res_range)
            };
            results.aps.atm_resid
                .iter()
                .filter(|&i| res_range.contains(&(results.residues[*i].nr)))    // 用户筛选用nr
                .map(|&i| results.residues[i].id)     // 索引用id
                .collect()
        },
        _ => {
            println!("Invalid selection");
            return
        }
    };
    
    // 分析输出
    println!("Writing binding energy terms...");
    let def_name = get_outfile(&format!("MMPBSA_{}_res.csv", sys_name));
    let mut energy_res = fs::File::create(wd.join(&def_name)).unwrap();
    energy_res.write_all("id,name,ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW\n".as_bytes()).unwrap();
    // 各项取平均输出
    let avg_dh_res: Array1<f64> = results.dh_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    let avg_mm_res: Array1<f64> = results.mm_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    let avg_pb_res: Array1<f64> = results.pb_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    let avg_sa_res: Array1<f64> = results.sa_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    let avg_elec_res: Array1<f64> = results.elec_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    let avg_vdw_res: Array1<f64> = results.vdw_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for (i, res) in results.residues.iter().enumerate() {
        if !target_residues.contains(&res.id) {
            continue;
        }
        write!(energy_res, "{},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n", 
            res.nr, res.name, avg_dh_res[i], avg_mm_res[i], avg_pb_res[i], avg_sa_res[i], avg_elec_res[i], avg_vdw_res[i])
            .expect("Error while writing residue-wised energy file");
    }
    println!("Binding energy terms have been writen to {}", &def_name);
}

fn get_residue_range(results: &Results, cutoff: f64) -> HashSet<usize> {
    let mut res_range: HashSet<usize> = HashSet::new();
    let total_frames = results.times.len() - 1;
    let ligand_x: Vec<f64> = results.ndx_lig.iter().map(|&a| results.coord[[total_frames, a, 0]]).collect();
    let ligand_y: Vec<f64> = results.ndx_lig.iter().map(|&a| results.coord[[total_frames, a, 1]]).collect();
    let ligand_z: Vec<f64> = results.ndx_lig.iter().map(|&a| results.coord[[total_frames, a, 2]]).collect();
    for res in &results.residues {
        let atoms_id: Vec<usize> = results.aps.atm_resid.iter()
            .enumerate()
            .filter(|&(_, &a)| a == res.id)
            .map(|(index, _)| index)
            .collect();
        for i in atoms_id {
            let x = results.coord[[total_frames, i, 0]];
            let y = results.coord[[total_frames, i, 1]];
            let z = results.coord[[total_frames, i, 2]];
            for j in 0..results.ndx_lig.len() {
                if (x - ligand_x[j]).powi(2) + (y - ligand_y[j]).powi(2) + (z - ligand_z[j]).powi(2) < cutoff.powi(2) {
                    res_range.insert(results.aps.atm_resid[i]);
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
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
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
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
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
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
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
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
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
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
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
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
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

fn get_time_index(ts: f64, results: &Results) -> usize {
    ((ts * 1000.0 - results.times[0]) / (results.times[1] - results.times[0])) as usize
}
