use std::fs::{self, File};
use std::io::{stdin, Write};
use std::path::Path;
use ndarray::{Array1, Array2};
use crate::settings::Settings;
use crate::get_input_selection;

pub struct Results {
    pub times: Array1<f64>,
    pub coord: Array2<f64>,
    pub residues: Array1<(i32, String)>,
    pub atm_name: Array1<String>,
    pub atm_resnum: Array1<usize>,
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
    pub fn new(times: Array1<f64>, coord: Array2<f64>, residues: Array1<(i32, String)>,
               atm_name: Array1<String>, atm_resnum: Array1<usize>, elec_res: Array2<f64>,
               vdw_res: Array2<f64>, pb_res: Array2<f64>, sa_res: Array2<f64>) -> Results {
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

        let mm_res = &elec_res + &vdw_res;
        let dh_res = &mm_res + &pb_res + &sa_res;

        Results {
            times,
            coord,
            residues,
            atm_name,
            atm_resnum,
            dh,
            mm,
            pb,
            sa,
            elec,
            vdw,
            dh_res,
            mm_res,
            pb_res,
            sa_res,
            elec_res,
            vdw_res,
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

pub fn analyze_controller(results: &Results, temperature: f64, sys_name: &String, wd: &Path, settings: &Settings) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!("-1 Write residue-wised bind energy to pdb file");
        println!(" 0 Return");
        println!(" 1 View binding energy terms summary");
        println!(" 2 View binding energy terms by trajectory");
        println!(" 3 View residue-wised binding energy summary");
        println!(" 4 View residue-wised binding energy by time: ΔH");
        println!(" 5 View residue-wised binding energy by time: ΔMM");
        println!(" 6 View residue-wised binding energy by time: ΔPB");
        println!(" 7 View residue-wised binding energy by time: ΔSA");
        println!(" 8 View residue-wised binding energy by time: Δelec");
        println!(" 9 View residue-wised binding energy by time: ΔvdW");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            -1 => write_energy_to_bf(results, wd, sys_name),
            0 => break,
            1 => analyze_summary(results, temperature, wd, sys_name, settings),
            2 => analyze_traj(results, wd, sys_name),
            3 => analyze_res_avg(results, wd, sys_name),
            4 => analyze_dh_res_traj(results, wd, sys_name),
            5 => analyze_mm_res_traj(results, wd, sys_name),
            6 => analyze_pb_res_traj(results, wd, sys_name),
            7 => analyze_sa_res_traj(results, wd, sys_name),
            8 => analyze_elec_res_traj(results, wd, sys_name),
            9 => analyze_vdw_res_traj(results, wd, sys_name),
            _ => println!("Invalid input")
        }
    }
}

fn write_energy_to_bf(results: &Results, wd: &Path, sys_name: &String) {
    let mut f = fs::File::create(wd.join(format!("binding_energy_{}.pdb", sys_name))).unwrap();
    let coord = &results.coord;
    f.write_all("REMARK  The B-factor column is filled with the INVERSED residue-wised binding energy (ΔH), in kcal/mol\n".as_bytes()).unwrap();
    for atom_id in 0..coord.shape()[0] {
        let res_id = results.atm_resnum[atom_id];
        let atom_name = results.atm_name[atom_id].as_str();
        write_atom_line(res_id, atom_id, atom_name, &results, coord[[atom_id, 0]], coord[[atom_id, 1]], coord[[atom_id, 2]], &mut f);
    }
    println!("Finished writing binding energy information to {}", format!("binding_energy_{}.pdb", sys_name));
}

fn write_atom_line(res_id: usize, atom_id: usize, atom_name: &str, results: &Results, x: f64, y: f64, z: f64, f: &mut File) {
    let str = format!("ATOM  {:5} {:<4} {:<3} A{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           {:<2}\n",
                              atom_id + 1, atom_name, results.residues[res_id].1, results.residues[res_id].0, x, y, z, 
                              -results.dh_res[[results.dh_res.shape()[0] - 1, res_id]] / 4.18, atom_name.get(0..1).unwrap());
    f.write_all(str.as_bytes()).unwrap();
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

    let f_name = get_outfile(&format!("MMPBSA_{}.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_sum = fs::File::create(wd.join(&f_name)).unwrap();
    energy_sum.write_all("Energy Term,value,info\n".as_bytes()).unwrap();
    energy_sum.write_all(format!("ΔH,{:.3},ΔH=ΔMM+ΔPB+ΔSA (kJ/mol)\n", dh_avg).as_bytes()).unwrap();
    energy_sum.write_all(format!("ΔMM,{:.3},ΔMM=Δelec+ΔvdW (kJ/mol)\n", mm_avg).as_bytes()).unwrap();
    energy_sum.write_all(format!("ΔPB,{:.3},(kJ/mol)\n", pb_avg).as_bytes()).unwrap();
    energy_sum.write_all(format!("ΔSA,{:.3},(kJ/mol)\n", sa_avg).as_bytes()).unwrap();
    energy_sum.write_all(b"\n").unwrap();
    energy_sum.write_all(format!("Δelec,{:.3},(kJ/mol)\n", elec_avg).as_bytes()).unwrap();
    energy_sum.write_all(format!("ΔvdW,{:.3},(kJ/mol)\n", vdw_avg).as_bytes()).unwrap();
    energy_sum.write_all(b"\n").unwrap();
    energy_sum.write_all(format!("TΔS,{:.3},(kJ/mol)\n", tds).as_bytes()).unwrap();
    energy_sum.write_all(format!("ΔG,{:.3},ΔG=ΔH-TΔS (kJ/mol)\n", dg).as_bytes()).unwrap();
    energy_sum.write_all(format!("Ki,{:.3e},Ki=exp(ΔG/RT) (nM)\n", ki).as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_traj.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_sum = fs::File::create(wd.join(&f_name)).unwrap();
    energy_sum.write_all("Time (ns),ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW,(kJ/mol)\n"
        .as_bytes()).unwrap();
    for i in 0..results.times.len() {
        energy_sum.write_all(format!("{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n",
                                    results.times[i] / 1000.0, results.dh[i],
                                    results.mm[i], results.pb[i], results.sa[i],
                                    results.elec[i], results.vdw[i]).as_bytes()).unwrap();
    }
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_res_avg(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Energy term (kJ/mol)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔH".as_bytes()).unwrap();
    let avg_dh_res: Array1<f64> = results.dh_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for dh in &avg_dh_res {
        energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔMM".as_bytes()).unwrap();
    let avg_mm_res: Array1<f64> = results.mm_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for mm in &avg_mm_res {
        energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔPB".as_bytes()).unwrap();
    let avg_pb_res: Array1<f64> = results.pb_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for pb in &avg_pb_res {
        energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔSA".as_bytes()).unwrap();
    let avg_sa_res: Array1<f64> = results.sa_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for sa in &avg_sa_res {
        energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔelec".as_bytes()).unwrap();
    let avg_elec_res: Array1<f64> = results.elec_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for elec in &avg_elec_res {
        energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔvdW".as_bytes()).unwrap();
    let avg_vdw_res: Array1<f64> = results.vdw_res.columns().into_iter().map(|col| col.sum() / results.times.len() as f64).collect();
    for vdw in &avg_vdw_res {
        energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_dh_res_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res_ΔH.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
        for dh in &results.dh_res.row(i) {
            energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_mm_res_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res_ΔMM.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
        for mm in &results.mm_res.row(i) {
            energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_pb_res_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res_ΔPB.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
        for pb in &results.pb_res.row(i) {
            energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_sa_res_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res_ΔSA.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
        for sa in &results.sa_res.row(i) {
            energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_elec_res_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res_Δelec.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
        for elec in &results.elec_res.row(i) {
            energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

fn analyze_vdw_res_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = get_outfile(&format!("MMPBSA_{}_res_ΔvdW.csv", sys_name));
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Time (ns)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    for i in 0..results.times.len() {
        energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
        for vdw in &results.vdw_res.row(i) {
            energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
        }
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

pub fn get_outfile<T: ToString + std::fmt::Display>(default_name: &T) -> String {
    println!("\nInput file name to write (default: {}):", default_name);
    let mut temp = String::new();
    stdin().read_line(&mut temp).unwrap();
    match temp.trim().is_empty() {
        true => default_name.to_string(),
        _ => temp.trim().to_string()
    }
}

pub fn get_infile<T: ToString>(default_name: &T) -> String {
    let mut temp = String::new();
    stdin().read_line(&mut temp).unwrap();
    match temp.trim().is_empty() {
        true => default_name.to_string(),
        _ => temp.trim().to_string()
    }
}