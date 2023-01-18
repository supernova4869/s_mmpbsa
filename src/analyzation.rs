use std::fs;
use std::io::{stdin, Write};
use std::path::Path;
use ndarray::Array1;
use crate::atom_property::AtomProperty;
use crate::get_input_selection;
use crate::mmpbsa::get_residues;
use crate::parse_tpr::TPR;

pub struct Results {
    pub times: Array1<f64>,
    pub residues: Array1<(i32, String)>,
    pub mm: Array1<f64>,
    pub pb: Array1<f64>,
    pub sa: Array1<f64>,
    pub elec: Array1<f64>,
    pub vdw: Array1<f64>,
    pub dh: Array1<f64>,
    pub dh_res: Array1<f64>,
    pub mm_res: Array1<f64>,
    pub elec_res: Array1<f64>,
    pub vdw_res: Array1<f64>,
    pub pb_res: Array1<f64>,
    pub sa_res: Array1<f64>,
}

impl Results {
    pub fn new(tpr: &TPR, times: Array1<f64>, ndx_com: &Vec<usize>,
               mm: Array1<f64>, pb: Array1<f64>, sa: Array1<f64>,
               elec: Array1<f64>, vdw: Array1<f64>, dh: Array1<f64>,
               dh_res: Array1<f64>, mm_res: Array1<f64>, elec_res: Array1<f64>,
               vdw_res: Array1<f64>, pb_res: Array1<f64>, sa_res: Array1<f64>) -> Results {

        // residues number and name
        let residues = get_residues(tpr, ndx_com);        

        Results {
            times,
            residues,
            mm,
            pb,
            sa,
            elec,
            vdw,
            dh,
            dh_res,
            mm_res,
            elec_res,
            vdw_res,
            pb_res,
            sa_res,
        }
    }

    // totally time average and ts
    fn summary(&self, temperature: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let rt2kj = 8.314462618 * temperature / 1e3;

        let dh_avg = self.dh.iter().sum::<f64>() / self.dh.len() as f64;
        let mm_avg = self.mm.iter().sum::<f64>() / self.mm.len() as f64;
        let elec_avg = self.elec.iter().sum::<f64>() / self.elec.len() as f64;
        let vdw_avg = self.vdw.iter().sum::<f64>() / self.vdw.len() as f64;
        let pb_avg = self.pb.iter().sum::<f64>() / self.pb.len() as f64;
        let sa_avg = self.sa.iter().sum::<f64>() / self.sa.len() as f64;

        let tds = self.mm.iter()
            .map(|&p| f64::exp((p - mm_avg) / rt2kj))
            .sum::<f64>() / self.mm.len() as f64;
        let tds = -rt2kj * tds.ln();
        let dg = dh_avg - tds;
        let ki = f64::exp(dg / rt2kj);
        return (dh_avg, mm_avg, pb_avg, sa_avg, elec_avg, vdw_avg, tds, dg, ki);
    }
}

pub fn analyze_controller(results: &Results, temperature: f64, sys_name: &String, wd: &Path) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!(" 0 Return");
        println!(" 1 View binding energy terms summary");
        println!(" 2 View binding energy terms by time");
        println!(" 3 View residue-wised binding energy summary");
        // println!(" 4 View residue-wised binding energy by time: dH");
        // println!(" 5 View residue-wised binding energy by time: MM");
        // println!(" 6 View residue-wised binding energy by time: PB");
        // println!(" 7 View residue-wised binding energy by time: SA");
        // println!(" 8 View residue-wised binding energy by time: elec");
        // println!(" 9 View residue-wised binding energy by time: vdW");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            0 => break,
            1 => {
                analyze_summary(results, temperature, wd, sys_name);
            }
            2 => {
                analyze_traj(results, wd, sys_name);
            }
            3 => {
                analyze_summary_res(results, wd, sys_name);
            }
            // 4 => {
            //     analyze_dh_res_traj(results, wd, sys_name);
            // }
            // 5 => {
            //     analyze_mm_res_traj(results, wd, sys_name);
            // }
            // 6 => {
            //     analyze_pb_res_traj(results, wd, sys_name);
            // }
            // 7 => {
            //     analyze_sa_res_traj(results, wd, sys_name);
            // }
            // 8 => {
            //     analyze_elec_res_traj(results, wd, sys_name);
            // }
            // 9 => {
            //     analyze_vdw_res_traj(results, wd, sys_name);
            // }
            _ => println!("Invalid input")
        }
    }
}

fn analyze_summary(results: &Results, temperature: f64, wd: &Path, sys_name: &String) {
    let (dh_avg, mm_avg, pb_avg, sa_avg, elec_avg,
        vdw_avg, tds, dg, ki) = results.summary(temperature);
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
    println!("Ki: {:.3}", ki);

    let f_name = format!("{}_MMPBSA.csv", sys_name);
    println!("\nWrite to file? [Y/n]");
    let mut temp = String::new();
    stdin().read_line(&mut temp).unwrap();
    if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
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
        energy_sum.write_all(format!("Ki,{:.3e},Ki=exp(ΔG/RT)\n", ki).as_bytes()).unwrap();
        println!("Binding energy terms have been writen to {}", &f_name);
    }
}

fn analyze_traj(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = format!("{}_MMPBSA_traj.csv", sys_name);
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

fn analyze_summary_res(results: &Results, wd: &Path, sys_name: &String) {
    let f_name = format!("{}_MMPBSA_res.csv", sys_name);
    println!("Writing binding energy terms...");
    let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
    energy_res.write_all("Energy term (kJ/mol)".as_bytes()).unwrap();
    for (i, res) in &results.residues {
        energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔH".as_bytes()).unwrap();
    for dh in &results.dh_res {
        energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔMM".as_bytes()).unwrap();
    for mm in &results.mm_res {
        energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔPB".as_bytes()).unwrap();
    for pb in &results.pb_res {
        energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔSA".as_bytes()).unwrap();
    for sa in &results.sa_res {
        energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔelec".as_bytes()).unwrap();
    for elec in &results.elec_res {
        energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
    }
    energy_res.write_all("\nΔvdW".as_bytes()).unwrap();
    for vdw in &results.vdw_res {
        energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
    }
    energy_res.write_all("\n".as_bytes()).unwrap();
    println!("Binding energy terms have been writen to {}", &f_name);
}

// 考虑将de_elec等传入, 提供每帧输出

// fn analyze_dh_res_traj(results: &Results, wd: &Path, sys_name: &String) {
//     let f_name = format!("{}_res_dH.csv", sys_name);
//     let mut temp = String::new();
//     stdin().read_line(&mut temp).unwrap();
//     if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
//         println!("Writing binding energy terms...");
//         let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
//         energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//         for (i, res) in &results.residues {
//             energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
//         }
//         for i in 0..results.times.len() {
//             energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
//             for dh in &results.dh_res {
//                 energy_res.write_all(format!(",{:.3}", dh).as_bytes()).unwrap();
//             }
//         }
//         energy_res.write_all("\n".as_bytes()).unwrap();
//         println!("Binding energy terms have been writen to {}", &f_name);
//     }
// }
//
// fn analyze_mm_res_traj(results: &Results, wd: &Path, sys_name: &String) {
//     let f_name = format!("{}_res_MM.csv", sys_name);
//     let mut temp = String::new();
//     stdin().read_line(&mut temp).unwrap();
//     if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
//         println!("Writing binding energy terms...");
//         let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
//         energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//         for (i, res) in &results.residues {
//             energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
//         }
//         for i in 0..results.times.len() {
//             energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
//             for mm in &results.mm_res {
//                 energy_res.write_all(format!(",{:.3}", mm).as_bytes()).unwrap();
//             }
//         }
//         energy_res.write_all("\n".as_bytes()).unwrap();
//         println!("Binding energy terms have been writen to {}", &f_name);
//     }
// }
//
// fn analyze_pb_res_traj(results: &Results, wd: &Path, sys_name: &String) {
//     let f_name = format!("{}_res_PB.csv", sys_name);
//     let mut temp = String::new();
//     stdin().read_line(&mut temp).unwrap();
//     if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
//         println!("Writing binding energy terms...");
//         let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
//         energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//         for (i, res) in &results.residues {
//             energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
//         }
//         for i in 0..results.times.len() {
//             energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
//             for pb in &results.pb_res {
//                 energy_res.write_all(format!(",{:.3}", pb).as_bytes()).unwrap();
//             }
//         }
//         energy_res.write_all("\n".as_bytes()).unwrap();
//         println!("Binding energy terms have been writen to {}", &f_name);
//     }
// }
//
// fn analyze_sa_res_traj(results: &Results, wd: &Path, sys_name: &String) {
//     let f_name = format!("{}_res_SA.csv", sys_name);
//     let mut temp = String::new();
//     stdin().read_line(&mut temp).unwrap();
//     if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
//         println!("Writing binding energy terms...");
//         let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
//         energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//         for (i, res) in &results.residues {
//             energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
//         }
//         for i in 0..results.times.len() {
//             energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
//             for sa in &results.sa_res {
//                 energy_res.write_all(format!(",{:.3}", sa).as_bytes()).unwrap();
//             }
//         }
//         energy_res.write_all("\n".as_bytes()).unwrap();
//         println!("Binding energy terms have been writen to {}", &f_name);
//     }
// }
//
// fn analyze_elec_res_traj(results: &Results, wd: &Path, sys_name: &String) {
//     let f_name = format!("{}_res_elec.csv", sys_name);
//     let mut temp = String::new();
//     stdin().read_line(&mut temp).unwrap();
//     if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
//         println!("Writing binding energy terms...");
//         let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
//         energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//         for (i, res) in &results.residues {
//             energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
//         }
//         for i in 0..results.times.len() {
//             energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
//             for elec in &results.elec_res {
//                 energy_res.write_all(format!(",{:.3}", elec).as_bytes()).unwrap();
//             }
//         }
//         energy_res.write_all("\n".as_bytes()).unwrap();
//         println!("Binding energy terms have been writen to {}", &f_name);
//     }
// }
//
// fn analyze_vdw_res_traj(results: &Results, wd: &Path, sys_name: &String) {
//     let f_name = format!("{}_res_vdW.csv", sys_name);
//     let mut temp = String::new();
//     stdin().read_line(&mut temp).unwrap();
//     if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
//         println!("Writing binding energy terms...");
//         let mut energy_res = fs::File::create(wd.join(&f_name)).unwrap();
//         energy_res.write_all("Time (ns)".as_bytes()).unwrap();
//         for (i, res) in &results.residues {
//             energy_res.write_all(format!(",{}#{}", i, res).as_bytes()).unwrap();
//         }
//         for i in 0..results.times.len() {
//             energy_res.write_all(format!("\n{}", results.times[i] / 1000.0).as_bytes()).unwrap();
//             for vdw in &results.vdw_res {
//                 energy_res.write_all(format!(",{:.3}", vdw).as_bytes()).unwrap();
//             }
//         }
//         energy_res.write_all("\n".as_bytes()).unwrap();
//         println!("Binding energy terms have been writen to {}", &f_name);
//     }
// }
