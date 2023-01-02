use std::fs;
use std::io::{stdin, Write};
use std::path::Path;
use ndarray::Array1;
use crate::get_input_selection;
use crate::parameters::Parameters;

pub struct Results {
    pub times: Array1<f64>,     // ps
    pub mm: Array1<f64>,
    pub pb: Array1<f64>,
    pub sa: Array1<f64>,
    pub cou: Array1<f64>,
    pub vdw: Array1<f64>,
    pub dh: Array1<f64>,
    pub dh_res: Array1<f64>,
    pub mm_res: Array1<f64>,
    pub cou_res: Array1<f64>,
    pub vdw_res: Array1<f64>,
    pub dpb_res: Array1<f64>,
    pub dsa_res: Array1<f64>,
    pub pb_res: Array1<f64>,
    pub sa_res: Array1<f64>,
}

impl Results {
    // totally time average and ts
    fn summary(&self, temperature: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let rt2kj = 8.314462618 * temperature / 1e3;

        let dh_avg = self.dh.iter().sum::<f64>() / self.dh.len() as f64;
        let mm_avg = self.mm.iter().sum::<f64>() / self.mm.len() as f64;
        let cou_avg = self.cou.iter().sum::<f64>() / self.cou.len() as f64;
        let vdw_avg = self.vdw.iter().sum::<f64>() / self.vdw.len() as f64;
        let pb_avg = self.pb.iter().sum::<f64>() / self.pb.len() as f64;
        let sa_avg = self.sa.iter().sum::<f64>() / self.sa.len() as f64;

        let tds = self.mm.iter()
            .map(|&p| f64::exp((p - mm_avg) / rt2kj))
            .sum::<f64>() / self.mm.len() as f64;
        let tds = -rt2kj * tds.ln();
        let dg = dh_avg - tds;
        let ki = f64::exp(dg / rt2kj);
        return (dh_avg, mm_avg, pb_avg, sa_avg, cou_avg, vdw_avg, tds, dg, ki);
    }

    fn get_traj_energy(&self) -> (Array1<f64>, Array1<f64>, Array1<f64>, Array1<f64>, Array1<f64>, Array1<f64>) {
        (self.dh.clone(), self.mm.clone(), self.pb.clone(), self.sa.clone(), self.cou.clone(), self.vdw.clone())
    }
}

pub fn analyze_controller(results: &Results, temperature: f64, sys_name: &String, wd: &Path) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!(" 0 Return");
        println!(" 1 View binding energy terms summary");
        println!(" 2 View binding energy terms during trajectory");
        println!(" 3 View residue-wised binding energy terms summary");
        println!(" 4 View residue-wised binding energy (term: MM)");
        println!(" 5 View residue-wised binding energy (term: PB)");
        println!(" 6 View residue-wised binding energy (term: SA)");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            0 => break,
            1 => {
                analyze_summary(results, temperature, wd, sys_name);
            }
            2 => {
                analyze_traj(results, wd, sys_name);
            }
            _ => println!("Coming")
        }
    }
}

fn analyze_summary(results: &Results, temperature: f64, wd: &Path, sys_name: &String) {
    let (dh_avg, mm_avg, pb_avg, sa_avg, cou_avg,
        vdw_avg, tds, dg, ki) = results.summary(temperature);
    println!("Energy terms:");
    println!("ΔH: {:.3} kJ/mol", dh_avg);
    println!("ΔMM: {:.3} kJ/mol", mm_avg);
    println!("ΔPB: {:.3} kJ/mol", pb_avg);
    println!("ΔSA: {:.3} kJ/mol", sa_avg);
    println!();
    println!("Δele: {:.3} kJ/mol", cou_avg);
    println!("Δvdw: {:.3} kJ/mol", vdw_avg);
    println!();
    println!("TΔS: {:.3} kJ/mol", tds);
    println!("ΔG: {:.3} kJ/mol", dg);
    println!("Ki: {:.3}", ki);

    let f_name = format!("{}_MMPBSA.csv", sys_name);
    println!("Write to file? [Y/n]");
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
        energy_sum.write_all(format!("Δelec,{:.3},(kJ/mol)\n", cou_avg).as_bytes()).unwrap();
        energy_sum.write_all(format!("ΔvdW,{:.3},(kJ/mol)\n", vdw_avg).as_bytes()).unwrap();
        energy_sum.write_all(b"\n").unwrap();
        energy_sum.write_all(format!("TΔS,{:.3},(kJ/mol)\n", tds).as_bytes()).unwrap();
        energy_sum.write_all(format!("ΔG,{:.3},ΔG=ΔH-TΔS (kJ/mol)\n", dg).as_bytes()).unwrap();
        energy_sum.write_all(format!("Ki,{:.3e},Ki=exp(ΔG/RT)\n", ki).as_bytes()).unwrap();
        println!("Binding energy terms have been writen to {}", &f_name);
    }
}

fn analyze_traj(results: &Results, wd: &Path, sys_name: &String) {
    let (dh, mm, cou,
        vdw, pb, sa) = &results.get_traj_energy();
    let f_name = format!("{}_MMPBSA_traj.csv", sys_name);
    println!("Write to file? [Y/n]");
    let mut temp = String::new();
    stdin().read_line(&mut temp).unwrap();
    if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
        println!("Writing binding energy terms...");
        let mut energy_sum = fs::File::create(wd.join(&f_name)).unwrap();
        energy_sum.write_all("Time (ns),ΔH,ΔMM,ΔPB,ΔSA,Δelec,ΔvdW,(kJ/mol)\n"
            .as_bytes()).unwrap();
        for i in 0..dh.len() {
            energy_sum.write_all(format!("{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}\n",
                                         results.times[i] / 1000.0, dh[i], mm[i], pb[i], sa[i],
                                         cou[i], vdw[i]).as_bytes()).unwrap();
        }
        println!("Binding energy terms have been writen to {}", &f_name);
    }
}