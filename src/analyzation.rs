use std::fs;
use std::io::{stdin, Write};
use std::path::Path;
use ndarray::Array1;
use crate::apbs_param::PBESet;
use crate::get_input_selection;

pub struct Results {
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
}

pub fn analyze_controller(wd: &Path, sys_name: &String, results: &Results, temperature: f64) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!(" 0 Return");
        println!(" 1 Output binding energy terms summary to csv file");
        println!(" 2 Output binding energy terms during trajectory to csv file");
        println!(" 3 Output residue-decomposed binding energy terms summary to csv file");
        println!(" 4 Output residue-decomposed binding energy (term: MM) to csv file");
        println!(" 5 Output residue-decomposed binding energy (term: PB) to csv file");
        println!(" 6 Output residue-decomposed binding energy (term: SA) to csv file");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            0 => break,
            1 => {
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
                println!("Write to {}? [Y/n]", f_name);
                let mut temp = String::new();
                stdin().read_line(&mut temp).unwrap();
                if temp.trim().is_empty() || temp.trim() == "Y" || temp.trim() == "y" {
                    println!("Writing binding energy terms...");
                    let mut energy_sum = fs::File::create(wd.join(&f_name)).unwrap();
                    energy_sum.write_all("Energy Term,kJ/mol,info\n".as_bytes()).unwrap();
                    energy_sum.write_all(format!("ΔH,{:.3},ΔH=ΔMM+ΔPB+ΔSA\n", dh_avg).as_bytes()).unwrap();
                    energy_sum.write_all(format!("ΔMM,{:.3},ΔMM=Δelectrostatic+Δvan der Waals\n", mm_avg).as_bytes()).unwrap();
                    energy_sum.write_all(format!("ΔPB,{:.3}\n", pb_avg).as_bytes()).unwrap();
                    energy_sum.write_all(format!("ΔSA,{:.3}\n", sa_avg).as_bytes()).unwrap();
                    energy_sum.write_all(b"\n").unwrap();
                    energy_sum.write_all(format!("Δelectrostatic,{:.3}\n", cou_avg).as_bytes()).unwrap();
                    energy_sum.write_all(format!("Δvan der Waals,{:.3}\n", vdw_avg).as_bytes()).unwrap();
                    energy_sum.write_all(b"\n").unwrap();
                    energy_sum.write_all(format!("TΔS,{:.3}\n", tds).as_bytes()).unwrap();
                    energy_sum.write_all(format!("ΔG,{:.3},ΔG=ΔH-TΔS\n", dg).as_bytes()).unwrap();
                    energy_sum.write_all(format!("Ki,{:.3e},Ki=exp(ΔG/RT)\n", ki).as_bytes()).unwrap();
                    println!("Binding energy terms have been writen to {}", &f_name);
                }
            }
            _ => println!("Coming")
        }
    }
}
