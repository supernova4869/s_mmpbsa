use std::fs;
use std::io::Write;
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
    pub fn new(total_frames: usize, total_res_num: usize) -> Results {
        let mm: Array1<f64> = Array1::zeros(total_frames);
        let pb: Array1<f64> = Array1::zeros(total_frames);
        let sa: Array1<f64> = Array1::zeros(total_frames);
        let cou: Array1<f64> = Array1::zeros(total_frames);
        let vdw: Array1<f64> = Array1::zeros(total_frames);
        let dh: Array1<f64> = Array1::zeros(total_frames);

        let dh_res: Array1<f64> = Array1::zeros(total_res_num);
        let mm_res: Array1<f64> = Array1::zeros(total_res_num);
        let cou_res: Array1<f64> = Array1::zeros(total_res_num);
        let vdw_res: Array1<f64> = Array1::zeros(total_res_num);
        let dpb_res: Array1<f64> = Array1::zeros(total_res_num);
        let dsa_res: Array1<f64> = Array1::zeros(total_res_num);

        let pb_res: Array1<f64> = Array1::zeros(total_res_num);
        let sa_res: Array1<f64> = Array1::zeros(total_res_num);

        Results {
            mm,
            pb,
            sa,
            cou,
            vdw,
            dh,
            dh_res,
            mm_res,
            cou_res,
            vdw_res,
            dpb_res,
            dsa_res,
            pb_res,
            sa_res,
        }
    }

    // totally time average and ts
    fn summary(&self, pbe_set: &PBESet) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let rt2kj = 8.314462618 * pbe_set.temp / 1e3;

        let dh_total = self.dh.iter().sum::<f64>() / self.dh.len() as f64;
        let mm_total = self.mm.iter().sum::<f64>() / self.mm.len() as f64;
        let cou_total = self.cou.iter().sum::<f64>() / self.cou.len() as f64;
        let vdw_total = self.vdw.iter().sum::<f64>() / self.vdw.len() as f64;
        let pb_total = self.pb.iter().sum::<f64>() / self.pb.len() as f64;
        let sa_total = self.sa.iter().sum::<f64>() / self.sa.len() as f64;

        let tds_total = self.mm.iter()
            .map(|&p| f64::exp((p - mm_total) / rt2kj))
            .sum::<f64>() / self.mm.len() as f64;
        let tds_total = -rt2kj * tds_total.ln();
        let dg_total = dh_total - tds_total;
        let ki = f64::exp(dg_total / rt2kj);
        return (dh_total, mm_total, pb_total, sa_total, cou_total, vdw_total, tds_total, dg_total, ki);
    }
}

pub fn analyze_controller(wd: &Path, sys_name: &String, results: &Results, pbe_set: &PBESet) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!(" 0 Return");
        println!(" 1 Output binding energy terms summary to csv file");
        println!(" 2 Output binding energy terms of each time point to csv file");
        println!(" 3 Output residue-decomposed binding energy terms summary to csv file");
        println!(" 4 Output residue-decomposed binding energy (term: MM) to csv file");
        println!(" 5 Output residue-decomposed binding energy (term: PB) to csv file");
        println!(" 6 Output residue-decomposed binding energy (term: SA) to csv file");
        let sel_fun: i32 = get_input_selection();
        match sel_fun {
            0 => break,
            1 => {
                println!("Writing binding energy terms...");
                let (dh_total, mm_total, pb_total, sa_total,
                    cou_total, vdw_total, tds_total, dg_total, ki) = results.summary(pbe_set);
                let f_name = format!("{}_MMPBSA.csv", sys_name);
                let mut energy_sum = fs::File::create(wd.join(&f_name)).unwrap();
                energy_sum.write_all("Energy Term,kJ/mol,info\n".as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔH,{:.3},ΔH=ΔMM+ΔPB+ΔSA\n", dh_total).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔMM,{:.3},ΔMM=Δelectrostatic+Δvan der Waals\n", mm_total).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔPB,{:.3}\n", pb_total).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔSA,{:.3}\n", sa_total).as_bytes()).unwrap();
                energy_sum.write_all(b"\n").unwrap();
                energy_sum.write_all(format!("Δelectrostatic,{:.3}\n", cou_total).as_bytes()).unwrap();
                energy_sum.write_all(format!("Δvan der Waals,{:.3}\n", vdw_total).as_bytes()).unwrap();
                energy_sum.write_all(b"\n").unwrap();
                energy_sum.write_all(format!("TΔS,{:.3}\n", tds_total).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔG,{:.3},ΔG=ΔH-TΔS\n", dg_total).as_bytes()).unwrap();
                energy_sum.write_all(format!("Ki,{:.3e}\n", ki).as_bytes()).unwrap();
                println!("Binding energy terms have been writen to {}", &f_name);
            }
            _ => println!("Coming")
        }
    }
}
