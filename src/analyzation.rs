use std::fs;
use std::io::Write;
use std::path::Path;
use crate::get_input_value;

pub fn analyze_controller(wd:&Path, sys_name: &String, results: (f64, f64, f64, f64, f64, f64, f64, f64, f64)) {
    loop {
        println!("\n                 ************ MM-PBSA analyzation ************");
        println!(" 0 Return");
        println!(" 1 Output binding energy terms summary to csv file");
        println!(" 2 Output binding energy terms of each time point to csv file");
        println!(" 3 Output residue-decomposed binding energy terms summary to csv file");
        println!(" 4 Output residue-decomposed binding energy (term: MM) to csv file");
        println!(" 5 Output residue-decomposed binding energy (term: PB) to csv file");
        println!(" 6 Output residue-decomposed binding energy (term: SA) to csv file");
        // COU和VDW都要体现出来, 想想怎样列
        // 得为结果专门建一个数据结构
        let sel_fun: i32 = get_input_value();
        match sel_fun {
            0 => break,
            1 => {
                println!("Writing binding energy terms...");
                let (dH, MM, PB, SA, COU, VDW, TdS, dG, Ki) = results;
                let f_name = format!("{}_MMPBSA.csv", sys_name);
                let mut energy_sum = fs::File::create(wd.join(&f_name)).unwrap();
                energy_sum.write_all("Energy Term,kJ/mol,info\n".as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔH,{:.3},ΔH=ΔMM+ΔPB+ΔSA\n", dH).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔMM,{:.3},ΔMM=Δelectrostatic+Δvan der Waals\n", MM).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔPB,{:.3}\n", PB).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔSA,{:.3}\n", SA).as_bytes()).unwrap();
                energy_sum.write_all(b"\n").unwrap();
                energy_sum.write_all(format!("Δelectrostatic,{:.3}\n", COU).as_bytes()).unwrap();
                energy_sum.write_all(format!("Δvan der Waals,{:.3}\n", VDW).as_bytes()).unwrap();
                energy_sum.write_all(b"\n").unwrap();
                energy_sum.write_all(format!("TΔS,{:.3}\n", TdS).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔG,{:.3},ΔG=ΔH-TΔS\n", dG).as_bytes()).unwrap();
                energy_sum.write_all(format!("Ki,{:.3e}\n", Ki).as_bytes()).unwrap();
                println!("Binding energy terms have been writen to {}", &f_name);
            },
            _ => println!("Coming")
        }
    }
}
