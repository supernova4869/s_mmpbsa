use std::fs;
use std::io::Write;
use crate::get_input_value;

pub fn analyze_controller(sys_name: &String, results: (f64, f64, f64, f64, f64, f64, f64, f64, f64)) {
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
                let mut energy_sum = fs::File::create(f_name.as_str()).unwrap();
                energy_sum.write_all(format!("Energy Term,kJ/mol,info\n").as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔH,{:.3},ΔH=ΔMM+ΔPB+ΔSA\n", dH).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔMM,{:.3},ΔMM=ΔCoulomb+Δvan der Waals\n", MM).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔPB,{:.3}\n", PB).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔSA,{:.3}\n", SA).as_bytes()).unwrap();
                energy_sum.write_all(b"\n").unwrap();
                energy_sum.write_all(format!("ΔCoulomb,{:.3}\n", COU).as_bytes()).unwrap();
                energy_sum.write_all(format!("Δvan der Waals,{:.3}\n", VDW).as_bytes()).unwrap();
                energy_sum.write_all(b"\n").unwrap();
                energy_sum.write_all(format!("TΔS,{:.3}\n", TdS).as_bytes()).unwrap();
                energy_sum.write_all(format!("ΔG,{:.3}\n", dG).as_bytes()).unwrap();
                energy_sum.write_all(format!("Ki,{:.3}\n", Ki).as_bytes()).unwrap();
                println!("Binding energy terms have been writen to {}", f_name);
            },
            _ => println!("Coming")
        }
    }
}
