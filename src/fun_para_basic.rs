use std::io::stdin;
use std::path::Path;
use crate::parameters::Parameters;
use crate::{get_input_value, convert_cur_dir, confirm_file_validity};
use crate::fun_para_trj::set_para_trj;
use crate::parse_tpr::TPR;

pub fn set_para_basic(trj: &String, tpr: &TPR, ndx: &String, wd: &Path, settings: &mut Parameters) {
    let mut trj = String::from(trj);
    let mut ndx = String::from(ndx);

    loop {
        println!("\n                 ************ MM/PB-SA Files ************");
        println!(" 0 Go to next step");
        println!(" 1 Assign trajectory file (xtc or trr), current: {}", match trj.len() {
            0 => "undefined",
            _ => trj.as_str()
        });
        println!(" 2 Assign index file (ndx), current: {}", match ndx.len() {
            0 => "undefined",
            _ => ndx.as_str()
        });
        println!(" 3 Exit program");
        let i = get_input_value();
        match i {
            0 => {
                if trj.len() == 0 {
                    println!("Trajectory file not assigned.");
                } else if ndx.len() == 0 {
                    // 可能要改, 以后不需要index也能算
                    println!("Index file not assigned.");
                } else {
                    // go to next step
                    set_para_trj(&trj, &tpr, &ndx, &wd, settings);
                }
            }
            1 => {
                println!("Input trajectory file path (if in the same directory with tpr, then simply input (e.g.) `?md.xtc`):");
                trj.clear();
                stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                trj = convert_cur_dir(&trj, &settings);
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr"], &settings);
            }
            2 => {
                println!("Input index file path (if in the same directory with tpr, then simply input (e.g.) `?index.ndx`):");
                ndx.clear();
                stdin().read_line(&mut ndx).expect("Failed while reading index file");
                ndx = convert_cur_dir(&ndx, &settings);
                ndx = confirm_file_validity(&mut ndx, vec!["ndx"], &settings);
            }
            3 => break,
            _ => println!("Error input.")
        };
    }
}