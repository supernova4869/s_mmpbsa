use std::path::Path;
use std::rc::Rc;
use xdrfile::{Frame, XTCTrajectory};
use crate::{get_input_selection, parameters::Parameters};
use crate::atom_radius::Radius;
use crate::fun_para_mmpbsa::set_para_mmpbsa;
use crate::index_parser::Index;
use crate::parse_tpr::TPR;

pub fn set_para_trj(trj: &String, tpr: &mut TPR, ndx: &String, wd: &Path, atom_radius: &Radius, settings: &mut Parameters) {
    let mut complex_grp: i32 = -1;
    let mut receptor_grp: i32 = -1;
    let mut ligand_grp: i32 = -1;
    println!("Reading trajectory...");
    let xtc = XTCTrajectory::open_read(trj).expect("Error reading trajectory");
    let frames: Vec<Rc<Frame>> = xtc.into_iter().map(|p| p.unwrap()).collect();
    let mut bt: f64 = frames[0].time as f64;
    let mut et: f64 = frames[frames.len() - 1].time as f64;
    let mut dt: f64 = (frames[1].time - frames[0].time) as f64;
    let index = Index::new(ndx);
    loop {
        println!("\n                 ************ Trajectory Parameters ************");
        println!("-10 Return");
        println!("  0 Go to next step");
        println!("  1 Select complex group, current:            {}", show_grp(complex_grp, &index));
        println!("  2 Select receptor groups, current:          {}", show_grp(receptor_grp, &index));
        println!("  3 Select ligand groups, current:            {}", show_grp(ligand_grp, &index));
        println!("  4 Set start time of analysis, current:      {} ns", bt / 1000.0);
        println!("  5 Set end time of analysis, current:        {} ns", et / 1000.0);
        println!("  6 Set time interval of analysis, current:   {} ns", dt / 1000.0);
        let i = get_input_selection();
        match i {
            -10 => return,
            0 => {
                set_para_mmpbsa(trj, tpr, ndx, wd,
                                complex_grp as usize,
                                receptor_grp as usize,
                                ligand_grp as usize,
                                bt, et, dt, atom_radius,
                                settings);
            }
            1 => {
                println!("Current groups:");
                index.list_groups();
                println!("Input complex group num:");
                complex_grp = get_input_selection();
            }
            2 => {
                println!("Current groups:");
                index.list_groups();
                println!("Input receptor group num:");
                receptor_grp = get_input_selection();
            }
            3 => {
                println!("Current groups:");
                index.list_groups();
                println!("Input ligand group num:");
                ligand_grp = get_input_selection();
            }
            4 => {
                println!("Input start time (ns), should be divisible of {} ns:", dt / 1000.0);
                let mut new_bt = get_input_selection::<f64>() * 1000.0;
                while (new_bt - bt) % dt != 0.0 || new_bt > frames[frames.len() - 1].time as f64 || new_bt < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_bt / 1000.0);
                    println!("Input start time (ns) again, should be divisible of {} ns:", dt / 1000.0);
                    new_bt = get_input_selection::<f64>() * 1000.0;
                }
                bt = new_bt;
            }
            5 => {
                println!("Input end time (ns), should be divisible of {} ns:", dt / 1000.0);
                let mut new_et = get_input_selection::<f64>() * 1000.0;
                while (new_et - et) % dt != 0.0 || new_et > frames[frames.len() - 1].time as f64 || new_et < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_et / 1000.0);
                    println!("Input end time (ns) again, should be divisible of {} ns:", dt / 1000.0);
                    new_et = get_input_selection::<f64>() * 1000.0;
                }
                et = new_et;
            }
            6 => {
                println!("Input interval time (ns), should be divisible of {} ns:", dt / 1000.0);
                let mut new_dt = get_input_selection::<f64>() * 1000.0;
                while new_dt % dt != 0.0 {
                    println!("The input {} ns is not a valid time step.", new_dt / 1000.0);
                    println!("Input interval time (ns) again, should be divisible of {} ns:", dt / 1000.0);
                    new_dt = get_input_selection::<f64>() * 1000.0;
                }
                dt = new_dt;
            }
            _ => println!("Invalid input")
        }
    }
}

fn show_grp(grp: i32, ndx: &Index) -> String {
    match grp {
        -1 => String::from("undefined"),
        _ => format!("{}): {}, {} atoms",
                     grp,
                     ndx.groups[grp as usize].name,
                     ndx.groups[grp as usize].indexes.len())
    }
}