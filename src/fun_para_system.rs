use core::f64;
use std::env;
use colored::*;
use ndarray::{Array1, Array3};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::path::Path;
use std::collections::BTreeSet;

use crate::parameters::Config;
use crate::parse_pdb::PDB;
use crate::parse_gro::GRO;
use crate::settings::Settings;
use crate::utils::{append_new_name, get_input_selection};
use crate::fun_para_mmpbsa::set_para_mmpbsa;
use crate::parse_ndx::{Index, IndexGroup};
use crate::parse_tpr::TPR;
use crate::atom_property::AtomProperties;
use crate::parse_tpr::Residue;
use crate::utils::{convert_tpr, trjconv};
use crate::read_xtc::read_xtc;

pub fn set_para_trj(trj: &String, tpr: &mut TPR, ndx_name: &String, config: &Option<Config>, 
                    wd: &Path, tpr_name: &str, settings: &mut Settings) {
    let mut receptor_grp: Option<usize> = None;
    let mut ligand_grp: Option<usize> = None;
    let mut bt: f64 = 0.0;                                  // ps
    let mut et: f64 = f64::INFINITY;                        // ps
    let mut dt = 1000.0;                               // ps
    let mut ie_multi = 10;                             // multipli
    println!("Reading {} file...", ndx_name);
    let ndx = Index::from(ndx_name);
    loop {
        println!("\n                 ************ Trajectory Parameters ************");
        println!("-10 Return");
        println!(" -1 Toggle whether to fix PBC conditions, current: {}", settings.fix_pbc);
        if receptor_grp.is_none() {
            println!("{}", "  0 Go to next step (incomplete)".red().bold());
        } else if ligand_grp.is_none() {
            println!("{}", "  0 Go to next step (solvation only)".yellow().bold());
        } else {
            println!("{}", "  0 Go to next step (complete)".green().bold());
        }
        println!("  1 Select receptor group, current:           {}", show_grp(receptor_grp, &ndx));
        println!("  2 Select ligand group, current:             {}", show_grp(ligand_grp, &ndx));
        println!("  3 Set start time to analyze, current:       {} ns", bt / 1000.0);
        println!("  4 Set end time to analyze, current:         {} ns", et / 1000.0);
        println!("  5 Set time interval, current:               {} ns", dt / 1000.0);
        println!("  6 Set IE sampling rate multiple, current:   {}", ie_multi);
        let i = get_input_selection();
        match i {
            Ok(-10) => return,
            Ok(-1) => {
                settings.fix_pbc = !settings.fix_pbc;
            }
            Ok(0) => {
                if let Some(receptor_grp) = receptor_grp {
                    prepare_system_tpr(receptor_grp, ligand_grp, trj, tpr, &ndx, tpr_name, ndx_name, bt, et, dt, 
                        dt / ie_multi as f64, config, wd, settings);
                } else {
                    println!("Please select receptor groups.");
                };
            }
            Ok(1) => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input receptor group num:");
                receptor_grp = get_input_selection().ok();
                // 处理输入，直到获得有效值或 None
                while let Some(grp_num) = receptor_grp {
                    if grp_num < ndx.groups.len() {
                        break; // 输入有效，退出循环
                    } else {
                        println!("Error with group {} (maximum {}), input again", grp_num, ndx.groups.len() - 1);
                        receptor_grp = get_input_selection().ok();
                    }
                }
            }
            Ok(2) => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input ligand group num (directly enter for nothing):");
                ligand_grp = get_input_selection().ok();
                // 处理输入，直到获得有效值或 None
                while let Some(grp_num) = ligand_grp {
                    if grp_num < ndx.groups.len() {
                        break; // 输入有效，退出循环
                    } else {
                        println!("Error with group {} (maximum {}), input again", grp_num, ndx.groups.len() - 1);
                        ligand_grp = get_input_selection().ok();
                    }
                }
            }
            Ok(3) => {
                println!("Input start time (in ns):");
                let mut new_bt = get_input_selection::<f64>().unwrap_or(0.0) * 1000.0;
                while new_bt < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_bt / 1000.0);
                    println!("Input start time (in ns) again:");
                    new_bt = get_input_selection::<f64>().unwrap_or(0.0) * 1000.0;
                }
                bt = new_bt;
            }
            Ok(4) => {
                println!("Input end time (in ns):");
                let mut new_et = get_input_selection::<f64>().unwrap_or(f64::INFINITY) * 1000.0;
                while new_et < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_et / 1000.0);
                    println!("Input end time (in ns) again:");
                    new_et = get_input_selection::<f64>().unwrap_or(f64::INFINITY) * 1000.0;
                }
                et = new_et;
            }
            Ok(5) => {
                println!("Input interval time (in ns) for MM-PBSA:");
                let mut new_dt = get_input_selection::<f64>().unwrap_or(1.0) * 1000.0;
                while new_dt < 0.0 {
                    println!("The input {} ns is not a valid interval time.", new_dt / 1000.0);
                    println!("Input interval time (in ns) again:");
                    new_dt = get_input_selection::<f64>().unwrap_or(1.0) * 1000.0;
                }
                dt = new_dt;
            }
            Ok(6) => {
                println!("Input interpolation multiplier (positive integer) for IE:");
                let mut new_ie_multi = get_input_selection::<usize>().unwrap_or(10);
                while new_ie_multi <= 0 {
                    println!("Must be positive integer, input again:");
                    new_ie_multi = get_input_selection::<usize>().unwrap_or(10);
                }
                ie_multi = new_ie_multi;
            }
            _ => println!("Invalid input")
        }
    }
}

// convert rec and lig to begin at 0 and continous
pub fn normalize_index(ndx_rec: &BTreeSet<usize>, ndx_lig: &Option<BTreeSet<usize>>) -> (BTreeSet<usize>, Option<BTreeSet<usize>>) {
    if let Some(ndx_lig) = ndx_lig {
        // 利用BTreeSet的union和有序特性
        let union_set: BTreeSet<usize> = ndx_rec.union(ndx_lig).cloned().collect();
        let index_map: Vec<usize> = union_set.iter().cloned().collect();
        
        let ndx_rec_norm = ndx_rec
            .iter()
            .map(|&i| index_map.binary_search(&i).unwrap())
            .collect();
        
        let ndx_lig_norm = ndx_lig
            .iter()
            .map(|&i| index_map.binary_search(&i).unwrap())
            .collect();
        
        (ndx_rec_norm, Some(ndx_lig_norm))
    } else {
        let continuous_set: BTreeSet<usize> = (0..ndx_rec.len()).collect();
        (continuous_set, None)
    }
}

pub fn get_residues_tpr(tpr: &TPR, ndx_com: &BTreeSet<usize>) -> Vec<Residue> {
    let mut residues: Vec<Residue> = Vec::with_capacity(ndx_com.len());
    let mut idx = 0;
    let mut resind_offset = 0;
    
    for mol in &tpr.molecule_types {
        let mol_type = &tpr.molecule_blocks[mol.molecule_type_id];
        for _ in 0..mol_type.molecules_num {
            for atom in &mol.atoms {
                idx += 1;
                if ndx_com.contains(&idx) && residues.len() <= atom.resind + resind_offset {
                    let mut cur_res = mol.residues[atom.resind].to_owned();
                    // 检查重复的nr
                    if !residues.is_empty() {
                        let last_nr = residues.last().unwrap().nr;
                        if residues.iter().any(|r| r.nr == cur_res.nr) {
                            cur_res.nr = last_nr + 1;
                        }
                    }
                    
                    residues.push(cur_res);
                }
            }
            resind_offset += mol.residues.len();
        }
    }
    residues
}

fn show_grp(grp_id: Option<usize>, ndx: &Index) -> String {
    if let Some(grp_id) = grp_id {
        if let Some(grp) = ndx.groups.get(grp_id) {
            format!("{}): {}", grp_id, grp)
        } else {
            String::from("undefined")
        }
    } else {
        String::from("undefined")
    }
}

fn prepare_system_tpr(receptor_grp: usize, ligand_grp: Option<usize>, 
                  trj: &String, tpr: &mut TPR, ndx: &Index, 
                  tpr_name: &str, ndx_name: &String, bt: f64, et: f64, 
                  dt: f64, dt_ie: f64, config: &Option<Config>,
                  wd: &Path, settings: &mut Settings) {
    // atom indexes
    println!("Preparing atom indexes...");
    let ndx_lig = match ligand_grp {
        Some(ligand_grp) => Some(&ndx.groups[ligand_grp].indexes),
        None => None
    };
    let ndx_rec = &ndx.groups[receptor_grp].indexes;
    let ndx_com: BTreeSet<usize> = match ndx_lig {
        Some(ndx_lig) => {
            ndx_rec.union(ndx_lig).cloned().collect()
        }
        None => ndx_rec.iter().cloned().collect()
    };

    // atom properties
    println!("Parsing atom properties...");
    let mut aps = AtomProperties::from_tpr(tpr, &ndx_com);
    println!("Collecting residues list...");
    let residues = get_residues_tpr(tpr, &ndx_com);

    if trj.ends_with("pdb") {
        let pdb = PDB::from(trj);
        let mut coordinates: Array3<f64> = Array3::zeros((pdb.models.len(), pdb.models[0].atoms.len(), 3));
        for (i, model) in pdb.models.iter().enumerate() {
            for (j, atom) in model.atoms.iter().enumerate() {
                coordinates[[i, j, 0]] = atom.x;
                coordinates[[i, j, 1]] = atom.y;
                coordinates[[i, j, 2]] = atom.z;
            }
        }
        let time_list = (0..coordinates.shape()[0]).map(|t| (t + 1) as f64 * 1000.0).collect();

        println!("Normalizing index...");
        let ndx_lig = match ligand_grp {
            Some(ligand_grp) => Some(ndx.groups[ligand_grp].indexes.clone()),
            None => None
        };
        let (ndx_rec, ndx_lig) = 
            normalize_index(&ndx.groups[receptor_grp].indexes, &ndx_lig);
        
        set_para_mmpbsa(&time_list, &time_list, &coordinates, tpr, &ndx, config,
            wd, &mut aps, &ndx_rec, &ndx_lig, receptor_grp, ligand_grp, &residues, settings);
    } else if trj.ends_with("gro") {
        let gro = GRO::from(trj);
        let mut coordinates: Array3<f64> = Array3::zeros((1, gro.atoms.len(), 3));
        for (i, atom) in gro.atoms.iter().enumerate() {
            coordinates[[0, i, 0]] = atom.x * 10.0;
            coordinates[[0, i, 1]] = atom.y * 10.0;
            coordinates[[0, i, 2]] = atom.z * 10.0;
        }
        let time_list = (0..coordinates.shape()[0]).map(|t| (t + 1) as f64 * 1000.0).collect();
        
        println!("Normalizing index...");
        let ndx_lig = match ligand_grp {
            Some(ligand_grp) => Some(ndx.groups[ligand_grp].indexes.clone()),
            None => None
        };
        let (ndx_rec, ndx_lig) = 
            normalize_index(&ndx.groups[receptor_grp].indexes, &ndx_lig);
        
        set_para_mmpbsa(&time_list, &time_list, &coordinates, tpr, &ndx, config,
            wd, &mut aps, &ndx_rec, &ndx_lig, receptor_grp, ligand_grp, &residues, settings);
    } else {
        // pre-treat trajectory
        let trj_mmpbsa = append_new_name(trj, ".xtc", "_MMPBSA_"); // get trj output file name
        let tpr_name = append_new_name(tpr_name, ".tpr", ""); // fuck the passed tpr name is dump
        
        // step 1: generate new index
        println!("Generating Index...");
        let ndx_whole = if let Some(ndx_lig) = &ndx_lig {
            Index::new(vec![
                IndexGroup::new("Complex", &ndx_rec.union(&ndx_lig).cloned().collect()), 
                IndexGroup::new("Receptor", &ndx_rec),
                IndexGroup::new("Ligand", &ndx_lig)
            ])
        } else {
            Index::new(vec![
                IndexGroup::new("Complex", &ndx_rec)
            ])
        };
        let ndx_mmpbsa = append_new_name(ndx_name, ".ndx", "_MMPBSA_");
        ndx_whole.to_ndx(&ndx_mmpbsa);
        
        // step 2: extract new trj with old tpr and new index
        println!("Extracting trajectory, be patient...\x1b[90m");   // turn gray
        // currently use smaller dt_ie
        // let trj_fullpath = Path::new(trj).to_str().unwrap();
        if settings.fix_pbc {
            trjconv(&vec!["Complex", "Complex", "Complex"], &env::current_dir().unwrap(), settings, &trj, &tpr_name, &ndx_mmpbsa, &trj_mmpbsa, 
                &["-b", &bt.to_string(), "-e", &et.to_string(), "-dt", &dt_ie.to_string(), "-pbc", "cluster", "-center"]);
        } else {
            trjconv(&vec!["Complex"], &env::current_dir().unwrap(), settings, &trj, &tpr_name, &ndx_mmpbsa, &trj_mmpbsa, 
                &["-b", &bt.to_string(), "-e", &et.to_string(), "-dt", &dt_ie.to_string()]);
        }
        
        // step 3: extract new tpr from old tpr
        let tpr_mmpbsa = append_new_name(&tpr_name, ".tpr", "_MMPBSA_"); // get extracted tpr file name
        convert_tpr(&vec!["Complex"], &env::current_dir().unwrap(), settings, &tpr_name, &ndx_mmpbsa, &tpr_mmpbsa);
        
        // step 4: generate new index with new tpr
        // must normalize index here after trajectory extracion, or the traj may contain less atoms
        println!("\x1b[0mNormalizing index...\x1b[90m");   // turn white
        let ndx_lig = match ligand_grp {
            Some(ligand_grp) => Some(ndx.groups[ligand_grp].indexes.clone()),
            None => None
        };
        let (ndx_rec, ndx_lig) = 
            normalize_index(&ndx.groups[receptor_grp].indexes, &ndx_lig);

        let ndx_whole = if let Some(ndx_lig) = &ndx_lig {
            Index::new(vec![
                IndexGroup::new("Complex", &ndx_rec.union(&ndx_lig).cloned().collect()), 
                IndexGroup::new("Receptor", &ndx_rec),
                IndexGroup::new("Ligand", &ndx_lig)
            ])
        } else {
            Index::new(vec![
                IndexGroup::new("Complex", &ndx_rec)
            ])
        };
        ndx_whole.to_ndx(&ndx_mmpbsa);
        
        // 需要处理一下atom_properties的id
        aps.atom_props.iter_mut().enumerate().for_each(|(i, ap)| ap.id = i);
        
        // 生成初始结构方便查看
        let init_struct = append_new_name(trj, "_struct.gro", "_MMPBSA_"); // get trj output file name
        trjconv(&vec!["Complex"], &env::current_dir().unwrap(), settings, &trj_mmpbsa, &tpr_name, &ndx_mmpbsa, &init_struct, &vec!["-dump", "0"]);
        
        // step 5: Read trajectory and get time and coordinate
        println!("\x1b[0mPreparing trajectories for IE calculation...");
        println!("Loading trajectory coordinates...");
        let time_box_info = read_xtc(&trj_mmpbsa);
        let (time_list_ie, coordinates_ie): (Vec<f64>, Vec<Vec<[f32; 3]>>) = time_box_info
            .par_iter()
            .map(|frame| (frame.0 as f64, frame.1.clone()))
            .unzip();

        let num_frames = time_box_info.len();
        let num_atoms = if num_frames > 0 { time_box_info[0].1.len() } else { 0 };
        let coordinates_ie = Array3::from_shape_fn((num_frames, num_atoms, 3), |(i, j, k)| {
            coordinates_ie[i][j][k] as f64 * 10.0
        });

        let start_time = time_list_ie.first().copied().unwrap_or(0.0);
        let end_time = time_list_ie.last().copied().unwrap_or(0.0);
        let num_time_points = ((end_time - start_time) / dt + 1.0) as usize;
        let time_list = Array1::linspace(start_time, end_time, num_time_points).to_vec();

        set_para_mmpbsa(&time_list, &time_list_ie, &coordinates_ie, tpr, &ndx, config,
            wd, &mut aps, &ndx_rec, &ndx_lig, receptor_grp, ligand_grp, &residues, settings);
    };
}
