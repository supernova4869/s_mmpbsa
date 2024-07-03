use std::cmp::Ordering;
use std::collections::HashSet;
use std::path::Path;
use std::rc::Rc;
use xdrfile::Frame;
use std::fs;

use crate::parse_pdbqt::PDBQT;
use crate::settings::Settings;
use crate::utils::{get_input_selection, append_new_name};
use crate::fun_para_mmpbsa::{set_para_mmpbsa, set_para_mmpbsa_pdbqt};
use crate::index_parser::{Index, IndexGroup};
use crate::parse_tpr::TPR;
use crate::atom_property::AtomProperty;
use crate::parse_tpr::Residue;
use indicatif::{ProgressBar, ProgressStyle};
use crate::utils::{convert_tpr, trjconv};

pub fn set_para_trj(trj: &String, tpr: &mut TPR, ndx_name: &String, wd: &Path, tpr_name: &str, settings: &mut Settings) {
    let mut receptor_grp: Option<usize> = None;
    let mut ligand_grp: Option<usize> = None;
    let mut bt: f64 = 0.0;                                  // ps
    let mut et: f64 = tpr.dt * tpr.nsteps as f64;           // ps
    let mut dt: f64 = tpr.dt * tpr.nstxout as f64;          // ps
    let unit_dt: f64 = tpr.dt * tpr.nstxout as f64;         // ps
    let ndx = Index::from(ndx_name);
    loop {
        println!("\n                 ************ Trajectory Parameters ************");
        println!("-10 Return");
        println!(" -1 Toggle whether to fix PBC conditions, current: {}", settings.fix_pbc);
        println!("  0 Go to next step");
        println!("  1 Select receptor groups, current:          {}", show_grp(receptor_grp, &ndx));
        println!("  2 Select ligand groups, current:            {}", show_grp(ligand_grp, &ndx));
        println!("  3 Set start time to analyze, current:       {} ns", bt / 1000.0);
        println!("  4 Set end time to analyze, current:         {} ns", et / 1000.0);
        println!("  5 Set time interval to analyze, current:    {} ns", dt / 1000.0);
        let i = get_input_selection();
        match i {
            -10 => return,
            -1 => {
                settings.fix_pbc = !settings.fix_pbc;
            }
            0 => {
                match receptor_grp {
                    Some(receptor_grp) => {
                        // atom indexes
                        println!("Preparing atom indexes...");
                        let ndx_lig = match ligand_grp {
                            Some(ligand_grp) => Some(&ndx.groups[ligand_grp].indexes),
                            None => None
                        };
                        let ndx_rec = &ndx.groups[receptor_grp].indexes;
                        let ndx_com = match ndx_lig {
                            Some(ndx_lig) => {
                                match ndx_lig[0] > ndx_rec[0] {
                                    true => {
                                        let mut ndx_com = ndx_rec.to_vec();
                                        ndx_com.extend(ndx_lig);
                                        ndx_com
                                    }
                                    false => {
                                        let mut ndx_com = ndx_lig.to_vec();
                                        ndx_com.extend(ndx_rec);
                                        ndx_com
                                    }
                                }
                            }
                            None => ndx_rec.to_vec()
                        };

                        // atom properties
                        println!("Parsing atom properties...");
                        let mut aps = AtomProperty::from_tpr(tpr, &ndx_com);
                        println!("Collecting residues list...");
                        let residues = get_residues_tpr(tpr, &ndx_com);

                        // pre-treat trajectory: fix pbc
                        println!("Extracting trajectory...");

                        let trj_whole = append_new_name(trj, "_1_whole.xtc", "_MMPBSA_"); // get trj output file name
                        let trj_center = append_new_name(trj, "_2_center.xtc", "_MMPBSA_");
                        let trj_cluster = append_new_name(trj, "_3_cluster.xtc", "_MMPBSA_");
                        let trj_mmpbsa = append_new_name(trj, "_4_pbc.xtc", "_MMPBSA_");
                        let tpr_name = append_new_name(tpr_name, ".tpr", "");       // fuck the tpr name is dump
                        
                        // add a Complex group to index file
                        let com_group = IndexGroup::new("Complex", &ndx_com);
                        let mut new_ndx = ndx.clone();
                        new_ndx.rm_group("Complex");
                        new_ndx.push(&com_group);
                        let ndx_whole = append_new_name(ndx_name, "_whole.ndx", "_MMPBSA_"); // get extracted index file name
                        new_ndx.to_ndx(&ndx_whole);
                        
                        // echo "Complex" | gmx trjconv -f md.xtc -s md.tpr -n index.idx -o md_trj_whole.xtc -pbc whole
                        trjconv("Complex", wd, settings, trj, &tpr_name, &ndx_whole, &trj_whole, &["-pbc", "whole"], settings.debug_mode);
                        // echo "Complex" | gmx convert-tpr -s md.tpr -n index.idx -o md_trj_com.tpr
                        let tpr_mmpbsa = append_new_name(&tpr_name, ".tpr", "_MMPBSA_"); // get extracted tpr file name
                        convert_tpr("Complex", wd, settings, &tpr_name, &ndx_whole, &tpr_mmpbsa, settings.debug_mode);
                        if !settings.debug_mode {
                            fs::remove_file(&ndx_whole).unwrap();
                        }

                        // Index normalization
                        let (ndx_com, ndx_rec, ndx_lig) = 
                            normalize_index(&ndx.groups[receptor_grp].indexes, match ligand_grp {
                                Some(ligand_grp) => Some(&ndx.groups[ligand_grp].indexes),
                                None => None
                            });

                        // extract index file
                        let ndx_mmpbsa = match ligand_grp {
                            Some(ligand_grp) => {
                                Index::new(vec![
                                    IndexGroup::new("Complex", &ndx_com), 
                                    IndexGroup::new(&ndx.groups[receptor_grp].name, &ndx_rec),
                                    IndexGroup::new(&ndx.groups[ligand_grp].name, &ndx_lig)
                                ])
                            },
                            None => {
                                Index::new(vec![IndexGroup::new("Receptor", &ndx_com)])
                            }
                        };
                        ndx_mmpbsa.to_ndx(Path::new(wd).join("_MMPBSA_index.ndx").to_str().unwrap());
                        let ndx_mmpbsa = Path::new(wd).join("_MMPBSA_index.ndx");
                        let ndx_mmpbsa = ndx_mmpbsa.to_str().unwrap();

                        let trj_mmpbsa = if settings.fix_pbc {
                            println!("Fixing PBC conditions...");
                            match ligand_grp {
                                Some(ligand_grp) => {
                                    println!("Fixing PBC 0/3...");
                                    // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $pdb    &>>$err -pbc mol -center
                                    trjconv(&(ndx.groups[ligand_grp].name.to_owned() + " Complex"),
                                        wd, settings, &trj_whole, &tpr_mmpbsa, &ndx_mmpbsa, &trj_center, &["-pbc", "mol", "-center"], settings.debug_mode);
                                    println!("Fixing PBC 1/3...");
                                    // echo -e "$com\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
                                    trjconv("Complex Complex",
                                        wd, settings, &trj_center, &tpr_mmpbsa, &ndx_mmpbsa, &trj_cluster, &["-pbc", "cluster"], settings.debug_mode);
                                    println!("Fixing PBC 2/3...");
                                    // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjcls -o $pdb    &>>$err -fit rot+trans
                                    trjconv("1 0",
                                        wd, settings, &trj_cluster, &tpr_mmpbsa, &ndx_mmpbsa, &trj_mmpbsa, &["-fit", "rot+trans"], settings.debug_mode);
                                    if !settings.debug_mode {
                                        fs::remove_file(&trj_center).unwrap();
                                        fs::remove_file(&trj_cluster).unwrap();
                                    }
                                },
                                None => {
                                    // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $trjcnt &>>$err -pbc mol -center
                                    println!("Fixing PBC 0/1...");
                                    trjconv("0 0 0", 
                                        wd, settings, &trj_whole, &tpr_mmpbsa, &ndx_mmpbsa, &trj_mmpbsa, &["-pbc", "mol", "-center", "-fit", "rot+trans"], settings.debug_mode);
                                }
                            }
                            if !settings.debug_mode {
                                fs::remove_file(&trj_whole).unwrap();
                            }
                            trj_mmpbsa
                        } else {
                            trj_whole
                        };
                        if !settings.debug_mode {
                            fs::remove_file(&tpr_mmpbsa).unwrap();
                            fs::remove_file(&ndx_mmpbsa).unwrap();
                        }
                        println!("Fixing PBC finished.");

                        set_para_mmpbsa(&trj_mmpbsa, tpr, &ndx, wd, &mut aps, 
                            &ndx_com,
                            &ndx_rec,
                            &ndx_lig,
                            receptor_grp,
                            ligand_grp,
                            bt, et, dt,
                            &residues,
                            settings);
                    }
                    _ => println!("Please select receptor groups.")
                }
            }
            1 => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input receptor group num:");
                receptor_grp = Some(get_input_selection());
            }
            2 => {
                println!("Current groups:");
                ndx.list_groups();
                println!("Input ligand group num (-1 for nothing):");
                ligand_grp = match get_input_selection() {
                    -1 => None,
                    i => Some(i as usize)
                };
            }
            3 => {
                println!("Input start time (ns), should be divisible of {} ps:", dt);
                let mut new_bt = get_input_selection::<f64>() * 1000.0;
                while new_bt * 1000.0 % dt != 0.0 || new_bt > tpr.nsteps as f64 * tpr.dt as f64 || new_bt < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_bt / 1000.0);
                    println!("Input start time (ns) again, should be divisible of {} fs:", dt);
                    new_bt = get_input_selection::<f64>() * 1000.0;
                }
                bt = new_bt;
            }
            4 => {
                println!("Input end time (ns), should be divisible of {} ps:", dt);
                let mut new_et = get_input_selection::<f64>() * 1000.0;
                while new_et * 1000.0 % dt != 0.0 || new_et > tpr.nsteps as f64 * tpr.dt as f64 || new_et < 0.0 {
                    println!("The input {} ns not a valid time in trajectory.", new_et / 1000.0);
                    println!("Input end time (ns) again, should be divisible of {} fs:", dt);
                    new_et = get_input_selection::<f64>() * 1000.0;
                }
                et = new_et;
            }
            5 => {
                println!("Input interval time (ns), should be divisible of {} ps:", unit_dt);
                let mut new_dt = get_input_selection::<f64>() * 1000.0;
                while new_dt * 1000.0 % unit_dt != 0.0 {
                    println!("The input {} ns is not a valid time step.", new_dt / 1000.0);
                    println!("Input interval time (ns) again, should be divisible of {} ps:", unit_dt);
                    new_dt = get_input_selection::<f64>() * 1000.0;
                }
                dt = new_dt;
            }
            _ => println!("Invalid input")
        }
    }
}

fn show_grp(grp: Option<usize>, ndx: &Index) -> String {
    match grp {
        None => String::from("undefined"),
        Some(grp) => format!("{}): {}, {} atoms",
                    grp,
                    ndx.groups[grp as usize].name,
                    ndx.groups[grp as usize].indexes.len())
    }
}

pub fn set_para_dock(receptor: &PDBQT, ligand: &PDBQT, wd: &Path, settings: &mut Settings) {
    let mut bf: i32 = 0;
    let mut ef: i32 = ligand.models.len() as i32 - 1;
    loop {
        println!("\n                 ************ Trajectory Parameters ************");
        println!("-10 Return");
        println!("  0 Go to next step");
        println!("  1 Set start model id to analyze, current:      {}", bf + 1);
        println!("  2 Set end model id to analyze, current:        {}", ef + 1);
        let i = get_input_selection();
        match i {
            -10 => return,
            0 => {
                // 建立frames
                let mut frames: Vec<Rc<Frame>> = vec![];
                for i in bf as usize..(ef as usize + 1) {
                    let mut frame = Frame::with_len(receptor.models[0].atoms.len() + ligand.models[0].atoms.len());
                    frame.step = i;
                    frame.time = i as f32 * 1000.0;
                    for (j, aj) in frame.coords.iter_mut().enumerate() {
                        if j < receptor.models[0].atoms.len() {
                            let cur_atom = &receptor.models[0].atoms[j];
                            aj[0] = cur_atom.x as f32 / 10.0;
                            aj[1] = cur_atom.y as f32 / 10.0;
                            aj[2] = cur_atom.z as f32 / 10.0;
                        } else {
                            let cur_atom = &ligand.models[i].atoms[j - receptor.models[0].atoms.len()];
                            aj[0] = cur_atom.x as f32 / 10.0;
                            aj[1] = cur_atom.y as f32 / 10.0;
                            aj[2] = cur_atom.z as f32 / 10.0;
                        }
                    }
                    frames.push(Rc::new(frame));
                }

                println!("Parsing atom properties...");
                let mut aps = AtomProperty::from_pdbqt(receptor, ligand);
                println!("Collecting residues list...");
                let residues = get_residues_pdbqt(receptor, ligand);
                let ndx_com = Vec::from_iter(0..(receptor.models[0].atoms.len() + ligand.models[0].atoms.len()));
                let ndx_rec = Vec::from_iter(0..receptor.models[0].atoms.len());
                let ndx_lig = Vec::from_iter(ndx_rec.len()..ndx_rec.len() + ligand.models[0].atoms.len());
                set_para_mmpbsa_pdbqt(&frames, &mut aps, bf as usize, ef as usize, 1, (ef - bf) as usize + 1, 298.15, 
                    &ndx_com, &ndx_rec, &ndx_lig, wd, &residues, settings);
            }
            1 => {
                println!("Input start model");
                let mut new_bf = get_input_selection::<i32>() - 1;
                while new_bf < 0 {
                    println!("The input {} not a valid model in ligands.", new_bf + 1);
                    println!("Input start model");
                    new_bf = get_input_selection::<i32>() - 1;
                }
                bf = new_bf as i32;
            }
            2 => {
                println!("Input end model");
                let mut new_ef = get_input_selection::<i32>() - 1;
                while new_ef >= ligand.models.len() as i32 {
                    println!("The input {} not a valid model in ligands.", new_ef + 1);
                    println!("Input end model");
                    new_ef = get_input_selection::<i32>() - 1;
                }
                ef = new_ef as i32;
            }
            _ => println!("Invalid input")
        }
    }
}

// convert rec and lig to begin at 0 and continous
pub fn normalize_index(ndx_rec: &Vec<usize>, ndx_lig: Option<&Vec<usize>>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let offset = match ndx_lig {
        Some(ndx_lig) => ndx_lig[0].min(ndx_rec[0]),
        None => ndx_rec[0]
    };
    let mut ndx_rec: Vec<usize> = ndx_rec.iter().map(|p| p - offset).collect();
    let mut ndx_lig = match ndx_lig {
        Some(ndx_lig) => ndx_lig.iter().map(|p| p - offset).collect(),
        None => ndx_rec.clone()
    };
    let ndx_com = match ndx_lig[0].cmp(&ndx_rec[0]) {
        Ordering::Greater => {
            ndx_lig = ndx_lig.iter().map(|p| p - ndx_lig[0] + ndx_rec.len()).collect();
            let mut ndx_com = ndx_rec.to_vec();
            ndx_com.extend(&ndx_lig);
            ndx_com
        }
        Ordering::Less => {
            ndx_rec = ndx_rec.iter().map(|p| p - ndx_rec[0] + ndx_lig.len()).collect();
            let mut ndx_com = ndx_lig.to_vec();
            ndx_com.extend(&ndx_rec);
            ndx_com
        }
        Ordering::Equal => Vec::from_iter(0..ndx_rec.len())
    };
    (ndx_com, ndx_rec, ndx_lig)
}

pub fn get_residues_tpr(tpr: &TPR, ndx_com: &Vec<usize>) -> Vec<Residue> {
    let mut residues: Vec<Residue> = vec![];
    let mut idx = 0;
    let mut resind_offset = 0;
    
    let pb = ProgressBar::new(tpr.n_atoms as u64);
    pb.set_style(ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:50.cyan/cyan} {percent}% {msg}").unwrap()
        .progress_chars("=>-"));
    for mol in &tpr.molecules {
        for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
            for atom in &mol.atoms {
                idx += 1;
                if ndx_com.contains(&idx) && residues.len() <= atom.resind + resind_offset {
                    residues.push(mol.residues[atom.resind].to_owned());
                }
                pb.inc(1);
                pb.set_message(format!("eta. {} s", pb.eta().as_secs()));
            }
            resind_offset += mol.residues.len();
        }
    }
    pb.finish();
    residues
}

pub fn get_residues_pdbqt(receptor: &PDBQT, ligand: &PDBQT) -> Vec<Residue> {
    let mut residues: Vec<Residue> = vec![];
    let receptor_resid: HashSet<i32> = HashSet::from_iter(receptor.models[0].atoms.iter().map(|a| a.resid));
    let mut receptor_resid: Vec<i32> = receptor_resid.into_iter().collect();
    receptor_resid.sort();
    let ligand_resid: HashSet<i32> = HashSet::from_iter(ligand.models[0].atoms.iter().map(|a| a.resid));
    let mut ligand_resid: Vec<i32> = ligand_resid.into_iter().collect();
    ligand_resid.sort();
    receptor_resid.extend(ligand_resid);
    let resid = receptor_resid;
    let mut atoms = receptor.models[0].atoms.to_vec();
    atoms.extend(ligand.models[0].atoms.to_vec());

    // residue names
    let mut resnames = Vec::new();
    let mut index = 0;
    for a in &atoms {
        if a.resid == resid[index] {
            resnames.push(a.resname.to_string());
            index += 1;
            if index == resid.len() {
                break;
            }
        }
    }
    for (index, resid) in resid.iter().enumerate() {
        residues.push(Residue::new(index, resnames[index].to_string(), *resid));
    }
    residues
}