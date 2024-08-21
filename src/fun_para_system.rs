use std::path::Path;
use std::fs;

use crate::settings::Settings;
use crate::utils::{append_new_name, get_input_selection, make_ndx, trajectory};
use crate::fun_para_mmpbsa::set_para_mmpbsa;
use crate::index_parser::{Index, IndexGroup};
use crate::parse_tpr::TPR;
use crate::atom_property::AtomProperties;
use crate::parse_tpr::Residue;
use crate::utils::{convert_tpr, trjconv};

pub fn set_para_trj(trj: &String, tpr: &mut TPR, ndx_name: &String, wd: &Path, tpr_name: &str, settings: &mut Settings) {
    let mut receptor_grp: Option<usize> = None;
    let mut ligand_grp: Option<usize> = None;
    let mut bt: f64 = 0.0;                                  // ps
    let mut et: f64 = tpr.dt * tpr.nsteps as f64;           // ps
    let mut dt = 1000.0;          // ps
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
                        let mut aps = AtomProperties::from_tpr(tpr, &ndx_com);
                        println!("Collecting residues list...");
                        let residues = get_residues_tpr(tpr, &ndx_com);

                        // pre-treat trajectory: fix pbc
                        let trj_whole = append_new_name(trj, "_1_whole.xtc", "_MMPBSA_"); // get trj output file name
                        let trj_center = append_new_name(trj, "_2_center.xtc", "_MMPBSA_");
                        let trj_cluster = append_new_name(trj, "_3_cluster.xtc", "_MMPBSA_");
                        let trj_mmpbsa = append_new_name(trj, "_4_pbc.xtc", "_MMPBSA_");
                        let tpr_name = append_new_name(tpr_name, ".tpr", "");       // fuck the tpr name is dump
                        
                        // step 1: generate new index
                        println!("Generating Index...");
                        // gmx make_ndx -f md.tpr -n index.idx -o md_trj_whole.xtc -pbc whole
                        let ndx_whole = append_new_name(ndx_name, "_whole.ndx", "_MMPBSA_"); // get extracted index file name
                        if let Some(ligand_grp) = ligand_grp {
                            make_ndx(&vec![
                                format!("{} | {}", receptor_grp, ligand_grp).as_str(),
                                format!("name {} Complex", ndx.groups.len()).as_str(),
                                format!("name {} Receptor", receptor_grp).as_str(),
                                format!("name {} Ligand", ligand_grp).as_str(),
                                "q"
                            ], wd, settings, &tpr_name, ndx_name, &ndx_whole);
                        } else {
                            make_ndx(&vec![
                                format!("{}", receptor_grp).as_str(),
                                format!("name {} Receptor", receptor_grp).as_str(),
                                "q"
                            ], wd, settings, &tpr_name, ndx_name, &ndx_whole);
                        }
                        
                        // step 2: extract new trj with old tpr
                        println!("Extracting trajectory, be patient...");
                        // echo "Complex" | gmx trjconv -f md.xtc -s md.tpr -n index.idx -o md_trj_whole.xtc -pbc whole
                        trjconv(&vec!["Complex"], wd, settings, trj, &tpr_name, &ndx_whole, &trj_whole, &["-pbc", "whole"]);
                        
                        // step 3: extract new tpr with old tpr
                        // echo "Complex" | gmx convert-tpr -s md.tpr -n index.idx -o md_trj_com.tpr
                        let tpr_mmpbsa = append_new_name(&tpr_name, ".tpr", "_MMPBSA_"); // get extracted tpr file name
                        convert_tpr(&vec!["Complex"], wd, settings, &tpr_name, &ndx_whole, &tpr_mmpbsa);
                        if !settings.debug_mode {
                            fs::remove_file(&ndx_whole).unwrap();
                        }
                        
                        // step 4: generate new index with new tpr
                        println!("Normalizing index...");
                        let (ndx_rec, ndx_lig) = 
                            normalize_index(&ndx.groups[receptor_grp].indexes, match ligand_grp {
                                Some(ligand_grp) => Some(&ndx.groups[ligand_grp].indexes),
                                None => None
                            });
                        // 需要处理一下atom_properties的id
                        aps.atom_props.iter_mut().enumerate().for_each(|(i, ap)| ap.id = i);
                        
                        // extract index file
                        let ndx_mmpbsa = match ligand_grp {
                            Some(_) => {
                                Index::new(vec![
                                    IndexGroup::new("Complex", &ndx_rec.iter().chain(ndx_lig.iter()).cloned().collect()), 
                                    IndexGroup::new("Receptor", &ndx_rec),
                                    IndexGroup::new("Ligand", &ndx_lig)
                                ])
                            },
                            None => {
                                Index::new(vec![IndexGroup::new("Receptor", &ndx_rec)])
                            }
                        };
                        ndx_mmpbsa.to_ndx(Path::new(wd).join("_MMPBSA_index.ndx").to_str().unwrap());
                        let ndx_mmpbsa = Path::new(wd).join("_MMPBSA_index.ndx");
                        let ndx_mmpbsa = ndx_mmpbsa.to_str().unwrap();

                        let trj_mmpbsa = if settings.fix_pbc {
                            match ligand_grp {
                                Some(_) => {
                                    println!("Fixing PBC conditions 0/3...");
                                    // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $pdb    &>>$err -pbc mol -center
                                    trjconv(&vec!["Ligand", "Complex"],
                                        wd, settings, &trj_whole, &tpr_mmpbsa, &ndx_mmpbsa, &trj_center, &["-pbc", "mol", "-center"]);
                                    println!("Fixing PBC conditions 1/3...");
                                    // echo -e "$com\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
                                    trjconv(&vec!["Complex", "Complex"],
                                        wd, settings, &trj_center, &tpr_mmpbsa, &ndx_mmpbsa, &trj_cluster, &["-pbc", "cluster"]);
                                    println!("Fixing PBC conditions 2/3...");
                                    // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjcls -o $pdb    &>>$err -fit rot+trans
                                    trjconv(&vec!["Ligand", "Complex"],
                                        wd, settings, &trj_cluster, &tpr_mmpbsa, &ndx_mmpbsa, &trj_mmpbsa, &["-fit", "rot+trans"]);
                                    if !settings.debug_mode {
                                        fs::remove_file(&trj_center).unwrap();
                                        fs::remove_file(&trj_cluster).unwrap();
                                    }
                                },
                                None => {
                                    // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $trjcnt &>>$err -pbc mol -center
                                    println!("Fixing PBC conditions 0/1...");
                                    trjconv(&vec!["Complex", "Complex", "Complex"], 
                                        wd, settings, &trj_whole, &tpr_mmpbsa, &ndx_mmpbsa, &trj_mmpbsa, &["-pbc", "mol", "-center", "-fit", "rot+trans"]);
                                }
                            }
                            if !settings.debug_mode {
                                fs::remove_file(&trj_whole).unwrap();
                            }
                            trj_mmpbsa
                        } else {
                            trj_whole
                        };
                        println!("Fixing PBC finished.");
                        println!("Loading trajectory file...");
                        trajectory(&vec!["Complex"], wd, settings, &trj_mmpbsa, &tpr_mmpbsa, &ndx_mmpbsa, "_MMPBSA_coord.xvg");
                        if !settings.debug_mode {
                            fs::remove_file(&tpr_mmpbsa).unwrap();
                            fs::remove_file(&ndx_mmpbsa).unwrap();
                        }

                        set_para_mmpbsa(&trj_mmpbsa, tpr, &ndx, wd, &mut aps, 
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

// convert rec and lig to begin at 0 and continous
pub fn normalize_index(ndx_rec: &Vec<usize>, ndx_lig: Option<&Vec<usize>>) -> (Vec<usize>, Vec<usize>) {
    if let Some(ndx_lig) = ndx_lig {
        let mut ndx_rec_norm = ndx_rec.clone();
        let mut ndx_lig_norm = ndx_lig.clone();
        let last_atom = ndx_rec.len() + ndx_lig.len() - 1;
        for cur_atom_id in 0..=last_atom {
            if !ndx_lig_norm.contains(&cur_atom_id) && !ndx_rec_norm.contains(&cur_atom_id) {
                let ndx_lig_norm2 = ndx_lig_norm.clone();
                let ndx_rec_norm2 = ndx_rec_norm.clone();
                let next_edge_id = ndx_lig_norm2.iter().find(|&&i| i > cur_atom_id).unwrap_or(&last_atom);
                let next_edge_id = next_edge_id.min(ndx_rec_norm2.iter().find(|&&i| i > cur_atom_id).unwrap_or(&last_atom));
                let offset = next_edge_id - cur_atom_id;
                ndx_lig_norm.iter_mut().for_each(|i| if *i > cur_atom_id { *i -= offset } );
                ndx_rec_norm.iter_mut().for_each(|i| if *i > cur_atom_id { *i -= offset } );
            }
        }
        (ndx_rec_norm, ndx_lig_norm)
    } else {
        ((0..ndx_rec.len()).collect(), (0..ndx_rec.len()).collect())
    }
}

pub fn get_residues_tpr(tpr: &TPR, ndx_com: &Vec<usize>) -> Vec<Residue> {
    let mut residues: Vec<Residue> = vec![];
    let mut idx = 0;
    let mut resind_offset = 0;
    
    for mol in &tpr.molecules {
        for _ in 0..tpr.molecule_types[mol.molecule_type_id].molecules_num {
            for atom in &mol.atoms {
                idx += 1;
                if ndx_com.contains(&idx) && residues.len() <= atom.resind + resind_offset {
                    residues.push(mol.residues[atom.resind].to_owned());
                }
            }
            resind_offset += mol.residues.len();
        }
    }
    residues
}
