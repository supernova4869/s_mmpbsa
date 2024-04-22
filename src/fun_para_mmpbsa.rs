use std::io::stdin;
use std::path::Path;
use crate::utils::{get_input, get_input_selection, append_new_name};
use crate::index_parser::{Index, IndexGroup};
use crate::settings::Settings;
use crate::apbs_param::{PBASet, PBESet};
use std::io::Write;
use std::fs::{File, self};
use std::cmp::Ordering;
use crate::atom_property::AtomProperty;
use crate::parse_tpr::TPR;
use crate::mmpbsa::{self, get_residues};
use crate::analyzation;
use crate::utils::{convert_tpr, trjconv};

pub fn set_para_mmpbsa(trj: &String, tpr: &mut TPR, ndx: &Index, wd: &Path,
                       tpr_name: &str, ndx_name: &str,
                       receptor_grp: usize, ligand_grp: Option<usize>,
                       bt: f64, et: f64, dt: f64, settings: &mut Settings) {
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
    let mut aps = AtomProperty::new(tpr, &ndx_com);
    println!("Collecting residues list...");
    let residues = get_residues(tpr, &ndx_com);

    // pre-treat trajectory: fix pbc
    println!("Extracting trajectory...");

    let trj_whole = append_new_name(trj, "_whole.xtc", "MMPBSA_"); // get trj output file name
    let trj_center = append_new_name(trj, "_center.xtc", "MMPBSA_");
    let trj_cluster = append_new_name(trj, "_cluster.xtc", "MMPBSA_");
    let trj_mmpbsa = append_new_name(trj, ".xtc", "MMPBSA_");
    let tpr_name = append_new_name(tpr_name, ".tpr", "");       // fuck the tpr name is dump
    
    // add a Complex group to index file
    let com_group = IndexGroup::new("Complex", &ndx_com);
    let mut new_ndx = ndx.clone();
    new_ndx.rm_group("Complex");
    new_ndx.push(&com_group);
    let ndx_whole = append_new_name(ndx_name, "_whole.ndx", "MMPBSA_"); // get extracted index file name
    new_ndx.to_ndx(&ndx_whole);
    
    // echo "Complex" | gmx trjconv -f md.xtc -s md.tpr -n index.idx -o md_trj_whole.xtc -pbc whole
    trjconv("Complex", wd, settings, trj, &tpr_name, &ndx_whole, &trj_whole, &["-pbc", "whole"], settings.debug_mode);
    // echo "Complex" | gmx convert-tpr -s md.tpr -n index.idx -o md_trj_com.tpr
    let tpr_mmpbsa = append_new_name(&tpr_name, ".tpr", "MMPBSA_"); // get extracted tpr file name
    convert_tpr("Complex", wd, settings, &tpr_name, &ndx_whole, &tpr_mmpbsa, settings.debug_mode);
    if !settings.debug_mode {
        fs::remove_file(&ndx_whole).unwrap();
    }

    // Index normalization
    let (ndx_com, ndx_rec, ndx_lig) = normalize_index(ndx_rec, ndx_lig);

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
    ndx_mmpbsa.to_ndx(Path::new(wd).join("MMPBSA_index.ndx").to_str().unwrap());
    let ndx_mmpbsa = Path::new(wd).join("MMPBSA_index.ndx");
    let ndx_mmpbsa = ndx_mmpbsa.to_str().unwrap();

    let trj_mmpbsa = if settings.fix_pbc {
        println!("Fixing PBC conditions...");
        match ligand_grp {
            Some(ligand_grp) => {
                // echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $pdb    &>>$err -pbc mol -center
                trjconv(&(ndx.groups[ligand_grp].name.to_owned() + " Complex"),
                    wd, settings, &trj_whole, &tpr_mmpbsa, &ndx_mmpbsa, &trj_center, &["-pbc", "mol", "-center"], settings.debug_mode);
                // echo -e "$com\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
                trjconv("Complex Complex",
                    wd, settings, &trj_center, &tpr_mmpbsa, &ndx_mmpbsa, &trj_cluster, &["-pbc", "cluster"], settings.debug_mode);
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
    
    // kinds of radius types
    let radius_types = vec!["ff", "amber", "Bondi", "mBondi", "mBondi2"];
    let mut pbe_set = PBESet::new(tpr.temp);
    let mut pba_set = PBASet::new(tpr.temp);
    loop {
        println!("\n                 ************ MM/PB-SA Parameters ************");
        println!("-10 Return");
        println!(" -3 Output PBSA parameters");
        println!(" -2 Output LJ parameters");
        println!(" -1 Output structural parameters");
        println!("  0 Start MM/PB-SA calculation");
        println!("  1 Toggle whether to use Debye-Huckel shielding method, current: {}", settings.use_dh);
        println!("  2 Toggle whether to use interaction entropy (IE) method, current: {}", settings.use_ts);
        println!("  3 Select atom radius type, current: {}", radius_types[settings.rad_type]);
        println!("  4 Input atom distance cutoff for MM calculation (A), current: {}", settings.r_cutoff);
        println!("  5 Input coarse grid expand factor (cfac), current: {}", settings.cfac);
        println!("  6 Input fine grid expand amount (fadd), current: {} A", settings.fadd);
        println!("  7 Input fine mesh spacing (df), current: {} A", settings.df);
        println!("  8 Prepare PB parameters for APBS");
        println!("  9 Prepare SA parameters for APBS");
        println!(" 10 Toggle whether to do alanine scanning, current: {}", settings.if_alanine_scanning);
        let i = get_input_selection();
        match i {
            -10 => return,
            -1 => {
                let mut paras = File::create(wd.join("paras_structure.txt")).unwrap();
                paras.write_all(format!("Receptor group: {}\n", 
                    ndx.groups[receptor_grp as usize].name).as_bytes()).unwrap();
                match ligand_grp {
                    Some(ligand_grp) => {
                        paras.write_all(format!("Ligand group: {}\n", 
                            ndx.groups[ligand_grp as usize].name).as_bytes()).unwrap();
                    }
                    None => {
                        paras.write_all("Ligand group: None\n".as_bytes()).unwrap();
                    }
                }
                paras.write_all(format!("Atom radius type: {}\n", radius_types[settings.rad_type]).as_bytes()).unwrap();
                paras.write_all(format!("Atoms:\n     id   name   type   charge   radius   resnum  resname\n").as_bytes()).unwrap();
                for idx in 0..ndx_com.len() {
                    paras.write_all(format!("{:7}{:>7}{:7}{:9.2}{:9.2}{:9}{:>9}\n", 
                        aps.atm_index[idx], aps.atm_name[idx], aps.atm_typeindex[idx], 
                        aps.atm_charge[idx], aps.atm_radius[idx], aps.atm_resid[idx] + 1, 
                        aps.atm_resname[idx]).as_bytes()).unwrap();
                }
                println!("Structural parameters have been written to paras_structure.txt");
            }
            -2 => {
                let mut paras = File::create(wd.join("paras_LJ.txt")).unwrap();
                paras.write_all(format!("Atom types num: {}\n", tpr.atom_types_num).as_bytes()).unwrap();
                paras.write_all("c6:\n".as_bytes()).unwrap();
                for i in 0..aps.c6.shape()[0] {
                    for j in 0..aps.c6.shape()[1] {
                        paras.write_all(format!("{:13.6E} ", aps.c6[[i, j]]).as_bytes()).unwrap();
                    }
                    paras.write_all("\n".as_bytes()).unwrap();
                }
                paras.write_all("c12:\n".as_bytes()).unwrap();
                for i in 0..aps.c12.shape()[0] {
                    for j in 0..aps.c12.shape()[1] {
                        paras.write_all(format!("{:13.6E} ", aps.c12[[i, j]]).as_bytes()).unwrap();
                    }
                    paras.write_all("\n".as_bytes()).unwrap();
                }
                println!("Forcefield parameters have been written to paras_LJ.txt");
            }
            -3 => {
                let mut paras = File::create(wd.join("paras_pbsa.txt")).unwrap();
                paras.write_all(format!("Use Debye-Huckel shielding method: {}\n", settings.use_dh).as_bytes()).unwrap();
                paras.write_all(format!("Use entropy contribution: {}\n", settings.use_ts).as_bytes()).unwrap();
                paras.write_all(format!("Atom radius type: {}\n", radius_types[settings.rad_type]).as_bytes()).unwrap();
                paras.write_all(format!("Atom distance cutoff for MM calculation (A): {}\n", settings.r_cutoff).as_bytes()).unwrap();
                paras.write_all(format!("Coarse grid expand factor (cfac): {}\n", settings.cfac).as_bytes()).unwrap();
                paras.write_all(format!("Fine grid expand amount (fadd): {} A\n", settings.fadd).as_bytes()).unwrap();
                paras.write_all(format!("Fine mesh spacing (df): {} A\n\n", settings.df).as_bytes()).unwrap();
                paras.write_all(format!("PB settings:\n{}\n\n", pbe_set).as_bytes()).unwrap();
                paras.write_all(format!("SA settings:\n{}\n", pba_set).as_bytes()).unwrap();
                println!("PBSA parameters have been written to paras_pbsa.txt");
            }
            0 => {
                println!("Applying {} radius...", radius_types[settings.rad_type]);
                aps.apply_radius(settings.rad_type, tpr, ndx_com.len(), &radius_types);

                // Temp directory for PBSA
                let mut sys_name = String::from("_system");
                println!("Input system name (default: {}):", sys_name);
                let mut input = String::new();
                stdin().read_line(&mut input).expect("Error input");
                if input.trim().len() != 0 {
                    sys_name = input.trim().to_string();
                }
                let temp_dir = wd.join(&sys_name);
                if let Some(_) = settings.apbs.as_ref() {
                    println!("Temporary files will be placed at {}/", temp_dir.display());
                    if !temp_dir.is_dir() {
                        fs::create_dir(&temp_dir).expect(format!("Failed to create temp directory: {}.", &sys_name).as_str());
                    } else {
                        println!("Directory {}/ not empty. Clear? [Y/n]", temp_dir.display());
                        let mut input = String::from("");
                        stdin().read_line(&mut input).expect("Get input error");
                        if input.trim().len() == 0 || input.trim() == "Y" || input.trim() == "y" {
                            fs::remove_dir_all(&temp_dir).expect("Remove dir failed");
                            fs::create_dir(&temp_dir).expect(format!("Failed to create temp directory: {}.", &sys_name).as_str());
                        }
                    }
                } else {
                    println!("Note: Since APBS not found, solvation energy will not be calculated.");
                };
                let results = mmpbsa::fun_mmpbsa_calculations(&trj_mmpbsa, &temp_dir, &sys_name, &aps,
                                                                &ndx_com, &ndx_rec, &ndx_lig, &residues,
                                                                bt, et, dt, &pbe_set, &pba_set, settings);
                // Clean trj
                if !settings.debug_mode {
                    fs::remove_file(&trj_mmpbsa).unwrap();
                }
                analyzation::analyze_controller(&results, pbe_set.temp, &sys_name, wd, ndx_com.len(), settings);
            }
            1 => {
                settings.use_dh = !settings.use_dh;
            }
            2 => {
                settings.use_ts = !settings.use_ts;
            }
            3 => {
                println!("Input atom radius type (default mBondi), Supported:{}", {
                    let mut s = String::new();
                    for (k, v) in radius_types.iter().enumerate() {
                        s.push_str(format!("\n{}):\t{}", k, v).as_str());
                    }
                    s
                });
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.rad_type = 3;
                } else {
                    let s = s.trim().parse().expect("Input not valid number.");
                    if s == 0 {
                        settings.rad_type = 0;
                    } else if s < radius_types.len() {
                        settings.rad_type = s;
                    } else {
                        println!("Radius type {} not supported. Will use mBondi instead.", radius_types[s]);
                        settings.rad_type = 3;
                    }
                }
            }
            4 => {
                println!("Input cutoff value (A), default 0 (inf):");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.r_cutoff = f64::INFINITY;
                } else {
                    settings.r_cutoff = s.trim().parse().expect("Input not valid number.");
                    if settings.r_cutoff == 0.0 {
                        settings.r_cutoff = f64::INFINITY;
                    }
                }
            }
            5 => {
                println!("Input coarse grid expand factor, default 3:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.cfac = 3.0;
                } else {
                    settings.cfac = s.trim().parse().expect("Input not valid number.");
                }
            }
            6 => {
                println!("Input fine grid expand amount (A), default 10:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.fadd = 10.0;
                } else {
                    settings.fadd = s.trim().parse().expect("Input not valid number.");
                }
            }
            7 => {
                println!("Input fine mesh spacing (A), default 0.5:");
                let mut s = String::new();
                stdin().read_line(&mut s).expect("Input error");
                if s.trim().is_empty() {
                    settings.df = 0.5;
                } else {
                    settings.df = s.trim().parse().expect("Input not valid number.");
                }
            }
            8 => {
                let pb_fpath = wd.join("PB_settings.yaml");
                pbe_set.save_params(&pb_fpath);
                println!("PB parameters have been wrote to {0}.\n\
                    Edit it and input its path to reload (default: {0}).", &pb_fpath.to_str().unwrap());
                let pb_fpath = get_input(pb_fpath.to_str().unwrap().to_string());
                pbe_set = PBESet::load_params(pb_fpath);
            }
            9 => {
                let sa_fpath = wd.join("SA_settings.yaml");
                pba_set.save_params(&sa_fpath);
                println!("SA parameters have been wrote to {0}.\n\
                    Edit it and input its path to reload (default: {0}).", &sa_fpath.to_str().unwrap());
                let sa_fpath = get_input(sa_fpath.to_str().unwrap().to_string());
                pba_set = PBASet::load_params(sa_fpath);
            }
            10 => {
                println!("We will proceed with the alanine scanning proposal\nput forward by the Chinese representative.");
                settings.if_alanine_scanning = !settings.if_alanine_scanning;
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
