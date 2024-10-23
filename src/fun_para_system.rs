use std::process::{exit, Command, Stdio};use std::env::{self, current_exe};
use std::path::Path;
use std::fs::{self, File};
use std::io::Write;

use crate::{dump_tpr, parse_mol2::MOL2};
use crate::parse_pdb::{PDBModel, PDB};
use crate::settings::Settings;
use crate::utils::{append_new_name, get_input_selection, make_ndx, multiwfn, sobtop, trajectory};
use crate::fun_para_mmpbsa::set_para_mmpbsa;
use crate::index_parser::{Index, IndexGroup};
use crate::parse_tpr::TPR;
use crate::atom_property::AtomProperties;
use crate::parse_tpr::Residue;
use crate::utils::{convert_tpr, convert_trj, trjconv, pdb2gmx, grompp};
use crate::parse_xvg::read_coord_xvg;

pub fn set_para_trj(trj: &String, tpr: &mut TPR, ndx_name: &String, wd: &Path, tpr_name: &str, settings: &mut Settings) {
    let mut receptor_grp: Option<usize> = None;
    let mut ligand_grp: Option<usize> = None;
    let mut bt: f64 = 0.0;                                  // ps
    let mut et: f64 = tpr.dt * tpr.nsteps as f64;           // ps
    let mut dt = 1000.0;                               // ps
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
                if let Some(receptor_grp) = receptor_grp {
                    prepare_system_tpr(receptor_grp, ligand_grp, trj, tpr, &ndx, tpr_name, ndx_name, bt, et, dt, wd, settings);
                } else {
                    println!("Please select receptor groups.");
                };
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

fn prepare_pymol_complex_pdb(rec_name: &str, lig_name: &str, temp_dir: &Path) -> (PDB, usize, usize) {
    println!("Preparing complex structures...");
    let rec_path = temp_dir.join(format!("MMPBSA_docking_{}_addH.pdb", rec_name));
    let lig_path = temp_dir.join(format!("MMPBSA_docking_{}.pdb", lig_name));
    let rec_pdb = PDB::from(rec_path.to_str().unwrap());
    let lig_pdb = PDB::from(lig_path.to_str().unwrap());
    let rec_atoms_num = rec_pdb.models[0].atoms.len();
    let lig_atoms_num = lig_pdb.models[0].atoms.len();
    let mut pdb: Vec<PDBModel> = vec![];
    for (i, m) in lig_pdb.models.iter().enumerate() {
        let mut rec = rec_pdb.models.get(i).unwrap_or(&rec_pdb.models[0]).clone();
        rec.push_atoms(&m.atoms);
        pdb.push(rec);
    }
    (PDB::new(&pdb), rec_atoms_num, lig_atoms_num)
}

pub fn set_para_trj_pdbqt(receptor_path: &String, ligand_path: &String, flex_path: &Option<String>,
                          ff: &String, method: &String, basis: &String, 
                          total_charge: i32, multiplicity: usize, 
                          wd: &Path, settings: &mut Settings) {
    let receptor_file_path = Path::new(receptor_path);
    let rec_name = receptor_file_path.file_stem().unwrap().to_str().unwrap();
    let ligand_file_path = Path::new(ligand_path);
    let lig_name = ligand_file_path.file_stem().unwrap().to_str().unwrap();
    let flex_name = if let Some(flex_path) = flex_path {
        let flex_file_path = Path::new(flex_path);
        Some(flex_file_path.file_stem().unwrap().to_str().unwrap())
    } else {
        None
    };
    let temp_dir = wd.join(format!("{}_{}", rec_name, lig_name));
    let temp_dir = Path::new(&temp_dir);
    if !temp_dir.is_dir() {
        fs::create_dir(temp_dir).unwrap();
    }
    
    // prepare pdbqt files
    pdbqt2pdb(rec_name, lig_name, &flex_name, temp_dir, settings);

    // fake tpr
    prepare_system_tpr_pdb(rec_name, lig_name, &flex_name, ff, method, basis, total_charge, multiplicity, temp_dir, settings);
    dump_tpr(&wd.join("md.tpr").display().to_string(), 
        &wd.join("md.dump").display().to_string(), 
        settings.gmx_path.as_ref().unwrap());
    let tpr = TPR::from(wd.join("md.dump").to_str().unwrap(), settings);

    // fake trj
    let (pdb, rec_atoms_num, lig_atoms_num) = prepare_pymol_complex_pdb(rec_name, lig_name, temp_dir);
    println!("Preparing docking multi comformations...");
    let trj_path = format!("MMPBSA_{}_{}.pdb", rec_name, lig_name);
    let trj_path = wd.join(&trj_path);
    pdb.to_pdb(&trj_path.to_str().unwrap());
    trajectory(&vec!["0"], wd, settings, &trj_path.to_str().unwrap(), "md.tpr", 
        &temp_dir.join("MMPBSA_index.ndx").to_str().unwrap(), wd.join("_MMPBSA_coord.xvg").to_str().unwrap());
    let (_, coordinates) = read_coord_xvg(wd.join("_MMPBSA_coord.xvg").to_str().unwrap());
    let time_list = (0..coordinates.shape()[0]).map(|t| (t + 1) as f64 * 1000.0).collect();

    // fake ndx
    let ndx_com = (0..(rec_atoms_num + lig_atoms_num)).collect();
    let ndx_rec: Vec<usize> = (0..rec_atoms_num).collect();
    let ndx_lig: Vec<usize> = (rec_atoms_num..(rec_atoms_num + lig_atoms_num)).collect();
    let ndx = Index::new(vec![IndexGroup::new("Receptor", &ndx_rec), IndexGroup::new("Ligand", &ndx_lig)]);

    println!("Parsing atom properties...");
    let mut aps = AtomProperties::from_tpr(&tpr, &ndx_com);
    println!("Collecting residues list...");
    let residues = get_residues_tpr(&tpr, &ndx_com);

    if !settings.debug_mode {
        fs::remove_dir_all(temp_dir).unwrap();
        println!("Removed temp directory.");
    }

    let mut bt: usize = 0;
    let mut et: usize = pdb.models.len() - 1;
    loop {
        println!("\n                 ************ Trajectory Parameters ************");
        println!("-10 Return");
        println!("  0 Go to next step");
        println!("  1 Set start pose to analyze, current:       {}", bt + 1);
        println!("  2 Set end pose to analyze, current:         {}", et + 1);
        let i = get_input_selection();
        match i {
            -10 => return,
            -1 => {
                settings.fix_pbc = !settings.fix_pbc;
            }
            0 => {
                set_para_mmpbsa(&time_list, &coordinates, &tpr, &ndx, wd, 
                    &mut aps, &ndx_rec, &ndx_lig, 0, Some(1), &residues, settings);
            }
            1 => {
                println!("Input start pose, should be integer:");
                let mut new_bt = get_input_selection::<usize>() - 1;
                while new_bt > pdb.models.len() {
                    println!("The input {} not a valid pose in trajectory.", new_bt);
                    println!("Input start pose again:");
                    new_bt = get_input_selection::<usize>() - 1;
                }
                bt = new_bt;
            }
            2 => {
                println!("Input end pose, should be integer:");
                let mut new_et = get_input_selection::<usize>() - 1;
                while new_et > pdb.models.len() {
                    println!("The input {} not a valid pose in trajectory.", new_et);
                    println!("Input end pose again:");
                    new_et = get_input_selection::<usize>() - 1;
                }
                et = new_et;
            }
            _ => println!("Invalid input")
        }
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
                let next_edge_lig = ndx_lig_norm2.iter().find(|&&i| i > cur_atom_id);
                let next_edge_rec = ndx_rec_norm2.iter().find(|&&i| i > cur_atom_id);
                let offset = if next_edge_lig.is_none() {
                    if next_edge_rec.is_none() {
                        0
                    } else {
                        next_edge_rec.unwrap() - cur_atom_id
                    }
                } else {
                    if next_edge_rec.is_none() {
                        next_edge_lig.unwrap() - cur_atom_id
                    } else {
                        next_edge_lig.unwrap().min(next_edge_rec.unwrap()) - cur_atom_id
                    }
                };
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
                  tpr_name: &str, ndx_name: &String, 
                  bt: f64, et: f64, dt: f64, 
                  wd: &Path, settings: &mut Settings) {
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
    let trj_mmpbsa = append_new_name(trj, ".xtc", "_MMPBSA_"); // get trj output file name
    let tpr_name = append_new_name(tpr_name, ".tpr", ""); // fuck the passed tpr name is dump
    
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
            // complex is receptor
            format!("name {} Complex", receptor_grp).as_str(),
            "q"
        ], wd, settings, &tpr_name, ndx_name, &ndx_whole);
    }
    
    // step 2: extract new trj with old tpr and new index
    println!("Extracting trajectory, be patient...");
    let (bt, et, dt) = (bt.to_string(), et.to_string(), dt.to_string());
    let other_params = vec![
        "-b", &bt,
        "-e", &et,
        "-dt", &dt
    ];
    trjconv(&vec!["Complex"], wd, settings, trj, &tpr_name, &ndx_whole, &trj_mmpbsa, &other_params);
    
    // step 3: extract new tpr from old tpr
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
            Index::new(vec![
                IndexGroup::new("Complex", &ndx_rec)
            ])
        }
    };
    ndx_mmpbsa.to_ndx(wd.join("_MMPBSA_index.ndx").to_str().unwrap());
    let ndx_mmpbsa = wd.join("_MMPBSA_index.ndx");
    let ndx_mmpbsa = ndx_mmpbsa.to_str().unwrap();
    // 在这里 remove pbc, convert-trj有bug, 不能处理不完整蛋白, 故先trjconv再convert-trj
    if settings.fix_pbc {
        let other_params = vec!["-rmpbc", "-select", "Complex"];
        let pbc_name = append_new_name(&trj_mmpbsa, "_pbc.xtc", "");
        convert_trj(&vec![], wd, settings, &trj_mmpbsa, &tpr_name, &ndx_mmpbsa, &pbc_name, &other_params);
        println!("Loading trajectory coordinates...");
        trajectory(&vec!["Complex"], wd, settings, &pbc_name, &tpr_mmpbsa, &ndx_mmpbsa, "_MMPBSA_coord.xvg");
    } else {
        println!("Loading trajectory coordinates...");
        trajectory(&vec!["Complex"], wd, settings, &trj_mmpbsa, &tpr_mmpbsa, &ndx_mmpbsa, "_MMPBSA_coord.xvg");
    }

    let (time_list, coordinates) = read_coord_xvg(wd.join("_MMPBSA_coord.xvg").to_str().unwrap());

    set_para_mmpbsa(&time_list, &coordinates, tpr, &ndx, wd, &mut aps, &ndx_rec, &ndx_lig, receptor_grp, ligand_grp, &residues, settings);
}

fn pdbqt2pdb(rec_name: &str, lig_name: &str, flex_name: &Option<&str>, temp_dir: &Path, settings: &Settings) {
    let out_rec_name = append_new_name(rec_name, ".pdb", "MMPBSA_docking_");
    let out_lig_name = append_new_name(lig_name, ".pdb", "MMPBSA_docking_");
    let out_flex_name = if let Some(flex_name) = flex_name {
        Some(append_new_name(flex_name, ".pdb", "MMPBSA_docking_"))
    } else {
        None
    };
    let pml_path = temp_dir.join("MMPBSA_docking.pml");
    let mut pml_file = fs::File::create(&pml_path).unwrap();
    writeln!(pml_file, "cmd.load(r\"../{}.pdbqt\", \"Protein\")", rec_name).unwrap();
    writeln!(pml_file, "cmd.save(r\"{}\", selection=\"(Protein)\", state=1)", &out_rec_name).unwrap();
    if !out_flex_name.is_none() {
        writeln!(pml_file, "cmd.load(r\"../{}.pdbqt\", \"Flex\")", flex_name.unwrap()).unwrap();
        writeln!(pml_file, "cmd.save(r\"{}\", selection=\"(Flex)\", state=0)", out_flex_name.unwrap()).unwrap();
    }
    writeln!(pml_file, "cmd.load(r\"../{}.pdbqt\", \"Ligand\")", lig_name).unwrap();
    writeln!(pml_file, "cmd.h_add(\"all\")").unwrap();
    writeln!(pml_file, "cmd.save(r\"LIG.mol2\", selection=\"(Ligand)\", state=1)").unwrap();
    writeln!(pml_file, "cmd.save(r\"{}\", selection=\"(Ligand)\", state=0)", &out_lig_name).unwrap();
    writeln!(pml_file, "quit").unwrap();
    println!("\nLoading docking results files with PyMOL...");
    let result = Command::new(settings.pymol_path.as_ref().unwrap())
        .args(vec!["-cq", pml_path.to_str().unwrap()])
        .current_dir(temp_dir)
        .stdout(Stdio::null())
        .spawn();
    match result {
        Ok(mut child) => {
            child.wait().ok();
        }
        Err(_) => {
            eprintln!("The configured PyMOL '{}' not found.", settings.pymol_path.as_ref().unwrap());
        }
    }
}

fn copy_dir(src: &Path, dest: &Path) {
    if !dest.is_dir() {
        fs::create_dir_all(dest).unwrap();
    }

    // 读取源目录中的所有条目
    for entry in fs::read_dir(src).unwrap() {
        let entry = entry.unwrap();
        let src_path = entry.path();
        let dest_path = dest.join(entry.file_name());

        if src_path.is_dir() {
            // 递归复制子目录
            copy_dir(&src_path, &dest_path);
        } else {
            // 复制文件
            fs::copy(&src_path, &dest_path).unwrap();
        }
    }
}

fn prepare_system_tpr_pdb(rec_name: &str, lig_name: &str, flex_name: &Option<&str>, ff: &String, method: &String, basis: &String, 
                          total_charge: i32, multiplicity: usize, temp_dir: &Path, settings: &Settings) {
    // prepare protein top
    let protein_name = format!("MMPBSA_docking_{}.pdb", rec_name);
    let ff_dir = env::current_exe().unwrap().parent().unwrap().join("include").join(ff.to_string() + &".ff/");
    let dest = temp_dir.join(ff.to_string() + &".ff/");
    fs::copy(env::current_exe().unwrap().parent().unwrap().join("include").join("residuetypes.dat"), 
        temp_dir.join("residuetypes.dat")).unwrap();
    fs::copy(env::current_exe().unwrap().parent().unwrap().join("include").join("elements.dat"), 
        temp_dir.join("elements.dat")).unwrap();
    fs::copy(env::current_exe().unwrap().parent().unwrap().join("include").join("xlateat.dat"), 
        temp_dir.join("xlateat.dat")).unwrap();
    fs::copy(env::current_exe().unwrap().parent().unwrap().join("include").join("specbond.dat"), 
        temp_dir.join("specbond.dat")).unwrap();
    println!("Copying {} to {}...", ff_dir.display(), dest.display());
    copy_dir(&ff_dir, &dest);

    let protein_out = if let Some(flex_name) = flex_name {
        // prepare protein
        let flex_name = format!("MMPBSA_docking_{}.pdb", flex_name);
        let new_pdb = combine_flex(&protein_name, &flex_name, temp_dir);
        // fuck pdb2gmx cannot prepare all models
        let mut total_pdb = vec![];
        for i in 0..new_pdb.models.len() {
            let complete_protein_name = format!("MMPBSA_docking_{}_{}.pdb", rec_name, i);
            let complete_protein_name = temp_dir.join(&complete_protein_name);
            let complete_protein_name = complete_protein_name.to_str().unwrap();
            new_pdb.models[i].to_pdb(complete_protein_name);
            pdb2gmx(&vec![], temp_dir, settings, &complete_protein_name, &complete_protein_name, ff, "spc");
            let mut complete_mdl = PDBModel::from(&fs::read_to_string(complete_protein_name).unwrap());
            complete_mdl.modelid = i as i32 + 1;
            total_pdb.push(complete_mdl);
            fs::remove_file(complete_protein_name).unwrap();
        }
        let total_pdb = PDB::new(&total_pdb);
        let protein_out_pdb = append_new_name(&protein_name, "_addH.pdb", "");
        println!("Preparing flexible residues...");
        total_pdb.to_pdb(temp_dir.join(&protein_out_pdb).to_str().unwrap());
        let protein_out = append_new_name(&protein_name, ".gro", "");
        pdb2gmx(&vec![], temp_dir, settings, temp_dir.join(&protein_out_pdb).to_str().unwrap(), &protein_out, ff, "spc");
        protein_out
    } else {
        let protein_out = append_new_name(&protein_name, ".gro", "");
        pdb2gmx(&vec![], temp_dir, settings, &protein_name, &protein_out, ff, "spc");
        let protein_out_pdb = append_new_name(&protein_name, "_addH.pdb", "");
        pdb2gmx(&vec![], temp_dir, settings, &protein_name, &protein_out_pdb, ff, "spc");
        protein_out
    };

    println!("Calculating ligand charge, be patient...");
    let ligand_name = "LIG.mol2";
    let ligand_path = temp_dir.join(&ligand_name);
    let ligand_path = ligand_path.to_str().unwrap().trim_start_matches(r"\\?\");
    calc_charge(lig_name, temp_dir, method, basis, total_charge, multiplicity, settings);

    println!("Preparing docking parameters...");
    // prepare ligand top
    let lig_gro_path = temp_dir.join(append_new_name(&ligand_name, ".gro", ""));
    let lig_gro_path = lig_gro_path.to_str().unwrap();
    println!("{}", lig_gro_path);
    let itp_path = temp_dir.join(append_new_name(&ligand_name, ".itp", ""));
    let itp_path = itp_path.to_str().unwrap();
    let top_path = temp_dir.join(append_new_name(&ligand_name, ".top", ""));
    let top_path = top_path.to_str().unwrap();
    sobtop(&vec!["7", "10", temp_dir.join("LIG.chg").to_str().unwrap(), "0", 
        "2", lig_gro_path, "1", "2", "4", top_path, itp_path, "0"], settings, ligand_path);

    // include ligand top into protein
    let protein_top = temp_dir.join("topol.top").display().to_string();
    let topol = fs::read_to_string(protein_top).unwrap();
    let mut top_contents: Vec<&str> = topol.split("\n").collect();
    let ln = top_contents.iter().enumerate().find_map(|(i, &t)| if t.starts_with("#include") {
        Some(i)
    } else {
        None
    }).unwrap();
    let itp_line = format!("#include \"{}\"", itp_path.trim_start_matches(r"\\?\"));
    top_contents.insert(ln + 1, itp_line.as_str());
    top_contents.insert(top_contents.len() - 1, "LIG                 1");
    let new_top = top_contents.join("\n");
    let mut new_top_file = File::create(temp_dir.join("topol.top")).unwrap();
    File::write_all(&mut new_top_file, new_top.as_bytes()).unwrap();

    // include ligand structure into protein
    let protein_gro = temp_dir.join(protein_out).display().to_string();
    let structure = fs::read_to_string(protein_gro).unwrap();
    let mut struct_contents: Vec<&str> = structure.split("\n").collect();
    let protein_atom_num: usize = struct_contents[1].trim().parse().unwrap();
    let ligand_gro = fs::read_to_string(temp_dir.join("LIG.gro")).unwrap();
    let ligand_gro: Vec<&str> = ligand_gro.split("\n").collect();
    let ligand_atoms_num: usize = ligand_gro[1].trim().parse().unwrap();
    let ligand_atoms_gro = ligand_gro[2..(ligand_atoms_num + 2)].to_vec().join("\n");
    struct_contents.insert(protein_atom_num + 2, &ligand_atoms_gro);
    let total_atoms_num = format!("{:5}", protein_atom_num + ligand_atoms_num);
    struct_contents[1] = total_atoms_num.as_str();
    let new_gro = struct_contents.join("\n");
    let complex_gro_name = format!("MMPBSA_{}_{}.gro", rec_name, lig_name);
    let complex_gro_path = temp_dir.join(complex_gro_name);
    let mut new_gro_file = File::create(&complex_gro_path).unwrap();
    File::write_all(&mut new_gro_file, new_gro.as_bytes()).unwrap();

    // make complex ndx
    make_ndx(&vec!["q"], temp_dir, settings, complex_gro_path.to_str().unwrap(), "", "MMPBSA_index.ndx");

    // grompp
    let mdp = current_exe().unwrap().parent().unwrap().join("include").join("md.mdp");
    grompp(&vec![], temp_dir, settings, mdp.to_str().unwrap(), complex_gro_path.to_str().unwrap(), "../md.tpr");
}

fn calc_charge(lig_name: &str, temp_dir: &Path, method: &String, basis: &String, total_charge: i32, multiplicity: usize, settings: &Settings) {
    let lig_file = format!("MMPBSA_docking_{}.pdb", lig_name) ;
    let lig_pdb = PDB::from(temp_dir.join(&lig_file).to_str().unwrap());
    let elements = lig_pdb.models[0].get_elements();
    let coord = lig_pdb.models[0].get_coordinates();
    
    if settings.chg_m.as_ref().unwrap().eq("acpype") {
        let amber_home = settings.amber_dir.as_ref().unwrap();
        // Add ENV Var
        env::set_var("AMBERHOME", amber_home);
        // Add PATH
        let path = env::var("PATH").unwrap();
        env::set_var("PATH", format!("{}:{}", path, Path::new(amber_home).join("bin").to_str().unwrap()));
        Command::new(Path::new(amber_home).join("bin").join("antechamber"))
            .args(vec!["-i", "LIG.mol2", 
                       "-fi", "mol2", 
                       "-o", "LIG_c.mol2", 
                       "-fo", "mol2", 
                       "-nc", total_charge.to_string().as_str(), 
                       "-m", multiplicity.to_string().as_str(), 
                       "-s", "2", 
                       "-df", "2", 
                       "-at", "amber", 
                       "-c", "bcc", 
                       "-ek", "maxcyc=0", 
                       "-pf", "y", 
                       "-gn", settings.nkernels.to_string().as_str(),
                       "-dr", if settings.debug_mode {"y"} else {"n"}
            ])
            .current_dir(temp_dir)
            .stdin(Stdio::inherit())
            .stdout(if settings.debug_mode { Stdio::inherit() } else { Stdio::null() })
            .stderr(Stdio::inherit())
            .status()
            .expect("Cannot properly run antechamber");
        // multiwfn(&vec!["100", "2", "3", "LIG.chg", "0", "q"], settings, "LIG_c.mol2", temp_dir);
        // fuck Multiwfn outputs chg with mass
        let new_lig = MOL2::from(temp_dir.join("LIG_c.mol2").to_str().unwrap());
        new_lig.to_chg(temp_dir.join("LIG.chg").to_str().unwrap());
    } else if settings.chg_m.as_ref().unwrap().eq("gaussian") {
        // write gjf file
        let level = format!("{}/{} em=GD3BJ", method, basis);
        let mut gjf = File::create(temp_dir.join("LIG.gjf")).unwrap();
        writeln!(&mut gjf, "%nproc={}", settings.nkernels * 2).unwrap();
        writeln!(&mut gjf, "%chk=LIG.chk").unwrap();
        writeln!(&mut gjf, "# {}", level).unwrap();
        writeln!(&mut gjf, "").unwrap();
        writeln!(&mut gjf, "{}", &lig_file).unwrap();
        writeln!(&mut gjf, "").unwrap();
        writeln!(&mut gjf, "{} {}", total_charge, multiplicity).unwrap();
        for (i, a) in coord.rows().into_iter().enumerate() {
            writeln!(&mut gjf, " {:2}{:27.8}{:14.8}{:14.8}", elements[i], a[0], a[1], a[2]).unwrap();
        }
        writeln!(&mut gjf, "").unwrap();

        let infile = File::open(temp_dir.join("LIG.gjf")).unwrap();
        let outfile = File::create(temp_dir.join("LIG.out")).unwrap();
        let gauss_dir = settings.gaussian_dir.as_ref().unwrap();
        // Add ENV Var
        env::set_var("GAUSS_EXEDIR", gauss_dir);
        // Add PATH
        let path = env::var("PATH").unwrap();
        env::set_var("PATH", format!("{}:{}", path, gauss_dir));

        let gaussian_status = Command::new(Path::new(gauss_dir).join(settings.gaussian_exe.as_ref().unwrap()).to_str().unwrap())
            .current_dir(temp_dir)
            .stdin(Stdio::from(infile))
            .stdout(Stdio::from(outfile))
            .stderr(Stdio::inherit())
            .status()
            .expect("Failed to start process");
        if gaussian_status.code() != Some(0) {
            println!("Gaussian not normally exited. Change calculation level.");
            exit(1);
        }
        Command::new(Path::new(gauss_dir).join("formchk").to_str().unwrap())
            .current_dir(temp_dir)
            .arg("LIG.chk")
            .stdout(Stdio::null())
            .stderr(Stdio::inherit())
            .status()
            .expect("Failed to start process");
        let fchk_path = if cfg!(windows) {
            temp_dir.join("LIG.fch")
        } else {
            temp_dir.join("LIG.fchk")
        };
        multiwfn(&vec!["7", "18", "1", "y", "0", "0", "q"], settings, 
                fchk_path.to_str().unwrap().trim_start_matches(r"\\?\"), 
                Path::new(temp_dir.to_str().unwrap().trim_start_matches(r"\\?\"))).unwrap();
    }
}

fn combine_flex(protein_name: &String, flex_name: &String, temp_dir: &Path) -> PDB {
    let protein = PDB::from(temp_dir.join(protein_name).to_str().unwrap());
    let flex = PDB::from(temp_dir.join(&flex_name).to_str().unwrap());
    let mut models: Vec<PDBModel> = vec![];
    for (i, flex_model) in flex.models.iter().enumerate() {
        let mut pro = protein.models[0].clone();
        pro.modelid = i as i32 + 1;
        for f_atom in &flex_model.atoms {
            let pos = find_res_by_name_chain(&pro, f_atom.resid, &f_atom.chainname) + 1;
            pro.insert_atoms(pos, f_atom);
        }
        models.push(pro);
    }
    PDB::new(&models)
}

fn find_res_by_name_chain(pro_mdl: &PDBModel, ref_resid: i32, chain_id: &String) -> usize {
    let cur_res = pro_mdl.atoms.iter().enumerate().filter_map(|(i, a)| 
        if a.resid == ref_resid && a.chainname.eq(chain_id) {
            Some((i, a))
        } else {
            None
        }).last().unwrap();
    cur_res.0
}