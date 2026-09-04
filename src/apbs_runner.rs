#![allow(dead_code)]
#![allow(unused_assignments)]

// APBS routines - Front-end routines
// Port of src/routines.c

use apbs_generic::error::{ApbsError, ApbsResult};
use apbs_generic::nosh::{NOsh, NOshCalc, NOshPrint};
use apbs_generic::valist::Valist;
use apbs_generic::vpbe::Vpbe;

use std::collections::HashMap;
use std::cell::RefCell;
// use std::fs;
// use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::OnceLock;

thread_local! {
    static STDOUT_CAPTURE: RefCell<Option<Vec<String>>> = const { RefCell::new(None) };
}

macro_rules! println {
    () => {
        std::println!()
    };
    ($($arg:tt)*) => {{
        let line = format!($($arg)*);
        let captured = STDOUT_CAPTURE.with(|capture| {
            let mut capture = capture.borrow_mut();
            if let Some(lines) = capture.as_mut() {
                lines.push(line.clone());
                true
            } else {
                false
            }
        });
        if !captured {
            std::println!("{}", line);
        }
    }};
}

// /// Atomic force components.
// /// Port of AtomForce from routines.h line 96.
// #[derive(Debug, Clone, Default)]
// pub struct AtomForce {
//     pub ib_force: [f64; 3],   // Ion-boundary force
//     pub qf_force: [f64; 3],   // Charge-field force
//     pub db_force: [f64; 3],   // Dielectric boundary force
//     pub sasa_force: [f64; 3], // SASA force (coupled to gamma)
//     pub sav_force: [f64; 3],  // SAV force (coupled to press)
//     pub wca_force: [f64; 3],  // WCA integral force (coupled to bconc)
// }

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

fn kt_to_kjmol(temp: f64) -> f64 {
    apbs_generic::vunit::KB * 1e-3 * apbs_generic::vunit::NA * temp
}

fn report_mg_total_energy(
    pbeparm: &apbs_generic::pbeparm::PBEparm,
    energy_kt: f64,
) {
    let conv = kt_to_kjmol(pbeparm.temp);
    println!("  Total electrostatic energy = {:1.12E} kJ/mol", energy_kt * conv);
}

fn report_mg_energy_components(
    pbeparm: &apbs_generic::pbeparm::PBEparm,
    vpmg: &apbs_mg::vpmg::Vpmg,
) {
    let conv = kt_to_kjmol(pbeparm.temp);
    let total = vpmg.energy(1);
    let qf = vpmg.qf_energy(1);
    let qm = vpmg.qm_energy(1);
    let diel = vpmg.diel_energy(1);

    println!("  Total electrostatic energy = {:1.12E} kJ/mol", total * conv);
    println!("  Fixed charge energy = {} kJ/mol", 0.5 * qf * conv);
    println!("  Mobile charge energy = {} kJ/mol", qm * conv);
    println!("  Dielectric energy = {} kJ/mol", diel * conv);
    println!("  Per-atom energies:");
    let num_atoms = vpmg.pbe.alist.number_atoms();
    for i in 0..num_atoms {
        println!(
            "      Atom {}:  {:1.12E} kJ/mol",
            i,
            0.5 * vpmg.qf_atom_energy(i) * conv
        );
    }
}

// fn ensure_bem_tool_path() {
//     use std::path::PathBuf;

//     let mut candidates = Vec::new();
//     if let Some(dir) = option_env!("APBS_BUNDLED_C_BACKEND_BINDIR") {
//         candidates.push(PathBuf::from(dir));
//     }
//     candidates.push(PathBuf::from("/tmp/apbs-cmake-build/_deps/tabi-build/bin"));
//     candidates.push(PathBuf::from("/tmp/apbs-cmake-build/temp"));

//     let Some(tool_dir) = candidates
//         .into_iter()
//         .find(|dir| dir.join("NanoShaper").is_file()) else {
//         return;
//     };

//     let current = std::env::var_os("PATH").unwrap_or_default();
//     let mut parts = std::env::split_paths(&current).collect::<Vec<_>>();
//     if parts.iter().any(|p| p == &tool_dir) {
//         return;
//     }
//     parts.insert(0, tool_dir);
//     if let Ok(joined) = std::env::join_paths(parts) {
//         std::env::set_var("PATH", joined);
//     }
// }

// fn write_valist_as_pqr(alist: &Valist, path: &Path) -> ApbsResult<()> {
//     let mut out = String::new();
//     for (idx, atom) in alist.atoms().iter().enumerate() {
//         let atom_name = if atom.atom_name.is_empty() { "X" } else { &atom.atom_name };
//         let res_name = if atom.res_name.is_empty() { "MOL" } else { &atom.res_name };
//         out.push_str(&format!(
//             "ATOM {:>6} {:<4} {:<4} {:>4} {:>12.6} {:>12.6} {:>12.6} {:>12.6} {:>12.6}\n",
//             idx + 1,
//             atom_name,
//             res_name,
//             1,
//             atom.position[0],
//             atom.position[1],
//             atom.position[2],
//             atom.charge,
//             atom.radius,
//         ));
//     }
//     fs::write(path, out).map_err(|e| ApbsError::Io(format!("Failed to write temporary PQR '{}': {}", path.display(), e)))
// }

// fn make_force_temp_pqr(base_name: &str, atom_index: usize, axis: usize, sign: &str) -> PathBuf {
//     let pid = std::process::id();
//     std::env::temp_dir().join(format!(
//         "apbs-rust-fem-force-{}-{}-a{}-d{}-{}.pqr",
//         pid, base_name, atom_index, axis, sign
//     ))
// }

enum ElecResult {
    Mg {
        energy: f64,
        vpmg: apbs_mg::vpmg::Vpmg,
    },
}

enum CalcEnergyResult {
    ElecMg {
        name: String,
        energy: f64,
        vpmg: apbs_mg::vpmg::Vpmg,
    },
    Apolar {
        name: String,
        energy: f64,
    },
}

struct CalcBlockResult {
    index: usize,
    log: String,
    result: Option<CalcEnergyResult>,
}

fn parallel_blocks_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("APBS_RUST_PARALLEL_BLOCKS")
            .map(|value| matches!(value.as_str(), "1" | "true" | "TRUE" | "yes" | "YES"))
            .unwrap_or(false)
    })
}

fn max_parallel_blocks() -> usize {
    static VALUE: OnceLock<usize> = OnceLock::new();
    *VALUE.get_or_init(|| {
        std::env::var("APBS_RUST_BLOCK_THREADS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .filter(|&value| value > 0)
            .unwrap_or(2)
    })
}

fn calc_can_run_parallel(calc: &NOshCalc) -> bool {
    match calc {
        NOshCalc::Elec(elec) => {
            elec.pbeparm.bcfl != apbs_generic::vhal::Vbcfl::Focus
                && elec.pbeparm.writes.is_empty()
                && !elec.pbeparm.writemat
        }
        NOshCalc::Apolar(_) => true,
        NOshCalc::Print(_) => false,
    }
}

fn capture_stdout<F>(f: F) -> (Option<CalcEnergyResult>, String)
where
    F: FnOnce() -> Option<CalcEnergyResult>,
{
    STDOUT_CAPTURE.with(|capture| {
        *capture.borrow_mut() = Some(Vec::new());
    });
    let result = f();
    let lines = STDOUT_CAPTURE.with(|capture| {
        capture
            .borrow_mut()
            .take()
            .unwrap_or_default()
            .join("\n")
    });
    let log = if lines.is_empty() {
        String::new()
    } else {
        format!("{}\n", lines)
    };
    (result, log)
}

/// Main APBS execution routine
pub fn run_apbs(input_file: &str) -> ApbsResult<()> {
    // Parse input file
    let mut nosh = NOsh::new();
    nosh.read(input_file)?;

    if parallel_blocks_enabled() {
        return run_apbs_parallel_segments(&nosh);
    }

    run_apbs_sequential(&nosh)
}

fn run_apbs_sequential(nosh: &NOsh) -> ApbsResult<()> {
    // Energy storage: name -> energy value (in kT for ELEC, kJ/mol for APOLAR)
    let mut energies: HashMap<String, f64> = HashMap::new();
    // Previous Vpmg for focus BC handoff
    let mut prev_vpmg: Option<apbs_mg::vpmg::Vpmg> = None;

    // Process each calculation
    for calc in &nosh.calcs {
        match calc {
            NOshCalc::Elec(elecalc) => {
                println!("Running ELEC calc: {}", elecalc.name);
                if debug_enabled() {
                    eprintln!("[DEBUG]   molecules: {:?}", elecalc.molecules);
                    eprintln!("[DEBUG]   bcfl: {:?}, pbetype: {:?}, dime: {:?}", elecalc.pbeparm.bcfl, elecalc.pbeparm.pbetype, elecalc.mgparm.dime);
                    eprintln!("[DEBUG]   mgtype: {:?}, glen: {:?}, cglen: {:?}, fglen: {:?}", elecalc.mgparm.r#type, elecalc.mgparm.glen, elecalc.mgparm.cglen, elecalc.mgparm.fglen);
                    eprintln!("[DEBUG]   nion: {}, ionq: {:?}, ionc: {:?}, ionr: {:?}", elecalc.pbeparm.nion, &elecalc.pbeparm.ionq[..3], &elecalc.pbeparm.ionc[..3], &elecalc.pbeparm.ionr[..3]);
                }
                let is_focus = elecalc.pbeparm.bcfl == apbs_generic::vhal::Vbcfl::Focus;
                let focus_ref = if is_focus { prev_vpmg.as_mut() } else { None };
                match run_elec(&nosh, elecalc, focus_ref) {
                    Ok(ElecResult::Mg { energy, vpmg }) => {
                        println!("Finished ELEC calc: {} ({:.6e} kT)", elecalc.name, energy);
                        energies.insert(elecalc.name.clone(), energy);
                        prev_vpmg = Some(vpmg);
                    }
                    Err(e) => {
                        eprintln!("ELEC calc failed: {} ({})", elecalc.name, e);
                        prev_vpmg = None;
                    }
                }
            }
            NOshCalc::Apolar(apolarcalc) => {
                println!("Running APOLAR calc: {}", apolarcalc.name);
                match run_apolar(&nosh, apolarcalc) {
                    Ok(energy) => {
                        println!("Finished APOLAR calc: {} ({:.6e} kJ/mol)", apolarcalc.name, energy);
                        energies.insert(apolarcalc.name.clone(), energy);
                    }
                    Err(e) => {
                        eprintln!("APOLAR calc failed: {} ({})", apolarcalc.name, e);
                    }
                }
            }
            NOshCalc::Print(_printcalc) => {
                // Print statements handled after all calcs
            }
        }
    }

    process_prints(nosh, &energies);
    Ok(())
}

fn run_apbs_parallel_segments(nosh: &NOsh) -> ApbsResult<()> {
    let mut energies: HashMap<String, f64> = HashMap::new();
    let mut prev_vpmg: Option<apbs_mg::vpmg::Vpmg> = None;
    let mut index = 0;

    while index < nosh.calcs.len() {
        if calc_can_run_parallel(&nosh.calcs[index]) {
            let start = index;
            while index < nosh.calcs.len() && calc_can_run_parallel(&nosh.calcs[index]) {
                index += 1;
            }
            run_parallel_segment(nosh, start, index, &mut energies, &mut prev_vpmg);
            continue;
        }

        run_one_calc_sequential(nosh, &nosh.calcs[index], &mut energies, &mut prev_vpmg);
        index += 1;
    }

    process_prints(nosh, &energies);
    Ok(())
}

fn run_parallel_segment(
    nosh: &NOsh,
    start: usize,
    end: usize,
    energies: &mut HashMap<String, f64>,
    prev_vpmg: &mut Option<apbs_mg::vpmg::Vpmg>,
) {
    let max_jobs = max_parallel_blocks();
    let mut all_results = Vec::with_capacity(end - start);

    for chunk_start in (start..end).step_by(max_jobs) {
        let chunk_end = (chunk_start + max_jobs).min(end);
        let mut chunk_results = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(chunk_end - chunk_start);
            for idx in chunk_start..chunk_end {
                handles.push(scope.spawn(move || run_one_calc_captured(nosh, idx)));
            }

            let mut results = Vec::with_capacity(handles.len());
            for handle in handles {
                match handle.join() {
                    Ok(result) => results.push(result),
                    Err(_) => {
                        eprintln!("APBS parallel block worker panicked; skipping failed block");
                    }
                }
            }
            results
        });
        all_results.append(&mut chunk_results);
    }

    all_results.sort_by_key(|result| result.index);
    for result in all_results {
        print!("{}", result.log);
        apply_calc_result(result.result, energies, prev_vpmg);
    }
}

fn run_one_calc_captured(nosh: &NOsh, index: usize) -> CalcBlockResult {
    let (result, log) = capture_stdout(|| execute_calc_without_focus(nosh, &nosh.calcs[index]));
    CalcBlockResult { index, log, result }
}

fn run_one_calc_sequential(
    nosh: &NOsh,
    calc: &NOshCalc,
    energies: &mut HashMap<String, f64>,
    prev_vpmg: &mut Option<apbs_mg::vpmg::Vpmg>,
) {
    match calc {
        NOshCalc::Elec(elecalc) => {
            println!("Running ELEC calc: {}", elecalc.name);
            if debug_enabled() {
                eprintln!("[DEBUG]   molecules: {:?}", elecalc.molecules);
                eprintln!("[DEBUG]   bcfl: {:?}, pbetype: {:?}, dime: {:?}", elecalc.pbeparm.bcfl, elecalc.pbeparm.pbetype, elecalc.mgparm.dime);
                eprintln!("[DEBUG]   mgtype: {:?}, glen: {:?}, cglen: {:?}, fglen: {:?}", elecalc.mgparm.r#type, elecalc.mgparm.glen, elecalc.mgparm.cglen, elecalc.mgparm.fglen);
                eprintln!("[DEBUG]   nion: {}, ionq: {:?}, ionc: {:?}, ionr: {:?}", elecalc.pbeparm.nion, &elecalc.pbeparm.ionq[..3], &elecalc.pbeparm.ionc[..3], &elecalc.pbeparm.ionr[..3]);
            }
            let is_focus = elecalc.pbeparm.bcfl == apbs_generic::vhal::Vbcfl::Focus;
            let focus_ref = if is_focus { prev_vpmg.as_mut() } else { None };
            let result = match run_elec(nosh, elecalc, focus_ref) {
                Ok(ElecResult::Mg { energy, vpmg }) => {
                    println!("Finished ELEC calc: {} ({:.6e} kT)", elecalc.name, energy);
                    Some(CalcEnergyResult::ElecMg {
                        name: elecalc.name.clone(),
                        energy,
                        vpmg,
                    })
                }
                Err(e) => {
                    eprintln!("ELEC calc failed: {} ({})", elecalc.name, e);
                    None
                }
            };
            apply_calc_result(result, energies, prev_vpmg);
        }
        NOshCalc::Apolar(apolarcalc) => {
            println!("Running APOLAR calc: {}", apolarcalc.name);
            let result = match run_apolar(nosh, apolarcalc) {
                Ok(energy) => {
                    println!("Finished APOLAR calc: {} ({:.6e} kJ/mol)", apolarcalc.name, energy);
                    Some(CalcEnergyResult::Apolar {
                        name: apolarcalc.name.clone(),
                        energy,
                    })
                }
                Err(e) => {
                    eprintln!("APOLAR calc failed: {} ({})", apolarcalc.name, e);
                    None
                }
            };
            apply_calc_result(result, energies, prev_vpmg);
        }
        NOshCalc::Print(_) => {}
    }
}

fn execute_calc_without_focus(nosh: &NOsh, calc: &NOshCalc) -> Option<CalcEnergyResult> {
    match calc {
        NOshCalc::Elec(elecalc) => {
            println!("Running ELEC calc: {}", elecalc.name);
            match run_elec(nosh, elecalc, None) {
                Ok(ElecResult::Mg { energy, vpmg }) => {
                    println!("Finished ELEC calc: {} ({:.6e} kT)", elecalc.name, energy);
                    Some(CalcEnergyResult::ElecMg {
                        name: elecalc.name.clone(),
                        energy,
                        vpmg,
                    })
                }
                Err(e) => {
                    eprintln!("ELEC calc failed: {} ({})", elecalc.name, e);
                    None
                }
            }
        }
        NOshCalc::Apolar(apolarcalc) => {
            println!("Running APOLAR calc: {}", apolarcalc.name);
            match run_apolar(nosh, apolarcalc) {
                Ok(energy) => {
                    println!("Finished APOLAR calc: {} ({:.6e} kJ/mol)", apolarcalc.name, energy);
                    Some(CalcEnergyResult::Apolar {
                        name: apolarcalc.name.clone(),
                        energy,
                    })
                }
                Err(e) => {
                    eprintln!("APOLAR calc failed: {} ({})", apolarcalc.name, e);
                    None
                }
            }
        }
        NOshCalc::Print(_) => None,
    }
}

fn apply_calc_result(
    result: Option<CalcEnergyResult>,
    energies: &mut HashMap<String, f64>,
    prev_vpmg: &mut Option<apbs_mg::vpmg::Vpmg>,
) {
    match result {
        Some(CalcEnergyResult::ElecMg { name, energy, vpmg }) => {
            energies.insert(name, energy);
            *prev_vpmg = Some(vpmg);
        }
        Some(CalcEnergyResult::Apolar { name, energy }) => {
            energies.insert(name, energy);
        }
        None => {
            *prev_vpmg = None;
        }
    }
}

fn process_prints(nosh: &NOsh, energies: &HashMap<String, f64>) {
    // Process PRINT statements
    if debug_enabled() {
        eprintln!("[DEBUG] Energies collected: {:?}", energies.keys().collect::<Vec<_>>());
    }
    // Conversion factor: kT -> kJ/mol = Vunit_kb * 1e-3 * Vunit_Na * temp
    // Vunit_kb = 1.3806581e-23 J/K, Vunit_Na = 6.0221367e+23 /mol
    // Using Rust constants: KB * 1e-3 * NA * temp
    for calc in &nosh.calcs {
        if let NOshCalc::Print(print_stmt) = calc {
            if debug_enabled() {
                eprintln!("[DEBUG] Processing PRINT: {:?}", print_stmt);
            }
            match print_stmt {
                NOshPrint::ElecEnergy { left, right } => {
                    if right.is_empty() {
                        // Single energy
                        if let Some(&energy) = energies.get(left) {
                            let temp = nosh.get_elec_temp(left).unwrap_or(298.15);
                            let energy_kjmol = kt_to_kjmol(temp) * energy;
                            println!("  Global net ELEC energy = {:.12E} kJ/mol", energy_kjmol);
                        } else {
                            eprintln!("Warning: calc '{}' not found for print", left);
                        }
                    } else {
                        // Subtraction: left - right
                        match (energies.get(left).copied(), energies.get(right).copied()) {
                            (Some(e_left), Some(e_right)) => {
                                let temp = nosh.get_elec_temp(left).unwrap_or(298.15);
                                let delta_kjmol = kt_to_kjmol(temp) * (e_left - e_right);
                                println!("  Global net ELEC energy = {:.12E} kJ/mol", delta_kjmol);
                            }
                            _ => {
                                eprintln!(
                                    "Warning: missing energy for print '{} - {}'; skipping",
                                    left, right
                                );
                            }
                        }
                    }
                }
                NOshPrint::ApolEnergy { name } => {
                    if let Some(&energy) = energies.get(name) {
                        // APOLAR energy is already in kJ/mol (gamma is in kJ/mol/A^2)
                        println!("  Global net APOL energy = {:.12E} kJ/mol", energy);
                    } else {
                        eprintln!("Warning: calc '{}' not found for print", name);
                    }
                }
                NOshPrint::Raw { print_type, .. } => {
                    eprintln!(
                        "Warning: unsupported PRINT statement was preserved but not executed: {}",
                        print_type
                    );
                }
            }
        }
    }
}

/// Run an ELEC (electrostatics) calculation
fn run_elec(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
    pmg_old: Option<&mut apbs_mg::vpmg::Vpmg>,
) -> ApbsResult<ElecResult> {
    // Load molecules
    let mut molecules = Vec::new();
    for mol_name in &elec.molecules {
        let mol_path = nosh.get_mol_path(mol_name)?;
        let mut alist = Valist::new();
        match elec.input_format {
            apbs_generic::nosh::NOshInputFormat::Pqr => {
                alist.read_pqr(None, &mol_path)?;
            }
            apbs_generic::nosh::NOshInputFormat::Pdb => {
                alist.read_pdb(None, &mol_path)?;
            }
            _ => {
                return Err(ApbsError::UnsupportedFormat(
                    "Only PQR and PDB input supported".to_string(),
                ));
            }
        }
        molecules.push(Arc::new(alist));
    }

    // Create PBE object
    let pbeparm = &elec.pbeparm;
    let mgparm = &elec.mgparm;

    // Create PBE
    let alist = if !molecules.is_empty() {
        molecules[0].clone()
    } else {
        return Err(ApbsError::InvalidParameter("No molecules loaded".to_string()));
    };

    let vpbe = Arc::new(Vpbe::new(
        alist.clone(),
        pbeparm.nion,
        &pbeparm.ionc,
        &pbeparm.ionr,
        &pbeparm.ionq,
        pbeparm.temp,
        pbeparm.pdie,
        pbeparm.sdie,
        pbeparm.srad,
        pbeparm.bcfl as i32,
        pbeparm.sdens,
        pbeparm.zmem,
        pbeparm.lmem,
        pbeparm.mdie,
        pbeparm.memv,
    )?);

    // Create MG solver
    let mut mgparm_clone = mgparm.clone();
    for d in 0..3 {
        if matches!(mgparm_clone.cmeth, apbs_generic::mgparm::MGparmCentMeth::Molecule) && mgparm_clone.center[d] == 0.0 {
            mgparm_clone.center[d] = vpbe.solute_center[d];
        }
        if matches!(mgparm_clone.ccmeth, apbs_generic::mgparm::MGparmCentMeth::Molecule) && mgparm_clone.ccenter[d] == 0.0 {
            mgparm_clone.ccenter[d] = vpbe.solute_center[d];
        }
        if matches!(mgparm_clone.fcmeth, apbs_generic::mgparm::MGparmCentMeth::Molecule) && mgparm_clone.fcenter[d] == 0.0 {
            mgparm_clone.fcenter[d] = vpbe.solute_center[d];
        }
    }
    // Compute grid spacing if not set (for mg-auto)
    for d in 0..3 {
        if mgparm_clone.grid[d] == 0.0 && mgparm_clone.glen[d] > 0.0 && mgparm_clone.dime[d] > 1 {
            mgparm_clone.grid[d] = mgparm_clone.glen[d] / (mgparm_clone.dime[d] - 1) as f64;
        }
    }
    let mut pmgp = apbs_mg::vpmgp::Vpmgp::new(&mgparm_clone);
    pmgp.bcfl = pbeparm.bcfl;
    // Set nonlinearity flag based on PBE type
    pmgp.nonlin = match pbeparm.pbetype {
        apbs_generic::vhal::VhalPBEType::NPBE => 1,
        apbs_generic::vhal::VhalPBEType::NRPBE => 1,
        _ => 0,
    };
    if debug_enabled() {
        eprintln!("[DEBUG]   Vpmgp: nx={}, ny={}, nz={}, hx={:.4}, hy={:.4}, hz={:.4}, xmin={:.2}, ymin={:.2}, zmin={:.2}, bcfl={:?}", pmgp.nx, pmgp.ny, pmgp.nz, pmgp.hx, pmgp.hy, pmgp.hzed, pmgp.xmin, pmgp.ymin, pmgp.zmin, pmgp.bcfl);
    }

    let focus_flag = if pmg_old.is_some() { 1 } else { 0 };
    let mut vpmg = apbs_mg::vpmg::Vpmg::new(
        pmgp,
        vpbe,
        focus_flag,
        pmg_old.as_deref(),
        pbeparm,
        pbeparm.calc_energy,
    );

    // Set partition mask to 1.0 everywhere (single partition)
    vpmg.unset_part();

    // Fill coefficients
    vpmg.fillco(
        pbeparm.srfm,
        pbeparm.swin,
        mgparm.chgm,
    )?;

    // Debug: check coefficient arrays
    let nf = (vpmg.pmgp.nx * vpmg.pmgp.ny * vpmg.pmgp.nz) as usize;
    let a1_nan = vpmg.a1cf[..nf].iter().any(|x| x.is_nan());
    let f_nan = vpmg.fcf[..nf].iter().any(|x| x.is_nan());
    let c_nan = vpmg.ccf[..nf].iter().any(|x| x.is_nan());
    let fcf_nonzero = vpmg.fcf[..nf].iter().filter(|x| x.abs() > 1.0e-30).count();
    let fcf_max = vpmg.fcf[..nf].iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let num_atoms = vpmg.pbe.alist.number_atoms();
    let charge_total: f64 = (0..num_atoms).map(|i| vpmg.pbe.alist.get_atom(i).charge).sum();
    if debug_enabled() {
        eprintln!("[DEBUG]   fillco done: atoms={}, charge_total={:.4e}, a1cf_nan={}, fcf_nan={}, ccf_nan={}", num_atoms, charge_total, a1_nan, f_nan, c_nan);
        eprintln!("[DEBUG]   fcf stats: max={:.4e}, nonzero={}/{}, a1cf[0]={:.4e}, fcf[0]={:.4e}, ccf[0]={:.4e}", fcf_max, fcf_nonzero, nf, vpmg.a1cf[0], vpmg.fcf[0], vpmg.ccf[0]);
        let num_atoms_dbg = vpmg.pbe.alist.number_atoms();
        if num_atoms_dbg > 0 {
            let a0 = vpmg.pbe.alist.get_atom(0);
            let a1 = vpmg.pbe.alist.get_atom(1);
            eprintln!("[DEBUG]   atom[0]: pos=({:.2},{:.2},{:.2}), chg={:.4e}, rad={:.2}", a0.position[0], a0.position[1], a0.position[2], a0.charge, a0.radius);
            eprintln!("[DEBUG]   atom[1]: pos=({:.2},{:.2},{:.2}), chg={:.4e}, rad={:.2}", a1.position[0], a1.position[1], a1.position[2], a1.charge, a1.radius);
            eprintln!("[DEBUG]   grid: xmin={:.2} xmax={:.2} ymin={:.2} ymax={:.2} zmin={:.2} zmax={:.2}", vpmg.pmgp.xmin, vpmg.pmgp.xmax, vpmg.pmgp.ymin, vpmg.pmgp.ymax, vpmg.pmgp.zmin, vpmg.pmgp.zmax);
            eprintln!("[DEBUG]   zmagic={:.6e}, hx={:.4e}, hy={:.4e}, hz={:.4e}", vpmg.pbe.zmagic, vpmg.pmgp.hx, vpmg.pmgp.hy, vpmg.pmgp.hzed);
            eprintln!("[DEBUG]   charge_meth={:?}, chgm passed to fillco", vpmg.charge_meth);
        }
    }

    // Solve
    vpmg.solve()?;

    // Debug: check solution
    let u_nan = vpmg.u[..nf].iter().any(|x| x.is_nan());
    if debug_enabled() {
        eprintln!("[DEBUG]   solve done: u_nan={}, u[0]={:.4e}", u_nan, vpmg.u[0]);
    }

    // Compute energy
    // For focusing: compute external energy contributions from coarse grid first
    if let Some(old) = pmg_old {
        vpmg.compute_ext_energy(old);
    }
    let energy = vpmg.energy(pbeparm.calc_energy as i32);
    match pbeparm.calc_energy {
        apbs_generic::pbeparm::PBEparmCalcEnergy::Total => {
            report_mg_total_energy(pbeparm, energy);
        }
        apbs_generic::pbeparm::PBEparmCalcEnergy::Comps => {
            report_mg_energy_components(pbeparm, &vpmg);
        }
        apbs_generic::pbeparm::PBEparmCalcEnergy::No => {}
    }

    // Write output
    for write_spec in &pbeparm.writes {
        write_output(&vpmg, write_spec)?;
    }

    Ok(ElecResult::Mg { energy, vpmg })
}

#[cfg(any())]
fn run_fem(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
) -> ApbsResult<f64> {
    if let Some(result) = run_legacy_solver_backend(nosh, elec, LegacySolverKind::Fem)? {
        return Ok(result.energy_kt);
    }
    Err(ApbsError::Fetk(format!(
        "FEM calculation '{}' is unavailable because the apbs-fem crate was removed; use a configured legacy backend",
        elec.name
    )))
}

#[cfg(any())]
fn run_fem(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
) -> ApbsResult<f64> {
    let femparm = elec.femparm.as_ref().ok_or_else(|| {
        ApbsError::Fetk("FEM calculation requested without FEM parameters".to_string())
    })?;
    femparm.check()?;

    if let Some(result) = run_legacy_solver_backend(nosh, elec, LegacySolverKind::Fem)? {
        return Ok(result.energy_kt);
    }
    if matches!(elec_backend_mode(), ElecBackendMode::Legacy) {
        return Err(ApbsError::Fetk(format!(
            "FEM calculation '{}' requested APBS_ELEC_BACKEND=legacy but no usable legacy backend was found",
            elec.name
        )));
    }

    if !apbs_fem::fe_enabled() {
        return Err(ApbsError::Fetk(format!(
            "FEM calculation '{}' parsed successfully, but FEtk support is not enabled in this build (mode={:?})",
            elec.name,
            elec_backend_mode()
        )));
    }

    let pbeparm = &elec.pbeparm;
    let molecules = load_input_molecules(nosh, elec)?;
    let mol_index = pbeparm
        .molid
        .map(|id| id.saturating_sub(1) as usize)
        .unwrap_or(0);
    let alist = molecules.get(mol_index).cloned().ok_or_else(|| {
        ApbsError::InvalidParameter(format!(
            "FEM calculation '{}' requested molid={}, but only {} molecule(s) were loaded",
            elec.name,
            pbeparm.molid.unwrap_or(1),
            molecules.len()
        ))
    })?;

    if pbeparm.use_diel_map {
        eprintln!(
            "Warning: FEM calc '{}' ignores external dielectric maps, matching APBS C initFE() behavior",
            elec.name
        );
    }
    if pbeparm.use_kappa_map {
        eprintln!(
            "Warning: FEM calc '{}' ignores external kappa maps, matching APBS C initFE() behavior",
            elec.name
        );
    }
    if pbeparm.use_charge_map {
        eprintln!(
            "Warning: FEM calc '{}' ignores external charge maps, matching APBS C initFE() behavior",
            elec.name
        );
    }
    let cube_center = alist.center;
    let mol_path = nosh.get_mol_path(&elec.molecules[mol_index])?;
    let (energy, fetk) = solve_fem_native(nosh, elec, pbeparm, femparm, &mol_path, cube_center)?;

    match pbeparm.calc_force {
        apbs_generic::pbeparm::PBEparmCalcForce::No => {}
        apbs_generic::pbeparm::PBEparmCalcForce::Total => {
            if elec.input_format != apbs_generic::nosh::NOshInputFormat::Pqr {
                eprintln!(
                    "Warning: FEM calc '{}' requested total force output, but native FEM numerical force currently requires PQR input; skipping",
                    elec.name
                );
            } else {
                let forces = compute_fem_forces_numerical(nosh, elec, pbeparm, femparm, &alist, &elec.name)?;
                let conv = apbs_generic::vunit::KB * pbeparm.temp * 1e-3 * apbs_generic::vunit::NA;
                let net = forces.iter().fold([0.0; 3], |mut acc, f| {
                    acc[0] += f[0];
                    acc[1] += f[1];
                    acc[2] += f[2];
                    acc
                });
                eprintln!("  Calculating FEM forces...");
                eprintln!("  Printing net FEM force for molecule {} (kJ/mol/A)", pbeparm.molid.unwrap_or(0));
                eprintln!("femF tot  {:4.3e}  {:4.3e}  {:4.3e}",
                    conv * net[0], conv * net[1], conv * net[2]);
            }
        }
        apbs_generic::pbeparm::PBEparmCalcForce::Comps => {
            if elec.input_format != apbs_generic::nosh::NOshInputFormat::Pqr {
                eprintln!(
                    "Warning: FEM calc '{}' requested per-atom force output, but native FEM numerical force currently requires PQR input; skipping",
                    elec.name
                );
            } else {
                let forces = compute_fem_forces_numerical(nosh, elec, pbeparm, femparm, &alist, &elec.name)?;
                let conv = apbs_generic::vunit::KB * pbeparm.temp * 1e-3 * apbs_generic::vunit::NA;
                eprintln!("  Calculating FEM forces...");
                eprintln!("  Printing per-atom FEM forces for molecule {} (kJ/mol/A)", pbeparm.molid.unwrap_or(0));
                for (iatom, f) in forces.iter().enumerate() {
                    eprintln!("femF tot {}  {:4.3e}  {:4.3e}  {:4.3e}",
                        iatom, conv * f[0], conv * f[1], conv * f[2]);
                }
            }
        }
    }

    for write_spec in &pbeparm.writes {
        match write_output_fem(&fetk, write_spec, &elec.name) {
            Ok(()) => {}
            Err(ApbsError::UnsupportedFormat(msg)) => {
                eprintln!("Warning: {}", msg);
            }
            Err(err) => return Err(err),
        }
    }

    Ok(energy)
}

#[cfg(any())]
fn solve_fem_native(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
    pbeparm: &apbs_generic::pbeparm::PBEparm,
    femparm: &apbs_generic::femparm::FEMparm,
    mol_path: &str,
    cube_center: [f64; 3],
) -> ApbsResult<(f64, apbs_fem::vfetk::Vfetk)> {
    let mut fetk = apbs_fem::vfetk::Vfetk::new();
    fetk.init(mol_path, elec.input_format, pbeparm, femparm)?;

    if femparm.use_mesh {
        let (mesh_format, mesh_relpath) = nosh.get_mesh_info(femparm.mesh_id)?;
        let mesh_path = nosh.get_mesh_path(femparm.mesh_id)?;
        if mesh_format != "mcsf" {
            return Err(ApbsError::Fetk(format!(
                "FEM calculation '{}' requests mesh format '{}', but only MCSF is a viable FE input target",
                elec.name, mesh_format
            )));
        }
        if !apbs_fem::fe_enabled() {
            return Err(ApbsError::Fetk(format!(
                "FEM calculation '{}' requests external MCSF mesh '{}' (resolved to '{}'), but FEtk support is not enabled in this build",
                elec.name, mesh_relpath, mesh_path
            )));
        }
        fetk.load_external_mesh_file(&mesh_path)?;
    } else {
        fetk.load_cube(cube_center, femparm.glen, false)?;
        fetk.refine()?;
        fetk.refine()?;
    }

    if pbeparm.calc_energy == apbs_generic::pbeparm::PBEparmCalcEnergy::Comps {
        eprintln!(
            "Warning: FEM calc '{}' requests component energy output; APBS C also does not provide verbose FEM energy decomposition",
            elec.name
        );
    }

    let mut energy = 0.0;
    let mut prev_energy: Option<f64> = None;
    let max_cycles = femparm.maxsolve.max(1);
    for cycle in 0..max_cycles {
        if cycle > 0 {
            fetk.refine()?;
        }
        fetk.solve()?;
        if debug_enabled() {
            eprintln!(
                "[DEBUG]   FEM hierarchy: level={}, num_levels={}",
                fetk.level(),
                fetk.num_levels()
            );
        }
        energy = fetk.energy()?;
        if let Some(prev) = prev_energy {
            let de = (energy - prev).abs();
            if femparm.target_res > 0.0 && de <= femparm.target_res {
                if debug_enabled() {
                    eprintln!(
                        "[DEBUG]   FEM converged at cycle {}: |dE|={:.6e} <= targetRes={:.6e}",
                        cycle + 1,
                        de,
                        femparm.target_res
                    );
                }
                break;
            }
        }
        prev_energy = Some(energy);
    }

    Ok((energy, fetk))
}

#[cfg(any())]
fn compute_fem_forces_numerical(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
    pbeparm: &apbs_generic::pbeparm::PBEparm,
    femparm: &apbs_generic::femparm::FEMparm,
    alist: &Valist,
    calc_name: &str,
) -> ApbsResult<Vec<[f64; 3]>> {
    const DELTA: f64 = 1.0e-3;

    let mut forces = vec![[0.0; 3]; alist.number_atoms()];
    for iatom in 0..alist.number_atoms() {
        for axis in 0..3 {
            let mut plus = alist.clone();
            plus.get_atom_mut(iatom).position[axis] += DELTA;
            plus.get_statistics();
            let plus_path = make_force_temp_pqr(calc_name, iatom, axis, "plus");
            write_valist_as_pqr(&plus, &plus_path)?;

            let mut minus = alist.clone();
            minus.get_atom_mut(iatom).position[axis] -= DELTA;
            minus.get_statistics();
            let minus_path = make_force_temp_pqr(calc_name, iatom, axis, "minus");
            write_valist_as_pqr(&minus, &minus_path)?;

            let plus_energy = solve_fem_native(
                nosh,
                elec,
                pbeparm,
                femparm,
                plus_path.to_str().ok_or_else(|| ApbsError::Fetk("Invalid temporary plus-path".to_string()))?,
                plus.center,
            )?.0;
            let minus_energy = solve_fem_native(
                nosh,
                elec,
                pbeparm,
                femparm,
                minus_path.to_str().ok_or_else(|| ApbsError::Fetk("Invalid temporary minus-path".to_string()))?,
                minus.center,
            )?.0;

            let _ = fs::remove_file(&plus_path);
            let _ = fs::remove_file(&minus_path);

            forces[iatom][axis] = -(plus_energy - minus_energy) / (2.0 * DELTA);
        }
    }

    Ok(forces)
}

#[cfg(any())]
fn write_output_fem(
    fetk: &apbs_fem::vfetk::Vfetk,
    write_spec: &apbs_generic::pbeparm::WriteSpec,
    calc_name: &str,
) -> ApbsResult<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use apbs_generic::vhal::{VdataFormat, VdataType};
    use std::fs::File;
    use std::io::Write;

    let c_style_path = match write_spec.writefmt {
        VdataFormat::DX => format!("{}.dx", write_spec.stem),
        VdataFormat::AVS => format!("{}.ucd", write_spec.stem),
        VdataFormat::DXBin => format!("{}.dxbin", write_spec.stem),
        _ => write_spec.stem.clone(),
    };

    match (write_spec.writetype, write_spec.writefmt) {
        (VdataType::Pot, VdataFormat::Flat | VdataFormat::GZ) => {
            let sol = fetk.solution()?;
            let file = File::create(&write_spec.stem)
                .map_err(|e| ApbsError::Io(format!("Failed to create FEM output '{}': {}", write_spec.stem, e)))?;
            let mut writer: Box<dyn Write> = match write_spec.writefmt {
                VdataFormat::Flat => Box::new(file),
                VdataFormat::GZ => Box::new(GzEncoder::new(file, Compression::default())),
                _ => unreachable!(),
            };
            writeln!(writer, "# APBS FEM flat solution for calc {}", calc_name)
                .map_err(|e| ApbsError::Io(format!("Write error for '{}': {}", write_spec.stem, e)))?;
            writeln!(writer, "# length {}", sol.len())
                .map_err(|e| ApbsError::Io(format!("Write error for '{}': {}", write_spec.stem, e)))?;
            for value in sol {
                writeln!(writer, "{:.16e}", value)
                    .map_err(|e| ApbsError::Io(format!("Write error for '{}': {}", write_spec.stem, e)))?;
            }
            Ok(())
        }
        (VdataType::Pot | VdataType::Smol | VdataType::Sspl | VdataType::Vdw | VdataType::Ivdw | VdataType::Ndens | VdataType::Qdens,
         VdataFormat::DX | VdataFormat::AVS) => fetk.write_data(write_spec.writetype, write_spec.writefmt, &c_style_path),
        (VdataType::Pot | VdataType::Smol | VdataType::Sspl | VdataType::Vdw | VdataType::Ivdw | VdataType::Ndens | VdataType::Qdens,
         unsupported_fmt) => Err(ApbsError::UnsupportedFormat(format!(
            "FEM calc '{}' does not support native {:?} output for write type {:?}; supported formats are dx/avs and pot also supports flat/gz",
            calc_name, unsupported_fmt, write_spec.writetype
        ))),
        (other, _) => Err(ApbsError::UnsupportedFormat(format!(
            "FEM calc '{}' does not support native write type {:?}; native FEtk supports pot/smol/sspl/vdw/ivdw/ndens/qdens",
            calc_name, other
        ))),
    }
}

#[cfg(any())]
fn load_input_molecules(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
) -> ApbsResult<Vec<Arc<Valist>>> {
    let mut molecules = Vec::new();
    for mol_name in &elec.molecules {
        let mol_path = nosh.get_mol_path(mol_name)?;
        let mut alist = Valist::new();
        match elec.input_format {
            apbs_generic::nosh::NOshInputFormat::Pqr => {
                alist.read_pqr(None, &mol_path)?;
            }
            apbs_generic::nosh::NOshInputFormat::Pdb => {
                alist.read_pdb(None, &mol_path)?;
            }
            _ => {
                return Err(ApbsError::UnsupportedFormat(
                    "Only PQR and PDB molecule input are currently supported".to_string(),
                ));
            }
        }
        molecules.push(Arc::new(alist));
    }
    if molecules.is_empty() {
        return Err(ApbsError::InvalidParameter("No molecules loaded".to_string()));
    }
    Ok(molecules)
}

#[cfg(any())]
fn run_bem(
    nosh: &NOsh,
    elec: &apbs_generic::nosh::NOshElec,
) -> ApbsResult<f64> {
    use std::ffi::CString;

    let bemparm = elec.bemparm.as_ref().ok_or_else(|| {
        ApbsError::Bem("BEM calculation requested without BEM parameters".to_string())
    })?;
    bemparm.check()?;

    if let Some(result) = run_legacy_solver_backend(nosh, elec, LegacySolverKind::Bem)? {
        let kt_factor =
            apbs_generic::vunit::KB * 1e-3 * apbs_generic::vunit::NA * elec.pbeparm.temp;
        println!(
            "  Global net ELEC energy = {:.12E} kJ/mol",
            result.energy_kt * kt_factor
        );
        if let Some(coulombic_kjmol) = result.coulombic_kjmol {
            println!(
                "  Global net COULOMBIC energy = {:.12E} kJ/mol",
                coulombic_kjmol
            );
        }
        return Ok(result.energy_kt);
    }
    if matches!(elec_backend_mode(), ElecBackendMode::Legacy) {
        return Err(ApbsError::Bem(format!(
            "BEM calculation '{}' requested APBS_ELEC_BACKEND=legacy but no usable legacy backend was found",
            elec.name
        )));
    }

    let pbeparm = &elec.pbeparm;
    ensure_bem_tool_path();
    let mol_index = pbeparm
        .molid
        .map(|id| id.saturating_sub(1) as usize)
        .unwrap_or(0);
    let mol_name = elec.molecules.get(mol_index).ok_or_else(|| {
        ApbsError::InvalidParameter(format!(
            "BEM calculation '{}' requested molid={}, but only {} molecule(s) were loaded",
            elec.name,
            pbeparm.molid.unwrap_or(1),
            elec.molecules.len()
        ))
    })?;
    let mol_path = nosh.get_mol_path(mol_name)?;
    let c_mol_path = CString::new(mol_path.clone())
        .map_err(|_| ApbsError::Bem(format!("Invalid molecule path for native BEM: {}", mol_path)))?;
    let input_format = match elec.input_format {
        apbs_generic::nosh::NOshInputFormat::Pqr => 0,
        apbs_generic::nosh::NOshInputFormat::Pdb => 1,
        other => {
            return Err(ApbsError::UnsupportedFormat(format!(
                "Native BEM molecule input only supports PQR/PDB, got {:?}",
                other
            )))
        }
    };

    if pbeparm.calc_force != apbs_generic::pbeparm::PBEparmCalcForce::No {
        eprintln!(
            "Warning: BEM calc '{}' requests force output, but upstream APBS/TABIPB does not expose a real BEM force decomposition here; returning energy only",
            elec.name
        );
    }
    if !pbeparm.writes.is_empty() {
        eprintln!(
            "Warning: BEM calc '{}' requests {} APBS write output(s), but upstream BEM uses TABIPB/NanoShaper outputs rather than APBS write-data handlers",
            elec.name,
            pbeparm.writes.len()
        );
    }
    if pbeparm.writemat {
        eprintln!(
            "Warning: BEM calc '{}' requested writemat '{}', but upstream APBS BEM writemat handler is a stub; ignoring",
            elec.name,
            pbeparm.writematstem
        );
    }
    if bemparm.outdata != 0 {
        eprintln!(
            "  Native BEM calc '{}' requested TABIPB outdata={} and may emit files such as 'output.vtk' in the working directory",
            elec.name,
            bemparm.outdata
        );
    }
    let ionic_strength = (0..pbeparm.nion)
        .map(|i| pbeparm.ionc[i] * pbeparm.ionq[i] * pbeparm.ionq[i])
        .sum::<f64>();

    let output = unsafe {
        crate::bem_ffi::ApbsRustBem_run_from_file(
            c_mol_path.as_ptr(),
            input_format,
            bemparm.mesh,
            pbeparm.sdens,
            pbeparm.srad,
            pbeparm.temp,
            pbeparm.pdie,
            pbeparm.sdie,
            ionic_strength,
            bemparm.tree_order,
            bemparm.tree_n0,
            bemparm.mac,
            bemparm.outdata,
        )
    };
    if output.ok_ != 1 {
        return Err(ApbsError::Bem(format!(
            "Native TABIPB-backed BEM solve failed for calc '{}' and molecule '{}'",
            elec.name, mol_name
        )));
    }

    println!(
        "  Global net ELEC energy = {:.12E} kJ/mol",
        output.solvation_energy_
    );
    println!(
        "  Global net COULOMBIC energy = {:.12E} kJ/mol",
        output.coulombic_energy_
    );

    let kt_factor = apbs_generic::vunit::KB * 1e-3 * apbs_generic::vunit::NA * pbeparm.temp;
    Ok(output.solvation_energy_ / kt_factor)
}

/// Run an APOLAR (non-polar) calculation.
/// Port of initAPOL from routines.c line 4917.
struct ApolarRunResult {
    energy: f64,
    sasa: f64,
    atomsasa: Vec<f64>,
}

fn run_apolar(
    nosh: &NOsh,
    apolar: &apbs_generic::nosh::NOshApolar,
) -> ApbsResult<f64> {
    Ok(run_apolar_full(nosh, apolar)?.energy)
}

/// Run an APOLAR (non-polar) calculation and return both the total energy and
/// the per-atom solvent accessible surface areas.
fn run_apolar_full(
    nosh: &NOsh,
    apolar: &apbs_generic::nosh::NOshApolar,
) -> ApbsResult<ApolarRunResult> {
    let mut apolparm = apolar.apolparm.clone();

    // Load molecules
    let mut molecules = Vec::new();
    for mol_name in &apolar.molecules {
        let mol_path = nosh.get_mol_path(mol_name)?;
        let mut alist = Valist::new();
        alist.read_pqr(None, &mol_path)?;
        molecules.push(Arc::new(alist));
    }

    let alist = if !molecules.is_empty() {
        molecules[0].clone()
    } else {
        return Err(ApbsError::InvalidParameter("No molecules loaded for APOLAR".to_string()));
    };

    // Determine solute bounding box
    let natoms = alist.number_atoms();
    let mut xmin = f64::MAX;
    let mut xmax = f64::MIN;
    let mut ymin = f64::MAX;
    let mut ymax = f64::MIN;
    let mut zmin = f64::MAX;
    let mut zmax = f64::MIN;

    for i in 0..natoms {
        let atom = alist.get_atom(i);
        let pos = atom.position;
        let rad = atom.radius;
        if pos[0] + rad > xmax { xmax = pos[0] + rad; }
        if pos[0] - rad < xmin { xmin = pos[0] - rad; }
        if pos[1] + rad > ymax { ymax = pos[1] + rad; }
        if pos[1] - rad < ymin { ymin = pos[1] - rad; }
        if pos[2] + rad > zmax { zmax = pos[2] + rad; }
        if pos[2] - rad < zmin { zmin = pos[2] - rad; }
    }

    let solute_xlen = xmax - xmin;
    let solute_ylen = ymax - ymin;
    let solute_zlen = zmax - zmin;

    // Set up hash table for cell list
    let srad = apolparm.srad;
    let srad_pad = srad + 2.0 * apolparm.dpos;

    // Create cell list and accessibility object
    let nhash = [
        (solute_xlen / 0.5).max(3.0).min(apbs_generic::vhal::MAX_HASH_DIM as f64) as usize,
        (solute_ylen / 0.5).max(3.0).min(apbs_generic::vhal::MAX_HASH_DIM as f64) as usize,
        (solute_zlen / 0.5).max(3.0).min(apbs_generic::vhal::MAX_HASH_DIM as f64) as usize,
    ];
    let clist = Arc::new(apbs_generic::vclist::Vclist::new_auto(
        &alist, srad_pad, nhash,
    )?);
    // let acc = apbs_generic::vacc::Vacc::new(&alist, &clist, apolparm.sdens);

    let numatoms = alist.number_atoms();
    let mut atomsasa = vec![0.0; numatoms];

    // Calculate energies if requested
    if apolparm.calc_energy != apbs_generic::apolparm::APOLparmCalcEnergy::No {
        // Per-atom SASA (needs &mut self for surface building)
        if apolparm.gamma.abs() > apbs_generic::vhal::VSMALL {
            // Use a mutable block for surface-building operations
            {
                let mut acc_mut = apbs_generic::vacc::Vacc::new(&alist, &clist, apolparm.sdens);
                apolparm.sasa = acc_mut.sasa(srad);
                for i in 0..numatoms {
                    atomsasa[i] = acc_mut.atom_sasa(i, srad);
                }
            }
        } else {
            apolparm.sasa = 0.0;
        }

        energy_apol(&apolparm, apolparm.sasa, &atomsasa);
    }

    let energy = apolparm.gamma * apolparm.sasa;

    Ok(ApolarRunResult {
        energy,
        sasa: apolparm.sasa,
        atomsasa,
    })
}

/// Report/print apolar energy components.
/// Port of energyAPOL from routines.c line 5149.
fn energy_apol(
    apolparm: &apbs_generic::apolparm::APOLparm,
    sasa: f64,
    atomsasa: &[f64],
) {
    if debug_enabled() {
        eprintln!("\nSolvent Accessible Surface Area (SASA) for each atom:");
        for (i, &sa) in atomsasa.iter().enumerate() {
            eprintln!("  SASA for atom {}: {:1.12E}", i, sa);
        }
    }
    println!("Total solvent accessible surface area: {} A^2", sasa);

    match apolparm.calc_energy {
        apbs_generic::apolparm::APOLparmCalcEnergy::No => {},
        apbs_generic::apolparm::APOLparmCalcEnergy::Comps => {
            println!("Per-atom surface energies:");
            for (i, &sa) in atomsasa.iter().enumerate() {
                let surface_energy = apolparm.gamma * sa;
                println!(
                    "    Atom {}: {:1.12E} kJ/mol",
                    i,
                    surface_energy
                );
            }
            println!("Total surface energy = {:1.12E} kJ/mol", apolparm.gamma * sasa);
        },
        apbs_generic::apolparm::APOLparmCalcEnergy::Total => {
            println!("Total surface area: {} A^2", sasa);
            println!("Total surface energy = {:1.12E} kJ/mol", apolparm.gamma * sasa);
        },
    }
}

/// Compute apolar forces.
/// Port of forceAPOL from routines.c line 5207.
#[cfg(any())]
fn force_apol(
    acc: &apbs_generic::vacc::Vacc,
    apolparm: &mut apbs_generic::apolparm::APOLparm,
    alist: &Valist,
    clist: &Arc<apbs_generic::vclist::Vclist>,
) -> Vec<AtomForce> {
    let srad = apolparm.srad;
    let press = apolparm.press;
    let gamma = apolparm.gamma;
    let offset = apolparm.dpos;
    let bconc = apolparm.bconc;
    let natom = alist.number_atoms();

    match apolparm.calc_force {
        apbs_generic::apolparm::APOLparmCalcForce::Total => {
            let mut atom_force = AtomForce::default();

            for i in 0..natom {
                let mut dsasa = [0.0f64; 3];
                let mut dsav = [0.0f64; 3];
                let mut force = [0.0f64; 3];

                if gamma.abs() > apbs_generic::vhal::VSMALL {
                    acc.atomd_sasa(offset, srad, i, &mut dsasa);
                }
                if press.abs() > apbs_generic::vhal::VSMALL {
                    acc.atomd_sav(srad, i, &mut dsav);
                }
                if bconc.abs() > apbs_generic::vhal::VSMALL {
                    let _ = acc.wca_force_atom(apolparm, i, &mut force);
                }

                for j in 0..3 {
                    atom_force.sasa_force[j] += dsasa[j];
                    atom_force.sav_force[j] += dsav[j];
                    atom_force.wca_force[j] += force[j];
                }
            }

            if debug_enabled() {
                eprintln!("  Printing net forces (kJ/mol/A)");
                eprintln!("  Legend:");
                eprintln!("    sasa  -- SASA force");
                eprintln!("    sav   -- SAV force");
                eprintln!("    wca   -- WCA force\n");
                eprintln!("  sasa  {:4.3e} {:4.3e} {:4.3e}",
                    atom_force.sasa_force[0], atom_force.sasa_force[1], atom_force.sasa_force[2]);
                eprintln!("  sav   {:4.3e} {:4.3e} {:4.3e}",
                    atom_force.sav_force[0], atom_force.sav_force[1], atom_force.sav_force[2]);
                eprintln!("  wca   {:4.3e} {:4.3e} {:4.3e}",
                    atom_force.wca_force[0], atom_force.wca_force[1], atom_force.wca_force[2]);
            }

            vec![atom_force]
        },
        apbs_generic::apolparm::APOLparmCalcForce::Comps => {
            let mut forces = Vec::with_capacity(natom);

            if debug_enabled() {
                eprintln!("  Printing per atom forces (kJ/mol/A)");
                eprintln!("  Legend:");
                eprintln!("    tot  n -- Total force for atom n");
                eprintln!("    sasa n -- SASA force for atom n");
                eprintln!("    sav  n -- SAV force for atom n");
                eprintln!("    wca  n -- WCA force for atom n\n");
            }

            for i in 0..natom {
                let mut dsasa = [0.0f64; 3];
                let mut dsav = [0.0f64; 3];
                let mut force = [0.0f64; 3];

                if gamma.abs() > apbs_generic::vhal::VSMALL {
                    acc.atomd_sasa(offset, srad, i, &mut dsasa);
                }
                if press.abs() > apbs_generic::vhal::VSMALL {
                    acc.atomd_sav(srad, i, &mut dsav);
                }
                if bconc.abs() > apbs_generic::vhal::VSMALL {
                    let _ = acc.wca_force_atom(apolparm, i, &mut force);
                }

                let xf = -(gamma * dsasa[0] + press * dsav[0] + bconc * force[0]);
                let yf = -(gamma * dsasa[1] + press * dsav[1] + bconc * force[1]);
                let zf = -(gamma * dsasa[2] + press * dsav[2] + bconc * force[2]);

                let mut atom_force = AtomForce::default();
                atom_force.sasa_force = dsasa;
                atom_force.sav_force = dsav;
                atom_force.wca_force = force;

                if debug_enabled() {
                    eprintln!("  tot  {} {:4.3e} {:4.3e} {:4.3e}", i, xf, yf, zf);
                    eprintln!("  sasa {} {:4.3e} {:4.3e} {:4.3e}",
                        i, atom_force.sasa_force[0], atom_force.sasa_force[1], atom_force.sasa_force[2]);
                    eprintln!("  sav  {} {:4.3e} {:4.3e} {:4.3e}",
                        i, atom_force.sav_force[0], atom_force.sav_force[1], atom_force.sav_force[2]);
                    eprintln!("  wca  {} {:4.3e} {:4.3e} {:4.3e}",
                        i, atom_force.wca_force[0], atom_force.wca_force[1], atom_force.wca_force[2]);
                }

                forces.push(atom_force);
            }

            forces
        },
        _ => Vec::new(),
    }
}

/// Write output to file
fn write_output(
    vpmg: &apbs_mg::vpmg::Vpmg,
    write_spec: &apbs_generic::pbeparm::WriteSpec,
) -> ApbsResult<()> {
    use apbs_generic::vhal::{VdataType, VdataFormat};

    match write_spec.writetype {
        VdataType::Pot
        | VdataType::AtomPot
        | VdataType::Charge
        | VdataType::Smol
        | VdataType::Sspl
        | VdataType::Vdw
        | VdataType::Ivdw
        | VdataType::Lap
        | VdataType::Edens
        | VdataType::Ndens
        | VdataType::Qdens
        | VdataType::DielX
        | VdataType::DielY
        | VdataType::DielZ
        | VdataType::Kappa => {
            let (origin_dx, origin_dy, origin_dz) = match write_spec.writetype {
                VdataType::DielX => (0.5 * vpmg.pmgp.hx, 0.0, 0.0),
                VdataType::DielY => (0.0, 0.5 * vpmg.pmgp.hy, 0.0),
                VdataType::DielZ => (0.0, 0.0, 0.5 * vpmg.pmgp.hzed),
                _ => (0.0, 0.0, 0.0),
            };
            let mut grid = apbs_mg::vgrid::Vgrid::new(
                vpmg.pmgp.nx as usize,
                vpmg.pmgp.ny as usize,
                vpmg.pmgp.nz as usize,
                vpmg.pmgp.hx,
                vpmg.pmgp.hy,
                vpmg.pmgp.hzed,
                vpmg.pmgp.xmin + origin_dx,
                vpmg.pmgp.ymin + origin_dy,
                vpmg.pmgp.zmin + origin_dz,
            );
            vpmg.fill_array(&mut grid.data, write_spec.writetype)?;
            let quantity = match write_spec.writetype {
                VdataType::Charge => "CHARGE DISTRIBUTION (e)",
                VdataType::Pot => "POTENTIAL (kT/e)",
                VdataType::AtomPot => "ATOM POTENTIALS",
                VdataType::Smol => "SOLVENT ACCESSIBILITY -- MOLECULAR",
                VdataType::Sspl => "SOLVENT ACCESSIBILITY -- SPLINE",
                VdataType::Vdw => "SOLVENT ACCESSIBILITY -- VAN DER WAALS",
                VdataType::Ivdw => "ION ACCESSIBILITY -- SPLINE",
                VdataType::Lap => "POTENTIAL LAPLACIAN (kT/e/A^2)",
                VdataType::Edens => "ENERGY DENSITY (kT/e/A)^2",
                VdataType::Ndens => "ION NUMBER DENSITY (M)",
                VdataType::Qdens => "ION CHARGE DENSITY (e_c * M)",
                VdataType::DielX => "X-SHIFTED DIELECTRIC MAP",
                VdataType::DielY => "Y-SHIFTED DIELECTRIC MAP",
                VdataType::DielZ => "Z-SHIFTED DIELECTRIC MAP",
                VdataType::Kappa => "KAPPA MAP",
            };

            match write_spec.writefmt {
                VdataFormat::DX => {
                    grid.write_dx(&format!("{}.dx", write_spec.stem), quantity)?;
                }
                VdataFormat::GZ => {
                    grid.write_gz(&format!("{}.dx.gz", write_spec.stem), quantity)?;
                }
                VdataFormat::DXBin => {
                    grid.write_dxbin(&format!("{}.dxbin", write_spec.stem), quantity)?;
                }
                VdataFormat::UHBD => {
                    grid.write_uhbd(&format!("{}.grd", write_spec.stem), quantity)?;
                }
                VdataFormat::Flat => {
                    let count = if write_spec.writetype == VdataType::AtomPot {
                        vpmg.pbe.alist.number_atoms()
                    } else {
                        grid.data.len()
                    };
                    grid.write_flat_values(&format!("{}.txt", write_spec.stem), quantity, count)?;
                }
                VdataFormat::AVS => {
                    return Err(ApbsError::UnsupportedFormat(
                        "AVS output is not supported for uniform MG meshes, matching C APBS behavior".to_string(),
                    ));
                }
                VdataFormat::MCSF => {
                    return Err(ApbsError::UnsupportedFormat(
                        "MCSF output is not supported for uniform MG meshes, matching C APBS behavior".to_string(),
                    ));
                }
            }
        }
    }

    Ok(())
}

// /// Compute forces on all atoms for an ELEC calculation.
// /// Port of forceMG from routines.c line 1803.
// ///
// /// For PCF_TOTAL: returns a single AtomForce with summed forces over all atoms.
// /// For PCF_COMPS: returns one AtomForce per atom.
// /// For PCF_NONE: returns empty vec.
// pub fn force_mg(
//     pmg: &apbs_mg::vpmg::Vpmg,
//     pbeparm: &apbs_generic::pbeparm::PBEparm,
//     mgparm: &apbs_generic::mgparm::MGparm,
// ) -> Vec<AtomForce> {
//     use apbs_generic::pbeparm::PBEparmCalcForce;

//     match pbeparm.calc_force {
//         PBEparmCalcForce::Total => {
//             let n = pmg.pbe.alist.number_atoms();
//             let mut result = AtomForce::default();
//             for j in 0..n {
//                 let qf = pmg.qf_force_full(j, mgparm.chgm);
//                 let ib = pmg.ib_force(j, pbeparm.srfm).unwrap_or([0.0; 3]);
//                 let db = pmg.db_force(j, pbeparm.srfm).unwrap_or([0.0; 3]);
//                 for k in 0..3 {
//                     result.qf_force[k] += qf[k];
//                     result.ib_force[k] += ib[k];
//                     result.db_force[k] += db[k];
//                 }
//             }
//             let conversion = apbs_generic::vunit::KB * pbeparm.temp * 1e-3 * apbs_generic::vunit::NA;
//             eprintln!("  Calculating forces...");
//             eprintln!("  Printing net forces for molecule {} (kJ/mol/A)", pbeparm.molid.unwrap_or(0));
//             eprintln!("  Legend:");
//             eprintln!("    qf  -- fixed charge force");
//             eprintln!("    db  -- dielectric boundary force");
//             eprintln!("    ib  -- ionic boundary force");
//             eprintln!("  qf  {:4.3e}  {:4.3e}  {:4.3e}",
//                 conversion * result.qf_force[0],
//                 conversion * result.qf_force[1],
//                 conversion * result.qf_force[2]);
//             eprintln!("  ib  {:4.3e}  {:4.3e}  {:4.3e}",
//                 conversion * result.ib_force[0],
//                 conversion * result.ib_force[1],
//                 conversion * result.ib_force[2]);
//             eprintln!("  db  {:4.3e}  {:4.3e}  {:4.3e}",
//                 conversion * result.db_force[0],
//                 conversion * result.db_force[1],
//                 conversion * result.db_force[2]);
//             vec![result]
//         }
//         PBEparmCalcForce::Comps => {
//             let n = pmg.pbe.alist.number_atoms();
//             let mut result = Vec::with_capacity(n);
//             eprintln!("  Calculating forces...");
//             eprintln!("  Printing per-atom forces for molecule {} (kJ/mol/A)", pbeparm.molid.unwrap_or(0));
//             eprintln!("  Legend:");
//             eprintln!("    tot n -- total force for atom n");
//             eprintln!("    qf  n -- fixed charge force for atom n");
//             eprintln!("    db  n -- dielectric boundary force for atom n");
//             eprintln!("    ib  n -- ionic boundary force for atom n");
//             let conversion = apbs_generic::vunit::KB * pbeparm.temp * 1e-3 * apbs_generic::vunit::NA;
//             for j in 0..n {
//                 let qf = pmg.qf_force_full(j, mgparm.chgm);
//                 let ib = pmg.ib_force(j, pbeparm.srfm).unwrap_or([0.0; 3]);
//                 let db = pmg.db_force(j, pbeparm.srfm).unwrap_or([0.0; 3]);
//                 let total = [
//                     qf[0] + ib[0] + db[0],
//                     qf[1] + ib[1] + db[1],
//                     qf[2] + ib[2] + db[2],
//                 ];
//                 eprintln!("mgF  tot {}  {:4.3e}  {:4.3e}  {:4.3e}",
//                     j, conversion * total[0], conversion * total[1], conversion * total[2]);
//                 eprintln!("mgF  qf  {}  {:4.3e}  {:4.3e}  {:4.3e}",
//                     j, conversion * qf[0], conversion * qf[1], conversion * qf[2]);
//                 eprintln!("mgF  ib  {}  {:4.3e}  {:4.3e}  {:4.3e}",
//                     j, conversion * ib[0], conversion * ib[1], conversion * ib[2]);
//                 eprintln!("mgF  db  {}  {:4.3e}  {:4.3e}  {:4.3e}",
//                     j, conversion * db[0], conversion * db[1], conversion * db[2]);
//                 result.push(AtomForce {
//                     qf_force: qf,
//                     ib_force: ib,
//                     db_force: db,
//                     ..Default::default()
//                 });
//             }
//             result
//         }
//         _ => Vec::new(),
//     }
// }

// /// Store per-atom energy decomposition.
// /// Port of storeAtomEnergy from routines.c line 2056.
// pub fn store_atom_energy(pmg: &apbs_mg::vpmg::Vpmg) -> Vec<f64> {
//     let n = pmg.pbe.alist.number_atoms();
//     let mut energies = Vec::with_capacity(n);
//     for i in 0..n {
//         energies.push(pmg.qf_atom_energy(i));
//     }
//     energies
// }

// /// Write calculation results to flat-text file.
// /// Port of writedataFlat from routines.c line 2075.
// pub fn writedata_flat(
//     fname: &str,
//     nosh: &NOsh,
//     energies: &HashMap<String, f64>,
// ) -> ApbsResult<()> {
//     use std::fs::File;
//     use std::io::Write;

//     let mut file = File::create(fname)
//         .map_err(|e| ApbsError::Io(format!("Failed to open {}: {}", fname, e)))?;

//     // Write timestamp
//     let now = std::time::SystemTime::now()
//         .duration_since(std::time::UNIX_EPOCH)
//         .unwrap_or_default();
//     writeln!(file, "Timestamp: {} seconds since epoch", now.as_secs())
//         .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//     // Write each ELEC calculation
//     for calc in &nosh.calcs {
//         if let NOshCalc::Elec(elec) = calc {
//             let pbeparm = &elec.pbeparm;
//             let mgparm = &elec.mgparm;
//             let conversion = apbs_generic::vunit::KB * pbeparm.temp * 1e-3 * apbs_generic::vunit::NA;

//             writeln!(file, "elec name {}", elec.name)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::mgparm::MGparmCalcType;
//             match mgparm.r#type {
//                 MGparmCalcType::Dummy => writeln!(file, "    mg-dummy"),
//                 MGparmCalcType::Manual => writeln!(file, "    mg-manual"),
//                 MGparmCalcType::Auto => writeln!(file, "    mg-auto"),
//                 MGparmCalcType::Parallel => writeln!(file, "    mg-para"),
//                 _ => writeln!(file, "    mg-unknown"),
//             }.map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             writeln!(file, "    mol {}", pbeparm.molid.unwrap_or(0))
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    dime {} {} {}", mgparm.dime[0], mgparm.dime[1], mgparm.dime[2])
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::vhal::VhalPBEType;
//             match pbeparm.pbetype {
//                 VhalPBEType::NPBE => writeln!(file, "    npbe"),
//                 VhalPBEType::LPBE => writeln!(file, "    lpbe"),
//                 _ => writeln!(file, "    unknown-pbe"),
//             }.map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             if pbeparm.nion > 0 {
//                 for i in 0..pbeparm.nion as usize {
//                     writeln!(file, "    ion {:4.3} {:4.3} {:4.3}",
//                         pbeparm.ionr[i], pbeparm.ionq[i], pbeparm.ionc[i])
//                         .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                 }
//             }

//             writeln!(file, "    pdie {:4.3}", pbeparm.pdie)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    sdie {:4.3}", pbeparm.sdie)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::vhal::VsurfMeth;
//             match pbeparm.srfm {
//                 VsurfMeth::Mol => { writeln!(file, "    srfm mol").map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?; }
//                 VsurfMeth::MolSmooth => { writeln!(file, "    srfm smol").map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?; }
//                 VsurfMeth::Spline => { writeln!(file, "    srfm spl2").map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?; }
//                 _ => {}
//             }
//             writeln!(file, "    srad {:4.3}", pbeparm.srad)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::vhal::Vbcfl;
//             match pbeparm.bcfl {
//                 Vbcfl::Zero => writeln!(file, "    bcfl zero"),
//                 Vbcfl::SDH => writeln!(file, "    bcfl sdh"),
//                 Vbcfl::MDH => writeln!(file, "    bcfl mdh"),
//                 Vbcfl::Focus => writeln!(file, "    bcfl focus"),
//                 Vbcfl::Map => writeln!(file, "    bcfl map"),
//                 Vbcfl::Mem => writeln!(file, "    bcfl mem"),
//                 _ => writeln!(file, "    bcfl unknown"),
//             }.map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             writeln!(file, "    temp {:4.3}", pbeparm.temp)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             if let Some(&energy) = energies.get(&elec.name) {
//                 use apbs_generic::pbeparm::PBEparmCalcEnergy;
//                 match pbeparm.calc_energy {
//                     PBEparmCalcEnergy::Total => {
//                         writeln!(file, "        totEnergy {:1.12E} kJ/mol", energy * conversion)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                     PBEparmCalcEnergy::Comps => {
//                         writeln!(file, "        totEnergy {:1.12E} kJ/mol", energy * conversion)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                     _ => {}
//                 }
//             }

//             writeln!(file, "    end")
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "end")
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//         }
//     }

//     // Handle PRINT statements
//     for calc in &nosh.calcs {
//         if let NOshCalc::Print(print_stmt) = calc {
//             match print_stmt {
//                 NOshPrint::ElecEnergy { left, right } => {
//                     if right.is_empty() {
//                         if let Some(&energy) = energies.get(left) {
//                             let temp = nosh.get_elec_temp(left).unwrap_or(298.15);
//                             let conversion = apbs_generic::vunit::KB * temp * 1e-3 * apbs_generic::vunit::NA;
//                             writeln!(file, "print energy {} ({}) end", left, left)
//                                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                             writeln!(file, "    globalEnergy {:1.12E} kJ/mol", energy * conversion)
//                                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         }
//                     } else {
//                         let e_left = energies.get(left).copied().unwrap_or(0.0);
//                         let e_right = energies.get(right).copied().unwrap_or(0.0);
//                         let temp = nosh.get_elec_temp(left).unwrap_or(298.15);
//                         let conversion = apbs_generic::vunit::KB * temp * 1e-3 * apbs_generic::vunit::NA;
//                         writeln!(file, "print energy {} ({}) - {} ({}) end", left, left, right, right)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "    globalEnergy {:1.12E} kJ/mol", (e_left - e_right) * conversion)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                 }
//                 NOshPrint::ApolEnergy { name } => {
//                     if let Some(&energy) = energies.get(name) {
//                         writeln!(file, "print APOL energy {} ({}) end", name, name)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "    globalEnergy {:1.12E} kJ/mol", energy)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                 }
//                 _ => {}
//             }
//         }
//     }

//     Ok(())
// }

/// Result of one ELEC or APOLAR calculation executed in-process.
#[derive(Debug, Clone)]
pub(crate) struct SolverCalcResult {
    pub name: String,
    pub kind: SolverCalcKind,
    /// Per-atom values. For ELEC they are kJ/mol energies; for APOLAR they are
    /// per-atom solvent accessible surface areas in A^2.
    pub per_atom: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SolverCalcKind {
    Elec,
    Apolar,
}

#[derive(Debug, Default)]
pub(crate) struct SolverRun {
    pub calcs: Vec<SolverCalcResult>,
    /// Captured APBS-style log of the run (useful for debug output).
    pub log: String,
}

/// mg-auto focusing levels are named `<calc>_<level>` (e.g. `rec_SOL_0`),
/// while the finest level keeps the plain `<calc>` name from the input file.
fn is_focus_level_name(name: &str) -> bool {
    name.rsplit_once('_').map_or(false, |(_, suffix)| {
        !suffix.is_empty() && suffix.bytes().all(|b| b.is_ascii_digit())
    })
}

/// Solve an APBS input file in-process (no sub-process) and return per-atom
/// energies for the finest focusing level of every ELEC/APOLAR calculation.
///
/// All normal APBS progress output is captured (and returned as `log`) so the
/// caller can decide whether to surface it, e.g. in debug mode.
pub(crate) fn run_apbs_in_process(input_file: &str) -> ApbsResult<SolverRun> {
    let mut nosh = NOsh::new();
    nosh.read(input_file)?;

    // Route every println! in this module into the capture buffer.
    STDOUT_CAPTURE.with(|capture| {
        *capture.borrow_mut() = Some(Vec::new());
    });

    let mut run = SolverRun::default();
    // Previous Vpmg for focus BC handoff between focusing levels.
    let mut prev_vpmg: Option<apbs_mg::vpmg::Vpmg> = None;

    for calc in &nosh.calcs {
        match calc {
            NOshCalc::Elec(elecalc) => {
                println!("Running ELEC calc: {}", elecalc.name);
                let is_focus = elecalc.pbeparm.bcfl == apbs_generic::vhal::Vbcfl::Focus;
                let focus_ref = if is_focus {
                    prev_vpmg.as_mut()
                } else {
                    None
                };
                match run_elec(&nosh, elecalc, focus_ref) {
                    Ok(ElecResult::Mg { energy, vpmg }) => {
                        println!("Finished ELEC calc: {} ({:.6e} kT)", elecalc.name, energy);
                        if !is_focus_level_name(&elecalc.name) {
                            run.calcs.push(SolverCalcResult {
                                name: elecalc.name.clone(),
                                kind: SolverCalcKind::Elec,
                                per_atom: per_atom_elec_energies(elecalc, &vpmg),
                            });
                        }
                        prev_vpmg = Some(vpmg);
                    }
                    Err(e) => {
                        let msg = format!("ELEC calc '{}' failed: {}", elecalc.name, e);
                        eprintln!("{}", msg);
                        prev_vpmg = None;
                        STDOUT_CAPTURE.with(|capture| {
                            *capture.borrow_mut() = None;
                        });
                        return Err(ApbsError::InvalidParameter(msg));
                    }
                }
            }
            NOshCalc::Apolar(apolarcalc) => {
                println!("Running APOLAR calc: {}", apolarcalc.name);
                match run_apolar_full(&nosh, apolarcalc) {
                    Ok(result) => {
                        println!(
                            "Finished APOLAR calc: {} ({:.6e} kJ/mol)",
                            apolarcalc.name, result.energy
                        );
                        run.calcs.push(SolverCalcResult {
                            name: apolarcalc.name.clone(),
                            kind: SolverCalcKind::Apolar,
                            per_atom: result.atomsasa,
                        });
                    }
                    Err(e) => {
                        let msg = format!("APOLAR calc '{}' failed: {}", apolarcalc.name, e);
                        eprintln!("{}", msg);
                        STDOUT_CAPTURE.with(|capture| {
                            *capture.borrow_mut() = None;
                        });
                        return Err(ApbsError::InvalidParameter(msg));
                    }
                }
            }
            NOshCalc::Print(_) => {}
        }
    }

    run.log = STDOUT_CAPTURE.with(|capture| {
        capture
            .borrow_mut()
            .take()
            .unwrap_or_default()
            .join("\n")
    });
    Ok(run)
}

fn per_atom_elec_energies(
    elec: &apbs_generic::nosh::NOshElec,
    vpmg: &apbs_mg::vpmg::Vpmg,
) -> Vec<f64> {
    let conv = kt_to_kjmol(elec.pbeparm.temp);
    let num_atoms = vpmg.pbe.alist.number_atoms();
    (0..num_atoms)
        .map(|i| 0.5 * vpmg.qf_atom_energy(i) * conv)
        .collect()
}


// /// Write calculation results to XML file.
// /// Port of writedataXML from routines.c line 2336.
// pub fn writedata_xml(
//     fname: &str,
//     nosh: &NOsh,
//     energies: &HashMap<String, f64>,
// ) -> ApbsResult<()> {
//     use std::fs::File;
//     use std::io::Write;

//     let mut file = File::create(fname)
//         .map_err(|e| ApbsError::Io(format!("Failed to open {}: {}", fname, e)))?;

//     writeln!(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
//         .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//     writeln!(file, "<APBS>")
//         .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//     let now = std::time::SystemTime::now()
//         .duration_since(std::time::UNIX_EPOCH)
//         .unwrap_or_default();
//     writeln!(file, "  <date>Timestamp: {} seconds since epoch</date>", now.as_secs())
//         .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//     for calc in &nosh.calcs {
//         if let NOshCalc::Elec(elec) = calc {
//             let pbeparm = &elec.pbeparm;
//             let mgparm = &elec.mgparm;
//             let conversion = apbs_generic::vunit::KB * pbeparm.temp * 1e-3 * apbs_generic::vunit::NA;

//             writeln!(file, "  <elec>")
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <name>{}</name>", elec.name)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::mgparm::MGparmCalcType;
//             let mgtype = match mgparm.r#type {
//                 MGparmCalcType::Dummy => "mg-dummy",
//                 MGparmCalcType::Manual => "mg-manual",
//                 MGparmCalcType::Auto => "mg-auto",
//                 MGparmCalcType::Parallel => "mg-para",
//                 _ => "unknown",
//             };
//             writeln!(file, "    <type>{}</type>", mgtype)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <molid>{}</molid>", pbeparm.molid.unwrap_or(0))
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <nx>{}</nx>", mgparm.dime[0])
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <ny>{}</ny>", mgparm.dime[1])
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <nz>{}</nz>", mgparm.dime[2])
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::vhal::VhalPBEType;
//             let pbetyp = match pbeparm.pbetype {
//                 VhalPBEType::NPBE => "npbe",
//                 VhalPBEType::LPBE => "lpbe",
//                 _ => "unknown",
//             };
//             writeln!(file, "    <pbe>{}</pbe>", pbetyp)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             for i in 0..pbeparm.nion as usize {
//                 writeln!(file, "    <ion>")
//                     .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                 writeln!(file, "      <radius>{}</radius>", pbeparm.ionr[i])
//                     .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                 writeln!(file, "      <charge>{}</charge>", pbeparm.ionq[i])
//                     .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                 writeln!(file, "      <concentration>{}</concentration>", pbeparm.ionc[i])
//                     .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                 writeln!(file, "    </ion>")
//                     .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             }

//             writeln!(file, "    <pdie>{}</pdie>", pbeparm.pdie)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <sdie>{}</sdie>", pbeparm.sdie)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::vhal::VsurfMeth;
//             let srfmtxt = match pbeparm.srfm {
//                 VsurfMeth::Mol => "mol",
//                 VsurfMeth::MolSmooth => "smol",
//                 VsurfMeth::Spline => "spl2",
//                 VsurfMeth::Spline3 => "spl3",
//                 VsurfMeth::Spline4 => "spl4",
//                 _ => "unknown",
//             };
//             writeln!(file, "    <srfm>{}</srfm>", srfmtxt)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <srad>{}</srad>", pbeparm.srad)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             use apbs_generic::vhal::Vbcfl;
//             let bcfltxt = match pbeparm.bcfl {
//                 Vbcfl::Zero => "zero",
//                 Vbcfl::SDH => "sdh",
//                 Vbcfl::MDH => "mdh",
//                 Vbcfl::Focus => "focus",
//                 Vbcfl::Map => "map",
//                 Vbcfl::Mem => "mem",
//                 _ => "unknown",
//             };
//             writeln!(file, "    <bcfl>{}</bcfl>", bcfltxt)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//             writeln!(file, "    <temp>{}</temp>", pbeparm.temp)
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//             if let Some(&energy) = energies.get(&elec.name) {
//                 use apbs_generic::pbeparm::PBEparmCalcEnergy;
//                 match pbeparm.calc_energy {
//                     PBEparmCalcEnergy::Total => {
//                         writeln!(file, "    <totEnergy>{:1.12E}</totEnergy>", energy * conversion)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                     PBEparmCalcEnergy::Comps => {
//                         writeln!(file, "    <totEnergy>{:1.12E}</totEnergy>", energy * conversion)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                     _ => {}
//                 }
//             }

//             writeln!(file, "  </elec>")
//                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//         }
//     }

//     // Handle PRINT statements
//     for calc in &nosh.calcs {
//         if let NOshCalc::Print(print_stmt) = calc {
//             match print_stmt {
//                 NOshPrint::ElecEnergy { left, right } => {
//                     if right.is_empty() {
//                         if let Some(&energy) = energies.get(left) {
//                             let temp = nosh.get_elec_temp(left).unwrap_or(298.15);
//                             let conversion = apbs_generic::vunit::KB * temp * 1e-3 * apbs_generic::vunit::NA;
//                             writeln!(file, "  <printEnergy>")
//                                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                             writeln!(file, "    <equation>{} ({}) end</equation>", left, left)
//                                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                             writeln!(file, "    <globalEnergy>{:1.12E}</globalEnergy>", energy * conversion)
//                                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                             writeln!(file, "  </printEnergy>")
//                                 .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         }
//                     } else {
//                         let e_left = energies.get(left).copied().unwrap_or(0.0);
//                         let e_right = energies.get(right).copied().unwrap_or(0.0);
//                         let temp = nosh.get_elec_temp(left).unwrap_or(298.15);
//                         let conversion = apbs_generic::vunit::KB * temp * 1e-3 * apbs_generic::vunit::NA;
//                         writeln!(file, "  <printEnergy>")
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "    <equation>{} ({}) - {} ({}) end</equation>", left, left, right, right)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "    <globalEnergy>{:1.12E}</globalEnergy>", (e_left - e_right) * conversion)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "  </printEnergy>")
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                 }
//                 NOshPrint::ApolEnergy { name } => {
//                     if let Some(&energy) = energies.get(name) {
//                         writeln!(file, "  <printEnergy>")
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "    <equation>APOL energy {} ({}) end</equation>", name, name)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "    <globalEnergy>{:1.12E}</globalEnergy>", energy)
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                         writeln!(file, "  </printEnergy>")
//                             .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;
//                     }
//                 }
//                 _ => {}
//             }
//         }
//     }

//     writeln!(file, "</APBS>")
//         .map_err(|e| ApbsError::Io(format!("Write error: {}", e)))?;

//     Ok(())
// }
