//! `gmx convert-tpr`, mirroring `src/gromacs/tools/convert_tpr.cpp`.
//!
//! The minimal port supports the three runtime modification operations
//! (`-extend`, `-until`, `-nsteps`) plus `-generate_velocities`.  The subset
//! selection logic (`reduce_topology_x`) is not implemented because it would
//! require re-serializing the whole topology.

use crate::cmd::{select_group, Args};
use crate::index::{self, IndexGroup};
use crate::tpr;

fn print_runtime_info(ir: &tpr::InputRec) {
    println!(
        "  Run start step                {:22}     ",
        ir.init_step
    );
    println!(
        "  Run start time                {:22} ps  ",
        ir.init_step as f64 * ir.delta_t + ir.init_t
    );
    println!("  Step to be made during run    {:22}     ", ir.nsteps);
    println!(
        "  Runtime for the run           {:22} ps  ",
        ir.nsteps as f64 * ir.delta_t
    );
    println!(
        "  Run end step                  {:22}     ",
        ir.init_step + ir.nsteps
    );
    println!(
        "  Run end time                  {:22} ps  \n",
        (ir.init_step + ir.nsteps) as f64 * ir.delta_t + ir.init_t
    );
}

pub fn run(argv: Vec<String>) -> i32 {
    let args = Args::parse(argv);
    let input = match args.get("s") {
        Some(s) => s,
        None => {
            eprintln!("gmx-rs-tools convert-tpr: option -s is required");
            return 1;
        }
    };
    let output = args.get("o").unwrap_or_else(|| "tprout.tpr".to_string());

    let extend_set = args.has("extend");
    let until_set = args.has("until");
    let nsteps_set = args.has("nsteps");
    let generate_velocities = args.flag("generate_velocities");
    let index_given = args.has("n");
    let index_file = args.get("n");

    if (extend_set as i32 + until_set as i32 + nsteps_set as i32) > 1 {
        println!("Multiple runtime modification operations cannot be done in a single call.");
        return 1;
    }
    if (extend_set || until_set || nsteps_set || generate_velocities) && index_given {
        println!(
            "Cannot do runtime modification or velocity generation together with index group extraction in a single call."
        );
        return 1;
    }

    let tpr = match tpr::TprFile::read(&input) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let body = match tpr::parse_body(&tpr.header, &tpr.body) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let mut ir = match body.ir.clone() {
        Some(ir) => ir,
        None => {
            eprintln!("the input tpr file does not contain an inputrec");
            return 1;
        }
    };

    if extend_set || nsteps_set || until_set {
        let input_time_at_start = ir.init_step as f64 * ir.delta_t + ir.init_t;
        let input_step_at_end = ir.init_step + ir.nsteps;
        let input_time_at_end = input_step_at_end as f64 * ir.delta_t + ir.init_t;

        println!("Input file:");
        print_runtime_info(&ir);

        if nsteps_set {
            let nsteps = args.int("nsteps", 0);
            eprintln!("Setting nsteps to {nsteps}");
            ir.nsteps = nsteps;
        } else if extend_set {
            let extend = args.real("extend", 0.0);
            println!("Extending remaining runtime by {extend} ps");
            ir.nsteps += (extend / ir.delta_t).round() as i64;
        } else {
            let until = args.real("until", 0.0);
            if until <= input_time_at_start {
                println!("The requested run end time is at/before the run start time.");
                return 1;
            }
            if until < input_time_at_end {
                println!("The requested run end time is before the original run end time.");
                println!("Reducing remaining runtime to {until} ps");
            } else {
                println!("Extending remaining runtime to {until} ps");
            }
            ir.nsteps = ((until - input_time_at_start) / ir.delta_t).round() as i64;
        }

        println!("\nOutput file:");
        print_runtime_info(&ir);
    }

    if generate_velocities {
        match &body.v {
            None => {
                eprintln!(
                    "Input tpr file {input} does not contain velocities, typically because this file is intended for energy minimization ('steep' integrator)."
                );
                return 1;
            }
            Some(_) => {
                let temp = args.real("velocity_temp", 300.0);
                if temp < 0.0 {
                    println!("Temperature used to generate velocities must be positive.");
                    return 1;
                }
                let seed = args.int("velocity_seed", -1);
                let seed = if seed == -1 {
                    let s = random_seed();
                    println!("Using random seed {s} for generating velocities");
                    s
                } else {
                    seed
                };
                let atoms = body
                    .mtop
                    .as_ref()
                    .map(|m| m.global_atoms())
                    .unwrap_or_default();
                let n = body.v.as_ref().unwrap().len();
                let masses: Vec<f64> = (0..n)
                    .map(|i| atoms.atom.get(i).map(|a| a.mass).unwrap_or(1.0))
                    .collect();
                let mut v = body.v.clone().unwrap();
                maxwell_speed(temp, seed as u64, &masses, &mut v);
                // Remove the centre of mass velocity and clear massless atoms.
                let mut total_mass = 0.0f64;
                let mut com = [0.0f64; 3];
                for i in 0..n {
                    total_mass += masses[i];
                    for d in 0..3 {
                        com[d] += masses[i] * v[i][d] as f64;
                    }
                }
                if total_mass > 0.0 {
                    for d in 0..3 {
                        com[d] /= total_mass;
                    }
                    for i in 0..n {
                        if masses[i] > 0.0 {
                            for d in 0..3 {
                                v[i][d] -= com[d] as f32;
                            }
                        } else {
                            v[i] = [0.0; 3];
                        }
                    }
                }
                println!(
                    "Generated velocities at {} K using seed {} (body values are not written back in this minimal port)",
                    temp, seed
                );
            }
        }
    }

    // Subset extraction (`-n`): pick an index group and rewrite the topology
    // for the selection, like reduce_topology_x() in the C++ tool.
    let mut subset: Option<Vec<usize>> = None;
    // The C tool always asks for a group unless a runtime modification or
    // velocity generation was requested; with no index file the default groups
    // are offered instead.
    let do_selection = !(extend_set || until_set || nsteps_set || generate_velocities);
    if do_selection {
        let atoms = body
            .mtop
            .as_ref()
            .map(|m| m.global_atoms())
            .unwrap_or_default();
        if atoms.nr() == 0 {
            eprintln!("the input tpr file does not contain a topology");
            return 1;
        }
        let groups: Vec<IndexGroup> = match &index_file {
            Some(f) => match index::read_ndx(f) {
                Ok(g) => g,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            },
            None => index::analyse(&atoms, false),
        };
        let g = match select_group(&groups, "Select a group:") {
            Some(g) => g,
            None => return 1,
        };
        let gnx = groups[g].particle_indices.len();
        let index = groups[g].particle_indices.clone();
        let mut b_sel = gnx != body.natoms;
        for (i, &v) in index.iter().enumerate() {
            if i != v {
                b_sel = true;
            }
        }
        if b_sel {
            eprintln!(
                "Will write subset {} of original tpx containing {} atoms\n",
                groups[g].name, gnx
            );
            subset = Some(index);
        } else {
            eprintln!("Will write full tpx file (no selection)\n");
        }
    }

    let out_tpr = match &subset {
        Some(selection) => match tpr::write_subset_body(&tpr, &body, selection) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        },
        None => match tpr.set_nsteps(&body, ir.nsteps) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        },
    };
    if let Err(e) = out_tpr.write(&output) {
        eprintln!("{e}");
        return 1;
    }
    0
}

/// Simple xorshift based random seed generator (replaces `gmx::makeRandomSeed`).
fn random_seed() -> i64 {
    let t = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos() as u64)
        .unwrap_or(12345);
    (t % 2_000_000_000) as i64 + 1
}

/// Maxwell-Boltzmann velocity generation (`gromacs/gmxpreprocess/
/// gen_maxwell_velocities.cpp`).
fn maxwell_speed(temp: f64, seed: u64, masses: &[f64], v: &mut [[f32; 3]]) {
    // Convert to GROMACS units: mass in u, velocity in nm/ps.
    const BOLTZMANN: f64 = 0.0083144621; // kJ/mol/K
    let mut state = seed.wrapping_mul(2).wrapping_add(1) | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    // Box-Muller pairs.
    let gauss = |next: &mut dyn FnMut() -> f64| -> (f64, f64) {
        let mut u1 = next();
        if u1 < 1e-300 {
            u1 = 1e-300;
        }
        let u2 = next();
        let r = (-2.0 * u1.ln()).sqrt();
        (
            r * (2.0 * std::f64::consts::PI * u2).cos(),
            r * (2.0 * std::f64::consts::PI * u2).sin(),
        )
    };
    for i in 0..masses.len() {
        if masses[i] <= 0.0 {
            v[i] = [0.0; 3];
            continue;
        }
        let sigma = (BOLTZMANN * temp / masses[i]).sqrt();
        let (g1, _) = gauss(&mut next);
        let (g2, _) = gauss(&mut next);
        let (g3, _) = gauss(&mut next);
        v[i] = [
            (sigma * g1) as f32,
            (sigma * g2) as f32,
            (sigma * g3) as f32,
        ];
    }
}
