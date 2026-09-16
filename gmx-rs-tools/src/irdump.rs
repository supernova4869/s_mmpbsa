//! Printing of a decoded inputrec, mirroring `pr_inputrec()` and the static
//! helpers in `mdtypes/inputrec.cpp`, plus `dumpKeyValueTree()` from
//! `utility/keyvaluetree.cpp` and the option defaults of the `applied-forces`
//! and `fast-multipole-method` MD modules (`mdrun/mdmodules.cpp`,
//! `applied_forces/*`, `fmm/*`).

use crate::cmd::{fmt_e, fmt_g};
use crate::enum_names as en;
use crate::enum_names::lookup;
use crate::ir::{
    AwhBiasParams, AwhParams, ExpandedVals, FepVals, InputRec, KvtValue, PullParams, Rot, SimTempVals,
    SwapCoords,
};
use std::fmt::Write;

const INDENT: usize = 3;

fn indent(out: &mut String, n: usize) {
    for _ in 0..n {
        out.push(' ');
    }
}

/// `pr_title()`: writes `title:` and returns the new indent.
fn title(out: &mut String, n: usize, text: &str) -> usize {
    indent(out, n);
    let _ = writeln!(out, "{text}:");
    n + INDENT
}

fn p_str(out: &mut String, n: usize, key: &str, value: &str) {
    indent(out, n);
    let _ = writeln!(out, "{key:<30} = {value}");
}

fn p_int(out: &mut String, n: usize, key: &str, value: i32) {
    indent(out, n);
    let _ = writeln!(out, "{key:<30} = {value}");
}

fn p_int64(out: &mut String, n: usize, key: &str, value: i64) {
    indent(out, n);
    let _ = writeln!(out, "{key:<30} = {value}");
}

fn p_real(out: &mut String, n: usize, key: &str, value: f64) {
    indent(out, n);
    let _ = writeln!(out, "{key:<30} = {}", fmt_g(value, 0, 6));
}

fn p_bool(out: &mut String, n: usize, key: &str, value: bool) {
    p_str(out, n, key, if value { "true" } else { "false" });
}

/// `pr_rvecs()`: `title (nx3):` followed by one `%12.5e` triple per row.
fn p_rvecs(out: &mut String, n: usize, key: &str, values: &[[f64; 3]]) {
    let inner = title_width(out, n, key, values.len(), true);
    for (i, v) in values.iter().enumerate() {
        indent(out, inner);
        let _ = writeln!(
            out,
            "{key}[{i:5}]={{{}, {}, {}}}",
            fmt_e(v[0], 12, 5),
            fmt_e(v[1], 12, 5),
            fmt_e(v[2], 12, 5)
        );
    }
}

/// `pr_title_nxn()` with `nx3` and the entries following.
fn title_width(out: &mut String, n: usize, key: &str, rows: usize, _dims: bool) -> usize {
    indent(out, n);
    let _ = writeln!(out, "{key} ({rows}x3):");
    n + INDENT
}

/// `pr_rvec()`: `title (n):` followed by `%12.5e` scalars.
fn p_rvec(out: &mut String, n: usize, key: &str, values: &[f64]) {
    if values.is_empty() {
        indent(out, n);
        let _ = writeln!(out, "{key}: not available");
        return;
    }
    indent(out, n);
    let _ = writeln!(out, "{key} ({}):", values.len());
    for (i, v) in values.iter().enumerate() {
        indent(out, n + INDENT);
        let _ = writeln!(out, "{key}[{i}]={}", fmt_e(*v, 12, 5));
    }
}

/// `pr_ivec()`.
fn p_ivec(out: &mut String, n: usize, key: &str, values: &[i32]) {
    indent(out, n);
    let _ = writeln!(out, "{key} ({}):", values.len());
    for (i, v) in values.iter().enumerate() {
        indent(out, n + INDENT);
        let _ = writeln!(out, "{key}[{i}]={v}");
    }
}

/// `pr_ivec_block()`: consecutive runs of at least three indices are
/// compressed into a range.
fn p_ivec_block(out: &mut String, n: usize, key: &str, values: &[i32]) {
    if values.is_empty() {
        indent(out, n);
        let _ = writeln!(out, "{key}: not available");
        return;
    }
    indent(out, n);
    let _ = writeln!(out, "{key} ({}):", values.len());
    let inner = n + INDENT;
    let mut i = 0usize;
    while i < values.len() {
        let mut j = i + 1;
        while j < values.len() && values[j] == values[j - 1] + 1 {
            j += 1;
        }
        if j - i < 3 {
            while i < j {
                indent(out, inner);
                let _ = writeln!(out, "{key}[{i}]={}", values[i]);
                i += 1;
            }
        } else {
            indent(out, inner);
            let _ = writeln!(
                out,
                "{key}[{},...,{}] = {{{},...,{}}}",
                i,
                j - 1,
                values[i],
                values[j - 1]
            );
            i = j;
        }
    }
}

/// `pr_matrix()` in the non-MDP format.
fn p_matrix(out: &mut String, n: usize, key: &str, m: &[[f64; 3]; 3]) {
    p_rvecs(out, n, key, m);
}

/// Prints the whole `inputrec:` block.
///
/// `original_inputrec` mirrors `gmx dump -orgir`: the parameters are printed
/// as they were read from the file instead of being adjusted with the defaults
/// of the MD modules of the installed GROMACS.
pub fn print_inputrec(out: &mut String, ir: &InputRec, original_inputrec: bool) {
    let indent = title(out, 0, "inputrec");

    p_str(out, indent, "integrator", lookup(en::INTEGRATION_ALGORITHM, ir.integrator as i64));
    p_real(out, indent, "tinit", ir.init_t);
    p_real(out, indent, "dt", ir.delta_t);
    p_int64(out, indent, "nsteps", ir.nsteps);
    p_int64(out, indent, "init-step", ir.init_step);
    p_int(out, indent, "simulation-part", ir.simulation_part);
    p_bool(out, indent, "mts", ir.use_mts);
    if ir.use_mts {
        for (index, level) in ir.mts_levels.iter().enumerate().skip(1) {
            let mut groups = String::new();
            for (i, name) in [
                "longrange-nonbonded",
                "nonbonded",
                "pair",
                "dihedral",
                "angle",
                "pull",
                "awh",
            ]
            .iter()
            .enumerate()
            {
                if level.force_groups & (1 << i) != 0 {
                    if !groups.is_empty() {
                        groups.push(' ');
                    }
                    groups.push_str(name);
                }
            }
            let key = format!("mts-level{}-forces", index + 1);
            p_str(out, indent, &key, &groups);
            let key = format!("mts-level{}-factor", index + 1);
            p_int(out, indent, &key, level.step_factor);
        }
    }
    p_real(out, indent, "mass-repartition-factor", ir.mass_repartition_factor);
    p_str(
        out,
        indent,
        "comm-mode",
        lookup(en::COM_REMOVAL_ALGORITHM, ir.comm_mode as i64),
    );
    p_int(out, indent, "nstcomm", ir.nstcomm);

    p_real(out, indent, "bd-fric", ir.bd_fric);
    p_int64(out, indent, "ld-seed", ir.ld_seed);

    p_real(out, indent, "emtol", ir.em_tol);
    p_real(out, indent, "emstep", ir.em_stepsize);
    p_int(out, indent, "niter", ir.niter);
    p_real(out, indent, "fcstep", ir.fc_stepsize);
    p_int(out, indent, "nstcgsteep", ir.nstcgsteep);
    p_int(out, indent, "nbfgscorr", ir.nbfgscorr);

    p_real(out, indent, "rtpi", ir.rtpi);

    p_int(out, indent, "nstxout", ir.nstxout);
    p_int(out, indent, "nstvout", ir.nstvout);
    p_int(out, indent, "nstfout", ir.nstfout);
    p_int(out, indent, "nstlog", ir.nstlog);
    p_int(out, indent, "nstcalcenergy", ir.nstcalcenergy);
    p_int(out, indent, "nstenergy", ir.nstenergy);
    p_int(out, indent, "nstxout-compressed", ir.nstxout_compressed);
    p_real(out, indent, "compressed-x-precision", ir.x_compression_precision);

    p_str(out, indent, "cutoff-scheme", lookup(en::CUTOFF_SCHEME, ir.cutoff_scheme as i64));
    p_int(out, indent, "nstlist", ir.nstlist);
    p_str(out, indent, "pbc", ir.pbc().name());
    p_bool(out, indent, "periodic-molecules", ir.b_periodic_mols);
    p_real(out, indent, "verlet-buffer-tolerance", ir.verletbuf_tol);
    p_real(
        out,
        indent,
        "verlet-buffer-pressure-tolerance",
        ir.verlet_buffer_pressure_tolerance,
    );
    p_real(out, indent, "rlist", ir.rlist);

    p_str(
        out,
        indent,
        "coulombtype",
        lookup(en::COULOMB_INTERACTION_TYPE, ir.coulombtype as i64),
    );
    p_str(
        out,
        indent,
        "coulomb-modifier",
        lookup(en::INTERACTION_MODIFIERS, ir.coulomb_modifier as i64),
    );
    p_real(out, indent, "rcoulomb-switch", ir.rcoulomb_switch);
    p_real(out, indent, "rcoulomb", ir.rcoulomb);
    if ir.epsilon_r != 0.0 {
        p_real(out, indent, "epsilon-r", ir.epsilon_r);
    } else {
        p_str(out, indent, "epsilon-r", "inf");
    }
    if ir.epsilon_rf != 0.0 {
        p_real(out, indent, "epsilon-rf", ir.epsilon_rf);
    } else {
        p_str(out, indent, "epsilon-rf", "inf");
    }
    p_str(out, indent, "vdw-type", lookup(en::VAN_DER_WAALS_TYPE, ir.vdwtype as i64));
    p_str(
        out,
        indent,
        "vdw-modifier",
        lookup(en::INTERACTION_MODIFIERS, ir.vdw_modifier as i64),
    );
    p_real(out, indent, "rvdw-switch", ir.rvdw_switch);
    p_real(out, indent, "rvdw", ir.rvdw);
    p_str(
        out,
        indent,
        "DispCorr",
        lookup(en::DISPERSION_CORRECTION_TYPE, ir.disp_corr as i64),
    );
    p_real(out, indent, "table-extension", ir.tabext);

    p_real(out, indent, "fourierspacing", ir.fourier_spacing);
    p_int(out, indent, "fourier-nx", ir.nkx);
    p_int(out, indent, "fourier-ny", ir.nky);
    p_int(out, indent, "fourier-nz", ir.nkz);
    p_int(out, indent, "pme-order", ir.pme_order);
    p_real(out, indent, "ewald-rtol", ir.ewald_rtol);
    p_real(out, indent, "ewald-rtol-lj", ir.ewald_rtol_lj);
    p_str(
        out,
        indent,
        "lj-pme-comb-rule",
        lookup(en::LONG_RANGE_VD_W, ir.ljpme_combination_rule as i64),
    );
    p_str(
        out,
        indent,
        "ewald-geometry",
        lookup(en::EWALD_GEOMETRY, ir.ewald_geometry as i64),
    );
    p_real(out, indent, "epsilon-surface", ir.epsilon_surface);

    p_str(
        out,
        indent,
        "ensemble-temperature-setting",
        lookup(en::ENSEMBLE_TEMPERATURE_SETTING, ir.ensemble_temperature_setting as i64),
    );
    if ir.ensemble_temperature_setting == 1 {
        p_real(out, indent, "ensemble-temperature", ir.ensemble_temperature);
    }
    p_str(out, indent, "tcoupl", lookup(en::TEMPERATURE_COUPLING, ir.tcoupl as i64));
    p_int(out, indent, "nsttcouple", ir.nsttcouple);
    p_int(out, indent, "nh-chain-length", ir.opts.nhchainlength);
    p_bool(out, indent, "print-nose-hoover-chain-variables", ir.b_print_nh_chains);

    p_str(out, indent, "pcoupl", lookup(en::PRESSURE_COUPLING, ir.pcoupl as i64));
    if ir.pcoupl != 0 {
        p_str(
            out,
            indent,
            "pcoupltype",
            lookup(en::PRESSURE_COUPLING_TYPE, ir.pcoupltype as i64),
        );
        p_int(out, indent, "nstpcouple", ir.nstpcouple);
        p_real(out, indent, "tau-p", ir.tau_p);
        p_matrix(out, indent, "compressibility", &ir.compress);
        p_matrix(out, indent, "ref-p", &ir.ref_p);
    }
    p_str(
        out,
        indent,
        "refcoord-scaling",
        lookup(en::REF_COORD_SCALING, ir.refcoord_scaling as i64),
    );

    // `prRVecs()` prints an empty `(0x3)` block when the vector was allocated
    // but cleared, and "not available" when it was never allocated.
    if ir.posres_com_allocated {
        p_rvecs(out, indent, "posres-com", &ir.posres_com);
        p_rvecs(out, indent, "posres-comB", &ir.posres_com_b);
    } else {
        out.push_str("   posres-com: not available\n");
        out.push_str("   posres-comB: not available\n");
    }

    p_bool(out, indent, "QMMM", ir.opts.b_qmmm);
    out.push_str("qm-opts:\n");
    p_int(out, indent, "ngQM", ir.opts.ngqm);

    p_str(
        out,
        indent,
        "constraint-algorithm",
        lookup(en::CONSTRAINT_ALGORITHM, ir.e_constr_alg as i64),
    );
    p_bool(out, indent, "continuation", ir.b_continuation);

    p_bool(out, indent, "Shake-SOR", ir.b_shake_sor);
    p_real(out, indent, "shake-tol", ir.shake_tol);
    p_int(out, indent, "lincs-order", ir.n_proj_order);
    p_int(out, indent, "lincs-iter", ir.n_lincs_iter);
    p_real(out, indent, "lincs-warnangle", ir.lincs_warn_angle);

    p_int(out, indent, "nwall", ir.nwall);
    p_str(out, indent, "wall-type", lookup(en::WALL_TYPE, ir.wall_type as i64));
    p_real(out, indent, "wall-r-linpot", ir.wall_r_linpot);
    p_int(out, indent, "wall-atomtype[0]", ir.wall_atomtype[0]);
    p_int(out, indent, "wall-atomtype[1]", ir.wall_atomtype[1]);
    p_real(out, indent, "wall-density[0]", ir.wall_density[0]);
    p_real(out, indent, "wall-density[1]", ir.wall_density[1]);
    p_real(out, indent, "wall-ewald-zfac", ir.wall_ewald_zfac);

    p_bool(out, indent, "pull", ir.b_pull);
    if ir.b_pull {
        print_pull(out, indent, &ir.pull);
    }

    p_bool(out, indent, "awh", ir.b_do_awh);
    if ir.b_do_awh {
        print_awh(out, indent, &ir.awh);
    }

    p_bool(out, indent, "rotation", ir.b_rot);
    if ir.b_rot {
        print_rot(out, indent, &ir.rot);
    }

    p_bool(out, indent, "interactiveMD", ir.b_imd);
    if ir.b_imd {
        print_imd(out, indent, &ir.imd.ind);
    }

    p_str(
        out,
        indent,
        "disre",
        lookup(en::DISTANCE_RESTRAINT_REFINEMENT, ir.e_disre as i64),
    );
    p_str(
        out,
        indent,
        "disre-weighting",
        lookup(en::DISTANCE_RESTRAINT_WEIGHTING, ir.e_disre_weighting as i64),
    );
    p_bool(out, indent, "disre-mixed", ir.b_disre_mixed);
    p_real(out, indent, "dr-fc", ir.dr_fc);
    p_real(out, indent, "dr-tau", ir.dr_tau);
    p_int(out, indent, "nstdisreout", ir.nstdisreout);

    p_real(out, indent, "orire-fc", ir.orires_fc);
    p_real(out, indent, "orire-tau", ir.orires_tau);
    p_int(out, indent, "nstorireout", ir.nstorireout);

    p_str(
        out,
        indent,
        "free-energy",
        lookup(en::FREE_ENERGY_PERTURBATION_TYPE, ir.efep as i64),
    );
    if ir.efep != 0 || ir.b_sim_temp {
        print_fepvals(out, indent, &ir.fepvals);
    }
    if ir.b_expanded {
        print_expandedvals(out, indent, &ir.expandedvals, ir.fepvals.n_lambda);
    }

    p_real(out, indent, "cos-acceleration", ir.cos_accel);
    p_matrix(out, indent, "deform", &ir.deform);

    p_bool(out, indent, "simulated-tempering", ir.b_sim_temp);
    if ir.b_sim_temp {
        print_simtempvals(out, indent, &ir.simtempvals, ir.fepvals.n_lambda);
    }

    p_str(
        out,
        indent,
        "swapcoords",
        lookup(en::SWAP_TYPE, ir.e_swap_coords as i64),
    );
    if ir.e_swap_coords != 0 {
        print_swap(out, indent, &ir.swap);
    }

    p_int(out, indent, "userint1", ir.userint1);
    p_int(out, indent, "userint2", ir.userint2);
    p_int(out, indent, "userint3", ir.userint3);
    p_int(out, indent, "userint4", ir.userint4);
    p_real(out, indent, "userreal1", ir.userreal1);
    p_real(out, indent, "userreal2", ir.userreal2);
    p_real(out, indent, "userreal3", ir.userreal3);
    p_real(out, indent, "userreal4", ir.userreal4);

    dump_params(out, indent, &ir.params, !original_inputrec);

    print_grp_opts(out, indent, "grpopts", ir);
}

/// `pr_grp_opts()` for the non-MDP format, with the inputs that come from the
/// inputrec itself (QMMM and the acceleration groups).
fn print_grp_opts(out: &mut String, n: usize, key: &str, ir: &InputRec) {
    let opts = &ir.opts;
    let _ = writeln!(out, "{key}:");

    indent(out, n);
    let _ = write!(out, "nrdf:");
    for v in &opts.nrdf {
        let _ = write!(out, "  {}", fmt_g(*v, 10, 6));
    }
    out.push('\n');

    indent(out, n);
    let _ = write!(out, "ref-t:");
    for v in &opts.ref_t {
        let _ = write!(out, "  {}", fmt_g(*v, 10, 6));
    }
    out.push('\n');

    indent(out, n);
    let _ = write!(out, "tau-t:");
    for v in &opts.tau_t {
        let _ = write!(out, "  {}", fmt_g(*v, 10, 6));
    }
    out.push('\n');

    let _ = write!(out, "annealing:");
    for v in &opts.annealing {
        let _ = write!(out, "  {}", pad_left(lookup(en::SIMULATED_ANNEALING, *v as i64), 10));
    }
    out.push('\n');

    let _ = write!(out, "annealing-npoints:");
    for v in &opts.anneal_npoints {
        let _ = write!(out, "  {}", pad_left(&v.to_string(), 10));
    }
    out.push('\n');

    for (i, npoints) in opts.anneal_npoints.iter().enumerate() {
        if *npoints > 0 {
            let _ = write!(out, "annealing-time [{i}]:\t");
            for v in opts.anneal_time.get(i).map(|v| v.as_slice()).unwrap_or(&[]) {
                let _ = write!(out, "  {:>10.1}", v);
            }
            out.push('\n');
            let _ = write!(out, "annealing-temp [{i}]:\t");
            for v in opts.anneal_temp.get(i).map(|v| v.as_slice()).unwrap_or(&[]) {
                let _ = write!(out, "  {:>10.1}", v);
            }
            out.push('\n');
        }
    }

    indent(out, n);
    let _ = write!(out, "acc:\t");
    for a in &opts.acceleration {
        for v in a {
            let _ = write!(out, "  {}", fmt_g(*v, 10, 6));
        }
    }
    out.push('\n');

    indent(out, n);
    let _ = write!(out, "nfreeze:");
    for f in &opts.nfreeze {
        for v in f {
            let _ = write!(out, "  {}", pad_left(if *v != 0 { "Y" } else { "N" }, 10));
        }
    }
    out.push('\n');

    for i in 0..opts.ngener {
        indent(out, n);
        let _ = write!(out, "energygrp-flags[{i:3}]:");
        for m in 0..opts.ngener {
            let idx = opts.ngener as usize * i as usize + m as usize;
            let v = opts.egp_flags.get(idx).copied().unwrap_or(0);
            let _ = write!(out, " {v}");
        }
        out.push('\n');
    }
}

/// `%*s` with a positive width.
fn pad_left(s: &str, width: usize) -> String {
    if s.len() < width {
        format!("{s:>width$}")
    } else {
        s.to_string()
    }
}

fn print_fepvals(out: &mut String, n: usize, fep: &FepVals) {
    p_real(out, n, "init-lambda", fep.init_lambda_without_states);
    p_int(out, n, "init-lambda-state", fep.init_fep_state);
    p_real(out, n, "delta-lambda", fep.delta_lambda);
    p_int(out, n, "nstdhdl", fep.nstdhdl);
    p_int(out, n, "n-lambdas", fep.n_lambda);
    if fep.n_lambda > 0 {
        indent(out, n);
        out.push_str("separate-dvdl:\n");
        for (i, name) in en::FREE_ENERGY_PERTURBATION_COUPLING_TYPE
            .iter()
            .enumerate()
        {
            let v = fep.separate_dvdl.get(i).copied().unwrap_or(false);
            let _ = write!(out, "{name:>18} = ");
            out.push_str(if v { "  TRUE" } else { "  FALSE" });
            out.push('\n');
        }
        out.push_str("all-lambdas:\n");
        for (i, name) in en::FREE_ENERGY_PERTURBATION_COUPLING_TYPE
            .iter()
            .enumerate()
        {
            let _ = write!(out, "{name:>18} = ");
            for j in 0..fep.n_lambda as usize {
                let v = fep.all_lambda[i].get(j).copied().unwrap_or(0.0);
                let _ = write!(out, "  {}", fmt_g(v, 10, 6));
            }
            out.push('\n');
        }
    }
    p_int(out, n, "calc-lambda-neighbors", fep.lambda_neighbors);
    p_str(
        out,
        n,
        "dhdl-print-energy",
        lookup(en::FREE_ENERGY_PRINT_ENERGY, fep.edhdl_print_energy as i64),
    );
    p_real(out, n, "sc-alpha", fep.sc_alpha);
    p_int(out, n, "sc-power", fep.sc_power);
    p_real(out, n, "sc-r-power", fep.sc_r_power);
    p_real(out, n, "sc-sigma", fep.sc_sigma);
    p_real(out, n, "sc-sigma-min", fep.sc_sigma_min);
    p_bool(out, n, "sc-coul", fep.b_sc_coul);
    p_int(out, n, "dh-hist-size", fep.dh_hist_size);
    p_real(out, n, "dh-hist-spacing", fep.dh_hist_spacing);
    p_str(
        out,
        n,
        "separate-dhdl-file",
        lookup(en::SEPARATE_DHDL_FILE, fep.separate_dhdl_file as i64),
    );
    p_str(
        out,
        n,
        "dhdl-derivatives",
        lookup(en::DH_DL_DERIVATIVE_CALCULATION, fep.dhdl_derivatives as i64),
    );
    p_str(
        out,
        n,
        "sc-function",
        lookup(en::SOFTCORE_TYPE, fep.softcore_function as i64),
    );
    p_real(out, n, "sc-gapsys-scale-linpoint-lj", fep.sc_gapsys_scale_linpoint_lj);
    p_real(out, n, "sc-gapsys-scale-linpoint-q", fep.sc_gapsys_scale_linpoint_q);
    p_real(out, n, "sc-gapsys-sigma-lj", fep.sc_gapsys_sigma_lj);
}

fn print_simtempvals(out: &mut String, n: usize, simtemp: &SimTempVals, n_lambda: i32) {
    p_str(
        out,
        n,
        "simulated-tempering-scaling",
        lookup(en::SIMULATED_TEMPERING, simtemp.scale as i64),
    );
    p_real(out, n, "sim-temp-low", simtemp.low);
    p_real(out, n, "sim-temp-high", simtemp.high);
    p_rvec(
        out,
        n,
        "simulated tempering temperatures",
        &simtemp.temperatures[..simtemp.temperatures.len().min(n_lambda.max(0) as usize)],
    );
}

fn print_expandedvals(out: &mut String, n: usize, expand: &ExpandedVals, n_lambda: i32) {
    p_int(out, n, "nstexpanded", expand.nstexpanded);
    p_str(
        out,
        n,
        "lmc-stats",
        lookup(en::LAMBDA_WEIGHT_CALCULATION, expand.elamstats as i64),
    );
    p_str(
        out,
        n,
        "lmc-move",
        lookup(en::LAMBDA_MOVE_CALCULATION, expand.elmcmove as i64),
    );
    p_str(
        out,
        n,
        "lmc-weights-equil",
        lookup(en::LAMBDA_WEIGHT_WILL_REACH_EQUILIBRIUM, expand.elmceq as i64),
    );
    match expand.elmceq {
        1 => p_int(out, n, "weight-equil-number-all-lambda", expand.equil_n_at_lam),
        2 => p_int(out, n, "weight-equil-number-samples", expand.equil_samples),
        3 => p_int(out, n, "weight-equil-number-steps", expand.equil_steps),
        4 => p_real(out, n, "weight-equil-wl-delta", expand.equil_wl_delta),
        5 => p_real(out, n, "weight-equil-count-ratio", expand.equil_ratio),
        _ => {}
    }
    p_int(out, n, "lmc-seed", expand.lmc_seed);
    p_real(out, n, "mc-temperature", expand.mc_temp);
    p_int(out, n, "lmc-repeats", expand.lmc_repeats);
    p_int(out, n, "lmc-gibbsdelta", expand.gibbsdeltalam);
    p_int(out, n, "lmc-forced-nstart", expand.lmc_forced_nstart);
    p_bool(out, n, "symmetrized-transition-matrix", expand.b_symmetrized_t_matrix);
    p_int(out, n, "nst-transition-matrix", expand.nst_tij);
    p_int(out, n, "mininum-var-min", expand.minvarmin);
    p_int(out, n, "weight-c-range", expand.c_range);
    p_real(out, n, "wl-scale", expand.wl_scale);
    p_real(out, n, "wl-ratio", expand.wl_ratio);
    p_real(out, n, "init-wl-delta", expand.init_wl_delta);
    p_bool(out, n, "wl-oneovert", expand.b_wl_oneovert);
    let n_lambda = n_lambda.max(0) as usize;
    p_rvec(
        out,
        n,
        "init-lambda-weights",
        &expand.init_lambda_weights[..expand.init_lambda_weights.len().min(n_lambda)],
    );
    p_rvec(
        out,
        n,
        "init-lambda-counts",
        &expand.init_lambda_counts[..expand.init_lambda_counts.len().min(n_lambda)],
    );
    p_rvec(
        out,
        n,
        "init-wl-histogram-counts",
        &expand.init_wl_histogram_counts[..expand.init_wl_histogram_counts.len().min(n_lambda)],
    );
}

fn print_pull(out: &mut String, n: usize, pull: &PullParams) {
    p_real(out, n, "pull-cylinder-r", pull.cylinder_r);
    p_real(out, n, "pull-constr-tol", pull.constr_tol);
    p_bool(out, n, "pull-print-COM", pull.b_print_com);
    p_bool(out, n, "pull-print-ref-value", pull.b_print_ref_value);
    p_bool(out, n, "pull-print-components", pull.b_print_comp);
    p_int(out, n, "pull-nstxout", pull.nstxout);
    p_int(out, n, "pull-nstfout", pull.nstfout);
    p_bool(
        out,
        n,
        "pull-pbc-ref-prev-step-com",
        pull.b_set_pbc_ref_to_prev_step_com,
    );
    p_bool(out, n, "pull-xout-average", pull.b_xout_average);
    p_bool(out, n, "pull-fout-average", pull.b_fout_average);
    p_int(out, n, "pull-ngroups", pull.ngroup);
    for (g, group) in pull.group.iter().enumerate() {
        indent(out, n);
        let _ = writeln!(out, "pull-group {g}:");
        let inner = n + 2;
        p_ivec_block(out, inner, "atom", &group.ind);
        p_rvec(out, inner, "weight", &group.weight);
        p_int(out, inner, "pbcatom", group.pbcatom);
    }
    p_int(out, n, "pull-ncoords", pull.ncoord);
    for (c, coord) in pull.coord.iter().enumerate() {
        indent(out, n);
        let _ = writeln!(out, "pull-coord {c}:");
        p_str(
            out,
            n,
            "type",
            lookup(en::PULLING_ALGORITHM, coord.etype as i64),
        );
        if coord.etype == 5 {
            p_str(out, n, "potential-provider", &coord.external_potential_provider);
        }
        p_str(
            out,
            n,
            "geometry",
            lookup(en::PULL_GROUP_GEOMETRY, coord.egeom as i64),
        );
        for (g, group) in coord.group.iter().enumerate() {
            p_int(out, n, &format!("group[{g}]"), *group);
        }
        p_ivec(out, n, "dim", &coord.dim);
        p_rvec(out, n, "origin", &coord.origin);
        p_rvec(out, n, "vec", &coord.vec);
        p_bool(out, n, "start", coord.b_start);
        p_real(out, n, "init", coord.init);
        p_real(out, n, "rate", coord.rate);
        p_real(out, n, "k", coord.k);
        p_real(out, n, "kB", coord.kb);
    }
}

fn print_awh_dim(out: &mut String, n: usize, prefix: &str, dim: &crate::ir::AwhDimParams) {
    indent(out, n);
    let n = n + 1;
    let _ = writeln!(out, "{prefix}:");
    p_str(
        out,
        n,
        "coord-provider",
        lookup(en::AWH_COORDINATE_PROVIDER_TYPE, dim.coord_provider as i64),
    );
    p_int(out, n, "coord-index", dim.coord_index + 1);
    p_real(out, n, "start", dim.origin);
    p_real(out, n, "end", dim.end);
    p_real(out, n, "period", dim.period);
    p_real(out, n, "force-constant", dim.force_constant);
    p_real(out, n, "diffusion", dim.diffusion);
    p_real(out, n, "cover-diameter", dim.cover_diameter);
}

fn print_awh_bias(out: &mut String, n: usize, bias: &AwhBiasParams, prefix: &str) {
    p_real(out, n, &format!("{prefix}-error-init"), bias.error_init);
    p_str(
        out,
        n,
        &format!("{prefix}-growth"),
        lookup(en::AWH_HISTOGRAM_GROWTH_TYPE, bias.e_growth as i64),
    );
    p_real(out, n, &format!("{prefix}-growth-factor"), bias.growth_factor);
    p_str(
        out,
        n,
        &format!("{prefix}-target"),
        lookup(en::AWH_TARGET_TYPE, bias.e_target as i64),
    );
    p_real(
        out,
        n,
        &format!("{prefix}-target-beta-scaling"),
        bias.target_beta_scaling,
    );
    p_real(out, n, &format!("{prefix}-target-cutoff"), bias.target_cutoff);
    p_bool(
        out,
        n,
        &format!("{prefix}-target-metric-scaling"),
        bias.scale_target_by_metric,
    );
    p_real(
        out,
        n,
        &format!("{prefix}-target-metric-scaling-limit"),
        bias.target_metric_scaling_limit,
    );
    p_bool(out, n, &format!("{prefix}-user-data"), bias.b_user_data);
    p_int(out, n, &format!("{prefix}-share-group"), bias.share_group);
    p_bool(
        out,
        n,
        &format!("{prefix}-equilibrate-histogram"),
        bias.equilibrate_histogram,
    );
    p_real(
        out,
        n,
        &format!("{prefix}-histogram-tolerance"),
        bias.histogram_tolerance,
    );
    p_int(out, n, &format!("{prefix}-ndim"), bias.dim_params.len() as i32);
    for (d, dim) in bias.dim_params.iter().enumerate() {
        print_awh_dim(out, n, &format!("{prefix}-dim{}", d + 1), dim);
    }
}

fn print_awh(out: &mut String, n: usize, awh: &AwhParams) {
    p_str(
        out,
        n,
        "awh-potential",
        lookup(en::AWH_POTENTIAL_TYPE, awh.potential as i64),
    );
    p_int(out, n, "awh-seed", awh.seed as i32);
    p_int(out, n, "awh-nstout", awh.nstout);
    p_int(out, n, "awh-nstsample", awh.nst_sample_coord);
    p_int(out, n, "awh-nsamples-update", awh.num_samples_update_free_energy);
    p_bool(out, n, "awh-share-bias-multisim", awh.share_bias_multisim);
    p_int(out, n, "awh-nbias", awh.bias.len() as i32);
    for (k, bias) in awh.bias.iter().enumerate() {
        print_awh_bias(out, n, bias, &format!("awh{}", k + 1));
    }
}

fn print_rot(out: &mut String, n: usize, rot: &Rot) {
    p_int(out, n, "rot-nstrout", rot.nstrout);
    p_int(out, n, "rot-nstsout", rot.nstsout);
    p_int(out, n, "rot-ngroups", rot.grp.len() as i32);
    for (g, grp) in rot.grp.iter().enumerate() {
        indent(out, n);
        let _ = writeln!(out, "rot-group {g}:");
        let n = n + 2;
        p_str(
            out,
            n,
            "rot-type",
            lookup(en::ENFORCED_ROTATION_GROUP_TYPE, grp.e_type as i64),
        );
        p_bool(out, n, "rot-massw", grp.b_mass_w != 0);
        p_ivec_block(out, n, "atom", &grp.ind);
        p_rvecs(out, n, "x-ref", &grp.x_ref_original);
        p_rvec(out, n, "rot-vec", &grp.input_vec);
        p_rvec(out, n, "rot-pivot", &grp.pivot);
        p_real(out, n, "rot-rate", grp.rate);
        p_real(out, n, "rot-k", grp.k);
        p_real(out, n, "rot-slab-dist", grp.slab_dist);
        p_real(out, n, "rot-min-gauss", grp.min_gaussian);
        p_real(out, n, "rot-eps", grp.eps);
        p_str(
            out,
            n,
            "rot-fit-method",
            lookup(en::ROTATION_GROUP_FITTING, grp.e_fittype as i64),
        );
        p_int(out, n, "rot-potfit-nstep", grp.pot_angle_nstep);
        p_real(out, n, "rot-potfit-step", grp.pot_angle_step);
    }
}

fn print_imd(out: &mut String, n: usize, ind: &[i32]) {
    p_int(out, n, "IMD-atoms", ind.len() as i32);
    p_ivec_block(out, n, "atom", ind);
}

fn print_swap(out: &mut String, n: usize, swap: &SwapCoords) {
    p_int(out, n, "swap-frequency", swap.nstswap);

    let split0 = swap.groups.first();
    p_bool(out, n, "massw_split0", swap.massw_split[0]);
    if let Some(g) = split0 {
        p_ivec_block(out, n, "split atoms group 0", &g.ind);
    }
    p_bool(out, n, "massw_split1", swap.massw_split[1]);
    if let Some(g) = swap.groups.get(1) {
        p_ivec_block(out, n, "split atoms group 1", &g.ind);
    }
    if let Some(g) = swap.groups.get(2) {
        let key = format!("solvent group {}", g.molname);
        p_ivec_block(out, n, &key, &g.ind);
    }
    for g in swap.groups.iter().skip(3) {
        let key = format!("ion group {}", g.molname);
        p_ivec_block(out, n, &key, &g.ind);
    }

    p_real(out, n, "cyl0-r", swap.cyl0r);
    p_real(out, n, "cyl0-up", swap.cyl0u);
    p_real(out, n, "cyl0-down", swap.cyl0l);
    p_real(out, n, "cyl1-r", swap.cyl1r);
    p_real(out, n, "cyl1-up", swap.cyl1u);
    p_real(out, n, "cyl1-down", swap.cyl1l);
    p_int(out, n, "coupl-steps", swap.n_average);
    for comp in 0..2 {
        let label = if comp == 0 { 'A' } else { 'B' };
        for g in swap.groups.iter().skip(3) {
            let key = format!("{}-in-{}", g.molname, label);
            p_int(out, n, &key, g.nmol_req[comp]);
        }
    }
    p_real(out, n, "threshold", swap.threshold);
    p_real(out, n, "bulk-offsetA", swap.bulk_offset[0]);
    p_real(out, n, "bulk-offsetB", swap.bulk_offset[1]);
}

// ---------------------------------------------------------------------------
// The module parameter key value tree
// ---------------------------------------------------------------------------

/// Defaults of the `applied-forces` and `fast-multipole-method` module options,
/// in the order in which `MDModules::makeModuleOptions()` registers them.
/// Each entry is the key path relative to the root of `ir->params` and the text
/// `simpleValueToString()` writes for the option's default value.
const MODULE_PARAMS: &[(&[&str], &str)] = &[
    (&["applied-forces", "electric-field", "x", "E0"], "0"),
    (&["applied-forces", "electric-field", "x", "omega"], "0"),
    (&["applied-forces", "electric-field", "x", "t0"], "0"),
    (&["applied-forces", "electric-field", "x", "sigma"], "0"),
    (&["applied-forces", "electric-field", "y", "E0"], "0"),
    (&["applied-forces", "electric-field", "y", "omega"], "0"),
    (&["applied-forces", "electric-field", "y", "t0"], "0"),
    (&["applied-forces", "electric-field", "y", "sigma"], "0"),
    (&["applied-forces", "electric-field", "z", "E0"], "0"),
    (&["applied-forces", "electric-field", "z", "omega"], "0"),
    (&["applied-forces", "electric-field", "z", "t0"], "0"),
    (&["applied-forces", "electric-field", "z", "sigma"], "0"),
    (&["applied-forces", "density-guided-simulation", "active"], "false"),
    (&["applied-forces", "density-guided-simulation", "group"], "protein"),
    (
        &["applied-forces", "density-guided-simulation", "similarity-measure"],
        "inner-product",
    ),
    (
        &["applied-forces", "density-guided-simulation", "atom-spreading-weight"],
        "unity",
    ),
    (
        &["applied-forces", "density-guided-simulation", "force-constant"],
        "1e+09",
    ),
    (
        &[
            "applied-forces",
            "density-guided-simulation",
            "gaussian-transform-spreading-width",
        ],
        "0.2",
    ),
    (
        &[
            "applied-forces",
            "density-guided-simulation",
            "gaussian-transform-spreading-range-in-multiples-of-width",
        ],
        "4",
    ),
    (
        &[
            "applied-forces",
            "density-guided-simulation",
            "reference-density-filename",
        ],
        "reference.mrc",
    ),
    (&["applied-forces", "density-guided-simulation", "nst"], "1"),
    (
        &["applied-forces", "density-guided-simulation", "normalize-densities"],
        "true",
    ),
    (
        &["applied-forces", "density-guided-simulation", "adaptive-force-scaling"],
        "false",
    ),
    (
        &[
            "applied-forces",
            "density-guided-simulation",
            "adaptive-force-scaling-time-constant",
        ],
        "4",
    ),
    (
        &["applied-forces", "density-guided-simulation", "shift-vector"],
        "",
    ),
    (
        &["applied-forces", "density-guided-simulation", "transformation-matrix"],
        "",
    ),
    (&["applied-forces", "qmmm-cp2k", "active"], "false"),
    (&["applied-forces", "qmmm-cp2k", "qmgroup"], "System"),
    (&["applied-forces", "qmmm-cp2k", "qmmethod"], "PBE"),
    (&["applied-forces", "qmmm-cp2k", "qmfilenames"], ""),
    (&["applied-forces", "qmmm-cp2k", "qmcharge"], "0"),
    (&["applied-forces", "qmmm-cp2k", "qmmultiplicity"], "1"),
    (&["applied-forces", "colvars", "active"], "false"),
    (&["applied-forces", "colvars", "configfile"], ""),
    (&["applied-forces", "colvars", "seed"], "-1"),
    (&["applied-forces", "nnpot", "active"], "false"),
    (&["applied-forces", "nnpot", "modelfile"], "model.pt"),
    (&["applied-forces", "nnpot", "input-group"], "System"),
    (&["applied-forces", "nnpot", "link-type"], "H"),
    (&["applied-forces", "nnpot", "link-distance"], "0.1"),
    (&["applied-forces", "nnpot", "pair-cutoff"], "0"),
    (&["applied-forces", "nnpot", "embedding"], "mechanical"),
    (&["applied-forces", "nnpot", "nnp-charge"], "0"),
    (&["applied-forces", "nnpot", "model-input1"], ""),
    (&["applied-forces", "nnpot", "model-input2"], ""),
    (&["applied-forces", "nnpot", "model-input3"], ""),
    (&["applied-forces", "nnpot", "model-input4"], ""),
    (&["applied-forces", "nnpot", "model-input5"], ""),
    (&["applied-forces", "nnpot", "model-input6"], ""),
    (&["applied-forces", "nnpot", "model-input7"], ""),
    (&["applied-forces", "nnpot", "model-input8"], ""),
    (&["applied-forces", "nnpot", "model-input9"], ""),
    (&["fast-multipole-method", "fmm", "backend"], "inactive"),
    (&["fast-multipole-method", "fmm", "exafmm-order"], "6"),
    (&["fast-multipole-method", "fmm", "exafmm-direct-range"], "2"),
    (&["fast-multipole-method", "fmm", "exafmm-direct-provider"], "GROMACS"),
    (&["fast-multipole-method", "fmm", "exafmm-tree-type"], "uniform"),
    (&["fast-multipole-method", "fmm", "exafmm-tree-depth"], "0"),
    (
        &["fast-multipole-method", "fmm", "exafmm-max-particles-per-cell"],
        "0",
    ),
    (&["fast-multipole-method", "fmm", "fmsolvr-order"], "8"),
    (&["fast-multipole-method", "fmm", "fmsolvr-direct-range"], "1"),
    (&["fast-multipole-method", "fmm", "fmsolvr-direct-provider"], "FMM"),
    (
        &["fast-multipole-method", "fmm", "fmsolvr-dipole-compensation"],
        "true",
    ),
    (&["fast-multipole-method", "fmm", "fmsolvr-tree-depth"], "3"),
    (&["fast-multipole-method", "fmm", "fmsolvr-sparse"], "false"),
];

/// Writes the `ir->params` key value tree the way
/// `MDModules().adjustInputrecBasedOnModules()` + `dumpKeyValueTree()` do:
/// the tree is rebuilt from the module options (which fixes the order), pulling
/// each value from the file when it is present there.
fn dump_params(out: &mut String, n: usize, params: &KvtValue, adjust: bool) {
    if !adjust {
        dump_kvt_object(out, n, params);
        return;
    }
    // Build the nested structure of the module option tree.
    let mut root: Vec<(String, Node)> = Vec::new();
    for (path, default) in MODULE_PARAMS {
        let mut level = &mut root;
        for (i, key) in path.iter().enumerate() {
            let last = i + 1 == path.len();
            let pos = match level.iter().position(|(k, _)| k == key) {
                Some(pos) => pos,
                None => {
                    level.push((
                        (*key).to_string(),
                        if last {
                            Node::Leaf(default.to_string())
                        } else {
                            Node::Object(Vec::new())
                        },
                    ));
                    level.len() - 1
                }
            };
            match &mut level[pos].1 {
                Node::Object(children) => level = children,
                Node::Leaf(_) => break,
            }
        }
    }
    // Overwrite leaves with the values stored in the file.
    for (path, _) in MODULE_PARAMS {
        if let Some(value) = lookup_kvt(params, path) {
            let text = value.to_text();
            if let Some(Node::Leaf(slot)) = node_at_mut(&mut root, path) {
                *slot = text;
            }
        }
    }
    dump_node(out, n, &root);
}

/// `dumpKeyValueTree()` on a tree read straight from the file.
fn dump_kvt_object(out: &mut String, n: usize, value: &KvtValue) {
    let props = match value.as_object() {
        Some(props) => props,
        None => return,
    };
    for (key, value) in props {
        match value {
            KvtValue::Object(_) => {
                indent(out, n);
                let _ = writeln!(out, "{key}:");
                dump_kvt_object(out, n + 2, value);
            }
            KvtValue::Array(elements) if elements.iter().all(|e| e.as_object().is_some()) => {
                indent(out, n);
                let _ = writeln!(out, "{key}:");
                for element in elements {
                    dump_kvt_object(out, n + 2, element);
                }
            }
            KvtValue::Array(elements) => {
                indent(out, n);
                let width = 33 - n;
                let _ = write!(out, "{key:<width$} = [");
                for element in elements {
                    let _ = write!(out, " {}", element.to_text());
                }
                out.push_str(" ]\n");
            }
            other => {
                indent(out, n);
                let width = 33 - n;
                let _ = write!(out, "{key:<width$} = {}\n", other.to_text());
            }
        }
    }
}

enum Node {
    Object(Vec<(String, Node)>),
    Leaf(String),
}

fn node_at_mut<'a>(nodes: &'a mut Vec<(String, Node)>, path: &[&str]) -> Option<&'a mut Node> {
    let (first, rest) = path.split_first()?;
    let pos = nodes.iter().position(|(k, _)| k == first)?;
    if rest.is_empty() {
        return Some(&mut nodes[pos].1);
    }
    match &mut nodes[pos].1 {
        Node::Object(children) => node_at_mut(children, rest),
        Node::Leaf(_) => None,
    }
}

fn lookup_kvt<'a>(value: &'a KvtValue, path: &[&str]) -> Option<&'a KvtValue> {
    let (first, rest) = path.split_first()?;
    let next = value.get(first)?;
    if rest.is_empty() {
        Some(next)
    } else {
        lookup_kvt(next, rest)
    }
}

fn dump_node(out: &mut String, n: usize, nodes: &[(String, Node)]) {
    for (key, node) in nodes {
        match node {
            Node::Object(children) => {
                indent(out, n);
                let _ = writeln!(out, "{key}:");
                dump_node(out, n + 2, children);
            }
            Node::Leaf(text) => {
                indent(out, n);
                let width = 33 - n;
                let _ = write!(out, "{key:<width$} = {text}\n");
            }
        }
    }
}
