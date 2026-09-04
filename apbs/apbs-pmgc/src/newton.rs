// APBS PMGC newton - Newton driver for nonlinear problems
// Port of pmgc/newtond.c

use rayon::prelude::*;
use std::sync::OnceLock;

const PAR_THRESHOLD: usize = 16_384;
const DEFAULT_NEWTON_ITMAX_CAP: i32 = 20;

fn configured_itmax(requested: i32) -> i32 {
    static VALUE: OnceLock<i32> = OnceLock::new();
    *VALUE.get_or_init(|| {
        if let Ok(value) = std::env::var("APBS_RUST_NEWTON_ITMAX") {
            return value.parse::<i32>().unwrap_or(requested).max(1);
        }
        if std::env::var("APBS_RUST_OVERRIDE_ITMAX").is_ok() {
            return requested.max(1);
        }
        requested.clamp(1, DEFAULT_NEWTON_ITMAX_CAP)
    })
}

/// Newton iteration driver for nonlinear PBE
pub fn newton(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc_all: &[f64], fc_all: &[f64],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    itmax: i32, errtol: f64,
    nlev: i32,
    pc: &[f64], iz: &[i32],
    nu1: i32, nu2: i32,
    omegal: f64, irite: i32,
    mgsolv: i32,
) {
    let n = nx * ny * nz;
    let itmax = configured_itmax(itmax);
    let ac_fine = &ac[0..4 * n];
    let cc_fine = &cc_all[0..n];
    let fc_fine = &fc_all[0..n];
    let mut j_cc = vec![0.0; n];
    let mut du = vec![0.0; n];
    // Multigrid V-cycle solves use the (fine-level) Jacobian as the diagonal
    // coefficient; all coarser levels keep the pre-built linear operator.
    let mut jac_cc_all = cc_all.to_vec();
    let mut rhs = vec![0.0; n];
    let fnorm = crate::blas::xnrm2(n, fc_fine, 0);

    for _iter in 0..itmax {
        // Compute nonlinear residual: r = f - N(u)
        compute_residual(nx, ny, nz, ac_fine, cc_fine, fc_fine, u, r);

        let rnorm = crate::blas::xnrm2(n, r, 0);
        // Use the same relative-residual test as the rest of the PMG driver
        // (and as the C Vnewton driver): stopping on the absolute rnorm kept
        // the loop at its 20-iteration cap even after convergence in
        // rnorm/||f|| was reached.
        let converged = if fnorm > 0.0 {
            rnorm < errtol * fnorm
        } else {
            rnorm < errtol
        };
        if converged {
            return;
        }

        // Compute Jacobian: J = dN/du = A + cc * cosh(u) on the diagonal.
        if n >= PAR_THRESHOLD {
            j_cc.par_iter_mut().enumerate().for_each(|(i, out)| {
                let u_val = u[i].clamp(crate::mgfas::SINH_MIN, crate::mgfas::SINH_MAX);
                *out = cc_fine[i] * u_val.cosh();
            });
        } else {
            for i in 0..n {
                let u_val = u[i].clamp(crate::mgfas::SINH_MIN, crate::mgfas::SINH_MAX);
                // smooth()/mresid use (o_c + cc_param) as diagonal,
                // so pass only the reaction part here.
                j_cc[i] = cc_fine[i] * u_val.cosh();
            }
        }
        jac_cc_all[..n].copy_from_slice(&j_cc);

        // Solve the linearized system J * du = r with one full multigrid V-cycle.
        if n >= PAR_THRESHOLD {
            du.par_iter_mut().for_each(|v| *v = 0.0);
        } else {
            du.fill(0.0);
        }
        rhs.copy_from_slice(r);
        crate::mgcs::mgcs(
            nlev,
            nx as i32, ny as i32, nz as i32,
            ipc, rpc,
            ac, &jac_cc_all, &rhs,
            0, 0, 0,
            pc, iz,
            &mut du, w1, w2, r,
            nu1, nu2, omegal, irite, mgsolv,
        );

        // Update: u = u + du
        if n >= PAR_THRESHOLD {
            u.par_iter_mut().zip(&du).for_each(|(ui, dui)| *ui += *dui);
        } else {
            for i in 0..n {
                u[i] += du[i];
            }
        }
    }
}

fn compute_residual(
    nx: usize, ny: usize, nz: usize,
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &[f64], r: &mut [f64],
) {
    let n = nx * ny * nz;
    let o_c = &ac[0..n];
    let o_e = &ac[n..2 * n];
    let o_n = &ac[2 * n..3 * n];
    let u_c = &ac[3 * n..4 * n];

    // r = f - Au
    crate::blas::mresid(nx, ny, nz, &[], &[], o_c, o_e, o_n, u_c, cc, u, fc, r);

    // mresid already contributes linear +cc*u on diagonal.
    // Add only nonlinear increment cc*(sinh(u)-u).
    if n >= PAR_THRESHOLD {
        r.par_iter_mut().enumerate().for_each(|(i, ri)| {
            let u_val = u[i].clamp(crate::mgfas::SINH_MIN, crate::mgfas::SINH_MAX);
            *ri -= cc[i] * (u_val.sinh() - u_val);
        });
    } else {
        for i in 0..n {
            let u_val = u[i].clamp(crate::mgfas::SINH_MIN, crate::mgfas::SINH_MAX);
            r[i] -= cc[i] * (u_val.sinh() - u_val);
        }
    }
}
