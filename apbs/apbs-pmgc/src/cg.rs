// APBS PMGC CG - Conjugate Gradient solver
// Port of pmgc/cgd.c

use crate::blas;
use std::sync::OnceLock;

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

/// Conjugate Gradient solver for Ax = b
/// Returns the number of iterations and whether convergence was achieved
pub fn cg(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    x: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, errtol: f64,
    rinf_norm: f64,
) {
    let n = nx * ny * nz;

    // Split ac into components: [o_c, o_e, o_n, u_c]
    let o_e = &ac[n..2 * n];
    let o_n = &ac[2 * n..3 * n];
    let u_c = &ac[3 * n..4 * n];

    // Debug: check fc and operator for NaN before residual
    let fc_nan = fc.iter().any(|x| x.is_nan());
    let fc_max = fc.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let o_c_nan = ac[0..n].iter().any(|x| x.is_nan());
    let cc_nan = cc.iter().any(|x| x.is_nan());
    if debug_enabled() {
        eprintln!("[DEBUG-CG] pre-resid: n={}, fc_nan={}, fc_max={:.4e}, o_c_nan={}, cc_nan={}", n, fc_nan, fc_max, o_c_nan, cc_nan);
    }

    // r = f - Ax
    crate::blas::mresid(nx, ny, nz, ipc, rpc, &ac[0..n], o_e, o_n, u_c, cc, x, fc, r);

    let r_nan = r.iter().any(|x| x.is_nan());
    let r_max = r.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    if debug_enabled() {
        eprintln!("[DEBUG-CG] post-resid: r_nan={}, r_max={:.4e}", r_nan, r_max);
    }

    // p = r (initial search direction)
    blas::xcopy(n, r, 0, w1, 0);

    // rr = r . r
    let mut rr = blas::xdot(n, r, 0, r, 0);
    let _rr0 = rr;

    if debug_enabled() {
        eprintln!("[DEBUG-CG] init: rr={:.4e}, rnorm={:.4e}, errtol={:.4e}, rinf_norm={:.4e}", rr, rr.sqrt(), errtol, rinf_norm);
    }

    if rr.sqrt() < errtol * rinf_norm {
        *iters = 0;
        return;
    }

    for iter in 0..itmax {
        *iters = iter + 1;

        // w2 = A * p
        crate::matvec::matvec(nx, ny, nz, ipc, rpc, ac, cc, fc, w1, w2);

        // alpha = rr / (p . Ap)
        let pap = blas::xdot(n, w1, 0, w2, 0);
        if debug_enabled() {
            eprintln!("[DEBUG-CG] iter={}: pap={:.4e}, rr={:.4e}", iter, pap, rr);
        }
        if pap.abs() < 1.0e-30 {
            if debug_enabled() {
                eprintln!("[DEBUG-CG] breaking: pap too small");
            }
            break;
        }
        let alpha = rr / pap;

        // x = x + alpha * p
        blas::xaxpy(n, alpha, w1, 0, x, 0);

        // r = r - alpha * Ap
        blas::xaxpy(n, -alpha, w2, 0, r, 0);

        let rr_new = blas::xdot(n, r, 0, r, 0);

        // Check convergence
        if rr_new.sqrt() < errtol * rinf_norm {
            break;
        }

        let beta = rr_new / rr;
        rr = rr_new;

        // p = r + beta * p
        for i in 0..n {
            w1[i] = r[i] + beta * w1[i];
        }
    }
}
