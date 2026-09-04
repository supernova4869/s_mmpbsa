// APBS PMGC mgcs - Linear V-cycle solver
// Port of pmgc/mgcsd.c

use std::sync::OnceLock;

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

fn force_alpha_one() -> bool {
    static FORCE: OnceLock<bool> = OnceLock::new();
    *FORCE.get_or_init(|| std::env::var("APBS_RUST_FORCE_ALPHA1").ok().as_deref() == Some("1"))
}

/// Linear V-cycle multigrid solver
pub fn mgcs(
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    ac_off: usize, cc_off: usize, fc_off: usize,
    pc: &[f64], iz: &[i32],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu1: i32, nu2: i32,
    omegal: f64,
    irite: i32,
    mgsolv: i32,
) {
    mgcs_level(
        0,
        nlev,
        nx,
        ny,
        nz,
        ipc,
        rpc,
        ac,
        cc,
        fc,
        ac_off,
        cc_off,
        fc_off,
        0,
        pc,
        iz,
        u,
        w1,
        w2,
        r,
        nu1,
        nu2,
        omegal,
        irite,
        mgsolv,
    );
}

fn mgcs_level(
    lev_idx: usize,
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    ac_off: usize, cc_off: usize, fc_off: usize,
    pc_off: usize,
    pc: &[f64], iz: &[i32],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu1: i32, nu2: i32,
    omegal: f64,
    irite: i32,
    mgsolv: i32,
) {
    let current_ipc = &ipc[lev_idx * 20..];
    let current_rpc = &rpc[lev_idx * 20..];
    let nf = (nx * ny * nz) as usize;
    let numdia = current_ipc[0];

    if nlev == 1 {
        solve_coarsest(
            nx as usize, ny as usize, nz as usize,
            current_ipc, current_rpc,
            &ac[ac_off..ac_off + 4 * nf],
            &cc[cc_off..cc_off + nf],
            &fc[fc_off..fc_off + nf],
            u,
            mgsolv,
    );
    return;
    }

    // Pre-smoothing (nu1 iterations, forward red-black order)
    smooth_on_level(
        nx as usize, ny as usize, nz as usize,
        current_ipc, current_rpc,
        &ac[ac_off..ac_off + 4 * nf],
        &cc[cc_off..cc_off + nf],
        &fc[fc_off..fc_off + nf],
        u, w1, w2, r,
        nu1, omegal, numdia, 0,
    );

    // Compute residual: r = f - Au
    crate::blas::mresid(
        nx as usize, ny as usize, nz as usize,
        current_ipc, current_rpc,
        &ac[ac_off..ac_off + nf],
        &ac[ac_off + nf..ac_off + 2 * nf],
        &ac[ac_off + 2 * nf..ac_off + 3 * nf],
        &ac[ac_off + 3 * nf..ac_off + 4 * nf],
        &cc[cc_off..cc_off + nf],
        u,
        &fc[fc_off..fc_off + nf],
        r,
    );

    // Make coarse grid
    let (nx_c, ny_c, nz_c) = crate::build_str::make_coarse(nx, ny, nz);
    let nc = (nx_c * ny_c * nz_c) as usize;
    let pc_len = 27 * nc;
    let pc_slice = if pc_off + pc_len <= pc.len() {
        &pc[pc_off..pc_off + pc_len]
    } else {
        &[]
    };

    // Restrict RESIDUAL to coarse grid as RHS
    let mut r_c = vec![0.0; nc];
    restrict_vec(
        nx as usize, ny as usize, nz as usize,
        nx_c as usize, ny_c as usize, nz_c as usize,
        r, &mut r_c, pc_slice,
    );

    let ac_off_c = ac_off + 4 * nf;
    let cc_off_c = cc_off + nf;
    let ac_coarse = &ac[ac_off_c..ac_off_c + 4 * nc];

    // Initialize coarse grid solution to zero
    let mut u_c_vec = vec![0.0; nc];

    // Recurse on coarse level
    mgcs_level(
        lev_idx + 1,
        nlev - 1,
        nx_c, ny_c, nz_c,
        ipc, rpc,
        ac, cc, &r_c,
        ac_off_c, cc_off_c, 0,
        pc_off + pc_len,
        pc, iz,
        &mut u_c_vec,
        w1, w2, r,
        nu1, nu2,
        omegal,
        irite,
        mgsolv,
    );

    // Prolongate correction to fine grid and add directly.
    let mut corr = vec![0.0; nf];
    prolongate_vec(
        nx as usize, ny as usize, nz as usize,
        nx_c as usize, ny_c as usize, nz_c as usize,
        &u_c_vec, &mut corr, pc_slice,
    );
    let alpha = if force_alpha_one() {
        1.0
    } else {
        correction_damping(
        nx_c as usize,
        ny_c as usize,
        nz_c as usize,
        &ipc[(lev_idx + 1) * 20..],
        &rpc[(lev_idx + 1) * 20..],
        ac_coarse,
        &cc[cc_off_c..cc_off_c + nc],
        &u_c_vec,
        &r_c,
    )
    };
    if debug_enabled() {
        let corr_norm = crate::blas::xnrm2(nf, &corr, 0);
        let coarse_norm = crate::blas::xnrm2(nc, &u_c_vec, 0);
        let rhs_norm = crate::blas::xnrm2(nc, &r_c, 0);
        eprintln!(
            "[DEBUG-MGCS] lev={} -> {}, alpha={:.4e}, coarse_norm={:.4e}, corr_norm={:.4e}, rhs_norm={:.4e}",
            lev_idx,
            lev_idx + 1,
            alpha,
            coarse_norm,
            corr_norm,
            rhs_norm
        );
    }
    for i in 0..nf {
        u[i] += alpha * corr[i];
    }

    // Post-smoothing (nu2 iterations, adjoint/reversed red-black order)
    smooth_on_level(
        nx as usize, ny as usize, nz as usize,
        current_ipc, current_rpc,
        &ac[ac_off..ac_off + 4 * nf],
        &cc[cc_off..cc_off + nf],
        &fc[fc_off..fc_off + nf],
        u, w1, w2, r,
        nu2, omegal, numdia, 1,
    );
}

/// Smooth on current level
fn smooth_on_level(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu: i32, omega: f64, numdia: i32, iadjoint: i32,
) {
    let n = nx * ny * nz;
    let o_c = &ac[0..n];
    let o_e = &ac[n..2*n];
    let o_n = &ac[2*n..3*n];
    let u_c = &ac[3*n..4*n];
    let mut iters = 0;
    crate::gs::gsrb(
        nx, ny, nz, ipc, rpc, o_c, cc, fc, o_e, o_n, u_c,
        u, w1, w2, r,
        nu, &mut iters, 0.0, omega, 0, iadjoint, numdia,
    );
}

/// Solve on coarsest grid directly using CG
fn solve_coarsest(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
    mgsolv: i32,
) {
    let nf = nx * ny * nz;
    u.fill(0.0);
    if mgsolv == 1 {
        solve_coarsest_direct(nx, ny, nz, ac, cc, fc, u);
    } else {
        let errtol = 1.0e-6;
        let mut iters = 0;
        let rinf = crate::blas::xnrm2(nf, fc, 0);
        let mut w1 = vec![0.0; nf];
        let mut w2 = vec![0.0; nf];
        let mut r = vec![0.0; nf];

        crate::cg::cg(
            nx, ny, nz, ipc, rpc, ac, cc, fc,
            u, &mut w1, &mut w2, &mut r,
            100, &mut iters, errtol, rinf,
        );
    }
}

fn solve_coarsest_direct(
    nx: usize, ny: usize, nz: usize,
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
) {
    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }

    let nf = nx * ny * nz;
    let nxny = nx * ny;
    let o_c = &ac[0..nf];
    let o_e = &ac[nf..2 * nf];
    let o_n = &ac[2 * nf..3 * nf];
    let u_c = &ac[3 * nf..4 * nf];

    let nix = nx - 2;
    let niy = ny - 2;
    let niz = nz - 2;
    let nint = nix * niy * niz;
    if nint == 0 {
        return;
    }

    let row_of = |i: usize, j: usize, k: usize| -> usize {
        (i - 1) + (j - 1) * nix + (k - 1) * nix * niy
    };

    let mut a = vec![0.0f64; nint * nint];
    let mut b = vec![0.0f64; nint];
    for k in 1..(nz - 1) {
        for j in 1..(ny - 1) {
            for i in 1..(nx - 1) {
                let ip = i + j * nx + k * nxny;
                let row = row_of(i, j, k);
                a[row * nint + row] = o_c[ip] + cc[ip];
                b[row] = fc[ip];

                if i > 1 {
                    let col = row_of(i - 1, j, k);
                    a[row * nint + col] = -o_e[ip - 1];
                }
                if i + 1 < nx - 1 {
                    let col = row_of(i + 1, j, k);
                    a[row * nint + col] = -o_e[ip];
                }
                if j > 1 {
                    let col = row_of(i, j - 1, k);
                    a[row * nint + col] = -o_n[ip - nx];
                }
                if j + 1 < ny - 1 {
                    let col = row_of(i, j + 1, k);
                    a[row * nint + col] = -o_n[ip];
                }
                if k > 1 {
                    let col = row_of(i, j, k - 1);
                    a[row * nint + col] = -u_c[ip - nxny];
                }
                if k + 1 < nz - 1 {
                    let col = row_of(i, j, k + 1);
                    a[row * nint + col] = -u_c[ip];
                }
            }
        }
    }

    for piv in 0..nint {
        let mut pivot_row = piv;
        let mut pivot_abs = a[piv * nint + piv].abs();
        for r in (piv + 1)..nint {
            let cand = a[r * nint + piv].abs();
            if cand > pivot_abs {
                pivot_abs = cand;
                pivot_row = r;
            }
        }
        if pivot_abs <= 1.0e-30 {
            return;
        }
        if pivot_row != piv {
            for c in piv..nint {
                a.swap(piv * nint + c, pivot_row * nint + c);
            }
            b.swap(piv, pivot_row);
        }
        let diag = a[piv * nint + piv];
        for r in (piv + 1)..nint {
            let factor = a[r * nint + piv] / diag;
            if factor == 0.0 {
                continue;
            }
            a[r * nint + piv] = 0.0;
            for c in (piv + 1)..nint {
                a[r * nint + c] -= factor * a[piv * nint + c];
            }
            b[r] -= factor * b[piv];
        }
    }

    let mut x = vec![0.0f64; nint];
    for r in (0..nint).rev() {
        let mut sum = b[r];
        for c in (r + 1)..nint {
            sum -= a[r * nint + c] * x[c];
        }
        let diag = a[r * nint + r];
        if diag.abs() <= 1.0e-30 {
            return;
        }
        x[r] = sum / diag;
    }

    for k in 1..(nz - 1) {
        for j in 1..(ny - 1) {
            for i in 1..(nx - 1) {
                let ip = i + j * nx + k * nxny;
                u[ip] = x[row_of(i, j, k)];
            }
        }
    }
}

/// Restriction: fine grid to coarse grid (simple averaging with clamping)
pub(crate) fn restrict_vec(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    fine: &[f64], coarse: &mut [f64], pc: &[f64],
) {
    crate::matvec::restrc(nxf, nyf, nzf, nxc, nyc, nzc, fine, coarse, pc);
}

/// Bilinear prolongation: coarse grid to fine grid
/// For 2:1 coarsening, each fine point gets interpolated from coarse neighbors.
pub(crate) fn prolongate_vec(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    coarse: &[f64], fine: &mut [f64], pc: &[f64],
) {
    crate::matvec::interp_pmg(nxc, nyc, nzc, nxf, nyf, nzf, coarse, fine, pc);
}

/// Extraction operator corresponding to APBS C `Vextrac`.
#[allow(dead_code)]
pub(crate) fn extract_vec(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    fine: &[f64], coarse: &mut [f64],
) {
    for kc in 1..nzc.saturating_sub(1) {
        let kf = 2 * kc;
        for jc in 1..nyc.saturating_sub(1) {
            let jf = 2 * jc;
            for ic in 1..nxc.saturating_sub(1) {
                let if_ = 2 * ic;
                if if_ < nxf && jf < nyf && kf < nzf {
                    let ipc = ic + jc * nxc + kc * nxc * nyc;
                    let ipf = if_ + jf * nxf + kf * nxf * nyf;
                    coarse[ipc] = fine[ipf];
                }
            }
        }
    }
}

fn correction_damping(
    nx: usize,
    ny: usize,
    nz: usize,
    ipc: &[i32],
    rpc: &[f64],
    ac: &[f64],
    cc: &[f64],
    corr: &[f64],
    rhs: &[f64],
) -> f64 {
    let nf = nx * ny * nz;
    let mut ac_corr = vec![0.0; nf];
    crate::matvec::matvec(nx, ny, nz, ipc, rpc, ac, cc, &vec![0.0; nf], corr, &mut ac_corr);
    let num = crate::blas::xdot(nf, corr, 0, rhs, 0);
    let den = crate::blas::xdot(nf, corr, 0, &ac_corr, 0);
    if den.is_finite() && den > 0.0 && num.is_finite() {
        num / den
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::{correction_damping, extract_vec, prolongate_vec, restrict_vec};

    #[test]
    fn restrict_vec_preserves_constant_field() {
        let fine = vec![3.25; 5 * 5 * 5];
        let mut coarse = vec![0.0; 3 * 3 * 3];
        restrict_vec(5, 5, 5, 3, 3, 3, &fine, &mut coarse, &[]);
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    let v = coarse[i + j * 3 + k * 9];
                    if i == 1 && j == 1 && k == 1 {
                        assert!((v - 3.25).abs() < 1.0e-12);
                    } else {
                        assert!(v.abs() < 1.0e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn prolongate_vec_preserves_constant_field() {
        let coarse = vec![1.75; 3 * 3 * 3];
        let mut fine = vec![0.0; 5 * 5 * 5];
        prolongate_vec(5, 5, 5, 3, 3, 3, &coarse, &mut fine, &[]);
        for k in 0..5 {
            for j in 0..5 {
                for i in 0..5 {
                    let v = fine[i + j * 5 + k * 25];
                    if i == 0 || i == 4 || j == 0 || j == 4 || k == 0 || k == 4 {
                        assert!(v.abs() < 1.0e-12);
                    } else {
                        assert!((v - 1.75).abs() < 1.0e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn prolongate_vec_copies_coincident_coarse_nodes() {
        let nxc = 3usize;
        let nyc = 3usize;
        let nzc = 3usize;
        let mut coarse = vec![0.0; nxc * nyc * nzc];
        let idx = 1 + 1 * nxc + 1 * nxc * nyc;
        coarse[idx] = 7.0;

        let mut fine = vec![0.0; 5 * 5 * 5];
        prolongate_vec(5, 5, 5, nxc, nyc, nzc, &coarse, &mut fine, &[]);

        let coincident = 2 + 2 * 5 + 2 * 5 * 5;
        assert!((fine[coincident] - 7.0).abs() < 1.0e-12);
    }

    #[test]
    fn extract_vec_copies_coincident_fine_nodes() {
        let nxf = 5usize;
        let nyf = 5usize;
        let nzf = 5usize;
        let mut fine = vec![0.0; nxf * nyf * nzf];
        let ipf = 2 + 2 * nxf + 2 * nxf * nyf;
        fine[ipf] = 9.5;
        let mut coarse = vec![0.0; 3 * 3 * 3];
        extract_vec(nxf, nyf, nzf, 3, 3, 3, &fine, &mut coarse);
        let ipc = 1 + 1 * 3 + 1 * 3 * 3;
        assert!((coarse[ipc] - 9.5).abs() < 1.0e-12);
    }

    #[test]
    fn correction_damping_falls_back_to_one_for_nonpositive_denominator() {
        let nx = 1usize;
        let ny = 1usize;
        let nz = 1usize;
        let ipc = [4i32];
        let rpc = [1.0f64, 1.0, 1.0];
        let ac = vec![0.0; 4];
        let cc = vec![0.0];
        let corr = vec![2.0];
        let residual = vec![3.0];
        let alpha = correction_damping(nx, ny, nz, &ipc, &rpc, &ac, &cc, &corr, &residual);
        assert!((alpha - 1.0).abs() < 1.0e-12);
    }
}
