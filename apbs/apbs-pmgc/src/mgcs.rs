// APBS PMGC mgcs - Linear multigrid V-cycle solver
// Port of pmgc/mgcsd.c (Vmvcs): outer iteration loop with the istop=1
// relative L1 residual stopping test, nu1/nu2 red-black Gauss-Seidel
// smoothing, pc-based restriction/interpolation, the Hackbusch/Reusken
// correction damping on every up-cycle level, and a direct (dpbfa/dpbsl)
// or CG solve on the coarsest level.

use crate::build_b::Banded;
use crate::matvec::{split_bands14, split_bands4, Bands14};

/// One multigrid level: operator, Helmholtz term, right-hand side and the
/// prolongation operator toward the next coarser level.
pub struct Level {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    /// 4*nf bands for the fine level (7-point), 14*nf for Galerkin levels.
    pub ac: Vec<f64>,
    pub cc: Vec<f64>,
    pub fc: Vec<f64>,
    pub numdia: i32,
    /// 27 * n_c interpolation weights toward level+1 (None on the coarsest).
    pub pc: Option<Vec<f64>>,
    /// Factored banded interior system (coarsest level, mgsolv == 1).
    pub banded: Option<Banded>,
}

impl Level {
    pub fn nf(&self) -> usize {
        self.nx * self.ny * self.nz
    }
}

pub fn zero_boundary(nx: usize, ny: usize, nz: usize, x: &mut [f64]) {
    let nxny = nx * ny;
    for k in 0..nz {
        for j in 0..ny {
            x[j * nx + k * nxny] = 0.0;
            x[nx - 1 + j * nx + k * nxny] = 0.0;
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            x[i + k * nxny] = 0.0;
            x[i + (ny - 1) * nx + k * nxny] = 0.0;
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            x[i + j * nx] = 0.0;
            x[i + j * nx + (nz - 1) * nxny] = 0.0;
        }
    }
}

/// Smooth one level via Vsmooth/Vgsrb. When `iresid` is set the residual is
/// returned in `r_out`.
#[allow(clippy::too_many_arguments)]
fn smooth_level(
    level: &Level,
    x: &mut [f64],
    fc: &[f64],
    nu: i32,
    iresid: i32,
    iadjoint: i32,
    r_out: Option<&mut [f64]>,
) {
    let nf = level.nf();
    let mut iters = 0i32;
    let mut dummy_r = vec![0.0f64; nf];
    let r: &mut [f64] = match r_out {
        Some(r) => r,
        None => &mut dummy_r,
    };
    let mut w1 = vec![0.0f64; nf];
    let mut w2 = vec![0.0f64; nf];
    crate::gs::gsrb(
        level.nx, level.ny, level.nz,
        &[], &[],
        &level.ac, &level.cc, fc,
        x, &mut w1, &mut w2, r,
        nu, &mut iters, 0.0, 0.0, iresid, iadjoint, level.numdia,
    );
}

/// Apply A * x -> y on one level (Vmatvec dispatch by stencil).
fn matvec_level(level: &Level, x: &[f64], y: &mut [f64]) {
    let nf = level.nf();
    if level.numdia == 4 {
        let _ = split_bands4(&level.ac, nf);
        crate::matvec::matvec(
            level.nx, level.ny, level.nz,
            &[], &[],
            &level.ac, &level.cc, &[],
            x, y,
        );
    } else {
        let b: Bands14 = split_bands14(&level.ac, nf).expect("bands14");
        crate::matvec::matvec27_c(level.nx, level.ny, level.nz, &b, &level.cc, x, y);
    }
}

/// Residual r = f - A x on the finest level (Vmresid).
fn mresid_level0(level: &Level, x: &[f64], r: &mut [f64]) {
    let nf = level.nf();
    crate::blas::mresid(
        level.nx, level.ny, level.nz,
        &[], &[],
        &level.ac[..nf],
        &level.ac[nf..2 * nf],
        &level.ac[2 * nf..3 * nf],
        &level.ac[3 * nf..4 * nf],
        &level.cc,
        x, &level.fc, r,
    );
}

/// Port of Vmvcs. `u` holds the finest-level solution (C's level-1 `x`).
/// Returns the number of outer iterations performed.
pub fn mvcs(
    u: &mut [f64],
    levels: &[Level],
    nu1: i32,
    nu2: i32,
    itmax: i32,
    errtol: f64,
    mgsolv: i32,
    epsiln: f64,
) -> i32 {
    let nlev = levels.len();
    assert!(nlev >= 1, "mvcs needs at least one level");
    let fine = &levels[0];

    // Per-level scratch (x for levels >= 1; w0/w1/w2 per level).
    let mut xs: Vec<Vec<f64>> = levels.iter().map(|l| vec![0.0f64; l.nf()]).collect();
    let mut w0s: Vec<Vec<f64>> = levels.iter().map(|l| vec![0.0f64; l.nf()]).collect();
    let mut w1s: Vec<Vec<f64>> = levels.iter().map(|l| vec![0.0f64; l.nf()]).collect();
    let mut w2s: Vec<Vec<f64>> = levels.iter().map(|l| vec![0.0f64; l.nf()]).collect();

    if nlev == 1 {
        // Vmvcs with nlev == 1: a single coarsest-style solve.
        solve_coarsest(levels, mgsolv, epsiln, &mut xs, &w0s[nlev - 1], &mut w1s[nlev - 1]);
        u.copy_from_slice(&xs[0]);
        return 1;
    }

    // Stopping-test denominator (istop == 1: L1 norm of the fine rhs).
    let mut rsden = crate::blas::xnrm1_grid(fine.nx, fine.ny, fine.nz, &fine.fc);
    if rsden == 0.0 {
        rsden = 1.0;
    }
    let mut iters = 0i32;

    loop {
        // nu1 pre-smoothings on the fine grid, with residual.
        smooth_level(&levels[0], u, &levels[0].fc, nu1, 1, 0, Some(&mut w1s[0]));
        w0s[0].copy_from_slice(&w1s[0]);

        // Down-cycle: restrict the residual and pre-smooth coarser levels.
        for lev in 1..nlev {
            let prev = lev - 1;
            zero_boundary(
                levels[prev].nx, levels[prev].ny, levels[prev].nz,
                &mut w1s[prev],
            );
            crate::matvec::restrc(
                levels[prev].nx, levels[prev].ny, levels[prev].nz,
                levels[lev].nx, levels[lev].ny, levels[lev].nz,
                &w1s[prev], &mut w0s[lev],
                levels[prev].pc.as_deref().expect("pc"),
            );
            if lev != nlev - 1 {
                xs[lev].fill(0.0);
                smooth_level(&levels[lev], &mut xs[lev], &w0s[lev], nu1, 1, 0, Some(&mut w1s[lev]));
            }
        }

        // Coarsest level solve.
        {
            let last = nlev - 1;
            solve_coarsest(levels, mgsolv, epsiln, &mut xs, &w0s[last], &mut w1s[last]);
        }

        // Up-cycle: interpolate the correction, apply the damping parameter,
        // correct and post-smooth.
        for lev in (0..nlev - 1).rev() {
            let coarser = lev + 1;
            crate::matvec::interp_pmg(
                levels[coarser].nx, levels[coarser].ny, levels[coarser].nz,
                levels[lev].nx, levels[lev].ny, levels[lev].nz,
                &xs[coarser], &mut w1s[lev],
                levels[lev].pc.as_deref().expect("pc"),
            );

            // Hackbusch/Reusken damping = CG steplength on the coarser level.
            matvec_level(&levels[coarser], &xs[coarser], &mut w2s[coarser]);
            let xnum = crate::blas::xdot_grid(
                levels[coarser].nx, levels[coarser].ny, levels[coarser].nz,
                &xs[coarser], &w0s[coarser],
            );
            let xden = crate::blas::xdot_grid(
                levels[coarser].nx, levels[coarser].ny, levels[coarser].nz,
                &xs[coarser], &w2s[coarser],
            );
            let xdamp = xnum / xden;

            if lev == 0 {
                crate::blas::xaxpy_grid(
                    levels[0].nx, levels[0].ny, levels[0].nz,
                    xdamp, &w1s[0], u,
                );
                smooth_level(&levels[0], u, &levels[0].fc, nu2, 0, 1, None);
            } else {
                crate::blas::xaxpy_grid(
                    levels[lev].nx, levels[lev].ny, levels[lev].nz,
                    xdamp, &w1s[lev], &mut xs[lev],
                );
                smooth_level(&levels[lev], &mut xs[lev], &w0s[lev], nu2, 0, 1, None);
            }
        }

        // Iteration complete: istop == 1 stopping test on the fine residual.
        iters += 1;
        mresid_level0(&levels[0], u, &mut w1s[0]);
        let rsnrm = crate::blas::xnrm1_grid(levels[0].nx, levels[0].ny, levels[0].nz, &w1s[0]);
        if !(iters < itmax && rsnrm / rsden > errtol) {
            break;
        }
    }

    iters
}

/// Coarsest-level solve (Vmvcs's mgsolv == 1 direct / mgsolv == 0 CG path).
/// The right-hand side is `w0_last` (the restricted residual); the result is
/// stored into xs[last].
fn solve_coarsest(
    levels: &[Level],
    mgsolv: i32,
    epsiln: f64,
    xs: &mut [Vec<f64>],
    w0_last: &[f64],
    w1_last: &mut [f64],
) {
    let last = levels.len() - 1;
    let level = &levels[last];
    let nf = level.nf();
    if mgsolv == 1 {
        let banded = level.banded.as_ref().expect("factored coarsest operator");
        w1_last.copy_from_slice(w0_last);
        crate::lapack::dpbsl(&banded.abd, banded.lda, banded.n, banded.m, w1_last);
        xs[last].copy_from_slice(w1_last);
        zero_boundary(level.nx, level.ny, level.nz, &mut xs[last]);
    } else {
        xs[last].fill(0.0);
        let mut iters = 0i32;
        let rinf = crate::blas::xnrm2(nf, w0_last, 0);
        let mut w1 = vec![0.0f64; nf];
        let mut w2 = vec![0.0f64; nf];
        let mut r = vec![0.0f64; nf];
        crate::cg::cg(
            level.nx, level.ny, level.nz,
            &[], &[],
            &level.ac, &level.cc, w0_last,
            &mut xs[last], &mut w1, &mut w2, &mut r,
            100, &mut iters, epsiln, rinf,
        );
    }
}

// ---------------------------------------------------------------------------
// Legacy single-V-cycle entry used by the Newton (NPBE) path. This is the
// previous recursive implementation with an injection-built hierarchy and
// CG-steplength correction damping; it is kept so the nonlinear branch
// behaves exactly as before.
// ---------------------------------------------------------------------------

/// Linear V-cycle multigrid solver (legacy interface)
#[allow(clippy::too_many_arguments)]
pub fn mgcs(
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    ac_off: usize, cc_off: usize, fc_off: usize,
    pc: &[f64], _iz: &[i32],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu1: i32, nu2: i32,
    omegal: f64,
    _irite: i32,
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
        u,
        w1,
        w2,
        r,
        nu1,
        nu2,
        omegal,
        mgsolv,
    );
}

#[allow(clippy::too_many_arguments)]
fn mgcs_level(
    lev_idx: usize,
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    ac_off: usize, cc_off: usize, fc_off: usize,
    pc_off: usize,
    pc: &[f64],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu1: i32, nu2: i32,
    omegal: f64,
    mgsolv: i32,
) {
    let nf = (nx * ny * nz) as usize;
    let numdia = ipc[lev_idx * 20];

    if nlev == 1 {
        legacy_solve_coarsest(
            nx as usize, ny as usize, nz as usize,
            &ac[ac_off..ac_off + 4 * nf],
            &cc[cc_off..cc_off + nf],
            &fc[fc_off..fc_off + nf],
            u,
            mgsolv,
        );
        return;
    }

    // Pre-smoothing (nu1 iterations, forward red-black order)
    legacy_smooth_on_level(
        nx as usize, ny as usize, nz as usize,
        &ac[ac_off..ac_off + 4 * nf],
        &cc[cc_off..cc_off + nf],
        &fc[fc_off..fc_off + nf],
        u, w1, w2, r,
        nu1, omegal, numdia, 0,
    );

    // Compute residual: r = f - Au
    crate::blas::mresid(
        nx as usize, ny as usize, nz as usize,
        ipc, rpc,
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
        pc,
        &mut u_c_vec,
        w1, w2, r,
        nu1, nu2,
        omegal,
        mgsolv,
    );

    // Prolongate correction to fine grid and add directly.
    let mut corr = vec![0.0; nf];
    prolongate_vec(
        nx as usize, ny as usize, nz as usize,
        nx_c as usize, ny_c as usize, nz_c as usize,
        &u_c_vec, &mut corr, pc_slice,
    );
    let alpha = correction_damping(
        nx_c as usize,
        ny_c as usize,
        nz_c as usize,
        &ipc[(lev_idx + 1) * 20..],
        &rpc[(lev_idx + 1) * 20..],
        ac_coarse,
        &cc[cc_off_c..cc_off_c + nc],
        &u_c_vec,
        &r_c,
    );
    for i in 0..nf {
        u[i] += alpha * corr[i];
    }

    // Post-smoothing (nu2 iterations, adjoint/reversed red-black order)
    legacy_smooth_on_level(
        nx as usize, ny as usize, nz as usize,
        &ac[ac_off..ac_off + 4 * nf],
        &cc[cc_off..cc_off + nf],
        &fc[fc_off..fc_off + nf],
        u, w1, w2, r,
        nu2, omegal, numdia, 1,
    );
}

/// Smooth on current level (legacy)
#[allow(clippy::too_many_arguments)]
fn legacy_smooth_on_level(
    nx: usize, ny: usize, nz: usize,
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu: i32, omega: f64, numdia: i32, iadjoint: i32,
) {
    let n = nx * ny * nz;
    let _ = n;
    let mut iters = 0;
    crate::gs::gsrb(
        nx, ny, nz, &[], &[],
        ac, cc, fc,
        u, w1, w2, r,
        nu, &mut iters, 0.0, omega, 0, iadjoint, numdia,
    );
}

/// Solve on coarsest grid directly (legacy)
fn legacy_solve_coarsest(
    nx: usize, ny: usize, nz: usize,
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
    mgsolv: i32,
) {
    let nf = nx * ny * nz;
    u.fill(0.0);
    if mgsolv == 1 {
        legacy_solve_coarsest_direct(nx, ny, nz, ac, cc, fc, u);
    } else {
        let errtol = 1.0e-6;
        let mut iters = 0;
        let rinf = crate::blas::xnrm2(nf, fc, 0);
        let mut w1 = vec![0.0; nf];
        let mut w2 = vec![0.0; nf];
        let mut r = vec![0.0; nf];
        crate::cg::cg(
            nx, ny, nz, &[], &[],
            ac, cc, fc,
            u, &mut w1, &mut w2, &mut r,
            100, &mut iters, errtol, rinf,
        );
    }
}

fn legacy_solve_coarsest_direct(
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
        for rr in (piv + 1)..nint {
            let cand = a[rr * nint + piv].abs();
            if cand > pivot_abs {
                pivot_abs = cand;
                pivot_row = rr;
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
        for rr in (piv + 1)..nint {
            let factor = a[rr * nint + piv] / diag;
            if factor == 0.0 {
                continue;
            }
            a[rr * nint + piv] = 0.0;
            for c in (piv + 1)..nint {
                a[rr * nint + c] -= factor * a[piv * nint + c];
            }
            b[rr] -= factor * b[piv];
        }
    }

    let mut x = vec![0.0f64; nint];
    for rr in (0..nint).rev() {
        let mut sum = b[rr];
        for c in (rr + 1)..nint {
            sum -= a[rr * nint + c] * x[c];
        }
        let diag = a[rr * nint + rr];
        if diag.abs() <= 1.0e-30 {
            return;
        }
        x[rr] = sum / diag;
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
    use super::{extract_vec, prolongate_vec, restrict_vec};

    #[test]
    fn restrict_vec_preserves_constant_field_without_pc() {
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
    fn prolongate_vec_preserves_constant_field_without_pc() {
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
}
