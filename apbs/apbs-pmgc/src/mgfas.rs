// APBS PMGC mgfas - Full Approximation Scheme (nonlinear V-cycle)
// Port of pmgc/mgfasd.c

/// Nonlinear FAS V-cycle for the Poisson-Boltzmann equation
pub fn mgfas(
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    pc: &[f64],
    iz: &[i32],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu1: i32, nu2: i32,
    omegan: f64,
    irite: i32,
) {
    mgfas_level(
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
        0,
        0,
        0,
        0,
        pc,
        iz,
        u,
        w1,
        w2,
        r,
        nu1,
        nu2,
        omegan,
        irite,
    );
}

fn mgfas_level(
    lev_idx: usize,
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    ac_off: usize, cc_off: usize, fc_off: usize,
    pc_off: usize,
    pc: &[f64],
    iz: &[i32],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu1: i32, nu2: i32,
    omegan: f64,
    irite: i32,
) {
    let current_ipc = &ipc[lev_idx * 20..];
    let current_rpc = &rpc[lev_idx * 20..];
    let nf = (nx * ny * nz) as usize;
    let ac_cur = &ac[ac_off..ac_off + 4 * nf];
    let cc_cur = &cc[cc_off..cc_off + nf];
    let fc_cur = &fc[fc_off..fc_off + nf];

    if nlev == 1 {
        solve_coarsest_nonlinear(
            nx as usize, ny as usize, nz as usize,
            current_ipc, current_rpc, ac_cur, cc_cur, fc_cur,
            u, w1, w2, r,
        );
        return;
    }

    // Pre-smoothing
    smooth_nonlinear(
        nx as usize, ny as usize, nz as usize,
        current_ipc, current_rpc, ac_cur, cc_cur, fc_cur,
        u, w1, w2, r,
        nu1, omegan,
    );

    // Compute defect: r = f - N(u)
    compute_nonlinear_residual(
        nx as usize, ny as usize, nz as usize,
        current_ipc, current_rpc, ac_cur, cc_cur, fc_cur,
        u, r,
    );

    // Restrict defect and solution to coarse grid
    let (nx_c, ny_c, nz_c) = crate::build_str::make_coarse(nx, ny, nz);
    let nc = (nx_c * ny_c * nz_c) as usize;
    let pc_len = 27 * nc;
    let pc_slice = if pc_off + pc_len <= pc.len() {
        &pc[pc_off..pc_off + pc_len]
    } else {
        &[]
    };

    let mut u_c = vec![0.0; nc];
    let mut r_c = vec![0.0; nc];
    let mut fc_c = vec![0.0; nc];
    let mut cc_c = vec![0.0; nc];

    restrict_fas(
        nx as usize, ny as usize, nz as usize,
        nx_c as usize, ny_c as usize, nz_c as usize,
        u, r, fc_cur, cc_cur,
        &mut u_c, &mut r_c, &mut fc_c, &mut cc_c, pc_slice,
    );
    let u_c_initial = u_c.clone();

    // FAS correction: modify coarse RHS
    // fc_c = r_c + N_c(u_c) - Ac * u_c
    // (simplified: use restricted fine-grid operator on coarse solution)
    let nf_c = nc;
    let mut ac_u_c = vec![0.0; nf_c];
    let ac_off_c = ac_off + 4 * nf;
    let ac_coarse = &ac[ac_off_c..ac_off_c + 4 * nf_c];
    apply_operator(
        nx_c as usize, ny_c as usize, nz_c as usize,
        ac_coarse, &u_c, &mut ac_u_c,
    );
    for i in 0..nc {
        fc_c[i] = r_c[i] + fc_c[i] - ac_u_c[i];
    }

    mgfas_level(
        lev_idx + 1,
        nlev - 1,
        nx_c, ny_c, nz_c,
        ipc, rpc,
        ac, &cc_c, &fc_c,
        ac_off_c, 0, 0,
        pc_off + pc_len,
        pc, iz,
        &mut u_c,
        w1, w2, r,
        nu1, nu2,
        omegan,
        irite,
    );

    // Prolongate coarse correction
    let mut coarse_delta = u_c.clone();
    for i in 0..nc {
        coarse_delta[i] -= u_c_initial[i];
    }

    let mut corr = vec![0.0; nf];
    prolongate_fas(
        nx as usize, ny as usize, nz as usize,
        nx_c as usize, ny_c as usize, nz_c as usize,
        &coarse_delta, &mut corr, pc_slice,
    );

    // Add correction
    for i in 0..nf {
        u[i] += corr[i];
    }

    // Post-smoothing
    smooth_nonlinear(
        nx as usize, ny as usize, nz as usize,
        current_ipc, current_rpc, ac_cur, cc_cur, fc_cur,
        u, w1, w2, r,
        nu2, omegan,
    );
}

/// Compute nonlinear residual r = f - N(u) where N includes PB nonlinearity
pub(crate) fn compute_nonlinear_residual(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &[f64], r: &mut [f64],
) {
    let nf = nx * ny * nz;
    let o_c = &ac[0..nf];
    let o_e = &ac[nf..2 * nf];
    let o_n = &ac[2 * nf..3 * nf];
    let u_c = &ac[3 * nf..4 * nf];

    // r = f - Au (linear part)
    crate::blas::mresid(nx, ny, nz, ipc, rpc, o_c, o_e, o_n, u_c, cc, u, fc, r);

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let ip = i + j * nx + k * nx * ny;
                // mresid already includes the linear +cc*u contribution on the diagonal.
                // For nonlinear PB, add only the nonlinear increment: cc*(sinh(u)-u).
                let u_val = u[ip].clamp(SINH_MIN, SINH_MAX);
                let sinh_minus_u = u_val.sinh() - u_val;
                r[ip] -= cc[ip] * sinh_minus_u;
            }
        }
    }
}

/// Smooth nonlinear equation
fn smooth_nonlinear(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    nu: i32, omega: f64,
) {
    let nf = nx * ny * nz;
    let numdia = ipc[0];

    // Use GSRB with Newton-like linearization
    for _ in 0..nu {
        // Compute nonlinear residual
        compute_nonlinear_residual(nx, ny, nz, ipc, rpc, ac, cc, fc, u, r);
        let rhs_nl = r.to_vec();

        // Linearize and solve for correction
        let mut du = vec![0.0; nf];
        crate::smooth::smooth(
            nx, ny, nz, ipc, rpc,
            ac, cc, &rhs_nl,
            &mut du, w1, w2, r,
            numdia, 1, omega,
            1, 0,
        );

        // Update solution
        for i in 0..nf {
            u[i] += du[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{compute_nonlinear_residual, restrict_fas};

    #[test]
    fn nonlinear_residual_uses_local_cc_not_rpc_zkappa2() {
        let nx = 3usize;
        let ny = 3usize;
        let nz = 3usize;
        let nf = nx * ny * nz;
        let ac = vec![0.0; 4 * nf];
        let mut cc = vec![0.0; nf];
        let fc = vec![0.0; nf];
        let mut u = vec![0.0; nf];
        let mut r = vec![0.0; nf];
        let ipc = vec![4];
        let rpc = vec![1.0, 1.0, 1.0, 999.0];
        let center = 1 + nx + nx * ny;
        cc[center] = 2.5;
        u[center] = 1.0;

        compute_nonlinear_residual(nx, ny, nz, &ipc, &rpc, &ac, &cc, &fc, &u, &mut r);

        let expected = -2.5 * 1.0f64.sinh();
        assert!((r[center] - expected).abs() < 1.0e-12);
    }

    #[test]
    fn restrict_fas_extracts_solution_on_coarse_grid() {
        let nxf = 5usize;
        let nyf = 5usize;
        let nzf = 5usize;
        let nxc = 3usize;
        let nyc = 3usize;
        let nzc = 3usize;
        let mut u_f = vec![0.0; nxf * nyf * nzf];
        u_f[2 + 2 * nxf + 2 * nxf * nyf] = 7.5;
        let zeros_f = vec![0.0; nxf * nyf * nzf];
        let mut u_c = vec![0.0; nxc * nyc * nzc];
        let mut r_c = vec![0.0; nxc * nyc * nzc];
        let mut fc_c = vec![0.0; nxc * nyc * nzc];
        let mut cc_c = vec![0.0; nxc * nyc * nzc];

        restrict_fas(
            nxf, nyf, nzf,
            nxc, nyc, nzc,
            &u_f, &zeros_f, &zeros_f, &zeros_f,
            &mut u_c, &mut r_c, &mut fc_c, &mut cc_c, &[],
        );

        let center = 1 + 1 * nxc + 1 * nxc * nyc;
        assert!((u_c[center] - 7.5).abs() < 1.0e-12);
    }
}

/// Solve on coarsest grid (nonlinear)
fn solve_coarsest_nonlinear(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    u: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
) {
    // Use Newton iteration on coarsest grid
    smooth_nonlinear(nx, ny, nz, ipc, rpc, ac, cc, fc, u, w1, w2, r, 100, 1.0);
}

/// Apply operator to vector
fn apply_operator(
    nx: usize, ny: usize, nz: usize,
    ac: &[f64], x: &[f64], y: &mut [f64],
) {
    let nf = nx * ny * nz;
    let o_c = &ac[0..nf];
    let o_e = &ac[nf..2 * nf];
    let o_n = &ac[2 * nf..3 * nf];
    let u_c = &ac[3 * nf..4 * nf];

    let nxny = nx * ny;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let ip = i + j * nx + k * nxny;
                let mut sum = o_c[ip] * x[ip];
                if i > 0 { sum -= o_e[ip - 1] * x[ip - 1]; }
                if i < nx - 1 { sum -= o_e[ip] * x[ip + 1]; }
                if j > 0 { sum -= o_n[ip - nx] * x[ip - nx]; }
                if j < ny - 1 { sum -= o_n[ip] * x[ip + nx]; }
                if k > 0 { sum -= u_c[ip - nxny] * x[ip - nxny]; }
                if k < nz - 1 { sum -= u_c[ip] * x[ip + nxny]; }
                y[ip] = sum;
            }
        }
    }
}

/// Restriction for FAS (restricts solution and defect)
fn restrict_fas(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    u_f: &[f64], r_f: &[f64], fc_f: &[f64], cc_f: &[f64],
    u_c: &mut [f64], r_c: &mut [f64], fc_c: &mut [f64], cc_c: &mut [f64],
    pc: &[f64],
) {
    crate::mgcs::extract_vec(nxf, nyf, nzf, nxc, nyc, nzc, u_f, u_c);
    crate::mgcs::restrict_vec(nxf, nyf, nzf, nxc, nyc, nzc, r_f, r_c, pc);
    crate::mgcs::restrict_vec(nxf, nyf, nzf, nxc, nyc, nzc, fc_f, fc_c, pc);
    crate::mgcs::restrict_vec(nxf, nyf, nzf, nxc, nyc, nzc, cc_f, cc_c, pc);
}

/// Prolongate for FAS
fn prolongate_fas(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    coarse: &[f64], fine: &mut [f64], pc: &[f64],
) {
    crate::mgcs::prolongate_vec(nxf, nyf, nzf, nxc, nyc, nzc, coarse, fine, pc);
}
pub(crate) const SINH_MIN: f64 = -85.0;
pub(crate) const SINH_MAX: f64 = 85.0;
