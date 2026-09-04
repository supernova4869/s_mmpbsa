// APBS PMGC mypdec - PDE-specific routines
// Port of pmgc/mypdec.c

/// Compute the Hackbusch-Reusken damping parameter for coarse-grid correction
pub fn hackbusch_reusken(
    nx: usize, ny: usize, nz: usize,
    o_c: &[f64], o_e: &[f64], o_n: &[f64], u_c: &[f64],
    cc: &[f64],
) -> f64 {
    // Estimate the spectral radius of the iteration matrix
    // using the formula: omega = 2 / (1 + sqrt(1 - rho^2))
    // where rho is the spectral radius of the Gauss-Seidel iteration matrix

    let nf = nx * ny * nz;

    // Compute ratio of off-diagonal to diagonal norms
    let mut diag_sum = 0.0;
    let mut off_sum = 0.0;

    for i in 0..nf {
        diag_sum += o_c[i].abs() + cc[i].abs();
        if i > 0 { off_sum += o_e[i - 1].abs(); }
        if i < nf - 1 { off_sum += o_e[i].abs(); }
        if i >= nx { off_sum += o_n[i - nx].abs(); }
        if i < nf - nx { off_sum += o_n[i].abs(); }
        if i >= nx * ny { off_sum += u_c[i - nx * ny].abs(); }
        if i < nf - nx * ny { off_sum += u_c[i].abs(); }
    }

    let rho = if diag_sum > 1.0e-20 {
        (off_sum / diag_sum).min(0.999)
    } else {
        0.5
    };

    // Optimal omega for SOR
    2.0 / (1.0 + (1.0 - rho * rho).sqrt())
}

/// Compute the residual for the nonlinear PB equation
pub fn pb_residual(
    nx: usize, ny: usize, nz: usize,
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
    u: &[f64], r: &mut [f64],
    hx: f64, hy: f64, hz: f64,
    _zkappa2: f64,
) {
    let nxny = nx * ny;

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let ip = i + j * nx + k * nxny;

                let mut lap = 0.0;

                // x-direction
                if i > 0 {
                    let a_w = 0.5 * (a1cf[ip - 1] + a1cf[ip]);
                    lap -= a_w * (u[ip] - u[ip - 1]) / (hx * hx);
                }
                if i < nx - 1 {
                    let a_e = 0.5 * (a1cf[ip] + a1cf[ip + 1]);
                    lap -= a_e * (u[ip] - u[ip + 1]) / (hx * hx);
                }

                // y-direction
                if j > 0 {
                    let a_s = 0.5 * (a2cf[ip - nx] + a2cf[ip]);
                    lap -= a_s * (u[ip] - u[ip - nx]) / (hy * hy);
                }
                if j < ny - 1 {
                    let a_n = 0.5 * (a2cf[ip] + a2cf[ip + nx]);
                    lap -= a_n * (u[ip] - u[ip + nx]) / (hy * hy);
                }

                // z-direction
                if k > 0 {
                    let a_d = 0.5 * (a3cf[ip - nxny] + a3cf[ip]);
                    lap -= a_d * (u[ip] - u[ip - nxny]) / (hz * hz);
                }
                if k < nz - 1 {
                    let a_u = 0.5 * (a3cf[ip] + a3cf[ip + nxny]);
                    lap -= a_u * (u[ip] - u[ip + nxny]) / (hz * hz);
                }

                // Reaction term: ccf * sinh(u)
                r[ip] = lap + ccf[ip] * u[ip].sinh() - fcf[ip];
            }
        }
    }
}
