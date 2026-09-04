// APBS PMGC buildA - Build operator A (finite volume discretization)
// Port of pmgc/buildAd.c

/// Build the finite volume discretization operator
pub fn build_a(
    nx: usize, ny: usize, nz: usize,
    ipkey: i32, mgdisc: i32, numdia: i32,
    ac: &mut [f64],   // operator storage: [o_c(nf), o_e(nf), o_n(nf), u_c(nf)]
    cc: &mut [f64], fc: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
) {
    if mgdisc == 0 {
        build_a_fv(
            nx, ny, nz, ipkey, numdia,
            ac, cc, fc,
            xf, yf, zf, gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf, ccf, fcf,
        );
    } else {
        build_a_fe(
            nx, ny, nz, ipkey, numdia,
            ac, cc, fc,
            xf, yf, zf, gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf, ccf, fcf,
        );
    }
}

/// Box method (finite volume) discretization.
///
/// This follows APBS C `VbuildA_fv`: only interior unknowns are assembled and
/// Dirichlet boundary values are folded into the RHS.
fn build_a_fv(
    nx: usize, ny: usize, nz: usize,
    _ipkey: i32, _numdia: i32,
    ac: &mut [f64], cc: &mut [f64], fc: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
) {
    let nxny = nx * ny;
    let nf = nx * ny * nz;

    ac[..4 * nf].fill(0.0);
    cc[..nf].fill(0.0);
    fc[..nf].fill(0.0);

    let ijkx = |j: usize, k: usize, face: usize, ny_: usize, nz_: usize| -> usize {
        face * ny_ * nz_ + k * ny_ + j
    };
    let ijky = |i: usize, k: usize, face: usize, nx_: usize, nz_: usize| -> usize {
        face * nx_ * nz_ + k * nx_ + i
    };
    let ijkz = |i: usize, j: usize, face: usize, nx_: usize, ny_: usize| -> usize {
        face * nx_ * ny_ + j * nx_ + i
    };

    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }

    for k in 1..(nz - 1) {
        let hzm1 = zf[k] - zf[k - 1];
        let hz = zf[k + 1] - zf[k];

        for j in 1..(ny - 1) {
            let hym1 = yf[j] - yf[j - 1];
            let hy = yf[j + 1] - yf[j];

            for i in 1..(nx - 1) {
                let hxm1 = xf[i] - xf[i - 1];
                let hx = xf[i + 1] - xf[i];
                let ip = i + j * nx + k * nxny;

                let coef_oe = (hym1 + hy) * (hzm1 + hz) / (4.0 * hx);
                let coef_oem1 = (hym1 + hy) * (hzm1 + hz) / (4.0 * hxm1);
                let coef_on = (hxm1 + hx) * (hzm1 + hz) / (4.0 * hy);
                let coef_onm1 = (hxm1 + hx) * (hzm1 + hz) / (4.0 * hym1);
                let coef_uc = (hxm1 + hx) * (hym1 + hy) / (4.0 * hz);
                let coef_ucm1 = (hxm1 + hx) * (hym1 + hy) / (4.0 * hzm1);
                let coef_fc = (hxm1 + hx) * (hym1 + hy) * (hzm1 + hz) / 8.0;

                fc[ip] = coef_fc * fcf[ip];
                cc[ip] = coef_fc * ccf[ip];

                ac[ip] = coef_oe * a1cf[ip]
                    + coef_oem1 * a1cf[ip - 1]
                    + coef_on * a2cf[ip]
                    + coef_onm1 * a2cf[ip - nx]
                    + coef_uc * a3cf[ip]
                    + coef_ucm1 * a3cf[ip - nxny];

                let east_is_interior = i < nx - 2;
                ac[nf + ip] = if east_is_interior { coef_oe * a1cf[ip] } else { 0.0 };
                if !east_is_interior {
                    fc[ip] += coef_oe * a1cf[ip] * gxcf[ijkx(j, k, 1, ny, nz)];
                }

                let north_is_interior = j < ny - 2;
                ac[2 * nf + ip] = if north_is_interior { coef_on * a2cf[ip] } else { 0.0 };
                if !north_is_interior {
                    fc[ip] += coef_on * a2cf[ip] * gycf[ijky(i, k, 1, nx, nz)];
                }

                let up_is_interior = k < nz - 2;
                ac[3 * nf + ip] = if up_is_interior { coef_uc * a3cf[ip] } else { 0.0 };
                if !up_is_interior {
                    fc[ip] += coef_uc * a3cf[ip] * gzcf[ijkz(i, j, 1, nx, ny)];
                }

                if i == 1 {
                    fc[ip] += coef_oem1 * a1cf[ip - 1] * gxcf[ijkx(j, k, 0, ny, nz)];
                }
                if j == 1 {
                    fc[ip] += coef_onm1 * a2cf[ip - nx] * gycf[ijky(i, k, 0, nx, nz)];
                }
                if k == 1 {
                    fc[ip] += coef_ucm1 * a3cf[ip - nxny] * gzcf[ijkz(i, j, 0, nx, ny)];
                }
            }
        }
    }
}

/// Finite element method discretization (delegates to FV for now)
fn build_a_fe(
    nx: usize, ny: usize, nz: usize,
    ipkey: i32, numdia: i32,
    ac: &mut [f64], cc: &mut [f64], fc: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
) {
    build_a_fv(
        nx, ny, nz, ipkey, numdia,
        ac, cc, fc,
        xf, yf, zf, gxcf, gycf, gzcf,
        a1cf, a2cf, a3cf, ccf, fcf,
    );
}
