// APBS PMGC buildG - Galerkin coarse-grid operator assembly
// Port of pmgc/buildGd.c (VbuildG_7 / VbuildG_27): A_c = P^T A_f P with the
// 27-component interpolation operator P (pc) built by VbuildP_trilin.
//
// The C versions hand-unroll every coarse stencil band; here the same
// Galerkin product is assembled algorithmically per coarse point, which is
// mathematically identical (P has exactly one 27-point support per coarse
// point, and each fine point belongs to at most two supports).

use crate::matvec::{Bands14, Bands4};

/// Prolongation plane order: matches the Vrestrc2/VinterpPMG2/VbuildPb_trilin
/// component order (oPC, oPN, oPS, oPE, oPW, oPNE, oPNW, oPSE, oPSW, u*, d*).
#[inline]
fn plane_index(d: (isize, isize, isize)) -> usize {
    match d {
        (0, 0, 0) => 0,
        (0, 1, 0) => 1,
        (0, -1, 0) => 2,
        (1, 0, 0) => 3,
        (-1, 0, 0) => 4,
        (1, 1, 0) => 5,
        (-1, 1, 0) => 6,
        (1, -1, 0) => 7,
        (-1, -1, 0) => 8,
        (0, 0, 1) => 9,
        (0, 1, 1) => 10,
        (0, -1, 1) => 11,
        (1, 0, 1) => 12,
        (-1, 0, 1) => 13,
        (1, 1, 1) => 14,
        (-1, 1, 1) => 15,
        (1, -1, 1) => 16,
        (-1, -1, 1) => 17,
        (0, 0, -1) => 18,
        (0, 1, -1) => 19,
        (0, -1, -1) => 20,
        (1, 0, -1) => 21,
        (-1, 0, -1) => 22,
        (1, 1, -1) => 23,
        (-1, 1, -1) => 24,
        (1, -1, -1) => 25,
        (-1, -1, -1) => 26,
        _ => unreachable!(),
    }
}

#[inline]
fn pc_weight(pc: &[f64], nc: usize, c: usize, d: (isize, isize, isize)) -> f64 {
    pc[plane_index(d) * nc + c]
}

/// Fine operator abstraction: 7-point (numdia 4) or 27-point (numdia 14).
pub enum FineOp<'a> {
    Seven(Bands4<'a>),
    TwentySeven(Bands14<'a>),
}

impl<'a> FineOp<'a> {
    #[inline]
    fn diag(&self, a: usize) -> f64 {
        match self {
            FineOp::Seven(b) => b.o_c[a],
            FineOp::TwentySeven(b) => b.o_c[a],
        }
    }

    /// Coefficient A[a][a+delta]; the caller guarantees the neighbor index
    /// stays inside the grid. Returns 0.0 for stencil entries that the
    /// fine discretization does not couple.
    #[inline]
    fn offdiag(&self, a: usize, delta: (isize, isize, isize), nx: usize, _ny: usize) -> f64 {
        let nxny = nx * _ny;
        match self {
            FineOp::Seven(b) => match delta {
                (1, 0, 0) => b.o_e[a],
                (-1, 0, 0) => b.o_e[a - 1],
                (0, 1, 0) => b.o_n[a],
                (0, -1, 0) => b.o_n[a - nx],
                (0, 0, 1) => b.u_c[a],
                (0, 0, -1) => b.u_c[a - nxny],
                _ => 0.0,
            },
            FineOp::TwentySeven(b) => match delta {
                (1, 0, 0) => b.o_e[a],
                (-1, 0, 0) => b.o_e[a - 1],
                (0, 1, 0) => b.o_n[a],
                (0, -1, 0) => b.o_n[a - nx],
                (1, 1, 0) => b.o_ne[a],
                (-1, 1, 0) => b.o_nw[a],
                (1, -1, 0) => b.o_nw[a - nx + 1],
                (-1, -1, 0) => b.o_ne[a - nx - 1],
                (0, 0, 1) => b.u_c[a],
                (0, 0, -1) => b.u_c[a - nxny],
                (1, 0, 1) => b.u_e[a],
                (-1, 0, 1) => b.u_w[a],
                (0, 1, 1) => b.u_n[a],
                (0, -1, 1) => b.u_s[a],
                (1, 1, 1) => b.u_ne[a],
                (-1, 1, 1) => b.u_nw[a],
                (1, -1, 1) => b.u_se[a],
                (-1, -1, 1) => b.u_sw[a],
                (1, 0, -1) => b.u_w[a - nxny + 1],
                (-1, 0, -1) => b.u_e[a - nxny - 1],
                (0, 1, -1) => b.u_s[a - nxny + nx],
                (0, -1, -1) => b.u_n[a - nxny - nx],
                (1, 1, -1) => b.u_sw[a - nxny + nx + 1],
                (-1, 1, -1) => b.u_se[a - nxny + nx - 1],
                (1, -1, -1) => b.u_nw[a - nxny + 1 - nx],
                (-1, -1, -1) => b.u_ne[a - nxny - nx - 1],
                _ => 0.0,
            },
        }
    }
}

/// Per-axis (coarse owner, offset) candidates for fine coordinate t.
#[inline]
fn owners_axis(t: usize) -> ((usize, isize), Option<(usize, isize)>) {
    if t % 2 == 0 {
        ((t / 2, 0), None)
    } else {
        (((t - 1) / 2, 1), Some(((t + 1) / 2, -1)))
    }
}

#[inline]
fn loc_index(e: (isize, isize, isize)) -> usize {
    ((e.0 + 1) + 3 * (e.1 + 1) + 9 * (e.2 + 1)) as usize
}

/// The 13 stored off-diagonal directions of the 14-band layout.
const STORED_DIRS: [(isize, isize, isize); 13] = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (-1, 1, 0),
    (1, 0, 1),
    (-1, 0, 1),
    (0, 1, 1),
    (0, -1, 1),
    (1, 1, 1),
    (-1, 1, 1),
    (1, -1, 1),
    (-1, -1, 1),
];

/// Galerkin coarse operator: ac_out (14 bands * nxc*nyc*nzc) = P^T A_f P.
/// Only interior coarse points are written; boundaries stay zero (only
/// interior rows are ever consumed by Vgsrb27x/Vmresid27/Vbuildband).
#[allow(clippy::too_many_arguments)]
pub fn build_galerkin(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    pc: &[f64],
    fine: FineOp,
    ac_out: &mut [f64],
) {
    let nc = nxc * nyc * nzc;
    assert!(ac_out.len() >= 14 * nc, "coarse operator buffer too small");
    ac_out[..14 * nc].fill(0.0);

    for kc in 1..nzc.saturating_sub(1) {
        for jc in 1..nyc.saturating_sub(1) {
            for ic in 1..nxc.saturating_sub(1) {
                let c = ic + jc * nxc + kc * nxc * nyc;
                let mut loc = [0.0f64; 27];

                for dz in -1isize..=1 {
                    for dy in -1isize..=1 {
                        for dx in -1isize..=1 {
                            let w = pc_weight(pc, nc, c, (dx, dy, dz));
                            if w == 0.0 {
                                continue;
                            }
                            let ax = 2 * ic as isize + dx;
                            let ay = 2 * jc as isize + dy;
                            let az = 2 * kc as isize + dz;
                            if ax < 0 || ay < 0 || az < 0
                                || ax >= nxf as isize || ay >= nyf as isize || az >= nzf as isize
                            {
                                continue;
                            }
                            let a =
                                (ax as usize) + (ay as usize) * nxf + (az as usize) * nxf * nyf;

                            // p == q terms of P^T A P: the fine diagonal A[a][a]
                            // couples the weights of EVERY coarse support that
                            // contains a (supports overlap at odd fine points),
                            // so distribute to each owner of a.
                            {
                                let (oa1, oa2) = owners_axis(ax as usize);
                                let (ob1, ob2) = owners_axis(ay as usize);
                                let (oc1, oc2) = owners_axis(az as usize);
                                let oxs = [Some(oa1), oa2];
                                let oys = [Some(ob1), ob2];
                                let ozs = [Some(oc1), oc2];
                                for &(c2x, dpx) in oxs.iter().flatten() {
                                    for &(c2y, dpy) in oys.iter().flatten() {
                                        for &(c2z, dpz) in ozs.iter().flatten() {
                                            if c2x >= nxc || c2y >= nyc || c2z >= nzc {
                                                continue;
                                            }
                                            let c2 = c2x + c2y * nxc + c2z * nxc * nyc;
                                            let e = (
                                                c2x as isize - ic as isize,
                                                c2y as isize - jc as isize,
                                                c2z as isize - kc as isize,
                                            );
                                            let w2 = pc_weight(pc, nc, c2, (dpx, dpy, dpz));
                                            loc[loc_index(e)] += w * fine.diag(a) * w2;
                                        }
                                    }
                                }
                            }

                            for (ddx, ddy, ddz) in OFFDIAG_DIRS {
                                let bx = ax + ddx;
                                let by = ay + ddy;
                                let bz = az + ddz;
                                if bx < 0 || by < 0 || bz < 0
                                    || bx >= nxf as isize
                                    || by >= nyf as isize
                                    || bz >= nzf as isize
                                {
                                    continue;
                                }
                                // Bands store |A| (matvec subtracts), so the
                                // true operator coefficient is the negated band.
                                let coef = -fine.offdiag(
                                    a,
                                    (ddx, ddy, ddz),
                                    nxf,
                                    nyf,
                                );
                                if coef == 0.0 {
                                    continue;
                                }

                                let (ox1, ox2) = owners_axis(bx as usize);
                                let (oy1, oy2) = owners_axis(by as usize);
                                let (oz1, oz2) = owners_axis(bz as usize);
                                let xs = [Some(ox1), ox2];
                                let ys = [Some(oy1), oy2];
                                let zs = [Some(oz1), oz2];
                                for &(cxx, dpx) in xs.iter().flatten() {
                                    for &(cyy, dpy) in ys.iter().flatten() {
                                        for &(czz, dpz) in zs.iter().flatten() {
                                            if cxx >= nxc || cyy >= nyc || czz >= nzc {
                                                continue;
                                            }
                                            let c2 = cxx
                                                + cyy * nxc
                                                + czz * nxc * nyc;
                                            let e = (
                                                cxx as isize - ic as isize,
                                                cyy as isize - jc as isize,
                                                czz as isize - kc as isize,
                                            );
                                            let w2 = pc_weight(pc, nc, c2, (dpx, dpy, dpz));
                                            loc[loc_index(e)] += w * coef * w2;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Store the 14 bands for this coarse point. The loc array
                // holds the true coarse matrix entries; off-diagonal bands are
                // stored negated (positive), matching the fine-band convention.
                let base = |band: usize| band * nc + c;
                ac_out[base(0)] = loc[loc_index((0, 0, 0))];
                for (band, dir) in STORED_DIRS.iter().enumerate() {
                    ac_out[base(band + 1)] = -loc[loc_index(*dir)];
                }
            }
        }
    }
}

/// All 26 off-stencil directions of the 27-point neighborhood.
const OFFDIAG_DIRS: [(isize, isize, isize); 26] = [
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
    (1, 1, 0),
    (-1, 1, 0),
    (1, -1, 0),
    (-1, -1, 0),
    (1, 0, 1),
    (-1, 0, 1),
    (0, 1, 1),
    (0, -1, 1),
    (1, 1, 1),
    (-1, 1, 1),
    (1, -1, 1),
    (-1, -1, 1),
    (1, 0, -1),
    (-1, 0, -1),
    (0, 1, -1),
    (0, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (1, -1, -1),
    (-1, -1, -1),
];

#[cfg(test)]
mod tests {
    use super::{build_galerkin, FineOp};
    use crate::matvec::{matvec27_c, split_bands14, split_bands4, Bands4};

    /// For a constant 7-point operator (oC=6, oE=oN=uC=-1) the Galerkin
    /// product with the constant trilinear P must reproduce the classic
    /// 27-point constant-coefficient stencil, and matvec27 must agree
    /// with the direct P^T A P application on a random vector.
    #[test]
    fn galerkin_constant_operator_matvec_consistency() {
        let (nxf, nyf, nzf) = (9usize, 9usize, 9usize);
        let (nxc, nyc, nzc) = (5usize, 5usize, 5usize);
        let nf = nxf * nyf * nzf;
        let nc = nxc * nyc * nzc;

        // Positive-convention bands (matvec subtracts off-diagonals).
        let mut ac_f = vec![0.0f64; 4 * nf];
        for i in 0..nf {
            ac_f[i] = 6.0;
            ac_f[nf + i] = 1.0;
            ac_f[2 * nf + i] = 1.0;
            ac_f[3 * nf + i] = 1.0;
        }
        let pc = crate::build_p::build_p_trilin_block(nxc, nyc, nzc);
        let mut ac_c = vec![0.0f64; 14 * nc];
        build_galerkin(
            nxf, nyf, nzf, nxc, nyc, nzc,
            &pc,
            FineOp::Seven(Bands4 {
                o_c: &ac_f[0..nf],
                o_e: &ac_f[nf..2 * nf],
                o_n: &ac_f[2 * nf..3 * nf],
                u_c: &ac_f[3 * nf..4 * nf],
            }),
            &mut ac_c,
        );

        // Random-ish fine vector, prolongate manually via interp_pmg.
        let mut xc = vec![0.0f64; nc];
        for k in 1..nzc - 1 {
            for j in 1..nyc - 1 {
                for i in 1..nxc - 1 {
                    let ip = i + j * nxc + k * nxc * nyc;
                    xc[ip] = ((ip * 37 + 11) % 17) as f64 - 8.0;
                }
            }
        }
        let mut xf = vec![0.0f64; nf];
        crate::matvec::interp_pmg(nxc, nyc, nzc, nxf, nyf, nzf, &xc, &mut xf, &pc);

        // Apply fine operator (with cc = 0) at interior points.
        let b4 = split_bands4(&ac_f, nf).unwrap();
        let mut af_x = vec![0.0f64; nf];
        for k in 1..nzf - 1 {
            for j in 1..nyf - 1 {
                for i in 1..nxf - 1 {
                    let ip = i + j * nxf + k * nxf * nyf;
                    af_x[ip] = 6.0 * xf[ip]
                        - b4.o_e[ip] * xf[ip + 1]
                        - b4.o_e[ip - 1] * xf[ip - 1]
                        - b4.o_n[ip] * xf[ip + nxf]
                        - b4.o_n[ip - nxf] * xf[ip - nxf]
                        - b4.u_c[ip] * xf[ip + nxf * nyf]
                        - b4.u_c[ip - nxf * nyf] * xf[ip - nxf * nyf];
                }
            }
        }

        // Restrict to coarse: expect equals A_c * xc at interior points.
        let mut af_c = vec![0.0f64; nc];
        crate::matvec::restrc(nxf, nyf, nzf, nxc, nyc, nzc, &af_x, &mut af_c, &pc);

        let b14 = split_bands14(&ac_c, nc).unwrap();
        let mut ac_x = vec![0.0f64; nc];
        matvec27_c(nxc, nyc, nzc, &b14, &vec![0.0; nc], &xc, &mut ac_x);

        for k in 1..nzc - 1 {
            for j in 1..nyc - 1 {
                for i in 1..nxc - 1 {
                    let ip = i + j * nxc + k * nxc * nyc;
                    let scale = af_c[ip].abs().max(1.0);
                    assert!(
                        (af_c[ip] - ac_x[ip]).abs() / scale < 1.0e-12,
                        "mismatch at {}: {} vs {}",
                        ip,
                        af_c[ip],
                        ac_x[ip]
                    );
                }
            }
        }
    }
}
