// APBS PMGC matvec - Matrix-vector operations
// Port of pmgc/matvecd.c

use rayon::prelude::*;

const PAR_THRESHOLD: usize = 16_384;

#[derive(Clone, Copy)]
pub struct Pc27<'a> {
    pub o_pc: &'a [f64],
    pub o_pn: &'a [f64],
    pub o_ps: &'a [f64],
    pub o_pe: &'a [f64],
    pub o_pw: &'a [f64],
    pub o_pne: &'a [f64],
    pub o_pnw: &'a [f64],
    pub o_pse: &'a [f64],
    pub o_psw: &'a [f64],
    pub u_pc: &'a [f64],
    pub u_pn: &'a [f64],
    pub u_ps: &'a [f64],
    pub u_pe: &'a [f64],
    pub u_pw: &'a [f64],
    pub u_pne: &'a [f64],
    pub u_pnw: &'a [f64],
    pub u_pse: &'a [f64],
    pub u_psw: &'a [f64],
    pub d_pc: &'a [f64],
    pub d_pn: &'a [f64],
    pub d_ps: &'a [f64],
    pub d_pe: &'a [f64],
    pub d_pw: &'a [f64],
    pub d_pne: &'a [f64],
    pub d_pnw: &'a [f64],
    pub d_pse: &'a [f64],
    pub d_psw: &'a [f64],
}

fn idx3(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
    i + j * nx + k * nx * ny
}

fn use_pc_transfer() -> bool {
    true
}

fn zero_dirichlet_boundary(nx: usize, ny: usize, nz: usize, x: &mut [f64]) {
    if nx == 0 || ny == 0 || nz == 0 {
        return;
    }
    for k in 0..nz {
        for j in 0..ny {
            x[idx3(0, j, k, nx, ny)] = 0.0;
            x[idx3(nx - 1, j, k, nx, ny)] = 0.0;
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            x[idx3(i, 0, k, nx, ny)] = 0.0;
            x[idx3(i, ny - 1, k, nx, ny)] = 0.0;
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            x[idx3(i, j, 0, nx, ny)] = 0.0;
            x[idx3(i, j, nz - 1, nx, ny)] = 0.0;
        }
    }
}

pub fn split_pc27<'a>(pc: &'a [f64], nxc: usize, nyc: usize, nzc: usize) -> Option<Pc27<'a>> {
    let n = nxc * nyc * nzc;
    if pc.len() < 27 * n {
        return None;
    }
    let mut rem = &pc[..27 * n];
    let take = |rem: &mut &'a [f64]| {
        let (head, tail) = rem.split_at(n);
        *rem = tail;
        head
    };
    Some(Pc27 {
        o_pc: take(&mut rem),
        o_pn: take(&mut rem),
        o_ps: take(&mut rem),
        o_pe: take(&mut rem),
        o_pw: take(&mut rem),
        o_pne: take(&mut rem),
        o_pnw: take(&mut rem),
        o_pse: take(&mut rem),
        o_psw: take(&mut rem),
        u_pc: take(&mut rem),
        u_pn: take(&mut rem),
        u_ps: take(&mut rem),
        u_pe: take(&mut rem),
        u_pw: take(&mut rem),
        u_pne: take(&mut rem),
        u_pnw: take(&mut rem),
        u_pse: take(&mut rem),
        u_psw: take(&mut rem),
        d_pc: take(&mut rem),
        d_pn: take(&mut rem),
        d_ps: take(&mut rem),
        d_pe: take(&mut rem),
        d_pw: take(&mut rem),
        d_pne: take(&mut rem),
        d_pnw: take(&mut rem),
        d_pse: take(&mut rem),
        d_psw: take(&mut rem),
    })
}

/// Matrix-vector product y = Ax using packed diagonal storage
pub fn matvec(
    nx: usize, ny: usize, nz: usize,
    _ipc: &[i32], _rpc: &[f64],
    o_c: &[f64], cc: &[f64], _fc: &[f64],
    x: &[f64], y: &mut [f64],
) {
    let n = nx * ny * nz;
    let nxny = nx * ny;

    // 7-point stencil: o_c is center, o_e is east face, o_n is north face, u_c is up face
    // ac layout: [o_c(nf), o_e(nf), o_n(nf), u_c(nf)]
    let o_e = &o_c[n..2 * n];
    let o_n = &o_c[2 * n..3 * n];
    let u_c = &o_c[3 * n..4 * n];

    y[..n].fill(0.0);
    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }

    if n >= PAR_THRESHOLD {
        y[..n].par_chunks_mut(nxny).enumerate().for_each(|(k, y_plane)| {
            if k == 0 || k + 1 == nz {
                return;
            }
            for j in 1..(ny - 1) {
                for i in 1..(nx - 1) {
                    let local = i + j * nx;
                    let ip = local + k * nxny;
                    let mut sum = o_c[ip] * x[ip];
                    sum -= o_e[ip - 1] * x[ip - 1];
                    sum -= o_e[ip] * x[ip + 1];
                    sum -= o_n[ip - nx] * x[ip - nx];
                    sum -= o_n[ip] * x[ip + nx];
                    sum -= u_c[ip - nxny] * x[ip - nxny];
                    sum -= u_c[ip] * x[ip + nxny];
                    sum += cc[ip] * x[ip];
                    y_plane[local] = sum;
                }
            }
        });
        return;
    }

    for k in 1..(nz - 1) {
        for j in 1..(ny - 1) {
            for i in 1..(nx - 1) {
                let ip = i + j * nx + k * nxny;
                let mut sum = o_c[ip] * x[ip];

                sum -= o_e[ip - 1] * x[ip - 1];
                sum -= o_e[ip] * x[ip + 1];
                sum -= o_n[ip - nx] * x[ip - nx];
                sum -= o_n[ip] * x[ip + nx];
                sum -= u_c[ip - nxny] * x[ip - nxny];
                sum -= u_c[ip] * x[ip + nxny];
                sum += cc[ip] * x[ip];
                y[ip] = sum;
            }
        }
    }
}

/// Restriction operator corresponding to APBS C `Vrestrc`.
///
/// The active Rust PMGC path now threads placeholder `pc` blocks through this
/// entry point. Interior points consume the provided 27-component weights;
/// boundaries still fall back to the older averaging behavior.
pub fn restrc(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    xin: &[f64], xout: &mut [f64],
    _pc: &[f64],
) {
    xout.fill(0.0);
    if use_pc_transfer() {
        if let Some(pc27) = split_pc27(_pc, nxc, nyc, nzc) {
            for kc in 1..nzc.saturating_sub(1) {
                for jc in 1..nyc.saturating_sub(1) {
                    for ic in 1..nxc.saturating_sub(1) {
                        let ipc = idx3(ic, jc, kc, nxc, nyc);
                        let if_ = 2 * ic;
                        let jf = 2 * jc;
                        let kf = 2 * kc;

                        let tmp_o = pc27.o_pc[ipc] * xin[idx3(if_, jf, kf, nxf, nyf)]
                            + pc27.o_pn[ipc] * xin[idx3(if_, jf + 1, kf, nxf, nyf)]
                            + pc27.o_ps[ipc] * xin[idx3(if_, jf - 1, kf, nxf, nyf)]
                            + pc27.o_pe[ipc] * xin[idx3(if_ + 1, jf, kf, nxf, nyf)]
                            + pc27.o_pw[ipc] * xin[idx3(if_ - 1, jf, kf, nxf, nyf)]
                            + pc27.o_pne[ipc] * xin[idx3(if_ + 1, jf + 1, kf, nxf, nyf)]
                            + pc27.o_pnw[ipc] * xin[idx3(if_ - 1, jf + 1, kf, nxf, nyf)]
                            + pc27.o_pse[ipc] * xin[idx3(if_ + 1, jf - 1, kf, nxf, nyf)]
                            + pc27.o_psw[ipc] * xin[idx3(if_ - 1, jf - 1, kf, nxf, nyf)];
                        let tmp_u = pc27.u_pc[ipc] * xin[idx3(if_, jf, kf + 1, nxf, nyf)]
                            + pc27.u_pn[ipc] * xin[idx3(if_, jf + 1, kf + 1, nxf, nyf)]
                            + pc27.u_ps[ipc] * xin[idx3(if_, jf - 1, kf + 1, nxf, nyf)]
                            + pc27.u_pe[ipc] * xin[idx3(if_ + 1, jf, kf + 1, nxf, nyf)]
                            + pc27.u_pw[ipc] * xin[idx3(if_ - 1, jf, kf + 1, nxf, nyf)]
                            + pc27.u_pne[ipc] * xin[idx3(if_ + 1, jf + 1, kf + 1, nxf, nyf)]
                            + pc27.u_pnw[ipc] * xin[idx3(if_ - 1, jf + 1, kf + 1, nxf, nyf)]
                            + pc27.u_pse[ipc] * xin[idx3(if_ + 1, jf - 1, kf + 1, nxf, nyf)]
                            + pc27.u_psw[ipc] * xin[idx3(if_ - 1, jf - 1, kf + 1, nxf, nyf)];
                        let tmp_d = pc27.d_pc[ipc] * xin[idx3(if_, jf, kf - 1, nxf, nyf)]
                            + pc27.d_pn[ipc] * xin[idx3(if_, jf + 1, kf - 1, nxf, nyf)]
                            + pc27.d_ps[ipc] * xin[idx3(if_, jf - 1, kf - 1, nxf, nyf)]
                            + pc27.d_pe[ipc] * xin[idx3(if_ + 1, jf, kf - 1, nxf, nyf)]
                            + pc27.d_pw[ipc] * xin[idx3(if_ - 1, jf, kf - 1, nxf, nyf)]
                            + pc27.d_pne[ipc] * xin[idx3(if_ + 1, jf + 1, kf - 1, nxf, nyf)]
                            + pc27.d_pnw[ipc] * xin[idx3(if_ - 1, jf + 1, kf - 1, nxf, nyf)]
                            + pc27.d_pse[ipc] * xin[idx3(if_ + 1, jf - 1, kf - 1, nxf, nyf)]
                            + pc27.d_psw[ipc] * xin[idx3(if_ - 1, jf - 1, kf - 1, nxf, nyf)];
                        xout[ipc] = tmp_o + tmp_u + tmp_d;
                    }
                }
            }
            zero_dirichlet_boundary(nxc, nyc, nzc, xout);
            return;
        }
    }

    // Legacy fallback path without pc weights.
    let nxnyc = nxc * nyc;
    for kc in 0..nzc {
        for jc in 0..nyc {
            for ic in 0..nxc {
                let ipc_idx = ic + jc * nxc + kc * nxnyc;
                let if_ = 2 * ic;
                let jf = 2 * jc;
                let kf = 2 * kc;
                let mut sum = 0.0;
                let mut count = 0;
                for dk in 0..2 {
                    for dj in 0..2 {
                        for di in 0..2 {
                            let fi = if_ + di;
                            let fj = jf + dj;
                            let fk = kf + dk;
                            if fi < nxf && fj < nyf && fk < nzf {
                                sum += xin[fi + fj * nxf + fk * nxf * nyf];
                                count += 1;
                            }
                        }
                    }
                }
                xout[ipc_idx] = if count > 0 { sum / count as f64 } else { 0.0 };
            }
        }
    }
    zero_dirichlet_boundary(nxc, nyc, nzc, xout);
}

/// Prolongation operator corresponding to APBS C `VinterpPMG`.
///
/// The active Rust PMGC path now threads placeholder `pc` blocks through this
/// entry point. Interior points consume the provided 27-component weights;
/// boundaries still fall back to the older trilinear behavior.
pub fn interp_pmg(
    nxc: usize, nyc: usize, nzc: usize,
    nxf: usize, nyf: usize, nzf: usize,
    xin: &[f64], xout: &mut [f64],
    _pc: &[f64],
) {
    xout.fill(0.0);
    if use_pc_transfer() {
        if let Some(pc27) = split_pc27(_pc, nxc, nyc, nzc) {
            for kc in 0..nzc.saturating_sub(1) {
                for jc in 0..nyc.saturating_sub(1) {
                    for ic in 0..nxc.saturating_sub(1) {
                        let ipc = idx3(ic, jc, kc, nxc, nyc);
                        let fi = 2 * ic;
                        let fj = 2 * jc;
                        let fk = 2 * kc;

                        xout[idx3(fi, fj, fk, nxf, nyf)] = xin[ipc];

                    xout[idx3(fi + 1, fj, fk, nxf, nyf)] =
                        pc27.o_pe[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                        + pc27.o_pw[idx3(ic + 1, jc, kc, nxc, nyc)] * xin[idx3(ic + 1, jc, kc, nxc, nyc)];
                    xout[idx3(fi, fj + 1, fk, nxf, nyf)] =
                        pc27.o_pn[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                        + pc27.o_ps[idx3(ic, jc + 1, kc, nxc, nyc)] * xin[idx3(ic, jc + 1, kc, nxc, nyc)];
                    xout[idx3(fi, fj, fk + 1, nxf, nyf)] =
                        pc27.u_pc[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                        + pc27.d_pc[idx3(ic, jc, kc + 1, nxc, nyc)] * xin[idx3(ic, jc, kc + 1, nxc, nyc)];

                    xout[idx3(fi + 1, fj + 1, fk, nxf, nyf)] =
                        pc27.o_pne[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                        + pc27.o_pnw[idx3(ic + 1, jc, kc, nxc, nyc)] * xin[idx3(ic + 1, jc, kc, nxc, nyc)]
                        + pc27.o_pse[idx3(ic, jc + 1, kc, nxc, nyc)] * xin[idx3(ic, jc + 1, kc, nxc, nyc)]
                        + pc27.o_psw[idx3(ic + 1, jc + 1, kc, nxc, nyc)] * xin[idx3(ic + 1, jc + 1, kc, nxc, nyc)];
                    xout[idx3(fi + 1, fj, fk + 1, nxf, nyf)] =
                        pc27.u_pe[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                        + pc27.u_pw[idx3(ic + 1, jc, kc, nxc, nyc)] * xin[idx3(ic + 1, jc, kc, nxc, nyc)]
                        + pc27.d_pe[idx3(ic, jc, kc + 1, nxc, nyc)] * xin[idx3(ic, jc, kc + 1, nxc, nyc)]
                        + pc27.d_pw[idx3(ic + 1, jc, kc + 1, nxc, nyc)] * xin[idx3(ic + 1, jc, kc + 1, nxc, nyc)];
                    xout[idx3(fi, fj + 1, fk + 1, nxf, nyf)] =
                        pc27.u_pn[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                        + pc27.u_ps[idx3(ic, jc + 1, kc, nxc, nyc)] * xin[idx3(ic, jc + 1, kc, nxc, nyc)]
                        + pc27.d_pn[idx3(ic, jc, kc + 1, nxc, nyc)] * xin[idx3(ic, jc, kc + 1, nxc, nyc)]
                        + pc27.d_ps[idx3(ic, jc + 1, kc + 1, nxc, nyc)] * xin[idx3(ic, jc + 1, kc + 1, nxc, nyc)];
                        xout[idx3(fi + 1, fj + 1, fk + 1, nxf, nyf)] =
                            pc27.u_pne[ipc] * xin[idx3(ic, jc, kc, nxc, nyc)]
                            + pc27.u_pnw[idx3(ic + 1, jc, kc, nxc, nyc)] * xin[idx3(ic + 1, jc, kc, nxc, nyc)]
                            + pc27.u_pse[idx3(ic, jc + 1, kc, nxc, nyc)] * xin[idx3(ic, jc + 1, kc, nxc, nyc)]
                            + pc27.u_psw[idx3(ic + 1, jc + 1, kc, nxc, nyc)] * xin[idx3(ic + 1, jc + 1, kc, nxc, nyc)]
                            + pc27.d_pne[idx3(ic, jc, kc + 1, nxc, nyc)] * xin[idx3(ic, jc, kc + 1, nxc, nyc)]
                            + pc27.d_pnw[idx3(ic + 1, jc, kc + 1, nxc, nyc)] * xin[idx3(ic + 1, jc, kc + 1, nxc, nyc)]
                            + pc27.d_pse[idx3(ic, jc + 1, kc + 1, nxc, nyc)] * xin[idx3(ic, jc + 1, kc + 1, nxc, nyc)]
                            + pc27.d_psw[idx3(ic + 1, jc + 1, kc + 1, nxc, nyc)] * xin[idx3(ic + 1, jc + 1, kc + 1, nxc, nyc)];
                    }
                }
            }
            zero_dirichlet_boundary(nxf, nyf, nzf, xout);
            return;
        }
    }

    // Legacy fallback path without pc weights.
    let nxnyf = nxf * nyf;
    let nxnyc = nxc * nyc;
    for kf in 0..nzf {
        for jf in 0..nyf {
            for if_ in 0..nxf {
                let fc_val = if_ as f64 / 2.0;
                let fj = jf as f64 / 2.0;
                let fk = kf as f64 / 2.0;
                let ic = (fc_val.floor() as usize).min(nxc - 1);
                let jc = (fj.floor() as usize).min(nyc - 1);
                let kc = (fk.floor() as usize).min(nzc - 1);
                let dx = fc_val - ic as f64;
                let dy = fj - jc as f64;
                let dz = fk - kc as f64;
                let w000 = (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
                let w100 = dx * (1.0 - dy) * (1.0 - dz);
                let w010 = (1.0 - dx) * dy * (1.0 - dz);
                let w110 = dx * dy * (1.0 - dz);
                let w001 = (1.0 - dx) * (1.0 - dy) * dz;
                let w101 = dx * (1.0 - dy) * dz;
                let w011 = (1.0 - dx) * dy * dz;
                let w111 = dx * dy * dz;
                let c000 = xin[ic + jc * nxc + kc * nxnyc];
                let c100 = xin[(ic + 1).min(nxc - 1) + jc * nxc + kc * nxnyc];
                let c010 = xin[ic + (jc + 1).min(nyc - 1) * nxc + kc * nxnyc];
                let c110 =
                    xin[(ic + 1).min(nxc - 1) + (jc + 1).min(nyc - 1) * nxc + kc * nxnyc];
                let c001 = xin[ic + jc * nxc + (kc + 1).min(nzc - 1) * nxnyc];
                let c101 =
                    xin[(ic + 1).min(nxc - 1) + jc * nxc + (kc + 1).min(nzc - 1) * nxnyc];
                let c011 =
                    xin[ic + (jc + 1).min(nyc - 1) * nxc + (kc + 1).min(nzc - 1) * nxnyc];
                let c111 = xin[(ic + 1).min(nxc - 1)
                    + (jc + 1).min(nyc - 1) * nxc
                    + (kc + 1).min(nzc - 1) * nxnyc];
                let ip = if_ + jf * nxf + kf * nxnyf;
                xout[ip] = w000 * c000
                    + w100 * c100
                    + w010 * c010
                    + w110 * c110
                    + w001 * c001
                    + w101 * c101
                    + w011 * c011
                    + w111 * c111;
            }
        }
    }
    zero_dirichlet_boundary(nxf, nyf, nzf, xout);
}

#[cfg(test)]
mod tests {
    use super::{idx3, interp_pmg, restrc, split_pc27};

    #[test]
    fn restrc_preserves_constant_field_without_pc() {
        let xin = vec![2.0; 5 * 5 * 5];
        let mut xout = vec![0.0; 3 * 3 * 3];
        restrc(5, 5, 5, 3, 3, 3, &xin, &mut xout, &[]);
        for k in 0..3 {
            for j in 0..3 {
                for i in 0..3 {
                    let v = xout[idx3(i, j, k, 3, 3)];
                    if i == 1 && j == 1 && k == 1 {
                        assert!((v - 2.0).abs() < 1.0e-12);
                    } else {
                        assert!(v.abs() < 1.0e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn interp_pmg_preserves_constant_field_without_pc() {
        let xin = vec![4.0; 3 * 3 * 3];
        let mut xout = vec![0.0; 5 * 5 * 5];
        interp_pmg(3, 3, 3, 5, 5, 5, &xin, &mut xout, &[]);
        for k in 0..5 {
            for j in 0..5 {
                for i in 0..5 {
                    let v = xout[idx3(i, j, k, 5, 5)];
                    if i == 0 || i == 4 || j == 0 || j == 4 || k == 0 || k == 4 {
                        assert!(v.abs() < 1.0e-12);
                    } else {
                        assert!((v - 4.0).abs() < 1.0e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn split_pc27_matches_apbs_component_order() {
        let n = 2usize;
        let pc: Vec<f64> = (0..27 * n).map(|i| i as f64).collect();
        let parts = split_pc27(&pc, 1, 1, n).expect("pc should split");

        assert_eq!(parts.o_pc, &[0.0, 1.0]);
        assert_eq!(parts.o_pn, &[2.0, 3.0]);
        assert_eq!(parts.o_ps, &[4.0, 5.0]);
        assert_eq!(parts.o_pe, &[6.0, 7.0]);
        assert_eq!(parts.o_pw, &[8.0, 9.0]);
        assert_eq!(parts.u_pc, &[18.0, 19.0]);
        assert_eq!(parts.d_pc, &[36.0, 37.0]);
        assert_eq!(parts.d_psw, &[52.0, 53.0]);
    }

    #[test]
    fn split_pc27_rejects_short_buffer() {
        assert!(split_pc27(&vec![0.0; 10], 1, 1, 1).is_none());
    }

    #[test]
    fn interp_pmg_uses_pc_weights_for_x_midpoint() {
        let nxc = 3usize;
        let nyc = 3usize;
        let nzc = 3usize;
        let n = nxc * nyc * nzc;
        let mut pc = vec![0.0; 27 * n];
        let center = idx3(1, 1, 1, nxc, nyc);
        let east = idx3(2, 1, 1, nxc, nyc);
        pc[3 * n + center] = 0.25; // oPE(center)
        pc[4 * n + east] = 0.75;   // oPW(east)

        let mut xin = vec![0.0; n];
        xin[center] = 2.0;
        xin[east] = 10.0;
        let mut xout = vec![0.0; 5 * 5 * 5];
        interp_pmg(nxc, nyc, nzc, 5, 5, 5, &xin, &mut xout, &pc);

        let midpoint = idx3(3, 2, 2, 5, 5);
        assert!((xout[midpoint] - (0.25 * 2.0 + 0.75 * 10.0)).abs() < 1.0e-12);
    }

    #[test]
    fn restrc_uses_pc_weights_for_center_sum() {
        let nxc = 3usize;
        let nyc = 3usize;
        let nzc = 3usize;
        let n = nxc * nyc * nzc;
        let mut pc = vec![0.0; 27 * n];
        let center = idx3(1, 1, 1, nxc, nyc);
        pc[center] = 1.0;         // oPC
        pc[9 * n + center] = 2.0; // uPC
        pc[18 * n + center] = 3.0; // dPC

        let mut xin = vec![0.0; 5 * 5 * 5];
        xin[idx3(2, 2, 2, 5, 5)] = 1.5;
        xin[idx3(2, 2, 3, 5, 5)] = 2.0;
        xin[idx3(2, 2, 1, 5, 5)] = 4.0;
        let mut xout = vec![0.0; n];
        restrc(5, 5, 5, nxc, nyc, nzc, &xin, &mut xout, &pc);

        assert!((xout[center] - (1.0 * 1.5 + 2.0 * 2.0 + 3.0 * 4.0)).abs() < 1.0e-12);
    }
}
