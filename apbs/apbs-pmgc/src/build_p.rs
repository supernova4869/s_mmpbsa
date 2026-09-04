// APBS PMGC buildP - Prolongation operator builders
// The 7-point path implements the APBS C interior op7 formulas with
// trilinear boundary fallback. The 27-point path is still not implemented.

fn idx3(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
    i + j * nx + k * nx * ny
}

fn set_plane(pc: &mut [f64], plane: usize, n: usize, value: f64) {
    let start = plane * n;
    let end = start + n;
    for v in &mut pc[start..end] {
        *v = value;
    }
}

pub fn build_p_trilin_block(nxc: usize, nyc: usize, nzc: usize) -> Vec<f64> {
    let n = nxc * nyc * nzc;
    let mut pc = vec![0.0; 27 * n];

    set_plane(&mut pc, 0, n, 1.0);      // oPC
    set_plane(&mut pc, 1, n, 0.5);      // oPN
    set_plane(&mut pc, 2, n, 0.5);      // oPS
    set_plane(&mut pc, 3, n, 0.5);      // oPE
    set_plane(&mut pc, 4, n, 0.5);      // oPW
    set_plane(&mut pc, 5, n, 0.25);     // oPNE
    set_plane(&mut pc, 6, n, 0.25);     // oPNW
    set_plane(&mut pc, 7, n, 0.25);     // oPSE
    set_plane(&mut pc, 8, n, 0.25);     // oPSW
    set_plane(&mut pc, 9, n, 0.5);      // uPC
    set_plane(&mut pc, 10, n, 0.25);    // uPN
    set_plane(&mut pc, 11, n, 0.25);    // uPS
    set_plane(&mut pc, 12, n, 0.25);    // uPE
    set_plane(&mut pc, 13, n, 0.25);    // uPW
    set_plane(&mut pc, 14, n, 0.125);   // uPNE
    set_plane(&mut pc, 15, n, 0.125);   // uPNW
    set_plane(&mut pc, 16, n, 0.125);   // uPSE
    set_plane(&mut pc, 17, n, 0.125);   // uPSW
    set_plane(&mut pc, 18, n, 0.5);     // dPC
    set_plane(&mut pc, 19, n, 0.25);    // dPN
    set_plane(&mut pc, 20, n, 0.25);    // dPS
    set_plane(&mut pc, 21, n, 0.25);    // dPE
    set_plane(&mut pc, 22, n, 0.25);    // dPW
    set_plane(&mut pc, 23, n, 0.125);   // dPNE
    set_plane(&mut pc, 24, n, 0.125);   // dPNW
    set_plane(&mut pc, 25, n, 0.125);   // dPSE
    set_plane(&mut pc, 26, n, 0.125);   // dPSW

    pc
}

fn assign_if_finite(dst: &mut f64, numer: f64, denom: f64) {
    if denom.abs() > 1.0e-30 {
        let val = numer / denom;
        if val.is_finite() {
            *dst = val;
        }
    }
}

pub fn build_p_op7_block(nxf: usize, nyf: usize, nzf: usize, ac: &[f64]) -> Vec<f64> {
    let nf = nxf * nyf * nzf;
    assert!(ac.len() >= 4 * nf, "build_p_op7_block requires 4*nf operator entries");

    let (nxc, nyc, nzc) = (nxf / 2 + 1, nyf / 2 + 1, nzf / 2 + 1);
    let nc = nxc * nyc * nzc;
    let mut pc = build_p_trilin_block(nxc, nyc, nzc);

    let (o_c, rest) = ac.split_at(nf);
    let (o_e, rest) = rest.split_at(nf);
    let (o_n, u_c) = rest.split_at(nf);

    let at = |plane: usize, idx: usize| -> usize { plane * nc + idx };

    for kc in 1..nzc.saturating_sub(1) {
        let k = 2 * kc;
        for jc in 1..nyc.saturating_sub(1) {
            let j = 2 * jc;
            for ic in 1..nxc.saturating_sub(1) {
                let i = 2 * ic;

                let ip = idx3(i, j, k, nxf, nyf);
                let im1 = i - 1;
                let ip1 = i + 1;
                let jm1 = j - 1;
                let jp1 = j + 1;
                let km1 = k - 1;
                let kp1 = k + 1;
                let c = idx3(ic, jc, kc, nxc, nyc);
                let o_pc_v = 1.0;
                let mut o_pn_v = pc[at(1, c)];
                let mut o_ps_v = pc[at(2, c)];
                let mut o_pe_v = pc[at(3, c)];
                let mut o_pw_v = pc[at(4, c)];
                let mut o_pne_v = pc[at(5, c)];
                let mut o_pnw_v = pc[at(6, c)];
                let mut o_pse_v = pc[at(7, c)];
                let mut o_psw_v = pc[at(8, c)];
                let mut u_pc_v = pc[at(9, c)];
                let mut u_pn_v = pc[at(10, c)];
                let mut u_ps_v = pc[at(11, c)];
                let mut u_pe_v = pc[at(12, c)];
                let mut u_pw_v = pc[at(13, c)];
                let mut u_pne_v = pc[at(14, c)];
                let mut u_pnw_v = pc[at(15, c)];
                let mut u_pse_v = pc[at(16, c)];
                let mut u_psw_v = pc[at(17, c)];
                let mut d_pc_v = pc[at(18, c)];
                let mut d_pn_v = pc[at(19, c)];
                let mut d_ps_v = pc[at(20, c)];
                let mut d_pe_v = pc[at(21, c)];
                let mut d_pw_v = pc[at(22, c)];
                let mut d_pne_v = pc[at(23, c)];
                let mut d_pnw_v = pc[at(24, c)];
                let mut d_pse_v = pc[at(25, c)];
                let mut d_psw_v = pc[at(26, c)];

                assign_if_finite(
                    &mut o_pn_v,
                    o_n[ip],
                    o_c[idx3(i, jp1, k, nxf, nyf)]
                        - o_e[idx3(im1, jp1, k, nxf, nyf)]
                        - o_e[idx3(i, jp1, k, nxf, nyf)]
                        - u_c[idx3(i, jp1, km1, nxf, nyf)]
                        - u_c[idx3(i, jp1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_ps_v,
                    o_n[idx3(i, jm1, k, nxf, nyf)],
                    o_c[idx3(i, jm1, k, nxf, nyf)]
                        - o_e[idx3(im1, jm1, k, nxf, nyf)]
                        - o_e[idx3(i, jm1, k, nxf, nyf)]
                        - u_c[idx3(i, jm1, km1, nxf, nyf)]
                        - u_c[idx3(i, jm1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_pe_v,
                    o_e[ip],
                    o_c[idx3(ip1, j, k, nxf, nyf)]
                        - u_c[idx3(ip1, j, km1, nxf, nyf)]
                        - u_c[idx3(ip1, j, k, nxf, nyf)]
                        - o_n[idx3(ip1, j, k, nxf, nyf)]
                        - o_n[idx3(ip1, jm1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_pw_v,
                    o_e[idx3(im1, j, k, nxf, nyf)],
                    o_c[idx3(im1, j, k, nxf, nyf)]
                        - u_c[idx3(im1, j, km1, nxf, nyf)]
                        - u_c[idx3(im1, j, k, nxf, nyf)]
                        - o_n[idx3(im1, j, k, nxf, nyf)]
                        - o_n[idx3(im1, jm1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_pne_v,
                    o_n[idx3(ip1, j, k, nxf, nyf)] * o_pe_v
                        + o_e[idx3(i, jp1, k, nxf, nyf)] * o_pn_v,
                    o_c[idx3(ip1, jp1, k, nxf, nyf)]
                        - u_c[idx3(ip1, jp1, km1, nxf, nyf)]
                        - u_c[idx3(ip1, jp1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_pnw_v,
                    o_n[idx3(im1, j, k, nxf, nyf)] * o_pw_v
                        + o_e[idx3(im1, jp1, k, nxf, nyf)] * o_pn_v,
                    o_c[idx3(im1, jp1, k, nxf, nyf)]
                        - u_c[idx3(im1, jp1, km1, nxf, nyf)]
                        - u_c[idx3(im1, jp1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_pse_v,
                    o_n[idx3(ip1, jm1, k, nxf, nyf)] * o_pe_v
                        + o_e[idx3(i, jm1, k, nxf, nyf)] * o_ps_v,
                    o_c[idx3(ip1, jm1, k, nxf, nyf)]
                        - u_c[idx3(ip1, jm1, km1, nxf, nyf)]
                        - u_c[idx3(ip1, jm1, k, nxf, nyf)],
                );
                assign_if_finite(
                    &mut o_psw_v,
                    o_n[idx3(im1, jm1, k, nxf, nyf)] * o_pw_v
                        + o_e[idx3(im1, jm1, k, nxf, nyf)] * o_ps_v,
                    o_c[idx3(im1, jm1, k, nxf, nyf)]
                        - u_c[idx3(im1, jm1, km1, nxf, nyf)]
                        - u_c[idx3(im1, jm1, k, nxf, nyf)],
                );

                assign_if_finite(
                    &mut d_pc_v,
                    u_c[idx3(i, j, km1, nxf, nyf)],
                    o_c[idx3(i, j, km1, nxf, nyf)]
                        - o_n[idx3(i, j, km1, nxf, nyf)]
                        - o_n[idx3(i, jm1, km1, nxf, nyf)]
                        - o_e[idx3(im1, j, km1, nxf, nyf)]
                        - o_e[idx3(i, j, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_pn_v,
                    o_n[idx3(i, j, km1, nxf, nyf)] * d_pc_v
                        + u_c[idx3(i, jp1, km1, nxf, nyf)] * o_pn_v,
                    o_c[idx3(i, jp1, km1, nxf, nyf)]
                        - o_e[idx3(im1, jp1, km1, nxf, nyf)]
                        - o_e[idx3(i, jp1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_ps_v,
                    o_n[idx3(i, jm1, km1, nxf, nyf)] * d_pc_v
                        + u_c[idx3(i, jm1, km1, nxf, nyf)] * o_ps_v,
                    o_c[idx3(i, jm1, km1, nxf, nyf)]
                        - o_e[idx3(im1, jm1, km1, nxf, nyf)]
                        - o_e[idx3(i, jm1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_pe_v,
                    u_c[idx3(ip1, j, km1, nxf, nyf)] * o_pe_v
                        + o_e[idx3(i, j, km1, nxf, nyf)] * d_pc_v,
                    o_c[idx3(ip1, j, km1, nxf, nyf)]
                        - o_n[idx3(ip1, j, km1, nxf, nyf)]
                        - o_n[idx3(ip1, jm1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_pw_v,
                    u_c[idx3(im1, j, km1, nxf, nyf)] * o_pw_v
                        + o_e[idx3(im1, j, km1, nxf, nyf)] * d_pc_v,
                    o_c[idx3(im1, j, km1, nxf, nyf)]
                        - o_n[idx3(im1, j, km1, nxf, nyf)]
                        - o_n[idx3(im1, jm1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_pne_v,
                    u_c[idx3(ip1, jp1, km1, nxf, nyf)] * o_pne_v
                        + o_e[idx3(i, jp1, km1, nxf, nyf)] * d_pn_v
                        + o_n[idx3(ip1, j, km1, nxf, nyf)] * d_pe_v,
                    o_c[idx3(ip1, jp1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_pnw_v,
                    u_c[idx3(im1, jp1, km1, nxf, nyf)] * o_pnw_v
                        + o_e[idx3(im1, jp1, km1, nxf, nyf)] * d_pn_v
                        + o_n[idx3(im1, j, km1, nxf, nyf)] * d_pw_v,
                    o_c[idx3(im1, jp1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_pse_v,
                    u_c[idx3(ip1, jm1, km1, nxf, nyf)] * o_pse_v
                        + o_e[idx3(i, jm1, km1, nxf, nyf)] * d_ps_v
                        + o_n[idx3(ip1, jm1, km1, nxf, nyf)] * d_pe_v,
                    o_c[idx3(ip1, jm1, km1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut d_psw_v,
                    u_c[idx3(im1, jm1, km1, nxf, nyf)] * o_psw_v
                        + o_e[idx3(im1, jm1, km1, nxf, nyf)] * d_ps_v
                        + o_n[idx3(im1, jm1, km1, nxf, nyf)] * d_pw_v,
                    o_c[idx3(im1, jm1, km1, nxf, nyf)],
                );

                assign_if_finite(
                    &mut u_pc_v,
                    u_c[ip],
                    o_c[idx3(i, j, kp1, nxf, nyf)]
                        - o_n[idx3(i, j, kp1, nxf, nyf)]
                        - o_n[idx3(i, jm1, kp1, nxf, nyf)]
                        - o_e[idx3(im1, j, kp1, nxf, nyf)]
                        - o_e[idx3(i, j, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_pn_v,
                    o_n[idx3(i, j, kp1, nxf, nyf)] * u_pc_v
                        + u_c[idx3(i, jp1, k, nxf, nyf)] * o_pn_v,
                    o_c[idx3(i, jp1, kp1, nxf, nyf)]
                        - o_e[idx3(im1, jp1, kp1, nxf, nyf)]
                        - o_e[idx3(i, jp1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_ps_v,
                    o_n[idx3(i, jm1, kp1, nxf, nyf)] * u_pc_v
                        + u_c[idx3(i, jm1, k, nxf, nyf)] * o_ps_v,
                    o_c[idx3(i, jm1, kp1, nxf, nyf)]
                        - o_e[idx3(im1, jm1, kp1, nxf, nyf)]
                        - o_e[idx3(i, jm1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_pe_v,
                    u_c[idx3(ip1, j, k, nxf, nyf)] * o_pe_v
                        + o_e[idx3(i, j, kp1, nxf, nyf)] * u_pc_v,
                    o_c[idx3(ip1, j, kp1, nxf, nyf)]
                        - o_n[idx3(ip1, j, kp1, nxf, nyf)]
                        - o_n[idx3(ip1, jm1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_pw_v,
                    u_c[idx3(im1, j, k, nxf, nyf)] * o_pw_v
                        + o_e[idx3(im1, j, kp1, nxf, nyf)] * u_pc_v,
                    o_c[idx3(im1, j, kp1, nxf, nyf)]
                        - o_n[idx3(im1, j, kp1, nxf, nyf)]
                        - o_n[idx3(im1, jm1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_pne_v,
                    u_c[idx3(ip1, jp1, k, nxf, nyf)] * o_pne_v
                        + o_e[idx3(i, jp1, kp1, nxf, nyf)] * u_pn_v
                        + o_n[idx3(ip1, j, kp1, nxf, nyf)] * u_pe_v,
                    o_c[idx3(ip1, jp1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_pnw_v,
                    u_c[idx3(im1, jp1, k, nxf, nyf)] * o_pnw_v
                        + o_e[idx3(im1, jp1, kp1, nxf, nyf)] * u_pn_v
                        + o_n[idx3(im1, j, kp1, nxf, nyf)] * u_pw_v,
                    o_c[idx3(im1, jp1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_pse_v,
                    u_c[idx3(ip1, jm1, k, nxf, nyf)] * o_pse_v
                        + o_e[idx3(i, jm1, kp1, nxf, nyf)] * u_ps_v
                        + o_n[idx3(ip1, jm1, kp1, nxf, nyf)] * u_pe_v,
                    o_c[idx3(ip1, jm1, kp1, nxf, nyf)],
                );
                assign_if_finite(
                    &mut u_psw_v,
                    u_c[idx3(im1, jm1, k, nxf, nyf)] * o_psw_v
                        + o_e[idx3(im1, jm1, kp1, nxf, nyf)] * u_ps_v
                        + o_n[idx3(im1, jm1, kp1, nxf, nyf)] * u_pw_v,
                    o_c[idx3(im1, jm1, kp1, nxf, nyf)],
                );

                pc[at(0, c)] = o_pc_v;
                pc[at(1, c)] = o_pn_v;
                pc[at(2, c)] = o_ps_v;
                pc[at(3, c)] = o_pe_v;
                pc[at(4, c)] = o_pw_v;
                pc[at(5, c)] = o_pne_v;
                pc[at(6, c)] = o_pnw_v;
                pc[at(7, c)] = o_pse_v;
                pc[at(8, c)] = o_psw_v;
                pc[at(9, c)] = u_pc_v;
                pc[at(10, c)] = u_pn_v;
                pc[at(11, c)] = u_ps_v;
                pc[at(12, c)] = u_pe_v;
                pc[at(13, c)] = u_pw_v;
                pc[at(14, c)] = u_pne_v;
                pc[at(15, c)] = u_pnw_v;
                pc[at(16, c)] = u_pse_v;
                pc[at(17, c)] = u_psw_v;
                pc[at(18, c)] = d_pc_v;
                pc[at(19, c)] = d_pn_v;
                pc[at(20, c)] = d_ps_v;
                pc[at(21, c)] = d_pe_v;
                pc[at(22, c)] = d_pw_v;
                pc[at(23, c)] = d_pne_v;
                pc[at(24, c)] = d_pnw_v;
                pc[at(25, c)] = d_pse_v;
                pc[at(26, c)] = d_psw_v;
            }
        }
    }

    pc
}

/// Build a dormant trilinear coarse-to-fine prolongation stencil.
///
/// The active Rust PMGC solve path does not consume this operator. It is only a
/// simple trilinear placeholder and does not implement the APBS C
/// operator-dependent `VbuildP_op7` / `VbuildP_op27` formulas.
pub fn build_p(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    o_pc: &mut [f64], o_pn: &mut [f64], _o_ps: &mut [f64],
    o_pe: &mut [f64], _o_pw: &mut [f64],
    o_pne: &mut [f64], o_pnw: &mut [f64],
    o_pse: &mut [f64], o_psw: &mut [f64],
    u_pc: &mut [f64], u_pn: &mut [f64], _u_ps: &mut [f64],
    u_pe: &mut [f64], _u_pw: &mut [f64],
    u_pne: &mut [f64], u_pnw: &mut [f64],
    u_pse: &mut [f64], u_psw: &mut [f64],
    d_pc: &mut [f64], d_pn: &mut [f64], _d_ps: &mut [f64],
    d_pe: &mut [f64], _d_pw: &mut [f64],
    d_pne: &mut [f64], d_pnw: &mut [f64],
    d_pse: &mut [f64], d_psw: &mut [f64],
) {
    // Trilinear interpolation weights
    // For each fine grid point, compute weights to the 8 surrounding coarse grid points

    let nxny_f = nxf * nyf;
    let _nxny_c = nxc * nyc;

    for kf in 0..nzf {
        for jf in 0..nyf {
            for ipf in 0..nxf {
                let ip = ipf + jf * nxf + kf * nxny_f;

                // Compute fractional coarse grid coordinates
                let fc = ipf as f64 / 2.0;
                let fj = jf as f64 / 2.0;
                let fk = kf as f64 / 2.0;

                let ic = fc.floor() as usize;
                let jc = fj.floor() as usize;
                let kc = fk.floor() as usize;

                let dx = fc - ic as f64;
                let dy = fj - jc as f64;
                let dz = fk - kc as f64;

                // Clamp to coarse grid
                let _ic = ic.min(nxc - 1);
                let _jc = jc.min(nyc - 1);
                let _kc = kc.min(nzc - 1);

                // Trilinear weights for 8 coarse neighbors
                let w000 = (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
                let w100 = dx * (1.0 - dy) * (1.0 - dz);
                let w010 = (1.0 - dx) * dy * (1.0 - dz);
                let w110 = dx * dy * (1.0 - dz);
                let w001 = (1.0 - dx) * (1.0 - dy) * dz;
                let w101 = dx * (1.0 - dy) * dz;
                let w011 = (1.0 - dx) * dy * dz;
                let w111 = dx * dy * dz;

                // Center weights
                o_pc[ip] = w000;
                o_pe[ip] = w100;
                o_pn[ip] = w010;
                u_pc[ip] = w001;
                u_pe[ip] = w101;
                u_pn[ip] = w011;

                // Diagonal weights (corner contributions)
                o_pne[ip] = w110;
                o_pnw[ip] = 0.0;
                o_pse[ip] = 0.0;
                o_psw[ip] = 0.0;
                u_pne[ip] = w111;
                u_pnw[ip] = w011;
                u_pse[ip] = w101;
                u_psw[ip] = w001;

                // Down-plane weights
                d_pc[ip] = w000;
                d_pe[ip] = w100;
                d_pn[ip] = w010;
                d_pne[ip] = w110;
                d_pnw[ip] = 0.0;
                d_pse[ip] = 0.0;
                d_psw[ip] = 0.0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{build_p, build_p_op7_block, build_p_trilin_block};
    use crate::build_a::build_a;
    use crate::matvec::{interp_pmg, restrc, split_pc27};

    #[test]
    fn build_p_trilinear_weights_partition_unity_at_cell_center() {
        let (nxf, nyf, nzf) = (3usize, 3usize, 3usize);
        let (nxc, nyc, nzc) = (2usize, 2usize, 2usize);
        let nf = nxf * nyf * nzf;
        let mut o_pc = vec![0.0; nf];
        let mut o_pn = vec![0.0; nf];
        let mut o_ps = vec![0.0; nf];
        let mut o_pe = vec![0.0; nf];
        let mut o_pw = vec![0.0; nf];
        let mut o_pne = vec![0.0; nf];
        let mut o_pnw = vec![0.0; nf];
        let mut o_pse = vec![0.0; nf];
        let mut o_psw = vec![0.0; nf];
        let mut u_pc = vec![0.0; nf];
        let mut u_pn = vec![0.0; nf];
        let mut u_ps = vec![0.0; nf];
        let mut u_pe = vec![0.0; nf];
        let mut u_pw = vec![0.0; nf];
        let mut u_pne = vec![0.0; nf];
        let mut u_pnw = vec![0.0; nf];
        let mut u_pse = vec![0.0; nf];
        let mut u_psw = vec![0.0; nf];
        let mut d_pc = vec![0.0; nf];
        let mut d_pn = vec![0.0; nf];
        let mut d_ps = vec![0.0; nf];
        let mut d_pe = vec![0.0; nf];
        let mut d_pw = vec![0.0; nf];
        let mut d_pne = vec![0.0; nf];
        let mut d_pnw = vec![0.0; nf];
        let mut d_pse = vec![0.0; nf];
        let mut d_psw = vec![0.0; nf];

        build_p(
            nxf, nyf, nzf, nxc, nyc, nzc,
            &mut o_pc, &mut o_pn, &mut o_ps, &mut o_pe, &mut o_pw,
            &mut o_pne, &mut o_pnw, &mut o_pse, &mut o_psw,
            &mut u_pc, &mut u_pn, &mut u_ps, &mut u_pe, &mut u_pw,
            &mut u_pne, &mut u_pnw, &mut u_pse, &mut u_psw,
            &mut d_pc, &mut d_pn, &mut d_ps, &mut d_pe, &mut d_pw,
            &mut d_pne, &mut d_pnw, &mut d_pse, &mut d_psw,
        );

        let center = 1 + 1 * nxf + 1 * nxf * nyf;
        let sum = o_pc[center]
            + o_pe[center]
            + o_pn[center]
            + o_pne[center]
            + u_pc[center]
            + u_pe[center]
            + u_pn[center]
            + u_pne[center];
        assert!((sum - 1.0).abs() < 1.0e-12);
        assert!((o_pc[center] - 0.125).abs() < 1.0e-12);
    }

    #[test]
    fn build_p_trilin_block_matches_apbs_c_constant_weights() {
        let pc = build_p_trilin_block(2, 2, 2);
        let parts = split_pc27(&pc, 2, 2, 2).expect("pc block");
        assert!(parts.o_pc.iter().all(|v| (*v - 1.0).abs() < 1.0e-12));
        assert!(parts.o_pe.iter().all(|v| (*v - 0.5).abs() < 1.0e-12));
        assert!(parts.o_pne.iter().all(|v| (*v - 0.25).abs() < 1.0e-12));
        assert!(parts.u_pc.iter().all(|v| (*v - 0.5).abs() < 1.0e-12));
        assert!(parts.d_psw.iter().all(|v| (*v - 0.125).abs() < 1.0e-12));
    }

    #[test]
    fn build_p_op7_matches_trilinear_on_uniform_stencil_interior() {
        let (nxf, nyf, nzf) = (5usize, 5usize, 5usize);
        let nf = nxf * nyf * nzf;
        let mut ac = vec![0.0; 4 * nf];
        ac[0..nf].fill(6.0);
        ac[nf..2 * nf].fill(1.0);
        ac[2 * nf..3 * nf].fill(1.0);
        ac[3 * nf..4 * nf].fill(1.0);

        let pc = build_p_op7_block(nxf, nyf, nzf, &ac);
        let parts = split_pc27(&pc, 3, 3, 3).expect("pc block");
        let c = 1 + 3 + 9;

        assert!((parts.o_pc[c] - 1.0).abs() < 1.0e-12);
        assert!((parts.o_pn[c] - 0.5).abs() < 1.0e-12);
        assert!((parts.o_ps[c] - 0.5).abs() < 1.0e-12);
        assert!((parts.o_pe[c] - 0.5).abs() < 1.0e-12);
        assert!((parts.o_pw[c] - 0.5).abs() < 1.0e-12);
        assert!((parts.o_pne[c] - 0.25).abs() < 1.0e-12);
        assert!((parts.o_pnw[c] - 0.25).abs() < 1.0e-12);
        assert!((parts.o_pse[c] - 0.25).abs() < 1.0e-12);
        assert!((parts.o_psw[c] - 0.25).abs() < 1.0e-12);
        assert!((parts.u_pc[c] - 0.5).abs() < 1.0e-12);
        assert!((parts.d_pc[c] - 0.5).abs() < 1.0e-12);
        assert!((parts.u_pne[c] - 0.125).abs() < 1.0e-12);
        assert!((parts.d_psw[c] - 0.125).abs() < 1.0e-12);
    }

    #[test]
    fn build_p_op7_interpolation_preserves_constant_fields_on_uniform_stencil() {
        let (nxf, nyf, nzf) = (5usize, 5usize, 5usize);
        let nf = nxf * nyf * nzf;
        let mut ac = vec![0.0; 4 * nf];
        ac[0..nf].fill(6.0);
        ac[nf..2 * nf].fill(1.0);
        ac[2 * nf..3 * nf].fill(1.0);
        ac[3 * nf..4 * nf].fill(1.0);
        let pc = build_p_op7_block(nxf, nyf, nzf, &ac);

        let coarse = vec![2.5; 3 * 3 * 3];
        let mut fine = vec![0.0; nf];
        interp_pmg(3, 3, 3, nxf, nyf, nzf, &coarse, &mut fine, &pc);
        for k in 0..nzf {
            for j in 0..nyf {
                for i in 0..nxf {
                    let v = fine[i + j * nxf + k * nxf * nyf];
                    assert!(v.is_finite());
                    if i == 0 || i + 1 == nxf || j == 0 || j + 1 == nyf || k == 0 || k + 1 == nzf {
                        assert!(v.abs() < 1.0e-12);
                    } else {
                        assert!((v - 2.5).abs() < 1.0e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn build_p_op7_restriction_stays_finite_on_uniform_stencil() {
        let (nxf, nyf, nzf) = (5usize, 5usize, 5usize);
        let nf = nxf * nyf * nzf;
        let mut ac = vec![0.0; 4 * nf];
        ac[0..nf].fill(6.0);
        ac[nf..2 * nf].fill(1.0);
        ac[2 * nf..3 * nf].fill(1.0);
        ac[3 * nf..4 * nf].fill(1.0);
        let pc = build_p_op7_block(nxf, nyf, nzf, &ac);
        let fine_const = vec![1.75; nf];
        let mut coarse_out = vec![0.0; 3 * 3 * 3];
        restrc(nxf, nyf, nzf, 3, 3, 3, &fine_const, &mut coarse_out, &pc);
        assert!(coarse_out.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn build_p_op7_from_built_operator_stays_finite() {
        let (nx, ny, nz) = (5usize, 5usize, 5usize);
        let nf = nx * ny * nz;
        let mut ac = vec![0.0; 4 * nf];
        let mut cc = vec![0.0; nf];
        let mut fc = vec![0.0; nf];
        let xf = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let yf = xf.clone();
        let zf = xf.clone();
        let gxcf = vec![0.0; 2 * ny * nz];
        let gycf = vec![0.0; 2 * nx * nz];
        let gzcf = vec![0.0; 2 * nx * ny];
        let mut a1cf = vec![1.0; nf];
        let mut a2cf = vec![1.0; nf];
        let mut a3cf = vec![1.0; nf];
        let ccf = vec![0.0; nf];
        let fcf = vec![0.0; nf];

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let ip = i + j * nx + k * nx * ny;
                    a1cf[ip] = 1.0 + 0.05 * i as f64;
                    a2cf[ip] = 1.0 + 0.05 * j as f64;
                    a3cf[ip] = 1.0 + 0.05 * k as f64;
                }
            }
        }

        build_a(
            nx, ny, nz, 0, 0, 4,
            &mut ac, &mut cc, &mut fc,
            &xf, &yf, &zf,
            &gxcf, &gycf, &gzcf,
            &a1cf, &a2cf, &a3cf,
            &ccf, &fcf,
        );

        let pc = build_p_op7_block(nx, ny, nz, &ac);
        assert!(pc.iter().all(|v| v.is_finite()));
        let max_abs = pc.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(max_abs < 10.0, "unexpected large op7 weight: {max_abs}");
    }
}
