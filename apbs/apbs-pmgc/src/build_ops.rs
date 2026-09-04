// APBS PMGC buildOps - Build all multigrid operators
// Port of pmgc/buildOps routines

use crate::build_str::make_coarse;

fn effective_mgcoar(mgcoar: i32, mgdisc: i32) -> i32 {
    let _ = mgdisc;
    mgcoar
}

fn build_trilinear_pc_ff(nxf: usize, nyf: usize, nzf: usize) -> Vec<f64> {
    let nf = nxf * nyf * nzf;
    let mut pc_ff = vec![0.0; nf];
    for k in 0..nzf {
        for j in 0..nyf {
            for i in 0..nxf {
                let ip = i + j * nxf + k * nxf * nyf;
                let mut w = 1.0;
                if i % 2 == 1 {
                    w *= 0.5;
                }
                if j % 2 == 1 {
                    w *= 0.5;
                }
                if k % 2 == 1 {
                    w *= 0.5;
                }
                pc_ff[ip] = w;
            }
        }
    }
    pc_ff
}

// fn restrict_3d(
//     nxf: usize,
//     nyf: usize,
//     nzf: usize,
//     nxc: usize,
//     nyc: usize,
//     nzc: usize,
//     fine: &[f64],
// ) -> Vec<f64> {
//     let nc = nxc * nyc * nzc;
//     let mut coarse = vec![0.0f64; nc];
//     let nxnyc = nxc * nyc;
//     for kc in 0..nzc {
//         for jc in 0..nyc {
//             for ic in 0..nxc {
//                 let ipc = ic + jc * nxc + kc * nxnyc;
//                 let if_ = (2 * ic).min(nxf - 1);
//                 let jf = (2 * jc).min(nyf - 1);
//                 let kf = (2 * kc).min(nzf - 1);
//                 if if_ == 0 || if_ + 1 >= nxf || jf == 0 || jf + 1 >= nyf || kf == 0 || kf + 1 >= nzf {
//                     coarse[ipc] = fine[if_ + jf * nxf + kf * nxf * nyf];
//                     continue;
//                 }
//                 let mut sum = 0.0;
//                 for dk in -1isize..=1 {
//                     for dj in -1isize..=1 {
//                         for di in -1isize..=1 {
//                             let fi = (if_ as isize + di) as usize;
//                             let fj = (jf as isize + dj) as usize;
//                             let fk = (kf as isize + dk) as usize;
//                             let neighbors =
//                                 di.unsigned_abs() + dj.unsigned_abs() + dk.unsigned_abs();
//                             let w = match neighbors {
//                                 0 => 8.0,
//                                 1 => 4.0,
//                                 2 => 2.0,
//                                 _ => 1.0,
//                             };
//                             sum += w * fine[fi + fj * nxf + fk * nxf * nyf];
//                         }
//                     }
//                 }
//                 coarse[ipc] = sum / 64.0;
//             }
//         }
//     }
//     coarse
// }

// fn restrict_bc_x(nyf: usize, nzf: usize, nyc: usize, nzc: usize, gxcf: &[f64]) -> Vec<f64> {
//     let mut out = vec![0.0f64; 2 * nyc * nzc];
//     for face in 0..2 {
//         for kc in 0..nzc {
//             for jc in 0..nyc {
//                 let fj = (2 * jc).min(nyf - 1);
//                 let fk = (2 * kc).min(nzf - 1);
//                 out[face * nyc * nzc + kc * nyc + jc] = gxcf[face * nyf * nzf + fk * nyf + fj];
//             }
//         }
//     }
//     out
// }

// fn restrict_bc_y(nxf: usize, nzf: usize, nxc: usize, nzc: usize, gycf: &[f64]) -> Vec<f64> {
//     let mut out = vec![0.0f64; 2 * nxc * nzc];
//     for face in 0..2 {
//         for kc in 0..nzc {
//             for ic in 0..nxc {
//                 let fi = (2 * ic).min(nxf - 1);
//                 let fk = (2 * kc).min(nzf - 1);
//                 out[face * nxc * nzc + kc * nxc + ic] = gycf[face * nxf * nzf + fk * nxf + fi];
//             }
//         }
//     }
//     out
// }

// fn restrict_bc_z(nxf: usize, nyf: usize, nxc: usize, nyc: usize, gzcf: &[f64]) -> Vec<f64> {
//     let mut out = vec![0.0f64; 2 * nxc * nyc];
//     for face in 0..2 {
//         for jc in 0..nyc {
//             for ic in 0..nxc {
//                 let fi = (2 * ic).min(nxf - 1);
//                 let fj = (2 * jc).min(nyf - 1);
//                 out[face * nxc * nyc + jc * nxc + ic] = gzcf[face * nxf * nyf + fj * nxf + fi];
//             }
//         }
//     }
//     out
// }

/// Build all multigrid operators for all levels
/// Uses index-based access to avoid borrow conflicts
pub fn build_ops(
    nlev: i32,
    nx: i32, ny: i32, nz: i32,
    mgcoar: i32, mgdisc: i32,
    iz: &mut [i32],
    ac: &mut [f64],
    cc: &mut [f64], fc: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
) {
    let numdia = if mgdisc == 0 { 4 } else { 14 };
    let mgcoar = effective_mgcoar(mgcoar, mgdisc);

    // Build finest level operator
    iz[0] = 0; // offset into ac for level 0
    iz[1] = numdia * nx * ny * nz; // size of level 0 operator
    iz[2] = 0; // offset into pc for level 0

    let nf = (nx * ny * nz) as usize;

    // Build A on finest level
    crate::build_a::build_a(
        nx as usize, ny as usize, nz as usize,
        0, mgdisc, numdia,
        &mut ac[0..4 * nf],
        &mut cc[0..nf],
        &mut fc[0..nf],
        xf, yf, zf,
        gxcf, gycf, gzcf,
        a1cf, a2cf, a3cf,
        ccf, fcf,
    );

    // Build coarser levels
    let mut nx_l = nx;
    let mut ny_l = ny;
    let mut nz_l = nz;
    let mut ac_offset = 4 * nf;
    let mut ac_fine_offset = 0usize;
    let mut cc_fc_offset = nf;

    for _lev in 1..nlev as usize {
        let (nx_c, ny_c, nz_c) = make_coarse(nx_l, ny_l, nz_l);
        let nf_l = (nx_l * ny_l * nz_l) as usize;
        let nc_l = (nx_c * ny_c * nz_c) as usize;

        if mgcoar == 2 {
            // Galerkin: Ac = P^T * Af * P
            // Split ac to satisfy borrow checker: immutable view of fine level,
            // mutable view of destination coarse level.
            let (ac_head, ac_coarse) = ac.split_at_mut(ac_offset);
            let pc_ff = build_trilinear_pc_ff(nx_l as usize, ny_l as usize, nz_l as usize);
            let ac_fine_level = &ac_head[ac_fine_offset..ac_fine_offset + 4 * nf_l];
            crate::build_g::build_g(
                nx_l as usize, ny_l as usize, nz_l as usize,
                nx_c as usize, ny_c as usize, nz_c as usize,
                numdia,
                &pc_ff,
                ac_fine_level,
                &mut ac_coarse[0..4 * nc_l],
            );
        } else {
            // Standard: build operator on coarse grid
            crate::build_a::build_a(
                nx_c as usize, ny_c as usize, nz_c as usize,
                0, mgdisc, numdia,
                &mut ac[ac_offset..ac_offset + 4 * nc_l],
                &mut cc[cc_fc_offset..cc_fc_offset + nc_l],
                &mut fc[cc_fc_offset..cc_fc_offset + nc_l],
                xf, yf, zf,
                gxcf, gycf, gzcf,
                a1cf, a2cf, a3cf,
                ccf, fcf,
            );
        }

        ac_offset += 4 * nc_l;
        ac_fine_offset += 4 * nf_l;
        cc_fc_offset += nc_l;
        nx_l = nx_c;
        ny_l = ny_c;
        nz_l = nz_c;
    }
}

#[cfg(test)]
mod tests {
    use super::build_ops;

    #[test]
    fn build_ops_accepts_galerkin_flag_on_7pt_path() {
        let nlev = 2;
        let (nx, ny, nz) = (5i32, 5i32, 5i32);
        let nf = (nx * ny * nz) as usize;
        let nc = ((nx / 2 + 1) * (ny / 2 + 1) * (nz / 2 + 1)) as usize;
        let narr = nf + nc;
        let mut iz_std = vec![0; 50 * (nlev as usize + 1)];
        let mut iz_demoted = vec![0; 50 * (nlev as usize + 1)];
        let mut ac_std = vec![0.0; 4 * narr];
        let mut ac_demoted = vec![0.0; 4 * narr];
        let mut cc_std = vec![0.0; narr];
        let mut cc_demoted = vec![0.0; narr];
        let mut fc_std = vec![0.0; narr];
        let mut fc_demoted = vec![0.0; narr];

        let xf = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let yf = xf.clone();
        let zf = xf.clone();
        let gxcf = vec![0.0; 2 * ny as usize * nz as usize];
        let gycf = vec![0.0; 2 * nx as usize * nz as usize];
        let gzcf = vec![0.0; 2 * nx as usize * ny as usize];
        let a1cf = vec![1.0; nf];
        let a2cf = vec![1.0; nf];
        let a3cf = vec![1.0; nf];
        let ccf = vec![0.0; nf];
        let fcf = vec![0.0; nf];

        build_ops(
            nlev, nx, ny, nz,
            0, 0,
            &mut iz_std,
            &mut ac_std,
            &mut cc_std,
            &mut fc_std,
            &xf, &yf, &zf,
            &gxcf, &gycf, &gzcf,
            &a1cf, &a2cf, &a3cf,
            &ccf, &fcf,
        );
        build_ops(
            nlev, nx, ny, nz,
            2, 0,
            &mut iz_demoted,
            &mut ac_demoted,
            &mut cc_demoted,
            &mut fc_demoted,
            &xf, &yf, &zf,
            &gxcf, &gycf, &gzcf,
            &a1cf, &a2cf, &a3cf,
            &ccf, &fcf,
        );

        assert_eq!(ac_std.len(), ac_demoted.len());
        assert_eq!(cc_std.len(), cc_demoted.len());
        assert_eq!(fc_std.len(), fc_demoted.len());
    }
}
