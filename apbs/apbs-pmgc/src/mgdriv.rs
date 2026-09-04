// APBS PMGC mgdriv - Top-level multigrid driver
// Port of pmgc/mgdrvd.c

use std::sync::OnceLock;

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

fn nonlinear_solver_mode() -> &'static str {
    static MODE: OnceLock<String> = OnceLock::new();
    MODE.get_or_init(|| {
        std::env::var("APBS_RUST_NONLIN_SOLVER")
            .unwrap_or_else(|_| "newton".to_string())
            .to_lowercase()
    })
}

fn pc_mode() -> String {
    std::env::var("APBS_RUST_PC_MODE").unwrap_or_else(|_| "op7".to_string())
}

fn force_mgsolv() -> Option<i32> {
    std::env::var("APBS_RUST_FORCE_MGSOLV")
        .ok()
        .and_then(|s| s.parse::<i32>().ok())
}

fn override_itmax(default: usize) -> usize {
    std::env::var("APBS_RUST_OVERRIDE_ITMAX")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn override_errtol(default: f64) -> f64 {
    std::env::var("APBS_RUST_OVERRIDE_ERRTOL")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn restrict_mode() -> String {
    std::env::var("APBS_RUST_RESTRICT_MODE").unwrap_or_else(|_| "inject".to_string())
}

/// Top-level multigrid solver driver
pub fn mgdriv(
    iparm: &mut [i32], rparm: &mut [f64],
    _iwork: &mut [i32], _rwork: &mut [f64],
    u: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64], tcf: &mut [f64],
) {
    // Extract parameters from iparm
    let nx = iparm[0] as usize;
    let ny = iparm[1] as usize;
    let nz = iparm[2] as usize;
    let nlev = iparm[3] as usize;
    let _mgkey = iparm[4];
    let nonlin = iparm[5];
    let _mgcoar = iparm[9];
    let mgdisc = iparm[10];
    let mgsolv = force_mgsolv().unwrap_or(iparm[11]);
    let nu1 = iparm[12];
    let nu2 = iparm[13];

    // Extract real parameters
    let omegal = rparm[0];
    let omegan = rparm[1];
    let errtol = override_errtol(rparm[2]);
    let itmax = override_itmax(iparm[7] as usize);

    let numdia = if mgdisc == 0 { 4 } else { 14 };
    let nf = nx * ny * nz;
    let mut total_op_size = 0usize;
    let mut level_sizes = Vec::new();
    let mut level_nx = nx;
    let mut level_ny = ny;
    let mut level_nz = nz;
    for _ in 0..nlev {
        let nf_l = level_nx * level_ny * level_nz;
        level_sizes.push((level_nx, level_ny, level_nz, nf_l));
        total_op_size += 4 * nf_l;
        level_nx = level_nx / 2 + 1;
        level_ny = level_ny / 2 + 1;
        level_nz = level_nz / 2 + 1;
    }
    let mut ac = vec![0.0f64; total_op_size];
    let mut narr_total = 0usize;
    for &(_, _, _, nf_l) in &level_sizes {
        narr_total += nf_l;
    }
    let mut cc_all = vec![0.0f64; narr_total];
    let mut fc_all = vec![0.0f64; narr_total];
    let mut ac_offset = 0usize;
    let mut cc_offset = 0usize;
    let mut fc_offset = 0usize;

    crate::build_a::build_a(
        nx, ny, nz, 0, mgdisc, numdia,
        &mut ac[ac_offset..ac_offset + 4 * nf],
        &mut cc_all[cc_offset..cc_offset + nf],
        &mut fc_all[fc_offset..fc_offset + nf],
        xf, yf, zf, gxcf, gycf, gzcf,
        a1cf, a2cf, a3cf, ccf, fcf,
    );
    ac_offset += 4 * nf;
    cc_offset += nf;
    fc_offset += nf;

    let mut cur_xf = xf.to_vec();
    let mut cur_yf = yf.to_vec();
    let mut cur_zf = zf.to_vec();
    let mut cur_a1 = a1cf.to_vec();
    let mut cur_a2 = a2cf.to_vec();
    let mut cur_a3 = a3cf.to_vec();
    let mut cur_cc = ccf.to_vec();
    let mut cur_fc = fcf.to_vec();
    let mut cur_gxcf = gxcf.to_vec();
    let mut cur_gycf = gycf.to_vec();
    let mut cur_gzcf = gzcf.to_vec();
    let mut cur_nx = nx;
    let mut cur_ny = ny;
    let mut cur_nz = nz;

    for _lev in 1..nlev {
        let (nx_c, ny_c, nz_c) = crate::build_str::make_coarse(
            cur_nx as i32, cur_ny as i32, cur_nz as i32,
        );
        let nx_c = nx_c as usize;
        let ny_c = ny_c as usize;
        let nz_c = nz_c as usize;
        let nf_c = nx_c * ny_c * nz_c;

        let mut xf_c = vec![0.0f64; nx_c];
        let mut yf_c = vec![0.0f64; ny_c];
        let mut zf_c = vec![0.0f64; nz_c];
        for i in 0..nx_c {
            xf_c[i] = cur_xf[(2 * i).min(cur_nx - 1)];
        }
        for j in 0..ny_c {
            yf_c[j] = cur_yf[(2 * j).min(cur_ny - 1)];
        }
        for k in 0..nz_c {
            zf_c[k] = cur_zf[(2 * k).min(cur_nz - 1)];
        }

        let a1_c = restrict_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_a1);
        let a2_c = restrict_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_a2);
        let a3_c = restrict_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_a3);
        let cc_c = restrict_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_cc);
        let fc_c = restrict_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_fc);
        let gxcf_c = restrict_bc_x(cur_ny, cur_nz, ny_c, nz_c, &cur_xf, &xf_c, &cur_gxcf);
        let gycf_c = restrict_bc_y(cur_nx, cur_nz, nx_c, nz_c, &cur_yf, &yf_c, &cur_gycf);
        let gzcf_c = restrict_bc_z(cur_nx, cur_ny, nx_c, ny_c, &cur_zf, &zf_c, &cur_gzcf);

        crate::build_a::build_a(
            nx_c, ny_c, nz_c, 0, mgdisc, numdia,
            &mut ac[ac_offset..ac_offset + 4 * nf_c],
            &mut cc_all[cc_offset..cc_offset + nf_c],
            &mut fc_all[fc_offset..fc_offset + nf_c],
            &xf_c, &yf_c, &zf_c,
            &gxcf_c, &gycf_c, &gzcf_c,
            &a1_c, &a2_c, &a3_c, &cc_c, &fc_c,
        );

        ac_offset += 4 * nf_c;
        cc_offset += nf_c;
        fc_offset += nf_c;
        cur_xf = xf_c;
        cur_yf = yf_c;
        cur_zf = zf_c;
        cur_a1 = a1_c;
        cur_a2 = a2_c;
        cur_a3 = a3_c;
        cur_cc = cc_c;
        cur_fc = fc_c;
        cur_gxcf = gxcf_c;
        cur_gycf = gycf_c;
        cur_gzcf = gzcf_c;
        cur_nx = nx_c;
        cur_ny = ny_c;
        cur_nz = nz_c;
    }

    // Solve
    let irite = iparm[16];
    let mut w1 = vec![0.0f64; nf];
    let mut w2 = vec![0.0f64; nf];
    let mut r = vec![0.0f64; nf];
    let iz = vec![0i32; 50 * (nlev + 1)];
    let mut ipc = vec![0i32; 20 * (nlev + 1)];
    let mut rpc = vec![0.0f64; 20 * (nlev + 1)];
    let mut pc = vec![0.0f64; 27 * narr_total]; // prolongation storage

    // Initialize ipc[0] = numdia for each level
    for lev in 0..nlev {
        ipc[lev * 20] = numdia;
    }

    // Initialize rpc with per-level grid spacings.
    // For geometric 2:1 coarsening, each coarser spacing doubles.
    {
        let hx0 = if nx > 1 { xf[1] - xf[0] } else { 1.0 };
        let hy0 = if ny > 1 { yf[1] - yf[0] } else { 1.0 };
        let hz0 = if nz > 1 { zf[1] - zf[0] } else { 1.0 };
        for lev in 0..nlev {
            let scale = (1usize << lev) as f64;
            rpc[lev * 20 + 0] = hx0 * scale;
            rpc[lev * 20 + 1] = hy0 * scale;
            rpc[lev * 20 + 2] = hz0 * scale;
            rpc[lev * 20 + 3] = 0.0; // zkappa2 (for nonlinear)
        }
    }

    // Build prolongation blocks for each fine->coarse mapping.
    // Default to operator-dependent op7 prolongation, with env fallback:
    // APBS_RUST_PC_MODE=trilin
    {
        let mut pc_off = 0usize;
        let mut lx = nx as usize;
        let mut ly = ny as usize;
        let mut lz = nz as usize;
        let mut ac_fine_off = 0usize;
        let mode = pc_mode();
        for _lev in 0..nlev.saturating_sub(1) {
            let nx_c = lx / 2 + 1;
            let ny_c = ly / 2 + 1;
            let nz_c = lz / 2 + 1;
            let nf_l = lx * ly * lz;
            let block = if mode.eq_ignore_ascii_case("trilin") {
                crate::build_p::build_p_trilin_block(nx_c, ny_c, nz_c)
            } else {
                let ac_slice = &ac[ac_fine_off..ac_fine_off + 4 * nf_l];
                crate::build_p::build_p_op7_block(lx, ly, lz, ac_slice)
            };
            let block_len = block.len();
            if pc_off + block_len <= pc.len() {
                pc[pc_off..pc_off + block_len].copy_from_slice(&block);
            }
            pc_off += block_len;
            ac_fine_off += 4 * nf_l;
            lx = nx_c;
            ly = ny_c;
            lz = nz_c;
        }
    }

    // Outer iteration loop: run V-cycles until convergence
    for outer in 0..itmax {
        if nonlin == 0 {
            // Linear solver (V-cycle)
            crate::mgcs::mgcs(
                nlev as i32,
                nx as i32, ny as i32, nz as i32,
                &ipc, &rpc, &ac, &cc_all, &fc_all,
                0, 0, 0,  // ac_off, cc_off, fc_off for level 0
                &pc, &iz,
                u, &mut w1, &mut w2, &mut r,
                nu1, nu2, omegal, irite, mgsolv,
            );
        } else {
            // Nonlinear solver: default to the stable Newton path.
            // Set APBS_RUST_NONLIN_SOLVER=mgfas to exercise the experimental FAS path.
            if nonlinear_solver_mode() == "mgfas" {
                crate::mgfas::mgfas(
                    nlev as i32,
                    nx as i32, ny as i32, nz as i32,
                    &ipc, &rpc, &ac, &cc_all, &fc_all, &pc, &iz,
                    u, &mut w1, &mut w2, &mut r,
                    nu1, nu2, omegan, irite,
                );
            } else {
                crate::newton::newton(
                    nx, ny, nz,
                    &ipc, &rpc,
                    &ac, &cc_all, &fc_all,
                    u,
                    &mut w1, &mut w2, &mut r,
                    itmax as i32, errtol,
                    nlev as i32,
                    &pc, &iz,
                    nu1, nu2,
                    omegal, irite, mgsolv,
                );
                break;
            }
        }

        // Check convergence using the same residual definition as the active solver.
        if nonlin == 0 {
            crate::blas::mresid(
                nx, ny, nz,
                &ipc, &rpc,
                &ac[0..nf], &ac[nf..2*nf], &ac[2*nf..3*nf], &ac[3*nf..4*nf],
                &cc_all[0..nf],
                u, &fc_all[0..nf], &mut r,
            );
        } else {
            crate::mgfas::compute_nonlinear_residual(
                nx, ny, nz,
                &ipc, &rpc,
                &ac[0..4 * nf],
                &cc_all[0..nf],
                &fc_all[0..nf],
                u,
                &mut r,
            );
        }
        let rnorm = crate::blas::xnrm2(nf, &r, 0);
        let fnorm = crate::blas::xnrm2(nf, &fc_all[0..nf], 0);
        let u_norm = crate::blas::xnrm2(nf, u, 0);
        let rel_err = if fnorm > 0.0 { rnorm / fnorm } else { rnorm };

        if debug_enabled() && (outer < 3 || outer % 10 == 0 || rel_err < errtol) {
            eprintln!("[DEBUG-MGDRV] outer={}, rel_err={:.4e}, rnorm={:.4e}, fnorm={:.4e}, u_norm={:.4e}", outer, rel_err, rnorm, fnorm, u_norm);
        }

        if rel_err < errtol {
            break;
        }
    }

    // Debug: check final solution
    let u_max = u.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let u_nonzero = u.iter().filter(|x| x.abs() > 1.0e-30).count();
    if debug_enabled() {
        eprintln!("[DEBUG-MGDRV-FINAL] u_max={:.4e}, u_nonzero={}/{}, u[0]={:.4e}, u[mid]={:.4e}", u_max, u_nonzero, nf, u[0], u[nf/2]);
    }

    // Copy solution to true solution array
    tcf[..nf].copy_from_slice(&u[..nf]);
}

/// Restrict a 3D field from fine to coarse grid using coincident-point injection.
/// This better matches PMGC's coordinate-consistent coarsening than box averaging.
fn restrict_3d(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    fine: &[f64],
) -> Vec<f64> {
    let use_fullweight = !restrict_mode().eq_ignore_ascii_case("inject");
    let nc = nxc * nyc * nzc;
    let mut coarse = vec![0.0f64; nc];
    let nxnyc = nxc * nyc;

    for kc in 0..nzc {
        for jc in 0..nyc {
            for ic in 0..nxc {
                let ipc = ic + jc * nxc + kc * nxnyc;
                let if_ = (2 * ic).min(nxf - 1);
                let jf = (2 * jc).min(nyf - 1);
                let kf = (2 * kc).min(nzf - 1);

                // Boundary points use direct injection; interior uses full-weighting.
                if !use_fullweight
                    || if_ == 0 || if_ + 1 >= nxf
                    || jf == 0 || jf + 1 >= nyf
                    || kf == 0 || kf + 1 >= nzf
                {
                    coarse[ipc] = fine[if_ + jf * nxf + kf * nxf * nyf];
                    continue;
                }

                // Standard 3D full-weighting:
                // center: 8/64, faces: 4/64, edges: 2/64, corners: 1/64.
                let mut sum = 0.0f64;
                for dk in -1isize..=1 {
                    for dj in -1isize..=1 {
                        for di in -1isize..=1 {
                            let fi = (if_ as isize + di) as usize;
                            let fj = (jf as isize + dj) as usize;
                            let fk = (kf as isize + dk) as usize;
                            let neighbors = di.unsigned_abs() + dj.unsigned_abs() + dk.unsigned_abs();
                            let w = match neighbors {
                                0 => 8.0,
                                1 => 4.0,
                                2 => 2.0,
                                _ => 1.0,
                            };
                            sum += w * fine[fi + fj * nxf + fk * nxf * nyf];
                        }
                    }
                }
                coarse[ipc] = sum / 64.0;
            }
        }
    }
    coarse
}

/// Restrict x-face BC array
fn restrict_bc_x(
    nyf: usize, nzf: usize,
    nyc: usize, nzc: usize,
    _xf_f: &[f64], _xf_c: &[f64],
    gxcf: &[f64],
) -> Vec<f64> {
    // BC array has shape [2, nzf, nyf] (face=0 and face=1)
    // Restrict to [2, nzc, nyc]
    let mut gxcf_c = vec![0.0f64; 2 * nzc * nyc];
    for face in 0..2 {
        for kc in 0..nzc {
            for jc in 0..nyc {
                let fj = (2 * jc).min(nyf - 1);
                let fk = (2 * kc).min(nzf - 1);
                gxcf_c[face * nyc * nzc + kc * nyc + jc] =
                    gxcf[face * nyf * nzf + fk * nyf + fj];
            }
        }
    }
    gxcf_c
}

/// Restrict y-face BC array
fn restrict_bc_y(
    nxf: usize, nzf: usize,
    nxc: usize, nzc: usize,
    _yf_f: &[f64], _yf_c: &[f64],
    gycf: &[f64],
) -> Vec<f64> {
    let mut gycf_c = vec![0.0f64; 2 * nxc * nzc];
    for face in 0..2 {
        for kc in 0..nzc {
            for ic in 0..nxc {
                let ffi = (2 * ic).min(nxf - 1);
                let ffk = (2 * kc).min(nzf - 1);
                gycf_c[face * nxc * nzc + kc * nxc + ic] =
                    gycf[face * nxf * nzf + ffk * nxf + ffi];
            }
        }
    }
    gycf_c
}

/// Restrict z-face BC array
fn restrict_bc_z(
    nxf: usize, nyf: usize,
    nxc: usize, nyc: usize,
    _zf_f: &[f64], _zf_c: &[f64],
    gzcf: &[f64],
) -> Vec<f64> {
    let mut gzcf_c = vec![0.0f64; 2 * nxc * nyc];
    for face in 0..2 {
        for jc in 0..nyc {
            for ic in 0..nxc {
                let ffi = (2 * ic).min(nxf - 1);
                let ffj = (2 * jc).min(nyf - 1);
                gzcf_c[face * nxc * nyc + jc * nxc + ic] =
                    gzcf[face * nxf * nyf + ffj * nxf + ffi];
            }
        }
    }
    gzcf_c
}
