// APBS PMGC mgdriv - Top-level multigrid driver
// Port of pmgc/mgdrvd.c (Vmgdriv).

use std::sync::OnceLock;
//
// Linear (LPBE) path: faithful to APBS's default configuration
// (mgcoar = 2 Galerkin coarsening, mgprol = 0 trilinear interpolation,
// mgsolv = 1 direct coarsest solve, istop = 1 relative L1 residual test,
// itmax = 200, errtol = 1e-6): the fine 7-point operator is coarsened by
// Galerkin product with the trilinear pc, cc/fc are restricted with the
// same pc, and the V-cycle is Vmvcs (see mgcs::mvcs).
//
// Nonlinear path (NPBE): Newton/FAS branch, unchanged.

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

/// Top-level multigrid solver driver
#[allow(clippy::too_many_arguments)]
pub fn mgdriv(
    iparm: &mut [i32], rparm: &mut [f64],
    _iwork: &mut [i32], _rwork: &mut [f64],
    u: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64], tcf: &mut [f64],
) {
    let nx = iparm[0] as usize;
    let ny = iparm[1] as usize;
    let nz = iparm[2] as usize;
    let nonlin = iparm[5];
    let mgdisc = iparm[10];

    if nonlin == 0 {
        linear_mgdriv(
            iparm, rparm, u,
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf, ccf, fcf,
        );
    } else {
        nonlinear_mgdriv(
            iparm, rparm, u,
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf, ccf, fcf,
        );
    }

    let _ = mgdisc;
    let nf = nx * ny * nz;
    // Copy solution to true solution array
    tcf[..nf].copy_from_slice(&u[..nf]);
}

/// Linear LPBE path: Galerkin hierarchy + Vmvcs.
#[allow(clippy::too_many_arguments)]
fn linear_mgdriv(
    iparm: &[i32], rparm: &[f64],
    u: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
) {
    let nx = iparm[0] as usize;
    let ny = iparm[1] as usize;
    let nz = iparm[2] as usize;
    let nlev = iparm[3] as usize;
    let mgdisc = iparm[10];
    let nu1 = iparm[12];
    let nu2 = iparm[13];
    let itmax = iparm[7];
    let errtol = rparm[2];
    // Vpmgp_ctor2: LPBE uses mgsolv = 1 (direct banded coarsest solve).
    let mgsolv = if iparm[11] == 0 { 1 } else { iparm[11] };
    let numdia = if mgdisc == 0 { 4 } else { 14 };

    let nf = nx * ny * nz;

    // Level grid sizes (top down).
    let mut sizes = Vec::with_capacity(nlev);
    {
        let (mut cx, mut cy, mut cz) = (nx, ny, nz);
        for _ in 0..nlev {
            sizes.push((cx, cy, cz));
            let (a, b, c) = crate::build_str::make_coarse(cx as i32, cy as i32, cz as i32);
            cx = a as usize;
            cy = b as usize;
            cz = c as usize;
        }
    }

    // Fine-level operator (VbuildA on the finest grid).
    let mut ac0 = vec![0.0f64; 4 * nf];
    let mut cc0 = vec![0.0f64; nf];
    let mut fc0 = vec![0.0f64; nf];
    crate::build_a::build_a(
        nx, ny, nz, 0, mgdisc, numdia,
        &mut ac0, &mut cc0, &mut fc0,
        xf, yf, zf, gxcf, gycf, gzcf,
        a1cf, a2cf, a3cf, ccf, fcf,
    );

    // Build the Galerkin hierarchy (VbuildP trilinear pc + VbuildG +
    // Vrestrc of cc/fc per level, as Vbuildops does for mgcoar = 2).
    let mut levels: Vec<crate::mgcs::Level> = Vec::with_capacity(nlev);
    levels.push(crate::mgcs::Level {
        nx, ny, nz,
        ac: ac0,
        cc: cc0,
        fc: fc0,
        numdia,
        pc: None,
        banded: None,
    });
    for lev in 1..nlev {
        let (pnx, pny, pnz) = sizes[lev - 1];
        let (cnx, cny, cnz) = sizes[lev];
        let npf = pnx * pny * pnz;
        let ncl = cnx * cny * cnz;
        let prev = &levels[lev - 1];

        let pc = crate::build_p::build_p_trilin_block(cnx, cny, cnz);
        let mut ac_c = vec![0.0f64; 14 * ncl];
        let fine_op = if prev.numdia == 4 {
            crate::build_g::FineOp::Seven(
                crate::matvec::split_bands4(&prev.ac, npf).expect("bands4"),
            )
        } else {
            crate::build_g::FineOp::TwentySeven(
                crate::matvec::split_bands14(&prev.ac, npf).expect("bands14"),
            )
        };
        crate::build_g::build_galerkin(pnx, pny, pnz, cnx, cny, cnz, &pc, fine_op, &mut ac_c);

        let mut cc_c = vec![0.0f64; ncl];
        let mut fc_c = vec![0.0f64; ncl];
        crate::matvec::restrc(pnx, pny, pnz, cnx, cny, cnz, &prev.cc, &mut cc_c, &pc);
        crate::matvec::restrc(pnx, pny, pnz, cnx, cny, cnz, &prev.fc, &mut fc_c, &pc);

        levels[lev - 1].pc = Some(pc);
        levels.push(crate::mgcs::Level {
            nx: cnx, ny: cny, nz: cnz,
            ac: ac_c,
            cc: cc_c,
            fc: fc_c,
            numdia: 14,
            pc: None,
            banded: None,
        });
    }

    // Factor the coarsest-level interior system (Vbuildband + Vdpbfa).
    // If the factorization fails, fall back to the iterative coarsest
    // solver exactly as Vbuildops does.
    let mut mgsolv_eff = mgsolv;
    {
        let last = levels.len() - 1;
        let l = &levels[last];
        let mut banded = crate::build_b::build_band(l.nx, l.ny, l.nz, &l.ac);
        let info = crate::lapack::dpbfa(&mut banded.abd, banded.lda, banded.n, banded.m);
        if info != 0 {
            mgsolv_eff = 0;
        } else {
            levels[last].banded = Some(banded);
        }
    }

    let epsiln = 2.2204460492503131e-16; // Vnm_epsmac
    let iters = crate::mgcs::mvcs(
        u, &levels,
        nu1, nu2,
        itmax, errtol,
        mgsolv_eff, epsiln,
    );

    if debug_enabled() {
        eprintln!("[DEBUG-MGDRV] linear mvcs iters={}", iters);
    }
}

/// Nonlinear path: Newton (default) or experimental FAS, unchanged from the
/// previous driver revision.
#[allow(clippy::too_many_arguments)]
fn nonlinear_mgdriv(
    iparm: &[i32], rparm: &[f64],
    u: &mut [f64],
    xf: &[f64], yf: &[f64], zf: &[f64],
    gxcf: &[f64], gycf: &[f64], gzcf: &[f64],
    a1cf: &[f64], a2cf: &[f64], a3cf: &[f64],
    ccf: &[f64], fcf: &[f64],
) {
    let nx = iparm[0] as usize;
    let ny = iparm[1] as usize;
    let nz = iparm[2] as usize;
    let nlev = iparm[3] as usize;
    let mgdisc = iparm[10];
    let mgsolv = if iparm[11] == 0 { 1 } else { iparm[11] };
    let nu1 = iparm[12];
    let nu2 = iparm[13];
    let omegal = rparm[0];
    let errtol = rparm[2];
    let itmax = iparm[7] as usize;
    let irite = iparm[16];
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

        let a1_c = inject_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_a1);
        let a2_c = inject_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_a2);
        let a3_c = inject_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_a3);
        let cc_c = inject_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_cc);
        let fc_c = inject_3d(cur_nx, cur_ny, cur_nz, nx_c, ny_c, nz_c, &cur_fc);
        let gxcf_c = inject_bc(cur_ny, cur_nz, ny_c, nz_c, &cur_gxcf);
        let gycf_c = inject_bc(cur_nx, cur_nz, nx_c, nz_c, &cur_gycf);
        let gzcf_c = inject_bc(cur_nx, cur_ny, nx_c, ny_c, &cur_gzcf);

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
    let mut w1 = vec![0.0f64; nf];
    let mut w2 = vec![0.0f64; nf];
    let mut r = vec![0.0f64; nf];
    let iz = vec![0i32; 50 * (nlev + 1)];
    let mut ipc = vec![0i32; 20 * (nlev + 1)];
    let mut rpc = vec![0.0f64; 20 * (nlev + 1)];
    let pc = vec![0.0f64; 27 * narr_total.max(1)];

    for lev in 0..nlev {
        ipc[lev * 20] = numdia;
    }
    {
        let hx0 = if nx > 1 { xf[1] - xf[0] } else { 1.0 };
        let hy0 = if ny > 1 { yf[1] - yf[0] } else { 1.0 };
        let hz0 = if nz > 1 { zf[1] - zf[0] } else { 1.0 };
        for lev in 0..nlev {
            let scale = (1usize << lev) as f64;
            rpc[lev * 20 + 0] = hx0 * scale;
            rpc[lev * 20 + 1] = hy0 * scale;
            rpc[lev * 20 + 2] = hz0 * scale;
            rpc[lev * 20 + 3] = 0.0;
        }
    }

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
}

/// Coincident-point injection restriction (used by the nonlinear branch's
/// standard-coarsening hierarchy, mirroring Vbuildcopy0).
fn inject_3d(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    fine: &[f64],
) -> Vec<f64> {
    let nc = nxc * nyc * nzc;
    let mut coarse = vec![0.0f64; nc];
    for kc in 0..nzc {
        for jc in 0..nyc {
            for ic in 0..nxc {
                let if_ = (2 * ic).min(nxf - 1);
                let jf = (2 * jc).min(nyf - 1);
                let kf = (2 * kc).min(nzf - 1);
                coarse[ic + jc * nxc + kc * nxc * nyc] =
                    fine[if_ + jf * nxf + kf * nxf * nyf];
            }
        }
    }
    coarse
}

fn inject_bc(
    nf1: usize, nf2: usize,
    nc1: usize, nc2: usize,
    g: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0f64; 2 * nc1 * nc2];
    for face in 0..2 {
        for k in 0..nc2 {
            for j in 0..nc1 {
                let fj = (2 * j).min(nf1 - 1);
                let fk = (2 * k).min(nf2 - 1);
                out[face * nc1 * nc2 + k * nc1 + j] = g[face * nf1 * nf2 + fk * nf1 + fj];
            }
        }
    }
    out
}
