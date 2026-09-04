// APBS PMGC pack - Pack/unpack routines for multigrid data
// Port of pmgc/mikpckd.c and pmgc/mlinpckd.c

/// Pack integer parameters for multigrid levels
pub fn mikpack(
    iparm: &[i32],
    iz: &mut [i32],
    _ipc: &mut [i32],
    nlev: i32,
) {
    // Copy key parameters into iz array for each level
    for lev in 0..=nlev as usize {
        let base = lev * 50;
        iz[base] = iparm[0]; // nx
        iz[base + 1] = iparm[1]; // ny
        iz[base + 2] = iparm[2]; // nz
        iz[base + 3] = iparm[3]; // numdia
        iz[base + 4] = iparm[4]; // mgcoar
        iz[base + 5] = iparm[5]; // mgdisc
        iz[base + 6] = iparm[6]; // mgsolv
    }
}

/// Pack real parameters for multigrid levels
pub fn mrkpack(
    rparm: &[f64],
    rpc: &mut [f64],
    nlev: i32,
) {
    for lev in 0..=nlev as usize {
        let base = lev * 20;
        rpc[base] = rparm[0]; // hx
        rpc[base + 1] = rparm[1]; // hy
        rpc[base + 2] = rparm[2]; // hz
        rpc[base + 3] = rparm[3]; // zkappa2
        rpc[base + 4] = rparm[4]; // omegal
        rpc[base + 5] = rparm[5]; // omegan
    }
}

/// Linear pack: pack fine grid data into coarse grid storage
pub fn mlinpck(
    nx_f: usize, ny_f: usize, _nz_f: usize,
    nx_c: usize, ny_c: usize, nz_c: usize,
    fine: &[f64], coarse: &mut [f64],
) {
    let nxny_f = nx_f * ny_f;
    let nxny_c = nx_c * ny_c;

    for kc in 0..nz_c {
        for jc in 0..ny_c {
            for ic in 0..nx_c {
                let ipc = ic + jc * nx_c + kc * nxny_c;
                let ipf = 2 * ic + 2 * jc * nx_f + 2 * kc * nxny_f;
                coarse[ipc] = fine[ipf];
            }
        }
    }
}

/// Linear unpack: unpack coarse grid data into fine grid storage
pub fn mlinupck(
    nx_f: usize, ny_f: usize, _nz_f: usize,
    nx_c: usize, ny_c: usize, nz_c: usize,
    coarse: &[f64], fine: &mut [f64],
) {
    let nxny_f = nx_f * ny_f;
    let nxny_c = nx_c * ny_c;

    for kc in 0..nz_c {
        for jc in 0..ny_c {
            for ic in 0..nx_c {
                let ipc = ic + jc * nx_c + kc * nxny_c;
                let ipf = 2 * ic + 2 * jc * nx_f + 2 * kc * nxny_f;
                fine[ipf] = coarse[ipc];
            }
        }
    }
}
