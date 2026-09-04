// APBS PMGC buildStr - Build multigrid hierarchy structure
// Port of pmgc/mgsubd.c

/// Compute sizes of work arrays for the multigrid solver
pub fn mgsz(
    mgcoar: i32, mgdisc: i32, mgsolv: i32,
    nx: i32, ny: i32, nz: i32, nlev: i32,
    nxc: &mut i32, nyc: &mut i32, nzc: &mut i32,
    nf: &mut i32, nc: &mut i32, narr: &mut i32, narrc: &mut i32,
    n_rpc: &mut i32, n_iz: &mut i32, n_ipc: &mut i32,
    iretot: &mut i32, iintot: &mut i32,
) {
    let numdia = if mgdisc == 0 { 4 } else { 14 };

    // Compute coarsest grid dimensions
    let mut nxc_l = nx;
    let mut nyc_l = ny;
    let mut nzc_l = nz;
    for _ in 1..nlev {
        nxc_l = nxc_l / 2 + 1;
        nyc_l = nyc_l / 2 + 1;
        nzc_l = nzc_l / 2 + 1;
    }
    *nxc = nxc_l;
    *nyc = nyc_l;
    *nzc = nzc_l;

    *nf = nx * ny * nz;
    *nc = nxc_l * nyc_l * nzc_l;

    // Total unknowns across all levels
    let mut narr_l = 0;
    let mut nx_l = nx;
    let mut ny_l = ny;
    let mut nz_l = nz;
    for _ in 0..nlev {
        narr_l += nx_l * ny_l * nz_l;
        nx_l = nx_l / 2 + 1;
        ny_l = ny_l / 2 + 1;
        nz_l = nz_l / 2 + 1;
    }
    *narr = narr_l;
    *narrc = narr_l - *nf;

    // Parameter arrays
    *n_rpc = (nlev + 1) * 20;
    *n_iz = (nlev + 1) * 50;
    *n_ipc = (nlev + 1) * 20;

    // Operator storage
    let effective_mgcoar = mgcoar;

    let nac = if effective_mgcoar == 2 {
        14 * (*narrc) + numdia * (*nf)
    } else {
        numdia * narr_l
    };

    // Prolongation storage
    let npc = 27 * (*narrc);

    // Real work array
    let nb = if mgsolv == 1 {
        let nb_base = (nx - 2) * (ny - 2);
        if mgdisc == 0 { nb_base + 1 } else { nb_base * 3 + 1 }
    } else {
        0
    };

    *iretot = narr_l + nac + npc + (*nc) + (*nf) + nb * (*nc);
    *iintot = (*n_iz) + (*n_ipc) + narr_l + nac + npc + (*nc) + (*nf);
}

/// Compute per-level prolongation (`pc`) offsets for the APBS PMGC hierarchy.
///
/// The returned vector has length `nlev + 1`. Entry `lev` (1-based in APBS C,
/// 0-based here) gives the starting offset, in scalar `f64` entries, of the
/// 27-point prolongation block mapping level `lev` to `lev + 1`.
///
/// This follows the APBS C `iz(11, lev)` convention from `mgsubd.c`:
/// level 0 starts at 0, and each subsequent level advances by `27 * n_coarse`.
pub fn prolongation_offsets(nx: i32, ny: i32, nz: i32, nlev: i32) -> Vec<usize> {
    let levels = (nlev as usize).saturating_add(1);
    let mut offsets = vec![0usize; levels];
    if nlev <= 1 {
        return offsets;
    }

    let mut nx_l = nx;
    let mut ny_l = ny;
    let mut nz_l = nz;

    for lev in 1..nlev as usize {
        let (nx_c, ny_c, nz_c) = make_coarse(nx_l, ny_l, nz_l);
        let n_coarse = (nx_c * ny_c * nz_c) as usize;
        offsets[lev] = offsets[lev - 1] + 27 * n_coarse;
        nx_l = nx_c;
        ny_l = ny_c;
        nz_l = nz_c;
    }

    offsets
}

/// Compute coarse grid dimensions
pub fn make_coarse(
    nx_old: i32, ny_old: i32, nz_old: i32,
) -> (i32, i32, i32) {
    (nx_old / 2 + 1, ny_old / 2 + 1, nz_old / 2 + 1)
}

#[cfg(test)]
mod tests {
    use super::{make_coarse, prolongation_offsets};

    #[test]
    fn prolongation_offsets_match_apbs_c_lagging_level_layout() {
        let nlev = 4;
        let offsets = prolongation_offsets(33, 33, 33, nlev);
        assert_eq!(offsets.len(), 5);
        assert_eq!(offsets[0], 0);

        let (nx1, ny1, nz1) = make_coarse(33, 33, 33);
        let n1 = (nx1 * ny1 * nz1) as usize;
        assert_eq!(offsets[1], 27 * n1);

        let (nx2, ny2, nz2) = make_coarse(nx1, ny1, nz1);
        let n2 = (nx2 * ny2 * nz2) as usize;
        assert_eq!(offsets[2], 27 * (n1 + n2));

        let (nx3, ny3, nz3) = make_coarse(nx2, ny2, nz2);
        let n3 = (nx3 * ny3 * nz3) as usize;
        assert_eq!(offsets[3], 27 * (n1 + n2 + n3));
    }
}
