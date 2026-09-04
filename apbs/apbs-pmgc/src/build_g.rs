// APBS PMGC buildG - placeholder Galerkin-style coarse-grid operator
// Porting note: this keeps the Rust PMGC 4-band storage layout
// [o_c, o_e, o_n, u_c] and performs weighted local coarsening.
// It is not APBS C-equivalent `P^T A P` assembly and is no longer used
// by the active 7-point solve path.

fn idx(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
    k * nx * ny + j * nx + i
}

fn trilinear_weight(i: usize, j: usize, k: usize) -> f64 {
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
    w
}

fn local_weight(pc_ff: &[f64], i: usize, j: usize, k: usize, nx: usize, ny: usize) -> f64 {
    let ip = idx(i, j, k, nx, ny);
    if ip < pc_ff.len() {
        pc_ff[ip]
    } else {
        trilinear_weight(i, j, k)
    }
}

fn weighted_avg_2x2x2(
    src: &[f64],
    pc_ff: &[f64],
    nxf: usize,
    nyf: usize,
    nzf: usize,
    if0: usize,
    jf0: usize,
    kf0: usize,
) -> f64 {
    let mut wsum = 0.0;
    let mut vsum = 0.0;
    for dk in 0..=1 {
        for dj in 0..=1 {
            for di in 0..=1 {
                let i = if0 + di;
                let j = jf0 + dj;
                let k = kf0 + dk;
                if i < nxf && j < nyf && k < nzf {
                    let w = local_weight(pc_ff, i, j, k, nxf, nyf);
                    let ip = idx(i, j, k, nxf, nyf);
                    wsum += w;
                    vsum += w * src[ip];
                }
            }
        }
    }
    if wsum > 0.0 { vsum / wsum } else { 0.0 }
}

/// Placeholder coarse-grid builder for the dormant 4-band Rust PMGC layout.
///
/// This is not APBS C-equivalent `P^T A P` assembly. It only supports the
/// current 4-band `[o_c, o_e, o_n, u_c]` storage shape and performs a local
/// weighted average on each band.
pub fn build_g(
    nxf: usize, nyf: usize, nzf: usize,
    nxc: usize, nyc: usize, nzc: usize,
    numdia: i32,
    pc_ff: &[f64],
    ac_ff: &[f64],
    ac: &mut [f64],
) {
    assert_eq!(
        numdia, 4,
        "build_g placeholder only supports the 4-band Rust PMGC layout; APBS C 14-band/27-point Galerkin assembly is not implemented"
    );
    let nf_f = nxf * nyf * nzf;
    let nf_c = nxc * nyc * nzc;
    if ac_ff.len() < 4 * nf_f || ac.len() < 4 * nf_c {
        return;
    }

    let (diag_f, rem_f) = ac_ff.split_at(nf_f);
    let (east_f, rem_f) = rem_f.split_at(nf_f);
    let (north_f, up_f) = rem_f.split_at(nf_f);

    let (diag_c, rem_c) = ac.split_at_mut(nf_c);
    let (east_c, rem_c) = rem_c.split_at_mut(nf_c);
    let (north_c, up_c) = rem_c.split_at_mut(nf_c);

    for kc in 0..nzc {
        for jc in 0..nyc {
            for ic in 0..nxc {
                let if0 = (2 * ic).min(nxf.saturating_sub(1));
                let jf0 = (2 * jc).min(nyf.saturating_sub(1));
                let kf0 = (2 * kc).min(nzf.saturating_sub(1));
                let ipc = idx(ic, jc, kc, nxc, nyc);

                diag_c[ipc] = weighted_avg_2x2x2(diag_f, pc_ff, nxf, nyf, nzf, if0, jf0, kf0);
                east_c[ipc] = weighted_avg_2x2x2(east_f, pc_ff, nxf, nyf, nzf, if0, jf0, kf0);
                north_c[ipc] = weighted_avg_2x2x2(north_f, pc_ff, nxf, nyf, nzf, if0, jf0, kf0);
                up_c[ipc] = weighted_avg_2x2x2(up_f, pc_ff, nxf, nyf, nzf, if0, jf0, kf0);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::build_g;

    #[test]
    fn build_g_preserves_constant_operator_bands() {
        let (nxf, nyf, nzf) = (5usize, 5usize, 5usize);
        let (nxc, nyc, nzc) = (3usize, 3usize, 3usize);
        let nf_f = nxf * nyf * nzf;
        let nf_c = nxc * nyc * nzc;

        let mut ac_ff = vec![0.0f64; 4 * nf_f];
        for i in 0..nf_f {
            ac_ff[i] = 10.0;
            ac_ff[nf_f + i] = 2.0;
            ac_ff[2 * nf_f + i] = 3.0;
            ac_ff[3 * nf_f + i] = 4.0;
        }
        let pc_ff = vec![1.0f64; nf_f];
        let mut ac_c = vec![0.0f64; 4 * nf_c];
        build_g(nxf, nyf, nzf, nxc, nyc, nzc, 4, &pc_ff, &ac_ff, &mut ac_c);

        for i in 0..nf_c {
            assert!((ac_c[i] - 10.0).abs() < 1.0e-12);
            assert!((ac_c[nf_c + i] - 2.0).abs() < 1.0e-12);
            assert!((ac_c[2 * nf_c + i] - 3.0).abs() < 1.0e-12);
            assert!((ac_c[3 * nf_c + i] - 4.0).abs() < 1.0e-12);
        }
    }

    #[test]
    #[should_panic(expected = "build_g placeholder only supports the 4-band Rust PMGC layout")]
    fn build_g_rejects_non_4band_numdia() {
        let mut ac_c = vec![0.0; 4];
        build_g(1, 1, 1, 1, 1, 1, 14, &[1.0], &[1.0, 2.0, 3.0, 4.0], &mut ac_c);
    }
}
