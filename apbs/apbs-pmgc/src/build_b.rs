// APBS PMGC buildB - Build banded matrix for the coarsest-grid direct solver
// Port of pmgc/buildBd.c (Vbuildband / Vbuildband1_7 / Vbuildband1_27).
//
// Rows/columns cover the interior points only (i, j, k in 1..n-1), row
// order (k, j, i). The band offsets below are exactly those of
// Vbuildband1_27; for a 7-point operator the extra bands are zero and the
// identical table reduces to Vbuildband1_7.
//
// Note (faithful to the C): the banded diagonal stores oC WITHOUT the
// Helmholtz term cc.

use crate::matvec::{split_bands14, split_bands4};

pub struct Banded {
    pub abd: Vec<f64>,
    pub n: usize,
    pub m: usize,
    pub lda: usize,
}

/// Assemble the banded interior system of the coarsest level operator
/// (`ac`: 4 bands for a 7-point operator, 14 bands for 27-point).
pub fn build_band(nx: usize, ny: usize, nz: usize, ac: &[f64]) -> Banded {
    let n_grid = nx * ny * nz;
    let nix = nx - 2;
    let niy = ny - 2;
    let n = nix * niy * (nz - 2);
    let m = nix * niy + nix + 1;
    let lda = m + 1;
    let nxy = nix * niy;

    let mut abd = vec![0.0f64; lda * n];

    let fourteen = ac.len() >= 14 * n_grid;
    let b14 = fourteen.then(|| split_bands14(ac, n_grid)).flatten();
    let b4 = (!fourteen).then(|| split_bands4(ac, n_grid)).flatten();

    let mut jj = 0usize;
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                jj += 1;
                let row = jj - 1;
                let ip = i + j * nx + k * nx * ny;
                let (oc, oe_im1, on_jm1, one, onw, uc_km1, ue, uw, un, us, une, unw, use_, usw) =
                    if let Some(b) = b14 {
                        (
                            b.o_c[ip],
                            b.o_e[ip - 1],
                            b.o_n[ip - nx],
                            b.o_ne[ip - nx],
                            b.o_nw[ip - nx],
                            b.u_c[ip - nx * ny],
                            b.u_e[ip - nx * ny],
                            b.u_w[ip - nx * ny],
                            b.u_n[ip - nx * ny],
                            b.u_s[ip - nx * ny],
                            b.u_ne[ip - nx * ny],
                            b.u_nw[ip - nx * ny],
                            b.u_se[ip - nx * ny],
                            b.u_sw[ip - nx * ny],
                        )
                    } else {
                        let b = b4.expect("operator bands");
                        (
                            b.o_c[ip],
                            b.o_e[ip - 1],
                            b.o_n[ip - nx],
                            0.0,
                            0.0,
                            b.u_c[ip - nx * ny],
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        )
                    };

                // A(row - l, row) lives at abd[row * lda + m - l].
                let mut put = |l: usize, val: f64| {
                    abd[row * lda + m - l] = val;
                };
                put(0, oc);
                put(1, -oe_im1);
                put(nix, -on_jm1);
                put(nix - 1, -one);
                put(nix + 1, -onw);
                put(nxy, -uc_km1);
                put(nxy - 1, -ue);
                put(nxy + 1, -uw);
                put(nxy - nix, -un);
                put(nxy + nix, -us);
                put(nxy - nix - 1, -une);
                put(nxy - nix + 1, -unw);
                put(nxy + nix - 1, -use_);
                put(nxy + nix + 1, -usw);
            }
        }
    }

    Banded { abd, n, m, lda }
}
