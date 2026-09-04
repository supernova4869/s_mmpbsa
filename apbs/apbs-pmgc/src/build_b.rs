// APBS PMGC buildB - Build banded matrix for coarse-grid solver
// Port of pmgc/buildBd.c

/// Build banded matrix representation for the symmetric banded linear solver
pub fn build_b(
    nx: usize, ny: usize, nz: usize,
    o_c: &[f64], o_e: &[f64], o_n: &[f64], u_c: &[f64],
    cc: &[f64],
    ab: &mut [f64],
    nband: usize,
) {
    let nxny = nx * ny;
    let _n = nx * ny * nz;

    // Clear banded matrix
    for v in ab.iter_mut() {
        *v = 0.0;
    }

    // Fill banded storage: ab[0..nband] is the diagonal band
    // For 7-point stencil, bandwidth = nx * ny (all z-neighbors are within this bandwidth)

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let ip = i + j * nx + k * nxny;

                // Diagonal
                let diag_idx = nband - 1 + ip * nband;
                if diag_idx < ab.len() {
                    ab[diag_idx] = o_c[ip] + cc[ip];
                }

                // x-neighbors (bandwidth offset 1)
                if i > 0 {
                    let idx = nband - 1 - 1 + ip * nband;
                    if idx < ab.len() {
                        ab[idx] = -o_e[ip - 1];
                    }
                }
                if i < nx - 1 {
                    let idx = nband - 1 + 1 + ip * nband;
                    if idx < ab.len() {
                        ab[idx] = -o_e[ip];
                    }
                }

                // y-neighbors (bandwidth offset nx)
                if j > 0 {
                    let idx = nband - 1 - nx as usize + ip * nband;
                    if idx < ab.len() {
                        ab[idx] = -o_n[ip - nx];
                    }
                }
                if j < ny - 1 {
                    let idx = nband - 1 + nx as usize + ip * nband;
                    if idx < ab.len() {
                        ab[idx] = -o_n[ip];
                    }
                }

                // z-neighbors (bandwidth offset nx*ny)
                if k > 0 {
                    let idx = nband - 1 - nxny as usize + ip * nband;
                    if idx < ab.len() {
                        ab[idx] = -u_c[ip - nxny];
                    }
                }
                if k < nz - 1 {
                    let idx = nband - 1 + nxny as usize + ip * nband;
                    if idx < ab.len() {
                        ab[idx] = -u_c[ip];
                    }
                }
            }
        }
    }
}
