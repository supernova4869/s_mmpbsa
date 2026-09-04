// APBS PMGC BLAS - BLAS-like vector operations
// Port of pmgc/blas operations

use rayon::prelude::*;

const PAR_THRESHOLD: usize = 16_384;

/// Dot product of two vectors
pub fn xdot(n: usize, x: &[f64], ix: usize, y: &[f64], iy: usize) -> f64 {
    if n >= PAR_THRESHOLD && ix + n <= x.len() && iy + n <= y.len() {
        return x[ix..ix + n]
            .par_iter()
            .zip(&y[iy..iy + n])
            .map(|(a, b)| a * b)
            .sum();
    }
    let mut sum = 0.0;
    let mut xi = ix;
    let mut yi = iy;
    for _ in 0..n {
        sum += x[xi] * y[yi];
        xi += 1;
        yi += 1;
    }
    sum
}

/// y := alpha*x + y (AXPY)
pub fn xaxpy(n: usize, alpha: f64, x: &[f64], ix: usize, y: &mut [f64], iy: usize) {
    if n >= PAR_THRESHOLD && ix + n <= x.len() && iy + n <= y.len() {
        y[iy..iy + n]
            .par_iter_mut()
            .zip(&x[ix..ix + n])
            .for_each(|(yd, xs)| *yd += alpha * *xs);
        return;
    }
    let mut xi = ix;
    let mut yi = iy;
    for _ in 0..n {
        y[yi] += alpha * x[xi];
        xi += 1;
        yi += 1;
    }
}

/// y := x (copy)
pub fn xcopy(n: usize, x: &[f64], ix: usize, y: &mut [f64], iy: usize) {
    if n >= PAR_THRESHOLD && ix + n <= x.len() && iy + n <= y.len() {
        y[iy..iy + n]
            .par_iter_mut()
            .zip(&x[ix..ix + n])
            .for_each(|(yd, xs)| *yd = *xs);
        return;
    }
    let mut xi = ix;
    let mut yi = iy;
    for _ in 0..n {
        y[yi] = x[xi];
        xi += 1;
        yi += 1;
    }
}

/// x := alpha*x (scale)
pub fn xscal(n: usize, alpha: f64, x: &mut [f64], ix: usize) {
    if n >= PAR_THRESHOLD && ix + n <= x.len() {
        x[ix..ix + n].par_iter_mut().for_each(|v| *v *= alpha);
        return;
    }
    let mut xi = ix;
    for _ in 0..n {
        x[xi] *= alpha;
        xi += 1;
    }
}

/// Euclidean norm
pub fn xnrm2(n: usize, x: &[f64], ix: usize) -> f64 {
    if n >= PAR_THRESHOLD && ix + n <= x.len() {
        return x[ix..ix + n].par_iter().map(|v| v * v).sum::<f64>().sqrt();
    }
    let mut sum = 0.0;
    let mut xi = ix;
    for _ in 0..n {
        sum += x[xi] * x[xi];
        xi += 1;
    }
    sum.sqrt()
}

/// Compute residual: r = f - Au
pub fn mresid(
    nx: usize, ny: usize, nz: usize,
    _ipc: &[i32], _rpc: &[f64],
    o_c: &[f64], o_e: &[f64], o_n: &[f64], u_c: &[f64],
    cc: &[f64],
    x: &[f64], fc: &[f64], r: &mut [f64],
) {
    let n = nx * ny * nz;
    let r_len = n.min(r.len());
    r[..r_len].fill(0.0);
    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }
    let nxny = nx * ny;
    if n >= PAR_THRESHOLD {
        r[..n].par_chunks_mut(nxny).enumerate().for_each(|(k, r_plane)| {
            if k == 0 || k + 1 == nz {
                return;
            }
            for j in 1..(ny - 1) {
                for i in 1..(nx - 1) {
                    let local = i + j * nx;
                    let ip = local + k * nxny;
                    let mut ax = o_c[ip] * x[ip];
                    ax -= o_e[ip - 1] * x[ip - 1];
                    ax -= o_e[ip] * x[ip + 1];
                    ax -= o_n[ip - nx] * x[ip - nx];
                    ax -= o_n[ip] * x[ip + nx];
                    ax -= u_c[ip - nxny] * x[ip - nxny];
                    ax -= u_c[ip] * x[ip + nxny];
                    ax += cc[ip] * x[ip];
                    r_plane[local] = fc[ip] - ax;
                }
            }
        });
        return;
    }
    for k in 1..(nz - 1) {
        for j in 1..(ny - 1) {
            for i in 1..(nx - 1) {
                let ip = i + j * nx + k * nxny;
                let mut ax = o_c[ip] * x[ip];
                ax -= o_e[ip - 1] * x[ip - 1];
                ax -= o_e[ip] * x[ip + 1];
                ax -= o_n[ip - nx] * x[ip - nx];
                ax -= o_n[ip] * x[ip + nx];
                ax -= u_c[ip - nxny] * x[ip - nxny];
                ax -= u_c[ip] * x[ip + nxny];
                ax += cc[ip] * x[ip];
                r[ip] = fc[ip] - ax;
            }
        }
    }
}
