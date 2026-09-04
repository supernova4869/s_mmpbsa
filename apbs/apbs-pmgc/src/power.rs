// APBS PMGC power - Power iteration for eigenvalue estimation
// Port of pmgc/powerd.c

/// Power iteration to estimate the spectral radius of the iteration matrix
pub fn power(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64],
    x: &mut [f64],
    w1: &mut [f64], _w2: &mut [f64],
    niter: i32,
) -> f64 {
    let n = nx * ny * nz;
    let mut lambda_old = 0.0;

    // Initialize x to random
    for i in 0..n {
        x[i] = 1.0;
    }

    for _ in 0..niter {
        // w1 = A * x
        crate::matvec::matvec(nx, ny, nz, ipc, rpc, ac, cc, &vec![0.0; n], x, w1);

        // lambda = x . w1 / x . x
        let xdotx = crate::blas::xdot(n, x, 0, x, 0);
        let xdotax = crate::blas::xdot(n, x, 0, w1, 0);

        if xdotx.abs() > 1.0e-30 {
            lambda_old = xdotax / xdotx;
        }

        // x = w1 / ||w1||
        let norm = crate::blas::xnrm2(n, w1, 0);
        if norm.abs() > 1.0e-30 {
            for i in 0..n {
                x[i] = w1[i] / norm;
            }
        }
    }

    lambda_old
}
