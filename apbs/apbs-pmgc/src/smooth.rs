// APBS PMGC smooth - Smoother dispatcher
// Port of pmgc/smoothd.c

/// Smooth the solution using Gauss-Seidel or CG
pub fn smooth(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    x: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    numdia: i32,
    nu: i32,       // number of smoothing iterations
    omega: f64,
    isolves: i32,  // solver type: 0=CG, 1=GSRB
    iadjoint: i32,
) {
    if isolves == 0 {
        // Conjugate Gradient
        let errtol = 1.0e-8;
        let mut iters = 0;
        let rinf_norm = crate::blas::xnrm2(nx * ny * nz, fc, 0);
        crate::cg::cg(
            nx, ny, nz, ipc, rpc, ac, cc, fc,
            x, w1, w2, r,
            nu, &mut iters, errtol, rinf_norm,
        );
    } else {
        // Gauss-Seidel Red-Black
        let errtol = 1.0e-8;
        let mut iters = 0;
        let o_c = &ac[0..nx * ny * nz];
        let o_e = &ac[nx * ny * nz..2 * nx * ny * nz];
        let o_n = &ac[2 * nx * ny * nz..3 * nx * ny * nz];
        let u_c = &ac[3 * nx * ny * nz..4 * nx * ny * nz];
        crate::gs::gsrb(
            nx, ny, nz, ipc, rpc, o_c, cc, fc, o_e, o_n, u_c,
            x, w1, w2, r,
            nu, &mut iters, errtol, omega, 0, iadjoint, numdia,
        );
    }
}
