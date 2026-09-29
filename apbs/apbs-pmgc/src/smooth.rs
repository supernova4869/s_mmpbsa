// APBS PMGC smooth - Smoother dispatcher
// Port of pmgc/smoothd.c (Vsmooth)

/// Smooth the solution using Gauss-Seidel or CG.
/// `meth` follows Vsmooth: 0 = wjac (unsupported), 1 = Vgsrb,
/// 4 = Vcghs. The 27-point (numdia 14) operator is only valid with Vgsrb.
pub fn smooth(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    x: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    numdia: i32,
    nu: i32,       // number of smoothing iterations
    omega: f64,
    meth: i32,
    iresid: i32,
    iadjoint: i32,
) {
    if meth == 4 {
        // Conjugate Gradient (used to solve the coarsest level when
        // mgsolv == 0, mirroring Vmvcs's mgsmoo_s = 4 call)
        let errtol = 1.0e-8;
        let mut iters = 0;
        let rinf_norm = crate::blas::xnrm2(nx * ny * nz, fc, 0);
        crate::cg::cg(
            nx, ny, nz, ipc, rpc, ac, cc, fc,
            x, w1, w2, r,
            nu, &mut iters, errtol, rinf_norm,
        );
    } else {
        // Gauss-Seidel Red-Black (Vgsrb)
        let mut iters = 0;
        crate::gs::gsrb(
            nx, ny, nz, ipc, rpc, ac, cc, fc,
            x, w1, w2, r,
            nu, &mut iters, 0.0, omega, iresid, iadjoint, numdia,
        );
    }
}
