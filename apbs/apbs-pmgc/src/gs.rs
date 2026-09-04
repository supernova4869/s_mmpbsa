// APBS PMGC Gauss-Seidel - Red-Black Gauss-Seidel smoother
// Port of pmgc/gsd.c

use rayon::prelude::*;

const PAR_THRESHOLD: usize = 16_384;

#[derive(Clone, Copy)]
struct GridPtr(*mut f64);

unsafe impl Send for GridPtr {}
unsafe impl Sync for GridPtr {}

impl GridPtr {
    #[inline]
    unsafe fn read(self, idx: usize) -> f64 {
        *self.0.add(idx)
    }

    #[inline]
    unsafe fn write(self, idx: usize, value: f64) {
        *self.0.add(idx) = value;
    }
}

#[inline]
fn update_color_7pt_parallel(
    nx: usize,
    ny: usize,
    nz: usize,
    nxny: usize,
    color: usize,
    iadjoint: i32,
    o_c: &[f64],
    cc: &[f64],
    fc: &[f64],
    o_e: &[f64],
    o_n: &[f64],
    u_c: &[f64],
    x: &mut [f64],
) {
    let xptr = GridPtr(x.as_mut_ptr());
    (1..(nz - 1)).into_par_iter().for_each(move |k| {
        for j in 1..(ny - 1) {
            let mut i = 1 + ((color + iadjoint as usize + j + k) & 1);
            while i < nx - 1 {
                let ip = i + j * nx + k * nxny;
                let diag = o_c[ip] + cc[ip];
                if diag.abs() > 1.0e-20 {
                    // Same-color grid points are not 7-point neighbors, so the
                    // parallel writes below do not race with any read in this color sweep.
                    unsafe {
                        let mut rhs = fc[ip];
                        rhs += o_e[ip - 1] * xptr.read(ip - 1);
                        rhs += o_e[ip] * xptr.read(ip + 1);
                        rhs += o_n[ip - nx] * xptr.read(ip - nx);
                        rhs += o_n[ip] * xptr.read(ip + nx);
                        rhs += u_c[ip - nxny] * xptr.read(ip - nxny);
                        rhs += u_c[ip] * xptr.read(ip + nxny);
                        xptr.write(ip, rhs / diag);
                    }
                }
                i += 2;
            }
        }
    });
}

/// Gauss-Seidel Red-Black smoother (dispatcher)
pub fn gsrb(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    o_c: &[f64], cc: &[f64], fc: &[f64],
    o_e: &[f64], o_n: &[f64], u_c: &[f64],
    x: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, errtol: f64,
    omega: f64, iresid: i32, iadjoint: i32,
    numdia: i32,
) {
    if numdia == 4 {
        gsrb7x(
            nx, ny, nz, ipc, rpc, o_c, cc, fc, o_e, o_n, u_c,
            x, w1, w2, r, itmax, iters, errtol, omega, iresid, iadjoint,
        );
    } else {
        gsrb27x(
            nx, ny, nz, ipc, rpc, o_c, cc, fc,
            o_e, o_n, u_c,
            x, w1, w2, r, itmax, iters, errtol, omega, iresid, iadjoint,
        );
    }
}

/// 7-point Gauss-Seidel Red-Black smoother
pub fn gsrb7x(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    o_c: &[f64], cc: &[f64], fc: &[f64],
    o_e: &[f64], o_n: &[f64], u_c: &[f64],
    x: &mut [f64],
    _w1: &mut [f64], _w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, errtol: f64,
    _omega: f64, iresid: i32, iadjoint: i32,
) {
    let nxny = nx * ny;
    let n = nx * ny * nz;
    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }
    let use_parallel = n >= PAR_THRESHOLD;

    for iter in 0..itmax {
        *iters = iter + 1;

        for color in 0..2 {
            if use_parallel {
                update_color_7pt_parallel(
                    nx, ny, nz, nxny, color, iadjoint,
                    o_c, cc, fc, o_e, o_n, u_c, x,
                );
            } else {
                for k in 1..(nz - 1) {
                    for j in 1..(ny - 1) {
                        let mut i = 1 + ((color + iadjoint as usize + j + k) & 1);
                        while i < nx - 1 {
                            let ip = i + j * nx + k * nxny;
                            let diag = o_c[ip] + cc[ip];
                            if diag.abs() > 1.0e-20 {
                                let mut rhs = fc[ip];
                                rhs += o_e[ip - 1] * x[ip - 1];
                                rhs += o_e[ip] * x[ip + 1];
                                rhs += o_n[ip - nx] * x[ip - nx];
                                rhs += o_n[ip] * x[ip + nx];
                                rhs += u_c[ip - nxny] * x[ip - nxny];
                                rhs += u_c[ip] * x[ip + nxny];
                                x[ip] = rhs / diag;
                            }
                            i += 2;
                        }
                    }
                }
            }
        }

        // Check convergence via residual if requested
        if iresid != 0 && errtol > 0.0 && iter > 0 && iter % 10 == 0 {
            crate::blas::mresid(
                nx, ny, nz, ipc, rpc,
                o_c, o_e, o_n, u_c, cc, x, fc, r,
            );
            let rnorm = crate::blas::xnrm2(nx * ny * nz, r, 0);
            if rnorm < errtol {
                break;
            }
        }
    }

    if iresid != 0 {
        crate::blas::mresid(
            nx, ny, nz, ipc, rpc,
            o_c, o_e, o_n, u_c, cc, x, fc, r,
        );
    }
}

/// 27-point Gauss-Seidel Red-Black smoother
pub fn gsrb27x(
    nx: usize, ny: usize, nz: usize,
    ipc: &[i32], rpc: &[f64],
    o_c: &[f64], cc: &[f64], fc: &[f64],
    o_e: &[f64], o_n: &[f64], u_c: &[f64],
    x: &mut [f64],
    _w1: &mut [f64], _w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, _errtol: f64,
    _omega: f64, iresid: i32, iadjoint: i32,
) {
    let nxny = nx * ny;
    let n = nx * ny * nz;
    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }
    let use_parallel = n >= PAR_THRESHOLD;

    for iter in 0..itmax {
        *iters = iter + 1;

        for color in 0..2 {
            if use_parallel {
                update_color_7pt_parallel(
                    nx, ny, nz, nxny, color, iadjoint,
                    o_c, cc, fc, o_e, o_n, u_c, x,
                );
            } else {
                for k in 1..(nz - 1) {
                    for j in 1..(ny - 1) {
                        let mut i = 1 + ((color + iadjoint as usize + j + k) & 1);
                        while i < nx - 1 {
                            let ip = i + j * nx + k * nxny;
                            let diag = o_c[ip] + cc[ip];
                            if diag.abs() > 1.0e-20 {
                                let mut rhs = fc[ip];
                                rhs += o_e[ip - 1] * x[ip - 1];
                                rhs += o_e[ip] * x[ip + 1];
                                rhs += o_n[ip - nx] * x[ip - nx];
                                rhs += o_n[ip] * x[ip + nx];
                                rhs += u_c[ip - nxny] * x[ip - nxny];
                                rhs += u_c[ip] * x[ip + nxny];
                                x[ip] = rhs / diag;
                            }
                            i += 2;
                        }
                    }
                }
            }
        }
    }

    if iresid != 0 {
        crate::blas::mresid(
            nx, ny, nz, ipc, rpc,
            o_c, o_e, o_n, u_c, cc, x, fc, r,
        );
    }
}
