// APBS PMGC Gauss-Seidel - Red-Black Gauss-Seidel smoother
// Port of pmgc/gsd.c (Vgsrb7x / Vgsrb27x) and pmgc/smoothd.c (Vsmooth)

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
                // Same-color grid points are not 7-point neighbors, so the
                // parallel writes below do not race with any read in this
                // color sweep (matches the OpenMP-parallel Vgsrb7x).
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
                i += 2;
            }
        }
    });
}

/// Gauss-Seidel Red-Black smoother (dispatcher).
/// Port of Vgsrb (gsd.c): `ac` holds 4 bands for numdia == 4
/// ([oC, oE, oN, uC]) or 14 bands for numdia == 14 (27-point Galerkin).
#[allow(clippy::too_many_arguments)]
pub fn gsrb(
    nx: usize, ny: usize, nz: usize,
    _ipc: &[i32], _rpc: &[f64],
    ac: &[f64], cc: &[f64], fc: &[f64],
    x: &mut [f64],
    w1: &mut [f64], w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, errtol: f64,
    omega: f64, iresid: i32, iadjoint: i32,
    numdia: i32,
) {
    let n = nx * ny * nz;
    if numdia == 4 {
        let b = crate::matvec::split_bands4(ac, n)
            .expect("7-point operator needs 4 bands");
        gsrb7x(
            nx, ny, nz, _ipc, _rpc, b, cc, fc,
            x, w1, w2, r, itmax, iters, errtol, omega, iresid, iadjoint,
        );
    } else {
        let b = crate::matvec::split_bands14(ac, n)
            .expect("27-point operator needs 14 bands");
        gsrb27x(
            nx, ny, nz, b, cc, fc,
            x, w1, w2, r, itmax, iters, errtol, omega, iresid, iadjoint,
        );
    }
}

/// 7-point Gauss-Seidel Red-Black smoother.
/// Port of Vgsrb7x (gsd.c): exactly `itmax` dual sweeps, residual returned
/// in `r` when `iresid` is set.
#[allow(clippy::too_many_arguments)]
pub fn gsrb7x(
    nx: usize, ny: usize, nz: usize,
    _ipc: &[i32], _rpc: &[f64],
    b: crate::matvec::Bands4, cc: &[f64], fc: &[f64],
    x: &mut [f64],
    _w1: &mut [f64], _w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, _errtol: f64,
    _omega: f64, iresid: i32, iadjoint: i32,
) {
    let nxny = nx * ny;
    let n = nx * ny * nz;
    let (o_c, o_e, o_n, u_c) = (b.o_c, b.o_e, b.o_n, b.u_c);
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
                            let mut rhs = fc[ip];
                            rhs += o_e[ip - 1] * x[ip - 1];
                            rhs += o_e[ip] * x[ip + 1];
                            rhs += o_n[ip - nx] * x[ip - nx];
                            rhs += o_n[ip] * x[ip + nx];
                            rhs += u_c[ip - nxny] * x[ip - nxny];
                            rhs += u_c[ip] * x[ip + nxny];
                            x[ip] = rhs / diag;
                            i += 2;
                        }
                    }
                }
            }
        }
    }

    if iresid != 0 {
        crate::matvec::mresid7_c(nx, ny, nz, &b, cc, fc, x, r);
    }
}

/// 27-point Gauss-Seidel Red-Black smoother.
/// Port of Vgsrb27x (gsd.c). The OpenMP pragma is commented out in the C
/// code because same-color 27-point neighbors exist: updates must run in
/// the exact sequential k/j/i order of the C loop.
#[allow(clippy::too_many_arguments)]
pub fn gsrb27x(
    nx: usize, ny: usize, nz: usize,
    b: crate::matvec::Bands14, cc: &[f64], fc: &[f64],
    x: &mut [f64],
    _w1: &mut [f64], _w2: &mut [f64], r: &mut [f64],
    itmax: i32, iters: &mut i32, _errtol: f64,
    _omega: f64, iresid: i32, iadjoint: i32,
) {
    let nxny = nx * ny;
    if nx < 3 || ny < 3 || nz < 3 {
        return;
    }

    for iter in 0..itmax {
        *iters = iter + 1;

        for color in 0..2 {
            for k in 1..(nz - 1) {
                for j in 1..(ny - 1) {
                    let mut i = 1 + ((color + iadjoint as usize + j + k) & 1);
                    while i < nx - 1 {
                        let ip = i + j * nx + k * nxny;
                        let tmp_o = b.o_n[ip] * x[ip + nx]
                            + b.o_n[ip - nx] * x[ip - nx]
                            + b.o_e[ip] * x[ip + 1]
                            + b.o_e[ip - 1] * x[ip - 1]
                            + b.o_ne[ip] * x[ip + nx + 1]
                            + b.o_nw[ip] * x[ip + nx - 1]
                            + b.o_nw[ip - nx + 1] * x[ip - nx + 1]
                            + b.o_ne[ip - nx - 1] * x[ip - nx - 1];
                        let tmp_u = b.u_c[ip] * x[ip + nxny]
                            + b.u_n[ip] * x[ip + nxny + nx]
                            + b.u_s[ip] * x[ip + nxny - nx]
                            + b.u_e[ip] * x[ip + nxny + 1]
                            + b.u_w[ip] * x[ip + nxny - 1]
                            + b.u_ne[ip] * x[ip + nxny + nx + 1]
                            + b.u_nw[ip] * x[ip + nxny + nx - 1]
                            + b.u_se[ip] * x[ip + nxny - nx + 1]
                            + b.u_sw[ip] * x[ip + nxny - nx - 1];
                        let tmp_d = b.u_c[ip - nxny] * x[ip - nxny]
                            + b.u_s[ip - nxny + nx] * x[ip - nxny + nx]
                            + b.u_n[ip - nxny - nx] * x[ip - nxny - nx]
                            + b.u_w[ip - nxny + 1] * x[ip - nxny + 1]
                            + b.u_e[ip - nxny - 1] * x[ip - nxny - 1]
                            + b.u_sw[ip - nxny + nx + 1] * x[ip - nxny + nx + 1]
                            + b.u_se[ip - nxny + nx - 1] * x[ip - nxny + nx - 1]
                            + b.u_nw[ip - nxny + 1 - nx] * x[ip - nxny + 1 - nx]
                            + b.u_ne[ip - nxny - nx - 1] * x[ip - nxny - nx - 1];
                        x[ip] = (fc[ip] + (tmp_o + tmp_u + tmp_d))
                            / (b.o_c[ip] + cc[ip]);
                        i += 2;
                    }
                }
            }
        }
    }

    if iresid != 0 {
        crate::matvec::mresid27_c(nx, ny, nz, &b, cc, fc, x, r);
    }
}
