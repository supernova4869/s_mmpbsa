// APBS PMGC gpu_mg - GPU implementation of the linear multigrid V-cycle
//
// Runs `mgcs::mvcs` (the LPBE Poisson-Boltzmann solver loop) on a wgpu device:
// per-level operators / right-hand sides / prolongation weights are uploaded
// once per solve and the red-black Gauss-Seidel sweeps, residuals,
// restriction, prolongation, operator applications and reductions all run as
// compute kernels.  The control flow (V-cycle order, Hackbusch/Reusken
// damping, stopping test) stays on the CPU, and the coarsest-level direct
// solve (dpbfa/dpbsl) stays on the CPU because it is tiny.
//
// Numerical notes:
// * Every kernel evaluates its stencil in exactly the operand order of the
//   CPU code, and the fine-level 7-point smoother keeps the CPU's red-black
//   pattern, so fine-level sweeps are bitwise identical.
// * The 27-point Galerkin levels cannot keep the sequential CPU visit order on
//   a GPU (same-color 27-point neighbors exist, as noted in the original
//   OpenMP comment); they use an 8-color ordering instead, which converges to
//   the same solution up to floating-point rounding.
// * Reductions (dot products, L1 norms) sum in a different order than the
//   CPU, which can flip the stopping test one iteration earlier or later.  The
//   resulting potentials agree with the CPU solve to around the solver
//   tolerance (errtol), far below MM-PBSA reporting precision.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::OnceLock;

use bytemuck::{Pod, Zeroable};
use wgpu::util::{BufferInitDescriptor, DeviceExt};

use crate::mgcs::Level;

static MG_GPU_REQUESTED: AtomicBool = AtomicBool::new(false);

/// Request GPU execution for the linear multigrid solve.  When no usable GPU
/// is found (or a solve fails) the driver falls back to the CPU path.
pub fn set_mg_gpu_requested(on: bool) {
    MG_GPU_REQUESTED.store(on, Ordering::Relaxed);
}

pub(crate) fn mg_gpu_requested() -> bool {
    MG_GPU_REQUESTED.load(Ordering::Relaxed)
}

/// Whether GPU execution was requested for the PB solver path (used by the
/// boundary-condition kernels too).
pub fn gpu_requested() -> bool {
    mg_gpu_requested()
}

/// Fine grids smaller than this stay on the CPU: upload/readback latency and
/// per-dispatch overhead dominate below roughly a 50^3 grid.
const MIN_FINE_POINTS: usize = 65536;

static SOLVER: OnceLock<Option<std::sync::Mutex<GpuMg>>> = OnceLock::new();

fn solver() -> Option<&'static std::sync::Mutex<GpuMg>> {
    SOLVER
        .get_or_init(|| match GpuMg::new() {
            Ok(s) => {
                let name = s.adapter_name.clone();
                println!("GPU enabled for the PB multigrid solver on: {}.", name);
                Some(std::sync::Mutex::new(s))
            }
            Err(e) => {
                println!(
                    "GPU requested for PB/SA but not usable ({}); the multigrid solver stays on CPU.",
                    e
                );
                None
            }
        })
        .as_ref()
}

/// Eagerly initializes the GPU solver (if requested) so the "GPU enabled"
/// message appears at startup instead of on the first multigrid solve.
pub fn warm_up() {
    if mg_gpu_requested() {
        let _ = solver();
    }
}

/// GPU dispatch of `mgcs::mvcs`.  Returns `None` — leaving `u` untouched —
/// whenever the GPU path does not apply or fails, so the caller can fall back
/// to the CPU implementation.
#[allow(clippy::too_many_arguments)]
pub fn mvcs_gpu(
    u: &mut [f64],
    levels: &[Level],
    nu1: i32,
    nu2: i32,
    itmax: i32,
    errtol: f64,
    mgsolv: i32,
) -> Option<i32> {
    if !mg_gpu_requested() || mgsolv != 1 || levels.len() < 2 {
        return None;
    }
    if levels[0].numdia != 4 {
        // Only the standard fine-grid 7-point operator is supported.
        return None;
    }
    if levels[0].nf() < MIN_FINE_POINTS {
        return None;
    }
    if levels[..levels.len() - 1].iter().any(|l| l.pc.is_none()) {
        // pc weights must exist on every level but the coarsest.
        return None;
    }
    let Ok(mut s) = solver()?.lock() else {
        return None;
    };
    match s.mvcs(u, levels, nu1, nu2, itmax, errtol) {
        Ok(iters) => Some(iters),
        Err(e) => {
            println!("GPU PB multigrid solve failed ({}); recomputing on CPU.", e);
            None
        }
    }
}

/// GPU port of `build_g::build_galerkin` for one hierarchy level: computes
/// the 14-band coarse operator `ac_out` from the packed fine operator and
/// the trilinear pc weights.  Returns false — leaving `ac_out` untouched —
/// when the GPU path does not apply or fails, so the caller can fall back
/// to the CPU implementation.
pub fn galerkin_gpu(
    fine_ac: &[f64],
    fine_numdia: u32,
    pc: &[f64],
    fine_dims: (usize, usize, usize),
    coarse_dims: (usize, usize, usize),
    ac_out: &mut [f64],
) -> bool {
    if !mg_gpu_requested() {
        return false;
    }
    let Some(s) = solver() else {
        return false;
    };
    let Ok(mut s) = s.lock() else {
        return false;
    };
    match s.galerkin(fine_ac, fine_numdia, pc, fine_dims, coarse_dims, ac_out) {
        Ok(()) => true,
        Err(e) => {
            static WARNED: OnceLock<()> = OnceLock::new();
            WARNED.get_or_init(|| {
                println!("GPU Galerkin build failed ({}); that step stays on CPU.", e);
            });
            false
        }
    }
}

/// GPU evaluation of the multi-ion Debye-Hueckel boundary condition
/// (`vpmg::bcfl_mdh`).  `atoms` holds (x, y, z, EC*q, size) per atom; the
/// coordinate vectors give the grid geometry.  Returns the three face arrays
/// (gxcf, gycf, gzcf), or None on any fallback condition.
pub fn dh_boundary_gpu(
    atoms: &[[f64; 5]],
    xf: &[f64],
    yf: &[f64],
    zf: &[f64],
    xkappa: f64,
    pre1: f64,
) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    if !gpu_requested() {
        return None;
    }
    let Ok(mut s) = solver()?.lock() else {
        return None;
    };
    match s.dh_boundary(atoms, xf, yf, zf, xkappa, pre1) {
        Ok(r) => Some(r),
        Err(e) => {
            static WARNED: OnceLock<()> = OnceLock::new();
            WARNED.get_or_init(|| {
                println!("GPU Debye-Hueckel boundary failed ({}); that step stays on CPU.", e);
            });
            None
        }
    }
}

// ---------------------------------------------------------------------------
// WGSL kernels
// ---------------------------------------------------------------------------

const WGSL_SHADER: &str = r#"
struct P {
    n: u32, nx: u32, ny: u32, nz: u32,
    nxf: u32, nyf: u32, nzf: u32, nc: u32,
    color: u32, iadjoint: u32, m: u32, pad0: u32,
    alpha: f64,
    pre1: f64,
}

@group(0) @binding(0) var<storage> p: P;
@group(0) @binding(1) var<storage> ac: array<f64>;
@group(0) @binding(2) var<storage> cc: array<f64>;
@group(0) @binding(3) var<storage> fc: array<f64>;
@group(0) @binding(4) var<storage> pc: array<f64>;
@group(0) @binding(5) var<storage, read_write> x: array<f64>;
@group(0) @binding(6) var<storage, read_write> aux: array<f64>;
@group(0) @binding(7) var<storage, read_write> outp: array<f64>;

// Band d of the packed operator: element d * n + point.
fn b(d: u32, ip: u32) -> f64 {
    return ac[d * p.n + ip];
}

fn coords(lin: u32) -> vec3<u32> {
    let i = lin % p.nx;
    let j = (lin / p.nx) % p.ny;
    let k = lin / (p.nx * p.ny);
    return vec3<u32>(i, j, k);
}

fn is_boundary(i: u32, j: u32, k: u32) -> bool {
    return i == 0u || j == 0u || k == 0u
        || i == p.nx - 1u || j == p.ny - 1u || k == p.nz - 1u;
}

// pc weight of `plane` at coarse point `ipc` times the fine value `fi`
// (used by the restriction kernel).
fn w(plane: u32, ipc: u32, fi: u32) -> f64 {
    return pc[plane * p.nc + ipc] * x[fi];
}

// Coarse vector value `off` elements past `cell` (used by interp).
fn xw(cell: u32, off: u32) -> f64 {
    return x[cell + off];
}

// Fine-grid index of (ifx, jfy, kfz) + (dx, dy, dz); the fine dims live in
// p.nxf / p.nyf.
fn fine_idx(ifx: u32, jfy: u32, kfz: u32, dx: i32, dy: i32, dz: i32) -> u32 {
    let i = u32(i32(ifx) + dx);
    let j = u32(i32(jfy) + dy);
    let k = u32(i32(kfz) + dz);
    return i + j * p.nxf + k * (p.nxf * p.nyf);
}

// One color pass of the 7-point red-black Gauss-Seidel sweep.  Keeps the
// CPU's parity pattern (i + j + k + color + iadjoint odd), so fine-level
// sweeps match the CPU bitwise.
@compute @workgroup_size(64)
fn gsrb7(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) { return; }
    if ((c.x + c.y + c.z + p.color + p.iadjoint) & 1u) != 1u { return; }
    let nx = p.nx;
    let nxny = nx * p.ny;
    let ip = lin;
    let diag = b(0u, ip) + cc[ip];
    var rhs = fc[ip];
    rhs += b(1u, ip - 1u) * x[ip - 1u];
    rhs += b(1u, ip) * x[ip + 1u];
    rhs += b(2u, ip - nx) * x[ip - nx];
    rhs += b(2u, ip) * x[ip + nx];
    rhs += b(3u, ip - nxny) * x[ip - nxny];
    rhs += b(3u, ip) * x[ip + nxny];
    x[ip] = rhs / diag;
}

// One color of the 8-color 27-point sweep (same-color 27-point neighbors make
// the CPU's sequential red-black order unparallelizable).  The per-point
// update is bitwise the CPU expression.
@compute @workgroup_size(64)
fn gsrb27(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) { return; }
    if ((c.x & 1u) + 2u * (c.y & 1u) + 4u * (c.z & 1u) != p.color) { return; }
    let nx = p.nx;
    let nxny = nx * p.ny;
    let ip = lin;
    let tmp_o = b(2u, ip) * x[ip + nx]
        + b(2u, ip - nx) * x[ip - nx]
        + b(1u, ip) * x[ip + 1u]
        + b(1u, ip - 1u) * x[ip - 1u]
        + b(4u, ip) * x[ip + nx + 1u]
        + b(5u, ip) * x[ip + nx - 1u]
        + b(5u, ip - nx + 1u) * x[ip - nx + 1u]
        + b(4u, ip - nx - 1u) * x[ip - nx - 1u];
    let tmp_u = b(3u, ip) * x[ip + nxny]
        + b(8u, ip) * x[ip + nxny + nx]
        + b(9u, ip) * x[ip + nxny - nx]
        + b(6u, ip) * x[ip + nxny + 1u]
        + b(7u, ip) * x[ip + nxny - 1u]
        + b(10u, ip) * x[ip + nxny + nx + 1u]
        + b(11u, ip) * x[ip + nxny + nx - 1u]
        + b(12u, ip) * x[ip + nxny - nx + 1u]
        + b(13u, ip) * x[ip + nxny - nx - 1u];
    let tmp_d = b(3u, ip - nxny) * x[ip - nxny]
        + b(9u, ip - nxny + nx) * x[ip - nxny + nx]
        + b(8u, ip - nxny - nx) * x[ip - nxny - nx]
        + b(7u, ip - nxny + 1u) * x[ip - nxny + 1u]
        + b(6u, ip - nxny - 1u) * x[ip - nxny - 1u]
        + b(13u, ip - nxny + nx + 1u) * x[ip - nxny + nx + 1u]
        + b(12u, ip - nxny + nx - 1u) * x[ip - nxny + nx - 1u]
        + b(11u, ip - nxny + 1u - nx) * x[ip - nxny + 1u - nx]
        + b(10u, ip - nxny - nx - 1u) * x[ip - nxny - nx - 1u];
    x[ip] = (fc[ip] + (tmp_o + tmp_u + tmp_d)) / (b(0u, ip) + cc[ip]);
}

// Residual aux = fc - A x, 7-point, interior only (boundary untouched),
// operand order of Vmresid7_1s.
@compute @workgroup_size(64)
fn mresid7(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) { return; }
    let nx = p.nx;
    let nxny = nx * p.ny;
    let ip = lin;
    aux[ip] = fc[ip]
        + b(2u, ip) * x[ip + nx]
        + b(2u, ip - nx) * x[ip - nx]
        + b(1u, ip) * x[ip + 1u]
        + b(1u, ip - 1u) * x[ip - 1u]
        + b(3u, ip - nxny) * x[ip - nxny]
        + b(3u, ip) * x[ip + nxny]
        - (b(0u, ip) + cc[ip]) * x[ip];
}

// Residual aux = fc - A x, 27-point, interior only, operand order of
// Vmresid27_1s.
@compute @workgroup_size(64)
fn mresid27(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) { return; }
    let nx = p.nx;
    let nxny = nx * p.ny;
    let ip = lin;
    let tmp_o = b(2u, ip) * x[ip + nx]
        + b(2u, ip - nx) * x[ip - nx]
        + b(1u, ip) * x[ip + 1u]
        + b(1u, ip - 1u) * x[ip - 1u]
        + b(4u, ip) * x[ip + nx + 1u]
        + b(5u, ip) * x[ip + nx - 1u]
        + b(5u, ip - nx + 1u) * x[ip - nx + 1u]
        + b(4u, ip - nx - 1u) * x[ip - nx - 1u];
    let tmp_u = b(3u, ip) * x[ip + nxny]
        + b(8u, ip) * x[ip + nxny + nx]
        + b(9u, ip) * x[ip + nxny - nx]
        + b(6u, ip) * x[ip + nxny + 1u]
        + b(7u, ip) * x[ip + nxny - 1u]
        + b(10u, ip) * x[ip + nxny + nx + 1u]
        + b(11u, ip) * x[ip + nxny + nx - 1u]
        + b(12u, ip) * x[ip + nxny - nx + 1u]
        + b(13u, ip) * x[ip + nxny - nx - 1u];
    let tmp_d = b(3u, ip - nxny) * x[ip - nxny]
        + b(9u, ip - nxny + nx) * x[ip - nxny + nx]
        + b(8u, ip - nxny - nx) * x[ip - nxny - nx]
        + b(7u, ip - nxny + 1u) * x[ip - nxny + 1u]
        + b(6u, ip - nxny - 1u) * x[ip - nxny - 1u]
        + b(13u, ip - nxny + nx + 1u) * x[ip - nxny + nx + 1u]
        + b(12u, ip - nxny + nx - 1u) * x[ip - nxny + nx - 1u]
        + b(11u, ip - nxny + 1u - nx) * x[ip - nxny + 1u - nx]
        + b(10u, ip - nxny - nx - 1u) * x[ip - nxny - nx - 1u];
    aux[ip] = fc[ip] + tmp_o + tmp_u + tmp_d - (b(0u, ip) + cc[ip]) * x[ip];
}

// aux = A x, 7-point; interior computed, boundary zeroed (matches the CPU's
// fill(0) + interior loop).
@compute @workgroup_size(64)
fn matvec7(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) {
        aux[lin] = 0.0;
        return;
    }
    let nx = p.nx;
    let nxny = nx * p.ny;
    let ip = lin;
    var sum = b(0u, ip) * x[ip];
    sum -= b(1u, ip - 1u) * x[ip - 1u];
    sum -= b(1u, ip) * x[ip + 1u];
    sum -= b(2u, ip - nx) * x[ip - nx];
    sum -= b(2u, ip) * x[ip + nx];
    sum -= b(3u, ip - nxny) * x[ip - nxny];
    sum -= b(3u, ip) * x[ip + nxny];
    sum += cc[ip] * x[ip];
    aux[ip] = sum;
}

// aux = A x, 27-point, interior only, operand order of Vmatvec27_1s.
@compute @workgroup_size(64)
fn matvec27(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) { return; }
    let nx = p.nx;
    let nxny = nx * p.ny;
    let ip = lin;
    let tmp_o = -(b(2u, ip) * x[ip + nx]
        + b(2u, ip - nx) * x[ip - nx]
        + b(1u, ip) * x[ip + 1u]
        + b(1u, ip - 1u) * x[ip - 1u]
        + b(4u, ip) * x[ip + nx + 1u]
        + b(5u, ip) * x[ip + nx - 1u]
        + b(5u, ip - nx + 1u) * x[ip - nx + 1u]
        + b(4u, ip - nx - 1u) * x[ip - nx - 1u]);
    let tmp_u = -(b(3u, ip) * x[ip + nxny]
        + b(8u, ip) * x[ip + nxny + nx]
        + b(9u, ip) * x[ip + nxny - nx]
        + b(6u, ip) * x[ip + nxny + 1u]
        + b(7u, ip) * x[ip + nxny - 1u]
        + b(10u, ip) * x[ip + nxny + nx + 1u]
        + b(11u, ip) * x[ip + nxny + nx - 1u]
        + b(12u, ip) * x[ip + nxny - nx + 1u]
        + b(13u, ip) * x[ip + nxny - nx - 1u]);
    let tmp_d = -(b(3u, ip - nxny) * x[ip - nxny]
        + b(9u, ip - nxny + nx) * x[ip - nxny + nx]
        + b(8u, ip - nxny - nx) * x[ip - nxny - nx]
        + b(7u, ip - nxny + 1u) * x[ip - nxny + 1u]
        + b(6u, ip - nxny - 1u) * x[ip - nxny - 1u]
        + b(13u, ip - nxny + nx + 1u) * x[ip - nxny + nx + 1u]
        + b(12u, ip - nxny + nx - 1u) * x[ip - nxny + nx - 1u]
        + b(11u, ip - nxny + 1u - nx) * x[ip - nxny + 1u - nx]
        + b(10u, ip - nxny - nx - 1u) * x[ip - nxny - nx - 1u]);
    aux[ip] = tmp_o + tmp_u + tmp_d + (b(0u, ip) + cc[ip]) * x[ip];
}

// Restriction with the 27-component pc weights: coarse output for every
// coarse point (boundary zeroed), operand order of Vrestrc2.  x is the fine
// vector; (nxf, nyf) are the fine dims; nc is the coarse point count.
@compute @workgroup_size(64)
fn restrc(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) {
        aux[lin] = 0.0;
        return;
    }
    let ipc = lin;
    let ifx = 2u * c.x;
    let jfy = 2u * c.y;
    let kfz = 2u * c.z;
    let f000 = fine_idx(ifx, jfy, kfz, 0, 0, 0);
    let f010 = fine_idx(ifx, jfy, kfz, 0, 1, 0);
    let f0m0 = fine_idx(ifx, jfy, kfz, 0, -1, 0);
    let f100 = fine_idx(ifx, jfy, kfz, 1, 0, 0);
    let fm00 = fine_idx(ifx, jfy, kfz, -1, 0, 0);
    let f110 = fine_idx(ifx, jfy, kfz, 1, 1, 0);
    let fm10 = fine_idx(ifx, jfy, kfz, -1, 1, 0);
    let f1m0 = fine_idx(ifx, jfy, kfz, 1, -1, 0);
    let fmm0 = fine_idx(ifx, jfy, kfz, -1, -1, 0);
    let f001 = fine_idx(ifx, jfy, kfz, 0, 0, 1);
    let f011 = fine_idx(ifx, jfy, kfz, 0, 1, 1);
    let f0m1 = fine_idx(ifx, jfy, kfz, 0, -1, 1);
    let f101 = fine_idx(ifx, jfy, kfz, 1, 0, 1);
    let fm01 = fine_idx(ifx, jfy, kfz, -1, 0, 1);
    let f111 = fine_idx(ifx, jfy, kfz, 1, 1, 1);
    let fm11 = fine_idx(ifx, jfy, kfz, -1, 1, 1);
    let f1m1 = fine_idx(ifx, jfy, kfz, 1, -1, 1);
    let fmm1 = fine_idx(ifx, jfy, kfz, -1, -1, 1);
    let f00m = fine_idx(ifx, jfy, kfz, 0, 0, -1);
    let f01m = fine_idx(ifx, jfy, kfz, 0, 1, -1);
    let f0mm = fine_idx(ifx, jfy, kfz, 0, -1, -1);
    let f10m = fine_idx(ifx, jfy, kfz, 1, 0, -1);
    let fm0m = fine_idx(ifx, jfy, kfz, -1, 0, -1);
    let f11m = fine_idx(ifx, jfy, kfz, 1, 1, -1);
    let fm1m = fine_idx(ifx, jfy, kfz, -1, 1, -1);
    let f1mm = fine_idx(ifx, jfy, kfz, 1, -1, -1);
    let fmmm = fine_idx(ifx, jfy, kfz, -1, -1, -1);
    let tmp_o = w(0u, ipc, f000)
        + w(1u, ipc, f010)
        + w(2u, ipc, f0m0)
        + w(3u, ipc, f100)
        + w(4u, ipc, fm00)
        + w(5u, ipc, f110)
        + w(6u, ipc, fm10)
        + w(7u, ipc, f1m0)
        + w(8u, ipc, fmm0);
    let tmp_u = w(9u, ipc, f001)
        + w(10u, ipc, f011)
        + w(11u, ipc, f0m1)
        + w(12u, ipc, f101)
        + w(13u, ipc, fm01)
        + w(14u, ipc, f111)
        + w(15u, ipc, fm11)
        + w(16u, ipc, f1m1)
        + w(17u, ipc, fmm1);
    let tmp_d = w(18u, ipc, f00m)
        + w(19u, ipc, f01m)
        + w(20u, ipc, f0mm)
        + w(21u, ipc, f10m)
        + w(22u, ipc, fm0m)
        + w(23u, ipc, f11m)
        + w(24u, ipc, fm1m)
        + w(25u, ipc, f1mm)
        + w(26u, ipc, fmmm);
    aux[ipc] = tmp_o + tmp_u + tmp_d;
}

// Prolongation with the pc weights: one thread per FINE point computes its
// value from the owning coarse cell exactly as the CPU's per-cell writes
// assign it.  x is the coarse vector; (nxf, nyf, nzf) are the coarse dims.
// Boundaries are zeroed.
@compute @workgroup_size(64)
fn interp(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) {
        aux[lin] = 0.0;
        return;
    }
    let ncx = p.nxf;
    let ncy = p.nyf;
    let ic = c.x >> 1u;
    let jc = c.y >> 1u;
    let kc = c.z >> 1u;
    // The CPU only visits cells 0..n-2 per axis; corners owned by the skipped
    // last cell stay zero.
    if (ic + 1u >= ncx || jc + 1u >= ncy || kc + 1u >= p.nzf) {
        aux[lin] = 0.0;
        return;
    }
    let cell = ic + jc * ncx + kc * ncx * ncy;
    let nc = p.nc;
    let nxny = ncx * ncy;
    let dx = c.x & 1u;
    let dy = c.y & 1u;
    let dz = c.z & 1u;
    var v: f64;
    if (dx == 0u && dy == 0u && dz == 0u) {
        v = xw(cell, 0u);
    } else if (dx == 1u && dy == 0u && dz == 0u) {
        v = pc[3u * nc + cell] * xw(cell, 0u) + pc[4u * nc + cell + 1u] * xw(cell, 1u);
    } else if (dx == 0u && dy == 1u && dz == 0u) {
        v = pc[1u * nc + cell] * xw(cell, 0u) + pc[2u * nc + cell + ncx] * xw(cell, ncx);
    } else if (dx == 0u && dy == 0u && dz == 1u) {
        v = pc[9u * nc + cell] * xw(cell, 0u) + pc[18u * nc + cell + nxny] * xw(cell, nxny);
    } else if (dx == 1u && dy == 1u && dz == 0u) {
        v = pc[5u * nc + cell] * xw(cell, 0u)
            + pc[6u * nc + cell + 1u] * xw(cell, 1u)
            + pc[7u * nc + cell + ncx] * xw(cell, ncx)
            + pc[8u * nc + cell + ncx + 1u] * xw(cell, ncx + 1u);
    } else if (dx == 1u && dy == 0u && dz == 1u) {
        v = pc[12u * nc + cell] * xw(cell, 0u)
            + pc[13u * nc + cell + 1u] * xw(cell, 1u)
            + pc[21u * nc + cell + nxny] * xw(cell, nxny)
            + pc[22u * nc + cell + nxny + 1u] * xw(cell, nxny + 1u);
    } else if (dx == 0u && dy == 1u && dz == 1u) {
        v = pc[10u * nc + cell] * xw(cell, 0u)
            + pc[11u * nc + cell + ncx] * xw(cell, ncx)
            + pc[19u * nc + cell + nxny] * xw(cell, nxny)
            + pc[20u * nc + cell + ncx + nxny] * xw(cell, ncx + nxny);
    } else {
        v = pc[14u * nc + cell] * xw(cell, 0u)
            + pc[15u * nc + cell + 1u] * xw(cell, 1u)
            + pc[16u * nc + cell + ncx] * xw(cell, ncx)
            + pc[17u * nc + cell + ncx + 1u] * xw(cell, ncx + 1u)
            + pc[23u * nc + cell + nxny] * xw(cell, nxny)
            + pc[24u * nc + cell + nxny + 1u] * xw(cell, nxny + 1u)
            + pc[25u * nc + cell + ncx + nxny] * xw(cell, ncx + nxny)
            + pc[26u * nc + cell + ncx + nxny + 1u] * xw(cell, ncx + nxny + 1u);
    }
    aux[lin] = v;
}

// x += alpha * aux over the interior (Vxaxpy); x is the accumulator.
@compute @workgroup_size(64)
fn axpy(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) { return; }
    x[lin] = x[lin] + p.alpha * aux[lin];
}

// Zero the six Dirichlet faces.
@compute @workgroup_size(64)
fn zero_boundary(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin >= p.n) { return; }
    let c = coords(lin);
    if (is_boundary(c.x, c.y, c.z)) {
        x[lin] = 0.0;
    }
}

@compute @workgroup_size(64)
fn fill_zero(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x < p.n) {
        x[gid.x] = 0.0;
    }
}

var<workgroup> wshm: array<f64, 256>;

// Interior-point L1 norm: outp[workgroup] = sum |x|, one partial per
// workgroup; the CPU sums the partials.
@compute @workgroup_size(256)
fn reduce_abs(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
) {
    var v = f64(0.0);
    if (gid.x < p.m) {
        let dix = p.nx - 2u;
        let diy = p.ny - 2u;
        let i = 1u + gid.x % dix;
        let j = 1u + (gid.x / dix) % diy;
        let k = 1u + gid.x / (dix * diy);
        v = abs(x[i + j * p.nx + k * p.nx * p.ny]);
    }
    wshm[lid.x] = v;
    reduce_workgroup(lid.x);
    if (lid.x == 0u) {
        outp[wid.x] = wshm[0];
    }
}

// Interior-point dot product of x and aux, same reduction scheme.
@compute @workgroup_size(256)
fn reduce_dot(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
) {
    var v = f64(0.0);
    if (gid.x < p.m) {
        let dix = p.nx - 2u;
        let diy = p.ny - 2u;
        let i = 1u + gid.x % dix;
        let j = 1u + (gid.x / dix) % diy;
        let k = 1u + gid.x / (dix * diy);
        let ip = i + j * p.nx + k * p.nx * p.ny;
        v = x[ip] * aux[ip];
    }
    wshm[lid.x] = v;
    reduce_workgroup(lid.x);
    if (lid.x == 0u) {
        outp[wid.x] = wshm[0];
    }
}

fn reduce_workgroup(lid: u32) {
    var sw = 256u;
    loop {
        sw = sw >> 1u;
        workgroupBarrier();
        if (lid < sw) {
            wshm[lid] = wshm[lid] + wshm[lid + sw];
        }
        if (sw == 1u) { break; }
    }
    workgroupBarrier();
}

// ---------------------------------------------------------------------------
// Galerkin coarse-operator assembly (A_c = P^T A_f P), the GPU counterpart
// of `build_g::build_galerkin`.  One thread per coarse point assembles the
// 27 local stencil entries and writes the 14 packed bands.  Accumulation
// order matches the CPU loops (dz/dy/dx support loop, x/y/z owner
// distribution), so results agree with the CPU bit-for-bit up to the last
// rounding of each sum.  p: n/nx/ny/nz = coarse grid, nxf/nyf/nzf = fine
// grid, color = fine numdia (4 or 14); ac is the packed fine operator and
// pc the 27 weight planes.
// ---------------------------------------------------------------------------

var<private> OFFDIAG: array<vec3<i32>, 26> = array<vec3<i32>, 26>(
    vec3<i32>(1, 0, 0),
    vec3<i32>(-1, 0, 0),
    vec3<i32>(0, 1, 0),
    vec3<i32>(0, -1, 0),
    vec3<i32>(0, 0, 1),
    vec3<i32>(0, 0, -1),
    vec3<i32>(1, 1, 0),
    vec3<i32>(-1, 1, 0),
    vec3<i32>(1, -1, 0),
    vec3<i32>(-1, -1, 0),
    vec3<i32>(1, 0, 1),
    vec3<i32>(-1, 0, 1),
    vec3<i32>(0, 1, 1),
    vec3<i32>(0, -1, 1),
    vec3<i32>(1, 1, 1),
    vec3<i32>(-1, 1, 1),
    vec3<i32>(1, -1, 1),
    vec3<i32>(-1, -1, 1),
    vec3<i32>(1, 0, -1),
    vec3<i32>(-1, 0, -1),
    vec3<i32>(0, 1, -1),
    vec3<i32>(0, -1, -1),
    vec3<i32>(1, 1, -1),
    vec3<i32>(-1, 1, -1),
    vec3<i32>(1, -1, -1),
    vec3<i32>(-1, -1, -1),
);

var<private> STORED: array<vec3<i32>, 13> = array<vec3<i32>, 13>(
    vec3<i32>(1, 0, 0),
    vec3<i32>(0, 1, 0),
    vec3<i32>(0, 0, 1),
    vec3<i32>(1, 1, 0),
    vec3<i32>(-1, 1, 0),
    vec3<i32>(1, 0, 1),
    vec3<i32>(-1, 0, 1),
    vec3<i32>(0, 1, 1),
    vec3<i32>(0, -1, 1),
    vec3<i32>(1, 1, 1),
    vec3<i32>(-1, 1, 1),
    vec3<i32>(1, -1, 1),
    vec3<i32>(-1, -1, 1),
);

// Plane index of the pc weight for offset (dx, dy, dz); mirrors
// build_g::plane_index.
fn plane27(dx: i32, dy: i32, dz: i32) -> u32 {
    var layer: u32;
    if (dz == 0) {
        layer = 0u;
    } else if (dz > 0) {
        layer = 9u;
    } else {
        layer = 18u;
    }
    var row: u32;
    if (dy == 0) {
        row = select(select(4u, 3u, dx == 1), 0u, dx == 0);
    } else if (dy == 1) {
        row = select(select(6u, 5u, dx == 1), 1u, dx == 0);
    } else {
        row = select(select(8u, 7u, dx == 1), 2u, dx == 0);
    }
    return layer + row;
}

fn loc_index(dx: i32, dy: i32, dz: i32) -> u32 {
    return u32((dx + 1) + 3 * (dy + 1) + 9 * (dz + 1));
}

// Owner candidates of a fine coordinate along one axis (t >= 0).
fn owner1(t: i32) -> i32 {
    if (t % 2 == 0) { return t / 2; }
    return (t - 1) / 2;
}
fn woff1(t: i32) -> i32 {
    if (t % 2 == 0) { return 0; }
    return 1;
}
fn has_owner2(t: i32) -> bool {
    return t % 2 != 0;
}
fn owner2(t: i32) -> i32 {
    return (t + 1) / 2;
}

// Fine operator band `band` at fine point `idx` (p.nxf/p.nyf/p.nzf are the
// fine dims).
fn fband(band: u32, idx: u32) -> f64 {
    let nf = p.nxf * p.nyf * p.nzf;
    return ac[band * nf + idx];
}

// Fine operator coefficient A[a][a+d] in the packed 14-band layout; mirrors
// build_g::FineOp::offdiag.
fn fine_coef(a: u32, dx: i32, dy: i32, dz: i32, is27: bool) -> f64 {
    let nx = p.nxf;
    let nxny = p.nxf * p.nyf;
    if (dx == 1 && dy == 0 && dz == 0) { return fband(1u, a); }
    if (dx == -1 && dy == 0 && dz == 0) { return fband(1u, a - 1u); }
    if (dx == 0 && dy == 1 && dz == 0) { return fband(2u, a); }
    if (dx == 0 && dy == -1 && dz == 0) { return fband(2u, a - nx); }
    if (!is27) {
        if (dx == 0 && dy == 0 && dz == 1) { return fband(3u, a); }
        if (dx == 0 && dy == 0 && dz == -1) { return fband(3u, a - nxny); }
        return f64(0.0);
    }
    if (dx == 0 && dy == 0 && dz == 1) { return fband(3u, a); }
    if (dx == 0 && dy == 0 && dz == -1) { return fband(3u, a - nxny); }
    if (dx == 1 && dy == 1 && dz == 0) { return fband(4u, a); }
    if (dx == -1 && dy == 1 && dz == 0) { return fband(5u, a); }
    if (dx == 1 && dy == -1 && dz == 0) { return fband(5u, a - nx + 1u); }
    if (dx == -1 && dy == -1 && dz == 0) { return fband(4u, a - nx - 1u); }
    if (dx == 1 && dy == 0 && dz == 1) { return fband(6u, a); }
    if (dx == -1 && dy == 0 && dz == 1) { return fband(7u, a); }
    if (dx == 0 && dy == 1 && dz == 1) { return fband(8u, a); }
    if (dx == 0 && dy == -1 && dz == 1) { return fband(9u, a); }
    if (dx == 1 && dy == 1 && dz == 1) { return fband(10u, a); }
    if (dx == -1 && dy == 1 && dz == 1) { return fband(11u, a); }
    if (dx == 1 && dy == -1 && dz == 1) { return fband(12u, a); }
    if (dx == -1 && dy == -1 && dz == 1) { return fband(13u, a); }
    if (dx == 1 && dy == 0 && dz == -1) { return fband(7u, a - nxny + 1u); }
    if (dx == -1 && dy == 0 && dz == -1) { return fband(6u, a - nxny - 1u); }
    if (dx == 0 && dy == 1 && dz == -1) { return fband(9u, a - nxny + nx); }
    if (dx == 0 && dy == -1 && dz == -1) { return fband(8u, a - nxny - nx); }
    if (dx == 1 && dy == 1 && dz == -1) { return fband(13u, a - nxny + nx + 1u); }
    if (dx == -1 && dy == 1 && dz == -1) { return fband(12u, a - nxny + nx - 1u); }
    if (dx == 1 && dy == -1 && dz == -1) { return fband(11u, a - nxny + 1u - nx); }
    if (dx == -1 && dy == -1 && dz == -1) { return fband(10u, a - nxny - nx - 1u); }
    return f64(0.0);
}

// The Galerkin assembly for one coarse point.  Written without `continue`
// (nested-loop continue handling has been fragile in naga's WGSL front end):
// every filter is a boolean that gates the accumulation instead.
@compute @workgroup_size(64)
fn galerkin(@builtin(global_invocation_id) gid: vec3<u32>) {
    let lin = gid.x;
    if (lin < p.n) {
        let c = coords(lin);
        let interior = !is_boundary(c.x, c.y, c.z);
        let nc = p.n;
        let ncx = p.nx;
        let ncy = p.ny;
        let ncz = p.nz;
        let is27 = p.color != 4u;
        let ic = i32(c.x);
        let jc = i32(c.y);
        let kc = i32(c.z);
        var loc: array<f64, 27>;
        for (var i = 0; i < 27; i++) {
            loc[i] = f64(0.0);
        }

        if (interior) {
            for (var dz = -1; dz <= 1; dz++) {
                for (var dy = -1; dy <= 1; dy++) {
                    for (var dx = -1; dx <= 1; dx++) {
                        let w = pc[plane27(dx, dy, dz) * nc + lin];
                        let ax = 2 * ic + dx;
                        let ay = 2 * jc + dy;
                        let az = 2 * kc + dz;
                        let in_grid = w != f64(0.0)
                            && ax >= 0 && ay >= 0 && az >= 0
                            && ax < i32(p.nxf) && ay < i32(p.nyf) && az < i32(p.nzf);
                        if (in_grid) {
                            let a = u32(ax) + u32(ay) * p.nxf + u32(az) * (p.nxf * p.nyf);

                            var xs: array<i32, 2>;
                            var xs2: array<bool, 2>;
                            var xw: array<i32, 2>;
                            xs[0] = owner1(ax);
                            xs[1] = owner2(ax);
                            xs2[0] = true;
                            xs2[1] = has_owner2(ax);
                            xw[0] = woff1(ax);
                            xw[1] = -1;
                            var ys: array<i32, 2>;
                            var ys2: array<bool, 2>;
                            var yw: array<i32, 2>;
                            ys[0] = owner1(ay);
                            ys[1] = owner2(ay);
                            ys2[0] = true;
                            ys2[1] = has_owner2(ay);
                            yw[0] = woff1(ay);
                            yw[1] = -1;
                            var zs: array<i32, 2>;
                            var zs2: array<bool, 2>;
                            var zw: array<i32, 2>;
                            zs[0] = owner1(az);
                            zs[1] = owner2(az);
                            zs2[0] = true;
                            zs2[1] = has_owner2(az);
                            zw[0] = woff1(az);
                            zw[1] = -1;

                            // p == q: the fine diagonal couples every support
                            // containing a.
                            let val0 = w * fband(0u, a);
                            for (var sx = 0; sx < 2; sx++) {
                                for (var sy = 0; sy < 2; sy++) {
                                    for (var sz = 0; sz < 2; sz++) {
                                        let ok = (sx == 0 || xs2[sx])
                                            && (sy == 0 || ys2[sy])
                                            && (sz == 0 || zs2[sz])
                                            && xs[sx] < i32(ncx)
                                            && ys[sy] < i32(ncy)
                                            && zs[sz] < i32(ncz);
                                        if (ok) {
                                            let c2 = u32(xs[sx]) + u32(ys[sy]) * ncx
                                                + u32(zs[sz]) * (ncx * ncy);
                                            let ex = xs[sx] - ic;
                                            let ey = ys[sy] - jc;
                                            let ez = zs[sz] - kc;
                                            let w2 = pc[plane27(xw[sx], yw[sy], zw[sz]) * nc + c2];
                                            loc[loc_index(ex, ey, ez)] += val0 * w2;
                                        }
                                    }
                                }
                            }

                            for (var k = 0u; k < 26u; k++) {
                                let dd = OFFDIAG[k];
                                let bx = ax + dd.x;
                                let by = ay + dd.y;
                                let bz = az + dd.z;
                                let b_in = bx >= 0 && by >= 0 && bz >= 0
                                    && bx < i32(p.nxf) && by < i32(p.nyf) && bz < i32(p.nzf);
                                let coef = -fine_coef(a, dd.x, dd.y, dd.z, is27);
                                if (b_in && coef != f64(0.0)) {
                                    let val = w * coef;
                                    var bxs: array<i32, 2>;
                                    var bxs2: array<bool, 2>;
                                    var bxw: array<i32, 2>;
                                    bxs[0] = owner1(bx);
                                    bxs[1] = owner2(bx);
                                    bxs2[0] = true;
                                    bxs2[1] = has_owner2(bx);
                                    bxw[0] = woff1(bx);
                                    bxw[1] = -1;
                                    var bys: array<i32, 2>;
                                    var bys2: array<bool, 2>;
                                    var byw: array<i32, 2>;
                                    bys[0] = owner1(by);
                                    bys[1] = owner2(by);
                                    bys2[0] = true;
                                    bys2[1] = has_owner2(by);
                                    byw[0] = woff1(by);
                                    byw[1] = -1;
                                    var bzs: array<i32, 2>;
                                    var bzs2: array<bool, 2>;
                                    var bzw: array<i32, 2>;
                                    bzs[0] = owner1(bz);
                                    bzs[1] = owner2(bz);
                                    bzs2[0] = true;
                                    bzs2[1] = has_owner2(bz);
                                    bzw[0] = woff1(bz);
                                    bzw[1] = -1;
                                    for (var sx = 0; sx < 2; sx++) {
                                        for (var sy = 0; sy < 2; sy++) {
                                            for (var sz = 0; sz < 2; sz++) {
                                                let ok = (sx == 0 || bxs2[sx])
                                                    && (sy == 0 || bys2[sy])
                                                    && (sz == 0 || bzs2[sz])
                                                    && bxs[sx] < i32(ncx)
                                                    && bys[sy] < i32(ncy)
                                                    && bzs[sz] < i32(ncz);
                                                if (ok) {
                                                    let c2 = u32(bxs[sx]) + u32(bys[sy]) * ncx
                                                        + u32(bzs[sz]) * (ncx * ncy);
                                                    let ex = bxs[sx] - ic;
                                                    let ey = bys[sy] - jc;
                                                    let ez = bzs[sz] - kc;
                                                    let w2 = pc[plane27(bxw[sx], byw[sy], bzw[sz]) * nc + c2];
                                                    loc[loc_index(ex, ey, ez)] += val * w2;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Store the 14 packed bands; off-diagonal bands negated to the
        // positive band convention, exactly as the CPU does.
        aux[lin] = loc[13u];
        for (var k = 0u; k < 13u; k++) {
            let d = STORED[k];
            aux[(k + 1u) * nc + lin] = -loc[loc_index(d.x, d.y, d.z)];
        }
    }
}

// ---------------------------------------------------------------------------
// Multi-ion Debye-Hueckel boundary conditions: one thread per boundary point
// sums the screened Coulomb potentials of all atoms.  p: n = total face
// points, nx/ny/nz = grid dims, nxf = atom count, alpha = xkappa,
// pre1 = EC/(4 pi eps0 eps_w kB T) * 1e10.  ac holds the atoms as
// (x, y, z, q*EC, size) quintuplets; pc packs the two face values per axis
// ([0] x-lo, [1] x-hi, [2] y-lo, [3] y-hi, [4] z-lo, [5] z-hi) followed by
// the xf, yf, zf coordinate vectors; x is the output buffer ordered as
// [x-faces | y-faces | z-faces].
// ---------------------------------------------------------------------------

@compute @workgroup_size(64)
fn dh_boundary(@builtin(global_invocation_id) gid: vec3<u32>) {
    let m = gid.x;
    if (m >= p.n) { return; }
    let nx = p.nx;
    let ny = p.ny;
    let nz = p.nz;
    let natoms = p.nxf;

    let x_face_pts = 2u * ny * nz;
    let y_face_pts = 2u * nx * nz;
    var gx: f64;
    var gy: f64;
    var gz: f64;
    if (m < x_face_pts) {
        let face = m / (ny * nz);
        let r = m % (ny * nz);
        let k = r / ny;
        let j = r % ny;
        gx = select(pc[0u], pc[1u], face == 1u);
        gy = pc[6u + nx + j];
        gz = pc[6u + nx + ny + k];
    } else if (m < x_face_pts + y_face_pts) {
        let m2 = m - x_face_pts;
        let face = m2 / (nx * nz);
        let r = m2 % (nx * nz);
        let k = r / nx;
        let i = r % nx;
        gx = pc[6u + i];
        gy = select(pc[2u], pc[3u], face == 1u);
        gz = pc[6u + nx + ny + k];
    } else {
        let m2 = m - x_face_pts - y_face_pts;
        let face = m2 / (nx * ny);
        let r = m2 % (nx * ny);
        let j = r / nx;
        let i = r % nx;
        gx = pc[6u + i];
        gy = pc[6u + nx + j];
        gz = select(pc[4u], pc[5u], face == 1u);
    }

    var acc = f64(0.0);
    for (var ia = 0u; ia < natoms; ia += 1u) {
        let base = ia * 5u;
        let dx = gx - ac[base + 0u];
        let dy = gy - ac[base + 1u];
        let dz = gz - ac[base + 2u];
        let dist = sqrt(dx * dx + dy * dy + dz * dz);
        let charge = ac[base + 3u];
        let size = ac[base + 4u];
        // Hardware exp: the ~5e-8 relative deviation from libm is far
        // below the multigrid stopping tolerance (1e-6) and costs a fraction
        // of the high-accuracy Taylor path.
        let screened = exp(-p.alpha * (dist - size));
        var term: f64;
        if (p.alpha > f64(0.0)) {
            term = p.pre1 * (charge / dist) * screened / (f64(1.0) + p.alpha * size);
        } else {
            term = p.pre1 * (charge / dist);
        }
        acc += term;
    }
    x[m] = acc;
}

"#;

// Rust driver
// ---------------------------------------------------------------------------

const WORKGROUP: u32 = 64;
const RED_WORKGROUP: u32 = 256;
const PARAM_SLOTS: usize = 4096;
const PARAM_SLOT_SIZE: usize = 256;

/// Mirror of the WGSL `P` struct (56 bytes, storage layout).
#[repr(C)]
#[derive(Copy, Clone, Pod, Zeroable, Default)]
struct Params {
    n: u32,
    nx: u32,
    ny: u32,
    nz: u32,
    nxf: u32,
    nyf: u32,
    nzf: u32,
    nc: u32,
    color: u32,
    iadjoint: u32,
    m: u32,
    pad0: u32,
    alpha: f64,
    pre1: f64,
}

macro_rules! pipelines {
    ($($field:ident => $entry:literal),* $(,)?) => {
        struct Pipelines {
            $($field: wgpu::ComputePipeline,)*
        }
        impl Pipelines {
            fn new(
                device: &wgpu::Device,
                bgl: &wgpu::BindGroupLayout,
                module: &wgpu::ShaderModule,
            ) -> Pipelines {
                let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: Some("PB multigrid pipeline layout"),
                    bind_group_layouts: &[bgl],
                    push_constant_ranges: &[],
                });
                Pipelines {
                    $($field: device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                        label: Some($entry),
                        layout: Some(&layout),
                        module,
                        entry_point: Some($entry),
                        compilation_options: Default::default(),
                        cache: None,
                    }),)*
                }
            }
        }
    };
}

pipelines! {
    dh_boundary => "dh_boundary",
    galerkin => "galerkin",
    gsrb7 => "gsrb7",
    gsrb27 => "gsrb27",
    mresid7 => "mresid7",
    mresid27 => "mresid27",
    matvec7 => "matvec7",
    matvec27 => "matvec27",
    restrc => "restrc",
    interp => "interp",
    axpy => "axpy",
    zero_boundary => "zero_boundary",
    fill_zero => "fill_zero",
    reduce_abs => "reduce_abs",
    reduce_dot => "reduce_dot",
}

pub struct GpuMg {
    device: wgpu::Device,
    queue: wgpu::Queue,
    bgl: wgpu::BindGroupLayout,
    pipelines: Pipelines,
    adapter_name: String,
    /// Ring of parameter blocks, one slot per dispatch (dynamic offsets).
    params_buf: wgpu::Buffer,
    params_stage: Vec<u8>,
    /// Reduction partials (one f64 per workgroup, two regions for the paired
    /// damping dot products) and their MAP_READ staging copy.
    partials: wgpu::Buffer,
    partials_read: wgpu::Buffer,
    max_partials: u32,
    /// Byte stride between the two partials regions (256-aligned, as storage
    /// buffer binding offsets require).
    partials_stride: u64,
    /// Placeholder for unused read-only slots.
    dummy: wgpu::Buffer,
    /// Placeholder for unused read-write slots (must not alias `dummy`).
    dummy_rw: wgpu::Buffer,
}

/// One segment of recorded GPU work between two synchronizations.
struct Seg {
    encoder: wgpu::CommandEncoder,
    slot: usize,
}

#[allow(clippy::upper_case_acronyms)]
struct LevelBufs {
    nx: u32,
    ny: u32,
    nz: u32,
    n: u32,
    numdia: u32,
    ac: wgpu::Buffer,
    cc: wgpu::Buffer,
    fc: wgpu::Buffer,
    pc: Option<wgpu::Buffer>,
    x: wgpu::Buffer,
    w0: wgpu::Buffer,
    w1: wgpu::Buffer,
    w2: wgpu::Buffer,
}

fn interior_count(nx: u32, ny: u32, nz: u32) -> u32 {
    if nx < 3 || ny < 3 || nz < 3 {
        0
    } else {
        (nx - 2) * (ny - 2) * (nz - 2)
    }
}

impl GpuMg {
    fn new() -> Result<GpuMg, String> {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .map_err(|e| format!("no GPU adapter found: {}", e))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err(format!(
                "adapter '{}' lacks 64-bit float shader support",
                adapter.get_info().name
            ));
        }
        // PB grids need storage buffers well above the 128 MB default limit;
        // request what the adapter can give.  Every buffer is still checked
        // against the final limits per solve, with CPU fallback otherwise.
        let mut limits = wgpu::Limits::default();
        limits.max_buffer_size = limits.max_buffer_size.max(adapter.limits().max_buffer_size);
        limits.max_storage_buffer_binding_size = limits
            .max_storage_buffer_binding_size
            .max(adapter.limits().max_storage_buffer_binding_size);
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("s_mmpbsa PB multigrid device"),
            required_features: wgpu::Features::SHADER_F64,
            required_limits: limits,
            memory_hints: Default::default(),
            trace: Default::default(),
        }))
        .map_err(|e| format!("cannot open GPU device: {}", e))?;

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("s_mmpbsa PB multigrid shader"),
            source: wgpu::ShaderSource::Wgsl(std::borrow::Cow::Borrowed(WGSL_SHADER)),
        });
        let storage_read = |binding: u32, dynamic: bool| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: true },
                has_dynamic_offset: dynamic,
                min_binding_size: None,
            },
            count: None,
        };
        let storage_rw = |binding: u32| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Storage { read_only: false },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("s_mmpbsa PB multigrid bind group layout"),
            entries: &[
                storage_read(0, true),  // params ring (dynamic offset)
                storage_read(1, false), // ac operator bands
                storage_read(2, false), // cc Helmholtz term
                storage_read(3, false), // fc right-hand side
                storage_read(4, false), // pc prolongation weights
                storage_rw(5),          // x (primary vector)
                storage_rw(6),          // aux (second vector / kernel output)
                storage_rw(7),          // reduction partials
            ],
        });
        let pipelines = Pipelines::new(&device, &bgl, &module);
        let info = adapter.get_info();
        let adapter_name = format!("{} ({:?})", info.name, info.backend);
        let dummy = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg dummy"),
            size: 256,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let params_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg params ring"),
            size: (PARAM_SLOTS * PARAM_SLOT_SIZE) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let partials = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg partials"),
            size: 8,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let partials_read = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg partials readback"),
            size: 8,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let dummy_rw = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg dummy rw"),
            size: 256,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        Ok(GpuMg {
            device,
            queue,
            bgl,
            pipelines,
            adapter_name,
            params_buf,
            params_stage: vec![0u8; PARAM_SLOTS * PARAM_SLOT_SIZE],
            partials,
            partials_read,
            max_partials: 1,
            partials_stride: 256,
            dummy,
            dummy_rw,
        })
    }

    fn ensure_partials(&mut self, partials: u32) {
        let need = partials.max(1);
        if need <= self.max_partials {
            return;
        }
        let stride = (((need as u64) * 8).div_ceil(256)) * 256;
        self.partials = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg partials"),
            size: stride * 2,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        self.partials_read = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg partials readback"),
            size: stride * 2,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.max_partials = need;
        self.partials_stride = stride;
    }

    /// The eight bind-group entries shared by every kernel: params ring
    /// (dynamic offset), packed operator, Helmholtz term, RHS, pc weights,
    /// primary vector, aux vector and reduction partials region.
    #[allow(clippy::too_many_arguments)]
    fn bg_entries<'a>(
        &'a self,
        ac: &'a wgpu::Buffer,
        cc: &'a wgpu::Buffer,
        fc: &'a wgpu::Buffer,
        pc: Option<&'a wgpu::Buffer>,
        x: &'a wgpu::Buffer,
        aux: &'a wgpu::Buffer,
        out_off: u64,
    ) -> [wgpu::BindGroupEntry<'a>; 8] {
        let dummy = &self.dummy;
        let pc = pc.unwrap_or(dummy);
        [
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                    buffer: &self.params_buf,
                    offset: 0,
                    size: std::num::NonZeroU64::new(PARAM_SLOT_SIZE as u64),
                }),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: ac.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: cc.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: fc.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 4,
                resource: pc.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 5,
                resource: x.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 6,
                resource: aux.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 7,
                resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                    buffer: &self.partials,
                    offset: out_off,
                    size: std::num::NonZeroU64::new((self.max_partials as u64) * 8),
                }),
            },
        ]
    }

    /// Bind group for one dispatch on a multigrid level: `fc` is the RHS the
    /// kernels should use (the level's own or the restricted residual) and
    /// `out_off` selects the reduction region in bytes.
    #[allow(clippy::too_many_arguments)]
    fn make_bg(
        &self,
        level: &LevelBufs,
        fc: &wgpu::Buffer,
        x: &wgpu::Buffer,
        aux: &wgpu::Buffer,
        out_off: u64,
    ) -> wgpu::BindGroup {
        let entries = self.bg_entries(
            &level.ac,
            &level.cc,
            fc,
            level.pc.as_ref(),
            x,
            aux,
            out_off,
        );
        self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("PB mg bind group"),
            layout: &self.bgl,
            entries: &entries,
        })
    }

    /// GPU DH boundary (see `dh_boundary_gpu`).
    #[allow(clippy::too_many_arguments)]
    fn dh_boundary(
        &mut self,
        atoms: &[[f64; 5]],
        xf: &[f64],
        yf: &[f64],
        zf: &[f64],
        xkappa: f64,
        pre1: f64,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), String> {
        let (nx, ny, nz) = (xf.len(), yf.len(), zf.len());
        let natoms = atoms.len();
        let n_x = 2 * ny * nz;
        let n_y = 2 * nx * nz;
        let n_z = 2 * nx * ny;
        let total = n_x + n_y + n_z;
        self.fits(atoms.len() * 5)?;

        let mut fine_pc = Vec::with_capacity(6 + nx + ny + nz);
        fine_pc.extend_from_slice(&[xf[0], xf[nx - 1], yf[0], yf[ny - 1], zf[0], zf[nz - 1]]);
        fine_pc.extend_from_slice(xf);
        fine_pc.extend_from_slice(yf);
        fine_pc.extend_from_slice(zf);
        self.fits(fine_pc.len())?;

        let atoms_buf = self.device.create_buffer_init(&BufferInitDescriptor {
            label: Some("PB mg dh atoms"),
            contents: bytemuck::cast_slice(atoms),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let pc_buf = self.device.create_buffer_init(&BufferInitDescriptor {
            label: Some("PB mg dh pc"),
            contents: bytemuck::cast_slice(&fine_pc),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let out = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg dh out"),
            size: (total * 8) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let dummy_ro = self.dummy.clone();
        let dummy_rw = self.dummy_rw.clone();
        let entries = self.bg_entries(
            &atoms_buf,
            &dummy_ro,
            &dummy_ro,
            Some(&pc_buf),
            &out,
            &dummy_rw,
            0,
        );
        let bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("PB mg dh bind group"),
            layout: &self.bgl,
            entries: &entries,
        });
        let params = Params {
            n: total as u32,
            nx: nx as u32,
            ny: ny as u32,
            nz: nz as u32,
            nxf: natoms as u32,
            nyf: 0,
            nzf: 0,
            nc: 0,
            color: 0,
            iadjoint: 0,
            m: 0,
            pad0: 0,
            alpha: xkappa,
            pre1,
        };
        let mut seg = self.begin_seg();
        let pipe = self.pipelines.dh_boundary.clone();
        self.dispatch(&mut seg, &pipe, &bg, params, (total as u32).div_ceil(WORKGROUP));
        self.flush(seg);

        let values = self.read_buffer(&out, total)?;
        let mut gx = vec![0.0f64; n_x];
        let mut gy = vec![0.0f64; n_y];
        let mut gz = vec![0.0f64; n_z];
        gx.copy_from_slice(&values[..n_x]);
        gy.copy_from_slice(&values[n_x..n_x + n_y]);
        gz.copy_from_slice(&values[n_x + n_y..]);
        Ok((gx, gy, gz))
    }

    /// Storage-buffer size check against the device limits.
    fn fits(&self, len: usize) -> Result<(), String> {
        let bytes = (len * 8) as u64;
        let lim = self.device.limits();
        if bytes > lim.max_buffer_size || bytes > (lim.max_storage_buffer_binding_size as u64) {
            return Err(format!("grid buffer of {} bytes exceeds GPU limits", bytes));
        }
        Ok(())
    }

    fn begin_seg(&self) -> Seg {
        Seg {
            encoder: self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("PB mg segment"),
                }),
            slot: 0,
        }
    }

    fn dispatch(
        &mut self,
        seg: &mut Seg,
        pipeline: &wgpu::ComputePipeline,
        bg: &wgpu::BindGroup,
        params: Params,
        workgroups: u32,
    ) {
        assert!(
            seg.slot < PARAM_SLOTS,
            "PB multigrid parameter ring exhausted"
        );
        let offset = seg.slot * PARAM_SLOT_SIZE;
        self.params_stage[offset..offset + std::mem::size_of::<Params>()]
            .copy_from_slice(bytemuck::bytes_of(&params));
        // One pass per dispatch: pass boundaries synchronize storage access
        // between kernels, which the V-cycle's producer/consumer chains need.
        {
            let mut pass = seg.encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("PB mg pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, bg, &[offset as u32]);
            pass.dispatch_workgroups(workgroups, 1, 1);
        }
        seg.slot += 1;
    }

    /// Uploads the used parameter-ring slots and submits the segment.
    fn flush(&mut self, seg: Seg) {
        let used = seg.slot;
        let encoder = seg.encoder;
        if used > 0 {
            self.queue.write_buffer(
                &self.params_buf,
                0,
                &self.params_stage[..used * PARAM_SLOT_SIZE],
            );
        }
        self.queue.submit(Some(encoder.finish()));
    }

    /// Reads back `n` f64 elements of a storage buffer.
    fn read_buffer(&self, buf: &wgpu::Buffer, n: usize) -> Result<Vec<f64>, String> {
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg readback"),
            size: (n * 8) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("PB mg readback copy"),
            });
        enc.copy_buffer_to_buffer(buf, 0, &staging, 0, (n * 8) as u64);
        self.queue.submit(Some(enc.finish()));
        map_read(&self.device, &staging, n * 8)
    }

    /// GPU Galerkin product for one hierarchy level (see `galerkin_gpu`).
    fn galerkin(
        &mut self,
        fine_ac: &[f64],
        fine_numdia: u32,
        pc: &[f64],
        fine_dims: (usize, usize, usize),
        coarse_dims: (usize, usize, usize),
        ac_out: &mut [f64],
    ) -> Result<(), String> {
        let (nxf, nyf, nzf) = fine_dims;
        let (nxc, nyc, nzc) = coarse_dims;
        let nf = nxf * nyf * nzf;
        let nc = nxc * nyc * nzc;
        if fine_ac.len() < fine_numdia as usize * nf || pc.len() < 27 * nc {
            return Err("galerkin input sizes do not match the grid dims".to_string());
        }
        self.fits(fine_ac.len())?;
        self.fits(pc.len())?;

        let fine_buf = self.device.create_buffer_init(&BufferInitDescriptor {
            label: Some("PB mg galerkin fine op"),
            contents: bytemuck::cast_slice(fine_ac),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let pc_buf = self.device.create_buffer_init(&BufferInitDescriptor {
            label: Some("PB mg galerkin pc"),
            contents: bytemuck::cast_slice(pc),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let out = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("PB mg galerkin out"),
            size: (14 * nc * 8) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: true,
        });
        {
            let mut view = out.slice(..).get_mapped_range_mut();
            view.fill(0);
        }
        out.unmap();

        let dummy_ro = self.dummy.clone();
        let dummy_rw = self.dummy_rw.clone();
        let entries = self.bg_entries(
            &fine_buf,
            &dummy_ro,
            &dummy_ro,
            Some(&pc_buf),
            &dummy_rw,
            &out,
            0,
        );
        let bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("PB mg galerkin bind group"),
            layout: &self.bgl,
            entries: &entries,
        });
        let params = Params {
            n: nc as u32,
            nx: nxc as u32,
            ny: nyc as u32,
            nz: nzc as u32,
            nxf: nxf as u32,
            nyf: nyf as u32,
            nzf: nzf as u32,
            nc: 0,
            color: fine_numdia,
            iadjoint: 0,
            m: 0,
            pad0: 0,
            alpha: 0.0,
            pre1: 0.0,
        };
        let mut seg = self.begin_seg();
        let pipe = self.pipelines.galerkin.clone();
        self.dispatch(
            &mut seg,
            &pipe,
            &bg,
            params,
            (nc as u32).div_ceil(WORKGROUP),
        );
        self.flush(seg);

        let values = self.read_buffer(&out, 14 * nc)?;
        ac_out[..14 * nc].copy_from_slice(&values);
        Ok(())
    }

    /// Reads and sums one reduction region (0 or 1) of `count` partials.
    fn read_partials(&self, region: u32, count: u32) -> Result<f64, String> {
        if count == 0 {
            return Ok(0.0);
        }
        let offset = (region as u64) * self.partials_stride;
        let bytes = (count as u64) * 8;
        let mut enc = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("PB mg partials copy"),
            });
        enc.copy_buffer_to_buffer(&self.partials, offset, &self.partials_read, 0, bytes);
        self.queue.submit(Some(enc.finish()));
        let slice = self.partials_read.slice(0..bytes);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.device
            .poll(wgpu::PollType::Wait)
            .map_err(|e| format!("GPU poll failed: {}", e))?;
        rx.recv()
            .map_err(|_| "GPU map callback was dropped".to_string())?
            .map_err(|e| format!("GPU reduction mapping failed: {}", e))?;
        let view = slice.get_mapped_range();
        let data: &[f64] = bytemuck::cast_slice(&view);
        let sum: f64 = data.iter().sum();
        drop(view);
        self.partials_read.unmap();
        Ok(sum)
    }

    // -- kernel wrappers ----------------------------------------------------

    #[allow(clippy::too_many_arguments)]
    fn smooth(
        &mut self,
        seg: &mut Seg,
        level: &LevelBufs,
        fc: &wgpu::Buffer,
        x: &wgpu::Buffer,
        resid: &wgpu::Buffer,
        nu: i32,
        iadjoint: u32,
    ) {
        let bg = self.make_bg(level, fc, x, resid, 0);
        let p0 = grid_params(level);
        let n = level.n;
        if level.numdia == 4 {
            for _ in 0..nu {
                for color in 0..2u32 {
                    let mut p = p0;
                    p.color = color;
                    p.iadjoint = iadjoint;
                    let pipe = self.pipelines.gsrb7.clone();
            self.dispatch(seg, &pipe, &bg, p, n.div_ceil(WORKGROUP));
                }
            }
        } else {
            for _ in 0..nu {
                let order: Vec<u32> = if iadjoint != 0 {
                    (0..8u32).rev().collect()
                } else {
                    (0..8u32).collect()
                };
                for color in order {
                    let mut p = p0;
                    p.color = color;
                    let pipe = self.pipelines.gsrb27.clone();
            self.dispatch(seg, &pipe, &bg, p, n.div_ceil(WORKGROUP));
                }
            }
        }
    }

    fn resid(
        &mut self,
        seg: &mut Seg,
        level: &LevelBufs,
        fc: &wgpu::Buffer,
        x: &wgpu::Buffer,
        target: &wgpu::Buffer,
    ) {
        let bg = self.make_bg(level, fc, x, target, 0);
        let p = grid_params(level);
        let pipe = if level.numdia == 4 {
            self.pipelines.mresid7.clone()
        } else {
            self.pipelines.mresid27.clone()
        };
        self.dispatch(seg, &pipe, &bg, p, level.n.div_ceil(WORKGROUP));
    }

    fn zero_boundary(&mut self, seg: &mut Seg, level: &LevelBufs, x: &wgpu::Buffer) {
        let bg = self.make_bg(level, &level.fc, x, &self.dummy_rw.clone(), 0);
        let pipe = self.pipelines.zero_boundary.clone();
            self.dispatch(seg, &pipe,
            &bg,
            grid_params(level),
            level.n.div_ceil(WORKGROUP),
        );
    }

    fn fill_zero(&mut self, seg: &mut Seg, level: &LevelBufs, x: &wgpu::Buffer) {
        let bg = self.make_bg(level, &level.fc, x, &self.dummy_rw.clone(), 0);
        let pipe = self.pipelines.fill_zero.clone();
            self.dispatch(seg, &pipe,
            &bg,
            grid_params(level),
            level.n.div_ceil(WORKGROUP),
        );
    }

    #[allow(clippy::too_many_arguments)]
    fn restrc(
        &mut self,
        seg: &mut Seg,
        fine: &LevelBufs,
        coarse: &LevelBufs,
        x: &wgpu::Buffer,
        target: &wgpu::Buffer,
    ) {
        fine.pc.as_ref().expect("restriction needs pc weights");
        let bg = self.make_bg(fine, &fine.fc, x, target, 0);
        let mut p = grid_params(coarse);
        p.nxf = fine.nx;
        p.nyf = fine.ny;
        p.nzf = fine.nz;
        p.nc = coarse.n;
        let pipe = self.pipelines.restrc.clone();
            self.dispatch(seg, &pipe, &bg, p, coarse.n.div_ceil(WORKGROUP));
    }

    #[allow(clippy::too_many_arguments)]
    fn interp(
        &mut self,
        seg: &mut Seg,
        fine: &LevelBufs,
        coarse: &LevelBufs,
        x: &wgpu::Buffer,
        target: &wgpu::Buffer,
    ) {
        fine.pc.as_ref().expect("prolongation needs pc weights");
        let bg = self.make_bg(fine, &fine.fc, x, target, 0);
        let mut p = grid_params(fine);
        p.nxf = coarse.nx;
        p.nyf = coarse.ny;
        p.nzf = coarse.nz;
        p.nc = coarse.n;
        let pipe = self.pipelines.interp.clone();
            self.dispatch(seg, &pipe, &bg, p, fine.n.div_ceil(WORKGROUP));
    }

    fn matvec(
        &mut self,
        seg: &mut Seg,
        level: &LevelBufs,
        x: &wgpu::Buffer,
        target: &wgpu::Buffer,
    ) {
        let bg = self.make_bg(level, &level.fc, x, target, 0);
        let p = grid_params(level);
        let pipe = if level.numdia == 4 {
            self.pipelines.matvec7.clone()
        } else {
            self.pipelines.matvec27.clone()
        };
        self.dispatch(seg, &pipe, &bg, p, level.n.div_ceil(WORKGROUP));
    }

    /// Records both damping dot products of one coarse level into the two
    /// partials regions and reads the sums back.
    #[allow(clippy::too_many_arguments)]
    fn damping_dots(
        &mut self,
        level: &LevelBufs,
        xs: &wgpu::Buffer,
        w0: &wgpu::Buffer,
        w2: &wgpu::Buffer,
    ) -> Result<(f64, f64), String> {
        let m = interior_count(level.nx, level.ny, level.nz);
        if m == 0 {
            return Ok((0.0, 0.0));
        }
        let wg = m.div_ceil(RED_WORKGROUP);
        let mut seg = self.begin_seg();
        let bg_num = self.make_bg(level, &level.fc, xs, w0, 0);
        let bg_den = self.make_bg(level, &level.fc, xs, w2, self.partials_stride);
        let p = grid_params(level);
        self.dispatch(&mut seg, &self.pipelines.reduce_dot.clone(), &bg_num, p, wg);
        self.dispatch(&mut seg, &self.pipelines.reduce_dot.clone(), &bg_den, p, wg);
        self.flush(seg);
        let xnum = self.read_partials(0, wg)?;
        let xden = self.read_partials(1, wg)?;
        Ok((xnum, xden))
    }

    // -- the V-cycle --------------------------------------------------------

    fn v_cycle(&mut self, bufs: &[LevelBufs], levels: &[Level], nu1: i32, nu2: i32) -> Result<(), String> {
        let nlev = bufs.len();

        // Pre-smoothing on the fine level + residual into w1[0], w0[0] = w1[0].
        {
            let f = &bufs[0];
            let mut seg = self.begin_seg();
            self.smooth(&mut seg, f, &f.fc, &f.x, &f.w1, nu1, 0);
            self.resid(&mut seg, f, &f.fc, &f.x, &f.w1);
            self.flush(seg);
            let mut enc = self
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("PB mg w0 copy"),
                });
            enc.copy_buffer_to_buffer(&f.w1, 0, &f.w0, 0, (f.n * 8) as u64);
            self.queue.submit(Some(enc.finish()));
        }

        // Down-cycle: restrict the residual, pre-smooth coarser levels.
        for lev in 1..nlev {
            let prev = lev - 1;
            let (p_buf, c_buf) = (&bufs[prev], &bufs[lev]);
            {
                let mut seg = self.begin_seg();
                self.zero_boundary(&mut seg, p_buf, &p_buf.w1);
                self.restrc(&mut seg, p_buf, c_buf, &p_buf.w1, &c_buf.w0);
                self.flush(seg);
            }
            if lev != nlev - 1 {
                let mut seg = self.begin_seg();
                self.fill_zero(&mut seg, c_buf, &c_buf.x);
                // The smoother reads its right-hand side through the fc
                // binding, which here is the restricted residual w0[lev].
                self.smooth(&mut seg, c_buf, &c_buf.w0, &c_buf.x, &c_buf.w1, nu1, 0);
                self.resid(&mut seg, c_buf, &c_buf.w0, &c_buf.x, &c_buf.w1);
                self.flush(seg);
            }
        }

        // Coarsest level: direct banded solve on the CPU (dpbsl), then upload.
        {
            let last = nlev - 1;
            let l = &bufs[last];
            let banded = levels[last]
                .banded
                .as_ref()
                .ok_or("coarsest level is not factored for the direct solve")?;
            let nf = l.n as usize;
            let mut w = self.read_buffer(&l.w0, nf)?;
            crate::lapack::dpbsl(&banded.abd, banded.lda, banded.n, banded.m, &mut w);
            self.queue.write_buffer(&l.x, 0, bytemuck::cast_slice(&w));
            let mut seg = self.begin_seg();
            self.zero_boundary(&mut seg, l, &l.x);
            self.flush(seg);
        }

        // Up-cycle: interpolate the correction, damp it with the CG
        // steplength of the coarser level, correct and post-smooth.
        for lev in (0..nlev - 1).rev() {
            let coarser = lev + 1;
            let (f_buf, c_buf) = (&bufs[lev], &bufs[coarser]);
            let (xnum, xden);
            {
                let mut seg = self.begin_seg();
                self.interp(&mut seg, f_buf, c_buf, &c_buf.x, &f_buf.w1);
                self.matvec(&mut seg, c_buf, &c_buf.x, &c_buf.w2);
                self.flush(seg);
            }
            (xnum, xden) = self.damping_dots(c_buf, &c_buf.x, &c_buf.w0, &c_buf.w2)?;
            let xdamp = xnum / xden;

            let mut seg = self.begin_seg();
            {
                let bg = self.make_bg(f_buf, &f_buf.fc, &f_buf.x, &f_buf.w1, 0);
                let mut p = grid_params(f_buf);
                p.alpha = xdamp;
                self.dispatch(&mut seg, &self.pipelines.axpy.clone(), &bg, p, f_buf.n.div_ceil(WORKGROUP));
            }
            // Post-smoothing: the fine level uses its own right-hand side,
            // coarser levels the restricted residual.
            let rhs = if lev == 0 { &f_buf.fc } else { &f_buf.w0 };
            self.smooth(&mut seg, f_buf, rhs, &f_buf.x, &self.dummy_rw.clone(), nu2, 1);
            self.flush(seg);
        }

        Ok(())
    }

    /// GPU port of `mgcs::mvcs` (mgsolv == 1 path).
    #[allow(clippy::too_many_arguments)]
    fn mvcs(
        &mut self,
        u: &mut [f64],
        levels: &[Level],
        nu1: i32,
        nu2: i32,
        itmax: i32,
        errtol: f64,
    ) -> Result<i32, String> {
        let nlev = levels.len();
        let fine = &levels[0];
        let nf = fine.nf();

        let max_interior = levels
            .iter()
            .map(|l| interior_count(l.nx as u32, l.ny as u32, l.nz as u32))
            .max()
            .unwrap_or(1);
        self.ensure_partials(max_interior.div_ceil(RED_WORKGROUP));

        // ---- upload all levels ----
        let fits = |len: usize| -> Result<(), String> {
            let bytes = (len * 8) as u64;
            let lim = self.device.limits();
            if bytes > lim.max_buffer_size || bytes > (lim.max_storage_buffer_binding_size as u64) {
                return Err(format!("grid buffer of {} bytes exceeds GPU limits", bytes));
            }
            Ok(())
        };
        let mut bufs: Vec<LevelBufs> = Vec::with_capacity(nlev);
        for (lev, l) in levels.iter().enumerate() {
            fits(l.ac.len())?;
            fits(l.cc.len())?;
            fits(l.fc.len())?;
            let pc_buf = match &l.pc {
                Some(pc) => {
                    fits(pc.len())?;
                    Some(self.device.create_buffer_init(&BufferInitDescriptor {
                        label: Some("PB mg pc"),
                        contents: bytemuck::cast_slice(pc),
                        usage: wgpu::BufferUsages::STORAGE,
                    }))
                }
                None => None,
            };
            let ac = self.device.create_buffer_init(&BufferInitDescriptor {
                label: Some("PB mg ac"),
                contents: bytemuck::cast_slice(&l.ac),
                usage: wgpu::BufferUsages::STORAGE,
            });
            let cc = self.device.create_buffer_init(&BufferInitDescriptor {
                label: Some("PB mg cc"),
                contents: bytemuck::cast_slice(&l.cc),
                usage: wgpu::BufferUsages::STORAGE,
            });
            let fc = self.device.create_buffer_init(&BufferInitDescriptor {
                label: Some("PB mg fc"),
                contents: bytemuck::cast_slice(&l.fc),
                usage: wgpu::BufferUsages::STORAGE,
            });
            let lnf = l.nf();
            let zeros = vec![0.0f64; lnf];
            let mk = || {
                self.device.create_buffer_init(&BufferInitDescriptor {
                    label: Some("PB mg vec"),
                    contents: bytemuck::cast_slice(&zeros),
                    // COPY_SRC: the coarsest w0 is read back for the direct
                    // solve; COPY_DST: the solved coarsest x is written back
                    // with queue.write_buffer.
                    usage: wgpu::BufferUsages::STORAGE
                        | wgpu::BufferUsages::COPY_SRC
                        | wgpu::BufferUsages::COPY_DST,
                })
            };
            let x = if lev == 0 {
                fits(u.len().min(lnf))?;
                self.device.create_buffer_init(&BufferInitDescriptor {
                    label: Some("PB mg u"),
                    contents: bytemuck::cast_slice(&u[..lnf]),
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                })
            } else {
                mk()
            };
            bufs.push(LevelBufs {
                nx: l.nx as u32,
                ny: l.ny as u32,
                nz: l.nz as u32,
                n: lnf as u32,
                numdia: l.numdia as u32,
                ac,
                cc,
                fc,
                pc: pc_buf,
                x,
                w0: mk(),
                w1: mk(),
                w2: mk(),
            });
        }

        // rsden = L1 norm of the fine right-hand side (istop == 1).
        let mut rsden = {
            let f = &bufs[0];
            let m = interior_count(f.nx, f.ny, f.nz);
            if m == 0 {
                0.0
            } else {
                let wg = m.div_ceil(RED_WORKGROUP);
                let mut seg = self.begin_seg();
                // The reduce kernel only reads x; unused slots must not alias
                // it or each other (read-write usage is exclusive).
                let bg = self.make_bg(f, &self.dummy.clone(), &f.fc, &self.dummy_rw.clone(), 0);
                let mut p = grid_params(f);
                p.m = m;
                self.dispatch(&mut seg, &self.pipelines.reduce_abs.clone(), &bg, p, wg);
                self.flush(seg);
                self.read_partials(0, wg)?
            }
        };
        if rsden == 0.0 {
            rsden = 1.0;
        }

        let mut iters = 0i32;
        loop {
            self.v_cycle(&bufs, levels, nu1, nu2)?;

            // istop == 1 stopping test on the fine residual.
            iters += 1;
            let f = &bufs[0];
            let m = interior_count(f.nx, f.ny, f.nz);
            let rsnrm = {
                let mut seg = self.begin_seg();
                self.resid(&mut seg, f, &f.fc, &f.x, &f.w1);
                self.flush(seg);
                if m == 0 {
                    0.0
                } else {
                    let wg = m.div_ceil(RED_WORKGROUP);
                    let mut seg = self.begin_seg();
                    let bg = self.make_bg(f, &self.dummy.clone(), &f.w1, &self.dummy_rw.clone(), 0);
                    let mut p = grid_params(f);
                    p.m = m;
                    self.dispatch(&mut seg, &self.pipelines.reduce_abs.clone(), &bg, p, wg);
                    self.flush(seg);
                    self.read_partials(0, wg)?
                }
            };
            if !(iters < itmax && rsnrm / rsden > errtol) {
                break;
            }
        }

        // Read back the fine solution.
        let solution = self.read_buffer(&bufs[0].x, nf)?;
        u[..nf].copy_from_slice(&solution[..nf]);

        Ok(iters)
    }
}

fn grid_params(l: &LevelBufs) -> Params {
    Params {
        n: l.n,
        nx: l.nx,
        ny: l.ny,
        nz: l.nz,
        nxf: 0,
        nyf: 0,
        nzf: 0,
        nc: 0,
        color: 0,
        iadjoint: 0,
        m: interior_count(l.nx, l.ny, l.nz),
        pad0: 0,
        alpha: 0.0,
        pre1: 0.0,
    }
}

/// Maps a freshly submitted buffer and returns its first `bytes` as f64s.
fn map_read(device: &wgpu::Device, buffer: &wgpu::Buffer, bytes: usize) -> Result<Vec<f64>, String> {
    let slice = buffer.slice(0..bytes as u64);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device
        .poll(wgpu::PollType::Wait)
        .map_err(|e| format!("GPU poll failed: {}", e))?;
    rx.recv()
        .map_err(|_| "GPU map callback was dropped".to_string())?
        .map_err(|e| format!("GPU output mapping failed: {}", e))?;
    let view = slice.get_mapped_range();
    let data: Vec<f64> = bytemuck::cast_slice(&view).to_vec();
    drop(view);
    buffer.unmap();
    Ok(data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mgcs::{self, Level};

    /// Builds the same Galerkin hierarchy `linear_mgdriv` builds for a
    /// synthetic 65^3 constant-coefficient problem and factors the coarsest
    /// level, so both the CPU and GPU V-cycles can run on identical data.
    fn build_problem() -> (Vec<Level>, usize) {
        let (nx, ny, nz) = (65usize, 65usize, 65usize);
        let nlev = 5usize;
        let nf = nx * ny * nz;
        let axis = |n: usize| (0..n).map(|i| i as f64 / (n - 1) as f64).collect::<Vec<f64>>();
        let (xf, yf, zf) = (axis(nx), axis(ny), axis(nz));
        let gxcf = vec![0.0; 2 * ny * nz];
        let gycf = vec![0.0; 2 * nx * nz];
        let gzcf = vec![0.0; 2 * nx * ny];
        let a1 = vec![1.0; nf];
        let a2 = vec![1.0; nf];
        let a3 = vec![1.0; nf];
        let ccf = vec![0.1; nf];
        let fcf: Vec<f64> = (0..nf).map(|i| (((i * 7) % 23) as f64) / 23.0 - 0.5).collect();

        let mut ac0 = vec![0.0f64; 4 * nf];
        let mut cc0 = vec![0.0f64; nf];
        let mut fc0 = vec![0.0f64; nf];
        crate::build_a::build_a(
            nx, ny, nz, 0, 0, 4,
            &mut ac0, &mut cc0, &mut fc0,
            &xf, &yf, &zf, &gxcf, &gycf, &gzcf,
            &a1, &a2, &a3, &ccf, &fcf,
        );
        let mut levels = vec![Level {
            nx, ny, nz,
            ac: ac0,
            cc: cc0,
            fc: fc0,
            numdia: 4,
            pc: None,
            banded: None,
        }];
        let mut size = (nx, ny, nz);
        for lev in 1..nlev {
            let (pnx, pny, pnz) = size;
            let (cnx, cny, cnz) = crate::build_str::make_coarse(
                pnx as i32, pny as i32, pnz as i32,
            );
            let (cnx, cny, cnz) = (cnx as usize, cny as usize, cnz as usize);
            let npf = pnx * pny * pnz;
            let ncl = cnx * cny * cnz;
            let prev = &levels[lev - 1];
            let pc = crate::build_p::build_p_trilin_block(cnx, cny, cnz);
            let mut ac_c = vec![0.0f64; 14 * ncl];
            let fine_op = crate::build_g::FineOp::Seven(
                crate::matvec::split_bands4(&prev.ac, npf).expect("bands4"),
            );
            crate::build_g::build_galerkin(pnx, pny, pnz, cnx, cny, cnz, &pc, fine_op, &mut ac_c);
            let mut cc_c = vec![0.0f64; ncl];
            let mut fc_c = vec![0.0f64; ncl];
            crate::matvec::restrc(pnx, pny, pnz, cnx, cny, cnz, &prev.cc, &mut cc_c, &pc);
            crate::matvec::restrc(pnx, pny, pnz, cnx, cny, cnz, &prev.fc, &mut fc_c, &pc);
            levels[lev - 1].pc = Some(pc);
            levels.push(Level {
                nx: cnx, ny: cny, nz: cnz,
                ac: ac_c,
                cc: cc_c,
                fc: fc_c,
                numdia: 14,
                pc: None,
                banded: None,
            });
            size = (cnx, cny, cnz);
        }
        let last = nlev - 1;
        let l = &levels[last];
        let mut banded = crate::build_b::build_band(l.nx, l.ny, l.nz, &l.ac);
        let info = crate::lapack::dpbfa(&mut banded.abd, banded.lda, banded.n, banded.m);
        assert_eq!(info, 0, "synthetic coarsest factorization failed");
        levels[last].banded = Some(banded);
        (levels, nf)
    }

    /// Exhaustively checks the WGSL plane27 mapping against the CPU
    /// plane_index table by evaluating it on the device.
    #[test]
    fn wgsl_plane27_matches_cpu_table() {
        if GpuMg::new().is_err() {
            println!("No f64-capable GPU; skipping the plane27 probe.");
            return;
        }
        let s = solver().unwrap();
        let s = s.lock().unwrap();
        let shader = r#"
        fn plane27_probe(dx: i32, dy: i32, dz: i32) -> u32 {
            var layer: u32;
            if (dz == 0) {
                layer = 0u;
            } else if (dz > 0) {
                layer = 9u;
            } else {
                layer = 18u;
            }
            var row: u32;
            if (dy == 0) {
                row = select(select(4u, 3u, dx == 1), 0u, dx == 0);
            } else if (dy == 1) {
                row = select(select(6u, 5u, dx == 1), 1u, dx == 0);
            } else {
                row = select(select(8u, 7u, dx == 1), 2u, dx == 0);
            }
            return layer + row;
        }
        @group(0) @binding(0) var<storage, read_write> out: array<f64>;
        @compute @workgroup_size(1)
        fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
            let dz: i32 = i32(gid.x / 9u) - 1;
            let rem = gid.x % 9u;
            let dy: i32 = i32(rem / 3u) - 1;
            let dx: i32 = i32(rem % 3u) - 1;
            out[gid.x] = f64(plane27_probe(dx, dy, dz));
        }
        "#;
        let module = s.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("plane27 probe"),
            source: wgpu::ShaderSource::Wgsl(std::borrow::Cow::Borrowed(shader)),
        });
        let layout = s
            .device
            .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: None,
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });
        let pipeline = s
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: None,
                layout: Some(&s
                    .device
                    .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                        label: None,
                        bind_group_layouts: &[&layout],
                        push_constant_ranges: &[],
                    })),
                module: &module,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            });
        let buf = s.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: 27 * 8,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let bg = s.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: buf.as_entire_binding(),
            }],
        });
        let mut enc = s
            .device
            .create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(27, 1, 1);
        }
        s.queue.submit(Some(enc.finish()));
        let values = s.read_buffer(&buf, 27).unwrap();
        // The CPU plane_index table, replicated (build_g::plane_index is
        // private).
        let plane_index = |d: (i32, i32, i32)| -> usize {
            match d {
                (0, 0, 0) => 0, (0, 1, 0) => 1, (0, -1, 0) => 2,
                (1, 0, 0) => 3, (-1, 0, 0) => 4,
                (1, 1, 0) => 5, (-1, 1, 0) => 6,
                (1, -1, 0) => 7, (-1, -1, 0) => 8,
                (0, 0, 1) => 9, (0, 1, 1) => 10, (0, -1, 1) => 11,
                (1, 0, 1) => 12, (-1, 0, 1) => 13,
                (1, 1, 1) => 14, (-1, 1, 1) => 15,
                (1, -1, 1) => 16, (-1, -1, 1) => 17,
                (0, 0, -1) => 18, (0, 1, -1) => 19, (0, -1, -1) => 20,
                (1, 0, -1) => 21, (-1, 0, -1) => 22,
                (1, 1, -1) => 23, (-1, 1, -1) => 24,
                (1, -1, -1) => 25, (-1, -1, -1) => 26,
                _ => unreachable!(),
            }
        };
        for (gid, v) in values.iter().enumerate() {
            let dz = gid as i32 / 9 - 1;
            let rem = gid % 9;
            let dy = rem as i32 / 3 - 1;
            let dx = (rem % 3) as i32 - 1;
            let want = plane_index((dx, dy, dz));
            println!("({dx},{dy},{dz}) -> {} expect {}", *v as u32, want);
            assert_eq!(*v, want as f64, "plane27({dx},{dy},{dz}) wrong");
        }
    }

    /// The GPU Debye-Hueckel boundary must agree with the CPU point-major
    /// formulation (the same formula, plain double-precision exp).
    #[test]
    fn gpu_dh_boundary_matches_cpu() {
        if GpuMg::new().is_err() {
            println!("No f64-capable GPU found; skipping the DH boundary test.");
            return;
        }
        let mut rng_state = 987u64;
        let mut rand = move || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((rng_state >> 33) as f64) / ((u64::MAX >> 33) as f64)
        };
        let (nx, ny, nz) = (25usize, 23usize, 21usize);
        let (xf, yf, zf) = (
            (0..nx).map(|i| i as f64 * 0.5).collect::<Vec<f64>>(),
            (0..ny).map(|i| i as f64 * 0.5).collect::<Vec<f64>>(),
            (0..nz).map(|i| i as f64 * 0.5).collect::<Vec<f64>>(),
        );
        let atoms: Vec<[f64; 5]> = (0..40)
            .map(|_| {
                [
                    rand() * 10.0,
                    rand() * 10.0,
                    rand() * 10.0,
                    0.2 + rand(),
                    rand() * 0.3,
                ]
            })
            .collect();
        let xkappa = 1.29f64;
        let pre1 = 42.0f64;

        let mut cpu_x = vec![0.0f64; 2 * ny * nz];
        let mut cpu_y = vec![0.0f64; 2 * nx * nz];
        let mut cpu_z = vec![0.0f64; 2 * nx * ny];
        // CPU reference: the bcfl_mdh formulation.
        for (m, out) in cpu_x.iter_mut().enumerate() {
            let face = m / (ny * nz);
            let r = m % (ny * nz);
            let (k, j) = (r / ny, r % ny);
            let gpos = [
                if face == 0 { xf[0] } else { xf[nx - 1] },
                yf[j],
                zf[k],
            ];
            let mut acc = 0.0;
            for atom in &atoms {
                let (ax, ay, az, q, size) = (atom[0], atom[1], atom[2], atom[3], atom[4]);
                let dist = ((gpos[0] - ax).powi(2) + (gpos[1] - ay).powi(2) + (gpos[2] - az).powi(2)).sqrt();
                acc += pre1 * (q / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
            }
            *out = acc;
        }
        // (y and z faces follow the same formula; compare all three via GPU.)

        set_mg_gpu_requested(true);
        let (gpu_x, gpu_y, gpu_z) = dh_boundary_gpu(&atoms, &xf, &yf, &zf, xkappa, pre1)
            .expect("GPU DH boundary must run on this machine");
        set_mg_gpu_requested(false);

        let cmp = |gpu: &[f64], cpu: &[f64], what: &str| {
            let scale = cpu.iter().fold(1e-30f64, |m, v| m.max(v.abs()));
            let maxdiff = gpu
                .iter()
                .zip(cpu)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            // The GPU kernel uses the hardware f64 exp (~5e-8 relative
            // deviation from libm); the summed boundary values inherit that,
            // which is far below the multigrid stopping tolerance.
            assert!(
                maxdiff / scale < 5.0e-7,
                "{what}: GPU vs CPU DH deviation {maxdiff} at scale {scale}"
            );
        };
        cmp(&gpu_x, &cpu_x, "x-faces");
        // y/z faces: regenerate the CPU reference the same way.
        for (m, out) in cpu_y.iter_mut().enumerate() {
            let face = m / (nx * nz);
            let r = m % (nx * nz);
            let (k, i2) = (r / nx, r % nx);
            let gpos = [xf[i2], if face == 0 { yf[0] } else { yf[ny - 1] }, zf[k]];
            let mut acc = 0.0;
            for atom in &atoms {
                let (ax, ay, az, q, size) = (atom[0], atom[1], atom[2], atom[3], atom[4]);
                let dist = ((gpos[0] - ax).powi(2) + (gpos[1] - ay).powi(2) + (gpos[2] - az).powi(2)).sqrt();
                acc += pre1 * (q / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
            }
            *out = acc;
        }
        cmp(&gpu_y, &cpu_y, "y-faces");
        for (m, out) in cpu_z.iter_mut().enumerate() {
            let face = m / (nx * ny);
            let r = m % (nx * ny);
            let (j, i2) = (r / nx, r % nx);
            let gpos = [xf[i2], yf[j], if face == 0 { zf[0] } else { zf[nz - 1] }];
            let mut acc = 0.0;
            for atom in &atoms {
                let (ax, ay, az, q, size) = (atom[0], atom[1], atom[2], atom[3], atom[4]);
                let dist = ((gpos[0] - ax).powi(2) + (gpos[1] - ay).powi(2) + (gpos[2] - az).powi(2)).sqrt();
                acc += pre1 * (q / dist) * (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
            }
            *out = acc;
        }
        cmp(&gpu_z, &cpu_z, "z-faces");
    }

    /// The GPU V-cycle must reproduce the CPU solve to solver-tolerance level.
    #[test]
    fn gpu_mvcs_matches_cpu_mvcs() {
        if GpuMg::new().is_err() {
            println!("No f64-capable GPU found; skipping the multigrid GPU/CPU test.");
            return;
        }
        let (levels, nf) = build_problem();

        // CPU reference (the flag is irrelevant for the direct mgcs::mvcs call).
        let mut u_cpu = vec![0.0f64; nf];
        let it_cpu = mgcs::mvcs(&mut u_cpu, &levels, 2, 2, 200, 1.0e-6, 1, 2.2204460492503131e-16);

        // GPU solve.
        set_mg_gpu_requested(true);
        let mut u_gpu = vec![0.0f64; nf];
        let iters = mvcs_gpu(&mut u_gpu, &levels, 2, 2, 200, 1.0e-6, 1)
            .expect("GPU multigrid must run on this machine");
        set_mg_gpu_requested(false);
        assert!(iters >= 1);
        println!("multigrid iterations: cpu={it_cpu} gpu={iters}");

        let scale = u_cpu.iter().fold(1e-30f64, |m, v| m.max(v.abs()));
        let maxdiff = u_gpu
            .iter()
            .zip(&u_cpu)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        // Reduction-order differences on the GPU can shift the relative-L1
        // stopping test by one outer iteration; the solutions are only defined
        // to the solver tolerance (errtol), so 10x errtol is the agreement
        // bound.  Everything below that is far under MM-PBSA reporting
        // precision.
        assert!(
            maxdiff / scale < 1.0e-5,
            "GPU vs CPU multigrid deviation {maxdiff} at scale {scale} (cpu {it_cpu} iters, gpu {iters} iters)"
        );
    }

    /// Rough timing at a production-like grid size.  Run explicitly with
    /// `cargo test -p apbs-pmgc gpu_bench -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn gpu_bench_production_grid() {
        if GpuMg::new().is_err() {
            println!("No f64-capable GPU; skipping the benchmark.");
            return;
        }
        let (levels, nf) = build_problem_sized(129, 5);
        let mut u_cpu = vec![0.0f64; nf];
        let t0 = std::time::Instant::now();
        let it_cpu = mgcs::mvcs(&mut u_cpu, &levels, 2, 2, 200, 1.0e-6, 1, 2.2204460492503131e-16);
        let t_cpu = t0.elapsed();
        set_mg_gpu_requested(true);
        let mut u_gpu = vec![0.0f64; nf];
        let t0 = std::time::Instant::now();
        let it_gpu = mvcs_gpu(&mut u_gpu, &levels, 2, 2, 200, 1.0e-6, 1).expect("GPU solve");
        let t_gpu = t0.elapsed();
        set_mg_gpu_requested(false);
        println!(
            "129^3 solve: cpu {it_cpu} iters in {t_cpu:?}, gpu {it_gpu} iters in {t_gpu:?}"
        );
    }

    fn build_problem_sized(n: usize, nlev: usize) -> (Vec<Level>, usize) {
        let (nx, ny, nz) = (n, n, n);
        let nf = nx * ny * nz;
        let axis = |m: usize| (0..m).map(|i| i as f64 / (m - 1) as f64).collect::<Vec<f64>>();
        let (xf, yf, zf) = (axis(nx), axis(ny), axis(nz));
        let gxcf = vec![0.0; 2 * ny * nz];
        let gycf = vec![0.0; 2 * nx * nz];
        let gzcf = vec![0.0; 2 * nx * ny];
        let a1 = vec![1.0; nf];
        let a2 = vec![1.0; nf];
        let a3 = vec![1.0; nf];
        let ccf = vec![1.0; nf];
        let fcf: Vec<f64> = (0..nf).map(|i| (((i * 7) % 23) as f64) / 23.0 - 0.5).collect();
        let mut ac0 = vec![0.0f64; 4 * nf];
        let mut cc0 = vec![0.0f64; nf];
        let mut fc0 = vec![0.0f64; nf];
        crate::build_a::build_a(
            nx, ny, nz, 0, 0, 4,
            &mut ac0, &mut cc0, &mut fc0,
            &xf, &yf, &zf, &gxcf, &gycf, &gzcf,
            &a1, &a2, &a3, &ccf, &fcf,
        );
        let mut levels = vec![Level {
            nx, ny, nz,
            ac: ac0,
            cc: cc0,
            fc: fc0,
            numdia: 4,
            pc: None,
            banded: None,
        }];
        let mut size = (nx, ny, nz);
        for lev in 1..nlev {
            let (pnx, pny, pnz) = size;
            let (cnx, cny, cnz) = crate::build_str::make_coarse(pnx as i32, pny as i32, pnz as i32);
            let (cnx, cny, cnz) = (cnx as usize, cny as usize, cnz as usize);
            let npf = pnx * pny * pnz;
            let ncl = cnx * cny * cnz;
            let prev = &levels[lev - 1];
            let pc = crate::build_p::build_p_trilin_block(cnx, cny, cnz);
            let mut ac_c = vec![0.0f64; 14 * ncl];
            let fine_op = crate::build_g::FineOp::Seven(
                crate::matvec::split_bands4(&prev.ac, npf).expect("bands4"),
            );
            crate::build_g::build_galerkin(pnx, pny, pnz, cnx, cny, cnz, &pc, fine_op, &mut ac_c);
            let mut cc_c = vec![0.0f64; ncl];
            let mut fc_c = vec![0.0f64; ncl];
            crate::matvec::restrc(pnx, pny, pnz, cnx, cny, cnz, &prev.cc, &mut cc_c, &pc);
            crate::matvec::restrc(pnx, pny, pnz, cnx, cny, cnz, &prev.fc, &mut fc_c, &pc);
            levels[lev - 1].pc = Some(pc);
            levels.push(Level {
                nx: cnx, ny: cny, nz: cnz,
                ac: ac_c,
                cc: cc_c,
                fc: fc_c,
                numdia: 14,
                pc: None,
                banded: None,
            });
            size = (cnx, cny, cnz);
        }
        let last = nlev - 1;
        let l = &levels[last];
        let mut banded = crate::build_b::build_band(l.nx, l.ny, l.nz, &l.ac);
        let info = crate::lapack::dpbfa(&mut banded.abd, banded.lda, banded.n, banded.m);
        assert_eq!(info, 0);
        levels[last].banded = Some(banded);
        (levels, nf)
    }

    /// The GPU Galerkin assembly must reproduce the CPU product for both
    /// fine-operator flavours (7-point from level 0, 27-point Galerkin
    /// levels) on a real trilinear pc.
    #[test]
    fn gpu_galerkin_matches_cpu() {
        if GpuMg::new().is_err() {
            println!("No f64-capable GPU found; skipping the Galerkin GPU/CPU test.");
            return;
        }
        let mut rng_state = 12345u64;
        let mut rand = move || {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((rng_state >> 33) as f64) / (u64::MAX >> 33) as f64
        };

        // First a hand-checkable constant-operator case: the CPU values for
        // 9^3 -> 5^3 are diag 6.75, single-axis bands 0.375, two-axis
        // 0.3125, three-axis 0.09375 at the fully interior point.
        {
            let (nxf, nyf, nzf) = (9usize, 9usize, 9usize);
            let (nxc, nyc, nzc) = (5usize, 5usize, 5usize);
            let nf = nxf * nyf * nzf;
            let nc = nxc * nyc * nzc;
            let fine_ac: Vec<f64> = (0..4 * nf)
                .map(|i| if i < nf { 6.0 } else { 1.0 })
                .collect();
            let pc = crate::build_p::build_p_trilin_block(nxc, nyc, nzc);
            let mut cpu = vec![0.0f64; 14 * nc];
            crate::build_g::build_galerkin(
                nxf, nyf, nzf, nxc, nyc, nzc, &pc,
                crate::build_g::FineOp::Seven(
                    crate::matvec::split_bands4(&fine_ac, nf).expect("bands4"),
                ),
                &mut cpu,
            );
            set_mg_gpu_requested(true);
            let mut gpu = vec![0.0f64; 14 * nc];
            assert!(galerkin_gpu(&fine_ac, 4, &pc, (nxf, nyf, nzf), (nxc, nyc, nzc), &mut gpu));
            set_mg_gpu_requested(false);
            let lin = 2 + 2 * nxc + 2 * nxc * nyc;
            for band in 0..4 {
                println!(
                    "constant case band {band}: gpu={} cpu={}",
                    gpu[band * nc + lin], cpu[band * nc + lin]
                );
            }
            let maxdiff = gpu.iter().zip(&cpu).map(|(a, b)| (a - b).abs()).fold(0.0f64, f64::max);
            println!("constant case maxdiff = {maxdiff}");
        }

        for (nf_axis, fine_numdia) in [(33usize, 4u32), (17usize, 14u32)] {
            let (nxf, nyf, nzf) = (nf_axis, nf_axis, nf_axis);
            let (nxc, nyc, nzc) = crate::build_str::make_coarse(
                nxf as i32, nyf as i32, nzf as i32,
            );
            let (nxc, nyc, nzc) = (nxc as usize, nyc as usize, nzc as usize);
            let nf = nxf * nyf * nzf;
            let nc = nxc * nyc * nzc;
            let nb = fine_numdia as usize;

            // Random positive-convention bands plus a strong diagonal.
            let fine_ac: Vec<f64> = (0..nb * nf)
                .map(|i| if i < nf { 6.0 + rand() } else { rand() })
                .collect();
            let pc = crate::build_p::build_p_trilin_block(nxc, nyc, nzc);

            let mut cpu = vec![0.0f64; 14 * nc];
            let fine_op = if fine_numdia == 4 {
                crate::build_g::FineOp::Seven(
                    crate::matvec::split_bands4(&fine_ac, nf).expect("bands4"),
                )
            } else {
                crate::build_g::FineOp::TwentySeven(
                    crate::matvec::split_bands14(&fine_ac, nf).expect("bands14"),
                )
            };
            crate::build_g::build_galerkin(
                nxf, nyf, nzf, nxc, nyc, nzc, &pc, fine_op, &mut cpu,
            );

            set_mg_gpu_requested(true);
            let mut gpu = vec![0.0f64; 14 * nc];
            assert!(
                galerkin_gpu(&fine_ac, fine_numdia, &pc, (nxf, nyf, nzf), (nxc, nyc, nzc), &mut gpu),
                "GPU Galerkin must run on this machine"
            );
            set_mg_gpu_requested(false);

            let scale = cpu.iter().fold(1e-30f64, |m, v| m.max(v.abs()));
            let maxdiff = gpu
                .iter()
                .zip(&cpu)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max);
            let mut shown = 0;
            for (band, (g, cval)) in gpu.iter().zip(&cpu).enumerate() {
                if shown < 8 && (g - cval).abs() > 1e-9 * scale.max(1.0) {
                    let (bi, idx) = (band / nc, band % nc);
                    let k = idx % nzc;
                    let j = (idx / nzc) % nyc;
                    let i = idx / (nzc * nyc);
                    println!(
                        "mismatch band {bi} at ({i},{j},{k}): gpu={g} cpu={cval} interior={}",
                        i > 0 && i < nxc - 1 && j > 0 && j < nyc - 1 && k > 0 && k < nzc - 1
                    );
                    shown += 1;
                }
            }
            assert!(
                maxdiff / scale < 1.0e-12,
                "Galerkin GPU vs CPU deviation {maxdiff} at scale {scale} (numdia {fine_numdia})"
            );
        }
    }

    /// An unusual but valid case: no requested GPU must leave `u` untouched
    /// and return None.
    #[test]
    fn mvcs_gpu_respects_the_request_flag() {
        let (levels, nf) = build_problem();
        set_mg_gpu_requested(false);
        let mut u = vec![7.5f64; nf];
        assert!(mvcs_gpu(&mut u, &levels, 2, 2, 60, 1.0e-6, 1).is_none());
        assert!(u.iter().all(|&v| v == 7.5));
    }
}
