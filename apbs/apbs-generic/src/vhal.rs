// APBS vhal.rs - Global enums, constants, and macros
// Port of src/generic/vhal.h

/// Spatial dimension (3D)
pub const VAPBS_DIM: usize = 3;

/// Maximum number of molecules in a PBE calculation
pub const MAXMOL: usize = 5;

/// Maximum number of ion species
pub const MAXION: usize = 10;

/// Maximum focusing levels
pub const MAXFOCUS: usize = 5;

/// Minimum multigrid levels
pub const VMGNLEV: usize = 4;

/// Maximum grid spacing reduction per focusing level
pub const VREDFRAC: f64 = 0.25;

/// Vertices per simplex (3D tetrahedra)
pub const VAPBS_NVS: usize = 4;

/// Face definitions for volume
pub const VAPBS_RIGHT: usize = 0;
pub const VAPBS_FRONT: usize = 1;
pub const VAPBS_LEFT: usize = 2;
pub const VAPBS_BACK: usize = 3;
pub const VAPBS_UP: usize = 4;
pub const VAPBS_DOWN: usize = 5;

/// Small number for grid comparisons
pub const VPMGSMALL: f64 = 1e-12;

/// Sinh chopping bounds
pub const SINH_MIN: f64 = -85.0;
pub const SINH_MAX: f64 = 85.0;

/// Timer IDs
pub const APBS_TIMER_WALL_CLOCK: usize = 26;
pub const APBS_TIMER_SETUP: usize = 27;
pub const APBS_TIMER_SOLVER: usize = 28;
pub const APBS_TIMER_ENERGY: usize = 29;
pub const APBS_TIMER_FORCE: usize = 30;

/// Maximum hash dimension for cell lists
pub const MAX_HASH_DIM: usize = 50;

/// Very small number
pub const VSMALL: f64 = 1.0e-14;

/// Very large number
pub const VLARGE: f64 = 1.0e14;

/// Max value for SMPBE
pub const VSMPBE_MAX: f64 = 100.0;

/// Return codes
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VrcCode {
    Warning = -1,
    Failure = 0,
    Success = 1,
}

impl VrcCode {
    pub fn is_success(self) -> bool {
        self as i32 > 0
    }
}

/// Solver method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VsolMeth {
    CGMG = 0,
    Newton = 1,
    MG = 2,
    CG = 3,
    SOR = 4,
    RBGS = 5,
    WJ = 6,
    Richardson = 7,
    CGMGAqua = 8,
    NewtonAqua = 9,
}

/// Surface method for molecular surface definition
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VsurfMeth {
    Mol = 0,
    MolSmooth = 1,
    Spline = 2,
    Spline3 = 3,
    Spline4 = 4,
    Sacc = 5,
}

/// PBE type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VhalPBEType {
    LPBE = 0,
    NPBE = 1,
    LRPBE = 2,
    NRPBE = 3,
    SMPBE = 4,
}

/// MG nonlinearity key
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VhalIPKEYType {
    SMPBE = -2,
    LPBE = -1,
    NPBE = 0,
}

/// Nonlinear type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VhalNONLINType {
    LPBE = 0,
    NPBE = 1,
    SMPBE = 2,
    LPBEAQUA = 3,
    NPBEAQUA = 4,
}

/// Boundary condition method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Vbcfl {
    Zero = 0,
    SDH = 1,
    MDH = 2,
    Unused = 3,
    Focus = 4,
    Mem = 5,
    Map = 6,
}

/// Charge discretization method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VchrgMeth {
    Tril = 0,
    Bspl2 = 1,
    Bspl4 = 2,
}

/// Charge source
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VchrgSrc {
    Charge = 0,
    Permanent = 1,
    Induced = 2,
    NLInduced = 3,
}

/// Output data type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VdataType {
    Charge = 0,
    Pot = 1,
    AtomPot = 2,
    Smol = 3,
    Sspl = 4,
    Vdw = 5,
    Ivdw = 6,
    Lap = 7,
    Edens = 8,
    Ndens = 9,
    Qdens = 10,
    DielX = 11,
    DielY = 12,
    DielZ = 13,
    Kappa = 14,
}

/// Output file format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VdataFormat {
    DX = 0,
    UHBD = 1,
    AVS = 2,
    MCSF = 3,
    GZ = 4,
    Flat = 5,
    DXBin = 6,
}

/// Output format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum VoutputFormat {
    Null = 0,
    Flat = 1,
}

/// 3D index macro: IJK(i,j,k) -> flat index
/// Corresponds to: ((k)*(nx)*(ny))+((j)*(nx))+(i)
#[inline(always)]
pub fn ijk(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
    k * nx * ny + j * nx + i
}

/// IJKx variant: ((i)*(ny)*(nz))+((k)*(ny))+(j)
#[inline(always)]
pub fn ijkx(j: usize, k: usize, i: usize, ny: usize, nz: usize) -> usize {
    i * ny * nz + k * ny + j
}

/// IJKy variant: ((j)*(nx)*(nz))+((k)*(nx))+(i)
#[inline(always)]
pub fn ijky(i: usize, k: usize, j: usize, nx: usize, nz: usize) -> usize {
    j * nx * nz + k * nx + i
}

/// IJKz variant: same as IJK
#[inline(always)]
pub fn ijkz(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
    k * nx * ny + j * nx + i
}

/// VFCHI macro for spline interface
#[inline(always)]
pub fn vfchi(iint: f64, iflt: f64) -> f64 {
    1.5 + iint - iflt
}

/// Cube
#[inline(always)]
pub fn vcub(x: f64) -> f64 {
    x * x * x
}
