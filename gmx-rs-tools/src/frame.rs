//! Shared trajectory data model.

/// A three component coordinate vector.
pub type Rvec = [f32; 3];
/// A 3x3 boxm matrix, stored as boxm[vector][component] like GROMACS does.
pub type Matrix = [[f32; 3]; 3];

/// Periodic boundary condition type, mirroring `PbcType` from
/// `gromacs/pbcutil/pbcenums.h`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PbcType {
    Unset,
    No,
    Xyz,
    XY,
    Screw,
}

impl PbcType {
    /// Values must match `gmx::PbcType` from `md_enums.h`:
    /// `Xyz = 0, No = 1, XY = 2, Screw = 3, Unset = 4`.
    pub fn from_int(v: i32) -> PbcType {
        match v {
            0 => PbcType::Xyz,
            1 => PbcType::No,
            2 => PbcType::XY,
            3 => PbcType::Screw,
            _ => PbcType::Unset,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            PbcType::Xyz => "xyz",
            PbcType::No => "no",
            PbcType::XY => "xy",
            PbcType::Screw => "screw",
            PbcType::Unset => "Unset",
        }
    }
}

/// A single atom as needed by the tools implemented here.
#[derive(Debug, Clone)]
pub struct Atom {
    pub name: String,
    pub atom_type: String,
    /// B-state atom type (`atomtypeB` in the tpr).
    pub atom_type_b: String,
    pub resind: i32,
    pub mass: f64,
    pub charge: f64,
    pub mass_b: f64,
    pub charge_b: f64,
    pub ptype: i32,
    pub atomnumber: i32,
    pub elem: String,
    pub type_id: u16,
    pub type_id_b: u16,
}

/// Residue information, mirroring `t_resinfo`.
#[derive(Debug, Clone)]
pub struct ResInfo {
    pub name: String,
    pub nr: i32,
    pub ic: u8,
    pub chainid: u8,
}

/// Minimal topology information (`t_atoms` plus the molecule bookkeeping that
/// `make_ndx` and `trjconv` need).
#[derive(Debug, Clone, Default)]
pub struct Atoms {
    pub atom: Vec<Atom>,
    pub resinfo: Vec<ResInfo>,
    pub name: String,
}

impl Atoms {
    pub fn nr(&self) -> usize {
        self.atom.len()
    }

    pub fn nres(&self) -> usize {
        self.resinfo.len()
    }
}

/// One frame of a trajectory.
#[derive(Debug, Clone, Default)]
pub struct Frame {
    pub natoms: usize,
    pub step: Option<i64>,
    pub time: Option<f64>,
    pub lambda: Option<f32>,
    pub boxm: Option<Matrix>,
    pub x: Option<Vec<Rvec>>,
    pub v: Option<Vec<Rvec>>,
    pub f: Option<Vec<Rvec>>,
    /// XTC compression precision (nm).
    pub prec: Option<f32>,
    pub pbc_type: PbcType,
    /// Per atom topology, only present for structure files.
    pub atoms: Option<Atoms>,
    pub title: String,
}

impl Default for PbcType {
    fn default() -> Self {
        PbcType::Unset
    }
}

impl Frame {
    pub fn new(natoms: usize) -> Frame {
        Frame {
            natoms,
            ..Default::default()
        }
    }
}

/// Whether a boxm matrix is triclinic (any off-diagonal element set).
pub fn is_triclinic(boxm: &Matrix) -> bool {
    boxm[0][1] != 0.0 || boxm[0][2] != 0.0 || boxm[1][0] != 0.0 || boxm[1][2] != 0.0
        || boxm[2][0] != 0.0 || boxm[2][1] != 0.0
}
