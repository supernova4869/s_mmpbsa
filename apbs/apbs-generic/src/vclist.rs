// APBS vclist.rs - Cell list for spatial hashing
// Port of src/generic/vclist.h / vclist.c

use crate::error::{ApbsError, ApbsResult};
use crate::valist::Valist;
use crate::vhal::VAPBS_DIM;

/// Domain mode for cell list construction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VclistDomainMode {
    /// Auto encompass entire molecule
    Auto,
    /// User-specified domain
    Manual,
}

/// A cell in the spatial hash grid containing atom indices
#[derive(Debug, Clone)]
pub struct VclistCell {
    /// Indices of atoms in this cell
    pub atoms: Vec<usize>,
}

impl VclistCell {
    pub fn new() -> Self {
        Self { atoms: Vec::new() }
    }

    pub fn with_capacity(n: usize) -> Self {
        Self { atoms: Vec::with_capacity(n) }
    }

    /// Number of atoms in this cell
    pub fn natoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get atom index by local index
    pub fn get_atom_index(&self, i: usize) -> usize {
        self.atoms[i]
    }
}

/// Atom cell list for spatial lookups
#[derive(Debug)]
pub struct Vclist {
    /// Reference to the atom list (borrowed)
    // Note: In C this is a pointer; in Rust we store the data differently
    /// Grid dimensions (nx, ny, nz)
    pub npts: [usize; 3],
    /// Total cells = nx*ny*nz
    pub n: usize,
    /// Max probe radius
    pub max_radius: f64,
    /// Cell array
    pub cells: Vec<VclistCell>,
    /// Grid lower corner
    pub lower_corner: [f64; 3],
    /// Grid upper corner
    pub upper_corner: [f64; 3],
    /// Grid spacings
    pub spacs: [f64; 3],
    /// Domain mode
    pub mode: VclistDomainMode,
}

impl Vclist {
    /// Create a cell list automatically encompassing the molecule
    pub fn new_auto(
        alist: &Valist,
        max_radius: f64,
        npts: [usize; 3],
    ) -> ApbsResult<Self> {
        // Compute bounding box
        let mut lower = [f64::MAX; 3];
        let mut upper = [f64::MIN; 3];
        let mut max_rad = 0.0_f64;
        let mut found_any = false;

        for atom in alist.atoms() {
            found_any = true;
            for d in 0..VAPBS_DIM {
                if atom.position[d] < lower[d] {
                    lower[d] = atom.position[d];
                }
                if atom.position[d] > upper[d] {
                    upper[d] = atom.position[d];
                }
            }
            if atom.radius > max_rad {
                max_rad = atom.radius;
            }
        }

        if !found_any {
            // Default bounding box for empty atom list
            lower = [-1.0; 3];
            upper = [1.0; 3];
        }

        // Expand bounding box by sqrt(2) * (max_rad + max_radius)
        let expand = std::f64::consts::SQRT_2 * (max_rad + max_radius);
        for d in 0..VAPBS_DIM {
            lower[d] -= expand;
            upper[d] += expand;
        }

        Self::new_manual(alist, max_radius, npts, lower, upper)
    }

    /// Create a cell list with manually specified domain
    pub fn new_manual(
        alist: &Valist,
        max_radius: f64,
        npts: [usize; 3],
        lower_corner: [f64; 3],
        upper_corner: [f64; 3],
    ) -> ApbsResult<Self> {
        // Validate
        for d in 0..VAPBS_DIM {
            if npts[d] < 3 {
                return Err(ApbsError::InvalidParameter(format!(
                    "npts[{}] = {} must be >= 3",
                    d, npts[d]
                )));
            }
            if upper_corner[d] <= lower_corner[d] {
                return Err(ApbsError::InvalidParameter(format!(
                    "upper_corner[{}] = {} must be > lower_corner[{}] = {}",
                    d, upper_corner[d], d, lower_corner[d]
                )));
            }
        }

        let n = npts[0] * npts[1] * npts[2];

        let mut spacs = [0.0f64; 3];
        for d in 0..VAPBS_DIM {
            spacs[d] = (upper_corner[d] - lower_corner[d]) / (npts[d] - 1) as f64;
        }

        // Create empty cells
        let mut cells: Vec<VclistCell> = Vec::with_capacity(n);
        for _ in 0..n {
            cells.push(VclistCell::new());
        }

        // Two-pass atom assignment: count then assign
        // Pass 1: count atoms per cell
        for atom in alist.atoms() {
            let (imin, imax) = Self::grid_span(
                atom.position, atom.radius, max_radius,
                lower_corner, spacs, npts,
            );
            for i in imin[0]..=imax[0] {
                for j in imin[1]..=imax[1] {
                    for k in imin[2]..=imax[2] {
                        let idx = Self::array_index(i, j, k, npts);
                        cells[idx].atoms.push(0); // placeholder
                    }
                }
            }
        }

        // Pass 2: assign atom indices
        // Reset cells
        for cell in &mut cells {
            cell.atoms.clear();
        }

        for (atom_idx, atom) in alist.atoms().iter().enumerate() {
            let (imin, imax) = Self::grid_span(
                atom.position, atom.radius, max_radius,
                lower_corner, spacs, npts,
            );
            for i in imin[0]..=imax[0] {
                for j in imin[1]..=imax[1] {
                    for k in imin[2]..=imax[2] {
                        let idx = Self::array_index(i, j, k, npts);
                        cells[idx].atoms.push(atom_idx);
                    }
                }
            }
        }

        Ok(Self {
            npts,
            n,
            max_radius,
            cells,
            lower_corner,
            upper_corner,
            spacs,
            mode: VclistDomainMode::Manual,
        })
    }

    /// Get cell for a given position
    pub fn get_cell(&self, pos: [f64; 3]) -> Option<&VclistCell> {
        let mut ic = [0usize; VAPBS_DIM];
        for d in 0..VAPBS_DIM {
            let idx = (pos[d] - self.lower_corner[d]) / self.spacs[d];
            if idx < 0.0 || idx >= self.npts[d] as f64 {
                return None;
            }
            ic[d] = idx as usize;
        }

        let array_idx = Self::array_index(ic[0], ic[1], ic[2], self.npts);
        Some(&self.cells[array_idx])
    }

    /// Get max probe radius
    pub fn max_radius(&self) -> f64 {
        self.max_radius
    }

    /// Compute min/max grid indices for an atom
    fn grid_span(
        pos: [f64; 3],
        atom_radius: f64,
        max_radius: f64,
        lower_corner: [f64; 3],
        spacs: [f64; 3],
        npts: [usize; 3],
    ) -> ([usize; 3], [usize; 3]) {
        let mut imin = [0usize; VAPBS_DIM];
        let mut imax = [0usize; VAPBS_DIM];

        for d in 0..VAPBS_DIM {
            let r = atom_radius + max_radius;
            let min_val = (pos[d] - r - lower_corner[d]) / spacs[d];
            let max_val = (pos[d] + r - lower_corner[d]) / spacs[d];

            imin[d] = if min_val < 0.0 { 0 } else { min_val as usize };
            imax[d] = ((max_val + 1.0) as usize).min(npts[d] - 1);
        }

        (imin, imax)
    }

    /// 3D to 1D array index
    #[inline]
    fn array_index(i: usize, j: usize, k: usize, npts: [usize; 3]) -> usize {
        npts[1] * npts[2] * i + npts[2] * j + k
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vatom::Vatom;

    fn make_test_valist() -> Valist {
        let mut valist = Valist::new();
        let mut a1 = Vatom::new();
        a1.set_position([0.0, 0.0, 0.0]);
        a1.set_radius(1.5);
        let mut a2 = Vatom::new();
        a2.set_position([3.0, 0.0, 0.0]);
        a2.set_radius(1.5);
        valist.atoms.push(a1);
        valist.atoms.push(a2);
        valist.number = 2;
        valist
    }

    #[test]
    fn test_vclist_auto() {
        let valist = make_test_valist();
        let clist = Vclist::new_auto(&valist, 1.4, [10, 10, 10]).unwrap();
        assert_eq!(clist.npts, [10, 10, 10]);
        assert_eq!(clist.n, 1000);
    }

    #[test]
    fn test_vclist_get_cell() {
        let valist = make_test_valist();
        let clist = Vclist::new_auto(&valist, 1.4, [10, 10, 10]).unwrap();

        // Origin should be in some cell
        let cell = clist.get_cell([0.0, 0.0, 0.0]);
        assert!(cell.is_some());
    }
}
