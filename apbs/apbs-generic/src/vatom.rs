// APBS vatom.rs - Atom data structure
// Port of src/generic/vatom.h / vatom.c

use std::fmt;

/// Maximum residue/atom name length
pub const VMAX_RECLEN: usize = 64;

/// Atom data structure for PDB/PQR interface
#[derive(Debug, Clone)]
pub struct Vatom {
    /// Atomic coordinates (Angstroms)
    pub position: [f64; 3],
    /// Atomic radius (Angstroms)
    pub radius: f64,
    /// Partial charge (elementary charges)
    pub charge: f64,
    /// Partition ID (for parallel decomposition)
    pub part_id: f64,
    /// VdW well depth for WCA calculations (kJ/mol)
    pub epsilon: f64,
    /// Unique non-negative atom index
    pub id: i32,
    /// Residue name from PDB/PQR
    pub res_name: String,
    /// Atom name from PDB/PQR
    pub atom_name: String,
}

impl Default for Vatom {
    fn default() -> Self {
        Self {
            position: [0.0; 3],
            radius: 0.0,
            charge: 0.0,
            part_id: -1.0,
            epsilon: 0.0,
            id: 0,
            res_name: String::new(),
            atom_name: String::new(),
        }
    }
}

impl Vatom {
    /// Create a new atom with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Get position as slice
    #[inline]
    pub fn position(&self) -> &[f64; 3] {
        &self.position
    }

    /// Set position from array
    #[inline]
    pub fn set_position(&mut self, pos: [f64; 3]) {
        self.position = pos;
    }

    #[inline]
    pub fn radius(&self) -> f64 {
        self.radius
    }

    #[inline]
    pub fn set_radius(&mut self, r: f64) {
        self.radius = r;
    }

    #[inline]
    pub fn charge(&self) -> f64 {
        self.charge
    }

    #[inline]
    pub fn set_charge(&mut self, q: f64) {
        self.charge = q;
    }

    #[inline]
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    #[inline]
    pub fn set_epsilon(&mut self, e: f64) {
        self.epsilon = e;
    }

    #[inline]
    pub fn part_id(&self) -> f64 {
        self.part_id
    }

    #[inline]
    pub fn set_part_id(&mut self, id: i32) {
        self.part_id = id as f64;
    }

    #[inline]
    pub fn atom_id(&self) -> i32 {
        self.id
    }

    #[inline]
    pub fn set_atom_id(&mut self, id: i32) {
        self.id = id;
    }

    #[inline]
    pub fn res_name(&self) -> &str {
        &self.res_name
    }

    pub fn set_res_name(&mut self, name: &str) {
        self.res_name = name.to_string();
    }

    #[inline]
    pub fn atom_name(&self) -> &str {
        &self.atom_name
    }

    pub fn set_atom_name(&mut self, name: &str) {
        self.atom_name = name.to_string();
    }

    /// Copy all fields from another atom
    pub fn copy_from(&mut self, src: &Vatom) {
        self.position = src.position;
        self.radius = src.radius;
        self.charge = src.charge;
        self.part_id = src.part_id;
        self.epsilon = src.epsilon;
        self.id = src.id;
        self.res_name = src.res_name.clone();
        self.atom_name = src.atom_name.clone();
    }

    /// Return memory usage in bytes
    pub fn mem_chk(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.res_name.len()
            + self.atom_name.len()
    }
}

impl fmt::Display for Vatom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Atom {} ({}) {}: q={:.4} r={:.4} [{:.3}, {:.3}, {:.3}]",
            self.id,
            self.atom_name,
            self.res_name,
            self.charge,
            self.radius,
            self.position[0],
            self.position[1],
            self.position[2],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_creation() {
        let mut atom = Vatom::new();
        atom.set_position([1.0, 2.0, 3.0]);
        atom.set_charge(-0.5);
        atom.set_radius(1.5);
        atom.set_atom_name("CA");
        atom.set_res_name("ALA");

        assert_eq!(atom.position(), &[1.0, 2.0, 3.0]);
        assert_eq!(atom.charge(), -0.5);
        assert_eq!(atom.radius(), 1.5);
        assert_eq!(atom.atom_name(), "CA");
        assert_eq!(atom.res_name(), "ALA");
    }

    #[test]
    fn test_atom_copy() {
        let mut a1 = Vatom::new();
        a1.set_charge(1.5);
        a1.set_atom_name("CB");

        let mut a2 = Vatom::new();
        a2.copy_from(&a1);

        assert_eq!(a2.charge(), 1.5);
        assert_eq!(a2.atom_name(), "CB");
    }
}
