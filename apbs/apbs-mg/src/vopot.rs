// APBS vopot.rs - Potential oracle
// Port of src/mg/vopot.h / vopot.c

use crate::vmgrid::Vmgrid;
use apbs_generic::error::{ApbsError, ApbsResult};
use apbs_generic::vhal::Vbcfl;
use apbs_generic::vpbe::Vpbe;
use apbs_generic::vunit::{EC, EPS0, KB};
use std::sync::Arc;

/// Potential oracle - provides potential at any point
pub struct Vopot {
    pub mgrid: Vmgrid,
    pub pbe: Option<Arc<Vpbe>>,
    pub bcfl: Vbcfl,
}

impl Vopot {
    pub fn new(mgrid: Vmgrid, bcfl: Vbcfl) -> Self {
        Self {
            mgrid,
            pbe: None,
            bcfl,
        }
    }

    pub fn with_pbe(mgrid: Vmgrid, pbe: Arc<Vpbe>, bcfl: Vbcfl) -> Self {
        Self {
            mgrid,
            pbe: Some(pbe),
            bcfl,
        }
    }

    fn sdh_offgrid_potential(&self, pbe: &Vpbe, pt: [f64; 3]) -> ApbsResult<f64> {
        let eps_w = pbe.solvent_diel;
        let xkappa = pbe.xkappa;
        let t = pbe.temperature;
        let size = pbe.solute_radius;
        let center = pbe.solute_center;
        let charge = EC * pbe.solute_charge;

        let dx = center[0] - pt[0];
        let dy = center[1] - pt[1];
        let dz = center[2] - pt[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-8);
        let mut val = charge / (4.0 * std::f64::consts::PI * EPS0 * eps_w * dist);
        if xkappa.abs() > 0.0 {
            val *= (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
        }
        Ok(val * EC / (KB * t))
    }

    fn mdh_offgrid_potential(&self, pbe: &Vpbe, pt: [f64; 3]) -> f64 {
        let eps_w = pbe.solvent_diel;
        let xkappa = pbe.xkappa;
        let t = pbe.temperature;
        let mut u = 0.0;
        for atom in pbe.alist.atoms() {
            let charge = EC * atom.charge;
            let size = atom.radius;
            let dx = atom.position[0] - pt[0];
            let dy = atom.position[1] - pt[1];
            let dz = atom.position[2] - pt[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-8);
            let mut val = charge / (4.0 * std::f64::consts::PI * EPS0 * eps_w * dist);
            if xkappa.abs() > 0.0 {
                val *= (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
            }
            u += val * EC / (KB * t);
        }
        u
    }

    /// Get potential at a point
    pub fn pot(&self, pt: [f64; 3]) -> ApbsResult<f64> {
        match self.mgrid.value(pt) {
            Ok(v) => Ok(v),
            Err(_) => {
                // Off-grid: apply boundary condition
                match self.bcfl {
                    Vbcfl::Zero => Ok(0.0),
                    Vbcfl::SDH => {
                        let pbe = self.pbe.as_ref().ok_or_else(|| {
                            ApbsError::InvalidParameter(
                                "Vopot SDH off-grid requires bound Vpbe".to_string(),
                            )
                        })?;
                        self.sdh_offgrid_potential(pbe, pt)
                    }
                    Vbcfl::MDH => {
                        let pbe = self.pbe.as_ref().ok_or_else(|| {
                            ApbsError::InvalidParameter(
                                "Vopot MDH off-grid requires bound Vpbe".to_string(),
                            )
                        })?;
                        Ok(self.mdh_offgrid_potential(pbe, pt))
                    }
                    Vbcfl::Mem => {
                        let pbe = self.pbe.as_ref().ok_or_else(|| {
                            ApbsError::InvalidParameter(
                                "Vopot MEM off-grid requires bound Vpbe".to_string(),
                            )
                        })?;
                        let base = self.mdh_offgrid_potential(pbe, pt);
                        Ok(base + pbe.z_mem * pbe.memv)
                    }
                    Vbcfl::Focus | Vbcfl::Map => Err(ApbsError::UnsupportedFormat(format!(
                        "vopot off-grid evaluation for boundary {:?} is invalid in APBS semantics",
                        self.bcfl
                    ))),
                    _ => Err(ApbsError::UnsupportedFormat(format!(
                        "vopot off-grid evaluation for boundary {:?} is not implemented",
                        self.bcfl
                    ))),
                }
            }
        }
    }

    /// Get gradient at a point
    pub fn gradient(&self, pt: [f64; 3]) -> ApbsResult<[f64; 3]> {
        if let Some(grid) = self.mgrid.get_grid_by_point(pt) {
            grid.gradient(pt)
        } else {
            match self.bcfl {
                Vbcfl::Zero => Ok([0.0; 3]),
                Vbcfl::SDH => {
                    let pbe = self.pbe.as_ref().ok_or_else(|| {
                        ApbsError::InvalidParameter(
                            "Vopot SDH off-grid gradient requires bound Vpbe".to_string(),
                        )
                    })?;
                    let eps_w = pbe.solvent_diel;
                    let xkappa = pbe.xkappa;
                    let t = pbe.temperature;
                    let size = pbe.solute_radius;
                    let center = pbe.solute_center;
                    let charge = EC * pbe.solute_charge;
                    let dx = center[0] - pt[0];
                    let dy = center[1] - pt[1];
                    let dz = center[2] - pt[2];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-8);
                    let mut val = charge / (4.0 * std::f64::consts::PI * EPS0 * eps_w);
                    if xkappa.abs() > 0.0 {
                        val *= (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
                    }
                    val *= EC / (KB * t);
                    Ok([
                        val * dx / dist * (-1.0 / (dist * dist) + xkappa / dist),
                        val * dy / dist * (-1.0 / (dist * dist) + xkappa / dist),
                        val * dz / dist * (-1.0 / (dist * dist) + xkappa / dist),
                    ])
                }
                Vbcfl::MDH | Vbcfl::Mem => {
                    let pbe = self.pbe.as_ref().ok_or_else(|| {
                        ApbsError::InvalidParameter(
                            format!("Vopot {:?} off-grid gradient requires bound Vpbe", self.bcfl),
                        )
                    })?;
                    let eps_w = pbe.solvent_diel;
                    let xkappa = pbe.xkappa;
                    let t = pbe.temperature;
                    let mut grad = [0.0; 3];
                    for atom in pbe.alist.atoms() {
                        let charge = EC * atom.charge;
                        let size = atom.radius;
                        let dx = atom.position[0] - pt[0];
                        let dy = atom.position[1] - pt[1];
                        let dz = atom.position[2] - pt[2];
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-8);
                        let mut val = charge / (4.0 * std::f64::consts::PI * EPS0 * eps_w);
                        if xkappa.abs() > 0.0 {
                            val *= (-xkappa * (dist - size)).exp() / (1.0 + xkappa * size);
                        }
                        val *= EC / (KB * t);
                        let factor = -1.0 / (dist * dist) + xkappa / dist;
                        grad[0] += val * dx / dist * factor;
                        grad[1] += val * dy / dist * factor;
                        grad[2] += val * dz / dist * factor;
                    }
                    Ok(grad)
                }
                Vbcfl::Focus | Vbcfl::Map => Err(ApbsError::UnsupportedFormat(format!(
                    "vopot off-grid gradient for boundary {:?} is invalid in APBS semantics",
                    self.bcfl
                ))),
                _ => Err(ApbsError::UnsupportedFormat(format!(
                    "vopot off-grid gradient for boundary {:?} is not implemented",
                    self.bcfl
                ))),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Vopot;
    use crate::vmgrid::Vmgrid;
    use apbs_generic::vhal::Vbcfl;

    #[test]
    fn offgrid_zero_returns_zero() {
        let vmgrid = Vmgrid::new();
        let vopot = Vopot::new(vmgrid, Vbcfl::Zero);
        let val = vopot.pot([100.0, 100.0, 100.0]).unwrap_or(1.0);
        assert_eq!(val, 0.0);
    }

    #[test]
    fn offgrid_sdh_without_pbe_errors() {
        let vmgrid = Vmgrid::new();
        let vopot = Vopot::new(vmgrid, Vbcfl::SDH);
        let err = vopot.pot([100.0, 100.0, 100.0]).err();
        assert!(err.is_some());
        let msg = format!("{}", err.unwrap_or_else(|| panic!("expected error")));
        assert!(msg.contains("requires bound Vpbe"));
    }

    #[test]
    fn offgrid_focus_is_explicit_error() {
        let vmgrid = Vmgrid::new();
        let vopot = Vopot::new(vmgrid, Vbcfl::Focus);
        let err = vopot.pot([100.0, 100.0, 100.0]).err();
        assert!(err.is_some());
        let msg = format!("{}", err.unwrap_or_else(|| panic!("expected error")));
        assert!(msg.contains("invalid in APBS semantics"));
    }

    #[test]
    fn offgrid_focus_gradient_is_explicit_error() {
        let vmgrid = Vmgrid::new();
        let vopot = Vopot::new(vmgrid, Vbcfl::Focus);
        let err = vopot.gradient([100.0, 100.0, 100.0]).err();
        assert!(err.is_some());
        let msg = format!("{}", err.unwrap_or_else(|| panic!("expected error")));
        assert!(msg.contains("invalid in APBS semantics"));
    }
}
