// APBS apolparm.rs - Non-polar (apolar) calculation parameters
// Port of src/generic/apolparm.h / apolparm.c

use crate::error::{ApbsError, ApbsResult};
use crate::vhal::VsurfMeth;

/// Apolar energy calculation mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum APOLparmCalcEnergy {
    No = 0,
    Total = 1,
    Comps = 2,
}

/// Apolar force calculation mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum APOLparmCalcForce {
    No = 0,
    Total = 1,
    Comps = 2,
}

/// Apolar do-calculation flag
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum APOLparmDoCalc {
    No = 0,
    Yes = 1,
    Error = 2,
}

/// Non-polar (apolar) calculation parameters
#[derive(Debug, Clone)]
pub struct APOLparm {
    pub parsed: bool,
    /// Whether user explicitly set the grid spacing
    pub setgrid: bool,
    /// Grid spacing (hx, hy, hz)
    pub grid: [f64; 3],
    /// Molecule ID
    pub molid: i32,
    /// Bulk concentration (M)
    pub bconc: f64,
    /// Vacc sphere density
    pub sdens: f64,
    /// Atom position offset (A)
    pub dpos: f64,
    /// Solvent pressure (kJ/(mol*A^3))
    pub press: f64,
    /// Surface method
    pub srfm: VsurfMeth,
    /// Solvent radius (A)
    pub srad: f64,
    /// Cubic spline window
    pub swin: f64,
    /// Temperature (K)
    pub temp: f64,
    /// Surface tension (kJ/(mol*A^2))
    pub gamma: f64,
    /// Energy calculation flag
    pub calc_energy: APOLparmCalcEnergy,
    /// Force calculation flag
    pub calc_force: APOLparmCalcForce,
    /// Water O LJ radius (A)
    pub watsigma: f64,
    /// Water O LJ well depth (kJ/mol)
    pub watepsilon: f64,
    /// SASA result
    pub sasa: f64,
    /// SAV result
    pub sav: f64,
    /// WCA energy result
    pub wca_energy: f64,
    /// Total forces (x, y, z)
    pub tot_force: [f64; 3],
}

impl Default for APOLparm {
    fn default() -> Self {
        Self {
            parsed: false,
            setgrid: false,
            grid: [0.5; 3],
            molid: 0,
            bconc: 0.0,
            sdens: 10.0,
            dpos: 0.2,
            press: 0.0,
            srfm: VsurfMeth::Mol,
            srad: 1.4,
            swin: 0.3,
            temp: 298.15,
            gamma: 0.072,
            calc_energy: APOLparmCalcEnergy::No,
            calc_force: APOLparmCalcForce::No,
            watsigma: 3.166,
            watepsilon: 0.650,
            sasa: 0.0,
            sav: 0.0,
            wca_energy: 0.0,
            tot_force: [0.0; 3],
        }
    }
}

impl APOLparm {
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse a token from the input file
    pub fn parse_token(&mut self, token: &str, value: &str) -> ApbsResult<()> {
        match token.to_ascii_lowercase().as_str() {
            "mol" => {
                if let Ok(id) = value.parse::<i32>() {
                    self.molid = id;
                }
            }
            "temp" => {
                self.temp = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid temp value: {}", value),
                })?;
            }
            "srfm" => {
                self.srfm = match value.to_ascii_lowercase().as_str() {
                    "mol" => VsurfMeth::Mol,
                    "smol" | "molsmooth" => VsurfMeth::MolSmooth,
                    "spline" => VsurfMeth::Spline,
                    "spline3" => VsurfMeth::Spline3,
                    "spline4" => VsurfMeth::Spline4,
                    "sacc" => VsurfMeth::Sacc,
                    _ => return Err(ApbsError::InvalidParameter(format!("Unknown srfm: {}", value))),
                };
            }
            "swin" => {
                self.swin = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid swin value: {}", value),
                })?;
            }
            "srad" => {
                self.srad = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid srad value: {}", value),
                })?;
            }
            "gamma" => {
                self.gamma = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid gamma value: {}", value),
                })?;
            }
            "press" => {
                self.press = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid press value: {}", value),
                })?;
            }
            "bconc" => {
                self.bconc = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid bconc value: {}", value),
                })?;
            }
            "sdens" => {
                self.sdens = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid sdens value: {}", value),
                })?;
            }
            "dpos" => {
                self.dpos = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid dpos value: {}", value),
                })?;
            }
            "grid" => {
                self.setgrid = true;
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() >= 3 {
                    for d in 0..3 {
                        self.grid[d] = vals[d].parse().map_err(|_| ApbsError::Parse {
                            line: 0,
                            message: format!("Invalid grid[{}]: {}", d, vals[d]),
                        })?;
                    }
                }
            }
            "calcenergy" => {
                self.calc_energy = match value.to_ascii_lowercase().as_str() {
                    "no" => APOLparmCalcEnergy::No,
                    "total" => APOLparmCalcEnergy::Total,
                    "comps" => APOLparmCalcEnergy::Comps,
                    _ => return Err(ApbsError::InvalidParameter(format!("Invalid calcenergy: {}", value))),
                };
            }
            "calcforce" => {
                self.calc_force = match value.to_ascii_lowercase().as_str() {
                    "no" => APOLparmCalcForce::No,
                    "total" => APOLparmCalcForce::Total,
                    "comps" => APOLparmCalcForce::Comps,
                    _ => return Err(ApbsError::InvalidParameter(format!("Invalid calcforce: {}", value))),
                };
            }
            _ => {}
        }
        Ok(())
    }

    pub fn check(&self) -> ApbsResult<()> {
        if self.srad <= 0.0 {
            return Err(ApbsError::InvalidParameter(
                "Solvent radius must be positive".to_string(),
            ));
        }
        if self.temp <= 0.0 {
            return Err(ApbsError::InvalidParameter(
                "Temperature must be positive".to_string(),
            ));
        }
        Ok(())
    }

    pub fn copy_from(&mut self, other: &APOLparm) {
        *self = other.clone();
    }
}
