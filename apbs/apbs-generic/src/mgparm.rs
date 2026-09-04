// APBS mgparm.rs - Multigrid parameters
// Port of src/generic/mgparm.h / mgparm.c

use crate::error::{ApbsError, ApbsResult};
use crate::vhal::VchrgMeth;

/// MG calculation type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MGparmCalcType {
    Manual = 0,
    Auto = 1,
    Parallel = 2,
    Dummy = 3,
    None = 4,
}

/// Grid centering method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MGparmCentMeth {
    Point = 0,
    Molecule = 1,
    Focus = 2,
}

/// Multigrid parameters
#[derive(Debug, Clone)]
pub struct MGparm {
    /// Calculation type
    pub r#type: MGparmCalcType,
    /// Has been parsed
    pub parsed: bool,
    /// Grid dimensions (nx, ny, nz)
    pub dime: [i32; 3],
    /// Charge discretization method
    pub chgm: VchrgMeth,
    /// Charge source
    pub chgs: crate::vhal::VchrgSrc,
    /// MG levels (deprecated)
    pub nlev: i32,
    /// Error tolerance
    pub etol: f64,
    /// Grid spacings (hx, hy, hz)
    pub grid: [f64; 3],
    /// Grid side lengths (glenx, gleny, glenz)
    pub glen: [f64; 3],
    /// Centering method
    pub cmeth: MGparmCentMeth,
    /// Grid center
    pub center: [f64; 3],
    /// Molecule index for centering
    pub centmol: i32,
    /// Coarse grid side lengths
    pub cglen: [f64; 3],
    /// Coarse centering method
    pub ccmeth: MGparmCentMeth,
    /// Coarse grid center
    pub ccenter: [f64; 3],
    /// Coarse grid molecule index
    pub ccentmol: i32,
    /// Fine grid side lengths
    pub fglen: [f64; 3],
    /// Fine centering method
    pub fcmeth: MGparmCentMeth,
    /// Fine grid center
    pub fcenter: [f64; 3],
    /// Fine grid molecule index
    pub fcentmol: i32,
    /// Disjoint partition center
    pub part_disj_center: [f64; 3],
    /// Disjoint partition lengths
    pub part_disj_length: [f64; 3],
    /// Boundary ownership flags
    pub part_disj_own_side: [i32; 6],
    /// Processor grid dimensions
    pub pdime: [i32; 3],
    /// Processor rank
    pub proc_rank: i32,
    /// Total processors
    pub proc_size: i32,
    /// Overlap fraction
    pub ofrac: f64,
    /// Async processor ID
    pub r#async: i32,
    /// Nonlinearity type
    pub nonlintype: i32,
    /// Solver method
    pub method: i32,
    /// Enable lpbe/aqua
    pub use_aqua: i32,
}

impl Default for MGparm {
    fn default() -> Self {
        Self {
            r#type: MGparmCalcType::None,
            parsed: false,
            dime: [0; 3],
            chgm: VchrgMeth::Tril,
            chgs: crate::vhal::VchrgSrc::Charge,
            nlev: 0,
            etol: 1.0e-9,
            grid: [0.0; 3],
            glen: [0.0; 3],
            cmeth: MGparmCentMeth::Molecule,
            center: [0.0; 3],
            centmol: 0,
            cglen: [0.0; 3],
            ccmeth: MGparmCentMeth::Molecule,
            ccenter: [0.0; 3],
            ccentmol: 0,
            fglen: [0.0; 3],
            fcmeth: MGparmCentMeth::Molecule,
            fcenter: [0.0; 3],
            fcentmol: 0,
            part_disj_center: [0.0; 3],
            part_disj_length: [0.0; 3],
            part_disj_own_side: [0; 6],
            pdime: [1, 1, 1],
            proc_rank: 0,
            proc_size: 1,
            ofrac: 0.0,
            r#async: -1,
            nonlintype: 0,
            method: 2,
            use_aqua: 0,
        }
    }
}

impl MGparm {
    pub fn new(calc_type: MGparmCalcType) -> Self {
        Self {
            r#type: calc_type,
            ..Default::default()
        }
    }

    pub fn get_nx(&self) -> i32 {
        self.dime[0]
    }

    pub fn get_ny(&self) -> i32 {
        self.dime[1]
    }

    pub fn get_nz(&self) -> i32 {
        self.dime[2]
    }

    pub fn get_hx(&self) -> f64 {
        self.grid[0]
    }

    pub fn get_hy(&self) -> f64 {
        self.grid[1]
    }

    pub fn get_hz(&self) -> f64 {
        self.grid[2]
    }

    pub fn check(&mut self) -> ApbsResult<()> {
        match self.r#type {
            MGparmCalcType::Manual => {
                for d in 0..3 {
                    if self.dime[d] <= 0 || (self.dime[d] as usize) % 2 == 0 {
                        return Err(ApbsError::InvalidParameter(format!(
                            "Grid dimension dime[{}] = {} must be odd and positive",
                            d, self.dime[d]
                        )));
                    }
                    if self.glen[d] <= 0.0 {
                        return Err(ApbsError::InvalidParameter(format!(
                            "Grid length glen[{}] = {} must be positive",
                            d, self.glen[d]
                        )));
                    }
                }
            }
            MGparmCalcType::Auto => {
                for d in 0..3 {
                    if self.dime[d] <= 0 || (self.dime[d] as usize) % 2 == 0 {
                        return Err(ApbsError::InvalidParameter(format!(
                            "Grid dimension dime[{}] = {} must be odd and positive",
                            d, self.dime[d]
                        )));
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }

    pub fn copy_from(&mut self, other: &MGparm) {
        *self = other.clone();
    }

    pub fn parse_token(&mut self, token: &str, value: &str) -> ApbsResult<()> {
        match token.to_ascii_lowercase().as_str() {
            "dime" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() >= 3 {
                    for d in 0..3 {
                        self.dime[d] = vals[d].parse().map_err(|_| ApbsError::Parse {
                            line: 0,
                            message: format!("Invalid dime[{}]: {}", d, vals[d]),
                        })?;
                    }
                }
            }
            "glen" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() >= 3 {
                    for d in 0..3 {
                        self.glen[d] = vals[d].parse().map_err(|_| ApbsError::Parse {
                            line: 0,
                            message: format!("Invalid glen[{}]: {}", d, vals[d]),
                        })?;
                    }
                }
            }
            "grid" => {
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
            "nlev" => {
                self.nlev = value.trim().parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid nlev: {}", value),
                })?;
            }
            "chgm" => {
                self.chgm = match value.to_ascii_lowercase().as_str() {
                    "tril" => VchrgMeth::Tril,
                    "bspl2" | "spl2" => VchrgMeth::Bspl2,
                    "bspl4" | "spl4" => VchrgMeth::Bspl4,
                    _ => return Err(ApbsError::InvalidParameter(format!("Unknown chgm: {}", value))),
                };
            }
            "cglen" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() >= 3 {
                    for d in 0..3 {
                        self.cglen[d] = vals[d].parse().map_err(|_| ApbsError::Parse {
                            line: 0,
                            message: format!("Invalid cglen[{}]: {}", d, vals[d]),
                        })?;
                    }
                }
            }
            "fglen" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() >= 3 {
                    for d in 0..3 {
                        self.fglen[d] = vals[d].parse().map_err(|_| ApbsError::Parse {
                            line: 0,
                            message: format!("Invalid fglen[{}]: {}", d, vals[d]),
                        })?;
                    }
                }
            }
            "gcent" | "cgcent" | "fgcent" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                let (target_center, target_method, target_mol, label) = match token {
                    "gcent" => (&mut self.center, &mut self.cmeth, &mut self.centmol, "gcent"),
                    "cgcent" => (&mut self.ccenter, &mut self.ccmeth, &mut self.ccentmol, "cgcent"),
                    "fgcent" => (&mut self.fcenter, &mut self.fcmeth, &mut self.fcentmol, "fgcent"),
                    _ => unreachable!(),
                };
                if vals.len() >= 2 && vals[0].eq_ignore_ascii_case("mol") {
                    *target_mol = vals[1].parse().map_err(|_| ApbsError::Parse {
                        line: 0,
                        message: format!("Invalid {} molecule index: {}", label, vals[1]),
                    })?;
                    *target_method = MGparmCentMeth::Molecule;
                } else if vals.len() >= 3 {
                    for d in 0..3 {
                        target_center[d] = vals[d].parse().map_err(|_| ApbsError::Parse {
                            line: 0,
                            message: format!("Invalid {}[{}]: {}", label, d, vals[d]),
                        })?;
                    }
                    *target_method = MGparmCentMeth::Point;
                }
            }
            _ => {
                return Err(ApbsError::InvalidParameter(format!(
                    "Unknown multigrid parameter: {}",
                    token
                )));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{MGparm, MGparmCalcType};

    #[test]
    fn parse_token_reads_nlev() {
        let mut mgparm = MGparm::new(MGparmCalcType::Manual);
        mgparm.parse_token("nlev", "3").expect("parse nlev");
        assert_eq!(mgparm.nlev, 3);
    }
}
