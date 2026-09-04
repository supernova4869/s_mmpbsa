// APBS pbeparm.rs - Generic PBE parameters
// Port of src/generic/pbeparm.h / pbeparm.c

use std::fmt;

use crate::error::{ApbsError, ApbsResult};
use crate::vhal::{Vbcfl, VdataFormat, VdataType, VsurfMeth, VhalPBEType, MAXION};

/// Maximum write statements per calculation
pub const PBEPARM_MAXWRITE: usize = 20;

/// Energy calculation mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PBEparmCalcEnergy {
    No = 0,
    Total = 1,
    Comps = 2,
}

/// Force calculation mode
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PBEparmCalcForce {
    No = 0,
    Total = 1,
    Comps = 2,
}

/// Write specification
#[derive(Debug, Clone)]
pub struct WriteSpec {
    pub stem: String,
    pub writetype: VdataType,
    pub writefmt: VdataFormat,
}

/// Write type for output (alias for VdataType usage in output context)
pub type WriteType = VdataType;

/// Write format for output
pub type WriteFormat = VdataFormat;

/// Generic PBE parameters
#[derive(Debug, Clone)]
pub struct PBEparm {
    /// Molecule ID
    pub molid: Option<i32>,
    /// External dielectric map
    pub use_diel_map: bool,
    pub diel_map_id: Option<i32>,
    /// External kappa map
    pub use_kappa_map: bool,
    pub kappa_map_id: Option<i32>,
    /// External potential map
    pub use_pot_map: bool,
    pub pot_map_id: Option<i32>,
    /// External charge map
    pub use_charge_map: bool,
    pub charge_map_id: Option<i32>,
    /// PBE version
    pub pbetype: VhalPBEType,
    /// Boundary condition method
    pub bcfl: Vbcfl,
    /// Number of counterion species
    pub nion: usize,
    /// Counterion charges (e)
    pub ionq: [f64; MAXION],
    /// Counterion concentrations (M)
    pub ionc: [f64; MAXION],
    /// Counterion radii (A)
    pub ionr: [f64; MAXION],
    /// Solute dielectric
    pub pdie: f64,
    /// Solvent dielectric
    pub sdie: f64,
    /// Vacc sphere density
    pub sdens: f64,
    /// Surface method
    pub srfm: VsurfMeth,
    /// Solvent radius
    pub srad: f64,
    /// Cubic spline window
    pub swin: f64,
    /// Temperature (K)
    pub temp: f64,
    /// SMPBE size parameter
    pub smsize: f64,
    /// SMPBE volume parameter
    pub smvolume: f64,
    /// Energy calculation flag
    pub calc_energy: PBEparmCalcEnergy,
    /// Force calculation flag
    pub calc_force: PBEparmCalcForce,
    /// Membrane bottom z (A)
    pub zmem: f64,
    /// Membrane width (A)
    pub lmem: f64,
    /// Membrane dielectric
    pub mdie: f64,
    /// Membrane potential
    pub memv: f64,
    /// Write statements
    pub writes: Vec<WriteSpec>,
    /// Write operator matrix
    pub writemat: bool,
    pub writematstem: String,
    pub writematflag: i32,
    /// PBAM 3D map
    pub pbam_3dmapstem: String,
    pub pbam_3dmapflag: bool,
    /// Whether fully parsed
    pub parsed: bool,
}

impl Default for PBEparm {
    fn default() -> Self {
        Self {
            molid: None,
            use_diel_map: false,
            diel_map_id: None,
            use_kappa_map: false,
            kappa_map_id: None,
            use_pot_map: false,
            pot_map_id: None,
            use_charge_map: false,
            charge_map_id: None,
            pbetype: VhalPBEType::NPBE,
            bcfl: Vbcfl::SDH,
            nion: 0,
            ionq: [0.0; MAXION],
            ionc: [0.0; MAXION],
            ionr: [0.0; MAXION],
            pdie: 2.0,
            sdie: 78.36,
            sdens: 10.0,
            srfm: VsurfMeth::Mol,
            srad: 1.4,
            swin: 0.3,
            temp: 298.15,
            smsize: 0.0,
            smvolume: 0.0,
            calc_energy: PBEparmCalcEnergy::No,
            calc_force: PBEparmCalcForce::No,
            zmem: 0.0,
            lmem: 0.0,
            mdie: 2.0,
            memv: 0.0,
            writes: Vec::new(),
            writemat: false,
            writematstem: String::new(),
            writematflag: 0,
            pbam_3dmapstem: String::new(),
            pbam_3dmapflag: false,
            parsed: false,
        }
    }
}

impl PBEparm {
    pub fn new() -> Self {
        Self::default()
    }

    /// Get ion charge
    pub fn get_ion_charge(&self, iion: usize) -> f64 {
        self.ionq[iion]
    }

    /// Get ion concentration
    pub fn get_ion_conc(&self, iion: usize) -> f64 {
        self.ionc[iion]
    }

    /// Get ion radius
    pub fn get_ion_radius(&self, iion: usize) -> f64 {
        self.ionr[iion]
    }

    /// Check parameters for consistency
    pub fn check(&self) -> ApbsResult<()> {
        if self.pdie <= 0.0 {
            return Err(ApbsError::InvalidParameter(
                "Solute dielectric must be positive".to_string(),
            ));
        }
        if self.sdie <= 0.0 {
            return Err(ApbsError::InvalidParameter(
                "Solvent dielectric must be positive".to_string(),
            ));
        }
        if self.temp <= 0.0 {
            return Err(ApbsError::InvalidParameter(
                "Temperature must be positive".to_string(),
            ));
        }
        if self.nion > MAXION {
            return Err(ApbsError::InvalidParameter(format!(
                "Number of ion species ({}) exceeds maximum ({})",
                self.nion, MAXION
            )));
        }
        Ok(())
    }

    /// Copy parameters from another PBEparm
    pub fn copy_from(&mut self, other: &PBEparm) {
        *self = other.clone();
    }

    /// Parse a token from the input file
    pub fn parse_token(&mut self, token: &str, value: &str) -> ApbsResult<()> {
        match token.to_ascii_lowercase().as_str() {
            "molid" => {
                self.molid = value.parse().ok();
            }
            "pdie" => {
                self.pdie = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid pdie value: {}", value),
                })?;
            }
            "sdie" => {
                self.sdie = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid sdie value: {}", value),
                })?;
            }
            "temp" => {
                self.temp = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid temp value: {}", value),
                })?;
            }
            "usemap" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() < 2 {
                    return Err(ApbsError::InvalidParameter(format!(
                        "usemap requires '<type> <id>', got '{}'",
                        value
                    )));
                }
                let map_id = vals[1].parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid usemap id: {}", vals[1]),
                })?;
                match vals[0].to_ascii_lowercase().as_str() {
                    "diel" => {
                        self.use_diel_map = true;
                        self.diel_map_id = Some(map_id);
                    }
                    "kappa" => {
                        self.use_kappa_map = true;
                        self.kappa_map_id = Some(map_id);
                    }
                    "pot" => {
                        self.use_pot_map = true;
                        self.pot_map_id = Some(map_id);
                    }
                    "charge" => {
                        self.use_charge_map = true;
                        self.charge_map_id = Some(map_id);
                    }
                    _ => {
                        return Err(ApbsError::InvalidParameter(format!(
                            "Invalid usemap type: {}",
                            vals[0]
                        )));
                    }
                }
            }
            "srad" => {
                self.srad = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid srad value: {}", value),
                })?;
            }
            "sdens" => {
                self.sdens = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid sdens value: {}", value),
                })?;
            }
            "swin" => {
                self.swin = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid swin value: {}", value),
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
            "bcfl" => {
                self.bcfl = match value.parse::<i32>() {
                    Ok(0) => Vbcfl::Zero,
                    Ok(1) => Vbcfl::SDH,
                    Ok(2) => Vbcfl::MDH,
                    Ok(4) => Vbcfl::Focus,
                    Ok(5) => Vbcfl::Mem,
                    Ok(6) => Vbcfl::Map,
                    _ => match value.to_ascii_lowercase().as_str() {
                        "zero" => Vbcfl::Zero,
                        "sdh" => Vbcfl::SDH,
                        "mdh" => Vbcfl::MDH,
                        "focus" => Vbcfl::Focus,
                        "mem" => Vbcfl::Mem,
                        "map" => Vbcfl::Map,
                        _ => return Err(ApbsError::InvalidParameter(format!("Invalid bcfl: {}", value))),
                    }
                };
            }
            "npbe" => {
                self.pbetype = VhalPBEType::NPBE;
            }
            "lpbe" => {
                self.pbetype = VhalPBEType::LPBE;
            }
            "nion" => {
                self.nion = value.parse().map_err(|_| ApbsError::Parse {
                    line: 0,
                    message: format!("Invalid nion value: {}", value),
                })?;
            }
            "ionq" | "ionc" | "ionr" => {
                // These are parsed elsewhere with index
            }
            "calcenergy" => {
                self.calc_energy = match value.to_ascii_lowercase().as_str() {
                    "no" => PBEparmCalcEnergy::No,
                    "total" => PBEparmCalcEnergy::Total,
                    "comps" => PBEparmCalcEnergy::Comps,
                    _ => return Err(ApbsError::InvalidParameter(format!("Invalid calcenergy: {}", value))),
                };
            }
            "calcforce" => {
                self.calc_force = match value.to_ascii_lowercase().as_str() {
                    "no" => PBEparmCalcForce::No,
                    "total" => PBEparmCalcForce::Total,
                    "comps" => PBEparmCalcForce::Comps,
                    _ => return Err(ApbsError::InvalidParameter(format!("Invalid calcforce: {}", value))),
                };
            }
            "write" => {
                let vals: Vec<&str> = value.split_whitespace().collect();
                if vals.len() < 3 {
                    return Err(ApbsError::InvalidParameter(format!(
                        "write requires '<type> <format> <stem>', got '{}'",
                        value
                    )));
                }
                let writetype = match vals[0].to_ascii_lowercase().as_str() {
                    "charge" => VdataType::Charge,
                    "pot" => VdataType::Pot,
                    "atmpot" | "atompot" => VdataType::AtomPot,
                    "smol" => VdataType::Smol,
                    "sspl" => VdataType::Sspl,
                    "vdw" => VdataType::Vdw,
                    "ivdw" => VdataType::Ivdw,
                    "lap" => VdataType::Lap,
                    "edens" => VdataType::Edens,
                    "ndens" => VdataType::Ndens,
                    "qdens" => VdataType::Qdens,
                    "dielx" => VdataType::DielX,
                    "diely" => VdataType::DielY,
                    "dielz" => VdataType::DielZ,
                    "kappa" => VdataType::Kappa,
                    _ => return Err(ApbsError::InvalidParameter(format!("Invalid write type: {}", vals[0]))),
                };
                let writefmt = match vals[1].to_ascii_lowercase().as_str() {
                    "dx" => VdataFormat::DX,
                    "uhbd" => VdataFormat::UHBD,
                    "avs" => VdataFormat::AVS,
                    "mcsf" => VdataFormat::MCSF,
                    "gz" => VdataFormat::GZ,
                    "flat" => VdataFormat::Flat,
                    "dxbin" => VdataFormat::DXBin,
                    _ => return Err(ApbsError::InvalidParameter(format!("Invalid write format: {}", vals[1]))),
                };
                self.writes.push(WriteSpec {
                    stem: vals[2].to_string(),
                    writetype,
                    writefmt,
                });
            }
            _ => {
                return Err(ApbsError::InvalidParameter(format!(
                    "Unknown PBE parameter: {}",
                    token
                )));
            }
        }
        Ok(())
    }
}

impl fmt::Display for PBEparm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "PBE: pbetype={:?}, bcfl={:?}, pdie={}, sdie={}, temp={}, srad={}, nion={}",
            self.pbetype, self.bcfl, self.pdie, self.sdie, self.temp, self.srad, self.nion
        )
    }
}
