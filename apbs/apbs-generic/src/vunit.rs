// APBS vunit.rs - Physical unit conversion constants
// Port of src/generic/vunit.h

/// Joules to calories
pub const J_TO_CAL: f64 = 4.184e+00;

/// Calories to joules
pub const CAL_TO_J: f64 = 2.3900574e-01;

/// AMU to kg
pub const AMU_TO_KG: f64 = 1.6605402e-27;

/// kg to AMU
pub const KG_TO_AMU: f64 = 6.0221367e+26;

/// Elementary charge to Coulombs
pub const EC_TO_C: f64 = 1.6021773e-19;

/// Coulombs to elementary charge
pub const C_TO_EC: f64 = 6.2415065e+18;

/// Electron charge (Coulombs)
pub const EC: f64 = 1.6021773e-19;

/// Boltzmann constant (J/K)
pub const KB: f64 = 1.3806581e-23;

/// Avogadro's number
pub const NA: f64 = 6.0221367e+23;

/// Pi
pub const PI: f64 = std::f64::consts::PI;

/// Vacuum permittivity (F/m)
pub const EPS0: f64 = 8.8541878e-12;

/// e_c^2/A in ESU => kcal/mol
pub const ESU_EC2A: f64 = 3.3206364e+02;

/// k_B in ESU => kcal/mol
pub const ESU_KB: f64 = 1.9871913e-03;

/// Water dielectric constant (at 25C)
pub const EPSILON_WATER: f64 = 78.36;

/// Temperature of water (K, at 25C)
pub const TEMP_WATER: f64 = 298.15;

/// Elementary charge in ESU (statcoulombs)
pub const EC_ESU: f64 = 4.803242384e-10;

/// Boltzmann constant in erg/K
pub const KB_ERG: f64 = 1.380662e-16;

/// Avogadro's number (ESU units version)
pub const NA_ESU: f64 = 6.022045e+23;
