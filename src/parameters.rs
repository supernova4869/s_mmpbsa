use std::{fmt, fs};
use std::fmt::Formatter;
use std::fs::File;
use std::io::Write;
use std::marker::Copy;
use std::path::Path;
use serde::{Serialize, Deserialize};

/// PB equation solved by APBS.  The tokens must match what the built-in
/// solver's input parser accepts (`lpbe`, `npbe`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PBSolver {
    /// Linear PB equation (recommended).
    Lpbe,
    /// Nonlinear PB equation, for highly charged systems.
    Npbe,
}

impl fmt::Display for PBSolver {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            PBSolver::Lpbe => write!(f, "lpbe"),
            PBSolver::Npbe => write!(f, "npbe"),
        }
    }
}

/// Boundary condition of the coarse-grid PB equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Bcfl {
    Zero,
    /// Single Debye-Hückel.
    Sdh,
    /// Multiple Debye-Hückel (recommended).
    Mdh,
    /// Focusing from a previous coarse solve.
    Focus,
    Mem,
    Map,
}

impl fmt::Display for Bcfl {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let token = match self {
            Bcfl::Zero => "zero",
            Bcfl::Sdh => "sdh",
            Bcfl::Mdh => "mdh",
            Bcfl::Focus => "focus",
            Bcfl::Mem => "mem",
            Bcfl::Map => "map",
        };
        write!(f, "{token}")
    }
}

/// Model used to build the dielectric / ion-accessibility surface.
///
/// Canonical tokens are the ones the built-in solver accepts; the serde
/// aliases keep the historic spellings of old settings files loadable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Srfm {
    Mol,
    #[serde(alias = "molsmooth")]
    Smol,
    #[serde(alias = "spl2")]
    Spline,
    #[serde(alias = "spl3")]
    Spline3,
    #[serde(alias = "spl4")]
    Spline4,
    Sacc,
}

impl fmt::Display for Srfm {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let token = match self {
            Srfm::Mol => "mol",
            Srfm::Smol => "smol",
            Srfm::Spline => "spline",
            Srfm::Spline3 => "spline3",
            Srfm::Spline4 => "spline4",
            Srfm::Sacc => "sacc",
        };
        write!(f, "{token}")
    }
}

/// Charge mapping onto the grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Chgm {
    Tril,
    #[serde(alias = "bspl2")]
    Spl2,
    #[serde(alias = "bspl4")]
    Spl4,
}

impl fmt::Display for Chgm {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let token = match self {
            Chgm::Tril => "tril",
            Chgm::Spl2 => "spl2",
            Chgm::Spl4 => "spl4",
        };
        write!(f, "{token}")
    }
}

/// Which energies the solver reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CalcEnergy {
    No,
    Total,
    /// Per-atom energy components (what s_mmpbsa needs).
    Comps,
}

impl fmt::Display for CalcEnergy {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let token = match self {
            CalcEnergy::No => "no",
            CalcEnergy::Total => "total",
            CalcEnergy::Comps => "comps",
        };
        write!(f, "{token}")
    }
}

/// Named PB parameter combinations offered by the interactive menu.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PbePreset {
    Default,
    FastScreening,
    FineMesh,
    Nonlinear,
}

impl PbePreset {
    /// `(preset, label, one-line description)` in menu order.
    pub fn all() -> [(PbePreset, &'static str, &'static str); 4] {
        [
            (PbePreset::Default, "default",
                "linear PB, mdh boundary, smol surface, spl4 charges, df = 0.5 A"),
            (PbePreset::FastScreening, "fast",
                "coarse mesh df = 1.0 A, for quickly scanning a system"),
            (PbePreset::FineMesh, "fine",
                "dense mesh df = 0.25 A, for production runs"),
            (PbePreset::Nonlinear, "npbe",
                "nonlinear PB equation, for highly charged systems"),
        ]
    }

    /// Overwrites the PB fields the preset owns; everything else (temperature,
    /// dielectrics, ions, grid padding) the user configured stays untouched.
    pub fn apply(&self, pbe: &mut PBESet) {
        match self {
            PbePreset::Default => {
                let (cfac, fadd, temp) = (pbe.cfac, pbe.fadd, pbe.temp);
                *pbe = PBESet::new(temp);
                pbe.cfac = cfac;
                pbe.fadd = fadd;
            }
            PbePreset::FastScreening => pbe.df = 1.0,
            PbePreset::FineMesh => pbe.df = 0.25,
            PbePreset::Nonlinear => pbe.pb_solver = PBSolver::Npbe,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct Config {
    pub program_set: ProgramSet,
    pub mm_set: MMSet,
    pub pbe_set: PBESet,
    pub pba_set: PBASet,
}

impl Config {
    pub fn new() -> Config {
        Config {
            program_set: ProgramSet::new(),
            mm_set: MMSet::new(),
            pbe_set: PBESet::new(298.15),
            pba_set: PBASet::new(298.15),
        }
    }

    pub fn load<T: AsRef<Path>>(file: T) -> Result<Config, serde_yaml::Error> {
        let config_set = fs::read_to_string(&file).expect("Read Config parameters file error.");
        serde_yaml::from_str(config_set.as_str())
    }

    /// Writes a commented template with the current values, so the file can be
    /// edited and reloaded (`s_mmpbsa -c config.yaml`).
    pub fn save<T: AsRef<Path>>(&self, file: T) {
        let mut text = String::from(
            "# s_mmpbsa configuration file.  Run it with:  s_mmpbsa -c config.yaml\n\
             # Empty strings mean \"not set\"; comments describe every field.\n");
        text += "program_set:\n";
        text += &yaml_field(2, "sys_name", &yaml_scalar(&self.program_set.sys_name), "system name used by the output files");
        text += &yaml_field(2, "calc_mm", &yaml_scalar(&self.program_set.calc_mm), "calculate the molecular mechanics term");
        text += &yaml_field(2, "calc_pbsa", &yaml_scalar(&self.program_set.calc_pbsa), "calculate the PB/SA term");
        text += &yaml_field(2, "debug", &yaml_scalar(&self.program_set.debug), "keep intermediate PB/SA files");
        text += &yaml_field(2, "trj", &yaml_scalar(&self.program_set.trj), "trajectory file (xtc/trr/gro/pdb)");
        text += &yaml_field(2, "tpr", &yaml_scalar(&self.program_set.tpr), "run input file");
        text += &yaml_field(2, "ndx", &yaml_scalar(&self.program_set.ndx), "index file (generated when missing)");
        text += &yaml_field(2, "ala_scan_range", &yaml_scalar(&self.program_set.ala_scan_range), "alanine-scanning residues, e.g. \"12,30-35\"; empty = none");
        text += "mm_set:\n";
        text += &yaml_field(2, "cutoff", &yaml_scalar(&self.mm_set.cutoff), "MM distance cutoff (A); .inf = no cutoff");
        text += &yaml_field(2, "electric_screening", &yaml_scalar(&self.mm_set.electric_screening), "electrostatic screening (J. Chem. Inf. Model. 2021, 61, 2454)");
        text += &yaml_field(2, "interaction_entropy", &yaml_scalar(&self.mm_set.interaction_entropy), "interaction entropy (-TDS); turning it off speeds the run up a lot");
        text += &yaml_field(2, "rec_grp", &yaml_scalar(&self.mm_set.rec_grp), "receptor group name of the index file");
        text += &yaml_field(2, "lig_grp", &yaml_scalar(&self.mm_set.lig_grp), "ligand group name; empty = solvation-only calculation");
        text += &yaml_field(2, "start_time", &yaml_scalar(&self.mm_set.start_time), "first time to analyze (ns)");
        text += &yaml_field(2, "end_time", &yaml_scalar(&self.mm_set.end_time), "last time to analyze (ns); .inf = until the end");
        text += &yaml_field(2, "dt", &yaml_scalar(&self.mm_set.dt), "time interval between analyzed frames (ns)");
        text += &yaml_field(2, "ie_multiple", &yaml_scalar(&self.mm_set.ie_multiple), "interaction-entropy sampling rate multiple");
        text += &yaml_field(2, "fix_pbc", &yaml_scalar(&self.mm_set.fix_pbc), "fix periodic boundary conditions first");
        text += &yaml_field(2, "radius_type", &yaml_scalar(&self.mm_set.radius_type), "atom radius set: ff | amber | Bondi | mBondi | mBondi2");
        text += "pbe_set:\n";
        text += &self.pbe_set.template_body(2);
        text += "pba_set:\n";
        text += &yaml_field(2, "temp", &yaml_scalar(&self.pba_set.temp), "temperature (K)");
        text += &yaml_field(2, "srfm", &yaml_scalar(&self.pba_set.srfm), "model for the solvent surface/volume (sacc recommended)");
        text += &yaml_field(2, "swin", &yaml_scalar(&self.pba_set.swin), "cubic spline window (A)");
        text += &yaml_field(2, "srad", &yaml_scalar(&self.pba_set.srad), "probe radius (A)");
        text += &yaml_field(2, "gamma", &yaml_scalar(&self.pba_set.gamma), "surface tension (kJ/mol/A^2)");
        text += &yaml_field(2, "press", &yaml_scalar(&self.pba_set.press), "pressure for the SAV term (kJ/mol/A^3)");
        text += &yaml_field(2, "bconc", &yaml_scalar(&self.pba_set.bconc), "solvent bulk density (A^-3)");
        text += &yaml_field(2, "sdens", &yaml_scalar(&self.pba_set.sdens), "surface point density");
        text += &yaml_field(2, "dpos", &yaml_scalar(&self.pba_set.dpos), "spacing of the volume-integration grid (A)");
        text += &yaml_field(2, "grid", &format!("[{}, {}, {}]", self.pba_set.grid.0, self.pba_set.grid.1, self.pba_set.grid.2), "explicit apolar integration grid spacing (A); keep [0.1, 0.1, 0.1]");
        text += &yaml_field(2, "calc_force", &yaml_scalar(&self.pba_set.calc_force), "compute apolar forces");
        text += &yaml_field(2, "calc_energy", &yaml_scalar(&self.pba_set.calc_energy), "apolar energy output: no | total | comps");

        let mut f = File::create(&file).expect("Save Config parameters error.");
        f.write_all(text.as_bytes()).expect("Save Config parameters error.");
    }
}

#[derive(Serialize, Deserialize)]
pub struct ProgramSet {
    pub sys_name: String,
    pub calc_mm: bool,
    pub calc_pbsa: bool,
    pub debug: bool,
    pub trj: String,
    pub tpr: String,
    pub ndx: String,
    pub ala_scan_range: String,
}

impl ProgramSet {
    pub fn new() -> ProgramSet {
        ProgramSet {
            sys_name: "system".to_string(),
            calc_mm: true,
            calc_pbsa: true,
            debug: false,
            trj: "md.xtc".to_string(),
            tpr: "md.tpr".to_string(),
            ndx: "index.ndx".to_string(),
            ala_scan_range: String::new(),
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct MMSet {
    pub cutoff: f64,
    pub electric_screening: bool,
    pub interaction_entropy: bool,
    pub rec_grp: String,
    pub lig_grp: String,
    pub start_time: f64,
    pub end_time: f64,
    pub dt: f64,
    pub ie_multiple: usize,
    pub fix_pbc: bool,
    pub radius_type: String,
}

impl MMSet {
    pub fn new() -> MMSet {
        MMSet {
            cutoff: f64::INFINITY,
            electric_screening: true,
            interaction_entropy: true,
            rec_grp: String::from("Protein"),
            lig_grp: String::new(),
            start_time: 0.0,
            end_time: f64::INFINITY,
            dt: 1.0,
            ie_multiple: 10,
            fix_pbc: true,
            radius_type: "mBondi".to_string(),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct PBESet {
    pub temp: f64,
    pub pdie: f64,
    pub sdie: f64,
    pub pb_solver: PBSolver,
    pub bcfl: Bcfl,
    pub srfm: Srfm,
    pub chgm: Chgm,
    pub swin: f64,
    pub srad: f64,
    pub sdens: i32,
    pub ions: Vec<Ion>,
    pub calc_energy: CalcEnergy,
    /// Grid knobs of the `-c config.yaml` path: they are synced into the
    /// program settings (`settings.cfac/fadd/df`) when a config is loaded.
    /// In the interactive menu they come from menu items 6/7/8 instead.
    pub cfac: f64,
    pub fadd: f64,
    pub df: f64,
}

impl PBESet {
    pub fn new(temp: f64) -> PBESet {
        return PBESet {
            temp,
            pdie: 2.0,
            sdie: 78.4,
            pb_solver: PBSolver::Lpbe,
            bcfl: Bcfl::Mdh,
            srfm: Srfm::Smol,
            chgm: Chgm::Spl4,
            swin: 0.3,
            srad: 1.4,
            sdens: 10,
            ions: vec![
                Ion { charge: 1.0, conc: 0.15, radius: 0.95 },
                Ion { charge: -1.0, conc: 0.15, radius: 1.81 },
            ],
            calc_energy: CalcEnergy::Comps,
            cfac: 1.5,
            fadd: 5.0,
            df: 0.5,
        };
    }

    pub fn from(pbe_set: &PBESet) -> PBESet {
        let mut ions: Vec<Ion> = vec![];
        for ion in &pbe_set.ions {
            ions.push(ion.clone());
        }
        let new_pbe_set = PBESet {
            temp: pbe_set.temp,
            pdie: pbe_set.pdie,
            sdie: pbe_set.sdie,
            pb_solver: pbe_set.pb_solver,
            bcfl: pbe_set.bcfl,
            srfm: pbe_set.srfm,
            chgm: pbe_set.chgm,
            swin: pbe_set.swin,
            srad: pbe_set.srad,
            sdens: pbe_set.sdens,
            ions,
            calc_energy: pbe_set.calc_energy,
            cfac: pbe_set.cfac,
            fadd: pbe_set.fadd,
            df: pbe_set.df,
        };
        return new_pbe_set;
    }

    pub fn load<T: AsRef<Path>>(file: T) -> Result<PBESet, serde_yaml::Error> {
        let pbe_set = fs::read_to_string(&file).expect("Read PB parameters file error.");
        serde_yaml::from_str(pbe_set.as_str())
    }

    /// Writes a commented template with the current values, so the file can be
    /// edited and reloaded from the interactive menu.
    pub fn save<T: AsRef<Path>>(&self, file: T) {
        let text = format!(
            "# PB (Poisson-Boltzmann) parameters for s_mmpbsa.\n\
             # Edit the values, save the file, then press ENTER in s_mmpbsa to reload.\n{}",
            self.template_body(0)
        );
        let mut f = File::create(&file).expect("Save PB parameters error.");
        f.write_all(text.as_bytes()).expect("Save PB parameters error.");
    }

    /// The commented `key: value` lines of the PB parameters, indented by
    /// `indent` spaces.  Shared by `PBESet::save` and the `config.yaml`
    /// template of `Config::save`.
    fn template_body(&self, indent: usize) -> String {
        let mut s = String::new();
        s += &yaml_field(indent, "temp", &self.temp.to_string(), "temperature (K)");
        s += &yaml_field(indent, "pdie", &self.pdie.to_string(), "solute dielectric constant");
        s += &yaml_field(indent, "sdie", &self.sdie.to_string(), "solvent dielectric constant (water: 78.4 at 298.15 K, vacuum: 1)");
        s += &yaml_field(indent, "pb_solver", &self.pb_solver.to_string(), "PB equation: lpbe (linear, recommended) | npbe (nonlinear)");
        s += &yaml_field(indent, "bcfl", &self.bcfl.to_string(), "boundary condition: zero | sdh | mdh (recommended) | focus | mem | map");
        s += &yaml_field(indent, "srfm", &self.srfm.to_string(), "surface definition: mol | smol (recommended) | spline | spline3 | spline4 | sacc");
        s += &yaml_field(indent, "chgm", &self.chgm.to_string(), "charge mapping: tril | spl2 | spl4 (recommended)");
        s += &yaml_field(indent, "swin", &self.swin.to_string(), "spline window (A), only for spline/spline3/spline4 surfaces");
        s += &yaml_field(indent, "srad", &self.srad.to_string(), "solvent probe radius (A)");
        s += &yaml_field(indent, "sdens", &self.sdens.to_string(), "surface density (grid points/A^2), for mol/sacc surfaces");
        s += &yaml_field(indent, "calc_energy", &self.calc_energy.to_string(), "energy output: no | total | comps (comps is what s_mmpbsa needs)");
        s += &yaml_field(indent, "ions", "", "mobile ions (charge in e, conc in M, radius in A), one block per species");
        for ion in &self.ions {
            let pad = " ".repeat(indent);
            s += &format!("{pad}- charge: {}\n", ion.charge);
            s += &format!("{pad}  conc: {}\n", ion.conc);
            s += &format!("{pad}  radius: {}\n", ion.radius);
        }
        s += &yaml_field(indent, "cfac", &self.cfac.to_string(), "coarse grid expansion factor (with -c config; interactive: menu 6)");
        s += &yaml_field(indent, "fadd", &self.fadd.to_string(), "fine grid padding (A) (with -c config; interactive: menu 7)");
        s += &yaml_field(indent, "df", &self.df.to_string(), "fine mesh spacing (A) (with -c config; interactive: menu 8)");
        s
    }
}

impl fmt::Display for PBESet {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut ions = String::new();
        for ion in &self.ions {
            ions.push_str(format!("  ion {}\n", ion).as_str());
        }
        write!(f,
            "  temp  {:7}  # Temperature\
            \n  pdie  {:7}  # Solute dielectric constant\
            \n  sdie  {:7}  # Solvent dielectric constant, vacuum 1, water 78.4 (298.15 K)\
            \n  \
            \n  {}           # PB equation solving method, lpbe(linear), npbe(nonlinear)\
            \n  bcfl  {:>7}  # Boundary conditions for coarse-grid PB equation, zero, sdh/mdh(single/multiple Debye-Huckel), focus, map\
            \n  srfm  {:>7}  # Model for constructing dielectric and ion boundaries, mol(molecular surface), smol(smooth molecular surface), spline/spline3/spline4\
            \n  chgm  {:>7}  # Charge mapping to grid points method, tril(trilinear interpolation), spl2/spl4(cubic/quartic B-spline discretization)\
            \n  swin  {:7}  # Cubic spline window value, only used for spline surfaces\
            \n  \
            \n  srad  {:7}  # Solvent probe radius\
            \n  sdens {:7}  # Surface density, grid points per A^2, not used when (srad=0) or (srfm=spline*)\
            \n  \
            \n  # Ion charge, concentration, radius\
            \n{}  \
            \n  calcenergy {}",
                self.temp, self.pdie, self.sdie,
                self.pb_solver, self.bcfl, self.srfm, self.chgm, self.swin,
                self.srad, self.sdens, ions, self.calc_energy
        )
    }
}

impl Clone for PBESet {
    fn clone(&self) -> PBESet {
        PBESet::from(self)
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct Ion {
    pub charge: f64,
    pub conc: f64,
    radius: f64,
}

impl fmt::Display for Ion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "charge {:2} conc {} radius {}", self.charge, self.conc, self.radius)
    }
}

impl Copy for Ion {}

impl Clone for Ion {
    fn clone(&self) -> Self {
        *self
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PBASet {
    temp: f64,
    srfm: String,
    swin: f64,
    srad: f64,
    pub gamma: f64,
    press: f64,
    pub bconc: f64,
    sdens: f64,
    dpos: f64,
    grid: (f64, f64, f64),
    calc_force: bool,
    calc_energy: String,
}

impl PBASet {
    pub fn new(temp: f64) -> Self {
        PBASet {
            temp,
            srfm: "sacc".to_string(),
            swin: 0.3,
            srad: 1.4,
            gamma: 0.0226778,
            press: 0.0,
            bconc: 3.84982,
            sdens: 10.0,
            dpos: 0.2,
            grid: (0.1, 0.1, 0.1),
            calc_force: false,
            calc_energy: "comps".to_string(),
        }
    }

    pub fn from(pba_set: &PBASet) -> PBASet {
        PBASet {
            temp: pba_set.temp,
            srfm: pba_set.srfm.to_string(),
            swin: pba_set.swin,
            srad: pba_set.srad,
            gamma: pba_set.gamma,
            press: pba_set.press,
            bconc: pba_set.bconc,
            sdens: pba_set.sdens,
            dpos: pba_set.dpos,
            grid: pba_set.grid,
            calc_force: pba_set.calc_force,
            calc_energy: pba_set.calc_energy.to_string(),
        }
    }

    pub fn load<T: AsRef<Path>>(file: T) -> Result<PBASet, serde_yaml::Error> {
        let pba_set = fs::read_to_string(&file).expect("Read SA parameters file error.");
        serde_yaml::from_str(pba_set.as_str())
    }

    pub fn save<T: AsRef<Path>>(&self, file: T) {
        let mut f = File::create(&file).expect("Save SA parameters error.");
        f.write_all(serde_yaml::to_string(self).unwrap().as_bytes()).expect("Save SA parameters error.");
    }
}

impl fmt::Display for PBASet {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f,
            "  temp  {:7}  # Temperature\
            \n  srfm  {:>7}  # Model for constructing solvent-related surface or volume\
            \n  swin  {:7}  # Cubic spline window (A), used to define spline surface\
            \n  \
            \n  srad  {:7}  # Probe radius (A)\
            \n  gamma {:7}  # Surface tension (kJ/mol-A^2)\
            \n  \
            \n  press {:7}  # Pressure (kJ/mol-A^3)\
            \n  bconc {:7}  # Solvent bulk density (A^3)\
            \n  sdens {:7}\
            \n  dpos  {:7}\
            \n  grid  {:7} {:5} {:5}\
            \n  \
            \n  calcforce  {}\
            \n  calcenergy {}", self.temp, self.srfm, self.swin,
                self.srad, self.gamma, self.press, self.bconc, self.sdens,
                self.dpos, self.grid.0, self.grid.1, self.grid.2,
                match self.calc_force {
                    true => "yes",
                    false => "no"
                }, self.calc_energy
        )
    }
}

impl Clone for PBASet {
    fn clone(&self) -> PBASet {
        PBASet::from(self)
    }
}

/// One `key: value  # comment` line of a YAML template.
fn yaml_field(indent: usize, key: &str, value: &str, comment: &str) -> String {
    let kv = format!("{}{}: {}", " ".repeat(indent), key, value);
    let trimmed = kv.trim_end();
    if comment.is_empty() {
        format!("{trimmed}\n")
    } else if trimmed.len() < 40 {
        format!("{trimmed:<40}# {comment}\n")
    } else {
        format!("{trimmed}  # {comment}\n")
    }
}

/// Serializes one scalar with proper YAML quoting.
fn yaml_scalar<T: serde::Serialize>(value: &T) -> String {
    serde_yaml::to_string(value).unwrap().trim_end().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pbe_from(text: &str) -> Result<PBESet, serde_yaml::Error> {
        serde_yaml::from_str(text)
    }

    #[test]
    fn yaml_roundtrip_preserves_pb_values() {
        let mut pbe = PBESet::new(310.15);
        pbe.df = 0.75;
        pbe.bcfl = Bcfl::Focus;
        let yaml = serde_yaml::to_string(&pbe).unwrap();
        let back: PBESet = serde_yaml::from_str(&yaml).unwrap();
        assert_eq!(back, pbe);
    }

    /// The generated template must deserialize back to the same parameters —
    /// a template that cannot be loaded would lock the user out.
    #[test]
    fn pb_template_roundtrips() {
        let mut pbe = PBESet::new(298.15);
        pbe.pb_solver = PBSolver::Npbe;
        pbe.srfm = Srfm::Spline4;
        let text = pbe.template_body(0);
        assert!(text.contains("pb_solver: npbe"));
        assert!(text.contains("srfm: spline4"));
        let back: PBESet = pbe_from(&text).unwrap();
        assert_eq!(back, pbe);
    }

    /// Old settings files spelled the surface/charge models spl2/spl4; they
    /// must keep loading and map onto the canonical solver tokens.
    #[test]
    fn legacy_surface_and_charge_names_are_accepted() {
        let text = "temp: 298.15\npdie: 2.0\nsdie: 78.4\npb_solver: lpbe\nbcfl: mdh\n\
                    srfm: spl4\nchgm: spl4\nswin: 0.3\nsrad: 1.4\nsdens: 10\ncalc_energy: comps\n\
                    ions:\n- charge: 1.0\n  conc: 0.15\n  radius: 0.95\n\
                    cfac: 1.5\nfadd: 5.0\ndf: 0.5\n";
        let pbe: PBESet = pbe_from(text).unwrap();
        assert_eq!(pbe.srfm, Srfm::Spline4);
        assert_eq!(pbe.chgm, Chgm::Spl4);
    }

    /// A misspelled choice must fail at load time with the list of valid
    /// values, not turn into a per-frame solver failure.
    #[test]
    fn misspelled_choice_is_rejected_with_valid_values() {
        let text = "temp: 298.15\npdie: 2.0\nsdie: 78.4\npb_solver: lpbe\nbcfl: mdhh\n\
                    srfm: smol\nchgm: spl4\nswin: 0.3\nsrad: 1.4\nsdens: 10\ncalc_energy: comps\n\
                    ions: []\ncfac: 1.5\nfadd: 5.0\ndf: 0.5\n";
        let err = pbe_from(text).unwrap_err().to_string();
        assert!(err.contains("unknown variant"), "{err}");
        assert!(err.contains("mdh"), "{err}");
    }

    #[test]
    fn presets_change_only_their_own_fields() {
        let mut pbe = PBESet::new(300.0);
        pbe.pdie = 4.0;
        PbePreset::FastScreening.apply(&mut pbe);
        assert_eq!(pbe.df, 1.0);
        assert_eq!(pbe.pdie, 4.0, "preset must not touch unrelated fields");

        PbePreset::Nonlinear.apply(&mut pbe);
        assert_eq!(pbe.pb_solver, PBSolver::Npbe);

        PbePreset::Default.apply(&mut pbe);
        assert_eq!(pbe.df, 0.5);
        assert_eq!(pbe.pb_solver, PBSolver::Lpbe);
        assert_eq!(pbe.temp, 300.0, "default preset keeps the temperature");
    }

    /// The whole `-c config.yaml` template must deserialize back into the
    /// identical configuration.
    #[test]
    fn config_save_writes_loadable_template() {
        let mut config = Config::new();
        config.program_set.sys_name = "kaguya".to_string();
        config.mm_set.rec_grp = "Protein".to_string();
        config.pbe_set.df = 0.4;
        let path = std::env::temp_dir().join("s_mmpbsa_config_template_test.yaml");
        config.save(&path);
        let text = fs::read_to_string(&path).unwrap();
        assert!(text.contains("sys_name: kaguya"));
        assert!(text.contains("# boundary condition"));
        let back = Config::load(&path).unwrap();
        fs::remove_file(&path).ok();
        assert_eq!(back.program_set.sys_name, "kaguya");
        assert_eq!(back.pbe_set.df, 0.4);
        assert_eq!(back.pbe_set, config.pbe_set);
        assert!((back.pba_set.gamma - config.pba_set.gamma).abs() < 1e-12);
        assert!((back.pba_set.bconc - config.pba_set.bconc).abs() < 1e-12);
        assert_eq!(back.mm_set.rec_grp, "Protein");
    }
}
