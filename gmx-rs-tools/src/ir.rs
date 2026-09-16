//! Full decoding of the `t_inputrec` section of a tpr file.
//!
//! This mirrors `do_inputrec()` from `fileio/tpxio.cpp` (plus the helpers it
//! calls: `do_fepvals`, `do_simtempvals`, `do_expandedvals`, `do_pull`,
//! `do_rot`, `do_imd`, `do_swapcoords_tpx`, `do_legacy_efield`) and the
//! `KeyValueTree` serialization from `serialization/keyvaluetreeserializer.cpp`
//! and `applied_forces/awh/read_params.cpp`.
//!
//! Only the fields that `gmx dump` prints are kept, but every field that is
//! stored in the file has to be decoded, since the following field starts where
//! the previous one ends.

use crate::tpr::{CReader, TpxHeader};
use crate::xdr::{Result, XdrError};

// ---------------------------------------------------------------------------
// tpx versions that changed the inputrec layout (see the `tpxv` enum).
// ---------------------------------------------------------------------------
pub const TPXV_V51: i32 = 51;
pub const TPXV_V53: i32 = 53;
pub const TPXV_V59: i32 = 59;
pub const TPXV_V60: i32 = 60;
pub const TPXV_V62: i32 = 62;
pub const TPXV_V64: i32 = 64;
pub const TPXV_V67: i32 = 67;
pub const TPXV_V69: i32 = 69;
pub const TPXV_V71: i32 = 71;
pub const TPXV_V73: i32 = 73;
pub const TPXV_V74: i32 = 74;
pub const TPXV_V77: i32 = 77;
pub const TPXV_V79: i32 = 79;
pub const TPXV_V81: i32 = 81;
pub const TPXV_V82: i32 = 82;
pub const TPXV_V83: i32 = 83;
pub const TPXV_V90: i32 = 90;
pub const TPXV_V92: i32 = 92;
pub const TPXV_V93: i32 = 93;
pub const TPXV_V94: i32 = 94;
pub const TPXV_V95: i32 = 95;
pub const TPXV_V96: i32 = 96; // ComputationalElectrophysiology
pub const TPXV_V97: i32 = 97; // Use64BitRandomSeed
pub const TPXV_V99: i32 = 99; // InteractiveMolecularDynamics
pub const TPXV_V100: i32 = 100; // RemoveObsoleteParameters1
pub const TPXV_V101: i32 = 101; // PullCoordTypeGeom
pub const TPXV_V102: i32 = 102; // PullGeomDirRel
pub const TPXV_V104: i32 = 104; // CompElWithSwapLayerOffset
pub const TPXV_V105: i32 = 105; // CompElPolyatomicIonsAndMultipleIonTypes
pub const TPXV_V106: i32 = 106; // RemoveAdress
pub const TPXV_V107: i32 = 107; // PullCoordNGroup
pub const TPXV_V108: i32 = 108; // RemoveTwinRange
pub const TPXV_V109: i32 = 109; // ReplacePullPrintCOM12
pub const TPXV_V110: i32 = 110; // PullExternalPotential
pub const TPXV_V111: i32 = 111; // GenericParamsForElectricField
pub const TPXV_V112: i32 = 112; // AcceleratedWeightHistogram
pub const TPXV_V113: i32 = 113; // RemoveImplicitSolvation
pub const TPXV_V114: i32 = 114; // PullPrevStepCOMAsReference
pub const TPXV_V116: i32 = 116; // PullAverage
pub const TPXV_V117: i32 = 117; // GenericInternalParameters
pub const TPXV_V122: i32 = 122; // MTS
pub const TPXV_V124: i32 = 124; // TransformationPullCoord
pub const TPXV_V125: i32 = 125; // SoftcoreGapsys
pub const TPXV_V129: i32 = 129; // EnsembleTemperature
pub const TPXV_V130: i32 = 130; // AwhGrowthFactor
pub const TPXV_V131: i32 = 131; // MassRepartitioning
pub const TPXV_V132: i32 = 132; // AwhTargetMetricScaling
pub const TPXV_V133: i32 = 133; // VerletBufferPressureTol
pub const TPXV_V135: i32 = 135; // RefScaleMultipleCOMs
pub const TPXV_V136: i32 = 136; // InputHistogramCounts
pub const TPXV_V138: i32 = 138; // AwhHistogramTolerance

/// A value of a serialized `gmx::KeyValueTreeObject`.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum KvtValue {
    #[default]
    Null,
    Object(Vec<(String, KvtValue)>),
    Array(Vec<KvtValue>),
    Str(String),
    Bool(bool),
    Char(i8),
    UChar(u8),
    Int(i32),
    Int64(i64),
    Float(f32),
    Double(f64),
}

impl KvtValue {
    /// `gmx::simpleValueToString()`: the text `dumpKeyValueTree` writes for a
    /// leaf value.
    pub fn to_text(&self) -> String {
        match self {
            KvtValue::Str(s) => s.clone(),
            KvtValue::Bool(b) => if *b { "true" } else { "false" }.to_string(),
            KvtValue::Char(c) => format!("0x{:x}", *c as u8),
            KvtValue::UChar(c) => format!("{c:x}"),
            KvtValue::Int(i) => format!("{i}"),
            KvtValue::Int64(i) => format!("{i}"),
            KvtValue::Float(f) => crate::cmd::fmt_g(*f as f64, 0, 6),
            KvtValue::Double(d) => crate::cmd::fmt_g(*d, 0, 6),
            KvtValue::Object(_) | KvtValue::Array(_) => String::new(),
            KvtValue::Null => String::new(),
        }
    }

    pub fn as_object(&self) -> Option<&Vec<(String, KvtValue)>> {
        match self {
            KvtValue::Object(o) => Some(o),
            _ => None,
        }
    }

    /// Looks up `key` in an object.
    pub fn get(&self, key: &str) -> Option<&KvtValue> {
        self.as_object()?
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v)
    }
}

/// Reads a serialized `KeyValueTreeObject` (a count followed by
/// `doString(key)`/value pairs).
fn read_kvt_object(r: &mut CReader) -> Result<KvtValue> {
    let count = r.int()?;
    if count < 0 {
        return Err(XdrError::Invalid("negative key value tree size".into()));
    }
    let mut props = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let key = r.string()?;
        let value = read_kvt_value(r)?;
        props.push((key, value));
    }
    Ok(KvtValue::Object(props))
}

fn read_kvt_value(r: &mut CReader) -> Result<KvtValue> {
    // `'O'`, `'A'`, `'s'`, `'b'`, `'c'`, `'u'`, `'i'`, `'l'`, `'f'`, `'d'`.
    let tag = r.uchar()?;
    let value = match tag {
        b'O' => read_kvt_object(r)?,
        b'A' => {
            let count = r.int()?;
            if count < 0 {
                return Err(XdrError::Invalid("negative key value array size".into()));
            }
            let mut values = Vec::with_capacity(count as usize);
            for _ in 0..count {
                values.push(read_kvt_value(r)?);
            }
            KvtValue::Array(values)
        }
        b's' => KvtValue::Str(r.string()?),
        b'b' => KvtValue::Bool(r.bool()?),
        b'c' => KvtValue::Char(r.uchar()? as i8),
        b'u' => KvtValue::UChar(r.uchar()?),
        b'i' => KvtValue::Int(r.int()?),
        b'l' => KvtValue::Int64(r.int64()?),
        b'f' => KvtValue::Float(r.float()?),
        b'd' => KvtValue::Double(r.double()?),
        other => {
            return Err(XdrError::Invalid(format!(
                "unknown key value tree type tag {other}"
            )))
        }
    };
    Ok(value)
}

#[derive(Debug, Clone, Default)]
pub struct MtsLevel {
    pub force_groups: i32,
    pub step_factor: i32,
}

#[derive(Debug, Clone, Default)]
pub struct FepVals {
    pub init_fep_state: i32,
    pub init_lambda_without_states: f64,
    pub delta_lambda: f64,
    pub n_lambda: i32,
    pub all_lambda: [Vec<f64>; 7],
    pub separate_dvdl: [bool; 7],
    pub sc_alpha: f64,
    pub sc_power: i32,
    pub sc_r_power: f64,
    pub sc_sigma: f64,
    pub sc_sigma_min: f64,
    pub b_sc_coul: bool,
    pub nstdhdl: i32,
    pub separate_dhdl_file: i32,
    pub dhdl_derivatives: i32,
    pub dh_hist_size: i32,
    pub dh_hist_spacing: f64,
    pub edhdl_print_energy: i32,
    pub softcore_function: i32,
    pub sc_gapsys_scale_linpoint_lj: f64,
    pub sc_gapsys_scale_linpoint_q: f64,
    pub sc_gapsys_sigma_lj: f64,
    pub lambda_neighbors: i32,
    pub lambda_start_n: i32,
    pub lambda_stop_n: i32,
}

#[derive(Debug, Clone, Default)]
pub struct SimTempVals {
    pub scale: i32,
    pub high: f64,
    pub low: f64,
    pub temperatures: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct ExpandedVals {
    pub init_lambda_weights: Vec<f64>,
    pub nstexpanded: i32,
    pub elmcmove: i32,
    pub elamstats: i32,
    pub lmc_repeats: i32,
    pub gibbsdeltalam: i32,
    pub lmc_forced_nstart: i32,
    pub lmc_seed: i32,
    pub mc_temp: f64,
    pub b_symmetrized_t_matrix: bool,
    pub nst_tij: i32,
    pub minvarmin: i32,
    pub c_range: i32,
    pub wl_scale: f64,
    pub wl_ratio: f64,
    pub init_wl_delta: f64,
    pub b_wl_oneovert: bool,
    pub elmceq: i32,
    pub equil_steps: i32,
    pub equil_samples: i32,
    pub equil_n_at_lam: i32,
    pub equil_wl_delta: f64,
    pub equil_ratio: f64,
    pub init_lambda_counts: Vec<f64>,
    pub init_wl_histogram_counts: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct PullGroup {
    pub ind: Vec<i32>,
    pub weight: Vec<f64>,
    pub pbcatom: i32,
}

#[derive(Debug, Clone, Default)]
pub struct PullCoord {
    pub etype: i32,
    pub external_potential_provider: String,
    pub egeom: i32,
    pub ngroup: i32,
    pub group: Vec<i32>,
    pub dim: [i32; 3],
    pub expression: String,
    pub origin: [f64; 3],
    pub vec: [f64; 3],
    pub b_start: bool,
    pub init: f64,
    pub rate: f64,
    pub k: f64,
    pub kb: f64,
}

#[derive(Debug, Clone, Default)]
pub struct PullParams {
    pub ngroup: i32,
    pub ncoord: i32,
    pub cylinder_r: f64,
    pub constr_tol: f64,
    pub b_print_com: bool,
    pub b_print_ref_value: bool,
    pub b_print_comp: bool,
    pub nstxout: i32,
    pub nstfout: i32,
    pub b_set_pbc_ref_to_prev_step_com: bool,
    pub group: Vec<PullGroup>,
    pub coord: Vec<PullCoord>,
    pub b_xout_average: bool,
    pub b_fout_average: bool,
}

#[derive(Debug, Clone, Default)]
pub struct AwhDimParams {
    pub coord_provider: i32,
    pub coord_index: i32,
    pub origin: f64,
    pub end: f64,
    pub period: f64,
    pub force_constant: f64,
    pub diffusion: f64,
    pub coord_value_init: f64,
    pub cover_diameter: f64,
}

#[derive(Debug, Clone, Default)]
pub struct AwhBiasParams {
    pub e_target: i32,
    pub target_beta_scaling: f64,
    pub target_cutoff: f64,
    pub e_growth: i32,
    pub growth_factor: f64,
    pub b_user_data: bool,
    pub scale_target_by_metric: bool,
    pub target_metric_scaling_limit: f64,
    pub error_init: f64,
    pub share_group: i32,
    pub equilibrate_histogram: bool,
    pub histogram_tolerance: f64,
    pub dim_params: Vec<AwhDimParams>,
}

#[derive(Debug, Clone, Default)]
pub struct AwhParams {
    pub nstout: i32,
    pub seed: i64,
    pub nst_sample_coord: i32,
    pub num_samples_update_free_energy: i32,
    pub potential: i32,
    pub share_bias_multisim: bool,
    pub bias: Vec<AwhBiasParams>,
}

#[derive(Debug, Clone, Default)]
pub struct RotGroup {
    pub e_type: i32,
    pub b_mass_w: i32,
    pub ind: Vec<i32>,
    pub x_ref_original: Vec<[f64; 3]>,
    pub input_vec: [f64; 3],
    pub pivot: [f64; 3],
    pub rate: f64,
    pub k: f64,
    pub slab_dist: f64,
    pub min_gaussian: f64,
    pub eps: f64,
    pub e_fittype: i32,
    pub pot_angle_nstep: i32,
    pub pot_angle_step: f64,
}

#[derive(Debug, Clone, Default)]
pub struct Rot {
    pub nstrout: i32,
    pub nstsout: i32,
    pub grp: Vec<RotGroup>,
}

#[derive(Debug, Clone, Default)]
pub struct Imd {
    pub ind: Vec<i32>,
}

#[derive(Debug, Clone, Default)]
pub struct SwapGroup {
    pub molname: String,
    pub ind: Vec<i32>,
    pub nmol_req: [i32; 2],
}

#[derive(Debug, Clone, Default)]
pub struct SwapCoords {
    pub groups: Vec<SwapGroup>,
    pub massw_split: [bool; 2],
    pub nstswap: i32,
    pub n_average: i32,
    pub threshold: f64,
    pub cyl0r: f64,
    pub cyl0u: f64,
    pub cyl0l: f64,
    pub cyl1r: f64,
    pub cyl1u: f64,
    pub cyl1l: f64,
    pub bulk_offset: [f64; 2],
}

#[derive(Debug, Clone, Default)]
pub struct GrpOpts {
    pub ngtc: i32,
    pub nhchainlength: i32,
    pub ngfrz: i32,
    pub ngener: i32,
    pub nrdf: Vec<f64>,
    pub ref_t: Vec<f64>,
    pub tau_t: Vec<f64>,
    pub nfreeze: Vec<[i32; 3]>,
    pub acceleration: Vec<[f64; 3]>,
    pub egp_flags: Vec<i32>,
    pub annealing: Vec<i32>,
    pub anneal_npoints: Vec<i32>,
    pub anneal_time: Vec<Vec<f64>>,
    pub anneal_temp: Vec<Vec<f64>>,
    pub ngqm: i32,
    pub b_qmmm: bool,
}

/// The decoded `t_inputrec`.
#[derive(Debug, Clone, Default)]
pub struct InputRec {
    pub pbc_type: i32,
    pub b_periodic_mols: bool,

    pub integrator: i32,
    pub nsteps: i64,
    pub init_step: i64,
    pub simulation_part: i32,
    pub use_mts: bool,
    pub mts_levels: Vec<MtsLevel>,
    pub mass_repartition_factor: f64,
    pub ensemble_temperature_setting: i32,
    pub ensemble_temperature: f64,
    pub nstcalcenergy: i32,
    pub cutoff_scheme: i32,
    pub nstlist: i32,
    pub rtpi: f64,
    pub nstcomm: i32,
    pub comm_mode: i32,
    pub nstcgsteep: i32,
    pub nbfgscorr: i32,
    pub nstlog: i32,
    pub nstxout: i32,
    pub nstvout: i32,
    pub nstfout: i32,
    pub nstenergy: i32,
    pub nstxout_compressed: i32,
    pub init_t: f64,
    pub delta_t: f64,
    pub x_compression_precision: f64,
    pub verletbuf_tol: f64,
    pub verlet_buffer_pressure_tolerance: f64,
    pub rlist: f64,
    pub coulombtype: i32,
    pub coulomb_modifier: i32,
    pub rcoulomb_switch: f64,
    pub rcoulomb: f64,
    pub vdwtype: i32,
    pub vdw_modifier: i32,
    pub rvdw_switch: f64,
    pub rvdw: f64,
    pub disp_corr: i32,
    pub epsilon_r: f64,
    pub epsilon_rf: f64,
    pub tabext: f64,
    pub fourier_spacing: f64,
    pub nkx: i32,
    pub nky: i32,
    pub nkz: i32,
    pub pme_order: i32,
    pub ewald_rtol: f64,
    pub ewald_rtol_lj: f64,
    pub ljpme_combination_rule: i32,
    pub ewald_geometry: i32,
    pub epsilon_surface: f64,
    pub b_continuation: bool,
    pub tcoupl: i32,
    pub b_print_nh_chains: bool,
    pub nsttcouple: i32,
    pub pcoupl: i32,
    pub pcoupltype: i32,
    pub nstpcouple: i32,
    pub tau_p: f64,
    pub ref_p: [[f64; 3]; 3],
    pub compress: [[f64; 3]; 3],
    pub refcoord_scaling: i32,
    pub posres_com: Vec<[f64; 3]>,
    pub posres_com_b: Vec<[f64; 3]>,
    /// Whether `posresCom` was allocated while reading.  GROMACS keeps the
    /// buffer when it clears an "old style" single COM group, which is visible
    /// in the dump as an empty `(0x3)` block instead of "not available".
    pub posres_com_allocated: bool,
    pub shake_tol: f64,
    pub efep: i32,
    pub fepvals: FepVals,
    pub b_sim_temp: bool,
    pub simtempvals: SimTempVals,
    pub b_expanded: bool,
    pub expandedvals: ExpandedVals,
    pub e_disre: i32,
    pub e_disre_weighting: i32,
    pub b_disre_mixed: bool,
    pub dr_fc: f64,
    pub dr_tau: f64,
    pub nstdisreout: i32,
    pub orires_fc: f64,
    pub orires_tau: f64,
    pub nstorireout: i32,
    pub em_stepsize: f64,
    pub em_tol: f64,
    pub b_shake_sor: bool,
    pub niter: i32,
    pub fc_stepsize: f64,
    pub e_constr_alg: i32,
    pub n_proj_order: i32,
    pub lincs_warn_angle: f64,
    pub n_lincs_iter: i32,
    pub bd_fric: f64,
    pub ld_seed: i64,
    pub deform: [[f64; 3]; 3],
    pub cos_accel: f64,
    pub userint1: i32,
    pub userint2: i32,
    pub userint3: i32,
    pub userint4: i32,
    pub userreal1: f64,
    pub userreal2: f64,
    pub userreal3: f64,
    pub userreal4: f64,
    pub b_pull: bool,
    pub pull: PullParams,
    pub b_do_awh: bool,
    pub awh: AwhParams,
    pub b_rot: bool,
    pub rot: Rot,
    pub b_imd: bool,
    pub imd: Imd,
    pub opts: GrpOpts,
    pub nwall: i32,
    pub wall_type: i32,
    pub wall_r_linpot: f64,
    pub wall_atomtype: [i32; 2],
    pub wall_density: [f64; 2],
    pub wall_ewald_zfac: f64,
    pub e_swap_coords: i32,
    pub swap: SwapCoords,
    /// Legated electric field parameters, used for tpr versions that predate
    /// the key value tree.
    pub legacy_efield: Vec<[f64; 4]>,
    /// `params`: the module parameters stored as a key value tree.
    pub params: KvtValue,
    pub internal_parameters: KvtValue,
}

impl InputRec {
    pub fn pbc(&self) -> crate::frame::PbcType {
        crate::frame::PbcType::from_int(self.pbc_type)
    }
}

/// Tracing helper, enabled with `GMXRS_TRACE`.
fn tr(r: &CReader, msg: &str) {
    if std::env::var_os("GMXRS_TRACE").is_some() {
        eprintln!("[ir] {msg}: pos={}", r.position());
    }
}

fn rvec(r: &mut CReader) -> Result<[f64; 3]> {
    Ok([r.real()?, r.real()?, r.real()?])
}

fn rvec_array(r: &mut CReader, n: usize) -> Result<Vec<[f64; 3]>> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(rvec(r)?);
    }
    Ok(v)
}

fn ivec(r: &mut CReader) -> Result<[i32; 3]> {
    Ok([r.int()?, r.int()?, r.int()?])
}

fn do_fepvals(r: &mut CReader, fep: &mut FepVals, version: i32) -> Result<()> {
    if version >= TPXV_V79 {
        fep.init_fep_state = r.int()?;
        fep.init_lambda_without_states = r.double()?;
        fep.delta_lambda = r.double()?;
    } else if version >= TPXV_V59 {
        fep.init_lambda_without_states = r.double()?;
        fep.delta_lambda = r.double()?;
    } else {
        fep.init_lambda_without_states = r.real()?;
        fep.delta_lambda = r.real()?;
    }
    if version >= TPXV_V79 {
        fep.n_lambda = r.int()?;
        for g in 0..7 {
            if fep.n_lambda > 0 {
                fep.all_lambda[g] = r.real_array(fep.n_lambda as usize)?;
                // `doBoolArray()` reads the whole `separate_dvdl` array on
                // every iteration of the loop above, so the file stores the
                // seven bools seven times over.
                for slot in fep.separate_dvdl.iter_mut() {
                    *slot = r.bool()?;
                }
            } else if fep.init_lambda_without_states >= 0.0 {
                fep.separate_dvdl[0] = true;
            }
        }
    } else if version >= TPXV_V64 {
        fep.n_lambda = r.int()?;
        let lambda = r.real_array(fep.n_lambda.max(0) as usize)?;
        for g in 0..7 {
            fep.all_lambda[g] = lambda.clone();
        }
        if fep.init_lambda_without_states >= 0.0 {
            fep.separate_dvdl[0] = true;
        }
    } else {
        fep.n_lambda = 0;
        if fep.init_lambda_without_states >= 0.0 {
            fep.separate_dvdl[0] = true;
        }
    }
    fep.sc_alpha = r.real()?;
    fep.sc_power = r.int()?;
    if version >= TPXV_V79 {
        fep.sc_r_power = r.real()?;
    } else {
        fep.sc_r_power = 6.0;
    }
    fep.sc_sigma = r.real()?;
    fep.sc_sigma_min = if version >= TPXV_V71 {
        fep.sc_sigma
    } else {
        0.0
    };
    if version >= TPXV_V79 {
        fep.b_sc_coul = r.bool()?;
    } else {
        fep.b_sc_coul = true;
    }
    if version >= TPXV_V64 {
        fep.nstdhdl = r.int()?;
    } else {
        fep.nstdhdl = 1;
    }
    if version >= TPXV_V73 {
        fep.separate_dhdl_file = r.int()?;
        fep.dhdl_derivatives = r.int()?;
    } else {
        fep.separate_dhdl_file = 1;
        fep.dhdl_derivatives = 1;
    }
    if version >= TPXV_V71 {
        fep.dh_hist_size = r.int()?;
        fep.dh_hist_spacing = r.double()?;
    } else {
        fep.dh_hist_size = 0;
        fep.dh_hist_spacing = 0.1;
    }
    if version >= TPXV_V79 {
        fep.edhdl_print_energy = r.int()?;
    } else {
        fep.edhdl_print_energy = 0;
    }
    if version >= TPXV_V125 {
        fep.softcore_function = r.int()?;
        fep.sc_gapsys_scale_linpoint_lj = r.real()?;
        fep.sc_gapsys_scale_linpoint_q = r.real()?;
        fep.sc_gapsys_sigma_lj = r.real()?;
    } else {
        fep.softcore_function = 0;
        fep.sc_gapsys_scale_linpoint_lj = 0.85;
        fep.sc_gapsys_scale_linpoint_q = 0.3;
        fep.sc_gapsys_sigma_lj = 0.3;
    }
    if (version >= TPXV_V83 && version < TPXV_V90) || version >= TPXV_V92 {
        fep.lambda_neighbors = r.int()?;
        if fep.lambda_neighbors >= 0
            && fep.init_fep_state >= 0
            && fep.init_lambda_without_states < 0.0
        {
            let mut start = fep.init_fep_state - fep.lambda_neighbors;
            let mut stop = fep.init_fep_state + fep.lambda_neighbors + 1;
            if start < 0 {
                start = 0;
            }
            if stop >= fep.n_lambda {
                stop = fep.n_lambda;
            }
            fep.lambda_start_n = start;
            fep.lambda_stop_n = stop;
        } else {
            fep.lambda_start_n = 0;
            fep.lambda_stop_n = fep.n_lambda;
        }
    } else {
        fep.lambda_start_n = 0;
        fep.lambda_stop_n = fep.n_lambda;
    }
    Ok(())
}

fn do_simtempvals(
    r: &mut CReader,
    simtemp: &mut SimTempVals,
    n_lambda: i32,
    version: i32,
) -> Result<()> {
    if version >= TPXV_V79 {
        simtemp.scale = r.int()?;
        simtemp.high = r.real()?;
        simtemp.low = r.real()?;
        if n_lambda > 0 {
            simtemp.temperatures = r.real_array(n_lambda as usize)?;
        }
    }
    Ok(())
}

fn do_expandedvals(
    r: &mut CReader,
    expand: &mut ExpandedVals,
    fep: &mut FepVals,
    version: i32,
) -> Result<()> {
    let n_lambda = fep.n_lambda;
    fep.lambda_start_n = 0;
    fep.lambda_stop_n = n_lambda;
    if version >= TPXV_V79 {
        if n_lambda > 0 {
            expand.init_lambda_weights = r.real_array(n_lambda as usize)?;
            if version < TPXV_V136 {
                let _ = r.bool()?; // former bInit_weights
            }
        }
        expand.nstexpanded = r.int()?;
        expand.elmcmove = r.int()?;
        expand.elamstats = r.int()?;
        expand.lmc_repeats = r.int()?;
        expand.gibbsdeltalam = r.int()?;
        expand.lmc_forced_nstart = r.int()?;
        expand.lmc_seed = r.int()?;
        expand.mc_temp = r.real()?;
        expand.b_symmetrized_t_matrix = r.bool()?;
        expand.nst_tij = r.int()?;
        expand.minvarmin = r.int()?;
        expand.c_range = r.int()?;
        expand.wl_scale = r.real()?;
        expand.wl_ratio = r.real()?;
        expand.init_wl_delta = r.real()?;
        expand.b_wl_oneovert = r.bool()?;
        expand.elmceq = r.int()?;
        expand.equil_steps = r.int()?;
        expand.equil_samples = r.int()?;
        expand.equil_n_at_lam = r.int()?;
        expand.equil_wl_delta = r.real()?;
        expand.equil_ratio = r.real()?;
    }
    if n_lambda > 0 {
        expand.init_lambda_counts = r.real_array(n_lambda as usize)?;
        expand.init_wl_histogram_counts = r.real_array(n_lambda as usize)?;
    }
    Ok(())
}

fn do_pullgrp_tpx_pre95(
    r: &mut CReader,
    pgrp: &mut PullGroup,
    pcrd: &mut PullCoord,
) -> Result<()> {
    let num_atoms = r.int()?;
    pgrp.ind = r.int_array(num_atoms.max(0) as usize)?;
    let num_weights = r.int()?;
    pgrp.weight = r.real_array(num_weights.max(0) as usize)?;
    pgrp.pbcatom = r.int()?;
    pcrd.vec = rvec(r)?;
    pcrd.origin = [0.0; 3];
    let tmp = rvec(r)?;
    pcrd.init = tmp[0];
    pcrd.rate = r.real()?;
    pcrd.k = r.real()?;
    pcrd.kb = r.real()?;
    Ok(())
}

fn do_pull_group(r: &mut CReader, pgrp: &mut PullGroup) -> Result<()> {
    let num_atoms = r.int()?;
    pgrp.ind = r.int_array(num_atoms.max(0) as usize)?;
    let num_weights = r.int()?;
    pgrp.weight = r.real_array(num_weights.max(0) as usize)?;
    pgrp.pbcatom = r.int()?;
    Ok(())
}

fn do_pull_coord(
    r: &mut CReader,
    pcrd: &mut PullCoord,
    version: i32,
    e_pull_old: i32,
    e_geom_old: i32,
    dim_old: [i32; 3],
) -> Result<()> {
    if version >= TPXV_V107 {
        pcrd.etype = r.int()?;
        if version >= TPXV_V110 {
            if pcrd.etype == 5 {
                pcrd.external_potential_provider = r.string()?;
            }
        }
        pcrd.egeom = r.int()?;
        pcrd.ngroup = r.int()?;
        if pcrd.ngroup <= 4 {
            pcrd.group = r.int_array(pcrd.ngroup.max(0) as usize)?;
        } else {
            // More groups than this code supports: the contents are skipped
            // and `ngroup` reset, exactly as GROMACS does.
            let _ = r.int_array(pcrd.ngroup.max(0) as usize)?;
            pcrd.ngroup = 0;
        }
        pcrd.dim = ivec(r)?;
        if version >= TPXV_V124 {
            pcrd.expression = r.string()?;
        }
    } else {
        pcrd.ngroup = 2;
        pcrd.group = vec![r.int()?, r.int()?];
        if version >= TPXV_V101 {
            let egeom = r.int()?;
            pcrd.etype = r.int()?;
            pcrd.egeom = egeom;
            pcrd.ngroup = if pcrd.egeom == 3 { 4 } else { 2 };
            if pcrd.ngroup == 4 {
                pcrd.group.push(r.int()?);
                pcrd.group.push(r.int()?);
            }
            pcrd.dim = ivec(r)?;
        } else {
            pcrd.etype = e_pull_old;
            pcrd.egeom = e_geom_old;
            pcrd.dim = dim_old;
        }
    }
    pcrd.origin = rvec(r)?;
    pcrd.vec = rvec(r)?;
    if version >= TPXV_V101 {
        pcrd.b_start = r.bool()?;
    }
    pcrd.init = r.real()?;
    pcrd.rate = r.real()?;
    pcrd.k = r.real()?;
    pcrd.kb = r.real()?;
    Ok(())
}

fn do_pull(r: &mut CReader, pull: &mut PullParams, version: i32, e_pull_old: i32) -> Result<()> {
    let mut e_geom_old = -1;
    let mut dim_old = [0i32; 3];

    if version >= TPXV_V95 {
        pull.ngroup = r.int()?;
    }
    pull.ncoord = r.int()?;
    if version < TPXV_V95 {
        pull.ngroup = pull.ncoord + 1;
    }
    if version < TPXV_V101 {
        e_geom_old = r.int()?;
        dim_old = ivec(r)?;
        let _ = r.real()?; // inner cylinder radius
    }
    pull.cylinder_r = r.real()?;
    pull.constr_tol = r.real()?;
    if version >= TPXV_V95 {
        pull.b_print_com = r.bool()?;
    }
    if version >= TPXV_V109 {
        pull.b_print_ref_value = r.bool()?;
        pull.b_print_comp = r.bool()?;
    } else if version >= TPXV_V101 {
        let _ = r.int()?; // former bPrintCOM2
        pull.b_print_ref_value = r.bool()?;
        pull.b_print_comp = r.bool()?;
    } else {
        pull.b_print_ref_value = false;
        pull.b_print_comp = true;
    }
    pull.nstxout = r.int()?;
    pull.nstfout = r.int()?;
    if version >= TPXV_V114 {
        pull.b_set_pbc_ref_to_prev_step_com = r.bool()?;
    }

    let ngroup = pull.ngroup.max(0) as usize;
    let ncoord = pull.ncoord.max(0) as usize;
    pull.group = vec![PullGroup::default(); ngroup];
    pull.coord = vec![PullCoord::default(); ncoord];

    if version < TPXV_V95 {
        let mut e_geom = e_geom_old;
        if e_geom == 0 {
            // `pull-geometry = position` (removed).
            return Err(XdrError::Invalid(
                "pull-geometry=position is no longer supported".into(),
            ));
        }
        if e_geom > 0 {
            e_geom = match e_geom {
                1 => 0,
                2 => 1,
                3 => 2,
                4 => 3,
                5 => 4,
                _ => e_geom - 1,
            };
        }
        for g in 0..ngroup {
            // A pull coordinate for group 0 is read and thrown away.
            let cidx = if g > 0 { g - 1 } else { 0 };
            let mut group = PullGroup::default();
            let mut coord = PullCoord::default();
            do_pullgrp_tpx_pre95(r, &mut group, &mut coord)?;
            pull.group[g] = group;
            if cidx < ncoord {
                pull.coord[cidx] = coord;
                if g > 0 {
                    pull.coord[cidx].group = vec![0, g as i32];
                }
            }
        }
        let _ = e_geom;
        pull.b_print_com = !pull.group[0].ind.is_empty();
    } else {
        for g in 0..ngroup {
            do_pull_group(r, &mut pull.group[g])?;
        }
        for c in 0..ncoord {
            let mut coord = std::mem::take(&mut pull.coord[c]);
            do_pull_coord(r, &mut coord, version, e_pull_old, e_geom_old, dim_old)?;
            pull.coord[c] = coord;
        }
    }
    if version >= TPXV_V116 {
        pull.b_xout_average = r.bool()?;
        pull.b_fout_average = r.bool()?;
    }
    Ok(())
}

fn do_rotgrp(r: &mut CReader, rotg: &mut RotGroup) -> Result<()> {
    rotg.e_type = r.int()?;
    rotg.b_mass_w = r.int()?;
    let nat = r.int()?;
    rotg.ind = r.int_array(nat.max(0) as usize)?;
    rotg.x_ref_original = rvec_array(r, nat.max(0) as usize)?;
    rotg.input_vec = rvec(r)?;
    rotg.pivot = rvec(r)?;
    rotg.rate = r.real()?;
    rotg.k = r.real()?;
    rotg.slab_dist = r.real()?;
    rotg.min_gaussian = r.real()?;
    rotg.eps = r.real()?;
    rotg.e_fittype = r.int()?;
    rotg.pot_angle_nstep = r.int()?;
    rotg.pot_angle_step = r.real()?;
    Ok(())
}

fn do_rot(r: &mut CReader, rot: &mut Rot) -> Result<()> {
    let num_groups = r.int()?;
    rot.nstrout = r.int()?;
    rot.nstsout = r.int()?;
    rot.grp = vec![RotGroup::default(); num_groups.max(0) as usize];
    for g in 0..rot.grp.len() {
        let mut grp = std::mem::take(&mut rot.grp[g]);
        do_rotgrp(r, &mut grp)?;
        rot.grp[g] = grp;
    }
    Ok(())
}

fn do_imd(r: &mut CReader, imd: &mut Imd) -> Result<()> {
    let nat = r.int()?;
    imd.ind = r.int_array(nat.max(0) as usize)?;
    Ok(())
}

fn do_swapgroup(r: &mut CReader, g: &mut SwapGroup) -> Result<()> {
    g.molname = r.string()?;
    let num_atoms = r.int()?;
    g.ind = r.int_array(num_atoms.max(0) as usize)?;
    g.nmol_req = [r.int()?, r.int()?];
    Ok(())
}

fn do_swapcoords_tpx(r: &mut CReader, swap: &mut SwapCoords, version: i32) -> Result<()> {
    if version >= TPXV_V105 {
        let num_groups = r.int()?;
        swap.groups = vec![SwapGroup::default(); num_groups.max(0) as usize];
        for g in 0..swap.groups.len() {
            let mut group = std::mem::take(&mut swap.groups[g]);
            do_swapgroup(r, &mut group)?;
            swap.groups[g] = group;
        }
        swap.massw_split[0] = r.bool()?;
        swap.massw_split[1] = r.bool()?;
        swap.nstswap = r.int()?;
        swap.n_average = r.int()?;
        swap.threshold = r.real()?;
        swap.cyl0r = r.real()?;
        swap.cyl0u = r.real()?;
        swap.cyl0l = r.real()?;
        swap.cyl1r = r.real()?;
        swap.cyl1u = r.real()?;
        swap.cyl1l = r.real()?;
    } else {
        // Old CompEl files always store split0, split1, solvent, anions and
        // cations in a fixed order.
        swap.groups = vec![SwapGroup::default(); 5];
        swap.groups[0].molname = "split0".into();
        swap.groups[1].molname = "split1".into();
        swap.groups[2].molname = "solvent".into();
        swap.groups[3].molname = "anions".into();
        swap.groups[4].molname = "cations".into();

        let num_atoms = r.int()?;
        swap.groups[3].ind = vec![0; num_atoms.max(0) as usize];
        let num_atoms = r.int()?;
        swap.groups[2].ind = vec![0; num_atoms.max(0) as usize];
        let num_atoms = r.int()?;
        swap.groups[0].ind = vec![0; num_atoms.max(0) as usize];
        swap.massw_split[0] = r.bool()?;
        let num_atoms = r.int()?;
        swap.groups[1].ind = vec![0; num_atoms.max(0) as usize];
        swap.massw_split[1] = r.bool()?;
        swap.nstswap = r.int()?;
        swap.n_average = r.int()?;
        swap.threshold = r.real()?;
        swap.cyl0r = r.real()?;
        swap.cyl0u = r.real()?;
        swap.cyl0l = r.real()?;
        swap.cyl1r = r.real()?;
        swap.cyl1u = r.real()?;
        swap.cyl1l = r.real()?;

        let n = swap.groups[3].ind.len();
        swap.groups[3].ind = r.int_array(n)?;
        let n = swap.groups[2].ind.len();
        swap.groups[2].ind = r.int_array(n)?;
        let n = swap.groups[0].ind.len();
        swap.groups[0].ind = r.int_array(n)?;
        let n = swap.groups[1].ind.len();
        swap.groups[1].ind = r.int_array(n)?;

        for j in 0..2 {
            swap.groups[3].nmol_req[j] = r.int()?;
            swap.groups[4].nmol_req[j] = r.int()?;
        }
    }
    if version >= TPXV_V104 {
        swap.bulk_offset[0] = r.real()?;
        swap.bulk_offset[1] = r.real()?;
    }
    Ok(())
}

fn do_legacy_efield(r: &mut CReader, out: &mut Vec<[f64; 4]>) -> Result<()> {
    for _ in 0..3 {
        let n = r.int()?;
        let nt = r.int()?;
        let aa = r.real_array(n.max(0) as usize)?;
        let phi = r.real_array(n.max(0) as usize)?;
        let at = r.real_array(nt.max(0) as usize)?;
        let phit = r.real_array(nt.max(0) as usize)?;
        if n > 0 {
            if n > 1 || nt > 1 {
                return Err(XdrError::Invalid(
                    "Can not handle tpr files with more than one electric field term per direction."
                        .into(),
                ));
            }
            out.push([
                aa.first().copied().unwrap_or(0.0),
                at.first().copied().unwrap_or(0.0),
                phi.first().copied().unwrap_or(0.0),
                phit.first().copied().unwrap_or(0.0),
            ]);
        }
    }
    Ok(())
}

fn do_awh_dim(r: &mut CReader) -> Result<AwhDimParams> {
    Ok(AwhDimParams {
        coord_provider: r.int()?,
        coord_index: r.int()?,
        origin: r.double()?,
        end: r.double()?,
        period: r.double()?,
        force_constant: r.double()?,
        diffusion: r.double()?,
        coord_value_init: r.double()?,
        cover_diameter: r.double()?,
    })
}

fn do_awh_bias(
    r: &mut CReader,
    without_growth_factor: bool,
    without_target_metric_scaling: bool,
    without_histogram_tolerance: bool,
) -> Result<AwhBiasParams> {
    let mut bias = AwhBiasParams::default();
    bias.e_target = r.int()?;
    bias.target_beta_scaling = r.double()?;
    bias.target_cutoff = r.double()?;
    bias.e_growth = r.int()?;
    if without_growth_factor {
        bias.growth_factor = 3.0;
    } else {
        bias.growth_factor = r.double()?;
    }
    bias.b_user_data = r.int()? != 0;
    if without_target_metric_scaling {
        bias.target_metric_scaling_limit = 10.0;
    } else {
        bias.scale_target_by_metric = r.bool()?;
        bias.target_metric_scaling_limit = r.double()?;
    }
    bias.error_init = r.double()?;
    let num_dimensions = r.int()?;
    bias.share_group = r.int()?;
    bias.equilibrate_histogram = r.bool()?;
    if without_histogram_tolerance {
        bias.histogram_tolerance = 0.2;
    } else {
        bias.histogram_tolerance = r.double()?;
    }
    for _ in 0..num_dimensions.max(0) {
        let dim = do_awh_dim(r)?;
        bias.dim_params.push(dim);
    }
    Ok(bias)
}

fn do_awh(r: &mut CReader, awh: &mut AwhParams, version: i32) -> Result<()> {
    let number_of_biases = r.int()?;
    awh.nstout = r.int()?;
    awh.seed = r.int64()?;
    awh.nst_sample_coord = r.int()?;
    awh.num_samples_update_free_energy = r.int()?;
    awh.potential = r.int()?;
    awh.share_bias_multisim = r.bool()?;
    for _ in 0..number_of_biases.max(0) {
        let bias = do_awh_bias(
            r,
            version < TPXV_V130,
            version < TPXV_V132,
            version < TPXV_V138,
        )?;
        awh.bias.push(bias);
    }
    Ok(())
}

/// Decodes the inputrec prefix and returns it together with the byte offset of
/// `nsteps`.
pub fn parse_inputrec(r: &mut CReader, header: &TpxHeader) -> Result<(InputRec, usize)> {
    let mut ir = InputRec::default();
    let version = header.file_version;
    if version >= TPXV_V53 {
        ir.pbc_type = r.int()?;
        ir.b_periodic_mols = r.bool()?;
    }

    ir.integrator = r.int()?;
    let nsteps_offset = r.position();
    if version >= TPXV_V62 {
        ir.nsteps = r.int64()?;
        ir.init_step = r.int64()?;
    } else {
        ir.nsteps = r.int()? as i64;
        ir.init_step = r.int()? as i64;
    }
    ir.simulation_part = r.int()?;
    if version >= TPXV_V122 {
        ir.use_mts = r.bool()?;
        let mut num_levels = 0i32;
        if ir.use_mts {
            num_levels = r.int()?;
        }
        for _ in 0..num_levels.max(0) {
            ir.mts_levels.push(MtsLevel {
                force_groups: r.int()?,
                step_factor: r.int()?,
            });
        }
    }
    ir.mass_repartition_factor = if version >= TPXV_V131 { r.real()? } else { 1.0 };
    if version >= TPXV_V129 {
        ir.ensemble_temperature_setting = r.int()?;
        ir.ensemble_temperature = r.real()?;
    }
    ir.nstcalcenergy = if version >= TPXV_V67 { r.int()? } else { 1 };
    if version >= TPXV_V81 {
        ir.cutoff_scheme = r.int()?;
        if version < TPXV_V94 {
            ir.cutoff_scheme = match ir.cutoff_scheme {
                0 => 1,
                1 => 0,
                other => other,
            };
        }
    } else {
        ir.cutoff_scheme = 1; // Group
    }
    let _ = r.int()?; // used to be ns_type
    ir.nstlist = r.int()?;
    let _ = r.int()?; // used to be ndelta
    ir.rtpi = r.real()?;
    ir.nstcomm = r.int()?;
    ir.comm_mode = r.int()?;
    if version < TPXV_V100 {
        let _ = r.int()?; // nstcheckpoint
    }
    ir.nstcgsteep = r.int()?;
    ir.nbfgscorr = r.int()?;
    ir.nstlog = r.int()?;
    ir.nstxout = r.int()?;
    ir.nstvout = r.int()?;
    ir.nstfout = r.int()?;
    ir.nstenergy = r.int()?;
    ir.nstxout_compressed = r.int()?;
    if version >= TPXV_V59 {
        ir.init_t = r.double()?;
        ir.delta_t = r.double()?;
    } else {
        ir.init_t = r.real()?;
        ir.delta_t = r.real()?;
    }
    ir.x_compression_precision = r.real()?;
    ir.verletbuf_tol = if version >= TPXV_V81 { r.real()? } else { 0.0 };
    ir.verlet_buffer_pressure_tolerance = if version >= TPXV_V133 { r.real()? } else { -1.0 };
    ir.rlist = r.real()?;
    if version >= TPXV_V67 && version < TPXV_V108 {
        let _ = r.real()?; // rlistlong
    }
    if version >= TPXV_V82 && version != TPXV_V90 {
        let _ = r.int()?; // nstcalclr
    }
    ir.coulombtype = r.int()?;
    ir.coulomb_modifier = if version >= TPXV_V81 {
        r.int()?
    } else if ir.cutoff_scheme == 0 {
        1 // PotShift for the Verlet scheme
    } else {
        2 // None
    };
    ir.rcoulomb_switch = r.real()?;
    ir.rcoulomb = r.real()?;
    ir.vdwtype = r.int()?;
    ir.vdw_modifier = if version >= TPXV_V81 {
        r.int()?
    } else if ir.cutoff_scheme == 0 {
        1
    } else {
        2
    };
    ir.rvdw_switch = r.real()?;
    ir.rvdw = r.real()?;
    ir.disp_corr = r.int()?;
    ir.epsilon_r = r.real()?;
    ir.epsilon_rf = r.real()?;
    ir.tabext = r.real()?;
    if version < TPXV_V113 {
        // Implicit solvent parameters of old files (read and ignored).
        let _ = r.int()?;
        let _ = r.int()?;
        let _ = r.real()?;
        let _ = r.real()?;
        let _ = r.int()?;
        let _ = r.real()?;
        let _ = r.real()?;
        let _ = r.real()?;
        let _ = r.real()?;
        if version >= TPXV_V60 {
            let _ = r.real()?;
            let _ = r.int()?;
        }
        let _ = r.real()?;
    }
    ir.fourier_spacing = if version >= TPXV_V81 { r.real()? } else { 0.0 };
    ir.nkx = r.int()?;
    ir.nky = r.int()?;
    ir.nkz = r.int()?;
    ir.pme_order = r.int()?;
    ir.ewald_rtol = r.real()?;
    ir.ewald_rtol_lj = if version >= TPXV_V93 {
        r.real()?
    } else {
        ir.ewald_rtol
    };
    ir.ewald_geometry = r.int()?;
    ir.epsilon_surface = r.real()?;
    if version < TPXV_V100 {
        let _ = r.bool()?; // bOptFFT
    }
    // `LongRangeVdW::Geometric` is the default for older files.
    ir.ljpme_combination_rule = if version >= TPXV_V93 { r.int()? } else { 0 };
    ir.b_continuation = r.bool()?;
    ir.tcoupl = r.int()?;
    ir.b_print_nh_chains = if version >= TPXV_V79 { r.bool()? } else { false };
    ir.nsttcouple = if version >= TPXV_V71 {
        r.int()?
    } else {
        ir.nstcalcenergy
    };
    ir.pcoupl = r.int()?;
    ir.pcoupltype = r.int()?;
    ir.nstpcouple = if version >= TPXV_V71 {
        r.int()?
    } else {
        ir.nstcalcenergy
    };
    ir.tau_p = r.real()?;
    ir.ref_p = [rvec(r)?, rvec(r)?, rvec(r)?];
    ir.compress = [rvec(r)?, rvec(r)?, rvec(r)?];
    ir.refcoord_scaling = r.int()?;

    let num_posres_com = if version >= TPXV_V135 {
        r.int()?
    } else {
        1
    };
    ir.posres_com_allocated = num_posres_com > 0;
    ir.posres_com = rvec_array(r, num_posres_com.max(0) as usize)?;
    ir.posres_com_b = rvec_array(r, num_posres_com.max(0) as usize)?;

    if version < TPXV_V79 {
        let _ = r.int()?; // andersen_seed
    }
    ir.shake_tol = r.real()?;
    ir.efep = r.int()?;
    do_fepvals(r, &mut ir.fepvals, version)?;
    tr(r, "fepvals done");
    if version >= TPXV_V79 {
        ir.b_sim_temp = r.bool()?;
    }
    if ir.b_sim_temp {
        let n_lambda = ir.fepvals.n_lambda;
        do_simtempvals(r, &mut ir.simtempvals, n_lambda, version)?;
    }
    if version >= TPXV_V79 {
        ir.b_expanded = r.bool()?;
    }
    if ir.b_expanded {
        let mut fep = std::mem::take(&mut ir.fepvals);
        do_expandedvals(r, &mut ir.expandedvals, &mut fep, version)?;
        ir.fepvals = fep;
    }
    ir.e_disre = r.int()?;
    ir.e_disre_weighting = r.int()?;
    ir.b_disre_mixed = r.bool()?;
    ir.dr_fc = r.real()?;
    ir.dr_tau = r.real()?;
    ir.nstdisreout = r.int()?;
    ir.orires_fc = r.real()?;
    ir.orires_tau = r.real()?;
    ir.nstorireout = r.int()?;
    if version < TPXV_V79 {
        let _ = r.real()?; // dihre_fc
    }
    ir.em_stepsize = r.real()?;
    ir.em_tol = r.real()?;
    ir.b_shake_sor = r.bool()?;
    ir.niter = r.int()?;
    ir.fc_stepsize = r.real()?;
    ir.e_constr_alg = r.int()?;
    ir.n_proj_order = r.int()?;
    ir.lincs_warn_angle = r.real()?;
    ir.n_lincs_iter = r.int()?;
    ir.bd_fric = r.real()?;
    ir.ld_seed = if version >= TPXV_V97 {
        r.int64()?
    } else {
        r.int()? as i64
    };
    ir.deform = [rvec(r)?, rvec(r)?, rvec(r)?];
    ir.cos_accel = r.real()?;
    ir.userint1 = r.int()?;
    ir.userint2 = r.int()?;
    ir.userint3 = r.int()?;
    ir.userint4 = r.int()?;
    ir.userreal1 = r.real()?;
    ir.userreal2 = r.real()?;
    ir.userreal3 = r.real()?;
    ir.userreal4 = r.real()?;

    if version >= TPXV_V77 && version < TPXV_V106 {
        let b_adress = r.bool()?;
        if b_adress {
            let _ = r.int()?;
            let _ = r.real()?;
            let _ = r.real()?;
            let _ = r.real()?;
            let _ = r.int()?;
            let _ = r.int()?;
            let _ = rvec(r)?;
            let num_thermo_force_groups = r.int()?;
            let _ = r.real()?;
            let num_energy_groups = r.int()?;
            let _ = r.int()?;
            if num_thermo_force_groups > 0 {
                let _ = r.int_array(num_thermo_force_groups as usize)?;
            }
            if num_energy_groups > 0 {
                let _ = r.int_array(num_energy_groups as usize)?;
            }
        }
    }

    // Pulling: before `tpxv_PullCoordTypeGeom` an enum stored the old
    // algorithm, which has to be shifted to the current enum values.
    let mut e_pull_old = 0;
    if version >= TPXV_V101 {
        ir.b_pull = r.bool()?;
    } else {
        let e = r.int()?;
        ir.b_pull = e != 0;
        e_pull_old = match e {
            0 => 0,
            1 => 0,
            2 => 1,
            3 => 2,
            4 => 3,
            5 => 4,
            6 => 5,
            _ => e,
        };
    }
    if ir.b_pull {
        let mut pull = std::mem::take(&mut ir.pull);
        do_pull(r, &mut pull, version, e_pull_old)?;
        ir.pull = pull;
    }

    if version >= TPXV_V112 {
        ir.b_do_awh = r.bool()?;
        if ir.b_do_awh {
            let mut awh = std::mem::take(&mut ir.awh);
            do_awh(r, &mut awh, version)?;
            ir.awh = awh;
        }
    }
    if version >= TPXV_V74 {
        ir.b_rot = r.bool()?;
        if ir.b_rot {
            let mut rot = std::mem::take(&mut ir.rot);
            do_rot(r, &mut rot)?;
            ir.rot = rot;
        }
    }
    if version >= TPXV_V99 {
        ir.b_imd = r.bool()?;
            if ir.b_imd {
            let mut imd = std::mem::take(&mut ir.imd);
            do_imd(r, &mut imd)?;
            ir.imd = imd;
        }
    }

    // --- grpopts -----------------------------------------------------------
    ir.opts.ngtc = r.int()?;
    ir.opts.nhchainlength = if version >= TPXV_V69 { r.int()? } else { 1 };
    let num_acceleration_groups = r.int()?;
    ir.opts.ngfrz = r.int()?;
    ir.opts.ngener = r.int()?;
    if ir.opts.ngtc > 0 {
        let n = ir.opts.ngtc as usize;
        ir.opts.nrdf = r.real_array(n)?;
        ir.opts.ref_t = r.real_array(n)?;
        ir.opts.tau_t = r.real_array(n)?;
    }
    if ir.opts.ngfrz > 0 {
        let mut nfreeze = Vec::with_capacity(ir.opts.ngfrz as usize);
        for _ in 0..ir.opts.ngfrz {
            nfreeze.push(ivec(r)?);
        }
        ir.opts.nfreeze = nfreeze;
    }
    if num_acceleration_groups > 0 {
        ir.opts.acceleration = rvec_array(r, num_acceleration_groups as usize)?;
    }
    let egp = (ir.opts.ngener as i64 * ir.opts.ngener as i64).max(0) as usize;
    ir.opts.egp_flags = r.int_array(egp)?;
    if ir.opts.ngtc > 0 {
        let n = ir.opts.ngtc as usize;
        ir.opts.annealing = r.int_array(n)?;
        ir.opts.anneal_npoints = r.int_array(n)?;
        for j in 0..n {
            let k = ir.opts.anneal_npoints[j].max(0) as usize;
            ir.opts.anneal_time.push(r.real_array(k)?);
            ir.opts.anneal_temp.push(r.real_array(k)?);
        }
    }

    // --- walls -------------------------------------------------------------
    ir.nwall = r.int()?;
    ir.wall_type = r.int()?;
    ir.wall_r_linpot = r.real()?;
    ir.wall_atomtype = [r.int()?, r.int()?];
    ir.wall_density = [r.real()?, r.real()?];
    ir.wall_ewald_zfac = r.real()?;

    if version < TPXV_V111 {
        do_legacy_efield(r, &mut ir.legacy_efield)?;
    }
    tr(r, "after walls and electric field");
    if version >= TPXV_V96 {
        ir.e_swap_coords = r.int()?;
        if ir.e_swap_coords != 0 {
            let mut swap = std::mem::take(&mut ir.swap);
            do_swapcoords_tpx(r, &mut swap, version)?;
            ir.swap = swap;
        }
    }

    // --- QMMM --------------------------------------------------------------
    ir.opts.b_qmmm = r.bool()?;
    let _unused_qmmm_scheme = r.int()?;
    let _unused_scalefactor = r.real()?;
    ir.opts.ngqm = r.int()?;
    if ir.opts.ngqm > 0 && ir.opts.b_qmmm {
        let n = ir.opts.ngqm as usize;
        let _ = r.int_array(4 * n)?;
        for _ in 0..n {
            let _ = r.bool()?;
        }
        let _ = r.int_array(2 * n)?;
        let _ = r.real_array(2 * n)?;
        let _ = r.int_array(3 * n)?;
    }

    tr(r, "after qmmm");
    if version >= TPXV_V111 {
        ir.params = read_kvt_object(r)?;
    } else {
        // `do_legacy_efield` builds this part of the tree.
        let mut electric_field = Vec::new();
        for (dim, name) in ["x", "y", "z"].iter().enumerate() {
            // A direction is only added when the file actually stores an
            // electric field term for it.
            if let Some(v) = ir.legacy_efield.get(dim) {
                // `addValue<real>()`: `real` is `float` in this build.
                let fields = vec![
                    ("E0".to_string(), KvtValue::Float(v[0] as f32)),
                    ("omega".to_string(), KvtValue::Float(v[1] as f32)),
                    ("t0".to_string(), KvtValue::Float(v[2] as f32)),
                    ("sigma".to_string(), KvtValue::Float(v[3] as f32)),
                ];
                electric_field.push((name.to_string(), KvtValue::Object(fields)));
            }
        }
        ir.params = KvtValue::Object(vec![(
            "applied-forces".to_string(),
            KvtValue::Object(vec![(
                "electric-field".to_string(),
                KvtValue::Object(electric_field),
            )]),
        )]);
    }
    if version >= TPXV_V117 {
        ir.internal_parameters = read_kvt_object(r)?;
    }

    if version < TPXV_V129 {
        let has_annealing = ir.opts.annealing.iter().any(|&a| a != 0);
        ir.ensemble_temperature_setting = if has_annealing || ir.b_sim_temp {
            2 // variable
        } else if integrator_has_reference_temperature(&ir) {
            1 // constant
        } else {
            0 // not available
        };
        if ir.ensemble_temperature_setting == 1 {
            ir.ensemble_temperature = ir.opts.ref_t.first().copied().unwrap_or(0.0);
        } else {
            ir.ensemble_temperature = -1.0;
        }
    }
    Ok((ir, nsteps_offset))
}

/// `integratorHasReferenceTemperature()` from `mdtypes/inputrec.cpp`.
fn integrator_has_reference_temperature(ir: &InputRec) -> bool {
    // EI_SD(9), BD(3), TPI(7) and TPIC(8), or any temperature coupling.
    ir.tcoupl != 0 || matches!(ir.integrator, 3 | 7 | 8 | 9)
}
