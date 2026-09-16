//! Run input (.tpr) file handling, mirroring `gromacs/fileio/tpxio.cpp`.
//!
//! Every layout GROMACS itself accepts (version 58 and later) is decoded.  The
//! file consists of an XDR header followed by a body which, in order, contains
//!
//! ```text
//! state   : boxm, box_rel, boxv, temperature coupling state
//! mtop    : symbol table, name, force field parameters, molecule types,
//!           molecule blocks, #atoms, groups, exclusions
//! state   : x, v
//! ir      : pbc type, periodic molecules, inputrec
//! ```
//!
//! Files before `tpxv_AddSizeField` (119) store the body as plain XDR, newer
//! ones as a single opaque blob written with the compact in-memory serializer;
//! `CReader`/`CWriter` switch between the two.
//!
//! The inputrec itself is decoded by [`crate::ir::parse_inputrec`], and the
//! topology printing lives in [`crate::tpdump`].

use crate::frame::{Atom, Atoms, Matrix, ResInfo, Rvec};
use crate::xdr::{Reader, Result, Writer, XdrError};

/// When `GMXRS_TRACE` is set, prints the reader position at each parsing stage.
fn trace(msg: &str, r: &CReader) {
    if std::env::var_os("GMXRS_TRACE").is_some() {
        eprintln!("[tpr] {msg}: pos={}", r.position());
    }
}

/// The tpr body is written with `gmx::InMemorySerializer`, which is *not* XDR:
/// values are stored in their native width (with the byte order swapped to big
/// endian on little endian hosts) and strings are a `uint64` length followed by
/// the raw characters without padding or terminator.
pub struct CReader<'a> {
    data: &'a [u8],
    pos: usize,
    double_precision: bool,
    /// True when the body was written through the XDR serializer (tpx files
    /// older than `tpxv_AddSizeField`), which stores strings as XDR strings
    /// and pads `bool`/`uchar`/`ushort` to four bytes.
    xdr: bool,
}

impl<'a> CReader<'a> {
    pub fn new(data: &'a [u8], double_precision: bool, xdr: bool) -> Self {
        CReader {
            data,
            pos: 0,
            double_precision,
            xdr,
        }
    }

    pub fn position(&self) -> usize {
        self.pos
    }

    pub fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    fn bytes(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return Err(XdrError::Truncated("tpr body"));
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    pub fn int(&mut self) -> Result<i32> {
        let b = self.bytes(4)?;
        Ok(i32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    pub fn int64(&mut self) -> Result<i64> {
        let b = self.bytes(8)?;
        Ok(i64::from_be_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    pub fn bool(&mut self) -> Result<bool> {
        if self.xdr {
            Ok(self.int()? != 0)
        } else {
            Ok(self.bytes(1)?[0] != 0)
        }
    }

    pub fn uchar(&mut self) -> Result<u8> {
        if self.xdr {
            Ok(self.int()? as u8)
        } else {
            Ok(self.bytes(1)?[0])
        }
    }

    pub fn ushort(&mut self) -> Result<u16> {
        if self.xdr {
            Ok(self.int()? as u16)
        } else {
            let b = self.bytes(2)?;
            Ok(u16::from_be_bytes([b[0], b[1]]))
        }
    }

    pub fn float(&mut self) -> Result<f32> {
        let b = self.bytes(4)?;
        Ok(f32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    pub fn double(&mut self) -> Result<f64> {
        let b = self.bytes(8)?;
        Ok(f64::from_be_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    pub fn real(&mut self) -> Result<f64> {
        if self.double_precision {
            // `real` is `float` in this build, so a double precision file is
            // truncated on read, exactly as GROMACS does.
            Ok(self.double()? as f32 as f64)
        } else {
            Ok(self.float()? as f64)
        }
    }

    pub fn string(&mut self) -> Result<String> {
        if self.xdr {
            // `XdrSerializer::doString()`: length including the null byte,
            // then the plain XDR string.
            let _len_plus_one = self.int()?;
            let len = self.int()?;
            if len < 0 {
                return Err(XdrError::Invalid("negative string length in tpr body".into()));
            }
            let len = len as usize;
            let raw = self.bytes(len)?;
            let pad = (4 - (len % 4)) % 4;
            let _ = self.bytes(pad)?;
            let end = raw.iter().position(|&c| c == 0).unwrap_or(raw.len());
            return Ok(String::from_utf8_lossy(&raw[..end]).into_owned());
        }
        let len = self.int64()?;
        if len < 0 {
            return Err(XdrError::Invalid("negative string length in tpr body".into()));
        }
        let raw = self.bytes(len as usize)?;
        Ok(String::from_utf8_lossy(raw).into_owned())
    }

    pub fn int_array(&mut self, n: usize) -> Result<Vec<i32>> {
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.int()?);
        }
        Ok(v)
    }

    pub fn real_array(&mut self, n: usize) -> Result<Vec<f64>> {
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            v.push(self.real()?);
        }
        Ok(v)
    }
}

/// First tpx version that stores the body as a single opaque blob.
pub const TPXV_ADD_SIZE_FIELD: i32 = 119;
/// `TpxGeneration::AddSizeField` (values start at `Initial = 26`), the first
/// generation that stores `sizeOfTprBody` in the header.
pub const TPX_GENERATION_ADD_SIZE_FIELD: i32 = 27;
/// Oldest tpx version this code accepts, matching `tpx_incompatible_version`.
pub const TPX_MIN_SUPPORTED_VERSION: i32 = 58;
/// Value of `tpx_version` (the version this code writes).
pub const TPX_VERSION: i32 = 138;
/// Value of `tpx_generation`.
pub const TPX_GENERATION: i32 = 29;

// Versions that changed a serialized structure; see the `tpxv` enum in
// `gmx-2026.3/src/gromacs/fileio/tpxio.cpp`.
const TPXV_V51: i32 = 51; // box_rel in the state
const TPXV_V56: i32 = 56; // extra (removed) matrix in the state
const TPXV_V60: i32 = 60; // implicit solvent atom type arrays
const TPXV_V62: i32 = 62; // nsteps/init_step became int64
const TPXV_V63: i32 = 63; // resinfo gained nr/ic
const TPXV_V65: i32 = 65; // cmap
const TPXV_V66: i32 = 66; // ffparams reppow
const TPXV_V69: i32 = 69; // single temperature coupling array in the state
const TPXV_V79: i32 = 79; // fep_state, Urey-Bradley/Morse B state
const TPXV_V82: i32 = 82; // dihedral restraints B state, nstcalclr
const TPXV_V103: i32 = 103; // intermolecular interactions
const TPXV_V113: i32 = 113; // implicit solvation removed
const TPXV_V120: i32 = 120; // intermolecular exclusion group
const TPXV_V127: i32 = 127; // Thole rfac removed
const TPXV_V128: i32 = 128; // atomtypes removed from the topology
const TPXV_V134: i32 = 134; // Martini bonded B-state parameters
const TPXV_V135: i32 = 135; // multiple COM groups for refcoord-scaling
/// Number of `InteractionFunction` values (including the energies).
pub const INTERACTION_FUNCTION_COUNT: usize = 95;
/// `ftupd[]` from `tpxio.cpp`: interaction function types that were inserted
/// into the enum later than the listed file version.
///
/// Reading an older file requires shifting all function type numbers that are
/// larger or equal to the inserted type, and the corresponding interaction
/// lists are not present in the file at all.
const FTUPD: &[(i32, i32)] = &[
    (70, 9),   // RestraintBonds
    (98, 12),  // RestrictedBendingPotential
    (76, 13),  // LinearAngles
    (98, 21),  // RestrictedTorsionPotential
    (98, 22),  // CombinedBendingTorsionPotential
    (65, 27),  // DihedralEnergyCorrectionMap
    (60, 28),  // GeneralizedBorn12PolarizationUnused
    (61, 29),  // GeneralizedBorn13PolarizationUnused
    (61, 30),  // GeneralizedBorn14PolarizationUnused
    (72, 31),  // GeneralizedBornPolarizationUnused
    (72, 32),  // NonpolarSolvationUnused
    (93, 46),  // LennardJonesReciprocalSpace
    (76, 51),  // AnharmonicPolarization
    (90, 53),  // FlatBottomedPositionRestraints
    (121, 65), // VirtualSite1
    (118, 67), // VirtualSite2FlexibleDistance
    (117, 76), // DensityFitting
    (137, 78), // NeuralNetworkPotentialEnergy
    (69, 84),  // VirialTemperatureUnused
    (66, 85),  // PressureDispersionCorrection
    (79, 90),  // dVCoulombdLambda
    (79, 91),  // dVvanderWaalsdLambda
    (79, 92),  // dVbondeddLambda
    (79, 93),  // dVrestraintdLambda
    (79, 94),  // dVtemperaturedLambda
];
/// `NR_RBDIHS` / `NR_CBTDIHS`.
const NR_DIHEDRAL_PARAMS: usize = 6;
/// Number of `SimulationAtomGroupType` values.
const SIMULATION_ATOM_GROUP_COUNT: usize = 10;

#[derive(Debug, Clone, Default)]
pub struct TpxHeader {
    pub version_string: String,
    pub precision: i32,
    pub is_double: bool,
    pub file_version: i32,
    pub file_generation: i32,
    pub tag: String,
    pub natoms: i32,
    pub ngtc: i32,
    pub fep_state: i32,
    pub lambda: f64,
    pub b_ir: bool,
    pub b_top: bool,
    pub b_x: bool,
    pub b_v: bool,
    pub b_f: bool,
    pub b_box: bool,
    pub size_of_tpr_body: i64,
}

/// One serialized force field parameter value (`t_iparams` is a union of
/// reals and ints, dumped in the order they appear in the file).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PVal {
    Real(f64),
    Int(i64),
}

/// `gmx_ffparams_t`, keeping the parameter values so that `dump` can print
/// them the way `pr_iparams()` does.
#[derive(Debug, Clone, Default)]
pub struct FfParams {
    pub atnr: i32,
    pub functype: Vec<i32>,
    pub reppow: f64,
    pub fudge_qq: f64,
    /// Values per `functype` entry, in serialization order.
    pub iparams: Vec<Vec<PVal>>,
}

impl FfParams {
    /// `(c6, c12)` of parameter entry `i` when it is a short range
    /// Lennard-Jones pair (`LJ_SR`).
    ///
    /// The first `atnr * atnr` entries of the parameter list are these pairs,
    /// in the same order as `nbfp[i * atnr + j]`.
    pub fn lj_sr(&self, i: usize) -> Option<(f64, f64)> {
        if self.functype.get(i).copied() != Some(37) {
            return None;
        }
        let p = self.iparams.get(i)?;
        match (p.first()?, p.get(1)?) {
            (PVal::Real(c6), PVal::Real(c12)) => Some((*c6, *c12)),
            _ => None,
        }
    }

    /// `atnr` as `usize`.
    pub fn atnr_usize(&self) -> usize {
        self.atnr.max(0) as usize
    }
}

/// A molecule type (`gmx_moltype_t`) reduced to what the tools need.
#[derive(Debug, Clone, Default)]
pub struct MolType {
    pub name: String,
    pub atoms: Atoms,
    pub excls: Vec<Vec<i32>>,
    /// Edges of the bonded interaction graph, i.e. the pairs of atoms that
    /// `mk_graph_moltype()` connects when making molecules whole.  Built from
    /// the interaction types that carry the `IF_CHEMBOND` flag.
    pub bonds: Vec<(u32, u32)>,
    /// All interaction lists of the molecule type, indexed by
    /// `InteractionFunction`.  Needed because virtual sites are not chemical
    /// bonds but still join different parts of the molecular graph.
    pub ilists: Vec<Vec<i32>>,
}

/// A molecule block (`gmx_molblock_t`).
#[derive(Debug, Clone, Default)]
pub struct MolBlock {
    pub moltype_index: i32,
    pub nmol: i32,
    /// Position restraint coordinates of the first molecule (`posres_xA/B`).
    pub posres_xa: Vec<[f32; 3]>,
    pub posres_xb: Vec<[f32; 3]>,
}

/// Global topology (`gmx_mtop_t`), partially decoded.
#[derive(Debug, Clone, Default)]
pub struct Mtop {
    pub symtab: Vec<String>,
    pub name: String,
    pub natoms: usize,
    pub moltypes: Vec<MolType>,
    pub molblocks: Vec<MolBlock>,
    pub ffparams: FfParams,
    pub b_intermolecular: bool,
    /// System wide interaction lists, printed when
    /// `bIntermolecularInteractions` is set.
    pub intermolecular_ilists: Vec<Vec<i32>>,
    /// Group indices for each `SimulationAtomGroupType`.
    pub groups: Vec<Vec<i32>>,
    pub group_names: Vec<String>,
    /// Per atom group numbers for each `SimulationAtomGroupType`.
    pub group_numbers: Vec<Vec<u8>>,
    /// `cmap_grid->grid_spacing` (0 when the file has no cmap).
    pub cmap_grid_spacing: i32,
    /// One `grid_spacing^2 * 4` block of dihedral energy correction values
    /// per grid.
    pub cmap_data: Vec<Vec<f32>>,
}

impl Mtop {
    /// Atom ranges of each individual molecule, in global atom order.
    pub fn molecule_ranges(&self) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut offset = 0usize;
        for mb in &self.molblocks {
            let n = self
                .moltypes
                .get(mb.moltype_index as usize)
                .map(|mt| mt.atoms.nr())
                .unwrap_or(0);
            for _ in 0..mb.nmol.max(0) {
                out.push((offset, offset + n));
                offset += n;
            }
        }
        out
    }

    /// Builds the global atom list by replaying the molecule blocks.
    ///
    /// Residues of molecules that contain at most
    /// `maxResiduesPerMoleculeToTriggerRenumber` residues are renumbered
    /// sequentially, exactly like `gmx_mtop_global_atoms()` does.
    pub fn global_atoms(&self) -> Atoms {
        let max_res_renum = if self.molblocks.len() == 1 && self.molblocks[0].nmol == 1 {
            0
        } else {
            1
        };
        let mut max_resnr = 0;
        for mt in &self.moltypes {
            if mt.atoms.nres() as i32 > max_res_renum {
                for ri in &mt.atoms.resinfo {
                    if ri.nr > max_resnr {
                        max_resnr = ri.nr;
                    }
                }
            }
        }

        let mut atoms = Atoms::default();
        atoms.name = self.name.clone();
        for mb in &self.molblocks {
            let mt = match self.moltypes.get(mb.moltype_index as usize) {
                Some(m) => m,
                None => continue,
            };
            let nres = mt.atoms.nres();
            let nmol = mb.nmol.max(0);
            let dest_res0 = atoms.resinfo.len();
            for _ in 0..nmol {
                let res_offset = atoms.resinfo.len() as i32;
                for ri in &mt.atoms.resinfo {
                    atoms.resinfo.push(ri.clone());
                }
                for a in &mt.atoms.atom {
                    let mut a = a.clone();
                    a.resind += res_offset;
                    atoms.atom.push(a);
                }
            }
            if (nres as i32) <= max_res_renum {
                for j in 0..nmol {
                    for l in 0..nres {
                        max_resnr += 1;
                        atoms.resinfo[dest_res0 + j as usize * nres + l].nr = max_resnr;
                    }
                }
            }
        }
        atoms
    }
}

pub use crate::ir::InputRec;

/// The decoded body of a tpr file.
#[derive(Debug, Clone, Default)]
pub struct TprBody {
    pub boxm: Option<Matrix>,
    pub box_rel: Option<Matrix>,
    pub boxv: Option<Matrix>,
    pub mtop: Option<Mtop>,
    pub x: Option<Vec<Rvec>>,
    pub v: Option<Vec<Rvec>>,
    pub ir: Option<InputRec>,
    /// Byte offset of the inputrec section inside the body blob.
    pub ir_offset: usize,
    /// Byte offset of `nsteps` inside the body blob.
    pub nsteps_offset: usize,
    /// Symbol table index of the system (and single molecule type) name.
    pub mtop_name_symidx: i32,
    /// Byte offsets of the parts of the body that have to be reassembled when
    /// writing a modified topology (e.g. `convert-tpr -n`).
    pub moltype_count_offset: usize,
    pub molblock_count_offset: usize,
    pub natoms_offset: usize,
    /// Offset of the `bIntermolecularInteractions` flag (right after natoms).
    pub mtop_tail_offset: usize,
    pub mtop_end_offset: usize,
    /// Number of atom entries in each of the x/v sections.
    pub natoms: usize,
}

/// A run input file: XDR header plus the raw body blob.
#[derive(Debug, Clone)]
pub struct TprFile {
    pub header: TpxHeader,
    pub body: Vec<u8>,
}

fn read_matrix(r: &mut CReader) -> Result<Matrix> {
    let mut m = [[0.0f32; 3]; 3];
    for row in m.iter_mut() {
        for c in row.iter_mut() {
            *c = r.real()? as f32;
        }
    }
    Ok(m)
}

impl TprFile {
    /// Reads header and body.
    pub fn read(path: &str) -> Result<TprFile> {
        let data = std::fs::read(path)
            .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
        let file = TprFile::from_bytes(&data)?;
        // `do_tpxheader()` reports this on stderr for every tool that reads a
        // run input file.
        eprintln!(
            "Reading file {}, {} ({} precision)",
            path,
            file.header.version_string,
            if file.header.is_double { "double" } else { "single" }
        );
        Ok(file)
    }

    pub fn from_bytes(data: &[u8]) -> Result<TprFile> {
        let mut r = Reader::new(data);
        let version_string = r.gmx_string()?;
        if !version_string.starts_with("VERSION") {
            return Err(XdrError::Invalid(
                "this file is from a GROMACS version which is older than 2.0".into(),
            ));
        }
        let precision = r.int()?;
        let is_double = precision == 8;

        let mut h = TpxHeader {
            version_string,
            precision,
            is_double,
            ..Default::default()
        };
        h.file_version = r.int()?;
        if h.file_version >= 77 && h.file_version <= 79 {
            h.tag = r.gmx_string()?;
        }
        h.file_generation = r.int()?;
        if h.file_version >= 81 {
            h.tag = r.gmx_string()?;
        }
        h.natoms = r.int()?;
        h.ngtc = r.int()?;
        if h.file_version < 62 {
            let _ = r.int()?;
            let _ = r.real(is_double)?;
        }
        if h.file_version >= 79 {
            h.fep_state = r.int()?;
        }
        h.lambda = r.real(is_double)?;
        h.b_ir = r.bool()?;
        h.b_top = r.bool()?;
        h.b_x = r.bool()?;
        h.b_v = r.bool()?;
        h.b_f = r.bool()?;
        h.b_box = r.bool()?;
        if h.file_version >= TPXV_ADD_SIZE_FIELD
            && h.file_generation >= TPX_GENERATION_ADD_SIZE_FIELD
        {
            h.size_of_tpr_body = r.int64()?;
        }

        let body_start = r.position();
        let body_len = if h.size_of_tpr_body > 0 {
            h.size_of_tpr_body as usize
        } else {
            data.len() - body_start
        };
        if body_start + body_len > data.len() {
            return Err(XdrError::Truncated("tpr body"));
        }
        let body = data[body_start..body_start + body_len].to_vec();

        Ok(TprFile { header: h, body })
    }

    /// Serializes header and body.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut w = Writer::new();
        w.gmx_string(&self.header.version_string);
        w.int(self.header.precision);
        w.int(self.header.file_version);
        if self.header.file_version >= 77 && self.header.file_version <= 79 {
            w.gmx_string(&self.header.tag);
        }
        w.int(self.header.file_generation);
        if self.header.file_version >= 81 {
            w.gmx_string(&self.header.tag);
        }
        w.int(self.header.natoms);
        w.int(self.header.ngtc);
        if self.header.file_version < 62 {
            w.int(0);
            w.real(0.0, self.header.is_double);
        }
        if self.header.file_version >= 79 {
            w.int(self.header.fep_state);
        }
        w.real(self.header.lambda, self.header.is_double);
        w.bool(self.header.b_ir);
        w.bool(self.header.b_top);
        w.bool(self.header.b_x);
        w.bool(self.header.b_v);
        w.bool(self.header.b_f);
        w.bool(self.header.b_box);
        if self.header.file_version >= TPXV_ADD_SIZE_FIELD
            && self.header.file_generation >= TPX_GENERATION_ADD_SIZE_FIELD
        {
            w.int64(self.body.len() as i64);
        }
        w.opaque(&self.body);
        w.into_vec()
    }

    pub fn write(&self, path: &str) -> Result<()> {
        std::fs::write(path, self.to_bytes())
            .map_err(|e| XdrError::Invalid(format!("cannot write {path}: {e}")))?;
        Ok(())
    }

    /// Patches `nsteps` in the raw body and returns the updated file.
    pub fn set_nsteps(&self, body: &TprBody, nsteps: i64) -> Result<TprFile> {
        let mut out = self.clone();
        if body.ir_offset == 0 {
            return Err(XdrError::Invalid(
                "inputrec was not found in this tpr file".into(),
            ));
        }
        let off = body.nsteps_offset;
        if self.header.file_version >= TPXV_V62 {
            out.body[off..off + 8].copy_from_slice(&nsteps.to_be_bytes());
        } else {
            // Before version 62 the step counts are stored as 32 bit integers.
            out.body[off..off + 4].copy_from_slice(&(nsteps as i32).to_be_bytes());
        }
        Ok(out)
    }
}

/// Decodes the tpr body.
/// Decodes the tpr body.
pub fn parse_body(header: &TpxHeader, body: &[u8]) -> Result<TprBody> {
    if header.file_version < TPX_MIN_SUPPORTED_VERSION {
        return Err(XdrError::Invalid(format!(
            "reading tpx file version {} is not supported (the oldest supported \
             version is {}, matching GROMACS' own limit)",
            header.file_version, TPX_MIN_SUPPORTED_VERSION
        )));
    }
    // The body is stored as one opaque blob (serialized with the compact
    // in-memory serializer) only from `tpxv_AddSizeField` and
    // `TpxGeneration::AddSizeField` on; older files are plain XDR.
    let xdr_body = !(header.file_version >= TPXV_ADD_SIZE_FIELD
        && header.file_generation >= TPX_GENERATION_ADD_SIZE_FIELD);
    let mut r = CReader::new(body, header.is_double, xdr_body);
    let mut out = TprBody::default();

    // --- state, first part -------------------------------------------------
    if header.b_box {
        out.boxm = Some(read_matrix(&mut r)?);
        if header.file_version >= TPXV_V51 {
            out.box_rel = Some(read_matrix(&mut r)?);
        }
        out.boxv = Some(read_matrix(&mut r)?);
        if header.file_version < TPXV_V56 {
            // Unused extra matrix in very old files.
            let _ = read_matrix(&mut r)?;
        }
    }
    if header.ngtc > 0 {
        if header.file_version < TPXV_V69 {
            // Obsolete Berendsen temperature coupling lambdas.
            let _ = r.real_array(header.ngtc as usize)?;
        }
        let _ = r.real_array(header.ngtc as usize)?;
    }

    // --- global topology ---------------------------------------------------
    if header.b_top {
        trace("mtop start", &r);
        let mtop = parse_mtop(&mut r, header, &mut out)?;
        out.mtop = Some(mtop);
        trace("mtop end", &r);
    }
    out.mtop_end_offset = r.position();

    // --- state, second part ------------------------------------------------
    if header.b_x {
        trace("x start", &r);
        let mut x = Vec::with_capacity(header.natoms as usize);
        for _ in 0..header.natoms {
            x.push([r.real()? as f32, r.real()? as f32, r.real()? as f32]);
        }
        out.x = Some(x);
        trace("x end", &r);
    }
    if header.b_v {
        trace("v start", &r);
        let mut v = Vec::with_capacity(header.natoms as usize);
        for _ in 0..header.natoms {
            v.push([r.real()? as f32, r.real()? as f32, r.real()? as f32]);
        }
        out.v = Some(v);
        trace("v end", &r);
    }

    // --- inputrec ----------------------------------------------------------
    if header.b_ir {
        out.ir_offset = r.position();
        let (ir, nsteps_offset) = crate::ir::parse_inputrec(&mut r, header)?;
        out.nsteps_offset = nsteps_offset;
        out.ir = Some(ir);
        trace("inputrec prefix end", &r);
    }
    do_tpx_finalize(header, &mut out);
    Ok(out)
}

/// `do_tpx_finalize()`: the parts of the post-processing that change what
/// `gmx dump` prints.
fn do_tpx_finalize(header: &TpxHeader, body: &mut TprBody) {
    if header.file_version >= TPXV_V135 {
        return;
    }
    let ir = match body.ir.as_mut() {
        Some(ir) => ir,
        None => return,
    };
    let mtop = match body.mtop.as_ref() {
        Some(mtop) => mtop,
        None => return,
    };
    // `gmx_mtop_ftype_count()` over the position restraint interaction types.
    let mut n_posres = 0usize;
    for moltype in &mtop.moltypes {
        for ftype in [52usize, 53] {
            let list = moltype.ilists.get(ftype).map(|l| l.len()).unwrap_or(0);
            n_posres += list / (1 + interaction_function_nratoms_pub(ftype));
        }
    }
    // `RefCoordScaling::Com`
    if n_posres == 0 || ir.refcoord_scaling != 2 {
        ir.posres_com.clear();
        ir.posres_com_b.clear();
    }
}

/// Reads `do_symtab`.
fn parse_symtab(r: &mut CReader) -> Result<Vec<String>> {
    let nr = r.int()?;
    if nr < 0 {
        return Err(XdrError::Invalid("negative symbol table size".into()));
    }
    let mut v = Vec::with_capacity(nr as usize);
    for _ in 0..nr {
        v.push(r.string()?);
    }
    Ok(v)
}

fn symstr(r: &mut CReader, symtab: &[String]) -> Result<String> {
    let idx = r.int()?;
    Ok(symtab.get(idx as usize).cloned().unwrap_or_default())
}

/// `tpx_make_chain_identifiers()` from `fileio/confio.cpp`.
///
/// Every molecule with at least 15 atoms is given a chain identifier ('A',
/// 'B', ... then lower case, then digits); if only one identifier was handed
/// out all of them are blanked again.  This is what makes `gmx trjconv -s x.tpr`
/// write chain IDs into PDB output, while `gmx make_ndx` (which reads the
/// topology without molecule information) does not.
pub fn assign_chain_ids(atoms: &mut Atoms, ranges: &[(usize, usize)]) {
    const CHAIN_MIN_ATOMS: usize = 15;
    let mut next_chain_id = b'A';
    let mut out_of_ids = false;
    for &(begin, end) in ranges {
        let chain_id = if end - begin >= CHAIN_MIN_ATOMS && !out_of_ids {
            let id = next_chain_id;
            if next_chain_id == b'Z' {
                next_chain_id = b'a';
            } else if next_chain_id == b'z' {
                next_chain_id = b'0';
            } else if next_chain_id == b'9' {
                out_of_ids = true;
            } else {
                next_chain_id += 1;
            }
            id
        } else {
            b' '
        };
        for a in begin..end {
            if let Some(atom) = atoms.atom.get(a) {
                if let Some(ri) = atoms.resinfo.get_mut(atom.resind as usize) {
                    ri.chainid = chain_id;
                }
            }
        }
    }
    if next_chain_id == b'B' {
        for ri in atoms.resinfo.iter_mut() {
            ri.chainid = b' ';
        }
    }
}

/// `atomicnumber_to_element()`: the same deliberately incomplete table that
/// GROMACS uses when reading tpr files.
fn atomicnumber_to_element(n: i32) -> &'static str {
    match n {
        1 => "H",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        11 => "Na",
        12 => "Mg",
        15 => "P",
        16 => "S",
        17 => "Cl",
        18 => "Ar",
        19 => "K",
        20 => "Ca",
        25 => "Mn",
        26 => "Fe",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        35 => "Br",
        47 => "Ag",
        _ => "",
    }
}

fn parse_atoms(r: &mut CReader, symtab: &[String], file_version: i32) -> Result<Atoms> {
    let nr = r.int()? as usize;
    let nres = r.int()? as usize;
    let mut atom = Vec::with_capacity(nr);
    for _ in 0..nr {
        let mass = r.real()?;
        let charge = r.real()?;
        let _mass_b = r.real()?;
        let _charge_b = r.real()?;
        let type_id = r.ushort()?;
        let _type_b = r.ushort()?;
        let _ptype = r.int()?;
        let resind = r.int()?;
        let atomnumber = r.int()?;
        atom.push(Atom {
            name: String::new(),
            atom_type: String::new(),
            atom_type_b: String::new(),
            resind,
            mass,
            charge,
            mass_b: _mass_b,
            charge_b: _charge_b,
            ptype: _ptype,
            atomnumber,
            elem: atomicnumber_to_element(atomnumber).to_string(),
            type_id,
            type_id_b: _type_b,
        });
    }
    for i in 0..nr {
        atom[i].name = symstr(r, symtab)?;
    }
    for i in 0..nr {
        atom[i].atom_type = symstr(r, symtab)?;
    }
    for i in 0..nr {
        atom[i].atom_type_b = symstr(r, symtab)?;
    }
    let mut resinfo = Vec::with_capacity(nres);
    for j in 0..nres {
        let name = symstr(r, symtab)?;
        let (nr_, ic) = if file_version >= TPXV_V63 {
            (r.int()?, r.uchar()?)
        } else {
            (j as i32 + 1, b' ')
        };
        resinfo.push(ResInfo {
            name,
            nr: nr_,
            ic,
            chainid: b' ',
        });
    }
    Ok(Atoms {
        atom,
        resinfo,
        name: String::new(),
    })
}

/// `do_iparams`: advances the reader over the parameters of one function type.
/// `do_iparams()`: reads the parameter values of one function type.
/// `do_iparams()`: reads the parameter values of one function type, in the
/// order they are serialized (reals and ints interleaved for some types).
fn parse_iparams(r: &mut CReader, ftype: i32, file_version: i32) -> Result<Vec<PVal>> {
    let mut vals: Vec<PVal> = Vec::new();
    fn reals(vals: &mut Vec<PVal>, r: &mut CReader, n: usize) -> Result<()> {
        for _ in 0..n {
            vals.push(PVal::Real(r.real()?));
        }
        Ok(())
    }
    fn ints(vals: &mut Vec<PVal>, r: &mut CReader, n: usize) -> Result<()> {
        for _ in 0..n {
            vals.push(PVal::Int(r.int()? as i64));
        }
        Ok(())
    }
    match ftype {
        // Bonds, GROMOS96Bonds, HarmonicPotential, Angles, GROMOS96Angles,
        // ImproperDihedrals.
        0 | 1 | 5 | 10 | 11 | 24 => reals(&mut vals, r, 4)?,
        // RestrictedBendingPotential: the B state was added in version 134.
        12 => reals(&mut vals, r, if file_version >= TPXV_V134 { 4 } else { 2 })?,
        13 => reals(&mut vals, r, 4)?, // LinearAngles
        6 => reals(&mut vals, r, 2)?,  // FENEBonds
        9 => reals(&mut vals, r, 8)?,  // RestraintBonds
        // Tabulated bonds/angles/dihedrals: kA, table, kB
        7 | 8 | 18 | 26 => {
            vals.push(PVal::Real(r.real()?));
            vals.push(PVal::Int(r.int()? as i64));
            vals.push(PVal::Real(r.real()?));
        }
        14 => reals(&mut vals, r, 3)?, // CrossBondBonds
        15 => reals(&mut vals, r, 4)?, // CrossBondAngles
        // UreyBradleyPotential: B state added in version 79.
        16 => reals(&mut vals, r, if file_version >= TPXV_V79 { 8 } else { 4 })?,
        17 => reals(&mut vals, r, 6)?, // QuarticAngles (theta + 5 coefficients)
        38 => reals(&mut vals, r, 3)?, // BuckinghamShortRange
        // MorsePotential: B state added in version 79.
        2 => reals(&mut vals, r, if file_version >= TPXV_V79 { 6 } else { 3 })?,
        3 => reals(&mut vals, r, 3)?, // CubicBonds
        4 => {}                       // ConnectBonds
        48 => reals(&mut vals, r, 1)?, // Polarization
        51 => reals(&mut vals, r, 3)?, // AnharmonicPolarization
        49 => reals(&mut vals, r, 6)?, // WaterPolarization
        // TholePolarization: the unused rfac parameter was removed in 127.
        50 => reals(&mut vals, r, if file_version < TPXV_V127 { 4 } else { 3 })?,
        37 => reals(&mut vals, r, 2)?, // LennardJonesShortRange
        33 => reals(&mut vals, r, 4)?, // LennardJones14
        35 => reals(&mut vals, r, 5)?, // LennardJonesCoulomb14Q
        36 => reals(&mut vals, r, 4)?, // LennardJonesCoulombNonBondedPairs
        // Proper/periodic improper dihedrals, angle restraints
        19 | 25 | 58 | 59 => {
            reals(&mut vals, r, 4)?;
            vals.push(PVal::Int(r.int()? as i64));
        }
        // RestrictedTorsionPotential: B state added in version 134.
        21 => reals(&mut vals, r, if file_version >= TPXV_V134 { 4 } else { 2 })?,
        54 => {
            // DistanceRestraints
            ints(&mut vals, r, 2)?;
            reals(&mut vals, r, 4)?;
        }
        56 => {
            // OrientationRestraints
            ints(&mut vals, r, 3)?;
            reals(&mut vals, r, 3)?;
        }
        // DihedralRestraints: two obsolete ints before version 82.
        60 => {
            if file_version >= TPXV_V82 {
                reals(&mut vals, r, 6)?;
            } else {
                ints(&mut vals, r, 2)?;
                reals(&mut vals, r, 3)?;
            }
        }
        52 => reals(&mut vals, r, 12)?, // PositionRestraints (4 rvecs)
        53 => {
            // FlatBottomedPositionRestraints
            vals.push(PVal::Int(r.int()? as i64));
            reals(&mut vals, r, 5)?;
        }
        // CombinedBendingTorsion: B state added in version 134.
        22 => reals(&mut vals, r, if file_version >= TPXV_V134 { 2 * NR_DIHEDRAL_PARAMS } else { NR_DIHEDRAL_PARAMS })?,
        20 | 23 => reals(&mut vals, r, 2 * NR_DIHEDRAL_PARAMS)?, // RB / Fourier dihedrals
        62 | 63 => reals(&mut vals, r, 2)?, // Constraints
        64 => reals(&mut vals, r, 2)?,      // SETTLE
        65 => {}                            // VirtualSite1
        66 | 67 => reals(&mut vals, r, 1)?, // VirtualSite2 variants
        // VirtualSite3, VirtualSite3FlexibleDistance,
        // VirtualSite3FlexibleAngleDistance: `a`, `b`
        68 | 69 | 70 => reals(&mut vals, r, 2)?,
        // VirtualSite3Outside, VirtualSite4FlexibleDistance,
        // VirtualSite4FlexibleDistanceNormalization: `a`, `b`, `c`
        71 | 72 | 73 => reals(&mut vals, r, 3)?,
        74 => {
            // VirtualSiteN
            vals.push(PVal::Int(r.int()? as i64));
            reals(&mut vals, r, 1)?;
        }
        28 | 29 | 30 => {
            // Implicit solvent parameters: read (and ignored) for old files.
            if file_version < 68 {
                reals(&mut vals, r, 4)?;
            }
            if file_version < TPXV_V113 {
                reals(&mut vals, r, 5)?;
            }
        }
        27 => {
            // DihedralEnergyCorrectionMap
            ints(&mut vals, r, 2)?;
        }
        other => {
            return Err(XdrError::Invalid(format!(
                "unknown interaction function type {other} in force field parameters"
            )))
        }
    }
    Ok(vals)
}

fn parse_ffparams(r: &mut CReader, file_version: i32) -> Result<FfParams> {
    let atnr = r.int()?;
    let num_types = r.int()?;
    if num_types < 0 {
        return Err(XdrError::Invalid("negative force field type count".into()));
    }
    let mut functype = r.int_array(num_types as usize)?;
    // `reppow` only exists in files from version 66 on.
    let reppow = if file_version >= TPXV_V66 {
        r.double()?
    } else {
        12.0
    };
    let fudge_qq = r.real()?;
    // Shift the function types of files written before a type was inserted
    // into the enum (see `ftupd[]` in tpxio.cpp).
    for ft in functype.iter_mut() {
        for (fvnr, ftype) in FTUPD {
            if file_version < *fvnr && *ft >= *ftype {
                *ft += 1;
            }
        }
    }
    let mut iparams = Vec::with_capacity(functype.len());
    for ft in &functype {
        iparams.push(parse_iparams(r, *ft, file_version)?);
    }
    Ok(FfParams {
        atnr,
        functype,
        reppow,
        fudge_qq,
        iparams,
    })
}

/// True when a list of this function type is absent from the file because the
/// type did not exist yet when the file was written.
fn ilist_is_absent(iftype: usize, file_version: i32) -> bool {
    FTUPD
        .iter()
        .any(|(fvnr, ftype)| file_version < *fvnr && iftype as i32 == *ftype)
}

/// Interaction types with the `IF_CHEMBOND` flag and their atom counts, taken
/// from the function table in `gromacs/topology/ifunc.cpp`.
const CHEMBOND_TYPES: &[(usize, usize)] = &[
    (0, 2),  // Bonds
    (1, 2),  // GROMOS96Bonds
    (2, 2),  // MorsePotential
    (3, 2),  // CubicBonds
    (4, 2),  // ConnectBonds
    (6, 2),  // FENEBonds
    (7, 2),  // TabulatedBonds
    (48, 2), // Polarization
    (51, 2), // AnharmonicPolarization
    (62, 2), // Constraints
    (64, 3), // SETTLE
];

fn parse_ilists(
    r: &mut CReader,
    file_version: i32,
    mut bonds: Option<&mut Vec<(u32, u32)>>,
    mut all: Option<&mut Vec<Vec<i32>>>,
) -> Result<()> {
    for iftype in 0..INTERACTION_FUNCTION_COUNT {
        if ilist_is_absent(iftype, file_version) {
            // `do_ilists()` clears the list without reading anything, so the
            // entries stay indexed by interaction function type.
            if let Some(all) = all.as_deref_mut() {
                all.push(Vec::new());
            }
            continue;
        }
        let nr = r.int()?;
        if nr < 0 {
            return Err(XdrError::Invalid("negative interaction list size".into()));
        }
        let atoms = r.int_array(nr as usize)?;
        if let Some(bonds) = bonds.as_deref_mut() {
            if let Some((_, nratoms)) = CHEMBOND_TYPES.iter().find(|(t, _)| *t == iftype).map(|(t, n)| (*t, *n)) {
                let mut i = 0usize;
                while i + nratoms + 1 <= atoms.len() {
                    if iftype == 64 {
                        // SETTLE: bond the first atom with all others.  The
                        // atoms start at i + 1 (i holds the parameter index).
                        for j in 1..nratoms {
                            if atoms[i + 1] >= 0 && atoms[i + j + 1] >= 0 {
                                bonds.push((atoms[i + 1] as u32, atoms[i + j + 1] as u32));
                            }
                        }
                    } else {
                        for j in 1..nratoms {
                            if atoms[i + j] >= 0 && atoms[i + j + 1] >= 0 {
                                bonds.push((atoms[i + j] as u32, atoms[i + j + 1] as u32));
                            }
                        }
                    }
                    i += nratoms + 1;
                }
            }
        }
        if let Some(all) = all.as_deref_mut() {
            all.push(atoms);
        }
    }
    Ok(())
}

fn parse_list_of_lists(r: &mut CReader) -> Result<Vec<Vec<i32>>> {
    let num_lists = r.int()? as usize;
    let num_elements = r.int()? as usize;
    let ranges = r.int_array(num_lists + 1)?;
    let elements = r.int_array(num_elements)?;
    let mut out = Vec::with_capacity(num_lists);
    for i in 0..num_lists {
        let a = ranges[i] as usize;
        let b = ranges[i + 1] as usize;
        out.push(elements[a..b].to_vec());
    }
    Ok(out)
}

fn parse_moltype(r: &mut CReader, symtab: &[String], file_version: i32) -> Result<MolType> {
    let name = symstr(r, symtab)?;
    let atoms = parse_atoms(r, symtab, file_version)?;
    let mut bonds = Vec::new();
    let mut ilists = Vec::new();
    parse_ilists(r, file_version, Some(&mut bonds), Some(&mut ilists))?;
    // charge groups (obsolete): int nr, int[nr + 1]
    let nr = r.int()?;
    let _ = r.int_array(nr as usize + 1)?;
    let excls = parse_list_of_lists(r)?;
    Ok(MolType {
        name,
        atoms,
        excls,
        bonds,
        ilists,
    })
}

fn parse_mtop(r: &mut CReader, header: &TpxHeader, out: &mut TprBody) -> Result<Mtop> {
    let symtab = parse_symtab(r)?;
    trace(&format!("symtab: {} strings", symtab.len()), r);
    out.mtop_name_symidx = r.int()?;
    let name = symtab
        .get(out.mtop_name_symidx as usize)
        .cloned()
        .unwrap_or_default();
    let ffparams = parse_ffparams(r, header.file_version)?;
    trace("ffparams done", r);

    out.moltype_count_offset = r.position();
    let nmoltype = r.int()?;
    trace(&format!("nmoltype={nmoltype}"), r);
    let mut moltypes = Vec::with_capacity(nmoltype.max(0) as usize);
    for _ in 0..nmoltype {
        moltypes.push(parse_moltype(r, &symtab, header.file_version)?);
        trace("moltype done", r);
    }
    let nmolblock = r.int()?;
    out.molblock_count_offset = r.position() - 4;
    trace(&format!("nmolblock={nmolblock}"), r);
    let mut molblocks = Vec::with_capacity(nmolblock.max(0) as usize);
    for _ in 0..nmolblock {
        let moltype_index = r.int()?;
        let nmol = r.int()?;
        let _num_atoms_per_molecule = r.int()?;
        let n_posres_a = r.int()?;
        let mut posres_xa = Vec::new();
        if n_posres_a > 0 {
            for _ in 0..n_posres_a {
                posres_xa.push([r.real()? as f32, r.real()? as f32, r.real()? as f32]);
            }
        }
        let n_posres_b = r.int()?;
        let mut posres_xb = Vec::new();
        if n_posres_b > 0 {
            for _ in 0..n_posres_b {
                posres_xb.push([r.real()? as f32, r.real()? as f32, r.real()? as f32]);
            }
        }
        molblocks.push(MolBlock {
            moltype_index,
            nmol,
            posres_xa,
            posres_xb,
        });
    }
    out.natoms_offset = r.position();
    let natoms = r.int()? as usize;
    trace(&format!("mtop natoms={natoms}"), r);
    out.natoms = natoms;
    out.mtop_tail_offset = r.position();

    let mut b_intermolecular = false;
    let mut intermolecular_ilists = Vec::new();
    if header.file_version >= TPXV_V103 {
        b_intermolecular = r.bool()?;
        if b_intermolecular {
            parse_ilists(r, header.file_version, None, Some(&mut intermolecular_ilists))?;
        }
    }

    // `do_atomtypes()` (removed from the format in version 128)
    if header.file_version < TPXV_V128 {
        let nr = r.int()?;
        if nr > 0 {
            if header.file_version < TPXV_V113 {
                let _ = r.real_array(nr as usize * 3)?;
            }
            let _ = r.int_array(nr as usize)?;
            if header.file_version >= TPXV_V60 && header.file_version < TPXV_V113 {
                let _ = r.real_array(nr as usize * 2)?;
            }
        }
    }

    // cmap
    let mut cmap_grid_spacing = 0i32;
    let mut cmap_data: Vec<Vec<f32>> = Vec::new();
    if header.file_version >= TPXV_V65 {
        let ngrid = r.int()?;
        let grid_spacing = r.int()?;
        cmap_grid_spacing = grid_spacing;
        let nelem = (grid_spacing * grid_spacing) as usize;
        for _ in 0..ngrid.max(0) {
            let values = r.real_array(nelem * 4)?;
            cmap_data.push(values.into_iter().map(|v| v as f32).collect());
        }
    }
    trace("cmap done", r);

    // groups
    let mut groups = Vec::with_capacity(SIMULATION_ATOM_GROUP_COUNT);
    for _ in 0..SIMULATION_ATOM_GROUP_COUNT {
        let size = r.int()?;
        let arr = if size > 0 {
            r.int_array(size as usize)?
        } else {
            Vec::new()
        };
        groups.push(arr);
    }
    let n_group_names = r.int()?;
    let mut group_names = Vec::with_capacity(n_group_names.max(0) as usize);
    for _ in 0..n_group_names {
        group_names.push(symstr(r, &symtab)?);
    }
    let mut group_numbers = Vec::with_capacity(SIMULATION_ATOM_GROUP_COUNT);
    for _ in 0..SIMULATION_ATOM_GROUP_COUNT {
        let n = r.int()?;
        let mut v = Vec::new();
        if n > 0 {
            for _ in 0..n {
                v.push(r.uchar()?);
            }
        }
        group_numbers.push(v);
    }
    trace("groups done", r);

    // intermolecular exclusion group (added in version 120)
    if header.file_version >= TPXV_V120 {
        let excl_size = r.int64()?;
        if excl_size > 0 {
            let _ = r.int_array(excl_size as usize)?;
        }
    }
    trace("mtop end", r);

    let _ = header;
    Ok(Mtop {
        symtab,
        name,
        natoms,
        moltypes,
        molblocks,
        ffparams,
        b_intermolecular,
        intermolecular_ilists,
        groups,
        group_names,
        group_numbers,
        cmap_grid_spacing,
        cmap_data,
    })
}

/// Serializes a symbol reference: `do_symstr()` writes the index of the string
/// in the symbol table.
/// Writer for the compact (`gmx::InMemorySerializer`) body format: values use
/// their native width, `bool`/`uchar` are one byte, `ushort` two bytes, and
/// nothing is padded.
pub struct CWriter {
    pub data: Vec<u8>,
    double_precision: bool,
    xdr: bool,
}

impl CWriter {
    pub fn new(double_precision: bool, xdr: bool) -> Self {
        CWriter {
            data: Vec::new(),
            double_precision,
            xdr,
        }
    }

    pub fn int(&mut self, v: i32) {
        self.data.extend_from_slice(&v.to_be_bytes());
    }

    pub fn bool(&mut self, v: bool) {
        if self.xdr {
            self.int(if v { 1 } else { 0 });
        } else {
            self.data.push(if v { 1 } else { 0 });
        }
    }

    pub fn uchar(&mut self, v: u8) {
        if self.xdr {
            self.int(v as i32);
        } else {
            self.data.push(v);
        }
    }

    pub fn ushort(&mut self, v: u16) {
        if self.xdr {
            self.int(v as i32);
        } else {
            self.data.extend_from_slice(&v.to_be_bytes());
        }
    }

    pub fn real(&mut self, v: f64) {
        if self.double_precision {
            self.data.extend_from_slice(&v.to_be_bytes());
        } else {
            self.data.extend_from_slice(&(v as f32).to_be_bytes());
        }
    }

    pub fn bytes(&mut self, b: &[u8]) {
        self.data.extend_from_slice(b);
    }

    pub fn into_vec(self) -> Vec<u8> {
        self.data
    }

    /// `doString()` in the XDR serializer: length, then the XDR string.
    pub fn string(&mut self, s: &str) {
        if self.xdr {
            self.int(s.len() as i32 + 1);
            self.int(s.len() as i32);
            self.data.extend_from_slice(s.as_bytes());
            for _ in 0..(4 - (s.len() % 4)) % 4 {
                self.data.push(0);
            }
        } else {
            self.data
                .extend_from_slice(&(s.len() as u64).to_be_bytes());
            self.data.extend_from_slice(s.as_bytes());
        }
    }
}

fn put_symstr(w: &mut CWriter, symtab: &std::collections::HashMap<String, i32>, s: &str) {
    w.int(symtab.get(s).copied().unwrap_or(0));
}

fn put_atoms(w: &mut CWriter, atoms: &Atoms, symtab: &std::collections::HashMap<String, i32>) {
    w.int(atoms.nr() as i32);
    w.int(atoms.resinfo.len() as i32);
    // do_atom(): m, q, mB, qB, type, typeB, ptype, resind, atomnumber
    for a in &atoms.atom {
        w.real(a.mass);
        w.real(a.charge);
        w.real(a.mass_b);
        w.real(a.charge_b);
        w.ushort(a.type_id);
        w.ushort(a.type_id_b);
        w.int(a.ptype);
        w.int(a.resind);
        w.int(a.atomnumber);
    }
    for a in &atoms.atom {
        put_symstr(w, symtab, &a.name);
    }
    for a in &atoms.atom {
        put_symstr(w, symtab, &a.atom_type);
    }
    for a in &atoms.atom {
        put_symstr(w, symtab, &a.atom_type_b);
    }
    // do_resinfo(): symstr(name), nr, ic
    for ri in &atoms.resinfo {
        put_symstr(w, symtab, &ri.name);
        w.int(ri.nr);
        w.uchar(if ri.ic == 0 { b' ' } else { ri.ic });
    }
}

fn put_ilists(w: &mut CWriter, ilists: &[Vec<i32>], file_version: i32) {
    for i in 0..INTERACTION_FUNCTION_COUNT {
        // Lists of interaction types that did not exist yet are not stored at
        // all in older files, so they must not be written either.
        if ilist_is_absent(i, file_version) {
            continue;
        }
        let list = ilists.get(i).map(|l| l.as_slice()).unwrap_or(&[]);
        w.int(list.len() as i32);
        for v in list {
            w.int(*v);
        }
    }
}

/// `reduce_ilist()`: keeps the interactions whose particles are all kept and
/// renumbers them.
fn reduce_ilist(list: &[i32], invindex: &[i32], b_keep: &[bool], nratoms: usize) -> Vec<i32> {
    let mut out = Vec::new();
    if list.is_empty() || nratoms == 0 {
        return out;
    }
    let mut i = 0usize;
    while i + nratoms < list.len() {
        let mut keep = true;
        for j in 0..nratoms {
            let a = list[i + 1 + j];
            if a < 0 || (a as usize) < b_keep.len() && !b_keep[a as usize] {
                keep = false;
            }
            if a < 0 || (a as usize) >= b_keep.len() {
                keep = false;
            }
        }
        if keep {
            out.push(list[i]);
            for j in 0..nratoms {
                out.push(invindex[list[i + 1 + j] as usize]);
            }
        }
        i += nratoms + 1;
    }
    out
}

/// `reduce_topology_x()` from `convert_tpr.cpp`, expressed as a rewrite of the
/// serialized body.
///
/// The force field parameters, the symbol table and everything after the
/// molecule blocks are copied verbatim; only the atom list, the interaction
/// lists, the exclusions and the molecule block are rebuilt for the selection.
pub fn write_subset_body(
    tpr: &TprFile,
    body: &TprBody,
    selection: &[usize],
) -> Result<TprFile> {
    let mtop = body
        .mtop
        .as_ref()
        .ok_or_else(|| XdrError::Invalid("the tpr file has no topology".into()))?;
    if body.moltype_count_offset == 0 {
        return Err(XdrError::Invalid(
            "this tpr file could not be decoded far enough for a subset".into(),
        ));
    }
    let natoms = body.natoms;
    let mut b_keep = vec![false; natoms];
    let mut invindex = vec![-1i32; natoms];
    for (i, &a) in selection.iter().enumerate() {
        if a >= natoms {
            return Err(XdrError::Invalid(format!(
                "index {a} is larger than the number of atoms in the tpr file ({natoms})"
            )));
        }
        b_keep[a] = true;
        invindex[a] = i as i32;
    }

    let symtab = mtop.symtab.clone();
    let symidx: std::collections::HashMap<String, i32> = symtab
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i as i32))
        .collect();

    // --- new atoms ---------------------------------------------------------
    let all_atoms = mtop.global_atoms();
    let mut atoms = Atoms {
        name: mtop.name.clone(),
        atom: Vec::with_capacity(selection.len()),
        resinfo: Vec::new(),
    };
    // reduce_atom(): the residue of the selected atom is copied whenever its
    // original `resind` differs from the *new* index of the previous atom.
    // Because the previous atom was already renumbered, this normally produces
    // one residue per atom for selections that do not start at the first
    // residue; the quirk is reproduced here for compatibility.
    let mut nr: i32 = -1;
    for (i, &a) in selection.iter().enumerate() {
        let mut atom = all_atoms.atom[a].clone();
        let orig_resind = atom.resind;
        if i == 0 || orig_resind != nr {
            nr += 1;
            atoms
                .resinfo
                .push(all_atoms.resinfo[orig_resind as usize].clone());
        }
        atom.resind = nr;
        atoms.atom.push(atom);
    }
    // `reduce_atom()` only rewrites `atom`, `atomname` and `resinfo`; the
    // `atomtype`/`atomtypeB` arrays keep the global order, so a selection that
    // is not a prefix of the system receives the types of the first `gnx`
    // global atoms.  This is reproduced verbatim (including for selections
    // where it yields mismatched types), matching `gmx convert-tpr`.
    for (i, atom) in atoms.atom.iter_mut().enumerate() {
        if let Some(src) = all_atoms.atom.get(i) {
            atom.atom_type = src.atom_type.clone();
            atom.atom_type_b = src.atom_type_b.clone();
        }
    }

    // --- new interaction lists --------------------------------------------
    // `reduce_topology_x()` first builds a local topology for the whole
    // system (`gmx_mtop_generate_local_top`), which concatenates the
    // interaction lists of every molecule with the atom indices offset, and
    // only then applies the selection.
    let mut global_ilists: Vec<Vec<i32>> = vec![Vec::new(); INTERACTION_FUNCTION_COUNT];
    let mut global_excls: Vec<Vec<i32>> = vec![Vec::new(); natoms];
    {
        let mut offset = 0usize;
        for mb in &mtop.molblocks {
            let Some(mt) = mtop.moltypes.get(mb.moltype_index as usize) else {
                continue;
            };
            let srcnr = mt.atoms.nr();
            for _ in 0..mb.nmol.max(0) {
                for (i, list) in mt.ilists.iter().enumerate() {
                    if i >= INTERACTION_FUNCTION_COUNT {
                        break;
                    }
                    let nratoms = interaction_function_nratoms(i);
                    let mut k = 0usize;
                    while nratoms > 0 && k + nratoms < list.len() {
                        global_ilists[i].push(list[k]);
                        for j in 0..nratoms {
                            global_ilists[i].push(list[k + 1 + j] + offset as i32);
                        }
                        k += nratoms + 1;
                    }
                }
                for (a, list) in mt.excls.iter().enumerate() {
                    if offset + a < natoms {
                        global_excls[offset + a] = list
                            .iter()
                            .filter(|c| **c >= 0)
                            .map(|c| *c + offset as i32)
                            .collect();
                    }
                }
                offset += srcnr;
            }
        }
    }

    let mut ilists: Vec<Vec<i32>> = Vec::with_capacity(INTERACTION_FUNCTION_COUNT);
    // `gen_local_top()` is called with mergeConstr = true, which concatenates
    // the "constraints without coupling" list onto the plain constraints list.
    let mut constraints: Vec<i32> = Vec::new();
    for i in 0..INTERACTION_FUNCTION_COUNT {
        let nratoms = interaction_function_nratoms(i);
        let reduced = reduce_ilist(&global_ilists[i], &invindex, &b_keep, nratoms);
        if i == 62 {
            constraints.extend_from_slice(&reduced);
            ilists.push(Vec::new());
        } else if i == 63 {
            constraints.extend_from_slice(&reduced);
            ilists.push(Vec::new());
        } else {
            ilists.push(reduced);
        }
    }
    ilists[62] = constraints;

    // --- new exclusions ----------------------------------------------------
    let mut excls: Vec<Vec<i32>> = Vec::new();
    for &a in selection {
        let mut list = Vec::new();
        if let Some(src) = global_excls.get(a) {
            for &j in src {
                if j >= 0 && (j as usize) < natoms && b_keep[j as usize] {
                    list.push(invindex[j as usize]);
                }
            }
        }
        excls.push(list);
    }

    // --- serialize ---------------------------------------------------------
    let dp = tpr.header.is_double;
    let xdr_body = !(tpr.header.file_version >= TPXV_ADD_SIZE_FIELD
        && tpr.header.file_generation >= TPX_GENERATION_ADD_SIZE_FIELD);
    let mut w = CWriter::new(dp, xdr_body);
    // state, first part, and the symbol table + name + force field parameters
    w.bytes(&tpr.body[..body.moltype_count_offset]);
    // a single molecule type containing the whole selection
    w.int(1);
    w.int(body.mtop_name_symidx);
    put_atoms(&mut w, &atoms, &symidx);
    put_ilists(&mut w, &ilists, tpr.header.file_version);
    // obsolete charge groups: one per atom
    w.int(atoms.nr() as i32);
    for i in 0..=atoms.nr() {
        w.int(i as i32);
    }
    // exclusions as a ListOfLists
    w.int(excls.len() as i32);
    let total: usize = excls.iter().map(|l| l.len()).sum();
    w.int(total as i32);
    let mut acc = 0i32;
    w.int(0);
    for l in &excls {
        acc += l.len() as i32;
        w.int(acc);
    }
    for l in &excls {
        for v in l {
            w.int(*v);
        }
    }
    // a single molecule block
    w.int(1);
    let mb = mtop.molblocks.first().cloned().unwrap_or_default();
    w.int(0);
    w.int(1);
    w.int(atoms.nr() as i32);
    w.int(mb.posres_xa.len() as i32);
    for v in &mb.posres_xa {
        for c in v {
            w.real(*c as f64);
        }
    }
    w.int(mb.posres_xb.len() as i32);
    for v in &mb.posres_xb {
        for c in v {
            w.real(*c as f64);
        }
    }
    w.int(atoms.nr() as i32);
    // everything after natoms (intermolecular interactions, cmap, groups,
    // exclusion group) is copied unchanged
    w.bytes(&tpr.body[body.mtop_tail_offset..body.mtop_end_offset]);

    // --- state, second part: the selected coordinates and velocities -------
    let stride = 3 * real_size(tpr.header.is_double);
    let mut pos = body.mtop_end_offset;
    if tpr.header.b_x {
        for &a in selection {
            let off = pos + a * stride;
            w.bytes(&tpr.body[off..off + stride]);
        }
        pos += natoms * stride;
    }
    if tpr.header.b_v {
        for &a in selection {
            let off = pos + a * stride;
            w.bytes(&tpr.body[off..off + stride]);
        }
    }
    // the inputrec is unchanged
    w.bytes(&tpr.body[body.ir_offset..]);

    let mut out = tpr.clone();
    out.header.natoms = atoms.nr() as i32;
    out.header.size_of_tpr_body = 0;
    out.body = w.into_vec();
    Ok(out)
}

fn real_size(double_precision: bool) -> usize {
    if double_precision {
        8
    } else {
        4
    }
}

/// Number of particles of each interaction function (`interaction_function[]`
/// in `topology/ifunc.cpp`).
pub fn interaction_function_nratoms_pub(iftype: usize) -> usize {
    interaction_function_nratoms(iftype)
}

fn interaction_function_nratoms(iftype: usize) -> usize {
    match iftype {
        0..=9 => 2,      // bonds
        10..=18 => 3,    // angles
        19..=26 => 4,    // dihedrals
        27 => 5,         // DihedralEnergyCorrectionMap
        33 | 35 | 36 | 37 | 38 => 2,
        48 => 2,         // Polarization
        49 => 5,         // WaterPolarization
        50 => 4,         // TholePolarization
        51 => 2,         // AnharmonicPolarization
        52 | 53 => 1,    // position restraints
        54 => 2,         // DistanceRestraints
        56 => 2,         // OrientationRestraints
        58 => 4,         // AngleRestraints
        59 => 2,         // AngleZAxisRestraints
        60 => 4,         // DihedralRestraints
        62 | 63 => 2,    // constraints
        64 => 3,         // SETTLE
        65 => 2,         // VirtualSite1
        66 | 67 => 3,    // VirtualSite2 variants
        68 | 69 | 70 | 71 => 4, // VirtualSite3 variants
        72 | 73 => 5,    // VirtualSite4 variants
        74 => 2,         // VirtualSiteN
        _ => 0,
    }
}
