//! `gmx trjconv`, mirroring the frame loop of
//! `src/gromacs/tools/trjconv.cpp`.

use std::fmt::Write as _;

use crate::cmd::Args;
use crate::frame::{Atoms, Frame, Matrix, Rvec};
use crate::index::{self, IndexGroup};
use crate::pbc::{self, ECenter};
use crate::progress::Progress;
use crate::tpr;
use crate::trx::{self, TrxFormat, TrxWriter};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PbcMode {
    None,
    Mol,
    Res,
    Atom,
    NoJump,
    Cluster,
    Whole,
}

impl PbcMode {
    fn parse(s: &str) -> PbcMode {
        match s {
            "mol" => PbcMode::Mol,
            "res" => PbcMode::Res,
            "atom" => PbcMode::Atom,
            "nojump" => PbcMode::NoJump,
            "cluster" => PbcMode::Cluster,
            "whole" => PbcMode::Whole,
            _ => PbcMode::None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum UrMode {
    Rect,
    Tric,
    Compact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FitMode {
    None,
    Fit,
    FitXy,
    Translation,
    TransXy,
    Progressive,
}

use crate::cmd::select_group as require_group;

/// The frame progress of a running conversion.
///
/// Progress is measured against the requested time window when one was given
/// (`-b`/`-e`), and against the position in the input file otherwise, which is
/// what "converting the whole trajectory" means.
struct FrameProgress {
    bar: Progress,
    time_window: Option<(f64, f64)>,
}

impl FrameProgress {
    fn new(time_window: Option<(f64, f64)>) -> Self {
        FrameProgress {
            bar: Progress::new(),
            time_window,
        }
    }

    /// Reports the frame that has just been read.
    ///
    /// The state is updated for every frame so that the final draw is exact;
    /// `Progress` throttles the redraws itself.
    fn tick(&mut self, source: &mut trx::FrameSource, frames: u64, time: f64) {
        let fraction = match self.time_window {
            // A window with an unbounded end (`-e inf`, which is how "the whole
            // trajectory" is passed in) carries no information about how far
            // along the conversion is, so the position in the file is used.
            Some((start, end)) if end.is_finite() && end > start => {
                Some(((time - start) / (end - start)).clamp(0.0, 1.0))
            }
            _ => source.read_progress().and_then(|p| p.fraction()),
        };
        self.bar.update(fraction, frames, time);
    }

    fn suspend(&mut self) {
        self.bar.pause();
    }

    fn finish(&mut self) {
        self.bar.finish();
    }
}

fn rmod(x: f64, first: f64, step: f64) -> bool {
    if step == 0.0 {
        return false;
    }
    let r = (x - first) % step;
    r.abs() < 0.5 * step * 1e-5 || (step - r.abs()).abs() < 0.5 * step * 1e-5
}

/// Virtual site interaction types and their atom counts (`def_vsite` entries
/// in `gromacs/topology/ifunc.cpp`), indexed by `InteractionFunction`.
const VSITE_TYPES: &[(usize, usize)] = &[
    (65, 2),  // VirtualSite1
    (66, 3),  // VirtualSite2
    (67, 3),  // VirtualSite2FlexibleDistance
    (68, 4),  // VirtualSite3
    (69, 4),  // VirtualSite3FlexibleDistance
    (70, 4),  // VirtualSite3FlexibleAngleDistance
    (71, 4),  // VirtualSite3Outside
    (72, 5),  // VirtualSite4FlexibleDistance
    (73, 5),  // VirtualSite4FlexibleDistanceNormalization
    (74, 2),  // VirtualSiteN
];

/// Prints the same warning GROMACS emits when molecules cannot be made whole
/// because no topology was supplied.
fn warn_no_topology() {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        eprintln!(
            "\nWARNING: If there are molecules in the input trajectory file\n\
             \x20        that are broken across periodic boundaries, they\n\
             \x20        cannot be made whole (or treated as whole) without\n\
             \x20        you providing a run input file.\n"
        );
    });
}

/// `mk_graph_moltype()`: the chemical bonds plus, when the bond graph falls
/// apart into several parts, the interactions that reconnect them (virtual
/// sites, position restraints, ...).  Virtual sites are the case that occurs
/// in practice, so those are the ones handled here.
fn build_moltype_adjacency(mt: &tpr::MolType) -> Vec<Vec<u32>> {
    let n = mt.atoms.nr();
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for (a, b) in &mt.bonds {
        let (a, b) = (*a as usize, *b as usize);
        if a < n && b < n && a != b {
            adj[a].push(b as u32);
            adj[b].push(a as u32);
        }
    }
    if mt.bonds.is_empty() {
        // No bonded interactions at all (e.g. a file without a topology): fall
        // back to the exclusion lists.
        for (i, list) in mt.excls.iter().enumerate() {
            for &c in list {
                if c >= 0 && (c as usize) < n && c as usize != i {
                    adj[i].push(c as u32);
                }
            }
        }
        return adj;
    }

    // Connected parts of the chemical bond graph.
    let mut part: Vec<usize> = (0..n).collect();
    fn find(part: &mut Vec<usize>, mut x: usize) -> usize {
        while part[x] != x {
            part[x] = part[part[x]];
            x = part[x];
        }
        x
    }
    for (a, b) in &mt.bonds {
        let (ra, rb) = (find(&mut part, *a as usize), find(&mut part, *b as usize));
        if ra != rb {
            part[ra] = rb;
        }
    }
    let parts_snapshot = part.clone();
    let root = |p: &Vec<usize>, mut x: usize| {
        while p[x] != x {
            x = p[x];
        }
        x
    };
    for (iftype, nratoms) in VSITE_TYPES {
        let Some(list) = mt.ilists.get(*iftype) else {
            continue;
        };
        let mut i = 0usize;
        while i + nratoms + 1 <= list.len() {
            for j in 1..*nratoms {
                let (a, b) = (list[i + j], list[i + j + 1]);
                if a >= 0 && b >= 0 && (a as usize) < n && (b as usize) < n {
                    let (a, b) = (a as usize, b as usize);
                    // Only join parts that are still separate.
                    if root(&parts_snapshot, a) != root(&parts_snapshot, b) {
                        adj[a].push(b as u32);
                        adj[b].push(a as u32);
                    }
                }
            }
            i += nratoms + 1;
        }
    }
    adj
}

fn make_whole_graph(
    x: &mut [Rvec],
    boxm: &Matrix,
    pbc_type: crate::frame::PbcType,
    start: usize,
    adjacency: &[Vec<u32>],
) {
    let natoms = adjacency.len();
    if natoms == 0 {
        return;
    }
    // Mirrors mk_mshift(): every atom accumulates an integer number of box
    // vectors (its "shift") and shift_self() applies them at the end.  Adding
    // an integer multiple of a box vector - including zero - is what GROMACS
    // does, so coordinates that need no shift are still written back
    // unchanged in value (and -0.0 becomes 0.0, exactly as in the original).
    let npbcdim = pbc::npbc_dims(pbc_type);
    let triclinic = crate::frame::is_triclinic(boxm);
    let hbox = [0.5 * boxm[0][0], 0.5 * boxm[1][1], 0.5 * boxm[2][2]];
    let mut visited = vec![false; natoms];
    let mut ishift = vec![[0i32; 3]; natoms];
    // `mk_mshift()` always expands the lowest indexed grey atom next, which
    // matters because an atom's shift is fixed by the atom that discovered it.
    let mut queue = std::collections::BinaryHeap::new();
    for seed in 0..natoms {
        if visited[seed] {
            continue;
        }
        visited[seed] = true;
        queue.push(std::cmp::Reverse(seed));
        while let Some(std::cmp::Reverse(prev)) = queue.pop() {
            for &nb in &adjacency[prev] {
                if nb as usize >= natoms {
                    continue;
                }
                let cur = nb as usize;
                if visited[cur] {
                    continue;
                }
                visited[cur] = true;
                let (gi, gj) = (start + cur, start + prev);
                // mk_1shift / mk_1shift_tric: the neighbour's shift index is
                // the parent's plus at most one step per box vector.
                let mut dx = [
                    x[gj][0] - x[gi][0],
                    x[gj][1] - x[gi][1],
                    x[gj][2] - x[gi][2],
                ];
                let pi = ishift[prev];
                let mut ci = [0i32; 3];
                if triclinic {
                    for m in (0..npbcdim).rev() {
                        if dx[m] < -hbox[m] {
                            ci[m] = pi[m] - 1;
                            for d in (0..m).rev() {
                                dx[d] += boxm[m][d];
                            }
                        } else if dx[m] >= hbox[m] {
                            ci[m] = pi[m] + 1;
                            for d in (0..m).rev() {
                                dx[d] -= boxm[m][d];
                            }
                        } else {
                            ci[m] = pi[m];
                        }
                    }
                } else {
                    for m in 0..npbcdim {
                        ci[m] = if dx[m] < -hbox[m] {
                            pi[m] - 1
                        } else if dx[m] >= hbox[m] {
                            pi[m] + 1
                        } else {
                            pi[m]
                        };
                    }
                }
                ishift[cur] = ci;
                queue.push(std::cmp::Reverse(cur));
            }
        }
    }
    // `shift_self()` applies the accumulated shifts to every atom, including
    // the ones that need no shift.
    for (i, is) in ishift.iter().enumerate() {
        let gi = start + i;
        if triclinic {
            x[gi][0] = x[gi][0]
                + is[0] as f32 * boxm[0][0]
                + is[1] as f32 * boxm[1][0]
                + is[2] as f32 * boxm[2][0];
            x[gi][1] = x[gi][1] + is[1] as f32 * boxm[1][1] + is[2] as f32 * boxm[2][1];
            x[gi][2] = x[gi][2] + is[2] as f32 * boxm[2][2];
        } else {
            for d in 0..3 {
                x[gi][d] = x[gi][d] + is[d] as f32 * boxm[d][d];
            }
        }
    }
}

/// Builds per-molecule atom lists from the topology (molecule blocks) or from
/// residue indices when only a structure file is available.
fn build_molecules(
    atoms: &Atoms,
    mtop: Option<&tpr::Mtop>,
    per_residue: bool,
) -> Vec<Vec<usize>> {
    if per_residue {
        let mut out: Vec<Vec<usize>> = vec![Vec::new(); atoms.nres()];
        for (i, a) in atoms.atom.iter().enumerate() {
            let r = a.resind as usize;
            if r < out.len() {
                out[r].push(i);
            }
        }
        out.retain(|m| !m.is_empty());
        return out;
    }
    if let Some(mtop) = mtop {
        let mut out = Vec::new();
        let mut offset = 0usize;
        for mb in &mtop.molblocks {
            let n = mtop
                .moltypes
                .get(mb.moltype_index as usize)
                .map(|m| m.atoms.nr())
                .unwrap_or(0);
            for _ in 0..mb.nmol.max(0) {
                out.push((offset..offset + n).collect());
                offset += n;
            }
        }
        if !out.is_empty() {
            return out;
        }
    }
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); atoms.nres()];
    for (i, a) in atoms.atom.iter().enumerate() {
        let r = a.resind as usize;
        if r < out.len() {
            out[r].push(i);
        }
    }
    out.retain(|m| !m.is_empty());
    out
}

pub fn run(argv: Vec<String>) -> i32 {
    // Optional phase timing, enabled with GMXRS_TIME=1.
    let b_time = std::env::var_os("GMXRS_TIME").is_some();
    let t_start = std::time::Instant::now();
    let mut t_tpr = std::time::Duration::ZERO;
    let mut t_read = std::time::Duration::ZERO;
    let mut t_process = std::time::Duration::ZERO;
    let mut t_write = std::time::Duration::ZERO;
    let mut n_frames = 0usize;
    let args = Args::parse(argv);
    let in_file = match args.get("f") {
        Some(f) => f,
        None => {
            eprintln!("gmx-rs-tools trjconv: option -f is required");
            return 1;
        }
    };
    let out_file = match args.get("o") {
        Some(o) => o,
        None => {
            eprintln!("gmx-rs-tools trjconv: option -o is required");
            return 1;
        }
    };
    let out_format = match trx::format_from_path(&out_file) {
        Some(f) => f,
        None => {
            eprintln!("Output file name '{out_file}' has an unsupported extension");
            return 1;
        }
    };
    // `fprintf(stderr, "Will write %s: %s\n", ...)`
    eprintln!(
        "Will write {}: {}",
        out_format.extension(),
        out_format.description()
    );

    // Command line options.  These have to be read before the input file is
    // opened because -b/-e select which frames are read at all.
    let skip_nr = args.int("skip", 1).max(1);
    let delta_t = args.real("dt", 0.0);
    let tzero_in = args.real("t0", 0.0);
    let b_set_time = args.has("t0");
    let timestep = args.real("timestep", 0.0);
    let b_timestep = args.has("timestep");
    let tdump = args.real("dump", -1.0);
    let b_tdump = args.has("dump");
    // -b/-e select the time range that is read from the input.
    let b_begin = args.has("b");
    let tbegin = args.real("b", 0.0);
    let b_end = args.has("e");
    let tend = args.real("e", f64::MAX);
    let ndec = args.int("ndec", 3);
    let b_set_prec = args.has("ndec");
    let b_vel = args.flag("vel");
    let b_force = args.flag("force");
    let b_center = args.flag("center");
    let b_sep = args.flag("sep");
    let nzero = args.int("nzero", 0);
    let b_set_box = args.has("box");
    let newbox = args.real3("box", [-1.0, -1.0, -1.0]);
    let b_trans = args.has("trans");
    let trans = args.real3("trans", [0.0, 0.0, 0.0]);
    let b_shift = args.has("shift");
    let shift = args.real3("shift", [0.0, 0.0, 0.0]);
    let b_round = args.flag("round");
    let pbc_mode = args
        .get("pbc")
        .map(|s| PbcMode::parse(&s))
        .unwrap_or(PbcMode::None);
    let ur_mode = match args.get("ur").as_deref() {
        Some("tric") => UrMode::Tric,
        Some("compact") => UrMode::Compact,
        _ => UrMode::Rect,
    };
    let ecenter = ECenter::from_name(&args.get("boxcenter").unwrap_or_else(|| "tric".into()));
    let fit_mode = match args.get("fit").as_deref() {
        Some("rot+trans") => FitMode::Fit,
        Some("rotxy+transxy") => FitMode::FitXy,
        Some("translation") => FitMode::Translation,
        Some("transxy") => FitMode::TransXy,
        Some("progressive") => FitMode::Progressive,
        _ => FitMode::None,
    };
    let b_fit = matches!(
        fit_mode,
        FitMode::Fit | FitMode::FitXy | FitMode::Progressive
    );
    let b_reset = matches!(
        fit_mode,
        FitMode::Fit
            | FitMode::FitXy
            | FitMode::Translation
            | FitMode::TransXy
            | FitMode::Progressive
    );
    let nfitdim = if matches!(fit_mode, FitMode::FitXy | FitMode::TransXy) {
        2
    } else {
        3
    };
    if b_fit && pbc_mode != PbcMode::None {
        eprintln!(
            "PBC condition treatment does not work together with rotational fit.\n\
             Please do the PBC condition treatment first and then run trjconv in a second step\n\
             for the rotational fit."
        );
        return 1;
    }

    // Frames are streamed one at a time so that trajectories larger than
    // memory can be converted.
    let (in_format, mut source) = match trx::FrameSource::open(&in_file) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    // A bounded `-e` gives an exact fraction to report; otherwise the bar
    // follows how much of the input file has been read.
    let mut progress = FrameProgress::new(if b_end { Some((tbegin, tend)) } else { None });
    // Frames before -b are skipped while reading, like read_first_frame().
    let mut first_frame = loop {
        match source.next_frame() {
            Ok(Some(f)) => {
                progress.tick(&mut source, 0, f.time.unwrap_or(0.0));
                if b_begin && f.time.unwrap_or(0.0) < tbegin {
                    continue;
                }
                break f;
            }
            Ok(None) => {
                eprintln!("Could not read a frame from {in_file}");
                return 1;
            }
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        }
    };
    if b_end && !b_begin {
        // Without -b the window starts at the first frame that is read.
        progress.time_window = Some((first_frame.time.unwrap_or(0.0), tend));
    }
    // The run input file is read and the groups are picked next; the bar has
    // to stop drawing before the first of those messages, otherwise they would
    // be printed on top of it.  It starts again for the frame loop.
    progress.suspend();
    // Topology: prefer the run input file.
    let mut mtop: Option<tpr::Mtop> = None;
    let mut atoms: Option<Atoms> = None;
    let mut base_title = first_frame.title.clone();
    // Reference coordinates for -fit / -pbc nojump: the structure file
    // coordinates (xp from read_tps_conf()).
    let mut structure_x: Option<Vec<Rvec>> = None;
    let mut structure_boxm: Option<Matrix> = None;
    // The pbc type of the run input file applies to every frame.
    let mut pbc_override: Option<crate::frame::PbcType> = None;
    // GROMACS falls back to "topol.tpr" for -s whenever the topology is
    // needed, and fails if that file cannot be read.
    let top_file = args.get("s").unwrap_or_else(|| "topol.tpr".to_string());
    let b_tps = args.has("s")
        || b_fit
        || b_reset
        || matches!(
            pbc_mode,
            PbcMode::Whole | PbcMode::Mol | PbcMode::Res | PbcMode::Cluster
        )
        || matches!(out_format, TrxFormat::Gro | TrxFormat::Pdb);
    if b_tps {
        let t0 = std::time::Instant::now();
        match tpr::TprFile::read(&top_file) {
            Ok(t) => match tpr::parse_body(&t.header, &t.body) {
                Ok(body) => {
                    if let Some(m) = body.mtop {
                        base_title = m.name.clone();
                        let mut global = m.global_atoms();
                        // read_tps_conf() adds chain identifiers to the
                        // topology, which is why PDB output from a tpr has
                        // chain letters.
                        tpr::assign_chain_ids(&mut global, &m.molecule_ranges());
                        atoms = Some(global);
                        mtop = Some(m);
                    }
                    structure_x = body.x.clone();
                    structure_boxm = body.boxm;
                    let pbc = body
                        .ir
                        .as_ref()
                        .map(|ir| ir.pbc())
                        .unwrap_or(first_frame.pbc_type);
                    first_frame.pbc_type = pbc;
                    pbc_override = Some(pbc);
                }
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            },
            Err(e) => {
                eprintln!("File input/output error:\n{top_file}\n({e})");
                return 1;
            }
        }
        t_tpr = t0.elapsed();
    } else if let Some(a) = first_frame.atoms.clone() {
        atoms = Some(a);
    }
    let atoms = match atoms {
        Some(a) => a,
        None => Atoms::default(),
    };
    // `-pbc mol` and `-pbc cluster` need the molecule descriptions, which only
    // a run input file provides ("There are no molecule descriptions. I need a
    // .tpr file for this pbc option.").
    if mtop.is_none() && matches!(pbc_mode, PbcMode::Mol | PbcMode::Cluster) {
        let name = if pbc_mode == PbcMode::Mol { "mol" } else { "cluster" };
        eprintln!("Option -pbc {name} requires a .tpr file for the -s option");
        return 1;
    }

    let mut prec = 1.0f32;
    for _ in 0..ndec {
        prec *= 10.0;
    }

    // Index groups.
    let index_file = args.get("n");
    let groups: Vec<IndexGroup> = match &index_file {
        Some(f) => match index::read_ndx(f) {
            Ok(g) => g,
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        },
        None => {
            if atoms.nr() == 0 {
                eprintln!("No index file specified and no topology available for default groups");
                return 1;
            }
            index::analyse(&atoms, false)
        }
    };

    // Fit group.
    // The bar is erased while the groups are being picked, so that its line
    // does not mix with the prompts.
    progress.suspend();
    let mut fit_index: Vec<usize> = Vec::new();
    if b_reset {
        println!("Select group for {} fit", if b_fit { "least squares" } else { "translational" });
        let g = match require_group(&groups, "Select a group:") {
            Some(g) => g,
            None => return 1,
        };
        fit_index = groups[g].particle_indices.clone();
        if b_fit && fit_index.len() < 2 {
            eprintln!("Need at least 2 atoms to fit!");
            return 1;
        }
    } else if pbc_mode == PbcMode::Cluster {
        println!("Select group for clustering");
        let g = match require_group(&groups, "Select a group:") {
            Some(g) => g,
            None => return 1,
        };
        fit_index = groups[g].particle_indices.clone();
    }
    // Centering group.
    let mut center_index: Vec<usize> = Vec::new();
    if b_center {
        println!("Select group for centering");
        let g = match require_group(&groups, "Select a group:") {
            Some(g) => g,
            None => return 1,
        };
        center_index = groups[g].particle_indices.clone();
    }
    println!("Select group for output");
    let gout = match require_group(&groups, "Select a group:") {
        Some(g) => g,
        None => return 1,
    };
    let out_index: Vec<usize> = groups[gout].particle_indices.clone();

    for &i in &out_index {
        if i >= first_frame.natoms {
            eprintln!(
                "Index[{i}] {} is larger than the number of atoms in the trajectory file ({}).",
                i + 1,
                first_frame.natoms
            );
            return 1;
        }
    }

    let mut w_rls = vec![0.0f64; atoms.nr().max(first_frame.natoms)];
    if b_reset {
        for &i in &fit_index {
            if i < atoms.nr() {
                w_rls[i] = atoms.atom[i].mass;
            }
        }
    }

    // Molecule instances (start atom + molecule type) mirror `top->mols`
    // together with `idef` and drive the "whole molecule" treatment.
    let mut mol_instances: Vec<(usize, usize)> = Vec::new();
    if let Some(m) = &mtop {
        let mut offset = 0usize;
        for mb in &m.molblocks {
            for _ in 0..mb.nmol.max(0) {
                mol_instances.push((offset, mb.moltype_index as usize));
                offset += m
                    .moltypes
                    .get(mb.moltype_index as usize)
                    .map(|mt| mt.atoms.nr())
                    .unwrap_or(0);
            }
        }
    }
    let residue_groups = build_molecules(&atoms, None, true);
    let molecule_atom_lists: Vec<Vec<usize>> = mol_instances
        .iter()
        .map(|(start, mt)| {
            let n = mtop
                .as_ref()
                .and_then(|m| m.moltypes.get(*mt))
                .map(|mt| mt.atoms.nr())
                .unwrap_or(0);
            (*start..*start + n).collect()
        })
        .collect();
    let masses: Vec<f64> = (0..atoms.nr()).map(|i| atoms.atom[i].mass).collect();

    // Bonded interaction graph per molecule type (mk_graph_moltype()).  Atoms
    // that are not part of any bonded interaction are chained to their
    // predecessor so that they still end up next to the rest of the molecule.
    let mol_adjacency: Vec<Vec<Vec<u32>>> = match &mtop {
        Some(m) => m.moltypes.iter().map(build_moltype_adjacency).collect(),
        None => Vec::new(),
    };

    // Makes a set of coordinates whole, following the bond graph of the
    // molecule types (the role of `gmx_rmpbc_apply()`).
    let unwrap = |x: &mut [Rvec], boxm: &Matrix, pbc_type: crate::frame::PbcType| {
        if let Some(m) = &mtop {
            for (start, mt) in &mol_instances {
                if let Some(adj) = mol_adjacency.get(*mt) {
                    let _ = m;
                    make_whole_graph(x, boxm, pbc_type, *start, adj);
                }
            }
        } else {
            // Without a topology there is no bond graph: GROMACS warns and
            // leaves the coordinates untouched.
            warn_no_topology();
        }
    };

    // Reference structure for fitting.
    let mut xp: Vec<Rvec> = structure_x
        .clone()
        .or_else(|| first_frame.x.clone())
        .unwrap_or_default();
    let mut x_shift = [0.0f32; 3];
    if b_reset {
        // bRmPBC is set for -fit as well, and the reference box is the one
        // from the structure file.
        let b_rmpbc = b_fit
            || matches!(pbc_mode, PbcMode::Whole | PbcMode::Mol | PbcMode::Res);
        if b_rmpbc {
            let box0 = structure_boxm
                .or(first_frame.boxm)
                .unwrap_or([[0.0; 3]; 3]);
            let mut xp2 = xp.clone();
            let pbc0 = structure_x
                .as_ref()
                .map(|_| first_frame.pbc_type)
                .unwrap_or(first_frame.pbc_type);
            unwrap(&mut xp2, &box0, pbc0);
            xp = xp2;
        }
        // x_shift records how far the reference structure was moved so that
        // the fitted coordinates can be translated back.
        if !out_index.is_empty() && out_index[0] < xp.len() {
            x_shift = xp[out_index[0]];
        }
        let n = atoms.nr().max(xp.len());
        pbc::reset_x_ndim(nfitdim, fit_index.len(), &fit_index, n, None, &mut xp, &w_rls);
        if !out_index.is_empty() && out_index[0] < xp.len() {
            for d in 0..3 {
                x_shift[d] -= xp[out_index[0]][d];
            }
        }
    }

    let mut out_index_for_writer = out_index.clone();
    if out_index_for_writer.is_empty() {
        out_index_for_writer = (0..first_frame.natoms).collect();
    }

    let mut writer = match TrxWriter::create(&out_file, out_format, prec) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    if !atoms.atom.is_empty() {
        writer.set_atoms(&atoms, &out_index_for_writer);
    }

    let mut tshift = 0.0f64;
    let mut tzero = 0.0f64;
    let first_time = first_frame.time.unwrap_or(0.0);
    if b_set_time {
        tshift = tzero_in - first_time;
    } else {
        tzero = first_time;
    }

    // -dump: pick the frame nearest to the requested time.  This follows the
    // look-ahead logic of the original: the decision is taken as soon as a
    // frame at or after the dump time is seen, choosing between that frame and
    // the previous one; if the trajectory ends earlier the last frame is used.
    let mut dump_frame: Option<Frame> = None;
    if b_tdump {
        let mut prev: Option<Frame> = None;
        let mut candidate = Some(first_frame.clone());
        loop {
            let fr = match candidate.take() {
                Some(f) => f,
                None => match source.next_frame() {
                    Ok(Some(f)) => f,
                    Ok(None) => break,
                    Err(e) => {
                        eprintln!("{e}");
                        return 1;
                    }
                },
            };
            let t = fr.time.unwrap_or(0.0);
            if b_end && t > tend {
                break;
            }
            // GROMACS compares the times as `real` (float), which matters for
            // the tie break between two equally distant frames.
            let t_f = t as f32;
            let tdump_f = tdump as f32;
            if t_f >= tdump_f {
                let chosen = match &prev {
                    Some(p) if (t_f - tdump_f) > (tdump_f - p.time.unwrap_or(0.0) as f32) => {
                        p.clone()
                    }
                    _ => fr,
                };
                dump_frame = Some(chosen);
                break;
            }
            prev = Some(fr);
        }
        if dump_frame.is_none() {
            // No frame reached the dump time: use the last one read.
            dump_frame = prev;
        }
    }

    let mut outframe = 0usize;
    let mut file_nr = 0usize;
    let mut last_written = (0usize, 0.0f64);
    let mut previous: Option<Vec<Rvec>> = None;

    // Base name for -sep / generated file names.
    let (out_base, out_ext) = if b_sep {
        match out_file.rfind('.') {
            Some(p) => (out_file[..p].to_string(), out_file[p + 1..].to_string()),
            None => {
                eprintln!("Output file name '{out_file}' does not contain a '.'");
                return 1;
            }
        }
    } else {
        (String::new(), String::new())
    };

    // Frame loop.  With -dump only the selected frame is processed.
    let mut frame = 0usize;
    let mut next_in: Option<Frame> = if b_tdump {
        dump_frame.clone()
    } else {
        Some(first_frame.clone())
    };
    while let Some(mut fr) = match next_in.take() {
        Some(f) => Some(f),
        None => {
            let t = std::time::Instant::now();
            let r = match source.next_frame() {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            t_read += t.elapsed();
            r
        }
    } {
        n_frames += 1;
        progress.tick(&mut source, n_frames as u64, fr.time.unwrap_or(0.0));
        let t_process_start = std::time::Instant::now();
        // -e stops reading (and writing) beyond the requested end time.
        if b_end && fr.time.unwrap_or(0.0) > tend {
            break;
        }
        if let Some(pbc) = pbc_override {
            fr.pbc_type = pbc;
        }
        let natoms = fr.natoms;
        if fr.x.is_none() {
            if b_tdump {
                break;
            }
            frame += 1;
            continue;
        }
        if b_set_box {
            let mut boxm = fr.boxm.unwrap_or([[0.0; 3]; 3]);
            for m in 0..3 {
                if newbox[m] >= 0.0 {
                    boxm[m][m] = newbox[m] as f32;
                }
            }
            fr.boxm = Some(boxm);
        }
        if b_trans {
            if let Some(x) = &mut fr.x {
                for xi in x.iter_mut() {
                    for d in 0..3 {
                        xi[d] += trans[d] as f32;
                    }
                }
            }
        }

        let boxm = fr.boxm.unwrap_or([[0.0; 3]; 3]);
        if pbc_mode == PbcMode::NoJump {
            if frame == 0 && structure_x.is_none() {
                previous = fr.x.clone();
            } else {
                // For the first frame the reference is the structure file, for
                // all later frames it is the previous (already unwrapped) one.
                let prev: Vec<Rvec> = if frame == 0 {
                    xp.clone()
                } else {
                    previous.clone().unwrap_or_default()
                };
                if let Some(x) = &mut fr.x {
                    for i in 0..natoms {
                        let dx = pbc::pbc_dx(fr.pbc_type, &boxm, &x[i], &prev[i]);
                        for d in 0..3 {
                            x[i][d] = prev[i][d] + dx[d];
                        }
                    }
                }
                previous = fr.x.clone();
            }
        } else if pbc_mode == PbcMode::Cluster {
            // `-pbc cluster` is applied to every frame, like the C tool.
            // If a structure file is given the first frame is additionally
            // unwrapped against it (bTPS case of the nojump branch).
            if let Some(x) = &mut fr.x {
                pbc::calc_pbc_cluster(
                    ecenter,
                    x,
                    &fit_index,
                    &molecule_atom_lists,
                    fr.pbc_type,
                    &boxm,
                );
            }
        }

        if let Some(x) = &mut fr.x {
            if std::env::var_os("GMXRS_TRACE").is_some() {
                eprintln!(
                    "[trjconv] frame {frame}: pbc={pbc_mode:?} natoms={} mol_instances={} res_groups={} atoms={}",
                    x.len(),
                    mol_instances.len(),
                    residue_groups.len(),
                    atoms.nr()
                );
            }
            // bRmPBC is set for -fit, -pbc whole, -pbc mol and -pbc res.
            let b_rmpbc =
                b_fit || matches!(pbc_mode, PbcMode::Whole | PbcMode::Mol | PbcMode::Res);
            if b_rmpbc {
                unwrap(x, &boxm, fr.pbc_type);
            }
            match pbc_mode {
                PbcMode::Atom => match ur_mode {
                    UrMode::Rect => pbc::put_atoms_in_box(fr.pbc_type, &boxm, x),
                    UrMode::Tric => pbc::put_atoms_in_triclinic_unitcell(ecenter, &boxm, x),
                    UrMode::Compact => {
                        pbc::put_atoms_in_compact_unitcell(fr.pbc_type, ecenter, &boxm, x)
                    }
                },
                PbcMode::Whole => {}
                PbcMode::Mol => {
                    for mol in &molecule_atom_lists {
                        put_group_com_in_box(
                            x, mol, &boxm, fr.pbc_type, ecenter, ur_mode, &masses,
                        );
                    }
                }
                PbcMode::Res => {
                    for res in &residue_groups {
                        put_group_com_in_box(
                            x, res, &boxm, fr.pbc_type, ecenter, ur_mode, &masses,
                        );
                    }
                }
                _ => {}
            }
        }

        if fit_mode == FitMode::Progressive || b_reset {
            if let Some(x) = &mut fr.x {
                let n = x.len();
                pbc::reset_x_ndim(nfitdim, fit_index.len(), &fit_index, n, None, x, &w_rls);
                if b_fit {
                    pbc::do_fit_ndim(nfitdim, n, &w_rls, &xp, x);
                }
                // Progressive fitting uses the fitted coordinates of this
                // frame as the reference for the next one.
                if fit_mode == FitMode::Progressive {
                    xp = x.clone();
                }
                if !b_center {
                    for xi in x.iter_mut() {
                        for d in 0..3 {
                            xi[d] += x_shift[d];
                        }
                    }
                }
            }
        }

        if b_center {
            if let Some(x) = &mut fr.x {
                pbc::center_x(ecenter, x, &boxm, &center_index);
            }
        }

        // Output frame time.
        let mut frout_time = fr.time.unwrap_or(0.0);
        if b_timestep {
            frout_time = tzero + frame as f64 * timestep;
        } else if b_set_time {
            frout_time += tshift;
        }

        let mut b_write = if b_tdump {
            true
        } else {
            frame % skip_nr as usize == 0
        };
        if b_write && delta_t != 0.0 {
            let t = if b_round {
                (frout_time + 0.5).floor()
            } else {
                frout_time
            };
            let (t0, dt) = if b_round {
                ((tzero + 0.5).floor(), (delta_t + 0.5).floor())
            } else {
                (tzero, delta_t)
            };
            b_write = rmod(t, t0, dt);
        }

        if !b_write {
            frame += 1;
            continue;
        }

        if b_tdump {
            // The leading newline of the original terminates its in-line
            // progress output, which this port does not print.
            eprintln!("Dumping frame at t= {frout_time} ps");
        }

        let mut outfr = fr.clone();
        outfr.time = Some(frout_time);
        let has_time = fr.time.is_some();
        let has_step = fr.step.is_some();
        if outfr.step.is_none() {
            outfr.step = Some(frame as i64);
        }
        outfr.natoms = out_index_for_writer.len();
        if let Some(v) = fr.v.clone() {
            outfr.v = if b_vel { Some(v) } else { None };
        } else {
            outfr.v = None;
        }
        if let Some(f) = fr.f.clone() {
            outfr.f = if b_force { Some(f) } else { None };
        } else {
            outfr.f = None;
        }
        if !b_set_prec && fr.prec.is_some() {
            // -ndec not given: keep the input XTC precision.
            // GROMACS copies `fr.prec` (the multiplication factor stored in
            // the compressed coordinate block) straight into the output
            // frame, so a trajectory is re-encoded at the precision it was
            // written with instead of the next power of ten above it.
            let p = fr.prec.unwrap();
            if let TrxWriter::Xtc { prec, .. } = &mut writer {
                *prec = p;
            }
        }
        if b_shift {
            if let Some(x) = &mut outfr.x {
                for i in 0..out_index_for_writer.len() {
                    for d in 0..3 {
                        x[i][d] += (outframe as f64 * shift[d]) as f32;
                    }
                }
            }
        }

        // GROMACS only writes " t= " and " step= " into the title when the
        // input frame actually carried a time / step (frout.bTime/bStep).
        let title = build_title(&base_title, has_time, outfr.step, has_step, frout_time);
        let title = if title.is_empty() {
            "Generated by trjconv".to_string()
        } else {
            title
        };

        let t_write_start = std::time::Instant::now();
        t_process += t_write_start.duration_since(t_process_start);
        let result = if b_sep {
            let name = mk_filenm(&out_base, &out_ext, nzero as usize, file_nr);
            let mut single = match TrxWriter::create(&name, out_format, prec) {
                Ok(w) => w,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            };
            if !atoms.atom.is_empty() {
                single.set_atoms(&atoms, &out_index_for_writer);
            }
            let r = single.write_frame(&outfr, &out_index_for_writer, &title);
            r.and_then(|_| single.finish())
        } else {
            writer.write_frame(&outfr, &out_index_for_writer, &title)
        };
        t_write += t_write_start.elapsed();
        if let Err(e) = result {
            eprintln!("{e}");
            return 1;
        }
        if b_sep {
            file_nr += 1;
        }

        last_written = (outframe, frout_time);
        outframe += 1;
        let _ = in_format;
        frame += 1;
        if b_tdump {
            break;
        }
    }

    // With -sep every frame is written to its own file, so the combined output
    // file is not created at all (matching gmx trjconv).
    if !b_sep {
        if let Err(e) = writer.finish() {
            eprintln!("{e}");
            return 1;
        }
    }
    progress.finish();
    eprintln!(
        "Last written: frame {:6} time {:8.3}",
        last_written.0, last_written.1
    );
    if outframe == 0 {
        eprintln!("WARNING no output, last frame read at t={}", last_written.1);
    }
    if b_time {
        eprintln!(
            "[time] frames={n_frames} tpr={:.3}s read={:.3}s process={:.3}s write={:.3}s total={:.3}s",
            t_tpr.as_secs_f64(),
            t_read.as_secs_f64(),
            t_process.as_secs_f64(),
            t_write.as_secs_f64(),
            t_start.elapsed().as_secs_f64()
        );
    }
    0
}

fn put_group_com_in_box(
    x: &mut [Rvec],
    group: &[usize],
    boxm: &Matrix,
    pbc_type: crate::frame::PbcType,
    ecenter: ECenter,
    ur_mode: UrMode,
    masses: &[f64],
) {
    if group.is_empty() {
        return;
    }
    // Mirrors put_molecule_com_in_box()/put_residue_com_in_box(): the
    // mass-weighted centre of mass is put inside the unit cell and the whole
    // group is translated by the resulting shift.
    // The C code accumulates `com` as `real` (float) and `mtot` as double, and
    // divides with `svmul(1.0 / mtot, com, com)`; the rounding is visible in
    // the last printed decimal of coordinates near zero.
    let mut com = [0.0f32; 3];
    let mut mtot = 0.0f64;
    for &i in group {
        let m = masses.get(i).copied().unwrap_or(1.0) as f32;
        for d in 0..3 {
            com[d] += m * x[i][d];
        }
        mtot += m as f64;
    }
    if mtot == 0.0 {
        return;
    }
    let inv = 1.0 / mtot;
    for d in 0..3 {
        com[d] = (com[d] as f64 * inv) as f32;
    }
    let mut newcom = com;
    match ur_mode {
        UrMode::Rect => pbc::put_atoms_in_box(pbc_type, boxm, std::slice::from_mut(&mut newcom)),
        UrMode::Tric => {
            pbc::put_atoms_in_triclinic_unitcell(ecenter, boxm, std::slice::from_mut(&mut newcom))
        }
        UrMode::Compact => pbc::put_atoms_in_compact_unitcell(
            pbc_type,
            ecenter,
            boxm,
            std::slice::from_mut(&mut newcom),
        ),
    }
    let shift = [
        newcom[0] - com[0],
        newcom[1] - com[1],
        newcom[2] - com[2],
    ];
    if shift[0] != 0.0 || shift[1] != 0.0 || shift[2] != 0.0 {
        for &i in group {
            for d in 0..3 {
                x[i][d] += shift[d];
            }
        }
    }
}


fn mk_filenm(base: &str, ext: &str, ndigit: usize, file_nr: usize) -> String {
    let nbuf = file_nr.to_string();
    let mut out = String::from(base);
    if nbuf.len() < ndigit {
        out.push_str(&"0".repeat(ndigit - nbuf.len()));
    }
    out.push_str(&nbuf);
    out.push('.');
    out.push_str(ext);
    out
}

fn build_title(
    base: &str,
    has_time: bool,
    step: Option<i64>,
    has_step: bool,
    time: f64,
) -> String {
    let mut title = base.to_string();
    if let Some(p) = title.find(" t= ") {
        title.truncate(p);
    }
    if let Some(p) = title.find(" step= ") {
        title.truncate(p);
    }
    let mut s = title;
    if has_time {
        let _ = write!(s, " t= {time:9.5}");
    }
    if has_step {
        if let Some(step) = step {
            let _ = write!(s, " step= {step}");
        }
    }
    s
}
