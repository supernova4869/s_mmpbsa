//! Printing of a decoded run input file, mirroring `pr_mtop()` and friends
//! (`topology/topology.cpp`, `topology/atoms.cpp`, `topology/idef.cpp`,
//! `topology/forcefieldparameters.cpp`, `utility/txtdump.cpp`).
//!
//! Every `printf` conversion is reproduced by a helper that matches the C
//! behaviour, in particular `%e` (two digit signed exponent) and `%g`.

use std::fmt::Write as _;

use crate::tparsenames::{GROUP_SHORT_NAMES, INTERACTION_NAMES, PARTICLE_TYPE_NAMES};
use crate::tpr::{FfParams, Mtop, PVal};

/// `INDENT` from `txtdump.h`.
pub const INDENT: usize = 3;

/// `%*e`
fn e(v: f64, width: usize, prec: usize) -> String {
    crate::cmd::fmt_e(v, width, prec)
}

/// `%*g`
fn g(v: f64, width: usize, prec: usize) -> String {
    crate::cmd::fmt_g(v, width, prec)
}

fn pad(out: &mut String, n: usize) {
    for _ in 0..n {
        out.push(' ');
    }
}

/// `pr_title()`
pub fn title_of(out: &mut String, indent: usize, s: &str) -> usize {
    pad(out, indent);
    let _ = writeln!(out, "{s}:");
    indent + INDENT
}

/// `pr_title_n()`
pub fn title_n(out: &mut String, indent: usize, s: &str, n: usize) -> usize {
    pad(out, indent);
    let _ = writeln!(out, "{s} ({n}):");
    indent + INDENT
}

/// `pr_int()`
pub fn pr_int(out: &mut String, indent: usize, name: &str, v: i64) {
    pad(out, indent);
    let _ = writeln!(out, "{name:<30} = {v}");
}

/// `pr_str()`
pub fn pr_str(out: &mut String, indent: usize, name: &str, v: &str) {
    pad(out, indent);
    let _ = writeln!(out, "{name:<30} = {v}");
}

/// `pr_real()` / `pr_double()`: `%-30s = %g`
pub fn pr_g(out: &mut String, indent: usize, name: &str, v: f64) {
    pad(out, indent);
    let _ = writeln!(out, "{name:<30} = {}", g(v, 0, 6));
}

/// `pr_rvecs()`: an unallocated array prints "not available".
pub fn pr_rvecs(out: &mut String, indent: usize, name: &str, values: &[[f32; 3]]) {
    if values.is_empty() {
        pad(out, indent);
        let _ = writeln!(out, "{name}: not available");
        return;
    }
    pad(out, indent);
    let _ = writeln!(out, "{name} ({}x3):", values.len());
    let ind = indent + INDENT;
    for (i, v) in values.iter().enumerate() {
        pad(out, ind);
        let _ = write!(out, "{name}[{i:5}]={{");
        for (j, c) in v.iter().enumerate() {
            if j != 0 {
                let _ = write!(out, ", ");
            }
            let _ = write!(out, "{}", e(*c as f64, 12, 5));
        }
        let _ = writeln!(out, "}}");
    }
}

fn interaction_name(i: i32) -> &'static str {
    INTERACTION_NAMES
        .get(i as usize)
        .map(|(short, _)| *short)
        .unwrap_or("")
}

fn interaction_longname(i: usize) -> &'static str {
    INTERACTION_NAMES.get(i).map(|(_, l)| *l).unwrap_or("")
}

/// `pr_iparams()`: the values of one force field parameter set.
pub fn pr_iparams(out: &mut String, ftype: i32, vals: &[PVal]) {
    let r = |i: usize| -> f64 {
        match vals.get(i) {
            Some(PVal::Real(v)) => *v,
            _ => 0.0,
        }
    };
    let iv = |i: usize| -> i64 {
        match vals.get(i) {
            Some(PVal::Int(v)) => *v,
            _ => 0,
        }
    };
    // `name=%15.8e` pairs, comma separated.
    let named = |out: &mut String, names: &[&str]| {
        for (i, n) in names.iter().enumerate() {
            if i > 0 {
                let _ = write!(out, ", ");
            }
            let _ = write!(out, "{n}={}", e(r(i), 15, 8));
        }
        let _ = writeln!(out);
    };
    // `printHarmonicInteraction()`: the A/B parameter pairs of the harmonic
    // potential are printed with `%12.5e` instead of `%15.8e`.
    let harmonic = |out: &mut String, rname: &str, krname: &str| {
        let _ = writeln!(
            out,
            "{rname}A={}, {krname}A={}, {rname}B={}, {krname}B={}",
            e(r(0), 12, 5),
            e(r(1), 12, 5),
            e(r(2), 12, 5),
            e(r(3), 12, 5)
        );
    };
    match ftype {
        // Bonds, GROMOS96Bonds, HarmonicPotential: b0/cb
        0 | 1 | 5 => harmonic(out, "b0", "cb"),
        // Angles, GROMOS96Angles: th/ct
        10 | 11 => harmonic(out, "th", "ct"),
        // ImproperDihedrals: xi/cx
        24 => harmonic(out, "xi", "cx"),
        // MorsePotential
        2 => named(out, &["b0A", "cbA", "betaA", "b0B", "cbB", "betaB"]),
        // CubicBonds
        3 => named(out, &["b0", "kb", "kcub"]),
        // ConnectBonds: `ensureEmptyLine()` on a fresh writer emits nothing
        // at all, so the next entry continues on the same line.
        4 => {}
        // FENEBonds
        6 => named(out, &["bm", "kb"]),
        // TabulatedBonds, TabulatedBondsNoCoupling, TabulatedAngles, TabulatedDihedrals
        7 | 8 | 18 | 26 => {
            let _ = writeln!(
                out,
                "tab={}, kA={}, kB={}",
                iv(1),
                e(r(0), 15, 8),
                e(r(2), 15, 8)
            );
        }
        // RestraintBonds
        9 => named(
            out,
            &["lowA", "up1A", "up2A", "kA", "lowB", "up1B", "up2B", "kB"],
        ),
        // RestrictedBendingPotential: costheta0/ktheta
        12 => harmonic(out, "costheta0", "ktheta"),
        // LinearAngles
        13 => named(out, &["klinA", "aA", "klinB", "aB"]),
        // CrossBondBonds
        14 => named(out, &["r1e", "r2e", "krr"]),
        // CrossBondAngles
        15 => named(out, &["r1e", "r2e", "r3e", "krt"]),
        // UreyBradleyPotential
        16 => named(
            out,
            &[
                "thetaA", "kthetaA", "r13A", "kUBA", "thetaB", "kthetaB", "r13B", "kUBB",
            ],
        ),
        // QuarticAngles
        17 => {
            let _ = write!(out, "theta={}", e(r(0), 15, 8));
            for i in 0..5 {
                let _ = write!(out, ", c{i}={}", e(r(1 + i), 15, 8));
            }
            let _ = writeln!(out);
        }
        // ProperDihedrals, PeriodicImproperDihedrals, AngleRestraints,
        // AngleZAxisRestraints
        19 | 25 | 58 | 59 => {
            let _ = writeln!(
                out,
                "phiA={}, cpA={}, phiB={}, cpB={}, mult={}",
                e(r(0), 15, 8),
                e(r(1), 15, 8),
                e(r(2), 15, 8),
                e(r(3), 15, 8),
                iv(4)
            );
        }
        // RyckaertBellemansDihedrals
        20 => {
            for i in 0..6 {
                let _ = write!(
                    out,
                    "{}rbcA[{i}]={}",
                    if i == 0 { "" } else { ", " },
                    e(r(i), 15, 8)
                );
            }
            let _ = writeln!(out);
            for i in 0..6 {
                let _ = write!(
                    out,
                    "{}rbcB[{i}]={}",
                    if i == 0 { "" } else { ", " },
                    e(r(6 + i), 15, 8)
                );
            }
            let _ = writeln!(out);
        }
        // RestrictedTorsionPotential
        21 => named(out, &["phiA", "cpA", "phiB", "cpB"]),
        // CombinedBendingTorsionPotential
        22 => {
            let _ = write!(out, "kphi={}", e(r(0), 15, 8));
            for i in 1..6 {
                let _ = write!(out, ", cbtcA[{}]={}", i - 1, e(r(i), 15, 8));
            }
            let _ = writeln!(out);
        }
        // FourierDihedrals: the OPLS constants are recovered from the
        // Ryckaert-Bellemans coefficients stored in the file.
        23 => {
            let mut va = [0.0f64; 4];
            let mut vb = [0.0f64; 4];
            for (i, v) in [&mut va, &mut vb].into_iter().enumerate() {
                let off = 6 * i;
                v[3] = -0.25 * r(off + 4);
                v[2] = -0.5 * r(off + 3);
                v[1] = 4.0 * v[3] - r(off + 2);
                v[0] = 3.0 * v[2] - 2.0 * r(off + 1);
            }
            for (i, v) in va.iter().enumerate() {
                let _ = write!(
                    out,
                    "{}FourA[{i}]={}",
                    if i == 0 { "" } else { ", " },
                    e(*v, 15, 8)
                );
            }
            let _ = writeln!(out);
            for (i, v) in vb.iter().enumerate() {
                let _ = write!(
                    out,
                    "{}FourB[{i}]={}",
                    if i == 0 { "" } else { ", " },
                    e(*v, 15, 8)
                );
            }
            let _ = writeln!(out);
        }
        // DihedralEnergyCorrectionMap
        27 => {
            let _ = writeln!(out, "cmapA={}, cmapB={}", iv(0), iv(1));
        }
        // GeneralizedBorn*PolarizationUnused: nothing is printed.
        28 | 29 | 30 | 31 | 32 => {}
        // LennardJones14
        33 => named(out, &["c6A", "c12A", "c6B", "c12B"]),
        // LennardJonesCoulomb14Q
        35 => named(out, &["fqq", "qi", "qj", "c6", "c12"]),
        // LennardJonesCoulombNonBondedPairs
        36 => named(out, &["qi", "qj", "c6", "c12"]),
        // LennardJonesShortRange
        37 => named(out, &["c6", "c12"]),
        // BuckinghamShortRange
        38 => named(out, &["a", "b", "c"]),
        // Polarization
        48 => named(out, &["alpha"]),
        // WaterPolarization: the distances use `%9.6f`
        49 => {
            let _ = writeln!(
                out,
                "al_x={}, al_y={}, al_z={}, rOH={}, rHH={}, rOD={}",
                e(r(0), 15, 8),
                e(r(1), 15, 8),
                e(r(2), 15, 8),
                format!("{:9.6}", r(3)),
                format!("{:9.6}", r(4)),
                format!("{:9.6}", r(5))
            );
        }
        // TholePolarization
        50 => {
            let n = if vals.len() >= 4 { 4 } else { 3 };
            let names = ["a", "alpha1", "alpha2", "rfac"];
            named(out, &names[..n]);
        }
        // AnharmonicPolarization: no commas in the format string
        51 => {
            let _ = writeln!(
                out,
                "alpha={} drcut={} khyp={}",
                e(r(0), 15, 8),
                e(r(1), 15, 8),
                e(r(2), 15, 8)
            );
        }
        // PositionRestraints
        52 => {
            let _ = writeln!(
                out,
                "pos0A=({},{},{}), fcA=({},{},{}), pos0B=({},{},{}), fcB=({},{},{})",
                e(r(0), 15, 8),
                e(r(1), 15, 8),
                e(r(2), 15, 8),
                e(r(3), 15, 8),
                e(r(4), 15, 8),
                e(r(5), 15, 8),
                e(r(6), 15, 8),
                e(r(7), 15, 8),
                e(r(8), 15, 8),
                e(r(9), 15, 8),
                e(r(10), 15, 8),
                e(r(11), 15, 8)
            );
        }
        // FlatBottomedPositionRestraints
        53 => {
            let _ = writeln!(
                out,
                "pos0=({},{},{}), geometry={}, r={}, k={}",
                e(r(1), 15, 8),
                e(r(2), 15, 8),
                e(r(3), 15, 8),
                iv(0),
                e(r(4), 15, 8),
                e(r(5), 15, 8)
            );
        }
        // DistanceRestraints
        54 => {
            let _ = writeln!(
                out,
                "label={:>4}, type={:>1}, low={}, up1={}, up2={}, fac={})",
                iv(0),
                iv(1),
                e(r(2), 15, 8),
                e(r(3), 15, 8),
                e(r(4), 15, 8),
                e(r(5), 15, 8)
            );
        }
        // OrientationRestraints
        56 => {
            let _ = writeln!(
                out,
                "ex={:>4}, label={}, power={:>4}, c={}, obs={}, kfac={})",
                iv(0),
                iv(1),
                iv(2),
                e(r(3), 15, 8),
                e(r(4), 15, 8),
                e(r(5), 15, 8)
            );
        }
        // DihedralRestraints
        60 => {
            if vals.len() >= 6 {
                named(out, &["phiA", "dphiA", "kfacA", "phiB", "dphiB", "kfacB"])
            } else {
                let _ = writeln!(
                    out,
                    "phiA={}, dphiA={}, kfacA={}",
                    e(r(0), 15, 8),
                    e(r(1), 15, 8),
                    e(r(2), 15, 8)
                );
            }
        }
        // Constraints, ConstraintsNoCoupling
        62 | 63 => named(out, &["dA", "dB"]),
        // SETTLE
        64 => named(out, &["doh", "dhh"]),
        // VirtualSite1: nothing is printed.
        65 => {}
        // VirtualSite2, VirtualSite2FlexibleDistance
        66 | 67 => named(out, &["a"]),
        // VirtualSite3 variants
        68 | 69 | 70 => named(out, &["a", "b"]),
        // VirtualSite3Outside, VirtualSite4FlexibleDistance(+Normalization)
        71 | 72 | 73 => named(out, &["a", "b", "c"]),
        // VirtualSiteN
        74 => {
            let _ = writeln!(out, "n={:>2}, a={}", iv(0), e(r(1), 15, 8));
        }
        _ => {
            let _ = writeln!(out);
        }
    }
}

/// `pr_ffparams()`
pub fn pr_ffparams(
    out: &mut String,
    indent: usize,
    title: &str,
    ff: &FfParams,
    show_numbers: bool,
    cmap_grid_spacing: i32,
    cmap_data: &[Vec<f32>],
) {
    let ind = title_of(out, indent, title);
    pad(out, ind);
    let _ = writeln!(out, "atnr={}", ff.atnr);
    pad(out, ind);
    let _ = writeln!(out, "ntypes={}", ff.functype.len());
    for (i, ft) in ff.functype.iter().enumerate() {
        pad(out, ind + INDENT);
        let _ = write!(
            out,
            "functype[{}]={}, ",
            if show_numbers { i as i64 } else { -1 },
            interaction_name(*ft)
        );
        pr_iparams(
            out,
            *ft,
            ff.iparams.get(i).map(|v| v.as_slice()).unwrap_or(&[]),
        );
    }
    pr_g(out, ind, "reppow", ff.reppow);
    pr_g(out, ind, "fudgeQQ", ff.fudge_qq);
    pr_cmap(out, title_cmap(), cmap_grid_spacing, cmap_data, show_numbers);
}

fn title_cmap() -> &'static str {
    "cmap"
}

/// `pr_cmap()`: the dihedral energy correction maps.
fn pr_cmap(
    out: &mut String,
    title: &str,
    grid_spacing: i32,
    cmap_data: &[Vec<f32>],
    show_numbers: bool,
) {
    let _ = writeln!(out, "{title}");
    let dx = if grid_spacing != 0 {
        360.0 / grid_spacing as f64
    } else {
        0.0
    };
    let nelem = (grid_spacing * grid_spacing).max(0) as usize;
    for (i, data) in cmap_data.iter().enumerate() {
        let mut idx = -180.0f64;
        let _ = writeln!(out, "{:>8} {:>8} {:>8} {:>8}", "V", "dVdx", "dVdy", "d2dV");
        let _ = writeln!(
            out,
            "grid[{:>3}]={{",
            if show_numbers { i as i64 } else { -1 }
        );
        for j in 0..nelem {
            if grid_spacing != 0 && j % grid_spacing as usize == 0 {
                let _ = writeln!(out, "{idx:8.1}");
                idx += dx;
            }
            let v = |k: usize| data.get(j * 4 + k).copied().unwrap_or(0.0) as f64;
            let _ = writeln!(
                out,
                "{:8.3} {:8.3} {:8.3} {:8.3}",
                v(0),
                v(1),
                v(2),
                v(3)
            );
        }
        out.push('\n');
    }
}

/// `pr_atom()`
fn pr_atom(out: &mut String, indent: usize, atoms: &crate::frame::Atoms) {
    let ind = title_n(out, indent, "atom", atoms.nr());
    for (i, a) in atoms.atom.iter().enumerate() {
        pad(out, ind);
        let ptype = PARTICLE_TYPE_NAMES
            .get(a.ptype as usize)
            .copied()
            .unwrap_or("Atom");
        let _ = write!(out, "atom[{i:6}]={{");
        let _ = write!(out, "type={:3}, typeB={:3}, ", a.type_id, a.type_id_b);
        let _ = write!(out, "ptype={ptype:>8}, ");
        let _ = write!(out, "m={}, ", e(a.mass, 12, 5));
        let _ = write!(out, "q={}, ", e(a.charge, 12, 5));
        let _ = write!(out, "mB={}, ", e(a.mass_b, 12, 5));
        let _ = write!(out, "qB={}, ", e(a.charge_b, 12, 5));
        let _ = write!(out, "resind={:5}, ", a.resind);
        let _ = writeln!(out, "atomnumber={:3}}}", a.atomnumber);
    }
}

/// `pr_strings()` for atom names.
fn pr_strings(out: &mut String, indent: usize, title: &str, names: &[String], show_numbers: bool) {
    let ind = title_n(out, indent, title, names.len());
    for (i, n) in names.iter().enumerate() {
        pad(out, ind);
        let idx = if show_numbers { i as i64 } else { -1 };
        let _ = writeln!(out, "{title}[{idx}]={{name=\"{n}\"}}");
    }
}

/// `pr_strings2()` for atom types (A and B state).
fn pr_strings2(
    out: &mut String,
    indent: usize,
    title: &str,
    names: &[(String, String)],
    show_numbers: bool,
) {
    let ind = title_n(out, indent, title, names.len());
    for (i, (a, b)) in names.iter().enumerate() {
        pad(out, ind);
        let idx = if show_numbers { i as i64 } else { -1 };
        let _ = writeln!(out, "{title}[{idx}]={{name=\"{a}\",nameB=\"{b}\"}}");
    }
}

/// `pr_resinfo()`
fn pr_resinfo(out: &mut String, indent: usize, atoms: &crate::frame::Atoms, show_numbers: bool) {
    let ind = title_n(out, indent, "residue", atoms.resinfo.len());
    for (i, ri) in atoms.resinfo.iter().enumerate() {
        pad(out, ind);
        let idx = if show_numbers { i as i64 } else { -1 };
        let ic = if ri.ic == 0 { ' ' } else { ri.ic as char };
        let _ = writeln!(
            out,
            "residue[{idx}]={{name=\"{}\", nr={}, ic='{ic}'}}",
            ri.name, ri.nr
        );
    }
}

/// `pr_atoms()`
pub fn pr_atoms(
    out: &mut String,
    indent: usize,
    title: &str,
    atoms: &crate::frame::Atoms,
    show_numbers: bool,
) {
    let ind = title_of(out, indent, title);
    pr_atom(out, ind, atoms);
    let names: Vec<String> = atoms.atom.iter().map(|a| a.name.clone()).collect();
    pr_strings(out, ind, "atom", &names, show_numbers);
    let types: Vec<(String, String)> = atoms
        .atom
        .iter()
        .map(|a| (a.atom_type.clone(), a.atom_type_b.clone()))
        .collect();
    pr_strings2(out, ind, "type", &types, show_numbers);
    pr_resinfo(out, ind, atoms, show_numbers);
}

/// `pr_listoflists()`
pub fn pr_listoflists(
    out: &mut String,
    indent: usize,
    title: &str,
    lists: &[Vec<i32>],
    show_numbers: bool,
) {
    let ind = title_of(out, indent, title);
    pad(out, ind);
    let _ = writeln!(out, "numLists={}", lists.len());
    pad(out, ind);
    let total: usize = lists.iter().map(|l| l.len()).sum();
    let _ = writeln!(out, "numElements={total}");
    for (i, list) in lists.iter().enumerate() {
        pad(out, ind);
        let prefix = if list.is_empty() {
            format!("{title}[{i}]={{")
        } else {
            let idx = if show_numbers { i as i64 } else { -1 };
            format!("{title}[{idx}][num={}]={{", list.len())
        };
        let _ = write!(out, "{prefix}");
        // `pr_listoflists()` tracks the current line length and wraps at
        // `USE_WIDTH` (= LINE_WIDTH - RMARGIN = 70) columns.
        let mut size = ind + prefix.len();
        for (j, v) in list.iter().enumerate() {
            if j != 0 {
                let _ = write!(out, ", ");
                size += 2;
            }
            if size > 70 {
                out.push('\n');
                pad(out, ind + INDENT);
                size = ind + INDENT;
            }
            let text = v.to_string();
            size += text.len();
            let _ = write!(out, "{text}");
        }
        let _ = writeln!(out, "}}");
    }
}

/// `pr_ilist()`
pub fn pr_ilist(
    out: &mut String,
    indent: usize,
    iftype: usize,
    functype: &[i32],
    iatoms: &[i32],
    show_numbers: bool,
    show_parameters: bool,
    iparams: &[Vec<PVal>],
) {
    let ind = title_of(out, indent, interaction_longname(iftype));
    pad(out, ind);
    let _ = writeln!(out, "nr: {}", iatoms.len());
    if iatoms.is_empty() {
        return;
    }
    pad(out, ind);
    let _ = writeln!(out, "iatoms:");
    let mut i = 0usize;
    let mut j = 0usize;
    while i < iatoms.len() {
        pad(out, ind + INDENT);
        let t = iatoms[i];
        let ftype = functype.get(t as usize).copied().unwrap_or(0);
        if show_numbers {
            let _ = write!(out, "{j} type={t} ");
        }
        j += 1;
        let _ = write!(out, "({})", interaction_name(ftype));
        let nratoms = crate::tpr::interaction_function_nratoms_pub(ftype as usize);
        for k in 0..nratoms {
            if let Some(a) = iatoms.get(i + 1 + k) {
                let _ = write!(out, " {a:3}");
            }
        }
        if show_parameters {
            let _ = write!(out, "  ");
            pr_iparams(
                out,
                ftype,
                iparams.get(t as usize).map(|v| v.as_slice()).unwrap_or(&[]),
            );
        }
        // `pr_iparams()` already ends with a newline, so the `fprintf(fp,
        // "\n")` of `printIlist()` leaves an empty line after each listed
        // interaction when the parameters are shown.
        let _ = writeln!(out);
        i += 1 + nratoms;
    }
}

/// `pr_groups()`
pub fn pr_groups(out: &mut String, indent: usize, mtop: &Mtop, show_numbers: bool) {
    for (i, group) in mtop.groups.iter().enumerate() {
        let short = GROUP_SHORT_NAMES.get(i).copied().unwrap_or("");
        let _ = write!(out, "grp[{short:<12}] nr={}, name=[", group.len());
        for entry in group {
            let name = mtop
                .group_names
                .get(*entry as usize)
                .cloned()
                .unwrap_or_default();
            let _ = write!(out, " {name}");
        }
        let _ = writeln!(out, "]");
    }
    pr_strings(out, indent, "grpname", &mtop.group_names, show_numbers);
    pad(out, indent);
    let _ = write!(out, "groups          ");
    for n in GROUP_SHORT_NAMES.iter() {
        let _ = write!(out, " {n:>5.5}");
    }
    let _ = writeln!(out);
    pad(out, indent);
    let _ = write!(out, "allocated       ");
    let mut nat_max = 0usize;
    for gn in mtop.group_numbers.iter() {
        let _ = write!(out, " {:5}", gn.len());
        nat_max = nat_max.max(gn.len());
    }
    let _ = writeln!(out);
    if nat_max == 0 {
        pad(out, indent);
        let _ = write!(out, "groupnr[{:>5}] =", "*");
        for _ in mtop.group_numbers.iter() {
            let _ = write!(out, "  {:3} ", 0);
        }
        let _ = writeln!(out);
    } else {
        for i in 0..nat_max {
            pad(out, indent);
            let _ = write!(out, "groupnr[{i:5}] =");
            for gn in mtop.group_numbers.iter() {
                let v = gn.get(i).copied().unwrap_or(0);
                let _ = write!(out, "  {v:3} ");
            }
            let _ = writeln!(out);
        }
    }
}

/// `pr_molblock()`
fn pr_molblock(out: &mut String, indent: usize, title: &str, mtop: &Mtop, n: usize) {
    let mb = &mtop.molblocks[n];
    let ind = title_n(out, indent, title, n);
    pad(out, ind);
    let name = mtop
        .moltypes
        .get(mb.moltype_index as usize)
        .map(|m| m.name.clone())
        .unwrap_or_default();
    let _ = writeln!(out, "{:<20} = {} \"{name}\"", "moltype", mb.moltype_index);
    pr_int(out, ind, "#molecules", mb.nmol as i64);
    pr_int(out, ind, "#posres_xA", mb.posres_xa.len() as i64);
    if !mb.posres_xa.is_empty() {
        pr_rvecs(out, ind, "posres_xA", &mb.posres_xa);
    }
    pr_int(out, ind, "#posres_xB", mb.posres_xb.len() as i64);
    if !mb.posres_xb.is_empty() {
        pr_rvecs(out, ind, "posres_xB", &mb.posres_xb);
    }
}

/// `pr_moltype()`
fn pr_moltype(
    out: &mut String,
    indent: usize,
    mtop: &Mtop,
    n: usize,
    show_numbers: bool,
    show_parameters: bool,
) {
    let mt = &mtop.moltypes[n];
    let ind = title_n(out, indent, "moltype", n);
    pad(out, ind);
    let _ = writeln!(out, "name=\"{}\"", mt.name);
    pr_atoms(out, ind, "atoms", &mt.atoms, show_numbers);
    pr_listoflists(out, ind, "excls", &mt.excls, show_numbers);
    for (j, list) in mt.ilists.iter().enumerate() {
        pr_ilist(
            out,
            ind,
            j,
            &mtop.ffparams.functype,
            list,
            show_numbers,
            show_parameters,
            &mtop.ffparams.iparams,
        );
    }
}

/// `pr_mtop()`
pub fn pr_mtop(
    out: &mut String,
    indent: usize,
    t: &str,
    mtop: &Mtop,
    show_numbers: bool,
    show_parameters: bool,
) {
    let ind = title_of(out, indent, t);
    pad(out, ind);
    let _ = writeln!(out, "name=\"{}\"", mtop.name);
    pr_int(out, ind, "#atoms", mtop.natoms as i64);
    pr_int(out, ind, "#molblock", mtop.molblocks.len() as i64);
    for i in 0..mtop.molblocks.len() {
        pr_molblock(out, ind, "molblock", mtop, i);
    }
    pr_str(
        out,
        ind,
        "bIntermolecularInteractions",
        if mtop.b_intermolecular { "true" } else { "false" },
    );
    if mtop.b_intermolecular {
        for (j, list) in mtop.intermolecular_ilists.iter().enumerate() {
            pr_ilist(
                out,
                ind,
                j,
                &mtop.ffparams.functype,
                list,
                show_numbers,
                show_parameters,
                &mtop.ffparams.iparams,
            );
        }
    }
    pr_ffparams(
        out,
        ind,
        "ffparams",
        &mtop.ffparams,
        show_numbers,
        mtop.cmap_grid_spacing,
        &mtop.cmap_data,
    );
    for i in 0..mtop.moltypes.len() {
        pr_moltype(out, ind, mtop, i, show_numbers, show_parameters);
    }
    pr_groups(out, ind, mtop, show_numbers);
}
