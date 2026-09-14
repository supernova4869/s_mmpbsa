//! `gmx dump`, mirroring `src/gromacs/tools/dump.cpp`.

use std::fmt::Write as _;

use crate::cmd::{fmt_e, fmt_g, Args};
use crate::frame::Frame;
use crate::index::strcasecmp;
use crate::tpr;
use crate::trx;

const INDENT: usize = 3;

fn print_rvecs(out: &mut String, indent: usize, title: &str, v: &[[f32; 3]]) {
    let _ = writeln!(out, "{}{} ({}x3):", " ".repeat(indent), title, v.len());
    for (i, row) in v.iter().enumerate() {
        let _ = write!(
            out,
            "{}{}[{:5}]={{",
            " ".repeat(indent + INDENT),
            title,
            i
        );
        for (j, c) in row.iter().enumerate() {
            if j != 0 {
                let _ = write!(out, ", ");
            }
            let _ = write!(out, "{}", fmt_e(*c as f64, 12, 5));
        }
        let _ = writeln!(out, "}}");
    }
}

fn dump_xtc(path: &str) -> i32 {
    // Frames are streamed so that trajectories larger than memory can be
    // dumped.
    let (_, mut source) = match trx::FrameSource::open(path) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let mut nframe = 0usize;
    loop {
        let fr = match source.next_frame() {
            Ok(Some(f)) => f,
            Ok(None) => break,
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        };
        let mut out = String::new();
        let _ = writeln!(out, "{path} frame {nframe}:");
        let _ = writeln!(
            out,
            "{:3}natoms={:10}  step={:10}  time={}  prec={}",
            "",
            fr.natoms,
            fr.step.unwrap_or(0),
            fmt_e(fr.time.unwrap_or(0.0), 12, 7),
            fmt_g(fr.prec.unwrap_or(0.0) as f64, 10, 6)
        );
        if let Some(boxm) = &fr.boxm {
            print_rvecs(&mut out, INDENT, "box", boxm);
        }
        if let Some(x) = &fr.x {
            print_rvecs(&mut out, INDENT, "x", x);
        }
        if let Some(v) = &fr.v {
            print_rvecs(&mut out, INDENT, "v", v);
        }
        if let Some(f) = &fr.f {
            print_rvecs(&mut out, INDENT, "f", f);
        }
        print!("{out}");
        use std::io::Write;
        let _ = std::io::stdout().flush();
        nframe += 1;
    }
    0
}

fn dump_trr(path: &str) -> i32 {
    let (_, mut source) = match trx::FrameSource::open(path) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let mut nframe = 0usize;
    loop {
        let fr = match source.next_frame() {
            Ok(Some(f)) => f,
            Ok(None) => break,
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        };
        let mut out = String::new();
        let _ = writeln!(out, "{path} frame {nframe}:");
        let _ = writeln!(
            out,
            "{:3}natoms={:10}  step={:10}  time={}  lambda={}",
            "",
            fr.natoms,
            fr.step.unwrap_or(0),
            fmt_e(fr.time.unwrap_or(0.0), 12, 7),
            fmt_g(fr.lambda.unwrap_or(0.0) as f64, 10, 6)
        );
        if let Some(boxm) = &fr.boxm {
            print_rvecs(&mut out, INDENT, "box", boxm);
        }
        if let Some(x) = &fr.x {
            print_rvecs(&mut out, INDENT, "x", x);
        }
        if let Some(v) = &fr.v {
            print_rvecs(&mut out, INDENT, "v", v);
        }
        if let Some(f) = &fr.f {
            print_rvecs(&mut out, INDENT, "f", f);
        }
        print!("{out}");
        use std::io::Write;
        let _ = std::io::stdout().flush();
        nframe += 1;
    }
    0
}

/// Prints a compact summary of the tpr topology (the C tool prints every
/// topology substructure with `pr_mtop`; here only the essentials are shown).
fn print_mtop(out: &mut String, mtop: &tpr::Mtop, indent: usize) {
    let pad = " ".repeat(indent);
    let _ = writeln!(out, "{pad}natoms = {}", mtop.natoms);
    let _ = writeln!(out, "{pad}molecule types ({})", mtop.moltypes.len());
    for (i, mt) in mtop.moltypes.iter().enumerate() {
        let _ = writeln!(
            out,
            "{pad}  moltype[{}] '{}' with {} atoms",
            i,
            mt.name,
            mt.atoms.nr()
        );
    }
    let _ = writeln!(out, "{pad}molecule blocks ({})", mtop.molblocks.len());
    for (i, mb) in mtop.molblocks.iter().enumerate() {
        let _ = writeln!(
            out,
            "{pad}  molblock[{}]: type {} with {} molecules",
            i, mb.moltype_index, mb.nmol
        );
    }
    let _ = writeln!(out, "{pad}groups ({})", mtop.group_names.len());
    for (i, g) in mtop.groups.iter().enumerate() {
        if !g.is_empty() {
            let _ = writeln!(out, "{pad}  group[{}] has {} entries", i, g.len());
        }
    }
}

fn group_short_name(i: usize) -> &'static str {
    // Mirrors the SimulationAtomGroupType names used by "Group statistics".
    const NAMES: [&str; 10] = [
        "Temp.",
        "Energy",
        "Acceleration",
        "Freeze",
        "User1",
        "User2",
        "VCM",
        "Compressed",
        "Or. res. fit",
        "QMMM",
    ];
    NAMES.get(i).copied().unwrap_or("Unknown")
}

fn dump_tpr(path: &str, show_numbers: bool, _show_params: bool) -> i32 {
    let tpr = match tpr::TprFile::read(path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let body = match tpr::parse_body(&tpr.header, &tpr.body) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };

    let h = &tpr.header;
    let mut out = String::new();
    let _ = writeln!(out, "{path}:");
    if let Some(ir) = &body.ir {
        let _ = writeln!(out, "inputrec:");
        let fields: [(&str, String); 24] = [
            ("integrator", format!("{}", integrator_name(ir.integrator))),
            ("tinit", fmt_g(ir.init_t, 0, 6)),
            ("dt", fmt_g(ir.delta_t, 0, 6)),
            ("nsteps", format!("{}", ir.nsteps)),
            ("init-step", format!("{}", ir.init_step)),
            ("simulation-part", format!("{}", ir.simulation_part)),
            ("mts", format!("{}", ir.use_mts)),
            ("mass-repartition-factor", fmt_g(ir.mass_repartition_factor, 0, 6)),
            ("comm-mode", format!("{}", comm_mode_name(ir.comm_mode))),
            ("nstcomm", format!("{}", ir.nstcomm)),
            ("rtpi", fmt_g(ir.rtpi, 0, 6)),
            ("nstxout", format!("{}", ir.nstxout)),
            ("nstvout", format!("{}", ir.nstvout)),
            ("nstfout", format!("{}", ir.nstfout)),
            ("nstlog", format!("{}", ir.nstlog)),
            ("nstcalcenergy", format!("{}", ir.nstcalcenergy)),
            ("nstenergy", format!("{}", ir.nstenergy)),
            ("nstxout-compressed", format!("{}", ir.nstxout_compressed)),
            ("compressed-x-precision", fmt_g(ir.x_compression_precision, 0, 6)),
            ("cutoff-scheme", format!("{}", cutoff_scheme_name(ir.cutoff_scheme))),
            ("nstlist", format!("{}", ir.nstlist)),
            ("pbc", format!("{}", ir.pbc().name())),
            ("rlist", fmt_g(ir.rlist, 0, 6)),
            ("rcoulomb", fmt_g(ir.rcoulomb, 0, 6)),
        ];
        for (name, value) in fields {
            let _ = writeln!(out, "   {:<30} = {}", name, value);
        }
    }

    let _ = writeln!(out, "header:");
    let _ = writeln!(out, "   bIr    = {}present", if h.b_ir { "" } else { "not " });
    let _ = writeln!(out, "   bBox   = {}present", if h.b_box { "" } else { "not " });
    let _ = writeln!(out, "   bTop   = {}present", if h.b_top { "" } else { "not " });
    let _ = writeln!(out, "   bX     = {}present", if h.b_x { "" } else { "not " });
    let _ = writeln!(out, "   bV     = {}present", if h.b_v { "" } else { "not " });
    let _ = writeln!(out, "   bF     = {}present", if h.b_f { "" } else { "not " });
    let _ = writeln!(out, "   natoms = {}", h.natoms);
    let _ = writeln!(out, "   lambda = {:e}", h.lambda);
    let _ = writeln!(out, "   buffer size = {}", h.size_of_tpr_body);

    if let Some(mtop) = &body.mtop {
        let _ = writeln!(out, "topology:");
        print_mtop(&mut out, mtop, INDENT);
    }
    if let Some(boxm) = &body.boxm {
        print_rvecs(&mut out, INDENT, "box", boxm);
    }
    if let Some(x) = &body.x {
        print_rvecs(&mut out, INDENT, "x", x);
    }
    if let Some(v) = &body.v {
        print_rvecs(&mut out, INDENT, "v", v);
    }
    let _ = show_numbers;
    print!("{out}");

    // Group statistics, mirroring the tail of list_tpr().
    if let Some(mtop) = &body.mtop {
        println!("Group statistics");
        for g in 0..mtop.groups.len().min(10) {
            let n = mtop.groups[g].len();
            if n == 0 {
                continue;
            }
            let mut counts = vec![0usize; n];
            let mut total = 0usize;
            for &gi in &mtop.group_numbers.get(g).cloned().unwrap_or_default() {
                if (gi as usize) < counts.len() {
                    counts[gi as usize] += 1;
                    total += 1;
                }
            }
            let mut line = String::new();
            let _ = write!(line, "{:<12}: ", group_short_name(g));
            if total > 0 {
                for c in &counts {
                    let _ = write!(line, "  {:5}", c);
                }
            }
            let _ = write!(line, "  (total {total} atoms)");
            println!("{line}");
        }
    }
    0
}

fn integrator_name(v: i32) -> &'static str {
    // IntegrationAlgorithm, see gromacs/mdtypes/md_enums.cpp
    const NAMES: [&str; 13] = [
        "md",
        "steep",
        "cg",
        "bd",
        "sd2 - removed",
        "nm",
        "l-bfgs",
        "tpi",
        "tpic",
        "sd",
        "md-vv",
        "md-vv-avek",
        "mimic",
    ];
    NAMES.get(v as usize).copied().unwrap_or("unknown")
}

fn comm_mode_name(v: i32) -> &'static str {
    match v {
        0 => "Linear",
        1 => "Angular",
        2 => "No",
        _ => "unknown",
    }
}

fn cutoff_scheme_name(v: i32) -> &'static str {
    // CutoffScheme: Verlet = 0, Group = 1
    match v {
        0 => "Verlet",
        1 => "Group",
        _ => "unknown",
    }
}

/// Dumps a topology file by echoing it with the preprocessor applied.  The C
/// tool runs the full C preprocessor; here only `#include` free files are
/// echoed verbatim.
fn dump_top(path: &str) -> i32 {
    match std::fs::read_to_string(path) {
        Ok(content) => {
            print!("{content}");
            if !content.ends_with('\n') {
                println!();
            }
            0
        }
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            1
        }
    }
}

/// Dumps a sparse matrix (.mtx) file, mirroring `list_mtx()`.
fn dump_mtx(path: &str) -> i32 {
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            return 1;
        }
    };
    let numbers: Vec<f64> = content
        .split_whitespace()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect();
    if numbers.len() < 2 {
        eprintln!("invalid matrix file {path}");
        return 1;
    }
    let nrow = numbers[0] as usize;
    let ncol = numbers[1] as usize;
    if numbers.len() < 2 + nrow * ncol {
        eprintln!("matrix file {path} is truncated");
        return 1;
    }
    let mut full = vec![0.0f64; nrow * ncol];
    for i in 0..nrow {
        for j in 0..ncol {
            full[i * ncol + j] = numbers[2 + i * ncol + j];
        }
    }
    println!("{nrow} {ncol}");
    for i in 0..nrow {
        let mut line = String::new();
        for j in 0..ncol {
            let _ = write!(line, " {}", fmt_g(full[i * ncol + j], 0, 6));
        }
        println!("{line}");
    }
    0
}

pub fn run(argv: Vec<String>) -> i32 {
    let args = Args::parse(argv);
    let _ = strcasecmp("a", "a");
    if let Some(f) = args.get("f") {
        return match trx::format_from_path(&f) {
            Some(trx::TrxFormat::Xtc) => dump_xtc(&f),
            Some(trx::TrxFormat::Trr) => dump_trr(&f),
            _ => {
                eprintln!(
                    "File {f} is of an unsupported type. Try using the command\n 'less {f}'\n"
                );
                1
            }
        };
    }
    if let Some(s) = args.get("s") {
        return dump_tpr(&s, args.flag("nr"), args.flag("param"));
    }
    if let Some(p) = args.get("p") {
        return dump_top(&p);
    }
    if let Some(m) = args.get("mtx") {
        return dump_mtx(&m);
    }
    if args.has("e") || args.has("cp") {
        eprintln!("energy and checkpoint dumps are not implemented in this minimal port");
        return 1;
    }
    eprintln!("gmx-rs-tools dump: no input file given (use -s, -f, -p or -mtx)");
    1
}

/// Convenience for tests: dump a set of frames into a string.
pub fn dump_frames_string(frames: &[Frame]) -> String {
    let mut out = String::new();
    for fr in frames {
        if let Some(x) = &fr.x {
            print_rvecs(&mut out, INDENT, "x", x);
        }
    }
    out
}
