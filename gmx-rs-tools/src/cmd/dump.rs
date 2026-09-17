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
            " ".repeat(indent + 3),
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
    let mut progress = crate::progress::Progress::new();
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
        progress.update_from(source.read_progress(), nframe as u64 + 1, fr.time.unwrap_or(0.0));
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
    progress.finish();
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
    let mut progress = crate::progress::Progress::new();
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
        progress.update_from(source.read_progress(), nframe as u64 + 1, fr.time.unwrap_or(0.0));
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
    progress.finish();
    0
}


fn group_short_name(i: usize) -> &'static str {
    // `shortName()` from `topology/topology.cpp`.
    crate::tparsenames::GROUP_SHORT_NAMES
        .get(i)
        .copied()
        .unwrap_or("Unknown")
}

fn dump_tpr(path: &str, show_numbers: bool, show_params: bool, original_inputrec: bool) -> i32 {
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
        crate::irdump::print_inputrec(&mut out, ir, original_inputrec);
    }

    let _ = writeln!(out, "header:");
    let _ = writeln!(out, "   bIr    = {}present", if h.b_ir { "" } else { "not " });
    let _ = writeln!(out, "   bBox   = {}present", if h.b_box { "" } else { "not " });
    let _ = writeln!(out, "   bTop   = {}present", if h.b_top { "" } else { "not " });
    let _ = writeln!(out, "   bX     = {}present", if h.b_x { "" } else { "not " });
    let _ = writeln!(out, "   bV     = {}present", if h.b_v { "" } else { "not " });
    let _ = writeln!(out, "   bF     = {}present", if h.b_f { "" } else { "not " });
    let _ = writeln!(out, "   natoms = {}", h.natoms);
    let _ = writeln!(out, "   lambda = {}", fmt_e(h.lambda, 0, 6));
    let _ = writeln!(out, "   buffer size = {}", h.size_of_tpr_body);

    if let Some(mtop) = &body.mtop {
        crate::tpdump::pr_mtop(&mut out, 0, "topology", mtop, show_numbers, show_params);
    }
    // `list_tpr()` prints every state matrix, using "not available" for the
    // ones the file does not contain.
    let zeros = [[0.0f32; 3]; 3];
    match &body.boxm {
        Some(m) => print_rvecs(&mut out, 0, "box", m),
        None => out.push_str("box: not available\n"),
    }
    match &body.box_rel {
        Some(m) => print_rvecs(&mut out, 0, "box_rel", m),
        None => out.push_str("box_rel: not available\n"),
    }
    match &body.boxv {
        Some(m) => print_rvecs(&mut out, 0, "boxv", m),
        None => out.push_str("boxv: not available\n"),
    }
    for name in ["pres_prev", "svir_prev", "fvir_prev"] {
        if h.b_box {
            print_rvecs(&mut out, 0, name, &zeros);
        } else {
            let _ = writeln!(out, "{name}: not available");
        }
    }
    // `nosehoover_xi` is only filled in when the file has no temperature
    // coupling state, in which case `do_tpx_finalize()` creates it from the
    // inputrec.
    let nh = match &body.ir {
        Some(ir) if h.ngtc == 0 => ir.opts.nhchainlength.max(0) * ir.opts.ngtc.max(0),
        _ => 0,
    };
    if nh > 0 {
        let _ = write!(out, "nosehoover_xi:\t");
        for _ in 0..nh {
            let _ = write!(out, "  {:>10}", "0");
        }
        out.push('\n');
    } else {
        out.push_str("nosehoover_xi: not available\n");
    }
    match &body.x {
        Some(x) => print_rvecs(&mut out, 0, "x", x),
        None => out.push_str("x: not available\n"),
    }
    match &body.v {
        Some(v) => print_rvecs(&mut out, 0, "v", v),
        None => out.push_str("v: not available\n"),
    }
    print!("{out}");

    // Group statistics, mirroring the tail of list_tpr().
    if let Some(mtop) = &body.mtop {
        println!("Group statistics");
        for g in 0..mtop.groups.len().min(10) {
            let n = mtop.groups[g].len();
            let mut counts = vec![0usize; n];
            let mut total = 0usize;
            let numbers = mtop.group_numbers.get(g);
            for i in 0..mtop.natoms {
                // `getGroupType()`: an empty array means group 0 for every
                // atom.
                let gi = match numbers {
                    Some(v) if !v.is_empty() => v.get(i).copied().unwrap_or(0) as usize,
                    _ => 0,
                };
                if gi < counts.len() {
                    counts[gi] += 1;
                    total += 1;
                }
            }
            let mut line = String::new();
            let _ = write!(line, "{:<12}: ", group_short_name(g));
            for c in &counts {
                let _ = write!(line, "  {:5}", c);
            }
            let _ = write!(line, "  (total {total} atoms)");
            println!("{line}");
        }
    }
    0
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
        if args.has("om") {
            eprintln!(
                "gmx-rs-tools dump: writing an mdp file with -om is not implemented in this \
                 minimal port"
            );
            return 1;
        }
        if args.flag("sys") {
            eprintln!(
                "gmx-rs-tools dump: -sys (whole system topology instead of per molecule type) \
                 is not implemented in this minimal port"
            );
            return 1;
        }
        return dump_tpr(
            &s,
            args.flag_default_true("nr"),
            args.flag("param"),
            args.flag("orgir"),
        );
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
