//! The GROMACS tools s_mmpbsa needs, run in-process.
//!
//! `s_mmpbsa` used to shell out to `gmx dump`, `gmx trjconv`,
//! `gmx convert-tpr` and `gmx make_ndx`.  Those four tools are now provided by
//! the vendored `gmx-rs-tools` crate (the Rust port in `gmx-rs-tools/`), which
//! is linked into the s_mmpbsa binary, so no external GROMACS program is
//! needed.
//!
//! `gmx-rs-tools` decodes every run input file GROMACS itself accepts (tpx
//! version >= 58, i.e. GROMACS 4.0 and later), so old files no longer need an
//! external program.  The GROMACS binary named by `gmx_path` in `settings.ini`
//! is only used as a fallback when the built-in reader cannot decode a file at
//! all.  That setting is off by default.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

use gmx_rs_tools::{cmd, index, tpr};

use crate::parse_tpr::TPR;
use crate::settings::Settings;

/// Loads the system description of a run input file, in-process whenever the
/// file can be decoded and through `gmx dump` otherwise.
pub fn load_tpr(tpr_path: &str, settings: &Settings) -> TPR {
    let reason = match TPR::from_run_input(tpr_path) {
        Ok(tpr) => return tpr,
        Err(e) => format!("the built-in reader could not read it ({})", e),
    };
    let gmx = match gmx_program(settings) {
        Some(gmx) => gmx,
        None => {
            println!("Error: {} cannot be read: {}.", tpr_path, reason);
            no_gmx_advice();
        }
    };
    println!("Note: {}: {}; falling back to {}.", tpr_path, reason, gmx);
    let dump_to = dump_path(tpr_path);
    dump_tpr(tpr_path, &dump_to, &gmx, settings);
    TPR::from(&dump_to)
}

/// `.dump` file name used by the `gmx dump` fallback, next to the working
/// directory like the previous releases did it.
fn dump_path(tpr_path: &str) -> String {
    let stem = Path::new(tpr_path).file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("md");
    std::env::current_dir().unwrap()
        .join(format!("{}.dump", stem))
        .display().to_string()
}

/// Generates the default index groups of a run input file.
///
/// This is `gmx make_ndx -f <f> -o <o>` followed by `q`.  Quitting the editor
/// right away writes the default groups unchanged, which is what the index
/// generator does directly.
pub fn make_ndx(options: &[&str], wd: &Path, settings: &Settings, f: &str, n: &str, o: &str) {
    let quits = options.iter().all(|opt| opt.eq_ignore_ascii_case("q"));
    if quits {
        match default_index_groups(f, n, o) {
            Ok(()) => return,
            Err(e) => {
                println!("Error: unable to generate {}:\n{}", o, e);
                if gmx_program(settings).is_none() {
                    return;
                }
            }
        }
    } else {
        let mut args: Vec<String> = vec!["-f".into(), f.into(), "-o".into(), o.into()];
        if !n.is_empty() {
            args.extend(["-n".to_string(), n.to_string()]);
        }
        cmd::set_scripted_input(options.iter().map(|s| s.to_string()));
        let code = cmd::make_ndx::run(args);
        cmd::clear_scripted_input();
        if code == 0 {
            return;
        }
        if gmx_program(settings).is_none() {
            println!("Error: the built-in make_ndx could not write {}.", o);
            return;
        }
    }
    let gmx = gmx_or_advice(settings, f);
    println!("Note: falling back to {}.", gmx);
    let args = match n.is_empty() {
        true => ["make_ndx", "-f", f, "-o", o, "-quiet"].to_vec(),
        false => ["make_ndx", "-f", f, "-n", n, "-o", o, "-quiet"].to_vec(),
    };
    run_gmx(&gmx, options, &args, wd, settings);
}

/// Writes `index::analyse()`'s default groups, the result of an editor session
/// that is quit immediately.
fn default_index_groups(tpr_path: &str, n: &str, o: &str) -> Result<(), String> {
    let file = tpr::TprFile::read(tpr_path).map_err(|e| e.to_string())?;
    let body = tpr::parse_body(&file.header, &file.body).map_err(|e| e.to_string())?;
    let mtop = body.mtop.ok_or("the run input file does not contain a topology")?;
    let groups = if n.is_empty() {
        // No old index file: the editor starts from the default groups.
        index::analyse(&mtop.global_atoms(), false)
    } else {
        index::read_ndx(n).map_err(|e| e.to_string())?
    };
    index::write_ndx(o, &groups, false, mtop.natoms).map_err(|e| e.to_string())
}

/// Extracts and manipulates a trajectory.
pub fn trjconv(options: &[&str], wd: &Path, settings: &Settings, f: &str, s: &str, n: &str,
               o: &str, others: &[&str]) {
    let args: Vec<String> = ["-f", f, "-s", s, "-n", n, "-o", o].iter()
        .map(|a| a.to_string())
        .chain(others.iter().map(|a| a.to_string()))
        .collect();
    cmd::set_scripted_input(options.iter().map(|s| s.to_string()));
    let code = cmd::trjconv::run(args);
    cmd::clear_scripted_input();
    if code == 0 {
        return;
    }
    if gmx_program(settings).is_none() {
        println!("Error: the built-in trjconv could not write {}.", o);
        return;
    }
    let gmx = gmx_or_advice(settings, s);
    println!("Note: falling back to {}.", gmx);
    let args: Vec<&str> = ["trjconv", "-f", f, "-s", s, "-n", n, "-o", o, "-quiet"].iter()
        .chain(others.iter()).cloned().collect();
    run_gmx(&gmx, options, &args, wd, settings);
}

/// Writes a subset run input file, `gmx convert-tpr -n`.
pub fn convert_tpr(options: &[&str], wd: &Path, settings: &Settings, s: &str, n: &str, o: &str) {
    let args: Vec<String> = ["-s", s, "-n", n, "-o", o].iter()
        .map(|a| a.to_string()).collect();
    cmd::set_scripted_input(options.iter().map(|s| s.to_string()));
    let code = cmd::convert_tpr::run(args);
    cmd::clear_scripted_input();
    if code == 0 {
        return;
    }
    if gmx_program(settings).is_none() {
        println!("Error: the built-in convert-tpr could not write {}.", o);
        return;
    }
    let gmx = gmx_or_advice(settings, s);
    println!("Note: falling back to {}.", gmx);
    let args = ["convert-tpr", "-s", s, "-n", n, "-o", o, "-quiet"];
    run_gmx(&gmx, options, &args, wd, settings);
}

// --- GROMACS binary fallback ------------------------------------------------
//
// Only used when the built-in reader cannot decode a run input file at all and
// `gmx_path` names a GROMACS program.

/// `gmx dump -s <tpr>`, written to `dump_to`.
fn dump_tpr(tpr_path: &str, dump_to: &str, gmx: &str, settings: &Settings) {
    if settings.debug_mode {
        println!("CMD: {} dump -s {}", gmx, tpr_path);
    }
    let dump = match Command::new(&gmx)
        .arg("dump").arg("-s").arg(tpr_path)
        .output()
    {
        Ok(dump) => dump,
        Err(e) => {
            println!("Error: cannot run `{} dump -s {}` ({}).", gmx, tpr_path, e);
            no_gmx_advice();
        }
    };
    if !dump.status.success() {
        println!("Error: `{} dump -s {}` failed.", gmx, tpr_path);
        no_gmx_advice();
    }
    let dump = String::from_utf8_lossy(&dump.stdout);
    let mut outfile = std::fs::File::create(dump_to)
        .unwrap_or_else(|_| panic!("Cannot create {}.", dump_to));
    outfile.write_all(dump.as_bytes())
        .unwrap_or_else(|_| panic!("Cannot write {}.", dump_to));
    println!("Dumped tpr file to {}", dump_to);
}

/// Runs a GROMACS tool, feeding `options` to its standard input.
fn run_gmx(gmx: &str, options: &[&str], args: &[&str], wd: &Path, settings: &Settings) {
    if settings.debug_mode {
        println!("CMD: {} {}", gmx, args.join(" "));
    }
    let mut child = match Command::new(&gmx)
        .args(args)
        .current_dir(wd)
        .stdin(Stdio::piped())
        .stdout(if settings.debug_mode { Stdio::inherit() } else { Stdio::null() })
        .spawn()
    {
        Ok(child) => child,
        Err(e) => {
            println!("Error: cannot run `{} {}` ({}).", gmx, args.join(" "), e);
            no_gmx_advice();
        }
    };
    if let Some(stdin) = child.stdin.as_mut() {
        for option in options {
            if settings.debug_mode {
                println!("Input: {}", option);
            }
            writeln!(stdin, "{}", option).unwrap();
        }
    }
    child.wait().expect("gmx failed.");
}

/// The GROMACS program named by `gmx_path`, when one is configured.
///
/// The setting is off by default: a missing or empty entry means that no
/// GROMACS program is used at all.
fn gmx_program(settings: &Settings) -> Option<String> {
    let path = settings.gmx_path.as_deref().map(str::trim).unwrap_or("");
    if path.is_empty() {
        None
    } else {
        Some(path.to_string())
    }
}

/// The configured GROMACS program, or an explanation of what to do without it.
fn gmx_or_advice(settings: &Settings, tpr_path: &str) -> String {
    match gmx_program(settings) {
        Some(gmx) => gmx,
        None => {
            println!(
                "Error: {} needs GROMACS, which is not configured.",
                tpr_path
            );
            no_gmx_advice();
        }
    }
}

/// Explains what to do when GROMACS is needed but not configured.
fn no_gmx_advice() -> ! {
    println!(
        "s_mmpbsa needs GROMACS only for run input files written before GROMACS 2021.\n\
         Either set `gmx_path` in settings.ini to the gmx program, or re-convert the run\n\
         input file once on a machine that has GROMACS with\n\
         `gmx convert-tpr -s old.tpr -o new.tpr` and use the converted file."
    );
    println!("Press ENTER to exit.");
    let _ = std::io::stdin().read_line(&mut String::new());
    std::process::exit(0);
}
