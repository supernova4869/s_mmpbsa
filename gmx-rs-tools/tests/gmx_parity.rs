//! Optional parity tests against a real GROMACS installation.
//!
//! These are skipped unless `GMXRS_PARITY_TESTDATA` points at a directory that
//! contains `topol.tpr` and `traj.xtc`, and `GMX_BIN` (or `gmx` on `PATH`) is
//! an executable GROMACS driver:
//!
//! ```text
//! GMXRS_PARITY_TESTDATA=/tmp/gmx-rs-tools/work GMX_BIN=/opt/gromacs/bin/gmx \
//!     cargo test --test gmx_parity -- --nocapture
//! ```

use std::path::PathBuf;
use std::process::Command;

fn gmx_bin() -> Option<String> {
    if let Ok(b) = std::env::var("GMX_BIN") {
        return Some(b);
    }
    for candidate in ["gmx", "/opt/gromacs/bin/gmx"] {
        if Command::new(candidate).arg("--version").output().is_ok() {
            return Some(candidate.to_string());
        }
    }
    None
}

fn testdata() -> Option<PathBuf> {
    let dir = std::env::var("GMXRS_PARITY_TESTDATA").ok()?;
    let p = PathBuf::from(dir);
    if p.join("topol.tpr").exists() && p.join("traj.xtc").exists() {
        Some(p)
    } else {
        None
    }
}

fn run(cmd: &mut Command, stdin: &str) -> Vec<u8> {
    use std::io::Write;
    let mut child = cmd
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("spawn");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(stdin.as_bytes())
        .unwrap();
    child.wait_with_output().expect("wait").stdout
}

#[test]
fn dump_xtc_matches_gmx() {
    let (Some(gmx), Some(dir)) = (gmx_bin(), testdata()) else {
        eprintln!("skipping: set GMXRS_PARITY_TESTDATA and GMX_BIN to run parity tests");
        return;
    };
    let xtc = dir.join("traj.xtc");
    let mine = run(
        Command::new(env!("CARGO_BIN_EXE_gmx-rs-tools")).arg("dump").arg("-f").arg(&xtc),
        "",
    );
    let reference = run(
        Command::new(&gmx).arg("dump").arg("-f").arg(&xtc),
        "",
    );
    assert_eq!(
        String::from_utf8_lossy(&mine),
        String::from_utf8_lossy(&reference),
        "gmx-rs-tools dump -f does not match gmx dump -f"
    );
}

#[test]
fn trjconv_gro_matches_gmx() {
    let (Some(gmx), Some(dir)) = (gmx_bin(), testdata()) else {
        eprintln!("skipping: set GMXRS_PARITY_TESTDATA and GMX_BIN to run parity tests");
        return;
    };
    let out_dir = std::env::temp_dir().join("gmx-rs-tools-parity");
    std::fs::create_dir_all(&out_dir).unwrap();
    let mine_out = out_dir.join("mine.gro");
    let ref_out = out_dir.join("ref.gro");
    let tpr = dir.join("topol.tpr");
    let xtc = dir.join("traj.xtc");

    run(
        Command::new(env!("CARGO_BIN_EXE_gmx-rs-tools"))
            .args(["trjconv", "-f"])
            .arg(&xtc)
            .arg("-s")
            .arg(&tpr)
            .args(["-t0", "0", "-o"])
            .arg(&mine_out),
        "0\n",
    );
    run(
        Command::new(&gmx)
            .args(["trjconv", "-f"])
            .arg(&xtc)
            .arg("-s")
            .arg(&tpr)
            .args(["-t0", "0", "-o"])
            .arg(&ref_out),
        "0\n",
    );
    assert_eq!(
        std::fs::read(&mine_out).unwrap(),
        std::fs::read(&ref_out).unwrap(),
        "gmx-rs-tools trjconv -o out.gro does not match gmx trjconv"
    );
}

#[test]
fn make_ndx_default_groups_match_gmx() {
    let (Some(gmx), Some(dir)) = (gmx_bin(), testdata()) else {
        eprintln!("skipping: set GMXRS_PARITY_TESTDATA and GMX_BIN to run parity tests");
        return;
    };
    let out_dir = std::env::temp_dir().join("gmx-rs-tools-parity");
    std::fs::create_dir_all(&out_dir).unwrap();
    let mine_out = out_dir.join("mine.ndx");
    let ref_out = out_dir.join("ref.ndx");
    let tpr = dir.join("topol.tpr");

    run(
        Command::new(env!("CARGO_BIN_EXE_gmx-rs-tools"))
            .args(["make_ndx", "-f"])
            .arg(&tpr)
            .arg("-o")
            .arg(&mine_out),
        "q\n",
    );
    run(
        Command::new(&gmx)
            .args(["make_ndx", "-f"])
            .arg(&tpr)
            .arg("-o")
            .arg(&ref_out),
        "q\n",
    );
    assert_eq!(
        std::fs::read(&mine_out).unwrap(),
        std::fs::read(&ref_out).unwrap(),
        "default index groups do not match"
    );
}

/// `convert-tpr` subset extraction: the resulting tpr must be equivalent to
/// the one GROMACS writes (compared through `gmx dump -s`, since the two
/// implementations serialize the file independently).
#[test]
fn convert_tpr_subset_matches_gmx() {
    let (Some(gmx), Some(dir)) = (gmx_bin(), testdata()) else {
        eprintln!("skipping: set GMXRS_PARITY_TESTDATA and GMX_BIN to run parity tests");
        return;
    };
    let out_dir = std::env::temp_dir().join("gmx-rs-tools-parity");
    std::fs::create_dir_all(&out_dir).unwrap();
    let mine_out = out_dir.join("mine_subset.tpr");
    let ref_out = out_dir.join("ref_subset.tpr");
    let tpr = dir.join("topol.tpr");

    // Group 1 of the default groups.
    run(
        Command::new(env!("CARGO_BIN_EXE_gmx-rs-tools"))
            .args(["convert-tpr", "-s"])
            .arg(&tpr)
            .arg("-o")
            .arg(&mine_out),
        "1\n",
    );
    run(
        Command::new(&gmx)
            .args(["convert-tpr", "-s"])
            .arg(&tpr)
            .arg("-o")
            .arg(&ref_out),
        "1\n",
    );
    assert_eq!(
        std::fs::metadata(&mine_out).unwrap().len(),
        std::fs::metadata(&ref_out).unwrap().len(),
        "subset tpr has a different size"
    );
    let mine_dump = run(
        Command::new(&gmx).args(["dump", "-s"]).arg(&mine_out),
        "",
    );
    let ref_dump = run(
        Command::new(&gmx).args(["dump", "-s"]).arg(&ref_out),
        "",
    );
    let strip = |v: &[u8]| {
        String::from_utf8_lossy(v)
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(
        strip(&mine_dump),
        strip(&ref_dump),
        "gmx dump -s of the subset tpr differs"
    );
}
