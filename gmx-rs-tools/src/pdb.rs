//! PDB coordinate files, mirroring the parts of `gromacs/fileio/pdbio.cpp`
//! used by `gmx trjconv` and `gmx make_ndx`.


use crate::frame::{Atoms, Atom, Frame, Matrix, PbcType, ResInfo};
use crate::xdr::{Result as XResult, XdrError};

fn sub(line: &str, from: usize, to: usize) -> String {
    let b = line.as_bytes();
    if from >= b.len() {
        return String::new();
    }
    let end = std::cmp::min(to, b.len());
    String::from_utf8_lossy(&b[from..end]).trim().to_string()
}

fn gmx_angle(a: [f32; 3], b: [f32; 3]) -> f64 {
    let dot = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) as f64;
    let na = ((a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) as f64).sqrt();
    let nb = ((b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) as f64).sqrt();
    if na == 0.0 || nb == 0.0 {
        return 90.0;
    }
    (dot / (na * nb)).clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI
}

/// `norm()` from `gromacs/utility/vec.h`, computed in single precision like
/// GROMACS does for a `real == float` build.
fn norm(a: [f32; 3]) -> f32 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Reads all frames of a PDB file. Only ATOM/HETATM/CRYST1/MODEL/TITLE are
/// interpreted; everything else is skipped.
pub fn read_all(path: &str) -> XResult<Vec<Frame>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;

    let mut frames: Vec<Frame> = Vec::new();
    let mut current_atoms = Atoms::default();
    let mut current_x: Vec<[f32; 3]> = Vec::new();
    let mut current_box = [[0.0f32; 3]; 3];
    let mut have_box = false;
    let mut title = String::new();
    let mut model_nr: Option<i64> = None;
    let mut started = false;

    let finish = |frames: &mut Vec<Frame>,
                  atoms: &mut Atoms,
                  x: &mut Vec<[f32; 3]>,
                  boxm: &mut Matrix,
                  have_box: bool,
                  title: &str,
                  model_nr: Option<i64>| {
        if atoms.atom.is_empty() {
            return;
        }
        let natoms = atoms.atom.len();
        let mut frame = Frame::new(natoms);
        frame.x = Some(x.clone());
        frame.boxm = if have_box { Some(*boxm) } else { None };
        frame.pbc_type = if have_box { PbcType::Xyz } else { PbcType::No };
        frame.atoms = Some(atoms.clone());
        frame.title = title.to_string();
        frame.step = model_nr;
        // GROMACS numbers models starting from 1.
        frame.time = None;
        frames.push(frame);
        atoms.atom.clear();
        atoms.resinfo.clear();
        x.clear();
        boxm.iter_mut().for_each(|r| r.iter_mut().for_each(|c| *c = 0.0));
    };

    for (line_no, line) in content.lines().enumerate() {
        let _ = line_no;
        if line.len() >= 6 {
            let rec = &line[0..6];
            if rec == "ENDMDL" {
                finish(
                    &mut frames,
                    &mut current_atoms,
                    &mut current_x,
                    &mut current_box,
                    have_box,
                    &title,
                    model_nr,
                );
                started = false;
                continue;
            }
        }
        if line.starts_with("MODEL") {
            let n: i64 = line[5..].trim().parse().unwrap_or(0);
            model_nr = Some(n);
            started = true;
            continue;
        }
        if line.starts_with("CRYST1") {
            let toks: Vec<&str> = line[6..].split_whitespace().collect();
            if toks.len() >= 6 {
                let a: f32 = toks[0].parse().unwrap_or(0.0);
                let b: f32 = toks[1].parse().unwrap_or(0.0);
                let c: f32 = toks[2].parse().unwrap_or(0.0);
                let alpha: f64 = toks[3].parse().unwrap_or(90.0);
                let beta: f64 = toks[4].parse().unwrap_or(90.0);
                let gamma: f64 = toks[5].parse().unwrap_or(90.0);
                let (fa, fb, fc) = (a as f64 * 0.1, b as f64 * 0.1, c as f64 * 0.1);
                current_box = [[0.0; 3]; 3];
                current_box[0][0] = fa as f32;
                if alpha != 90.0 || beta != 90.0 || gamma != 90.0 {
                    let cosa = if alpha != 90.0 {
                        alpha.to_radians().cos()
                    } else {
                        0.0
                    };
                    let cosb = if beta != 90.0 {
                        beta.to_radians().cos()
                    } else {
                        0.0
                    };
                    let (cosg, sing) = if gamma != 90.0 {
                        (gamma.to_radians().cos(), gamma.to_radians().sin())
                    } else {
                        (0.0, 1.0)
                    };
                    let zx = fc * cosb;
                    let zy = fc * (cosa - cosb * cosg) / sing;
                    let zz = (fc * fc - zx * zx - zy * zy).sqrt();
                    current_box[1][0] = (fb * cosg) as f32;
                    current_box[1][1] = (fb * sing) as f32;
                    current_box[2][0] = zx as f32;
                    current_box[2][1] = zy as f32;
                    current_box[2][2] = zz as f32;
                } else {
                    current_box[1][1] = fb as f32;
                    current_box[2][2] = fc as f32;
                }
                have_box = true;
                started = true;
            }
            continue;
        }
        if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            let atomname = sub(line, 12, 16);
            let resname = sub(line, 17, 20);
            let chainid = line.as_bytes().get(21).copied().unwrap_or(b' ');
            let resnr: i32 = sub(line, 22, 26).parse().unwrap_or(0);
            let resic = line.as_bytes().get(26).copied().unwrap_or(b' ');
            let xc: f32 = sub(line, 30, 38).parse().unwrap_or(0.0);
            let yc: f32 = sub(line, 38, 46).parse().unwrap_or(0.0);
            let zc: f32 = sub(line, 46, 54).parse().unwrap_or(0.0);
            let elem = sub(line, 76, 78);

            let newres = current_atoms.resinfo.last().map_or(true, |ri| {
                ri.nr != resnr || ri.ic != resic || ri.name != resname
            });
            let resind = if newres {
                current_atoms.resinfo.push(ResInfo {
                    name: resname.clone(),
                    nr: resnr,
                    ic: resic,
                    chainid,
                });
                current_atoms.resinfo.len() as i32 - 1
            } else {
                current_atoms.resinfo.len() as i32 - 1
            };
            current_atoms.atom.push(Atom {
                name: atomname.clone(),
                atom_type: atomname,
                atom_type_b: String::new(),
                resind,
                mass: 0.0,
                charge: 0.0,
                mass_b: 0.0,
                charge_b: 0.0,
                ptype: 0,
                atomnumber: 0,
                elem,
                type_id: 0,
                type_id_b: 0,
            });
            current_x.push([xc / 10.0, yc / 10.0, zc / 10.0]);
            started = true;
            continue;
        }
        if line.starts_with("TITLE") || line.starts_with("HEADER") {
            let rest = line[6..].trim();
            if !rest.is_empty() {
                title = rest.split("      ").next().unwrap_or(rest).trim().to_string();
                started = true;
            }
            continue;
        }
        if line.starts_with("TER") {
            started = true;
        }
    }
    if !current_atoms.atom.is_empty() {
        finish(
            &mut frames,
            &mut current_atoms,
            &mut current_x,
            &mut current_box,
            have_box,
            &title,
            model_nr,
        );
    }
    let _ = started;
    Ok(frames)
}

fn write_pdb_box<W: std::io::Write>(out: &mut W, pbc_type: PbcType, boxm: &Matrix) {
    let alpha = if norm(boxm[1]) * norm(boxm[2]) != 0.0 {
        gmx_angle(boxm[1], boxm[2])
    } else {
        90.0
    };
    let beta = if norm(boxm[0]) * norm(boxm[2]) != 0.0 {
        gmx_angle(boxm[0], boxm[2])
    } else {
        90.0
    };
    let gamma = if norm(boxm[0]) * norm(boxm[1]) != 0.0 {
        gmx_angle(boxm[0], boxm[1])
    } else {
        90.0
    };
    let _ = writeln!(out, "REMARK    THIS IS A SIMULATION BOX");
    let space_group = if pbc_type != PbcType::Screw { "P 1" } else { "P 21 1 1" };
    let a: f32 = if pbc_type != PbcType::Screw { 10.0 } else { 20.0 };
    let _ = writeln!(
        out,
        "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:<11}{:4}",
        (a * norm(boxm[0])) as f64,
        (10.0 * norm(boxm[1])) as f64,
        (10.0 * norm(boxm[2])) as f64,
        alpha,
        beta,
        gamma,
        space_group,
        1
    );
}

/// Builds the fixed prefix and suffix of every atom line.
///
/// The coordinates are the only part of a PDB atom record that changes between
/// frames, so the rest is formatted once per conversion.
pub fn atom_fields(atoms: &Atoms, index: &[usize]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut prefixes = Vec::with_capacity(index.len());
    let mut suffixes = Vec::with_capacity(index.len());
    for &i in index {
        let resind = atoms.atom[i].resind as usize;
        let resnm = atoms
            .resinfo
            .get(resind)
            .map(|r| r.name.clone())
            .unwrap_or_default();
        let nm = atoms.atom[i].name.clone();
        let resnr = atoms.resinfo.get(resind).map(|r| r.nr).unwrap_or(0);
        let resic = atoms.resinfo.get(resind).map(|r| r.ic).unwrap_or(b' ');
        let ch = {
            let c = atoms.resinfo.get(resind).map(|r| r.chainid).unwrap_or(b' ');
            if c == 0 {
                b' '
            } else {
                c
            }
        };
        let elem = atoms.atom[i].elem.clone();
        let start_in_col13 = if !elem.is_empty() && elem.len() >= 2 && nm.len() >= 2 {
            nm[0..2].eq_ignore_ascii_case(&elem[0..2])
        } else {
            nm.chars().count() >= 4
        };
        let mut tmp_atomname = Vec::with_capacity(5);
        if !start_in_col13 {
            tmp_atomname.push(b' ');
        }
        tmp_atomname.extend_from_slice(&nm.as_bytes()[..nm.len().min(4)]);

        let mut p = Vec::with_capacity(32);
        p.extend_from_slice(b"ATOM  ");
        push_int_padded(&mut p, (i as i32 + 1) % 100000, 5);
        p.push(b' ');
        // "%-4.4s"
        for c in tmp_atomname.iter().take(4) {
            p.push(*c);
        }
        for _ in tmp_atomname.len().min(4)..4 {
            p.push(b' ');
        }
        p.push(b' '); // alternate location
        // "%4.4s" of `resnm + " "`: the C conversion right justifies, so a
        // name shorter than three characters is preceded by a space (this is
        // how `NA` ends up as " NA " and `SOL` as "SOL ").
        let printed = if resnm.len() < 4 { resnm.len() + 1 } else { 4 };
        for _ in printed..4 {
            p.push(b' ');
        }
        for c in resnm.as_bytes().iter().take(printed) {
            p.push(*c);
        }
        if resnm.len() < 4 {
            p.push(b' ');
        }
        p.push(ch);
        push_int_padded(&mut p, resnr % 10000, 4);
        p.push(if resic == 0 { b' ' } else { resic });
        p.extend_from_slice(b"   ");
        prefixes.push(p);

        let mut s = Vec::with_capacity(16);
        // occupancy, b-factor and element: "%6.2f%6.2f          %2s"
        s.extend_from_slice(b"  1.00  0.00          ");
        let e = elem.as_bytes();
        for _ in e.len()..2 {
            s.push(b' ');
        }
        s.extend_from_slice(&e[..e.len().min(2)]);
        suffixes.push(s);
    }
    (prefixes, suffixes)
}

/// `%*d`: right justified integer (the C conversion used for PDB columns).
fn push_int_padded(out: &mut Vec<u8>, value: i32, width: usize) {
    let mut buf = [0u8; 12];
    let mut n = 0usize;
    let neg = value < 0;
    let mut v = value.unsigned_abs();
    if v == 0 {
        buf[n] = b'0';
        n += 1;
    }
    while v > 0 {
        buf[n] = b'0' + (v % 10) as u8;
        v /= 10;
        n += 1;
    }
    let digits = n + usize::from(neg);
    for _ in digits..width {
        out.push(b' ');
    }
    if neg {
        out.push(b'-');
    }
    while n > 0 {
        n -= 1;
        out.push(buf[n]);
    }
}

/// Writes one frame, mirroring `write_pdbfile()` with the default (non
/// standard-compliant) settings used by `gmx trjconv`.
pub fn write_frame<W: std::io::Write>(
    out: &mut W,
    title: &str,
    prefixes: &[Vec<u8>],
    suffixes: &[Vec<u8>],
    x: &[[f32; 3]],
    index: &[usize],
    pbc_type: PbcType,
    boxm: &Matrix,
    model_nr: i32,
) {
    let _ = writeln!(out, "REMARK    GENERATED BY TRJCONV");
    let _ = writeln!(
        out,
        "TITLE     {}",
        if title.is_empty() { "GROMACS" } else { title }
    );
    write_pdb_box(out, pbc_type, boxm);
    let _ = writeln!(out, "MODEL {:>8}", if model_nr > 0 { model_nr } else { 1 });

    for (i, prefix) in prefixes.iter().enumerate() {
        let ai = if index.len() == prefixes.len() { index[i] } else { i };
        let _ = out.write_all(prefix);
        let _ = write!(
            out,
            "{:>8.3}{:>8.3}{:>8.3}",
            10.0 * x[ai][0] as f64,
            10.0 * x[ai][1] as f64,
            10.0 * x[ai][2] as f64
        );
        let _ = out.write_all(&suffixes[i]);
        let _ = writeln!(out);
    }
    let _ = writeln!(out, "TER");
    let _ = writeln!(out, "ENDMDL");
}
