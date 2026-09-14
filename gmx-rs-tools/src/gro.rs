//! GROMOS-87 (.gro) coordinate files, mirroring `gromacs/fileio/groio.cpp`.


use crate::frame::{Atoms, Frame, Matrix, ResInfo, Atom};
use crate::xdr::{Result as XResult, XdrError};

/// Parses one fixed-format float out of `chunk`, returning `None` on failure.
fn parse_chunk(chunk: &str) -> Option<f64> {
    let trimmed = chunk.trim();
    if trimmed.is_empty() {
        return None;
    }
    trimmed.parse::<f64>().ok()
}

struct GroFrameRaw {
    title: String,
    natoms: usize,
    names: Vec<(String, String, u8, i32)>, // atomname, resname, ic, resnr
    x: Vec<[f32; 3]>,
    v: Option<Vec<[f32; 3]>>,
    boxm: Matrix,
    have_box: bool,
}

/// Reads one frame starting at line index `*pos`, advancing it.
fn read_one(lines: &[String], pos: &mut usize) -> XResult<Option<GroFrameRaw>> {
    if *pos >= lines.len() {
        return Ok(None);
    }
    let title = lines[*pos].clone();
    *pos += 1;
    if *pos >= lines.len() {
        return Ok(None);
    }
    let natoms_line = &lines[*pos];
    *pos += 1;
    let natoms: i64 = natoms_line
        .trim()
        .parse()
        .map_err(|_| XdrError::Invalid("gro file does not have the number of atoms on the second line".into()))?;
    let natoms = natoms as usize;

    let mut names = Vec::with_capacity(natoms);
    let mut x = Vec::with_capacity(natoms);
    let mut v = Vec::with_capacity(natoms);
    let mut have_v = false;
    let mut ddist = 0usize;
    let mut first = true;

    for i in 0..natoms {
        if *pos >= lines.len() {
            return Err(XdrError::Invalid(format!(
                "Unexpected end of file in gro file at line {}",
                i + 2
            )));
        }
        let line = lines[*pos].clone();
        *pos += 1;
        if line.len() < 39 {
            return Err(XdrError::Invalid(format!(
                "Invalid line in gro file for atom {}:\n{}",
                i + 1,
                line
            )));
        }
        if first {
            first = false;
            let bytes = line.as_bytes();
            let p1 = bytes
                .iter()
                .position(|&c| c == b'.')
                .ok_or_else(|| XdrError::Invalid("A coordinate in gro file does not contain a '.'".into()))?;
            let p2 = bytes[p1 + 1..]
                .iter()
                .position(|&c| c == b'.')
                .map(|p| p + p1 + 1)
                .ok_or_else(|| XdrError::Invalid("A coordinate in gro file does not contain a '.'".into()))?;
            ddist = p2 - p1;
            if ddist < 5 {
                ddist = 5;
            }
        }

        let resnr: i32 = line[0..5].trim().parse().unwrap_or(0);
        let resname: String = line[5..10].trim().to_string();
        let atomname: String = line[10..15].trim().to_string();
        names.push((atomname, resname, b' ', resnr));

        let bytes = line.as_bytes();
        let mut ptr = 20usize;
        let mut vals = [0f64; 3];
        for slot in vals.iter_mut() {
            let end = std::cmp::min(ptr + ddist, bytes.len());
            let chunk = std::str::from_utf8(&bytes[ptr..end]).unwrap_or("");
            ptr = end;
            *slot = parse_chunk(chunk).unwrap_or(0.0);
        }
        x.push([vals[0] as f32, vals[1] as f32, vals[2] as f32]);

        // Velocities are optional and only present on longer lines.
        let mut vv = [0f32; 3];
        let mut got_v = false;
        for slot in vv.iter_mut() {
            if ptr >= bytes.len() {
                break;
            }
            let end = std::cmp::min(ptr + ddist, bytes.len());
            let chunk = std::str::from_utf8(&bytes[ptr..end]).unwrap_or("");
            ptr = end;
            match parse_chunk(chunk) {
                Some(val) => {
                    *slot = val as f32;
                    got_v = true;
                }
                None => *slot = 0.0,
            }
        }
        if got_v {
            have_v = true;
        }
        v.push(vv);
    }

    let mut boxm = [[0.0f32; 3]; 3];
    if *pos < lines.len() {
        let boxline = &lines[*pos];
        *pos += 1;
        let vals: Vec<f64> = boxline
            .split_whitespace()
            .filter_map(|t| t.parse::<f64>().ok())
            .collect();
        if vals.len() < 3 {
            return Err(XdrError::Invalid("Bad box in gro file".into()));
        }
        boxm[0][0] = vals[0] as f32;
        boxm[1][1] = vals[1] as f32;
        boxm[2][2] = vals[2] as f32;
        if vals.len() >= 9 {
            boxm[0][1] = vals[3] as f32;
            boxm[0][2] = vals[4] as f32;
            boxm[1][0] = vals[5] as f32;
            boxm[1][2] = vals[6] as f32;
            boxm[2][0] = vals[7] as f32;
            boxm[2][1] = vals[8] as f32;
        }
    }

    Ok(Some(GroFrameRaw {
        title,
        natoms,
        names,
        x,
        v: if have_v { Some(v) } else { None },
        boxm,
        have_box: true,
    }))
}

/// Reads every frame of a .gro file.
pub fn read_all(path: &str) -> XResult<Vec<Frame>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
    let lines: Vec<String> = content.lines().map(|l| l.to_string()).collect();

    let mut frames = Vec::new();
    let mut pos = 0usize;
    while let Some(raw) = read_one(&lines, &mut pos)? {
        let mut frame = Frame::new(raw.natoms);
        frame.x = Some(raw.x.clone());
        frame.v = raw.v.clone();
        if raw.have_box {
            frame.boxm = Some(raw.boxm);
        }
        frame.prec = Some(1.0);
        frame.pbc_type = crate::frame::PbcType::Unset;
        frame.title = raw.title.clone();
        frame.time = parse_title_time(&raw.title);
        frame.step = parse_title_step(&raw.title);

        // Build the atom/residue bookkeeping the same way get_w_conf() does.
        let mut atoms = Atoms::default();
        atoms.name = raw.title.clone();
        let mut oldres = -1i32;
        let mut newres = -1i32;
        let mut oldresname = String::new();
        for (atomname, resname, ic, resnr) in raw.names.iter() {
            if newres == -1 || *resnr != oldres || *resname != oldresname {
                oldres = *resnr;
                newres += 1;
                atoms.resinfo.push(ResInfo {
                    name: resname.clone(),
                    nr: *resnr,
                    ic: *ic,
                    chainid: b' ',
                });
            }
            oldresname = resname.clone();
            atoms.atom.push(Atom {
                name: atomname.clone(),
                atom_type: atomname.clone(),
                atom_type_b: String::new(),
                resind: newres,
                mass: 0.0,
                charge: 0.0,
                mass_b: 0.0,
                charge_b: 0.0,
                ptype: 0,
                atomnumber: 0,
                elem: String::new(),
                type_id: 0,
                type_id_b: 0,
            });
        }
        frame.atoms = Some(atoms);
        frames.push(frame);
    }
    Ok(frames)
}

fn parse_title_time(title: &str) -> Option<f64> {
    let idx = title.find("t=")?;
    let rest = &title[idx + 2..];
    let token: String = rest
        .trim_start()
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '.' || *c == '-' || *c == '+' || *c == 'e' || *c == 'E')
        .collect();
    token.parse::<f64>().ok()
}

fn parse_title_step(title: &str) -> Option<i64> {
    let idx = title.find("step=")?;
    let rest = &title[idx + 5..];
    rest.trim_start()
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+')
        .collect::<String>()
        .parse::<i64>()
        .ok()
}

/// Writes the boxm line exactly like `write_hconf_box()`.
pub fn write_box<W: std::io::Write>(out: &mut W, boxm: &Matrix) {
    let triclinic = boxm[0][1] != 0.0
        || boxm[0][2] != 0.0
        || boxm[1][0] != 0.0
        || boxm[1][2] != 0.0
        || boxm[2][0] != 0.0
        || boxm[2][1] != 0.0;
    if triclinic {
        let _ = writeln!(
            out,
            "{:10.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5}",
            boxm[0][0], boxm[1][1], boxm[2][2], boxm[0][1], boxm[0][2], boxm[1][0], boxm[1][2],
            boxm[2][0], boxm[2][1]
        );
    } else {
        let _ = writeln!(out, "{:10.5} {:9.5} {:9.5}", boxm[0][0], boxm[1][1], boxm[2][2]);
    }
}

/// Precomputes the fixed part of every atom line.  `write_hconf_indexed_p()`
/// writes residue number, residue name, atom name and atom number, none of
/// which change between frames, so they only need to be formatted once.
pub fn atom_prefixes(atoms: &Atoms, index: &[usize]) -> Vec<Vec<u8>> {
    let mut prefixes = Vec::with_capacity(index.len());
    for &ai in index {
        let resind = atoms.atom[ai].resind as usize;
        let (resnm, resnr) = if resind < atoms.resinfo.len() {
            (atoms.resinfo[resind].name.clone(), atoms.resinfo[resind].nr)
        } else {
            (" ??? ".to_string(), resind as i32 + 1)
        };
        let nm = atoms.atom[ai].name.clone();
        let mut s = Vec::with_capacity(20);
        // "%5d%-5.5s%5.5s%5d"
        push_int_padded(&mut s, resnr % 100000, 5);
        push_field(&mut s, &resnm, 5);
        push_field_right(&mut s, &nm, 5);
        push_int_padded(&mut s, (ai as i32 + 1) % 100000, 5);
        prefixes.push(s);
    }
    prefixes
}

/// `%*d`: right justified integer, equivalent to the C conversions GROMACS
/// uses for the fixed columns of a .gro file.
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

/// `%-5.5s`: at most five characters, padded on the right.
fn push_field(out: &mut Vec<u8>, s: &str, width: usize) {
    let bytes = s.as_bytes();
    let take = bytes.len().min(width);
    out.extend_from_slice(&bytes[..take]);
    for _ in take..width {
        out.push(b' ');
    }
}

/// `%5.5s`: at most five characters, padded on the left.
fn push_field_right(out: &mut Vec<u8>, s: &str, width: usize) {
    let bytes = s.as_bytes();
    let take = bytes.len().min(width);
    if bytes.len() < width {
        for _ in bytes.len()..width {
            out.push(b' ');
        }
    }
    out.extend_from_slice(&bytes[..take]);
}

/// Writes one frame in `.gro` format, mirroring `write_hconf_indexed_p()`.
pub fn write_frame<W: std::io::Write>(
    out: &mut W,
    title: &str,
    prefixes: &[Vec<u8>],
    x: &[[f32; 3]],
    index: &[usize],
    v: Option<&[[f32; 3]]>,
    boxm: &Matrix,
) {
    let _ = writeln!(
        out,
        "{}",
        if title.is_empty() {
            "GROningen MAchine for Chemical Simulation"
        } else {
            title
        }
    );
    let _ = writeln!(out, "{:5}", prefixes.len());
    for (i, prefix) in prefixes.iter().enumerate() {
        let ai = if index.len() == prefixes.len() {
            index[i]
        } else {
            i
        };
        let _ = out.write_all(prefix);
        let c = x[ai];
        match v {
            Some(v) => {
                let _ = writeln!(
                    out,
                    "{:8.3}{:8.3}{:8.3}{:8.4}{:8.4}{:8.4}",
                    c[0], c[1], c[2], v[ai][0], v[ai][1], v[ai][2]
                );
            }
            None => {
                let _ = writeln!(out, "{:8.3}{:8.3}{:8.3}", c[0], c[1], c[2]);
            }
        }
    }
    write_box(out, boxm);
}
