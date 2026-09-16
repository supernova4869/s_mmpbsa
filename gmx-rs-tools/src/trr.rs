//! TRR trajectory format, mirroring `gromacs/fileio/trrio.cpp`.
//!
//! ```text
//! int    magic      1993
//! string "GMX_trn_file"
//! int    ir_size, e_size, box_size, vir_size, pres_size, top_size, sym_size
//! int    x_size, v_size, f_size, natoms
//! int    step, nre
//! real   t, lambda
//! <boxm, x, v, f as real vectors>
//! ```

use crate::frame::Frame;
use crate::xdr::{Reader, Result, Writer, XdrError};

pub const TRR_MAGIC: i32 = 1993;
const TRR_VERSION: &str = "GMX_trn_file";

/// Writes a single TRR frame in single precision.
pub fn write_frame(w: &mut Writer, frame: &Frame) {
    let natoms = frame.x.as_ref().map(|x| x.len()).unwrap_or(frame.natoms);
    let real_size = 4usize;

    w.int(TRR_MAGIC);
    w.gmx_string(TRR_VERSION);
    w.int(0); // ir_size
    w.int(0); // e_size
    w.int(if frame.boxm.is_some() { 9 * real_size as i32 } else { 0 });
    w.int(0); // vir_size
    w.int(0); // pres_size
    w.int(0); // top_size
    w.int(0); // sym_size
    w.int(if frame.x.is_some() {
        (natoms * 3 * real_size) as i32
    } else {
        0
    });
    w.int(if frame.v.is_some() {
        (natoms * 3 * real_size) as i32
    } else {
        0
    });
    w.int(if frame.f.is_some() {
        (natoms * 3 * real_size) as i32
    } else {
        0
    });
    w.int(natoms as i32);
    w.int(frame.step.unwrap_or(0) as i32);
    w.int(0); // nre
    w.float(frame.time.unwrap_or(0.0) as f32);
    w.float(frame.lambda.unwrap_or(0.0));

    if let Some(boxm) = &frame.boxm {
        for row in boxm.iter() {
            for c in row.iter() {
                w.float(*c);
            }
        }
    }
    if let Some(x) = &frame.x {
        for c in x {
            for v in c {
                w.float(*v);
            }
        }
    }
    if let Some(v) = &frame.v {
        for c in v {
            for e in c {
                w.float(*e);
            }
        }
    }
    if let Some(f) = &frame.f {
        for c in f {
            for e in c {
                w.float(*e);
            }
        }
    }
}

fn n_float_size(box_size: i32, x_size: i32, v_size: i32, f_size: i32, natoms: i32) -> Result<i32> {
    let nflsize = if box_size != 0 {
        box_size / 9
    } else if x_size != 0 {
        x_size / (natoms * 3)
    } else if v_size != 0 {
        v_size / (natoms * 3)
    } else if f_size != 0 {
        f_size / (natoms * 3)
    } else {
        return Err(XdrError::Invalid("Can not determine precision of trr file".into()));
    };
    if nflsize != 4 && nflsize != 8 {
        return Err(XdrError::Invalid(format!("Float size {nflsize}. Maybe different CPU?")));
    }
    Ok(nflsize)
}

/// Reads one TRR frame; returns `None` at clean end of file.
pub fn read_frame(r: &mut Reader) -> Result<Option<Frame>> {
    if r.remaining() < 4 {
        return Ok(None);
    }
    let magic = r.int()?;
    if magic != TRR_MAGIC {
        return Err(XdrError::Invalid(
            "Failed to find GROMACS magic number in trr frame header, so this is not a trr file!"
                .into(),
        ));
    }
    let _version = r.gmx_string()?;
    let ir_size = r.int()?;
    let e_size = r.int()?;
    let box_size = r.int()?;
    let vir_size = r.int()?;
    let pres_size = r.int()?;
    let top_size = r.int()?;
    let sym_size = r.int()?;
    let x_size = r.int()?;
    let v_size = r.int()?;
    let f_size = r.int()?;
    let natoms = r.int()? as usize;
    let step = r.int()? as i64;
    let _nre = r.int()?;

    let double_precision = n_float_size(box_size, x_size, v_size, f_size, natoms as i32)? == 8;
    let t = r.real(double_precision)?;
    let lambda = r.real(double_precision)? as f32;

    if ir_size != 0 || e_size != 0 || top_size != 0 || sym_size != 0 {
        return Err(XdrError::Invalid(
            "trr file contains inputrec/energies/topology/symbol table".into(),
        ));
    }

    let mut frame = Frame::new(natoms);
    frame.step = Some(step);
    frame.time = Some(t);
    frame.lambda = Some(lambda);

    if box_size != 0 {
        let mut boxm = [[0.0f32; 3]; 3];
        for row in boxm.iter_mut() {
            for c in row.iter_mut() {
                *c = r.real(double_precision)? as f32;
            }
        }
        frame.boxm = Some(boxm);
    }
    if vir_size != 0 {
        for _ in 0..9 {
            let _ = r.real(double_precision)?;
        }
    }
    if pres_size != 0 {
        for _ in 0..9 {
            let _ = r.real(double_precision)?;
        }
    }
    if x_size != 0 {
        let mut x = Vec::with_capacity(natoms);
        for _ in 0..natoms {
            x.push([
                r.real(double_precision)? as f32,
                r.real(double_precision)? as f32,
                r.real(double_precision)? as f32,
            ]);
        }
        frame.x = Some(x);
    }
    if v_size != 0 {
        let mut v = Vec::with_capacity(natoms);
        for _ in 0..natoms {
            v.push([
                r.real(double_precision)? as f32,
                r.real(double_precision)? as f32,
                r.real(double_precision)? as f32,
            ]);
        }
        frame.v = Some(v);
    }
    if f_size != 0 {
        let mut f = Vec::with_capacity(natoms);
        for _ in 0..natoms {
            f.push([
                r.real(double_precision)? as f32,
                r.real(double_precision)? as f32,
                r.real(double_precision)? as f32,
            ]);
        }
        frame.f = Some(f);
    }

    Ok(Some(frame))
}
/// Reads the raw bytes of one TRR frame from a stream.
///
/// All section sizes are part of the header, so the exact frame length is known
/// before the payload is read.  Returns `Ok(None)` at a clean end of file.
/// Reads the fixed part of one frame: the frame header with the sizes of all
/// sections.
///
/// Returns the bytes read together with the number of payload bytes that
/// follow them, or `None` at a clean end of file.  Splitting a frame this way
/// lets [`read_frame_bytes`] decode it and [`frame_count`] skip the payload of
/// a trajectory that is never decoded.
fn read_frame_prefix<R: std::io::Read>(r: &mut R) -> Result<Option<(Vec<u8>, usize)>> {
    use crate::xdr::{read_exact, read_or_eof, xdr_pad};

    let mut buf: Vec<u8> = Vec::with_capacity(4096);
    let mut magic_bytes = [0u8; 4];
    if !read_or_eof(r, &mut magic_bytes)? {
        return Ok(None);
    }
    buf.extend_from_slice(&magic_bytes);
    let magic = i32::from_be_bytes(magic_bytes);
    if magic != TRR_MAGIC {
        return Err(XdrError::Invalid(
            "Failed to find GROMACS magic number in trr frame header, so this is not a trr file!"
                .into(),
        ));
    }

    // Version string: len + 1, len, characters, padding.
    let mut int4 = [0u8; 4];
    read_exact(r, &mut int4)?;
    buf.extend_from_slice(&int4);
    read_exact(r, &mut int4)?;
    buf.extend_from_slice(&int4);
    let len = i32::from_be_bytes(int4).max(0) as usize;
    let mut s = vec![0u8; len + xdr_pad(len)];
    read_exact(r, &mut s)?;
    buf.extend_from_slice(&s);

    // The seven section sizes plus the x/v/f sizes and the atom count.
    let mut sizes = [0u8; 11 * 4];
    read_exact(r, &mut sizes)?;
    buf.extend_from_slice(&sizes);
    let get = |i: usize| {
        i32::from_be_bytes([
            sizes[i * 4],
            sizes[i * 4 + 1],
            sizes[i * 4 + 2],
            sizes[i * 4 + 3],
        ])
    };
    let float_size = n_float_size(get(2), get(7), get(8), get(9), get(10))? as usize;

    // step, nre, then t and lambda (whose width follows the file precision).
    let mut tail = [0u8; 8];
    read_exact(r, &mut tail)?;
    buf.extend_from_slice(&tail);
    let mut reals = vec![0u8; 2 * float_size];
    read_exact(r, &mut reals)?;
    buf.extend_from_slice(&reals);

    // Payload: box, virial, pressure, x, v, f.
    let payload: usize = (1..=6).map(|i| get(i).max(0) as usize).sum::<usize>()
        + get(7).max(0) as usize
        + get(8).max(0) as usize
        + get(9).max(0) as usize;
    Ok(Some((buf, payload)))
}

pub fn read_frame_bytes<R: std::io::Read>(r: &mut R) -> Result<Option<Vec<u8>>> {
    use crate::xdr::read_exact;

    let Some((mut buf, payload)) = read_frame_prefix(r)? else {
        return Ok(None);
    };
    let mut data = vec![0u8; payload];
    read_exact(r, &mut data)?;
    buf.extend_from_slice(&data);
    Ok(Some(buf))
}

/// Counts the frames of a TRR file.
///
/// Only the fixed part of every frame is read; the coordinates, velocities and
/// forces are skipped with a seek, so a trajectory much larger than memory is
/// counted without decoding a single frame.
pub fn frame_count(path: &str) -> Result<usize> {
    use std::io::{BufReader, Seek, SeekFrom};

    let file = std::fs::File::open(path)
        .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
    let mut r = BufReader::with_capacity(1 << 16, file);
    let mut frames = 0usize;
    while let Some((_, payload)) = read_frame_prefix(&mut r)? {
        r.seek(SeekFrom::Current(payload as i64))
            .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
        frames += 1;
    }
    Ok(frames)
}

/// Decodes a frame produced by [`read_frame_bytes`].
pub fn decode_frame(bytes: &[u8]) -> Result<Frame> {
    let mut r = Reader::new(bytes);
    read_frame(&mut r)?.ok_or_else(|| XdrError::Invalid("empty trr frame".into()))
}
