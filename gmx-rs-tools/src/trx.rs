//! Trajectory format dispatch, mirroring the file type detection of
//! `gromacs/fileio/trxio.cpp` and `gromacs/fileio/filetypes.cpp`.

use std::io::Write as _;
use crate::frame::{Atoms, Frame};
use crate::xdr::{Result, XdrError};
use crate::{gro, pdb, trr, xtc};

/// A trajectory frame reduced to the coordinates of the selected atoms.
///
/// This is what [`CoordinateReader`] yields; it is deliberately plain data so
/// that callers can use it without knowing anything about the file formats.
#[derive(Debug, Clone, Default)]
pub struct CoordFrame {
    /// Index of the frame in the input file (0 based).
    pub frame: usize,
    pub step: Option<i64>,
    pub time: Option<f64>,
    pub boxm: Option<crate::frame::Matrix>,
    /// Coordinates of the selected atoms, in selection order.
    pub x: Vec<crate::frame::Rvec>,
    /// Velocities of the selected atoms, when the file contains them.
    pub v: Option<Vec<crate::frame::Rvec>>,
    /// Forces of the selected atoms, when the file contains them.
    pub f: Option<Vec<crate::frame::Rvec>>,
}

/// Which frames [`CoordinateReader`] should return.
#[derive(Debug, Clone, Copy)]
pub struct FrameRange {
    /// First time to return (frames before it are skipped).
    pub begin: Option<f64>,
    /// Last time to return; reading stops after it.
    pub end: Option<f64>,
    /// Return every `skip`-th frame of the input (1 = all).
    pub skip: usize,
    /// Stop after this many returned frames.
    pub max_frames: Option<usize>,
}

impl Default for FrameRange {
    fn default() -> Self {
        FrameRange {
            begin: None,
            end: None,
            skip: 1,
            max_frames: None,
        }
    }
}

/// Streams coordinates out of a trajectory.
///
/// ```no_run
/// use gmx_rs_tools::trx::{CoordinateReader, FrameRange};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut reader = CoordinateReader::open("traj.xtc", None, FrameRange::default())?.1;
/// while let Some(frame) = reader.next_frame()? {
///     println!("t = {:?} with {} atoms", frame.time, frame.x.len());
/// }
/// # Ok(())
/// # }
/// ```
pub struct CoordinateReader {
    source: FrameSource,
    /// Atom indices to return; `None` means every atom in the file.
    selection: Option<Vec<usize>>,
    range: FrameRange,
    /// Index of the next frame in the input file.
    next_frame: usize,
    /// Number of frames already returned.
    returned: usize,
    /// Number of frames seen since the first one at or after `begin`; `-skip`
    /// counts from there (this is what `gmx trjconv` does, since its frame
    /// counter starts after the begin-time skip).
    since_begin: usize,
}

impl CoordinateReader {
    /// Opens a trajectory; returns the detected format and the reader.
    pub fn open(
        path: &str,
        selection: Option<&[usize]>,
        range: FrameRange,
    ) -> Result<(TrxFormat, CoordinateReader)> {
        let (format, source) = FrameSource::open(path)?;
        let reader = CoordinateReader {
            source,
            selection: selection.map(|s| s.to_vec()),
            range,
            next_frame: 0,
            returned: 0,
            since_begin: 0,
        };
        Ok((format, reader))
    }

    /// Returns the next requested frame, or `None` at the end of the file.
    pub fn next_frame(&mut self) -> Result<Option<CoordFrame>> {
        if let Some(max) = self.range.max_frames {
            if self.returned >= max {
                return Ok(None);
            }
        }
        let skip = self.range.skip.max(1);
        loop {
            let Some(frame) = self.source.next_frame()? else {
                return Ok(None);
            };
            let index = self.next_frame;
            self.next_frame += 1;
            let time = frame.time;
            if let Some(begin) = self.range.begin {
                if time.unwrap_or(0.0) < begin {
                    continue;
                }
            }
            if let Some(end) = self.range.end {
                if time.unwrap_or(0.0) > end {
                    return Ok(None);
                }
            }
            let k = self.since_begin;
            self.since_begin += 1;
            if k % skip != 0 {
                continue;
            }
            let pick = |v: &Vec<crate::frame::Rvec>| -> Vec<crate::frame::Rvec> {
                match &self.selection {
                    Some(sel) => sel.iter().filter(|i| **i < v.len()).map(|i| v[*i]).collect(),
                    None => v.clone(),
                }
            };
            let out = CoordFrame {
                frame: index,
                step: frame.step,
                time,
                boxm: frame.boxm,
                x: frame.x.as_ref().map(pick).unwrap_or_default(),
                v: frame.v.as_ref().map(pick),
                f: frame.f.as_ref().map(pick),
            };
            self.returned += 1;
            return Ok(Some(out));
        }
    }
}

/// Convenience wrapper: reads every requested frame into memory.
pub fn read_coordinates(
    path: &str,
    selection: Option<&[usize]>,
    range: FrameRange,
) -> Result<Vec<CoordFrame>> {
    let mut reader = CoordinateReader::open(path, selection, range)?.1;
    let mut out = Vec::new();
    while let Some(frame) = reader.next_frame()? {
        out.push(frame);
    }
    Ok(out)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrxFormat {
    Xtc,
    Trr,
    Gro,
    Pdb,
}

impl TrxFormat {
    pub fn extension(self) -> &'static str {
        match self {
            TrxFormat::Xtc => "xtc",
            TrxFormat::Trr => "trr",
            TrxFormat::Gro => "gro",
            TrxFormat::Pdb => "pdb",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            TrxFormat::Xtc => "Compressed trajectory (portable xdr format)",
            TrxFormat::Trr => "Trajectory in portable xdr format",
            TrxFormat::Gro => "Coordinate file in Gromos-87 format",
            TrxFormat::Pdb => "Protein data bank file",
        }
    }

    /// True when the format can carry velocities.
    pub fn can_have_velocities(self) -> bool {
        matches!(self, TrxFormat::Trr | TrxFormat::Gro | TrxFormat::Pdb)
    }
}

/// Detects a trajectory/structure format from the file name extension,
/// mirroring `fn2ftp()`.
pub fn format_from_path(path: &str) -> Option<TrxFormat> {
    let lower = path.to_ascii_lowercase();
    let ext = lower.rsplit('.').next()?;
    match ext {
        "xtc" => Some(TrxFormat::Xtc),
        "trr" => Some(TrxFormat::Trr),
        "gro" => Some(TrxFormat::Gro),
        "pdb" | "ent" => Some(TrxFormat::Pdb),
        _ => None,
    }
}

/// Reads every frame of a trajectory or structure file.
pub fn read_frames(path: &str) -> Result<(TrxFormat, Vec<Frame>)> {
    let (format, mut source) = FrameSource::open(path)?;
    let mut frames = Vec::new();
    while let Some(f) = source.next_frame()? {
        frames.push(f);
    }
    Ok((format, frames))
}

/// Streams frames one at a time so that trajectories much larger than memory
/// can be processed.
pub enum FrameSource {
    Xtc(std::io::BufReader<std::fs::File>),
    Trr(std::io::BufReader<std::fs::File>),
    Memory(std::vec::IntoIter<Frame>),
}

impl FrameSource {
    /// Opens a trajectory or structure file and detects its format.
    pub fn open(path: &str) -> Result<(TrxFormat, FrameSource)> {
        let format = format_from_path(path).ok_or_else(|| {
            XdrError::Invalid(format!(
                "File {path} is not a supported trajectory or structure file"
            ))
        })?;
        let source = match format {
            TrxFormat::Xtc => FrameSource::Xtc(std::io::BufReader::with_capacity(
                1 << 20,
                std::fs::File::open(path)
                    .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?,
            )),
            TrxFormat::Trr => FrameSource::Trr(std::io::BufReader::with_capacity(
                1 << 20,
                std::fs::File::open(path)
                    .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?,
            )),
            TrxFormat::Gro => FrameSource::Memory(gro::read_all(path)?.into_iter()),
            TrxFormat::Pdb => FrameSource::Memory(pdb::read_all(path)?.into_iter()),
        };
        Ok((format, source))
    }

    /// Returns the next frame, or `None` at the end of the file.
    pub fn next_frame(&mut self) -> Result<Option<Frame>> {
        match self {
            FrameSource::Xtc(r) => match xtc::read_frame_bytes(r)? {
                Some(bytes) => {
                    let mut rd = crate::xdr::Reader::new(&bytes);
                    Ok(xtc::read_frame(&mut rd)?)
                }
                None => Ok(None),
            },
            FrameSource::Trr(r) => match trr::read_frame_bytes(r)? {
                Some(bytes) => Ok(Some(trr::decode_frame(&bytes)?)),
                None => Ok(None),
            },
            FrameSource::Memory(it) => Ok(it.next()),
        }
    }
}

/// Reads every frame of a trajectory or structure file (eager version, used by
/// the tests and by callers that need random access).
pub fn read_frames_eager(path: &str) -> Result<(TrxFormat, Vec<Frame>)> {
    let format = format_from_path(path).ok_or_else(|| {
        XdrError::Invalid(format!(
            "File {path} is not a supported trajectory or structure file"
        ))
    })?;
    match format {
        TrxFormat::Xtc => {
            let data = std::fs::read(path)
                .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
            let mut r = crate::xdr::Reader::new(&data);
            let mut frames = Vec::new();
            while let Some(f) = xtc::read_frame(&mut r)? {
                frames.push(f);
            }
            Ok((format, frames))
        }
        TrxFormat::Trr => {
            let data = std::fs::read(path)
                .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
            let mut r = crate::xdr::Reader::new(&data);
            let mut frames = Vec::new();
            while let Some(f) = trr::read_frame(&mut r)? {
                frames.push(f);
            }
            Ok((format, frames))
        }
        TrxFormat::Gro => Ok((format, gro::read_all(path)?)),
        TrxFormat::Pdb => Ok((format, pdb::read_all(path)?)),
    }
}

/// Reads only the first frame (used to obtain the atom count).
pub fn read_first_frame(path: &str) -> Result<Frame> {
    let (_, frames) = read_frames(path)?;
    frames
        .into_iter()
        .next()
        .ok_or_else(|| XdrError::Invalid(format!("could not read a frame from {path}")))
}

/// Streaming writer for trajectory output.
pub enum TrxWriter {
    Xtc {
        path: String,
        out: std::io::BufWriter<std::fs::File>,
        prec: f32,
    },
    Trr {
        path: String,
        out: std::io::BufWriter<std::fs::File>,
    },
    Gro {
        path: String,
        out: std::io::BufWriter<std::fs::File>,
        prefixes: Vec<Vec<u8>>,
    },
    Pdb {
        path: String,
        out: std::io::BufWriter<std::fs::File>,
        prefixes: Vec<Vec<u8>>,
        suffixes: Vec<Vec<u8>>,
        model: i32,
    },
}

/// Opens the output file with a buffer large enough for a whole frame.
fn create_output(path: &str) -> Result<std::io::BufWriter<std::fs::File>> {
    let f = std::fs::File::create(path)
        .map_err(|e| XdrError::Invalid(format!("cannot write {path}: {e}")))?;
    Ok(std::io::BufWriter::with_capacity(1 << 20, f))
}

impl TrxWriter {
    pub fn create(path: &str, format: TrxFormat, prec: f32) -> Result<TrxWriter> {
        match format {
            TrxFormat::Xtc => Ok(TrxWriter::Xtc {
                path: path.to_string(),
                out: create_output(path)?,
                prec,
            }),
            TrxFormat::Trr => Ok(TrxWriter::Trr {
                path: path.to_string(),
                out: create_output(path)?,
            }),
            TrxFormat::Gro => Ok(TrxWriter::Gro {
                path: path.to_string(),
                out: create_output(path)?,
                prefixes: Vec::new(),
            }),
            TrxFormat::Pdb => Ok(TrxWriter::Pdb {
                path: path.to_string(),
                out: create_output(path)?,
                prefixes: Vec::new(),
                suffixes: Vec::new(),
                model: 0,
            }),
        }
    }

    /// Sets the topology used for text output formats.
    pub fn set_atoms(&mut self, atoms: &Atoms, index: &[usize]) {
        match self {
            TrxWriter::Gro { prefixes, .. } => *prefixes = gro::atom_prefixes(atoms, index),
            TrxWriter::Pdb {
                prefixes,
                suffixes,
                ..
            } => {
                let (p, s) = pdb::atom_fields(atoms, index);
                *prefixes = p;
                *suffixes = s;
            }
            _ => {}
        }
    }

    pub fn write_frame(&mut self, frame: &Frame, index: &[usize], title: &str) -> Result<()> {
        match self {
            TrxWriter::Xtc { out, prec, .. } => {
                let mut w = crate::xdr::Writer::new();
                xtc::write_frame(&mut w, frame, *prec);
                out.write_all(&w.data)
                    .map_err(|e| XdrError::Invalid(format!("write error: {e}")))?;
                Ok(())
            }
            TrxWriter::Trr { out, .. } => {
                let mut w = crate::xdr::Writer::new();
                trr::write_frame(&mut w, frame);
                out.write_all(&w.data)
                    .map_err(|e| XdrError::Invalid(format!("write error: {e}")))?;
                Ok(())
            }
            TrxWriter::Gro { out, prefixes, .. } => {
                let x = frame.x.clone().unwrap_or_default();
                gro::write_frame(
                    out,
                    title,
                    prefixes,
                    &x,
                    index,
                    frame.v.as_deref(),
                    &frame.boxm.unwrap_or([[0.0; 3]; 3]),
                );
                Ok(())
            }
            TrxWriter::Pdb {
                out,
                prefixes,
                suffixes,
                model,
                ..
            } => {
                let x = frame.x.clone().unwrap_or_default();
                *model += 1;
                pdb::write_frame(
                    out,
                    title,
                    prefixes,
                    suffixes,
                    &x,
                    index,
                    frame.pbc_type,
                    &frame.boxm.unwrap_or([[0.0; 3]; 3]),
                    *model,
                );
                Ok(())
            }
        }
    }

    /// Writes the accumulated output to `path` (binary formats use the path
    /// provided at creation time).
    pub fn finish(self) -> Result<()> {
        match self {
            TrxWriter::Xtc { path, mut out, .. }
            | TrxWriter::Trr { path, mut out }
            | TrxWriter::Gro { path, mut out, .. }
            | TrxWriter::Pdb { path, mut out, .. } => {
                out.flush()
                    .map_err(|e| XdrError::Invalid(format!("cannot write {path}: {e}")))
            }
        }
    }
}
