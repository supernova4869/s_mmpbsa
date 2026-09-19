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

    /// Position in the input, forwarded to [`FrameSource::read_progress`].
    pub fn read_progress(&mut self) -> Option<ReadProgress> {
        self.source.read_progress()
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

/// Number of frames in a trajectory or structure file.
///
/// `xtc` and `trr` files are counted by walking their frame headers, which
/// skips the payload of every frame; structure files are read.  The count is
/// what the progress bars use as their total.
pub fn frame_count(path: &str) -> Result<usize> {
    match format_from_path(path) {
        Some(TrxFormat::Xtc) => xtc::frame_count(path),
        Some(TrxFormat::Trr) => trr::frame_count(path),
        Some(TrxFormat::Gro) => Ok(gro::read_all(path)?.len()),
        Some(TrxFormat::Pdb) => Ok(pdb::read_all(path)?.len()),
        None => Err(XdrError::Invalid(format!(
            "File {path} is not a supported trajectory or structure file"
        ))),
    }
}

/// Streams frames one at a time so that trajectories much larger than memory
/// can be processed.
pub enum FrameSource {
    Xtc(std::io::BufReader<std::fs::File>),
    Trr(std::io::BufReader<std::fs::File>),
    Memory {
        iter: std::vec::IntoIter<Frame>,
        total: usize,
    },
}

/// How far a [`FrameSource`] has progressed through its input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadProgress {
    /// Bytes consumed and the size of the input file.
    Bytes(u64, u64),
    /// Frames read and the number of frames in the file.
    Frames(u64, u64),
}

impl ReadProgress {
    /// Amount of the input consumed so far.
    pub fn done(self) -> u64 {
        match self {
            ReadProgress::Bytes(done, _) => done,
            ReadProgress::Frames(done, _) => done,
        }
    }

    /// Number of frames in the input, when the progress counts frames rather
    /// than bytes of it.
    pub fn known_total(self) -> Option<u64> {
        match self {
            ReadProgress::Bytes(..) => None,
            ReadProgress::Frames(_, total) => Some(total),
        }
    }

    /// Fraction of the input that has been consumed, when it is known.
    pub fn fraction(self) -> Option<f64> {
        let (done, total) = match self {
            ReadProgress::Bytes(done, total) => (done, total),
            ReadProgress::Frames(done, total) => (done, total),
        };
        if total == 0 {
            None
        } else {
            Some(done as f64 / total as f64)
        }
    }
}

impl FrameSource {
    /// Opens a trajectory or structure file and detects its format.
    pub fn open(path: &str) -> Result<(TrxFormat, FrameSource)> {
        let format = format_from_path(path).ok_or_else(|| {
            XdrError::Invalid(format!(
                "File {path} is not a supported trajectory or structure file"
            ))
        })?;
        let open_reader = || -> Result<std::io::BufReader<std::fs::File>> {
            Ok(std::io::BufReader::with_capacity(
                1 << 20,
                std::fs::File::open(path)
                    .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?,
            ))
        };
        let source = match format {
            TrxFormat::Xtc => FrameSource::Xtc(open_reader()?),
            TrxFormat::Trr => FrameSource::Trr(open_reader()?),
            TrxFormat::Gro => {
                let frames = gro::read_all(path)?;
                let total = frames.len();
                FrameSource::Memory {
                    iter: frames.into_iter(),
                    total,
                }
            }
            TrxFormat::Pdb => {
                let frames = pdb::read_all(path)?;
                let total = frames.len();
                FrameSource::Memory {
                    iter: frames.into_iter(),
                    total,
                }
            }
        };
        Ok((format, source))
    }

    /// Position in the input, for the progress bar.
    ///
    /// Binary trajectories report the number of bytes consumed of the file
    /// size, files that are read into memory report the frame count.
    pub fn read_progress(&mut self) -> Option<ReadProgress> {
        match self {
            FrameSource::Xtc(reader) | FrameSource::Trr(reader) => {
                use std::io::Seek as _;
                let total = reader.get_ref().metadata().ok()?.len();
                let consumed = reader.stream_position().ok()?;
                Some(ReadProgress::Bytes(consumed, total))
            }
            FrameSource::Memory { iter, total } => Some(ReadProgress::Frames(
                (*total - iter.len()) as u64,
                *total as u64,
            )),
        }
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
            FrameSource::Memory { iter, .. } => Ok(iter.next()),
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
/// A frame ready for serial encoding by [`TrxWriter`].
enum Pending {
    Xtc { frame: Frame, prec: f32 },
    Trr { frame: Frame },
    Gro { frame: Frame, title: String },
    Pdb { frame: Frame, title: String, model: i32 },
}

/// Writes a trajectory or a structure file.
///
/// Each frame is encoded and written immediately on the calling thread.
pub struct TrxWriter {
    path: String,
    format: TrxFormat,
    out: std::io::BufWriter<std::fs::File>,
    /// Precision (multiplication factor) of the XTC output.
    prec: f32,
    /// Atoms written to the output.
    index: Vec<usize>,
    /// Precomputed per atom fields of the text output formats.
    prefixes: Vec<Vec<u8>>,
    suffixes: Vec<Vec<u8>>,
    /// Model number of the next PDB frame.
    model: i32,
}

/// Opens the output file with a buffer large enough for a whole frame.
fn create_output(path: &str) -> Result<std::io::BufWriter<std::fs::File>> {
    let f = std::fs::File::create(path)
        .map_err(|e| XdrError::Invalid(format!("cannot write {path}: {e}")))?;
    Ok(std::io::BufWriter::with_capacity(1 << 20, f))
}

/// Selects the requested atoms from a frame before binary output.
fn select_frame(frame: &Frame, index: &[usize]) -> Frame {
    // An empty index means "write every atom"; avoid a needless copy.
    if index.is_empty() {
        return frame.clone();
    }
    let select = |data: &Option<Vec<crate::frame::Rvec>>| {
        data.as_ref().map(|values| {
            index.iter()
                .filter(|&&i| i < values.len())
                .map(|&i| values[i])
                .collect::<Vec<_>>()
        })
    };
    Frame {
        natoms: index.len(),
        x: select(&frame.x),
        v: select(&frame.v),
        f: select(&frame.f),
        ..frame.clone()
    }
}

impl TrxWriter {
    pub fn create(path: &str, format: TrxFormat, prec: f32) -> Result<TrxWriter> {
        Ok(TrxWriter {
            path: path.to_string(),
            format,
            out: create_output(path)?,
            prec,
            index: Vec::new(),
            prefixes: Vec::new(),
            suffixes: Vec::new(),
            model: 0,
        })
    }

    /// Sets the topology used for text output formats and the atoms written.
    pub fn set_atoms(&mut self, atoms: &Atoms, index: &[usize]) {
        match self.format {
            TrxFormat::Gro => self.prefixes = gro::atom_prefixes(atoms, index),
            TrxFormat::Pdb => {
                let (p, s) = pdb::atom_fields(atoms, index);
                self.prefixes = p;
                self.suffixes = s;
            }
            _ => {}
        }
        self.index = index.to_vec();
    }

    /// Sets the atoms written to the output; needed when there is no topology
    /// to build the fields of the text formats from (`set_atoms` does it for
    /// the formats that need one).
    pub fn set_index(&mut self, index: &[usize]) {
        self.index = index.to_vec();
    }

    /// Sets the precision of the XTC output for the frames written from now on.
    pub fn set_precision(&mut self, prec: f32) {
        self.prec = prec;
    }

    /// Encodes and writes one frame immediately.
    pub fn write_frame(&mut self, frame: Frame, title: &str) -> Result<()> {
        let job = match self.format {
            TrxFormat::Xtc => Pending::Xtc {
                frame,
                prec: self.prec,
            },
            TrxFormat::Trr => Pending::Trr { frame },
            TrxFormat::Gro => Pending::Gro {
                frame,
                title: title.to_string(),
            },
            TrxFormat::Pdb => {
                self.model += 1;
                Pending::Pdb {
                    frame,
                    title: title.to_string(),
                    model: self.model,
                }
            }
        };
        let bytes = self.encode(&job)?;
        self.out.write_all(&bytes).map_err(|e| {
            XdrError::Invalid(format!("cannot write {}: {e}", self.path))
        })?;
        Ok(())
    }

    /// Encodes one frame.
    fn encode(&self, job: &Pending) -> Result<Vec<u8>> {
        match job {
            Pending::Xtc { frame, prec } => {
                let frame = select_frame(frame, &self.index);
                let mut w = crate::xdr::Writer::new();
                xtc::write_frame(&mut w, &frame, *prec);
                Ok(w.data)
            }
            Pending::Trr { frame } => {
                let frame = select_frame(frame, &self.index);
                let mut w = crate::xdr::Writer::new();
                trr::write_frame(&mut w, &frame);
                Ok(w.data)
            }
            Pending::Gro { frame, title } => {
                let mut buf: Vec<u8> = Vec::new();
                gro::write_frame(
                    &mut buf,
                    title,
                    &self.prefixes,
                    frame.x.as_deref().unwrap_or(&[]),
                    &self.index,
                    frame.v.as_deref(),
                    &frame.boxm.unwrap_or([[0.0; 3]; 3]),
                );
                Ok(buf)
            }
            Pending::Pdb {
                frame,
                title,
                model,
            } => {
                let mut buf: Vec<u8> = Vec::new();
                pdb::write_frame(
                    &mut buf,
                    title,
                    &self.prefixes,
                    &self.suffixes,
                    frame.x.as_deref().unwrap_or(&[]),
                    &self.index,
                    frame.pbc_type,
                    &frame.boxm.unwrap_or([[0.0; 3]; 3]),
                    *model,
                );
                Ok(buf)
            }
        }
    }

    /// Writes the accumulated output to `path` (binary formats use the path
    /// provided at creation time).
    pub fn finish(mut self) -> Result<()> {
        self.out
            .flush()
            .map_err(|e| XdrError::Invalid(format!("cannot write {}: {e}", self.path)))
    }
}
