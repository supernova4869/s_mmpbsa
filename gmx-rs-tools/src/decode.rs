//! Parallel decoding of trajectory frames.
//!
//! `trjconv`, `dump -f` and `coords` stream frames, which means that a frame
//! has to be read, decoded, processed and written before the next one is
//! read.  The reading is I/O bound, but the compressed XTC (and TRR) codec is
//! pure CPU work, and on a large trajectory it is a large part of the runtime.
//!
//! [`ParallelSource`] keeps the frames in order while moving that CPU work off
//! the critical path: one thread reads the raw frames (so the disk is busy
//! while the pool decodes), and the frames of a batch are decoded in parallel.

use std::io::{BufReader, Seek};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::Arc;

use rayon::prelude::*;

use crate::frame::Frame;
use crate::trx::{ReadProgress, TrxFormat};
use crate::xdr::{Result, XdrError};
use crate::{trr, xtc};

/// How many raw frames the reader thread hands to the pool at once.
const DECODE_BATCH: usize = 16;
/// How many decoded frames the reader thread may run ahead by.
///
/// This bounds the memory of the pipeline: at most `DECODE_DEPTH` frames plus
/// one batch of raw frame bytes are in flight.
const DECODE_DEPTH: usize = 32;

/// A trajectory whose frames are decoded by a worker pool.
pub struct ParallelSource {
    frames: Receiver<Result<Frame>>,
    /// Bytes read by the reader thread, for the progress bars.
    position: Arc<AtomicU64>,
    size: u64,
    /// Kept so that the reader stops when the source is dropped.
    thread: Option<std::thread::JoinHandle<()>>,
}

impl ParallelSource {
    /// Starts the reader thread for an XTC or TRR file.
    pub fn open(path: &str) -> Result<ParallelSource> {
        let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        let format = crate::trx::format_from_path(path)
            .ok_or_else(|| XdrError::Invalid(format!("File {path} is not a trajectory")))?;
        let (tx, rx) = sync_channel(DECODE_DEPTH);
        let position = Arc::new(AtomicU64::new(0));
        let reader_position = Arc::clone(&position);
        let reader_path = path.to_string();
        let thread = std::thread::Builder::new()
            .name("gmx-rs-tools-reader".to_string())
            .spawn(move || read_frames(&reader_path, format, tx, reader_position))
            .map_err(|e| XdrError::Invalid(format!("cannot start the reader thread: {e}")))?;
        Ok(ParallelSource {
            frames: rx,
            position,
            size,
            thread: Some(thread),
        })
    }

    /// Returns the next frame of the input, in input order.
    pub fn next_frame(&mut self) -> Result<Option<Frame>> {
        match self.frames.recv() {
            Ok(frame) => frame.map(Some),
            // The reader thread has finished, so this is a clean end of file.
            Err(_) => Ok(None),
        }
    }

    pub fn read_progress(&mut self) -> ReadProgress {
        ReadProgress::Bytes(self.position.load(Ordering::Relaxed), self.size)
    }
}

impl Drop for ParallelSource {
    fn drop(&mut self) {
        // The reader thread blocks in `send` while the queue is full, and a
        // tool that stops early (`-e`, `-dump`) never drains it, so the
        // receiving end has to be closed before joining.
        let (_tx, closed) = std::sync::mpsc::channel();
        drop(std::mem::replace(&mut self.frames, closed));
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

/// Reads raw frames and decodes them in parallel, sending them in order.
fn read_frames(
    path: &str,
    format: TrxFormat,
    tx: SyncSender<Result<Frame>>,
    position: Arc<AtomicU64>,
) {
    let file = match std::fs::File::open(Path::new(path)) {
        Ok(f) => f,
        Err(e) => {
            let _ = tx.send(Err(cannot_read(path, e)));
            return;
        }
    };
    let mut r = BufReader::with_capacity(1 << 20, file);
    let read_raw: fn(&mut BufReader<std::fs::File>) -> Result<Option<Vec<u8>>> = match format {
        TrxFormat::Xtc => xtc::read_frame_bytes,
        TrxFormat::Trr => trr::read_frame_bytes,
        _ => unreachable!("only streamed formats use the parallel reader"),
    };
    let decode: fn(&[u8]) -> Result<Frame> = match format {
        TrxFormat::Xtc => decode_xtc,
        TrxFormat::Trr => trr::decode_frame,
        _ => unreachable!("only streamed formats use the parallel reader"),
    };

    loop {
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(DECODE_BATCH);
        let mut eof = false;
        for _ in 0..DECODE_BATCH {
            match read_raw(&mut r) {
                Ok(Some(raw)) => batch.push(raw),
                Ok(None) => {
                    eof = true;
                    break;
                }
                Err(e) => {
                    let _ = tx.send(Err(e));
                    return;
                }
            }
        }
        if let Ok(pos) = r.stream_position() {
            position.store(pos, Ordering::Relaxed);
        }
        if batch.len() > 1 {
            for frame in batch.par_iter().map(|raw| decode(raw)).collect::<Vec<_>>() {
                if tx.send(frame).is_err() {
                    return;
                }
            }
        } else if let Some(raw) = batch.first() {
            if tx.send(decode(raw)).is_err() {
                return;
            }
        }
        if eof {
            return;
        }
    }
}

fn decode_xtc(raw: &[u8]) -> Result<Frame> {
    let mut r = crate::xdr::Reader::new(raw);
    xtc::read_frame(&mut r)?.ok_or_else(|| XdrError::Truncated("xtc frame"))
}

fn cannot_read(path: &str, e: std::io::Error) -> XdrError {
    XdrError::Invalid(format!("cannot read {path}: {e}"))
}
