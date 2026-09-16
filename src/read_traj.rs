//! Reads the processed trajectory that the MM-PBSA calculation runs on.
//!
//! Coordinates used to be decoded with the `xdrfile` crate, i.e. the C library
//! bundled with the GROMACS sources.  They are now read through the vendored
//! `gmx-rs-tools` reader, the same streaming reader that `trjconv` and
//! `dump -f` use, so trajectories are handled by the code that is already
//! linked into the binary (and the crate's C dependency is gone).
//!
//! `read_traj` returns the time and the coordinates of every frame.  The
//! processed trajectory is always XTC, but the reader accepts the other
//! GROMACS trajectory and structure formats (`trr`, `gro`, `pdb`) as well.

use gmx_rs_tools::progress::Progress;
use gmx_rs_tools::trx::{self, CoordinateReader, FrameRange};

pub fn read_traj(trj: &str) -> Vec<(f64, Vec<[f32; 3]>)> {
    let (_, mut reader) = CoordinateReader::open(trj, None, FrameRange::default())
        .unwrap_or_else(|e| panic!("Cannot open trajectory {}: {}", trj, e));

    let mut progress = Progress::new();
    // The bar counts frames, so the number of frames of the trajectory is
    // needed before reading it.  Counting walks the frame headers of the file,
    // which is only worth it when the bar is actually drawn.
    let total = if progress.is_enabled() {
        trx::frame_count(trj).unwrap_or(0) as u64
    } else {
        0
    };
    let mut frame_data: Vec<(f64, Vec<[f32; 3]>)> = Vec::new();
    loop {
        let frame = match reader.next_frame() {
            Ok(Some(frame)) => frame,
            Ok(None) => break,
            Err(e) => panic!("Cannot read trajectory {}: {}", trj, e),
        };
        let time = frame.time.unwrap_or(0.0);
        let frames = frame_data.len() as u64 + 1;
        progress.update(frames, total, time);
        frame_data.push((time, frame.x));
    }
    progress.pause();

    println!("Finished reading trajectory with {} frames.", frame_data.len());

    frame_data
}
