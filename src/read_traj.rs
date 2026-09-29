//! Reads the processed trajectory that the MM-PBSA calculation runs on.
//!
//! Coordinates used to be decoded with the `xdrfile` crate, i.e. the C library
//! bundled with the GROMACS sources.  They are now read through the vendored
//! `gmx-rs-tools` reader, the same streaming reader that `trjconv` and
//! `dump -f` use, so trajectories are handled by the code that is already
//! linked into the binary (and the crate's C dependency is gone).
//!
//! Frames are converted to double precision as they stream in, so the
//! single-precision frames of the reader never accumulate: peak memory is one
//! f64 copy of the trajectory instead of an f32 copy plus two more
//! conversions downstream.

use gmx_rs_tools::progress::Progress;
use gmx_rs_tools::trx::{CoordinateReader, FrameRange};
use ndarray::Array3;

/// Reads every frame of the trajectory.
///
/// Returns the frame times (ps) and the coordinates as one
/// `(n_frames, n_atoms, 3)` array in GROMACS length units (nm).
pub fn read_traj(trj: &str) -> (Vec<f64>, Array3<f64>) {
    let (_, mut reader) = CoordinateReader::open(trj, None, FrameRange::default())
        .unwrap_or_else(|e| panic!("Cannot open trajectory {}: {}", trj, e));

    let mut progress = Progress::new();
    let mut times: Vec<f64> = Vec::new();
    let mut coords: Vec<f64> = Vec::new();
    let mut n_atoms: Option<usize> = None;
    loop {
        let frame = match reader.next_frame() {
            Ok(Some(frame)) => frame,
            Ok(None) => break,
            Err(e) => panic!("Cannot read trajectory {}: {}", trj, e),
        };
        match n_atoms {
            None => n_atoms = Some(frame.x.len()),
            Some(n) if n != frame.x.len() => panic!(
                "Trajectory {}: frame {} has {} atoms, expected {}",
                trj,
                times.len(),
                frame.x.len(),
                n
            ),
            Some(_) => {}
        }
        let time = frame.time.unwrap_or(0.0);
        times.push(time);
        coords.extend(frame.x.iter().flat_map(|c| c.iter().map(|&v| v as f64)));
        let frames = times.len() as u64;
        progress.update_from(reader.read_progress(), frames, time);
    }
    progress.pause();

    let n_atoms = n_atoms.unwrap_or(0);
    let coordinates = Array3::from_shape_vec((times.len(), n_atoms, 3), coords)
        .expect("Trajectory frames do not have a consistent atom count");
    println!("Finished reading trajectory with {} frames.", times.len());

    (times, coordinates)
}
