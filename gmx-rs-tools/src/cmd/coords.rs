//! `gmx-rs-tools coords` — read coordinates straight out of a trajectory.
//!
//! This is a thin front end over the trajectory readers used by `gmx-rs-tools dump`
//! and `gmx-rs-tools trjconv`: frames are streamed, an index group can be used to
//! select atoms, and a time window / frame stride can be applied, which makes
//! the output usable for very large trajectories.

use std::io::Write as _;

use crate::cmd::{select_group, Args};
use crate::index;
use crate::tpr;
use crate::trx::{self, FrameRange};

/// The frame progress of a running extraction.
///
/// Progress is measured against the requested time window when one was given
/// (`-b`/`-e`), and against the position in the input file otherwise.
struct FrameProgress {
    bar: crate::progress::Progress,
    time_window: Option<(f64, f64)>,
}

impl FrameProgress {
    fn new(time_window: Option<(f64, f64)>) -> Self {
        FrameProgress {
            bar: crate::progress::Progress::new(),
            time_window,
        }
    }

    /// Reports the frame that has just been read.
    fn tick(&mut self, reader: &mut trx::CoordinateReader, frames: u64, time: f64) {
        let fraction = match self.time_window {
            Some((start, end)) if end > start => {
                Some(((time - start) / (end - start)).clamp(0.0, 1.0))
            }
            _ => reader.read_progress().and_then(|p| p.fraction()),
        };
        self.bar.update(fraction, frames, time);
    }

    fn finish(&mut self) {
        self.bar.finish();
    }
}

pub fn run(argv: Vec<String>) -> i32 {
    let args = Args::parse(argv);
    let in_file = match args.get("f") {
        Some(f) => f,
        None => {
            eprintln!("gmx-rs-tools coords: option -f is required");
            return 1;
        }
    };
    let range = FrameRange {
        begin: args.get("b").and_then(|v| v.parse().ok()),
        end: args.get("e").and_then(|v| v.parse().ok()),
        skip: args.int("skip", 1).max(1) as usize,
        max_frames: args.get("nmax").and_then(|v| v.parse().ok()),
    };
    let ndec = args.int("ndec", 6).max(0) as usize;
    let with_velocities = args.flag("vel");

    // Atom selection.  A group can be picked from an index file (-n) or from
    // the default groups of a topology (-s); without either, every atom is
    // written.
    let mut selection: Option<Vec<usize>> = None;
    let need_groups = args.has("n") || args.has("s") || args.has("sel");
    if need_groups {
        let atoms = match args.get("s") {
            Some(s) => match tpr::TprFile::read(&s) {
                Ok(t) => match tpr::parse_body(&t.header, &t.body) {
                    Ok(body) => body
                        .mtop
                        .map(|m| m.global_atoms())
                        .unwrap_or_default(),
                    Err(e) => {
                        eprintln!("{e}");
                        return 1;
                    }
                },
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            },
            None => Default::default(),
        };
        let groups = match args.get("n") {
            Some(f) => match index::read_ndx(&f) {
                Ok(g) => g,
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            },
            None => {
                if atoms.nr() == 0 {
                    eprintln!("No index file specified and no topology available for default groups");
                    return 1;
                }
                index::analyse(&atoms, false)
            }
        };
        // `-sel` takes the group directly (script friendly); otherwise the
        // group is asked for on stdin like the other tools.
        let g = match args.get("sel") {
            Some(spec) => {
                let g = match spec.parse::<usize>() {
                    Ok(v) if v < groups.len() => v as i32,
                    _ => index::find_group(&spec, &groups),
                };
                if g < 0 {
                    eprintln!("No such group '{spec}'");
                    return 1;
                }
                g as usize
            }
            None => match select_group(&groups, "Select a group:") {
                Some(g) => g,
                None => return 1,
            },
        };
        selection = Some(groups[g].particle_indices.clone());
    }

    let out: Box<dyn std::io::Write> = match args.get("o") {
        Some(path) => match std::fs::File::create(&path) {
            Ok(f) => Box::new(std::io::BufWriter::with_capacity(1 << 20, f)),
            Err(e) => {
                eprintln!("cannot write {path}: {e}");
                return 1;
            }
        },
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    let mut out = out;

    let mut reader = match trx::CoordinateReader::open(&in_file, selection.as_deref(), range) {
        Ok((_, r)) => r,
        Err(e) => {
            eprintln!("{e}");
            return 1;
        }
    };
    let header = if with_velocities {
        "# frame time step index x y z vx vy vz"
    } else {
        "# frame time step index x y z"
    };
    if writeln!(out, "{header}").is_err() {
        return 1;
    }
    let mut nframes = 0usize;
    let mut natoms_out = 0usize;
    // With `-e` the bar follows the requested time window, otherwise how much
    // of the input file has been read.
    let mut progress = FrameProgress::new(match (range.begin, range.end) {
        (_, Some(end)) => Some((range.begin.unwrap_or(0.0), end)),
        _ => None,
    });
    loop {
        let frame = match reader.next_frame() {
            Ok(Some(f)) => f,
            Ok(None) => break,
            Err(e) => {
                eprintln!("{e}");
                return 1;
            }
        };
        progress.tick(&mut reader, nframes as u64 + 1, frame.time.unwrap_or(0.0));
        let time = frame.time.unwrap_or(0.0);
        let step = frame.step.unwrap_or(0);
        if natoms_out == 0 {
            natoms_out = frame.x.len();
        }
        for (i, x) in frame.x.iter().enumerate() {
            let index = selection
                .as_ref()
                .and_then(|s| s.get(i))
                .map(|v| v + 1)
                .unwrap_or(i + 1);
            let line = if with_velocities {
                let v = frame
                    .v
                    .as_ref()
                    .and_then(|v| v.get(i))
                    .copied()
                    .unwrap_or([0.0; 3]);
                format!(
                    "{} {time:.6} {step} {index} {:.ndec$} {:.ndec$} {:.ndec$} {:.ndec$} {:.ndec$} {:.ndec$}\n",
                    frame.frame,
                    x[0] as f64,
                    x[1] as f64,
                    x[2] as f64,
                    v[0] as f64,
                    v[1] as f64,
                    v[2] as f64,
                )
            } else {
                format!(
                    "{} {time:.6} {step} {index} {:.ndec$} {:.ndec$} {:.ndec$}\n",
                    frame.frame, x[0] as f64, x[1] as f64, x[2] as f64
                )
            };
            if out.write_all(line.as_bytes()).is_err() {
                return 1;
            }
        }
        nframes += 1;
    }
    progress.finish();
    if out.flush().is_err() {
        return 1;
    }
    eprintln!(
        "Read {nframes} frame{} with {natoms_out} atoms each",
        if nframes == 1 { "" } else { "s" }
    );
    0
}
