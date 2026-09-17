//! Frame progress for the tools that stream trajectories.
//!
//! The bar is drawn with [`indicatif`].  GROMACS prints an in-line
//! `\rReading frame ...` line instead; a bar carries the same information and
//! is easier to read.
//!
//! The bar is written to standard error so that it never mixes with the data a
//! tool writes to standard output, and `indicatif` hides it automatically when
//! that stream is not a terminal, which keeps pipes and log files clean.
//! [`set_enabled`] switches the bar off for the whole program, or forces it on
//! (drawing on `/dev/tty` when standard error is redirected); nothing else,
//! in particular no environment variable, changes whether it is drawn.

use std::io::IsTerminal;
use std::sync::atomic::{AtomicU8, Ordering};

use indicatif::{ProgressBar, ProgressDrawTarget};

use crate::trx::ReadProgress;

/// Redraw rate of the bar, in hertz.
const REFRESH_HZ: u8 = 10;
/// Length the byte fraction is scaled to when the frame total is unknown.
const FRACTION_LENGTH: u64 = 1000;

/// Program wide switch: 0 = not set, 1 = forced on, 2 = forced off.
static ENABLED: AtomicU8 = AtomicU8::new(0);

/// Switches the progress bar on or off for the whole program.
///
/// The default follows the terminal: the bar is drawn while standard error is
/// a terminal and hidden when the output is piped or logged.  A host program
/// that links this crate and prints its own progress can call
/// `set_enabled(false)` once at start-up to keep the embedded tools quiet, or
/// `set_enabled(true)` to keep the bar even when standard error is redirected.
pub fn set_enabled(enabled: bool) {
    ENABLED.store(if enabled { 1 } else { 2 }, Ordering::Relaxed);
}

/// Whether the bars are drawn.
///
/// This is what [`set_enabled`] installed; without it the terminal decides,
/// see [`Progress::new`].
pub fn is_enabled() -> bool {
    enabled()
}

fn enabled() -> bool {
    match ENABLED.load(Ordering::Relaxed) {
        1 => true,
        2 => false,
        // `indicatif` only draws when standard error is a terminal, so the
        // decision can be left to it as well.
        _ => std::io::stderr().is_terminal(),
    }
}

/// Draws the frame progress of a running conversion.
pub struct Progress {
    bar: Option<ProgressBar>,
    enabled: bool,
    spinner: bool,
    /// Requested time window (`-b`/`-e`), when the caller gave a bounded one.
    window: Option<(f64, f64)>,
}

impl Progress {
    /// Creates a bar that follows [`set_enabled`] and the terminal.
    pub fn new() -> Self {
        Progress {
            bar: None,
            enabled: enabled(),
            spinner: false,
            window: None,
        }
    }

    /// Makes the bar follow a requested time window instead of the position in
    /// the input.
    ///
    /// A conversion that reads a window out of a much larger trajectory would
    /// otherwise show almost no progress: the frames before the window are read
    /// and discarded, so its position in the file says very little.
    pub fn set_time_window(&mut self, window: Option<(f64, f64)>) {
        self.window = window;
    }

    /// Reports progress: `frames` frames have been processed out of `total`
    /// and the current frame is at `time`.
    ///
    /// A `total` of zero means that the number of frames is not known; the bar
    /// then follows `fraction`, how much of the input has been read, and the
    /// message carries the frame counter alone.  Counting the frames of an
    /// `xtc`/`trr` input would cost a pass over the whole file, which is not
    /// worth it in front of a conversion of that same file.
    pub fn update(&mut self, frames: u64, total: u64, fraction: Option<f64>, time: f64) {
        if !self.enabled {
            return;
        }
        let spinner = total == 0 && fraction.is_none();
        let message = if total > 0 {
            format!("frame {frames:>7}/{total:<7} t={time:>10.3} ps")
        } else {
            format!("frame {frames:>7}  t={time:>10.3} ps")
        };
        match &self.bar {
            Some(bar) => {
                if self.spinner != spinner {
                    self.spinner = spinner;
                    apply_style(bar, spinner);
                }
                set_position(bar, frames, total, fraction);
                bar.set_message(message);
                if spinner {
                    bar.tick();
                }
            }
            None => {
                // The bar is configured while hidden and gets its draw target
                // last, so that it is drawn with the state of the current frame
                // instead of an empty bar.
                let bar = ProgressBar::hidden();
                apply_style(&bar, spinner);
                set_position(&bar, frames, total, fraction);
                bar.set_message(message);
                bar.set_draw_target(draw_target());
                self.bar = Some(bar);
                self.spinner = spinner;
            }
        }
    }

    /// Reports progress from a [`ReadProgress`], which either knows the number
    /// of frames of the input (`gro`, `pdb`) or how much of a streamed file has
    /// been read.
    pub fn update_from(&mut self, read: Option<ReadProgress>, frames: u64, time: f64) {
        let total = read.and_then(|p| p.known_total()).unwrap_or(0);
        let fraction = match self.window {
            Some((begin, end)) if end.is_finite() && end > begin => {
                Some(((time - begin) / (end - begin)).clamp(0.0, 1.0))
            }
            _ => read.and_then(|p| p.fraction()),
        };
        self.update(frames, total, fraction, time);
    }

    /// True when the bar is actually drawn.
    ///
    /// Callers that have to do extra work for the bar (counting the frames of
    /// the input, which is a pass over the file) can use this to skip it when
    /// nothing is drawn anyway, the way it happens when the output is piped.
    pub fn is_enabled(&self) -> bool {
        self.enabled
    }

    /// Removes the bar so that other messages get a clean line.
    ///
    /// The next [`Progress::update`] draws a fresh bar.
    pub fn pause(&mut self) {
        if let Some(bar) = self.bar.take() {
            bar.finish_and_clear();
        }
    }

    /// Shows the final state and leaves it on screen.
    pub fn finish(&mut self) {
        if let Some(bar) = self.bar.take() {
            bar.finish();
        }
    }
}

impl Drop for Progress {
    fn drop(&mut self) {
        self.finish();
    }
}

/// Applies the shared bar style, determinate or counters-only.
fn apply_style(bar: &ProgressBar, spinner: bool) {
    if spinner {
        crate::utils::set_spinner_style(bar);
    } else {
        // The colourless variant: the tools keep the colour of the output
        // around them (s_mmpbsa prints grey sections around these bars).
        crate::utils::set_style_plain(bar);
    }
}

/// Draw target of the bar: the terminal of standard error, or `/dev/tty`
/// when the bar was forced on and standard error is redirected.
fn draw_target() -> ProgressDrawTarget {
    terminal_target().unwrap_or_else(|| ProgressDrawTarget::stderr_with_hz(REFRESH_HZ))
}

/// Scales the bar to the frame counter, or to the byte fraction when the
/// number of frames is not known.
fn set_position(bar: &ProgressBar, frames: u64, total: u64, fraction: Option<f64>) {
    if total > 0 {
        bar.set_length(total);
        bar.set_position(frames.min(total));
    } else if let Some(fraction) = fraction {
        bar.set_length(FRACTION_LENGTH);
        bar.set_position((fraction.clamp(0.0, 1.0) * FRACTION_LENGTH as f64).round() as u64);
    }
}

/// Terminal used when [`set_enabled(true)`](set_enabled) forced the bar on.
///
/// `indicatif` skips drawing when standard error is not a terminal, so the bar
/// is drawn on `/dev/tty` in that case: the messages can be piped or logged
/// while the bar is still visible on the terminal.
fn terminal_target() -> Option<ProgressDrawTarget> {
    if ENABLED.load(Ordering::Relaxed) != 1 {
        return None;
    }
    #[cfg(unix)]
    {
        let reader = std::fs::File::open("/dev/null").ok()?;
        let writer = std::fs::File::options().write(true).open("/dev/tty").ok()?;
        Some(ProgressDrawTarget::term(
            console::Term::read_write_pair(reader, writer),
            REFRESH_HZ,
        ))
    }
    #[cfg(not(unix))]
    {
        None
    }
}
