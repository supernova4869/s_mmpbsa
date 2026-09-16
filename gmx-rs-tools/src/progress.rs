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

/// Redraw rate of the bar, in hertz.
const REFRESH_HZ: u8 = 10;

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
}

impl Progress {
    /// Creates a bar that follows [`set_enabled`] and the terminal.
    pub fn new() -> Self {
        Progress {
            bar: None,
            enabled: enabled(),
            spinner: false,
        }
    }

    /// Reports progress: `frames` frames have been processed out of `total`
    /// and the current frame is at `time`.
    ///
    /// The `pos/len` fields of the bar are the frame counter and the number of
    /// frames in the input, so a `total` of zero (unknown, e.g. because the
    /// count was not taken) draws the frame counter of the message only.
    pub fn update(&mut self, frames: u64, total: u64, time: f64) {
        if !self.enabled {
            return;
        }
        let spinner = total == 0;
        let message = if spinner {
            format!("frame {frames:>7}  t={time:>10.3} ps")
        } else {
            format!("t={time:>10.3} ps")
        };
        match &self.bar {
            Some(bar) => {
                if self.spinner != spinner {
                    self.spinner = spinner;
                    apply_style(bar, spinner);
                }
                set_position(bar, frames, total);
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
                set_position(&bar, frames, total);
                bar.set_message(message);
                bar.set_draw_target(draw_target());
                self.bar = Some(bar);
                self.spinner = spinner;
            }
        }
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

/// Scales the bar to the frame counter.
fn set_position(bar: &ProgressBar, frames: u64, total: u64) {
    if total > 0 {
        bar.set_length(total);
        bar.set_position(frames.min(total));
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
