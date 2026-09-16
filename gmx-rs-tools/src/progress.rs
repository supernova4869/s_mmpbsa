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
/// Length the fraction is scaled to; `utils::set_style` prints it as the
/// `pos/len` fields of the bar, so this reads as a percentage.
const LENGTH: u64 = 100;

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

/// Draws the frame progress of a running conversion.
pub struct Progress {
    bar: Option<ProgressBar>,
    enabled: bool,
    spinner: bool,
}

impl Progress {
    /// Creates a bar that follows [`set_enabled`] and the terminal.
    pub fn new() -> Self {
        let enabled = match ENABLED.load(Ordering::Relaxed) {
            1 => true,
            2 => false,
            // `indicatif` only draws when standard error is a terminal, so the
            // decision can be left to it as well.
            _ => std::io::stderr().is_terminal(),
        };
        Progress {
            bar: None,
            enabled,
            spinner: false,
        }
    }

    /// Reports progress: `fraction` of the requested work is done, `frames`
    /// frames have been processed and the current frame is at `time`.
    ///
    /// A `None` fraction (the total is unknown) draws the counters only.
    pub fn update(&mut self, fraction: Option<f64>, frames: u64, time: f64) {
        if !self.enabled {
            return;
        }
        let message = format!("frame {frames:>7}  t={time:>10.3} ps");
        let spinner = fraction.is_none();
        match &self.bar {
            Some(bar) => {
                if self.spinner != spinner {
                    self.spinner = spinner;
                    apply_style(bar, spinner);
                }
                set_position(bar, fraction);
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
                set_position(&bar, fraction);
                bar.set_message(message);
                bar.set_draw_target(draw_target());
                self.bar = Some(bar);
                self.spinner = spinner;
            }
        }
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

/// Scales `fraction` to the length of the bar.
fn set_position(bar: &ProgressBar, fraction: Option<f64>) {
    if let Some(fraction) = fraction {
        bar.set_length(LENGTH);
        bar.set_position((fraction.clamp(0.0, 1.0) * LENGTH as f64).round() as u64);
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
