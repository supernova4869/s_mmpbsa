//! Small helpers shared by the tools.

use indicatif::{ProgressBar, ProgressStyle};

/// Applies the progress bar style used throughout s_mmpbsa to `pb`.
///
/// The bars of this crate are drawn with the same template as the ones
/// s_mmpbsa prints itself, so that the progress of a whole run looks uniform.
pub fn set_style(pb: &ProgressBar) {
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:50.cyan/cyan} {pos}/{len} {msg}",
        )
        .expect("valid progress bar template")
        .progress_chars("=>-"),
    );
}

/// Same layout as [`set_style`], but without the colour of the bar.
///
/// The tools of this crate draw their bars with this style so that they keep
/// the colour of the surrounding output: s_mmpbsa switches the terminal to
/// grey while it extracts trajectories, and a bar that carried its own colour
/// would stand out of it.
pub fn set_style_plain(pb: &ProgressBar) {
    pb.set_style(
        ProgressStyle::with_template("[{elapsed_precise}] {bar:50} {pos}/{len} {msg}")
            .expect("valid progress bar template")
            .progress_chars("=>-"),
    );
}

/// Style for a bar whose total is unknown: the same layout without the bar.
pub fn set_spinner_style(pb: &ProgressBar) {
    pb.set_style(
        ProgressStyle::with_template("[{elapsed_precise}] {spinner} {msg}")
            .expect("valid progress bar template"),
    );
}
