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

/// Style of the tools of this crate: the same layout as [`set_style`] without
/// the colour and without the counters.
///
/// The colour is left to the surrounding output, since s_mmpbsa switches the
/// terminal to grey while it extracts trajectories and a bar that carried its
/// own colour would stand out of it.  The frame counter is part of the message
/// instead of the `pos`/`len` fields, because the number of frames of a
/// streamed trajectory is not always known (see `progress::update`).
pub fn set_style_plain(pb: &ProgressBar) {
    pb.set_style(
        ProgressStyle::with_template("[{elapsed_precise}] {bar:50} {msg}")
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
