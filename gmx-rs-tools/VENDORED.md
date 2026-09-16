# gmx-rs-tools (vendored)

This directory is a copy of https://github.com/supernova4869/gmx-rs-tools, the
Rust port of the GROMACS tools `dump`, `trjconv`, `convert-tpr` and `make_ndx`
that s_mmpbsa links in place of calling the GROMACS binary.

The crate is a workspace member and is used through
`s_mmpbsa/src/gmx.rs`:

| s_mmpbsa used to call | in-process equivalent |
| --- | --- |
| `gmx dump -s md.tpr` | `gmx_rs_tools::tpr::{TprFile, parse_body}` |
| `gmx trjconv ...` | `gmx_rs_tools::cmd::trjconv::run` |
| `gmx convert-tpr ...` | `gmx_rs_tools::cmd::convert_tpr::run` |
| `gmx make_ndx -f md.tpr -o index.ndx` | `gmx_rs_tools::cmd::make_ndx::run`, or `index::analyse` + `index::write_ndx` for the default groups |

## Local changes

The tree currently **tracks upstream verbatim**: there is no `local.patch`, so
the copy is exactly what https://github.com/supernova4869/gmx-rs-tools
contains.  The changes s_mmpbsa needed while it was being wired in have been
folded into the upstream project, and the file is only re-created when a local
change has to be carried for a while:

* `tpr.rs` — `Mtop::ffparams` and the `nbfp` accessors (`FfParams::lj_sr()`,
  `FfParams::atnr_usize()`) are part of the crate itself.
* `cmd/mod.rs` — `set_scripted_input` / `clear_scripted_input` let another
  program answer the interactive questions (index group selections, `make_ndx`
  editor commands) of the command front ends; `select_group` and
  `read_input_line` read from that script when one is installed and from
  standard input otherwise.
* `cmd/make_ndx.rs` — the editor reads its input through
  `cmd::read_input_line` and stays silent when a script is installed.
* `cmd/trjconv.rs` — the output XTC precision is copied from the input frame,
  like GROMACS does, instead of being rounded up to the next power of ten; a
  time window with an unbounded end (`-e inf`, which is what "the whole
  trajectory" is passed in) falls back to the position in the file for the
  progress bar instead of pinning it at 0%.
* `utils.rs` — `set_style` (the cyan style s_mmpbsa also applies to its own
  bars), `set_style_plain` (the same layout without the colour, used by the
  tools so that their bars keep the colour of the surrounding output) and
  `set_spinner_style`.
* `progress.rs` — the bars follow `progress::set_enabled` and the terminal
  only: the `GMXRS_PROGRESS` environment variable was dropped, and the styles
  come from `utils` instead of being built in this module.

To carry a new local change, write it into this directory and record the
difference against the checkout, e.g.:

```bash
for f in src/cmd/trjconv.rs src/progress.rs; do
    diff -uN --label "a/$f" --label "b/$f" \
        /path/to/gmx-rs-tools/"$f" gmx-rs-tools/"$f"
done > gmx-rs-tools/local.patch
```

## Refreshing

```bash
scripts/sync_gmx_rs_tools.sh /path/to/gmx-rs-tools
```

The script copies `src/`, `tests/`, `scripts/`, `Cargo.toml`, `README.md` and
`NOTES.md` from the checkout and re-applies `local.patch` when one exists; if
the patch no longer applies, the upstream files changed in the same places and
the changes have to be reconciled by hand.
