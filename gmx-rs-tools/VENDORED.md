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

`local.patch` holds the changes s_mmpbsa needs on top of the crate as it is
published upstream:

* `tpr.rs` — `Mtop::ffparams` keeps the force field parameters (the `nbfp`
  array) that the reader previously skipped, so the Lennard-Jones `c6`/`c12`
  pairs can be read straight out of the run input file.
* `cmd/mod.rs` — `set_scripted_input` / `clear_scripted_input` let another
  program answer the interactive questions (index group selections, `make_ndx`
  editor commands) of the command front ends; `select_group` and
  `read_input_line` read from that script when one is installed and from
  standard input otherwise.
* `cmd/make_ndx.rs` — the editor reads its input through
  `cmd::read_input_line` and stays silent when a script is installed.
* `cmd/trjconv.rs` — the output XTC precision is copied from the input frame,
  like GROMACS does, instead of being rounded up to the next power of ten.

## Refreshing

```bash
scripts/sync_gmx_rs_tools.sh /path/to/gmx-rs-tools
```

The script copies `src/`, `tests/`, `scripts/`, `Cargo.toml`, `README.md` and
`NOTES.md` from the checkout and re-applies `local.patch`; if the patch no
longer applies, the upstream files changed in the same places and the changes
have to be reconciled by hand.
