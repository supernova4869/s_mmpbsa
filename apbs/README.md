# Rust APBS solver (apbs-rs) - vendored sources

These crates (`apbs-generic`, `apbs-pmgc`, `apbs-mg`) are vendored from the
Rust port of APBS, https://github.com/supernova4869/apbs-rs (version 3.4.1,
BSD-2-Clause).

They are compiled directly into the `s_mmpbsa` binary. The in-process driver is
`../src/apbs_runner.rs` (copied from `apbs-core/src/routines.rs` of apbs-rs).

To refresh these sources from a local apbs-rs checkout, run:

```bash
scripts/sync_apbs_rs.sh /path/to/apbs-rs
```
