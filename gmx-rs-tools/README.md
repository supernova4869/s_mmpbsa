# gmx-rs-tools — a minimal Rust port of four GROMACS tools

This crate extracts the core logic of four GROMACS 2026.3 command line tools and
re-implements it in Rust, without depending on any part of the GROMACS C++
code base:

| GROMACS command | C++ source | Rust |
| --- | --- | --- |
| `gmx dump` | `src/gromacs/tools/dump.cpp` | `src/cmd/dump.rs` |
| `gmx trjconv` | `src/gromacs/tools/trjconv.cpp` | `src/cmd/trjconv.rs` |
| `gmx convert-tpr` | `src/gromacs/tools/convert_tpr.cpp` | `src/cmd/convert_tpr.rs` |
| `gmx make_ndx` | `src/gromacs/tools/make_ndx.cpp` | `src/cmd/make_ndx.rs` |

The file readers are also usable as a library, so coordinates can be read
directly from a trajectory (see below).

The file formats these tools rely on are implemented from scratch as well:
XDR, the compressed XTC coordinates (`libxdrf.cpp`), TRR, GRO, PDB, index files
and the run input (TPR) container.  See [NOTES.md](NOTES.md) for the
function-by-function mapping to the C++ code and for the known simplifications.

`dump -s` (and the run-input handling of the other tools) decodes **every tpx
version GROMACS itself accepts**, i.e. version 58 and later; both body
encodings are supported (plain XDR before `tpxv_AddSizeField`, the compact
in-memory serializer from GROMACS 2021 on).

## Building and running

The crate uses `indicatif` (and, through it, `console`) for the progress bar of
the streaming tools; everything else is implemented from scratch.  Both crates
are in the local registry, so the build works offline:

```console
$ cd rust/gmx-rs-tools
$ cargo build --release
$ ./target/release/gmx-rs-tools            # usage
```

The binary is a multicall driver, like `gmx`:

```console
$ gmx-rs-tools dump -f traj.xtc
$ gmx-rs-tools dump -s topol.tpr
$ gmx-rs-tools trjconv -f traj.xtc -s topol.tpr -o out.gro -t0 0 -pbc mol -ur compact
$ gmx-rs-tools convert-tpr -s topol.tpr -o longer.tpr -nsteps 10000
$ gmx-rs-tools make_ndx -f topol.tpr -o index.ndx
$ gmx-rs-tools coords -f traj.xtc -s topol.tpr -sel Protein -b 10 -e 20
```

Group selections are read from standard input exactly like GROMACS does, so
existing scripts keep working (`printf '1\n0\n' | gmx-rs-tools trjconv ... -fit
rot+trans`).

## Implemented logic

**`coords`** (an addition on top of the extracted logic, not a GROMACS
command) reads coordinates straight out of a trajectory:

```console
$ gmx-rs-tools coords -f traj.xtc -s topol.tpr -sel Protein -b 10 -e 20 -skip 5 -o coords.txt
$ head -3 coords.txt
# frame time step index x y z
250 10.000000 5000 1 0.083000 1.261000 0.615000
250 10.000000 5000 2 0.046000 1.229000 0.551000
```

One line per atom (`frame`, `time`, `step`, 1-based atom index, coordinates),
with `-vel` adding velocities.  Atom selection uses the same index groups as
the other tools: `-n index.ndx` plus a group, `-sel <name|number>` to name one
directly, or `-s topol.tpr` for the default groups.  `-b`/`-e`/`-skip` behave
exactly like `gmx trjconv`'s and frames are streamed, so a window can be pulled
out of a 25 GB trajectory without loading it.

While frames are being read, `trjconv`, `dump -f` and `coords` draw a progress
bar on standard error:

```text
[00:00:24] ===============================>------------------ 64/100 frame       6  t=  5000.000 ps
```

The bar (`indicatif`) only appears when standard error is a terminal, so pipes
and logs stay clean.  The style comes from `utils::set_style`, the same
template s_mmpbsa prints its own bars with.  A host program that links the
crate can control the bars with `progress::set_enabled`: `false` switches them
off and `true` forces them on (they are then drawn on `/dev/tty`).

## Reading coordinates

The same capability is available as a library API:

```rust
use gmx_rs_tools::trx::{CoordinateReader, FrameRange};

// None = every atom; a slice of atom indices selects a subset.
let mut reader = CoordinateReader::open("traj.xtc", None, FrameRange {
    begin: Some(10.0),
    end: Some(20.0),
    skip: 5,
    max_frames: None,
})?.1;

while let Some(frame) = reader.next_frame()? {
    println!("frame {} t={:?}: {} atoms", frame.frame, frame.time, frame.x.len());
    let first = frame.x[0];          // [f32; 3], nm
    let _ = first;
}
```

`CoordinateReader` streams frames (constant memory), `read_coordinates()`
collects them into a `Vec<CoordFrame>` when that is more convenient, and
`trx::FrameSource` gives the full frames (including topology for structure
files) if the coordinates are not the only thing needed.

**`dump`** prints XTC/TRR frames (byte-identical formatting to GROMACS,
including the `%12.5e`/`%g` field widths) and the complete `gmx dump -s`
listing of a run input file: the full inputrec (including the force field
parameter blocks, free energy, pull, AWH, enforced rotation, IMD, walls, swap
ions, QMMM, the `applied-forces`/`fast-multipole-method` module parameter tree
and `grpopts`), the header, the whole topology (`pr_mtop()`, molecule by
molecule), the state matrices, the coordinates/velocities and the group
statistics.  `-nr`/`-nonr`, `-param` and `-orgir` behave like the original.

**`trjconv`** implements the frame loop of the original: `-skip`, `-dt`,
`-round`, `-dump`, `-t0`, `-timestep`, `-pbc none|atom|mol|res|whole|nojump|cluster`,
`-ur rect|tric|compact`, `-center`, `-boxcenter`, `-box`, `-trans`, `-shift`,
`-fit none|rot+trans|rotxy+transxy|translation|transxy|progressive`, `-ndec`,
`-vel`, `-force`, `-sep`, `-nzero`, `-b`/`-e` time ranges, plus index group
selection (`-n`).

Like the original, `-s` defaults to `topol.tpr` and `-n` to `index.ndx` (the
latter only when it exists), so the tool can be run from a simulation directory
without repeating those options.  If `topol.tpr` is needed but missing, the run
stops with GROMACS' "File input/output error" message rather than silently
writing a topology-less output file.

Trajectories are read frame by frame, so `trjconv` and `dump` handle
trajectories far larger than memory: converting a time window out of a 25 GB
trajectory only reads the frames inside that window.

**`convert-tpr`** implements `-nsteps`, `-extend` and `-until` by patching the
`nsteps` field of the serialized inputrec, and reports the runtime information
before and after the change.  The output is the input file with only the
`nsteps` bytes modified.  Like the original it also always asks for an index
group (default groups when no `-n` is given) and can write a subset tpr
(`-n index.ndx`, or just picking a group such as `Protein` or `SOL`):

```console
$ printf '1\n' | gmx-rs-tools convert-tpr -s md.tpr -o protein.tpr
Will write subset Protein of original tpx containing 3638 atoms
```

The subset file is rebuilt with a freshly serialized topology (atoms,
interaction lists, exclusions, molecule block) while the force field
parameters, symbol table, groups and inputrec are copied through unchanged.

**`make_ndx`** implements the interactive editor: `a`, `t`, `r`, `ri`, `chain`,
`res`, `name`, `del`, `keep`, `case`, `splitres`, `splitat`, `splitch`, `l`,
`h`, `!`, `&`, `|`, quoted group names, and the automatic generation of the
default groups (`System`, the ten protein groups, `Water`/`SOL`, `Ion`, one
group per remaining residue name, `Water_and_ions`).

## Verification

### Run input files of every version

`gmx dump -s` was compared **line by line** against `gmx` for **all 1341 `.tpr`
files on the test machine**, covering tpx versions 73 to 138 / GROMACS 4.5.5 to
2026.3, in four modes (default, `-orgir yes`, `-param` and `-nonr`): all four
modes are byte-identical for 1340 files.  The one exception has an interaction
list entry with type index `-1`, which makes `gmx dump -s -param` index
`functype[-1]`/`iparams[-1]` and print uninitialized memory; this port prints
zeros there (see [NOTES.md](NOTES.md)).

| GROMACS version | tpx version | files | result |
| --- | --- | --- | --- |
| 4.5.5 / 4.5.5-dev | 73 | 4 | identical |
| 2016 / 2016-dev / 2017-dev | 110 / 111 | 3 | identical |
| 2019.6 | 116 | 91 | identical |
| 2021.3 / 2021.5 | 122 | 272 | identical |
| 2022.4 / 2022.5 | 127 | 13 | identical |
| 2023-rc1 / 2023.1 / 2023.2 / 2023.3 | 129 | 75 | identical |
| 2024.1 / 2024.2 | 133 | 147 | identical |
| 2025.2 | 137 | 668 | identical |
| 2026.1 / 2026.3 | 138 | 68 | identical |

The five run input files that ship in the GROMACS source tree (including the
4.5.5 and double-precision 2016 ones) are parsed by a regression test in
`tests/roundtrip.rs`, and `make_ndx`, `trjconv` (all `-pbc` modes) and
`convert-tpr` were checked against `gmx` on old-version files as well.

### Command behaviour

Every implementation was compared against the matching GROMACS 2026.3 binary
(`gmx`, built from this source tree) on a 21 frame, 923 atom trajectory plus its
`topol.tpr`:

| Case | Result |
| --- | --- |
| `dump -f traj.xtc` | byte identical (19530 lines) |
| `dump -f traj.trr` | byte identical |
| `trjconv -o out.xtc` | byte identical (XTC compression is bit exact) |
| `trjconv -o out.trr` | byte identical |
| `trjconv -o out.gro` (default, `-t0`, `-skip`, `-dt`, `-dump`, `-sep`) | byte identical |
| `trjconv -o out.pdb` | byte identical |
| `trjconv -pbc atom/mol/res/whole/nojump`, `-ur rect/tric/compact` | byte identical |
| `trjconv -pbc cluster` (all `-ur`/`-boxcenter` combinations) | byte identical |
| `trjconv -center`, `-trans`, `-shift`, `-ndec`, `-fit translation` | byte identical |
| `trjconv -b/-e` time windows | byte identical |
| `trjconv -fit rot+trans` / `rotxy+transxy` / `progressive` | within 1 ULP (≤ 12 of 19383 GRO lines differ in the last printable digit) |
| `convert-tpr -nsteps 1000` | only the 2 bytes of `nsteps` differ from the input |
| `make_ndx` (default groups and scripted editor sessions) | byte identical |

A production system was used as a second check: a 5168736 byte `md.tpr`
(145382 atoms, 5 molecule types, written by GROMACS 2025.2) and its 25 GB
trajectory.

| Case | Result |
| --- | --- |
| `dump -s md.tpr` | byte identical (whole dump, including the 48908 line topology) |
| `dump -f md.xtc` (first 2000 lines) | byte identical |
| `trjconv -dump 0.2 -o out.gro` | byte identical (0.5 s, single frame read) |
| `trjconv -b 2 -e 6 -pbc mol -o out.gro` | byte identical (19.6 MB, only the frames in the window are read) |
| `make_ndx -f md.tpr` (48908 lines) | byte identical |
| `make_ndx` scripted session | byte identical |
| `convert-tpr -nsteps 1000` | only the 4 bytes of `nsteps` differ |
| `convert-tpr` subset tprs (all 8 default groups, and groups read from an index file) | same size and identical `gmx dump -s` output; a subset tpr also drives `gmx trjconv` to identical results |

All `trjconv` PBC modes (`atom`, `mol`, `res`, `whole`, `nojump` with `-ur
rect/tric/compact`), `-pbc cluster`, the output formats (`gro`, `xtc`, `trr`,
`pdb`), `-skip`,
`-fit translation`, `-b`/`-e` and `-dump` were compared on this system as well
and are byte identical.

## Performance

`trjconv` streams frames and writes the output through a buffered writer, and
the fixed columns of `gro`/`pdb` atom records are formatted once per conversion
instead of once per frame.  On the 145382 atom system (21 frames):

| Case | gmx-rs-tools | gmx |
| --- | --- | --- |
| `-o out.gro` | 0.55 s | 0.75 s |
| `-o out.xtc` | 0.17 s | 0.14 s |
| `-pbc mol -o out.gro` | 0.58 s | 0.81 s |
| `-pbc mol -ur compact -o out.gro` | 0.59 s | 0.82 s |

`gmx-rs-tools dump -f` is still about 1.5x slower than `gmx dump -f` because the
per-line formatting dominates; it does print the identical text.

The comparison is automated:

```console
$ GMXRS_PARITY_TESTDATA=/path/with/topol.tpr+traj.xtc \
  GMX_BIN=/opt/gromacs/bin/gmx cargo test --test gmx_parity
```

`scripts/make_testdata.sh` builds such a directory from the systems that ship
in the GROMACS source tree.

The format level tests (`cargo test`) are self contained: they generate their
own XTC/TRR/GRO/index data and check round trips, so they run without GROMACS.

## Layout

```text
src/xdr.rs     XDR primitives and the GROMACS string encoding
src/xtc.rs     XTC reader/writer including the compressed coordinate codec
src/trr.rs     TRR reader/writer
src/gro.rs     GROMOS-87 coordinate files
src/pdb.rs     PDB coordinate files
src/tpr.rs     TPR container, header, topology, state and inputrec
src/index.rs   index files, group matching and default group generation
src/pbc.rs     periodicity, centering and the quaternion least squares fit
src/trx.rs     trajectory format dispatch
src/cmd/*.rs   the four command line front ends
```
