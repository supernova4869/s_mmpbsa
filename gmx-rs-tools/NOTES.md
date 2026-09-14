# Where the logic comes from

The table below maps every non-trivial part of `gmx-rs-tools` to the GROMACS 2026.3
source it was extracted from.  Line references are to this repository.

## `gmx dump` — `src/gromacs/tools/dump.cpp`

| Rust | GROMACS |
| --- | --- |
| `cmd/dump.rs::dump_xtc` | `list_xtc()` (dump.cpp:228) and `pr_rvecs`/`pr_title` from `utility/txtdump.cpp` |
| `cmd/dump.rs::dump_trr` | `list_trr()` (dump.cpp:169) |
| `cmd/dump.rs::dump_tpr` | `list_tpr()` (dump.cpp:80) |
| `cmd/dump.rs::print_mtop` | `pr_mtop()` (a summary instead of the full dump) |
| `cmd/dump.rs::group_short_name` | `shortName(SimulationAtomGroupType)` |
| `cmd/dump.rs::dump_top`, `dump_mtx` | `list_top()`, `list_mtx()` |

The frame printing uses the same indentation (`INDENT = 3`), the same field
widths (`%10d`, `%12.7e`, `%12.5e`, `%10g`) and therefore produces identical
output for XTC and TRR input.

## `gmx trjconv` — `src/gromacs/tools/trjconv.cpp`

| Rust | GROMACS |
| --- | --- |
| `cmd/trjconv.rs::run` frame loop | `gmx_trjconv()` main loop (trjconv.cpp:1100-1600) |
| `cmd/trjconv.rs::select_group` | `qgroup()`/`rd_groups()` (`topology/index.cpp`) |
| `cmd/trjconv.rs::mk_filenm` | `mk_filenm()` (trjconv.cpp:90) |
| `cmd/trjconv.rs::rmod` | `bRmod()` |
| `cmd/trjconv.rs::build_title` | the `top_title`/`t=`/`step=` assembly in the write switch |
| `cmd/trjconv.rs::make_whole_graph` | `mk_multishell`/`mk_1shift`/`shift_self` (`pbcutil/mshift.cpp`) and `gmx_rmpbc_apply()` (`pbcutil/rmpbc.cpp`) |
| `cmd/trjconv.rs::build_moltype_adjacency` | `mk_graph_moltype()`/`mk_igraph()` (`pbcutil/mshift.cpp`), including the second phase that reconnects separate parts through virtual sites |
| `cmd/trjconv.rs::VSITE_TYPES` | the `def_vsite(...)` entries of `topology/ifunc.cpp` |
| `cmd/trjconv.rs::put_group_com_in_box` | `put_molecule_com_in_box()`/`put_residue_com_in_box()` (`pbcutil/pbcmethods.cpp`) |
| `pbc.rs::put_atoms_in_box` | `putAtomsInBoxTemplated` (pbc.cpp:1278) |
| `pbc.rs::put_atoms_in_triclinic_unitcell` | same name (pbc.cpp:1543) |
| `pbc.rs::put_atoms_in_compact_unitcell` | same name (pbc.cpp:1611) |
| `pbc.rs::calc_pbc_cluster` | `calc_pbc_cluster()` (pbcmethods.cpp) |
| `pbc.rs::pbc_dx_gmx`, `triclinic_shift_vectors`, `max_cutoff2` | `pbc_dx()` and `low_set_pbc()` (pbc.cpp:667, 460-600) |
| `pbc.rs::center_x` | `center_x()` (pbcmethods.cpp:409) |
| `pbc.rs::reset_x_ndim` | `reset_x_ndim()` (`math/do_fit.cpp:297`) |
| `pbc.rs::calc_fit_r`, `jacobi` | `calc_fit_R()` (`math/do_fit.cpp`) and `jacobi()` (`math/nrjac.cpp:69`) |
| `tpr.rs::assign_chain_ids` | `tpx_make_chain_identifiers()`/`ChainIdFiller` (`fileio/confio.cpp:255-337`) |
| `tpr.rs::FTUPD` | `ftupd[]` (`fileio/tpxio.cpp:280`) |
| `tpr.rs::TPX_GENERATION_ADD_SIZE_FIELD` | `TpxGeneration::AddSizeField` (`fileio/tpxio.cpp:238`) |

The `-pbc nojump` branch, the `-dump` nearest-frame selection, the
`-skip`/`-dt`/`-timestep` bookkeeping and the `-sep` file naming all follow the
original control flow, including the special case that the first frame is
compared against the structure file coordinates.

## `gmx convert-tpr` — `src/gromacs/tools/convert_tpr.cpp`

| Rust | GROMACS |
| --- | --- |
| `cmd/convert_tpr.rs::run` | `ConvertTpr::run()` (convert_tpr.cpp:376) |
| `cmd/convert_tpr.rs::print_runtime_info` | `print_runtime_info()` (convert_tpr.cpp:337) |
| `cmd/convert_tpr.rs::maxwell_speed` | `maxwell_speed()` (`gmxpreprocess/gen_maxwell_velocities.cpp`) |
| `tpr.rs::TprFile::set_nsteps` | the `ir->nsteps = ...` assignment plus `write_tpx_state()` |
| `tpr.rs::write_subset_body` | `reduce_topology_x()`, `reduce_atom()`, `reduce_rvec()`, `reduce_ilist()`, `reduce_listoflists()` and `bKeepIt()`/`invind()` (convert_tpr.cpp:120-250) |
| `tpr.rs::reduce_ilist` | `reduce_ilist()` (convert_tpr.cpp:198) |
| `tpr.rs::put_atoms`, `put_ilists` | `do_atoms()`/`do_ilists()` from `fileio/tpxio.cpp`, written in the compact `InMemorySerializer` format |

While the C++ tool re-serializes the complete run input file, this port keeps
the body blob and patches the eight bytes of the `int64 nsteps` field, which
produces the same file with a much smaller amount of code.

## `gmx make_ndx` — `src/gromacs/tools/make_ndx.cpp`

| Rust | GROMACS |
| --- | --- |
| `cmd/make_ndx.rs::Editor::run_editor` | `edit_index()` (make_ndx.cpp:1075) |
| `cmd/make_ndx.rs::Editor::parse_entry` | `parse_entry()` (make_ndx.cpp:960) |
| `cmd/make_ndx.rs::Editor::select_atomnumbers` | `select_atomnumbers()` |
| `cmd/make_ndx.rs::Editor::select_residuenumbers` | `select_residuenumbers()` and `select_residueindices()` |
| `cmd/make_ndx.rs::Editor::select_by_name` | `select_atomnames()`, `select_residuenames()`, `select_chainnames()`, `comp_name()` |
| `cmd/make_ndx.rs::Editor::split_group`, `split_chain` | `split_group()`, `split_chain()` |
| `cmd/make_ndx.rs::{or_groups, and_groups}` | `or_groups()`, `and_groups()` |
| `index.rs::analyse`, `analyse_prot`, `analyse_other` | `analyse()` (index.cpp:544), `analyse_prot()` (index.cpp:345), `analyse_other()` (index.cpp:211) |
| `index.rs::read_ndx`, `write_ndx` | `init_index()`, `write_index()` |
| `index.rs::find_group`, `strcasecmp_min` | `findGroupTemplated()`, `gmx_strcasecmp_min()` |
| `index.rs::ResidueTypeMap` | `residueTypeMapFromLibraryFile()`/`typeOfNamedDatabaseResidue()` (`topology/residuetypes.cpp`) |

## File formats

| Rust | GROMACS |
| --- | --- |
| `xdr.rs` | `fileio/xdr_serializer.cpp` (including the extra `len + 1` of `doString`) |
| `xtc.rs` | `fileio/xtcio.cpp` and `fileio/libxdrf.cpp` (`sendbits`, `receivebits`, `sendints`, `receiveints`, `sizeofint`, `sizeofints`, `xdr3dfcoord`) |
| `trr.rs` | `fileio/trrio.cpp` |
| `gro.rs` | `fileio/groio.cpp` (`get_w_conf`, `write_hconf_indexed_p`, `write_hconf_box`) |
| `pdb.rs` | `fileio/pdbio.cpp` (`read_atom`, `read_cryst1`, `write_pdbfile`, `gmx_fprintf_pdb_atomline`, `gmx_write_pdb_box`) |
| `tpr.rs` | `fileio/tpxio.cpp`: `do_tpxheader`, `do_tpx_state_first`, `do_tpx_state_second`, `do_mtop`, `do_ffparams`, `do_iparams`, `do_atom`, `do_inputrec` (prefix), `atomicnumber_to_element` |
| `tpr.rs::Mtop::global_atoms` | `gmx_mtop_global_atoms()`/`atomcat()` (`topology/mtop_util.cpp`) and `gmx_mtop_t::finalize()` (`topology/topology.cpp`) |

The TPR body is *not* XDR: `write_tpx_state()` serializes it with
`gmx::InMemorySerializer`, which stores values in their native width (swapped to
big endian on little endian hosts), encodes strings as a `uint64` length
followed by the raw characters, and never pads.  `tpr.rs::CReader` implements
exactly that encoding.

# Known limitations

`gmx-rs-tools coords` is not part of GROMACS; it is a thin front end over
`trx::CoordinateReader`, which is the same streaming reader used by
`gmx dump -f` and `trjconv`.  Its `-b`/`-e`/`-skip` semantics were checked
against `gmx trjconv` (the frame counter starts at the first frame at or after
`-b`), and the coordinates match `gmx dump -f` exactly.

These are deliberate simplifications of the minimal port:

* Only TPR files with `fileVersion >= 119` (`tpxv_AddSizeField`, GROMACS 2021
  and later) are decoded; older layouts are rejected with a clear message.
* The header is decoded for every generation, including the `sizeOfTprBody`
  field that was added in `TpxGeneration::AddSizeField` (= 27, i.e. GROMACS
  2022 and later).  The `ftupd[]` table is applied so that files written before
  an interaction function type existed shift their function types back into the
  current numbering (and their interaction lists are skipped, as they are
  absent from those files).
* The inputrec is decoded up to `epsilon_surface`.  `convert-tpr` only needs
  `nsteps`, `init_step`, `init_t` and `delta_t`, and `dump` prints the fields it
  decoded; the remaining inputrec fields stay unparsed inside the body blob.
* `convert-tpr` reproduces two quirks of `reduce_topology_x()` that make the
  subset file differ from a "clean" rewrite: (1) `reduce_atom()` compares the
  original residue index of each atom with the *already renumbered* index of
  the previous one, so a selection that does not start at the first residue
  ends up with one residue per atom; (2) the `atomtype`/`atomtypeB` arrays are
  never permuted, so the selection receives the atom types of the first `n`
  atoms of the system.  Both are reproduced so that the output matches
  `gmx convert-tpr` byte for byte.
* `dump` prints a compact topology summary where `pr_mtop()` prints every
  topology sub-structure, and `-om` (write an mdp file) is not implemented.
* `dump -p` echoes the topology file instead of running the C preprocessor on
  it, and `dump -e`/`-cp` (energy/checkpoint files) are not implemented.
* `trjconv` does not implement `-sub` (removed from GROMACS as well),
  `-cluster`, `-drop`/`-dropunder`/`-dropover`, `-conect`, `-split`, `-exec`,
  `-fr` and TNG input/output.
* `-fit rot+trans` accumulates the rotation in single precision like GROMACS;
  a handful of output values can still land on the other side of the last
  printed decimal (≤ 12 of 19383 lines in the reference test case).
* The same applies to `-fit translation` on some systems (1 of 19446 lines),
  where the printed value is `-0.000` instead of `0.000`.
* `make_ndx` reads `residuetypes.dat` from `GMXDATA`/`GMXLIB` or the usual
  install locations and falls back to a built-in table of common residue names
  when the file is not available.
* Only the interaction types GROMACS itself uses in `mk_graph_moltype()` for
  reconnecting separate parts are modelled, i.e. the virtual sites; position
  restraints and similar types are ignored during molecule unwrapping.
