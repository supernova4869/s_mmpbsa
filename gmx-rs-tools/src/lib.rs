//! Minimal Rust re-implementation of a subset of the GROMACS command line tools
//! `gmx dump`, `gmx trjconv`, `gmx convert-tpr` and `gmx make_ndx`.
//!
//! The modules in this crate mirror the corresponding C++ sources of GROMACS
//! 2026.3 (see `NOTES.md` for the file-by-file mapping).  Only the parts of the
//! logic needed by the four tools above are implemented.

pub mod cmd;
pub mod frame;
pub mod gro;
pub mod index;
pub mod pbc;
pub mod pdb;
pub mod trr;
pub mod trx;
pub mod tpr;
pub mod xdr;
pub mod xtc;
