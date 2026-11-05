# s_mmpbsa
s_mmpbsa: Supernova's tool of binding free energy calculation for Gromacs trajectory, using molecular mechanics Poisson-Boltzmann surface area (MM-PBSA) method.

The s_mmpbsa program follows LGPL license, and can be freely used for academic purposes.

## Documentation
The detailed documentations could be found here. 
https://s-mmpbsa.readthedocs.io/en/latest/

## Introduction
MM-PBSA method is the most popular method to rapidly calculate binding free energy, especially for biological systems. However, as a widely-used MD program, Gromacs has not officially support MM-PBSA calculation. Although there have been numerous programs that can calculate binding free energy with GROMACS trajectory, most of them have some limitations at different aspects, for example: (1) Difficult to use and install (2) Lack of support for new version of Gromacs (3) Too slow (4) Not cross-platform. Instead, s_mmpbsa provides a convenient interface (like [Multiwfn](http://sobereva.com/multiwfn/)) to calculate binding free energy for GROMACS trajectory. 

## Features of s_mmpbsa
- Open source and freely available.
- Less need for running environment preparation: PB/SA and the Gromacs tools s_mmpbsa needs have been ported to Rust and are compiled into the executable, so no Gromacs or APBS installation is required. A Python environment is only needed for plotting.
- In contrast to other programs, s_mmpbsa is developed with Rust.
- Interactive operation, no need to write parameter files. Also, user can write shell script to invoke s_mmpbsa for batch use.
- Optional GPU acceleration of the MM pairwise kernel and the PB multigrid solver (`use_gpu = "y"` in settings.ini, or `use_gpu` in the config.yaml); falls back to the CPU automatically when no suitable GPU is found. The MM part (interaction entropy first, the per-frame sums taken from it where the trajectories share frames) always runs for the whole trajectory before the PB/SA terms, which are solved frame by frame.
- The GPU float precision is selectable with `gpu_precision` in settings.ini: `"fp64"` (default) uses double precision and requires the `SHADER_F64` GPU feature (unavailable on some GPUs, notably several AMD models), while `"fp32"` uses single precision, works on every GPU including AMD ones (where it is the only usable GPU mode), and clamps the PB stopping tolerance to `1e-5` (still well below the default `errtol`). When an AMD GPU is detected at calculation start, the precision is switched to `"fp32"` automatically with a notice. On the real-world test system the fp32 PB energies deviate from fp64 by at most ~0.03 kJ/mol per frame (~0.01%); the GPU fp32 kernels can be somewhat slower than fp64 on some NVIDIA drivers, but remain far faster than the CPU path.
- With `pb_fast_path` on (the default), all PB frames are solved on one shared grid — the elementwise maximum of the per-frame molecule bounds — and each frame's converged coarse- and fine-grid potentials warm-start the next frame's multigrid solve as exact initial guesses (`APBS_RUST_WARM_START=0` disables the warm start). The converged energies agree with the cold-start run to solver tolerance (~1e-6 relative) while needing fewer V-cycles per solve. With `pb_fast_path` off, every frame is gridded from its own bounds, as before. The coefficient fill and the GPU solve paths are further optimized for this shape of run: the surface-sphere marking buckets spheres per grid slab, the flat O(grid) fill loops and the harmonic-mean dielectric smoothing run in parallel, the prolongation-weight blocks (grid-geometry constants) are cached across solves, the multigrid work vectors are zeroed on-device instead of uploaded, and every solve reuses pooled GPU buffers, so consecutive frames neither reallocate nor re-transfer unchanged data. These optimizations (and the shared grid above) are controlled by `pb_fast_path` in settings.ini / config.yaml (default `"y"`); `"n"` falls back to per-frame grids, the plain serial fill and fresh per-solve GPU buffers.
- Considers electric screening effect, as "J. Chem. Inf. Model. 2021, 61, 2454".
- Considers conformational entropy, as "J. Chem. Phys. 2017, 146, 124124".
- Could store analysis results for further reproducable analysis.

## Main function
- Binding energy calculation from MD simulation.
- Alanine scanning of protein-ligand complex.

## Requirement

### Basic requirements
- To use plotting backend, you need to run `sudo apt install libfontconfig1-dev` on Ubuntu.
- Gromacs: not needed. The trajectory and run input file handling (`gmx dump`, `gmx trjconv`, `gmx convert-tpr`, `gmx make_ndx`) is performed in-process; see the "Built-in Gromacs tools" section below.
- PyMOL (optional) to plot the B-factor colored structure.

APBS is currently not needed. PB/SA calculations have been performed in-process totally by the Rust-ported APBS solver (apbs-rs), which has already been compiled into the s_mmpbsa binary.

## Usage
Although s_mmpbsa supports fixing PBC conditions to trajectory `.MMPBSA_[name].xtc`, it is still recommended to comfirm that the trajectory has been correct, using visualization software such as [VMD](http://www.ks.uiuc.edu/Research/vmd/):

```
vmd .MMPBSA_[name]_struct.gro .MMPBSA_[name].xtc
```

### MD Binding energy calculation:
``` bash
# Firstly, add s_mmpbsa folder to $PATH.
# Start s_mmpbsa, and input as follow (support # comments, but not recommended and usually no need to input with comments)
s_mmpbsa -f md_pbc.xtc -s md.tpr -n index.ndx
0 # go to next step (Trajectory Parameters)
1 # select receptor group
[protein group number]
2 # select ligand group
[ligand group number]
0 # go to next step (MM-PBSA Parameters)
# Other options usually no need to change. The PB and SA parameters could be modified by 10 and 11 (a commented PB_settings.yaml / SA_settings.yaml is written, edit it and press ENTER to reload; misspelled values are rejected with the list of valid choices). Menu item 9 loads a named PB parameter preset (default / fast / fine / npbe)
0 # go to next step (start calculation)
[return] # use default name "system" or input specific name
# Wait for calculation finish
-1 # write pdb file with residue-wised INVERSED binding energy filled in B-factor column
# input the time point
1 # view summary
2 # output energy by time
3 # output energy by residue
# input the time point (default average)
# input the residues range (default all)
0 # exit s_mmpbsa program
```

### Alanine scanning:
```bash
# At MM-PBSA Parameters page
2 # Do alanine scanning
1 # Select mutation residues by layers (ACS Catal. 2024, 14, 15, 11447–11456)
0 # Start calculation
...
```
The results will contain energy terms of both wild type and mutants.

### Use Analysis mode:
```bash
# Firstly, add s_mmpbsa folder to $PATH.
# Start s_mmpbsa, and input as follow (support # comments, but not recommended and usually no need to input with comments)
s_mmpbsa -a MMPBSA_[name]_WT.sm # analysis mode
... (same as above)
```

The data was generated with .csv format and plotted as figures; The pdb files with inversed binding energy filled in B-factors column are drawn as png figures by PyMOL (if usable).

## Download
Release file: https://github.com/supernova4869/s_mmpbsa/releases/latest, where "s_mmpbsa.exe" and "s_mmpbsa" are s_mmpbsa executable files on Windows and Linux operation systems, respectively.

## Citation
The [introduction paper of s_mmpbsa](https://www.mdpi.com/1420-3049/31/10/1683) has been published.

The developer would be deeply grateful if you use s_mmpbsa and appropriately cite it in your research as follows:

```
@article{s_mmpbsa,
   author = {Zhang, Jiaxing and Gu, Tao and Li, Chuanxi and Qi, Wei},
   title = {s_mmpbsa: A Lite and Cross-Platform MM-PBSA Program},
   journal = {Molecules},
   volume = {31},
   number = {10},
   pages = {1683},
   ISSN = {1420-3049},
   DOI = {10.3390/molecules31101683},
   year = {2026},
   type = {Journal Article}
}
```

## About developer
Dr. Jiaxing Zhang (Contact: zhangjiaxing7137@tju.edu.cn, Tianjin University)

If you encountered any difficulty while using s_mmpbsa, or you found any bugs, or you have any suggestion on improving s_mmpbsa, please E-mail me or join my QQ group 864191465 to describe.

## Built-in PBSA solver
- The Poisson-Boltzmann and surface-area calculations are powered by the Rust-ported APBS solver from https://github.com/supernova4869/apbs-rs, vendored under `apbs/` and compiled directly into the s_mmpbsa executable. No external APBS program or sub-process is used.
- To refresh the vendored solver code from a local apbs-rs checkout, run `scripts/sync_apbs_rs.sh`.

## Built-in Gromacs tools
- `gmx dump`, `gmx trjconv`, `gmx convert-tpr` and `gmx make_ndx` are replaced by the Rust port, vendored under `gmx-rs-tools/` and compiled directly into the s_mmpbsa executable. `src/gmx.rs` is the code that drives it, so no external Gromacs program or sub-process is used.
