# s_mmpbsa
s_mmpbsa: Supernova's tool of binding free energy calculation for Gromacs trajectory, using molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method.

The s_mmpbsa program can be freely used for both academic and commerical purposes.

## Introduction

MM/PB-SA method is the most popular method to rapidly calculate binding free energy, especially for biological systems. However, as a widely-used MD program, Gromacs has not officially support MM/PB-SA calculation. Although there have been numerous programs that can calculate binding free energy with GROMACS trajectory, most of them have some limitations at different aspects, for example: (1) Difficult to use and install (2) Lack of support for new version of Gromacs (3) Too slow (4) Not cross-platform. Instead, s_mmpbsa provides a convinent interface (like [Multiwfn](http://sobereva.com/multiwfn/)) to calculate binding free energy for GROMACS trajectory. 

## Features of s_mmpbsa

- Open source and freely available.
- No need for preparing running environment (e.g., for Python), only needs Gromacs program when running on linux system. (In contrast to other programs such as gmx_MMPBSA.py, s_mmpbsa is developed with Rust).
- Interactive operation, no need to write parameter files. Also, user can write shell script to invoke s_mmpbsa for batch use.
- Very fast. Due to the efficency of rust program.
- Considers electric screening effect, as "CHIN. PHYS. LETT. 2021, 38(1), 018701" describes.

## Usage

### Preparation

Before calculation, you should fix the periodic boundary conditions (PBC) of the system.

If studying the protein-ligand system (e.g., enzyme-substrate), then it is better to fix PBC with `cluster` option, i.e.:

```bash
# build group for protein + ligand
gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc cluster -center -dt 1000 -o md_pbc.xtc
# select the protein + ligand group
```

If studying the system of two phase (e.g., solvent extraction), just fix PBC with `mol` option to ensure the structure completion of each molecule, i.e.:

```bash
# build group for two phase
gmx trjconv -f md.xtc -s md.tpr -n index.ndx -pbc mol -dt 1000 -o md_pbc.xtc
# select the system group
```

Better to check if the trajectory PBC has been totally fixed by xtc visualization software, such as [VMD](http://www.ks.uiuc.edu/Research/vmd/).

### Calculation

Then start MM/PB-SA calculation.

``` bash
# Firstly, add s_mmpbsa folder to $PATH.
# Start s_mmpbsa, and input as follow (do not include comments)
md.tpr
1 # load xtc file
md_pbc.xtc # if not PBC-fixed, click "return" and use default md.xtc
2 # load ndx file
[return] # default index.ndx
0 # go to next step (Trajectory Parameters)
1 # select receptor group
[protein group number]
2 # select ligand group
[ligand group number]
5 # set time interval, usually analysis per 1 ns
1
0 # go to next step (MM/PB-SA Parameters)
# Usually no need to change. The PB and SA parameters could be modified by 8 and 9
0 # go to next step (start calculation)
[return] # use default system name or input your name
# Wait for calculation finish
1 # view summary
[return] # output energy summary with default file name or input your name
10 # output other infomation by 2-9
-1 # view pdb file with residue-wised INVERSED binding energy filled in B-factor column
0 # return
-10 # return
-10 # return
-10 # exit s_mmpbsa program
```

## Download
Release file: https://github.com/supernova4869/s_mmpbsa/releases, where "s_mmpbsa.exe" and "s_mmpbsa" are s_mmpbsa executable files on Windows and Linux operation systems, respectively.

## Citation
The `s_mmpbsa` program should be properly cited if the work will be published. 

Currently, s_mmpbsa is still in-develop. If you want to utilize s_mmpbsa in your work, please cite the program as following:

```
Jiaxing Zhang, s_mmpbsa, Version [current version], https://github.com/supernova4869/s_mmpbsa (accessed on yy-mm-dd)
```

After the detailed paper about s_mmpbsa is published (if fortunately), please cite the corresponding paper instead of the web page here.

## About developer
Dr. Jiaxing Zhang (Contact: Jiaxing_Zhang@outlook.com, Tianjin University)

If you encountered any difficulty while using s_mmpbsa, or you found any bugs, or you have any suggestion on improving s_mmpbsa, please E-mail me or join my QQ group 864191465 to describe.

## New Folder (?
- Multi-threading
- Automatic Periodic boundary conditions (PBC) fix
- Add supporting of other PBSA solvers, e.g., AFMPB, Delphi2
