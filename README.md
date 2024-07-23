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
- Supports molecular docking results (by vina or DSDP) analysis.
- Can perform alanine scanning based on both MD and molecular docking results.

## Usage
For Ubuntu system, maybe user should run the following commands to avoid `cc` error.
```
sudo apt update
sudo apt install build-essential
```

Currently s_mmpbsa has supported to fix PBC conditions and write trajectory to `MMPBSA_[name].xtc`. However, it is still better to re-check if the trajectory PBC has been totally fixed by xtc visualization software, such as [VMD](http://www.ks.uiuc.edu/Research/vmd/).

``` bash
# Firstly, add s_mmpbsa folder to $PATH.
# Start s_mmpbsa, and input as follow (support # comments, but not recommended to input comments)
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
2 # output energy by time
[return]
3 # output energy by residue
1 # write residues within 3 A (also try other options)
[return]
4 # output other infomations
-1 # view pdb file with residue-wised INVERSED binding energy filled in B-factor column
0 # exit s_mmpbsa program
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
Dr. Jiaxing Zhang (Contact: zhangjiaxing7137@tju.edu.cn, Tianjin University)

If you encountered any difficulty while using s_mmpbsa, or you found any bugs, or you have any suggestion on improving s_mmpbsa, please E-mail me or join my QQ group 864191465 to describe.

## New Folder (?
- Add support of other PBSA solvers, e.g., AFMPB, Delphi2, and also built-in LPBE solver
- Add figures: plotting data, 2D interaction plots like Ligplot+, 3D plots scripts with PyMOL
- Multi-threading
