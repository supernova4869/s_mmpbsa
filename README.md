# s_mmpbsa
s_mmpbsa: Supernova's tool of binding free energy calculation for Gromacs trajectory, using molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method.

The s_mmpbsa program follows LGPL license, and can be freely used for academic purposes.

## Introduction
MM/PB-SA method is the most popular method to rapidly calculate binding free energy, especially for biological systems. However, as a widely-used MD program, Gromacs has not officially support MM/PB-SA calculation. Although there have been numerous programs that can calculate binding free energy with GROMACS trajectory, most of them have some limitations at different aspects, for example: (1) Difficult to use and install (2) Lack of support for new version of Gromacs (3) Too slow (4) Not cross-platform. Instead, s_mmpbsa provides a convinent interface (like [Multiwfn](http://sobereva.com/multiwfn/)) to calculate binding free energy for GROMACS trajectory. 

## Features of s_mmpbsa
- Open source and freely available.
- Less need for running environment preparation, only needs Gromacs program when running on linux system, and Python environment for plotting. (In contrast to other programs such as gmx_MMPBSA.py, s_mmpbsa is developed with Rust).
- Interactive operation, no need to write parameter files. Also, user can write shell script to invoke s_mmpbsa for batch use.
- Very fast. Due to the efficency of rust program.
- Considers electric screening effect, as "J. Chem. Inf. Model. 2021, 61, 2454".
- Considers interaction entropy, as "J. Chem. Phys. 2017, 146, 124124".
- Supports molecular docking results (by Autodock vina or DSDP) analysis (next small version).
- Can perform alanine scanning based on both MD and molecular docking results (next small version).
- Could store analyzation results for further reproducable analyzation.

## Requirement
The matplotlib python package is essential during analyzation if plotting figures.
On Debian/Ubuntu/Linux, run:
```
sudo apt -y install python3-matplotlib build-essential
```
On CentOS/Rocky, run:
```
sudo dnf -y install python3-matplotlib
```

## Usage
Although s_mmpbsa supports fixing PBC conditions to trajectory `_MMPBSA_[name].xtc`, it is still recommended to comfirm that the trajectory has been correct, using xtc visualization software such as [VMD](http://www.ks.uiuc.edu/Research/vmd/).

``` bash
# Typical calculation mode
# Firstly, add s_mmpbsa folder to $PATH.
# Start s_mmpbsa, and input as follow (support # comments, but not recommended and usually no need to input with comments)
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
10 # Do Alanine scanning
1 # select residues within 4 A
# Other options usually no need to change. The PB and SA parameters could be modified by 8 and 9
0 # go to next step (start calculation)
[return] # use default system name or input your name
# Wait for calculation finish
-1 # write pdb file with residue-wised INVERSED binding energy filled in B-factor column
 # input the time point (default average)
1 # view summary
2 # output energy by time
3 # output energy by residue
 # input the time point (default average)
1 # write residues within 3 A (also try other options)
4 # output energy by ligand atoms
0 # exit s_mmpbsa program
```

```bash
# Analyzation mode
# Firstly, add s_mmpbsa folder to $PATH.
# Start s_mmpbsa, and input as follow (support # comments, but not recommended and usually no need to input with comments)
a # analyzation mode
[return] # same path as last opened tpr file
[return] # set temperature 298.15 K
[return] # set system name, same as the previous run
-1 # write pdb file with residue-wised INVERSED binding energy filled in B-factor column
 # input the time point (default average)
1 # view summary
2 # output energy by time
3 # output energy by residue
 # input the time point (default average)
1 # write residues within 3 A (also try other options)
4 # output energy by ligand atoms
0 # exit s_mmpbsa program
```

The data was generated with .csv format and the pml files could be loaded by PyMOL.
```bash
pymol MMPBSA_binding_energy__system_avg.pml
```

## Download
Release file: https://github.com/supernova4869/s_mmpbsa/releases, where "s_mmpbsa.exe" and "s_mmpbsa" are s_mmpbsa executable files on Windows and Linux operation systems, respectively.

## Citation
The `s_mmpbsa` program should be properly cited if the work will be published. 

Currently, s_mmpbsa is still in-develop. If you want to utilize s_mmpbsa in your work, please cite the program as following:

```
Jiaxing Zhang, s_mmpbsa, Version [your version], https://github.com/supernova4869/s_mmpbsa (accessed on yy-mm-dd)
```

After the detailed paper about s_mmpbsa is published (if fortunately), please cite the corresponding paper instead of the web page here.

## About developer
Dr. Jiaxing Zhang (Contact: zhangjiaxing7137@tju.edu.cn, Tianjin University)

If you encountered any difficulty while using s_mmpbsa, or you found any bugs, or you have any suggestion on improving s_mmpbsa, please E-mail me or join my QQ group 864191465 to describe.

## New Folder (?
- Add support of other PBSA solvers, e.g., Delphi2, and also built-in LPBE solver
