# s_mmpbsa
s_mmpbsa: Supernova's tool of binding free energy calculation for Gromacs trajectory, using molecular mechanics Poisson-Boltzmann surface area (MM/PB-SA) method.

The s_mmpbsa program follows LGPL license, and can be freely used for academic purposes.

## Introduction
MM/PB-SA method is the most popular method to rapidly calculate binding free energy, especially for biological systems. However, as a widely-used MD program, Gromacs has not officially support MM/PB-SA calculation. Although there have been numerous programs that can calculate binding free energy with GROMACS trajectory, most of them have some limitations at different aspects, for example: (1) Difficult to use and install (2) Lack of support for new version of Gromacs (3) Too slow (4) Not cross-platform. Instead, s_mmpbsa provides a convenient interface (like [Multiwfn](http://sobereva.com/multiwfn/)) to calculate binding free energy for GROMACS trajectory. 

## Features of s_mmpbsa
- Open source and freely available.
- Less need for running environment preparation, only needs Gromacs program when running on linux system, and Python environment for plotting.
- In contrast to other programs, s_mmpbsa is developed with Rust.
- Interactive operation, no need to write parameter files. Also, user can write shell script to invoke s_mmpbsa for batch use.
- Considers electric screening effect, as "J. Chem. Inf. Model. 2021, 61, 2454".
- Considers conformational entropy, as "J. Chem. Phys. 2017, 146, 124124".
- Could store analyzation results for further reproducable analyzation.

## Main function
- Binding energy calculation from MD simulation.
- Molecular docking results rescoring.
- Alanine scanning of protein-ligand complex.

## Requirement

### Basic requirements
- Gromacs: The gromacs program is needed.
- Matplotlib: (Optional) The matplotlib python package is essential during analyzation if plotting figures.
- APBS: (Optional) The default PBSA kernel (already built-in, but it is also supported to use other version of APBS programs).

On Debian/Ubuntu/Linux, run:
```
sudo apt -y install python3-matplotlib build-essential python-pip
```
On CentOS/Rocky, run:
```
sudo dnf -y install python3-matplotlib python-pip
```

### Special requirements of molecular docking rescoring:
- PyMOL is an optional software to plot the B-factor colored structure.
- Gaussian is an optional software to do DFT calculations for RESP atom charge calculation.
- Multiwfn is an optional program (already built-in) to fit RESP atom charge.

## Usage
Although s_mmpbsa supports fixing PBC conditions to trajectory `_MMPBSA_[name].xtc`, it is still recommended to comfirm that the trajectory has been correct, using xtc visualization software such as [VMD](http://www.ks.uiuc.edu/Research/vmd/).

### MD Binding energy calculation:
``` bash
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

### Docking Rescoring function:
``` bash
receptor.pdbqt
2 # load ligand file
[return] # default DSDP.pdbqt
3 # load flexible residues file (flexible docking only)
[flexible residues file name]
# other functions are used to prepare the complex system
# if the ligand is charged, do NOT forget to change option 8
0 # go to next step
1 # select start model
[start model number]
2 # select end model
[end model number]
0 # go to next step (MM/PB-SA Parameters)
# Other options usually no need to change. The PB and SA parameters could be modified by 8 and 9
0 # go to next step (start calculation)
[return] # use default system name or input your name
# Wait for calculation finish
-1 # write pdb file with residue-wised INVERSED binding energy filled in B-factor column
 # input the time point (default average)
1 # view summary (Here the ΔG and TΔS values are useless)
2 # output energy by time
3 # output energy by residue
 # input the time point (default average)
1 # write residues within 3 A (also try other options)
4 # output energy by ligand atoms
0 # exit s_mmpbsa program
```

### Alanine scanning:
```bash
# At MM/PB-SA Parameters page
2 # Do alanine scanning
1 # Select mutation residues by layers (ACS Catal. 2024, 14, 15, 11447–11456)
0 # Start calculation
...
```
The results will contain energy terms of both wild type and mutants.

### Use Analyzation mode:
```bash
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

The data was generated with .csv format and plotted as figures; The pdb files with inversed binding energy filled in B-factors column are drawn as png figures by PyMOL (if usable).

## Download
Release file: https://github.com/supernova4869/s_mmpbsa/releases/latest, where "s_mmpbsa.exe" and "s_mmpbsa" are s_mmpbsa executable files on Windows and Linux operation systems, respectively.

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
