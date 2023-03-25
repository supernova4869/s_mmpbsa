# super_mmpbsa
super_mmpbsa: Supernova's tool of binding free energy calculation for Gromacs trajectory, using molecular mechanics Poisson-Boltzmann surface area (MM-PBSA) method.

Super_mmpbsa can be freely used for both academic and commerical purposes.

## Introduction

MM/PB-SA method is the most popular method to calculate binding free energy, especially for biological systems. However, as a widely-used MD program, Gromacs has not officially support MM/PB-SA calculation. Although there have been numerous programs that can calculate binding free energy with GROMACS trajectory, all of them have some limitations at different aspects, for example: (1) Difficult to use and install (2) Not support newer version of Gromacs (3) Too slow (4) Not cross-platform. Instead, super_mmpbsa provides a convinent interface (like [Multiwfn](http://sobereva.com/multiwfn/)) to calculate binding free energy for GROMACS trajectory. 

## Features of super_mmpbsa

- Freely available.
- No need for preparing runtime, only needs Gromacs program on linux system. (in contrast to gmx_MMPBSA.py, super_mmpbsa is developed with in Rust).
- Interactive operation, no need to write parameter files. Also, user can write shell script to invoke super_mmpbsa for batch use.
- Very fast. Due to the efficency of rust program.
- Considers electric screening effect, as "CHIN. PHYS. LETT. 2021, 38(1), 018701" describes.

## Download
Release file: https://github.com/supernovaZhangJiaXing/super_mmpbsa/releases, where "super_mmpbsa.exe" and "super_mmpbsa" are super_mmpbsa executable files on Windows and Linux operation systems, respectively.

## Citation
Super_mmpbsa should be properly cited if the work will be published. 

Currently, super_mmpbsa is still in-develop. If you want to utilize super_mmpbsa in your work, please cite the program as following:

```
Jiaxing Zhang, Super_mmpbsa, Version [current version], https://github.com/supernovaZhangJiaXing/super_mmpbsa (accessed on yy-mm-dd)
```

After the detailed paper about super_mmpbsa is published (if fortunately), please cite the corresponding paper instead of the web page here.

## About developer
Dr. Jiaxing Zhang (Contact: Jiaxing_Zhang@outlook.com, Tian Jin University)

If you encountered any difficulty while using super_mmpbsa, or you found some bugs, or you have any suggestion on improving super_mmpbsa, please send E-mail to me or join my QQ group 864191465 to describe.

## New Folder (?
- Multi-threading
- Automatic Periodic boundary conditions (PBC) fix
- Add supporting of other PBSA solvers, e.g., AFMPB, Delphi2
