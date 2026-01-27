====
Introduction
====

This document introduces the background, principles and basic concepts of s_mmpbsa to help users understand the working principle and application scenarios of this tool.

Introduction to MM-PBSA Method
----------------

**MM-PBSA** (Molecular Mechanics/Poisson-Boltzmann Surface Area) is a widely used method for calculating binding free energy of biomolecules. This method combines molecular mechanics (MM) and continuum solvent model (PB-SA), which can quickly and accurately predict the interaction strength between biomolecules.

The basic principle of the MM-PBSA method is to evaluate the binding strength between molecules by calculating the free energy change before and after binding. Specifically, the binding free energy (ΔG) can be expressed as:

.. math::

   \Delta G_{binding} = \Delta G_{complex} - (\Delta G_{receptor} + \Delta G_{ligand})

Where the free energy (G) of each molecule consists of the following components:

.. math::

   G = E_{MM} + G_{solv} - T\Delta S

- **E_{MM}**: Molecular mechanics energy, including bond energy, angle energy, dihedral energy, van der Waals energy and electrostatic energy
- **G_{solv}**: Solvation free energy, including polar solvation energy (calculated via Poisson-Boltzmann equation) and non-polar solvation energy (calculated via surface area)
- **TΔS**: Entropy contribution term, usually calculated via normal mode analysis

Why Choose s_mmpbsa?
-------------------

Although Gromacs is a widely used molecular dynamics simulation software, it does not officially support MM-PBSA calculations. There are many MM-PBSA tools on the market that can handle Gromacs trajectories, but most of them have the following limitations:

1. **Complex installation and usage**
2. **Not supporting new versions of Gromacs**
3. **Slow calculation speed**
4. **Not cross-platform**

In contrast, s_mmpbsa offers the following advantages:

- **Simple and easy to use**: Interactive operation interface, no need to write complex parameter files
- **Efficient calculation**: Developed in Rust language, with fast calculation speed
- **Cross-platform**: Supports Windows and Linux operating systems
- **Rich features**: Supports charge screening effect and conformational entropy calculation
- **Easy integration**: Can be called via scripts and supports batch processing

s_mmpbsa's Basic Workflow
--------------------- 

s_mmpbsa's workflow mainly includes the following steps:

1. **Input processing**: Read Gromacs tpr, xtc and ndx files
2. **Trajectory processing**: Process molecular dynamics trajectories, including extracting coordinates and handling periodic boundary conditions
3. **MM calculation**: Calculate molecular mechanics energy (bond energy, van der Waals energy, electrostatic energy, etc.)
4. **PB-SA calculation**: Calculate solvation free energy (polar and non-polar)
5. **Entropy calculation**: Calculate conformational entropy contribution (optional)
6. **Result analysis**: Generate binding free energy reports and visualization results

Application Scenarios
--------

s_mmpbsa is suitable for the following research scenarios:

1. **Drug design**: Evaluate the binding strength between drug molecules and targets, guide drug optimization
2. **Protein-protein interactions**: Study the stability and interaction mechanism of protein complexes
3. **Enzyme-substrate interactions**: Analyze binding free energy changes in enzyme-catalyzed reactions
4. **Mutation effect prediction**: Evaluate the contribution of key residues in proteins to binding through alanine scanning
5. **Molecular docking result validation**: Provide more accurate binding energy predictions for molecular docking results

Theoretical Innovations
----------

s_mmpbsa has made the following improvements on the basis of the traditional MM-PBSA method:

1. **Charge screening effect**: Considered the charge screening effect in biomolecular environments, improving the accuracy of polar interaction calculations (Reference: J. Chem. Inf. Model. 2021, 61, 2454)

2. **Conformational entropy calculation**: Implemented an efficient conformational entropy calculation method, providing more comprehensive thermodynamic information for binding free energy prediction (Reference: J. Chem. Phys. 2017, 146, 124124)

3. **Parallel computing optimization**: Through the parallel features of the Rust language, significantly improved computing efficiency, especially for large biomolecular systems

4. **Result visualization**: Provided rich result analysis and visualization functions, facilitating users to understand and interpret calculation results

License Information
----------

s_mmpbsa follows the LGPL license and can be used free of charge for academic research purposes. If you use s_mmpbsa in a commercial environment, please ensure that you comply with the relevant provisions of the LGPL license.

Citing s_mmpbsa
-----------

If you use s_mmpbsa in your research work, please cite it in the following format:

.. code-block:: text

   Jiaxing Zhang, s_mmpbsa, Version [your version], https://github.com/supernova4869/s_mmpbsa (accessed on yy-mm-dd)

We are preparing a detailed academic paper on s_mmpbsa, and please cite the corresponding paper after its publication.