========
Usage
========

This document details the various functions and usage methods of s_mmpbsa to help you deeply understand and fully utilize this tool for binding free energy calculations.

Command Line Parameters
----------

s_mmpbsa supports the following command line parameters:

.. code-block:: bash
   
   # Basic usage
   s_mmpbsa [options] [tpr file]
   
   # Options
   --version   Display version information

Interactive Command Line Interface
----------------

s_mmpbsa's interactive command line interface is divided into the following main parts:

1. **File loading**: Load tpr, xtc and ndx files
2. **Trajectory parameter setting**: Set receptor, ligand groups, time interval, etc.
3. **MM/PB-SA parameter setting**: Set various parameters related to calculation
4. **Calculation execution**: Execute MM/PB-SA calculation
5. **Result analysis**: Analyze and visualize calculation results

The operation methods and parameter settings of each part are described in detail below.

File Loading
--------

After starting s_mmpbsa, you first need to load the necessary input files:

- **tpr file**: Gromacs topology parameter file, containing system topology information and atomic parameters
- **xtc file**: Gromacs trajectory file, containing system coordinate information
- **ndx file**: Gromacs index file, containing system grouping information

Example operations for loading files:

.. code-block:: bash
   
   # Input tpr file path (can be absolute or relative path)
   md.tpr
   
   # Choose to load xtc file
   1
   md_centered.xtc  # Input xtc file path, press Enter to use default md.xtc
   
   # Choose to load ndx file
   2
   index.ndx  # Input ndx file path, press Enter to use default index.ndx
   
   # Proceed to next step
   0

Trajectory Parameter Setting
------------

After loading files, you need to set trajectory parameters, mainly including:

- **Receptor group**: Choose which group to use as the receptor
- **Ligand group**: Choose which group to use as the ligand
- **Time interval**: Set the time interval for analysis
- **Skipped frames**: Set the number of frames to skip before starting analysis
- **End analysis time point**: Set the time point to end analysis

Example operations for setting trajectory parameters:

.. code-block:: bash
   
   # Select receptor group
   1
   1  # Input the number of the receptor group, for example 1 represents Protein group
   
   # Select ligand group
   2
   13  # Input the number of the ligand group, for example 13 represents ligand
   
   # Set start time
   3
   0  # Input the number of frames to skip, default is 0
   
   # Set end time
   4
   0  # Input end time point, 0 means analyze until the end of the trajectory
   
   # Set time interval (unit: ns)
   5
   1  # Input time interval, for example 1 means analyze once every 1ns
   
   # Proceed to next step
   0

MM/PB-SA Parameter Setting
--------------

Next, set the relevant parameters for MM/PB-SA calculation. s_mmpbsa provides multiple parameter setting options:

.. code-block:: bash
   
   # Display current parameter settings
   1
   
   # Set temperature (unit: K)
   2
   298.15  # Input temperature, default is 298.15K
   
   # Set NaCl concentration (unit: mol/L)
   3
   0.15  # Input KCl concentration, default is 0.15mol/L
   
   # Set salt bridge search distance (unit: Å)
   4
   4.0  # Input salt bridge search distance, default is 4.0Å
   
   # Set hydrogen bond search distance (unit: Å)
   5
   3.5  # Input hydrogen bond search distance, default is 3.5Å
   
   # Set van der Waals cutoff distance (unit: Å)
   6
   14.0  # Input van der Waals cutoff distance, default is 14.0Å
   
   # Set number of parallel cores for MM calculation
   7
   4  # Input number of parallel cores, default is the number of system CPU cores
   
   # Set PB parameters
   8
   # Enter PB parameter setting submenu (see below for details)
   
   # Set SA parameters
   9
   # Enter SA parameter setting submenu (see below for details)
   
   # Proceed to next step
   0

PB Parameter Setting
~~~~~~~~~

In the PB parameter setting submenu, you can set the following parameters:

.. code-block:: bash
   
   # Display current PB parameter settings
   1
   
   # Set solvent dielectric constant
   2
   78.54  # Input solvent dielectric constant, default is 78.54
   
   # Set solute dielectric constant
   3
   1.0  # Input solute dielectric constant, default is 1.0
   
   # Set grid spacing (unit: Å)
   4
   0.5  # Input grid spacing, default is 0.5Å
   
   # Set APBS executable file path
   5
   /usr/local/bin/apbs  # Input APBS executable file path, press Enter to use built-in path
   
   # Return to previous menu
   0

SA Parameter Setting
~~~~~~~~~

In the SA parameter setting submenu, you can set the following parameters:

.. code-block:: bash
   
   # Display current SA parameter settings
   1
   
   # Set surface tension (unit: kJ/(mol·Å²))
   2
   0.0379  # Input surface tension, default is 0.0379 kJ/(mol·Å²)
   
   # Set non-polar solvation parameter (unit: kJ/(mol·Å³))
   3
   0.0  # Input non-polar solvation parameter, default is 0.0 kJ/(mol·Å³)
   
   # Set SAS calculation method (0: Shrake-Rupley, 1: MSMS)
   4
   0  # Input SAS calculation method, default is 0 (Shrake-Rupley)
   
   # Return to previous menu
   0

Executing Calculation
--------

After setting up, start executing the MM/PB-SA calculation. Before calculation, you need to input the system name:

.. code-block:: bash
   
   # Input system name
   system  # Input system name, default is system
   
   # A progress bar and current energy value will be displayed during calculation
   # After calculation is complete, you will automatically enter analysis mode

Result Analysis
--------

After calculation is complete, you can analyze the results. s_mmpbsa provides multiple analysis options:

.. code-block:: bash
   
   # Generate pdb file containing residue binding energy information
   -1
   0  # Input time point, 0 means average value
   
   # View result summary
   1
   
   # Output energy change data over time
   2
   
   # Output residue binding energy at a specific time point
   3
   0  # Input time point, 0 means average value
   1  # Choose to output residues within 3Å range (0: all residues, 1: 3Å range, 2: 5Å range, 3: 10Å range)
   
   # Output binding energy of ligand atoms
   4
   
   # Output hydrogen bond and salt bridge information
   5
   
   # Output interaction energy matrix
   6
   
   # Exit program
   0

Using Analysis Mode
------------

s_mmpbsa also provides a special analysis mode, which can re-analyze already calculated results without re-calculating:

.. code-block:: bash
   
   # Start analysis mode
   s_mmpbsa
   a  # Input 'a' in the interactive interface to enter analysis mode
   
   # Input working directory path
   ./results  # Input the directory containing .sm result files, default is current directory
   
   # Input temperature (unit: K)
   298.15  # Input temperature, default is 298.15K
   
   # Input system name
   system  # Input the system name used during previous calculation, default is system
   
   # Subsequent analysis operations are the same as after normal calculation completion

Alanine Scanning
----------

s_mmpbsa also supports alanine scanning function, which can systematically mutate protein residues to alanine and calculate the binding free energy change before and after mutation.

Steps for performing alanine scanning:

.. code-block:: bash
   
   # Prepare working directory
   mkdir -p ala_scan
   cd ala_scan
   
   # Copy necessary input files
   cp ../md.tpr ../md_centered.xtc ../index.ndx .
   
   # Execute alanine scanning
   s_mmpbsa md.tpr
   1
   md_centered.xtc
   2
   index.ndx
   0
   1
   1  # Select receptor group (protein)
   2
   13  # Select ligand group
   5
   1
   0
   0
   system  # System name
   -1  # Generate pdb file (optional)
   0  # Exit analysis
   a  # Enter analysis mode
   .  # Use current directory
   298.15  # Temperature
   system  # System name
   0  # Exit analysis
   
   # Now you can view the results of alanine scanning

Interpretation of Calculation Results
------------

s_mmpbsa's calculation results mainly include the following energy terms:

- **ΔG_bind**: Total binding free energy
- **ΔH**: Enthalpy change
- **TΔS**: Entropy contribution (Note: s_mmpbsa does not directly calculate entropy at present, this value is usually set to 0 or estimated through other methods)
- **ΔE_vdw**: Van der Waals interaction energy
- **ΔE_elec**: Electrostatic interaction energy
- **ΔG_polar**: Polar solvation free energy
- **ΔG_nonpolar**: Non-polar solvation free energy

A larger negative value of binding free energy indicates stronger binding. Usually, ΔG_bind < -10 kJ/mol indicates strong binding.

Notes
--------

When using s_mmpbsa, you need to pay attention to the following points:

1. **Trajectory quality**: Ensure good trajectory quality, with correct PBC handling, centering and fitting operations.

2. **Index file**: Ensure that the index file contains correct receptor and ligand groups.

3. **Parameter selection**: For different systems, parameters may need to be adjusted to obtain more accurate results.

4. **Parallel computing**: Setting an appropriate nkernels value in settings.ini can utilize multi-core CPU to accelerate calculation.

5. **Result verification**: It is recommended to compare with experimental data or results from other calculation methods to verify the reliability of calculation results.

More Information
--------

- :doc:`quick_start_en`：Quick Start Guide
- :doc:`installation_en`：Installation Instructions
- :doc:`api`：API Documentation (for developers)
- :doc:`faq`：Frequently Asked Questions