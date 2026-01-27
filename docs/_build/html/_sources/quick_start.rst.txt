========
Quick Start
========

This document provides the basic usage process of s_mmpbsa to help you quickly get started with this tool for binding free energy calculations.

Starting s_mmpbsa
------------

After installation, you can start s_mmpbsa in the following ways:

.. code-block:: bash
   
   # Run directly in the command line
   s_mmpbsa
   
   # Or specify the tpr file path directly
   s_mmpbsa md.tpr

After starting, you will see the welcome message of s_mmpbsa and then enter the interactive interface.

Basic Workflow
------------

s_mmpbsa's basic workflow includes the following steps:

1. Load input files (tpr, xtc and ndx files)
2. Set trajectory parameters (select receptor and ligand groups)
3. Set MM-PBSA parameters
4. Execute calculation
5. Analyze results

We will detail the operation method of each step below.

Loading Input Files
------------

After starting s_mmpbsa, you first need to load the necessary input files:

.. code-block:: bash
   
   # Input tpr file path
   md.tpr
   
   # Load xtc file (option 1)
   1
   md_centered.xtc  # If the trajectory has been processed with PBC, you can input it directly; otherwise press Enter to use the default md.xtc
   
   # Load ndx file (option 2)
   2
   index.ndx  # Press Enter to use the default index.ndx
   
   # Proceed to the next step (option 0)
   0

Setting Trajectory Parameters
------------

After loading input files, you need to set trajectory parameters, mainly selecting receptor and ligand groups:

.. code-block:: bash
   
   # Select receptor group (option 1)
   1
   [Select the number of the receptor group, for example 1 represents Protein]
   
   # Select ligand group (option 2)
   2
   [Select the number of the ligand group, for example 13 represents ligand]
   
   # Set time interval (option 5), usually analyze once every 1ns
   5
   1  # Time interval, unit is ns
   
   # Proceed to the next step (option 0)
   0

Setting MM-PBSA Parameters
--------------

Next, set the relevant parameters for MM-PBSA calculation:

.. code-block:: bash
   
   # Under normal circumstances, you can use the default parameters
   # If you need to modify PB parameters, you can select option 8
   # If you need to modify SA parameters, you can select option 9
   
   # Proceed to the next step (option 0)
   0

Executing Calculation
--------

After setting up, start executing the calculation:

.. code-block:: bash
   
   # Input system name (default is system)
   [Press Enter to use the default name or input a custom name]
   
   # Wait for calculation to complete
   # A progress bar and current energy value will be displayed during calculation

Analyzing Results
--------

After calculation is complete, you can analyze the results:

.. code-block:: bash
   
   # Generate pdb file containing residue binding energy information (option -1)
   -1
   [Press Enter to use the default time point (average value) or input a specific time point]
   
   # View result summary (option 1)
   1
   
   # Output energy change data over time (option 2)
   2
   
   # Output residue binding energy at a specific time point (option 3)
   3
   [Press Enter to use the default time point (average value) or input a specific time point]
   1  # Select to output residues within 3Å range
   
   # Output binding energy of ligand atoms (option 4)
   4
   
   # Exit program (option 0)
   0

Using Analysis Mode
------------

s_mmpbsa also provides a special analysis mode, which can re-analyze already calculated results:

.. code-block:: bash
   
   # Start analysis mode
   s_mmpbsa
   a  # Input 'a' to enter analysis mode
   
   # Input working directory path (default is current directory)
   [Press Enter to use current directory or input the directory containing .sm result files]
   
   # Input temperature (default is 298.15K)
   [Press Enter to use default temperature or input custom temperature]
   
   # Input system name (default is system)
   [Press Enter to use default name or input the system name used during previous calculation]
   
   # Subsequent analysis operations are the same as after normal calculation completion

Example: Calculating Protein-Ligand Binding Energy
-------------------------

The following is a complete example of calculating protein-ligand binding energy:

.. code-block:: bash
   
   # Start s_mmpbsa and load files
   s_mmpbsa
   md.tpr
   1
   md_centered.xtc
   2
   index.ndx
   0
   
   # Set trajectory parameters
   1
   1  # Assume 1 is the Protein group
   2
   13  # Assume 13 is the ligand group
   5
   1
   0
   
   # Set MM-PBSA parameters (use default values)
   0
   
   # Execute calculation
   protein_ligand  # System name
   
   # Analyze results
   -1
   1
   2
   3
   1
   4
   0

Usage Tips
--------

1. **Trajectory preparation**: Before calculation, it is recommended to use Gromacs' trjconv tool to process the trajectory, including removing PBC, centering and fitting operations, to obtain better calculation results.

2. **Index file**: Ensure that the index file contains correct receptor and ligand groups. If there is no ready-made index file, you can use Gromacs' make_ndx tool to create one.

3. **Time interval**: For long MD simulations, you can appropriately increase the time interval to reduce the calculation amount. Usually analyzing once every 1-2ns can obtain good statistical results.

4. **Parallel computing**: Setting an appropriate nkernels value in settings.ini can utilize multi-core CPU to accelerate calculation.

5. **Result visualization**: The generated pdb file can be opened with software such as PyMOL to view the distribution of residue binding energy (colored by B factor).

Frequently Asked Questions
----------

### How to handle large systems?

For large systems, you can try the following optimization measures:

- Increase the time interval to reduce the number of analyzed frames
- Increase the nkernels value to utilize more CPU cores
- Use smaller cutoff distance (by modifying the r_cutoff parameter)

### How to improve calculation accuracy?

Methods to improve calculation accuracy include:

- Ensure good trajectory quality and correct PBC handling
- Increase the number of sampling points, i.e., reduce the time interval
- Adjust PB parameters, such as grid size, solvent dielectric constant, etc.

### How to interpret calculation results?

A larger negative value of binding free energy indicates stronger binding. Usually, the calculation results will give the following energy terms:

- ΔG_bind: Total binding free energy
- ΔH: Enthalpy change
- TΔS: Entropy contribution
- ΔE_vdw: Van der Waals interaction energy
- ΔE_elec: Electrostatic interaction energy
- ΔG_polar: Polar solvation free energy
- ΔG_nonpolar: Non-polar solvation free energy

For more detailed information, please refer to the :doc:`usage` chapter.