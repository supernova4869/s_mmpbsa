sphinx-quickstart on 2025-09-01.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to s_mmpbsa Documentation
======================

**s_mmpbsa** is an efficient tool for calculating biomolecular binding free energy, specifically designed for analyzing Gromacs trajectories using the Molecular Mechanics Poisson-Boltzmann Surface Area (MM/PB-SA) method.

.. toctree::
   :maxdepth: 2
   :caption: Content Overview

   introduction_en
   installation_en
   quick_start_en
   usage_en
   features
   examples
   api_reference
   faq_en
   contribution
   authors

Introduction
----

**s_mmpbsa** provides a user-friendly interface (similar to `Multiwfn`) for calculating binding free energy from Gromacs trajectories. Compared to other similar tools, it offers advantages such as easy installation, efficient operation, cross-platform compatibility, and supports advanced features like charge screening effects and conformational entropy calculation.

Main Features
--------

- **MD Simulation Binding Energy Calculation**: Calculate binding free energy between biomolecules from molecular dynamics simulation results
- **Molecular Docking Result Rescoring**: Provide more accurate binding energy predictions for molecular docking results
- **Protein-Ligand Complex Alanine Scanning**: Analyze the contribution of key residues in proteins to binding

Features
--------

- Open source and free, under the LGPL license
- Minimal environment dependencies, only requiring Gromacs on Linux systems; Python environment needed for plotting functionality
- Developed in Rust language for excellent performance
- Interactive operation, no need to write parameter files
- Considers charge screening effects, as described in literature [J. Chem. Inf. Model. 2021, 61, 2454]
- Considers conformational entropy, as described in literature [J. Chem. Phys. 2017, 146, 124124]
- Can store analysis results for further reproducible analysis

Getting Started
--------

Please refer to the :doc:`installation_en` chapter to install s_mmpbsa, and then check the :doc:`quick_start_en` chapter to learn the basic usage流程.

For detailed usage instructions, please refer to the :doc:`usage_en` chapter, which contains usage methods and examples of various functions.

Getting Help
--------

If you encounter any problems during use or have any suggestions for improvement, please contact the developer or join the QQ group:

- **Developer**: Dr. Jiaxing Zhang (zhangjiaxing7137@tju.edu.cn, Tianjin University)
- **QQ Group**: 864191465

Citation
----

If you use s_mmpbsa in your research work, please cite it in the following format:

.. code-block:: text

   Jiaxing Zhang, s_mmpbsa, Version [your version], https://github.com/supernova4869/s_mmpbsa (accessed on yy-mm-dd)

.. note::
   When a detailed paper on s_mmpbsa is published, please cite the corresponding paper instead of this webpage.

Index and Tables
----------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`