====
Installation
====

This document introduces the installation steps, environment requirements and configuration methods of s_mmpbsa to help users quickly set up the running environment.

System Requirements
--------

### Operating System

s_mmpbsa supports the following operating systems:

- **Windows**: Windows 7/8/10/11
- **Linux**: Ubuntu, Debian, CentOS, Rocky and other mainstream Linux distributions

### Hardware Requirements

- **Processor**: Multi-core processor (4 cores or more recommended)
- **Memory**: At least 4GB RAM (8GB or more recommended for large systems)
- **Disk space**: At least 500MB available space

Software Dependencies
--------

### Basic Dependencies

s_mmpbsa's core functionality does not require any external program:

- **Gromacs**: Not required. The Gromacs tools s_mmpbsa needs (`gmx dump`, `gmx trjconv`, `gmx convert-tpr` and `gmx make_ndx`) are provided by the Rust port (gmx-rs-tools) compiled into the s_mmpbsa executable. Run input files written by Gromacs 2021 and later are read directly; for files from older Gromacs versions an installed `gmx` program is used, which can be avoided by re-converting the file once with `gmx convert-tpr -s old.tpr -o new.tpr`.

Optional dependencies:

- **Matplotlib**: Used to draw result charts. If you need to use s_mmpbsa's analysis and plotting functions, it needs to be installed.
- **APBS**: Not required. Poisson-Boltzmann and surface-area calculations are performed by the Rust-ported APBS solver (apbs-rs), compiled directly into the s_mmpbsa executable.
- **PyMOL**: Used to draw B-factor colored structures.

Installation Methods
--------

### Windows System

1. **Download s_mmpbsa**
   
   Download the latest version of the Windows executable file from the GitHub release page:
   
   .. code-block:: powershell
      
      # Download s_mmpbsa.exe from GitHub
      # Visit: https://github.com/supernova4869/s_mmpbsa/releases/latest
   
2. **Add to System Path**
   
   Add the folder containing s_mmpbsa.exe to the system environment variable PATH so that s_mmpbsa can be run from any location.
   
3. **Install Optional Dependencies (If Analysis Features Are Needed)**
   
   .. code-block:: powershell
      
      # Install matplotlib
      pip install matplotlib

### Linux System

1. **Download s_mmpbsa**
   
   .. code-block:: bash
      
      # Download the latest version from GitHub
      wget https://github.com/supernova4869/s_mmpbsa/releases/latest/download/s_mmpbsa
      
      # Add execution permission
      chmod +x s_mmpbsa
   
2. **Add to System Path**
   
   .. code-block:: bash
      
      # Add the folder containing s_mmpbsa to the system environment variable PATH so that s_mmpbsa can be run from any location.
      export PATH=$PATH:/path/to/s_mmpbsa/
   
3. **Install Necessary Dependencies**
   
   .. code-block:: bash
      
      # Install matplotlib and other necessary Python packages
      
      #### Ubuntu/Debian systems
      sudo apt -y install python3-matplotlib build-essential python-pip
      #### CentOS/Rocky systems
      sudo dnf -y install python3-matplotlib python-pip

Verifying Installation
--------

After installation is complete, you can verify whether s_mmpbsa is installed correctly by the following methods:

.. code-block:: bash
   
   # Run in command line
   s_mmpbsa --version
   
   # Or run s_mmpbsa directly
   s_mmpbsa

If the installation is successful, you will see the welcome message and version number of s_mmpbsa.

Configuring s_mmpbsa
-----------

s_mmpbsa's configuration file is `settings.ini`, which contains various setting parameters of the program. You can modify these parameters as needed to optimize program performance or adjust calculation settings.

### Configuration File Location

- The configuration file is usually located in the directory where the s_mmpbsa executable file is located
- The program will check the location of `settings.ini` when it starts, with priority: current directory > program directory.
- If `settings.ini` is not found, the program will use default settings.

### Main Configuration Parameters

The configuration file contains the following main parameters:

- **gmx_path**: Path to the Gromacs program. It is empty by default, which disables it: leave it empty unless you have a run input (tpr) file written before Gromacs 2021, which the built-in reader does not support. Such a file can also be converted once with `gmx convert-tpr -s old.tpr -o new.tpr`, after which the setting is not needed at all.
- **nkernels**: Number of cores used for parallel computing
- **debug_mode**: Whether to enable debug mode (y/n). When enabled, intermediate files will not be deleted.
- **r_cutoff**: Cutoff distance for non-bonded interactions. 0 means no cutoff.
- **elec_screen**: Electrostatic shielding method setting. 0 means not using electrostatic shielding. 1 means using Debye-Hückel shielding.

Frequently Asked Questions
----------

### Gromacs Not Found

s_mmpbsa only needs a Gromacs program for run input (tpr) files written before Gromacs 2021, which the built-in reader cannot decode. The `gmx_path` setting is empty (turned off) by default, so either set it correctly in settings.ini:

1. `gmx_path` is set to the Gromacs program, e.g. `gmx_path = "gmx"` when `gmx` is on the PATH, or an absolute path such as `gmx_path = "/opt/gromacs/bin/gmx"`
2. The program can be started (check the path, if it is an absolute one)

or re-convert the run input file once with `gmx convert-tpr -s old.tpr -o new.tpr` (using any machine that has Gromacs) and use the converted file, which removes the need for Gromacs entirely.

### APBS Related Errors

The built-in APBS solver is part of the s_mmpbsa binary; if a PB/SA calculation
fails, check the reported solver error (the `.apbs`/`.pqr` intermediate files are
kept when `debug_mode` is enabled) and the PB grid settings (`cfac`, `fadd`, `df`).

### Python/matplotlib Related Errors

If Python or matplotlib related errors occur when using analysis functions, please ensure that:

1. The correct version of Python is installed (Python 3.6 or higher recommended)
2. The matplotlib package is installed

### Performance Issues

If the calculation speed is slow, you can try:

1. Increase the value of the nkernels parameter in settings.ini to utilize more CPU cores
2. For large systems, consider increasing the calculation time interval (i.e., reducing the number of analyzed frames)

Getting Help
--------

If you encounter any problems during the installation process, you can:

- Check the issues page in the GitHub repository: https://github.com/supernova4869/s_mmpbsa/issues
- Contact the developer: zhangjiaxing7137@tju.edu.cn
- Join the QQ group: 864191465
