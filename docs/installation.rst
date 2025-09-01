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

s_mmpbsa's core functionality requires the following software:

- **Gromacs**: Used to process molecular dynamics trajectory files. s_mmpbsa has a built-in Gromacs, but also supports using external Gromacs programs. Multiple versions of Gromacs are supported, but it is recommended to use newer versions for the best compatibility.

Optional dependencies:

- **Matplotlib**: Used to draw result charts. If you need to use s_mmpbsa's analysis and plotting functions, it needs to be installed.
- **APBS**: Used to calculate Poisson-Boltzmann surface area. s_mmpbsa has a built-in APBS kernel, but also supports using external APBS programs.
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

- **gmx_path**: Path to the Gromacs program. If it is "built-in", the program will use the gmx program in /programs/gmx/.
- **apbs_path**: Path to the APBS program. If it is "built-in", the program will use the apbs program in /programs/apbs/.
- **nkernels**: Number of cores used for parallel computing
- **debug_mode**: Whether to enable debug mode (y/n). When enabled, intermediate files will not be deleted.
- **r_cutoff**: Cutoff distance for non-bonded interactions. 0 means no cutoff.
- **elec_screen**: Electrostatic shielding method setting. 0 means not using electrostatic shielding. 1 means using Debye-Hückel shielding.

Frequently Asked Questions
----------

### Gromacs Not Found

If s_mmpbsa cannot find the Gromacs program, please ensure that:

1. Gromacs is installed correctly
2. The directory containing the Gromacs executable file has been added to the system environment variable PATH
3. The gmx_path parameter is set correctly in settings.ini

### APBS Related Errors

If there are issues with the built-in APBS kernel, you can try:

1. Ensure that the built-in APBS program has executable permissions
2. Install an external APBS program
3. Set the apbs_path parameter in settings.ini to point to the external APBS program

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