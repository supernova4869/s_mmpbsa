=======
API Documentation
=======

This document provides detailed information about the internal structure and main functions of s_mmpbsa for developers, helping you understand, extend or modify the functionality of s_mmpbsa.

Project Structure
--------

s_mmpbsa's project structure is as follows:

.. code-block:: bash
   
   s_mmpbsa/
   ├── src/
   │   ├── main.rs           # Program entry point
   │   ├── mmpbsa.rs         # Main implementation of MM/PB-SA calculations
   │   ├── analyzation.rs    # Result analysis functionality
   │   ├── parse_tpr.rs      # TPR file parsing
   │   ├── parse_ndx.rs      # NDX file parsing
   │   ├── parse_xtc.rs      # XTC file parsing
   │   ├── pdb.rs            # PDB file processing
   │   ├── pbsa.rs           # PB/SA calculations
   │   ├── mm.rs             # MM calculations
   │   ├── utils.rs          # Utility functions
   │   └── settings.rs       # Settings management
   ├── examples/             # Example files
   ├── docs/                 # Documentation
   ├── Cargo.toml            # Rust dependency management
   └── README.md             # Project description

Main Modules
--------

The main modules of s_mmpbsa and their functions are described in detail below.

main Module
--------

The main module is the program's entry point, responsible for handling command line parameters, initializing the environment and coordinating the work of other modules.

**Main Functions**:
- Parse command line parameters
- Initialize program environment
- Load input files
- Coordinate MM/PB-SA calculations
- Provide interactive command line interface

**Core Functions**:

.. code-block:: rust
   
   // Program entry point
   fn main() {}
   
   // Display welcome message
   fn welcome() {}
   
   // Verify file validity
   fn confirm_file_validity(path: &str) -> bool {}
   
   // Get built-in program path
   fn get_built_in_gmx() -> String {}

mmpbsa Module
----------

The mmpbsa module implements the core functionality of MM/PB-SA calculations, including energy calculation, alanine scanning, etc.

**Main Functions**:
- Execute MM/PB-SA calculations
- Implement alanine scanning
- Manage temporary files
- Coordinate MM and PB/SA calculations

**Core Functions**:

.. code-block:: rust
   
   // Execute MM/PB-SA calculations
   pub fn fun_mmpbsa_calculations(tpr_path: &str, ...) -> Result<SMResult, Box<dyn Error>> {}
   
   // Implement alanine mutation
   pub fn ala_mutate(tpr_path: &str, ...) -> Result<(), Box<dyn Error>> {}
   
   // Set progress bar style
   fn set_style() -> indicatif::ProgressStyle {}
   
   // Calculate MM/PB-SA energy
   fn calculate_mmpbsa(...) -> Result<SMResult, Box<dyn Error>> {}
   
   // Calculate MM energy
   fn calc_mm(...) -> Result<(Array1<f64>, Array1<f64>), Box<dyn Error>> {}
   
   // Calculate PB/SA energy
   fn calc_pbsa(...) -> Result<(Array1<f64>, Array1<f64>), Box<dyn Error>> {}

analyzation Module
---------------

The analyzation module implements result analysis functionality, including processing, visualization and export of results.

**Main Functions**:
- Process MM/PB-SA calculation results
- Provide result visualization
- Export result data
- Support various analysis operations

**Core Data Structures and Functions**:

.. code-block:: rust
   
   // Data structure for storing MM/PB-SA calculation results
   pub struct SMResult {
       pub dh: Array1<f64>,          // Enthalpy change
       pub mm: Array1<f64>,          // MM energy
       pub pb: Array1<f64>,          // PB energy
       pub sa: Array1<f64>,          // SA energy
       pub time: Array1<f64>,        // Time points
       pub residues: Vec<String>,    // Residue names
       pub res_energy: Array2<f64>,  // Residue energies
       // ... other fields
   }
   
   // Main controller for analysis functionality
   pub fn analyze_controller(result: &SMResult, ...) -> Result<(), Box<dyn Error>> {}
   
   // Get time range
   pub fn get_time_range(result: &SMResult) -> (f64, f64) {}
   
   // Get index corresponding to time point
   pub fn get_time_index(result: &SMResult, time: f64) -> usize {}

parse_tpr Module
------------

The parse_tpr module is responsible for parsing Gromacs TPR files and extracting system topology information and atomic parameters.

**Main Functions**:
- Parse TPR file format
- Extract atomic types, charges, masses, etc.
- Build system topology structure
- Provide access interface to topology data

parse_ndx Module
------------

The parse_ndx module is responsible for parsing Gromacs NDX files and extracting system grouping information.

**Main Functions**:
- Parse NDX file format
- Extract group names and atomic indices
- Provide access interface to grouping data

parse_xtc Module
------------

The parse_xtc module is responsible for parsing Gromacs XTC files and extracting system coordinate information.

**Main Functions**:
- Parse XTC file format
- Extract atomic coordinate data
- Support random access of trajectories
- Handle large trajectory files

pdb Module
--------

The pdb module is responsible for handling PDB files, including reading, modifying and writing PDB files.

**Main Functions**:
- Read PDB files
- Modify atomic coordinates and properties in PDB files
- Write PDB files
- Support encoding energy information into PDB files

pbsa Module
--------

The pbsa module implements the functionality of PB and SA energy calculations, including calling external programs (such as APBS) for calculations.

**Main Functions**:
- Prepare input files for PB calculations
- Call APBS for PB calculations
- Calculate SA energy
- Process PB/SA calculation results

mm Module
--------

The mm module implements the functionality of MM energy calculations, including van der Waals and electrostatic interaction calculations.

**Main Functions**:
- Calculate van der Waals interaction energy
- Calculate electrostatic interaction energy
- Implement distance cutoff optimization
- Support parallel computing

utils Module
--------

The utils module provides various general utility functions for use by other modules.

**Main Functions**:
- File operations
- String processing
- Mathematical calculations
- System calls

settings Module
------------

The settings module is responsible for managing program settings, including reading, modifying and saving settings.

**Main Functions**:
- Read settings.ini file
- Provide access interface to settings
- Save setting changes
- Manage program path configuration

Key Data Structures
------------

SMResult Structure
~~~~~~~~~~~~~

The SMResult structure is the core data structure of s_mmpbsa, used to store the results of MM/PB-SA calculations.

**Main Fields**:
- **dh**: Enthalpy change array
- **mm**: MM energy array
- **pb**: PB energy array
- **sa**: SA energy array
- **time**: Time point array
- **residues**: Residue name list
- **res_energy**: Residue energy matrix
- **atom_energy**: Atomic energy matrix

**Main Methods**:
- **new()**: Create SMResult instance
- **to_bin()**: Serialize results to binary file
- **from_bin()**: Deserialize results from binary file

Using s_mmpbsa as a Library
-----------------

s_mmpbsa can also be used as a Rust library for other Rust programs to call its functions.

**Example Code**:

.. code-block:: rust
   
   use s_mmpbsa::mmpbsa::fun_mmpbsa_calculations;
   use s_mmpbsa::analyzation::SMResult;
   
   fn main() -> Result<(), Box<dyn std::error::Error>> {
       // Set calculation parameters
       let tpr_path = "path/to/md.tpr";
       let xtc_path = "path/to/md_xtc.xtc";
       let ndx_path = "path/to/index.ndx";
       let rec_group = 0;  // Receptor group index
       let lig_group = 1;  // Ligand group index
       let time_interval = 1.0;  // Time interval (ns)
       let temp = 298.15;  // Temperature (K)
       let conc = 0.15;    // Salt concentration (mol/L)
       
       // Execute MM/PB-SA calculation
       let result = fun_mmpbsa_calculations(
           tpr_path,
           xtc_path,
           ndx_path,
           rec_group,
           lig_group,
           time_interval,
           temp,
           conc,
       )?;
       
       // Process calculation results
       println!("Average binding energy: {:.2} kJ/mol", result.dh.mean().unwrap());
       
       // Save results to file
       result.to_bin("result.sm")?;
       
       Ok(())
   }

Extending s_mmpbsa
-----------

If you want to extend the functionality of s_mmpbsa, you can consider the following aspects:

1. **Add new energy calculation methods**: You can add new energy calculation methods in the mm module and pbsa module.

2. **Support new input file formats**: You can add support for new file formats in the parse_tpr, parse_ndx and parse_xtc modules.

3. **Enhance analysis functionality**: You can add new analysis methods and visualization functions in the analyzation module.

4. **Optimize performance**: You can optimize the calculation core to improve calculation speed and memory usage efficiency.

5. **Add new solvation models**: You can add support for other solvation models, such as GB model, 3D-RISM, etc.

Contribution Guidelines
--------

If you want to contribute to the s_mmpbsa project, please follow these steps:

1. Fork the GitHub repository
2. Create your feature branch
3. Commit your changes
4. Push to your branch
5. Create a new Pull Request

Before submitting code, please ensure that your code complies with the project's coding standards and passes all tests.

More Information
--------

- :doc:`usage`：Usage Guide
- :doc:`installation`：Installation Instructions
- :doc:`faq`：Frequently Asked Questions