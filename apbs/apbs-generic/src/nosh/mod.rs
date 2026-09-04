// APBS nosh - Input file parser
// Port of src/generic/nosh.h / nosh.c

pub mod lexer;
pub mod parse_read;
pub mod parse_elec;
pub mod parse_apolar;
pub mod parse_print;
pub mod setup;

use crate::error::{ApbsError, ApbsResult};
use crate::pbeparm::PBEparm;
use crate::mgparm::MGparm;
use crate::apolparm::APOLparm;

/// Input file format
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NOshInputFormat {
    Pqr,
    Pdb,
    Xml,
    Flat,
}

/// Calculation type
#[derive(Debug, Clone)]
pub enum NOshCalc {
    Elec(NOshElec),
    Apolar(NOshApolar),
    Print(NOshPrint),
}

/// ELECTROSTATICS calculation
#[derive(Debug, Clone)]
pub struct NOshElec {
    pub name: String,
    pub molecules: Vec<String>,
    pub input_format: NOshInputFormat,
    pub pbeparm: PBEparm,
    pub mgparm: MGparm,
}

/// APOLAR calculation
#[derive(Debug, Clone)]
pub struct NOshApolar {
    pub name: String,
    pub molecules: Vec<String>,
    pub apolparm: APOLparm,
}

/// PRINT statement
#[derive(Debug, Clone)]
pub enum NOshPrint {
    ElecEnergy { left: String, right: String },
    ApolEnergy { name: String },
    Raw { keyword: String, calc_name: String, print_type: String },
}

/// Input file parser
pub struct NOsh {
    pub calcs: Vec<NOshCalc>,
    pub mol_paths: Vec<(String, String)>, // (name, path)
    pub mesh_paths: Vec<(String, String)>, // (format, path)
    /// Directory containing the input file
    pub input_dir: String,
}

impl NOsh {
    pub fn new() -> Self {
        Self {
            calcs: Vec::new(),
            mol_paths: Vec::new(),
            mesh_paths: Vec::new(),
            input_dir: String::new(),
        }
    }

    /// Read and parse input file
    pub fn read(&mut self, filename: &str) -> ApbsResult<()> {
        // Store input file directory for resolving relative paths
        if let Some(dir) = std::path::Path::new(filename).parent() {
            self.input_dir = dir.to_string_lossy().to_string();
        }

        let content = std::fs::read_to_string(filename)
            .map_err(|e| ApbsError::Io(format!("{}: {}", filename, e)))?;

        let cleaned_lines: Vec<String> = content
            .lines()
            .filter_map(|line| {
                let without_comment = line.split('#').next().unwrap_or("").trim();
                if without_comment.is_empty() {
                    None
                } else {
                    Some(without_comment.to_string())
                }
            })
            .collect();
        let lines: Vec<&str> = cleaned_lines.iter().map(String::as_str).collect();

        let mut i = 0;
        while i < lines.len() {
            let trimmed = lines[i].trim();
            let upper = trimmed.to_uppercase();

            if upper.starts_with("READ") {
                i = parse_read::parse_read(&lines, i, self)?;
            } else if upper.starts_with("ELEC") {
                i = parse_elec::parse_elec(&lines, i, self)?;
            } else if upper.starts_with("APOLAR") {
                i = parse_apolar::parse_apolar(&lines, i, self)?;
            } else if upper.starts_with("PRINT") {
                i = parse_print::parse_print(&lines, i, self)?;
            } else if upper == "QUIT" {
                i += 1;
            } else {
                let keyword = trimmed
                    .split_whitespace()
                    .next()
                    .unwrap_or(trimmed);
                return Err(ApbsError::UnsupportedFormat(format!(
                    "Unsupported top-level APBS token: {}",
                    keyword
                )));
            }
        }

        // Post-parse setup
        setup::post_parse(self)?;

        Ok(())
    }

    /// Get molecule path by name
    pub fn get_mol_path(&self, name: &str) -> ApbsResult<String> {
        // First try exact name match
        for (n, path) in &self.mol_paths {
            if n == name {
                return Ok(self.resolve_path(path));
            }
        }
        // Then try index-based lookup (1-based)
        if let Ok(idx) = name.parse::<usize>() {
            if idx >= 1 && idx <= self.mol_paths.len() {
                return Ok(self.resolve_path(&self.mol_paths[idx - 1].1));
            }
        }
        Err(ApbsError::FileNotFound(format!("Molecule '{}' not found", name)))
    }

    pub fn get_mesh_path(&self, idx_1based: i32) -> ApbsResult<String> {
        let (_, path) = self.get_mesh_info(idx_1based)?;
        Ok(self.resolve_path(&path))
    }

    pub fn get_mesh_info(&self, idx_1based: i32) -> ApbsResult<(String, String)> {
        let idx = usize::try_from(idx_1based.saturating_sub(1)).map_err(|_| {
            ApbsError::FileNotFound(format!("Mesh '{}' not found", idx_1based))
        })?;
        if let Some((format, path)) = self.mesh_paths.get(idx) {
            Ok((format.clone(), path.clone()))
        } else {
            Err(ApbsError::FileNotFound(format!("Mesh '{}' not found", idx_1based)))
        }
    }

    /// Resolve a path relative to the input file directory
    fn resolve_path(&self, path: &str) -> String {
        let p = std::path::Path::new(path);
        if p.is_absolute() {
            path.to_string()
        } else if !self.input_dir.is_empty() {
            std::path::Path::new(&self.input_dir)
                .join(path)
                .to_string_lossy()
                .to_string()
        } else {
            path.to_string()
        }
    }

    /// Get temperature for an ELEC calculation by name
    pub fn get_elec_temp(&self, name: &str) -> Option<f64> {
        for calc in &self.calcs {
            if let NOshCalc::Elec(elec) = calc {
                if elec.name == name {
                    return Some(elec.pbeparm.temp);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::NOsh;

    #[test]
    fn unknown_top_level_token_errors() {
        let input = std::env::temp_dir().join("apbs-unknown-top-level.in");
        std::fs::write(&input, "frobnicate\nquit\n").expect("write temp input");

        let mut nosh = NOsh::new();
        let err = nosh.read(input.to_str().expect("temp path")).expect_err("expected parse failure");
        let msg = format!("{}", err);
        assert!(msg.contains("Unsupported top-level APBS token"));

        let _ = std::fs::remove_file(input);
    }
}
