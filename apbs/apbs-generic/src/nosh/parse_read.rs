// APBS nosh/parse_read - READ statement parser

use crate::error::ApbsResult;
use super::{NOsh, NOshInputFormat};

/// Parse a READ statement block
pub fn parse_read(lines: &[&str], start: usize, nosh: &mut NOsh) -> ApbsResult<usize> {
    let mut i = start;

    // Find matching END
    while i < lines.len() {
        let trimmed = lines[i].trim().to_uppercase();
        if trimmed.starts_with("END") {
            return Ok(i + 1);
        }

        // Parse: mol pqr filename  OR  mesh mcsf filename
        let parts: Vec<&str> = lines[i].trim().split_whitespace().collect();
        if parts.len() >= 3 {
            let record_type = parts[0].to_ascii_lowercase();
            let format_str = parts[1].to_uppercase();
            let filename = parts[2].to_string();

            if record_type == "mol" {
                let _format = match format_str.as_str() {
                    "PQR" => NOshInputFormat::Pqr,
                    "PDB" => NOshInputFormat::Pdb,
                    "XML" => NOshInputFormat::Xml,
                    "FLAT" => NOshInputFormat::Flat,
                    _ => NOshInputFormat::Pqr,
                };

                nosh.mol_paths.push((record_type, filename));
            } else if record_type == "mesh" {
                nosh.mesh_paths.push((format_str.to_ascii_lowercase(), filename));
            }
        }

        i += 1;
    }

    Ok(i)
}
