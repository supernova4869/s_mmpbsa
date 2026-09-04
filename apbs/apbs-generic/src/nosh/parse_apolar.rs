// APBS nosh/parse_apolar - APOLAR statement parser

use crate::error::ApbsResult;
use crate::apolparm::APOLparm;
use super::{NOsh, NOshCalc, NOshApolar};

/// Parse an APOLAR statement block
pub fn parse_apolar(lines: &[&str], start: usize, nosh: &mut NOsh) -> ApbsResult<usize> {
    let mut i = start;
    let mut name = String::new();
    let mut molecules = Vec::new();
    let mut apolparm = APOLparm::new();

    // Parse APOLAR line: "APOLAR name VALUE" or "APOLAR VALUE"
    let first_line: Vec<&str> = lines[i].trim().split_whitespace().collect();
    if first_line.len() > 2 && first_line[1].to_ascii_lowercase() == "name" {
        name = first_line[2].to_string();
    } else if first_line.len() > 1 {
        name = first_line[1].to_string();
    }
    i += 1;

    while i < lines.len() {
        let trimmed = lines[i].trim().to_uppercase();

        if trimmed.starts_with("END") {
            nosh.calcs.push(NOshCalc::Apolar(NOshApolar {
                name,
                molecules,
                apolparm,
            }));
            return Ok(i + 1);
        }

        // Parse tokens
        let parts: Vec<&str> = lines[i].trim().split_whitespace().collect();
        if !parts.is_empty() {
            let token = parts[0].to_uppercase();
            let value = if parts.len() > 1 { parts[1..].join(" ") } else { String::new() };

            match token.as_str() {
                "MOL" => {
                    if !value.is_empty() {
                        molecules.push(value.clone());
                    }
                    apolparm.parse_token("mol", &value)?;
                }
                _ => {
                    apolparm.parse_token(&token.to_lowercase(), &value)?;
                }
            }
        }

        i += 1;
    }

    // If we reach here without END, still push the calc
    nosh.calcs.push(NOshCalc::Apolar(NOshApolar {
        name,
        molecules,
        apolparm,
    }));

    Ok(i)
}
