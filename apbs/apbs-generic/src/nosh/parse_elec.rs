// APBS nosh/parse_elec - ELEC statement parser

use crate::error::{ApbsError, ApbsResult};
use crate::pbeparm::PBEparm;
use crate::mgparm::{MGparm, MGparmCalcType};
use super::{NOsh, NOshCalc, NOshElec, NOshInputFormat};
use std::sync::OnceLock;

fn debug_enabled() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("APBS_RUST_DEBUG").is_some())
}

/// Parse an ELEC statement block
pub fn parse_elec(lines: &[&str], start: usize, nosh: &mut NOsh) -> ApbsResult<usize> {
    let mut i = start;
    let mut name = String::new();
    let mut molecules = Vec::new();
    let input_format = NOshInputFormat::Pqr;
    let mut pbeparm = PBEparm::new();
    let mut mgparm = MGparm::new(MGparmCalcType::Manual);

    // Parse ELEC line: "ELEC name VALUE" or "ELEC VALUE" or just "ELEC"
    let first_line: Vec<&str> = lines[i].trim().split_whitespace().collect();
    if first_line.len() > 2 && first_line[1].to_ascii_lowercase() == "name" {
        // "ELEC name 5.0ns_com" format
        name = first_line[2].to_string();
    } else if first_line.len() > 1 {
        name = first_line[1].to_string();
    }
    i += 1;

    while i < lines.len() {
        let trimmed = lines[i].trim().to_uppercase();

        if trimmed.starts_with("END") {
            mgparm.check()?;
            nosh.calcs.push(NOshCalc::Elec(NOshElec {
                name: name.clone(),
                molecules: molecules.clone(),
                input_format,
                pbeparm: pbeparm.clone(),
                mgparm: mgparm.clone(),
            }));
            return Ok(i + 1);
        }

        // Split line into tokens
        let parts: Vec<&str> = lines[i].trim().split_whitespace().collect();
        if parts.is_empty() {
            i += 1;
            continue;
        }

        let token = parts[0].to_uppercase();
        let value = if parts.len() > 1 {
            parts[1..].join(" ")
        } else {
            String::new()
        };
        if debug_enabled() {
            eprintln!("[DEBUG-PARSE] ELEC token='{}'", token);
        }

        match token.as_str() {
            // Molecule reference: "mol 1", "mol 2", "mol 3"
            "MOL" => {
                if parts.len() > 1 {
                    molecules.push(parts[1].to_string());
                    // Also set molid on pbeparm (1-based in input, store as-is)
                    if let Ok(id) = parts[1].parse::<i32>() {
                        pbeparm.molid = Some(id);
                    }
                }
            }
            // MG-auto calculation type
            "MG-AUTO" => {
                mgparm.r#type = MGparmCalcType::Auto;
            }
            "MG-MANUAL" => {
                mgparm.r#type = MGparmCalcType::Manual;
            }
            // MG parameters: dime, glen, gcent, grid, chgm, cglen, fglen, cgcent, fgcent
            "DIME" | "GLEN" | "GCENT" | "GRID" | "CHGM" | "CGLEN" | "FGLEN" | "CGCENT" | "FGCENT" | "NLEV" => {
                mgparm.parse_token(&token.to_lowercase(), &value)?;
            }
            // PBE type
            "NPBE" => {
                pbeparm.pbetype = crate::vhal::VhalPBEType::NPBE;
            }
            "LPBE" => {
                pbeparm.pbetype = crate::vhal::VhalPBEType::LPBE;
            }
            // Ion line: "ion charge Q conc C radius R"
            "ION" => {
                parse_ion_line(&parts, &mut pbeparm)?;
            }
            "CHARGE" | "CONC" | "RADIUS" => {
                update_last_ion_field(&token, &parts, &mut pbeparm)?;
            }
            // All other tokens → try pbeparm first, then mgparm
            _ => {
                let token_lower = token.to_lowercase();
                // Some tokens belong to pbeparm, some to mgparm
                let pbeparm_tokens = [
                    "pdie", "sdie", "temp", "srad", "sdens", "swin", "srfm",
                    "bcfl", "nion", "calcenergy", "calcforce", "zmem", "lmem",
                    "mdie", "memv", "smsize", "smvolume", "usemap", "write",
                ];
                if pbeparm_tokens.contains(&token_lower.as_str()) {
                    pbeparm.parse_token(&token_lower, &value)?;
                } else {
                    return Err(ApbsError::InvalidParameter(format!(
                        "Unsupported ELEC token: {}",
                        parts[0]
                    )));
                }
            }
        }

        i += 1;
    }

    // If we reach here without END, still push the calc
    mgparm.check()?;
    nosh.calcs.push(NOshCalc::Elec(NOshElec {
        name,
        molecules,
        input_format,
        pbeparm,
        mgparm,
    }));

    Ok(i)
}

/// Parse an ion line: "ion charge Q conc C radius R"
fn parse_ion_line(parts: &[&str], pbeparm: &mut PBEparm) -> ApbsResult<()> {
    // Expected format: ion charge Q conc C radius R
    let mut charge = 0.0;
    let mut conc = 0.0;
    let mut radius = 0.0;

    let mut j = 1; // skip "ion"
    while j < parts.len() {
        match parts[j].to_ascii_lowercase().as_str() {
            "charge" => {
                j += 1;
                if j < parts.len() {
                    charge = parts[j].parse().map_err(|_| crate::error::ApbsError::Parse {
                        line: 0,
                        message: format!("Invalid ion charge: {}", parts[j]),
                    })?;
                }
            }
            "conc" => {
                j += 1;
                if j < parts.len() {
                    conc = parts[j].parse().map_err(|_| crate::error::ApbsError::Parse {
                        line: 0,
                        message: format!("Invalid ion conc: {}", parts[j]),
                    })?;
                }
            }
            "radius" => {
                j += 1;
                if j < parts.len() {
                    radius = parts[j].parse().map_err(|_| crate::error::ApbsError::Parse {
                        line: 0,
                        message: format!("Invalid ion radius: {}", parts[j]),
                    })?;
                }
            }
            _ => {}
        }
        j += 1;
    }

    // Store ion parameters
    if pbeparm.nion < crate::vhal::MAXION {
        let idx = pbeparm.nion;
        pbeparm.ionq[idx] = charge;
        pbeparm.ionc[idx] = conc;
        pbeparm.ionr[idx] = radius;
        pbeparm.nion += 1;
    }

    Ok(())
}

fn update_last_ion_field(token: &str, parts: &[&str], pbeparm: &mut PBEparm) -> ApbsResult<()> {
    if pbeparm.nion == 0 {
        return Err(ApbsError::InvalidParameter(format!(
            "Ion continuation '{}' found before any ion declaration",
            token
        )));
    }
    if parts.len() < 2 {
        return Err(ApbsError::InvalidParameter(format!(
            "Ion continuation '{}' requires a value",
            token
        )));
    }

    let idx = pbeparm.nion - 1;
    match token {
        "CHARGE" => {
            pbeparm.ionq[idx] = parts[1].parse().map_err(|_| ApbsError::Parse {
                line: 0,
                message: format!("Invalid ion charge: {}", parts[1]),
            })?;
        }
        "CONC" => {
            pbeparm.ionc[idx] = parts[1].parse().map_err(|_| ApbsError::Parse {
                line: 0,
                message: format!("Invalid ion conc: {}", parts[1]),
            })?;
        }
        "RADIUS" => {
            pbeparm.ionr[idx] = parts[1].parse().map_err(|_| ApbsError::Parse {
                line: 0,
                message: format!("Invalid ion radius: {}", parts[1]),
            })?;
        }
        _ => {
            return Err(ApbsError::InvalidParameter(format!(
                "Unsupported ion continuation token: {}",
                token
            )));
        }
    }

    Ok(())
}
