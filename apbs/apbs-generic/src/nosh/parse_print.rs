// APBS nosh/parse_print - PRINT statement parser

use crate::error::ApbsResult;
use super::{NOsh, NOshCalc, NOshPrint};

/// Parse PRINT statements.
/// In APBS, each "print" line is a complete statement (not a block).
/// This function processes all consecutive print lines and returns.
pub fn parse_print(lines: &[&str], start: usize, nosh: &mut NOsh) -> ApbsResult<usize> {
    let mut i = start;

    while i < lines.len() {
        let trimmed = lines[i].trim().to_uppercase();

        // Stop at non-PRINT lines (ELEC, APOLAR, READ, etc.)
        if !trimmed.starts_with("PRINT") {
            break;
        }

        let parts: Vec<&str> = lines[i].trim().split_whitespace().collect();
        if parts.len() < 2 {
            nosh.calcs.push(NOshCalc::Print(NOshPrint::Raw {
                keyword: "print".to_string(),
                calc_name: String::new(),
                print_type: lines[i].trim().to_string(),
            }));
            i += 1;
            continue;
        }

        // The second token is the keyword (elecEnergy, apolEnergy, etc.)
        let keyword = parts[1].to_lowercase();

        match keyword.as_str() {
            "elecenergy" => {
                let print_stmt = parse_elec_energy(&parts);
                nosh.calcs.push(NOshCalc::Print(print_stmt));
            }
            "apolenergy" | "apol" => {
                let print_stmt = parse_apole_energy(&parts);
                nosh.calcs.push(NOshCalc::Print(print_stmt));
            }
            _ => {
                nosh.calcs.push(NOshCalc::Print(NOshPrint::Raw {
                    keyword: "print".to_string(),
                    calc_name: String::new(),
                    print_type: lines[i].trim().to_string(),
                }));
            }
        }

        i += 1;
    }

    Ok(i)
}

/// Parse elecEnergy print statement
/// Format: print elecEnergy NAME1 - NAME2 end
/// or: print elecEnergy NAME end (just one energy)
fn parse_elec_energy(parts: &[&str]) -> NOshPrint {
    // Find the "-" separator (skip "print" and "elecEnergy")
    let mut dash_idx = None;
    for (j, p) in parts.iter().enumerate() {
        if *p == "-" && j > 1 {
            dash_idx = Some(j);
            break;
        }
    }

    if let Some(dash_idx) = dash_idx {
        // Two-energy subtraction: print elecEnergy NAME1 - NAME2 end
        let left = parts[dash_idx - 1].to_string();
        // Collect right name (skip "end" tokens)
        let right_parts: Vec<&str> = parts[(dash_idx + 1)..]
            .iter()
            .filter(|p| p.to_uppercase() != "END")
            .copied()
            .collect();
        let right = right_parts.join(" ");
        NOshPrint::ElecEnergy { left, right }
    } else {
        // Single energy: print elecEnergy NAME end
        let name = parts.get(2)
            .filter(|s| s.to_uppercase() != "END")
            .map(|s| s.to_string())
            .unwrap_or_default();
        NOshPrint::ElecEnergy { left: name, right: String::new() }
    }
}

/// Parse apolEnergy print statement
/// Format: print apolEnergy NAME end
fn parse_apole_energy(parts: &[&str]) -> NOshPrint {
    let name = parts.get(2)
        .filter(|s| s.to_uppercase() != "END")
        .map(|s| s.to_string())
        .unwrap_or_default();
    NOshPrint::ApolEnergy { name }
}

#[cfg(test)]
mod tests {
    use super::parse_print;
    use crate::nosh::{NOsh, NOshCalc, NOshPrint};

    #[test]
    fn unknown_print_is_preserved_as_raw() {
        let mut nosh = NOsh::new();
        let lines = vec!["print mysteryEnergy foo end"];
        let next = parse_print(&lines, 0, &mut nosh).expect("parse should succeed");
        assert_eq!(next, 1);
        assert_eq!(nosh.calcs.len(), 1);
        match &nosh.calcs[0] {
            NOshCalc::Print(NOshPrint::Raw { print_type, .. }) => {
                assert_eq!(print_type, "print mysteryEnergy foo end");
            }
            other => panic!("expected raw print, got {:?}", other),
        }
    }
}
