// APBS nosh/setup - Post-parse setup (auto-focusing, etc.)

use crate::error::ApbsResult;
use crate::error::ApbsError;
use crate::vhal::{VREDFRAC, MAXFOCUS, Vbcfl};
use crate::mgparm::{MGparmCalcType, MGparmCentMeth};
use super::{NOsh, NOshCalc, NOshElec};

/// Post-parse setup routines
pub fn post_parse(nosh: &mut NOsh) -> ApbsResult<()> {
    setup_calc_mgauto(nosh)?;
    Ok(())
}

/// Setup mg-auto focusing calculations
/// Port of NOsh_setupCalcMGAUTO from nosh.c
fn setup_calc_mgauto(nosh: &mut NOsh) -> ApbsResult<()> {
    let mut new_calcs = Vec::new();

    for calc in &nosh.calcs {
        if let NOshCalc::Elec(elec) = calc {
            if elec.mgparm.r#type == MGparmCalcType::Auto {
                let focusing_calcs = setup_mgauto_focusing(elec)?;
                new_calcs.extend(focusing_calcs);
            } else {
                new_calcs.push(calc.clone());
            }
        } else {
            new_calcs.push(calc.clone());
        }
    }

    nosh.calcs = new_calcs;
    Ok(())
}

/// Create focusing hierarchy for mg-auto
/// Port of NOsh_setupCalcMGAUTO focusing logic from nosh.c
fn setup_mgauto_focusing(elec: &NOshElec) -> ApbsResult<Vec<NOshCalc>> {
    let mgparm = &elec.mgparm;
    let pbeparm = &elec.pbeparm;

    let dime = mgparm.dime;
    let fglen = mgparm.fglen;
    let cglen = mgparm.cglen;
    let fcenter = mgparm.fcenter;
    let ccenter = mgparm.ccenter;
    let user_bcfl = pbeparm.bcfl;

    // Compute grid spacings for coarse and fine
    // cgrid[j] = cglen[j] / (dime[j] - 1), fgrid[j] = fglen[j] / (dime[j] - 1)
    let mut cgrid = [0.0f64; 3];
    let mut fgrid = [0.0f64; 3];
    for d in 0..3 {
        if dime[d] > 1 {
            cgrid[d] = cglen[d] / (dime[d] - 1) as f64;
            fgrid[d] = fglen[d] / (dime[d] - 1) as f64;
        }
    }

    // Compute number of focusing levels per dimension
    // Port of nosh.c lines 1893-1900
    let mut tnfocus = [0i32; 3];
    for d in 0..3 {
        if cgrid[d] > 0.0 && fgrid[d] > 0.0 {
            let ratio = fgrid[d] / cgrid[d];
            if ratio < VREDFRAC {
                // Need more levels to avoid reducing by more than VREDFRAC per level
                tnfocus[d] = (ratio.ln() / VREDFRAC.ln()).ceil() as i32 + 1;
            } else {
                // At least 2 levels (coarse + fine)
                tnfocus[d] = 2;
            }
        }
    }

    let nfocus = *tnfocus.iter().max().unwrap_or(&0);
    let nfocus = nfocus.max(1).min(MAXFOCUS as i32);

    if nfocus <= 1 {
        // No focusing needed - single level with manual calc type
        let mut new_elec = elec.clone();
        new_elec.mgparm.r#type = MGparmCalcType::Manual;
        for d in 0..3 {
            new_elec.mgparm.center[d] = fcenter[d];
            if dime[d] > 1 {
                new_elec.mgparm.glen[d] = fglen[d];
                new_elec.mgparm.grid[d] = fglen[d] / (dime[d] - 1) as f64;
            }
        }
        return Ok(vec![NOshCalc::Elec(new_elec)]);
    }

    // Compute reduction ratio per dimension: redrat[j] = (fgrid[j]/cgrid[j])^(1/(nfocus-1))
    // Port of nosh.c line 1905
    let mut redrat = [0.0f64; 3];
    for d in 0..3 {
        if cgrid[d] > 0.0 && fgrid[d] > 0.0 {
            let ratio = fgrid[d] / cgrid[d];
            redrat[d] = ratio.powf(1.0 / (nfocus - 1) as f64);
        } else {
            redrat[d] = 1.0;
        }
    }

    let mut calcs = Vec::new();
    let mut prev_grid = cgrid;
    let mut prev_glen = cglen;

    for ilevel in 0..nfocus {
        let mut new_mgparm = mgparm.clone();
        let mut new_pbeparm = pbeparm.clone();

        if ilevel == 0 {
            // Level 0 (coarsest): use coarse grid parameters
            for d in 0..3 {
                new_mgparm.glen[d] = cglen[d];
                new_mgparm.grid[d] = cgrid[d];
                new_mgparm.center[d] = ccenter[d];
            }
            new_mgparm.cmeth = mgparm.ccmeth;
            new_pbeparm.bcfl = user_bcfl;
        } else {
            // Intermediate and fine levels: apply reduction ratio
            for d in 0..3 {
                new_mgparm.grid[d] = prev_grid[d] * redrat[d];
                new_mgparm.glen[d] = prev_glen[d] * redrat[d];
                new_mgparm.center[d] = fcenter[d];
            }
            new_mgparm.cmeth = MGparmCentMeth::Focus;

            // All levels beyond coarsest get Focus BCs (port of nosh.c line 2078)
            new_pbeparm.bcfl = Vbcfl::Focus;

            if ilevel == nfocus - 1 {
                // Finest level: use fine centering method
                new_mgparm.cmeth = mgparm.fcmeth;
            }
        }

        new_mgparm.r#type = MGparmCalcType::Manual;

        // Mesh repositioning: ensure coarse grid fully contains fine grid
        // Port of nosh.c lines 1981-2074
        if ilevel > 0 {
            for d in 0..3 {
                let half_len = new_mgparm.glen[d] / 2.0;
                let xmin_level = new_mgparm.center[d] - half_len;

                // Get previous level bounds
                let prev_half = prev_glen[d] / 2.0;
                let prev_center = if ilevel == 1 {
                    ccenter[d]
                } else {
                    calcs
                        .last()
                        .and_then(|calc| match calc {
                            NOshCalc::Elec(prev) => Some(prev.mgparm.center[d]),
                            _ => None,
                        })
                        .unwrap_or(fcenter[d])
                };
                let prev_xmin = prev_center - prev_half;
                let prev_xmax = prev_center + prev_half;

                // If lower boundary extends past previous, shift center up
                if xmin_level < prev_xmin {
                    let delta = prev_xmin - xmin_level;
                    if ilevel == nfocus - 1 {
                        return Err(ApbsError::InvalidParameter(format!(
                            "Finest focusing mesh falls below coarser mesh in dimension {} by {}",
                            d, delta
                        )));
                    }
                    new_mgparm.center[d] += delta;
                }
                // If upper boundary extends past previous, shift center down
                let new_xmax = new_mgparm.center[d] + half_len;
                if new_xmax > prev_xmax {
                    let delta = new_xmax - prev_xmax;
                    if ilevel == nfocus - 1 {
                        return Err(ApbsError::InvalidParameter(format!(
                            "Finest focusing mesh falls above coarser mesh in dimension {} by {}",
                            d, delta
                        )));
                    }
                    new_mgparm.center[d] -= delta;
                }
            }
        }

        let level_name = if ilevel == 0 {
            format!("{}_{}", elec.name, ilevel)
        } else if ilevel == nfocus - 1 {
            elec.name.clone()
        } else {
            format!("{}_{}", elec.name, ilevel)
        };

        calcs.push(NOshCalc::Elec(NOshElec {
            name: level_name,
            molecules: elec.molecules.clone(),
            input_format: elec.input_format,
            pbeparm: new_pbeparm,
            mgparm: new_mgparm.clone(),
        }));

        prev_grid = new_mgparm.grid;
        prev_glen = new_mgparm.glen;
    }

    Ok(calcs)
}
