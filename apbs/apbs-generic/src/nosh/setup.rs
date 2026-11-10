// APBS nosh/setup - Post-parse setup (auto-focusing, etc.)

use std::sync::atomic::{AtomicBool, Ordering};

use crate::error::ApbsResult;
use crate::error::ApbsError;
use crate::vhal::{VREDFRAC, MAXFOCUS, Vbcfl};
use crate::mgparm::{MGparmCalcType, MGparmCentMeth};
use super::{NOsh, NOshCalc, NOshElec};

/// Variable-dime focusing: give every non-finest focusing level its own
/// dimensions (derived from its domain length and grid spacing) instead of
/// carrying the fine level's dime.  With the classic constant-dime scheme
/// the coarse solve costs the same as the fine solve; with spacing doubled
/// per focusing step the coarse chain costs ~1/8 of the fine level.  The
/// fine level (domain, spacing and dims) is untouched.
static VARIABLE_DIME_FOCUSING: AtomicBool = AtomicBool::new(false);

/// Enables/disables the variable-dime focusing hierarchy (see
/// [`VARIABLE_DIME_FOCUSING`]).
pub fn set_focusing_variable_dime(on: bool) {
    VARIABLE_DIME_FOCUSING.store(on, Ordering::Relaxed);
}

fn focusing_variable_dime() -> bool {
    VARIABLE_DIME_FOCUSING.load(Ordering::Relaxed)
}

/// Halves a level's dimension count while keeping the multigrid chain
/// healthy: a dime of n spans n-1 points, and the internal hierarchy
/// needs (n-1) to stay divisible by 4 (two coarsening steps at least).
/// Halving the span exactly, `nc - 1 = (nf - 1) / 2`, preserves every
/// remaining factor of two (fine dims are odd multiples of 32 + 1, so the
/// halved span is a multiple of 16); when the halved span lands odd, drop
/// another factor of two rather than ship an unhealthy dime.
fn halved_dime(nf: i32) -> i32 {
    let mut span = (nf - 1) / 2;
    while span % 2 != 0 && span > 32 {
        span /= 2;
    }
    span + 1
}

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

    // Variable-dime focusing: the coarse chain runs at (about) double the
    // fine spacing, so its levels carry ~1/8 of the fine level's points.
    // The coarse dime is derived from the coarse domain length and that
    // spacing; the fine level keeps the user dime exactly.
    let var_dime = focusing_variable_dime();
    let coarse_dime = [
        halved_dime(dime[0]),
        halved_dime(dime[1]),
        halved_dime(dime[2]),
    ];
    // The exact coarse spacing that follows from the odd dime count.
    let coarse_grid = [
        cglen[0] / (coarse_dime[0] - 1) as f64,
        cglen[1] / (coarse_dime[1] - 1) as f64,
        cglen[2] / (coarse_dime[2] - 1) as f64,
    ];

    // Compute reduction ratio per dimension: redrat[j] = (fgrid[j]/cgrid[j])^(1/(nfocus-1))
    // Port of nosh.c line 1905.  With variable dims the chain runs between
    // the coarse and fine spacings derived above.
    let mut redrat = [0.0f64; 3];
    for d in 0..3 {
        let base_grid = if var_dime { coarse_grid[d] } else { cgrid[d] };
        if base_grid > 0.0 && fgrid[d] > 0.0 {
            let ratio = fgrid[d] / base_grid;
            redrat[d] = ratio.powf(1.0 / (nfocus - 1) as f64);
        } else {
            redrat[d] = 1.0;
        }
    }

    let mut calcs = Vec::new();
    let mut prev_grid = if var_dime { coarse_grid } else { cgrid };
    let mut prev_glen = cglen;

    for ilevel in 0..nfocus {
        let mut new_mgparm = mgparm.clone();
        let mut new_pbeparm = pbeparm.clone();

        if ilevel == 0 {
            // Level 0 (coarsest): use coarse grid parameters
            for d in 0..3 {
                new_mgparm.glen[d] = cglen[d];
                new_mgparm.grid[d] = if var_dime { coarse_grid[d] } else { cgrid[d] };
                new_mgparm.center[d] = ccenter[d];
            }
            if var_dime {
                new_mgparm.dime = coarse_dime;
            }
            new_mgparm.cmeth = mgparm.ccmeth;
            new_pbeparm.bcfl = user_bcfl;
        } else if ilevel == nfocus - 1 {
            // Finest level: the user's dime, spacing and centering, exactly
            // as the input requested.
            for d in 0..3 {
                new_mgparm.grid[d] = fgrid[d];
                new_mgparm.glen[d] = fglen[d];
                new_mgparm.center[d] = fcenter[d];
            }
            new_mgparm.cmeth = mgparm.fcmeth;
            new_pbeparm.bcfl = Vbcfl::Focus;
        } else {
            // Intermediate levels: apply reduction ratio
            for d in 0..3 {
                new_mgparm.grid[d] = prev_grid[d] * redrat[d];
                new_mgparm.glen[d] = prev_glen[d] * redrat[d];
                new_mgparm.center[d] = fcenter[d];
            }
            new_mgparm.cmeth = MGparmCentMeth::Focus;
            new_pbeparm.bcfl = Vbcfl::Focus;
        }

        new_mgparm.r#type = MGparmCalcType::Manual;

        // Variable-dime intermediate levels carry the coarse chain's dims,
        // derived from their own (already reduced) domain and spacing.
        if var_dime && ilevel > 0 && ilevel < nfocus - 1 {
            // Intermediate levels share the coarse chain's dims (their
            // domain shrinks with the same ratio as their spacing).
            new_mgparm.dime = coarse_dime;
            let mut lv_grid = [0.0f64; 3];
            for d in 0..3 {
                lv_grid[d] = new_mgparm.glen[d] / (coarse_dime[d] - 1) as f64;
            }
            new_mgparm.grid = lv_grid;
        }

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mgparm::MGparm;
    use crate::nosh::NOshInputFormat;
    use crate::pbeparm::PBEparm;

    /// A typical mg-auto ELEC: dime 129, cglen ~1.38x fglen, like the
    /// s_mmpbsa default cfac/fadd geometry.
    fn sample_elec() -> NOshElec {
        let mut mgparm = MGparm::new(MGparmCalcType::Auto);
        mgparm.dime = [129, 97, 129];
        mgparm.glen = [63.72, 48.30, 64.32];
        mgparm.cglen = [88.08, 69.78, 82.17];
        mgparm.fglen = [63.72, 48.30, 64.32];
        mgparm.grid = [0.0; 3];
        NOshElec {
            name: "sys_0ns_rec_SOL".to_string(),
            molecules: vec!["1".to_string()],
            input_format: NOshInputFormat::Pqr,
            pbeparm: PBEparm::new(),
            mgparm,
        }
    }

    fn level_dimes(calcs: &[NOshCalc]) -> Vec<[i32; 3]> {
        calcs.iter().map(|calc| match calc {
            NOshCalc::Elec(e) => e.mgparm.dime,
            _ => panic!("expected ELEC calc"),
        }).collect()
    }

    /// Both focusing modes, in one test: the coarse-dime switch is a
    /// process-global flag, so the states must be exercised serially.
    #[test]
    fn focusing_dime_modes() {
        // Default (constant-dime): every level carries the input dime.
        set_focusing_variable_dime(false);
        let dimes = level_dimes(&setup_mgauto_focusing(&sample_elec()).unwrap());
        assert_eq!(dimes.len(), 2, "this geometry yields coarse+fine");
        assert!(dimes.iter().all(|&d| d == [129, 97, 129]),
            "constant dime across levels: {dimes:?}");

        // Variable-dime: the coarse level runs at ~2x the fine spacing
        // (dime ~halved per axis), stays odd, and the fine level keeps the
        // user dims exactly.
        set_focusing_variable_dime(true);
        let calcs = setup_mgauto_focusing(&sample_elec()).unwrap();
        set_focusing_variable_dime(false);
        let dimes = level_dimes(&calcs);
        assert_eq!(dimes.len(), 2);
        let coarse = dimes[0];
        let fine = dimes[1];
        assert_eq!(fine, [129, 97, 129], "fine level unchanged");
        for d in 0..3 {
            assert!(coarse[d] < fine[d], "coarse dim {d} reduced: {coarse:?} vs {fine:?}");
            assert!(coarse[d] % 2 == 1, "coarse dim {d} odd: {coarse:?}");
            // The coarse span keeps enough powers of two for the internal
            // multigrid chain (two coarsening steps at least).
            assert!((coarse[d] - 1) % 4 == 0, "coarse span {d} % 4: {coarse:?}");
        }
        // The coarse spacing is the doubled fine spacing (within rounding).
        let e0 = match &calcs[0] { NOshCalc::Elec(e) => e, _ => unreachable!() };
        for d in 0..3 {
            let hc = e0.mgparm.grid[d];
            let hf = e0.mgparm.fglen[d] / (fine[d] - 1) as f64;
            // Spacing ratio = 2x the domain ratio cglen/fglen (~1.4 for
            // the default cfac/fadd), so expect roughly 2.5-3 here.
            let ratio = hc / hf;
            assert!(ratio > 1.8 && ratio < 4.0,
                "coarse/fine spacing ratio dim {d}: {ratio}");
            // And the derived grid exactly covers cglen with coarse dime.
            let span = (coarse[d] - 1) as f64 * hc;
            assert!((span - e0.mgparm.cglen[d]).abs() < 1e-9,
                "coarse span dim {d}: {span} vs cglen {}",
                e0.mgparm.cglen[d]);
        }
        // The generated mgparms must pass the Manual check (odd dims etc.).
        for calc in &calcs {
            match calc {
                NOshCalc::Elec(e) => e.mgparm.clone()
                    .check()
                    .unwrap_or_else(|err| panic!("level check failed: {err}")),
                _ => unreachable!(),
            }
        }
    }
}
