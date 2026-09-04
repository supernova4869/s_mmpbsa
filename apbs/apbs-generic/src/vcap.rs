// APBS vcap.rs - Capped exponential/hyperbolic functions
// Port of src/generic/vcap.h

use std::f64;

/// Maximum argument for exp/sinh/cosh (prevents overflow)
const EXPMAX: f64 = 85.00;
/// Minimum argument for exp/sinh/cosh (prevents underflow)
const EXPMIN: f64 = -85.00;

/// Capped exponential function.
/// Returns exp(x) clamped to prevent overflow.
/// Sets `ichop` to 1 if the function was capped, 0 otherwise.
#[inline]
pub fn vcap_exp(x: f64, ichop: &mut bool) -> f64 {
    if x > EXPMAX {
        *ichop = true;
        EXPMAX.exp()
    } else if x < EXPMIN {
        *ichop = true;
        EXPMIN.exp()
    } else {
        *ichop = false;
        x.exp()
    }
}

/// Capped hyperbolic sine function.
/// Returns sinh(x) clamped to prevent overflow.
/// Sets `ichop` to 1 if the function was capped, 0 otherwise.
#[inline]
pub fn vcap_sinh(x: f64, ichop: &mut bool) -> f64 {
    if x > EXPMAX {
        *ichop = true;
        EXPMAX.sinh()
    } else if x < EXPMIN {
        *ichop = true;
        EXPMIN.sinh()
    } else {
        *ichop = false;
        x.sinh()
    }
}

/// Capped hyperbolic cosine function.
/// Returns cosh(x) clamped to prevent overflow.
/// Sets `ichop` to 1 if the function was capped, 0 otherwise.
#[inline]
pub fn vcap_cosh(x: f64, ichop: &mut bool) -> f64 {
    if x > EXPMAX || x < EXPMIN {
        *ichop = true;
        EXPMAX.cosh()
    } else {
        *ichop = false;
        x.cosh()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exp_normal() {
        let mut chop = false;
        let result = vcap_exp(1.0, &mut chop);
        assert!(!chop);
        assert!((result - 1.0_f64.exp()).abs() < 1e-12);
    }

    #[test]
    fn test_exp_overflow() {
        let mut chop = false;
        let result = vcap_exp(100.0, &mut chop);
        assert!(chop);
        assert!((result - EXPMAX.exp()).abs() < 1e-12);
    }

    #[test]
    fn test_exp_underflow() {
        let mut chop = false;
        let result = vcap_exp(-100.0, &mut chop);
        assert!(chop);
        assert!((result - EXPMIN.exp()).abs() < 1e-12);
    }

    #[test]
    fn test_sinh_normal() {
        let mut chop = false;
        let result = vcap_sinh(1.0, &mut chop);
        assert!(!chop);
        assert!((result - 1.0_f64.sinh()).abs() < 1e-12);
    }

    #[test]
    fn test_cosh_normal() {
        let mut chop = false;
        let result = vcap_cosh(1.0, &mut chop);
        assert!(!chop);
        assert!((result - 1.0_f64.cosh()).abs() < 1e-12);
    }
}
