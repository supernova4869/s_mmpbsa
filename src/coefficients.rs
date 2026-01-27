use crate::parameters::PBESet;

pub struct Coefficients {
    pub f: f64,
    pub lambda_d: f64, // nm
    pub pdie: f64,
}

impl Coefficients {
    pub fn new(pbe_set: &PBESet) -> Coefficients {
        let f = 138.935457520287;  // electric conversion factor, unit: kJ mol^−1 nm e^−2
        let eps0 = 8.854187812800001e-12;  // ε0, unit: F/m
        let kb = 1.380649e-23;  // Boltzmann constant kB, unit: J/K
        let na = 6.02214076e+23;  // NA, unit: mol^-1
        let e_charge = 1.602176634e-19;  // elementary charge e, unit: C
        let ion_strength: f64 = pbe_set.ions.iter()
            .map(|ion| ion.conc * 1e3 * na * ion.charge * ion.charge * e_charge * e_charge).sum();
        let lambda_d = (eps0 * pbe_set.sdie * kb * pbe_set.temp / ion_strength).sqrt() * 1e9;  // nm
        return Coefficients { f, lambda_d, pdie: pbe_set.pdie };
    }
}

// J. Chem. Inf. Model. 2021, 61, 2454
pub fn screening_method(r: f64, coeff: &Coefficients, sm: usize) -> f64 {
    match sm {
        0 => 1.0,
        1 => (-r / coeff.lambda_d).exp(),
        2 => if r > coeff.lambda_d {
            ((coeff.lambda_d - r) / coeff.lambda_d).exp()
        } else {
            1.0
        },
        _ => 1.0
    }
}