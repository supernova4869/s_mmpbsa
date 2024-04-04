use crate::apbs_param::PBESet;

pub struct Coefficients {
    pub kj_elec: f64,
    pub kap: f64,
    pub pdie: f64,
}

impl Coefficients {
    pub fn new(pbe_set: &PBESet) -> Coefficients {
        let eps0 = 8.854187812800001e-12;
        let kj_elec = 1389.35457520287;
        let kb = 1.380649e-23;
        let na = 6.02214076e+23;
        let qe = 1.602176634e-19;
        let ion_strength: f64 = pbe_set.ions.iter()
            .map(|ion| ion.charge * ion.charge * ion.conc).sum();
        let kap = 1e-10 / f64::sqrt(eps0 * kb * pbe_set.temp * pbe_set.sdie / (ion_strength * qe * qe * na * 1e3));
        let pdie = pbe_set.pdie;
        return Coefficients { kj_elec, kap, pdie };
    }
}
