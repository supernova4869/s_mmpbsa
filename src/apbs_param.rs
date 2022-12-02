use std::fmt;
use std::fmt::Formatter;

pub struct PBESet {
    pub temp: f64,
    pub pdie: f64,
    pub sdie: f64,
    pb_solver: String,
    bcfl: String,
    srfm: String,
    chgm: String,
    swin: f64,
    srad: f64,
    sdens: i32,
    pub ions: Vec<Ion>,
    calc_force: bool,
    calc_energy: String,
}

impl PBESet {
    pub fn new() -> PBESet {
        return PBESet {
            temp: 298.15,
            pdie: 2.0,
            sdie: 78.54,
            pb_solver: "npbe".to_string(),
            bcfl: "mdh".to_string(),
            srfm: "smol".to_string(),
            chgm: "spl4".to_string(),
            swin: 0.3,
            srad: 1.4,
            sdens: 10,
            ions: vec![
                Ion { charge: 1.0, conc: 0.15, radius: 0.95 },
                Ion { charge: -1.0, conc: 0.15, radius: 1.81 },
            ],
            calc_force: false,
            calc_energy: "comps".to_string(),
        };
    }
}

impl fmt::Display for PBESet {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut ions = String::new();
        for ion in &self.ions {
            ions.push_str(format!("  ion {}\n", ion).as_str());
        }
        write!(f, "  temp  {}\
        \n  pdie  {}\
        \n  sdie  {}\
        \n  \
        \n  {}\
        \n  bcfl  {}\
        \n  srfm  {}\
        \n  chgm  {}\
        \n  swin  {}\
        \n  \
        \n  srad  {}\
        \n  sdens {}\
        \n  \
        \n{}  \
        \n  calcforce  {}\
        \n  calcenergy {}", self.temp, self.pdie, self.sdie,
               self.pb_solver, self.bcfl, self.srfm, self.chgm, self.swin,
               self.srad, self.sdens, ions, match self.calc_force {
                true => "yes",
                false => "no"
            }, self.calc_energy)
    }
}

pub struct PBASet {
    temp: f64,
    srfm: String,
    swin: f64,
    srad: f64,
    pub gamma: f64,
    press: f64,
    bconc: f64,
    sdens: f64,
    dpos: f64,
    grid: (f64, f64, f64),
    calc_force: bool,
    calc_energy: String,
}

pub struct Ion {
    pub charge: f64,
    pub conc: f64,
    radius: f64,
}

impl fmt::Display for Ion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "charge  {} conc {} radius {}", self.charge, self.conc, self.radius)
    }
}

impl PBASet {
    pub fn new() -> Self {
        PBASet {
            temp: 298.15,
            srfm: "sacc".to_string(),
            swin: 0.3,
            srad: 1.4,
            gamma: 1.0,
            press: 0.0,
            bconc: 0.0,
            sdens: 10.0,
            dpos: 0.2,
            grid: (0.1, 0.1, 0.1),
            calc_force: false,
            calc_energy: "total".to_string(),
        }
    }
}

impl fmt::Display for PBASet {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "  temp  {}\
        \n  srfm  {}\
        \n  swin  {}\
        \n  \
        \n  srad  {}\
        \n  gamma  {}\
        \n  \
        \n  press  {}\
        \n  bconc  {}\
        \n  sdens  {}\
        \n  dpos  {}\
        \n  grid  {} {} {}\
        \n  \
        \n  calcforce  {}\
        \n  calcenergy {}", self.temp, self.srfm, self.swin,
               self.srad, self.gamma, self.press, self.bconc, self.sdens,
               self.dpos, self.grid.0, self.grid.1, self.grid.2,
               match self.calc_force {
                   true => "yes",
                   false => "no"
               }, self.calc_energy)
    }
}
