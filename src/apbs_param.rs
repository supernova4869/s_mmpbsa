use std::{fmt, fs};
use std::fmt::Formatter;
use std::fs::File;
use std::io::Write;
use std::marker::Copy;
use std::path::Path;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct PBESet {
    pub temp: f64,
    pub pdie: f64,
    pub sdie: f64,
    pub pb_solver: String,
    pub bcfl: String,
    pub srfm: String,
    pub chgm: String,
    pub swin: f64,
    pub srad: f64,
    pub sdens: i32,
    pub ions: Vec<Ion>,
    pub calc_force: bool,
    pub calc_energy: String,
}

impl PBESet {
    pub fn new(temp: f64) -> PBESet {
        return PBESet {
            temp,
            pdie: 2.0,
            sdie: 78.4,
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

    pub fn from(pbe_set: &PBESet) -> PBESet {
        let mut ions: Vec<Ion> = vec![];
        for ion in &pbe_set.ions {
            ions.push(ion.clone());
        }
        let new_pbe_set = PBESet {
            temp: pbe_set.temp,
            pdie: pbe_set.pdie,
            sdie: pbe_set.sdie,
            pb_solver: pbe_set.pb_solver.clone(),
            bcfl: pbe_set.bcfl.clone(),
            srfm: pbe_set.srfm.clone(),
            chgm: pbe_set.chgm.clone(),
            swin: pbe_set.swin.clone(),
            srad: pbe_set.srad.clone(),
            sdens: pbe_set.sdens.clone(),
            ions,
            calc_force: pbe_set.calc_force.clone(),
            calc_energy: pbe_set.calc_energy.clone(),
        };
        return new_pbe_set;
    }

    pub fn load_params<T: AsRef<Path>>(file: T) -> PBESet {
        let pbe_set = fs::read_to_string(&file).expect("Read PB parameters error.");
        let pbe_set: PBESet = serde_json::from_str(pbe_set.as_str()).expect("Read PB parameters error.");
        pbe_set
    }

    pub fn save_params<T: AsRef<Path>>(&self, file: T) {
        let mut f = File::create(&file).expect("Save SA parameters error.");
        f.write_all(serde_json::to_string_pretty(self).unwrap().as_bytes())
            .expect("Save PB parameters error.");
    }
}

impl fmt::Display for PBESet {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let mut ions = String::new();
        for ion in &self.ions {
            ions.push_str(format!("  ion {}\n", ion).as_str());
        }
        write!(f,
            "  temp  {:7}  # 温度\
            \n  pdie  {:7}  # 溶质介电常数\
            \n  sdie  {:7}  # 溶剂介电常数, 真空1, 水78.4 (298.15 K)\
            \n  \
            \n  {}           # PB方程求解方法, lpbe(线性), npbe(非线性), smbpe(大小修正)\
            \n  bcfl  {:>7}  # 粗略格点PB方程的边界条件, zero, sdh/mdh(single/multiple Debye-Huckel), focus, map\
            \n  srfm  {:>7}  # 构建介质和离子边界的模型, mol(分子表面), smol(平滑分子表面), spl2/4(三次样条/7阶多项式)\
            \n  chgm  {:>7}  # 电荷映射到格点的方法, spl0/2/4, 三线性插值, 立方/四次B样条离散\
            \n  swin  {:7}  # 立方样条的窗口值, 仅用于 srfm=spl2/4\
            \n  \
            \n  srad  {:7}  # 溶剂探测半径\
            \n  sdens {:7}  # 表面密度, 每A^2的格点数, (srad=0)或(srfm=spl2/4)时不使用\
            \n  \
            \n  # 离子电荷, 浓度, 半径\
            \n{}  \
            \n  calcforce  {}\
            \n  calcenergy {}", self.temp, self.pdie, self.sdie,
                self.pb_solver, self.bcfl, self.srfm, self.chgm, self.swin,
                self.srad, self.sdens, ions, match self.calc_force {
                    true => "yes",
                    false => "no"
                }, self.calc_energy
        )
    }
}

impl Clone for PBESet {
    fn clone(&self) -> PBESet {
        PBESet::from(self)
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Ion {
    pub charge: f64,
    pub conc: f64,
    radius: f64,
}

impl fmt::Display for Ion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "charge {:2} conc {} radius {}", self.charge, self.conc, self.radius)
    }
}

impl Copy for Ion {}

impl Clone for Ion {
    fn clone(&self) -> Self {
        *self
    }
}

#[derive(Serialize, Deserialize, Debug)]
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

impl PBASet {
    pub fn new(temp: f64) -> Self {
        PBASet {
            temp,
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

    pub fn from(pba_set: &PBASet) -> PBASet {
        PBASet {
            temp: pba_set.temp,
            srfm: pba_set.srfm.to_string(),
            swin: pba_set.swin,
            srad: pba_set.srad,
            gamma: pba_set.gamma,
            press: pba_set.press,
            bconc: pba_set.bconc,
            sdens: pba_set.sdens,
            dpos: pba_set.dpos,
            grid: pba_set.grid,
            calc_force: pba_set.calc_force,
            calc_energy: pba_set.calc_energy.to_string(),
        }
    }

    pub fn load_params<T: AsRef<Path>>(file: T) -> PBASet {
        let pba_set = fs::read_to_string(&file).expect("Read SA parameters error.");
        let pba_set: PBASet = serde_json::from_str(pba_set.as_str()).expect("Read SA parameters error.");
        pba_set
    }

    pub fn save_params<T: AsRef<Path>>(&self, file: T) {
        let mut f = File::create(&file).expect("Save SA parameters error.");
        f.write_all(serde_json::to_string_pretty(self).unwrap().as_bytes())
            .expect("Save SA parameters error.");
    }
}

impl fmt::Display for PBASet {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f,
            "  temp  {:7}  # 温度\
            \n  srfm  {:>7}  # 构建溶剂相关表面或体积的模型\
            \n  swin  {:7}  # 立方样条窗口(A), 用于定义样条表面\
            \n  \
            \n  srad  {:7}  # 探测半径(A)\
            \n  gamma {:7}  # 表面张力(kJ/mol-A^2)\
            \n  \
            \n  press {:7}  # 压力(kJ/mol-A^3)\
            \n  bconc {:7}  # 溶剂本体密度(A^3)\
            \n  sdens {:7}\
            \n  dpos  {:7}\
            \n  grid  {:7} {:5} {:5}\
            \n  \
            \n  calcforce  {}\
            \n  calcenergy {}", self.temp, self.srfm, self.swin,
                self.srad, self.gamma, self.press, self.bconc, self.sdens,
                self.dpos, self.grid.0, self.grid.1, self.grid.2,
                match self.calc_force {
                    true => "yes",
                    false => "no"
                }, self.calc_energy
        )
    }
}

impl Clone for PBASet {
    fn clone(&self) -> PBASet {
        PBASet::from(self)
    }
}
