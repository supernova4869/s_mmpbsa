use std::fs::{File, read_to_string};
use std::io::{Write, BufWriter};
use std::path::Path;
use std::fmt::{self, Formatter};

// 定义结构体
#[derive(Debug)]
#[allow(dead_code)]
pub struct GRO {
    pub title: String,
    pub atom_num: usize,
    pub atoms: Vec<GroAtom>,
    pub box_vec: Vec<f64>,
}

#[derive(Debug)]
pub struct GroAtom {
    pub resid: i32,
    pub resname: String,
    pub atname: String,
    pub atid: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

// 公开API
impl GRO {
    pub fn from<P: AsRef<Path>>(filename: P) -> GRO {
        let content = read_to_string(filename).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        
        // 处理title
        let title = lines[0].to_string();
        
        // 处理原子数量
        let atom_num: usize = lines[1].trim().parse().unwrap();
        
        // 处理原子坐标
        let mut atoms = Vec::with_capacity(atom_num);
        for i in 2..(2 + atom_num) {
            atoms.push(GroAtom::from(lines[i]));
        }
        
        // 处理盒子
        let box_line = lines[2 + atom_num];
        
        let box_vec: Vec<f64> = box_line
            .split_whitespace()
            .map(|s| s.parse::<f64>().unwrap())
            .collect();
        
        GRO {
            title,
            atom_num,
            atoms,
            box_vec
        }
    }

    #[allow(dead_code)]
    pub fn to_gro(&self, gro_file: &str) {
        let file = File::create(gro_file).unwrap();
        let mut writer = BufWriter::new(file);
        
        writeln!(writer, "{}, Created by s_mmpbsa (https://github.com/supernova4869/s_mmpbsa)", self.title).unwrap();
        writeln!(writer, "{}", self.atom_num).unwrap();
        
        for atom in &self.atoms {
            writeln!(writer, "{}", atom).unwrap();
        }
        
        for coord in &self.box_vec {
            write!(writer, "{:10.5}", coord).unwrap();
        }
        writeln!(writer).unwrap();
    }
}

impl GroAtom {
    fn from(atomline: &str) -> GroAtom {
        let resid = atomline[0..5].trim().parse::<i32>().unwrap();
        let resname = atomline[5..8].trim().to_string();
        let atname = atomline[8..15].trim().to_string();
        let atid = atomline[15..20].trim().parse::<i32>().unwrap();
        
        let x = atomline[20..28].trim().parse::<f64>().unwrap();
        let y = atomline[28..36].trim().parse::<f64>().unwrap();
        let z = atomline[36..44].trim().parse::<f64>().unwrap();
        
        let (vx, vy, vz) = if atomline.len() > 44 {
            let vx = if atomline.len() >= 52 {
                atomline[44..52].trim().parse::<f64>().unwrap_or(0.0)
            } else { 0.0 };
            let vy = if atomline.len() >= 60 {
                atomline[52..60].trim().parse::<f64>().unwrap_or(0.0)
            } else { 0.0 };
            let vz = if atomline.len() >= 68 {
                atomline[60..68].trim().parse::<f64>().unwrap_or(0.0)
            } else { 0.0 };
            (vx, vy, vz)
        } else {
            (0.0, 0.0, 0.0)
        };
        
        GroAtom {
            resid,
            resname,
            atname,
            atid,
            x,
            y,
            z,
            vx,
            vy,
            vz,
        }
    }
}

impl fmt::Display for GroAtom {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:5}{:4} {:6}{:5}{:8.3}{:8.3}{:8.3}{:8.4}{:8.4}{:8.4}",
            self.resid,
            self.resname,
            self.atname,
            self.atid,
            self.x,
            self.y,
            self.z,
            self.vx,
            self.vy,
            self.vz
        )
    }
}
