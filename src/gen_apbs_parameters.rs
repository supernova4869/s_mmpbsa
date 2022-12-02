use crate::apbs_param::*;
use crate::Parameters;

pub fn dim_apbs(file: &str, mol_index: i32, min_x: f64, max_x: f64, min_y: f64, max_y: f64, min_z: f64, max_z: f64,
            settings: &Parameters, pbe_set: &PBESet, pbe_set0: &PBESet, pba_set: &PBASet) -> String {
    let cfac = settings.cfac;
    let fadd = settings.fadd;
    let df = settings.df;
    let grid_type = settings.grid_type;

    // convert to A
    let min_x = min_x * 10.0;
    let min_y = min_y * 10.0;
    let min_z = min_z * 10.0;
    let max_x = max_x * 10.0;
    let max_y = max_y * 10.0;
    let max_z = max_z * 10.0;

    let x_len = (max_x - min_x).max(0.1);
    let x_center = (max_x + min_x) / 2.0;
    let y_len = (max_y - min_y).max(0.1);
    let y_center = (max_y + min_y) / 2.0;
    let z_len = (max_z - min_z).max(0.1);
    let z_center = (max_z + min_z) / 2.0;

    let mut c_x = 0.0;
    let mut c_y = 0.0;
    let mut c_z = 0.0;
    let mut f_x = 0.0;
    let mut f_y = 0.0;
    let mut f_z = 0.0;
    let mut n_x = 0;
    let mut n_y = 0;
    let mut n_z = 0;

    let split_level = 4;
    let t: i32 = 2_i32.pow(split_level + 1);

    match grid_type {
        0 => {
            let fpre = 1;
            let cfac = 1.7;
            f_x = x_len + 2.0 * fadd;
            c_x = f_x * cfac;
            n_x = t * ((f_x / (t as f64 * df)).round() as i32 + 1 + fpre) + 1;
            f_y = y_len + 2.0 * fadd;
            c_y = f_y * cfac;
            n_y = t * ((f_y / (t as f64 * df)).round() as i32 + 1 + fpre) + 1;
            f_z = z_len + 2.0 * fadd;
            c_z = f_z * cfac;
            n_z = t * ((f_z / (t as f64 * df)).round() as i32 + 1 + fpre) + 1;
        }
        _ => {
            c_x = x_len * cfac;
            f_x = (x_len + fadd).min(c_x);
            c_y = y_len * cfac;
            f_y = (y_len + fadd).min(c_y);
            c_z = z_len * cfac;
            f_z = (z_len + fadd).min(c_z);
            n_x = (f_x / df).round() as i32 - 1;
            n_x = ((n_x as f64 / t as f64).round() as i32 * t + 1).max(33);
            n_y = (f_y / df).round() as i32 - 1;
            n_y = ((n_y as f64 / t as f64).round() as i32 * t + 1).max(33);
            n_z = (f_z / df).round() as i32 - 1;
            n_z = ((n_z as f64 / t as f64).round() as i32 * t + 1).max(33);
        }
    }
    let mg_set = "mg-auto";

    let xyz_set = format!("  {mg_set}\n  mol {mol_index}\
        \n  dime   {n_x}  {n_y}  {n_z}\
        \n  cglen  {c_x:.3}  {c_y:.3}  {c_z:.3}\
        \n  fglen  {f_x:.3}  {f_y:.3}  {f_z:.3}\
        \n  fgcent {x_center:.3}  {y_center:.3}  {z_center:.3}\
        \n  cgcent {x_center:.3}  {y_center:.3}  {z_center:.3}\n");

    return format!("\nELEC name {}\n\
    {} \n\
    {} \n\
    end\n\n\
    ELEC name {}_VAC\n\
    {}\n\
    {} \n\
    end\n\n\
    APOLAR name {}_SAS\n  \
    mol {}\n{}\n\
    end\n\n\
    print elecEnergy {} - {}_VAC end\n\
    print apolEnergy {}_SAS end\n\n", file, xyz_set, pbe_set.to_string(), file,
                   xyz_set, pbe_set0.to_string(), file, mol_index, pba_set.to_string(),
                   file, file, file);
}