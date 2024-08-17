use std::path::Path;
use std::fs;
use ndarray::Array3;

pub fn read_coord_xvg(wd: &Path, coord_fname: &str) -> (Vec<f64>, Array3<f64>) {
    let xvg = fs::read_to_string(wd.join(coord_fname).to_str().unwrap()).unwrap();
    let xvg: Vec<&str> = xvg.split("\n").filter_map(|s| 
        if !s.is_empty() && !s.trim().starts_with("@") && !s.trim().starts_with("#") { 
            Some(s)
        } else {
            None
        }
    ).collect();
    let mut time_list = vec![];
    let mut coordinates = vec![];
    for coord in &xvg {
        let mut coord_ts: Vec<f64> = coord.split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        time_list.push(coord_ts.remove(0));
        coordinates.append(&mut coord_ts);
    }
    (time_list, Array3::from_shape_vec((xvg.len(), coordinates.len() / xvg.len() / 3, 3), coordinates).unwrap() * 10.0)
}