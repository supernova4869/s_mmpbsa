use std::collections::HashMap;
use std::path::Path;
use std::io::{self, Write};
use std::str::FromStr;
use std::fmt::Debug;
use std::process::{Command, ExitStatus, Stdio};
use ndarray::{Array2, Axis};
use crate::parse_tpr::Residue;
use crate::settings::Settings;

pub fn range2list(range_str: &str) -> Vec<i32> {
    let mut selection_range: Vec<i32> = vec![];
    if range_str.trim().is_empty() {
        return selection_range;
    }
    let range_str = range_str.replace(" ", "");
    let sub_selections: Vec<&str> = range_str.split(',').collect();
    for s in sub_selections {
        if s.contains('-') {
            let r: Vec<&str> = s.split('-').collect();
            let l: i32 = r[0].parse().expect("Range lower bound not int");
            let u: i32 = r[1].parse().expect("Range upper bound not int");
            selection_range.append(&mut (l..=u).collect());
        } else {
            let s: i32 = s.parse().expect("Index not int");
            selection_range.push(s);
        }
    }
    return selection_range
}

pub fn get_input<T: FromStr>(default: T) -> T where <T as FromStr>::Err: Debug {
    loop {
        let mut inp: String = String::new();
        io::stdin().read_line(&mut inp).expect("Failed to read line");
        let inp: Vec<&str> = inp.split("#").collect();
        let inp = inp[0].trim();
        match inp.is_empty() {
            true => return default,
            false => match inp.parse() {
                Ok(a) => return a,
                Err(_) => {
                    continue
                }
            }
        }
    }
}

pub fn get_input_selection<T: FromStr>() -> T {
    loop {
        let mut input = String::from("");
        io::stdin().read_line(&mut input).expect("Error input.");
        let input: Vec<&str> = input.split("#").collect();
        let input = input[0].trim();
        match input.parse() {
            Ok(num) => return num,
            Err(_) => {
                continue;
            }
        };
    }
}

pub fn append_new_name(origin_name: &str, append_name: &str, prefix: &str) -> String {
    let file_path = Path::new(origin_name);
    let file_stem = file_path.file_stem().unwrap();
    let new_name = file_path.parent().unwrap().join(prefix.to_string() + file_stem.to_str().unwrap() + append_name);
    new_name.to_str().unwrap().to_string()
}

fn cmd_options(settings: &Settings, cmd: &str, options: &Vec<&str>, args: &[&str], wd: &Path) -> Result<ExitStatus, std::io::Error> {
    if settings.debug_mode {
        println!("CMD: {} {}", cmd.to_string(), args.join(" "));
    }
    let mut child = Command::new(cmd)
        .args(args)
        .current_dir(wd)
        .stdin(Stdio::piped())  // 开启标准输入管道
        .stdout(if settings.debug_mode { Stdio::inherit() } else { Stdio::null() })  // 将标准输出继承自父进程
        .stderr(if settings.debug_mode { Stdio::inherit() } else { Stdio::null() })  // 将标准输出继承自父进程
        .spawn()
        .expect("Failed to start process");

    // 获取stdin的可写句柄
    if let Some(stdin) = child.stdin.as_mut() {
        options.iter().for_each(|s| {
            if settings.debug_mode {
                println!("Input: {}", s);
            }
            writeln!(stdin, "{}", s).unwrap()
        });
    }

    // 等待子进程完成
    child.wait()
}

pub fn pdb2gmx(options: &Vec<&str>, wd: &Path, settings: &Settings, f: &str, o: &str, ff: &str, water: &str) {
    let args = ["pdb2gmx", "-f", f, "-o", o, "-ff", ff, "-water", water, "-ignh"];
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn grompp(options: &Vec<&str>, wd: &Path, settings: &Settings, f: &str, c: &str, o: &str) {
    let args = ["grompp", "-f", f, "-c", c, "-o", o, "-maxwarn", "10"];
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn convert_tpr(options: &Vec<&str>, wd: &Path, settings: &Settings, s: &str, n: &str, o: &str) {
    let args = ["convert-tpr", "-s", s, "-n", n, "-o", o];
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn convert_trj(options: &Vec<&str>, wd: &Path, settings: &Settings, f: &str, s: &str, n: &str, o: &str, others: &[&str]) {
    let args: Vec<&str> = ["convert-trj", "-f", f, "-s", s, "-n", n, "-o", o].iter().chain(others.iter()).cloned().collect();
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn trjconv(options: &Vec<&str>, wd: &Path, settings: &Settings, f: &str, s: &str, n: &str, o: &str, others: &[&str]) {
    let args: Vec<&str> = ["trjconv", "-f", f, "-s", s, "-n", n, "-o", o].iter().chain(others.iter()).cloned().collect();
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn make_ndx(options: &Vec<&str>, wd: &Path, settings: &Settings, f: &str, n: &str, o: &str) {
    let args = match n.is_empty() {
        true => ["make_ndx", "-f", f, "-o", o].to_vec(),
        false => ["make_ndx", "-f", f, "-n", n, "-o", o].to_vec()
    };
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn trajectory(options: &Vec<&str>, wd: &Path, settings: &Settings, f: &str, s: &str, n: &str, ox: &str) {
    let args = ["trajectory", "-f", f, "-s", s, "-n", n, "-ox", ox].to_vec();
    cmd_options(settings, settings.gmx_path.as_ref().unwrap(), options, &args, wd).unwrap();
}

pub fn sobtop(options: &Vec<&str>, settings: &Settings, infile: &str) {
    let args = vec![infile];
    let sobtop_dir = Path::new(settings.sobtop_dir.as_ref().unwrap());
    // fuck, sobtop must be used at its own directory
    cmd_options(settings, sobtop_dir.join("sobtop").to_str().unwrap(), options, &args, &sobtop_dir).unwrap();
}

pub fn multiwfn(options: &Vec<&str>, settings: &Settings, infile: &str, wd: &Path) -> Result<ExitStatus, std::io::Error> {
    let args = vec![infile];
    let multiwfn_dir = Path::new(settings.multiwfn_dir.as_ref().unwrap());
    cmd_options(settings, multiwfn_dir.join("Multiwfn").to_str().unwrap(), options, &args, wd)
}

pub fn resname_3to1(name: &str) -> Option<String> {
    let mut resname_map: HashMap<&str, &str> = HashMap::new();
    resname_map.insert("ALA", "A");
    resname_map.insert("CYS", "C");
    resname_map.insert("ASP", "D");
    resname_map.insert("ASH", "D");
    resname_map.insert("GLU", "E");
    resname_map.insert("GLH", "E");
    resname_map.insert("PHE", "F");
    resname_map.insert("GLY", "G");
    resname_map.insert("HIS", "H");
    resname_map.insert("HID", "H");
    resname_map.insert("HIE", "H");
    resname_map.insert("HIP", "H");
    resname_map.insert("ILE", "I");
    resname_map.insert("LYS", "K");
    resname_map.insert("LEU", "L");
    resname_map.insert("MET", "M");
    resname_map.insert("ASN", "N");
    resname_map.insert("PRO", "P");
    resname_map.insert("GLN", "Q");
    resname_map.insert("ARG", "R");
    resname_map.insert("SER", "S");
    resname_map.insert("THR", "T");
    resname_map.insert("VAL", "V");
    resname_map.insert("TRP", "W");
    resname_map.insert("TYR", "Y");
    if let Some(s) = resname_map.get(name) {
        Some(s.to_string())
    } else {
        None
    }
}

pub fn get_residue_range_ca(coord: &Array2<f64>, ref_ids: &Vec<usize>, cutoff: f64, 
        atom_res: &Vec<usize>, atom_names: &Vec<String>, residues: &Vec<Residue>) -> Vec<usize> {
    let mut res_range: Vec<usize> = vec![];
    let ligand_coord = coord.select(Axis(0), ref_ids);
    for res in residues {
        let cur_res_ca_id: Vec<usize> = atom_names.iter().enumerate()
            .filter_map(|(id, name)| if atom_res[id] == res.id && name.eq("CA") {
                Some(id)
            } else {
                None
            })
            .collect();
        let receptor_residue_coord = coord.select(Axis(0), &cur_res_ca_id);
        for ri in ligand_coord.rows() {
            if (&receptor_residue_coord - &ri).map(|&i| i.powi(2)).sum() <= cutoff.powi(2) {
                res_range.push(res.id);
                break;
            }
        }
    }
    res_range
}
