use std::fs;
use std::process::Command;
use std::io::{stdin, Write};
use std::path::Path;

fn main() {
    welcome();
    let mut tpr = String::new();
    stdin().read_line(&mut tpr).expect("Failed to read tpr file.");
    let tpr = tpr.trim();
    load_tpr(tpr);
}

fn welcome() {
    println!("SuperMmpbsa: Supernova's tool of calculating binding free energy using\n\
molecular mechanics Poisson–Boltzmann surface area (MM-PBSA) method.\n\
Website: currently not available.\n\
Developed by Jiaxing Zhang (zhangjiaxing7137@tju.edu.cn), Tian Jin University.\n\
Version 0.1, first release: 2022-Oct-17\n\n\
Input path of .tpr file, e.g. D:\\ZhangYang\\md.tpr");
}

fn load_tpr(tpr:&str) {
    // settings = get_settings();
    // gmx = settings["environments"]["gmx"];
    let mdout = Command::new("gmx").arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let mdout = String::from_utf8(mdout.stdout).expect("Convert dump results failed.");
    let outpath = Path::new(tpr).parent().expect("").join("_mdout.mdp");
    let mut outfile = fs::File::create(outpath).unwrap();
    outfile.write(mdout.as_bytes()).unwrap();
}