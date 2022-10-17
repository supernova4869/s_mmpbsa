use std::fs;
use std::io::{Read, stdin, Write};
use std::path::Path;
use std::process::Command;

fn main() {
    welcome();
    let mut tpr = String::new();
    stdin().read_line(&mut tpr).expect("Failed to read tpr file.");
    let tpr = tpr.trim();
    let wd = Path::new(tpr).parent().expect(""); // working directory (path of tpr location)
    println!("Currently working at {}", wd.display());
    dump_tpr(tpr, wd);
    let ligand_grp = -1;
    let receptor_grp = -1;
    let complex_grp = -1;
    let trj = "";
    let ndx = "";
    // ndx = IndexParser.Index(Vector{IndexParser.IndexGroup}())
    let use_dh = true;
    let use_ts = true;
    loop {
        println!("\n                 ************ SuperMMPBSA functions ************");
        println!("-2 Toggle whether to use entropy contribution, current: {}", use_ts);
        println!("-1 Toggle whether to use Debye-Huckel shielding method, current: {}", use_dh);
        println!(" 0 Do MM-PBSA calculations now!");
        // println!(" 1 Assign trajectory file (xtc or trr), current: {}", isempty(trjfilename) ? "undefined" : trjfilename);
        // println!(" 2 Assign index file (ndx), current: $(isempty(ndxfilename) ? "undefined" : ndxfilename)");
        // println!(" 3 Select ligand groups, current: $(ligand_grp != -1 ? "$ligand_grp " * ndx.groups[ligand_grp + 1].name : "undefined")");
        // println!(" 4 Select receptor groups, current: $(receptor_grp != -1 ? "$receptor_grp " * ndx.groups[receptor_grp + 1].name : "undefined")");
        // println!(" 5 Select complex groups, current: $(complex_grp != -1 ? "$complex_grp " * ndx.groups[complex_grp + 1].name : "undefined")");
        println!(" 6 Exit program");
        let mut i:String = String::from("");
        stdin().read_line(&mut i).expect("Error input");
        let i:i32 = i.trim().parse().expect("Error input");
        match i {
            -2..=5 => println!("OK"),
            6 => break,
            _ => println!("Error input.")
        };
    }
}

fn welcome() {
    println!("SuperMMPBSA: Supernova's tool of calculating binding free energy using\n\
molecular mechanics Poisson–Boltzmann surface area (MM-PBSA) method.\n\
Website: https://github.com/supernovaZhangJiaXing/super_mmpbsa\n\
Developed by Jiaxing Zhang (zhangjiaxing7137@tju.edu.cn), Tian Jin University.\n\
Version 0.1, first release: 2022-Oct-17\n\n\
Hint: you can directly load .tpr file by command: `SuperMMPBSA WangBingBing.tpr`\n\
Input path of .tpr file, e.g. D:/Study/ZhangYang.tpr");
}

fn dump_tpr(tpr:&str, wd:&Path) {
    // gmx = settings["environments"]["gmx"];
    let tpr_dump = Command::new("gmx").arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(wd.join("_mdout.mdp")).unwrap();
    outfile.write(tpr_dump.as_bytes()).unwrap();
    println!("Finished reading tpr file, md parameters dumped to {}", wd.join("_mdout.mdp").display());
}