mod index_parser;

use std::fs;
use std::env;
use std::io::{Read, stdin, Write};
use std::path::Path;
use std::process::Command;

fn main() {
    //parameters
    let gmx = "gmx";
    let args: Vec<String> = env::args().collect();
    let mut tpr = String::new();
    let mut trj = String::from("");
    let mut ndx = String::from("");
    let use_dh = true;
    let use_ts = true;

    // start workflow
    welcome();
    match args.len() {
        1 => {
            println!("Input path of .tpr file, e.g. D:/Study/ZhangYang.tpr");
            stdin().read_line(&mut tpr).expect("Failed to read tpr file.");
        },
        2 => tpr = args[1].to_string(),
        _ => {
            for i in 1..args.len() {
                match args[i].as_str() {
                    "-f" => { trj = args[i + 1].to_string() },
                    "-s" => { tpr = args[i + 1].to_string() },
                    "-n" => { ndx = args[i + 1].to_string() },
                    _ => {
                        if i % 2 == 1 {
                            println!("Omitted invalid option: {}", args[i])
                        }
                    }
                }
            }
        },
    }
    let mut tpr_path = Path::new(tpr.trim());
    while !tpr_path.is_file() {
        println!("Not valid tpr file. Input path of .tpr file again.");
        stdin().read_line(&mut tpr).expect("Failed to read tpr file.");
        tpr_path = Path::new(tpr.trim());
    }
    // working directory (path of tpr location)
    let wd = tpr_path.parent().expect("Failed getting parent directory.");
    println!("Currently working at path: {}", wd.display());
    dump_tpr(&tpr, wd, gmx);
    loop {
        println!("\n                 ************ SuperMMPBSA functions ************");
        println!("-2 Toggle whether to use entropy contribution, current: {}", use_ts);
        println!("-1 Toggle whether to use Debye-Huckel shielding method, current: {}", use_dh);
        println!(" 0 Do MM-PBSA calculations now!");
        println!(" 1 Assign trajectory file (xtc or trr), current: {}", match trj.len() {
            0 => "undefined",
            _ => trj.as_str()
        });
        println!(" 2 Assign index file (ndx), current: {}", match ndx.len() {
            0 => "undefined",
            _ => ndx.as_str()
        });
        println!(" 3 Exit program");
        let mut i:String = String::from("");
        stdin().read_line(&mut i).expect("Error input");
        let i:i32 = i.trim().parse().expect("Error input");
        match i {
            -2 => { let use_dh = !use_dh; },
            -1 => { let use_ts = !use_ts; },
            0 => {
                if trj.len() == 0 {
                    println!("Trajectory file not assigned.");
                } else if ndx.len() == 0 {
                    // 可能要改, 以后不需要index也能算
                    println!("Index file not assigned.");
                } else {
                    calc_mmpbsa(&trj, &tpr, &ndx, use_dh, use_ts);
                }
            },
            1 => {
                println!("Input trajectory file path:");
                stdin().read_line(&mut trj).expect("Failed while reading trajectory");
                let p = Path::new(trj.trim());
                let ext = p.extension().unwrap().to_str().expect("");
                if !p.is_file() {
                    println!("Error, not file.");
                    trj = "".to_string();
                } else if !(ext.eq("xtc")) && !(ext.eq("trr")) {
                    println!("Error, not xtc/trr extension.");
                    trj = "".to_string();
                } else {
                    trj = trj.trim().to_string();
                }
            },
            2 => {
                println!("Input index file path:");
                stdin().read_line(&mut ndx).expect("Failed while reading index");
                let p = Path::new(ndx.trim());
                let ext = p.extension().unwrap().to_str().expect("");
                if !p.is_file() {
                    println!("Error, not file.");
                    ndx = "".to_string();
                } else if !(ext == "ndx") {
                    println!("Error, not ndx extension.");
                    ndx = "".to_string();
                } else {
                    ndx = ndx.trim().to_string();
                }
            },
            3 => break,
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
Usage 1: run `SuperMMPBSA` and follow the prompts.\n\
Usage 2: run `SuperMMPBSA WangBingBing.tpr` to directly load WangBingBing.tpr.\n\
Usage 3: run `SuperMMPBSA -f md.xtc -s md.tpr -n index.ndx` to assign all needed files.\n");
}

fn dump_tpr(tpr:&String, wd:&Path, gmx:&str) {
    // gmx = settings["environments"]["gmx"];
    let tpr_dump = Command::new(gmx).arg("dump").arg("-s").arg(tpr).output().expect("gmx dump failed.");
    let tpr_dump = String::from_utf8(tpr_dump.stdout).expect("Getting dump output failed.");
    let mut outfile = fs::File::create(wd.join("_mdout.mdp")).unwrap();
    outfile.write(tpr_dump.as_bytes()).unwrap();
    println!("Finished loading tpr file, md parameters dumped to {}", wd.join("_mdout.mdp").display());
}

fn calc_mmpbsa(trj:&String, tpr:&String, ndx:&String, use_dh:bool, use_ts:bool) {
    // TODO: 选组
    // let mut ligand_grp = -1;
    // let mut receptor_grp = -1;
    // let mut complex_grp = -1;
    let ndx = index_parser::Index::new(ndx);
    println!("{:?}", ndx.groups[0].indexes);
    // ndx = index_parser.Index(Vector{IndexParser.IndexGroup}())
    // 这部分留到第二步, 因为后面可能要修改选择原子的规则
    // println!(" 3 Select ligand groups, current: {}", match ligand_grp {
    //     -1 => "undefined",
    //     _ => format("{} {}", ligand_grp, ndx.groups[ligand_grp + 1].name)
    // };
    // println!(" 4 Select receptor groups, current: $(receptor_grp != -1 ? "$receptor_grp " * ndx.groups[receptor_grp + 1].name : "undefined")");
    // println!(" 5 Select complex groups, current: $(complex_grp != -1 ? "$complex_grp " * ndx.groups[complex_grp + 1].name : "undefined")");
    // println!("Select groups and do calculations.");
}