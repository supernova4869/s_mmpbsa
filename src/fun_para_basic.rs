use std::io::stdin;
use std::path::Path;
use crate::settings::Settings;
use crate::utils::get_input_selection;
use crate::{convert_cur_dir, confirm_file_validity, check_apbs};
use crate::fun_para_trj::set_para_trj;
use crate::parse_tpr::TPR;

pub fn set_para_basic(trj: &String, tpr: &mut TPR, ndx: &String, wd: &Path, settings: &mut Settings) {
    let mut trj = String::from(trj);
    let mut ndx = String::from(ndx);

    loop {
        println!("\n                 ************ MM/PB-SA Files ************");
        println!("-10 Exit program");
        println!(" -1 Set apbs path, current: {}", match &settings.apbs {
            Some(s) => s.to_string(),
            None => String::from("Not set")
        });
        println!("  0 Go to next step");
        println!("  1 Assign trajectory file (xtc or trr), current: {}", match trj.len() {
            0 => "undefined",
            _ => trj.as_str()
        });
        println!("  2 Assign index file (ndx), current: {}", match ndx.len() {
            0 => "undefined",
            _ => ndx.as_str()
        });
        let i = get_input_selection();
        match i {
            -1 => {
                println!("Input APBS path (if empty, means do not do PBSA calculation):");
                let s: String = get_input_selection();
                match check_apbs(Some(s)) {
                    Some(s) => settings.apbs = Some(s),
                    None => settings.apbs = None
                }
            }
            0 => {
                if trj.len() == 0 {
                    println!("Trajectory file not assigned.");
                } else if ndx.len() == 0 {
                    // 可能要改, 以后不需要index也能算
                    println!("Index file not assigned.");
                } else {
                    let trjwho = trj.to_string() + "trj_whole.xtc";

                    // echo "$com" | $trjconv -s $tpr -n $idx -f $trj -o $trjwho -pbc whole &>>$err
                    // [[ ! -f "$trjwho" ]] && { echo -e "!!! ERROR !!! gmx trjconv Failed ! Check $err to figure out why.\n"; exit; }

                    // echo "$com" | $converttpr -s $tpr -n $idx -o $tpx &>>$err

                    // awk >$idx -v ndx=$ndx -v pro=$pro -v lig=$lig -v withLig=$withLig '
                    //     BEGIN { RS="["
                    //         while(getline < ndx) {
                    //             gsub(" ","", $1); gsub("\t","", $1)
                    //             if($1==pro)    { ipro=$3; npro=NF-2; }
                    //             if($1==pro"]") { ipro=$2; npro=NF-1; }
                    //             if(withLig) {
                    //                 if($1==lig)    { ilig=$3; nlig=NF-2; }
                    //                 if($1==lig"]") { ilig=$2; nlig=NF-1; }
                    //             }
                    //         }
                    //         if(ipro<ilig) { dpro=0; dlig=npro }
                    //         else          { dlig=0; dpro=nlig }
                    //         print "[ "pro" ]";       for(i=1; i<=npro; i++)      printf "%d%s", i+dpro, i%15==0?"\n":" "; print ""
                    //         print "[ "lig" ]";       for(i=1; i<=nlig; i++)      printf "%d%s", i+dlig, i%15==0?"\n":" "; print ""
                    //         print "[ "pro"_"lig" ]"; for(i=1; i<=npro+nlig; i++) printf "%d%s", i,      i%15==0?"\n":" "; print ""
                    //     } '

                    // trjcnt=_$pid~trj_cnt.xtc; trjcls=_$pid~trj_cls.xtc
                    // if [[ $withLig -eq 1 ]]; then
                    // # usful for single protein and ligand
                    //     echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $trjcnt &>>$err -pbc mol -center
                    //     echo -e "$com\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
                    //     echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjcls -o $pdb    &>>$err -fit rot+trans
                    // else
                    //     echo -e "$lig\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $pdb    &>>$err -pbc mol -center
                    // fi
                    // echo -e ">> 1. preprocessing trajectory: OK !\n"

                    // go to next step
                    set_para_trj(&trj, tpr, &ndx, &wd, settings);
                }
            }
            1 => {
                println!("Input trajectory file path, default: ?md.xtc (if in the same directory with tpr, then simply input (e.g.) `?md.xtc`):");
                trj.clear();
                stdin().read_line(&mut trj).expect("Failed while reading trajectory file");
                if trj.trim().is_empty() {
                    trj = "?md.xtc".to_string();
                }
                trj = convert_cur_dir(&trj, &settings);
                trj = confirm_file_validity(&mut trj, vec!["xtc", "trr"], &settings);
            }
            2 => {
                println!("Input index file path, default: ?index.ndx (if in the same directory with tpr, then simply input (e.g.) `?index.ndx`):");
                ndx.clear();
                stdin().read_line(&mut ndx).expect("Failed while reading index file");
                if ndx.trim().is_empty() {
                    ndx = "?index.ndx".to_string();
                }
                ndx = convert_cur_dir(&ndx, &settings);
                ndx = confirm_file_validity(&mut ndx, vec!["ndx"], &settings);
            }
            -10 => break,
            _ => println!("Error input.")
        };
    }
}