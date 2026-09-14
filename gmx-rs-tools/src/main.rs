//! `gmx-rs-tools`: a multicall binary exposing the four extracted tools.

use std::process::ExitCode;

fn main() -> ExitCode {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        print_usage();
        return ExitCode::from(1);
    }
    let command = args.remove(0);
    let code = match command.as_str() {
        "dump" => gmx_rs_tools::cmd::dump::run(args),
        "coords" => gmx_rs_tools::cmd::coords::run(args),
        "trjconv" => gmx_rs_tools::cmd::trjconv::run(args),
        "convert-tpr" | "convert_tpr" => gmx_rs_tools::cmd::convert_tpr::run(args),
        "make_ndx" | "make-ndx" => gmx_rs_tools::cmd::make_ndx::run(args),
        "-h" | "--help" | "help" => {
            print_usage();
            0
        }
        other => {
            eprintln!("gmx-rs-tools: unknown command '{other}'");
            print_usage();
            1
        }
    };
    ExitCode::from(code as u8)
}

fn print_usage() {
    eprintln!(
        "gmx-rs-tools - minimal Rust port of selected GROMACS tools\n\n\
         Usage: gmx-rs-tools <command> [options]\n\n\
         Commands:\n\
         \x20 dump         Make binary files human readable\n\
         \x20 coords       Read coordinates out of a trajectory\n\
         \x20 trjconv      Convert and manipulate trajectories\n\
         \x20 convert-tpr  Make a modified run-input file\n\
         \x20 make_ndx     Make index files\n"
    );
}
