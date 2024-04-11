use std::path::Path;
use std::io;
use std::str::FromStr;
use std::fmt::Debug;

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
    let mut inp: String = String::new();
    io::stdin().read_line(&mut inp).expect("Failed to read line");
    let inp: Vec<&str> = inp.split("#").collect();
    let inp = inp[0].trim();
    match inp.is_empty() {
        true => default,
        false => inp.parse().expect("Failed to parse input")
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
                println!("Error input, input again.");
                continue;
            }
        };
    }
}

pub fn get_outfile<T: ToString + std::fmt::Display>(default_name: &T) -> String {
    println!("\nInput file name to write (default: {}):", default_name);
    let mut temp = String::new();
    io::stdin().read_line(&mut temp).unwrap();
    match temp.trim().is_empty() {
        true => default_name.to_string(),
        _ => temp.trim().to_string()
    }
}

pub fn append_new_name(origin_name: &str, append_name: &str) -> String {
    let file_path = Path::new(origin_name);
    let file_stem = file_path.file_stem().unwrap();
    let new_name = file_path.parent().unwrap().join(file_stem.to_str().unwrap().to_string() + append_name);
    new_name.to_str().unwrap().to_string()
}
