use std::fs;
use regex::Regex;

pub struct IndexGroup {
    pub name:String,
    pub indexes:Vec<i32>
}

pub struct Index {
    pub groups:Vec<IndexGroup>
}

impl Index {
    pub fn new(index_file:&String) -> Index {
        let ndx = fs::read_to_string(index_file).expect("Failed reading index file");
        let re = Regex::new(r"\[\s*(.+?)\s*]").unwrap();
        let mut group_names: Vec<String> = vec![];          // name of each group
        for cap in re.captures_iter(&ndx) {
            group_names.push((&cap[1]).to_string());
        }
        let mut group_atoms:Vec<Vec<i32>> = vec![];         // atoms of each group
        let ndx = re.replace_all(&ndx, "[]");
        let ndx = Regex::new(r"\s+").unwrap().replace_all(&ndx, " ");
        let ndx: Vec<&str> = ndx.split("[]").collect();
        for i in 0..ndx.len() {
            if ndx[i].len() != 0 {
                let atom_list:Vec<&str> = ndx[i].trim().split(" ").collect();
                let mut at_list:Vec<i32> = vec![];
                for atom in atom_list {
                    at_list.push(atom.parse().expect("Failed to get index atom"));
                }
                group_atoms.push(at_list);
            }
        }
        let mut groups: Vec<IndexGroup> = vec![];
        while group_names.len() > 0 {
            groups.insert(0, IndexGroup {
                name: group_names.pop().expect("Error reading index file."),
                indexes: group_atoms.pop().expect("Error reading index file.")
            });
        }
        let index = Index{ groups: groups };
        return index;
    }
}