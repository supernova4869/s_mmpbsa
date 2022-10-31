use std::fs;
use regex::Regex;

pub struct IndexGroup {
    pub name: String,
    pub indexes: Vec<i32>,      // starts at 0
}

pub struct Index {
    pub groups: Vec<IndexGroup>,
}

impl Index {
    pub fn new(index_file: &String) -> Index {
        let ndx = fs::read_to_string(index_file).expect("Failed reading index file");
        let re = Regex::new(r"\[\s*(.+?)\s*]").unwrap();
        let mut group_names: Vec<String> = vec![];          // name of each group
        for cap in re.captures_iter(&ndx) {
            group_names.push((&cap[1]).to_string());
        }
        let mut group_atoms: Vec<Vec<i32>> = vec![];         // atoms of each group
        let ndx = re.replace_all(&ndx, "[]");
        let ndx = Regex::new(r"\s+").unwrap().replace_all(&ndx, " ");
        let ndx: Vec<&str> = ndx.split("[]").collect();
        for i in 0..ndx.len() {
            if ndx[i].len() != 0 {
                let atom_list: Vec<&str> = ndx[i].trim().split(" ").collect();
                let mut at_list: Vec<i32> = vec![];
                for atom in atom_list {
                    let at: i32 = atom.parse().expect("Failed to get index atom");
                    // gromacs index file starts at 1, but we need 0 as index slice
                    at_list.push(at - 1);
                }
                group_atoms.push(at_list);
            }
        }
        let mut groups: Vec<IndexGroup> = vec![];
        while group_names.len() > 0 {
            groups.insert(0, IndexGroup {
                name: group_names.pop().expect("Error reading index file."),
                indexes: group_atoms.pop().expect("Error reading index file."),
            });
        }
        let index = Index { groups };
        return index;
    }

    pub fn list_groups(&self) {
        for i in 0..self.groups.len() {
            println!("{:>3}): {:<-15}{:>10} atoms", i, self.groups[i].name, self.groups[i].indexes.len())
        }
    }
}