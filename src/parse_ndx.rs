use std::{fs::{self, File}, io::Write};
use regex::Regex;
use std::fmt::Formatter;
use std::fmt;
use std::collections::BTreeSet;

#[derive(Clone)]
pub struct IndexGroup {
    pub name: String,
    pub indexes: BTreeSet<usize>,      // starts at 0
}

impl IndexGroup {
    pub fn new(name: &str, list: &BTreeSet<usize>) -> IndexGroup {
        IndexGroup { name: name.to_string(), indexes: list.to_owned() }
    }
}

impl fmt::Display for IndexGroup {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}, {} atoms", 
            self.name,
            self.indexes.len())
    }
}

#[derive(Clone)]
pub struct Index {
    pub groups: Vec<IndexGroup>,
}

impl Index {
    pub fn new(groups: Vec<IndexGroup>) -> Index {
        Index{ groups }
    }

    pub fn from(index_file: &String) -> Index {
        let ndx = fs::read_to_string(index_file).expect("Failed reading index file");
        let re = Regex::new(r"\[\s*(.+?)\s*]").unwrap();
        let mut group_names: Vec<String> = vec![];          // name of each group
        for cap in re.captures_iter(&ndx) {
            group_names.push((&cap[1]).to_string());
        }
        let mut group_atoms: Vec<BTreeSet<usize>> = vec![];         // atoms of each group
        let ndx = re.replace_all(&ndx, "[]");
        let ndx = Regex::new(r"\s+").unwrap().replace_all(&ndx, " ");
        let ndx: Vec<&str> = ndx.split("[]").collect();
        for i in 0..ndx.len() {
            if !ndx[i].trim().is_empty() {
                let atom_list: BTreeSet<usize> = ndx[i].trim().split(" ").map(|a| a.parse().unwrap()).collect();
                // let atom_list: BTreeSet<usize> = ;
                group_atoms.push(atom_list.iter().map(|a| a - 1).collect());
            }
        }
        let mut groups: Vec<IndexGroup> = vec![];
        while group_names.len() > 0 {
            groups.insert(0, IndexGroup {
                name: group_names.pop().unwrap(),
                indexes: group_atoms.pop().unwrap(),
            });
        }
        return Index { groups }
    }

    pub fn list_groups(&self) {
        for i in 0..self.groups.len() {
            println!("{:>3}): {:<-15}{:>10} atoms", i, self.groups[i].name, self.groups[i].indexes.len())
        }
    }

    #[allow(dead_code)]
    pub fn push(&mut self, ng: &IndexGroup) {
        self.groups.push(ng.to_owned());
    }

    pub fn to_ndx(&self, file_name: &str) {
        let mut f = File::create(file_name).unwrap();
        for ig in &self.groups {
            writeln!(f, "[ {} ]", ig.name).unwrap();
            let vec: Vec<_> = ig.indexes.iter().collect();
            for chunk in vec.chunks(10) {
                for &num in chunk {
                    write!(f, "{:7} ", num + 1).unwrap();
                }
                writeln!(f).unwrap();
            }
        }
    }

    #[allow(dead_code)]
    pub fn rm_group(&mut self, name: &str) {
        self.groups.retain(|g| g.name.ne(name));
    }
}