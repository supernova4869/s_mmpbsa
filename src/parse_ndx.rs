use std::{fs::File, io::Write, path::Path};
use std::io::{BufRead, BufReader};
use std::fmt::Formatter;
use std::fmt;
use std::collections::BTreeSet;

#[derive(Clone)]
pub struct IndexGroup {
    pub name: String,
    pub indexes: BTreeSet<usize>,      // starts at 0
}

impl IndexGroup {
    pub fn new(name: &str, index_set: &BTreeSet<usize>) -> IndexGroup {
        IndexGroup { name: name.to_string(), indexes: index_set.to_owned() }
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
        let index_file = File::open(index_file).unwrap();
        let reader = BufReader::new(index_file);

        let mut groups = Vec::new();
        let mut current_name: Option<String> = None;
        let mut current_indexes = BTreeSet::new();

        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();

            // 跳过空行
            if line.is_empty() {
                continue;
            }

            // 检查是不是组标题，格式 [ GroupName ]
            if line.starts_with('[') && line.ends_with(']') {
                // 保存上一个组
                if let Some(name) = current_name.take() {
                    groups.push(IndexGroup {
                        name,
                        indexes: std::mem::take(&mut current_indexes),
                    });
                }

                // 提取组名，去掉前后的括号和空格
                let name = line[1..line.len() - 1].trim().to_string();
                current_name = Some(name);
            } else {
                // 数字行，解析所有数字（可能跨多行）
                if current_name.is_some() {
                    for num_str in line.split_whitespace() {
                        if let Ok(num) = num_str.parse::<usize>() {
                            // 转换成 0-based
                            if num > 0 {
                                current_indexes.insert(num - 1);
                            }
                        }
                    }
                }
            }
        }

        // 保存最后一组
        if let Some(name) = current_name {
            groups.push(IndexGroup {
                name,
                indexes: current_indexes,
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

    pub fn to_ndx<P: AsRef<Path>>(&self, file_name: P) {
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