use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, stdin};
use std::path::Path;
use std::fmt::{self, Formatter};
use std::hash::{Hash, Hasher};

use regex::Regex;

use crate::parse_pdb::PDB;
use crate::utils;

// 定义ITP文件中的拓扑项类型
#[derive(Debug)]
#[derive(Clone)]
pub enum TopologySection {
    Atomtypes(HashSet<Atomtype>),
    Atoms(Vec<Atom>),
    Bonds(Vec<Bond>),
    Angles(Vec<Angle>),
    Dihedrals(Vec<Dihedral>),
    CMAP(Vec<CMAP>),
    Others(String, Vec<String>),
}

// 原子类型
#[derive(Debug, Clone)]
pub struct Atomtype {
    name: String,
    id: i32,
    mass: f64,
    charge: f64,
    ptype: char,
    sigma: f64,
    epsilon: f64,
}

// 实现 PartialEq
impl PartialEq for Atomtype {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

// 实现 Eq（Eq 没有方法，只是一个标记 trait）
impl Eq for Atomtype {}

// 实现 Hash
impl Hash for Atomtype {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

impl fmt::Display for Atomtype {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:>6}{:8} {:12.6}{:12.6} {:4}{:18.6E}{:18.6E}",
            self.name,
            self.id,
            self.mass,
            self.charge,
            self.ptype,
            self.sigma,
            self.epsilon,
        )
    }
}

// 原子结构
#[derive(Debug, Clone)]
pub struct Atom {
    index: u32,
    name: String,
    type_name: String,
    res_index: u32,
    res_name: String,
    charge: f64,
    mass: f64,
    other_fields: Vec<String>,
}

// 键结构
#[derive(Debug, Clone)]
pub struct Bond {
    atom1: u32,
    atom2: u32,
    type_: String,
    other_fields: Vec<String>,
}

// 角度结构
#[derive(Debug, Clone)]
pub struct Angle {
    atom1: u32,
    atom2: u32,
    atom3: u32,
    type_: String,
    other_fields: Vec<String>,
}

// 二面角结构
#[derive(Debug, Clone)]
pub struct Dihedral {
    atom1: u32,
    atom2: u32,
    atom3: u32,
    atom4: u32,
    type_: String,
    other_fields: Vec<String>,
}

// CMAP结构
#[derive(Debug, Clone)]
pub struct CMAP {
    atom1: u32,
    atom2: u32,
    atom3: u32,
    atom4: u32,
    atom5: u32,
    type_: String,
}

// 定义拓扑结构
#[derive(Debug)]
pub struct Topology {
    pub sections: Vec<TopologySection>,
}

impl Topology {
    // 创建新的拓扑结构
    pub fn from<P: AsRef<Path>>(input_path: P) -> Self {
        let sections = read_itp_file(input_path).expect("Failed to read input file");
        Topology {
            sections,
        }
    }

    // 交互式排序原子
    #[allow(dead_code)]
    pub fn sort_atoms_interactive(&mut self, output_path: &String) {
        // 添加循环直到用户输入有效选项
        let m: i32 = loop {
            println!(" 0 Input new orders directly.");
            println!(" 1 Load from atoms order txt file.");
            println!(" 2 Reorder from reference pdb file.");
            print!("Please select an option: ");
            io::stdout().flush().unwrap(); // 确保提示信息被立即显示
            let mut i = String::new();
            stdin().read_line(&mut i).expect("Failed to read input");

            match i.trim().parse::<i32>() {
                Ok(num) if (num == 0 || num == 1 || num == 2) => break num,
                Ok(_) => println!("Invalid option. Please enter 0 or 1."),
                Err(_) => println!("Invalid input. Please enter a number.")
            }
        };

        let order = match m {
            0 => {
                println!("INFO: if you want the initial sequence [1 2 3 4 5] to be resorted into [1 4 5 2 3],");
                println!("then simply input \"1 4 5 2 3\" or \"1, 4-5, 2-3\" are both OK (no quotes).");
                println!("Input the target sequence of topol items:");
                io::stdout().flush().unwrap(); // 确保提示信息被立即显示
                let mut input = String::new();
                stdin().read_line(&mut input).expect("Failed to read input");
                utils::range2list(&input).iter().map(|a| *a as u32).collect()
            },
            1 => {
                println!("INFO: if you want the initial sequence [1 2 3 4 5] to be resorted into [1 4 5 2 3],");
                println!("then simply input \"1 4 5 2 3\" or \"1, 4-5, 2-3\" are both OK (no quotes).");
                println!("Input the path of sequence file (containing the above information):");
                let mut file_path = String::new();
                stdin().read_line(&mut file_path).expect("Failed to read input");
                let file_path = file_path.trim();

                // 使用PathBuf确保跨平台路径兼容
                let path = std::path::PathBuf::from(file_path);

                match std::fs::read_to_string(path) {
                    Ok(content) => {
                        println!("{}", content);
                        utils::range2list(&content).iter().map(|a| *a as u32).collect()
                    },
                    Err(e) => {
                        println!("Failed to read sequence file: {}", e);
                        return;
                    }
                }
            },
            2 => {
                // 根据pdb生成新的原子顺序
                println!("Input the path of reference pdb file:");
                let mut file_path = String::new();
                stdin().read_line(&mut file_path).expect("Failed to read input");
                let file_path = file_path.trim();

                // 使用PathBuf确保跨平台路径兼容
                let path = std::path::PathBuf::from(file_path);

                // 解析PDB文件
                let pdb = PDB::from(path);

                // 获取第一个模型的原子
                let model = &pdb.models[0];
                // 同时获取PDB原子的名称和残基名称
                let pdb_atoms_info: Vec<(String, String)> = model.atoms.iter()
                    .map(|a| (a.atname.trim().to_string(), a.resname.trim().to_string()))
                    .collect();

                // 查找拓扑文件中对应的原子索引
                let mut order = Vec::new();
                if let Some(TopologySection::Atoms(topol_atoms)) = self.sections.iter().find(|s| matches!(s, TopologySection::Atoms(_))) {
                    for (pdb_atom_name, pdb_res_name) in &pdb_atoms_info {
                        // 在拓扑原子中查找同时匹配原子名称和残基名称的原子
                        if let Some(topol_atom) = topol_atoms.iter().find(|a| a.name == *pdb_atom_name && a.res_name == *pdb_res_name) {
                            order.push(topol_atom.index);
                        } else {
                            println!("Warning: Atom '{}' with residue '{}' from PDB not found in topology.", pdb_atom_name, pdb_res_name);
                        }
                    }
                }

                println!("Generated order from PDB: {:?}", order);
                order
            }
            _ => {
                println!("Error input.");
                return;
            }
        };

        if !order.is_empty() {
            // 更新所有依赖于原子索引的部分
            self.sections = reorder_atoms(self.sections.clone(), &order);
            println!("Sort topol atoms finished.");
            let output_path_buf = std::path::PathBuf::from(&output_path);
            if let Some(parent_dir) = output_path_buf.parent() {
                if !parent_dir.exists() {
                    std::fs::create_dir_all(parent_dir).expect("Failed to create output directory");
                    println!("Created output directory: {:?}", parent_dir);
                }
            }

            // 写入输出文件
            self.to_itp(output_path, true).expect("Failed to write output file");
            println!("Output written to {}", output_path);
        }
    }

    pub fn to_itp<P: AsRef<Path>>(&self, path: P, with_type: bool) -> io::Result<()> {
        let mut file = File::create(path)?;

        for section in &self.sections {
            match section {
                TopologySection::Atomtypes(atomtypes) => {
                    if with_type {
                        writeln!(file, "[ atomtypes ]")?;
                        writeln!(file, "; name   at.num      mass       charge   ptype     sigma (nm)    epsilon (kJ/mol)")?;
                        //; name   at.num      mass       charge   ptype     sigma (nm)    epsilon (kJ/mol)
                        //    n3       7    14.006703    0.000000    A      3.249999E-01    7.112800E-01
                        for atomtype in atomtypes {
                            writeln!(
                                file,
                                "{:>6}{:8} {:12.6}{:12.6} {:4}{:18.6E}{:18.6E}",
                                atomtype.name,
                                atomtype.id,
                                atomtype.mass,
                                atomtype.charge,
                                atomtype.ptype,
                                atomtype.sigma,
                                atomtype.epsilon,
                            )?;
                        }
                        writeln!(file)?;
                    }
                }, 
                TopologySection::Atoms(atoms) => {
                    writeln!(file, "[ atoms ]")?;
                    writeln!(file, "; Index   type   residue  resname   atom     cgnr        charge        mass  ...")?;
                    //  Index   type   residue  resname   atom        cgnr     charge       mass
                    //  1     c          1      MOL     C1            1    0.50802181   12.010736
                    for atom in atoms {
                        writeln!(
                            file,
                            "{:7}{:>7} {:9} {:>8} {:>6} {:8}{:14.8}{:12.6} {}",
                            atom.index,
                            atom.type_name,
                            atom.res_index,
                            atom.res_name,
                            atom.name,
                            atom.index,
                            atom.charge,
                            atom.mass,
                            atom.other_fields.join(" ")
                        )?;
                    }
                    writeln!(file)?;
                }, 
                TopologySection::Bonds(bonds) => {
                    writeln!(file, "[ bonds ]")?;
                    writeln!(file, "; ai  aj  type  ...")?;
                    for bond in bonds {
                        writeln!(
                            file,
                            "{:>6} {:>6} {:>8}   {}",
                            bond.atom1,
                            bond.atom2,
                            bond.type_,
                            bond.other_fields.join(" ")
                        )?;
                    }
                    writeln!(file)?;
                },
                TopologySection::Angles(angles) => {
                    writeln!(file, "[ angles ]")?;
                    writeln!(file, "; ai  aj  ak  type  ...")?;
                    for angle in angles {
                        writeln!(
                            file,
                            "{:>6} {:>6} {:>6} {:>8}   {}",
                            angle.atom1,
                            angle.atom2,
                            angle.atom3,
                            angle.type_,
                            angle.other_fields.join(" ")
                        )?;
                    }
                    writeln!(file)?;
                },
                TopologySection::Dihedrals(dihedrals) => {
                    writeln!(file, "[ dihedrals ]")?;
                    writeln!(file, "; ai  aj  ak  al  type  ...")?;
                    for dihedral in dihedrals {
                        writeln!(
                            file,
                            "{:>6} {:>6} {:>6} {:>6} {:>8}   {}",
                            dihedral.atom1,
                            dihedral.atom2,
                            dihedral.atom3,
                            dihedral.atom4,
                            dihedral.type_,
                            dihedral.other_fields.join(" ")
                        )?;
                    }
                    writeln!(file)?;
                },
                TopologySection::CMAP(cmaps) => {
                    writeln!(file, "[ cmap ]")?;
                    writeln!(file, "; ai  aj  ak  al  type  ...")?;
                    for cmap in cmaps {
                        writeln!(
                            file,
                            "{:>6} {:>6} {:>6} {:>6} {:>6} {:>8}",
                            cmap.atom1,
                            cmap.atom2,
                            cmap.atom3,
                            cmap.atom4,
                            cmap.atom5,
                            cmap.type_,
                        )?;
                    }
                    writeln!(file)?;
                },
                TopologySection::Others(name, content) => {
                    writeln!(file, "[ {} ]", name)?;
                    for line in content {
                        writeln!(file, "{}", line)?;
                    }
                    writeln!(file)?;
                },
            }
        }

        Ok(())
    }
}

// 读取ITP文件并解析为拓扑项
fn read_itp_file<P: AsRef<Path>>(path: P) -> io::Result<Vec<TopologySection>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sections = Vec::new();
    let mut current_section: Option<(String, Vec<String>)> = None;

    let re = Regex::new(";.*").unwrap();
    for line in reader.lines() {
        let line = line?;
        let line = re.replace_all(&line, "");
        let line = line.trim();

        // 跳过空行和注释
        if line.is_empty() || line.starts_with(';') {
            continue;
        }

        // 检查是否是新的部分
        if line.starts_with('[') && line.ends_with(']') {
            // 如果有当前部分，先保存
            if let Some((name, content)) = current_section.take() {
                sections.push(parse_section(&name, content));
            }

            // 开始新的部分
            let section_name = line[1..line.len()-1].trim();
            current_section = Some((section_name.to_string(), Vec::new()));
        } else if let Some((_, ref mut content)) = current_section {
            // 添加到当前部分
            content.push(line.to_string());
        }
    }

    // 保存最后一个部分
    if let Some((name, content)) = current_section {
        sections.push(parse_section(&name, content));
    }

    Ok(sections)
}

// 解析拓扑部分
fn parse_section(name: &str, content: Vec<String>) -> TopologySection {
    // 标准化部分名称（小写并去除空白）
    let normalized_name = name.to_lowercase().trim().to_string();
    
    match normalized_name.as_str() {
        "atomtypes" => {
            let mut atomtypes: HashSet<Atomtype> = HashSet::new();
            for item in &content {
                atomtypes.insert(parse_atomtype_line(item).unwrap());
            }
            TopologySection::Atomtypes(atomtypes)
        },
        "atoms" => {
            let atoms: Vec<Atom> = content
                .iter()
                .filter_map(|line| parse_atom_line(line))
                .collect();
            TopologySection::Atoms(atoms)
        },
        "bonds" => {
            let bonds: Vec<Bond> = content
                .iter()
                .filter_map(|line| parse_bond_line(line))
                .collect();
            TopologySection::Bonds(bonds)
        },
        "angles" => {
            let angles: Vec<Angle> = content
                .iter()
                .filter_map(|line| parse_angle_line(line))
                .collect();
            TopologySection::Angles(angles)
        },
        "dihedrals" => {
            let dihedrals: Vec<Dihedral> = content
                .iter()
                .filter_map(|line| parse_dihedral_line(line))
                .collect();
            TopologySection::Dihedrals(dihedrals)
        },
        "cmap" => {
            let cmap: Vec<CMAP> = content
                .iter()
                .filter_map(|line| parse_cmap_line(line))
                .collect();
            TopologySection::CMAP(cmap)
        },
        _ => TopologySection::Others(name.to_string(), content),
    }
}

// 解析原子类型行
fn parse_atomtype_line(line: &str) -> Option<Atomtype> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 7 {
        return None;
    }

    Some(Atomtype {
        name: parts[0].to_string(),
        id: parts[1].parse().ok()?,
        mass: parts[2].parse().ok()?,
        charge: parts[3].parse().ok()?,
        ptype: parts[4].chars().next()?,
        sigma: parts[5].parse().ok()?,
        epsilon: parts[6].parse().ok()?
    })
}

// 解析原子行 - 适配输入文件的列顺序
fn parse_atom_line(line: &str) -> Option<Atom> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 8 {
        return None;
    }

    Some(Atom {
        index: parts[0].parse().ok()?,
        type_name: parts[1].to_string(),
        res_index: parts[2].parse().ok()?,
        res_name: parts[3].to_string(),
        name: parts[4].to_string(),
        charge: parts[6].parse().ok()?,
        mass: parts[7].parse().ok()?,
        other_fields: parts[8..].iter().map(|s| s.to_string()).collect(),
    })
}

// 解析键行
fn parse_bond_line(line: &str) -> Option<Bond> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 3 {
        return None;
    }

    Some(Bond {
        atom1: parts[0].parse().ok()?,
        atom2: parts[1].parse().ok()?,
        type_: parts[2].to_string(),
        other_fields: parts[3..].iter().map(|s| s.to_string()).collect(),
    })
}

// 解析角度行
fn parse_angle_line(line: &str) -> Option<Angle> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 4 {
        return None;
    }

    Some(Angle {
        atom1: parts[0].parse().ok()?,
        atom2: parts[1].parse().ok()?,
        atom3: parts[2].parse().ok()?,
        type_: parts[3].to_string(),
        other_fields: parts[4..].iter().map(|s| s.to_string()).collect(),
    })
}

// 解析二面角行
fn parse_dihedral_line(line: &str) -> Option<Dihedral> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 5 {
        return None;
    }

    Some(Dihedral {
        atom1: parts[0].parse().ok()?,
        atom2: parts[1].parse().ok()?,
        atom3: parts[2].parse().ok()?,
        atom4: parts[3].parse().ok()?,
        type_: parts[4].to_string(),
        other_fields: parts[5..].iter().map(|s| s.to_string()).collect(),
    })
}

// 解析二面角行
fn parse_cmap_line(line: &str) -> Option<CMAP> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 5 {
        return None;
    }

    Some(CMAP {
        atom1: parts[0].parse().ok()?,
        atom2: parts[1].parse().ok()?,
        atom3: parts[2].parse().ok()?,
        atom4: parts[3].parse().ok()?,
        atom5: parts[4].parse().ok()?,
        type_: parts[5].to_string(),
    })
}

// 重排原子
pub fn reorder_atoms(sections: Vec<TopologySection>, order: &[u32]) -> Vec<TopologySection> {
    // 创建原子索引映射：旧索引 -> 新索引
    let mut index_map: Vec<u32> = vec![0; order.len() + 1];  // +1 因为原子索引从1开始
    for (new_idx, &old_idx) in order.iter().enumerate() {
        index_map[old_idx as usize] = (new_idx + 1) as u32;  // 新索引也从1开始
    }

    let mut reordered = Vec::new();

    for section in sections {
        match section {
            TopologySection::Atoms(atoms) => {
                // 重排原子 - 确保所有原子都被保留并正确写入
                let mut reordered_atoms = Vec::with_capacity(atoms.len());
                let mut remaining_atoms: Vec<Atom> = atoms.clone();

                // 首先按照指定顺序添加原子
                for &index in order {
                    if let Some(pos) = remaining_atoms.iter().position(|atom| atom.index == index) {
                        let mut atom = remaining_atoms.remove(pos);
                        atom.index = (reordered_atoms.len() + 1) as u32;
                        reordered_atoms.push(atom);
                    } else {
                        // 原子未找到，跳过
                    }
                }

                // 然后添加任何未在order中指定的原子
                for mut atom in remaining_atoms {
                    atom.index = (reordered_atoms.len() + 1) as u32;
                    reordered_atoms.push(atom);
                }

                reordered.push(TopologySection::Atoms(reordered_atoms));
            },
            TopologySection::Bonds(bonds) => {
                // 更新键中的原子索引
                let new_bonds: Vec<Bond> = bonds
                    .into_iter()
                    .map(|bond| {
                        let mut new_bond = bond.clone();
                        new_bond.atom1 = index_map[bond.atom1 as usize];
                        new_bond.atom2 = index_map[bond.atom2 as usize];
                        new_bond
                    })
                    .collect();

                reordered.push(TopologySection::Bonds(new_bonds));
            },
            TopologySection::Angles(angles) => {
                // 更新角度中的原子索引
                let new_angles: Vec<Angle> = angles
                    .into_iter()
                    .map(|angle| {
                        let mut new_angle = angle.clone();
                        new_angle.atom1 = index_map[angle.atom1 as usize];
                        new_angle.atom2 = index_map[angle.atom2 as usize];
                        new_angle.atom3 = index_map[angle.atom3 as usize];
                        new_angle
                    })
                    .collect();

                reordered.push(TopologySection::Angles(new_angles));
            },
            TopologySection::Dihedrals(dihedrals) => {
                // 更新二面角中的原子索引
                let new_dihedrals: Vec<Dihedral> = dihedrals
                    .into_iter()
                    .map(|dihedral| {
                        let mut new_dihedral = dihedral.clone();
                        new_dihedral.atom1 = index_map[dihedral.atom1 as usize];
                        new_dihedral.atom2 = index_map[dihedral.atom2 as usize];
                        new_dihedral.atom3 = index_map[dihedral.atom3 as usize];
                        new_dihedral.atom4 = index_map[dihedral.atom4 as usize];
                        new_dihedral
                    })
                    .collect();

                reordered.push(TopologySection::Dihedrals(new_dihedrals));
            },
            TopologySection::CMAP(cmaps) => {
                // 更新二面角中的原子索引
                let new_cmaps: Vec<CMAP> = cmaps
                    .into_iter()
                    .map(|cmap| {
                        let mut new_cmap = cmap.clone();
                        new_cmap.atom1 = index_map[cmap.atom1 as usize];
                        new_cmap.atom2 = index_map[cmap.atom2 as usize];
                        new_cmap.atom3 = index_map[cmap.atom3 as usize];
                        new_cmap.atom4 = index_map[cmap.atom4 as usize];
                        new_cmap.atom5 = index_map[cmap.atom5 as usize];
                        new_cmap
                    })
                    .collect();

                reordered.push(TopologySection::CMAP(new_cmaps));
            },
            other => reordered.push(other),
        }
    }

    reordered
}

