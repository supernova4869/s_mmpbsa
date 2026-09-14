//! `gmx make_ndx`, mirroring the interactive editor of
//! `src/gromacs/tools/make_ndx.cpp`.

use std::fmt::Write as _;

use crate::cmd::Args;
use crate::frame::Atoms;
use crate::index::{self, IndexGroup};
use crate::tpr;
use crate::trx;

const MAXNAMES: usize = 1024;
const NOTSET: i32 = -92637;

struct Editor {
    atoms: Option<Atoms>,
    x: Vec<[f32; 3]>,
    groups: Vec<IndexGroup>,
    natoms: usize,
    case_sensitive: bool,
    verbose: bool,
}

impl Editor {
    fn have_atoms(&self, s: &str) -> bool {
        if self.atoms.is_none() {
            println!("Can not process '{s}' without atom info, use option -f");
            false
        } else {
            true
        }
    }

    // --- parser helpers -------------------------------------------------

    /// `parse_int`: returns (value, new index) or None.
    fn parse_int(s: &[u8], i: usize) -> Option<(i32, usize)> {
        let mut j = i;
        while j < s.len() && s[j] == b' ' {
            j += 1;
        }
        let start = j;
        if j < s.len() && s[j].is_ascii_digit() {
            let mut nr: i32 = 0;
            while j < s.len() && s[j].is_ascii_digit() {
                nr = nr * 10 + (s[j] - b'0') as i32;
                j += 1;
            }
            if j < s.len() && s[j].is_ascii_alphabetic() {
                j += 1;
            }
            if j >= s.len() || !s[j].is_ascii_alphanumeric() {
                Some((nr, j))
            } else {
                let _ = start;
                None
            }
        } else {
            None
        }
    }

    /// `parse_int_char`: returns (value, insertion code, new index).
    fn parse_int_char(s: &[u8], i: usize) -> Option<(i32, u8, usize)> {
        let mut j = i;
        while j < s.len() && s[j] == b' ' {
            j += 1;
        }
        if j < s.len() && s[j].is_ascii_digit() {
            let mut nr: i32 = 0;
            while j < s.len() && s[j].is_ascii_digit() {
                nr = nr * 10 + (s[j] - b'0') as i32;
                j += 1;
            }
            let mut c = b' ';
            if j < s.len() && s[j].is_ascii_alphabetic() {
                c = s[j];
                j += 1;
            }
            if j >= s.len() || !s[j].is_ascii_alphanumeric() {
                Some((nr, c, j))
            } else {
                None
            }
        } else {
            None
        }
    }

    fn is_name_char(c: u8) -> bool {
        c != 0 && c != b' ' && c != b'!' && c != b'&' && c != b'|'
    }

    /// `parse_names`
    fn parse_names(&self, s: &[u8], i: usize) -> Option<(Vec<String>, usize)> {
        let mut j = i;
        let mut names = Vec::new();
        while j < s.len() && (Self::is_name_char(s[j]) || s[j] == b' ') {
            if Self::is_name_char(s[j]) {
                if names.len() >= MAXNAMES {
                    return None;
                }
                let mut name = String::new();
                while j < s.len() && Self::is_name_char(s[j]) {
                    name.push(s[j] as char);
                    j += 1;
                }
                if !self.case_sensitive {
                    name = name.to_uppercase();
                }
                names.push(name);
            } else {
                j += 1;
            }
        }
        if names.is_empty() {
            None
        } else {
            Some((names, j))
        }
    }

    /// `parse_string`: a quoted group name.
    fn parse_string(&self, s: &[u8], i: usize) -> Option<(i32, usize)> {
        let mut j = i;
        while j < s.len() && s[j] == b' ' {
            j += 1;
        }
        if j < s.len() && s[j] == b'"' {
            j += 1;
            let mut name = String::new();
            while j < s.len() && s[j] != b'"' {
                name.push(s[j] as char);
                j += 1;
            }
            if j < s.len() {
                j += 1;
            }
            return Some((index::find_group(&name, &self.groups), j));
        }
        None
    }

    fn comp_name(&self, name: &str, search: &str) -> bool {
        let mut matches = true;
        let nb = name.as_bytes();
        let sb = search.as_bytes();
        let mut i = 0usize;
        while i < nb.len() && i < sb.len() && matches {
            if sb[i] == b'?' {
                i += 1;
                continue;
            } else if sb[i] == b'*' {
                if i + 1 < sb.len() {
                    println!("WARNING: Currently '*' is only supported at the end of an expression");
                }
                return i + 1 >= sb.len();
            }
            matches = if self.case_sensitive {
                nb[i] == sb[i]
            } else {
                nb[i].to_ascii_uppercase() == sb[i].to_ascii_uppercase()
            };
            i += 1;
        }
        matches
            && (nb.len() == i || (i == nb.len() && sb.len() == i))
            && (i == sb.len() || (sb.len() == i + 1 && sb.get(i) == Some(&b'*')))
    }

    // --- selectors -------------------------------------------------------

    fn select_atomnumbers(
        &self,
        s: &[u8],
        n1: i32,
        i: usize,
    ) -> (Vec<usize>, String, usize) {
        let atoms = self.atoms.as_ref().unwrap();
        let mut index = Vec::new();
        let mut j = i;
        while j < s.len() && s[j] == b' ' {
            j += 1;
        }
        let mut gname = String::new();
        if j < s.len() && s[j] == b'-' {
            j += 1;
            let (up, nj) = Self::parse_int(s, j).unwrap_or((0, j));
            j = nj;
            if n1 < 1 || n1 as usize > atoms.nr() || up < 1 || up as usize > atoms.nr() {
                println!("Invalid atom range");
            } else {
                for k in (n1 - 1)..=(up - 1) {
                    index.push(k as usize);
                }
                println!(
                    "Found {} atom{} in range {}-{}",
                    index.len(),
                    if index.len() == 1 { "" } else { "s" },
                    n1,
                    up
                );
                gname = if n1 == up {
                    format!("a_{n1}")
                } else {
                    format!("a_{n1}-{up}")
                };
            }
        } else {
            let mut cur = n1;
            gname.push('a');
            loop {
                if cur - 1 >= 0 && (cur - 1) < atoms.nr() as i32 {
                    index.push((cur - 1) as usize);
                    let _ = write!(gname, "_{cur}");
                } else {
                    println!("Invalid atom number {cur}");
                    index.clear();
                }
                match Self::parse_int(s, j) {
                    Some((v, nj)) => {
                        cur = v;
                        j = nj;
                    }
                    None => break,
                }
                if index.is_empty() {
                    break;
                }
            }
        }
        (index, gname, j)
    }

    fn select_residuenumbers(
        &self,
        s: &[u8],
        n1: i32,
        c: u8,
        i: usize,
        by_index: bool,
    ) -> (Vec<usize>, String, usize) {
        let atoms = self.atoms.as_ref().unwrap();
        let mut index = Vec::new();
        let mut j = i;
        while j < s.len() && s[j] == b' ' {
            j += 1;
        }
        let res_matches = |i: usize, j: i32, c: u8| -> bool {
            let ri = &atoms.resinfo[atoms.atom[i].resind as usize];
            if by_index {
                atoms.atom[i].resind + 1 == j && (c == b' ' || ri.ic == c)
            } else {
                ri.nr == j && (c == b' ' || ri.ic == c)
            }
        };
        if j < s.len() && s[j] == b'-' {
            if c != b' ' {
                println!("Error: residue insertion codes can not be used with residue range selection");
                return (index, String::new(), j);
            }
            j += 1;
            let (up, nj) = Self::parse_int(s, j).unwrap_or((0, j));
            j = nj;
            for i in 0..atoms.nr() {
                for k in n1..=up {
                    if res_matches(i, k, c) {
                        index.push(i);
                    }
                }
            }
            println!(
                "Found {} atom{} with res.nr. in range {}-{}",
                index.len(),
                if index.len() == 1 { "" } else { "s" },
                n1,
                up
            );
            let gname = if n1 == up {
                format!("r_{n1}")
            } else {
                format!("r_{n1}-{up}")
            };
            (index, gname, j)
        } else {
            let mut cur = n1;
            let mut cc = c;
            let mut gname = String::from("r");
            loop {
                for i in 0..atoms.nr() {
                    if res_matches(i, cur, cc) {
                        index.push(i);
                    }
                }
                let _ = write!(gname, "_{cur}");
                match Self::parse_int_char(s, j) {
                    Some((v, nc, nj)) => {
                        cur = v;
                        cc = nc;
                        j = nj;
                    }
                    None => break,
                }
            }
            (index, gname, j)
        }
    }

    fn select_by_name(
        &self,
        names: &[String],
        kind: NameKind,
    ) -> Vec<usize> {
        let atoms = self.atoms.as_ref().unwrap();
        let mut index = Vec::new();
        for i in 0..atoms.nr() {
            let name = match kind {
                NameKind::Atom => atoms.atom[i].name.clone(),
                NameKind::Type => atoms.atom[i].atom_type.clone(),
                NameKind::Residue => {
                    atoms.resinfo[atoms.atom[i].resind as usize].name.clone()
                }
                NameKind::Chain => {
                    let c = atoms.resinfo[atoms.atom[i].resind as usize].chainid;
                    (c as char).to_string()
                }
            };
            if names.iter().any(|n| self.comp_name(&name, n)) {
                index.push(i);
            }
        }
        index
    }

    fn list_residues(&self) {
        let atoms = self.atoms.as_ref().unwrap();
        if atoms.nr() == 0 {
            return;
        }
        let mut start = atoms.atom[0].resind;
        let mut prev = start;
        for i in 0..atoms.nr() {
            let resind = atoms.atom[i].resind;
            if resind != prev || i == atoms.nr() - 1 {
                let diff = atoms.resinfo[resind as usize].name
                    != atoms.resinfo[start as usize].name;
                if diff || i == atoms.nr() - 1 {
                    let end = if diff { prev } else { resind };
                    if end < start + 3 {
                        let mut line = String::new();
                        for j in start..=end {
                            let _ = write!(
                                line,
                                "{:4} {:<5}",
                                j + 1,
                                atoms.resinfo[j as usize].name
                            );
                        }
                        println!("{line}");
                    } else {
                        println!(
                            " {:4} - {:4} {:<5}  ",
                            start + 1,
                            end + 1,
                            atoms.resinfo[start as usize].name
                        );
                    }
                    start = resind;
                }
            }
            prev = resind;
        }
    }

    fn split_group(&mut self, sel: usize, by_atom: bool) {
        let atoms = match &self.atoms {
            Some(a) => a.clone(),
            None => return,
        };
        let name_to_split = self.groups[sel].name.clone();
        println!(
            "Splitting group {} '{}' into {}",
            sel,
            name_to_split,
            if by_atom { "atoms" } else { "residues" }
        );
        let group = self.groups[sel].particle_indices.clone();
        let mut prev_atom: i32 = -1;
        let mut new_groups: Vec<IndexGroup> = Vec::new();
        for a in group {
            let resind = atoms.atom[a].resind;
            let name = atoms.resinfo[resind as usize].name.clone();
            if by_atom
                || prev_atom == -1
                || atoms.atom[prev_atom as usize].resind != resind
            {
                let gname = if by_atom {
                    format!("{}_{}_{}", name_to_split, atoms.atom[a].name, a + 1)
                } else {
                    format!(
                        "{}_{}_{}",
                        name_to_split, name, atoms.resinfo[resind as usize].nr
                    )
                };
                new_groups.push(IndexGroup::new(&gname, Vec::new()));
            }
            new_groups.last_mut().unwrap().particle_indices.push(a);
            prev_atom = a as i32;
        }
        self.groups.extend(new_groups);
    }

    #[allow(unused_assignments)]
    fn split_chain(&mut self, sel: usize) {
        let atoms = match &self.atoms {
            Some(a) => a.clone(),
            None => return,
        };
        let natoms = atoms.nr();
        let mut chains: Vec<(usize, usize)> = Vec::new();
        let mut ca_start = 0usize;
        while ca_start < natoms {
            while ca_start < natoms && atoms.atom[ca_start].name != "CA" {
                ca_start += 1;
            }
            if ca_start >= natoms {
                break;
            }
            let mut start = ca_start;
            while start > 0 && atoms.atom[start - 1].resind == atoms.atom[ca_start].resind {
                start -= 1;
            }
            let mut i = ca_start;
            let ca_end;
            loop {
                let mut e = i;
                loop {
                    i += 1;
                    if i >= natoms || atoms.atom[i].name == "CA" {
                        break;
                    }
                }
                if i >= natoms {
                    ca_end = e;
                    break;
                }
                let dx = [
                    self.x[e][0] - self.x[i][0],
                    self.x[e][1] - self.x[i][1],
                    self.x[e][2] - self.x[i][2],
                ];
                let norm = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                if norm >= 0.45 {
                    ca_end = e;
                    break;
                }
                e = i;
            }
            let mut end = ca_end;
            while end + 1 < natoms
                && atoms.atom[end + 1].resind == atoms.atom[ca_end].resind
            {
                end += 1;
            }
            chains.push((start, end));
            ca_start = end + 1;
        }
        if chains.len() == 1 {
            println!("Found 1 chain, will not split");
        } else {
            println!("Found {} chains", chains.len());
        }
        for (j, (s, e)) in chains.iter().enumerate() {
            println!("{}:{:6} atoms ({} to {})", j + 1, e - s + 1, s + 1, e + 1);
        }
        if chains.len() > 1 {
            for (j, (s, e)) in chains.iter().enumerate() {
                let particles: Vec<usize> = self.groups[sel]
                    .particle_indices
                    .iter()
                    .copied()
                    .filter(|a| a >= s && a <= e)
                    .collect();
                if !particles.is_empty() {
                    let name = format!("{}_chain{}", self.groups[sel].name, j + 1);
                    self.groups.push(IndexGroup::new(&name, particles));
                }
            }
        }
    }

    /// `parse_entry`: returns (group, name, index list) and the new position.
    #[allow(unused_assignments)]
    fn parse_entry(&self, s: &[u8], i: usize) -> Option<(Vec<usize>, String, usize, bool)> {
        let mut j = i;
        while j < s.len() && s[j] == b' ' {
            j += 1;
        }
        let mut compl = false;
        if j < s.len() && s[j] == b'!' {
            compl = true;
            j += 1;
            while j < s.len() && s[j] == b' ' {
                j += 1;
            }
        }
        let ostring = String::from_utf8_lossy(&s[j..]).to_string();
        let mut index: Vec<usize> = Vec::new();
        let mut gname = String::new();

        if let Some((v, nj)) = Self::parse_int(s, j) {
            if v >= 0 && (v as usize) < self.groups.len() {
                index = self.groups[v as usize].particle_indices.clone();
                gname = self.groups[v as usize].name.clone();
                println!("Copied index group {v} '{gname}'");
                j = nj;
            } else {
                println!("Group {v} does not exist");
                return None;
            }
        } else if let Some((v, nj)) = self.parse_string(s, j) {
            if v >= 0 && (v as usize) < self.groups.len() {
                index = self.groups[v as usize].particle_indices.clone();
                gname = self.groups[v as usize].name.clone();
                println!("Copied index group {v} '{gname}'");
                j = nj;
            } else {
                return None;
            }
        } else if j < s.len() && s[j] == b'a' {
            j += 1;
            if !self.have_atoms(&ostring) {
                return None;
            }
            if let Some((v, nj)) = Self::parse_int(s, j) {
                let (idx, g, njj) = self.select_atomnumbers(s, v, nj);
                index = idx;
                gname = g;
                j = njj;
            } else if let Some((names, nj)) = self.parse_names(s, j) {
                index = self.select_by_name(&names, NameKind::Atom);
                println!(
                    "Found {} atoms with name{} {}",
                    index.len(),
                    if names.len() == 1 { "" } else { "s" },
                    names.join(" ")
                );
                gname = names.join("_");
                j = nj;
            } else {
                return None;
            }
        } else if j < s.len() && s[j] == b't' {
            j += 1;
            if !self.have_atoms(&ostring) {
                return None;
            }
            let (names, nj) = self.parse_names(s, j)?;
            index = self.select_by_name(&names, NameKind::Type);
            println!(
                "Found {} atoms with type{} {}",
                index.len(),
                if names.len() == 1 { "" } else { "s" },
                names.join(" ")
            );
            gname = names.join("_");
            j = nj;
        } else if s[j..].starts_with(b"res") {
            j += 3;
            if !self.have_atoms(&ostring) {
                return None;
            }
            let (v, nj) = Self::parse_int(s, j)?;
            j = nj;
            if v < 0 || v as usize >= self.groups.len() {
                return None;
            }
            let atoms = self.atoms.as_ref().unwrap();
            let selected = &self.groups[v as usize];
            for resnr in &selected.particle_indices {
                if *resnr >= atoms.nres() {
                    println!(
                        "Index {} contains number>nres ({}>{})",
                        selected.name,
                        resnr + 1,
                        atoms.nres()
                    );
                    return None;
                }
            }
            for k in 0..atoms.nr() {
                let resnr = atoms.resinfo[atoms.atom[k].resind as usize].nr;
                if selected.particle_indices.iter().any(|r| *r as i32 + 1 == resnr) {
                    index.push(k);
                }
            }
            println!(
                "Found {} atoms in {} residues from group {}",
                index.len(),
                selected.particle_indices.len(),
                selected.name
            );
            gname = format!("atom_{}", selected.name);
        } else if s[j..].starts_with(b"ri") {
            j += 2;
            if !self.have_atoms(&ostring) {
                return None;
            }
            let (v, c, nj) = Self::parse_int_char(s, j)?;
            let (idx, g, njj) = self.select_residuenumbers(s, v, c, nj, true);
            index = idx;
            gname = g;
            j = njj;
        } else if j < s.len() && s[j] == b'r' {
            j += 1;
            if !self.have_atoms(&ostring) {
                return None;
            }
            if let Some((v, c, nj)) = Self::parse_int_char(s, j) {
                let (idx, g, njj) = self.select_residuenumbers(s, v, c, nj, false);
                index = idx;
                gname = g;
                j = njj;
            } else if let Some((names, nj)) = self.parse_names(s, j) {
                index = self.select_by_name(&names, NameKind::Residue);
                println!(
                    "Found {} atoms with residue name{} {}",
                    index.len(),
                    if names.len() == 1 { "" } else { "s" },
                    names.join(" ")
                );
                gname = names.join("_");
                j = nj;
            } else {
                return None;
            }
        } else if s[j..].starts_with(b"chain") {
            j += 5;
            if !self.have_atoms(&ostring) {
                return None;
            }
            let (names, nj) = self.parse_names(s, j)?;
            index = self.select_by_name(&names, NameKind::Chain);
            println!(
                "Found {} atoms with chain identifier{} {}",
                index.len(),
                if names.len() == 1 { "" } else { "s" },
                names.join(" ")
            );
            gname = format!("ch{}", names.join(""));
            j = nj;
        } else {
            return None;
        }

        if compl {
            let mut keep = vec![false; self.natoms];
            for &k in &index {
                if k < keep.len() {
                    keep[k] = true;
                }
            }
            index = (0..self.natoms).filter(|k| !keep[*k]).collect();
            gname = format!("!{gname}");
            println!("Complemented group: {} atoms", index.len());
        }

        Some((index, gname, j, true))
    }

    fn remove_group(&mut self, first: i32, last: i32) {
        for j in 0..=(last - first) {
            if first < 0 || first as usize >= self.groups.len() {
                println!("Group {} does not exist", first + j);
            } else {
                println!(
                    "Removed group {} '{}'",
                    first + j,
                    self.groups[first as usize].name
                );
                self.groups.remove(first as usize);
            }
        }
    }

    fn print_groups(&self, all: bool, only: Option<usize>) {
        println!();
        let (i0, i1) = if all {
            (0, self.groups.len())
        } else {
            let n = only.unwrap_or(0);
            (n, (n + 1).min(self.groups.len()))
        };
        for i in i0..i1 {
            println!(
                "{:3} {:<20}: {:5} atoms",
                i,
                self.groups[i].name,
                self.groups[i].particle_indices.len()
            );
        }
    }

    fn run_editor(&mut self, script: &str) {
        let lines: Vec<&str> = if script.is_empty() {
            Vec::new()
        } else {
            script.lines().collect()
        };
        let mut li = 0usize;
        let mut print_once = !crate::cmd::is_scripted();
        let mut newgroup: i32 = NOTSET;
        loop {
            if self.verbose || print_once || newgroup != NOTSET {
                self.print_groups(self.verbose || print_once, Some(newgroup.max(0) as usize));
                newgroup = NOTSET;
            }
            if self.verbose || print_once {
                self.print_menu();
                print_once = false;
            }
            if !crate::cmd::is_scripted() {
                println!("\n> ");
            }
            let line = if li < lines.len() {
                let l = lines[li];
                li += 1;
                l.to_string()
            } else {
                match crate::cmd::read_input_line() {
                    Some(buf) => buf,
                    None => break,
                }
            };
            let trimmed = line.trim_end_matches(['\n', '\r']);
            let mut string = trimmed;
            while string.starts_with(' ') {
                string = &string[1..];
            }
            let bytes = string.as_bytes();
            let is_quit = string.starts_with('q');

            if string.is_empty() {
                print_once = true;
            } else if string.starts_with('h') {
                self.print_help();
            } else if string.starts_with("del") {
                let s = &bytes[3..];
                if let Some((nr, i)) = Self::parse_int(s, 0) {
                    let mut j = i;
                    while j < s.len() && s[j] == b' ' {
                        j += 1;
                    }
                    let nr2 = if j < s.len() && s[j] == b'-' {
                        Self::parse_int(s, j + 1).map(|(v, _)| v).unwrap_or(nr)
                    } else {
                        nr
                    };
                    self.remove_group(nr, nr2);
                }
            } else if string.starts_with("keep") {
                let s = &bytes[4..];
                if let Some((nr, _)) = Self::parse_int(s, 0) {
                    self.remove_group(nr + 1, self.groups.len() as i32 - 1);
                    self.remove_group(0, nr - 1);
                }
            } else if string.starts_with("name") {
                let s = &bytes[4..];
                if let Some((nr, i)) = Self::parse_int(s, 0) {
                    if nr >= 0 && (nr as usize) < self.groups.len() {
                        let rest = String::from_utf8_lossy(&s[i..]).to_string();
                        if let Some(tok) = rest.split_whitespace().next() {
                            self.groups[nr as usize].name = tok.to_string();
                        }
                    }
                }
            } else if string.starts_with("case") {
                self.case_sensitive = !self.case_sensitive;
                println!(
                    "Switched to case {}",
                    if self.case_sensitive { "sensitive" } else { "insensitive" }
                );
            } else if string.starts_with('v') {
                self.verbose = !self.verbose;
                println!("Turned verbose {}", if self.verbose { "on" } else { "off" });
            } else if string.starts_with('l') {
                if self.have_atoms(string) {
                    self.list_residues();
                }
            } else if string.starts_with("splitch") {
                if let Some((sel, _)) = Self::parse_int(&bytes[7..], 0) {
                    if sel >= 0 && (sel as usize) < self.groups.len() {
                        self.split_chain(sel as usize);
                    }
                }
            } else if string.starts_with("splitres") {
                if let Some((sel, _)) = Self::parse_int(&bytes[8..], 0) {
                    if sel >= 0 && (sel as usize) < self.groups.len() {
                        self.split_group(sel as usize, false);
                    }
                }
            } else if string.starts_with("splitat") {
                if let Some((sel, _)) = Self::parse_int(&bytes[7..], 0) {
                    if sel >= 0 && (sel as usize) < self.groups.len() {
                        self.split_group(sel as usize, true);
                    }
                }
            } else if !is_quit {
                let mut pos;
                if let Some((mut nr, mut gname, p, _)) = self.parse_entry(bytes, 0) {
                    pos = p;
                    loop {
                        while pos < bytes.len() && bytes[pos] == b' ' {
                            pos += 1;
                        }
                        let (b_and, b_or) = if pos < bytes.len() && bytes[pos] == b'&' {
                            (true, false)
                        } else if pos < bytes.len() && bytes[pos] == b'|' {
                            (false, true)
                        } else {
                            (false, false)
                        };
                        if !b_and && !b_or {
                            break;
                        }
                        pos += 1;
                        if let Some((nr2, gname2, p2, _)) = self.parse_entry(bytes, pos) {
                            pos = p2;
                            if b_or {
                                nr = or_groups(&nr, &nr2);
                                gname = format!("{gname}_{gname2}");
                            } else {
                                nr = and_groups(&nr, &nr2);
                                gname = format!("{gname}_&_{gname2}");
                            }
                        } else {
                            break;
                        }
                    }
                    while pos < bytes.len() && bytes[pos] == b' ' {
                        pos += 1;
                    }
                    if pos < bytes.len() {
                        println!(
                            "\nSyntax error: \"{}\"",
                            String::from_utf8_lossy(&bytes[pos..])
                        );
                    } else if !nr.is_empty() {
                        self.groups.push(IndexGroup::new(&gname, nr));
                    } else {
                        println!("Group is empty");
                    }
                }
            }
            if is_quit {
                break;
            }
        }
    }

    fn print_menu(&self) {
        println!();
        println!(
            " nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups"
        );
        println!(" 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues");
        println!(" 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help");
        println!(" 'r': residue              'res' nr         'chain' char");
        println!(
            " \"name\": group             'case': case {}         'q': save and quit",
            if self.case_sensitive { "sensitive  " } else { "insensitive" }
        );
        println!(" 'ri': residue index");
    }

    fn print_help(&self) {
        println!(" nr                : selects an index group by number or quoted string.");
        println!(" 'a' nr1 [nr2 ...] : selects atoms, atom numbering starts at 1.");
        println!(" 'a' nr1 - nr2     : selects atoms in the range from nr1 to nr2.");
        println!(" 'a' name1[*] [name2[*] ...] : selects atoms by name(s).");
        println!(" 't' type1[*] [type2[*] ...] : as 'a', but for type, run input file required.");
        println!(" 'r' nr1[ic1] [nr2[ic2] ...] : selects residues by number and insertion code.");
        println!(" 'r' nr1 - nr2               : selects residues in the range from nr1 to nr2.");
        println!(" 'r' name1[*] [name2[*] ...] : as 'a', but for residue names.");
        println!(" 'ri' nr1 - nr2              : selects residue indices, 1-indexed.");
        println!(" 'chain' ch1 [ch2 ...]       : selects atoms by chain identifier(s).");
        println!(" !                 : takes the complement of a group.");
        println!(" & |               : AND and OR, processed left to right.");
        println!(" 'name' nr name    : rename group nr to name.");
        println!(" 'del' nr1 [- nr2] : deletes one group or groups in a range.");
        println!(" 'keep' nr         : deletes all groups except nr.");
        println!(" 'case'            : make all name compares case (in)sensitive.");
        println!(" 'splitch' nr      : split group into chains using CA distances.");
        println!(" 'splitres' nr     : split group into residues.");
        println!(" 'splitat' nr      : split group into atoms.");
        println!(" 'res' nr          : interpret numbers in group as residue numbers");
        println!(" Enter             : list the currently defined groups and commands");
        println!(" 'l'               : list the residues.");
        println!(" 'h'               : show this help.");
        println!(" 'q'               : save and quit.");
    }
}

enum NameKind {
    Atom,
    Type,
    Residue,
    Chain,
}

fn or_groups(nr1: &[usize], nr2: &[usize]) -> Vec<usize> {
    let a: Vec<usize> = nr1.to_vec();
    let b: Vec<usize> = nr2.to_vec();
    let mut not_incr = false;
    let mut max = 0usize;
    for (i, v) in a.iter().enumerate() {
        if i > 0 && *v <= max {
            not_incr = true;
        }
        max = *v;
    }
    for (i, v) in b.iter().enumerate() {
        if i > 0 && *v <= max {
            not_incr = true;
        }
        max = *v;
    }
    if not_incr {
        println!("One of your groups is not ascending");
        return Vec::new();
    }
    let mut out = Vec::new();
    let (mut i1, mut i2) = (0usize, 0usize);
    while i1 < a.len() || i2 < b.len() {
        if i2 == b.len() || (i1 < a.len() && a[i1] < b[i2]) {
            out.push(a[i1]);
            i1 += 1;
        } else {
            if i2 < b.len() && (i1 == a.len() || a[i1] > b[i2]) {
                out.push(b[i2]);
            }
            i2 += 1;
        }
    }
    println!(
        "Merged two groups with OR: {} {} -> {}",
        nr1.len(),
        nr2.len(),
        out.len()
    );
    out
}

fn and_groups(nr1: &[usize], nr2: &[usize]) -> Vec<usize> {
    let mut out = Vec::new();
    for a in nr1 {
        for b in nr2 {
            if a == b {
                out.push(*a);
            }
        }
    }
    println!(
        "Merged two groups with AND: {} {} -> {}",
        nr1.len(),
        nr2.len(),
        out.len()
    );
    out
}

/// Reads the index groups of the input files, newest first.
fn load_groups(editor: &mut Editor, files: &[String]) {
    let mut collected: Vec<IndexGroup> = Vec::new();
    for f in files {
        match index::read_ndx(f) {
            Ok(g) => {
                let mut newv = g;
                newv.extend(collected);
                collected = newv;
            }
            Err(e) => eprintln!("{e}"),
        }
    }
    editor.groups = collected;
}

pub fn run(argv: Vec<String>) -> i32 {
    let args = Args::parse(argv);
    let f = args.get("f");
    let n_files = args.get_all("n");
    let o = args.get("o").unwrap_or_else(|| "index.ndx".to_string());
    let duplicate = args.flag("twin");
    let verbose = args.flag("verbose");
    let natoms_opt = args.get("natoms").and_then(|s| s.parse::<usize>().ok());

    if f.is_none() && n_files.is_empty() {
        eprintln!("No input files (structure or index)");
        return 1;
    }

    let mut editor = Editor {
        atoms: None,
        x: Vec::new(),
        groups: Vec::new(),
        natoms: 0,
        case_sensitive: false,
        verbose,
    };

    let mut natoms_from_structure = None;
    if let Some(path) = &f {
        eprintln!("\nReading structure file");
        let ext = crate::trx::format_from_path(path);
        match ext {
            Some(crate::trx::TrxFormat::Gro) | Some(crate::trx::TrxFormat::Pdb) => {
                match trx::read_first_frame(path) {
                    Ok(fr) => {
                        editor.atoms = fr.atoms.clone();
                        editor.x = fr.x.clone().unwrap_or_default();
                        natoms_from_structure = Some(fr.natoms);
                    }
                    Err(e) => {
                        eprintln!("{e}");
                        return 1;
                    }
                }
            }
            _ => match tpr::TprFile::read(path) {
                Ok(t) => match tpr::parse_body(&t.header, &t.body) {
                    Ok(body) => {
                        if let Some(mtop) = &body.mtop {
                            editor.atoms = Some(mtop.global_atoms());
                            editor.x = body.x.clone().unwrap_or_default();
                            natoms_from_structure = Some(mtop.natoms);
                        }
                    }
                    Err(e) => {
                        eprintln!("{e}");
                        return 1;
                    }
                },
                Err(e) => {
                    eprintln!("{e}");
                    return 1;
                }
            },
        }
    }

    println!("Going to read {} old index file(s)", n_files.len());
    if !n_files.is_empty() {
        load_groups(&mut editor, &n_files);
    } else if let Some(atoms) = editor.atoms.clone() {
        editor.groups = index::analyse(&atoms, true);
    }

    editor.natoms = match natoms_opt.or(natoms_from_structure) {
        Some(n) => n,
        None => {
            let mut max = -1i64;
            for g in &editor.groups {
                for a in &g.particle_indices {
                    max = max.max(*a as i64);
                }
            }
            let n = (max + 1) as usize;
            println!("Deducing {n} atoms in the system from indices in the index file");
            n
        }
    };

    editor.run_editor("");

    if let Err(e) = index::write_ndx(&o, &editor.groups, duplicate, editor.natoms) {
        eprintln!("{e}");
        return 1;
    }
    0
}
