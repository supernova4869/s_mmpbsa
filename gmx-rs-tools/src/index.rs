//! Index (.ndx) file handling and default group generation, mirroring
//! `gromacs/topology/index.cpp`.

use std::collections::HashMap;
use std::fmt::Write as _;

use crate::frame::Atoms;
use crate::xdr::{Result, XdrError};

/// One index group, equivalent to GROMACS' `IndexGroup`.
#[derive(Debug, Clone, Default)]
pub struct IndexGroup {
    pub name: String,
    pub particle_indices: Vec<usize>,
}

impl IndexGroup {
    pub fn new(name: &str, idx: Vec<usize>) -> IndexGroup {
        IndexGroup {
            name: name.to_string(),
            particle_indices: idx,
        }
    }
}

/// Case insensitive comparison ignoring `-` and `_` (`gmx_strcasecmp_min`).
pub fn strcasecmp_min(a: &str, b: &str) -> i32 {
    let mut ai = a.chars();
    let mut bi = b.chars();
    loop {
        let mut ch1;
        loop {
            ch1 = match ai.next() {
                Some(c) => c,
                None => '\0',
            };
            if ch1 != '-' && ch1 != '_' {
                break;
            }
        }
        let mut ch2;
        loop {
            ch2 = match bi.next() {
                Some(c) => c,
                None => '\0',
            };
            if ch2 != '-' && ch2 != '_' {
                break;
            }
        }
        let c1 = ch1.to_ascii_uppercase();
        let c2 = ch2.to_ascii_uppercase();
        if c1 != c2 {
            return c1 as i32 - c2 as i32;
        }
        if c1 == '\0' {
            return 0;
        }
    }
}

/// `gmx_strncasecmp_min`: like above but only compares `n` input characters.
pub fn strncasecmp_min(a: &str, b: &str, n: usize) -> i32 {
    let mut ai = a.chars();
    let mut bi = b.chars();
    let mut ca = 0usize;
    let mut cb = 0usize;
    loop {
        let mut ch1;
        loop {
            ch1 = match ai.next() {
                Some(c) => {
                    ca += 1;
                    c
                }
                None => '\0',
            };
            if ch1 != '-' && ch1 != '_' {
                break;
            }
        }
        let mut ch2;
        loop {
            ch2 = match bi.next() {
                Some(c) => {
                    cb += 1;
                    c
                }
                None => '\0',
            };
            if ch2 != '-' && ch2 != '_' {
                break;
            }
        }
        let c1 = ch1.to_ascii_uppercase();
        let c2 = ch2.to_ascii_uppercase();
        if c1 != c2 {
            return c1 as i32 - c2 as i32;
        }
        if c1 == '\0' || ca >= n || cb >= n {
            return 0;
        }
    }
}

/// Plain case insensitive comparison (`gmx_strcasecmp`).
pub fn strcasecmp(a: &str, b: &str) -> i32 {
    let mut ai = a.chars().map(|c| c.to_ascii_uppercase());
    let mut bi = b.chars().map(|c| c.to_ascii_uppercase());
    loop {
        let c1 = ai.next().unwrap_or('\0');
        let c2 = bi.next().unwrap_or('\0');
        if c1 != c2 {
            return c1 as i32 - c2 as i32;
        }
        if c1 == '\0' {
            return 0;
        }
    }
}

/// `minstring()`: `-` becomes `_`.
fn minstring(s: &str) -> String {
    s.chars().map(|c| if c == '-' { '_' } else { c }).collect()
}

/// `findGroupTemplated()`: whole name, then prefix, then substring match.
pub fn find_group(s: &str, groups: &[IndexGroup]) -> i32 {
    let mut aa: i32 = -1;
    let mut multiple = false;
    let n = s.len();
    for (i, g) in groups.iter().enumerate() {
        if strcasecmp_min(s, &g.name) == 0 {
            if aa != -1 {
                multiple = true;
            }
            aa = i as i32;
        }
    }
    if aa == -1 {
        for (i, g) in groups.iter().enumerate() {
            if strncasecmp_min(s, &g.name, n) == 0 {
                if aa != -1 {
                    multiple = true;
                }
                aa = i as i32;
            }
        }
    }
    if aa == -1 {
        let key = minstring(&s.to_uppercase());
        for (i, g) in groups.iter().enumerate() {
            let name = minstring(&g.name.to_uppercase());
            if name.contains(&key) {
                if aa != -1 {
                    multiple = true;
                }
                aa = i as i32;
            }
        }
    }
    if multiple {
        println!("Error: Multiple groups '{s}' selected");
        return -1;
    }
    aa
}

/// Reads an index file (`init_index()`).
pub fn read_ndx(path: &str) -> Result<Vec<IndexGroup>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| XdrError::Invalid(format!("cannot read {path}: {e}")))?;
    let mut groups: Vec<IndexGroup> = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim_start();
        if trimmed.starts_with('[') {
            let name = trimmed
                .trim_start_matches('[')
                .trim_end_matches(']')
                .trim()
                .to_string();
            groups.push(IndexGroup::new(&name, Vec::new()));
        } else if !trimmed.is_empty() {
            if groups.is_empty() {
                return Err(XdrError::Invalid(
                    "The first header of your indexfile is invalid".into(),
                ));
            }
            let g = groups.last_mut().unwrap();
            for tok in trimmed.split_whitespace() {
                if let Ok(v) = tok.parse::<i64>() {
                    g.particle_indices.push((v - 1).max(0) as usize);
                }
            }
        }
    }
    Ok(groups)
}

/// Writes an index file (`write_index()`).
pub fn write_ndx(
    path: &str,
    groups: &[IndexGroup],
    duplicate: bool,
    num_atoms: usize,
) -> Result<()> {
    let mut out = String::new();
    for g in groups {
        let _ = write!(out, "[ {} ]", g.name);
        let mut k = 0;
        for &pi in &g.particle_indices {
            let sep = if k % 15 == 0 { '\n' } else { ' ' };
            let _ = write!(out, "{sep}{:4}", pi + 1);
            k += 1;
        }
        let _ = writeln!(out);
    }
    if duplicate {
        eprintln!(
            "Duplicating the whole system with an atom offset of {num_atoms} atoms."
        );
        for g in groups {
            let _ = write!(out, "[ {}_copy ]", g.name);
            let mut k = 0;
            for &pi in &g.particle_indices {
                let sep = if k % 15 == 0 { '\n' } else { ' ' };
                let _ = write!(out, "{sep}{:4}", pi + 1 + num_atoms);
                k += 1;
            }
            let _ = writeln!(out);
        }
    }
    std::fs::write(path, out).map_err(|e| XdrError::Invalid(format!("cannot write {path}: {e}")))?;
    Ok(())
}

/// Maps residue names onto molecule categories using `residuetypes.dat`.
pub struct ResidueTypeMap {
    map: HashMap<String, String>,
    known: Vec<(String, String)>,
}

impl ResidueTypeMap {
    /// Loads the map from `residuetypes.dat`, searching the usual GROMACS
    /// locations.  Falls back to a small built-in table when the file cannot
    /// be found so that the tool stays usable stand-alone.
    pub fn load() -> ResidueTypeMap {
        let mut candidates: Vec<String> = Vec::new();
        for var in ["GMXDATA", "GMXLIB"] {
            if let Ok(v) = std::env::var(var) {
                candidates.push(format!("{v}/top/residuetypes.dat"));
                candidates.push(format!("{v}/residuetypes.dat"));
            }
        }
        candidates.push("share/top/residuetypes.dat".to_string());
        candidates.push("/usr/share/gromacs/top/residuetypes.dat".to_string());
        candidates.push("/opt/gromacs/share/gromacs/top/residuetypes.dat".to_string());

        let mut map = HashMap::new();
        for c in &candidates {
            if let Ok(content) = std::fs::read_to_string(c) {
                for line in content.lines() {
                    let line = line.trim();
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }
                    let toks: Vec<&str> = line.split_whitespace().collect();
                    if toks.len() >= 2 {
                        map.insert(toks[0].to_uppercase(), toks[1].to_string());
                    }
                }
                break;
            }
        }
        if map.is_empty() {
            for (res, cat) in BUILTIN_RESIDUE_TYPES {
                map.insert(res.to_string(), cat.to_string());
            }
        }
        let mut known: Vec<(String, String)> = map
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        known.sort();
        ResidueTypeMap { map, known }
    }

    /// `typeOfNamedDatabaseResidue()`: the category, or "Other" when unknown.
    pub fn category(&self, resname: &str) -> String {
        match self.map.get(&resname.to_uppercase()) {
            Some(c) => c.clone(),
            None if resname.is_empty() => "Other".to_string(),
            None => "Other".to_string(),
        }
    }

    pub fn entries(&self) -> &[(String, String)] {
        &self.known
    }
}

const BUILTIN_RESIDUE_TYPES: &[(&str, &str)] = &[
    ("ALA", "Protein"),
    ("ARG", "Protein"),
    ("ASN", "Protein"),
    ("ASP", "Protein"),
    ("CYS", "Protein"),
    ("GLN", "Protein"),
    ("GLU", "Protein"),
    ("GLY", "Protein"),
    ("HIS", "Protein"),
    ("HID", "Protein"),
    ("HIE", "Protein"),
    ("HIP", "Protein"),
    ("HISE", "Protein"),
    ("HISD", "Protein"),
    ("ILE", "Protein"),
    ("LEU", "Protein"),
    ("LYS", "Protein"),
    ("MET", "Protein"),
    ("PHE", "Protein"),
    ("PRO", "Protein"),
    ("SER", "Protein"),
    ("THR", "Protein"),
    ("TRP", "Protein"),
    ("TYR", "Protein"),
    ("VAL", "Protein"),
    ("SOL", "Water"),
    ("WAT", "Water"),
    ("HOH", "Water"),
    ("TIP3", "Water"),
    ("TIP4", "Water"),
    ("TIP5", "Water"),
    ("SPC", "Water"),
    ("NA", "Ion"),
    ("CL", "Ion"),
    ("K", "Ion"),
    ("MG", "Ion"),
    ("CA", "Ion"),
    ("ZN", "Ion"),
    ("FE", "Ion"),
    ("F", "Ion"),
    ("BR", "Ion"),
    ("IOD", "Ion"),
];

fn mk_aid(atoms: &Atoms, restype: &[String], typestring: &str, b_match: bool) -> Vec<usize> {
    let mut a = Vec::new();
    for i in 0..atoms.nr() {
        let resind = atoms.atom[i].resind as usize;
        let mut res = strcasecmp(&restype[resind], typestring) == 0;
        if !b_match {
            res = !res;
        }
        if res {
            a.push(i);
        }
    }
    a
}

struct ProtGroupDef {
    atomnames: &'static [&'static str],
    group_name: &'static str,
    take_complement: bool,
    wholename: i32,
    compareto: i32,
}

const PROTEIN_GROUPS: &[ProtGroupDef] = &[
    ProtGroupDef { atomnames: &[], group_name: "Protein", take_complement: true, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["H", "HN"], group_name: "Protein-H", take_complement: true, wholename: 0, compareto: -1 },
    ProtGroupDef { atomnames: &["CA"], group_name: "C-alpha", take_complement: false, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["N", "CA", "C"], group_name: "Backbone", take_complement: false, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT"], group_name: "MainChain", take_complement: false, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["N", "CA", "CB", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT"], group_name: "MainChain+Cb", take_complement: false, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT", "H1", "H2", "H3", "H", "HN"], group_name: "MainChain+H", take_complement: false, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT", "H1", "H2", "H3", "H", "HN"], group_name: "SideChain", take_complement: true, wholename: -1, compareto: -1 },
    ProtGroupDef { atomnames: &["N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT", "H1", "H2", "H3", "H", "HN"], group_name: "SideChain-H", take_complement: true, wholename: 11, compareto: -1 },
    ProtGroupDef { atomnames: &["MN1", "MN2", "MCB1", "MCB2", "MCG1", "MCG2", "MCD1", "MCD2", "MCE1", "MCE2", "MNZ1", "MNZ2"], group_name: "Prot-Masses", take_complement: true, wholename: -1, compareto: 0 },
];

fn group_matches_atomname(def: &ProtGroupDef, atomname: &str) -> bool {
    // Note: an empty atom name list means "no match"; groups that should
    // contain everything use `take_complement` (this mirrors the C loop that
    // simply does not execute for such entries).
    let atnm: &str = {
        let s = atomname.trim_start_matches(|c: char| c.is_ascii_digit());
        s
    };
    for (j, name) in def.atomnames.iter().enumerate() {
        if def.wholename == -1 || (j as i32) < def.wholename {
            if strcasecmp(name, atnm) == 0 {
                return true;
            }
        } else if atnm.len() >= name.len()
            && strcasecmp(name, &atnm[..name.len()]) == 0
        {
            return true;
        }
    }
    false
}

fn analyse_prot(restype: &[String], atoms: &Atoms, groups: &mut Vec<IndexGroup>) {
    let is_protein = |resind: usize| strcasecmp(&restype[resind], "Protein") == 0;
    for def in PROTEIN_GROUPS {
        let mut aid: Vec<usize> = Vec::new();
        for n in 0..atoms.nr() {
            let resind = atoms.atom[n].resind as usize;
            if is_protein(resind) {
                let matches = group_matches_atomname(def, &atoms.atom[n].name);
                if def.take_complement != matches {
                    aid.push(n);
                }
            }
        }
        let skip = if def.compareto == -1 {
            false
        } else {
            let other = &groups[groups.len() - 1 - def.compareto as usize].particle_indices;
            other == &aid
        };
        if !skip {
            groups.push(IndexGroup::new(def.group_name, aid));
        }
    }
}

fn analyse_other(restype: &[String], atoms: &Atoms, groups: &mut Vec<IndexGroup>) {
    let special = |r: &str| {
        strcasecmp(r, "Protein") == 0
            || strcasecmp(r, "DNA") == 0
            || strcasecmp(r, "RNA") == 0
            || strcasecmp(r, "Water") == 0
    };
    let mut restp: Vec<String> = Vec::new();
    for k in 0..atoms.nr() {
        let resind = atoms.atom[k].resind as usize;
        if !special(&restype[resind]) {
            let rname = atoms.resinfo[resind].name.clone();
            if !restp.contains(&rname) {
                restp.push(rname);
            }
        }
    }
    for rname in restp {
        let mut aid = Vec::new();
        for j in 0..atoms.nr() {
            let resind = atoms.atom[j].resind as usize;
            if atoms.resinfo[resind].name == rname {
                aid.push(j);
            }
        }
        groups.push(IndexGroup::new(&rname, aid));
    }
}

/// `analyse()`: builds the default index groups for a topology.
pub fn analyse(atoms: &Atoms, verbose: bool) -> Vec<IndexGroup> {
    let residuetypes = ResidueTypeMap::load();
    let mut groups: Vec<IndexGroup> = Vec::new();
    let aid_all: Vec<usize> = (0..atoms.nr()).collect();
    groups.push(IndexGroup::new("System", aid_all));

    let restype: Vec<String> = (0..atoms.nres())
        .map(|i| residuetypes.category(&atoms.resinfo[i].name))
        .collect();

    if verbose {
        println!("Analysing residue names:");
        let mut counts: Vec<(String, usize)> = Vec::new();
        for rt in &restype {
            match counts.iter_mut().find(|(c, _)| c == rt) {
                Some((_, n)) => *n += 1,
                None => counts.push((rt.clone(), 1)),
            }
        }
        for (c, n) in counts {
            println!("There are: {n:5} {c:>10} residues");
        }
    }

    let mut categories: Vec<String> = Vec::new();
    for rt in &restype {
        if !categories.contains(rt) {
            categories.push(rt.clone());
        }
    }

    let mut have_analysed_other = false;
    for cat in &categories {
        let aid_category = mk_aid(atoms, &restype, cat, true);
        if strcasecmp(cat, "Protein") == 0 && !aid_category.is_empty() {
            if verbose {
                println!("Analysing Protein...");
            }
            analyse_prot(&restype, atoms, &mut groups);
            let aid_non_protein = mk_aid(atoms, &restype, "Protein", false);
            if !aid_non_protein.is_empty() && aid_non_protein.len() < atoms.nr() {
                groups.push(IndexGroup::new("non-Protein", aid_non_protein));
            }
        } else if strcasecmp(cat, "Water") == 0 && !aid_category.is_empty() {
            groups.push(IndexGroup::new(cat, aid_category.clone()));
            groups.push(IndexGroup::new("SOL", aid_category));
            let aid_solvent = mk_aid(atoms, &restype, "Water", false);
            if !aid_solvent.is_empty() && aid_solvent.len() < atoms.nr() {
                groups.push(IndexGroup::new("non-Water", aid_solvent));
            }
        } else if strcasecmp(cat, "Ion") == 0 && !aid_category.is_empty() {
            groups.push(IndexGroup::new(cat, aid_category));
        } else if !aid_category.is_empty() && !have_analysed_other {
            groups.push(IndexGroup::new(cat, aid_category));
            analyse_other(&restype, atoms, &mut groups);
            have_analysed_other = true;
        }
    }

    let mut iwater = None;
    let mut iion = None;
    let mut nwater = 0;
    let mut nion = 0;
    for (i, g) in groups.iter().enumerate() {
        if strcasecmp(&g.name, "Water") == 0 {
            iwater = Some(i);
            nwater = g.particle_indices.len();
        } else if strcasecmp(&g.name, "Ion") == 0 {
            iion = Some(i);
            nion = g.particle_indices.len();
        }
    }
    if nwater > 0 && nion > 0 {
        let mut a: Vec<usize> = Vec::new();
        a.extend_from_slice(&groups[iwater.unwrap()].particle_indices);
        a.extend_from_slice(&groups[iion.unwrap()].particle_indices);
        groups.push(IndexGroup::new("Water_and_ions", a));
    }

    groups
}
