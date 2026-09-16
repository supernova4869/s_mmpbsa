//! Command line front ends mirroring the four extracted GROMACS tools.

pub mod coords;
pub mod convert_tpr;
pub mod dump;
pub mod make_ndx;
pub mod trjconv;

use std::cell::RefCell;
use std::collections::VecDeque;
use std::fmt::Write as _;

use crate::index::IndexGroup;

thread_local! {
    /// Lines that the tools read instead of standard input.  Installed with
    /// [`set_scripted_input`] when the front ends are driven from another
    /// program instead of a terminal.
    static SCRIPTED_INPUT: RefCell<Option<VecDeque<String>>> =
        const { RefCell::new(None) };
}

/// Answers the interactive questions of the tools (index group selections and
/// `make_ndx` editor commands) from a list of lines instead of standard input.
///
/// This exists so that a program can call [`trjconv::run`],
/// [`convert_tpr::run`] or [`make_ndx::run`] in-process and still choose the
/// groups; the interactive behaviour is unchanged when no script is
/// installed.  The script is per thread, and it is consumed in order, one line
/// per question.  Running out of lines behaves like end of input.
///
/// ```no_run
/// gmx_rs_tools::cmd::set_scripted_input(["Protein", "System"]);
/// gmx_rs_tools::cmd::trjconv::run(vec!["-f".into(), "traj.xtc".into(), "-o".into(), "out.gro".into()]);
/// gmx_rs_tools::cmd::clear_scripted_input();
/// ```
pub fn set_scripted_input<I, S>(lines: I)
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let queue: VecDeque<String> = lines.into_iter().map(Into::into).collect();
    SCRIPTED_INPUT.with(|s| *s.borrow_mut() = Some(queue));
}

/// Removes the script installed by [`set_scripted_input`], restoring the
/// interactive standard input.
pub fn clear_scripted_input() {
    SCRIPTED_INPUT.with(|s| *s.borrow_mut() = None);
}

/// True while the tools are driven by a script instead of a terminal.
pub fn is_scripted() -> bool {
    SCRIPTED_INPUT.with(|s| s.borrow().is_some())
}

/// Reads one line of input: the next scripted line when a script is installed,
/// otherwise a line from standard input.  Returns `None` at the end of the
/// input.
pub fn read_input_line() -> Option<String> {
    let scripted = SCRIPTED_INPUT.with(|s| s.borrow_mut().as_mut().and_then(|q| q.pop_front()));
    if let Some(line) = scripted {
        return Some(line);
    }
    if is_scripted() {
        // The script ran out: report end of input instead of blocking.
        return None;
    }
    let mut buf = String::new();
    match std::io::stdin().read_line(&mut buf) {
        Ok(0) | Err(_) => None,
        Ok(_) => Some(buf),
    }
}

/// Asks the user to pick an index group, matching `qgroup()`/`rd_groups()`
/// from `gromacs/topology/index.cpp`.
///
/// Returns `None` when standard input is exhausted, which GROMACS treats as a
/// fatal "Cannot read from input" error.
pub fn select_group(groups: &[IndexGroup], prompt: &str) -> Option<usize> {
    use std::io::Write as _;
    if groups.is_empty() {
        eprintln!("Error: no groups in indexfile");
        return None;
    }
    if !is_scripted() {
        for (i, g) in groups.iter().enumerate() {
            eprintln!(
                "Group {:5} ({:15}) has {:5} elements",
                i,
                g.name,
                g.particle_indices.len()
            );
        }
    }
    if groups.len() == 1 {
        // `qgroup()`: `fprintf(stderr, "There is one group in the index\n")`.
        if !is_scripted() {
            eprintln!("There is one group in the index");
        }
        return Some(0);
    }
    loop {
        if !is_scripted() {
            eprint!("{prompt} ");
            let _ = std::io::stderr().flush();
        }
        let line = match read_input_line() {
            Some(line) => line,
            None => {
                eprintln!("Cannot read from input");
                return None;
            }
        };
        let s = line.trim();
        if s.is_empty() {
            continue;
        }
        if let Ok(v) = s.parse::<usize>() {
            if v < groups.len() {
                // `qgroup()`: `printf("Selected %d: '%s'\n", ...)` on stdout.
                println!("Selected {}: '{}'", v, groups[v].name);
                return Some(v);
            }
        } else {
            let g = crate::index::find_group(s, groups);
            if g >= 0 {
                println!("Selected {}: '{}'", g, groups[g as usize].name);
                return Some(g as usize);
            }
        }
        println!("Error: No such group '{s}'");
    }
}

/// Very small gmx-style command line parser: options start with `-` and take
/// the following token(s) as values unless the next token also looks like an
/// option.
#[derive(Debug, Default)]
pub struct Args {
    pub opts: Vec<(String, Vec<String>)>,
}

fn looks_like_option(tok: &str) -> bool {
    if !tok.starts_with('-') {
        return false;
    }
    // Negative numbers are values, not options.
    tok[1..]
        .chars()
        .next()
        .map(|c| !(c.is_ascii_digit() || c == '.'))
        .unwrap_or(false)
}

impl Args {
    pub fn parse<I: IntoIterator<Item = String>>(argv: I) -> Args {
        let mut args = Args::default();
        let mut it = argv.into_iter().peekable();
        while let Some(tok) = it.next() {
            if looks_like_option(&tok) {
                let name = tok.trim_start_matches('-').to_string();
                let mut values = Vec::new();
                while let Some(next) = it.peek() {
                    if looks_like_option(next) || next.starts_with('-') {
                        break;
                    }
                    values.push(it.next().unwrap());
                }
                args.opts.push((name, values));
            }
        }
        args
    }

    pub fn get(&self, name: &str) -> Option<String> {
        self.opts
            .iter()
            .find(|(n, _)| n == name)
            .and_then(|(_, v)| v.first().cloned())
    }

    pub fn get_all(&self, name: &str) -> Vec<String> {
        self.opts
            .iter()
            .filter(|(n, _)| n == name)
            .flat_map(|(_, v)| v.clone())
            .collect()
    }

    pub fn has(&self, name: &str) -> bool {
        self.opts.iter().any(|(n, _)| n == name)
    }

    /// gmx style boolean option: `-sep` or `-sep yes` / `-sep no`.
    pub fn flag(&self, name: &str) -> bool {
        match self.opts.iter().find(|(n, _)| n == name) {
            None => false,
            Some((_, v)) => match v.first() {
                None => true,
                Some(s) => !matches!(s.as_str(), "no" | "false" | "0"),
            },
        }
    }

    /// gmx style boolean option with a `true` default: `-nr`/`-nonr` or
    /// `-nr yes|no`, as in `BooleanOption("nr").defaultValue(true)`.
    pub fn flag_default_true(&self, name: &str) -> bool {
        let mut value = true;
        if self.has(&format!("no{name}")) {
            value = false;
        }
        if let Some((_, v)) = self.opts.iter().find(|(n, _)| n == name) {
            value = match v.first() {
                None => true,
                Some(s) => !matches!(s.as_str(), "no" | "false" | "0"),
            };
        }
        value
    }

    pub fn int(&self, name: &str, default: i64) -> i64 {
        self.get(name)
            .and_then(|s| s.parse().ok())
            .unwrap_or(default)
    }

    pub fn real(&self, name: &str, default: f64) -> f64 {
        self.get(name)
            .and_then(|s| s.parse().ok())
            .unwrap_or(default)
    }

    pub fn real3(&self, name: &str, default: [f64; 3]) -> [f64; 3] {
        let vals = self.get_all(name);
        if vals.len() >= 3 {
            [
                vals[0].parse().unwrap_or(default[0]),
                vals[1].parse().unwrap_or(default[1]),
                vals[2].parse().unwrap_or(default[2]),
            ]
        } else {
            default
        }
    }
}

/// C `%g` style formatting.
pub fn fmt_g(v: f64, width: usize, precision: usize) -> String {
    let precision = if precision == 0 { 1 } else { precision };
    let s = if v == 0.0 {
        "0".to_string()
    } else {
        // `%g` picks the `%e` style when the exponent of the value rounded to
        // `precision` significant digits is below -4 or not below the
        // precision.  Using the rounded exponent matters: 9.9999997e-05, the
        // single precision representation of 1e-4, has exponent -4 and is
        // therefore printed as 0.0001.
        let exp = {
            let s = format!("{:.*e}", precision - 1, v);
            match s.find('e') {
                Some(pos) => s[pos + 1..].parse::<i32>().unwrap_or(0),
                None => 0,
            }
        };
        if exp < -4 || exp >= precision as i32 {
            let mut s = format!("{:.*e}", precision - 1, v);
            // Normalise the exponent to at least two digits with a sign.
            if let Some(pos) = s.find('e') {
                let (mant, e) = s.split_at(pos);
                let e = &e[1..];
                let (sign, digits) = if let Some(d) = e.strip_prefix('-') {
                    ("-", d)
                } else {
                    ("+", e.strip_prefix('+').unwrap_or(e))
                };
                s = format!("{mant}e{sign}{digits:0>2}");
            }
            if s.contains("e") {
                let pos = s.find('e').unwrap();
                let (mant, e) = s.split_at(pos);
                let mant = mant.trim_end_matches('0').trim_end_matches('.');
                s = format!("{mant}{e}");
            }
            s
        } else {
            let decimals = (precision as i32 - 1 - exp).max(0) as usize;
            let s = format!("{:.*}", decimals, v);
            if s.contains('.') {
                s.trim_end_matches('0').trim_end_matches('.').to_string()
            } else {
                s
            }
        }
    };
    if s.len() < width {
        let mut out = String::new();
        let _ = write!(out, "{:>width$}", s, width = width);
        out
    } else {
        s
    }
}

/// C `%e` style formatting with a two digit exponent.
pub fn fmt_e(v: f64, width: usize, precision: usize) -> String {
    let mut s = format!("{:.*e}", precision, v);
    if let Some(pos) = s.find('e') {
        let (mant, e) = s.split_at(pos);
        let e = &e[1..];
        let (sign, digits) = if let Some(d) = e.strip_prefix('-') {
            ("-", d)
        } else {
            ("+", e.strip_prefix('+').unwrap_or(e))
        };
        s = format!("{mant}e{sign}{digits:0>2}");
    }
    if s.len() < width {
        format!("{:>width$}", s, width = width)
    } else {
        s
    }
}

pub fn indent_str(n: usize) -> String {
    " ".repeat(n)
}
