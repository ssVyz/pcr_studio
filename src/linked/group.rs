//! Groups picked FASTQ files into read sets by their names: R1/R2 files of
//! one sample become pairs (several lanes together), the chunk files of a
//! Nanopore run one set.
//!
//! Recognised names (extension `.fastq`/`.fq`, optionally `.gz`):
//! - Illumina: `Sample_S1_L001_R1_001` (lane `L001` and the `_001` after
//!   `R1`/`R2` are not part of the sample name),
//! - `_R1`/`_R2` anywhere in the name, or a final `_1`/`_2` (`SRR123_1`)
//!   when exactly these two files exist,
//! - Nanopore (MinKNOW/Guppy, with a `pass`, `fail` or `runid` part):
//!   trailing chunk numbers are dropped (`FAX123_pass_barcode01_ab12cd34_0`).
//!
//! Other files are each a sample of their own.

use super::{Layout, Role};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

/// Whether the file name says FASTQ (`.fastq`, `.fq`, optionally `.gz`).
pub fn is_fastq(path: &Path) -> bool {
    let name = path.file_name().map(|s| s.to_string_lossy().to_lowercase()).unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    name.ends_with(".fastq") || name.ends_with(".fq")
}

/// FASTQ files in `dir` and its subfolders (up to three levels), sorted.
pub fn fastq_files_in(dir: &Path) -> Vec<PathBuf> {
    fn walk(dir: &Path, depth: usize, out: &mut Vec<PathBuf>) {
        let Ok(entries) = std::fs::read_dir(dir) else { return };
        for e in entries.flatten() {
            let p = e.path();
            if p.is_dir() {
                if depth < 3 {
                    walk(&p, depth + 1, out);
                }
            } else if is_fastq(&p) {
                out.push(p);
            }
        }
    }
    let mut out = Vec::new();
    walk(dir, 0, &mut out);
    out.sort();
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Grouping {
    /// One read set per sample (pairs and lanes together).
    BySample,
    /// One read set per file, or per R1/R2 pair.
    Separate,
    /// One read set with all files.
    Together,
}

impl Grouping {
    pub const ALL: [Grouping; 3] = [Grouping::BySample, Grouping::Separate, Grouping::Together];
}

impl std::fmt::Display for Grouping {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Grouping::BySample => "One read set per sample",
            Grouping::Separate => "One read set per file or pair",
            Grouping::Together => "All files in one read set",
        })
    }
}

/// A read set to create from picked files. Paired: R1, R2 alternate.
#[derive(Debug, Clone, PartialEq)]
pub struct Proposal {
    pub name: String,
    pub layout: Layout,
    pub files: Vec<(PathBuf, Role)>,
}

/// A picked file and what its name says.
struct Parsed {
    path: PathBuf,
    dir: PathBuf,
    stem: String,
    /// Name parts between `_` and `.` (byte ranges in `stem`).
    tokens: Vec<(usize, usize)>,
    /// R1/R2 and the index of the token saying so.
    mate: Option<(Role, usize)>,
}

impl Parsed {
    fn new(path: &Path) -> Parsed {
        let stem = super::stem(path);
        let mut tokens = Vec::new();
        let mut start = 0;
        for (i, c) in stem.char_indices() {
            if c == '_' || c == '.' {
                tokens.push((start, i));
                start = i + 1;
            }
        }
        tokens.push((start, stem.len()));
        let dir = path.parent().map(Path::to_path_buf).unwrap_or_default();
        let mut p = Parsed { path: path.to_path_buf(), dir, stem, tokens, mate: None };
        p.mate = (0..p.tokens.len()).rev().find_map(|i| match p.token(i).to_ascii_uppercase().as_str() {
            "R1" => Some((Role::R1, i)),
            "R2" => Some((Role::R2, i)),
            _ => None,
        });
        p
    }

    fn token(&self, i: usize) -> &str {
        let (a, b) = self.tokens[i];
        &self.stem[a..b]
    }

    fn last(&self) -> usize {
        self.tokens.len() - 1
    }

    /// Name before the last token (for `_1`/`_2` mates).
    fn prefix(&self) -> &str {
        &self.stem[..self.tokens[self.last()].0]
    }

    /// The name without the tokens in `drop`, keeping the separators of the rest.
    fn without(&self, drop: &HashSet<usize>) -> String {
        let mut out = String::new();
        for (i, &(a, b)) in self.tokens.iter().enumerate() {
            if drop.contains(&i) {
                continue;
            }
            if !out.is_empty() {
                out.push_str(&self.stem[a - 1..a]);
            }
            out.push_str(&self.stem[a..b]);
        }
        out
    }

    /// Name with the mate token blanked out: equal for the two files of a pair.
    fn mate_key(&self) -> String {
        match self.mate {
            Some((_, i)) => {
                let (a, b) = self.tokens[i];
                format!("{}\u{0}{}", &self.stem[..a], &self.stem[b..])
            }
            None => self.stem.clone(),
        }
    }

    /// Sample name: without mate, lane and Illumina's `_001`, and without the
    /// chunk numbers of Nanopore files.
    fn sample(&self) -> String {
        let is_num = |i: usize| !self.token(i).is_empty() && self.token(i).bytes().all(|b| b.is_ascii_digit());
        let mut drop = HashSet::new();
        if let Some((_, m)) = self.mate {
            drop.insert(m);
            if m < self.last() && is_num(m + 1) && self.token(m + 1).len() == 3 {
                drop.insert(m + 1);
            }
        }
        for i in 0..self.tokens.len() {
            let t = self.token(i);
            if t.len() == 4 && t.starts_with('L') && t[1..].bytes().all(|b| b.is_ascii_digit()) {
                drop.insert(i);
            }
        }
        let nanopore = (0..self.tokens.len()).any(|i| matches!(self.token(i), "pass" | "fail" | "runid"));
        if nanopore {
            let mut i = self.last();
            while i > 0 && is_num(i) {
                drop.insert(i);
                i -= 1;
            }
        }
        let name = self.without(&drop);
        if name.is_empty() { self.stem.clone() } else { name }
    }
}

/// Final `_1`/`_2` mates: a last token "1" or "2" counts when, among the
/// files with the same folder and name before it, exactly "1" and "2" exist.
fn number_mates(files: &mut [Parsed]) {
    let mut endings: HashMap<(PathBuf, String), HashSet<String>> = HashMap::new();
    for f in files.iter().filter(|f| f.mate.is_none() && f.tokens.len() > 1) {
        let last = f.token(f.last());
        if last.bytes().all(|b| b.is_ascii_digit()) {
            endings.entry((f.dir.clone(), f.prefix().to_string())).or_default().insert(last.to_string());
        }
    }
    let pair: HashSet<String> = ["1", "2"].iter().map(|s| s.to_string()).collect();
    for f in files.iter_mut().filter(|f| f.mate.is_none() && f.tokens.len() > 1) {
        if endings.get(&(f.dir.clone(), f.prefix().to_string())) == Some(&pair) {
            let last = f.last();
            f.mate = Some((if f.token(last) == "1" { Role::R1 } else { Role::R2 }, last));
        }
    }
}

/// The files of one group as a read set: paired if every file has its mate
/// (pairs ordered by name, i.e. by lane), single-end otherwise.
fn make_set(name: String, files: Vec<&Parsed>, pair: bool) -> Proposal {
    if pair && !files.is_empty() && files.iter().all(|f| f.mate.is_some()) {
        let mut by_key: BTreeMap<String, (Option<&Parsed>, Option<&Parsed>)> = BTreeMap::new();
        for &f in &files {
            let e = by_key.entry(f.mate_key()).or_default();
            match f.mate {
                Some((Role::R1, _)) if e.0.is_none() => e.0 = Some(f),
                Some((Role::R2, _)) if e.1.is_none() => e.1 = Some(f),
                _ => e.0 = None, // a duplicate: no clean pairs
            }
        }
        if by_key.values().all(|(a, b)| a.is_some() && b.is_some()) && by_key.len() * 2 == files.len() {
            let files = by_key
                .into_values()
                .flat_map(|(a, b)| [(a.unwrap().path.clone(), Role::R1), (b.unwrap().path.clone(), Role::R2)])
                .collect();
            return Proposal { name, layout: Layout::Paired, files };
        }
    }
    let mut paths: Vec<PathBuf> = files.iter().map(|f| f.path.clone()).collect();
    paths.sort();
    Proposal { name, layout: Layout::Single, files: paths.into_iter().map(|p| (p, Role::Single)).collect() }
}

/// Longest common start of the names, without trailing separators.
fn common_name(names: &[String]) -> String {
    let Some(first) = names.first() else { return String::new() };
    let mut len = first.len();
    for n in &names[1..] {
        len = len.min(first.bytes().zip(n.bytes()).take_while(|(a, b)| a == b).count());
    }
    while !first.is_char_boundary(len) {
        len -= 1;
    }
    first[..len].trim_end_matches(['_', '.', '-', ' ']).to_string()
}

/// Proposes read sets for the picked files. With `pair` off, R1/R2 in names
/// is ignored (every file is single-end).
pub fn propose(paths: &[PathBuf], pair: bool, grouping: Grouping) -> Vec<Proposal> {
    let mut seen = HashSet::new();
    let mut files: Vec<Parsed> = paths.iter().filter(|p| seen.insert((*p).clone())).map(|p| Parsed::new(p)).collect();
    if pair {
        number_mates(&mut files);
    } else {
        for f in &mut files {
            f.mate = None;
        }
    }
    match grouping {
        Grouping::Together => {
            if files.is_empty() {
                return Vec::new();
            }
            let samples: Vec<String> = files.iter().map(Parsed::sample).collect();
            let name = match common_name(&samples) {
                n if n.is_empty() => "Read set".to_string(),
                n => n,
            };
            vec![make_set(name, files.iter().collect(), pair)]
        }
        Grouping::BySample | Grouping::Separate => {
            let mut groups: BTreeMap<(PathBuf, String), Vec<&Parsed>> = BTreeMap::new();
            for f in &files {
                let key = if grouping == Grouping::BySample { f.sample() } else { f.mate_key() };
                groups.entry((f.dir.clone(), key)).or_default().push(f);
            }
            let mut out = Vec::new();
            for ((_, key), members) in groups {
                let name = if grouping == Grouping::BySample { key } else { members[0].sample() };
                let set = make_set(name, members.clone(), pair);
                if grouping == Grouping::Separate && set.layout == Layout::Single && members.len() > 1 {
                    // Mates without a partner: one set per file.
                    out.extend(members.iter().map(|&f| make_set(f.stem.clone(), vec![f], false)));
                } else {
                    out.push(set);
                }
            }
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn paths(names: &[&str]) -> Vec<PathBuf> {
        names.iter().map(|n| PathBuf::from(format!("/data/{n}"))).collect()
    }

    fn summary(p: &[Proposal]) -> Vec<(String, Layout, Vec<String>)> {
        p.iter()
            .map(|s| {
                let files = s.files.iter().map(|(f, r)| format!("{}:{}", r.label(), f.file_name().unwrap().to_string_lossy())).collect();
                (s.name.clone(), s.layout, files)
            })
            .collect()
    }

    #[test]
    fn illumina_samples_and_lanes() {
        let files = paths(&[
            "S1_S1_L002_R2_001.fastq.gz",
            "S1_S1_L001_R1_001.fastq.gz",
            "S1_S1_L001_R2_001.fastq.gz",
            "S1_S1_L002_R1_001.fastq.gz",
            "S2_S2_L001_R1_001.fastq.gz",
            "S2_S2_L001_R2_001.fastq.gz",
        ]);
        let p = propose(&files, true, Grouping::BySample);
        assert_eq!(
            summary(&p),
            [
                (
                    "S1_S1".to_string(),
                    Layout::Paired,
                    vec![
                        "R1:S1_S1_L001_R1_001.fastq.gz".to_string(),
                        "R2:S1_S1_L001_R2_001.fastq.gz".into(),
                        "R1:S1_S1_L002_R1_001.fastq.gz".into(),
                        "R2:S1_S1_L002_R2_001.fastq.gz".into()
                    ]
                ),
                ("S2_S2".to_string(), Layout::Paired, vec!["R1:S2_S2_L001_R1_001.fastq.gz".to_string(), "R2:S2_S2_L001_R2_001.fastq.gz".into()]),
            ]
        );
        // One set per pair; pairing off; everything together.
        assert_eq!(propose(&files, true, Grouping::Separate).len(), 3);
        // Without pairing, R1 and R2 are separate samples (lanes still together).
        let single = propose(&files, false, Grouping::BySample);
        let names: Vec<(&str, usize)> = single.iter().map(|s| (s.name.as_str(), s.files.len())).collect();
        assert_eq!(names, [("S1_S1_R1_001", 2), ("S1_S1_R2_001", 2), ("S2_S2_R1_001", 1), ("S2_S2_R2_001", 1)]);
        assert!(single.iter().all(|s| s.layout == Layout::Single && s.files[0].1 == Role::Single));
        let all = propose(&files, true, Grouping::Together);
        assert_eq!((all.len(), all[0].layout, all[0].files.len(), all[0].name.as_str()), (1, Layout::Paired, 6, "S"));
    }

    #[test]
    fn sra_nanopore_and_others() {
        let p = propose(&paths(&["SRR123_2.fastq.gz", "SRR123_1.fastq.gz"]), true, Grouping::BySample);
        assert_eq!(summary(&p), [("SRR123".to_string(), Layout::Paired, vec!["R1:SRR123_1.fastq.gz".to_string(), "R2:SRR123_2.fastq.gz".into()])]);
        // Numbered files that are not pairs, and chunks of a Nanopore run.
        let p = propose(
            &paths(&["sample_1.fq", "sample_2.fq", "sample_3.fq", "FAX1_pass_barcode01_ab12_0.fastq.gz", "FAX1_pass_barcode01_ab12_1.fastq.gz"]),
            true,
            Grouping::BySample,
        );
        let names: Vec<(&str, Layout, usize)> = p.iter().map(|s| (s.name.as_str(), s.layout, s.files.len())).collect();
        assert_eq!(
            names,
            [
                ("FAX1_pass_barcode01_ab12", Layout::Single, 2),
                ("sample_1", Layout::Single, 1),
                ("sample_2", Layout::Single, 1),
                ("sample_3", Layout::Single, 1)
            ]
        );
        // An R1 without its R2 is single-end.
        let p = propose(&paths(&["x_R1.fastq", "y_R1.fastq", "y_R2.fastq"]), true, Grouping::BySample);
        let names: Vec<(&str, Layout, usize)> = p.iter().map(|s| (s.name.as_str(), s.layout, s.files.len())).collect();
        assert_eq!(names, [("x", Layout::Single, 1), ("y", Layout::Paired, 2)]);
        assert!(is_fastq(Path::new("a.FQ.gz")) && is_fastq(Path::new("b.fastq")) && !is_fastq(Path::new("c.fasta.gz")));
        assert_eq!(common_name(&["ab_1".into(), "ab_2".into()]), "ab");
    }
}
