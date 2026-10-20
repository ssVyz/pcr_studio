//! primersearch (sister repo `C:\repos\primersearch`): finds a small set of
//! primers that together cover every sequence of an aligned set. PCR Studio
//! runs it on the selected columns and keeps its text output as a library
//! document.
//!
//! Invocation follows the tool's README ("Programmatic / subprocess use"):
//! `--no-config --silent --progress jsonl`, the text report on stdout (`-o -`),
//! progress as JSON lines on stderr. Every parameter is passed explicitly.

use crate::fasta::Record;
use crate::model::{Document, RefCoords};
use crate::seq::{self, GAP};
use crate::ui::canvas::group_digits;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::ffi::OsString;
use std::path::Path;

/// How primer variants are formed (`--mode`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mode {
    /// Exact-match primers only.
    #[default]
    NoAmbiguities,
    /// IUPAC codes absorb related variants, up to a budget.
    Incremental,
    /// Exhaustive search for the best set of n oligos, allowing mismatches.
    OptimizeByMismatch,
}

impl Mode {
    pub const ALL: [Mode; 3] = [Mode::NoAmbiguities, Mode::Incremental, Mode::OptimizeByMismatch];

    fn cli(self) -> &'static str {
        match self {
            Mode::NoAmbiguities => "no-ambiguities",
            Mode::Incremental => "incremental",
            Mode::OptimizeByMismatch => "optimize-by-mismatch",
        }
    }

    /// Whether `--exclude-n` / `--only-twofold` apply.
    pub fn uses_iupac(self) -> bool {
        self != Mode::NoAmbiguities
    }
}

impl std::fmt::Display for Mode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Mode::NoAmbiguities => "No ambiguities (exact match)",
            Mode::Incremental => "Incremental (IUPAC codes)",
            Mode::OptimizeByMismatch => "Optimize by mismatch",
        })
    }
}

/// Which sequences count in optimize-by-mismatch (`--mismatch-mode`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Criterion {
    #[default]
    AtMost,
    Exactly,
}

impl Criterion {
    pub const ALL: [Criterion; 2] = [Criterion::AtMost, Criterion::Exactly];

    fn cli(self) -> &'static str {
        match self {
            Criterion::AtMost => "lower-or-equal",
            Criterion::Exactly => "exact",
        }
    }
}

impl std::fmt::Display for Criterion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Criterion::AtMost => "At most the mismatches",
            Criterion::Exactly => "Exactly the mismatches",
        })
    }
}

/// Run parameters; the defaults are the tool's own.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Params {
    pub mode: Mode,
    /// The selection is the primer site: only generate the variants that
    /// cover it (`--fixed`; the Tm threshold is not enforced).
    pub fixed: bool,
    pub reverse: bool,
    /// Minimum Tm (°C).
    pub tm: f64,
    /// 3'-terminal bases that must be perfectly conserved.
    pub three_prime: usize,
    /// Incremental: coverage (%) at which more ambiguity codes are allowed.
    pub target: f64,
    /// Incremental: ambiguity codes per primer.
    pub max_amb: usize,
    /// Incremental: seeds per range (0 = no cap).
    pub max_seeds: usize,
    /// No N, and only two-fold codes (incremental, optimize by mismatch).
    pub exclude_n: bool,
    pub only_twofold: bool,
    /// Optimize by mismatch: coverage criterion, set size, mismatches,
    /// ambiguity codes per oligo and the search limits (0 = no limit).
    pub criterion: Criterion,
    pub n_oligos: usize,
    pub mismatches: usize,
    pub ambiguities: usize,
    pub max_candidates: u64,
    pub max_work: u64,
    /// Primer sequences in groups of three bases.
    pub spacer: bool,
    /// Worker threads; 0 = all cores.
    pub threads: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            mode: Mode::NoAmbiguities,
            fixed: false,
            reverse: false,
            tm: 55.0,
            three_prime: 0,
            target: 50.0,
            max_amb: 1,
            max_seeds: 50,
            exclude_n: false,
            only_twofold: false,
            criterion: Criterion::AtMost,
            n_oligos: 1,
            mismatches: 0,
            ambiguities: 0,
            max_candidates: 500_000_000,
            max_work: 2_000_000_000,
            spacer: true,
            threads: 0,
        }
    }
}

/// Reaction conditions in primersearch's units (oligo µM, salts mM).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Conditions {
    pub oligo_um: f64,
    pub na_mm: f64,
    pub mg_mm: f64,
    pub dntp_mm: f64,
}

/// Arguments of a run on `input`, starting with the input path.
pub fn args(p: &Params, c: &Conditions, input: &Path) -> Vec<OsString> {
    let mut a: Vec<String> = ["--format", "text", "-o", "-", "--no-config", "--silent", "--progress", "jsonl"].map(String::from).into();
    let mut flag = |name: &str, value: String| {
        a.push(name.into());
        a.push(value);
    };
    flag("--tm", p.tm.to_string());
    flag("--oligo", c.oligo_um.to_string());
    flag("--na", c.na_mm.to_string());
    flag("--mg", c.mg_mm.to_string());
    flag("--dntp", c.dntp_mm.to_string());
    flag("--mode", p.mode.cli().into());
    if p.three_prime > 0 {
        flag("--three-prime", p.three_prime.to_string());
    }
    match p.mode {
        Mode::NoAmbiguities => {}
        Mode::Incremental => {
            flag("--target", p.target.to_string());
            flag("--max-amb", p.max_amb.to_string());
            flag("--max-seeds", p.max_seeds.to_string());
        }
        Mode::OptimizeByMismatch => {
            flag("--mismatch-mode", p.criterion.cli().into());
            flag("--n-oligos", p.n_oligos.to_string());
            flag("--mismatches", p.mismatches.to_string());
            flag("--ambiguities", p.ambiguities.to_string());
            flag("--max-candidates", p.max_candidates.to_string());
            flag("--max-work", p.max_work.to_string());
        }
    }
    if p.threads > 0 {
        flag("-j", p.threads.to_string());
    }
    let switches = [
        (p.fixed, "--fixed"),
        (p.reverse, "--rev"),
        (p.mode.uses_iupac() && p.exclude_n, "--exclude-n"),
        (p.mode.uses_iupac() && p.only_twofold, "--only-twofold"),
        (true, if p.spacer { "--spacer" } else { "--no-spacer" }),
    ];
    a.extend(switches.into_iter().filter(|(on, _)| *on).map(|(_, s)| s.to_string()));
    std::iter::once(input.as_os_str().to_owned()).chain(a.into_iter().map(OsString::from)).collect()
}

/// The command line for the report, with a placeholder for the temporary input.
pub fn command_line(args: &[OsString]) -> String {
    let rest = args.iter().skip(1).map(|a| {
        let a = a.to_string_lossy();
        if a.contains(' ') { format!("\"{a}\"") } else { a.into_owned() }
    });
    std::iter::once("primersearch <input FASTA>".to_string()).chain(rest).collect::<Vec<_>>().join(" ")
}

/// Progress bar position (`None`: no estimate) and status text of a
/// `{"type":"progress"}` event. Greedy rounds restart their `pct` every
/// round, so the share of covered sequences is used there; the set search
/// cannot predict its total work and reports how much of `max_work` it used.
pub fn progress(ev: &Value) -> Option<(Option<f32>, String)> {
    if ev["type"] != "progress" {
        return None;
    }
    let phase = ev["phase"].as_str().unwrap_or("");
    let n = |k: &str| ev[k].as_u64();
    let fraction = match (phase, n("covered"), n("total")) {
        ("set_search", _, _) => None,
        ("round" | "fixed", Some(c), Some(t)) if t > 0 => Some(c as f32 / t as f32),
        _ => Some(ev["pct"].as_f64().unwrap_or(0.0) as f32 / 100.0),
    };
    let message = match (phase, n("candidates"), n("evaluations")) {
        ("set_search", Some(c), Some(e)) => {
            let mut m = format!("Searching sets over {} candidates · {} evaluations", group_digits(c as usize), group_digits(e as usize));
            if let Some(w) = n("max_work").filter(|&w| w > 0) {
                m.push_str(&format!(" ({:.1}% of max work)", e as f64 / w as f64 * 100.0));
            }
            m
        }
        _ => {
            let raw = ev["message"].as_str().unwrap_or("");
            let mut m = match raw.strip_prefix("Optimize by mismatch: ") {
                Some(rest) => rest.chars().take(1).flat_map(char::to_uppercase).chain(rest.chars().skip(1)).collect(),
                None => raw.to_string(),
            };
            if let ("reduce", Some(g)) = (phase, n("generated")) {
                m.push_str(&format!(" · {} candidates generated", group_digits(g as usize)));
            }
            m
        }
    };
    Some((fraction.map(|f| f.clamp(0.0, 1.0)), message))
}

/// The selected columns of the chosen rows, as passed to primersearch.
pub struct Input {
    pub slice: Document,
    /// Alignment column of every slice column.
    pub columns: Vec<usize>,
    /// Number of chosen rows (some may have no base in the selection).
    pub chosen: usize,
}

/// What primersearch will make of the input: it drops every sequence with a
/// gap (also no coverage) or an ambiguity code in the slice.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct InputStats {
    pub chosen: usize,
    /// Rows with at least one base in the selection (written to the input).
    pub written: usize,
    pub usable: usize,
    pub gaps: usize,
    pub ambiguous: usize,
    pub columns: usize,
}

impl Input {
    pub fn new(doc: &Document, rows: &[usize], c0: usize, c1: usize, only: Option<&[bool]>) -> Input {
        let (slice, columns) = doc.slice_rows(rows, c0, c1, only);
        Input { slice, columns, chosen: rows.len() }
    }

    /// Slice row `i` over all columns; gaps and no coverage as `-`.
    fn gapped(&self, i: usize) -> Vec<u8> {
        let r = &self.slice.rows[i];
        (0..self.columns.len()).map(|c| r.at(c).unwrap_or(GAP)).collect()
    }

    pub fn stats(&self) -> InputStats {
        let mut s = InputStats { chosen: self.chosen, written: self.slice.rows.len(), columns: self.columns.len(), ..Default::default() };
        for i in 0..self.slice.rows.len() {
            let g = self.gapped(i);
            if g.contains(&GAP) {
                s.gaps += 1;
            } else if !g.iter().all(|&b| seq::is_acgt(b)) {
                s.ambiguous += 1;
            } else {
                s.usable += 1;
            }
        }
        s
    }

    pub fn records(&self) -> Vec<Record> {
        (0..self.slice.rows.len())
            .map(|i| Record { name: self.slice.rows[i].name.clone(), description: String::new(), seq: self.gapped(i) })
            .collect()
    }
}

const MAX_STRETCHES: usize = 12;

/// Where the positions of the output lie in the alignment: one line per
/// stretch of consecutive reference positions (or alignment columns).
pub fn position_map(columns: &[usize], rc: Option<&RefCoords>) -> Vec<String> {
    let base = |c: usize| rc.is_some_and(|rc| rc.is_base.get(c).copied().unwrap_or(false));
    let follows = |prev: usize, c: usize| match rc {
        Some(rc) if base(prev) && base(c) => rc.pos_at[c] == rc.pos_at[prev] + 1,
        Some(_) => !base(prev) && !base(c) && c == prev + 1,
        None => c == prev + 1,
    };
    let mut stretches: Vec<(usize, usize)> = Vec::new();
    let mut start = 0;
    for i in 1..=columns.len() {
        if i == columns.len() || !follows(columns[i - 1], columns[i]) {
            stretches.push((start, i));
            start = i;
        }
    }
    let mut lines: Vec<String> = stretches
        .iter()
        .take(MAX_STRETCHES)
        .map(|&(a, b)| {
            let (c0, c1) = (columns[a], columns[b - 1] + 1);
            let cols = format!("columns {}–{}", group_digits(c0 + 1), group_digits(c1));
            let place = match rc {
                Some(rc) if base(c0) => {
                    let (x, y) = (rc.pos_at[c0] as usize, rc.pos_at[c1 - 1] as usize);
                    format!("ref {}–{} ({cols})", group_digits(x), group_digits(y))
                }
                Some(_) => format!("{cols} (insertion relative to the reference)"),
                None => cols,
            };
            format!("  {:>13}  →  {place}", format!("{}–{}", a + 1, b))
        })
        .collect();
    if stretches.len() > MAX_STRETCHES {
        lines.push(format!("  … {} more stretches", stretches.len() - MAX_STRETCHES));
    }
    lines
}

/// "38 s", "4 min 12 s", "1 h 05 min".
pub fn duration(seconds: f64) -> String {
    let s = seconds.max(0.0);
    if s < 10.0 {
        format!("{s:.1} s")
    } else if s < 60.0 {
        format!("{s:.0} s")
    } else if s < 3600.0 {
        format!("{} min {:02} s", (s / 60.0) as u64, (s % 60.0) as u64)
    } else {
        format!("{} h {:02} min", (s / 3600.0) as u64, ((s % 3600.0) / 60.0) as u64)
    }
}

/// Describes the run above the tool's own output.
pub struct Report<'a> {
    pub document: &'a str,
    /// Selected range, e.g. "ref 19,100–19,400 (columns 20,101–20,450)".
    pub selection: &'a str,
    /// The columns with gaps in the reference were left out.
    pub collapsed: bool,
    /// Which rows were chosen, e.g. "All sequences, without the reference".
    pub rows: &'a str,
    pub stats: InputStats,
    pub version: &'a str,
    pub path: &'a Path,
    pub command: &'a str,
    pub started: &'a str,
    pub finished: &'a str,
    pub seconds: f64,
    pub positions: &'a [String],
}

impl Report<'_> {
    /// The report header followed by the tool's output.
    pub fn text(&self, output: &str) -> String {
        let s = &self.stats;
        let d = group_digits;
        let mut selection = self.selection.to_string();
        if self.collapsed {
            selection.push_str("; columns with gaps in the reference left out");
        }
        let mut sequences = format!("{}: {} chosen", self.rows, d(s.chosen));
        if s.written < s.chosen {
            sequences.push_str(&format!(", {} without bases in the selection left out", d(s.chosen - s.written)));
        }
        let mut usable = d(s.usable);
        if s.gaps + s.ambiguous > 0 {
            usable.push_str(&format!(
                " of {} (primersearch drops {} with gaps or no coverage and {} with N/ambiguity codes)",
                d(s.written),
                d(s.gaps),
                d(s.ambiguous)
            ));
        }
        let fields = [
            ("Document", self.document.to_string()),
            ("Selection", selection),
            ("Columns", d(s.columns)),
            ("Sequences", sequences),
            ("Usable", usable),
            ("Tool", format!("primersearch {} ({})", self.version, self.path.display())),
            ("Command", self.command.to_string()),
            ("Started", self.started.to_string()),
            ("Finished", format!("{} ({})", self.finished, duration(self.seconds))),
        ];
        let title = "PCR STUDIO · PRIMERSEARCH RUN";
        let mut out = format!("{title}\n{}\n", "=".repeat(title.chars().count()));
        for (k, v) in fields {
            out.push_str(&format!("{:<12} {v}\n", format!("{k}:")));
        }
        out.push_str("\nPositions in the output count the columns passed to primersearch:\n");
        for l in self.positions {
            out.push_str(l);
            out.push('\n');
        }
        out.push('\n');
        out.push_str(output.trim_end());
        out.push('\n');
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Row;
    use serde_json::json;

    fn strs(a: &[OsString]) -> Vec<String> {
        a.iter().map(|s| s.to_string_lossy().into_owned()).collect()
    }

    fn doc(rows: &[&str]) -> Document {
        Document::new(rows.iter().enumerate().map(|(i, s)| Row::from_gapped(format!("s{i}"), String::new(), s.as_bytes().to_vec())).collect())
    }

    #[test]
    fn arguments() {
        let c = Conditions { oligo_um: 0.2, na_mm: 50.0, mg_mm: 3.0, dntp_mm: 0.8 };
        let a = strs(&args(&Params::default(), &c, Path::new("in.fasta")));
        assert_eq!(
            a.join(" "),
            "in.fasta --format text -o - --no-config --silent --progress jsonl --tm 55 --oligo 0.2 --na 50 --mg 3 --dntp 0.8 \
             --mode no-ambiguities --spacer"
        );
        // IUPAC restrictions only with ambiguity modes; mode options only for their mode.
        let p = Params { mode: Mode::Incremental, exclude_n: true, reverse: true, three_prime: 5, threads: 4, spacer: false, ..Params::default() };
        let a = strs(&args(&p, &c, Path::new("in.fasta"))).join(" ");
        assert!(a.contains("--mode incremental --three-prime 5 --target 50 --max-amb 1 --max-seeds 50 -j 4 --rev --exclude-n --no-spacer"), "{a}");
        assert!(!a.contains("--n-oligos"));
        let p = Params { mode: Mode::OptimizeByMismatch, fixed: true, criterion: Criterion::Exactly, n_oligos: 2, mismatches: 1, ..Params::default() };
        let a = strs(&args(&p, &c, Path::new("in.fasta"))).join(" ");
        assert!(a.contains("--mismatch-mode exact --n-oligos 2 --mismatches 1 --ambiguities 0 --max-candidates 500000000 --max-work 2000000000 --fixed"), "{a}");
        let p = Params { exclude_n: true, only_twofold: true, ..Params::default() };
        assert!(!strs(&args(&p, &c, Path::new("x"))).iter().any(|s| s == "--exclude-n" || s == "--only-twofold"));
        let line = command_line(&args(&Params::default(), &c, Path::new("C:/tmp/in.fasta")));
        assert!(line.starts_with("primersearch <input FASTA> --format text"));
    }

    #[test]
    fn progress_events() {
        let round = json!({"type":"progress","message":"Round 2","pct":10.0,"phase":"round","round":2,"covered":30,"total":120});
        assert_eq!(progress(&round), Some((Some(0.25), "Round 2".into())));
        let search = json!({"type":"progress","message":"x","pct":70.0,"phase":"set_search","candidates":1700,"evaluations":5000,"max_work":10000});
        assert_eq!(progress(&search), Some((None, "Searching sets over 1,700 candidates · 5,000 evaluations (50.0% of max work)".into())));
        let reduce = json!({"type":"progress","message":"Optimize by mismatch: reducing","pct":60.0,"phase":"reduce","generated":12345});
        assert_eq!(progress(&reduce), Some((Some(0.6), "Reducing · 12,345 candidates generated".into())));
        assert_eq!(progress(&json!({"type":"error","message":"boom"})), None);
    }

    #[test]
    fn input_and_positions() {
        // Reference AC-GTA: column 2 is an insertion in the other rows.
        let d = doc(&["AC-GTA", "ACTGTA", "ACTGNA", "AC-GTA", "--TG--", "------A"]);
        let rc = RefCoords::new(&d.rows[0], d.width);
        let rows = [1, 2, 3, 4, 5];
        // All columns: row 3 has a gap in the insertion, row 4 is partial, row 2 has an N.
        let input = Input::new(&d, &rows, 0, 6, None);
        assert_eq!(input.stats(), InputStats { chosen: 5, written: 4, usable: 1, gaps: 2, ambiguous: 1, columns: 6 });
        assert_eq!(input.records()[3].seq, b"--TG--");
        let map = position_map(&input.columns, Some(&rc));
        assert_eq!(map.len(), 3);
        assert!(map[0].ends_with("1–2  →  ref 1–2 (columns 1–2)"), "{map:?}");
        assert!(map[1].ends_with("3–3  →  columns 3–3 (insertion relative to the reference)"));
        assert!(map[2].ends_with("4–6  →  ref 3–5 (columns 4–6)"));
        // Collapsed gaps: the insertion column is left out; one stretch.
        let input = Input::new(&d, &rows, 0, 6, Some(&rc.is_base));
        assert_eq!(input.stats().usable, 2);
        assert_eq!(position_map(&input.columns, Some(&rc)), vec![format!("  {:>13}  →  ref 1–5 (columns 1–6)", "1–5")]);
        assert_eq!(position_map(&[4, 5, 9], None).len(), 2);
    }

    /// Runs the bundled primersearch (`assets/tools`) through the runner:
    /// `cargo test -- --ignored bundled_primersearch`.
    #[test]
    #[ignore]
    fn bundled_primersearch() {
        use crate::tools::{self, PRIMERSEARCH, ToolConfig, ToolStatus, run};
        use std::sync::atomic::AtomicBool;
        let ToolStatus::Ready { path, .. } = tools::check(&PRIMERSEARCH, &ToolConfig::default()) else {
            panic!("primersearch is not in assets/tools");
        };
        // 40 sequences; every fourth has a SNP, one has a gap.
        let base = b"ATGGCGTACGTTAGCATCGATCGGATCCTAGGCTAGCTAGGATCGATCGTAGCTAGCTAG";
        let rows = (0..40)
            .map(|i| {
                let mut s = base.to_vec();
                if i % 4 == 0 {
                    s[30] = b'A';
                }
                if i == 7 {
                    s[10] = GAP;
                }
                Row::from_gapped(format!("s{i}"), String::new(), s)
            })
            .collect();
        let input = Input::new(&Document::new(rows), &(0..40).collect::<Vec<_>>(), 0, 60, None);
        assert_eq!(input.stats().usable, 39);
        let file = run::TempFile::new("fasta");
        crate::fasta::write_fasta(&file.path, &input.records()).unwrap();
        let c = Conditions { oligo_um: 0.2, na_mm: 50.0, mg_mm: 3.0, dntp_mm: 0.8 };
        let go = |p: &Params, input: &Path, cancel: bool| {
            let mut cmd = run::command(&path);
            cmd.args(args(p, &c, input));
            let events = std::sync::Mutex::new(Vec::new());
            let out = run::run(cmd, &run::Handle::default(), &AtomicBool::new(cancel), &|ev| events.lock().unwrap().push(ev.clone())).unwrap();
            (out, events.into_inner().unwrap())
        };
        let (out, events) = go(&Params::default(), &file.path, false);
        assert_eq!(out.code, Some(0), "{}", out.failure());
        assert!(out.stdout.contains("PRIMER SEARCH RESULTS"), "{}", out.stdout);
        assert!(out.stdout.contains("Valid sequences:            39"), "{}", out.stdout);
        assert!(events.iter().any(|e| progress(e).is_some()), "{events:?}");
        let p = Params { mode: Mode::OptimizeByMismatch, fixed: true, mismatches: 1, ..Params::default() };
        let (out, _) = go(&p, &file.path, false);
        assert_eq!(out.code, Some(0), "{}", out.failure());
        // Failures come back with the tool's message; cancelling stops the run.
        let (out, _) = go(&Params::default(), Path::new("Z:/missing/input.fasta"), false);
        assert_eq!(out.code, Some(1));
        assert!(out.failure().contains("failed to read"), "{}", out.failure());
        assert!(go(&Params::default(), &file.path, true).0.cancelled);
    }

    #[test]
    fn report_text() {
        assert_eq!((duration(4.24), duration(38.0), duration(252.0), duration(3900.0)), ("4.2 s".into(), "38 s".into(), "4 min 12 s".into(), "1 h 05 min".into()));
        let positions = vec!["  1–20  →  ref 1–20".to_string()];
        let r = Report {
            document: "Doc",
            selection: "ref 1–20",
            collapsed: true,
            rows: "All sequences",
            stats: InputStats { chosen: 10, written: 9, usable: 7, gaps: 1, ambiguous: 1, columns: 20 },
            version: "0.1.1",
            path: Path::new("primersearch.exe"),
            command: "primersearch <input FASTA>",
            started: "2026-10-20 10:00:00",
            finished: "2026-10-20 10:04:12",
            seconds: 252.0,
            positions: &positions,
        };
        let t = r.text("RESULTS\n\n");
        assert!(t.contains("Selection:   ref 1–20; columns with gaps in the reference left out\n"), "{t}");
        assert!(t.contains("Usable:      7 of 9 (primersearch drops 1 with gaps or no coverage and 1 with N/ambiguity codes)\n"));
        assert!(t.contains("Sequences:   All sequences: 10 chosen, 1 without bases in the selection left out
"), "{t}");
        assert!(t.contains("Finished:    2026-10-20 10:04:12 (4 min 12 s)\n"));
        assert!(t.ends_with("  1–20  →  ref 1–20\n\nRESULTS\n"));
    }
}
