//! simpletrim: quality check and trimming of FASTQ reads. PCR Studio runs
//! `simpletrim analyze` on one file, or one lane's R1/R2 pair, of a linked
//! read set.
//!
//! The run is in text mode: the report then comes on stdout as simpletrim
//! writes it (with `--json` it would be JSON), and stderr carries plain
//! progress lines and `error: …` lines.

use super::run::Output;
use crate::linked::{LinkedFile, LinkedSet, Role};
use crate::ui::canvas::group_digits;
use std::ffi::OsString;
use std::path::{Path, PathBuf};

/// Reads (pairs) `simpletrim analyze` inspects unless told to read everything.
pub const DEFAULT_MAX_READS: u64 = 1_000_000;

/// What one analyze run reads: a single file, or the R1 and R2 files of a lane.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnalyzeInput {
    /// "Lane 2 of 4", "File 3 of 12"; empty for a set of one file or pair.
    pub part: String,
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
}

impl AnalyzeInput {
    pub fn files(&self) -> impl Iterator<Item = &Path> {
        std::iter::once(self.r1.as_path()).chain(self.r2.as_deref())
    }

    /// "Lane 2 of 4: x_R1.fastq.gz + x_R2.fastq.gz" (pick list entry).
    pub fn label(&self) -> String {
        let name = |p: &Path| p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| p.display().to_string());
        let files = match &self.r2 {
            Some(r2) => format!("{} + {}", name(&self.r1), name(r2)),
            None => name(&self.r1),
        };
        if self.part.is_empty() { files } else { format!("{}: {files}", self.part) }
    }
}

impl std::fmt::Display for AnalyzeInput {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.label())
    }
}

/// The runs a read set offers: each lane's R1/R2 pair of a paired set, each
/// file of a single-end set.
pub fn inputs(set: &LinkedSet) -> Result<Vec<AnalyzeInput>, String> {
    let sources = crate::linked::scan::sources(set)?;
    let n = sources.len();
    let word = if sources.first().is_some_and(|s| s.len() == 2) { "Lane" } else { "File" };
    sources
        .iter()
        .enumerate()
        .map(|(i, files)| {
            let part = if n > 1 { format!("{word} {} of {n}", i + 1) } else { String::new() };
            match files {
                [f] => Ok(AnalyzeInput { part, r1: f.path.clone(), r2: None }),
                [a, b] if a.role == Role::R1 && b.role == Role::R2 => Ok(AnalyzeInput { part, r1: a.path.clone(), r2: Some(b.path.clone()) }),
                _ => Err("The R1 and R2 files of the read set are not in order".into()),
            }
        })
        .collect()
}

/// The linked files an input reads (to check them before the run).
pub fn linked_files<'a>(set: &'a LinkedSet, input: &AnalyzeInput) -> Vec<&'a LinkedFile> {
    set.files.iter().filter(|f| input.files().any(|p| p == f.path)).collect()
}

/// Arguments of `simpletrim analyze`; `all_reads` lifts the read and base limits.
pub fn analyze_args(input: &AnalyzeInput, all_reads: bool, threads: usize) -> Vec<OsString> {
    let mut a: Vec<OsString> = vec!["analyze".into(), "--r1".into(), input.r1.clone().into()];
    if let Some(r2) = &input.r2 {
        a.extend(["--r2".into(), r2.clone().into()]);
    }
    if all_reads {
        a.extend(["--max-reads", "0", "--max-bases", "0"].map(OsString::from));
    }
    a.extend(["-w".into(), threads.max(1).to_string().into()]);
    a
}

/// Worker threads for a run: all cores.
pub fn threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4)
}

/// The command line for the report.
pub fn command_line(args: &[OsString]) -> String {
    std::iter::once("simpletrim".to_string())
        .chain(args.iter().map(|a| {
            let a = a.to_string_lossy();
            if a.contains(' ') { format!("\"{a}\"") } else { a.into_owned() }
        }))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Status text of a stderr line ("Reading up to 1,000,000 pairs ...",
/// "Detecting adapters ..."); `None` for blank and error lines.
pub fn status_line(line: &str) -> Option<String> {
    let line = line.trim();
    if line.is_empty() || line.starts_with("error:") {
        return None;
    }
    Some(line.trim_end_matches(" ...").trim_end_matches("...").to_string())
}

/// Reads (pairs) done according to a progress line ("read 2,000,000 pairs",
/// "read 300,000 reads").
pub fn reads_done(line: &str) -> Option<u64> {
    let mut words = line.trim().strip_prefix("read ")?.split_whitespace();
    let n = words.next()?.replace(',', "").parse().ok()?;
    matches!(words.next(), Some("pairs" | "reads")).then_some(n)
}

/// Why a run failed: simpletrim's `error: …` lines, else whatever it wrote.
pub fn error_text(out: &Output) -> String {
    let errors: Vec<&str> = out.stderr.lines().filter_map(|l| l.trim().strip_prefix("error:")).map(str::trim).collect();
    if errors.is_empty() { out.failure() } else { errors.join("\n") }
}

/// Describes the run above simpletrim's report.
pub struct Report<'a> {
    pub document: &'a str,
    /// "FASTQ · paired-end · 2 lanes".
    pub kind: &'a str,
    pub input: &'a AnalyzeInput,
    pub all_reads: bool,
    pub version: &'a str,
    pub path: &'a Path,
    pub command: &'a str,
    pub started: &'a str,
    pub finished: &'a str,
    pub seconds: f64,
}

impl Report<'_> {
    /// The report header followed by simpletrim's report.
    pub fn text(&self, output: &str) -> String {
        let unit = if self.input.r2.is_some() { "pairs" } else { "reads" };
        let reads = if self.all_reads {
            format!("all {unit}")
        } else {
            format!("up to {} {unit} from the start of the file (long reads: up to 500 Mb)", group_digits(DEFAULT_MAX_READS as usize))
        };
        let mut fields = vec![("Document", format!("{} ({})", self.document, self.kind))];
        if !self.input.part.is_empty() {
            fields.push(("Analyzed", self.input.part.to_lowercase()));
        }
        match &self.input.r2 {
            Some(r2) => fields.extend([("R1", self.input.r1.display().to_string()), ("R2", r2.display().to_string())]),
            None => fields.push(("File", self.input.r1.display().to_string())),
        }
        fields.extend([
            ("Reads", reads),
            ("Tool", format!("simpletrim {} ({})", self.version, self.path.display())),
            ("Command", self.command.to_string()),
            ("Started", self.started.to_string()),
            ("Finished", format!("{} ({})", self.finished, super::primersearch::duration(self.seconds))),
        ]);
        let title = "PCR STUDIO · SIMPLETRIM ANALYZE";
        let mut out = format!("{title}\n{}\n", "=".repeat(title.chars().count()));
        for (k, v) in fields {
            out.push_str(&format!("{:<10} {v}\n", format!("{k}:")));
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
    use crate::linked::{Format, Layout, Summary, Technology};

    fn file(path: &str, role: Role) -> LinkedFile {
        LinkedFile { path: path.into(), role, size: 1, modified: None, fingerprint: 0, content_hash: 0, records: 10, bases: 1000 }
    }

    fn set(layout: Layout, files: Vec<LinkedFile>) -> LinkedSet {
        LinkedSet { format: Format::Fastq, layout, technology: Technology::Illumina, files, summary: Summary::default() }
    }

    fn strs(a: &[OsString]) -> Vec<String> {
        a.iter().map(|s| s.to_string_lossy().into_owned()).collect()
    }

    #[test]
    fn inputs_of_read_sets() {
        let one = set(Layout::Single, vec![file("C:/r/a.fastq.gz", Role::Single)]);
        let i = inputs(&one).unwrap();
        assert_eq!(i, [AnalyzeInput { part: String::new(), r1: "C:/r/a.fastq.gz".into(), r2: None }]);
        assert_eq!(i[0].label(), "a.fastq.gz");

        let lanes = set(
            Layout::Paired,
            vec![
                file("C:/r/s_L001_R1.fq.gz", Role::R1),
                file("C:/r/s_L001_R2.fq.gz", Role::R2),
                file("C:/r/s_L002_R1.fq.gz", Role::R1),
                file("C:/r/s_L002_R2.fq.gz", Role::R2),
            ],
        );
        let i = inputs(&lanes).unwrap();
        assert_eq!(i.len(), 2);
        assert_eq!(i[1].label(), "Lane 2 of 2: s_L002_R1.fq.gz + s_L002_R2.fq.gz");
        assert_eq!(linked_files(&lanes, &i[1]).iter().map(|f| f.role).collect::<Vec<_>>(), [Role::R1, Role::R2]);

        let files = set(Layout::Single, vec![file("a.fq", Role::Single), file("b.fq", Role::Single), file("c.fq", Role::Single)]);
        assert_eq!(inputs(&files).unwrap()[2].label(), "File 3 of 3: c.fq");

        let swapped = set(Layout::Paired, vec![file("x_R2.fq", Role::R2), file("x_R1.fq", Role::R1)]);
        assert!(inputs(&swapped).is_err());
        assert!(inputs(&set(Layout::Paired, vec![file("x_R1.fq", Role::R1)])).is_err());
    }

    #[test]
    fn arguments_and_command_line() {
        let pair = AnalyzeInput { part: String::new(), r1: "C:/My reads/x_R1.fq.gz".into(), r2: Some("C:/My reads/x_R2.fq.gz".into()) };
        let a = analyze_args(&pair, false, 8);
        assert_eq!(strs(&a), ["analyze", "--r1", "C:/My reads/x_R1.fq.gz", "--r2", "C:/My reads/x_R2.fq.gz", "-w", "8"]);
        assert_eq!(command_line(&a), "simpletrim analyze --r1 \"C:/My reads/x_R1.fq.gz\" --r2 \"C:/My reads/x_R2.fq.gz\" -w 8");
        let single = AnalyzeInput { part: String::new(), r1: "x.fq".into(), r2: None };
        assert_eq!(strs(&analyze_args(&single, true, 0)), ["analyze", "--r1", "x.fq", "--max-reads", "0", "--max-bases", "0", "-w", "1"]);
    }

    #[test]
    fn stderr_lines() {
        assert_eq!(status_line("Reading up to 1,000,000 pairs ...").as_deref(), Some("Reading up to 1,000,000 pairs"));
        assert_eq!(status_line("read 2,000,000 pairs").as_deref(), Some("read 2,000,000 pairs"));
        assert_eq!(status_line("error: x contains no reads"), None);
        assert_eq!(status_line("  "), None);
        assert_eq!(reads_done("read 2,000,000 pairs"), Some(2_000_000));
        assert_eq!(reads_done("read 300,000 reads"), Some(300_000));
        assert_eq!(reads_done("Reading up to 1,000,000 pairs ..."), None);
        assert_eq!(reads_done("read 12 bases"), None);
        let out = Output {
            stderr: "Reading up to 1,000,000 pairs ...\nerror: R1 (bad.fastq): line 1: expected a FASTQ header".into(),
            code: Some(1),
            ..Default::default()
        };
        assert_eq!(error_text(&out), "R1 (bad.fastq): line 1: expected a FASTQ header");
        let out = Output { stderr: "something odd".into(), code: Some(3), ..Default::default() };
        assert_eq!(error_text(&out), "something odd");
    }

    /// Runs the bundled simpletrim (skipped when `assets/tools` has none):
    /// the text report on stdout, progress lines live, errors from `error:` lines.
    #[test]
    fn bundled_simpletrim_analyzes() {
        use super::super::run::{Handle, Line, command, run};
        use std::sync::Mutex;
        use std::sync::atomic::AtomicBool;
        let Ok(exe) = super::super::locate(&super::super::SIMPLETRIM, &Default::default()) else {
            eprintln!("simpletrim is not in assets/tools; skipped");
            return;
        };
        let dir = std::env::temp_dir().join(format!("pcr_studio_simpletrim_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut fastq = String::new();
        let mut x: u64 = 7;
        for i in 0..500 {
            let seq: String = (0..100)
                .map(|_| {
                    x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
                    b"ACGT"[(x >> 62) as usize] as char
                })
                .collect();
            fastq.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(100)));
        }
        let good = dir.join("reads.fastq");
        let bad = dir.join("bad.fastq");
        std::fs::write(&good, fastq).unwrap();
        std::fs::write(&bad, "not a fastq\n").unwrap();
        let analyze = |path: &Path| {
            let input = AnalyzeInput { part: String::new(), r1: path.to_path_buf(), r2: None };
            let mut cmd = command(&exe);
            cmd.args(analyze_args(&input, false, 2));
            let lines = Mutex::new(Vec::new());
            let on_line = |l: Line| {
                if let Line::Text(t) = l {
                    lines.lock().unwrap().extend(status_line(t));
                }
            };
            let out = run(cmd, &Handle::default(), &AtomicBool::new(false), &on_line).unwrap();
            (out, lines.into_inner().unwrap())
        };
        let (out, lines) = analyze(&good);
        assert_eq!(out.code, Some(0), "{}", out.stderr);
        assert!(out.stdout.starts_with("simpletrim analyze: 500 reads inspected"), "{}", out.stdout);
        assert!(lines.iter().any(|l| l == "Reading up to 1,000,000 reads"), "{lines:?}");
        let (out, _) = analyze(&bad);
        assert_eq!(out.code, Some(1));
        assert!(error_text(&out).contains("expected a FASTQ header"), "{}", error_text(&out));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn report_header() {
        let input = AnalyzeInput { part: "Lane 2 of 4".into(), r1: "C:/r/x_R1.fq.gz".into(), r2: Some("C:/r/x_R2.fq.gz".into()) };
        let r = Report {
            document: "Run 7",
            kind: "FASTQ · paired-end · 4 lanes",
            input: &input,
            all_reads: false,
            version: "0.0.5",
            path: Path::new("C:/tools/simpletrim.exe"),
            command: "simpletrim analyze …",
            started: "2026-11-11 10:00:00",
            finished: "2026-11-11 10:00:02",
            seconds: 2.0,
        };
        let t = r.text("simpletrim analyze: 1,000,000 read pairs inspected\n\n");
        assert!(t.starts_with("PCR STUDIO · SIMPLETRIM ANALYZE\n===="));
        assert!(t.contains("Document:  Run 7 (FASTQ · paired-end · 4 lanes)\n"));
        assert!(t.contains("Analyzed:  lane 2 of 4\nR1:        C:/r/x_R1.fq.gz\nR2:        C:/r/x_R2.fq.gz\n"));
        assert!(t.contains("Reads:     up to 1,000,000 pairs from the start"));
        assert!(t.contains("Finished:  2026-11-11 10:00:02 (2.0 s)\n\nsimpletrim analyze: 1,000,000 read pairs inspected\n"));
        assert!(t.ends_with("inspected\n"));
    }
}
