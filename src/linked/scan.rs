//! Reading linked files: one pass over every file when it is linked (totals,
//! length distribution, qualities, checks that R1 and R2 belong together),
//! and the first records for the viewer.

use super::fastq;
use super::{Format, Layout, LinkedFile, LinkedSet, Role, Summary, Technology};
use crate::fasta;
use crate::model::{DocKind, Document, Row};
use crate::seq;
use bio::io::fasta::FastaRead;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use xxhash_rust::xxh3::Xxh3;

/// What one pass over a file found.
#[derive(Debug, Clone, Default)]
pub struct FileScan {
    pub size: u64,
    pub modified: Option<u64>,
    pub fingerprint: u64,
    pub content_hash: u64,
    pub records: u64,
    /// Number of records per length (without gaps).
    pub lengths: BTreeMap<u64, u64>,
    /// Shortest and longest record including gaps.
    pub min_width: u64,
    pub max_width: u64,
    /// A, C, G, T, other bases (N, ambiguity codes).
    pub counts: [u64; 5],
    pub gapped: bool,
    /// Sum of the Phred qualities of all bases (FASTQ).
    pub quality_sum: u64,
    /// Hash of all read names in order, without `/1` and `/2` (equal for the two files of a pair).
    pub names: u64,
    /// Headers of the first records (technology detection).
    pub heads: Vec<String>,
}

impl FileScan {
    pub fn bases(&self) -> u64 {
        self.counts.iter().sum()
    }

    fn add_record(&mut self, seq: &[u8]) {
        let mut len = 0u64;
        for &b in seq {
            match CLASS[b as usize] {
                c @ 0..=4 => {
                    self.counts[c as usize] += 1;
                    len += 1;
                }
                GAP_CLASS => self.gapped = true,
                _ => {}
            }
        }
        *self.lengths.entry(len).or_default() += 1;
        let width = seq.len() as u64;
        self.min_width = if self.records == 0 { width } else { self.min_width.min(width) };
        self.max_width = self.max_width.max(width);
        self.records += 1;
    }
}

const GAP_CLASS: u8 = 5;
const IGNORED: u8 = 6;

/// Base class of every byte: 0-3 A/C/G/T (and U), 4 other residues, 5 gaps, 6 ignored.
static CLASS: [u8; 256] = {
    let mut t = [4u8; 256];
    let mut i = 0;
    while i < 256 {
        let b = (i as u8).to_ascii_uppercase();
        t[i] = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            b'-' | b'.' | b'~' => GAP_CLASS,
            b' ' | b'\t' | b'\r' | b'\n' | b'*' | b'0'..=b'9' => IGNORED,
            _ => 4,
        };
        i += 1;
    }
    t
};

/// Read name without a trailing `/1` or `/2` (old Illumina mate suffixes).
fn mate_name(id: &[u8]) -> &[u8] {
    match id {
        [rest @ .., b'/', b'1' | b'2'] => rest,
        _ => id,
    }
}

const HEADS: usize = 50;
/// Records between progress reports and cancel checks.
const REPORT_EVERY: u64 = 4096;

/// Reads a whole file once. `read` accumulates the bytes read from disk (all
/// files of a run share it); `progress` gets its value now and then.
pub fn scan_file(path: &Path, format: Format, read: &Arc<AtomicU64>, cancel: &AtomicBool, progress: &(dyn Fn(u64) + Sync)) -> Result<FileScan, String> {
    let meta = std::fs::metadata(path).map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    let (size, fingerprint) = super::fingerprint(path).map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    let hash = Arc::new(Mutex::new(Xxh3::new()));
    let mut reader = super::open(path, read.clone(), Some(hash.clone()))?;
    let first = reader.fill_buf().map_err(|e| format!("Cannot read {}: {e}", path.display()))?.iter().find(|b| !b.is_ascii_whitespace()).copied();
    match (format, first) {
        (Format::Fasta, Some(b'@')) => return Err(format!("{} is a FASTQ file: use File → Import FASTQ (linked)…", path.display())),
        (Format::Fastq, Some(b'>')) => return Err(format!("{} is a FASTA file: use File → Import FASTA…", path.display())),
        _ => {}
    }
    let mut s = FileScan { size, modified: super::modified_secs(&meta), fingerprint, ..Default::default() };
    let cancelled = || Err::<FileScan, String>("Cancelled".into());
    match format {
        Format::Fastq => {
            let mut r = fastq::Reader::new(reader);
            let mut rec = fastq::Record::default();
            let mut names = Xxh3::new();
            while r.next(&mut rec).map_err(|e| format!("{} is not a valid FASTQ file ({e})", path.display()))? {
                s.add_record(&rec.seq);
                s.quality_sum += rec.qual.iter().map(|&q| q.saturating_sub(33) as u64).sum::<u64>();
                names.update(mate_name(rec.id()));
                names.update(b"\n");
                if s.heads.len() < HEADS {
                    s.heads.push(String::from_utf8_lossy(&rec.head).to_string());
                }
                if s.records.is_multiple_of(REPORT_EVERY) {
                    if cancel.load(Ordering::Relaxed) {
                        return cancelled();
                    }
                    progress(read.load(Ordering::Relaxed));
                }
            }
            s.names = names.digest();
        }
        Format::Fasta => {
            let mut r = bio::io::fasta::Reader::from_bufread(reader);
            let mut rec = bio::io::fasta::Record::new();
            loop {
                r.read(&mut rec).map_err(|e| format!("{} is not a valid FASTA file ({e})", path.display()))?;
                if rec.is_empty() {
                    break;
                }
                s.add_record(rec.seq());
                if s.heads.len() < HEADS {
                    s.heads.push(rec.id().to_string());
                }
                if s.records.is_multiple_of(REPORT_EVERY) {
                    if cancel.load(Ordering::Relaxed) {
                        return cancelled();
                    }
                    progress(read.load(Ordering::Relaxed));
                }
            }
        }
    }
    if s.records == 0 {
        return Err(format!("{} contains no {} records", path.display(), format.label()));
    }
    s.content_hash = hash.lock().map(|h| h.digest()).unwrap_or(0);
    Ok(s)
}

/// Scans several files in parallel; results in the order of `files`.
pub fn scan_files(files: &[(PathBuf, Format)], cancel: &AtomicBool, progress: &(dyn Fn(u64) + Sync)) -> Vec<Result<FileScan, String>> {
    let read = Arc::new(AtomicU64::new(0));
    files.par_iter().map(|(p, f)| scan_file(p, *f, &read, cancel, progress)).collect()
}

/// N50: the length such that records at least this long hold half of all bases.
fn n50(lengths: &BTreeMap<u64, u64>) -> u64 {
    let total: u64 = lengths.iter().map(|(l, n)| l * n).sum();
    let mut acc = 0;
    for (&len, &n) in lengths.iter().rev() {
        acc += len * n;
        if acc * 2 >= total {
            return len;
        }
    }
    0
}

/// Totals over the scans of a set.
pub fn summarize(format: Format, scans: &[FileScan]) -> Summary {
    let mut lengths: BTreeMap<u64, u64> = BTreeMap::new();
    let mut counts = [0u64; 5];
    let mut quality = 0u64;
    for s in scans {
        for (&l, &n) in &s.lengths {
            *lengths.entry(l).or_default() += n;
        }
        for (c, n) in counts.iter_mut().zip(s.counts) {
            *c += n;
        }
        quality += s.quality_sum;
    }
    let records: u64 = scans.iter().map(|s| s.records).sum();
    let bases: u64 = counts.iter().sum();
    let acgt: u64 = counts[..4].iter().sum();
    let pct = |n: u64, of: u64| if of == 0 { 0.0 } else { 100.0 * n as f64 / of as f64 };
    let min_width = scans.iter().map(|s| s.min_width).min().unwrap_or(0);
    let width = scans.iter().map(|s| s.max_width).max().unwrap_or(0);
    let gapped = scans.iter().any(|s| s.gapped);
    Summary {
        records,
        bases,
        min_len: lengths.keys().next().copied().unwrap_or(0),
        max_len: lengths.keys().next_back().copied().unwrap_or(0),
        n50: n50(&lengths),
        width,
        mean_quality: (format == Format::Fastq && bases > 0).then(|| quality as f64 / bases as f64),
        gc_percent: pct(counts[1] + counts[2], acgt),
        n_percent: pct(counts[4], bases),
        aligned: format == Format::Fasta && (gapped || (records > 1 && min_width == width)),
    }
}

/// Whether `s` looks like a UUID (Nanopore read names).
fn is_uuid(s: &str) -> bool {
    s.len() == 36 && s.char_indices().all(|(i, c)| if matches!(i, 8 | 13 | 18 | 23) { c == '-' } else { c.is_ascii_hexdigit() })
}

/// Sequencing technology from the read headers, else from the read lengths.
pub fn detect_technology(heads: &[String], summary: &Summary) -> Technology {
    let id = |h: &String| h.split_ascii_whitespace().next().unwrap_or("").to_string();
    if heads.iter().any(|h| h.contains("runid=") || h.contains("start_time=") || is_uuid(&id(h))) {
        return Technology::Nanopore;
    }
    // Illumina: instrument:run:flowcell:lane:tile:x:y (older: 4 colons).
    if !heads.is_empty() && heads.iter().all(|h| id(h).matches(':').count() >= 4) {
        return Technology::Illumina;
    }
    match summary {
        s if s.records == 0 => Technology::Unknown,
        s if s.max_len <= 600 => Technology::Illumina,
        s if s.n50 >= 1000 => Technology::Nanopore,
        _ => Technology::Unknown,
    }
}

/// Builds a linked set from its files and their scans; paired files must have
/// the same reads (count and names) in R1 and R2.
pub fn build_set(format: Format, layout: Layout, files: &[(PathBuf, Role)], scans: Vec<FileScan>) -> Result<LinkedSet, String> {
    if layout == Layout::Paired {
        for (pair, s) in files.chunks(2).zip(scans.chunks(2)) {
            let [(p1, _), (p2, _)] = pair else { return Err("Paired files must come in R1/R2 pairs".into()) };
            let name = |p: &PathBuf| p.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
            if s[0].records != s[1].records {
                return Err(format!(
                    "{} and {} have different numbers of reads ({} and {}): they are not the two files of one paired run.",
                    name(p1),
                    name(p2),
                    s[0].records,
                    s[1].records
                ));
            }
            if s[0].names != s[1].names {
                return Err(format!(
                    "The read names in {} and {} do not match: they are not the two files of one paired run.",
                    name(p1),
                    name(p2)
                ));
            }
        }
    }
    let summary = summarize(format, &scans);
    let technology = match format {
        Format::Fastq => detect_technology(&scans[0].heads, &summary),
        Format::Fasta => Technology::Unknown,
    };
    let files = files
        .iter()
        .zip(&scans)
        .map(|((path, role), s)| LinkedFile {
            path: path.clone(),
            role: *role,
            size: s.size,
            modified: s.modified,
            fingerprint: s.fingerprint,
            content_hash: s.content_hash,
            records: s.records,
            bases: s.bases(),
        })
        .collect();
    Ok(LinkedSet { format, layout, technology, files, summary })
}

/// Reads shown in the viewer for a FASTQ set (pairs: half as many pairs).
pub const PREVIEW_READS: usize = 1000;
/// FASTA files are shown completely up to these limits.
const PREVIEW_SEQUENCES: usize = 100_000;
const PREVIEW_BASES: u64 = 100_000_000;

/// The records of a linked set the viewer shows.
pub struct Preview {
    pub doc: Document,
    pub kind: DocKind,
    /// All records of the set are in `doc`.
    pub complete: bool,
}

pub fn preview(set: &LinkedSet) -> Result<Preview, String> {
    let rows = match set.format {
        Format::Fastq => preview_fastq(set)?,
        Format::Fasta => {
            let path = &set.files.first().ok_or("No file is linked")?.path;
            let reader = super::open_plain(path)?;
            let (rows, _, _) =
                fasta::read_rows(reader, PREVIEW_SEQUENCES, PREVIEW_BASES, &mut || {}).map_err(|e| format!("Invalid FASTA in {}: {e}", path.display()))?;
            rows
        }
    };
    let complete = rows.len() as u64 == set.summary.records;
    let kind = match set.format {
        Format::Fastq => DocKind::Sequences,
        Format::Fasta => fasta::alignment_kind(rows.iter().any(|r| r.data.contains(&seq::GAP)), &rows.iter().map(Row::end).collect::<Vec<_>>()),
    };
    Ok(Preview { doc: Document::new(rows), kind, complete })
}

fn read_row(rec: &fastq::Record) -> Row {
    Row::from_gapped(
        String::from_utf8_lossy(rec.id()).to_string(),
        String::from_utf8_lossy(rec.desc()).to_string(),
        seq::normalize(&rec.seq),
    )
}

/// First reads: of the first R1/R2 pair alternating (metadata "Mate"), else
/// from the files in order (metadata "File" with several files).
fn preview_fastq(set: &LinkedSet) -> Result<Vec<Row>, String> {
    let mut rows = Vec::new();
    let invalid = |p: &Path, e: String| format!("{} is not a valid FASTQ file ({e})", p.display());
    if set.layout == Layout::Paired {
        let [f1, f2, ..] = set.files.as_slice() else { return Err("A paired read set needs R1 and R2".into()) };
        let mut r1 = fastq::Reader::new(super::open_plain(&f1.path)?);
        let mut r2 = fastq::Reader::new(super::open_plain(&f2.path)?);
        let (mut a, mut b) = (fastq::Record::default(), fastq::Record::default());
        while rows.len() < PREVIEW_READS && r1.next(&mut a).map_err(|e| invalid(&f1.path, e))? && r2.next(&mut b).map_err(|e| invalid(&f2.path, e))? {
            for (rec, mate) in [(&a, "R1"), (&b, "R2")] {
                let mut row = read_row(rec);
                row.set_meta("Mate", mate.into());
                rows.push(row);
            }
        }
        return Ok(rows);
    }
    let mut rec = fastq::Record::default();
    for f in &set.files {
        let mut r = fastq::Reader::new(super::open_plain(&f.path)?);
        while rows.len() < PREVIEW_READS && r.next(&mut rec).map_err(|e| invalid(&f.path, e))? {
            let mut row = read_row(&rec);
            if set.files.len() > 1 {
                row.set_meta("File", f.file_name());
            }
            rows.push(row);
        }
        if rows.len() >= PREVIEW_READS {
            break;
        }
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::super::tests::{scratch, write_gz};
    use super::*;

    fn fastq(names: &[&str], seqs: &[&str]) -> String {
        names.iter().zip(seqs).map(|(n, s)| format!("@{n}\n{s}\n+\n{}\n", "I".repeat(s.len()))).collect()
    }

    fn scan(files: &[(PathBuf, Format)]) -> Vec<FileScan> {
        let cancel = AtomicBool::new(false);
        scan_files(files, &cancel, &|_| {}).into_iter().map(Result::unwrap).collect()
    }

    #[test]
    fn scans_paired_fastq() {
        let dir = scratch("scan_pairs");
        let r1 = dir.join("s_R1.fastq.gz");
        let r2 = dir.join("s_R2.fastq");
        write_gz(&r1, fastq(&["M1:1:FC:1:1:1:1 1:N:0:1", "M1:1:FC:1:1:1:2 1:N:0:1"], &["ACGTN", "GGCC"]).as_bytes());
        std::fs::write(&r2, fastq(&["M1:1:FC:1:1:1:1 2:N:0:1", "M1:1:FC:1:1:1:2 2:N:0:1"], &["TTTT", "CCGGA"])).unwrap();
        let files = vec![(r1.clone(), Role::R1), (r2.clone(), Role::R2)];
        let scans = scan(&[(r1.clone(), Format::Fastq), (r2.clone(), Format::Fastq)]);
        assert_eq!((scans[0].records, scans[0].bases(), scans[0].counts[4]), (2, 9, 1));
        assert_eq!(scans[0].size, std::fs::metadata(&r1).unwrap().len());
        let set = build_set(Format::Fastq, Layout::Paired, &files, scans).unwrap();
        let s = &set.summary;
        assert_eq!((s.records, s.bases, s.min_len, s.max_len, s.width), (4, 18, 4, 5, 5));
        assert_eq!(s.mean_quality, Some(40.0));
        assert!((s.gc_percent - 100.0 * 10.0 / 17.0).abs() < 1e-9);
        assert_eq!(set.technology, Technology::Illumina);
        assert_eq!(set.files[1].records, 2);
        // The preview alternates R1 and R2.
        let p = preview(&set).unwrap();
        assert!(p.complete && p.kind == DocKind::Sequences);
        let mates: Vec<&str> = p.doc.rows.iter().map(|r| r.meta_value("Mate").unwrap()).collect();
        assert_eq!(mates, ["R1", "R2", "R1", "R2"]);
        assert_eq!(p.doc.rows[1].data, b"TTTT");
        // R2 from another run: different names.
        std::fs::write(&r2, fastq(&["x1", "x2"], &["TTTT", "CCGG"])).unwrap();
        let scans = scan(&[(r1.clone(), Format::Fastq), (r2.clone(), Format::Fastq)]);
        assert!(build_set(Format::Fastq, Layout::Paired, &files, scans).unwrap_err().contains("names"));
        std::fs::write(&r2, fastq(&["M1:1:FC:1:1:1:1"], &["TTTT"])).unwrap();
        let scans = scan(&[(r1.clone(), Format::Fastq), (r2.clone(), Format::Fastq)]);
        assert!(build_set(Format::Fastq, Layout::Paired, &files, scans).unwrap_err().contains("numbers of reads"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scans_fasta_and_rejects_wrong_formats() {
        let dir = scratch("scan_fasta");
        let fa = dir.join("aln.fasta");
        std::fs::write(&fa, ">a\nAC-GT\n>b\nACCGT\n").unwrap();
        let scans = scan(&[(fa.clone(), Format::Fasta)]);
        let set = build_set(Format::Fasta, Layout::Single, &[(fa.clone(), Role::Single)], scans).unwrap();
        let s = &set.summary;
        assert_eq!((s.records, s.bases, s.min_len, s.max_len, s.width, s.aligned, s.mean_quality), (2, 9, 4, 5, 5, true, None));
        let p = preview(&set).unwrap();
        assert!(p.complete && p.kind == DocKind::Alignment);
        let cancel = AtomicBool::new(false);
        let as_fastq = scan_files(&[(fa.clone(), Format::Fastq)], &cancel, &|_| {});
        assert!(as_fastq[0].as_ref().unwrap_err().contains("is a FASTA file"));
        let broken = dir.join("broken.fastq");
        std::fs::write(&broken, "@r\nACGT\n+\nII\n").unwrap();
        assert!(scan_files(&[(broken, Format::Fastq)], &cancel, &|_| {})[0].as_ref().unwrap_err().contains("not a valid FASTQ"));
        let empty = dir.join("empty.fastq");
        std::fs::write(&empty, "").unwrap();
        assert!(scan_files(&[(empty, Format::Fastq)], &cancel, &|_| {})[0].as_ref().unwrap_err().contains("no FASTQ records"));
        cancel.store(true, Ordering::Relaxed);
        let many = dir.join("many.fastq");
        let names: Vec<String> = (0..5000).map(|i| format!("r{i}")).collect();
        let refs: Vec<&str> = names.iter().map(String::as_str).collect();
        std::fs::write(&many, fastq(&refs, &vec!["ACGT"; 5000])).unwrap();
        assert_eq!(scan_files(&[(many, Format::Fastq)], &cancel, &|_| {})[0].as_ref().unwrap_err(), "Cancelled");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn technology_and_n50() {
        let summary = |max_len, n50| Summary { records: 10, max_len, n50, ..Default::default() };
        let heads = |h: &[&str]| h.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        assert_eq!(
            detect_technology(&heads(&["0a1b2c3d-1234-4abc-8def-0123456789ab runid=abc ch=12"]), &summary(5000, 3000)),
            Technology::Nanopore
        );
        assert_eq!(detect_technology(&heads(&["A00123:8:H5:2:1101:1000:2000 1:N:0:ACGT"]), &summary(151, 151)), Technology::Illumina);
        assert_eq!(detect_technology(&heads(&["SRR1.1 length=150"]), &summary(150, 150)), Technology::Illumina);
        assert_eq!(detect_technology(&heads(&["read7"]), &summary(40_000, 8000)), Technology::Nanopore);
        assert_eq!(detect_technology(&heads(&["read7"]), &summary(900, 700)), Technology::Unknown);
        let lengths: BTreeMap<u64, u64> = [(2, 3), (10, 1), (5, 2)].into_iter().collect();
        // 26 bases; 10 + 5 + 5 = 20 ≥ 13.
        assert_eq!(n50(&lengths), 5);
        assert_eq!(mate_name(b"r1/2"), b"r1");
    }
}
