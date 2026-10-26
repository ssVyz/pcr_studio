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
use std::collections::{BTreeMap, HashSet};
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
    /// Phred qualities of the bases of each row (FASTQ; empty for FASTA).
    pub qual: Vec<Vec<u8>>,
}

/// The first records: FASTQ reads (pairs alternating R1, R2), or a FASTA file
/// completely up to a size limit.
pub fn preview(set: &LinkedSet) -> Result<Preview, String> {
    let (rows, qual) = match set.format {
        Format::Fastq => {
            let first: Vec<u64> = (0..sample_size(set).min(units(set)?)).collect();
            fastq_rows(set, &first, &AtomicBool::new(false), &|_| {})?
        }
        Format::Fasta => {
            let path = &set.files.first().ok_or("No file is linked")?.path;
            let reader = super::open_plain(path)?;
            let (rows, _, _) =
                fasta::read_rows(reader, PREVIEW_SEQUENCES, PREVIEW_BASES, &mut || {}).map_err(|e| format!("Invalid FASTA in {}: {e}", path.display()))?;
            (rows, Vec::new())
        }
    };
    let complete = rows.len() as u64 == set.summary.records;
    let kind = match set.format {
        Format::Fastq => DocKind::Sequences,
        Format::Fasta => fasta::alignment_kind(rows.iter().any(|r| r.data.contains(&seq::GAP)), &rows.iter().map(Row::end).collect::<Vec<_>>()),
    };
    Ok(Preview { doc: Document::new(rows), kind, complete, qual })
}

/// A random sample of the reads of a FASTQ set, as many as [`preview`] shows,
/// drawn evenly from all files and shown in file order. Reads the files (in
/// parallel) up to the last drawn read; `progress` gets the share of their
/// bytes read so far.
pub fn sample(set: &LinkedSet, seed: u64, cancel: &AtomicBool, progress: &(dyn Fn(f32) + Sync)) -> Result<Preview, String> {
    if set.format != Format::Fastq {
        return Err("Only FASTQ reads can be sampled".into());
    }
    let picks = draw(units(set)?, sample_size(set), seed);
    let (rows, qual) = fastq_rows(set, &picks, cancel, progress)?;
    let complete = rows.len() as u64 == set.summary.records;
    Ok(Preview { doc: Document::new(rows), kind: DocKind::Sequences, complete, qual })
}

/// Reads (pairs for paired sets) the viewer shows at once.
fn sample_size(set: &LinkedSet) -> u64 {
    match set.layout {
        Layout::Paired => PREVIEW_READS as u64 / 2,
        Layout::Single => PREVIEW_READS as u64,
    }
}

/// Where the reads of a FASTQ set come from: single files, or the R1 and R2
/// files of each lane.
fn sources(set: &LinkedSet) -> Result<Vec<&[LinkedFile]>, String> {
    match set.layout {
        Layout::Paired if set.files.len() < 2 || !set.files.len().is_multiple_of(2) => Err("A paired read set needs R1 and R2 files".into()),
        Layout::Paired => Ok(set.files.chunks(2).collect()),
        Layout::Single => Ok(set.files.chunks(1).collect()),
    }
}

/// Number of reads (pairs) of a FASTQ set.
fn units(set: &LinkedSet) -> Result<u64, String> {
    Ok(sources(set)?.iter().map(|files| files[0].records).sum())
}

/// SplitMix64: small and good enough to draw a sample.
struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A number below `n`.
    fn below(&mut self, n: u64) -> u64 {
        ((self.next() as u128 * n as u128) >> 64) as u64
    }
}

/// `k` different numbers below `n` (all of them if `n <= k`), each equally
/// likely (Floyd's algorithm), ascending.
fn draw(n: u64, k: u64, seed: u64) -> Vec<u64> {
    if n <= k {
        return (0..n).collect();
    }
    let mut rng = SplitMix(seed);
    let mut picked = HashSet::with_capacity(k as usize);
    for j in n - k..n {
        let t = rng.below(j + 1);
        if !picked.insert(t) {
            picked.insert(j);
        }
    }
    let mut picks: Vec<u64> = picked.into_iter().collect();
    picks.sort_unstable();
    picks
}

/// A read as a viewer row, with the Phred qualities (Phred+33) of its bases.
fn read_row(rec: &fastq::Record) -> (Row, Vec<u8>) {
    let mut data = Vec::with_capacity(rec.seq.len());
    let mut qual = Vec::with_capacity(rec.seq.len());
    for (&b, &q) in rec.seq.iter().zip(&rec.qual) {
        if let Some(b) = seq::normalize_byte(b) {
            data.push(b);
            qual.push(q.saturating_sub(33));
        }
    }
    let row = Row::from_gapped(String::from_utf8_lossy(rec.id()).to_string(), String::from_utf8_lossy(rec.desc()).to_string(), data);
    // Leading and trailing gaps are not part of the row.
    let qual = qual[row.start..row.start + row.data.len()].to_vec();
    (row, qual)
}

/// The reads with the given numbers (ascending; counted over all files, pairs
/// for paired sets) as rows with their qualities, pairs alternating R1, R2
/// (metadata "Mate"); metadata "File" with several files or lanes. The files
/// (lanes) are read in parallel, each up to its last wanted read; `progress`
/// gets the share of their bytes read so far.
fn fastq_rows(set: &LinkedSet, picks: &[u64], cancel: &AtomicBool, progress: &(dyn Fn(f32) + Sync)) -> Result<(Vec<Row>, Vec<Vec<u8>>), String> {
    let sources = sources(set)?;
    let several = sources.len() > 1;
    // The wanted reads of each source, counted within it.
    let mut jobs: Vec<(&[LinkedFile], Vec<u64>)> = Vec::new();
    let (mut offset, mut rest) = (0, picks);
    for files in sources {
        let n = files[0].records;
        let (mine, after) = rest.split_at(rest.partition_point(|&p| p < offset + n));
        if !mine.is_empty() {
            jobs.push((files, mine.iter().map(|p| p - offset).collect()));
        }
        (offset, rest) = (offset + n, after);
    }
    let total = jobs.iter().flat_map(|(files, _)| files.iter()).map(|f| f.size).sum::<u64>().max(1);
    let read = Arc::new(AtomicU64::new(0));
    let report = || progress(read.load(Ordering::Relaxed) as f32 / total as f32);
    let parts = jobs.par_iter().map(|(files, wanted)| read_source(files, wanted, several, &read, cancel, &report)).collect::<Result<Vec<_>, String>>()?;
    let (mut rows, mut qual) = (Vec::new(), Vec::new());
    for (r, q) in parts {
        rows.extend(r);
        qual.extend(q);
    }
    Ok((rows, qual))
}

/// The wanted reads (ascending numbers) of one file or R1/R2 pair of files.
fn read_source(
    files: &[LinkedFile],
    wanted: &[u64],
    several: bool,
    read: &Arc<AtomicU64>,
    cancel: &AtomicBool,
    report: &(dyn Fn() + Sync),
) -> Result<(Vec<Row>, Vec<Vec<u8>>), String> {
    let mut readers = files.iter().map(|f| super::open(&f.path, read.clone(), None).map(fastq::Reader::new)).collect::<Result<Vec<_>, String>>()?;
    let mut recs = vec![fastq::Record::default(); files.len()];
    let (mut rows, mut qual) = (Vec::new(), Vec::new());
    let (mut i, mut w) = (0u64, 0usize);
    'reads: while w < wanted.len() {
        for ((r, rec), f) in readers.iter_mut().zip(recs.iter_mut()).zip(files) {
            if !r.next(rec).map_err(|e| format!("{} is not a valid FASTQ file ({e})", f.path.display()))? {
                break 'reads;
            }
        }
        if i == wanted[w] {
            w += 1;
            for (mate, (rec, f)) in recs.iter().zip(files).enumerate() {
                let (mut row, q) = read_row(rec);
                if files.len() == 2 {
                    row.set_meta("Mate", ["R1", "R2"][mate].into());
                }
                if several {
                    row.set_meta("File", f.file_name());
                }
                rows.push(row);
                qual.push(q);
            }
        }
        i += 1;
        if i.is_multiple_of(REPORT_EVERY) {
            if cancel.load(Ordering::Relaxed) {
                return Err("Cancelled".into());
            }
            report();
        }
    }
    Ok((rows, qual))
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
    fn qualities_follow_the_bases() {
        // '1' is dropped, the leading and trailing gaps are not part of the row.
        let rec = fastq::Record { head: b"r x".to_vec(), seq: b".AC1G-".to_vec(), qual: b"!#+5?I".to_vec() };
        let (row, qual) = read_row(&rec);
        assert_eq!((row.start, row.data.as_slice(), qual.as_slice()), (1, &b"ACG"[..], &[2, 10, 30][..]));
        let rec = fastq::Record { head: b"e".to_vec(), seq: b"--".to_vec(), qual: b"II".to_vec() };
        assert!(read_row(&rec).1.is_empty());
    }

    #[test]
    fn draws_distinct_numbers_evenly() {
        assert_eq!(draw(3, 5, 1), [0, 1, 2]);
        let a = draw(10_000, 1000, 7);
        assert_eq!(a.len(), 1000);
        assert!(a.windows(2).all(|w| w[0] < w[1]) && *a.last().unwrap() < 10_000);
        assert_eq!(a, draw(10_000, 1000, 7));
        assert_ne!(a, draw(10_000, 1000, 8));
        // About half of them in each half.
        let low = a.iter().filter(|&&x| x < 5000).count();
        assert!((400..600).contains(&low), "{low}");
    }

    #[test]
    fn samples_reads_from_all_lanes() {
        let dir = scratch("sample");
        let mut files = Vec::new();
        let mut paths = Vec::new();
        for (lane, n) in [(1, 600), (2, 700)] {
            let names: Vec<String> = (0..n).map(|i| format!("L{lane}:{i}")).collect();
            let refs: Vec<&str> = names.iter().map(String::as_str).collect();
            for (role, seq) in [(Role::R1, "ACGT"), (Role::R2, "GGCCA")] {
                let p = dir.join(format!("s_L00{lane}_{}.fastq", role.label()));
                std::fs::write(&p, fastq(&refs, &vec![seq; n])).unwrap();
                files.push((p.clone(), role));
                paths.push((p, Format::Fastq));
            }
        }
        let set = build_set(Format::Fastq, Layout::Paired, &files, scan(&paths)).unwrap();
        let cancel = AtomicBool::new(false);
        let p = sample(&set, 42, &cancel, &|_| {}).unwrap();
        let rows = &p.doc.rows;
        assert!(!p.complete && rows.len() == PREVIEW_READS && p.qual.len() == rows.len());
        // Mates stay together, in file order, from both lanes.
        for pair in rows.chunks(2) {
            assert_eq!((pair[0].meta_value("Mate"), pair[1].meta_value("Mate")), (Some("R1"), Some("R2")));
            assert_eq!(pair[0].name, pair[1].name);
            assert_eq!(pair[1].data, b"GGCCA");
        }
        assert!(rows.iter().all(|r| r.meta_value("File").is_some()));
        assert!(rows.iter().any(|r| r.name.starts_with("L1:")) && rows.iter().any(|r| r.name.starts_with("L2:")));
        assert!(p.qual.iter().zip(rows).all(|(q, r)| q.len() == r.data.len() && q.iter().all(|&x| x == 40)));
        let other = sample(&set, 43, &cancel, &|_| {}).unwrap();
        assert_ne!(rows.iter().map(|r| &r.name).collect::<Vec<_>>(), other.doc.rows.iter().map(|r| &r.name).collect::<Vec<_>>());
        // The first reads: all from the first lane.
        let first = preview(&set).unwrap();
        assert_eq!((first.doc.rows[0].name.as_str(), first.doc.rows[999].name.as_str()), ("L1:0", "L1:499"));
        // Single files: the whole set when it is small, and cancelling.
        let one = &files[..1];
        let set = build_set(Format::Fastq, Layout::Single, &[(one[0].0.clone(), Role::Single)], scan(&paths[..1])).unwrap();
        let p = sample(&set, 1, &cancel, &|_| {}).unwrap();
        assert!(p.complete && p.doc.rows.len() == 600 && p.doc.rows[0].meta.is_empty());
        let many = dir.join("many.fastq");
        let names: Vec<String> = (0..10_000).map(|i| format!("r{i}")).collect();
        let refs: Vec<&str> = names.iter().map(String::as_str).collect();
        std::fs::write(&many, fastq(&refs, &vec!["ACGT"; 10_000])).unwrap();
        let set = build_set(Format::Fastq, Layout::Single, &[(many.clone(), Role::Single)], scan(&[(many, Format::Fastq)])).unwrap();
        cancel.store(true, Ordering::Relaxed);
        assert_eq!(sample(&set, 1, &cancel, &|_| {}).err().as_deref(), Some("Cancelled"));
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
