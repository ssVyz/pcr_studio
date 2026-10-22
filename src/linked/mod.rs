//! Linked files: library documents that point to sequence files on disk
//! (FASTQ read sets, FASTA files) instead of holding the sequences.
//!
//! PCR Studio never changes a linked file. Each file is remembered by path,
//! size, modification time and a fingerprint (hash of the size and the first
//! and last MiB). A moved or deleted file shows as missing, a modified one as
//! changed; Locate… links a moved file again after checking it is the same
//! file.

pub mod fastq;
pub mod group;
pub mod scan;

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use xxhash_rust::xxh3::Xxh3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Format {
    Fastq,
    Fasta,
}

impl Format {
    pub fn label(self) -> &'static str {
        match self {
            Format::Fastq => "FASTQ",
            Format::Fasta => "FASTA",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Layout {
    Single,
    /// R1/R2 files; the files alternate R1, R2 (one pair per lane).
    Paired,
}

/// Sequencing technology of a read set (picks QC and mapping presets later).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Technology {
    #[default]
    Unknown,
    Illumina,
    Nanopore,
}

impl Technology {
    pub const ALL: [Technology; 3] = [Technology::Unknown, Technology::Illumina, Technology::Nanopore];
}

impl std::fmt::Display for Technology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Technology::Unknown => "Not detected",
            Technology::Illumina => "Illumina (short reads)",
            Technology::Nanopore => "Nanopore (long reads)",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Role {
    Single,
    R1,
    R2,
}

impl Role {
    pub fn label(self) -> &'static str {
        match self {
            Role::Single => "File",
            Role::R1 => "R1",
            Role::R2 => "R2",
        }
    }
}

/// One linked file and what identifies it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LinkedFile {
    pub path: PathBuf,
    pub role: Role,
    pub size: u64,
    /// Modification time (seconds since 1970) when the file was last checked.
    pub modified: Option<u64>,
    /// Hash of the size and the first and last MiB ([`fingerprint`]).
    pub fingerprint: u64,
    /// Hash of the whole file as stored (compressed files: the compressed bytes).
    pub content_hash: u64,
    pub records: u64,
    pub bases: u64,
}

impl LinkedFile {
    pub fn file_name(&self) -> String {
        self.path.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| self.path.display().to_string())
    }
}

/// Totals over all files of a set, recorded when it is linked.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Summary {
    pub records: u64,
    /// Bases without gaps.
    pub bases: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub n50: u64,
    /// Longest record including gaps (alignment width for aligned FASTA).
    pub width: u64,
    /// Mean Phred quality of all bases (FASTQ).
    pub mean_quality: Option<f64>,
    /// G+C share of the A/C/G/T bases, percent.
    pub gc_percent: f64,
    /// Share of N and other ambiguity codes among all bases, percent.
    pub n_percent: f64,
    /// FASTA with gaps, or with several sequences of one length (an alignment).
    pub aligned: bool,
}

/// What a linked document points to.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LinkedSet {
    pub format: Format,
    pub layout: Layout,
    pub technology: Technology,
    pub files: Vec<LinkedFile>,
    pub summary: Summary,
}

impl LinkedSet {
    /// "FASTQ · paired-end", "FASTQ · 3 files", "FASTA".
    pub fn kind_label(&self) -> String {
        match (self.format, self.layout) {
            (Format::Fastq, Layout::Paired) if self.files.len() > 2 => format!("FASTQ · paired-end · {} lanes", self.files.len() / 2),
            (Format::Fastq, Layout::Paired) => "FASTQ · paired-end".into(),
            (Format::Fastq, Layout::Single) if self.files.len() > 1 => format!("FASTQ · {} files", self.files.len()),
            (Format::Fastq, Layout::Single) => "FASTQ".into(),
            (Format::Fasta, _) => "FASTA".into(),
        }
    }

    /// What one record is called.
    pub fn unit(&self) -> &'static str {
        match self.format {
            Format::Fastq => "reads",
            Format::Fasta => "sequences",
        }
    }
}

/// State of a linked file on disk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileState {
    Ok,
    /// Not found or not reachable (e.g. a disconnected drive).
    Missing,
    /// Found, but not the linked file any more.
    Changed,
}

impl FileState {
    pub fn label(self) -> &'static str {
        match self {
            FileState::Ok => "OK",
            FileState::Missing => "Missing",
            FileState::Changed => "Changed",
        }
    }
}

/// Number of missing and changed files of a set.
pub fn problems(states: &[FileState]) -> (usize, usize) {
    let missing = states.iter().filter(|&&s| s == FileState::Missing).count();
    let changed = states.iter().filter(|&&s| s == FileState::Changed).count();
    (missing, changed)
}

/// "1 file missing", "2 files changed", "1 file missing, 1 changed"; `None` if all are fine.
pub fn problem_label(states: &[FileState]) -> Option<String> {
    let files = |n: usize| if n == 1 { "1 file".to_string() } else { format!("{n} files") };
    match problems(states) {
        (0, 0) => None,
        (m, 0) => Some(format!("{} missing", files(m))),
        (0, c) => Some(format!("{} changed", files(c))),
        (m, c) => Some(format!("{} missing, {c} changed", files(m))),
    }
}

const PROBE: u64 = 1 << 20;

/// Hash of a file's size and its first and last MiB: cheap to compute, and
/// different for a different file (gzip files end with a checksum).
/// Returns (size, fingerprint).
pub fn fingerprint(path: &Path) -> std::io::Result<(u64, u64)> {
    let mut f = File::open(path)?;
    let size = f.metadata()?.len();
    let mut h = Xxh3::new();
    h.update(&size.to_le_bytes());
    let head = size.min(PROBE);
    let mut buf = vec![0u8; head as usize];
    f.read_exact(&mut buf)?;
    h.update(&buf);
    if size > PROBE {
        let tail = (size - PROBE).min(PROBE);
        f.seek(SeekFrom::End(-(tail as i64)))?;
        buf.resize(tail as usize, 0);
        f.read_exact(&mut buf)?;
        h.update(&buf);
    }
    Ok((size, h.digest()))
}

/// Modification time in seconds since 1970.
pub fn modified_secs(meta: &std::fs::Metadata) -> Option<u64> {
    meta.modified().ok()?.duration_since(std::time::UNIX_EPOCH).ok().map(|d| d.as_secs())
}

/// Checks a linked file. Size and modification time unchanged: OK without
/// reading it; a new modification time with the same size: the fingerprint
/// decides. The second value is the new modification time of a file that is
/// still the same (to remember, so the next check is quick again).
pub fn check(f: &LinkedFile) -> (FileState, Option<u64>) {
    let Ok(meta) = std::fs::metadata(&f.path) else { return (FileState::Missing, None) };
    if !meta.is_file() {
        return (FileState::Missing, None);
    }
    if meta.len() != f.size {
        return (FileState::Changed, None);
    }
    let modified = modified_secs(&meta);
    if modified == f.modified {
        return (FileState::Ok, None);
    }
    match fingerprint(&f.path) {
        Ok((_, fp)) if fp == f.fingerprint => (FileState::Ok, modified),
        Ok(_) => (FileState::Changed, None),
        Err(_) => (FileState::Missing, None),
    }
}

/// Checks all files of a set; remembers new modification times of unchanged
/// files in `set`. Returns the states and whether `set` was updated.
pub fn check_set(set: &mut LinkedSet) -> (Vec<FileState>, bool) {
    let mut updated = false;
    let states = set
        .files
        .iter_mut()
        .map(|f| {
            let (state, modified) = check(f);
            if let Some(m) = modified {
                f.modified = Some(m);
                updated = true;
            }
            state
        })
        .collect();
    (states, updated)
}

/// Whether `path` is the linked file `f` (same size and fingerprint).
/// Returns its modification time if so.
pub fn same_file(f: &LinkedFile, path: &Path) -> Result<Option<u64>, String> {
    let meta = std::fs::metadata(path).map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    if meta.len() != f.size {
        return Err(format!(
            "{} is not the linked file {}: it has a different size ({} instead of {} bytes).",
            path.display(),
            f.file_name(),
            meta.len(),
            f.size
        ));
    }
    let (_, fp) = fingerprint(path).map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    if fp != f.fingerprint {
        return Err(format!("{} is not the linked file {}: its content differs.", path.display(), f.file_name()));
    }
    Ok(modified_secs(&meta))
}

/// Links file `idx` of `set` to `path` after checking it is the same file.
/// Other missing or changed files of the set that were in the same folder are
/// looked for (by name) in the new folder and linked too if they match.
/// Returns how many files were linked.
pub fn relink(set: &mut LinkedSet, idx: usize, path: &Path) -> Result<usize, String> {
    let file = set.files.get(idx).ok_or("No such linked file")?;
    let modified = same_file(file, path)?;
    let old_dir = file.path.parent().map(Path::to_path_buf);
    let new_dir = path.parent().map(Path::to_path_buf);
    let f = &mut set.files[idx];
    f.path = path.to_path_buf();
    f.modified = modified;
    let mut n = 1;
    if let (Some(old_dir), Some(new_dir)) = (old_dir, new_dir) {
        n += relink_folder(set, &old_dir, &new_dir);
    }
    Ok(n)
}

/// Links the files of `set` that were in `old_dir` and are not OK to files of
/// the same name in `new_dir`, where those are the same files. Returns how many.
pub fn relink_folder(set: &mut LinkedSet, old_dir: &Path, new_dir: &Path) -> usize {
    let mut n = 0;
    for f in set.files.iter_mut() {
        if f.path.parent() != Some(old_dir) || check(f).0 == FileState::Ok {
            continue;
        }
        let Some(name) = f.path.file_name() else { continue };
        let candidate = new_dir.join(name);
        if let Ok(modified) = same_file(f, &candidate) {
            f.path = candidate;
            f.modified = modified;
            n += 1;
        }
    }
    n
}

/// `v` with three significant digits and the first fitting unit (steps of 1000).
fn scaled(v: u64, base: &str, units: [&str; 4]) -> String {
    if v < 1000 {
        return format!("{v} {base}");
    }
    let mut x = v as f64;
    let mut unit = "";
    for u in units {
        x /= 1000.0;
        unit = u;
        if x < 999.5 {
            break;
        }
    }
    match x {
        x if x < 9.995 => format!("{x:.2} {unit}"),
        x if x < 99.95 => format!("{x:.1} {unit}"),
        x => format!("{x:.0} {unit}"),
    }
}

/// File size: "830 bytes", "12.4 MB", "1.27 GB".
pub fn size_label(bytes: u64) -> String {
    scaled(bytes, "bytes", ["KB", "MB", "GB", "TB"])
}

/// Number of bases: "950 bp", "35.2 kb", "1.24 Gb".
pub fn bases_label(bases: u64) -> String {
    scaled(bases, "bp", ["kb", "Mb", "Gb", "Tb"])
}

/// File name without `.gz` and the sequence extension (`reads.fastq.gz` → `reads`).
pub fn stem(path: &Path) -> String {
    let name = path.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    let name = name.strip_suffix(".gz").or_else(|| name.strip_suffix(".GZ")).unwrap_or(&name);
    match name.rsplit_once('.') {
        Some((s, _)) if !s.is_empty() => s.to_string(),
        _ => name.to_string(),
    }
}

/// Reader that counts the bytes it passes on (progress) and optionally hashes them.
struct Tap<R> {
    inner: R,
    read: Arc<AtomicU64>,
    hash: Option<Arc<Mutex<Xxh3>>>,
}

impl<R: Read> Read for Tap<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.read.fetch_add(n as u64, Ordering::Relaxed);
        if let Some(h) = &self.hash
            && let Ok(mut h) = h.lock()
        {
            h.update(&buf[..n]);
        }
        Ok(n)
    }
}

const BUF: usize = 256 * 1024;

/// Opens a sequence file for reading; gzip files (also multi-member ones as
/// written by bgzip and sequencers) are recognised by their first bytes and
/// decompressed. The bytes read from disk are added to `read`, and hashed into
/// `hash` if given.
pub fn open(path: &Path, read: Arc<AtomicU64>, hash: Option<Arc<Mutex<Xxh3>>>) -> Result<Box<dyn BufRead + Send>, String> {
    let file = File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
    let mut raw = BufReader::with_capacity(BUF, Tap { inner: file, read, hash });
    let magic = raw.fill_buf().map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    if magic.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::with_capacity(BUF, flate2::bufread::MultiGzDecoder::new(raw))))
    } else {
        Ok(Box::new(raw))
    }
}

/// Opens a sequence file without progress or hashing.
pub fn open_plain(path: &Path) -> Result<Box<dyn BufRead + Send>, String> {
    open(path, Arc::new(AtomicU64::new(0)), None)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::io::Write;

    /// A fresh scratch folder for one test.
    pub fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("pcr_studio_linked_{}_{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    pub fn write_gz(path: &Path, content: &[u8]) {
        let mut e = flate2::write::GzEncoder::new(File::create(path).unwrap(), flate2::Compression::fast());
        e.write_all(content).unwrap();
        e.finish().unwrap();
    }

    fn linked(path: &Path) -> LinkedFile {
        let meta = std::fs::metadata(path).unwrap();
        let (size, fingerprint) = super::fingerprint(path).unwrap();
        LinkedFile {
            path: path.to_path_buf(),
            role: Role::Single,
            size,
            modified: modified_secs(&meta),
            fingerprint,
            content_hash: 0,
            records: 0,
            bases: 0,
        }
    }

    #[test]
    fn opens_plain_and_gzip() {
        let dir = scratch("open");
        let plain = dir.join("a.fastq");
        std::fs::write(&plain, b"@r\nACGT\n+\nIIII\n").unwrap();
        let gz = dir.join("a.fastq.gz");
        // Two gzip members, as bgzip writes them.
        let mut bytes = Vec::new();
        for part in [&b"@r\nAC"[..], &b"GT\n+\nIIII\n"[..]] {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            e.write_all(part).unwrap();
            bytes.extend(e.finish().unwrap());
        }
        std::fs::write(&gz, &bytes).unwrap();
        for p in [&plain, &gz] {
            let read = Arc::new(AtomicU64::new(0));
            let mut s = String::new();
            open(p, read.clone(), None).unwrap().read_to_string(&mut s).unwrap();
            assert_eq!(s, "@r\nACGT\n+\nIIII\n");
            assert_eq!(read.load(Ordering::Relaxed), std::fs::metadata(p).unwrap().len());
        }
        assert_eq!(stem(&gz), "a");
        assert_eq!(stem(Path::new("x/run.1.fq")), "run.1");
        assert_eq!(
            [size_label(830), size_label(12_400_000), size_label(1_274_000_000), bases_label(999_600)],
            ["830 bytes", "12.4 MB", "1.27 GB", "1.00 Mb"]
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn states_and_relinking() {
        let dir = scratch("relink");
        let a = dir.join("r_R1.fastq");
        let b = dir.join("r_R2.fastq");
        let big: Vec<u8> = (0..3 * PROBE as usize).map(|i| b"ACGT"[i % 4]).collect();
        std::fs::write(&a, &big).unwrap();
        std::fs::write(&b, b"second file").unwrap();
        let mut set = LinkedSet {
            format: Format::Fastq,
            layout: Layout::Paired,
            technology: Technology::Unknown,
            files: vec![linked(&a), linked(&b)],
            summary: Summary::default(),
        };
        assert_eq!(check_set(&mut set).0, [FileState::Ok, FileState::Ok]);
        // A different modification time alone: still the same file.
        set.files[1].modified = Some(1);
        let (states, updated) = check_set(&mut set);
        assert_eq!(states, [FileState::Ok, FileState::Ok]);
        assert!(updated && set.files[1].modified != Some(1));
        // Same size, other content: changed.
        set.files[1].modified = Some(1);
        std::fs::write(&b, b"Second file").unwrap();
        assert_eq!(check_set(&mut set).0[1], FileState::Changed);
        std::fs::write(&b, b"second file").unwrap();
        // Both move to another folder: locating one finds the other.
        let moved = dir.join("moved");
        std::fs::create_dir_all(&moved).unwrap();
        std::fs::rename(&a, moved.join("r_R1.fastq")).unwrap();
        std::fs::rename(&b, moved.join("r_R2.fastq")).unwrap();
        assert_eq!(check_set(&mut set).0, [FileState::Missing, FileState::Missing]);
        assert_eq!(problem_label(&check_set(&mut set).0).as_deref(), Some("2 files missing"));
        std::fs::write(dir.join("other.fastq"), b"something else").unwrap();
        assert!(relink(&mut set, 0, &dir.join("other.fastq")).is_err());
        assert_eq!(relink(&mut set, 0, &moved.join("r_R1.fastq")).unwrap(), 2);
        assert_eq!(check_set(&mut set).0, [FileState::Ok, FileState::Ok]);
        assert_eq!(set.files[1].path, moved.join("r_R2.fastq"));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
