//! Core data model: sequence rows, documents, annotations and per-column statistics.

use crate::seq::{self, GAP};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// One sequence placed in alignment coordinates.
///
/// A row covers the columns `start .. start + data.len()`. Inside that range a
/// `-` is a real gap; outside it the row has no coverage (drawn blank, ignored
/// by consensus and primer evaluation).
#[derive(Debug, Clone, Default)]
pub struct Row {
    pub name: String,
    pub description: String,
    pub start: usize,
    pub data: Vec<u8>,
    pub meta: Vec<(String, String)>,
}

impl Row {
    /// Creates a row from gapped data, turning leading/trailing gaps into
    /// "no coverage" so partial sequences are not counted as deletions.
    pub fn from_gapped(name: String, description: String, data: Vec<u8>) -> Row {
        let first = data.iter().position(|&b| b != GAP);
        let (start, data) = match first {
            None => (0, Vec::new()),
            Some(first) => {
                let last = data.iter().rposition(|&b| b != GAP).unwrap();
                (first, data[first..=last].to_vec())
            }
        };
        Row { name, description, start, data, meta: Vec::new() }
    }

    #[inline]
    pub fn end(&self) -> usize {
        self.start + self.data.len()
    }

    /// Character at an alignment column, `None` if the row does not cover it.
    #[inline]
    pub fn at(&self, col: usize) -> Option<u8> {
        if col >= self.start && col < self.end() { Some(self.data[col - self.start]) } else { None }
    }

    pub fn ungapped(&self) -> Vec<u8> {
        seq::ungap(&self.data)
    }

    pub fn ungapped_len(&self) -> usize {
        self.data.iter().filter(|&&b| b != GAP).count()
    }

    pub fn meta_value(&self, key: &str) -> Option<&str> {
        self.meta.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    pub fn set_meta(&mut self, key: &str, value: String) {
        if let Some(entry) = self.meta.iter_mut().find(|(k, _)| k == key) {
            entry.1 = value;
        } else {
            self.meta.push((key.to_string(), value));
        }
    }

    /// Row contents over `[c0, c1)`, blank (no coverage) columns as `' '`.
    pub fn slice_with_blanks(&self, c0: usize, c1: usize) -> Vec<u8> {
        (c0..c1).map(|c| self.at(c).unwrap_or(b' ')).collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DocKind {
    /// Plain sequence list (e.g. an unaligned FASTA import).
    Sequences,
    /// Multiple alignment (aligned FASTA import).
    Alignment,
    /// Result of map-to-reference.
    Contig,
}

impl DocKind {
    pub fn as_str(self) -> &'static str {
        match self {
            DocKind::Sequences => "sequences",
            DocKind::Alignment => "alignment",
            DocKind::Contig => "contig",
        }
    }

    pub fn parse(s: &str) -> DocKind {
        match s {
            "alignment" => DocKind::Alignment,
            "contig" => DocKind::Contig,
            _ => DocKind::Sequences,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            DocKind::Sequences => "Sequence list",
            DocKind::Alignment => "Alignment",
            DocKind::Contig => "Contig",
        }
    }
}

/// The sequence data of a document. Immutable once created; shared via `Arc`.
#[derive(Debug, Clone)]
pub struct Document {
    pub rows: Vec<Row>,
    pub width: usize,
}

impl Document {
    pub fn new(rows: Vec<Row>) -> Document {
        let width = rows.iter().map(|r| r.end()).max().unwrap_or(0);
        Document { rows, width }
    }

    /// Sorted list of all metadata keys present in any row.
    pub fn meta_keys(&self) -> Vec<String> {
        let mut keys: Vec<String> = Vec::new();
        for r in &self.rows {
            for (k, _) in &r.meta {
                if !keys.contains(k) {
                    keys.push(k.clone());
                }
            }
        }
        keys
    }
}

/// Summary information about a stored document.
#[derive(Debug, Clone)]
pub struct DocInfo {
    pub id: i64,
    pub name: String,
    pub kind: DocKind,
    pub n_rows: usize,
    pub width: usize,
    pub reference: Option<usize>,
    pub created: String,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnnotationKind {
    ForwardPrimer,
    ReversePrimer,
    Probe,
    Region,
}

impl AnnotationKind {
    pub const ALL: [AnnotationKind; 4] = [
        AnnotationKind::ForwardPrimer,
        AnnotationKind::ReversePrimer,
        AnnotationKind::Probe,
        AnnotationKind::Region,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            AnnotationKind::ForwardPrimer => "primer_fwd",
            AnnotationKind::ReversePrimer => "primer_rev",
            AnnotationKind::Probe => "probe",
            AnnotationKind::Region => "region",
        }
    }

    pub fn parse(s: &str) -> AnnotationKind {
        match s {
            "primer_fwd" => AnnotationKind::ForwardPrimer,
            "primer_rev" => AnnotationKind::ReversePrimer,
            "probe" => AnnotationKind::Probe,
            _ => AnnotationKind::Region,
        }
    }

    /// Direction of the arrow drawn for this annotation: +1 right, -1 left, 0 none.
    pub fn direction(self) -> i8 {
        match self {
            AnnotationKind::ForwardPrimer => 1,
            AnnotationKind::ReversePrimer => -1,
            AnnotationKind::Probe | AnnotationKind::Region => 0,
        }
    }
}

impl std::fmt::Display for AnnotationKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            AnnotationKind::ForwardPrimer => "Forward primer",
            AnnotationKind::ReversePrimer => "Reverse primer",
            AnnotationKind::Probe => "Probe",
            AnnotationKind::Region => "Region",
        })
    }
}

/// An annotation on alignment columns `[start, end)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    pub id: i64,
    pub name: String,
    pub kind: AnnotationKind,
    pub start: usize,
    pub end: usize,
    pub note: String,
}

/// Per-column residue counts over all rows except the reference.
/// Order: A, C, G, T, gap, other (N / ambiguity codes).
#[derive(Debug, Clone, Default)]
pub struct ColumnStats {
    pub counts: Vec<[u32; 6]>,
}

impl ColumnStats {
    pub fn compute(doc: &Document, exclude: Option<usize>) -> ColumnStats {
        let width = doc.width;
        let counts = doc
            .rows
            .par_iter()
            .enumerate()
            .fold(
                || vec![[0u32; 6]; width],
                |mut acc, (i, row)| {
                    if Some(i) != exclude {
                        for (off, &b) in row.data.iter().enumerate() {
                            acc[row.start + off][seq::stat_index(b)] += 1;
                        }
                    }
                    acc
                },
            )
            .reduce(
                || vec![[0u32; 6]; width],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(b.iter()) {
                        for k in 0..6 {
                            x[k] += y[k];
                        }
                    }
                    a
                },
            );
        ColumnStats { counts }
    }

    #[inline]
    pub fn coverage(&self, col: usize) -> u32 {
        self.counts[col].iter().sum()
    }

    /// Mean pairwise identity of the column (gaps count as a residue),
    /// `None` without at least one covering sequence.
    pub fn identity(&self, col: usize) -> Option<f32> {
        let c = &self.counts[col];
        let n: u64 = c.iter().map(|&x| x as u64).sum();
        match n {
            0 => None,
            1 => Some(1.0),
            _ => {
                let same: u64 = c[..5].iter().map(|&x| (x as u64) * (x as u64).saturating_sub(1)).sum();
                Some(same as f32 / (n * (n - 1)) as f32)
            }
        }
    }

    /// Share of covering sequences that carry the most common residue
    /// (gaps count as a residue, N/ambiguity codes are ignored).
    pub fn conservation(&self, col: usize) -> Option<f32> {
        let c = &self.counts[col];
        let n: u32 = c[..5].iter().sum();
        if n == 0 {
            return None;
        }
        Some(*c[..5].iter().max().unwrap() as f32 / n as f32)
    }

    pub fn conservation_track(&self) -> Vec<Option<f32>> {
        (0..self.counts.len()).map(|c| self.conservation(c)).collect()
    }

    pub fn identity_track(&self) -> Vec<Option<f32>> {
        (0..self.counts.len()).map(|c| self.identity(c)).collect()
    }

    /// Consensus sequence (one char per column) at the given threshold.
    pub fn consensus(&self, threshold: ConsensusThreshold, no_coverage: u8) -> Vec<u8> {
        self.counts.iter().map(|c| consensus_char(c, threshold, no_coverage)).collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusThreshold {
    Majority,
    Percent(u8),
}

impl ConsensusThreshold {
    pub const ALL: [ConsensusThreshold; 8] = [
        ConsensusThreshold::Majority,
        ConsensusThreshold::Percent(50),
        ConsensusThreshold::Percent(60),
        ConsensusThreshold::Percent(75),
        ConsensusThreshold::Percent(90),
        ConsensusThreshold::Percent(95),
        ConsensusThreshold::Percent(99),
        ConsensusThreshold::Percent(100),
    ];
}

impl std::fmt::Display for ConsensusThreshold {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConsensusThreshold::Majority => f.write_str("Highest (majority)"),
            ConsensusThreshold::Percent(p) => write!(f, "{p}% identity"),
        }
    }
}

/// Consensus rule: the smallest set of residues (by decreasing frequency)
/// whose combined share reaches the threshold is called as its IUPAC code.
/// A gap wins only when it is the most frequent state. `N`/ambiguity codes in
/// the reads are ignored for calling.
pub fn consensus_char(c: &[u32; 6], threshold: ConsensusThreshold, no_coverage: u8) -> u8 {
    let total: u32 = c[..5].iter().sum();
    if total == 0 {
        return if c[5] > 0 { b'N' } else { no_coverage };
    }
    let mut order: [usize; 5] = [0, 1, 2, 3, 4];
    order.sort_by(|&a, &b| c[b].cmp(&c[a]).then(a.cmp(&b)));
    if order[0] == 4 {
        return GAP;
    }
    let bits = [seq::A, seq::C, seq::G, seq::T];
    match threshold {
        ConsensusThreshold::Majority => {
            let top = c[order[0]];
            let m = (0..4).filter(|&i| c[i] == top).fold(0u8, |m, i| m | bits[i]);
            seq::code(m)
        }
        ConsensusThreshold::Percent(p) => {
            let need = (total as u64 * p as u64).div_ceil(100);
            let mut acc = 0u64;
            let mut m = 0u8;
            for &i in &order {
                if acc >= need || c[i] == 0 {
                    break;
                }
                acc += c[i] as u64;
                if i < 4 {
                    m |= bits[i];
                }
            }
            if m == 0 { GAP } else { seq::code(m) }
        }
    }
}

/// Maps alignment columns to ungapped reference coordinates (1-based).
#[derive(Debug, Clone, Default)]
pub struct RefCoords {
    /// For each column: 1-based reference position of the last reference base
    /// at or before this column (0 before the first base).
    pub pos_at: Vec<u32>,
    /// Whether the reference has a base (not gap / no coverage) in the column.
    pub is_base: Vec<bool>,
    /// Column of each reference base (index = ref pos - 1).
    pub col_of: Vec<u32>,
}

impl RefCoords {
    pub fn new(row: &Row, width: usize) -> RefCoords {
        let mut pos_at = Vec::with_capacity(width);
        let mut is_base = Vec::with_capacity(width);
        let mut col_of = Vec::new();
        let mut p = 0u32;
        for c in 0..width {
            let base = matches!(row.at(c), Some(b) if b != GAP);
            if base {
                p += 1;
                col_of.push(c as u32);
            }
            pos_at.push(p);
            is_base.push(base);
        }
        RefCoords { pos_at, is_base, col_of }
    }

    /// Column of a 1-based reference position (clamped).
    pub fn column_of(&self, pos: usize) -> Option<usize> {
        if self.col_of.is_empty() {
            return None;
        }
        let i = pos.clamp(1, self.col_of.len()) - 1;
        Some(self.col_of[i] as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc(rows: &[&str]) -> Document {
        Document::new(
            rows.iter()
                .enumerate()
                .map(|(i, s)| Row::from_gapped(format!("s{i}"), String::new(), s.as_bytes().to_vec()))
                .collect(),
        )
    }

    #[test]
    fn terminal_gaps_are_no_coverage() {
        let r = Row::from_gapped("x".into(), String::new(), b"--AC-G--".to_vec());
        assert_eq!(r.start, 2);
        assert_eq!(r.data, b"AC-G");
        assert_eq!(r.at(0), None);
        assert_eq!(r.at(4), Some(b'-'));
        assert_eq!(r.end(), 6);
    }

    #[test]
    fn stats_and_consensus() {
        let d = doc(&["ACGT", "ACGA", "ACCA", "A-CA"]);
        let s = ColumnStats::compute(&d, None);
        assert_eq!(s.counts[0], [4, 0, 0, 0, 0, 0]);
        let cons = s.consensus(ConsensusThreshold::Majority, b'?');
        assert_eq!(&cons, b"ACSA");
        let c100 = s.consensus(ConsensusThreshold::Percent(100), b'?');
        assert_eq!(&c100, b"ACSW");
        assert_eq!(s.identity(0), Some(1.0));
        assert_eq!(s.conservation(3), Some(0.75));
        assert!(s.identity(3).unwrap() < 1.0);
    }

    #[test]
    fn stats_exclude_reference() {
        let d = doc(&["TTTT", "AAAA"]);
        let s = ColumnStats::compute(&d, Some(0));
        assert_eq!(s.counts[0], [1, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn ref_coords() {
        let r = Row::from_gapped("r".into(), String::new(), b"AC--GT".to_vec());
        let rc = RefCoords::new(&r, 6);
        assert_eq!(rc.pos_at, vec![1, 2, 2, 2, 3, 4]);
        assert_eq!(rc.column_of(3), Some(4));
    }
}
