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
    /// Plain text (e.g. the output of an external tool); no sequences.
    Text,
}

impl DocKind {
    pub fn as_str(self) -> &'static str {
        match self {
            DocKind::Sequences => "sequences",
            DocKind::Alignment => "alignment",
            DocKind::Contig => "contig",
            DocKind::Text => "text",
        }
    }

    pub fn parse(s: &str) -> DocKind {
        match s {
            "alignment" => DocKind::Alignment,
            "contig" => DocKind::Contig,
            "text" => DocKind::Text,
            _ => DocKind::Sequences,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            DocKind::Sequences => "Sequence list",
            DocKind::Alignment => "Alignment",
            DocKind::Contig => "Contig",
            DocKind::Text => "Text",
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

    /// Columns `[c0, c1)` of all rows except `exclude` (the reference), keeping
    /// names, descriptions and metadata. With `only` (one flag per alignment
    /// column) just the flagged columns are taken. Rows without any base in the
    /// slice and columns that are empty in every remaining row are dropped;
    /// terminal gaps become "no coverage" as on import.
    pub fn slice_columns(&self, c0: usize, c1: usize, exclude: Option<usize>, only: Option<&[bool]>) -> Document {
        let rows: Vec<usize> = (0..self.rows.len()).filter(|&i| Some(i) != exclude).collect();
        self.slice_rows(&rows, c0, c1, only).0
    }

    /// Like `slice_columns`, for the rows `rows` (in that order). Also returns
    /// the alignment column of every column of the slice.
    pub fn slice_rows(&self, rows: &[usize], c0: usize, c1: usize, only: Option<&[bool]>) -> (Document, Vec<usize>) {
        let c1 = c1.min(self.width);
        let c0 = c0.min(c1);
        let picked: Vec<(&Row, Vec<u8>)> = rows
            .iter()
            .filter_map(|&i| self.rows.get(i))
            .map(|r| (r, pick_columns(r.slice_with_blanks(c0, c1), c0, only)))
            .filter(|(_, d)| d.iter().any(|&b| b != GAP && b != b' '))
            .collect();
        let n = picked.first().map_or(0, |(_, d)| d.len());
        let keep: Vec<bool> = (0..n).map(|i| picked.iter().any(|(_, d)| d[i] != GAP && d[i] != b' ')).collect();
        let columns: Vec<usize> = (c0..c1)
            .filter(|&c| only.is_none_or(|m| m.get(c).copied().unwrap_or(false)))
            .zip(&keep)
            .filter(|(_, k)| **k)
            .map(|(c, _)| c)
            .collect();
        let out = picked
            .into_iter()
            .map(|(r, d)| {
                let data: Vec<u8> = d.iter().zip(&keep).filter(|(_, k)| **k).map(|(&b, _)| if b == b' ' { GAP } else { b }).collect();
                let mut row = Row::from_gapped(r.name.clone(), r.description.clone(), data);
                row.meta = r.meta.clone();
                row
            })
            .collect();
        (Document::new(out), columns)
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
    /// Library folder, `None` = top level.
    pub folder: Option<i64>,
}

/// A library folder.
#[derive(Debug, Clone, PartialEq)]
pub struct Folder {
    pub id: i64,
    pub name: String,
    pub parent: Option<i64>,
    pub collapsed: bool,
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

    /// Most frequent state of the column (A/C/G/T or gap; N ignored),
    /// `None` without coverage.
    pub fn majority_state(&self, col: usize) -> Option<u8> {
        let c = &self.counts[col];
        let (i, &n) = c[..5].iter().enumerate().max_by_key(|&(i, &n)| (n, 5 - i))?;
        (n > 0).then_some(seq::STAT_CHARS[i])
    }

    /// Number of covering sequences (N/ambiguity codes excluded) and how many
    /// of them carry a base compatible with `cmp`.
    fn matching(&self, col: usize, cmp: Option<u8>) -> (u64, u64) {
        let c = &self.counts[col];
        let n: u64 = c[..5].iter().map(|&x| x as u64).sum();
        let mask = cmp.map(seq::mask).unwrap_or(0);
        let bits = [seq::A, seq::C, seq::G, seq::T];
        let m: u64 = (0..4).filter(|&i| mask & bits[i] != 0).map(|i| c[i] as u64).sum();
        (n, m)
    }

    /// Mean pairwise identity of the column. Gaps are never identical (gap/gap
    /// pairs do not count), and a column whose comparison residue `cmp`
    /// (consensus or reference) is a gap or missing scores 0.
    /// `None` without covering sequences.
    pub fn identity(&self, col: usize, cmp: Option<u8>) -> Option<f32> {
        let c = &self.counts[col];
        let n: u64 = c[..5].iter().map(|&x| x as u64).sum();
        if n == 0 {
            return None;
        }
        if !cmp.is_some_and(|b| seq::mask(b) != 0) {
            return Some(0.0);
        }
        if n == 1 {
            return Some(if c[4] == 0 { 1.0 } else { 0.0 });
        }
        let same: u64 = c[..4].iter().map(|&x| (x as u64) * (x as u64).saturating_sub(1)).sum();
        Some(same as f32 / (n * (n - 1)) as f32)
    }

    /// Share of covering sequences identical to the comparison residue `cmp`
    /// (consensus or reference). Gaps never match; a gap or missing comparison
    /// residue scores 0. `None` without covering sequences.
    pub fn conservation(&self, col: usize, cmp: Option<u8>) -> Option<f32> {
        let (n, m) = self.matching(col, cmp);
        (n > 0).then(|| m as f32 / n as f32)
    }

    /// (identity, conservation) tracks against a per-column comparison residue.
    pub fn graph_tracks(&self, cmp: impl Fn(usize) -> Option<u8>) -> (Vec<Option<f32>>, Vec<Option<f32>>) {
        (0..self.counts.len())
            .map(|c| {
                let r = cmp(c);
                (self.identity(c, r), self.conservation(c, r))
            })
            .unzip()
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

    /// 1-based reference range covered by columns `[c0, c1)`; `None` if they
    /// contain no reference base.
    pub fn range(&self, c0: usize, c1: usize) -> Option<(u32, u32)> {
        let c1 = c1.min(self.pos_at.len());
        if c0 >= c1 {
            return None;
        }
        let first = (c0..c1).find(|&c| self.is_base[c]).map(|c| self.pos_at[c])?;
        Some((first, self.pos_at[c1 - 1]))
    }

    /// Columns `[c0, c1)` from the first to the last base of the 1-based
    /// reference positions `p0..=p1`; `None` outside the reference.
    pub fn columns(&self, p0: usize, p1: usize) -> Option<(usize, usize)> {
        let a = *self.col_of.get(p0.checked_sub(1)?)? as usize;
        let b = *self.col_of.get(p1.checked_sub(1)?)? as usize;
        (p0 <= p1).then_some((a, b + 1))
    }
}

/// Keeps the entries of `data` (columns `c0..`) whose column is flagged in
/// `only`; everything without `only`.
pub fn pick_columns(data: Vec<u8>, c0: usize, only: Option<&[bool]>) -> Vec<u8> {
    match only {
        Some(m) => data.into_iter().enumerate().filter(|(i, _)| m.get(c0 + i).copied().unwrap_or(false)).map(|(_, b)| b).collect(),
        None => data,
    }
}

/// The alignment columns shown in the viewer, numbered from 0 as display
/// columns: all of them, or a subset (gaps collapsed: only the columns where
/// the reference has a base). Positions between columns are "boundaries";
/// boundary `c` lies before column `c`.
#[derive(Debug, Clone, Default)]
pub struct ColumnMap {
    width: usize,
    subset: Option<Subset>,
}

#[derive(Debug, Clone)]
struct Subset {
    /// Alignment column of each display column.
    shown: Vec<u32>,
    /// Number of shown columns before each alignment column (len = width + 1).
    before: Vec<u32>,
}

impl ColumnMap {
    pub fn all(width: usize) -> ColumnMap {
        ColumnMap { width, subset: None }
    }

    /// Only the flagged columns; all of them if none is flagged.
    pub fn only(flags: &[bool]) -> ColumnMap {
        let shown: Vec<u32> = (0..flags.len()).filter(|&c| flags[c]).map(|c| c as u32).collect();
        if shown.is_empty() || shown.len() == flags.len() {
            return ColumnMap::all(flags.len());
        }
        let mut before = Vec::with_capacity(flags.len() + 1);
        let mut n = 0u32;
        for &f in flags {
            before.push(n);
            n += f as u32;
        }
        before.push(n);
        ColumnMap { width: flags.len(), subset: Some(Subset { shown, before }) }
    }

    /// Whether some columns are hidden.
    pub fn is_subset(&self) -> bool {
        self.subset.is_some()
    }

    /// Number of display columns.
    pub fn len(&self) -> usize {
        self.subset.as_ref().map_or(self.width, |s| s.shown.len())
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Number of hidden columns.
    pub fn hidden(&self) -> usize {
        self.width - self.len()
    }

    /// Alignment column of display column `d` (clamped to the last one).
    pub fn col(&self, d: usize) -> usize {
        match &self.subset {
            Some(s) => s.shown[d.min(s.shown.len() - 1)] as usize,
            None => d.min(self.width.saturating_sub(1)),
        }
    }

    /// Display boundary of alignment boundary `c`: the number of shown columns before it.
    pub fn disp(&self, c: usize) -> usize {
        match &self.subset {
            Some(s) => s.before[c.min(self.width)] as usize,
            None => c.min(self.width),
        }
    }

    /// Hidden columns in `[c0, c1)`.
    pub fn hidden_in(&self, c0: usize, c1: usize) -> usize {
        (c1.min(self.width).saturating_sub(c0)).saturating_sub(self.disp(c1) - self.disp(c0).min(self.disp(c1)))
    }

    /// Whether columns are hidden between display columns `d - 1` and `d`.
    pub fn hidden_before(&self, d: usize) -> bool {
        d > 0 && d < self.len() && self.col(d) > self.col(d - 1) + 1
    }

    /// Range of an annotation whose start (`left`) or end edge is dragged to
    /// display boundary `d` while the other end stays at alignment boundary
    /// `fixed`. The moved edge snaps to the shown column next to it (hidden
    /// columns there stay outside), and the range keeps at least one column.
    pub fn resize(&self, fixed: usize, left: bool, d: usize) -> (usize, usize) {
        if left {
            let d = d.min(self.disp(fixed).saturating_sub(1));
            let start = if d < self.len() { self.col(d) } else { self.width };
            (start.min(fixed.saturating_sub(1)), fixed)
        } else {
            let d = d.max(self.disp(fixed) + 1).min(self.len());
            let end = if d > 0 { self.col(d - 1) + 1 } else { 0 };
            (fixed, end.max(fixed + 1))
        }
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
        assert_eq!(s.identity(0, s.majority_state(0)), Some(1.0));
        assert!(s.identity(3, s.majority_state(3)).unwrap() < 1.0);
        assert_eq!(s.conservation(3, s.majority_state(3)), Some(0.75));
        // Against a different comparison residue only matching sequences count.
        assert_eq!(s.conservation(3, Some(b'T')), Some(0.25));
    }

    #[test]
    fn gaps_are_never_identical() {
        // Column 1 is mostly gap: consensus state is a gap -> 0, not "conserved".
        let d = doc(&["AAA", "A-A", "A-A", "A-A", "ACA"]);
        let s = ColumnStats::compute(&d, None);
        assert_eq!(s.majority_state(1), Some(seq::GAP));
        assert_eq!(s.identity(1, s.majority_state(1)), Some(0.0));
        assert_eq!(s.conservation(1, s.majority_state(1)), Some(0.0));
        // Compared to a reference base, only sequences carrying it count.
        assert_eq!(s.conservation(1, Some(b'A')), Some(0.2));
        // Gap/gap pairs no longer count as identical pairs.
        let id = s.identity(1, Some(b'A')).unwrap();
        assert!(id < 0.01, "{id}");
        // Reference without coverage at a column: nothing is identical.
        assert_eq!(s.conservation(0, None), Some(0.0));
    }

    #[test]
    fn stats_exclude_reference() {
        let d = doc(&["TTTT", "AAAA"]);
        let s = ColumnStats::compute(&d, Some(0));
        assert_eq!(s.counts[0], [1, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn slices_columns() {
        //            ref       insertion column 2 only in the reference
        let d = doc(&["ACGTAC", "AC-TAC", "--GTA-", "A-----"]);
        let s = d.slice_columns(1, 5, Some(0), None);
        // Row 3 has no base in 1..5 and is dropped; column 2 is kept (row 2 has G).
        assert_eq!(s.rows.len(), 2);
        assert_eq!(s.rows[0].name, "s1");
        assert_eq!(s.rows[0].data, b"C-TA");
        assert_eq!((s.rows[1].start, s.rows[1].data.as_slice()), (1, &b"GTA"[..]));
        // A column that is a gap in every remaining row disappears.
        let d = doc(&["ACGT", "A-GT", "A-GT"]);
        let s = d.slice_columns(0, 4, Some(0), None);
        assert_eq!(s.width, 3);
        assert_eq!(s.rows[0].data, b"AGT");
        // Only flagged columns (the reference's bases): insertions in other rows
        // go, and a row with bases only there is dropped.
        let d = doc(&["AC--GT", "ACTTGT", "--TT--", "AC--GA"]);
        let rc = RefCoords::new(&d.rows[0], d.width);
        let s = d.slice_columns(0, 6, Some(0), Some(&rc.is_base));
        assert_eq!(s.rows.len(), 2);
        assert_eq!(s.rows[0].data, b"ACGT");
        assert_eq!(s.rows[1].data, b"ACGA");
        assert_eq!(pick_columns(b"CTTG".to_vec(), 1, Some(&rc.is_base)), b"CG");
        // Chosen rows in the given order, with the alignment column of each slice column.
        let (s, cols) = d.slice_rows(&[3, 1], 0, 6, Some(&rc.is_base));
        assert_eq!(s.rows.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(), ["s3", "s1"]);
        assert_eq!(cols, [0, 1, 4, 5]);
        let d = doc(&["ACGT", "A-GT", "A-GT"]);
        assert_eq!(d.slice_rows(&[1, 2], 0, 4, None).1, [0, 2, 3]);
        assert_eq!(pick_columns(b"CTTG".to_vec(), 1, None), b"CTTG");
    }

    #[test]
    fn column_map() {
        let all = ColumnMap::all(5);
        assert!(!all.is_subset());
        assert_eq!((all.len(), all.col(3), all.disp(3), all.col(9), all.disp(9)), (5, 3, 3, 4, 5));
        assert_eq!(all.hidden_in(0, 5), 0);
        // Reference AC--GT-: columns 2, 3 and 6 are hidden.
        let r = Row::from_gapped("r".into(), String::new(), b"AC--GT".to_vec());
        let m = ColumnMap::only(&RefCoords::new(&r, 7).is_base);
        assert!(m.is_subset());
        assert_eq!((m.len(), m.hidden()), (4, 3));
        assert_eq!((0..4).map(|d| m.col(d)).collect::<Vec<_>>(), vec![0, 1, 4, 5]);
        assert_eq!((0..=7).map(|c| m.disp(c)).collect::<Vec<_>>(), vec![0, 1, 2, 2, 2, 3, 4, 4]);
        assert_eq!(m.hidden_in(1, 5), 2);
        assert_eq!(m.hidden_in(0, 2), 0);
        assert!(m.hidden_before(2) && !m.hidden_before(1) && !m.hidden_before(3));
        // Resizing snaps to shown columns and keeps hidden ones at the moved edge outside.
        assert_eq!(m.resize(6, true, 2), (4, 6));
        assert_eq!(m.resize(0, false, 2), (0, 2));
        assert_eq!(m.resize(0, false, 4), (0, 6));
        // At least one shown column.
        assert_eq!(m.resize(5, true, 4), (4, 5));
        assert_eq!(m.resize(1, false, 0), (1, 2));
        // A fixed start inside hidden columns: they stay in, the next shown column is the minimum.
        assert_eq!(m.resize(2, false, 0), (2, 5));
        // Without hidden columns this is plain column arithmetic.
        let a = ColumnMap::all(40);
        assert_eq!(a.resize(20, true, 5), (5, 20));
        assert_eq!(a.resize(20, true, 25), (19, 20));
        assert_eq!(a.resize(10, false, 30), (10, 30));
        assert_eq!(a.resize(10, false, 3), (10, 11));
        assert_eq!(a.resize(1, true, 0), (0, 1));
        // Nothing (or everything) flagged: all columns are shown.
        assert!(!ColumnMap::only(&[false, false]).is_subset());
        assert!(!ColumnMap::only(&[true, true]).is_subset());
    }

    #[test]
    fn ref_coords() {
        let r = Row::from_gapped("r".into(), String::new(), b"AC--GT".to_vec());
        let rc = RefCoords::new(&r, 6);
        assert_eq!(rc.pos_at, vec![1, 2, 2, 2, 3, 4]);
        assert_eq!(rc.column_of(3), Some(4));
        // Ranges: insertion columns map to the preceding base.
        assert_eq!(rc.range(1, 5), Some((2, 3)));
        assert_eq!(rc.range(2, 4), None);
        assert_eq!(rc.columns(2, 3), Some((1, 5)));
        assert_eq!(rc.columns(1, 4), Some((0, 6)));
        assert_eq!(rc.columns(0, 2), None);
        assert_eq!(rc.columns(3, 5), None);
        assert_eq!(rc.columns(3, 2), None);
    }
}
