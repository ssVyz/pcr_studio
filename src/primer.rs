//! Evaluates an oligo (primer/probe) taken from a column selection against
//! every sequence of the alignment: mismatch distribution, 3' end mismatches,
//! the most common target variants and a per-group inclusivity breakdown.

use crate::model::Document;
use crate::seq::{self, GAP};
use crate::thermo::{self, TmConditions};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OligoSource {
    Consensus,
    Reference,
}

impl std::fmt::Display for OligoSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            OligoSource::Consensus => "Consensus",
            OligoSource::Reference => "Reference",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    Reverse,
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Orientation::Forward => "Forward (5'→3' left to right)",
            Orientation::Reverse => "Reverse (reverse complement)",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PrimerOptions {
    pub source: OligoSource,
    pub orientation: Orientation,
    /// Number of 3'-terminal bases where a mismatch counts as "3' mismatch".
    pub three_prime_window: usize,
    pub tm: TmConditions,
}

impl Default for PrimerOptions {
    fn default() -> Self {
        PrimerOptions {
            source: OligoSource::Consensus,
            orientation: Orientation::Forward,
            three_prime_window: 5,
            tm: TmConditions::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Variant {
    /// Target sequence over the selection, oriented like the oligo, `.` = match.
    pub display: String,
    pub count: usize,
    pub mismatches: usize,
}

#[derive(Debug, Clone)]
pub struct GroupStat {
    pub label: String,
    pub sequences: usize,
    pub perfect: usize,
    pub one_mismatch: usize,
    pub three_prime: usize,
    pub no_coverage: usize,
}

#[derive(Debug, Clone)]
pub struct PrimerReport {
    pub col_start: usize,
    pub col_end: usize,
    /// Oligo 5'->3'.
    pub oligo: Vec<u8>,
    pub tm: Option<(f64, f64)>,
    pub gc: f64,
    /// Sequences covering the whole selection.
    pub sequences: usize,
    pub no_coverage: usize,
    /// Sequences with 0, 1, 2 and >= 3 mismatches.
    pub mismatch_hist: [usize; 4],
    /// Sequences with at least one mismatch in the 3' window.
    pub three_prime: usize,
    /// Covering sequences with N/ambiguity codes inside the selection
    /// (ambiguous positions are counted as compatible, not as mismatches).
    pub with_ambiguity: usize,
    /// Mismatch count per oligo position (5'->3').
    pub per_position: Vec<usize>,
    pub variants: Vec<Variant>,
    pub groups: Vec<GroupStat>,
    pub group_key: Option<String>,
}

impl PrimerReport {
    pub fn percent(&self, n: usize) -> f64 {
        if self.sequences == 0 { 0.0 } else { 100.0 * n as f64 / self.sequences as f64 }
    }
}

struct RowEval {
    covered: bool,
    mismatches: usize,
    three_prime: bool,
    ambiguous: bool,
}

/// Evaluates columns `[c0, c1)`. `template` is the full-width consensus or
/// reference row used as the oligo source.
pub fn evaluate(
    doc: &Document,
    reference: Option<usize>,
    template: &[u8],
    c0: usize,
    c1: usize,
    opts: &PrimerOptions,
    group_key: Option<&str>,
) -> PrimerReport {
    let c1 = c1.min(doc.width).max(c0);
    let tmpl: Vec<u8> = (c0..c1).map(|c| template.get(c).copied().unwrap_or(GAP)).collect();
    let reverse = opts.orientation == Orientation::Reverse;
    let forward_oligo = seq::ungap(&tmpl);
    let oligo = if reverse { seq::revcomp(&forward_oligo) } else { forward_oligo.clone() };
    let olen = forward_oligo.len();

    // Oligo position (0-based, 5'->3') for each selected column.
    let mut col_pos: Vec<usize> = Vec::with_capacity(tmpl.len());
    let mut idx = 0usize;
    for &b in &tmpl {
        let fwd = if b == GAP { idx.saturating_sub(1) } else { idx };
        let fwd = fwd.min(olen.saturating_sub(1));
        col_pos.push(if reverse { olen.saturating_sub(1) - fwd } else { fwd });
        if b != GAP {
            idx += 1;
        }
    }
    let window = opts.three_prime_window;
    let mut per_position = vec![0usize; olen];
    let mut variants: HashMap<Vec<u8>, (usize, usize)> = HashMap::new();

    let mut evals: Vec<Option<RowEval>> = Vec::with_capacity(doc.rows.len());
    for (i, row) in doc.rows.iter().enumerate() {
        if Some(i) == reference {
            evals.push(None);
            continue;
        }
        if tmpl.is_empty() || row.start > c0 || row.end() < c1 {
            evals.push(Some(RowEval { covered: false, mismatches: 0, three_prime: false, ambiguous: false }));
            continue;
        }
        let slice = &row.data[c0 - row.start..c1 - row.start];
        let mut mm = 0usize;
        let mut three = false;
        let mut amb = false;
        for (j, (&o, &r)) in tmpl.iter().zip(slice.iter()).enumerate() {
            let mismatch = if o == GAP && r == GAP {
                false
            } else if o == GAP || r == GAP {
                true
            } else if seq::is_acgt(r) {
                !seq::compatible(o, r)
            } else {
                amb = true;
                !seq::compatible(o, r)
            };
            if mismatch {
                mm += 1;
                if olen > 0 {
                    let p = col_pos[j];
                    per_position[p] += 1;
                    if olen - 1 - p < window {
                        three = true;
                    }
                }
            }
        }
        let e = variants.entry(slice.to_vec()).or_insert((0, mm));
        e.0 += 1;
        evals.push(Some(RowEval { covered: true, mismatches: mm, three_prime: three, ambiguous: amb }));
    }

    let mut report = PrimerReport {
        col_start: c0,
        col_end: c1,
        tm: thermo::tm_range(&oligo, &opts.tm),
        gc: seq::gc_percent(&oligo),
        oligo,
        sequences: 0,
        no_coverage: 0,
        mismatch_hist: [0; 4],
        three_prime: 0,
        with_ambiguity: 0,
        per_position,
        variants: Vec::new(),
        groups: Vec::new(),
        group_key: group_key.map(|s| s.to_string()),
    };
    for e in evals.iter().flatten() {
        if !e.covered {
            report.no_coverage += 1;
            continue;
        }
        report.sequences += 1;
        report.mismatch_hist[e.mismatches.min(3)] += 1;
        report.three_prime += e.three_prime as usize;
        report.with_ambiguity += e.ambiguous as usize;
    }

    // Most common target variants, oriented like the oligo.
    let mut vs: Vec<(Vec<u8>, usize, usize)> = variants.into_iter().map(|(k, (n, mm))| (k, n, mm)).collect();
    vs.sort_by(|a, b| b.1.cmp(&a.1).then(a.2.cmp(&b.2)).then(a.0.cmp(&b.0)));
    for (slice, count, mismatches) in vs.into_iter().take(12) {
        let (t, s) = if reverse { (seq::revcomp(&tmpl), seq::revcomp(&slice)) } else { (tmpl.clone(), slice) };
        let display: String = t
            .iter()
            .zip(s.iter())
            .filter(|(o, r)| !(**o == GAP && **r == GAP))
            .map(|(&o, &r)| if o == r && o != GAP { '.' } else { r as char })
            .collect();
        report.variants.push(Variant { display, count, mismatches });
    }

    if let Some(key) = group_key {
        let mut groups: Vec<GroupStat> = Vec::new();
        let mut by_label: HashMap<String, usize> = HashMap::new();
        for (i, e) in evals.iter().enumerate() {
            let Some(e) = e else { continue };
            let label = doc.rows[i].meta_value(key).unwrap_or("(none)").to_string();
            let gi = *by_label.entry(label.clone()).or_insert_with(|| {
                groups.push(GroupStat { label, sequences: 0, perfect: 0, one_mismatch: 0, three_prime: 0, no_coverage: 0 });
                groups.len() - 1
            });
            let g = &mut groups[gi];
            if !e.covered {
                g.no_coverage += 1;
                continue;
            }
            g.sequences += 1;
            g.perfect += (e.mismatches == 0) as usize;
            g.one_mismatch += (e.mismatches == 1) as usize;
            g.three_prime += e.three_prime as usize;
        }
        groups.sort_by(|a, b| b.sequences.cmp(&a.sequences).then(a.label.cmp(&b.label)));
        report.groups = groups;
    }
    report
}

/// Plain-text version of the report (for the clipboard).
pub fn report_text(r: &PrimerReport, ref_range: Option<(u32, u32)>) -> String {
    let mut s = String::new();
    s.push_str(&format!("5'-{}-3'\n", String::from_utf8_lossy(&r.oligo)));
    s.push_str(&format!("Length\t{}\n", r.oligo.len()));
    if let Some((a, b)) = ref_range {
        s.push_str(&format!("Reference position\t{a}-{b}\n"));
    }
    s.push_str(&format!("Tm\t{}\n", format_tm(r.tm)));
    s.push_str(&format!("GC\t{:.0}%\n", r.gc));
    s.push_str(&format!("Sequences\t{}\n", r.sequences));
    let labels = ["Perfect match", "1 mismatch", "2 mismatches", "3+ mismatches"];
    for (l, &n) in labels.iter().zip(r.mismatch_hist.iter()) {
        s.push_str(&format!("{l}\t{n}\t{:.2}%\n", r.percent(n)));
    }
    s.push_str(&format!("3' mismatch\t{}\t{:.2}%\n", r.three_prime, r.percent(r.three_prime)));
    s.push_str(&format!("No coverage\t{}\n", r.no_coverage));
    s
}

pub fn format_tm(tm: Option<(f64, f64)>) -> String {
    match tm {
        None => "–".into(),
        Some((a, b)) if (a - b).abs() < 0.05 => format!("{a:.1} °C"),
        Some((a, b)) => format!("{a:.1}–{b:.1} °C"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Row;

    fn doc(rows: &[&str]) -> Document {
        Document::new(
            rows.iter()
                .enumerate()
                .map(|(i, s)| {
                    let mut r = Row::from_gapped(format!("s{i}"), String::new(), s.as_bytes().to_vec());
                    r.set_meta("clade", if i % 2 == 0 { "A".into() } else { "B".into() });
                    r
                })
                .collect(),
        )
    }

    #[test]
    fn counts_mismatches() {
        //            reference  perfect   1mm(3')   2mm       no cov   gap
        let d = doc(&["ACGTACGT", "ACGTACGT", "ACGTACGA", "TCGAACGT", "----ACGT", "ACG-ACGT"]);
        let template = d.rows[0].data.clone();
        let opts = PrimerOptions::default();
        let r = evaluate(&d, Some(0), &template, 0, 8, &opts, Some("clade"));
        assert_eq!(r.oligo, b"ACGTACGT");
        assert_eq!(r.sequences, 4);
        assert_eq!(r.no_coverage, 1);
        assert_eq!(r.mismatch_hist, [1, 2, 1, 0]);
        // Default window = last 5 bases: the mismatches at position 4 of rows 3 and 5 count too.
        assert_eq!(r.three_prime, 3);
        let narrow = evaluate(&d, Some(0), &template, 0, 8, &PrimerOptions { three_prime_window: 2, ..opts }, None);
        assert_eq!(narrow.three_prime, 1);
        assert_eq!(r.per_position[7], 1);
        assert_eq!(r.variants[0].count, 1);
        assert!(r.groups.len() == 2);
    }

    #[test]
    fn reverse_orientation() {
        let d = doc(&["ACGTACGT", "TCGTACGT"]);
        let template = d.rows[0].data.clone();
        let opts = PrimerOptions { orientation: Orientation::Reverse, ..Default::default() };
        let r = evaluate(&d, Some(0), &template, 0, 8, &opts, None);
        assert_eq!(r.oligo, seq::revcomp(b"ACGTACGT"));
        // Mismatch at the left end of the alignment = 3' end of a reverse primer.
        assert_eq!(r.three_prime, 1);
        assert_eq!(*r.per_position.last().unwrap(), 1);
    }
}
