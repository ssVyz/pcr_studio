//! Motif search (IUPAC aware, optional mismatches, both strands) over the
//! reference, the consensus or all sequences. Hits are reported in alignment
//! columns; matches may span gaps.

use crate::model::{Document, Row};
use crate::seq::{self, GAP};
use rayon::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scope {
    Reference,
    Consensus,
    AllSequences,
}

impl std::fmt::Display for Scope {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Scope::Reference => "Reference",
            Scope::Consensus => "Consensus",
            Scope::AllSequences => "All sequences",
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    /// Row index, `None` for the consensus.
    pub row: Option<usize>,
    pub c0: usize,
    /// Exclusive end column.
    pub c1: usize,
    pub reverse: bool,
    pub mismatches: u8,
}

fn search_one(data: &[u8], start: usize, motifs: &[(Vec<u8>, bool)], max_mm: usize, row: Option<usize>, out: &mut Vec<Hit>, limit: usize) {
    let mut bases = Vec::with_capacity(data.len());
    let mut cols = Vec::with_capacity(data.len());
    for (i, &b) in data.iter().enumerate() {
        if b != GAP && b != b' ' {
            bases.push(seq::mask(b));
            cols.push(start + i);
        }
    }
    for (motif, reverse) in motifs {
        let m: Vec<u8> = motif.iter().map(|&b| seq::mask(b)).collect();
        if m.is_empty() || bases.len() < m.len() {
            continue;
        }
        'pos: for i in 0..=bases.len() - m.len() {
            let mut mm = 0usize;
            for j in 0..m.len() {
                if bases[i + j] & m[j] == 0 {
                    mm += 1;
                    if mm > max_mm {
                        continue 'pos;
                    }
                }
            }
            out.push(Hit { row, c0: cols[i], c1: cols[i + m.len() - 1] + 1, reverse: *reverse, mismatches: mm as u8 });
            if out.len() >= limit {
                return;
            }
        }
    }
}

/// Searches `motif`. `consensus` is the full-width consensus row.
pub fn search(
    doc: &Document,
    reference: Option<usize>,
    consensus: &[u8],
    motif: &str,
    scope: Scope,
    max_mismatches: usize,
    both_strands: bool,
    limit: usize,
) -> Vec<Hit> {
    let motif = seq::normalize(motif.as_bytes());
    let motif = seq::ungap(&motif);
    if motif.is_empty() {
        return Vec::new();
    }
    let mut motifs = vec![(motif.clone(), false)];
    let rc = seq::revcomp(&motif);
    if both_strands && rc != motif {
        motifs.push((rc, true));
    }
    let mut hits = Vec::new();
    match scope {
        Scope::Consensus => search_one(consensus, 0, &motifs, max_mismatches, None, &mut hits, limit),
        Scope::Reference => {
            if let Some(r) = reference {
                let row: &Row = &doc.rows[r];
                search_one(&row.data, row.start, &motifs, max_mismatches, Some(r), &mut hits, limit);
            }
        }
        Scope::AllSequences => {
            let per_row: Vec<Vec<Hit>> = doc
                .rows
                .par_iter()
                .enumerate()
                .map(|(i, row)| {
                    let mut v = Vec::new();
                    search_one(&row.data, row.start, &motifs, max_mismatches, Some(i), &mut v, limit);
                    v
                })
                .collect();
            for v in per_row {
                hits.extend(v);
                if hits.len() >= limit {
                    hits.truncate(limit);
                    break;
                }
            }
        }
    }
    hits.sort_by_key(|h| (h.c0, h.row, h.reverse));
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_across_gaps_and_strands() {
        let doc = Document::new(vec![
            Row::from_gapped("r".into(), String::new(), b"AACG-TTGGG".to_vec()),
            Row::from_gapped("s".into(), String::new(), b"AACCCTTGGG".to_vec()),
        ]);
        let hits = search(&doc, Some(0), b"", "CGTT", Scope::Reference, 0, false, 100);
        assert_eq!(hits, vec![Hit { row: Some(0), c0: 2, c1: 7, reverse: false, mismatches: 0 }]);
        // AACG is the reverse complement of CGTT.
        let hits = search(&doc, Some(0), b"", "CGTT", Scope::AllSequences, 0, true, 100);
        assert!(hits.iter().any(|h| h.reverse && h.c0 == 0));
        let hits = search(&doc, Some(0), b"", "CCCN", Scope::AllSequences, 0, false, 100);
        assert_eq!(hits.len(), 1);
        let hits = search(&doc, Some(0), b"", "CGTA", Scope::Reference, 1, false, 100);
        assert_eq!(hits.len(), 1);
    }
}
