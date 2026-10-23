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

/// Bases of a full-width row (gaps and blank columns skipped) with their
/// alignment columns, prepared for placing oligos.
#[derive(Debug, Clone, Default)]
pub struct Template {
    masks: Vec<u8>,
    cols: Vec<usize>,
}

impl Template {
    pub fn new(row: &[u8]) -> Template {
        let (masks, cols) = row.iter().enumerate().filter(|&(_, &b)| b != GAP && b != b' ').map(|(c, &b)| (seq::mask(b), c)).unzip();
        Template { masks, cols }
    }
}

/// Best position of an oligo on a template.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Placement {
    pub c0: usize,
    /// Exclusive end column.
    pub c1: usize,
    /// The reverse complement matches (minus strand).
    pub reverse: bool,
    pub mismatches: usize,
    /// Other positions (either strand) with as few mismatches.
    pub others: usize,
}

const N_MASK: u8 = seq::A | seq::C | seq::G | seq::T;

/// Position of `oligo` (5'→3') with the fewest mismatches on either strand
/// of `t`, if that is at most `max_mismatches`. IUPAC codes in the oligo
/// match the bases they include; an N in the template only matches an N, so
/// runs of N do not attract oligos. Matches may span gaps. Of equally good
/// positions the leftmost wins, plus strand first.
pub fn best_placement(t: &Template, oligo: &[u8], max_mismatches: usize) -> Option<Placement> {
    let fwd = seq::ungap(oligo);
    let m = fwd.len();
    if m == 0 || t.masks.len() < m {
        return None;
    }
    let rc = seq::revcomp(&fwd);
    let mut strands = vec![(fwd, false)];
    if rc != strands[0].0 {
        strands.push((rc, true));
    }
    let strands: Vec<(Vec<u8>, bool)> = strands.into_iter().map(|(s, r)| (s.iter().map(|&b| seq::mask(b)).collect(), r)).collect();
    let differs = |tb: u8, ob: u8| tb & ob == 0 || (tb == N_MASK && ob != N_MASK);
    let mut best: Option<Placement> = None;
    let mut bound = max_mismatches;
    for i in 0..=t.masks.len() - m {
        let window = &t.masks[i..i + m];
        for (motif, reverse) in &strands {
            let mut mm = 0usize;
            for (&tb, &ob) in window.iter().zip(motif) {
                if differs(tb, ob) {
                    mm += 1;
                    if mm > bound {
                        break;
                    }
                }
            }
            if mm > bound {
                continue;
            }
            match &mut best {
                Some(b) if b.mismatches == mm => b.others += 1,
                _ => {
                    best = Some(Placement { c0: t.cols[i], c1: t.cols[i + m - 1] + 1, reverse: *reverse, mismatches: mm, others: 0 });
                    bound = mm;
                }
            }
        }
    }
    best
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

    #[test]
    fn places_oligos_on_both_strands() {
        //                       0         1         2
        //                       01234567890123456789012345
        let t = Template::new(b"  GATTACAGG--CCTTAGACCNNNN");
        // Plus strand, across the gap columns 11-12.
        let p = best_placement(&t, b"CAGGCCTT", 0).unwrap();
        assert_eq!((p.c0, p.c1, p.reverse, p.mismatches, p.others), (7, 17, false, 0, 0));
        // Reverse complement of GATTAC: minus strand.
        let p = best_placement(&t, b"GTAATC", 0).unwrap();
        assert_eq!((p.c0, p.c1, p.reverse), (2, 8, true));
        // Degenerate oligo (R = A/G, Y = C/T).
        assert_eq!(best_placement(&t, b"RATTAY", 0).map(|p| (p.c0, p.reverse)), Some((2, false)));
        // One mismatch: only found when allowed.
        assert_eq!(best_placement(&t, b"GATCAC", 0), None);
        assert_eq!(best_placement(&t, b"GATCAC", 1).map(|p| (p.c0, p.mismatches)), Some((2, 1)));
        // Template N runs count as mismatches, N in the oligo matches anything.
        assert_eq!(best_placement(&t, b"ACCAAAA", 3), None);
        assert_eq!(best_placement(&t, b"ACCNNNN", 0).map(|p| p.c0), Some(19));
    }

    #[test]
    fn counts_equally_good_positions() {
        let t = Template::new(b"ACGTTTACGTTTTTCGTA");
        let p = best_placement(&t, b"ACGT", 0).unwrap();
        // ACGT is its own reverse complement: searched once, found twice.
        assert_eq!((p.c0, p.reverse, p.others), (0, false, 1));
        // Windows with one mismatch come first; the exact match at 13 replaces them.
        let p = best_placement(&t, b"TCGT", 1).unwrap();
        assert_eq!((p.c0, p.mismatches, p.others), (13, 0, 0));
        assert_eq!(best_placement(&t, b"", 3), None);
        assert_eq!(best_placement(&t, &[b'A'; 40], 40), None);
    }
}
