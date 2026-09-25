//! Merges independent query-to-target alignments into one contig
//! (a target-anchored multiple alignment).
//!
//! Every target base gets a column. Insertions relative to the target are
//! pooled per target position into a block as wide as the longest insertion
//! there; inside a block insertions are left-justified, except insertions
//! before a query's first aligned base (overhangs), which are right-justified
//! so they stay attached to the query. Queries spanning a block without an
//! insertion of their own get gaps there.

use super::align::PairAln;
use crate::seq::{self, GAP};

pub struct Layout {
    /// Width of the insertion block before target position `p` (`p == len`: after the end).
    pub ins_width: Vec<u32>,
    /// Column of target base `p`; `col_of[len]` is the total width.
    pub col_of: Vec<usize>,
}

impl Layout {
    pub fn new<'a>(target_len: usize, alns: impl Iterator<Item = &'a PairAln>) -> Layout {
        let mut ins_width = vec![0u32; target_len + 1];
        for a in alns {
            for (pos, bases) in &a.ins {
                let w = &mut ins_width[*pos as usize];
                *w = (*w).max(bases.len() as u32);
            }
        }
        let mut col_of = Vec::with_capacity(target_len + 1);
        let mut acc = 0usize;
        for (p, &w) in ins_width.iter().enumerate() {
            acc += w as usize;
            col_of.push(p + acc);
        }
        Layout { ins_width, col_of }
    }

    pub fn width(&self) -> usize {
        *self.col_of.last().unwrap()
    }

    /// Renders an alignment into (start column, gapped data).
    pub fn render(&self, a: &PairAln) -> (usize, Vec<u8>) {
        let ts = a.t_start;
        let te = a.t_end();
        let mut ins = a.ins.iter().peekable();
        let mut data = Vec::with_capacity(a.aligned.len() + 16);
        let mut start = self.col_of[ts];
        if let Some((pos, bases)) = ins.peek()
            && *pos as usize == ts {
                start -= bases.len();
                data.extend_from_slice(bases);
                ins.next();
            }
        for t in ts..te {
            if t > ts {
                let w = self.ins_width[t] as usize;
                if w > 0 {
                    let mut own = 0;
                    if let Some((pos, bases)) = ins.peek()
                        && *pos as usize == t {
                            data.extend_from_slice(bases);
                            own = bases.len();
                            ins.next();
                        }
                    data.extend(std::iter::repeat_n(GAP, w - own));
                }
            }
            data.push(a.aligned[t - ts]);
        }
        if let Some((pos, bases)) = ins.peek()
            && *pos as usize == te {
                data.extend_from_slice(bases);
            }
        (start, data)
    }
}

/// Derives the next mapping target from a rendered contig: the per-column
/// majority of the queries (gap-majority columns are dropped; columns without
/// query coverage keep the reference base). Also returns the reference
/// aligned to the new target.
pub fn consensus_target(width: usize, rows: &[(usize, Vec<u8>)], reference: &(usize, Vec<u8>)) -> (Vec<u8>, PairAln) {
    let mut counts = vec![[0u32; 6]; width];
    for (start, data) in rows {
        for (i, &b) in data.iter().enumerate() {
            counts[start + i][seq::stat_index(b)] += 1;
        }
    }
    let ref_at = |c: usize| -> Option<u8> {
        if c >= reference.0 && c < reference.0 + reference.1.len() { Some(reference.1[c - reference.0]) } else { None }
    };
    let mut target = Vec::with_capacity(width);
    let mut kept: Vec<Option<usize>> = Vec::with_capacity(width); // column -> target index
    for (c, cnt) in counts.iter().enumerate() {
        let r = ref_at(c);
        let cov: u32 = cnt[..5].iter().sum();
        let call = if cov == 0 {
            match r {
                Some(b) if b != GAP => Some(b),
                _ if cnt[5] > 0 => Some(b'N'),
                _ => None,
            }
        } else {
            let best = (0..5).max_by_key(|&i| (cnt[i], r.map(|b| seq::stat_index(b) == i).unwrap_or(false), 5 - i)).unwrap();
            if best == 4 { None } else { Some(seq::STAT_CHARS[best]) }
        };
        match call {
            Some(b) => {
                kept.push(Some(target.len()));
                target.push(b);
            }
            None => kept.push(None),
        }
    }
    // Reference aligned to the new target.
    let mut ref_aln = PairAln::default();
    let mut started = false;
    let mut pending: Vec<u8> = Vec::new();
    let mut next_t = 0usize;
    let (rs, re) = (reference.0, reference.0 + reference.1.len());
    for c in rs..re {
        let r = reference.1[c - rs];
        match kept[c] {
            Some(j) => {
                if !started {
                    ref_aln.t_start = j;
                    started = true;
                }
                if !pending.is_empty() {
                    ref_aln.ins.push((j as u32, std::mem::take(&mut pending)));
                }
                ref_aln.aligned.push(r);
                next_t = j + 1;
            }
            None => {
                if r != GAP {
                    pending.push(r);
                }
            }
        }
    }
    if !started {
        ref_aln.t_start = 0;
    }
    if !pending.is_empty() {
        ref_aln.ins.push((next_t as u32, pending));
    }
    (target, ref_aln)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(t_start: usize, aligned: &str, ins: &[(u32, &str)]) -> PairAln {
        PairAln {
            t_start,
            aligned: aligned.as_bytes().to_vec(),
            ins: ins.iter().map(|(p, s)| (*p, s.as_bytes().to_vec())).collect(),
            ..Default::default()
        }
    }

    #[test]
    fn merges_insertions() {
        // target ACGTACGT
        let r = aln(0, "ACGTACGT", &[]);
        let a = aln(0, "ACGTACGT", &[(4, "TT")]);
        let b = aln(2, "GTA-GT", &[(4, "C")]);
        let c = aln(5, "CGT", &[(5, "GG")]); // leading overhang-style insertion
        let lay = Layout::new(8, [&r, &a, &b, &c].into_iter());
        // Blocks: 2 columns before target 4, 2 columns before target 5.
        assert_eq!(lay.width(), 12);
        let s = |x: (usize, Vec<u8>)| (x.0, String::from_utf8(x.1).unwrap());
        assert_eq!(s(lay.render(&r)), (0, "ACGT--A--CGT".into()));
        assert_eq!(s(lay.render(&a)), (0, "ACGTTTA--CGT".into()));
        assert_eq!(s(lay.render(&b)), (2, "GTC-A---GT".into()));
        // Right-justified leading insertion directly before its first base.
        assert_eq!(s(lay.render(&c)), (7, "GGCGT".into()));
    }

    #[test]
    fn consensus_target_drops_gap_columns() {
        let reference = (0usize, b"ACGTA".to_vec());
        let rows = vec![(0usize, b"ACCTA".to_vec()), (0, b"AC-TA".to_vec()), (0, b"A--TA".to_vec())];
        let (t, ref_aln) = consensus_target(5, &rows, &reference);
        // column 2: C, -, - -> gap majority -> dropped; column 1: C,C,- -> C
        assert_eq!(t, b"ACTA".to_vec());
        assert_eq!(ref_aln.t_start, 0);
        assert_eq!(ref_aln.aligned, b"ACTA".to_vec());
        assert_eq!(ref_aln.ins, vec![(2, b"G".to_vec())]);
    }
}
