//! Pairwise alignment of one query against the mapping target.
//!
//! Seed-and-extend in the spirit of the Geneious mapper: k-mer hits between
//! query and target are chained with rust-bio's sparse DP (`sdpkpp`), the
//! chain is turned into exact-match segments, and everything between and
//! around the segments is filled with rust-bio's dynamic-programming aligner.
//! The query is aligned end to end (global on the query, local on the target);
//! optionally, query parts that hang over the ends of the target are clipped.

use super::index::KmerIndex;
use crate::seq;
use bio::alignment::AlignmentOperation as Op;
use bio::alignment::pairwise::{Aligner, MIN_SCORE, Scoring};
use bio::alignment::sparse;

#[derive(Debug, Clone, Copy)]
pub struct AlignParams {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// Clip query bases that extend past the ends of the target.
    pub clip_overhangs: bool,
    /// Largest DP matrix (cells) filled by the exact aligner.
    pub max_cells: usize,
}

impl Default for AlignParams {
    fn default() -> Self {
        AlignParams {
            match_score: 2,
            mismatch: -3,
            gap_open: -6,
            gap_extend: -1,
            clip_overhangs: true,
            max_cells: 6_000_000,
        }
    }
}

const CLIP_PENALTY: i32 = -24;

/// A query aligned to target positions `t_start .. t_start + aligned.len()`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PairAln {
    pub t_start: usize,
    /// Query character (or `-`) aligned to each target position.
    pub aligned: Vec<u8>,
    /// Insertions: (target position they precede, inserted bases).
    /// Position `t_start` = before the first aligned base, `t_end()` = after the last.
    pub ins: Vec<(u32, Vec<u8>)>,
    pub reverse: bool,
    pub matches: u32,
    pub mismatches: u32,
    pub gap_cols: u32,
    pub clipped_5: u32,
    pub clipped_3: u32,
    pub score: i64,
}

impl PairAln {
    pub fn t_end(&self) -> usize {
        self.t_start + self.aligned.len()
    }

    pub fn identity(&self) -> f32 {
        let d = self.matches + self.mismatches + self.gap_cols;
        if d == 0 { 0.0 } else { self.matches as f32 / d as f32 }
    }

    /// Query bases placed in the alignment (not clipped).
    pub fn aligned_query_len(&self) -> usize {
        self.aligned.iter().filter(|&&b| b != seq::GAP).count() + self.ins.iter().map(|(_, v)| v.len()).sum::<usize>()
    }

    /// Recomputes match statistics against the target.
    pub fn compute_stats(&mut self, target: &[u8]) {
        let (mut m, mut mm, mut g) = (0u32, 0u32, 0u32);
        for (i, &q) in self.aligned.iter().enumerate() {
            let t = target[self.t_start + i];
            if q == seq::GAP {
                g += 1;
            } else if q == t && seq::is_acgt(q) {
                m += 1;
            } else if !seq::compatible(q, t) {
                mm += 1;
            }
        }
        let t_end = self.t_end();
        for (pos, bases) in &self.ins {
            let overhang = (*pos as usize == self.t_start && self.t_start == 0) || (*pos as usize == t_end && t_end == target.len());
            if !overhang {
                g += bases.len() as u32;
            }
        }
        self.matches = m;
        self.mismatches = mm;
        self.gap_cols = g;
        self.score = m as i64 * 2 - mm as i64 * 3 - g as i64 * 2;
    }
}

/// Incrementally assembles a [`PairAln`].
struct Builder {
    t_start: Option<usize>,
    t_cur: usize,
    aligned: Vec<u8>,
    ins: Vec<(u32, Vec<u8>)>,
    pending: Vec<u8>,
    clipped_5: u32,
    clipped_3: u32,
}

impl Builder {
    fn new(t_cur: usize) -> Builder {
        Builder { t_start: None, t_cur, aligned: Vec::new(), ins: Vec::new(), pending: Vec::new(), clipped_5: 0, clipped_3: 0 }
    }

    fn begin(&mut self) {
        if self.t_start.is_none() {
            self.t_start = Some(self.t_cur);
        }
        if !self.pending.is_empty() {
            self.ins.push((self.t_cur as u32, std::mem::take(&mut self.pending)));
        }
    }

    fn m(&mut self, q: u8) {
        self.begin();
        self.aligned.push(q);
        self.t_cur += 1;
    }

    fn d(&mut self) {
        self.begin();
        self.aligned.push(seq::GAP);
        self.t_cur += 1;
    }

    fn i(&mut self, q: u8) {
        self.pending.push(q);
    }

    fn skip_target(&mut self, n: usize) {
        if self.t_start.is_none() && self.pending.is_empty() {
            self.t_cur += n;
        }
    }

    fn finish(mut self, reverse: bool, target: &[u8]) -> Option<PairAln> {
        if self.t_start.is_none() {
            self.t_start = Some(self.t_cur);
        }
        if !self.pending.is_empty() {
            self.ins.push((self.t_cur as u32, std::mem::take(&mut self.pending)));
        }
        if self.aligned.is_empty() {
            return None;
        }
        let mut aln = PairAln {
            t_start: self.t_start.unwrap(),
            aligned: self.aligned,
            ins: self.ins,
            reverse,
            clipped_5: self.clipped_5,
            clipped_3: self.clipped_3,
            ..Default::default()
        };
        aln.compute_stats(target);
        Some(aln)
    }
}

#[derive(Debug, Clone, Copy)]
struct Seg {
    q: usize,
    t: usize,
    len: usize,
}

/// Turns a chain of k-mer hits into non-overlapping exact-match segments.
fn chain_segments(anchors: &[(u32, u32)], path: &[usize], k: usize) -> Vec<Seg> {
    let mut segs: Vec<Seg> = Vec::new();
    for &i in path {
        let (x, y) = (anchors[i].0 as usize, anchors[i].1 as usize);
        if let Some(last) = segs.last_mut() {
            let (qe, te) = (last.q + last.len, last.t + last.len);
            let same_diag = y as i64 - x as i64 == last.t as i64 - last.q as i64;
            if same_diag && x <= qe {
                if x + k > qe {
                    last.len = x + k - last.q;
                }
                continue;
            }
            let shift = qe.saturating_sub(x).max(te.saturating_sub(y));
            if shift >= k {
                continue;
            }
            segs.push(Seg { q: x + shift, t: y + shift, len: k - shift });
        } else {
            segs.push(Seg { q: x, t: y, len: k });
        }
    }
    segs
}

fn chain(anchors: &[(u32, u32)], k: usize) -> (u32, Vec<usize>) {
    if anchors.is_empty() {
        return (0, Vec::new());
    }
    let r = sparse::sdpkpp(anchors, k, 1, -2, -1);
    (r.score, r.path)
}

struct Ctx<'a> {
    p: &'a AlignParams,
}

impl Ctx<'_> {
    fn scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32 + use<>> {
        let (ms, mm) = (self.p.match_score, self.p.mismatch);
        let f = move |a: u8, b: u8| {
            if a == b && seq::is_acgt(a) {
                ms
            } else if seq::compatible(a, b) {
                0
            } else {
                mm
            }
        };
        Scoring::new(self.p.gap_open, self.p.gap_extend, f)
    }

    fn apply(&self, b: &mut Builder, q: &[u8], ops: &[Op], prefix_side: bool) {
        let mut qi = 0usize;
        for op in ops {
            match *op {
                Op::Match | Op::Subst => {
                    b.m(q[qi]);
                    qi += 1;
                }
                Op::Del => b.d(),
                Op::Ins => {
                    b.i(q[qi]);
                    qi += 1;
                }
                Op::Xclip(n) => {
                    if n == 0 {
                        continue;
                    }
                    if prefix_side && qi == 0 {
                        b.clipped_5 += n as u32;
                    } else {
                        b.clipped_3 += n as u32;
                    }
                    qi += n;
                }
                Op::Yclip(n) => b.skip_target(n),
            }
        }
    }

    /// Overhanging query bases: clipped or inserted.
    fn overhang(&self, b: &mut Builder, q: &[u8], five_prime: bool) {
        if self.p.clip_overhangs {
            if five_prime {
                b.clipped_5 += q.len() as u32;
            } else {
                b.clipped_3 += q.len() as u32;
            }
        } else {
            for &x in q {
                b.i(x);
            }
        }
    }

    /// Aligns the query part before the first anchor; the alignment must end
    /// exactly at target position `ts`.
    fn prefix(&self, b: &mut Builder, q: &[u8], t: &[u8], ts: usize) {
        b.t_cur = ts;
        if q.is_empty() {
            return;
        }
        let slack = q.len() / 4 + 32;
        let w0 = ts.saturating_sub(q.len() + slack);
        let window = &t[w0..ts];
        if window.is_empty() {
            self.overhang(b, q, true);
            b.t_cur = ts;
            return;
        }
        if q.len() * window.len() > self.p.max_cells {
            let n = q.len().min(window.len());
            let extra = q.len() - n;
            self.overhang(b, &q[..extra], true);
            b.t_cur = ts - n;
            for &x in &q[extra..] {
                b.m(x);
            }
            return;
        }
        let mut scoring = self.scoring().yclip_prefix(0);
        if self.p.clip_overhangs && w0 == 0 {
            scoring = scoring.xclip_prefix(CLIP_PENALTY);
        }
        let mut aligner = Aligner::with_scoring(scoring);
        let aln = aligner.custom(q, window);
        b.t_cur = w0;
        self.apply(b, q, &aln.operations, true);
        debug_assert_eq!(b.t_cur, ts);
        b.t_cur = ts;
    }

    /// Aligns the query part after the last anchor, starting at target `te`.
    fn suffix(&self, b: &mut Builder, q: &[u8], t: &[u8], te: usize) {
        if q.is_empty() {
            return;
        }
        let slack = q.len() / 4 + 32;
        let w1 = (te + q.len() + slack).min(t.len());
        let window = &t[te..w1];
        if window.is_empty() {
            self.overhang(b, q, false);
            return;
        }
        if q.len() * window.len() > self.p.max_cells {
            let n = q.len().min(window.len());
            for &x in &q[..n] {
                b.m(x);
            }
            self.overhang(b, &q[n..], false);
            return;
        }
        let mut scoring = self.scoring().yclip_suffix(0);
        if self.p.clip_overhangs && w1 == t.len() {
            scoring = scoring.xclip_suffix(CLIP_PENALTY);
        }
        let mut aligner = Aligner::with_scoring(scoring);
        let aln = aligner.custom(q, window);
        self.apply(b, q, &aln.operations, false);
    }

    /// Globally aligns a gap between two anchored segments.
    fn fill(&self, b: &mut Builder, q: &[u8], t: &[u8], k: usize, depth: usize) {
        if q.is_empty() {
            for _ in t {
                b.d();
            }
            return;
        }
        if t.is_empty() {
            for &x in q {
                b.i(x);
            }
            return;
        }
        if q.len() * t.len() <= self.p.max_cells {
            let mut aligner = Aligner::with_scoring(self.scoring());
            let aln = aligner.global(q, t);
            self.apply(b, q, &aln.operations, false);
            return;
        }
        // Too large for exact DP: re-seed inside the gap with shorter words.
        let sub_k = k.saturating_sub(3).max(8);
        if depth < 3 {
            let idx = KmerIndex::new(t, sub_k, 8);
            let anchors = idx.anchors(q);
            let (score, path) = chain(&anchors, sub_k);
            if score > 0 {
                let segs = chain_segments(&anchors, &path, sub_k);
                let (mut qp, mut tp) = (0usize, 0usize);
                for s in &segs {
                    self.fill(b, &q[qp..s.q], &t[tp..s.t], sub_k, depth + 1);
                    for &x in &q[s.q..s.q + s.len] {
                        b.m(x);
                    }
                    qp = s.q + s.len;
                    tp = s.t + s.len;
                }
                self.fill(b, &q[qp..], &t[tp..], sub_k, depth + 1);
                return;
            }
        }
        // Last resort: ungapped diagonal, remaining length as one indel.
        let n = q.len().min(t.len());
        for &x in &q[..n] {
            b.m(x);
        }
        for &x in &q[n..] {
            b.i(x);
        }
        for _ in n..t.len() {
            b.d();
        }
    }
}

/// Aligns `q` (one orientation) along an anchor chain.
fn align_chain(q: &[u8], t: &[u8], anchors: &[(u32, u32)], path: &[usize], k: usize, reverse: bool, p: &AlignParams) -> Option<PairAln> {
    let segs = chain_segments(anchors, path, k);
    let first = *segs.first()?;
    let ctx = Ctx { p };
    let mut b = Builder::new(first.t);
    ctx.prefix(&mut b, &q[..first.q], t, first.t);
    let (mut qp, mut tp) = (first.q, first.t);
    for s in &segs {
        ctx.fill(&mut b, &q[qp..s.q], &t[tp..s.t], k, 0);
        for &x in &q[s.q..s.q + s.len] {
            b.m(x);
        }
        qp = s.q + s.len;
        tp = s.t + s.len;
    }
    ctx.suffix(&mut b, &q[qp..], t, tp);
    b.finish(reverse, t)
}

/// Unseeded fallback for short queries: exact semi-global DP over the whole target.
fn align_full(q: &[u8], t: &[u8], reverse: bool, p: &AlignParams) -> Option<PairAln> {
    let ctx = Ctx { p };
    let mut scoring = ctx.scoring().yclip(0);
    if p.clip_overhangs {
        scoring = scoring.xclip(CLIP_PENALTY);
    }
    let mut aligner = Aligner::with_scoring(scoring);
    let aln = aligner.custom(q, t);
    if aln.score == MIN_SCORE {
        return None;
    }
    let mut b = Builder::new(0);
    ctx.apply(&mut b, q, &aln.operations, true);
    b.finish(reverse, t)
}

/// Maps a query to the target in the better of both orientations.
pub fn map_query(q: &[u8], t: &[u8], idx: &KmerIndex, both_strands: bool, p: &AlignParams) -> Option<PairAln> {
    let k = idx.k;
    let fwd_anchors = idx.anchors(q);
    let rc = seq::revcomp(q);
    let rev_anchors = if both_strands { idx.anchors(&rc) } else { Vec::new() };
    let (fs, fpath) = chain(&fwd_anchors, k);
    let (rs, rpath) = chain(&rev_anchors, k);

    let mut candidates: Vec<PairAln> = Vec::new();
    let best = fs.max(rs);
    if best > 0 {
        if fs > 0 && fs * 2 >= best {
            candidates.extend(align_chain(q, t, &fwd_anchors, &fpath, k, false, p));
        }
        if rs > 0 && rs * 2 >= best {
            candidates.extend(align_chain(&rc, t, &rev_anchors, &rpath, k, true, p));
        }
    } else if q.len() * t.len() <= p.max_cells * 4 {
        candidates.extend(align_full(q, t, false, p));
        if both_strands {
            candidates.extend(align_full(&rc, t, true, p));
        }
    }
    candidates.into_iter().max_by_key(|a| (a.score, !a.reverse))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    fn render(a: &PairAln) -> String {
        String::from_utf8(a.aligned.clone()).unwrap()
    }

    #[test]
    fn exact_substring() {
        let t = pseudo_random(2000, 1);
        let q = t[500..800].to_vec();
        let idx = KmerIndex::new(&t, 13, 20);
        let a = map_query(&q, &t, &idx, true, &AlignParams::default()).unwrap();
        assert_eq!(a.t_start, 500);
        assert_eq!(a.aligned, q);
        assert!(!a.reverse);
        assert_eq!(a.mismatches, 0);
        assert!(a.ins.is_empty());
    }

    #[test]
    fn reverse_complement() {
        let t = pseudo_random(2000, 2);
        let q = seq::revcomp(&t[100..400]);
        let idx = KmerIndex::new(&t, 13, 20);
        let a = map_query(&q, &t, &idx, true, &AlignParams::default()).unwrap();
        assert!(a.reverse);
        assert_eq!(a.t_start, 100);
        assert_eq!(a.aligned, t[100..400].to_vec());
    }

    #[test]
    fn snp_insertion_deletion() {
        let t = pseudo_random(3000, 3);
        let mut q = t[1000..1600].to_vec();
        q[100] = if q[100] == b'A' { b'C' } else { b'A' }; // SNP
        q.drain(300..303); // 3 bp deletion in query
        q.splice(450..450, b"GGGTT".iter().copied()); // 5 bp insertion
        let idx = KmerIndex::new(&t, 13, 20);
        let a = map_query(&q, &t, &idx, true, &AlignParams::default()).unwrap();
        assert_eq!(a.t_start, 1000);
        assert_eq!(a.t_end(), 1600);
        assert_eq!(a.mismatches, 1);
        let dels = a.aligned.iter().filter(|&&b| b == seq::GAP).count();
        assert_eq!(dels, 3);
        let ins: usize = a.ins.iter().map(|(_, v)| v.len()).sum();
        assert_eq!(ins, 5);
        assert!(render(&a).len() == 600);
    }

    #[test]
    fn overhang_is_clipped_or_inserted() {
        let t = pseudo_random(1000, 4);
        let mut q = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec();
        q.extend_from_slice(&t[..300]);
        let idx = KmerIndex::new(&t, 13, 20);
        let clipped = map_query(&q, &t, &idx, false, &AlignParams::default()).unwrap();
        assert_eq!(clipped.t_start, 0);
        assert_eq!(clipped.clipped_5, 50);
        assert_eq!(clipped.mismatches, 0);
        let p = AlignParams { clip_overhangs: false, ..Default::default() };
        let kept = map_query(&q, &t, &idx, false, &p).unwrap();
        assert_eq!(kept.ins.len(), 1);
        assert_eq!(kept.ins[0].0, 0);
        assert_eq!(kept.ins[0].1.len(), 50);
    }

    #[test]
    fn divergent_long_sequence() {
        // ~10% divergence plus several indels over a genome-sized sequence.
        let t = pseudo_random(20_000, 5);
        let mut q = Vec::new();
        let noise = pseudo_random(20_000, 99);
        for (i, &b) in t.iter().enumerate() {
            if i % 997 == 0 {
                continue; // deletion
            }
            if noise[i] == b'A' && i % 3 == 0 {
                q.push(if b == b'G' { b'T' } else { b'G' });
            } else {
                q.push(b);
            }
            if i % 1511 == 0 {
                q.extend_from_slice(b"AC");
            }
        }
        let idx = KmerIndex::new(&t, 11, 20);
        let a = map_query(&q, &t, &idx, true, &AlignParams::default()).unwrap();
        assert!(a.identity() > 0.85, "identity {}", a.identity());
        // Position 0 is deleted in the query, so it starts at target position 1.
        assert_eq!(a.t_start, 1);
        assert_eq!(a.t_end(), t.len());
        assert_eq!(a.aligned_query_len(), q.len());
    }

    #[test]
    fn unrelated_short_query_fallback() {
        let t = pseudo_random(500, 6);
        let q = b"ACGTAC".to_vec();
        let idx = KmerIndex::new(&t, 13, 20);
        // Too short for any seed: falls back to full DP and still returns something.
        let a = map_query(&q, &t, &idx, false, &AlignParams::default());
        assert!(a.is_some());
    }
}
