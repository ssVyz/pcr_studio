//! k-mer index of the mapping target ("index word length" in Geneious terms).

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// Multiplicative hasher for 2-bit encoded k-mers (keys are already well mixed
/// after multiplication; SipHash would dominate the lookup cost).
#[derive(Default)]
pub struct KmerHasher(u64);

impl Hasher for KmerHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(8) ^ b as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        }
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = (x ^ (x >> 29)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

pub type KmerMap<V> = HashMap<u64, V, BuildHasherDefault<KmerHasher>>;

#[inline]
fn two_bit(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Calls `f(position, code)` for every k-mer made only of A/C/G/T.
pub fn for_each_kmer(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64)) {
    if k == 0 || k > 31 || seq.len() < k {
        return;
    }
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let mut code = 0u64;
    let mut valid = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        match two_bit(b) {
            Some(v) => {
                code = ((code << 2) | v) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= k {
            f(i + 1 - k, code);
        }
    }
}

pub struct KmerIndex {
    pub k: usize,
    table: KmerMap<Vec<u32>>,
}

impl KmerIndex {
    /// Indexes all k-mers of `target`; words occurring more than `max_occ`
    /// times are dropped (repeat masking).
    pub fn new(target: &[u8], k: usize, max_occ: usize) -> KmerIndex {
        let mut table: KmerMap<Vec<u32>> = KmerMap::default();
        for_each_kmer(target, k, |pos, code| table.entry(code).or_default().push(pos as u32));
        if max_occ > 0 {
            table.retain(|_, v| v.len() <= max_occ);
        }
        KmerIndex { k, table }
    }

    /// Sorted, de-duplicated (query position, target position) k-mer hits.
    pub fn anchors(&self, query: &[u8]) -> Vec<(u32, u32)> {
        let mut out = Vec::new();
        for_each_kmer(query, self.k, |qpos, code| {
            if let Some(hits) = self.table.get(&code) {
                for &t in hits {
                    out.push((qpos as u32, t));
                }
            }
        });
        out.sort_unstable();
        out.dedup();
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_anchors() {
        let idx = KmerIndex::new(b"GATTATT", 2, 0);
        let a = idx.anchors(b"TTAT");
        // TT@2,5; TA@3; AT@1,4 (0-based) -- the white-paper example.
        assert_eq!(a, vec![(0, 2), (0, 5), (1, 3), (2, 1), (2, 4)]);
    }

    #[test]
    fn skips_ambiguous() {
        let mut seen = Vec::new();
        for_each_kmer(b"ACNGTA", 2, |p, _| seen.push(p));
        assert_eq!(seen, vec![0, 3, 4]);
    }
}
