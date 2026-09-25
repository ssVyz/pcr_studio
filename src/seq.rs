//! Nucleotide helpers: normalization, IUPAC ambiguity sets, reverse complement.

/// Gap character used everywhere inside the program.
pub const GAP: u8 = b'-';

/// Bit masks for the four bases. IUPAC codes are unions of these.
pub const A: u8 = 1;
pub const C: u8 = 2;
pub const G: u8 = 4;
pub const T: u8 = 8;

/// Normalizes one input character to the internal alphabet:
/// uppercase IUPAC nucleotides plus `-` for gaps. `U` becomes `T`,
/// `.` and `~` become gaps, `?` and anything unknown become `N`.
/// Returns `None` for whitespace and digits (they are dropped).
pub fn normalize_byte(b: u8) -> Option<u8> {
    let u = b.to_ascii_uppercase();
    match u {
        b'A' | b'C' | b'G' | b'T' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D'
        | b'H' | b'V' | b'N' => Some(u),
        b'U' => Some(b'T'),
        b'-' | b'.' | b'~' => Some(GAP),
        b'?' | b'X' => Some(b'N'),
        _ if u.is_ascii_whitespace() || u.is_ascii_digit() || u == b'*' => None,
        _ => Some(b'N'),
    }
}

pub fn normalize(seq: &[u8]) -> Vec<u8> {
    seq.iter().filter_map(|&b| normalize_byte(b)).collect()
}

/// IUPAC code -> base mask. Gaps and unknown characters map to 0.
#[inline]
pub fn mask(b: u8) -> u8 {
    match b {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' | b'U' => T,
        b'R' => A | G,
        b'Y' => C | T,
        b'S' => G | C,
        b'W' => A | T,
        b'K' => G | T,
        b'M' => A | C,
        b'B' => C | G | T,
        b'D' => A | G | T,
        b'H' => A | C | T,
        b'V' => A | C | G,
        b'N' => A | C | G | T,
        _ => 0,
    }
}

/// Base mask -> IUPAC code (0 -> `-`).
pub fn code(mask: u8) -> u8 {
    const CODES: [u8; 16] = [
        b'-', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B',
        b'N',
    ];
    CODES[(mask & 15) as usize]
}

/// True for the four unambiguous bases.
#[inline]
pub fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Index used in column statistics: A,C,G,T -> 0..4, gap -> 4, anything else -> 5.
#[inline]
pub fn stat_index(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        GAP => 4,
        _ => 5,
    }
}

pub const STAT_CHARS: [u8; 6] = [b'A', b'C', b'G', b'T', GAP, b'N'];

/// Two IUPAC characters are compatible when their base sets intersect.
#[inline]
pub fn compatible(a: u8, b: u8) -> bool {
    mask(a) & mask(b) != 0
}

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other, // S, W, N, gap
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Removes gap characters.
pub fn ungap(seq: &[u8]) -> Vec<u8> {
    seq.iter().copied().filter(|&b| b != GAP).collect()
}

/// GC content in percent of the non-gap positions. S counts as GC, W as AT,
/// other ambiguity codes contribute their GC fraction.
pub fn gc_percent(seq: &[u8]) -> f64 {
    let mut gc = 0.0;
    let mut n = 0.0;
    for &b in seq {
        let m = mask(b);
        if m == 0 {
            continue;
        }
        let total = m.count_ones() as f64;
        let gcs = (m & (C | G)).count_ones() as f64;
        gc += gcs / total;
        n += 1.0;
    }
    if n == 0.0 { 0.0 } else { 100.0 * gc / n }
}

/// Expands an IUPAC oligo into all concrete sequences, up to `limit` variants.
/// Returns `None` if the expansion would exceed the limit.
pub fn expand_ambiguities(seq: &[u8], limit: usize) -> Option<Vec<Vec<u8>>> {
    let mut out: Vec<Vec<u8>> = vec![Vec::with_capacity(seq.len())];
    for &b in seq {
        let m = mask(b);
        let bases: Vec<u8> = [b'A', b'C', b'G', b'T']
            .into_iter()
            .zip([A, C, G, T])
            .filter(|(_, bit)| m & bit != 0)
            .map(|(c, _)| c)
            .collect();
        if bases.is_empty() {
            continue;
        }
        if out.len() * bases.len() > limit {
            return None;
        }
        if bases.len() == 1 {
            for s in &mut out {
                s.push(bases[0]);
            }
        } else {
            let mut next = Vec::with_capacity(out.len() * bases.len());
            for s in &out {
                for &nb in &bases {
                    let mut v = s.clone();
                    v.push(nb);
                    next.push(v);
                }
            }
            out = next;
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalization() {
        assert_eq!(normalize(b"acgu. \n~?x-12"), b"ACGT--NN-".to_vec());
    }

    #[test]
    fn iupac_roundtrip() {
        for c in b"ACGTRYSWKMBDHVN" {
            assert_eq!(code(mask(*c)), *c);
        }
        assert!(compatible(b'R', b'A'));
        assert!(!compatible(b'R', b'C'));
    }

    #[test]
    fn revcomp_works() {
        assert_eq!(revcomp(b"ACGTRN-"), b"-NYACGT".to_vec());
    }

    #[test]
    fn gc_content() {
        assert!((gc_percent(b"GGCCAATT") - 50.0).abs() < 1e-9);
        assert!((gc_percent(b"SSWW") - 50.0).abs() < 1e-9);
    }

    #[test]
    fn expansion() {
        let v = expand_ambiguities(b"ARY", 16).unwrap();
        assert_eq!(v.len(), 4);
        assert!(expand_ambiguities(b"NNNNN", 100).is_none());
    }
}
