//! Oligonucleotide melting temperature (nearest-neighbor model).
//!
//! SantaLucia (1998) unified nearest-neighbor parameters with the SantaLucia
//! salt correction; divalent cations are converted to a monovalent equivalent
//! with the von Ahsen et al. (2001) formula `120 * sqrt([Mg] - [dNTP])`.
//! These are the same model choices as primer3's defaults
//! (`santalucia_auto` + `santalucia` salt correction).

use crate::seq;
use serde::{Deserialize, Serialize};

/// Reaction conditions used for Tm calculations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TmConditions {
    /// Oligo concentration in nM.
    pub oligo_nm: f64,
    /// Monovalent cation concentration (Na+/K+) in mM.
    pub monovalent_mm: f64,
    /// Divalent cation concentration (Mg2+) in mM.
    pub divalent_mm: f64,
    /// dNTP concentration in mM.
    pub dntp_mm: f64,
}

impl Default for TmConditions {
    fn default() -> Self {
        TmConditions { oligo_nm: 50.0, monovalent_mm: 50.0, divalent_mm: 1.5, dntp_mm: 0.6 }
    }
}

const R: f64 = 1.987;
const KELVIN: f64 = 273.15;

/// Nearest-neighbor stack (dH kcal/mol, dS cal/K/mol) for the 5'->3'
/// dinucleotide `a b` on the top strand. Bases indexed A=0 C=1 G=2 T=3.
const NN: [[(f64, f64); 4]; 4] = [
    // A?        AA             AC             AG             AT
    [(-7.9, -22.2), (-8.4, -22.4), (-7.8, -21.0), (-7.2, -20.4)],
    // C?        CA             CC             CG             CT
    [(-8.5, -22.7), (-8.0, -19.9), (-10.6, -27.2), (-7.8, -21.0)],
    // G?        GA             GC             GG             GT
    [(-8.2, -22.2), (-9.8, -24.4), (-8.0, -19.9), (-8.4, -22.4)],
    // T?        TA             TC             TG             TT
    [(-7.2, -21.3), (-8.2, -22.2), (-8.5, -22.7), (-7.9, -22.2)],
];

fn idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Monovalent-equivalent concentration (mM) of the divalent cations.
pub fn divalent_to_monovalent(divalent: f64, dntp: f64) -> f64 {
    if divalent <= 0.0 {
        return 0.0;
    }
    let free = (divalent - dntp.max(0.0)).max(0.0);
    120.0 * free.sqrt()
}

fn is_self_complementary(s: &[u8]) -> bool {
    s == seq::revcomp(s).as_slice()
}

/// Melting temperature (°C) of an unambiguous oligo (A/C/G/T only).
/// Returns `None` for sequences shorter than 2 bases or with other characters.
pub fn tm_exact(oligo: &[u8], cond: &TmConditions) -> Option<f64> {
    if oligo.len() < 2 {
        return None;
    }
    let ids: Vec<usize> = oligo.iter().map(|&b| idx(b)).collect::<Option<_>>()?;
    let mut dh = 0.0; // kcal/mol
    let mut ds = 0.0; // cal/K/mol
    for w in ids.windows(2) {
        let (h, s) = NN[w[0]][w[1]];
        dh += h;
        ds += s;
    }
    // Initiation, applied per terminal base pair.
    for &end in [ids[0], ids[ids.len() - 1]].iter() {
        if end == 0 || end == 3 {
            dh += 2.3;
            ds += 4.1;
        } else {
            dh += 0.1;
            ds += -2.8;
        }
    }
    let symmetric = is_self_complementary(oligo);
    if symmetric {
        ds += -1.4;
    }
    let na = cond.monovalent_mm + divalent_to_monovalent(cond.divalent_mm, cond.dntp_mm);
    let n = oligo.len() as f64;
    ds += 0.368 * (n - 1.0) * (na.max(1e-6) / 1000.0).ln();
    let ct = if symmetric { cond.oligo_nm / 1e9 } else { cond.oligo_nm / 4e9 };
    Some(dh * 1000.0 / (ds + R * ct.ln()) - KELVIN)
}

/// Tm of an oligo that may contain IUPAC ambiguity codes: the minimum and
/// maximum over all concrete variants (at most 1024 are evaluated).
/// Gaps are ignored.
pub fn tm_range(oligo: &[u8], cond: &TmConditions) -> Option<(f64, f64)> {
    let clean = seq::ungap(oligo);
    let variants = seq::expand_ambiguities(&clean, 1024)?;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for v in &variants {
        let t = tm_exact(v, cond)?;
        lo = lo.min(t);
        hi = hi.max(t);
    }
    if lo.is_finite() { Some((lo, hi)) } else { None }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values computed with primer3's oligotm() (santalucia_auto,
    // santalucia salt correction) at 50 nM oligo, 50 mM Na, 1.5 mM Mg, 0.6 mM dNTP.
    const PRIMER3_REFERENCE: &[(&str, f64)] = &[
        ("TGCTCGATGCTAGCTAGCTC", 59.404_411_691_062_14),
        ("CCCCCATCCGATCAGGGGG", 63.840_326_030_193_52),
    ];

    /// (sequence, Tm at defaults, Tm at 250 nM / 50 mM / 3 mM Mg / 0.8 mM dNTP,
    /// Tm at 200 nM / 100 mM / no Mg), all from primer3's oligotm().
    const PRIMER3_CONDITIONS: &[(&str, f64, f64, f64)] = &[
        ("ACGTACGTACGTACGTAC", 53.713452, 57.655866, 53.466352),
        ("GGGGGGGGGGCCCCCCCCCC", 79.425446, 83.940241, 78.861769),
        ("ATATATATATATATATATAT", 30.791977, 34.690355, 30.305271),
        ("AGCTTGCAAGCTTGCA", 54.485419, 58.909188, 54.510770),
        ("TTGACCTGATCGGAATTCCA", 56.535221, 60.553669, 56.033170),
        ("GCGCAAATTTGCGCAT", 55.077538, 59.309059, 55.101802),
        ("CAGTCAGTTTGGGACCATAGG", 58.008750, 61.888345, 57.410889),
        ("ACCGTTAGCATTGACCAGTTGCAAT", 63.739508, 67.342774, 62.805888),
    ];

    #[test]
    fn matches_primer3() {
        let cond = TmConditions::default();
        for (s, expected) in PRIMER3_REFERENCE {
            let tm = tm_exact(s.as_bytes(), &cond).unwrap();
            assert!((tm - expected).abs() < 0.01, "{s}: {tm} vs {expected}");
        }
        let c2 = TmConditions { oligo_nm: 250.0, monovalent_mm: 50.0, divalent_mm: 3.0, dntp_mm: 0.8 };
        let c3 = TmConditions { oligo_nm: 200.0, monovalent_mm: 100.0, divalent_mm: 0.0, dntp_mm: 0.0 };
        for (s, a, b, d) in PRIMER3_CONDITIONS {
            for (cond, expected) in [(cond, *a), (c2, *b), (c3, *d)] {
                let tm = tm_exact(s.as_bytes(), &cond).unwrap();
                assert!((tm - expected).abs() < 0.01, "{s} {cond:?}: {tm} vs {expected}");
            }
        }
    }

    #[test]
    fn ambiguity_range() {
        let cond = TmConditions::default();
        let (lo, hi) = tm_range(b"TGCTCGATGCTAGCTAGCTS", &cond).unwrap();
        assert!(lo < hi);
        assert!(tm_exact(b"ACGN", &cond).is_none());
    }
}
