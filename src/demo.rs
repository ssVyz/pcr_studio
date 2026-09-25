//! Synthetic demo data: a reference genome and many related genomes from a
//! few clades (clade-specific substitutions and indels, per-sequence noise,
//! partial sequences, Ns and reverse-complemented entries).

use crate::model::Row;
use crate::seq;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn chance(&mut self, p: f64) -> bool {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64 <= p
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
}

enum Edit {
    Sub(usize, u8),
    Del(usize, usize),
    Ins(usize, Vec<u8>),
}

fn apply(reference: &[u8], edits: &[Edit]) -> Vec<u8> {
    let mut subs = std::collections::HashMap::new();
    let mut dels = std::collections::HashMap::new();
    let mut ins: std::collections::HashMap<usize, &Vec<u8>> = std::collections::HashMap::new();
    for e in edits {
        match e {
            Edit::Sub(p, b) => {
                subs.insert(*p, *b);
            }
            Edit::Del(p, l) => {
                dels.insert(*p, *l);
            }
            Edit::Ins(p, s) => {
                ins.insert(*p, s);
            }
        }
    }
    let mut out = Vec::with_capacity(reference.len() + 64);
    let mut i = 0;
    while i < reference.len() {
        if let Some(s) = ins.get(&i) {
            out.extend_from_slice(s);
        }
        if let Some(&l) = dels.get(&i) {
            i += l;
            continue;
        }
        out.push(*subs.get(&i).unwrap_or(&reference[i]));
        i += 1;
    }
    out
}

/// Returns (reference, genomes). Genome metadata: clade, country, year.
pub fn generate(len: usize, n: usize, seed: u64) -> (Row, Vec<Row>) {
    let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    // Mildly structured reference: GC-rich and AT-rich stretches.
    let mut reference = Vec::with_capacity(len);
    let mut gc_rich = false;
    for i in 0..len {
        if i % 800 == 0 {
            gc_rich = rng.chance(0.5);
        }
        let b = if gc_rich && rng.chance(0.3) { if rng.chance(0.5) { b'G' } else { b'C' } } else { rng.base() };
        reference.push(b);
    }
    let clades = ["A", "B", "C", "D", "E", "F"];
    let clade_weights = [0.35, 0.25, 0.15, 0.12, 0.08, 0.05];
    let countries = ["DE", "US", "BR", "IN", "ZA", "CN", "FR", "KE"];
    let mut clade_edits: Vec<Vec<Edit>> = Vec::new();
    for (ci, _) in clades.iter().enumerate() {
        let mut e = Vec::new();
        let divergence = 0.004 + 0.004 * ci as f64;
        for p in 0..len {
            if rng.chance(divergence) {
                let mut b = rng.base();
                while b == reference[p] {
                    b = rng.base();
                }
                e.push(Edit::Sub(p, b));
            }
        }
        for _ in 0..(2 + ci) {
            let p = 500 + rng.below((len - 1000) as u64) as usize;
            if rng.chance(0.5) {
                e.push(Edit::Del(p, 1 + rng.below(9) as usize));
            } else {
                let l = 1 + rng.below(12) as usize;
                e.push(Edit::Ins(p, (0..l).map(|_| rng.base()).collect()));
            }
        }
        clade_edits.push(e);
    }
    let mut genomes = Vec::with_capacity(n);
    for i in 0..n {
        let r = rng.below(1000) as f64 / 1000.0;
        let mut acc = 0.0;
        let mut ci = clades.len() - 1;
        for (k, w) in clade_weights.iter().enumerate() {
            acc += w;
            if r < acc {
                ci = k;
                break;
            }
        }
        let mut g = apply(&reference, &clade_edits[ci]);
        for b in g.iter_mut() {
            if rng.chance(0.0008) {
                *b = rng.base();
            }
        }
        if rng.chance(0.03) {
            let p = rng.below(g.len() as u64) as usize;
            let l = (20 + rng.below(200) as usize).min(g.len() - p);
            for b in &mut g[p..p + l] {
                *b = b'N';
            }
        }
        if rng.chance(0.15) {
            let start = rng.below(g.len() as u64 / 3) as usize;
            let end = g.len() - rng.below(g.len() as u64 / 3) as usize;
            g = g[start..end.max(start + 1000).min(g.len())].to_vec();
        }
        if rng.chance(0.1) {
            g = seq::revcomp(&g);
        }
        let country = countries[rng.below(countries.len() as u64) as usize];
        let year = 2005 + rng.below(20);
        let name = format!("GEN{:05}.1", i + 1);
        let mut row = Row { name, description: format!("[clade={}] [country={country}] [year={year}]", clades[ci]), start: 0, data: g, meta: Vec::new() };
        row.set_meta("clade", clades[ci].to_string());
        row.set_meta("country", country.to_string());
        row.set_meta("year", year.to_string());
        genomes.push(row);
    }
    let mut r = Row { name: "REF_genome".into(), description: "synthetic reference".into(), start: 0, data: reference, meta: Vec::new() };
    r.set_meta("clade", "reference".into());
    (r, genomes)
}

#[cfg(test)]
mod tests {
    /// Writes demo FASTA files: `cargo test --release -- --ignored write_demo_fasta`
    /// (directory from PCR_STUDIO_DEMO_DIR, default: system temp).
    #[test]
    #[ignore]
    fn write_demo_fasta() {
        let dir = std::env::var("PCR_STUDIO_DEMO_DIR").map(std::path::PathBuf::from).unwrap_or_else(|_| std::env::temp_dir());
        let (r, g) = super::generate(35_000, 2_000, 7);
        let rec = |r: &crate::model::Row| crate::fasta::Record { name: r.name.clone(), description: r.description.clone(), seq: r.data.clone() };
        crate::fasta::write_fasta(&dir.join("demo_reference.fasta"), &[rec(&r)]).unwrap();
        let recs: Vec<_> = g.iter().map(rec).collect();
        crate::fasta::write_fasta(&dir.join("demo_genomes.fasta"), &recs).unwrap();
        println!("wrote demo files to {}", dir.display());
    }

    #[test]
    fn deterministic() {
        let (r1, g1) = super::generate(3000, 20, 1);
        let (r2, g2) = super::generate(3000, 20, 1);
        assert_eq!(r1.data, r2.data);
        assert_eq!(g1[5].data, g2[5].data);
        assert_eq!(g1.len(), 20);
    }
}
