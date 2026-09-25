//! Map to reference: aligns many sequences to one reference and merges the
//! alignments into a contig, optionally refined by iterating against the
//! consensus ("fine tuning"). See `context/geneious_map_to_reference_context.md`.

pub mod align;
pub mod contig;
pub mod index;

use crate::model::{Document, Row};
use crate::seq;
use align::{AlignParams, PairAln};
use index::KmerIndex;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Sensitivity {
    Fastest,
    Medium,
    High,
    Highest,
    Custom,
}

impl Sensitivity {
    pub const ALL: [Sensitivity; 5] =
        [Sensitivity::Fastest, Sensitivity::Medium, Sensitivity::High, Sensitivity::Highest, Sensitivity::Custom];

    /// (index word length, ignore words repeated more than, min identity %)
    pub fn preset(self) -> Option<(usize, usize, f32)> {
        match self {
            Sensitivity::Fastest => Some((15, 12, 90.0)),
            Sensitivity::Medium => Some((13, 20, 80.0)),
            Sensitivity::High => Some((11, 30, 70.0)),
            Sensitivity::Highest => Some((9, 50, 60.0)),
            Sensitivity::Custom => None,
        }
    }
}

impl std::fmt::Display for Sensitivity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Sensitivity::Fastest => "Low sensitivity / fastest",
            Sensitivity::Medium => "Medium sensitivity / fast",
            Sensitivity::High => "High sensitivity / slower",
            Sensitivity::Highest => "Highest sensitivity / slow",
            Sensitivity::Custom => "Custom sensitivity",
        })
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapParams {
    pub sensitivity: Sensitivity,
    /// Index word length (k).
    pub word_len: usize,
    /// Ignore words occurring more often than this in the target (0 = keep all).
    pub max_word_repeats: usize,
    /// Minimum identity (%) over the aligned part for a sequence to be mapped.
    pub min_identity: f32,
    /// Minimum number of aligned query bases.
    pub min_aligned: usize,
    /// Fine tuning: additional iterations against the consensus.
    pub fine_tuning: usize,
    pub both_strands: bool,
    /// Clip query parts extending past the reference ends.
    pub trim_to_reference: bool,
}

impl Default for MapParams {
    fn default() -> Self {
        let (k, rep, id) = Sensitivity::Medium.preset().unwrap();
        MapParams {
            sensitivity: Sensitivity::Medium,
            word_len: k,
            max_word_repeats: rep,
            min_identity: id,
            min_aligned: 30,
            fine_tuning: 0,
            both_strands: true,
            trim_to_reference: true,
        }
    }
}

/// Progress/cancellation hooks for a mapping run.
pub struct Control<'a> {
    pub progress: &'a (dyn Fn(f32, &str) + Sync),
    pub cancel: &'a AtomicBool,
}

#[derive(Debug, Clone, Default)]
pub struct MapReport {
    pub queries: usize,
    pub mapped: usize,
    pub reverse: usize,
    pub iterations: usize,
    pub mean_identity: f32,
    pub seconds: f32,
}

pub struct MapOutput {
    /// Contig: row 0 is the reference, followed by the mapped sequences.
    pub contig: Document,
    pub unmapped: Vec<Row>,
    pub report: MapReport,
}

fn align_params(p: &MapParams) -> AlignParams {
    AlignParams { clip_overhangs: p.trim_to_reference, ..AlignParams::default() }
}

/// Maps all queries to `target`; `None` entries were not mapped.
fn map_round(
    target: &[u8],
    queries: &[Vec<u8>],
    p: &MapParams,
    ctl: &Control,
    round: usize,
    rounds: usize,
) -> Result<Vec<Option<PairAln>>, String> {
    let k = p.word_len.clamp(6, 31);
    let index = KmerIndex::new(target, k, p.max_word_repeats);
    let ap = align_params(p);
    let done = AtomicUsize::new(0);
    let n = queries.len().max(1);
    let label = if rounds > 1 { format!("Mapping (iteration {} of {rounds})", round + 1) } else { "Mapping".to_string() };
    let results: Vec<Option<PairAln>> = queries
        .par_iter()
        .map(|q| {
            if ctl.cancel.load(Ordering::Relaxed) {
                return None;
            }
            let r = align::map_query(q, target, &index, p.both_strands, &ap).filter(|a| {
                a.identity() * 100.0 >= p.min_identity && a.aligned_query_len() >= p.min_aligned.min(q.len())
            });
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d.is_multiple_of(64) || d == n {
                (ctl.progress)((round as f32 + d as f32 / n as f32) / rounds as f32, &label);
            }
            r
        })
        .collect();
    if ctl.cancel.load(Ordering::Relaxed) {
        return Err("Mapping cancelled".into());
    }
    Ok(results)
}

/// Runs map to reference. `reference` and `queries` may be gapped; gaps are removed.
pub fn map_to_reference(reference: &Row, queries: &[&Row], p: &MapParams, ctl: &Control) -> Result<MapOutput, String> {
    let t0 = std::time::Instant::now();
    let ref_seq = reference.ungapped();
    if ref_seq.len() < p.word_len {
        return Err("The reference sequence is too short".into());
    }
    let query_seqs: Vec<Vec<u8>> = queries.iter().map(|q| q.ungapped()).collect();
    let rounds = 1 + p.fine_tuning;

    let mut target = ref_seq.clone();
    let mut ref_aln = PairAln { t_start: 0, aligned: ref_seq.clone(), ..Default::default() };
    let mut results = map_round(&target, &query_seqs, p, ctl, 0, rounds)?;
    let mut iterations = 1;

    for round in 1..rounds {
        let layout = contig::Layout::new(target.len(), results.iter().flatten().chain(std::iter::once(&ref_aln)));
        let rendered: Vec<(usize, Vec<u8>)> = results.iter().flatten().map(|a| layout.render(a)).collect();
        let ref_rendered = layout.render(&ref_aln);
        let (next, next_ref) = contig::consensus_target(layout.width(), &rendered, &ref_rendered);
        if next == target {
            break;
        }
        (ctl.progress)(round as f32 / rounds as f32, "Re-mapping to consensus");
        let next_results = map_round(&next, &query_seqs, p, ctl, round, rounds)?;
        target = next;
        ref_aln = next_ref;
        results = next_results;
        iterations += 1;
    }

    (ctl.progress)(0.99, "Building contig");
    let layout = contig::Layout::new(target.len(), results.iter().flatten().chain(std::iter::once(&ref_aln)));
    let (rs, rdata) = layout.render(&ref_aln);
    let mut rows = Vec::with_capacity(queries.len() + 1);
    let mut ref_row = Row { name: reference.name.clone(), description: reference.description.clone(), start: rs, data: rdata, meta: reference.meta.clone() };
    ref_row.set_meta("Role", "Reference".into());
    rows.push(ref_row);

    let mut unmapped = Vec::new();
    let mut report = MapReport { queries: queries.len(), iterations, ..Default::default() };
    let mut id_sum = 0.0f64;
    for (q, r) in queries.iter().zip(results.iter()) {
        match r {
            Some(a) => {
                let (start, data) = layout.render(a);
                let mut row = Row { name: q.name.clone(), description: q.description.clone(), start, data, meta: q.meta.clone() };
                row.set_meta("Direction", if a.reverse { "Reverse".into() } else { "Forward".into() });
                row.set_meta("Identity %", format!("{:.2}", a.identity() * 100.0));
                row.set_meta("Mismatches", a.mismatches.to_string());
                if a.clipped_5 + a.clipped_3 > 0 {
                    row.set_meta("Trimmed bp", (a.clipped_5 + a.clipped_3).to_string());
                }
                report.mapped += 1;
                report.reverse += a.reverse as usize;
                id_sum += a.identity() as f64;
                rows.push(row);
            }
            None => unmapped.push(Row {
                name: q.name.clone(),
                description: q.description.clone(),
                start: 0,
                data: seq::ungap(&q.data),
                meta: q.meta.clone(),
            }),
        }
    }
    report.mean_identity = if report.mapped > 0 { (id_sum / report.mapped as f64 * 100.0) as f32 } else { 0.0 };
    report.seconds = t0.elapsed().as_secs_f32();
    Ok(MapOutput { contig: Document::new(rows), unmapped, report })
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

    fn row(name: &str, data: Vec<u8>) -> Row {
        Row { name: name.into(), data, ..Default::default() }
    }

    #[test]
    fn maps_and_builds_contig() {
        let reference = pseudo_random(3000, 11);
        let mut q1 = reference[200..900].to_vec();
        q1.splice(300..300, b"AAAA".iter().copied());
        let q2 = seq::revcomp(&reference[500..1500]);
        let q3 = pseudo_random(400, 77); // unrelated
        let owned = vec![row("q1", q1), row("q2", q2), row("q3", q3)];
        let queries: Vec<&Row> = owned.iter().collect();
        let cancel = AtomicBool::new(false);
        let ctl = Control { progress: &|_, _| {}, cancel: &cancel };
        for fine_tuning in [0, 2] {
            let p = MapParams { fine_tuning, ..Default::default() };
            let out = map_to_reference(&row("ref", reference.clone()), &queries, &p, &ctl).unwrap();
            assert_eq!(out.report.mapped, 2);
            assert_eq!(out.unmapped.len(), 1);
            let doc = &out.contig;
            assert_eq!(doc.rows.len(), 3);
            assert_eq!(doc.width, 3004);
            let r = &doc.rows[0];
            assert_eq!(r.ungapped(), reference);
            // q1's insertion shows as a gap block in the reference row.
            assert_eq!(r.data.iter().filter(|&&b| b == seq::GAP).count(), 4);
            assert_eq!(doc.rows[2].meta_value("Direction"), Some("Reverse"));
            // Every column of the reverse read agrees with the reference.
            let q2 = &doc.rows[2];
            for c in q2.start..q2.end() {
                let rb = r.at(c).unwrap();
                let qb = q2.at(c).unwrap();
                assert!(rb == qb || rb == seq::GAP, "column {c}");
            }
        }
    }

    #[test]
    fn fine_tuning_remaps_to_consensus() {
        // Every sample shares a SNP and an insertion relative to the reference,
        // so the consensus differs and a second round runs against it.
        let reference = pseudo_random(2000, 31);
        let mut sample = reference.clone();
        sample[700] = if sample[700] == b'A' { b'C' } else { b'A' };
        sample.splice(1200..1200, b"TTGCA".iter().copied());
        let owned: Vec<Row> = (0..6).map(|i| row(&format!("s{i}"), sample[i * 100..1500 + i * 50].to_vec())).collect();
        let queries: Vec<&Row> = owned.iter().collect();
        let cancel = AtomicBool::new(false);
        let ctl = Control { progress: &|_, _| {}, cancel: &cancel };
        let p = MapParams { fine_tuning: 3, ..Default::default() };
        let out = map_to_reference(&row("ref", reference.clone()), &queries, &p, &ctl).unwrap();
        assert_eq!(out.report.mapped, 6);
        assert_eq!(out.report.iterations, 2, "second round runs, third converges");
        let doc = &out.contig;
        assert_eq!(doc.rows[0].ungapped(), reference);
        // The insertion occupies the same five columns in every sample.
        let gap_cols: Vec<usize> = (0..doc.width).filter(|&c| doc.rows[0].at(c) == Some(seq::GAP)).collect();
        assert_eq!(gap_cols.len(), 5);
        for r in &doc.rows[1..] {
            let ins: Vec<u8> = gap_cols.iter().filter_map(|&c| r.at(c)).collect();
            if !ins.is_empty() {
                assert_eq!(ins, b"TTGCA".to_vec());
            }
            assert_eq!(r.ungapped(), sample[..].windows(r.ungapped().len()).find(|w| *w == r.ungapped().as_slice()).unwrap().to_vec());
        }
    }

    /// Scale check: `cargo test --release -- --ignored bench_genomes --nocapture`
    #[test]
    #[ignore]
    fn bench_genomes() {
        let reference = pseudo_random(35_000, 21);
        let noise = pseudo_random(35_000 * 8, 5);
        let queries: Vec<Row> = (0..2000)
            .map(|i| {
                let mut q = Vec::with_capacity(35_100);
                for (j, &b) in reference.iter().enumerate() {
                    let n = noise[(j * 7 + i * 13) % noise.len()];
                    if (j + i) % 2003 == 0 {
                        continue;
                    }
                    q.push(if n == b'A' && (j + i) % 11 == 0 { if b == b'C' { b'T' } else { b'C' } } else { b });
                    if (j + 3 * i) % 4001 == 0 {
                        q.push(b'G');
                    }
                }
                row(&format!("g{i}"), q)
            })
            .collect();
        let cancel = AtomicBool::new(false);
        let ctl = Control { progress: &|_, _| {}, cancel: &cancel };
        let t = std::time::Instant::now();
        let refs: Vec<&Row> = queries.iter().collect();
        let out = map_to_reference(&row("ref", reference), &refs, &MapParams::default(), &ctl).unwrap();
        println!("mapped {} of {} in {:.1}s, width {}, mean identity {:.2}", out.report.mapped, queries.len(), t.elapsed().as_secs_f32(), out.contig.width, out.report.mean_identity);
        assert_eq!(out.report.mapped, 2000);

        // Fine tuning on a subset: re-mapping to the consensus must keep everything mapped.
        let subset: Vec<&Row> = queries.iter().take(300).collect();
        let t = std::time::Instant::now();
        let p = MapParams { fine_tuning: 2, ..MapParams::default() };
        let ft = map_to_reference(&out.contig.rows[0], &subset, &p, &ctl).unwrap();
        println!("fine tuning: {} iterations, {} mapped in {:.1}s", ft.report.iterations, ft.report.mapped, t.elapsed().as_secs_f32());
        assert_eq!(ft.report.mapped, 300);
        assert_eq!(ft.contig.rows[0].ungapped(), out.contig.rows[0].ungapped());
    }
}
