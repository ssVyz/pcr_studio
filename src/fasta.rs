//! FASTA import/export and sequence metadata helpers.

use crate::model::{DocKind, Document, Row};
use crate::seq;
use bio::io::fasta::FastaRead;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

/// Result of reading a FASTA file.
pub struct Imported {
    pub name: String,
    pub kind: DocKind,
    pub document: Document,
}

/// Reads a FASTA file (also gzip-compressed). Headers are split into id and
/// description; `key=value` tokens (optionally in square brackets) in the
/// description become row metadata. Gapped files, or files whose sequences all
/// have the same length, are imported as alignments.
pub fn read_fasta(path: &Path, progress: &dyn Fn(f32)) -> Result<Imported, String> {
    let total = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0).max(1);
    let read = Arc::new(AtomicU64::new(0));
    let mut reader = crate::linked::open(path, read.clone(), None)?;
    let first = reader.fill_buf().map_err(|e| format!("Cannot read {}: {e}", path.display()))?.iter().find(|b| !b.is_ascii_whitespace()).copied();
    if first == Some(b'@') {
        return Err(format!("{} is a FASTQ file. FASTQ files are linked, not imported: use File → Import FASTQ (linked)…", path.display()));
    }
    let mut last_report = 0u64;
    let (rows, any_gap, _) = read_rows(reader, usize::MAX, u64::MAX, &mut || {
        let done = read.load(Ordering::Relaxed);
        if done > last_report {
            last_report = done;
            progress(done as f32 / total as f32);
        }
    })
    .map_err(|e| format!("Invalid FASTA in {}: {e}", path.display()))?;
    if rows.is_empty() {
        return Err(format!("{} contains no FASTA records", path.display()));
    }
    let kind = alignment_kind(any_gap, &rows.iter().map(Row::end).collect::<Vec<_>>());
    Ok(Imported { name: crate::linked::stem(path), kind, document: Document::new(rows) })
}

/// Reads FASTA records as rows (see [`read_fasta`]) until `max_rows` rows or
/// at least `max_bases` bases are read; `tick` runs every 256 records.
/// Returns the rows, whether any has a gap, and whether records are left.
pub fn read_rows(reader: impl BufRead, max_rows: usize, max_bases: u64, tick: &mut dyn FnMut()) -> Result<(Vec<Row>, bool, bool), String> {
    let mut reader = bio::io::fasta::Reader::from_bufread(reader);
    let mut rec = bio::io::fasta::Record::new();
    let mut rows = Vec::new();
    let mut any_gap = false;
    let mut bases = 0u64;
    loop {
        if rows.len() >= max_rows || bases >= max_bases {
            reader.read(&mut rec).map_err(|e| e.to_string())?;
            return Ok((rows, any_gap, !rec.is_empty()));
        }
        reader.read(&mut rec).map_err(|e| e.to_string())?;
        if rec.is_empty() {
            return Ok((rows, any_gap, false));
        }
        let data = seq::normalize(rec.seq());
        any_gap |= data.contains(&seq::GAP);
        let desc = rec.desc().unwrap_or("").to_string();
        let mut row = Row::from_gapped(rec.id().to_string(), desc.clone(), data);
        row.meta = parse_description_meta(&desc);
        bases += row.ungapped_len() as u64;
        rows.push(row);
        if rows.len() % 256 == 0 {
            tick();
        }
    }
}

/// Aligned FASTA: gapped, or several sequences that all share one length.
pub fn alignment_kind(any_gap: bool, lengths: &[usize]) -> DocKind {
    let same_len = lengths.len() > 1 && lengths.iter().all(|&l| l == lengths[0]);
    if any_gap || same_len { DocKind::Alignment } else { DocKind::Sequences }
}

/// Document kind for pasted records, decided as on FASTA import.
pub fn pasted_kind(records: &[PastedRecord]) -> DocKind {
    let any_gap = records.iter().any(|r| r.seq.contains(&seq::GAP));
    alignment_kind(any_gap, &records.iter().map(|r| r.seq.len()).collect::<Vec<_>>())
}

/// One record of pasted or typed sequence text.
#[derive(Debug, Clone, PartialEq)]
pub struct PastedRecord {
    /// FASTA header id (first word after `>`), `None` without a header.
    pub id: Option<String>,
    /// Rest of the header line.
    pub description: String,
    /// Sequence in the internal alphabet.
    pub seq: Vec<u8>,
    /// Characters that are not nucleotide codes and became `N`.
    pub unknown: usize,
}

/// Parses pasted text: a plain sequence or FASTA records. Sequence characters
/// are normalized as on FASTA import: whitespace and digits (GenBank/EMBL
/// listings) are dropped, `U` becomes `T`, gaps are kept and anything unknown
/// becomes `N`. Text before the first header is a record without a header.
pub fn parse_pasted(text: &str) -> Vec<PastedRecord> {
    let mut out: Vec<PastedRecord> = Vec::new();
    let empty = || PastedRecord { id: None, description: String::new(), seq: Vec::new(), unknown: 0 };
    for line in text.lines() {
        if let Some(header) = line.trim_start().strip_prefix('>') {
            let header = header.trim();
            let (id, desc) = header.split_once(char::is_whitespace).unwrap_or((header, ""));
            out.push(PastedRecord { id: (!id.is_empty()).then(|| id.to_string()), description: desc.trim().to_string(), ..empty() });
            continue;
        }
        let mut chars = line.chars().filter(|c| !c.is_whitespace()).peekable();
        if chars.peek().is_none() {
            continue;
        }
        if out.is_empty() {
            out.push(empty());
        }
        let rec = out.last_mut().unwrap();
        for c in chars {
            // Per character, so non-ASCII text becomes one N, not one per byte.
            let Some(b) = (if c.is_ascii() { seq::normalize_byte(c as u8) } else { Some(b'N') }) else { continue };
            if b == b'N' && !c.eq_ignore_ascii_case(&'n') {
                rec.unknown += 1;
            }
            rec.seq.push(b);
        }
    }
    out
}

/// Builds a document from pasted records. A single record gets `name` and
/// `description`; several records keep their header names ("Sequence N"
/// without a header) and descriptions. `key=value` tokens in descriptions
/// become metadata, terminal gaps "no coverage", as on FASTA import.
pub fn document_from_pasted(records: &[PastedRecord], name: &str, description: &str) -> (DocKind, Document) {
    let single = records.len() == 1;
    let rows: Vec<Row> = records
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let (name, desc) = if single {
                (name.to_string(), description.to_string())
            } else {
                (r.id.clone().unwrap_or_else(|| format!("Sequence {}", i + 1)), r.description.clone())
            };
            let mut row = Row::from_gapped(name, desc.clone(), r.seq.clone());
            row.meta = parse_description_meta(&desc);
            row
        })
        .collect();
    (pasted_kind(records), Document::new(rows))
}

/// Extracts `key=value` pairs from a FASTA description, e.g.
/// `[organism=Human adenovirus 5] [country=USA]` or `host=human year=2019`.
pub fn parse_description_meta(desc: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    if desc.contains('[') {
        let mut rest = desc;
        while let Some(open) = rest.find('[') {
            let Some(close) = rest[open..].find(']') else { break };
            let inner = &rest[open + 1..open + close];
            if let Some((k, v)) = inner.split_once('=') {
                let k = k.trim();
                if !k.is_empty() {
                    out.push((k.to_string(), v.trim().to_string()));
                }
            }
            rest = &rest[open + close + 1..];
        }
    } else {
        for tok in desc.split_whitespace() {
            if let Some((k, v)) = tok.split_once('=')
                && !k.is_empty() && !v.is_empty() {
                    out.push((k.to_string(), v.to_string()));
                }
        }
    }
    out
}

/// One FASTA record to be written.
pub struct Record {
    pub name: String,
    pub description: String,
    pub seq: Vec<u8>,
}

/// Writes records with bio's FASTA writer.
pub fn write_fasta(path: &Path, records: &[Record]) -> Result<(), String> {
    let file = std::fs::File::create(path).map_err(|e| format!("Cannot create {}: {e}", path.display()))?;
    let mut writer = bio::io::fasta::Writer::new(BufWriter::new(file));
    for r in records {
        let desc = if r.description.is_empty() { None } else { Some(r.description.as_str()) };
        writer
            .write(&r.name, desc, &r.seq)
            .map_err(|e| format!("Write error: {e}"))?;
    }
    writer.flush().map_err(|e| format!("Write error: {e}"))
}

/// Parsed metadata table: first column = sequence name, other columns = keys.
pub struct MetaTable {
    pub keys: Vec<String>,
    pub rows: Vec<(String, Vec<String>)>,
}

/// Reads a CSV/TSV/semicolon-separated table with a header line.
pub fn read_meta_table(path: &Path) -> Result<MetaTable, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = loop {
        match lines.next() {
            None => return Err("The table is empty".into()),
            Some(l) => {
                let l = l.map_err(|e| e.to_string())?;
                if !l.trim().is_empty() {
                    break l;
                }
            }
        }
    };
    let header = header.trim_start_matches('\u{feff}').to_string();
    let delim = if header.contains('\t') {
        '\t'
    } else if header.matches(';').count() > header.matches(',').count() {
        ';'
    } else {
        ','
    };
    let head = split_delimited(&header, delim);
    if head.len() < 2 {
        return Err("The table needs a name column and at least one metadata column".into());
    }
    let keys: Vec<String> = head[1..].iter().map(|k| k.trim().to_string()).collect();
    let mut rows = Vec::new();
    for l in lines {
        let l = l.map_err(|e| e.to_string())?;
        if l.trim().is_empty() {
            continue;
        }
        let fields = split_delimited(&l, delim);
        let name = fields[0].trim().to_string();
        let mut vals: Vec<String> = fields[1..].iter().map(|v| v.trim().to_string()).collect();
        vals.resize(keys.len(), String::new());
        rows.push((name, vals));
    }
    Ok(MetaTable { keys, rows })
}

/// Splits one line, honoring double quotes.
fn split_delimited(line: &str, delim: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '"' {
            if quoted && chars.peek() == Some(&'"') {
                cur.push('"');
                chars.next();
            } else {
                quoted = !quoted;
            }
        } else if c == delim && !quoted {
            out.push(std::mem::take(&mut cur));
        } else {
            cur.push(c);
        }
    }
    out.push(cur);
    out
}

/// Splits sequence names at a delimiter into metadata fields `field_1..n`.
pub fn name_fields(name: &str, delim: char) -> Vec<String> {
    name.split(delim).map(|s| s.to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("pcr_studio_test_{}_{name}", std::process::id()))
    }

    #[test]
    fn roundtrip_aligned() {
        let p = tmp("aln.fasta");
        std::fs::write(&p, ">a [host=human] [year=2019]\n--ACGT\nAC\n>b host=pig\nACACGTAC\n").unwrap();
        let imp = read_fasta(&p, &|_| {}).unwrap();
        assert_eq!(imp.kind, DocKind::Alignment);
        let r = &imp.document.rows[0];
        assert_eq!(r.start, 2);
        assert_eq!(r.data, b"ACGTAC");
        assert_eq!(r.meta_value("host"), Some("human"));
        assert_eq!(imp.document.rows[1].meta_value("host"), Some("pig"));
        let out = tmp("out.fasta");
        write_fasta(&out, &[Record { name: "x".into(), description: "d".into(), seq: b"AC-GT".to_vec() }]).unwrap();
        let text = std::fs::read_to_string(&out).unwrap();
        assert_eq!(text, ">x d\nAC-GT\n");
        let _ = std::fs::remove_file(p);
        let _ = std::fs::remove_file(out);
    }

    #[test]
    fn unaligned_detection() {
        let p = tmp("seqs.fasta");
        std::fs::write(&p, ">a\nACGT\n>b\nACGTTT\n").unwrap();
        let imp = read_fasta(&p, &|_| {}).unwrap();
        assert_eq!(imp.kind, DocKind::Sequences);
        let _ = std::fs::remove_file(p);
    }

    #[test]
    fn gzip_limits_and_fastq() {
        let p = tmp("seqs.fasta.gz");
        let mut e = flate2::write::GzEncoder::new(std::fs::File::create(&p).unwrap(), flate2::Compression::fast());
        std::io::Write::write_all(&mut e, b">a\nACGT\n>b\nACGTTT\n>c\nA\n").unwrap();
        e.finish().unwrap();
        let imp = read_fasta(&p, &|_| {}).unwrap();
        assert!(imp.name.ends_with("_seqs"));
        assert_eq!(imp.document.rows[1].data, b"ACGTTT");
        let rows = |max_rows, max_bases| read_rows(crate::linked::open_plain(&p).unwrap(), max_rows, max_bases, &mut || {}).unwrap();
        let (r, _, more) = rows(2, u64::MAX);
        assert_eq!((r.len(), more), (2, true));
        // Stops once the bases reach the limit (4 + 6 ≥ 5).
        let (r, _, more) = rows(10, 5);
        assert_eq!((r.len(), more), (2, true));
        let (r, any_gap, more) = rows(3, u64::MAX);
        assert_eq!((r.len(), any_gap, more), (3, false, false));
        let q = tmp("reads.fastq");
        std::fs::write(&q, "@r\nACGT\n+\nIIII\n").unwrap();
        assert!(read_fasta(&q, &|_| {}).err().unwrap().contains("is a FASTQ file"));
        let _ = std::fs::remove_file(p);
        let _ = std::fs::remove_file(q);
    }

    #[test]
    fn pasted_plain_and_genbank() {
        // GenBank ORIGIN listing: numbers and spaces are ignored, case and U normalized.
        let recs = parse_pasted("        1 acgtacgtac gtnnu\n       16 ACG\u{a0}T\r\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, None);
        assert_eq!(recs[0].seq, b"ACGTACGTACGTNNTACGT");
        assert_eq!(recs[0].unknown, 0);
        // Unknown characters become N and are counted; non-ASCII counts once.
        let recs = parse_pasted("AC-GJ?–T");
        assert_eq!(recs[0].seq, b"AC-GNNNT");
        assert_eq!(recs[0].unknown, 3);
        assert!(parse_pasted(" \n\n").is_empty());
    }

    #[test]
    fn pasted_fasta_records() {
        let recs = parse_pasted("\n>s1 [host=human] demo\nACGT\nAC\n> \nGG\n>s3\n");
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].id.as_deref(), Some("s1"));
        assert_eq!(recs[0].description, "[host=human] demo");
        assert_eq!(recs[0].seq, b"ACGTAC");
        assert_eq!(recs[1].id, None);
        assert!(recs[2].seq.is_empty());
        let (kind, doc) = document_from_pasted(&recs[..2], "ignored", "");
        assert_eq!(kind, DocKind::Sequences);
        assert_eq!(doc.rows[0].name, "s1");
        assert_eq!(doc.rows[0].meta_value("host"), Some("human"));
        assert_eq!(doc.rows[1].name, "Sequence 2");
        // A single record takes the given name and description.
        let (kind, doc) = document_from_pasted(&recs[..1], "My seq", "year=2020");
        assert_eq!(kind, DocKind::Sequences);
        assert_eq!((doc.rows[0].name.as_str(), doc.rows[0].meta_value("year")), ("My seq", Some("2020")));
        // Gaps or equal lengths make an alignment, as on import.
        let (kind, doc) = document_from_pasted(&parse_pasted("--AC-GT"), "g", "");
        assert_eq!((kind, doc.rows[0].start, doc.width), (DocKind::Alignment, 2, 7));
        let (kind, _) = document_from_pasted(&parse_pasted(">a\nACGT\n>b\nTTTT"), "x", "");
        assert_eq!(kind, DocKind::Alignment);
    }

    #[test]
    fn table_parsing() {
        let p = tmp("meta.csv");
        std::fs::write(&p, "name,genotype,\"place, country\"\ns1,A,\"Hamburg, DE\"\ns2,B\n").unwrap();
        let t = read_meta_table(&p).unwrap();
        assert_eq!(t.keys, vec!["genotype", "place, country"]);
        assert_eq!(t.rows[0].1[1], "Hamburg, DE");
        assert_eq!(t.rows[1].1, vec!["B".to_string(), String::new()]);
        let _ = std::fs::remove_file(p);
    }
}
