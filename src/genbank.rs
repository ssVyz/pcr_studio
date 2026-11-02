//! GenBank flat files (NCBI nucleotide records): header, features with their
//! locations and qualifiers, and the sequence.
//!
//! Only what PCR Studio uses is read: LOCUS (name, length), DEFINITION,
//! ACCESSION, VERSION, the feature table and ORIGIN; other sections are
//! skipped. Files may hold several records.

/// One record (`LOCUS` … `//`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Record {
    pub locus: String,
    /// Length from the LOCUS line (bp); without one, the length of the sequence.
    pub length: usize,
    pub definition: String,
    /// Accession with version ("MN908947.3"), else the accession, else the locus name.
    pub version: String,
    pub features: Vec<Feature>,
    /// Upper case; empty if the record has no ORIGIN section (or only its header was read).
    pub sequence: Vec<u8>,
}

/// A feature of the feature table.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Feature {
    /// Feature key: "gene", "CDS", "mat_peptide", …
    pub key: String,
    /// Location as written, without white space ("complement(join(1..20,30..40))").
    pub location: String,
    /// Qualifiers in record order, values without quotes ("gene" → "S").
    pub qualifiers: Vec<(String, String)>,
}

/// One part of a location: 1-based positions `start..=end` of the record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
    pub complement: bool,
}

/// A parsed location.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Location {
    /// Parts on this record, sorted by position.
    pub spans: Vec<Span>,
    /// Parts on other records (`J00194.1:100..202`), which are left out.
    pub remote: usize,
}

impl Location {
    /// +1 on the plus strand, -1 on the minus strand (every part
    /// complemented), 0 when mixed or without parts.
    pub fn strand(&self) -> i8 {
        match (self.spans.iter().all(|s| s.complement), self.spans.iter().any(|s| s.complement)) {
            _ if self.spans.is_empty() => 0,
            (true, _) => -1,
            (false, false) => 1,
            (false, true) => 0,
        }
    }
}

impl Feature {
    pub fn qualifier(&self, key: &str) -> Option<&str> {
        self.qualifiers.iter().find(|(k, v)| k == key && !v.trim().is_empty()).map(|(_, v)| v.as_str())
    }

    /// A short name from the qualifiers: the gene for genes, the product for
    /// coding sequences, peptides and RNAs, otherwise the first descriptive
    /// qualifier; the feature key when there is none.
    pub fn name(&self) -> String {
        const PRODUCT_FIRST: [&str; 13] = [
            "CDS",
            "mat_peptide",
            "sig_peptide",
            "transit_peptide",
            "propeptide",
            "mRNA",
            "rRNA",
            "tRNA",
            "ncRNA",
            "misc_RNA",
            "tmRNA",
            "precursor_RNA",
            "prim_transcript",
        ];
        let order: &[&str] = match self.key.as_str() {
            "gene" => &["gene", "locus_tag", "gene_synonym", "note"],
            "source" => &["organism", "strain", "isolate", "note"],
            k if PRODUCT_FIRST.contains(&k) => &["product", "gene", "locus_tag", "note"],
            _ => &[
                "gene",
                "product",
                "label",
                "standard_name",
                "regulatory_class",
                "rpt_type",
                "mobile_element_type",
                "bound_moiety",
                "note",
                "locus_tag",
            ],
        };
        order.iter().find_map(|k| self.qualifier(k)).map_or_else(|| self.key.clone(), |v| shorten(v, 60))
    }
}

/// First line of `s`, at most `max` characters (with “…” when cut).
fn shorten(s: &str, max: usize) -> String {
    let s = s.lines().next().unwrap_or("").trim();
    if s.chars().count() <= max {
        return s.to_string();
    }
    let mut t: String = s.chars().take(max - 1).collect();
    t.push('…');
    t
}

/// Where the parser is in a record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Section {
    Definition,
    Features,
    Origin,
    Other,
}

/// Collects the lines of one record.
struct Builder {
    rec: Record,
    accession: String,
    section: Section,
    feature: Option<Feature>,
    /// The last qualifier's value is a quoted text that continues on the next line.
    open_quote: bool,
}

impl Builder {
    fn new(locus_line: &str) -> Builder {
        let toks: Vec<&str> = locus_line.split_whitespace().skip(1).collect();
        let length = toks.windows(2).find(|w| w[1] == "bp" || w[1] == "aa").and_then(|w| w[0].parse().ok()).unwrap_or(0);
        let rec = Record { locus: toks.first().map(|s| s.to_string()).unwrap_or_default(), length, ..Record::default() };
        Builder { rec, accession: String::new(), section: Section::Other, feature: None, open_quote: false }
    }

    fn flush_feature(&mut self) {
        if let Some(mut f) = self.feature.take() {
            for (_, v) in &mut f.qualifiers {
                *v = unquote(v);
            }
            self.rec.features.push(f);
        }
        self.open_quote = false;
    }

    fn finish(mut self) -> Record {
        self.flush_feature();
        let mut r = self.rec;
        if r.version.is_empty() {
            r.version = if self.accession.is_empty() { r.locus.clone() } else { self.accession };
        }
        if r.length == 0 {
            r.length = r.sequence.len();
        }
        r
    }

    /// A line starting in column 1 (a new section).
    fn keyword(&mut self, line: &str) {
        self.flush_feature();
        let (kw, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        let rest = rest.trim();
        self.section = match kw {
            "DEFINITION" => {
                self.rec.definition = rest.to_string();
                Section::Definition
            }
            "FEATURES" => Section::Features,
            "ORIGIN" => Section::Origin,
            "ACCESSION" => {
                self.accession = rest.split_whitespace().next().unwrap_or("").to_string();
                Section::Other
            }
            "VERSION" => {
                self.rec.version = rest.split_whitespace().next().unwrap_or("").to_string();
                Section::Other
            }
            _ => Section::Other,
        };
    }

    /// An indented line of the current section.
    fn continuation(&mut self, line: &str) {
        let content = line.trim();
        match self.section {
            Section::Definition => {
                if !content.is_empty() {
                    self.rec.definition.push(' ');
                    self.rec.definition.push_str(content);
                }
            }
            Section::Origin => self.rec.sequence.extend(content.bytes().filter(u8::is_ascii_alphabetic).map(|b| b.to_ascii_uppercase())),
            Section::Features => self.feature_line(line, content),
            Section::Other => {}
        }
    }

    fn feature_line(&mut self, line: &str, content: &str) {
        if content.is_empty() {
            return;
        }
        let indent = line.len() - line.trim_start().len();
        // Feature keys start in column 6, qualifiers and continued locations in column 22.
        if !self.open_quote && indent < 21 {
            self.flush_feature();
            let (key, location) = content.split_once(char::is_whitespace).unwrap_or((content, ""));
            self.feature = Some(Feature { key: key.to_string(), location: strip_space(location), qualifiers: Vec::new() });
            return;
        }
        let Some(f) = self.feature.as_mut() else { return };
        if !self.open_quote && content.starts_with('/') {
            let (key, value) = match content[1..].split_once('=') {
                Some((k, v)) => (k.to_string(), v.to_string()),
                None => (content[1..].to_string(), String::new()),
            };
            self.open_quote = value.starts_with('"') && value.matches('"').count() % 2 == 1;
            f.qualifiers.push((key, value));
        } else if let Some((key, value)) = f.qualifiers.last_mut() {
            // Translations are wrapped without spaces; free text with one.
            if key != "translation" {
                value.push(' ');
            }
            value.push_str(content);
            if self.open_quote {
                self.open_quote = value.matches('"').count() % 2 == 1;
            }
        } else {
            f.location.push_str(&strip_space(content));
        }
    }
}

fn strip_space(s: &str) -> String {
    s.chars().filter(|c| !c.is_whitespace()).collect()
}

/// A qualifier value without its quotes (`""` inside stands for `"`).
fn unquote(v: &str) -> String {
    let v = v.trim();
    match v.strip_prefix('"') {
        Some(inner) => inner.strip_suffix('"').unwrap_or(inner).replace("\"\"", "\""),
        None => v.to_string(),
    }
}

/// Reads all records of a GenBank flat file. Text before the first LOCUS
/// line is ignored; a missing `//` ends the record at the next LOCUS line or
/// the end of the text.
pub fn parse(text: &str) -> Result<Vec<Record>, String> {
    let mut records = Vec::new();
    let mut cur: Option<Builder> = None;
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if line.starts_with("LOCUS") {
            if let Some(b) = cur.take() {
                records.push(b.finish());
            }
            cur = Some(Builder::new(line));
            continue;
        }
        let Some(b) = cur.as_mut() else { continue };
        if line.starts_with("//") {
            records.push(cur.take().unwrap().finish());
        } else if line.starts_with(|c: char| !c.is_whitespace()) {
            b.keyword(line);
        } else {
            b.continuation(line);
        }
    }
    if let Some(b) = cur {
        records.push(b.finish());
    }
    if records.is_empty() {
        return Err("This is not a GenBank record (no LOCUS line).".into());
    }
    Ok(records)
}

/// Parses a feature location: single bases, ranges (`<1..>200`),
/// `complement(…)`, `join(…)`, `order(…)`, sites between bases (`12^13`) and
/// parts on other records (counted, not kept).
pub fn parse_location(s: &str) -> Result<Location, String> {
    let bytes: Vec<u8> = s.bytes().filter(|b| !b.is_ascii_whitespace()).collect();
    let mut p = LocParser { s: &bytes, i: 0, loc: Location::default() };
    p.loc(false)?;
    if p.i != bytes.len() {
        return Err(format!("unexpected text at “{}”", String::from_utf8_lossy(&bytes[p.i..])));
    }
    p.loc.spans.sort_by_key(|s| (s.start, s.end));
    Ok(p.loc)
}

struct LocParser<'a> {
    s: &'a [u8],
    i: usize,
    loc: Location,
}

impl LocParser<'_> {
    fn eat(&mut self, token: &str) -> bool {
        if self.s[self.i..].starts_with(token.as_bytes()) {
            self.i += token.len();
            true
        } else {
            false
        }
    }

    fn expect(&mut self, token: &str) -> Result<(), String> {
        if self.eat(token) { Ok(()) } else { Err(format!("“{token}” expected")) }
    }

    fn loc(&mut self, complement: bool) -> Result<(), String> {
        if self.eat("complement(") {
            self.loc(!complement)?;
            return self.expect(")");
        }
        for op in ["join(", "order(", "bond(", "one-of("] {
            if self.eat(op) {
                loop {
                    self.loc(complement)?;
                    if !self.eat(",") {
                        break;
                    }
                }
                return self.expect(")");
            }
        }
        if self.eat("gap(") {
            while self.i < self.s.len() && self.s[self.i] != b')' {
                self.i += 1;
            }
            return self.expect(")");
        }
        self.span(complement)
    }

    fn span(&mut self, complement: bool) -> Result<(), String> {
        // A part on another record: "accession.version:range".
        let end = self.s[self.i..].iter().position(|&b| matches!(b, b',' | b')')).map_or(self.s.len(), |p| self.i + p);
        if self.s[self.i..end].contains(&b':') {
            self.i = end;
            self.loc.remote += 1;
            return Ok(());
        }
        let a = self.position(false)?;
        let b = if self.eat("..") {
            self.position(true)?
        } else if self.eat("^") {
            // A site between two bases: shown on the first.
            self.number()?;
            a
        } else if self.s.get(self.i) == Some(&b'.') && self.s.get(self.i + 1).is_some_and(u8::is_ascii_digit) {
            // Old notation "102.110": a single base somewhere in the range.
            self.i += 1;
            self.number()?
        } else {
            a
        };
        if a == 0 || b < a {
            return Err(format!("invalid range {a}..{b}"));
        }
        self.loc.spans.push(Span { start: a, end: b, complement });
        Ok(())
    }

    /// A position with an optional `<` / `>` (partial) or `(a.b)` (one of a
    /// range: its first base at a start, its last at an end).
    fn position(&mut self, end: bool) -> Result<usize, String> {
        let _ = self.eat("<") || self.eat(">");
        if self.eat("(") {
            let a = self.number()?;
            self.expect(".")?;
            let b = self.number()?;
            self.expect(")")?;
            return Ok(if end { b } else { a });
        }
        self.number()
    }

    fn number(&mut self) -> Result<usize, String> {
        let start = self.i;
        while self.i < self.s.len() && self.s[self.i].is_ascii_digit() {
            self.i += 1;
        }
        std::str::from_utf8(&self.s[start..self.i])
            .ok()
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| format!("position expected at “{}”", String::from_utf8_lossy(&self.s[start..])))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
LOCUS       TEST01                   120 bp    RNA     linear   VRL 18-MAR-2020
DEFINITION  Test virus isolate X, complete
            genome.
ACCESSION   TEST01
VERSION     TEST01.2
KEYWORDS    .
FEATURES             Location/Qualifiers
     source          1..120
                     /organism=\"Test virus\"
                     /mol_type=\"genomic RNA\"
     5'UTR           1..10
     gene            11..100
                     /gene=\"ORF1\"
     CDS             join(11..40,
                     42..100)
                     /gene=\"ORF1\"
                     /codon_start=1
                     /product=\"polyprotein with a
                     long name\"
                     /note=\"a \"\"quoted\"\" word\"
                     /translation=\"MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEAR
                     QHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAAL\"
     mat_peptide     complement(<101..>110)
                     /product=\"leader\"
                     /pseudo
     misc_feature    J00194.1:1..10
ORIGIN
        1 attaaaggtt tataccttcc caggtaacaa accaaccaac tttcgatctc ttgtagatct
       61 gttctctaaa cgaactttaa aatctgtgtg gctgtcactc ggctgcatgc ttagtgcact
//
LOCUS       SECOND
ORIGIN
        1 acgtn
//
";

    #[test]
    fn parses_records() {
        let recs = parse(SAMPLE).unwrap();
        assert_eq!(recs.len(), 2);
        let r = &recs[0];
        assert_eq!((r.locus.as_str(), r.length, r.version.as_str()), ("TEST01", 120, "TEST01.2"));
        assert_eq!(r.definition, "Test virus isolate X, complete genome.");
        assert_eq!(r.sequence.len(), 120);
        assert!(r.sequence.starts_with(b"ATTAAAGGTT"));
        let keys: Vec<&str> = r.features.iter().map(|f| f.key.as_str()).collect();
        assert_eq!(keys, ["source", "5'UTR", "gene", "CDS", "mat_peptide", "misc_feature"]);
        let cds = &r.features[3];
        assert_eq!(cds.location, "join(11..40,42..100)");
        assert_eq!(cds.qualifier("product"), Some("polyprotein with a long name"));
        assert_eq!(cds.qualifier("codon_start"), Some("1"));
        assert_eq!(cds.qualifier("note"), Some("a \"quoted\" word"));
        assert!(cds.qualifier("translation").unwrap().starts_with("MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDG"));
        assert_eq!(r.features[4].qualifiers[1], ("pseudo".to_string(), String::new()));
        // Names: gene for genes, product for CDS/peptides, organism for the source.
        let names: Vec<String> = r.features.iter().map(Feature::name).collect();
        assert_eq!(names, ["Test virus", "5'UTR", "ORF1", "polyprotein with a long name", "leader", "misc_feature"]);
        // Without VERSION/ACCESSION the locus name; without a length the sequence length.
        assert_eq!((recs[1].version.as_str(), recs[1].length, recs[1].sequence.as_slice()), ("SECOND", 5, &b"ACGTN"[..]));
        assert!(parse("nothing here").is_err());
    }

    #[test]
    fn header_only() {
        // An interrupted download: the header still gives name and length.
        let r = &parse("LOCUS       NC_045512   29903 bp ss-RNA linear VRL\nDEFINITION  SARS-CoV-2.\nVERSION     NC_045512.2\n").unwrap()[0];
        assert_eq!((r.length, r.version.as_str(), r.features.len()), (29903, "NC_045512.2", 0));
    }

    #[test]
    fn locations() {
        let span = |start, end, complement| Span { start, end, complement };
        let l = parse_location("join(11..40,42..100)").unwrap();
        assert_eq!((l.spans.clone(), l.strand()), (vec![span(11, 40, false), span(42, 100, false)], 1));
        let l = parse_location("complement(join(5..9, 1..3))").unwrap();
        assert_eq!((l.spans.clone(), l.strand()), (vec![span(1, 3, true), span(5, 9, true)], -1));
        let l = parse_location("join(complement(20..30),complement(1..10))").unwrap();
        assert_eq!(l.strand(), -1);
        assert_eq!(parse_location("<1..>200").unwrap().spans, [span(1, 200, false)]);
        assert_eq!(parse_location("467").unwrap().spans, [span(467, 467, false)]);
        assert_eq!(parse_location("12^13").unwrap().spans, [span(12, 12, false)]);
        assert_eq!(parse_location("(102.110)..(200.210)").unwrap().spans, [span(102, 210, false)]);
        assert_eq!(parse_location("102.110").unwrap().spans, [span(102, 110, false)]);
        let l = parse_location("join(1..5,J00194.1:100..202,complement(8..9))").unwrap();
        assert_eq!((l.spans.len(), l.remote, l.strand()), (2, 1, 0));
        let l = parse_location("J00194.1:100..202").unwrap();
        assert_eq!((l.spans.len(), l.remote, l.strand()), (0, 1, 0));
        for bad in ["", "join(1..5", "5..1", "abc", "1..5)"] {
            assert!(parse_location(bad).is_err(), "{bad}");
        }
    }
}
