//! Minimal FASTQ reader working on bytes: four-line records, but also
//! sequence and quality wrapped over several lines, and empty reads (as some
//! trimmers write them).

use std::io::BufRead;

#[derive(Debug, Clone, Default)]
pub struct Record {
    /// Header line without the leading `@`.
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl Record {
    /// Read name: the header up to the first whitespace.
    pub fn id(&self) -> &[u8] {
        let end = self.head.iter().position(|b| b.is_ascii_whitespace()).unwrap_or(self.head.len());
        &self.head[..end]
    }

    /// Header after the read name.
    pub fn desc(&self) -> &[u8] {
        let id = self.id().len();
        self.head[id..].trim_ascii()
    }
}

pub struct Reader<R> {
    inner: R,
    line: Vec<u8>,
    line_no: u64,
}

impl<R: BufRead> Reader<R> {
    pub fn new(inner: R) -> Reader<R> {
        Reader { inner, line: Vec::new(), line_no: 0 }
    }

    /// Reads the next line without its line break; `false` at the end.
    fn read_line(&mut self) -> Result<bool, String> {
        self.line.clear();
        let n = self.inner.read_until(b'\n', &mut self.line).map_err(|e| e.to_string())?;
        if n == 0 {
            return Ok(false);
        }
        self.line_no += 1;
        while matches!(self.line.last(), Some(b'\n' | b'\r')) {
            self.line.pop();
        }
        Ok(true)
    }

    /// Reads the next record into `rec`; `false` at the end of the file.
    pub fn next(&mut self, rec: &mut Record) -> Result<bool, String> {
        // Blank lines between records are skipped.
        loop {
            if !self.read_line()? {
                return Ok(false);
            }
            if !self.line.is_empty() {
                break;
            }
        }
        if self.line[0] != b'@' {
            return Err(format!("line {}: expected a FASTQ record header starting with '@'", self.line_no));
        }
        let start = self.line_no;
        let incomplete = || format!("line {start}: the FASTQ record is incomplete");
        rec.head.clear();
        rec.head.extend_from_slice(&self.line[1..]);
        rec.seq.clear();
        rec.qual.clear();
        loop {
            if !self.read_line()? {
                return Err(incomplete());
            }
            if self.line.first() == Some(&b'+') {
                break;
            }
            rec.seq.extend_from_slice(&self.line);
        }
        // At least one quality line, also for an empty read.
        let mut first = true;
        while first || rec.qual.len() < rec.seq.len() {
            first = false;
            if !self.read_line()? {
                if rec.seq.is_empty() {
                    break;
                }
                return Err(incomplete());
            }
            rec.qual.extend_from_slice(&self.line);
        }
        if rec.qual.len() != rec.seq.len() {
            return Err(format!("line {start}: the read has {} bases but {} quality values", rec.seq.len(), rec.qual.len()));
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_all(s: &str) -> Result<Vec<Record>, String> {
        let mut r = Reader::new(s.as_bytes());
        let mut out = Vec::new();
        let mut rec = Record::default();
        while r.next(&mut rec)? {
            out.push(rec.clone());
        }
        Ok(out)
    }

    #[test]
    fn parses_records() {
        let recs = read_all("@r1 1:N:0:1\nACGT\n+\nIIII\n\n@r2\r\nAC\r\nGT\r\n+r2\r\nII\r\nI#\r\n@empty\n\n+\n\n@last\nA\n+\n5").unwrap();
        assert_eq!(recs.len(), 4);
        assert_eq!((recs[0].id(), recs[0].desc()), (&b"r1"[..], &b"1:N:0:1"[..]));
        assert_eq!((recs[1].seq.as_slice(), recs[1].qual.as_slice()), (&b"ACGT"[..], &b"III#"[..]));
        assert!(recs[2].seq.is_empty() && recs[2].qual.is_empty());
        assert_eq!(recs[3].qual, b"5");
        // An empty read at the very end without its empty quality line.
        assert_eq!(read_all("@e\n\n+").unwrap().len(), 1);
    }

    #[test]
    fn rejects_broken_records() {
        assert!(read_all(">fasta\nACGT\n").unwrap_err().contains("'@'"));
        assert!(read_all("@r\nACGT\n+\nII").unwrap_err().contains("incomplete"));
        assert!(read_all("@r\nACGT\n").unwrap_err().contains("incomplete"));
        assert!(read_all("@r\nAC\n+\nIII\n").unwrap_err().contains("2 bases but 3"));
    }
}
