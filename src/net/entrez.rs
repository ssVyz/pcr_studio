//! NCBI Entrez E-utilities (<https://www.ncbi.nlm.nih.gov/books/NBK25497/>).
//!
//! NCBI asks every program to send its name (`tool`) and a contact email
//! with each request, and to make at most 3 requests per second (10 with an
//! API key). Requests here always carry both and are spaced accordingly; the
//! email is required.

use super::{NetConfig, Response, get, read_text};
use std::io::{BufRead, BufReader, Read};
use std::sync::Mutex;
use std::time::{Duration, Instant};

const BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
const TOOL: &str = "PCR_Studio";

/// Who asks: the user's email (required by NCBI) and an optional API key.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Credentials {
    pub email: String,
    pub api_key: String,
}

/// A plausible email address (name@domain.tld, no spaces).
pub fn valid_email(s: &str) -> bool {
    let s = s.trim();
    let Some((user, domain)) = s.split_once('@') else { return false };
    !user.is_empty()
        && !s.contains(char::is_whitespace)
        && !domain.contains('@')
        && domain.contains('.')
        && !domain.starts_with('.')
        && !domain.ends_with('.')
}

/// An accession, with or without version: letters, digits, `_` and `.`,
/// with at least one digit (MN908947.3, NC_045512, NZ_CP012345.1).
pub fn valid_accession(s: &str) -> bool {
    let s = s.trim();
    (1..=40).contains(&s.len()) && s.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.') && s.bytes().any(|b| b.is_ascii_digit())
}

/// Waits until the next request is within NCBI's rate limit.
fn throttle(api_key: bool) {
    static NEXT: Mutex<Option<Instant>> = Mutex::new(None);
    let gap = Duration::from_millis(if api_key { 110 } else { 350 });
    let mut next = NEXT.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(t) = *next {
        let now = Instant::now();
        if t > now {
            std::thread::sleep(t - now);
        }
    }
    *next = Some(Instant::now() + gap);
}

/// One E-utility request (`efetch.fcgi`, …) with `tool`, `email` and the API key added.
fn call(net: &NetConfig, who: &Credentials, utility: &str, params: &[(&str, &str)]) -> Result<Response, String> {
    let email = who.email.trim();
    if !valid_email(email) {
        return Err("Enter your email address in Settings (NCBI) first: NCBI requires it with every request.".into());
    }
    let key = who.api_key.trim();
    throttle(!key.is_empty());
    let mut query = params.to_vec();
    query.extend([("tool", TOOL), ("email", email)]);
    if !key.is_empty() {
        query.push(("api_key", key));
    }
    get(net, &format!("{BASE}{utility}"), &query)
}

/// The gist of an error answer (plain text, JSON or XML) as a message.
fn error_message(status: u16, body: &mut dyn Read) -> String {
    if status == 429 {
        return "NCBI received too many requests; wait a moment and try again.".into();
    }
    let gist = gist(&read_text(body, 4096).unwrap_or_default());
    if gist.is_empty() { format!("NCBI answered with HTTP status {status}.") } else { format!("NCBI answered with HTTP status {status}: {gist}") }
}

/// Text without markup, on one line, at most 200 characters.
fn gist(s: &str) -> String {
    let mut out = String::new();
    let mut in_tag = false;
    for c in s.chars() {
        match c {
            '<' => in_tag = true,
            '>' if in_tag => {
                in_tag = false;
                out.push(' ');
            }
            _ if !in_tag => out.push(c),
            _ => {}
        }
    }
    let one_line = out.split_whitespace().collect::<Vec<_>>().join(" ");
    if one_line.chars().count() > 200 { format!("{}…", one_line.chars().take(199).collect::<String>()) } else { one_line }
}

/// Downloads the GenBank flat file of nucleotide record `accession`
/// (`efetch`, `db=nuccore`, `rettype=gb`). `accept_length` sees the length
/// on the LOCUS line first; if it returns false, the download stops after
/// the header (up to the feature table), which is returned with `false`.
pub fn fetch_genbank(net: &NetConfig, who: &Credentials, accession: &str, accept_length: impl Fn(usize) -> bool) -> Result<(String, bool), String> {
    let acc = accession.trim();
    if !valid_accession(acc) {
        return Err(format!("“{acc}” is not an accession number (like MN908947.3 or NC_045512)."));
    }
    let mut resp = call(net, who, "efetch.fcgi", &[("db", "nuccore"), ("id", acc), ("rettype", "gb"), ("retmode", "text")])?;
    if resp.status != 200 {
        return Err(error_message(resp.status, &mut resp.body));
    }
    let mut reader = BufReader::new(resp.body);
    let read_err = |e: std::io::Error| format!("Reading the record from NCBI failed: {e}");
    // The header, up to the feature table.
    let mut header = Vec::new();
    loop {
        let start = header.len();
        if reader.read_until(b'\n', &mut header).map_err(read_err)? == 0 {
            break;
        }
        let line = &header[start..];
        if line.starts_with(b"FEATURES") || line.starts_with(b"ORIGIN") || line.starts_with(b"//") || header.len() > 1 << 20 {
            break;
        }
    }
    let mut text = String::from_utf8_lossy(&header).into_owned();
    if !text.trim_start().starts_with("LOCUS") {
        let gist = gist(&text);
        return Err(if gist.is_empty() {
            format!("NCBI returned no nucleotide record for “{acc}”.")
        } else {
            format!("NCBI returned no GenBank record for “{acc}”: {gist}")
        });
    }
    let length = crate::genbank::parse(&text).ok().and_then(|r| r.first().map(|r| r.length)).unwrap_or(0);
    if !accept_length(length) {
        return Ok((text, false));
    }
    text.push_str(&read_text(&mut reader, u64::MAX)?);
    Ok((text, true))
}

/// Checks that NCBI can be reached with these settings (one small request);
/// returns a description of the connection.
pub fn test_connection(net: &NetConfig, who: &Credentials) -> Result<String, String> {
    let route = match super::proxy_in_use(net)? {
        Some(p) => format!("through the proxy {p}"),
        None => "directly (no proxy)".to_string(),
    };
    let mut resp = call(net, who, "einfo.fcgi", &[("retmode", "json")])?;
    if resp.status != 200 {
        return Err(error_message(resp.status, &mut resp.body));
    }
    let text = read_text(&mut resp.body, 1 << 20)?;
    if !text.contains("dblist") {
        return Err(format!("Unexpected answer from NCBI: {}", gist(&text)));
    }
    Ok(format!("Connected to NCBI E-utilities {route}."))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checks() {
        for ok in ["a@b.de", " name.surname@uni-hamburg.de "] {
            assert!(valid_email(ok), "{ok}");
        }
        for bad in ["", "a@b", "@b.de", "a b@c.de", "a@@b.de", "a@.de", "a@b.", "ab.de"] {
            assert!(!valid_email(bad), "{bad}");
        }
        for ok in ["MN908947.3", "NC_045512", "NZ_CP012345.1", "AB123456"] {
            assert!(valid_accession(ok), "{ok}");
        }
        for bad in ["", "ABC", "MN908947.3&rettype=fasta", "a b1", "x/1"] {
            assert!(!valid_accession(bad), "{bad}");
        }
        assert_eq!(gist("<html><body><h1>Error</h1>\n  <p>Bad   request</p></body></html>"), "Error Bad request");
        // Nothing is sent without an email address.
        let err = fetch_genbank(&NetConfig::default(), &Credentials::default(), "NC_045512.2", |_| true).unwrap_err();
        assert!(err.contains("email"), "{err}");
    }

    /// Fetches a record from NCBI (network). Run with the email NCBI should see:
    /// `PCR_STUDIO_NCBI_EMAIL=you@example.org cargo test -- --ignored entrez`
    #[test]
    #[ignore]
    fn live_fetch() {
        let email = std::env::var("PCR_STUDIO_NCBI_EMAIL").expect("set PCR_STUDIO_NCBI_EMAIL");
        let who = Credentials { email, api_key: String::new() };
        let net = NetConfig::default();
        println!("{}", test_connection(&net, &who).unwrap());
        let (text, complete) = fetch_genbank(&net, &who, "NC_045512.2", |len| len == 29903).unwrap();
        assert!(complete);
        let rec = &crate::genbank::parse(&text).unwrap()[0];
        assert_eq!((rec.version.as_str(), rec.length, rec.sequence.len()), ("NC_045512.2", 29903, 29903));
        assert!(rec.features.iter().any(|f| f.key == "CDS" && f.qualifier("gene") == Some("S")));
        // A length that does not match stops after the header.
        let (header, complete) = fetch_genbank(&net, &who, "NC_045512.2", |len| len == 100).unwrap();
        assert!(!complete && header.len() < 5000, "{}", header.len());
        assert!(fetch_genbank(&net, &who, "XX_0000000000", |_| true).is_err());
    }
}
