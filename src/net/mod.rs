//! Network access. Everything PCR Studio sends over the internet goes
//! through this module; nothing outside `net` opens a connection.
//!
//! - [`entrez`]: NCBI E-utilities (GenBank records by accession).
//!
//! Requests block, so callers run them on background threads. HTTPS
//! certificates are checked against the operating system's trust store, so
//! networks that inspect HTTPS with their own (installed) certificate work.
//! The proxy comes from the settings, else from the system: the
//! `ALL_PROXY` / `HTTPS_PROXY` / `HTTP_PROXY` variables or, on Windows, the
//! proxy server of the Internet settings (automatic configuration scripts
//! are not read).

pub mod entrez;

use std::io::Read;
use std::time::Duration;

/// Connection settings.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct NetConfig {
    /// Proxy ("host:port" or "http://host:port"); empty: the system's.
    pub proxy: String,
}

const USER_AGENT: &str = concat!("PCR-Studio/", env!("CARGO_PKG_VERSION"));

/// Largest response body read (bytes).
const MAX_BODY: u64 = 256 * 1024 * 1024;

/// An HTTP response: its status and a reader of its body (size-limited).
pub struct Response {
    pub status: u16,
    pub body: Box<dyn Read + Send>,
}

/// The proxy URL of the settings (a missing scheme means `http://`).
fn proxy_url(p: &str) -> String {
    if p.contains("://") { p.to_string() } else { format!("http://{p}") }
}

fn proxy(cfg: &NetConfig) -> Result<Option<ureq::Proxy>, String> {
    match cfg.proxy.trim() {
        "" => Ok(ureq::Proxy::try_from_env()),
        p => ureq::Proxy::new(&proxy_url(p)).map(Some).map_err(|e| format!("The proxy “{p}” in Settings is not valid ({e}); example: proxy.example.org:8080")),
    }
}

/// The proxy requests go through ("host:port"), `None` for a direct connection.
pub fn proxy_in_use(cfg: &NetConfig) -> Result<Option<String>, String> {
    Ok(proxy(cfg)?.map(|p| format!("{}:{}", p.host(), p.port())))
}

fn agent(cfg: &NetConfig) -> Result<ureq::Agent, String> {
    let tls = ureq::tls::TlsConfig::builder().root_certs(ureq::tls::RootCerts::PlatformVerifier).build();
    Ok(ureq::Agent::config_builder()
        .proxy(proxy(cfg)?)
        .tls_config(tls)
        .timeout_connect(Some(Duration::from_secs(20)))
        .timeout_global(Some(Duration::from_secs(180)))
        .http_status_as_error(false)
        .user_agent(USER_AGENT)
        .build()
        .new_agent())
}

/// Host name of an URL, for messages.
fn host_of(url: &str) -> &str {
    let rest = url.split_once("://").map_or(url, |(_, r)| r);
    rest.split(['/', '?', ':']).next().unwrap_or(rest)
}

/// GET `url` with query parameters (they are percent-encoded here).
pub fn get(cfg: &NetConfig, url: &str, query: &[(&str, &str)]) -> Result<Response, String> {
    let host = host_of(url);
    let mut req = agent(cfg)?.get(url);
    for (k, v) in query {
        req = req.query(*k, *v);
    }
    let resp = req.call().map_err(|e| describe(&e, host, cfg))?;
    let status = resp.status().as_u16();
    let body = resp.into_body().into_with_config().limit(MAX_BODY).reader();
    Ok(Response { status, body: Box::new(body) })
}

/// A connection error as a message with a hint at what to check.
fn describe(e: &ureq::Error, host: &str, cfg: &NetConfig) -> String {
    let hint = if cfg.proxy.trim().is_empty() {
        "If your network needs a proxy, enter it in Settings (NCBI)."
    } else {
        "Check the proxy in Settings (NCBI)."
    };
    match e {
        ureq::Error::HostNotFound => format!("{host} was not found: no internet connection, or the network needs a proxy. {hint}"),
        ureq::Error::Timeout(_) => format!("The connection to {host} timed out. {hint}"),
        ureq::Error::ConnectionFailed | ureq::Error::Io(_) => format!("Could not connect to {host} ({e}). {hint}"),
        ureq::Error::Tls(_) | ureq::Error::Rustls(_) | ureq::Error::Pem(_) => format!(
            "The secure connection to {host} failed ({e}). If your network inspects HTTPS traffic, its certificate must be installed in Windows."
        ),
        ureq::Error::InvalidProxyUrl => "The proxy in Settings is not valid; example: proxy.example.org:8080".into(),
        _ => format!("The request to {host} failed: {e}"),
    }
}

/// Reads up to `max` bytes of a body as text (invalid UTF-8 replaced).
pub fn read_text(body: &mut dyn Read, max: u64) -> Result<String, String> {
    let mut bytes = Vec::new();
    body.take(max).read_to_end(&mut bytes).map_err(|e| format!("Reading the answer failed: {e}"))?;
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers() {
        assert_eq!(host_of("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?x=1"), "eutils.ncbi.nlm.nih.gov");
        assert_eq!(host_of("http://localhost:8080/"), "localhost");
        assert_eq!(proxy_url("proxy.example.org:8080"), "http://proxy.example.org:8080");
        assert_eq!(proxy_url("socks5://p:1080"), "socks5://p:1080");
        let cfg = NetConfig { proxy: "proxy.example.org:3128".into() };
        assert_eq!(proxy_in_use(&cfg).unwrap().as_deref(), Some("proxy.example.org:3128"));
        assert!(agent(&cfg).is_ok());
    }
}
