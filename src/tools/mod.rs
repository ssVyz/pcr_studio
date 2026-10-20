//! External command-line tools (e.g. primersearch): where their binaries are,
//! whether they work, and running them. A tool that is not available only
//! disables the features that need it.
//!
//! Each tool's location is set in the settings: bundled with PCR Studio (the
//! default), found on `PATH`, or a custom path.

pub mod primersearch;
pub mod run;

use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// A tool PCR Studio knows how to drive.
#[derive(Debug)]
pub struct ToolSpec {
    /// Key in the settings.
    pub id: &'static str,
    pub name: &'static str,
    /// Executable name without the platform suffix (`.exe`).
    pub exe: &'static str,
    /// Oldest version with everything PCR Studio uses.
    pub min_version: (u32, u32, u32),
    /// What PCR Studio uses the tool for (shown in the settings).
    pub purpose: &'static str,
}

pub const PRIMERSEARCH: ToolSpec = ToolSpec {
    id: "primersearch",
    name: "primersearch",
    exe: "primersearch",
    // 0.1.1 added `--progress jsonl`.
    min_version: (0, 1, 1),
    purpose: "Finds primer sets covering the selected columns (Tools menu).",
};

/// All tools, in the order the settings list them.
pub const ALL: [&ToolSpec; 1] = [&PRIMERSEARCH];

/// Where a tool's binary comes from.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ToolSource {
    /// PCR Studio's own copy in `assets/tools`.
    #[default]
    Bundled,
    /// The first match on `PATH`.
    Path,
    /// A path set by the user.
    Custom,
}

impl ToolSource {
    pub const ALL: [ToolSource; 3] = [ToolSource::Bundled, ToolSource::Path, ToolSource::Custom];
}

impl std::fmt::Display for ToolSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            ToolSource::Bundled => "Bundled with PCR Studio",
            ToolSource::Path => "From PATH",
            ToolSource::Custom => "Custom location",
        })
    }
}

/// Location setting of one tool.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default)]
pub struct ToolConfig {
    pub source: ToolSource,
    /// Used with `ToolSource::Custom`.
    pub custom_path: String,
}

/// Result of checking a tool.
#[derive(Debug, Clone, PartialEq)]
pub enum ToolStatus {
    Checking,
    Ready { path: PathBuf, version: String },
    Unavailable(String),
}

impl ToolStatus {
    pub fn is_ready(&self) -> bool {
        matches!(self, ToolStatus::Ready { .. })
    }
}

fn file_name(spec: &ToolSpec) -> String {
    format!("{}{}", spec.exe, std::env::consts::EXE_SUFFIX)
}

/// Folders searched for bundled tools: `assets/tools` next to the executable,
/// then (for `cargo run` and copies of the executable) the one in the source
/// tree the program was built from.
pub fn bundled_dirs() -> Vec<PathBuf> {
    let mut dirs = Vec::new();
    if let Some(dir) = std::env::current_exe().ok().and_then(|p| p.parent().map(Path::to_path_buf)) {
        dirs.push(dir.join("assets").join("tools"));
    }
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("assets").join("tools");
    if !dirs.contains(&source) {
        dirs.push(source);
    }
    dirs
}

/// The binary the setting points to, or why there is none.
pub fn locate(spec: &ToolSpec, cfg: &ToolConfig) -> Result<PathBuf, String> {
    let name = file_name(spec);
    match cfg.source {
        ToolSource::Bundled => {
            let dirs = bundled_dirs();
            dirs.iter().map(|d| d.join(&name)).find(|p| p.is_file()).ok_or_else(|| {
                let places: Vec<String> = dirs.iter().map(|d| d.display().to_string()).collect();
                format!("{name} is not in {}", places.join(" or "))
            })
        }
        ToolSource::Path => std::env::var_os("PATH")
            .iter()
            .flat_map(std::env::split_paths)
            .map(|d| d.join(&name))
            .find(|p| p.is_file())
            .ok_or_else(|| format!("{name} was not found on PATH")),
        ToolSource::Custom => {
            let p = PathBuf::from(cfg.custom_path.trim());
            if cfg.custom_path.trim().is_empty() {
                Err("No path set".into())
            } else if p.is_file() {
                Ok(p)
            } else {
                Err(format!("{} does not exist", p.display()))
            }
        }
    }
}

/// Locates the tool and asks it for its version (`--version`). Blocks for up
/// to a few seconds; call it off the UI thread.
pub fn check(spec: &ToolSpec, cfg: &ToolConfig) -> ToolStatus {
    let path = match locate(spec, cfg) {
        Ok(p) => p,
        Err(e) => return ToolStatus::Unavailable(e),
    };
    match probe_version(&path) {
        Ok(text) => match parse_version(&text) {
            Some(v) if v >= spec.min_version => ToolStatus::Ready { path, version: format!("{}.{}.{}", v.0, v.1, v.2) },
            Some(v) => {
                let (a, b, c) = spec.min_version;
                ToolStatus::Unavailable(format!("version {}.{}.{} is too old; PCR Studio needs {a}.{b}.{c} or newer", v.0, v.1, v.2))
            }
            None => ToolStatus::Unavailable(format!("unexpected answer to --version: {}", text.trim())),
        },
        Err(e) => ToolStatus::Unavailable(e),
    }
}

const PROBE_TIMEOUT: Duration = Duration::from_secs(10);

/// Output of `<path> --version`.
fn probe_version(path: &Path) -> Result<String, String> {
    use std::io::Read;
    let mut cmd = run::command(path);
    cmd.arg("--version").stdin(std::process::Stdio::null()).stdout(std::process::Stdio::piped()).stderr(std::process::Stdio::null());
    let mut child = cmd.spawn().map_err(|e| format!("cannot start {}: {e}", path.display()))?;
    let started = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(s)) => break s,
            Ok(None) if started.elapsed() > PROBE_TIMEOUT => {
                let _ = child.kill();
                let _ = child.wait();
                return Err("no answer to --version".into());
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(20)),
            Err(e) => return Err(e.to_string()),
        }
    };
    let mut out = String::new();
    if let Some(mut o) = child.stdout.take() {
        let _ = o.read_to_string(&mut out);
    }
    if !status.success() {
        return Err(format!("--version failed ({status})"));
    }
    Ok(out)
}

/// `major.minor.patch` from the last word of the first line, e.g.
/// "primersearch 0.1.1".
fn parse_version(text: &str) -> Option<(u32, u32, u32)> {
    let word = text.lines().next()?.split_whitespace().last()?;
    let mut parts = word.trim_start_matches('v').split('.').map(|p| p.parse::<u32>().ok());
    let v = (parts.next()??, parts.next()??, parts.next().unwrap_or(Some(0))?);
    Some(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn versions() {
        assert_eq!(parse_version("primersearch 0.1.1\n"), Some((0, 1, 1)));
        assert_eq!(parse_version("tool v2.10\n"), Some((2, 10, 0)));
        assert_eq!(parse_version("primersearch\n"), None);
        assert_eq!(parse_version(""), None);
        assert!((0, 1, 1) >= PRIMERSEARCH.min_version && (0, 1, 0) < PRIMERSEARCH.min_version);
    }

    #[test]
    fn locating() {
        let custom = |p: &str| ToolConfig { source: ToolSource::Custom, custom_path: p.into() };
        assert!(locate(&PRIMERSEARCH, &custom("")).is_err());
        assert!(locate(&PRIMERSEARCH, &custom("Z:/no/such/tool.exe")).is_err());
        let exe = std::env::current_exe().unwrap();
        assert_eq!(locate(&PRIMERSEARCH, &custom(exe.to_str().unwrap())), Ok(exe));
        assert!(matches!(check(&PRIMERSEARCH, &custom("")), ToolStatus::Unavailable(_)));
    }
}
