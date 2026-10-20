//! Runs a tool as a child process: no console window, stdout collected,
//! JSON-lines events from stderr (`--progress jsonl` convention), cancellable.

use serde_json::Value;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Duration;

/// A command for `program` that opens no console window on Windows.
pub fn command(program: &Path) -> Command {
    #[allow(unused_mut)]
    let mut cmd = Command::new(program);
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }
    cmd
}

/// A file in the temp folder for a tool's input, deleted when dropped.
pub struct TempFile {
    pub path: PathBuf,
}

impl TempFile {
    pub fn new(extension: &str) -> TempFile {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, Ordering::Relaxed);
        TempFile { path: std::env::temp_dir().join(format!("pcr_studio_{}_{n}.{extension}", std::process::id())) }
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

/// The running process, so it can be killed from another thread (e.g. when
/// PCR Studio exits during a run).
#[derive(Default)]
pub struct Handle {
    child: Mutex<Option<Child>>,
}

impl Handle {
    pub fn kill(&self) {
        if let Ok(mut c) = self.child.lock()
            && let Some(child) = c.as_mut()
        {
            let _ = child.kill();
        }
    }
}

/// What a finished (or killed) run produced.
#[derive(Debug, Default)]
pub struct Output {
    pub stdout: String,
    /// Plain-text stderr lines (everything that is not a JSON event).
    pub stderr: String,
    /// Message of a `{"type":"error"}` event.
    pub error: Option<String>,
    /// Exit code (`None` if there is none, e.g. a killed process elsewhere than on Windows).
    pub code: Option<i32>,
    pub cancelled: bool,
}

impl Output {
    /// The best description of a failure: the error event, else stderr.
    pub fn failure(&self) -> String {
        match (&self.error, self.stderr.trim(), self.code) {
            (Some(e), _, _) => e.clone(),
            (None, s, _) if !s.is_empty() => s.to_string(),
            (None, _, Some(c)) => format!("exit code {c}"),
            (None, _, None) => "the process was terminated".into(),
        }
    }
}

/// Runs `cmd` to completion. Every JSON object on stderr goes to `on_event`
/// (from a reader thread); setting `cancel` kills the process.
pub fn run(mut cmd: Command, handle: &Handle, cancel: &AtomicBool, on_event: &(dyn Fn(&Value) + Sync)) -> Result<Output, String> {
    cmd.stdin(Stdio::null()).stdout(Stdio::piped()).stderr(Stdio::piped());
    let program = cmd.get_program().to_string_lossy().into_owned();
    let mut child = cmd.spawn().map_err(|e| format!("Cannot start {program}: {e}"))?;
    let stdout = child.stdout.take();
    let stderr = child.stderr.take();
    *handle.child.lock().map_err(|e| e.to_string())? = Some(child);

    let result = std::thread::scope(|s| {
        // Both pipes are drained continuously so a full pipe never blocks the tool.
        let out = s.spawn(move || {
            let mut buf = Vec::new();
            if let Some(mut o) = stdout {
                let _ = o.read_to_end(&mut buf);
            }
            String::from_utf8_lossy(&buf).into_owned()
        });
        let err = s.spawn(move || {
            let (mut plain, mut error) = (Vec::new(), None);
            if let Some(e) = stderr {
                for line in BufReader::new(e).split(b'\n').map_while(Result::ok) {
                    let line = String::from_utf8_lossy(&line).trim().to_string();
                    match serde_json::from_str::<Value>(&line).ok().filter(Value::is_object) {
                        Some(ev) if ev["type"] == "error" => error = ev["message"].as_str().map(|m| m.trim().to_string()),
                        Some(ev) => on_event(&ev),
                        None if !line.is_empty() => plain.push(line),
                        None => {}
                    }
                }
            }
            (plain, error)
        });
        let mut cancelled = false;
        let status = loop {
            if cancel.load(Ordering::Relaxed) && !cancelled {
                cancelled = true;
                handle.kill();
            }
            let polled = match handle.child.lock() {
                Ok(mut c) => c.as_mut().map(|ch| ch.try_wait()),
                Err(e) => break Err(e.to_string()),
            };
            match polled {
                Some(Ok(Some(st))) => break Ok(st),
                Some(Ok(None)) => std::thread::sleep(Duration::from_millis(50)),
                Some(Err(e)) => break Err(e.to_string()),
                None => break Err("the process handle was lost".into()),
            }
        };
        let stdout = out.join().unwrap_or_default();
        let (plain, error) = err.join().unwrap_or_default();
        status.map(|st| {
            // The error event repeats the plain "error: …" line.
            let plain: Vec<String> = plain.into_iter().filter(|l| error.as_ref().is_none_or(|e| *l != format!("error: {e}"))).collect();
            Output { stdout, stderr: plain.join("\n"), error, code: st.code(), cancelled }
        })
    });
    if let Ok(mut c) = handle.child.lock() {
        *c = None;
    }
    result
}
