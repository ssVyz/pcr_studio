//! Background jobs: work runs on a std thread, the UI polls progress and
//! receives the result as a message.

use iced::Task;
use iced::futures::channel::oneshot;
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::{Arc, Mutex};

/// A value that can travel inside a `Clone` message and be taken exactly once.
pub struct Payload<T>(Arc<Mutex<Option<T>>>);

impl<T> Clone for Payload<T> {
    fn clone(&self) -> Self {
        Payload(self.0.clone())
    }
}

impl<T> std::fmt::Debug for Payload<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Payload")
    }
}

impl<T> Payload<T> {
    pub fn new(v: T) -> Self {
        Payload(Arc::new(Mutex::new(Some(v))))
    }
    pub fn take(&self) -> Option<T> {
        self.0.lock().ok()?.take()
    }
}

#[derive(Default)]
pub struct Progress {
    fraction: AtomicU32,
    message: Mutex<String>,
    pub cancel: AtomicBool,
}

impl Progress {
    pub fn set(&self, fraction: f32, message: &str) {
        self.fraction.store(fraction.clamp(0.0, 1.0).to_bits(), Ordering::Relaxed);
        if let Ok(mut m) = self.message.lock()
            && m.as_str() != message {
                *m = message.to_string();
            }
    }
    pub fn set_fraction(&self, fraction: f32) {
        self.fraction.store(fraction.clamp(0.0, 1.0).to_bits(), Ordering::Relaxed);
    }
    pub fn fraction(&self) -> f32 {
        f32::from_bits(self.fraction.load(Ordering::Relaxed))
    }
    pub fn message(&self) -> String {
        self.message.lock().map(|m| m.clone()).unwrap_or_default()
    }
}

pub struct Job {
    pub title: String,
    pub progress: Arc<Progress>,
    pub cancellable: bool,
}

/// Runs `work` on a new thread; `done` turns its result into a message
/// (`None` if the work panicked).
pub fn spawn<T, M>(
    title: &str,
    cancellable: bool,
    work: impl FnOnce(&Progress) -> T + Send + 'static,
    done: impl FnOnce(Option<T>) -> M + Send + 'static,
) -> (Job, Task<M>)
where
    T: Send + 'static,
    M: Send + 'static,
{
    let progress = Arc::new(Progress::default());
    progress.set(0.0, title);
    let p = progress.clone();
    let (tx, rx) = oneshot::channel();
    std::thread::Builder::new()
        .name(title.to_string())
        .spawn(move || {
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(&p))).ok();
            let _ = tx.send(r);
        })
        .expect("failed to spawn worker thread");
    let task = Task::perform(rx, move |r| done(r.ok().flatten()));
    (Job { title: title.to_string(), progress, cancellable }, task)
}
