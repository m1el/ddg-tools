//! A progress counter and a cancel flag for long builds and evaluations (SPEC.md §18, WP12/WP13
//! items; usable from Python and Blender).
//!
//! Put a [`Progress`] in [`crate::BuildOptions::progress`] or [`crate::EvalOptions::progress`] and
//! keep a clone. Another thread polls [`Progress::stage`], [`Progress::done`] and
//! [`Progress::total`], and may call [`Progress::cancel`]; the running call then returns
//! [`PatError::Cancelled`] at the next check (between chunks of a few hundred neighborhoods or
//! queries; the work already done is discarded). Counting is by atomics, so it changes no result:
//! a run that is not cancelled returns the same bits with or without a [`Progress`].

use crate::PatError;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};

#[derive(Debug, Default)]
struct Inner {
    cancel: AtomicBool,
    done: AtomicU64,
    total: AtomicU64,
    stage: AtomicUsize,
}

/// The stages a build or evaluation reports, in the order of [`Progress::stage`].
pub const STAGES: [&str; 6] = ["idle", "neighborhoods", "network", "fit", "tree", "evaluate"];

/// A shared progress counter and cancel flag (cheap to clone; clones share the state).
///
/// Two options that hold clones of the same `Progress` compare equal.
#[derive(Debug, Clone, Default)]
pub struct Progress(Arc<Inner>);

impl PartialEq for Progress {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl Progress {
    /// A new counter, not cancelled, at stage `"idle"`.
    pub fn new() -> Progress {
        Progress::default()
    }

    /// Ask the running call to stop. Sticky until [`Progress::reset`].
    pub fn cancel(&self) {
        self.0.cancel.store(true, Ordering::Relaxed);
    }

    /// Whether [`Progress::cancel`] was called.
    pub fn is_cancelled(&self) -> bool {
        self.0.cancel.load(Ordering::Relaxed)
    }

    /// Clear the cancel flag and the counters (to reuse the handle for the next call).
    pub fn reset(&self) {
        self.0.cancel.store(false, Ordering::Relaxed);
        self.0.done.store(0, Ordering::Relaxed);
        self.0.total.store(0, Ordering::Relaxed);
        self.0.stage.store(0, Ordering::Relaxed);
    }

    /// Items finished in the current stage.
    pub fn done(&self) -> u64 {
        self.0.done.load(Ordering::Relaxed)
    }

    /// Items in the current stage (0 before the first stage starts).
    pub fn total(&self) -> u64 {
        self.0.total.load(Ordering::Relaxed)
    }

    /// The current stage: one of [`STAGES`] (`"neighborhoods"`, `"network"`, `"fit"`, `"tree"`,
    /// `"evaluate"`, or `"idle"`).
    pub fn stage(&self) -> &'static str {
        STAGES[self.0.stage.load(Ordering::Relaxed).min(STAGES.len() - 1)]
    }

    /// `done / total` of the current stage in `[0, 1]` (0 when the stage has no items yet).
    pub fn fraction(&self) -> f64 {
        let t = self.total();
        if t == 0 { 0.0 } else { (self.done() as f64 / t as f64).min(1.0) }
    }

    /// Start `stage` (an index into [`STAGES`]) with `total` items.
    pub(crate) fn begin(&self, stage: usize, total: u64) {
        self.0.done.store(0, Ordering::Relaxed);
        self.0.total.store(total, Ordering::Relaxed);
        self.0.stage.store(stage, Ordering::Relaxed);
    }

    /// `n` more items are finished.
    pub(crate) fn advance(&self, n: u64) {
        self.0.done.fetch_add(n, Ordering::Relaxed);
    }

    /// `Err(Cancelled)` if the flag is set.
    pub(crate) fn check(&self) -> Result<(), PatError> {
        if self.is_cancelled() { Err(PatError::Cancelled) } else { Ok(()) }
    }
}

/// `p.check()` for an optional handle.
pub(crate) fn check(p: Option<&Progress>) -> Result<(), PatError> {
    p.map_or(Ok(()), Progress::check)
}

pub(crate) const NEIGHBORHOODS: usize = 1;
pub(crate) const NETWORK: usize = 2;
pub(crate) const FIT: usize = 3;
pub(crate) const TREE: usize = 4;
pub(crate) const EVALUATE: usize = 5;
