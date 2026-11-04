//! Explicit thread pools (SPEC.md §5.5): every parallel stage runs on its own pool, never the
//! global one.

use crate::PatError;

/// A rayon pool with `threads` workers (`None`: one per logical core). `Some(0)` is an error.
pub fn make_pool(threads: Option<usize>) -> Result<rayon::ThreadPool, PatError> {
    let n = match threads {
        Some(0) => return Err(PatError::InvalidOption { what: "threads = Some(0)".into() }),
        Some(t) => t,
        None => std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .thread_name(|i| format!("pat-{i}"))
        .build()
        .map_err(|e| PatError::ThreadPool { reason: e.to_string() })
}
