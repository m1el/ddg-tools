//! The crate's error type.
//!
//! Every error names the element that caused it (an index, a tensor name, an option). Index lists
//! are kept in full; their `Display` shows the first few.

use std::fmt;

/// Errors of every stage (validation, weights, prediction, fitting, evaluation).
#[derive(Debug, Clone, PartialEq)]
pub enum PatError {
    /// `positions` and `normals` have different lengths.
    LengthMismatch {
        /// Number of positions.
        positions: usize,
        /// Number of normals.
        normals: usize,
    },
    /// The cloud has fewer points than the stage needs.
    TooFewPoints {
        /// Number of points given.
        n: usize,
        /// Number needed.
        need: usize,
    },
    /// A coordinate is NaN or infinite.
    NonFinite {
        /// Index of the element (point, query, coefficient row).
        index: usize,
        /// What the element is (`"position"`, `"normal"`, `"query"`, …).
        what: &'static str,
    },
    /// A normal has length zero.
    ZeroNormal {
        /// Index of the point.
        index: usize,
    },
    /// A normal is not of unit length under [`crate::NormalPolicy::RequireUnit`].
    NonUnitNormal {
        /// Index of the point.
        index: usize,
        /// Its length.
        norm: f64,
    },
    /// The weights file is malformed (magic, missing tensor, shape, non-finite value).
    Weights {
        /// What is wrong, naming the tensor.
        reason: String,
    },
    /// The GPU was requested but is not available (feature, driver or device).
    Gpu {
        /// What is missing.
        reason: String,
    },
    /// The network produced non-finite outputs for these neighborhoods.
    NonFiniteOutput {
        /// Indices of the center points.
        indices: Vec<u32>,
    },
    /// Fitting produced non-finite tori for these points.
    NonFiniteTorus {
        /// Indices of the points.
        indices: Vec<u32>,
    },
    /// An option is out of its domain.
    InvalidOption {
        /// The option and why.
        what: String,
    },
    /// Two inputs disagree in size (coefficients vs. cloud, neighborhoods vs. cloud, buffers).
    SizeMismatch {
        /// Which inputs and their sizes.
        what: String,
    },
    /// The thread pool could not be created.
    ThreadPool {
        /// The reason given by rayon.
        reason: String,
    },
    /// A serialized file is malformed.
    Format {
        /// What is wrong.
        reason: String,
    },
    /// A query's neighbor set degenerates (PerQuery λ with `σ_x = 0`: every member at the query).
    DegenerateQuery {
        /// Index of the query.
        index: usize,
        /// Why.
        reason: String,
    },
    /// The network framework failed.
    Backend {
        /// What failed.
        reason: String,
    },
}

fn list(ix: &[u32]) -> String {
    const SHOW: usize = 10;
    let head: Vec<String> = ix.iter().take(SHOW).map(|i| i.to_string()).collect();
    if ix.len() > SHOW {
        format!("{} (and {} more)", head.join(", "), ix.len() - SHOW)
    } else {
        head.join(", ")
    }
}

impl fmt::Display for PatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PatError::LengthMismatch { positions, normals } => {
                write!(f, "{positions} positions but {normals} normals")
            }
            PatError::TooFewPoints { n, need } => {
                write!(f, "the cloud has {n} points, at least {need} are needed")
            }
            PatError::NonFinite { index, what } => write!(f, "{what} {index} is not finite"),
            PatError::ZeroNormal { index } => write!(f, "normal {index} has length zero"),
            PatError::NonUnitNormal { index, norm } => write!(
                f,
                "normal {index} has length {norm} (unit normals are required; use NormalPolicy::Normalize to opt in)"
            ),
            PatError::Weights { reason } => write!(f, "weights: {reason}"),
            PatError::Gpu { reason } => write!(f, "GPU: {reason}"),
            PatError::NonFiniteOutput { indices } => write!(
                f,
                "the network gave non-finite coefficients for {} neighborhoods: {}",
                indices.len(),
                list(indices)
            ),
            PatError::NonFiniteTorus { indices } => write!(
                f,
                "fitting gave non-finite tori at {} points: {}",
                indices.len(),
                list(indices)
            ),
            PatError::InvalidOption { what } => write!(f, "invalid option: {what}"),
            PatError::SizeMismatch { what } => write!(f, "size mismatch: {what}"),
            PatError::ThreadPool { reason } => write!(f, "thread pool: {reason}"),
            PatError::Format { reason } => write!(f, "format: {reason}"),
            PatError::DegenerateQuery { index, reason } => write!(f, "query {index}: {reason}"),
            PatError::Backend { reason } => write!(f, "network backend: {reason}"),
        }
    }
}

impl std::error::Error for PatError {}
