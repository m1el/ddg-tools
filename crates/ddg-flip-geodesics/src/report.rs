//! Options and report of the shortening (SPEC §4.2–§4.4, F13, F14). The shortening itself
//! (Alg. 1, Alg. 2) is work package F2; these types fix what it takes and returns.

use crate::error::FlipError;

/// Options of the shortening (Alg. 2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShortenOptions {
    /// A joint is straight when `α^min > π − straight_tolerance` (radians); default `1e-5`
    /// (`gc: fg.h:161`, F3).
    pub straight_tolerance: f64,
    /// Iteration limit; `None` is `100 (|E| + m₀)` with `m₀` the initial number of segments
    /// (F13).
    pub max_iterations: Option<usize>,
    /// Stop once the length is below this fraction of the initial length (Sec. 6.1, e.g.
    /// `0.9`; `gc: fg.cpp:967-973`).
    pub stop_length_ratio: Option<f64>,
    /// Stop once every queued angle is at least this (Sec. 4.2 end).
    pub stop_min_angle: Option<f64>,
    /// At most this many final sweeps that find work (F5); default `4`.
    pub max_sweeps: usize,
    /// Run [`Network::check`](crate::Network::check) after every iteration (tests; costly).
    pub check_invariants: bool,
}

impl Default for ShortenOptions {
    fn default() -> Self {
        ShortenOptions {
            straight_tolerance: 1e-5,
            max_iterations: None,
            stop_length_ratio: None,
            stop_min_angle: None,
            max_sweeps: 4,
            check_invariants: false,
        }
    }
}

impl ShortenOptions {
    /// Checks the ranges.
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`] naming the option.
    pub fn validate(&self) -> Result<(), FlipError> {
        let bad = |what: String| Err(FlipError::InvalidOptions { what });
        if !(self.straight_tolerance.is_finite() && self.straight_tolerance >= 0.0) {
            return bad(format!(
                "straight_tolerance must be finite and >= 0, got {}",
                self.straight_tolerance
            ));
        }
        if let Some(r) = self.stop_length_ratio
            && !(r.is_finite() && r > 0.0)
        {
            return bad(format!("stop_length_ratio must be finite and > 0, got {r}"));
        }
        if let Some(a) = self.stop_min_angle
            && !a.is_finite()
        {
            return bad(format!("stop_min_angle must be finite, got {a}"));
        }
        Ok(())
    }
}

/// Why the shortening stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShortenStatus {
    /// A final sweep found no joint to shorten.
    Geodesic,
    /// `max_iterations` was reached.
    IterationLimit,
    /// The length fell below `stop_length_ratio` times the initial length.
    LengthRatio,
    /// Every queued angle reached `stop_min_angle`.
    MinAngle,
    /// `max_sweeps` sweeps found work and the last one left joints queued.
    SweepLimit,
}

/// What the shortening did.
#[derive(Debug, Clone, PartialEq)]
pub struct ShortenReport {
    /// Why it stopped.
    pub status: ShortenStatus,
    /// Joints shortened (Alg. 1 calls).
    pub iterations: usize,
    /// Edge flips made.
    pub flips: usize,
    /// Flips skipped because the triangulation refused them although `β < π` (F6).
    pub refused_flips: usize,
    /// Replacements skipped by the length guard (F6).
    pub stalled: usize,
    /// Final sweeps run (F5).
    pub sweeps: usize,
    /// Single-segment loops moved across a face (App. B); the only steps that lengthen.
    pub loop_perturbations: usize,
    /// Wedges with a twice-incident edge met by Alg. 1 (App. A; SPEC §7 layer 3).
    pub twice_incident_wedges: usize,
    /// Paths that contracted to nothing and were removed (F12).
    pub contracted: Vec<usize>,
    /// Length before.
    pub initial_length: f64,
    /// Length after (the running sum, F14).
    pub final_length: f64,
}
