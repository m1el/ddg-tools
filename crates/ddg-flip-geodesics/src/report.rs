//! Options and report of the shortening (SPEC §4.2–§4.4, F13, F14), and what one step of it
//! did ([`Shortening`], [`ShortenStep`]).

use crate::error::FlipError;
use crate::network::{Joint, SegmentId};
use crate::wedge::Side;

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
    /// A final sweep found no joint to shorten, and every joint side below `π −
    /// straight_tolerance` is pinned or blocked by a path end: the final state certifies (SPEC
    /// §7 layer 1a) for any certificate tolerance `≥ straight_tolerance`.
    Geodesic,
    /// A final sweep found no joint it could shorten, but some joint side is still below
    /// `π − straight_tolerance`: one that stalled (the outer arc was not shorter, typically
    /// because the triangulation refused a flip, F6) and has not changed since, or one blocked
    /// by another part of the network. The final state does not certify at
    /// `straight_tolerance` (it may at the certificate's default `1e-4`, when the side is bent
    /// by less).
    Stalled,
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
    /// Out-and-back joints contracted (`deg W = 0`, SPEC §4.1 step 3, F12, F26).
    pub contractions: usize,
    /// Alg. 1 calls whose wedge had, at some point of the flipping, an edge incident on it
    /// twice (a self-edge at `b` with both ends inside the wedge; App. A; SPEC §7 layer 3).
    pub twice_incident_wedges: usize,
    /// Alg. 1 calls whose outer arc had, at some point of the flipping, an interior vertex of
    /// degree one (App. A, Thm A.5).
    pub degree_one_wedges: usize,
    /// Queue entries skipped because their segment is gone or their angle changed (F5).
    pub stale_entries: usize,
    /// Queue entries dropped because their wedge was not flexible (F5).
    pub blocked_entries: usize,
    /// Paths that contracted to nothing (F12), ascending. They keep their index and have no
    /// segment.
    pub contracted: Vec<usize>,
    /// Length before.
    pub initial_length: f64,
    /// Length after (the running sum, F14).
    pub final_length: f64,
}

/// What one shortening step did to a joint (one Alg. 1 call, or the App. B move of a
/// single-segment loop).
#[derive(Debug, Clone, PartialEq)]
pub enum Shortening {
    /// The joint was replaced by the outer arc of its wedge (Alg. 1, `deg W ≥ 1`): two
    /// segments of total length `old_length` by `new_segments` of total length
    /// `new_length < old_length`.
    Replaced {
        /// The new segments, in path order.
        new_segments: Vec<SegmentId>,
        /// `ℓ_ab + ℓ_bc`.
        old_length: f64,
        /// The length of the outer arc, summed in path order.
        new_length: f64,
    },
    /// An out-and-back joint (`deg W = 0`) was contracted: both segments removed (SPEC §4.1
    /// step 3, F12, F26).
    Contracted {
        /// `ℓ_ab + ℓ_bc`, the length removed.
        old_length: f64,
        /// Whether the path has no segment left.
        path_contracted: bool,
    },
    /// A single-segment loop was moved across the face on the side of its wedge (SPEC §4.3,
    /// App. B). The length grows.
    LoopMoved {
        /// The new segments, in path order.
        new_segments: Vec<SegmentId>,
        /// The length of the loop before.
        old_length: f64,
        /// The length after.
        new_length: f64,
    },
    /// The outer arc was not shorter (refused flips, F6): the path was left as it was; the
    /// flips made are kept.
    Stalled {
        /// `ℓ_ab + ℓ_bc`.
        old_length: f64,
        /// The length of the outer arc.
        new_length: f64,
    },
}

impl Shortening {
    /// The change of the network's length (negative when it shortened).
    pub fn length_change(&self) -> f64 {
        match self {
            Shortening::Replaced {
                old_length,
                new_length,
                ..
            }
            | Shortening::LoopMoved {
                old_length,
                new_length,
                ..
            } => new_length - old_length,
            Shortening::Contracted { old_length, .. } => -old_length,
            Shortening::Stalled { .. } => 0.0,
        }
    }
}

/// One step of the shortening, as passed to the observer of
/// [`GeodesicNetwork::shorten_with`](crate::GeodesicNetwork::shorten_with).
#[derive(Debug, Clone, PartialEq)]
pub struct ShortenStep {
    /// The iteration (from `1`).
    pub iteration: usize,
    /// The joint shortened.
    pub joint: Joint,
    /// The side of its wedge.
    pub side: Side,
    /// The joint's vertex.
    pub vertex: usize,
    /// The wedge angle on that side before the step.
    pub angle: f64,
    /// Flips made by the step.
    pub flips: usize,
    /// Flips the triangulation refused in the step.
    pub refused_flips: usize,
    /// What happened.
    pub shortening: Shortening,
    /// The network's length after the step (the running sum, F14).
    pub length: f64,
}
