//! The intrinsic optimal Delaunay triangulation (iODT, paper Sec. 4.3).
//!
//! An optimal Delaunay triangulation (Chen and Xu 2004) improves a triangulation by moving its
//! vertices, not only by adding vertices: each vertex is moved towards the mean of the
//! circumcentres of its triangles, weighted by their areas (Chen and Holst 2011, Eq. 4.13),
//! and the triangulation is flipped back to Delaunay. On an intrinsic triangulation only the
//! inserted vertices can move (moving an input vertex would change the surface), so new
//! vertices are made by splitting the edges longer than a target length (Tournois et al.
//! 2009), as the paper does. The update at a vertex is computed in its own frame, as the mean
//! of the vectors to the circumcentres in the planar layout of each triangle (the paper:
//! "rather than averaging these locations, we simply average the vectors").
//!
//! One iteration ([`IntrinsicTriangulation::optimal_delaunay`]):
//!
//! 1. split, at its midpoint, every edge longer than the target length (longest first);
//! 2. flip to Delaunay;
//! 3. for every movable vertex in the order of the ids: compute the update with the current
//!    triangulation (Gauss–Seidel), take the step `α` of the options; if the vertex would leave
//!    its star, halve the step until it does not (Chen and Holst's line search `α = 1, ½, ¼,
//!    …`; at most [`OptimalDelaunayOptions::max_halvings`] times, then the vertex stays);
//! 4. flip to Delaunay.
//!
//! It stops when an iteration split nothing and moved no vertex by more than the motion
//! tolerance (relative to the mean length of the vertex's edges), or at the iteration limit,
//! or when the limit on inserted vertices is reached. Everything is deterministic. The
//! decisions are SPEC N21 to N23.

use crate::delaunay::{DelaunayOptions, DelaunayReport};
use crate::error::{IntrinsicError, MoveRefusal};
use crate::operators::ShapeStatistics;
use crate::refine::InsertionLimit;
use crate::triangulation::IntrinsicTriangulation;

/// Which edges [`IntrinsicTriangulation::optimal_delaunay`] splits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TargetLength {
    /// None: only the vertices that are already inserted are moved.
    NoSplits,
    /// Edges longer than this length are split.
    Absolute(f64),
    /// Edges longer than this multiple of the mean edge length of the input mesh are split.
    MeanInputEdge(f64),
}

/// What happens to inserted vertices on the boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryVertices {
    /// They stay where they are (Chen and Holst 2011 move interior vertices only).
    Fixed,
    /// They move along the boundary by the component of their update along it.
    AlongBoundary,
}

/// Which centres of the triangles the update averages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OdtCenters {
    /// The circumcentres of all triangles (paper Sec. 4.3).
    Circumcenters,
    /// The circumcentres, except for triangles with a vertex on the boundary, which contribute
    /// their barycentre (Chen and Holst 2011, Sec. 4: a stretched triangle at the boundary
    /// cannot be flipped away). The same as `Circumcenters` on a closed surface.
    BarycentersAtBoundary,
}

/// Options of [`IntrinsicTriangulation::optimal_delaunay`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimalDelaunayOptions {
    /// Which edges are split at the start of every iteration.
    pub target_length: TargetLength,
    /// The largest number of iterations.
    pub max_iterations: usize,
    /// The run has converged when an iteration split no edge and moved no vertex by more than
    /// this fraction of the mean length of its edges (measured before the move).
    pub motion_tolerance: f64,
    /// The step `α` of the update (Chen and Holst 2011, Eq. 4.13), in `(0, 1]`: `1` moves a
    /// vertex all the way to the weighted mean of the circumcentres.
    pub step: f64,
    /// How often the step of one vertex is halved when the move would leave its star (or make
    /// a triangle too thin) before the vertex is left where it is for this iteration.
    pub max_halvings: usize,
    /// Inserted vertices on the boundary.
    pub boundary: BoundaryVertices,
    /// The centres the update averages.
    pub centers: OdtCenters,
    /// The limit on the vertices the splits insert; reaching it ends the run (after that
    /// iteration's moves and flips) and is reported, not an error.
    pub max_insertions: InsertionLimit,
    /// The Delaunay flips.
    pub delaunay: DelaunayOptions,
}

impl Default for OptimalDelaunayOptions {
    /// Split edges longer than the mean edge length of the input, at most 50 iterations, a
    /// motion tolerance of `1e-3`, the full step with up to 8 halvings, boundary vertices
    /// fixed, barycentres for the triangles at the boundary, at most 10 insertions per vertex
    /// (1000 at least), the default Delaunay flips. SPEC N21 says why.
    fn default() -> Self {
        OptimalDelaunayOptions {
            target_length: TargetLength::MeanInputEdge(1.0),
            max_iterations: 50,
            motion_tolerance: 1e-3,
            step: 1.0,
            max_halvings: 8,
            boundary: BoundaryVertices::Fixed,
            centers: OdtCenters::BarycentersAtBoundary,
            max_insertions: InsertionLimit::PerVertex(10),
            delaunay: DelaunayOptions::default(),
        }
    }
}

/// How [`IntrinsicTriangulation::optimal_delaunay`] ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptimalDelaunayStatus {
    /// The last iteration split no edge and moved no vertex by more than the tolerance.
    Converged,
    /// The limit on iterations was reached.
    IterationLimit,
    /// The limit on inserted vertices was reached.
    InsertionLimit,
}

/// What one iteration did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OdtIteration {
    /// Edges split.
    pub n_split: usize,
    /// Splits refused (a new triangle would have been too thin, N10), not counted in `n_split`.
    pub n_split_refused: usize,
    /// Vertices moved.
    pub n_moved: usize,
    /// Of those, the vertices moved with a halved step.
    pub n_halved: usize,
    /// Vertices that could not be moved even with the smallest step.
    pub n_not_moved: usize,
    /// Vertices skipped because an edge joins them to themselves.
    pub n_skipped: usize,
    /// Flips (after the splits and after the moves).
    pub n_flips: usize,
    /// The longest move.
    pub max_motion: f64,
    /// The longest move relative to the mean length of the vertex's edges before it.
    pub max_relative_motion: f64,
}

/// What [`IntrinsicTriangulation::optimal_delaunay`] did and what it left.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimalDelaunayReport {
    /// How it ended.
    pub status: OptimalDelaunayStatus,
    /// Number of iterations.
    pub n_iterations: usize,
    /// Vertices added (edge splits).
    pub n_inserted: usize,
    /// The target length that was used (`None` for [`TargetLength::NoSplits`]).
    pub target_length: Option<f64>,
    /// Every iteration, in order.
    pub iterations: Vec<OdtIteration>,
    /// The largest motion of the last iteration.
    pub last_max_motion: f64,
    /// The largest relative motion of the last iteration.
    pub last_max_relative_motion: f64,
    /// Flips in all, the first flip to Delaunay included.
    pub n_flips: usize,
    /// The report of the last flip to Delaunay.
    pub delaunay: DelaunayReport,
    /// Angles and areas of the result.
    pub statistics: ShapeStatistics,
}

impl IntrinsicTriangulation {
    /// The intrinsic optimal Delaunay triangulation (paper Sec. 4.3; see the
    /// [module documentation](self) for one iteration and SPEC N21 to N23 for the rules):
    /// splits edges longer than the target length, moves the inserted vertices towards the
    /// area-weighted mean of the circumcentres of their triangles, flips to Delaunay, and
    /// repeats. Input vertices never move. The result is a Delaunay triangulation of the same
    /// surface (up to the flips the Delaunay test refuses, reported in
    /// [`OptimalDelaunayReport::delaunay`]).
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidOptions`] for a target length that is not positive and finite,
    /// a step outside `(0, 1]`, a negative or non-finite motion tolerance, or invalid Delaunay
    /// options; [`IntrinsicError::MoveNotLocated`] and the errors of the flips and splits if
    /// the triangulation turns out to be inconsistent (nothing is repaired). Reaching a limit
    /// is not an error.
    pub fn optimal_delaunay(
        &mut self,
        options: &OptimalDelaunayOptions,
    ) -> Result<OptimalDelaunayReport, IntrinsicError> {
        let invalid = |reason: String| Err(IntrinsicError::InvalidOptions { reason });
        let target = match options.target_length {
            TargetLength::NoSplits => None,
            TargetLength::Absolute(h) => Some(h),
            TargetLength::MeanInputEdge(k) => {
                let lengths = self.input.geometry().edge_lengths();
                Some(k * lengths.iter().sum::<f64>() / lengths.len() as f64)
            }
        };
        if let Some(h) = target
            && !(h.is_finite() && h > 0.0)
        {
            return invalid(format!(
                "the target length of the optimal Delaunay triangulation must be positive and \
                 finite, got {h}"
            ));
        }
        if !(options.step > 0.0 && options.step <= 1.0) {
            return invalid(format!(
                "the step of the optimal Delaunay update must be in (0, 1], got {}",
                options.step
            ));
        }
        if !(options.motion_tolerance.is_finite() && options.motion_tolerance >= 0.0) {
            return invalid(format!(
                "the motion tolerance must be finite and not negative, got {}",
                options.motion_tolerance
            ));
        }
        let first = self.flip_to_delaunay(&options.delaunay)?;
        let mut n_flips = first.n_flips;
        let mut last_delaunay = first;
        let budget = match options.max_insertions {
            InsertionLimit::PerVertex(k) => k.saturating_mul(self.n_vertices()).max(1000),
            InsertionLimit::Vertices(n) => n,
            InsertionLimit::Unlimited => usize::MAX,
        };
        let barycenters = options.centers == OdtCenters::BarycentersAtBoundary;
        let mut n_inserted = 0;
        let mut iterations = Vec::new();
        let mut status = OptimalDelaunayStatus::IterationLimit;
        for _ in 0..options.max_iterations {
            let mut it = OdtIteration {
                n_split: 0,
                n_split_refused: 0,
                n_moved: 0,
                n_halved: 0,
                n_not_moved: 0,
                n_skipped: 0,
                n_flips: 0,
                max_motion: 0.0,
                max_relative_motion: 0.0,
            };
            let mut limit_reached = false;
            // 1. Splits, longest edge first.
            if let Some(h) = target {
                let mut long: Vec<(f64, u32)> = (0..self.n_edges())
                    .filter(|&e| self.edge_lengths[e] > h)
                    .map(|e| (self.edge_lengths[e], e as u32))
                    .collect();
                long.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
                for (_, e) in long {
                    if n_inserted >= budget {
                        limit_reached = true;
                        break;
                    }
                    match self.split_edge(e as usize, 0.5) {
                        Ok(_) => {
                            it.n_split += 1;
                            n_inserted += 1;
                        }
                        Err(IntrinsicError::InsertionRefused { .. }) => it.n_split_refused += 1,
                        Err(other) => return Err(other),
                    }
                }
                if it.n_split > 0 {
                    let report = self.flip_to_delaunay(&options.delaunay)?;
                    it.n_flips += report.n_flips;
                }
            }
            // 2. Moves, Gauss–Seidel in the order of the ids.
            for v in self.n_input_vertices()..self.n_vertices() {
                let boundary = self.is_boundary_vertex(v);
                if boundary && options.boundary == BoundaryVertices::Fixed {
                    continue;
                }
                let goal = match self.odt_target(v, barycenters) {
                    Ok(goal) => goal,
                    Err(IntrinsicError::MoveRefused {
                        reason: MoveRefusal::SelfEdge { .. },
                        ..
                    }) => {
                        it.n_skipped += 1;
                        continue;
                    }
                    Err(other) => return Err(other),
                };
                let goal = if boundary { [goal[0], 0.0] } else { goal };
                let length = goal[0].hypot(goal[1]) * options.step;
                if length == 0.0 {
                    continue;
                }
                let (mut sum, mut count) = (0.0, 0usize);
                for h in self.outgoing_halfedges(v) {
                    sum += self.edge_lengths[self.halfedge_edge[h] as usize];
                    count += 1;
                }
                let mean_edge = sum / count as f64;
                let mut alpha = options.step;
                let mut moved = None;
                for halving in 0..=options.max_halvings {
                    match self.move_vertex(v, [alpha * goal[0], alpha * goal[1]]) {
                        Ok(m) => {
                            moved = Some((m.distance, halving > 0));
                            break;
                        }
                        Err(IntrinsicError::MoveRefused {
                            reason:
                                MoveRefusal::NotInStar { .. }
                                | MoveRefusal::TooThin { .. }
                                | MoveRefusal::OnInputVertex { .. },
                            ..
                        }) => alpha *= 0.5,
                        Err(other) => return Err(other),
                    }
                }
                match moved {
                    Some((distance, halved)) => {
                        it.n_moved += 1;
                        it.n_halved += usize::from(halved);
                        it.max_motion = it.max_motion.max(distance);
                        it.max_relative_motion = it.max_relative_motion.max(distance / mean_edge);
                    }
                    None => it.n_not_moved += 1,
                }
            }
            // 3. Back to Delaunay.
            let report = self.flip_to_delaunay(&options.delaunay)?;
            it.n_flips += report.n_flips;
            last_delaunay = report;
            n_flips += it.n_flips;
            iterations.push(it);
            if limit_reached {
                status = OptimalDelaunayStatus::InsertionLimit;
                break;
            }
            if it.n_split == 0 && it.max_relative_motion <= options.motion_tolerance {
                status = OptimalDelaunayStatus::Converged;
                break;
            }
        }
        let last = iterations.last().copied();
        Ok(OptimalDelaunayReport {
            status,
            n_iterations: iterations.len(),
            n_inserted,
            target_length: target,
            last_max_motion: last.map_or(0.0, |it| it.max_motion),
            last_max_relative_motion: last.map_or(0.0, |it| it.max_relative_motion),
            iterations,
            n_flips,
            delaunay: last_delaunay,
            statistics: self.shape_statistics(),
        })
    }
}
