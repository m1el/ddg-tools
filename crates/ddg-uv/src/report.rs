//! What a run of the pipeline says about itself.

use ddg_optimal_cones::Status;

use crate::measure::AreaDistortion;
use crate::pack::PackReport;
use crate::feedback::ChartFeedback;
use crate::spikes::SpikeReport;

/// What was done about triangles without area.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct DegenerateReport {
    /// Number of triangles of the given triangulation whose edge lengths do not form a
    /// triangle of positive area.
    pub n_degenerate_triangles: usize,
    /// Whether the edge lengths were changed (`epsilon > 0`).
    pub mollified: bool,
    /// The amount that was added to the length of **every** edge of the mesh.
    pub epsilon: f64,
    /// The absolute margin of the mollification: its relative margin times the mean edge
    /// length. Zero if mollification did not run.
    pub delta: f64,
    /// Number of triangles that missed the margin (at least `n_degenerate_triangles` when
    /// mollification ran).
    pub n_triangles_below_margin: usize,
    /// The polygon of the triangle that determined `epsilon`.
    pub worst_polygon: Option<u32>,
}

/// What happened to the seams.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SeamReport {
    /// Number of distinct seam edges in the input.
    pub n_given: usize,
    /// Number of them that were cut: interior polygon edges.
    pub n_cut: usize,
    /// Number of them that lie on the boundary of the mesh already (nothing to cut).
    pub n_on_boundary: usize,
    /// Number of them that are not an edge of any polygon (ignored, with a warning).
    pub n_not_edges: usize,
    /// Number of edges that were cut in addition, to give a seam of a single edge two sides.
    pub n_extension_edges: usize,
    /// Number of connected pieces of the surface after the seams were cut.
    pub n_pieces: usize,
}

/// What validation found and did.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ValidationReport {
    /// Vertices that no polygon uses. They are ignored.
    pub n_loose_vertices: usize,
    /// Polygons whose winding was reversed to orient their component consistently.
    pub n_flipped_polygons: usize,
    /// Triangles without area.
    pub degenerate: DegenerateReport,
    /// Seams.
    pub seams: SeamReport,
}

/// How a closed piece was handled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClosedHandling {
    /// The piece has a boundary (of the mesh, or made by seams).
    NotClosed,
    /// The closed problem was solved as it stands.
    RawOptimum,
    /// Free cone sites were used.
    FreeSites {
        /// Number of sites the pipeline chose.
        n_automatic: usize,
        /// Number of pinned vertices of the caller that served as sites.
        n_pinned: usize,
        /// Whether this is the fallback of [`crate::ClosedSurfaces::Auto`]: the number of
        /// cones of the raw optimum that was discarded.
        fallback_from_cones: Option<usize>,
    },
    /// A slit was cut.
    Slit {
        /// Number of edges of the slit.
        n_edges: usize,
        /// Its length on the surface.
        length: f64,
    },
}

/// What cluster merging did on a piece.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MergeOutcome {
    /// Merging was not requested.
    NotRequested,
    /// No cluster was merged (there was none, or none within the size limit and the bounds).
    NothingToMerge,
    /// Clusters were merged and the merged cones are used.
    Merged,
    /// Clusters were merged, the result failed the guard, and the unmerged optimum is used.
    Rejected,
}

/// Summary of the cone solve of one piece.
#[derive(Debug, Clone, PartialEq)]
pub struct ConeSolveReport {
    /// Whether the solver certified its optimum.
    pub status: Status,
    /// The largest residual of the optimality system (state, adjoint, cone, bound, dual
    /// feasibility, Gauss–Bonnet) in the solver's certificate.
    pub max_residual: f64,
    /// The certificate's duality gap.
    pub duality_gap: f64,
    /// The certificate's complementarity gap (non-negative part of the gap).
    pub complementarity_gap: f64,
    /// The certificate's bound on what inexact residuals add to or take from the gap.
    pub gap_uncertainty: f64,
    /// The optimal objective `½ Σ m u² + λ Σ wR |µ|`.
    pub objective: f64,
    /// Newton steps of the regularized stage, of the exact stage, and safeguard rounds.
    pub newton_iterations: [usize; 3],
    /// Number of factorizations of the Newton matrix.
    pub factorizations: usize,
    /// Number of cones of the optimal measure, before merging.
    pub n_cones_optimal: usize,
    /// Number of edge-connected clusters of cones in the optimal measure.
    pub n_clusters: usize,
    /// What merging did.
    pub merge: MergeOutcome,
    /// The solver's distortion `½ Σ m u²` of the optimal measure.
    pub distortion_optimal: f64,
    /// The solver's distortion of the merged cones (equal to `distortion_optimal` if nothing
    /// was merged). Reported also when the merged cones were rejected.
    pub distortion_merged: f64,
    /// The solver's distortion of the cones that are used.
    pub distortion: f64,
    /// If the angles were re-optimized: the distortion before that step, and the number of
    /// cones before it (cones whose angle became zero are dropped).
    pub before_reoptimization: Option<(f64, usize)>,
    /// Wall-clock seconds inside the solver (both solves under
    /// [`crate::ClosedSurfaces::Auto`] with fallback).
    pub seconds: f64,
    /// The solver's own message, if it had one.
    pub message: Option<String>,
}

/// Statistics of a set of cones.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ConeStats {
    /// Number of cones.
    pub n_cones: usize,
    /// Number of cones with a negative angle.
    pub n_negative: usize,
    /// `Σ |µ|`.
    pub total_angle: f64,
    /// `Σ µ`.
    pub angle_sum: f64,
    /// The smallest cone angle (0 if there are no cones).
    pub min_angle: f64,
    /// The largest cone angle (0 if there are no cones).
    pub max_angle: f64,
}

impl ConeStats {
    /// The statistics of a list of cone angles.
    pub fn of(angles: impl Iterator<Item = f64>) -> ConeStats {
        let mut stats = ConeStats::default();
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for angle in angles {
            stats.n_cones += 1;
            stats.n_negative += usize::from(angle < 0.0);
            stats.total_angle += angle.abs();
            stats.angle_sum += angle;
            lo = lo.min(angle);
            hi = hi.max(angle);
        }
        if stats.n_cones > 0 {
            stats.min_angle = lo;
            stats.max_angle = hi;
        }
        stats
    }
}

/// What the flattener reported about a chart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlattenReport {
    /// Number of edges the cut planner added (cuts through cones, handle loops), without
    /// seams and slits.
    pub n_cut_edges: usize,
    /// Their total length on the surface.
    pub cut_length: f64,
    /// Number of triangles whose image is mirrored.
    pub flipped_triangles: usize,
    /// Number of triangles whose image has no area.
    pub degenerate_triangles: usize,
    /// Boundary vertices of the chart with more than a full turn around them (an exact sign of
    /// local overlap; it happens at cones of negative angle).
    pub wedges_over_2pi: usize,
    /// Pairs of boundary edges of the chart that cross each other.
    pub boundary_crossings: usize,
    /// The fraction of the chart's footprint in the plane that is covered more than once,
    /// estimated by sampling; `None` if the estimate is switched off.
    pub overlapped_fraction: Option<f64>,
    /// The largest quasi-conformal error `σ_max / σ_min` of a triangle (1 is conformal).
    pub quasi_conformal_max: f64,
    /// Its area-weighted mean.
    pub quasi_conformal_mean: f64,
    /// The largest relative change of a boundary length that closing the boundary curve
    /// needed.
    pub closure_length_change: f64,
    /// The largest relative difference of the lengths of the two sides of a cut edge.
    pub seam_length_mismatch: f64,
    /// The largest `|ℓ' − ℓ*| / ℓ*` over the boundary edges of the piece (the mesh's boundary,
    /// seams and slits), where `ℓ'` is the length of the edge in the layout and
    /// `ℓ* = e^{(u_i + u_j)/2} ℓ` the length the scale factors ask for. With optimized
    /// boundary values this says how well the flattening realizes the solver's values; with
    /// `u = 0` it says how well the lengths are kept. The flattener changes boundary lengths a
    /// little to close the boundary curve ([`FlattenReport::closure_length_change`]).
    pub boundary_length_error: f64,
    /// Area of the chart on the surface.
    pub surface_area: f64,
    /// Signed area of the chart in the plane, before packing.
    pub plane_area: f64,
    /// Area distortion of the chart by the measure of the dataset's README (see
    /// [`crate::measure`]), computed from the unpacked layout.
    pub area_distortion: AreaDistortion,
}

/// One chart: a connected piece of the surface after the seams were cut, flattened as one.
#[derive(Debug, Clone, PartialEq)]
pub struct ChartReport {
    /// Lowest-numbered polygon of the chart.
    pub first_polygon: u32,
    /// Number of polygons.
    pub n_polygons: usize,
    /// Number of triangles.
    pub n_triangles: usize,
    /// Genus of the piece.
    pub genus: usize,
    /// Number of boundary loops of the piece after the seams were cut (before any slit).
    pub n_boundary_loops: usize,
    /// How the piece was handled if it is closed.
    pub closed: ClosedHandling,
    /// Whether cone angles were restricted to positive values on this piece.
    pub positive_cones_only: bool,
    /// Whether the chart was solved and flattened with optimized boundary values
    /// ([`crate::BoundaryValues::Optimized`]): `false` if the option is off, the chart is
    /// closed or has no interior vertex, or the flattener had to drop the values (a warning
    /// says so).
    pub boundary_optimized: bool,
    /// The cone solve, or `None` if the piece has no interior vertex or no cone was possible
    /// (a closed piece of genus 1 with positive cones only).
    pub solve: Option<ConeSolveReport>,
    /// The spike check, or `None` if it was switched off, skipped, or the chart has no interior
    /// vertex.
    pub spikes: Option<SpikeReport>,
    /// The feedback loop of the cone solver, or `None` if it was switched off
    /// ([`crate::UnwrapOptions::cone_feedback`]) or the chart was not solved.
    pub feedback: Option<ChartFeedback>,
    /// The cones of the chart.
    pub cones: ConeStats,
    /// The flattening.
    pub flatten: FlattenReport,
}

/// Wall-clock seconds per stage.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Timings {
    /// Stage 1: checking the arrays, building and orienting the triangle mesh, edge lengths
    /// and geometry.
    pub validate: f64,
    /// Stage 2: cutting the seams and splitting into pieces.
    pub seams: f64,
    /// Stage 3a: cones (the cone solver, slits, the merge guard).
    pub cones: f64,
    /// Of that, inside the cone solver.
    pub cone_solver: f64,
    /// Of that, the spike check.
    pub spike_check: f64,
    /// Stage 3b: cutting through the cones and flattening.
    pub flatten: f64,
    /// Stage 3c: the feedback loop of the cone solver (its comparison, extra solves and
    /// flattenings); 0 if it is off.
    pub cone_feedback: f64,
    /// Stage 4: packing.
    pub pack: f64,
    /// The final checks and measurements.
    pub verify: f64,
    /// Everything.
    pub total: f64,
}

/// How well the two sides of the cut edges agree in length in the layout (before packing, which
/// scales all charts alike and changes no ratio).
///
/// An edge that was cut has one polygon on each side, in a chart of its own or in the same
/// chart, and each side has its own image. The relative mismatch of an edge is
/// `|ℓ_a − ℓ_b| / max(ℓ_a, ℓ_b)` for the lengths `ℓ_a`, `ℓ_b` of the two images.
///
/// With [`crate::BoundaryValues::Zero`] the scale factor is zero on a seam, so both sides keep
/// the length of the edge and the mismatch is what closing the boundary curve of each chart
/// changed (a fraction of a percent on the dataset). With
/// [`crate::BoundaryValues::Optimized`] the two pieces have different boundary values and the
/// mismatch is the difference of `e^{(u_i + u_j)/2}` between them; a texture painted across the
/// seam is stretched differently on its two sides by that much.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SeamLengths {
    /// Number of seam edges (given by the caller, and the extension of one-edge seams) that
    /// have a polygon on each side.
    pub n_seam_edges: usize,
    /// The largest mismatch over those, 0 if there are none.
    pub max_seam_mismatch: f64,
    /// The seam edge with the largest mismatch, as vertices of the caller's mesh `[min, max]`.
    pub worst_seam_edge: Option<[u32; 2]>,
    /// The mean mismatch over the seam edges, 0 if there are none.
    pub mean_seam_mismatch: f64,
    /// Number of the other cut edges (automatic slits, cuts through cones and around handles)
    /// that have a polygon on each side.
    pub n_other_cut_edges: usize,
    /// The largest mismatch over those. The flattener gives both sides of such an edge the same
    /// length, so this is at rounding level (but not for the slit of a closed surface, whose
    /// boundary is closed with the same adjustments as any other).
    pub max_other_mismatch: f64,
}

/// What happened to one chart that overlapped itself ([`crate::SplitOptions`]).
#[derive(Debug, Clone, PartialEq)]
pub struct SplitInfo {
    /// Lowest polygon of the chart that was cut (in the numbering of the run it was cut in).
    pub first_polygon: u32,
    /// The round (0 = the first re-run) in which the cut was made.
    pub round: usize,
    /// Number of edges of the cut; 0 if no cut separating the two layers was found.
    pub n_cut_edges: usize,
    /// Length of the cut on the surface.
    pub cut_length: f64,
    /// Overlapped fraction of the chart before the cut.
    pub overlap_before: f64,
    /// Overlapped area over footprint of the parts after the cut (0 if there was no cut).
    pub overlap_after: f64,
    /// Number of charts the parts form after the cut (those of the run after it).
    pub n_parts: usize,
    /// What became of the cut.
    pub outcome: SplitOutcome,
}

/// See [`SplitInfo`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitOutcome {
    /// The cut is part of the result.
    Kept,
    /// The cut did not lower the overlapped area enough, or mirrored a triangle: undone.
    Rejected,
    /// No cut was found (the overlapping parts are neighbours, or no finite cut separates
    /// them because of diagonals of the triangulation).
    NoCut,
}

/// What the automatic splitting of overlapping charts did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SplitReport {
    /// Runs of the pipeline, the first included.
    pub runs: usize,
    /// Cuts that are part of the result.
    pub n_kept: usize,
    /// Cuts that were undone.
    pub n_rejected: usize,
    /// Charts for which no cut was found.
    pub n_no_cut: usize,
    /// Whether overlapping charts remain that were not split because the budget
    /// (`max_splits`) or the rounds ran out.
    pub exhausted: bool,
    /// Total length of the kept cuts.
    pub cut_length: f64,
    /// The area-weighted overlapped fraction of the layout before any split, if estimated.
    pub overlap_before: Option<f64>,
    /// The largest overlapped fraction of a chart before any split, if estimated.
    pub max_chart_overlap_before: Option<f64>,
    /// One entry per cut tried, in order.
    pub splits: Vec<SplitInfo>,
}

/// The report of [`crate::unwrap`].
#[derive(Debug, Clone, PartialEq)]
pub struct Report {
    /// Number of vertices of the input (including loose ones).
    pub n_vertices: usize,
    /// Number of polygons.
    pub n_polygons: usize,
    /// Number of polygon corners.
    pub n_corners: usize,
    /// Number of triangles.
    pub n_triangles: usize,
    /// Stage 1 and 2.
    pub validation: ValidationReport,
    /// Steering: how many pinned, prescribed and forbidden vertices were applied (they are
    /// interior vertices of their piece) and how many had no effect because they lie on a
    /// boundary or a seam, where there are no cones, or are loose.
    pub steering: SteeringReport,
    /// One entry per chart, in chart order.
    pub charts: Vec<ChartReport>,
    /// The cones of the whole map.
    pub cones: ConeStats,
    /// Number of edges that were cut, seams included.
    pub n_cut_edges: usize,
    /// Area distortion of the whole map by the measure of the dataset's README, with one
    /// scale for all charts (see [`crate::measure`]).
    pub area_distortion: AreaDistortion,
    /// Totals over the charts: flipped triangles, triangles without area in the plane.
    pub flipped_triangles: usize,
    /// See `flipped_triangles`.
    pub degenerate_triangles: usize,
    /// The largest overlapped fraction of a chart, if it was estimated.
    pub max_overlapped_fraction: Option<f64>,
    /// The area-weighted overlapped fraction of the whole layout, if it was estimated: the sum
    /// of the overlapped areas of the charts over the sum of their footprints.
    pub overlapped_fraction: Option<f64>,
    /// Whether the two sides of the cut edges have the same length in the layout.
    pub seam_lengths: SeamLengths,
    /// Packing, if it ran.
    pub packing: Option<PackReport>,
    /// The automatic splitting of overlapping charts, if it was switched on.
    pub splitting: Option<SplitReport>,
    /// Where the time went.
    pub timings: Timings,
    /// Things the caller should know that did not stop the unwrap.
    pub warnings: Vec<String>,
}

/// How the steering was applied.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SteeringReport {
    /// Pinned vertices that became free cone sites.
    pub n_pinned_applied: usize,
    /// Pinned vertices without effect (on a boundary, on a seam or slit, or loose).
    pub n_pinned_without_effect: usize,
    /// Forbidden vertices that were excluded from carrying a cone.
    pub n_forbidden_applied: usize,
    /// Forbidden vertices without effect (no cone is possible there anyway).
    pub n_forbidden_without_effect: usize,
    /// Prescribed cones that were imposed.
    pub n_prescribed_applied: usize,
    /// The vertices of prescribed cones that could not be imposed, ascending: they lie on a
    /// boundary or seam (where the curvature is not a cone angle) or are loose.
    pub prescribed_not_applied: Vec<u32>,
}

/// A plain value of the flattened report ([`Report::summary`]).
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// A flag.
    Bool(bool),
    /// A count.
    Int(i64),
    /// A measurement.
    Float(f64),
    /// A name or status.
    Text(String),
}

impl Report {
    /// The report as an ordered list of named plain values, for bindings whose reports are
    /// dictionaries. Per-chart details are reduced to totals and extremes; keys that start
    /// with `seconds_` are wall-clock times. Warnings are not included (see
    /// [`Report::warnings`]).
    ///
    /// The keys are chosen not to collide with those the Python binding adds itself
    /// (`n_vertices`, `n_charts`, `n_cones`, `n_cut_edges`, `n_degenerate_triangles`,
    /// `n_flipped_uv_triangles`, `measured_*`, `seconds_pipeline`, `seconds_validate`, …);
    /// `n_unreferenced_vertices`
    /// and `n_seam_edges_ignored` keep the names the binding's stand-in pipeline used. The
    /// two overlap figures are left out when the estimate is switched off, and the packing
    /// figures when packing is. No value is NaN.
    pub fn summary(&self) -> Vec<(&'static str, Value)> {
        let int = |n: usize| Value::Int(n as i64);
        let v = &self.validation;
        let solves: Vec<&ConeSolveReport> = self
            .charts
            .iter()
            .filter_map(|c| c.solve.as_ref())
            .collect();
        let fold_max = |values: &mut dyn Iterator<Item = f64>| values.fold(0.0f64, f64::max);
        let mut out = vec![
            ("n_unreferenced_vertices", int(v.n_loose_vertices)),
            ("n_flipped_faces", int(v.n_flipped_polygons)),
            (
                "n_zero_area_triangles",
                int(v.degenerate.n_degenerate_triangles),
            ),
            ("mollified", Value::Bool(v.degenerate.mollified)),
            ("mollify_epsilon", Value::Float(v.degenerate.epsilon)),
            ("n_seam_edges_cut", int(v.seams.n_cut)),
            ("n_seam_edges_on_boundary", int(v.seams.n_on_boundary)),
            ("n_seam_edges_ignored", int(v.seams.n_not_edges)),
            ("n_seam_extension_edges", int(v.seams.n_extension_edges)),
            ("n_pinned_applied", int(self.steering.n_pinned_applied)),
            (
                "n_prescribed_applied",
                int(self.steering.n_prescribed_applied),
            ),
            (
                "n_prescribed_not_applied",
                int(self.steering.prescribed_not_applied.len()),
            ),
            (
                "n_forbidden_applied",
                int(self.steering.n_forbidden_applied),
            ),
            ("n_closed_pieces", int(self.n_closed_pieces())),
            ("n_slits", int(self.n_slits())),
            ("n_negative_cones", int(self.cones.n_negative)),
            ("total_cone_angle", Value::Float(self.cones.total_angle)),
            ("max_cone_angle", Value::Float(self.cones.max_angle)),
            ("min_cone_angle", Value::Float(self.cones.min_angle)),
            (
                "cones_certified",
                Value::Bool(solves.iter().all(|s| s.status == Status::Converged)),
            ),
            (
                "cone_solver_max_residual",
                Value::Float(fold_max(&mut solves.iter().map(|s| s.max_residual))),
            ),
            (
                "cone_solver_max_gap",
                Value::Float(fold_max(
                    &mut solves
                        .iter()
                        .map(|s| s.complementarity_gap.abs() + s.gap_uncertainty),
                )),
            ),
            (
                "cone_solver_factorizations",
                int(solves.iter().map(|s| s.factorizations).sum()),
            ),
            (
                "n_merges_rejected",
                int(solves
                    .iter()
                    .filter(|s| s.merge == MergeOutcome::Rejected)
                    .count()),
            ),
            (
                "area_distortion",
                Value::Float(self.area_distortion.scale_free),
            ),
            (
                "area_distortion_unscaled",
                Value::Float(self.area_distortion.raw),
            ),
            ("flatten_flipped_triangles", int(self.flipped_triangles)),
            (
                "flatten_collapsed_triangles",
                int(self.degenerate_triangles),
            ),
            (
                "max_quasi_conformal_error",
                Value::Float(fold_max(
                    &mut self.charts.iter().map(|c| c.flatten.quasi_conformal_max),
                )),
            ),
        ];
        out.push((
            "boundary_values_optimized",
            Value::Bool(self.charts.iter().any(|c| c.boundary_optimized)),
        ));
        out.push((
            "max_boundary_length_error",
            Value::Float(fold_max(
                &mut self.charts.iter().map(|c| c.flatten.boundary_length_error),
            )),
        ));
        out.push((
            "max_seam_length_mismatch",
            Value::Float(self.seam_lengths.max_seam_mismatch),
        ));
        out.push((
            "mean_seam_length_mismatch",
            Value::Float(self.seam_lengths.mean_seam_mismatch),
        ));
        out.push((
            "max_other_cut_length_mismatch",
            Value::Float(self.seam_lengths.max_other_mismatch),
        ));
        // Present only if the spike check ran.
        let checks: Vec<&SpikeReport> = self
            .charts
            .iter()
            .filter_map(|c| c.spikes.as_ref())
            .collect();
        if !checks.is_empty() {
            let spikes: Vec<u32> = {
                let mut all: Vec<u32> = checks
                    .iter()
                    .flat_map(|s| s.vertices.iter().copied())
                    .collect();
                all.sort_unstable();
                all
            };
            out.push(("n_spike_vertices", int(spikes.len())));
            out.push((
                "cone_free_distortion",
                Value::Float(fold_max(&mut checks.iter().map(|s| s.cone_free_distortion))),
            ));
            if !spikes.is_empty() {
                out.push(("spike_vertices", Value::Text(spike_list(&spikes))));
            }
        }
        // Present only if the feedback loop ran.
        let loops: Vec<&ChartFeedback> = self
            .charts
            .iter()
            .filter_map(|c| c.feedback.as_ref())
            .collect();
        if !loops.is_empty() {
            let mut reweighted: Vec<u32> = loops
                .iter()
                .flat_map(|f| f.reweighted.iter().copied())
                .collect();
            reweighted.sort_unstable();
            out.push((
                "cone_feedback_rounds",
                int(loops.iter().map(|f| f.report.rounds.len()).max().unwrap_or(0)),
            ));
            out.push(("n_cone_feedback_reweighted", int(reweighted.len())));
            out.push((
                "n_cone_feedback_pending",
                int(loops.iter().map(|f| f.pending.len()).sum()),
            ));
            out.push((
                "cone_feedback_max_mismatch",
                Value::Float(fold_max(&mut loops.iter().filter_map(|f| {
                    f.report.rounds.first().map(|r| r.max_abs_eligible)
                }))),
            ));
            out.push((
                "cone_feedback_final_mismatch",
                Value::Float(fold_max(&mut loops.iter().filter_map(|f| {
                    f.report.rounds.last().map(|r| r.max_abs_eligible)
                }))),
            ));
            if !reweighted.is_empty() {
                out.push((
                    "cone_feedback_vertices",
                    Value::Text(spike_list(&reweighted)),
                ));
            }
        }
        // Present only if the overlap was estimated: a NaN would break comparisons of reports.
        if let Some(fraction) = self.overlapped_fraction {
            out.push(("overlapped_fraction", Value::Float(fraction)));
        }
        if let Some(fraction) = self.max_overlapped_fraction {
            out.push(("max_chart_overlapped_fraction", Value::Float(fraction)));
        }
        if let Some(split) = &self.splitting {
            out.push(("n_splits", int(split.n_kept)));
            out.push(("n_splits_rejected", int(split.n_rejected)));
            out.push(("n_split_no_cut", int(split.n_no_cut)));
            out.push(("split_runs", int(split.runs)));
            out.push(("split_exhausted", Value::Bool(split.exhausted)));
            out.push(("split_cut_length", Value::Float(split.cut_length)));
            if let Some(x) = split.overlap_before {
                out.push(("overlapped_fraction_before_split", Value::Float(x)));
            }
            if let Some(x) = split.max_chart_overlap_before {
                out.push(("max_chart_overlapped_fraction_before_split", Value::Float(x)));
            }
        }
        if let Some(pack) = &self.packing {
            out.push(("packed", Value::Bool(true)));
            out.push(("pack_scale", Value::Float(pack.scale)));
            out.push(("pack_margin", Value::Float(pack.margin)));
            out.push(("pack_fill", Value::Float(pack.fill)));
        } else {
            out.push(("packed", Value::Bool(false)));
        }
        let t = &self.timings;
        out.extend([
            ("seconds_stage_validate", Value::Float(t.validate)),
            ("seconds_stage_seams", Value::Float(t.seams)),
            ("seconds_stage_cones", Value::Float(t.cones)),
            ("seconds_cone_solver", Value::Float(t.cone_solver)),
            ("seconds_spike_check", Value::Float(t.spike_check)),
            ("seconds_stage_flatten", Value::Float(t.flatten)),
            ("seconds_cone_feedback", Value::Float(t.cone_feedback)),
            ("seconds_stage_pack", Value::Float(t.pack)),
            ("seconds_stage_verify", Value::Float(t.verify)),
            ("seconds_unwrap", Value::Float(t.total)),
        ]);
        out
    }

    /// Number of pieces that were closed after the seams were cut.
    pub fn n_closed_pieces(&self) -> usize {
        self.charts
            .iter()
            .filter(|c| c.closed != ClosedHandling::NotClosed)
            .count()
    }

    /// Number of pieces that received an automatic slit.
    pub fn n_slits(&self) -> usize {
        self.charts
            .iter()
            .filter(|c| matches!(c.closed, ClosedHandling::Slit { .. }))
            .count()
    }
}

/// At most 64 vertex numbers, comma separated, then the number of those left out.
fn spike_list(vertices: &[u32]) -> String {
    const LISTED: usize = 64;
    let mut text = vertices
        .iter()
        .take(LISTED)
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    if vertices.len() > LISTED {
        text.push_str(&format!(", and {} more", vertices.len() - LISTED));
    }
    text
}
