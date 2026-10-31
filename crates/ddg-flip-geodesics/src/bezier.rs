//! Geodesic Bézier curves by midpoint subdivision (paper Sec. 6.2, Fig. 23, after Morera,
//! Carvalho and Velho 2008; SPEC §6.1, F19).
//!
//! The curve is one open path whose pinned vertices are the control points. A **round** of
//! the subdivision on the control points `P₀ … P_d` of a piece is de Casteljau's algorithm at
//! `½` with geodesics in place of straight lines (the reference's recursion,
//! `gc: fg.cpp:1068-1248`):
//!
//! 1. Shorten (Alg. 2, [`GeodesicNetwork::shorten`]): every part of the path between two
//!    consecutive pinned points is a geodesic.
//! 2. In the active range (at first `P₀ … P_d`), for each part of length `L` insert its
//!    midpoint: walk its segments to the one where the arc length `L/2` falls and split that
//!    edge there, at `t = (L/2 − s)/ℓ` along the segment (`s` the length before it; the split
//!    parameter of SPEC §8, `fo-bezier-midpoint-split`, to be proved with the Lean work of F5);
//!    pin the new vertex. If the point is within `snap_tolerance · ℓ` of either end of the
//!    edge, that vertex is pinned instead (F19; the reference checks only the far end,
//!    `gc: fg.cpp:1195-1229`). Unpin the old control points inside the active range (not its
//!    first and last; `gc: fg.cpp:1164-1183`), make the new points the active range, shorten,
//!    and repeat until the active range has one part; its midpoint is `B(½)`.
//! 3. Recurse on `P₀ … B(½)` and `B(½) … P_d` with one round less.
//!
//! After `r` rounds the pinned points are the control polygons of the `2^r` pieces, in curve
//! order (`d 2^r + 1` points); the path through them converges to the geodesic Bézier curve
//! (Morera et al. 2008, `fo-geodesic-bezier`, cited). In the plane they are the points of
//! planar de Casteljau subdivision (`fo-de-casteljau-midpoint`; `tests/bezier.rs`).
//!
//! **Refused splits.** The triangulation may refuse a midpoint's split. When a new triangle
//! would be too thin (IT N5; a geodesic across a flat region is often one long edge between
//! needles), the triangulation is flipped to Delaunay with the edges of the curve held fixed
//! and the split is made again ([`BezierOptions::delaunay_on_thin_split`]). When the integer
//! representation cannot place the vertex because it lies within rounding of an input edge
//! (it places such vertices strictly inside input faces), the point is moved by `2e-9 ℓ`
//! along the edge, towards its middle; this is decided before the split from the crossings of the edge, which
//! both representations compute ([`BezierOptions::unlocated_nudge`]). Both are counted in the report; the reference has
//! neither refusal.
//!
//! **Bookkeeping.** The control points are kept as a list in curve order and checked, after
//! every change, against the pinned vertices met along the path: a curve that passes through
//! a control point a second time (pinning is per vertex, F10) cannot be subdivided and is an
//! error, never repaired (the reference assumes a simple curve, `gc: fg.cpp:1095`). Splits go
//! through [`GeodesicNetwork::split_edge`], which keeps the place of every segment on the
//! split edge.
//!
//! **Representations.** Every decision (the part lengths, the segment, `t`, the snap) is made
//! from the edge lengths and the connectivity alone, and splits give bitwise the same lengths
//! on both representations (IT N10, R4), so signposts and integer coordinates make the same
//! splits, flips and segments (SPEC §7 test 8). Only the place of a new vertex on the input
//! mesh is computed differently.

use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, InsertRefusal, IntrinsicError, Triangulation,
};

use crate::error::FlipError;
use crate::init::waypoint_path;
use crate::network::GeodesicNetwork;
use crate::report::{ShortenOptions, ShortenReport, ShortenStatus};

/// More rounds than this are refused: `d 2^r + 1` control points.
pub const MAX_BEZIER_ROUNDS: usize = 24;

/// Options of [`geodesic_bezier`] and [`GeodesicNetwork::bezier_subdivide`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BezierOptions {
    /// The options of every shortening (Alg. 2).
    pub shorten: ShortenOptions,
    /// A midpoint within `snap_tolerance · ℓ` of an end of the edge it falls on (`ℓ` that
    /// edge's length) is that end (F19; `gc: fg.cpp:1186`, `useVertexEPS`). Default `1e-4`;
    /// in `[0, ½)`.
    pub snap_tolerance: f64,
    /// When the triangulation refuses a midpoint's split because a new triangle would be too
    /// thin (IT N5, [`InsertRefusal::TooThin`]): flip the triangulation to Delaunay with every
    /// edge that carries a segment held fixed (the iCDT of [`crate::cdt`], the marks cleared
    /// again), which leaves the curve as it is and fattens the triangles at the edge, then
    /// split again (once). Without it the refusal is the error. Default on. The reference has
    /// no such refusal; a geodesic across a flat region is often one long edge between
    /// needles.
    pub delaunay_on_thin_split: Option<DelaunayOptions>,
    /// A midpoint within `unlocated_nudge · ℓ` of a point where the split edge crosses an input
    /// edge is moved along the edge by `2 unlocated_nudge · ℓ`, towards the middle of the
    /// edge. The integer representation places a new vertex strictly inside
    /// an input face and refuses a point on an input edge ([`InsertRefusal::NotLocated`],
    /// IC-SPEC §4.4), which lattice-like meshes produce often; the crossings come from
    /// [`Triangulation::trace_edge_on_input`] on both representations, so both move the same
    /// points. A split still refused that way is made again once at `t ± unlocated_nudge`.
    /// Default `1e-9`, far below the snap tolerance; `0` turns both off. In `[0, 1e-3]`.
    pub unlocated_nudge: f64,
}

impl Default for BezierOptions {
    fn default() -> Self {
        BezierOptions {
            shorten: ShortenOptions::default(),
            snap_tolerance: 1e-4,
            delaunay_on_thin_split: Some(DelaunayOptions::default()),
            unlocated_nudge: 1e-9,
        }
    }
}

impl BezierOptions {
    /// Checks the ranges.
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`] naming the option.
    pub fn validate(&self) -> Result<(), FlipError> {
        self.shorten.validate()?;
        let s = self.snap_tolerance;
        if !(s.is_finite() && (0.0..0.5).contains(&s)) {
            return Err(FlipError::InvalidOptions {
                what: format!("snap_tolerance must be in [0, 0.5), got {s}"),
            });
        }
        let x = self.unlocated_nudge;
        if !(x.is_finite() && (0.0..=1e-3).contains(&x)) {
            return Err(FlipError::InvalidOptions {
                what: format!("unlocated_nudge must be in [0, 1e-3], got {x}"),
            });
        }
        Ok(())
    }
}

/// How a midpoint was made.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MidpointKind {
    /// The edge was split there.
    Split,
    /// The point was within the snap tolerance of the end where the segment starts (in path
    /// order); that vertex was pinned.
    SnappedToTail,
    /// The point was within the snap tolerance of the end where the segment ends; that
    /// vertex was pinned.
    SnappedToHead,
}

/// One midpoint of the subdivision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BezierMidpoint {
    /// The round it belongs to (`1` for the first round, i.e. `B(½)` of the whole curve and
    /// the points made on the way to it).
    pub round: usize,
    /// The vertex pinned for it (the new vertex of a split).
    pub vertex: u32,
    /// The edge it fell on (before the split: the edge that was split).
    pub edge: u32,
    /// The parameter of the point along that edge from the tail of its halfedge (the split
    /// parameter for [`MidpointKind::Split`]).
    pub t: f64,
    /// Split or snapped.
    pub kind: MidpointKind,
    /// The length `L` of the geodesic part it is the midpoint of.
    pub part_length: f64,
}

/// What the subdivision did.
#[derive(Debug, Clone, PartialEq)]
pub struct BezierReport {
    /// The pinned points at the end, in curve order: the control polygons of the `2^r`
    /// pieces, consecutive pieces sharing their end point (`d 2^r + 1` vertices).
    pub control_points: Vec<usize>,
    /// Every midpoint, in the order they were made.
    pub midpoints: Vec<BezierMidpoint>,
    /// Edge splits made.
    pub n_splits: usize,
    /// Midpoints snapped to an existing vertex.
    pub n_snapped: usize,
    /// Splits refused as too thin and made after constrained Delaunay flipping
    /// ([`BezierOptions::delaunay_on_thin_split`]).
    pub n_delaunay_retries: usize,
    /// Flips made by those Delaunay flippings.
    pub delaunay_flips: usize,
    /// Midpoints moved off a crossing with an input edge before the split, or after a split
    /// refused as not locatable ([`BezierOptions::unlocated_nudge`]).
    pub n_nudged: usize,
    /// Shortenings run (Alg. 2 calls).
    pub n_shortenings: usize,
    /// Shortenings that did not end [`ShortenStatus::Geodesic`].
    pub n_not_geodesic: usize,
    /// The status of every shortening that did not end `Geodesic`, in order (empty when all
    /// did).
    pub other_statuses: Vec<ShortenStatus>,
    /// Iterations (Alg. 1 calls) over all shortenings.
    pub iterations: usize,
    /// Flips over all shortenings.
    pub flips: usize,
    /// Refused flips over all shortenings.
    pub refused_flips: usize,
    /// Stalled joints over all shortenings.
    pub stalled: usize,
    /// The length of the path before the first shortening.
    pub initial_length: f64,
    /// The length at the end.
    pub final_length: f64,
}

impl BezierReport {
    fn new(initial_length: f64) -> BezierReport {
        BezierReport {
            control_points: Vec::new(),
            midpoints: Vec::new(),
            n_splits: 0,
            n_snapped: 0,
            n_delaunay_retries: 0,
            delaunay_flips: 0,
            n_nudged: 0,
            n_shortenings: 0,
            n_not_geodesic: 0,
            other_statuses: Vec::new(),
            iterations: 0,
            flips: 0,
            refused_flips: 0,
            stalled: 0,
            initial_length,
            final_length: initial_length,
        }
    }

    fn add(&mut self, r: &ShortenReport) {
        self.n_shortenings += 1;
        if r.status != ShortenStatus::Geodesic {
            self.n_not_geodesic += 1;
            self.other_statuses.push(r.status);
        }
        self.iterations += r.iterations;
        self.flips += r.flips;
        self.refused_flips += r.refused_flips;
        self.stalled += r.stalled;
    }
}

/// A geodesic Bézier curve with control points `control_vertices` (`P₀ … P_d`, `d ≥ 1`) on
/// `tri`, subdivided `rounds` times (SPEC §6.1): the initial path is [`waypoint_path`] through
/// the control points with every control point pinned, then
/// [`GeodesicNetwork::bezier_subdivide`]. The network returned has the curve as its one path;
/// its pinned vertices are [`BezierReport::control_points`].
///
/// # Errors
/// [`FlipError::InvalidOptions`] for options out of range or more than
/// [`MAX_BEZIER_ROUNDS`] rounds; the errors of [`waypoint_path`] (fewer than two control
/// points, a vertex out of range, two consecutive equal ones, not connected); the errors of
/// [`GeodesicNetwork::new`] (a path that crosses itself); and those of
/// [`GeodesicNetwork::bezier_subdivide`].
pub fn geodesic_bezier<T: Triangulation>(
    tri: T,
    control_vertices: &[usize],
    rounds: usize,
    options: &BezierOptions,
) -> Result<(GeodesicNetwork<T>, BezierReport), FlipError> {
    check_rounds(rounds)?;
    options.validate()?;
    let initial = waypoint_path(&tri, control_vertices, false, true)?;
    let mut network = GeodesicNetwork::new(tri, &initial.paths, &initial.pinned, None)?;
    let report = network.bezier_subdivide(rounds, options)?;
    Ok((network, report))
}

fn check_rounds(rounds: usize) -> Result<(), FlipError> {
    if rounds > MAX_BEZIER_ROUNDS {
        return Err(FlipError::InvalidOptions {
            what: format!("at most {MAX_BEZIER_ROUNDS} Bézier rounds, got {rounds}"),
        });
    }
    Ok(())
}

/// The curve as it is: its segments' lengths and vertices in path order.
struct Curve {
    /// `n + 1` vertices for `n` segments.
    vertices: Vec<usize>,
    /// Per segment: edge, forward.
    segments: Vec<(usize, bool)>,
}

impl<T: Triangulation> GeodesicNetwork<T> {
    /// Subdivides the network's curve as a geodesic Bézier curve `rounds` times (SPEC §6.1, the
    /// [module documentation](crate::bezier)): the network must be one open path, and its
    /// pinned vertices met along it, in order and including both ends, are the control points
    /// `P₀ … P_d`. The first step is a shortening, so the path need not be geodesic yet;
    /// `rounds = 0` only shortens.
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`] for options out of range or more than
    /// [`MAX_BEZIER_ROUNDS`] rounds, [`FlipError::Bezier`] for a network that is not one open
    /// path or whose control points repeat (nothing is changed in these cases);
    /// [`FlipError::Bezier`] when a shortening leads the path through a control point a
    /// second time; the errors of [`shorten`](Self::shorten) and of
    /// [`split_edge`](Self::split_edge) (a split the triangulation refuses). After an error
    /// that is not a refused input the network is valid but partly subdivided.
    pub fn bezier_subdivide(
        &mut self,
        rounds: usize,
        options: &BezierOptions,
    ) -> Result<BezierReport, FlipError> {
        check_rounds(rounds)?;
        options.validate()?;
        let net = &self.network;
        if net.n_paths() != 1 || net.is_closed(0) || net.n_path_segments(0) == 0 {
            return Err(FlipError::Bezier {
                reason: format!(
                    "the network must be one open path; it has {} paths{}",
                    net.n_paths(),
                    if net.n_paths() == 1 {
                        if net.is_closed(0) {
                            ", closed"
                        } else {
                            ", contracted"
                        }
                    } else {
                        ""
                    }
                ),
            });
        }
        let curve = self.curve();
        let mut controls: Vec<usize> = curve
            .vertices
            .iter()
            .copied()
            .filter(|&v| net.is_pinned(v))
            .collect();
        // A control point met twice (other than P₀ = P_d) cannot be unpinned in one place
        // only.
        let n = controls.len();
        let mut sorted = controls[1..n - 1].to_vec();
        sorted.push(controls[0]);
        if controls[n - 1] != controls[0] {
            sorted.push(controls[n - 1]);
        }
        sorted.sort_unstable();
        if let Some(w) = sorted.windows(2).find(|w| w[0] == w[1]) {
            return Err(FlipError::Bezier {
                reason: format!(
                    "the path meets the pinned vertex {} more than once (pinning is per \
                     vertex, so its control points must be distinct, except the first may be \
                     the last)",
                    w[0]
                ),
            });
        }

        let mut report = BezierReport::new(self.length);
        let r = self.shorten(&options.shorten)?;
        report.add(&r);
        self.check_controls(&controls)?;
        let last = controls.len() - 1;
        self.subdivide_range(rounds, 1, 0, last, &mut controls, options, &mut report)?;
        report.final_length = self.length;
        report.control_points = controls;
        Ok(report)
    }

    /// The vertices and segments of path 0 in order.
    fn curve(&self) -> Curve {
        let complex = self.tri.complex();
        let path = self.network.path(0);
        Curve {
            vertices: path.vertices(complex),
            segments: path
                .segments
                .iter()
                .map(|s| (s.edge as usize, s.forward))
                .collect(),
        }
    }

    /// The positions in the curve's vertex list of its control points, checking that the
    /// pinned vertices met along the path are exactly `controls`, in order.
    fn check_controls(&self, controls: &[usize]) -> Result<Vec<usize>, FlipError> {
        let curve = self.curve();
        let mut at = Vec::with_capacity(controls.len());
        for (i, &v) in curve.vertices.iter().enumerate() {
            if !self.network.is_pinned(v) {
                continue;
            }
            match controls.get(at.len()) {
                Some(&c) if c == v => at.push(i),
                expected => {
                    return Err(FlipError::Bezier {
                        reason: format!(
                            "the path meets the pinned vertex {v} where control point {} \
                             ({expected:?}) was expected: the curve passes through a control \
                             point a second time",
                            at.len()
                        ),
                    });
                }
            }
        }
        if at.len() != controls.len() {
            return Err(FlipError::Bezier {
                reason: format!(
                    "the path meets {} of its {} control points",
                    at.len(),
                    controls.len()
                ),
            });
        }
        Ok(at)
    }

    /// The reference's `bezierSubdivideRecursive` on the control points `controls[a ..= b]`
    /// (`gc: fg.cpp:1084-1248`). Returns the new index of the last one.
    #[allow(clippy::too_many_arguments)]
    fn subdivide_range(
        &mut self,
        rounds: usize,
        round: usize,
        a: usize,
        b: usize,
        controls: &mut Vec<usize>,
        options: &BezierOptions,
        report: &mut BezierReport,
    ) -> Result<usize, FlipError> {
        if rounds == 0 {
            return Ok(b);
        }
        let len_before = controls.len();
        let (mut lo, mut hi) = (a, b);
        let mid = loop {
            let parts = hi - lo;
            // The midpoint of every part, in curve order; each is inserted into `controls`
            // right away, so part `j` lies between `controls[lo + 2j]` and
            // `controls[lo + 2j + 1]` at the time it is split.
            for j in 0..parts {
                let v = self.insert_midpoint(lo + 2 * j, round, controls, options, report)?;
                controls.insert(lo + 2 * j + 1, v);
            }
            hi = lo + 2 * parts;
            if parts == 1 {
                let r = self.shorten(&options.shorten)?;
                report.add(&r);
                self.check_controls(controls)?;
                break lo + 1;
            }
            // Unpin the old control points inside the active range (`lo + 2, lo + 4, …,
            // hi − 2`); the range keeps its first and last.
            let complex = self.tri.complex();
            for k in (1..parts).rev() {
                let v = controls.remove(lo + 2 * k);
                self.network.unpin(complex, v)?;
            }
            hi -= parts - 1;
            let r = self.shorten(&options.shorten)?;
            report.add(&r);
            self.check_controls(controls)?;
            // The new points are the active range.
            lo += 1;
            hi -= 1;
        };
        let b = b + (controls.len() - len_before);
        let mid_after =
            self.subdivide_range(rounds - 1, round + 1, a, mid, controls, options, report)?;
        let b = b + (mid_after - mid);
        self.subdivide_range(
            rounds - 1,
            round + 1,
            mid_after,
            b,
            controls,
            options,
            report,
        )
    }

    /// The parameter to split edge `e` at instead of `t` when the point `t ℓ` is within
    /// `nudge · ℓ` of a crossing of `e` with an input edge ([`EdgeTrace`]'s `distance`s, which
    /// both representations compute to rounding): `t ± 2 nudge`, towards the middle of `e`
    /// (computed from `t` alone, so both representations get the same bits; the point is then
    /// between `nudge ℓ` and `3 nudge ℓ` from the crossing). `None` when no crossing is that
    /// close, `nudge` is `0`, or the edge cannot be traced (a signpost trace that fails IT N8;
    /// the split decides then).
    ///
    /// [`EdgeTrace`]: ddg_intrinsic::EdgeTrace
    fn off_crossings(&self, e: usize, t: f64, nudge: f64) -> Option<f64> {
        if nudge <= 0.0 {
            return None;
        }
        let l = self.tri.edge_length(e);
        let at = t * l;
        let trace = self
            .tri
            .trace_edge_on_input(e, &EdgeTraceOptions::default())
            .ok()?;
        trace
            .crossings
            .iter()
            .any(|c| (at - c.distance).abs() <= nudge * l)
            .then(|| towards_middle(t, 2.0 * nudge))
    }

    /// Finds the midpoint of the geodesic part from `controls[i]` to `controls[i + 1]`, splits
    /// its edge there (or snaps to an end), pins it and returns its vertex.
    fn insert_midpoint(
        &mut self,
        i: usize,
        round: usize,
        controls: &[usize],
        options: &BezierOptions,
        report: &mut BezierReport,
    ) -> Result<usize, FlipError> {
        let at = self.check_controls(controls)?;
        let curve = self.curve();
        let part = &curve.segments[at[i]..at[i + 1]];
        let lengths = self.tri.edge_lengths();
        let total = part
            .iter()
            .map(|&(e, _)| lengths[e])
            .fold(0.0, |s, l| s + l);
        let half = 0.5 * total;
        // The first segment whose end reaches `L/2` (the last one at worst).
        let mut before = 0.0;
        let mut k = 0;
        while k + 1 < part.len() {
            let next = before + lengths[part[k].0];
            if next >= half {
                break;
            }
            before = next;
            k += 1;
        }
        let (e, forward) = part[k];
        let l = lengths[e];
        let u = (half - before) / l;
        let mut t = if forward { u } else { 1.0 - u };
        let snap = options.snap_tolerance;
        let (vertex, kind) = if u <= snap {
            (curve.vertices[at[i] + k], MidpointKind::SnappedToTail)
        } else if u >= 1.0 - snap {
            (curve.vertices[at[i] + k + 1], MidpointKind::SnappedToHead)
        } else {
            if let Some(moved) = self.off_crossings(e, t, options.unlocated_nudge) {
                t = moved;
                report.n_nudged += 1;
            }
            let mut split = self.split_edge(e, t);
            if let (
                Err(FlipError::Intrinsic(IntrinsicError::InsertionRefused { reason, .. })),
                Some(delaunay),
            ) = (&split, &options.delaunay_on_thin_split)
                && matches!(reason, InsertRefusal::TooThin { .. })
            {
                // Edges that carry segments keep their ids and directions, so `e` and `t`
                // still name the point.
                let cdt = self.constrained_delaunay(delaunay, false)?;
                report.n_delaunay_retries += 1;
                report.delaunay_flips += cdt.delaunay.n_flips;
                split = self.split_edge(e, t);
            }
            if let Err(FlipError::Intrinsic(IntrinsicError::InsertionRefused { reason, .. })) =
                &split
                && matches!(reason, InsertRefusal::NotLocated { .. })
                && options.unlocated_nudge > 0.0
            {
                // Still within rounding of an input edge (not caught by `off_crossings`): move
                // it along the edge.
                t = towards_middle(t, options.unlocated_nudge);
                report.n_nudged += 1;
                split = self.split_edge(e, t);
            }
            (split?.insertion.vertex as usize, MidpointKind::Split)
        };
        if kind == MidpointKind::Split {
            report.n_splits += 1;
        } else {
            report.n_snapped += 1;
        }
        self.network.pin(vertex);
        report.midpoints.push(BezierMidpoint {
            round,
            vertex: vertex as u32,
            edge: e as u32,
            t,
            kind,
            part_length: total,
        });
        Ok(vertex)
    }
}

/// `t` moved by `delta` towards `½`.
fn towards_middle(t: f64, delta: f64) -> f64 {
    if t < 0.5 { t + delta } else { t - delta }
}
