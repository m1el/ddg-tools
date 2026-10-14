//! Tracing queries (paper Sec. 3.2.2, Alg. 2): follow a straight path over a triangulation.
//!
//! # Method
//!
//! The triangles along the path are unfolded one after the other into a single plane in which
//! the start point is the origin and the path is the positive x-axis. Entering a triangle
//! through an edge, its third vertex is placed from the two end points of that edge and the
//! three edge lengths ([`unfold_apex`]); the path leaves through the edge on whose other side
//! of the axis that vertex lies, i.e. the decision is the sign of one coordinate. The crossing
//! point is the intersection of the axis with that edge ([`axis_crossing`]). The direction of
//! the path is never transformed from one triangle to the next, so the only quantity that
//! accumulates rounding error is the position of the unfolded vertices.
//!
//! This is the unfolding the paper describes; the authors' implementation instead carries a
//! direction vector from face to face in barycentric coordinates and projects it back into
//! the valid range where rounding pushes it out. Nothing is projected or perturbed here
//! (SPEC N3). The only values moved into their range are reported in [`Trace::clamp`].
//!
//! # Paths through vertices
//!
//! A path that runs exactly along an edge is supported: it ends on that edge or at its far
//! vertex. A path that arrives at a vertex with length left over stops there
//! ([`TraceStop::Vertex`]); it is not continued, because the continuation through a cone
//! point is a matter of convention that no operation of this crate needs.

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::corner_angle;
use ddg_core::topology::NONE;

use crate::error::IntrinsicError;
use crate::formulas::{angle_scale, axis_crossing, barycentric, unfold_apex};
use crate::point::{SurfacePoint, point_defect};
use crate::surface::{Surface, edge_parameter, face_lengths, is_boundary_vertex, next, prev};

/// Options of a tracing query.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TraceOptions {
    /// Whether [`Trace::crossings`] is filled.
    pub record_crossings: bool,
    /// The query fails with [`IntrinsicError::CrossingLimit`] when the path has crossed this
    /// many edges. `None` stands for `100 · (number of faces of the surface traced over) +
    /// 100`. In exact arithmetic every path ends by itself; the limit bounds the work if
    /// rounding should ever keep one from advancing.
    pub max_crossings: Option<usize>,
}

/// Why a traced path ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceStop {
    /// The requested length was travelled.
    Completed,
    /// The path left the surface through a boundary edge, or the direction at a boundary
    /// vertex pointed out of the surface. [`Trace::length`] is the length travelled.
    Boundary,
    /// The path arrived at a vertex with length left over (see the
    /// [module documentation](self)). [`Trace::length`] is the length travelled.
    Vertex,
}

/// A point where a traced path crosses an edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing {
    /// The halfedge through which the path leaves a face: the path goes from the face of this
    /// halfedge to the face of its twin, from left to right of the halfedge's direction.
    pub halfedge: u32,
    /// The edge that is crossed.
    pub edge: u32,
    /// The parameter of the crossing point along the edge, in `[0, 1]`, in the direction of
    /// the edge (see [`SurfacePoint::Edge`]).
    pub t: f64,
    /// The distance of the crossing point from the start of the path, along the path.
    /// Distances do not decrease from one crossing to the next, up to rounding.
    pub distance: f64,
}

/// The result of a tracing query.
#[derive(Debug, Clone, PartialEq)]
pub struct Trace {
    /// Where the path ended. A point of a face in general; a point of an edge or a vertex
    /// when the path ran along an edge; a point of a boundary edge when it left the surface.
    pub end: SurfacePoint,
    /// Why it ended.
    pub stop: TraceStop,
    /// The length travelled: the requested length when `stop` is [`TraceStop::Completed`].
    pub length: f64,
    /// The direction that points back along the path at its end, in the frame of `end` (see
    /// [`SurfacePoint`]). When a vertex is inserted at the end of a traced path this is the
    /// signpost of its halfedge back to where the path came from (paper Sec. 3.2.3). For a
    /// path of length zero (a direction that leaves the surface at a boundary vertex) it is
    /// the direction that was asked for.
    pub back_direction: f64,
    /// The face in which the path leaves its start, or `None` if it runs along an edge.
    pub first_face: Option<u32>,
    /// The edge the path runs along, if it does.
    pub along_edge: Option<u32>,
    /// Number of edges crossed.
    pub n_crossings: usize,
    /// The crossings in order, if [`TraceOptions::record_crossings`] was set; empty otherwise.
    pub crossings: Vec<Crossing>,
    /// The largest amount by which a value was moved into its valid range: a barycentric
    /// coordinate of the end point that came out negative, or an angle inside the corner of a
    /// vertex that came out beyond the corner, in the units of that value. Zero in most
    /// queries and of the order of the rounding error otherwise.
    pub clamp: f64,
}

/// A point the path is expected to end at, and how far from it the end may be.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Target {
    pub(crate) point: SurfacePoint,
    /// Absolute tolerance on the distance between the target and the end of the path.
    pub(crate) tolerance: f64,
}

/// A trace and how it relates to its target.
#[derive(Debug, Clone)]
pub(crate) struct Outcome {
    pub(crate) trace: Trace,
    /// The distance between the target and the end of the path, if they met within the
    /// tolerance. Then `trace.end` is the target.
    pub(crate) mismatch: Option<f64>,
    /// The smallest such distance over all faces of the path that contain the target.
    pub(crate) closest: Option<f64>,
}

/// Reduces an angle to `[0, 2π)`.
#[inline]
fn wrap(angle: f64) -> f64 {
    let a = angle % TAU;
    let a = if a < 0.0 { a + TAU } else { a };
    if a >= TAU { 0.0 } else { a }
}

#[inline]
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

/// Where in its first face a path starts, for recognizing a path along an edge of that face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Origin {
    /// At the corner in this slot.
    Corner(usize),
    /// On the side in this slot.
    Side(usize),
    /// Inside the face, or in a face the path has entered.
    Elsewhere,
}

/// How the path leaves the face it is in.
#[derive(Debug, Clone, Copy)]
enum Exit {
    /// Through the halfedge in this slot of the face.
    Edge(usize),
    /// Through the vertex in this slot, which lies exactly on the path.
    Vertex(usize),
}

struct Walker<'a, S: Surface> {
    s: &'a S,
    start: SurfacePoint,
    direction: f64,
    length: f64,
    /// The path is followed while crossings lie before this abscissa: `length`, plus the
    /// tolerance of the target if there is one.
    reach: f64,
    record: bool,
    max_crossings: usize,
    target: Option<Target>,
    crossings: Vec<Crossing>,
    n_crossings: usize,
    clamp: f64,
    closest: Option<f64>,
    first_face: Option<u32>,
}

/// A tracing query without a target, with the public options.
pub(crate) fn trace_query<S: Surface>(
    s: &S,
    start: SurfacePoint,
    direction: f64,
    length: f64,
    options: &TraceOptions,
) -> Result<Trace, IntrinsicError> {
    run(s, start, direction, length, options, None).map(|outcome| outcome.trace)
}

/// The tracing kernel: follows the path from `start` in `direction` for `length` over `s`.
/// With a target, the path ends as soon as the target lies in the current face within its
/// tolerance of the path's end, and is otherwise followed until the tolerance is used up.
pub(crate) fn run<S: Surface>(
    s: &S,
    start: SurfacePoint,
    direction: f64,
    length: f64,
    options: &TraceOptions,
    target: Option<Target>,
) -> Result<Outcome, IntrinsicError> {
    if let Some(reason) = point_defect(&start, s.n_vertices(), s.n_edges(), s.n_faces()) {
        return Err(IntrinsicError::InvalidPoint {
            point: start,
            reason,
        });
    }
    if !(direction.is_finite() && (0.0..TAU).contains(&direction)) {
        return Err(IntrinsicError::InvalidDirection { direction });
    }
    if !(length.is_finite() && length > 0.0) {
        return Err(IntrinsicError::InvalidLength { length });
    }
    let mut walker = Walker {
        s,
        start,
        direction,
        length,
        reach: length + target.map_or(0.0, |t| t.tolerance),
        record: options.record_crossings,
        max_crossings: options
            .max_crossings
            .unwrap_or_else(|| s.n_faces().saturating_mul(100).saturating_add(100)),
        target,
        crossings: Vec::new(),
        n_crossings: 0,
        clamp: 0.0,
        closest: None,
        first_face: None,
    };
    match start {
        SurfacePoint::Vertex(v) => walker.start_at_vertex(v as usize),
        SurfacePoint::Edge { edge, t } => walker.start_on_edge(edge as usize, t),
        SurfacePoint::Face { face, bary } => walker.start_in_face(face as usize, bary),
    }
}

impl<S: Surface> Walker<'_, S> {
    fn degenerate(&self, face: usize) -> IntrinsicError {
        IntrinsicError::DegenerateFace {
            surface: S::NAME,
            face: face as u32,
            lengths: face_lengths(self.s, face),
        }
    }

    fn note_clamp(&mut self, amount: f64) {
        if amount > self.clamp {
            self.clamp = amount;
        }
    }

    fn outcome(
        &mut self,
        end: SurfacePoint,
        stop: TraceStop,
        length: f64,
        back_direction: f64,
        along_edge: Option<u32>,
        mismatch: Option<f64>,
    ) -> Outcome {
        Outcome {
            trace: Trace {
                end,
                stop,
                length,
                back_direction,
                first_face: self.first_face,
                along_edge,
                n_crossings: self.n_crossings,
                crossings: std::mem::take(&mut self.crossings),
                clamp: self.clamp,
            },
            mismatch,
            closest: self.closest,
        }
    }

    // ----- starts -----

    /// Start at a vertex in a normalized direction (paper Alg. 2).
    fn start_at_vertex(&mut self, v: usize) -> Result<Outcome, IntrinsicError> {
        let s = self.s;
        let boundary = is_boundary_vertex(s, v);
        if boundary && self.direction > PI {
            let (end, back) = (self.start, self.direction);
            return Ok(self.outcome(end, TraceStop::Boundary, 0.0, back, None, None));
        }
        // The outgoing halfedge with the largest signpost not above the direction opens the
        // wedge (the corner of a face) that contains the direction.
        let h0 = s.outgoing(v) as usize;
        let mut below: Option<(usize, f64)> = None;
        let mut last = (h0, s.signpost(h0));
        let mut h = h0;
        loop {
            let phi = s.signpost(h);
            if phi <= self.direction && below.is_none_or(|(_, best)| phi >= best) {
                below = Some((h, phi));
            }
            if phi >= last.1 {
                last = (h, phi);
            }
            let t = s.twin(prev(h));
            if t == NONE || t as usize == h0 {
                break;
            }
            h = t as usize;
        }
        // At a boundary vertex the direction π is the boundary edge that arrives (SPEC N2). It
        // has no halfedge leaving the vertex, so no signpost equals it; it is the far side of
        // the last wedge. Without this the angle inside that wedge comes out a rounding error
        // short of the corner in half the cases, and the path runs beside the edge.
        if boundary && self.direction == PI {
            return Ok(self.along(prev(last.0), false, 1.0));
        }
        // No signpost at or below the direction: the wedge is the one that contains angle
        // zero, opened by the halfedge with the largest signpost.
        let (h, base) = below.unwrap_or((last.0, last.1 - TAU));
        let (face, k) = (h / 3, h % 3);
        let l = face_lengths(s, face);
        let corner = corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3])
            .ok_or_else(|| self.degenerate(face))?;
        // Back from normalized to actual angles inside the corner.
        let mut delta = (self.direction - base) / angle_scale(s.angle_sum(v), boundary);
        if delta > corner {
            self.note_clamp(delta - corner);
            delta = corner;
        }
        let along_far_side = prev(h);
        if delta <= 0.0 {
            return Ok(self.along(h, true, 0.0));
        }
        if delta >= corner {
            return Ok(self.along(along_far_side, false, 1.0));
        }
        let mut q = [[0.0; 2]; 3];
        q[(k + 1) % 3] = [l[k] * delta.cos(), -l[k] * delta.sin()];
        let apex = unfold_apex(q[k], q[(k + 1) % 3], l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .ok_or_else(|| self.degenerate(face))?;
        q[(k + 2) % 3] = apex;
        self.first_face = Some(face as u32);
        if let Some(done) = self.meet_target(face, &q, Origin::Corner(k))? {
            return Ok(done);
        }
        // Directions within rounding of a side of the corner run along that side.
        if q[(k + 1) % 3][1] >= 0.0 {
            self.first_face = None;
            return Ok(self.along(h, true, 0.0));
        }
        if apex[1] <= 0.0 {
            self.note_clamp(-apex[1] / l[(k + 2) % 3]);
            self.first_face = None;
            return Ok(self.along(along_far_side, false, 1.0));
        }
        self.walk(face, q, Exit::Edge((k + 1) % 3))
    }

    /// Start at a point of an edge, in a direction relative to the edge.
    fn start_on_edge(&mut self, e: usize, t: f64) -> Result<Outcome, IntrinsicError> {
        let s = self.s;
        if t == 0.0 || t == 1.0 {
            return Err(IntrinsicError::InvalidPoint {
                point: self.start,
                reason: "a tracing query cannot start at the end of an edge; start at the vertex",
            });
        }
        let h0 = s.edge_halfedge(e) as usize;
        let (forward, backward) = if s.is_forward(h0) {
            (h0 as u32, s.twin(h0))
        } else {
            (s.twin(h0), h0 as u32)
        };
        let psi = self.direction;
        if psi == 0.0 || psi == PI {
            // Along the edge, in its direction or against it.
            let with_edge = psi == 0.0;
            return Ok(if forward != NONE {
                self.along(forward as usize, with_edge, t)
            } else {
                self.along(backward as usize, !with_edge, 1.0 - t)
            });
        }
        // The face on the side the direction points to, entered through its halfedge `h` at
        // parameter `th`, the direction making the angle `beta ∈ (0, π)` with `h`.
        let (h, th, beta) = if psi < PI {
            (forward, t, psi)
        } else {
            (backward, 1.0 - t, psi - PI)
        };
        if h == NONE {
            let (end, back) = (self.start, self.direction);
            return Ok(self.outcome(end, TraceStop::Boundary, 0.0, back, None, None));
        }
        let h = h as usize;
        let (face, k) = (h / 3, h % 3);
        let l = face_lengths(s, face);
        // In the frame of the path the halfedge points in the direction −beta.
        let u = [beta.cos(), -beta.sin()];
        let mut q = [[0.0; 2]; 3];
        q[k] = [-th * l[k] * u[0], -th * l[k] * u[1]];
        q[(k + 1) % 3] = [(1.0 - th) * l[k] * u[0], (1.0 - th) * l[k] * u[1]];
        let apex = unfold_apex(q[k], q[(k + 1) % 3], l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .ok_or_else(|| self.degenerate(face))?;
        q[(k + 2) % 3] = apex;
        self.first_face = Some(face as u32);
        if let Some(done) = self.meet_target(face, &q, Origin::Side(k))? {
            return Ok(done);
        }
        let exit = Self::exit_after_entry(k, apex);
        self.walk(face, q, exit)
    }

    /// Start at a point of a face, in a direction relative to the face's first halfedge.
    fn start_in_face(&mut self, face: usize, bary: [f64; 3]) -> Result<Outcome, IntrinsicError> {
        if bary.contains(&1.0) {
            return Err(IntrinsicError::InvalidPoint {
                point: self.start,
                reason: "a tracing query cannot start at the corner of a face; start at the vertex",
            });
        }
        let l = face_lengths(self.s, face);
        // The face in its own frame: corner 0 at the origin, corner 1 on the x-axis.
        let v1 = [l[0], 0.0];
        let v2 =
            unfold_apex([0.0, 0.0], v1, l[0], l[1], l[2]).ok_or_else(|| self.degenerate(face))?;
        let p = [
            bary[1] * v1[0] + bary[2] * v2[0],
            bary[1] * v1[1] + bary[2] * v2[1],
        ];
        // Rotate by −direction about the start point.
        let (sin, cos) = self.direction.sin_cos();
        let q = [[0.0, 0.0], v1, v2].map(|v| {
            let d = sub(v, p);
            [d[0] * cos + d[1] * sin, -d[0] * sin + d[1] * cos]
        });
        self.first_face = Some(face as u32);
        if let Some(done) = self.meet_target(face, &q, Origin::Elsewhere)? {
            return Ok(done);
        }
        // The path leaves through the edge that crosses the axis upwards.
        for k in 0..3 {
            if q[k][1] < 0.0 && q[(k + 1) % 3][1] > 0.0 {
                return self.walk(face, q, Exit::Edge(k));
            }
        }
        // A corner lies exactly on the line of the path.
        for k in 0..3 {
            if q[k][1] != 0.0 {
                continue;
            }
            let (a, b) = (q[(k + 1) % 3], q[(k + 2) % 3]);
            if q[k][0] > 0.0 && a[1] * b[1] < 0.0 {
                return self.walk(face, q, Exit::Vertex(k));
            }
        }
        Err(IntrinsicError::DirectionLeavesFace {
            point: self.start,
            direction: self.direction,
        })
    }

    /// The exit from a face entered through its halfedge in slot `k`: the tail of that
    /// halfedge is above the path, its head below, and the third vertex decides.
    #[inline]
    fn exit_after_entry(k: usize, apex: [f64; 2]) -> Exit {
        if apex[1] > 0.0 {
            Exit::Edge((k + 1) % 3)
        } else if apex[1] < 0.0 {
            Exit::Edge((k + 2) % 3)
        } else {
            Exit::Vertex((k + 2) % 3)
        }
    }

    // ----- the walk -----

    /// Follows the path from face to face. `q` holds the unfolded corners of `face`.
    fn walk(
        &mut self,
        mut face: usize,
        mut q: [[f64; 2]; 3],
        mut exit: Exit,
    ) -> Result<Outcome, IntrinsicError> {
        let s = self.s;
        loop {
            let slot = match exit {
                Exit::Edge(slot) => slot,
                Exit::Vertex(slot) => {
                    let x = q[slot][0];
                    if x > self.length {
                        return Ok(self.end_in_face(face, &q));
                    }
                    let stop = if x == self.length {
                        TraceStop::Completed
                    } else {
                        TraceStop::Vertex
                    };
                    let back = self.back_direction_at_corner(face, slot, &q)?;
                    let end = SurfacePoint::Vertex(s.tail(3 * face + slot));
                    return Ok(self.outcome(end, stop, x, back, None, None));
                }
            };
            let (tail, head) = (q[slot], q[(slot + 1) % 3]);
            let (t, x) = axis_crossing(tail, head);
            if x >= self.reach {
                return Ok(self.end_in_face(face, &q));
            }
            let h_out = 3 * face + slot;
            let edge = s.edge(h_out);
            let t_edge = edge_parameter(s, h_out, t);
            let h_in = s.twin(h_out);
            if h_in == NONE {
                // The direction of travel makes the angle atan2(−u_y, u_x) with the halfedge.
                let u = sub(head, tail);
                let along = if s.is_forward(h_out) { 0.0 } else { PI };
                let back = wrap((-u[1]).atan2(u[0]) + along + PI);
                let end = SurfacePoint::Edge { edge, t: t_edge };
                return Ok(self.outcome(end, TraceStop::Boundary, x, back, None, None));
            }
            if self.n_crossings >= self.max_crossings {
                return Err(IntrinsicError::CrossingLimit {
                    start: self.start,
                    limit: self.max_crossings,
                    travelled: x,
                });
            }
            self.n_crossings += 1;
            if self.record {
                self.crossings.push(Crossing {
                    halfedge: h_out as u32,
                    edge,
                    t: t_edge,
                    distance: x,
                });
            }
            let h_in = h_in as usize;
            let k = h_in % 3;
            face = h_in / 3;
            let l = face_lengths(s, face);
            let apex = unfold_apex(head, tail, l[k], l[(k + 1) % 3], l[(k + 2) % 3])
                .ok_or_else(|| self.degenerate(face))?;
            q[k] = head;
            q[(k + 1) % 3] = tail;
            q[(k + 2) % 3] = apex;
            if let Some(done) = self.meet_target(face, &q, Origin::Elsewhere)? {
                return Ok(done);
            }
            exit = Self::exit_after_entry(k, apex);
        }
    }

    /// The path ends inside `face`: its end point in barycentric coordinates.
    fn end_in_face(&mut self, face: usize, q: &[[f64; 2]; 3]) -> Outcome {
        let mut bary = barycentric(q, [self.length, 0.0]);
        let mut moved = 0.0f64;
        for b in &mut bary {
            if *b < 0.0 {
                moved = moved.max(-*b);
                *b = 0.0;
            }
        }
        if moved > 0.0 {
            self.note_clamp(moved);
            let sum = bary[0] + bary[1] + bary[2];
            bary = bary.map(|b| b / sum);
        }
        let back = Self::back_direction_in_face(q);
        let end = SurfacePoint::Face {
            face: face as u32,
            bary,
        };
        let length = self.length;
        self.outcome(end, TraceStop::Completed, length, back, None, None)
    }

    /// The direction back along the path in the frame of the face (its first halfedge).
    fn back_direction_in_face(q: &[[f64; 2]; 3]) -> f64 {
        let u = sub(q[1], q[0]);
        wrap((-u[1]).atan2(u[0]) + PI)
    }

    /// The direction back along the path at the vertex in `slot` of `face`, as a normalized
    /// angle in the frame of that vertex: the signpost of the halfedge leaving the corner plus
    /// the normalized angle from that halfedge to the path.
    fn back_direction_at_corner(
        &mut self,
        face: usize,
        slot: usize,
        q: &[[f64; 2]; 3],
    ) -> Result<f64, IntrinsicError> {
        let s = self.s;
        let h = 3 * face + slot;
        let v = s.tail(h) as usize;
        let boundary = is_boundary_vertex(s, v);
        let l = face_lengths(s, face);
        let corner = corner_angle(l[(slot + 1) % 3], l[slot], l[(slot + 2) % 3])
            .ok_or_else(|| self.degenerate(face))?;
        // Angle from the halfedge leaving the corner to the direction (−1, 0).
        let u = sub(q[(slot + 1) % 3], q[slot]);
        let mut omega = u[1].atan2(-u[0]);
        if omega < 0.0 {
            self.note_clamp(-omega);
            omega = 0.0;
        } else if omega > corner {
            self.note_clamp(omega - corner);
            omega = corner;
        }
        let phi = s.signpost(h) + angle_scale(s.angle_sum(v), boundary) * omega;
        Ok(if boundary { phi.min(PI) } else { wrap(phi) })
    }

    /// If the target lies in `face` within its tolerance of the end of the path, the path
    /// ends there. `origin` says where the path starts if this is its first face: a corner
    /// it starts from is not a place to arrive at, and a target on an edge of the face that
    /// the start lies on means that the path runs along that edge.
    fn meet_target(
        &mut self,
        face: usize,
        q: &[[f64; 2]; 3],
        origin: Origin,
    ) -> Result<Option<Outcome>, IntrinsicError> {
        let s = self.s;
        let Some(target) = self.target else {
            return Ok(None);
        };
        let end = [self.length, 0.0];
        let distance = |p: [f64; 2]| (p[0] - end[0]).hypot(p[1] - end[1]);
        let mut best: Option<(f64, usize)> = None;
        let mut consider = |d: f64, slot: usize| {
            if best.is_none_or(|(b, _)| d < b) {
                best = Some((d, slot));
            }
        };
        match target.point {
            SurfacePoint::Vertex(v) => {
                for (k, &corner) in q.iter().enumerate() {
                    if s.tail(3 * face + k) == v && origin != Origin::Corner(k) {
                        consider(distance(corner), k);
                    }
                }
            }
            SurfacePoint::Edge { edge, t } => {
                for k in 0..3 {
                    if s.edge(3 * face + k) == edge {
                        let th = edge_parameter(s, 3 * face + k, t);
                        let (a, b) = (q[k], q[(k + 1) % 3]);
                        consider(
                            distance([a[0] + th * (b[0] - a[0]), a[1] + th * (b[1] - a[1])]),
                            k,
                        );
                    }
                }
            }
            SurfacePoint::Face { face: f, bary } => {
                if f as usize == face {
                    consider(
                        distance([
                            bary[0] * q[0][0] + bary[1] * q[1][0] + bary[2] * q[2][0],
                            bary[0] * q[0][1] + bary[1] * q[1][1] + bary[2] * q[2][1],
                        ]),
                        0,
                    );
                }
            }
        }
        let Some((d, slot)) = best else {
            return Ok(None);
        };
        if self.closest.is_none_or(|c| d < c) {
            self.closest = Some(d);
        }
        if d > target.tolerance {
            return Ok(None);
        }
        let back = match target.point {
            SurfacePoint::Vertex(_) => self.back_direction_at_corner(face, slot, q)?,
            SurfacePoint::Edge { .. } => {
                let u = sub(q[(slot + 1) % 3], q[slot]);
                let along = if s.is_forward(3 * face + slot) {
                    0.0
                } else {
                    PI
                };
                wrap((-u[1]).atan2(u[0]) + along + PI)
            }
            SurfacePoint::Face { .. } => Self::back_direction_in_face(q),
        };
        // A target on an edge of the face that the start lies on too: the path runs along
        // that edge. `slot` is the corner of a vertex target, the side of an edge target.
        let on_vertex = matches!(target.point, SurfacePoint::Vertex(_));
        let on_edge = matches!(target.point, SurfacePoint::Edge { .. });
        let along_edge = match origin {
            // The sides at corner o are o (leaving) and o + 2 (arriving).
            Origin::Corner(o) if on_vertex => Some(if slot == (o + 1) % 3 { o } else { slot }),
            Origin::Corner(o) if on_edge && (slot == o || slot == (o + 2) % 3) => Some(slot),
            // The ends of side k are the corners k and k + 1.
            Origin::Side(k) if on_vertex && (slot == k || slot == (k + 1) % 3) => Some(k),
            Origin::Side(k) if on_edge && slot == k => Some(k),
            _ => None,
        }
        .map(|side| s.edge(3 * face + side));
        if along_edge.is_some() {
            self.first_face = None;
        }
        let length = self.length;
        Ok(Some(self.outcome(
            target.point,
            TraceStop::Completed,
            length,
            back,
            along_edge,
            Some(d),
        )))
    }

    // ----- paths along an edge -----

    /// The path runs along halfedge `h` from its parameter `t0`, in its direction (`forward`)
    /// or against it.
    fn along(&mut self, h: usize, forward: bool, t0: f64) -> Outcome {
        let s = self.s;
        let edge = s.edge(h);
        let len = s.length(edge as usize);
        // Distance to the vertex the path is heading for.
        let to_vertex = if forward { (1.0 - t0) * len } else { t0 * len };
        let far = if forward { s.tail(next(h)) } else { s.tail(h) };
        // The direction back at that vertex: the signpost of the halfedge that leaves it along
        // this edge. A boundary halfedge has no twin; its head sees it at the angle π.
        let back_at_far = if forward {
            match s.twin(h) {
                NONE => PI,
                t => s.signpost(t as usize),
            }
        } else {
            s.signpost(h)
        };
        // In the frame of the edge the path points along it (0) or against it (π).
        let with_edge = forward == s.is_forward(h);
        let back_on_edge = if with_edge { PI } else { 0.0 };
        let length = self.length;
        if let Some(target) = self.target {
            let at = match target.point {
                SurfacePoint::Vertex(v) if v == far => Some(to_vertex),
                SurfacePoint::Edge { edge: e, t } if e == edge => {
                    let th = edge_parameter(s, h, t);
                    let d = if forward {
                        (th - t0) * len
                    } else {
                        (t0 - th) * len
                    };
                    (d >= 0.0).then_some(d)
                }
                _ => None,
            };
            if let Some(at) = at {
                let d = (at - length).abs();
                if self.closest.is_none_or(|c| d < c) {
                    self.closest = Some(d);
                }
                if d <= target.tolerance {
                    let back = match target.point {
                        SurfacePoint::Vertex(_) => back_at_far,
                        _ => back_on_edge,
                    };
                    return self.outcome(
                        target.point,
                        TraceStop::Completed,
                        length,
                        back,
                        Some(edge),
                        Some(d),
                    );
                }
            }
        }
        if length < to_vertex {
            let th = if forward {
                t0 + length / len
            } else {
                t0 - length / len
            };
            let end = SurfacePoint::Edge {
                edge,
                t: edge_parameter(s, h, th),
            };
            return self.outcome(
                end,
                TraceStop::Completed,
                length,
                back_on_edge,
                Some(edge),
                None,
            );
        }
        let stop = if length == to_vertex {
            TraceStop::Completed
        } else {
            TraceStop::Vertex
        };
        let end = SurfacePoint::Vertex(far);
        self.outcome(end, stop, to_vertex, back_at_far, Some(edge), None)
    }
}

#[cfg(test)]
mod tests {
    //! Targets that are not vertices. Vertices on input edges and inside input faces do not
    //! exist before vertex insertion does; the kernel is written for them and tested here.

    use ddg_core::Geometry;
    use ddg_core::generators::grid_disk;

    use super::*;
    use crate::InputSurface;

    /// A flat 4 × 4 grid, sheared a little so that nothing is aligned, with its positions.
    fn surface() -> (InputSurface, Vec<[f64; 2]>) {
        let mut shape = grid_disk(4, 4).unwrap();
        for p in &mut shape.positions {
            p[0] += 0.3 * p[1];
        }
        let mesh = shape.topology().unwrap();
        let geometry = Geometry::from_positions(&mesh, &shape.positions).unwrap();
        let xy = shape.positions.iter().map(|p| [p[0], p[1]]).collect();
        (InputSurface::new(&mesh, &geometry).unwrap(), xy)
    }

    /// The position of a point of an edge or a face, and the angle in the plane of direction
    /// zero of its frame.
    fn place(s: &InputSurface, xy: &[[f64; 2]], p: SurfacePoint) -> ([f64; 2], f64) {
        let towards = |a: [f64; 2], b: [f64; 2]| (b[1] - a[1]).atan2(b[0] - a[0]);
        match p {
            SurfacePoint::Edge { edge, t } => {
                let [a, b] = s.mesh().edge_vertices(edge as usize).map(|v| xy[v]);
                (
                    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])],
                    towards(a, b),
                )
            }
            SurfacePoint::Face { face, bary } => {
                let q = s.mesh().faces()[face as usize].map(|v| xy[v as usize]);
                (
                    [0, 1].map(|k| bary[0] * q[0][k] + bary[1] * q[1][k] + bary[2] * q[2][k]),
                    towards(q[0], q[1]),
                )
            }
            SurfacePoint::Vertex(_) => unreachable!("not used here"),
        }
    }

    fn interior_edge(s: &InputSurface, n: usize) -> u32 {
        (0..s.mesh().n_edges())
            .filter(|&e| !s.mesh().is_boundary_edge(e))
            .nth(n)
            .unwrap() as u32
    }

    #[test]
    fn a_path_ends_at_a_target_on_an_edge_or_in_a_face() {
        let (s, xy) = surface();
        let options = TraceOptions {
            record_crossings: true,
            max_crossings: None,
        };
        let points = [
            SurfacePoint::Face {
                face: 2,
                bary: [0.2, 0.3, 0.5],
            },
            SurfacePoint::Edge {
                edge: interior_edge(&s, 3),
                t: 0.35,
            },
            SurfacePoint::Face {
                face: 27,
                bary: [0.6, 0.1, 0.3],
            },
            SurfacePoint::Edge {
                edge: interior_edge(&s, 30),
                t: 0.8,
            },
            SurfacePoint::Face {
                face: 13,
                bary: [0.1, 0.45, 0.45],
            },
        ];
        let mut with_crossings = 0;
        for &start in &points {
            for &end in &points {
                if start == end {
                    continue;
                }
                let ((a, zero), (b, zero_at_end)) = (place(&s, &xy, start), place(&s, &xy, end));
                let length = (b[0] - a[0]).hypot(b[1] - a[1]);
                let world = (b[1] - a[1]).atan2(b[0] - a[0]);
                let direction = wrap(world - zero);
                let target = Target {
                    point: end,
                    tolerance: 1e-12,
                };
                let outcome = run(&s, start, direction, length, &options, Some(target)).unwrap();
                let context = format!("{start} to {end}: {outcome:?}");
                let mismatch = outcome.mismatch.unwrap_or_else(|| panic!("{context}"));
                assert!(mismatch < 1e-14, "{context}");
                assert_eq!(outcome.closest, Some(mismatch), "{context}");
                let trace = &outcome.trace;
                assert_eq!(trace.end, end, "{context}");
                assert_eq!(
                    (trace.stop, trace.length),
                    (TraceStop::Completed, length),
                    "{context}"
                );
                // The way back, in the frame of the edge or the face of the target.
                let back = wrap(world + PI - zero_at_end);
                let apart = (trace.back_direction - back).abs();
                assert!(apart.min(TAU - apart) < 1e-12, "{context}: expected {back}");
                // The same path without a target crosses the same edges, and one more if it
                // ends on an edge and rounding puts its end just beyond.
                let free = run(&s, start, direction, length, &options, None)
                    .unwrap()
                    .trace;
                let n = trace.n_crossings;
                assert!(
                    free.n_crossings == n || free.n_crossings == n + 1,
                    "{context}"
                );
                assert_eq!(trace.crossings[..], free.crossings[..n], "{context}");
                if free.n_crossings > n {
                    let SurfacePoint::Edge { edge, .. } = end else {
                        panic!("{context}");
                    };
                    assert_eq!(free.crossings[n].edge, edge, "{context}");
                }
                with_crossings += usize::from(n > 0);

                // A target that is somewhere else on the way is not met, and the path says by
                // how much: 1e-3 too short.
                let short = run(&s, start, direction, length - 1e-3, &options, Some(target));
                let short = short.unwrap();
                assert_eq!(short.mismatch, None, "{context}");
                let closest = short.closest.unwrap();
                assert!((closest - 1e-3).abs() < 1e-12, "{context}: {closest}");
            }
        }
        assert!(with_crossings >= 15);
    }

    #[test]
    fn a_path_along_an_edge_ends_at_a_target_on_that_edge() {
        let (s, _) = surface();
        let options = TraceOptions::default();
        let mesh = s.mesh();
        for v in [6usize, 12, 18] {
            for h in mesh.outgoing_halfedges(v) {
                let e = mesh.edge(h);
                let len = s.geometry().edge_lengths()[e];
                let forward = mesh.edge_vertices(e)[0] == v;
                // 0.6 of the way from v, as a parameter of the edge.
                let t = if forward { 0.6 } else { 0.4 };
                let target = Target {
                    point: SurfacePoint::Edge { edge: e as u32, t },
                    tolerance: 1e-12 * len,
                };
                let start = SurfacePoint::Vertex(v as u32);
                let hit = run(&s, start, s.signpost(h), 0.6 * len, &options, Some(target)).unwrap();
                assert!(hit.mismatch.unwrap() < 1e-15, "{hit:?}");
                assert_eq!(hit.trace.end, target.point);
                assert_eq!(hit.trace.along_edge, Some(e as u32));
                assert_eq!(hit.trace.back_direction, if forward { PI } else { 0.0 });
                // Half the way only: 0.1 of the edge short of the target.
                let miss =
                    run(&s, start, s.signpost(h), 0.5 * len, &options, Some(target)).unwrap();
                assert_eq!(miss.mismatch, None);
                assert!((miss.closest.unwrap() - 0.1 * len).abs() < 1e-15);
                assert!(matches!(miss.trace.end, SurfacePoint::Edge { .. }));
            }
        }
    }
}
