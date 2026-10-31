//! Reading the integers, geometry: the positions of the crossings (SPEC §2.3, §3.1, §3.2; work
//! package WP2).
//!
//! The connectivity of every curve is decided by the integers ([`crate::trace`]); floating
//! point decides only where along the edges the crossings lie. Both directions unfold a strip of
//! triangles into the plane and intersect one straight segment with its interior edges:
//!
//! * **An input edge over the triangulation** (Supp. Alg. 2 lines 6–17, Fig. 8):
//!   [`IntegerTriangulation::extract_geometric_curve`] (`ExtractCurve`),
//!   [`IntegerTriangulation::extract_geometric_crossing`] (`ExtractGeometricCrossing`),
//!   [`IntegerTriangulation::segment_geometry`] (the strip of one [`CurveSegment`]) and
//!   [`IntegerTriangulation::extract_edge_geometry`] (all pieces of an input edge, Alg. 4). The
//!   strip is made of faces of the triangulation, laid out with its lengths; the parameter
//!   along a piece of a split input edge is mapped to the whole input edge with the places of
//!   the piece's ends (`uStart`, `uEnd`, Alg. 2 lines 10–16; SPEC I17).
//! * **An edge of the triangulation over the input mesh** (`ExtractIntrinsicEdge`, SPEC §3.1,
//!   FlipOut contract R2): [`IntegerTriangulation::extract_intrinsic_edge`]. The input edge of
//!   every crossing is identified from the integers (Alg. 1 and 3), the strip is made of input
//!   faces, laid out with the input lengths.
//!
//! The strip intersection (formula id `ic-strip-intersection`, pending in `CITED.md`) evaluates
//! two formulas of the signpost crate: every face of the strip is unfolded across the edge it
//! is entered by (`it-unfold-apex`), and in the frame in which the straight segment is the
//! positive x-axis, each interior edge is cut where its end points change sides
//! (`it-ray-edge-intersection`). With exact lengths the segment lies inside the strip (the strip
//! is the face sequence of the geodesic), so every parameter lies in `[0, 1]`. **SPEC I5:** a
//! parameter that falls outside by rounding is clamped and counted ([`StripReport`]), never
//! used to change connectivity.
//!
//! **Lengths.** The straight segment's length in the unfolded strip is compared with the length
//! it should have (the intrinsic edge's length, or the piece of the input edge between its
//! ends); the difference is reported ([`StripReport::max_relative_mismatch`],
//! [`EdgeTrace::mismatch`]). It is not a failure: the connectivity is exact whatever the
//! rounding of the positions.

use std::collections::HashMap;

use ddg_intrinsic::formulas::{axis_crossing, triangle_layout, unfold_apex};
use ddg_intrinsic::{
    Crossing as InputCrossing, EdgeTrace, EdgeTraceOptions, IntrinsicError, SurfacePoint,
};

use crate::certify::MAX_CERTIFIED_CROSSINGS;
use crate::coords::IntegerTriangulation;
use crate::error::{IntegerError, Place, ValidityItem};
use crate::trace::{Crossing, CurveSegment, InputEdgePath, PathPiece, TracedInputEdge, trace_with};

/// A geometric crossing `z = (ij→, p, ab, u, v)` (paper Sec. 2.5): a combinatorial crossing of
/// a halfedge of the triangulation with an input edge, and where it lies on both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometricCrossing {
    /// The combinatorial crossing `(ij→, p)`, oriented in the direction in which the curve was
    /// followed (the curve leaves the face of `ij→` through it).
    pub crossing: Crossing,
    /// The input edge `ab`.
    pub input_edge: u32,
    /// The parameter of the point along the input edge **in the direction of the edge** (from
    /// its lower-numbered to its higher-numbered vertex, as [`SurfacePoint::Edge`]). The paper
    /// measures `u` along the halfedge in the direction of the curve; that is `u` or `1 − u`.
    pub u: f64,
    /// The parameter of the point along the halfedge `ij→` of the crossing, from its tail `i`.
    pub v: f64,
    /// The sine of the angle at which the two edges cross, in `(0, 1]` (`0` only for a clamped
    /// crossing of parallel lines). It conditions `u` and `v`: moving either edge by `δ`
    /// normal to itself moves the crossing by `δ / sin` along the other, so positions computed
    /// in different ways agree to their rounding divided by this.
    pub sin_angle: f64,
}

impl GeometricCrossing {
    /// The point on the input mesh: `Edge { input_edge, u }`.
    pub fn on_input(&self) -> SurfacePoint {
        SurfacePoint::Edge {
            edge: self.input_edge,
            t: self.u,
        }
    }

    /// The point on the triangulation `tri`: a point of the edge of `ij→`, with the parameter
    /// in the direction of that edge (its
    /// [`edge_halfedge`](ddg_intrinsic::Complex::edge_halfedge)).
    pub fn on_triangulation(&self, tri: &IntegerTriangulation) -> SurfacePoint {
        let h = self.crossing.halfedge as usize;
        let c = tri.complex();
        SurfacePoint::Edge {
            edge: c.edge(h) as u32,
            t: if c.is_forward(h) {
                self.v
            } else {
                1.0 - self.v
            },
        }
    }
}

/// What the floating point of strip intersections did (SPEC I5).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StripReport {
    /// Number of geometric crossings computed.
    pub n_crossings: usize,
    /// Crossings at which a parameter (along the crossed edge, or along the straight segment)
    /// fell outside `[0, 1]` by rounding and was clamped (I5, `clamped_crossings`). Zero with
    /// exact arithmetic.
    pub clamped_crossings: usize,
    /// The largest amount by which a parameter was moved into `[0, 1]`.
    pub max_clamp: f64,
    /// The largest difference between the length of a straight segment in its unfolded strip
    /// and the length it should have, relative to the latter: the intrinsic edge's length for
    /// [`extract_intrinsic_edge`](IntegerTriangulation::extract_intrinsic_edge), the length of
    /// the piece of the input edge between the segment's ends for the curves of input edges.
    pub max_relative_mismatch: f64,
    /// Crossings visited by the walks that identify input edges (Alg. 1) in
    /// [`extract_intrinsic_edge`](IntegerTriangulation::extract_intrinsic_edge); zero
    /// elsewhere. The cost of an extraction is linear in this plus the number of crossings.
    pub walk_steps: u64,
}

impl StripReport {
    /// Adds the counts of `other` and keeps the larger maxima.
    pub fn merge(&mut self, other: &StripReport) {
        self.n_crossings += other.n_crossings;
        self.clamped_crossings += other.clamped_crossings;
        self.max_clamp = self.max_clamp.max(other.max_clamp);
        self.max_relative_mismatch = self.max_relative_mismatch.max(other.max_relative_mismatch);
        self.walk_steps += other.walk_steps;
    }

    fn note_hit(&mut self, hit: &Hit) {
        self.n_crossings += 1;
        if hit.clamp > 0.0 {
            self.clamped_crossings += 1;
            self.max_clamp = self.max_clamp.max(hit.clamp);
        }
    }

    fn note_mismatch(&mut self, relative: f64) {
        // NaN (a zero expected length) counts as the largest mismatch.
        if relative.is_nan() || relative > self.max_relative_mismatch {
            self.max_relative_mismatch = if relative.is_nan() {
                f64::INFINITY
            } else {
                relative
            };
        }
    }
}

/// The result of [`IntegerTriangulation::extract_geometric_curve`] (Supp. Alg. 2): the curve
/// segment through a crossing, the input edge it belongs to, and its geometric crossings.
#[derive(Debug, Clone, PartialEq)]
pub struct GeometricCurve {
    /// The combinatorial segment (WP1's [`CurveSegment`], unchanged).
    pub segment: CurveSegment,
    /// The input edge it belongs to (Alg. 3).
    pub traced: TracedInputEdge,
    /// The parameter of the segment's start vertex along the input edge, in the direction of
    /// the edge (`uStart`; `0` or `1` at an input vertex, `q⁰` at an inserted one).
    pub u_start: f64,
    /// The same for the segment's end vertex (`uEnd`).
    pub u_end: f64,
    /// The geometric crossings in the order of [`CurveSegment::crossings`].
    pub crossings: Vec<GeometricCrossing>,
    /// Clamps and the length mismatch of the strip.
    pub report: StripReport,
}

/// The result of [`IntegerTriangulation::extract_edge_geometry`] (Supp. Alg. 4 with the
/// geometry of Alg. 2): the path of an input edge with the positions of its crossings.
#[derive(Debug, Clone, PartialEq)]
pub struct InputEdgeGeometry {
    /// The combinatorial path (WP1's [`InputEdgePath`], unchanged).
    pub path: InputEdgePath,
    /// Per vertex of the path: its parameter along the input edge, in the direction of the edge
    /// (`0` or `1` at the input vertices, `q⁰` at the inserted ones).
    pub vertex_u: Vec<f64>,
    /// Per piece of the path: its geometric crossings in order (empty for a piece along an edge
    /// of the triangulation).
    pub crossings: Vec<Vec<GeometricCrossing>>,
    /// Clamps and the largest length mismatch over the pieces.
    pub report: StripReport,
}

// ----- the strip -----

type Corners = [[f64; 2]; 3];

/// Why a strip could not be laid out.
enum StripFault {
    /// A face whose lengths do not form a triangle.
    Degenerate { face: usize, lengths: [f64; 3] },
    /// Exit `i` is not a halfedge of the face the strip has reached.
    NotInFace(usize),
    /// Exit `i` has no twin.
    Boundary(usize),
}

/// A strip of triangles unfolded into the plane.
struct Strip {
    /// The corners of the first face, by slot.
    first: Corners,
    /// The corners of the last face, by slot.
    last: Corners,
    /// The last face.
    last_face: usize,
    /// Per exit halfedge: its tail and head.
    sides: Vec<[[f64; 2]; 2]>,
}

/// Lays out the strip that starts in `first_face` and leaves each face through the next halfedge
/// of `exits` into the face of its twin. The first face is placed by
/// [`triangle_layout`]; every next face is unfolded across the halfedge it is entered by.
///
/// [formula: it-unfold-apex]
fn lay_out_strip(
    face_lengths: &dyn Fn(usize) -> [f64; 3],
    twin: &dyn Fn(usize) -> Option<usize>,
    first_face: usize,
    exits: &[usize],
) -> Result<Strip, StripFault> {
    let l = face_lengths(first_face);
    let first = triangle_layout(l).ok_or(StripFault::Degenerate {
        face: first_face,
        lengths: l,
    })?;
    let (mut q, mut face) = (first, first_face);
    let mut sides = Vec::with_capacity(exits.len());
    for (i, &h) in exits.iter().enumerate() {
        if h / 3 != face {
            return Err(StripFault::NotInFace(i));
        }
        let slot = h % 3;
        let (tail, head) = (q[slot], q[(slot + 1) % 3]);
        sides.push([tail, head]);
        let g = twin(h).ok_or(StripFault::Boundary(i))?;
        face = g / 3;
        let k = g % 3;
        let l = face_lengths(face);
        let apex = unfold_apex(head, tail, l[k], l[(k + 1) % 3], l[(k + 2) % 3])
            .ok_or(StripFault::Degenerate { face, lengths: l })?;
        q[k] = head;
        q[(k + 1) % 3] = tail;
        q[(k + 2) % 3] = apex;
    }
    Ok(Strip {
        first,
        last: q,
        last_face: face,
        sides,
    })
}

/// Where the straight segment meets one interior edge of the strip.
struct Hit {
    /// Parameter along the crossed halfedge from its tail, in `[0, 1]`.
    t: f64,
    /// Parameter along the segment, in `[0, 1]`.
    s: f64,
    /// Distance from the start of the segment (`s` times its length).
    x: f64,
    /// How far a parameter was moved into `[0, 1]` (0 if neither was).
    clamp: f64,
    /// The sine of the angle between the segment and the halfedge.
    sin_angle: f64,
}

/// The straight segment from `p` to `q` in the plane of a strip.
struct Chord {
    p: [f64; 2],
    d: [f64; 2],
    len: f64,
}

impl Chord {
    fn new(p: [f64; 2], q: [f64; 2]) -> Option<Chord> {
        let d = [q[0] - p[0], q[1] - p[1]];
        let len = d[0].hypot(d[1]);
        (len > 0.0 && len.is_finite()).then_some(Chord { p, d, len })
    }

    /// `x` in the frame in which the segment is the positive x-axis from the origin.
    fn frame(&self, x: [f64; 2]) -> [f64; 2] {
        let v = [x[0] - self.p[0], x[1] - self.p[1]];
        [
            (v[0] * self.d[0] + v[1] * self.d[1]) / self.len,
            (self.d[0] * v[1] - self.d[1] * v[0]) / self.len,
        ]
    }

    /// Where the segment crosses the halfedge from `tail` to `head`, which it crosses from its
    /// left to its right (the tail lies below the segment's line, the head above). Parameters
    /// outside `[0, 1]` are clamped and the amount is reported (SPEC I5); a halfedge parallel to
    /// the segment (no finite parameter) is cut at its middle, reported as a clamp of `0.5`.
    ///
    /// [formula: it-ray-edge-intersection]
    fn hit(&self, tail: [f64; 2], head: [f64; 2]) -> Hit {
        let (a, b) = (self.frame(tail), self.frame(head));
        let (raw, _) = axis_crossing(a, b);
        let (t, mut clamp) = if !raw.is_finite() {
            (0.5, 0.5)
        } else if raw < 0.0 {
            (0.0, -raw)
        } else if raw > 1.0 {
            (1.0, raw - 1.0)
        } else {
            (raw, 0.0)
        };
        let x_raw = a[0] + t * (b[0] - a[0]);
        let s_raw = x_raw / self.len;
        let s = s_raw.clamp(0.0, 1.0);
        clamp = clamp.max((s_raw - s).abs());
        let x = if s == s_raw { x_raw } else { s * self.len };
        let sin_angle = ((b[1] - a[1]) / (b[0] - a[0]).hypot(b[1] - a[1])).abs();
        Hit {
            t,
            s,
            x,
            clamp,
            sin_angle: if sin_angle.is_finite() {
                sin_angle
            } else {
                0.0
            },
        }
    }
}

/// The point at parameter `t` (from the tail) of the side from `a` to `b`.
fn lerp(a: [f64; 2], b: [f64; 2], t: f64) -> [f64; 2] {
    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
}

fn invalid_edge(e: usize, detail: String) -> IntegerError {
    IntegerError::invalid(ValidityItem::Extraction, Place::Edge(e as u32), detail)
}

/// Sentinel that stops a walk of [`trace_with`] early; it never leaves this module.
fn stop_walk() -> IntegerError {
    IntegerError::Overflow { edge: u32::MAX }
}

/// How the path of an intrinsic edge leaves its tail, or arrives at its head, on the input mesh.
#[derive(Clone, Copy)]
enum EndPlace {
    /// At input vertex `v`, inside the corner of the input face of the outgoing input halfedge
    /// `out` (the corner at `v` between `out` and the next input halfedge counter-clockwise).
    Corner { out: usize },
    /// At the point `t` (in the edge's direction) of input edge `edge`.
    OnEdge { edge: usize, t: f64 },
    /// Inside the input face of the strip's first (or last) face, at these barycentric
    /// coordinates.
    InFace { bary: [f64; 3] },
}

/// The input faces crossed by an intrinsic edge, as decided by the integers.
struct InputStrip {
    first_face: usize,
    /// The input halfedges through which the path leaves each face.
    exits: Vec<usize>,
    start: EndPlace,
    end: EndPlace,
}

impl IntegerTriangulation {
    // ----- input edges over the triangulation (Alg. 2) -----

    /// The parameter, in the direction of input edge `input_edge`, of vertex `v` of the
    /// triangulation that lies on it: `0` or `1` for its end points, `q⁰` for an inserted
    /// vertex on it.
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] naming the vertex if it does not lie on the edge.
    pub fn parameter_on_input_edge(
        &self,
        v: usize,
        input_edge: usize,
    ) -> Result<f64, IntegerError> {
        let mesh = self.input.mesh();
        let bad = |detail: String| {
            IntegerError::invalid(ValidityItem::Extraction, Place::Vertex(v as u32), detail)
        };
        if input_edge >= mesh.n_edges() || v >= self.complex.n_vertices() {
            return Err(bad(format!(
                "vertex {v} or input edge {input_edge} does not exist"
            )));
        }
        let [x, y] = mesh.edge_vertices(input_edge);
        match self.locations[v] {
            SurfacePoint::Vertex(w) if w as usize == x => Ok(0.0),
            SurfacePoint::Vertex(w) if w as usize == y => Ok(1.0),
            SurfacePoint::Edge { edge, t } if edge as usize == input_edge => Ok(t),
            other => Err(bad(format!(
                "the curve of input edge {input_edge} ends here, but the vertex lies at {other}"
            ))),
        }
    }

    /// The geometric crossings of one curve segment (Supp. Alg. 2 lines 6–17): the strip of
    /// faces of the triangulation that `segment` crosses is laid out with its edge lengths, and
    /// the straight segment from the start vertex to the end vertex is intersected with every
    /// crossed halfedge. `u_start` and `u_end` are the parameters of the segment's ends along
    /// `input_edge` (in the edge's direction); the parameter along the segment is mapped to
    /// `u = (1 − ũ) u_start + ũ u_end` (lines 10–16). The report's mismatch compares the
    /// segment's length in the strip with `|u_end − u_start|` times the input edge's length.
    ///
    /// Strip intersection: formula id `ic-strip-intersection` (pending), evaluated with
    /// `it-unfold-apex` and `it-ray-edge-intersection`.
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] naming the crossing when the segment is not a chain
    /// of faces (a crossing not in the face the previous one entered, a boundary halfedge, an
    /// end corner in another face); [`IntegerError::IndexOutOfRange`];
    /// [`IntegerError::Intrinsic`] (a degenerate face, or a segment of length zero in the
    /// plane).
    pub fn segment_geometry(
        &self,
        segment: &CurveSegment,
        input_edge: usize,
        u_start: f64,
        u_end: f64,
    ) -> Result<(Vec<GeometricCrossing>, StripReport), IntegerError> {
        let Some(&first) = segment.crossings.first() else {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Vertex(segment.start_vertex),
                "a curve segment without crossings has no strip".to_string(),
            ));
        };
        let mesh = self.input.mesh();
        if input_edge >= mesh.n_edges() {
            return Err(IntegerError::IndexOutOfRange {
                what: "input edge",
                index: input_edge,
                bound: mesh.n_edges(),
            });
        }
        let nh = self.complex.n_halfedges();
        for z in &segment.crossings {
            self.check_crossing(*z)?;
        }
        for (what, h) in [
            ("start corner", segment.start as usize),
            ("end corner", segment.end.corner as usize),
        ] {
            if h >= nh {
                return Err(IntegerError::IndexOutOfRange {
                    what,
                    index: h,
                    bound: nh,
                });
            }
        }
        if segment.start as usize / 3 != first.halfedge as usize / 3 {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                first.place(),
                format!(
                    "the start corner {} is not in the face of the first crossing",
                    segment.start
                ),
            ));
        }
        let exits: Vec<usize> = segment
            .crossings
            .iter()
            .map(|z| z.halfedge as usize)
            .collect();
        let strip = lay_out_strip(
            &|f| self.face_lengths(f),
            &|h| self.complex.twin(h),
            first.halfedge as usize / 3,
            &exits,
        )
        .map_err(|fault| match fault {
            StripFault::Degenerate { face, lengths } => {
                IntegerError::Intrinsic(IntrinsicError::DegenerateFace {
                    surface: "intrinsic triangulation",
                    face: face as u32,
                    lengths,
                })
            }
            StripFault::NotInFace(i) => IntegerError::invalid(
                ValidityItem::Extraction,
                segment.crossings[i].place(),
                "the crossing is not in the face the curve entered".to_string(),
            ),
            StripFault::Boundary(i) => IntegerError::invalid(
                ValidityItem::Extraction,
                segment.crossings[i].place(),
                "the curve leaves the surface through a boundary halfedge".to_string(),
            ),
        })?;
        if segment.end.corner as usize / 3 != strip.last_face {
            return Err(IntegerError::invalid(
                ValidityItem::Extraction,
                Place::Halfedge(segment.end.corner),
                format!(
                    "the end corner is not in face {}, which the last crossing enters",
                    strip.last_face
                ),
            ));
        }
        let p = strip.first[segment.start as usize % 3];
        let q = strip.last[segment.end.corner as usize % 3];
        let chord = Chord::new(p, q).ok_or_else(|| {
            IntegerError::Intrinsic(IntrinsicError::DegenerateFace {
                surface: "intrinsic triangulation",
                face: strip.last_face as u32,
                lengths: self.face_lengths(strip.last_face),
            })
        })?;
        let mut report = StripReport::default();
        let expected = (u_end - u_start).abs() * self.input.geometry().edge_lengths()[input_edge];
        report.note_mismatch((chord.len - expected).abs() / expected);
        let mut out = Vec::with_capacity(exits.len());
        for (z, side) in segment.crossings.iter().zip(&strip.sides) {
            let hit = chord.hit(side[0], side[1]);
            report.note_hit(&hit);
            out.push(GeometricCrossing {
                crossing: *z,
                input_edge: input_edge as u32,
                u: (1.0 - hit.s) * u_start + hit.s * u_end,
                v: hit.t,
                sin_angle: hit.sin_angle,
            });
        }
        Ok((out, report))
    }

    /// Supp. Alg. 2, `ExtractCurve(ζ)`, with its geometry: the curve segment through `start`
    /// (both directions, [`extract_curve`](Self::extract_curve)), the input edge it belongs to
    /// ([`identify_traced_edge`](Self::identify_traced_edge), walking back past inserted
    /// vertices), the parameters of its ends along that input edge, and the geometric crossings
    /// ([`segment_geometry`](Self::segment_geometry)).
    ///
    /// # Errors
    /// Those of the steps.
    pub fn extract_geometric_curve(&self, start: Crossing) -> Result<GeometricCurve, IntegerError> {
        let segment = self.extract_curve(start)?;
        let traced = self.identify_traced_edge(&segment)?;
        let input_edge = traced.input_edge as usize;
        let u_start = self.parameter_on_input_edge(segment.start_vertex as usize, input_edge)?;
        let u_end = self.parameter_on_input_edge(segment.end.vertex as usize, input_edge)?;
        let (crossings, report) = self.segment_geometry(&segment, input_edge, u_start, u_end)?;
        Ok(GeometricCurve {
            segment,
            traced,
            u_start,
            u_end,
            crossings,
            report,
        })
    }

    /// `ExtractGeometricCrossing(ζ)` (paper Sec. 3.1):
    /// [`extract_geometric_curve`](Self::extract_geometric_curve) through `z`, returning only the
    /// crossing `z`. Its cost is that of the whole curve segment.
    ///
    /// # Errors
    /// As [`extract_geometric_curve`](Self::extract_geometric_curve).
    pub fn extract_geometric_crossing(
        &self,
        z: Crossing,
    ) -> Result<(GeometricCrossing, StripReport), IntegerError> {
        let curve = self.extract_geometric_curve(z)?;
        let found = curve
            .crossings
            .iter()
            .find(|c| c.crossing == z)
            .copied()
            .ok_or_else(|| {
                IntegerError::invalid(
                    ValidityItem::Extraction,
                    z.place(),
                    "the curve through the crossing does not contain it".to_string(),
                )
            })?;
        Ok((found, curve.report))
    }

    /// Supp. Alg. 4 with the geometry of Alg. 2: the path of input halfedge `input_halfedge`
    /// over the triangulation ([`extract_edge`](Self::extract_edge)) and the positions of all
    /// its crossings, each piece between two vertices of the triangulation on the input edge
    /// laid out in its own strip and mapped to the whole input edge by the places of its ends
    /// (SPEC I17).
    ///
    /// # Errors
    /// Those of [`extract_edge`](Self::extract_edge) and
    /// [`segment_geometry`](Self::segment_geometry).
    pub fn extract_edge_geometry(
        &self,
        input_halfedge: usize,
    ) -> Result<InputEdgeGeometry, IntegerError> {
        let path = self.extract_edge(input_halfedge)?;
        let input_edge = path.input_edge as usize;
        let vertex_u = path
            .vertices
            .iter()
            .map(|&v| self.parameter_on_input_edge(v as usize, input_edge))
            .collect::<Result<Vec<_>, _>>()?;
        let mut report = StripReport::default();
        let mut crossings = Vec::with_capacity(path.pieces.len());
        for (i, piece) in path.pieces.iter().enumerate() {
            match piece {
                PathPiece::Along { .. } => crossings.push(Vec::new()),
                PathPiece::Crossings(segment) => {
                    let (c, r) =
                        self.segment_geometry(segment, input_edge, vertex_u[i], vertex_u[i + 1])?;
                    report.merge(&r);
                    crossings.push(c);
                }
            }
        }
        Ok(InputEdgeGeometry {
            path,
            vertex_u,
            crossings,
            report,
        })
    }

    // ----- edges of the triangulation over the input mesh (ExtractIntrinsicEdge) -----

    /// `ExtractIntrinsicEdge(e)` (SPEC §3.1, FlipOut contract R2): the path of edge `e` over the
    /// input mesh, from the tail of its [`edge_halfedge`](ddg_intrinsic::Complex::edge_halfedge)
    /// to its head, as the signpost representation's
    /// [`trace_edge`](ddg_intrinsic::IntrinsicTriangulation::trace_edge) returns it.
    ///
    /// * `n_e ≥ 0`: for each of the `n_e` crossings of `e` in order from its tail, the input edge
    ///   crossing it is identified exactly from the integers: the curve is followed from the
    ///   crossing to a vertex (Alg. 1) and identified there (Alg. 3, walking back past inserted
    ///   vertices on input edges). The walk stops early at a crossing of `e` that an earlier
    ///   walk identified, and every crossing of `e` that a walk passes is identified with it, so
    ///   each stretch of a curve is walked once ([`StripReport::walk_steps`]). The input faces of
    ///   the path follow exactly: the face in which the edge leaves its tail (from the
    ///   roundabout at an input vertex, from `q⁰` at an inserted one), then across each
    ///   identified input edge in turn; consecutive input edges that do not bound a common input
    ///   face, or a path that does not arrive at the head's place (its roundabout or `q⁰`), are
    ///   an error naming the edge. The positions are those of the straight segment from the
    ///   tail's place to the head's place through this strip of input faces, unfolded with the
    ///   input lengths (`ic-extract-intrinsic-edge`, pending; evaluated with `it-unfold-apex`
    ///   and `it-ray-edge-intersection`); clamps per SPEC I5.
    /// * `n_e = −1`: the edge runs along an input edge (named by the roundabout at an input
    ///   vertex end or by `q⁰` at an inserted end); no crossings.
    ///
    /// The returned [`EdgeTrace`] follows the signpost crate's conventions: a crossing's
    /// `halfedge` is the input halfedge through which the path leaves an input face, `t` is
    /// measured in the direction of the input edge, `distance` along the straight segment from
    /// the tail (in the plane of the unfolded input faces), `first_face` is the input face in
    /// which the path leaves its tail. `mismatch` is the difference between the segment's
    /// length in the unfolded strip and the edge's length (for an edge along an input edge, the
    /// difference from the length of the piece of the input edge between its ends); `clamp` the
    /// largest clamp. Unlike a signpost trace, neither is a reason to fail: the connectivity is
    /// exact. `options.tolerance` is checked for being finite and not negative, as the signpost
    /// crate does, and is not used otherwise; `options.max_crossings` is not used (the work is
    /// bounded by the integers).
    ///
    /// # Errors
    /// [`IntegerError::IndexOutOfRange`]; [`IntegerError::Intrinsic`] for invalid options or a
    /// degenerate input face; [`IntegerError::InvalidCoordinates`] naming the edge, a crossing,
    /// a vertex or a face when the integers are not those of a realization.
    pub fn extract_intrinsic_edge(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<EdgeTrace, IntegerError> {
        self.extract_intrinsic_edge_with_report(e, options)
            .map(|(trace, _)| trace)
    }

    /// [`extract_intrinsic_edge`](Self::extract_intrinsic_edge) with its report: clamps (SPEC
    /// I5), the relative length mismatch and the work done.
    ///
    /// # Errors
    /// As [`extract_intrinsic_edge`](Self::extract_intrinsic_edge).
    pub fn extract_intrinsic_edge_with_report(
        &self,
        e: usize,
        options: &EdgeTraceOptions,
    ) -> Result<(EdgeTrace, StripReport), IntegerError> {
        let ne = self.complex.n_edges();
        if e >= ne {
            return Err(IntegerError::IndexOutOfRange {
                what: "edge",
                index: e,
                bound: ne,
            });
        }
        if !(options.tolerance.is_finite() && options.tolerance >= 0.0) {
            return Err(IntegerError::Intrinsic(IntrinsicError::InvalidOptions {
                reason: format!(
                    "the tolerance of an edge trace must be finite and not negative, got {}",
                    options.tolerance
                ),
            }));
        }
        let h = self.complex.edge_halfedge(e);
        let (a, b) = (self.complex.from_vertex(h), self.complex.to_vertex(h));
        let vertices = [a as u32, b as u32];
        let length = self.lengths[e];
        let n = self.normal[e];
        let mut report = StripReport::default();
        if n == -1 {
            let input_edge = self.input_edge_along(e)?;
            let piece = (self.parameter_on_input_edge(b, input_edge)?
                - self.parameter_on_input_edge(a, input_edge)?)
            .abs()
                * self.input.geometry().edge_lengths()[input_edge];
            let mismatch = (piece - length).abs();
            report.note_mismatch(mismatch / length);
            let trace = EdgeTrace {
                edge: e as u32,
                vertices,
                crossings: Vec::new(),
                along_input_edge: Some(input_edge as u32),
                first_face: None,
                mismatch,
                clamp: 0.0,
            };
            return Ok((trace, report));
        }
        if n < -1 {
            return Err(IntegerError::invalid(
                ValidityItem::Ranges,
                Place::Edge(e as u32),
                format!("normal coordinate {n} is below −1"),
            ));
        }
        let (labels, steps) = self.identify_crossings_of(e)?;
        report.walk_steps = steps;
        let strip = self.input_strip_of(e, &labels)?;
        let mesh = self.input.mesh();
        let input_lengths = self.input.geometry().edge_lengths();
        let laid = lay_out_strip(
            &|f| [0, 1, 2].map(|k| input_lengths[mesh.edge(3 * f + k)]),
            &|g| mesh.twin(g),
            strip.first_face,
            &strip.exits,
        )
        .map_err(|fault| match fault {
            StripFault::Degenerate { face, lengths } => {
                IntegerError::Intrinsic(IntrinsicError::DegenerateFace {
                    surface: "input mesh",
                    face: face as u32,
                    lengths,
                })
            }
            // The exits were chosen inside the faces they leave, and with twins.
            StripFault::NotInFace(i) | StripFault::Boundary(i) => invalid_edge(
                e,
                format!("crossing {i}: the strip of input faces is broken"),
            ),
        })?;
        let place = |end: EndPlace, corners: &Corners, face: usize| -> [f64; 2] {
            match end {
                EndPlace::Corner { out } => corners[out % 3],
                EndPlace::OnEdge { edge, t } => {
                    let s = (0..3)
                        .find(|&s| mesh.edge(3 * face + s) == edge)
                        .expect("the side was found when the strip was built");
                    let g = 3 * face + s;
                    let forward = mesh.from_vertex(g) == mesh.edge_vertices(edge)[0];
                    let th = if forward { t } else { 1.0 - t };
                    lerp(corners[s], corners[(s + 1) % 3], th)
                }
                EndPlace::InFace { bary, .. } => [0, 1].map(|c| {
                    bary[0] * corners[0][c] + bary[1] * corners[1][c] + bary[2] * corners[2][c]
                }),
            }
        };
        let p = place(strip.start, &laid.first, strip.first_face);
        let q = place(strip.end, &laid.last, laid.last_face);
        let chord = Chord::new(p, q).ok_or_else(|| {
            invalid_edge(
                e,
                "the tail and the head are the same point of the unfolded input faces".into(),
            )
        })?;
        let mismatch = (chord.len - length).abs();
        report.note_mismatch(mismatch / length);
        let mut crossings = Vec::with_capacity(strip.exits.len());
        let mut max_clamp = 0.0f64;
        for (&g, side) in strip.exits.iter().zip(&laid.sides) {
            let hit = chord.hit(side[0], side[1]);
            report.note_hit(&hit);
            max_clamp = max_clamp.max(hit.clamp);
            let edge = mesh.edge(g);
            let forward = mesh.from_vertex(g) == mesh.edge_vertices(edge)[0];
            crossings.push(InputCrossing {
                halfedge: g as u32,
                edge: edge as u32,
                t: if forward { hit.t } else { 1.0 - hit.t },
                distance: hit.x,
            });
        }
        let trace = EdgeTrace {
            edge: e as u32,
            vertices,
            crossings,
            along_input_edge: None,
            first_face: Some(strip.first_face as u32),
            mismatch,
            clamp: max_clamp,
        };
        Ok((trace, report))
    }

    /// The input edge along which edge `e` (with `n = −1`) runs, named by the roundabouts of its
    /// halfedges with tails in `V⁰` and by `q⁰` of its inserted ends; all must agree.
    fn input_edge_along(&self, e: usize) -> Result<usize, IntegerError> {
        let mesh = self.input.mesh();
        let h = self.complex.edge_halfedge(e);
        let ends = [
            (self.complex.from_vertex(h), Some(h)),
            (self.complex.to_vertex(h), self.complex.twin(h)),
        ];
        let mut named: Option<usize> = None;
        for (v, g) in ends {
            let this = if self.is_input_vertex(v) {
                match g {
                    Some(g) => {
                        let r = self.roundabout(g).ok_or_else(|| {
                            invalid_edge(
                                e,
                                format!(
                                    "halfedge {g} leaves input vertex {v} without a roundabout"
                                ),
                            )
                        })?;
                        if r >= self.local.degree(v) {
                            return Err(invalid_edge(
                                e,
                                format!("the roundabout {r} of halfedge {g} is not below deg₀"),
                            ));
                        }
                        Some(self.local.edge(v, r))
                    }
                    // The head of a boundary edge: named by the other end.
                    None => None,
                }
            } else {
                match self.locations[v] {
                    SurfacePoint::Edge { edge, .. } => Some(edge as usize),
                    other => {
                        return Err(invalid_edge(
                            e,
                            format!(
                                "the edge runs along an input edge, but its end {v} lies at {other}"
                            ),
                        ));
                    }
                }
            };
            match (named, this) {
                (Some(x), Some(y)) if x != y => {
                    return Err(IntegerError::invalid(
                        ValidityItem::AlongInputEdge,
                        Place::Edge(e as u32),
                        format!("its ends name different input edges, {x} and {y}"),
                    ));
                }
                (None, y) => named = y,
                _ => {}
            }
        }
        let input_edge = named
            .ok_or_else(|| invalid_edge(e, "no end of the edge names an input edge".into()))?;
        if input_edge >= mesh.n_edges() {
            return Err(invalid_edge(
                e,
                format!("input edge {input_edge} does not exist"),
            ));
        }
        Ok(input_edge)
    }

    /// For every crossing `p` of edge `e` (from the tail of its `edge_halfedge`), the input edge
    /// that crosses it, and the number of crossings visited.
    ///
    /// The crossings are identified in order. A walk from crossing `p` stops at a crossing of
    /// `e` whose input edge is known (an earlier index, or a later one an earlier walk passed);
    /// the later crossings of `e` it passes are remembered in a map, so memory grows with the
    /// work done, not with a (possibly corrupt) `n_e`.
    fn identify_crossings_of(&self, e: usize) -> Result<(Vec<u32>, u64), IntegerError> {
        let h = self.complex.edge_halfedge(e);
        let twin = self.complex.twin(h);
        let n = self.normal[e].max(0);
        if n as u64 > MAX_CERTIFIED_CROSSINGS {
            return Err(IntegerError::TooManyCrossings {
                total: n as u128,
                limit: MAX_CERTIFIED_CROSSINGS,
            });
        }
        let along = |z: Crossing| -> Option<i64> {
            let g = z.halfedge as usize;
            if g == h {
                Some(z.index)
            } else if Some(g) == twin {
                Some(n - 1 - z.index)
            } else {
                None
            }
        };
        /// What is known of a crossing of `e` beyond the one being identified.
        #[derive(Clone, Copy)]
        enum Later {
            Pending,
            Known(u32),
        }
        let mut labels: Vec<u32> = Vec::with_capacity((n as usize).min(1 << 16));
        let mut later: HashMap<i64, Later> = HashMap::new();
        let bound = Some(self.walk_bound());
        let mut steps = 0u64;
        for p in 0..n {
            if let Some(Later::Known(label)) = later.remove(&p) {
                labels.push(label);
                continue;
            }
            let mut pending: Vec<i64> = Vec::new();
            let mut found: Option<u32> = None;
            let mut closed: Option<Crossing> = None;
            let mut at_start = true;
            let walked = trace_with(self, Crossing::new(h, p), bound, &mut |z| {
                steps += 1;
                let Some(q) = along(z) else {
                    return Ok(());
                };
                if q == p {
                    // The start; a second visit means the curve is closed.
                    if std::mem::replace(&mut at_start, false) {
                        return Ok(());
                    }
                    closed = Some(z);
                    return Err(stop_walk());
                }
                if q < p {
                    found = Some(labels[q as usize]);
                    return Err(stop_walk());
                }
                match later.get(&q) {
                    Some(Later::Known(label)) => {
                        found = Some(*label);
                        Err(stop_walk())
                    }
                    Some(Later::Pending) => {
                        closed = Some(z);
                        Err(stop_walk())
                    }
                    None => {
                        later.insert(q, Later::Pending);
                        pending.push(q);
                        Ok(())
                    }
                }
            });
            let label = match (walked, found, closed) {
                (_, Some(label), _) => label,
                (_, _, Some(z)) => {
                    return Err(IntegerError::invalid(
                        ValidityItem::Extraction,
                        z.place(),
                        format!(
                            "crossing visited twice: the curve through crossing {p} of edge {e} \
                             is closed"
                        ),
                    ));
                }
                (Err(err), ..) => return Err(err),
                (Ok((arrival, _)), ..) => {
                    // Identify from the end the walk reached: the curve leaves that vertex by
                    // `arrival.back`.
                    let segment = CurveSegment {
                        start: arrival.corner,
                        start_vertex: arrival.vertex,
                        crossings: vec![arrival.back],
                        end: arrival,
                    };
                    self.identify_traced_edge(&segment)?.input_edge
                }
            };
            for q in pending {
                later.insert(q, Later::Known(label));
            }
            labels.push(label);
        }
        Ok((labels, steps))
    }

    /// The face of the input mesh in which halfedge `g` of the triangulation leaves its tail
    /// `v ∈ V⁰` strictly inside a corner (its edge does not run along an input edge): with
    /// roundabout `r`, the corner between the input halfedges with local indices `r − 1` and
    /// `r`. Returns the outgoing input halfedge with index `r − 1`; its face is the face.
    fn wedge_at(&self, e: usize, v: usize, g: usize) -> Result<usize, IntegerError> {
        let mesh = self.input.mesh();
        let r = self.roundabout(g).ok_or_else(|| {
            invalid_edge(
                e,
                format!("halfedge {g} leaves input vertex {v} without a roundabout"),
            )
        })?;
        let deg = self.local.degree(v);
        let k = if self.local.is_boundary(v) {
            // Directions at a boundary vertex lie between index 0 and deg₀ − 1, which has no
            // outgoing input halfedge.
            if r == 0 || r >= deg {
                return Err(invalid_edge(
                    e,
                    format!(
                        "halfedge {g} leaves boundary vertex {v} with roundabout {r}, which is \
                         not a direction inside the surface (deg₀ = {deg})"
                    ),
                ));
            }
            r - 1
        } else {
            if r >= deg {
                return Err(invalid_edge(
                    e,
                    format!("the roundabout {r} of halfedge {g} is not below deg₀ = {deg}"),
                ));
            }
            (r + deg - 1) % deg
        };
        let edge = self.local.edge(v, k);
        let x = mesh.edge_halfedge(edge);
        let out = if mesh.from_vertex(x) == v {
            Some(x)
        } else {
            mesh.twin(x).filter(|&t| mesh.from_vertex(t) == v)
        };
        out.ok_or_else(|| {
            invalid_edge(
                e,
                format!(
                    "input edge {edge} (direction {k} at vertex {v}) has no halfedge leaving it"
                ),
            )
        })
    }

    /// The strip of input faces of edge `e`, given the input edge of every crossing.
    fn input_strip_of(&self, e: usize, labels: &[u32]) -> Result<InputStrip, IntegerError> {
        let mesh = self.input.mesh();
        let h = self.complex.edge_halfedge(e);
        let (a, b) = (self.complex.from_vertex(h), self.complex.to_vertex(h));
        // The candidate places of the tail: one, or two for a vertex on an interior input edge.
        let starts: Vec<(usize, EndPlace)> = if self.is_input_vertex(a) {
            let out = self.wedge_at(e, a, h)?;
            vec![(out / 3, EndPlace::Corner { out })]
        } else {
            match self.locations[a] {
                SurfacePoint::Edge { edge, t } => {
                    let edge = edge as usize;
                    if edge >= mesh.n_edges() {
                        return Err(invalid_edge(
                            e,
                            format!("its tail lies on input edge {edge}, which does not exist"),
                        ));
                    }
                    let x = mesh.edge_halfedge(edge);
                    std::iter::once(x)
                        .chain(mesh.twin(x))
                        .map(|g| (g / 3, EndPlace::OnEdge { edge, t }))
                        .collect()
                }
                SurfacePoint::Face { face, bary } if (face as usize) < mesh.n_faces() => {
                    vec![(face as usize, EndPlace::InFace { bary })]
                }
                other => {
                    return Err(invalid_edge(
                        e,
                        format!("its tail {a} is an inserted vertex at {other}"),
                    ));
                }
            }
        };
        let mut first_error = None;
        for (face, start) in starts {
            match self.follow_input_strip(e, a, b, face, start, labels) {
                Ok(strip) => return Ok(strip),
                Err(err) => {
                    first_error.get_or_insert(err);
                }
            }
        }
        Err(first_error.expect("there is at least one candidate start"))
    }

    /// Follows the identified input edges from the start face; checks the arrival.
    fn follow_input_strip(
        &self,
        e: usize,
        a: usize,
        b: usize,
        first_face: usize,
        start: EndPlace,
        labels: &[u32],
    ) -> Result<InputStrip, IntegerError> {
        let mesh = self.input.mesh();
        let side_edges = |f: usize| [0, 1, 2].map(|k| mesh.edge(3 * f + k));
        let mut exits = Vec::with_capacity(labels.len());
        let mut face = first_face;
        // The halfedges of the current face through which the path may leave it.
        let mut allowed: Vec<usize> = match start {
            // Leaving a corner, the path crosses the side opposite it.
            EndPlace::Corner { out } => vec![mesh.next(out)],
            EndPlace::OnEdge { edge, .. } => (0..3)
                .map(|k| 3 * face + k)
                .filter(|&g| mesh.edge(g) != edge)
                .collect(),
            EndPlace::InFace { .. } => (0..3).map(|k| 3 * face + k).collect(),
        };
        if let EndPlace::OnEdge { edge, .. } = start
            && allowed.len() != 2
        {
            return Err(invalid_edge(
                e,
                format!(
                    "its tail lies on input edge {edge}, which is not a side of input face {face}"
                ),
            ));
        }
        let mut entry: Option<usize> = None;
        for (i, &label) in labels.iter().enumerate() {
            let label = label as usize;
            let Some(&g) = allowed.iter().find(|&&g| mesh.edge(g) == label) else {
                return Err(invalid_edge(
                    e,
                    if i == 0 {
                        format!(
                            "its first crossing is with input edge {label}, which does not bound \
                             input face {face} (sides {:?}) where the edge leaves its tail {a} \
                             on that side",
                            side_edges(face)
                        )
                    } else {
                        format!(
                            "crossings {} and {i} are with input edges {} and {label}, which do \
                             not bound a common input face (face {face} has sides {:?})",
                            i - 1,
                            labels[i - 1],
                            side_edges(face)
                        )
                    },
                ));
            };
            let Some(t) = mesh.twin(g) else {
                return Err(invalid_edge(
                    e,
                    format!("crossing {i} is with the boundary input edge {label}"),
                ));
            };
            exits.push(g);
            entry = Some(t);
            face = t / 3;
            allowed = vec![mesh.next(t), mesh.prev(t)];
        }
        // The arrival at the head.
        let end = if self.is_input_vertex(b) {
            let Some(g) = self.complex.twin(self.complex.edge_halfedge(e)) else {
                return Err(invalid_edge(
                    e,
                    "a boundary edge that does not run along an input edge".into(),
                ));
            };
            let out = self.wedge_at(e, b, g)?;
            let opposite_ok = match (entry, start) {
                (Some(t), _) => mesh.next(out) == t,
                // In one face: from another corner, or from a point of the side opposite b.
                (None, EndPlace::Corner { out: o }) => o != out,
                (None, EndPlace::OnEdge { edge, .. }) => mesh.edge(mesh.next(out)) == edge,
                (None, EndPlace::InFace { .. }) => true,
            };
            if out / 3 != face || !opposite_ok {
                return Err(invalid_edge(
                    e,
                    format!(
                        "the path reaches input face {face}, but the roundabout at its head {b} \
                         gives the corner of input halfedge {out} in face {}",
                        out / 3
                    ),
                ));
            }
            EndPlace::Corner { out }
        } else {
            match self.locations[b] {
                SurfacePoint::Edge { edge, t } => {
                    let edge = edge as usize;
                    let on_side = (0..3)
                        .map(|k| 3 * face + k)
                        .any(|g| Some(g) != entry && mesh.edge(g) == edge);
                    // In one face: not on the side the tail lies on, nor on a side at the
                    // tail's corner (the edge would run along that side).
                    let along_start = entry.is_none()
                        && match start {
                            EndPlace::OnEdge { edge: s, .. } => s == edge,
                            EndPlace::Corner { out } => mesh.edge(mesh.next(out)) != edge,
                            EndPlace::InFace { .. } => false,
                        };
                    if !on_side || along_start {
                        return Err(invalid_edge(
                            e,
                            format!(
                                "the path reaches input face {face} (sides {:?}), but its head \
                                 {b} lies on input edge {edge}",
                                side_edges(face)
                            ),
                        ));
                    }
                    EndPlace::OnEdge { edge, t }
                }
                SurfacePoint::Face { face: f, bary } if f as usize == face => {
                    EndPlace::InFace { bary }
                }
                other => {
                    return Err(invalid_edge(
                        e,
                        format!(
                            "the path reaches input face {face}, but its head {b} lies at {other}"
                        ),
                    ));
                }
            }
        };
        Ok(InputStrip {
            first_face,
            exits,
            start,
            end,
        })
    }
}
