//! Extraction of the paths on the input mesh `T⁰` and the certificate computed there (SPEC §5,
//! §7 layer 1b).
//!
//! [`GeodesicNetwork::extract`] returns one [`PathOnInput`] per path: the polyline on `T⁰`
//! that the edge path of `T` is. Each segment is the trace of its edge
//! ([`Triangulation::trace_edge_on_input`]: for signposts the walk over `T⁰` checked by IT N8,
//! for integer coordinates `ExtractIntrinsicEdge`, exact in the crossed edges), reversed for a
//! backward segment; the traces of the distinct edges are computed in parallel and assembled
//! in path order, so the result does not depend on the number of threads (F24). A trace that
//! fails is an error naming the path, segment and edge (F18); nothing is repaired.
//!
//! The polyline of a closed path starts and ends at the same point (the first point is
//! repeated), so a path of `n` segments has `P + 1` points for `P` pieces either way, and the
//! joint between the last and the first segment is the first point.
//!
//! # The certificate on `T⁰`
//!
//! [`InputCertificate`] is computed from the polylines and the *input* mesh alone (its edge
//! lengths and corner angles), independently of `T`; every piece is laid out in the plane in
//! its input face:
//!
//! * **Crossings.** At every crossing of an input edge, the pieces before and after, the
//!   second face unfolded across the edge (the unfolding of the diamond of `it-unfold-diamond`
//!   is used as a construction; no formula of the Lean tree is evaluated), make an angle
//!   `π`. The largest deviation is reported (`max_crossing_deviation`, in radians) and must
//!   be at most `crossing_tolerance`.
//! * **Joints.** At every joint, which is an input vertex, the two sides measured on `T⁰` are
//!   the input corner angles between the faces plus the partial angles of the first and last
//!   piece. Eq. 2 (`fo-locally-shortest` in the SPEC; not tagged, F5 has not put it in Lean
//!   yet) must hold for both sides within `input_tolerance` unless the joint is exempt on `T`
//!   (pinned, or every bent side blocked by a path end, [`classify`]). Independently of
//!   that, the sides measured here are compared with the wedge angles `T` computed
//!   (`max_wedge_mismatch`).
//! * **Lengths.** The length of the polyline equals the length of the edge path on `T`
//!   within `length_tolerance` relative.

use std::f64::consts::PI;

use ddg_core::TriMesh;
use ddg_intrinsic::{EdgeTrace, EdgeTraceOptions, PathOnInput, SurfacePoint, Triangulation};
use rayon::prelude::*;

use crate::certify::{JointStatus, classify};
use crate::error::FlipError;
use crate::network::{GeodesicNetwork, Joint, Network};
use crate::wedge::{JointGeometry, Side};

/// Options of [`GeodesicNetwork::extract`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtractOptions {
    /// The options of the trace of an edge over the input mesh (IT `EdgeTraceOptions`; for
    /// integer coordinates only validated).
    pub trace: EdgeTraceOptions,
    /// Compute the certificate on `T⁰` ([`Extraction::certificate`]); default on.
    pub certify_input: bool,
    /// Largest deviation from `π` (radians) of the angle at a crossing of an input edge;
    /// default `1e-6`.
    pub crossing_tolerance: f64,
    /// A joint that is not exempt must have both sides `≥ π − input_tolerance` on `T⁰`
    /// (radians); default `1e-4`, the certificate tolerance of `T` (F3).
    pub input_tolerance: f64,
    /// Largest relative difference between the length of a polyline and of its edge path on
    /// `T`; default `1e-9`.
    pub length_tolerance: f64,
}

impl Default for ExtractOptions {
    fn default() -> Self {
        ExtractOptions {
            trace: EdgeTraceOptions::default(),
            certify_input: true,
            crossing_tolerance: 1e-6,
            input_tolerance: 1e-4,
            length_tolerance: 1e-9,
        }
    }
}

impl ExtractOptions {
    /// Checks the ranges.
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`] naming the option.
    pub fn validate(&self) -> Result<(), FlipError> {
        for (name, x) in [
            ("crossing_tolerance", self.crossing_tolerance),
            ("input_tolerance", self.input_tolerance),
            ("length_tolerance", self.length_tolerance),
            ("trace.tolerance", self.trace.tolerance),
        ] {
            if !(x.is_finite() && x >= 0.0) {
                return Err(FlipError::InvalidOptions {
                    what: format!("{name} must be finite and >= 0, got {x}"),
                });
            }
        }
        Ok(())
    }
}

/// What the extraction of one edge did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeExtraction {
    /// The edge of `T`.
    pub edge: u32,
    /// Number of input edges its trace crosses (`0` for an edge along an input edge).
    pub n_crossings: u32,
    /// The mismatch of the trace at the far end (absolute): for signposts the distance of IT
    /// N8, for integer coordinates the difference between the segment's length in the unfolded
    /// strip and the edge's length (`length_mismatch`, SPEC §5).
    pub mismatch: f64,
    /// The largest clamp of a position (IC I5; zero for signposts).
    pub clamp: f64,
    /// The length of the edge in `T`.
    pub length: f64,
}

/// The counts of an extraction.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtractReport {
    /// The distinct edges that carry a segment, ascending, with what their extraction did.
    pub edges: Vec<EdgeExtraction>,
    /// Number of pieces over all polylines.
    pub n_pieces: usize,
    /// Number of crossings of input edges over all polylines (a crossing used by two segments
    /// on one edge counts twice).
    pub n_crossings: usize,
    /// The largest `mismatch / length` over the edges.
    pub max_relative_mismatch: f64,
    /// The largest clamp over the edges.
    pub max_clamp: f64,
}

/// The certificate on the input mesh (module documentation).
#[derive(Debug, Clone, PartialEq)]
pub struct InputCertificate {
    /// Whether nothing failed: no fault, and every maximum within its tolerance.
    pub passed: bool,
    /// The tolerances used: `(crossing_tolerance, input_tolerance, length_tolerance)`.
    pub tolerances: [f64; 3],
    /// Crossings of input edges checked.
    pub n_crossings: usize,
    /// The largest deviation from `π` of the angle at a crossing (radians); `0` without
    /// crossings.
    pub max_crossing_deviation: f64,
    /// Joints checked against Eq. 2 (not exempt).
    pub n_joints: usize,
    /// Joints exempt (pinned, or blocked by a path end).
    pub n_exempt: usize,
    /// The smallest side over the joints checked on `T⁰`; `+∞` if there are none. Must be at
    /// least `π − input_tolerance`.
    pub min_joint_angle: f64,
    /// The largest difference between a wedge angle measured on `T⁰` and the one computed on
    /// `T`, over all joints and both finite sides (radians).
    pub max_wedge_mismatch: f64,
    /// The largest relative difference between the length of a polyline and of its edge path
    /// on `T`.
    pub max_length_error: f64,
    /// Per path: the relative difference of the lengths.
    pub length_errors: Vec<f64>,
    /// Number of inconsistencies found (a point not in the face of its piece, two pieces that
    /// do not meet across the crossed edge, a degenerate input face, ...).
    pub n_faults: usize,
    /// The first of them, in words naming the path and point.
    pub faults: Vec<String>,
}

/// The paths on the input mesh with their reports ([`GeodesicNetwork::extract_with_report`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Extraction {
    /// One polyline per path, in path order (a contracted path has none: it is an empty
    /// polyline `points = []`).
    pub paths: Vec<PathOnInput>,
    /// Per path and segment: the index in `paths[p].points` of the point where the segment
    /// starts (so segment `i` runs through points `starts[i] ..= starts[i + 1]`; the last
    /// segment ends at the last point).
    pub segment_starts: Vec<Vec<usize>>,
    /// What the extraction did.
    pub report: ExtractReport,
    /// The certificate on `T⁰`, if asked for.
    pub certificate: Option<InputCertificate>,
}

// ----- the input mesh in the plane -----

/// The input mesh with its lengths: lays out faces in the plane.
struct Input<'a> {
    mesh: &'a TriMesh,
    lengths: &'a [f64],
    corner_angles: &'a [[f64; 3]],
}

type Layout = [[f64; 2]; 3];

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn norm(a: [f64; 2]) -> f64 {
    a[0].hypot(a[1])
}

impl Input<'_> {
    fn face_lengths(&self, f: usize) -> [f64; 3] {
        [0, 1, 2].map(|k| self.lengths[self.mesh.edge(3 * f + k)])
    }

    /// The corners of face `f` in the plane, counter-clockwise: corner `k` is the tail of
    /// halfedge `3f + k`.
    fn layout(&self, f: usize) -> Result<Layout, String> {
        let [l0, l1, l2] = self.face_lengths(f);
        let x = (l0 * l0 + l2 * l2 - l1 * l1) / (2.0 * l0);
        let y2 = l2 * l2 - x * x;
        if !(y2 > 0.0 && l0 > 0.0) {
            return Err(format!(
                "input face {f} is degenerate (lengths {l0}, {l1}, {l2})"
            ));
        }
        Ok([[0.0, 0.0], [l0, 0.0], [x, y2.sqrt()]])
    }

    fn corner_of(&self, f: usize, v: usize) -> Option<usize> {
        self.mesh.faces()[f].iter().position(|&x| x as usize == v)
    }

    /// The place of `point` in face `f` laid out as `layout`.
    fn locate(&self, f: usize, layout: &Layout, point: &SurfacePoint) -> Option<[f64; 2]> {
        match *point {
            SurfacePoint::Vertex(v) => self.corner_of(f, v as usize).map(|k| layout[k]),
            SurfacePoint::Edge { edge, t } => {
                let [lo, hi] = self.mesh.edge_vertices(edge as usize);
                let (a, b) = (self.corner_of(f, lo)?, self.corner_of(f, hi)?);
                // The edge must be a side of the face (two vertices of a face can also be joined
                // by another edge).
                (0..3).find(|&s| self.mesh.edge(3 * f + s) == edge as usize)?;
                Some([
                    (1.0 - t) * layout[a][0] + t * layout[b][0],
                    (1.0 - t) * layout[a][1] + t * layout[b][1],
                ])
            }
            SurfacePoint::Face { face, bary } => (face as usize == f).then(|| {
                [
                    bary[0] * layout[0][0] + bary[1] * layout[1][0] + bary[2] * layout[2][0],
                    bary[0] * layout[0][1] + bary[1] * layout[1][1] + bary[2] * layout[2][1],
                ]
            }),
        }
    }

    /// Face `b` laid out in the plane of face `a`'s `layout`, across `edge`, on the other side
    /// of it from `a`.
    fn unfold(&self, a: usize, layout: &Layout, b: usize, edge: usize) -> Option<Layout> {
        let ka = (0..3).find(|&k| self.mesh.edge(3 * a + k) == edge)?;
        let (u, w) = (self.mesh.faces()[a][ka] as usize, self.mesh.faces()[a][(ka + 1) % 3] as usize);
        let (pu, pw) = (layout[ka], layout[(ka + 1) % 3]);
        let apex_a = layout[(ka + 2) % 3];
        let (ku, kw) = (self.corner_of(b, u)?, self.corner_of(b, w)?);
        (0..3).find(|&k| self.mesh.edge(3 * b + k) == edge)?;
        let kc = 3 - ku - kw;
        // The side between corners i and j is the halfedge slot s with {s, s+1} = {i, j}.
        let side = |i: usize, j: usize| if (i + 1) % 3 == j { i } else { j };
        let d_u = self.lengths[self.mesh.edge(3 * b + side(kc, ku))];
        let d_w = self.lengths[self.mesh.edge(3 * b + side(kc, kw))];
        let e = sub(pw, pu);
        let l = norm(e);
        let dir = [e[0] / l, e[1] / l];
        let x = (d_u * d_u - d_w * d_w + l * l) / (2.0 * l);
        let y = (d_u * d_u - x * x).max(0.0).sqrt();
        let left = [-dir[1], dir[0]];
        let side_a = cross(dir, sub(apex_a, pu)).signum();
        let apex = [
            pu[0] + x * dir[0] - side_a * y * left[0],
            pu[1] + x * dir[1] - side_a * y * left[1],
        ];
        let mut out = [[0.0; 2]; 3];
        out[ku] = pu;
        out[kw] = pw;
        out[kc] = apex;
        Some(out)
    }
}

// ----- traces to pieces -----

/// One edge's path over `T⁰`, in the direction of its `edge_halfedge`.
struct EdgePieces {
    points: Vec<SurfacePoint>,
    faces: Vec<u32>,
    info: EdgeExtraction,
}

fn pieces_of_edge<T: Triangulation>(
    tri: &T,
    e: usize,
    options: &EdgeTraceOptions,
) -> Result<EdgePieces, FlipError> {
    let trace: EdgeTrace = tri.trace_edge_on_input(e, options)?;
    let mesh = tri.input().mesh();
    let length = tri.edge_length(e);
    let [a, b] = trace.vertices;
    let (from, to) = (tri.vertex_location(a as usize), tri.vertex_location(b as usize));
    let info = EdgeExtraction {
        edge: e as u32,
        n_crossings: trace.crossings.len() as u32,
        mismatch: trace.mismatch,
        clamp: trace.clamp,
        length,
    };
    if let Some(input_edge) = trace.along_input_edge {
        return Ok(EdgePieces {
            points: vec![from, to],
            faces: vec![(mesh.edge_halfedge(input_edge as usize) / 3) as u32],
            info,
        });
    }
    let internal = |what: String| FlipError::Internal { what };
    let mut face = trace
        .first_face
        .ok_or_else(|| internal(format!("the trace of edge {e} has no first input face")))?;
    let mut points = Vec::with_capacity(trace.crossings.len() + 2);
    let mut faces = Vec::with_capacity(trace.crossings.len() + 1);
    points.push(from);
    for crossing in &trace.crossings {
        faces.push(face);
        points.push(SurfacePoint::Edge {
            edge: crossing.edge,
            t: crossing.t,
        });
        face = mesh
            .twin(crossing.halfedge as usize)
            .map(|t| (t / 3) as u32)
            .ok_or_else(|| {
                internal(format!(
                    "the trace of edge {e} crosses boundary halfedge {} of the input mesh",
                    crossing.halfedge
                ))
            })?;
    }
    faces.push(face);
    points.push(to);
    Ok(EdgePieces {
        points,
        faces,
        info,
    })
}

// ----- the certificate -----

/// The certificate's findings for one path, merged in path order.
#[derive(Default)]
struct PathCertificate {
    n_crossings: usize,
    max_crossing_deviation: f64,
    n_joints: usize,
    n_exempt: usize,
    min_joint_angle: f64,
    max_wedge_mismatch: f64,
    length_error: f64,
    faults: Vec<String>,
    n_faults: usize,
}

const MAX_FAULTS: usize = 16;

impl PathCertificate {
    fn fault(&mut self, what: String) {
        self.n_faults += 1;
        if self.faults.len() < MAX_FAULTS {
            self.faults.push(what);
        }
    }
}

/// The position along the fan of the corners of vertex `v` of the direction from `v` towards
/// `towards` inside face `f` (laid out as `layout`): the sum of the corner angles before the
/// corner of `f`, plus the angle from the outgoing edge of that corner.
fn fan_position(
    input: &Input<'_>,
    v: usize,
    f: usize,
    layout: &Layout,
    towards: [f64; 2],
) -> Result<f64, String> {
    let k = input
        .corner_of(f, v)
        .ok_or_else(|| format!("vertex {v} is not a corner of input face {f}"))?;
    let h = 3 * f + k;
    let origin = layout[k];
    let u = sub(layout[(k + 1) % 3], origin);
    let w = sub(towards, origin);
    let theta = input.corner_angles[f][k];
    let phi = cross(u, w).atan2(dot(u, w));
    if !(-1e-6..=theta + 1e-6).contains(&phi) {
        return Err(format!(
            "the piece leaves vertex {v} outside input face {f} (angle {phi} from its edge, \
             corner {theta})"
        ));
    }
    let mut before = 0.0;
    for g in input.mesh.outgoing_halfedges(v) {
        if g == h {
            return Ok(before + phi.clamp(0.0, theta));
        }
        before += input.corner_angles[g / 3][g % 3];
    }
    Err(format!("corner {h} is not in the fan of vertex {v}"))
}

fn vertex_angle_sum(input: &Input<'_>, v: usize) -> f64 {
    input
        .mesh
        .outgoing_halfedges(v)
        .map(|g| input.corner_angles[g / 3][g % 3])
        .fold(0.0, |s, x| s + x)
}

/// Certifies path `p` (polyline `path` with segment starts `starts`) on the input mesh.
#[allow(clippy::too_many_arguments)]
fn certify_path<T: Triangulation>(
    tri: &T,
    network: &Network,
    input: &Input<'_>,
    p: usize,
    path: &PathOnInput,
    starts: &[usize],
    edge_path_length: f64,
    options: &ExtractOptions,
) -> PathCertificate {
    let mut cert = PathCertificate {
        min_joint_angle: f64::INFINITY,
        ..PathCertificate::default()
    };
    let n_points = path.points.len();
    if n_points < 2 {
        return cert;
    }
    let n_pieces = n_points - 1;
    // Layouts of the faces of the pieces, once each.
    let mut layouts: Vec<Option<Layout>> = Vec::with_capacity(n_pieces);
    for &f in &path.faces {
        match input.layout(f as usize) {
            Ok(l) => layouts.push(Some(l)),
            Err(what) => {
                cert.fault(format!("path {p}: {what}"));
                layouts.push(None);
            }
        }
    }
    let place = |piece: usize, point: usize| -> Option<[f64; 2]> {
        let layout = layouts[piece].as_ref()?;
        input.locate(path.faces[piece] as usize, layout, &path.points[point])
    };
    // Lengths.
    let mut length = 0.0;
    for i in 0..n_pieces {
        match (place(i, i), place(i, i + 1)) {
            (Some(a), Some(b)) => length += norm(sub(a, b)),
            _ => {
                if layouts[i].is_some() {
                    cert.fault(format!(
                        "path {p}, piece {i}: an end is not a point of input face {}",
                        path.faces[i]
                    ));
                }
            }
        }
    }
    cert.length_error = (length - edge_path_length).abs() / edge_path_length;
    // Crossings of input edges.
    for i in 1..n_pieces {
        let SurfacePoint::Edge { edge, .. } = path.points[i] else {
            continue;
        };
        let (fa, fb) = (path.faces[i - 1] as usize, path.faces[i] as usize);
        let (Some(la), Some(_)) = (layouts[i - 1].as_ref(), layouts[i].as_ref()) else {
            continue;
        };
        let Some(lb) = input.unfold(fa, la, fb, edge as usize) else {
            cert.fault(format!(
                "path {p}, point {i}: input faces {fa} and {fb} of the pieces around the \
                 crossing of input edge {edge} do not share that edge"
            ));
            continue;
        };
        let at = input.locate(fa, la, &path.points[i]);
        let before = input.locate(fa, la, &path.points[i - 1]);
        let after = input.locate(fb, &lb, &path.points[i + 1]);
        let (Some(at), Some(before), Some(after)) = (at, before, after) else {
            cert.fault(format!(
                "path {p}, point {i}: a neighbour of the crossing of input edge {edge} is not \
                 in the face of its piece"
            ));
            continue;
        };
        let (u, w) = (sub(before, at), sub(after, at));
        let angle = cross(u, w).abs().atan2(dot(u, w));
        cert.n_crossings += 1;
        cert.max_crossing_deviation = cert.max_crossing_deviation.max((PI - angle).abs());
    }
    // Joints.
    let ids = network.path_segment_ids(p);
    let n = ids.len();
    let closed = network.is_closed(p);
    let first = if closed { 0 } else { 1 };
    for i in first..n {
        let joint = Joint {
            path: p,
            incoming: ids[(i + n - 1) % n],
            outgoing: ids[i],
        };
        let j = starts[i];
        // The pieces around the joint and the points beyond them.
        let (piece_in, piece_out, before, after) = if j == 0 {
            (n_pieces - 1, 0, n_points - 2, 1)
        } else {
            (j - 1, j, j - 1, j + 1)
        };
        let SurfacePoint::Vertex(v) = path.points[j] else {
            continue;
        };
        let v = v as usize;
        let (Some(l_in), Some(l_out)) = (layouts[piece_in].as_ref(), layouts[piece_out].as_ref())
        else {
            continue;
        };
        let (f_in, f_out) = (path.faces[piece_in] as usize, path.faces[piece_out] as usize);
        let (Some(towards_a), Some(towards_c)) = (
            input.locate(f_in, l_in, &path.points[before]),
            input.locate(f_out, l_out, &path.points[after]),
        ) else {
            cert.fault(format!(
                "path {p}, joint {i} at vertex {v}: a neighbouring point is not in the face of its piece"
            ));
            continue;
        };
        let positions = fan_position(input, v, f_in, l_in, towards_a)
            .and_then(|a| fan_position(input, v, f_out, l_out, towards_c).map(|c| (a, c)));
        let (pos_in, pos_out) = match positions {
            Ok(x) => x,
            Err(what) => {
                cert.fault(format!("path {p}, joint {i}: {what}"));
                continue;
            }
        };
        // The left wedge sweeps counter-clockwise from bc to ba.
        let theta = vertex_angle_sum(input, v);
        let (left, right) = if input.mesh.is_boundary_vertex(v) {
            if pos_in >= pos_out {
                (pos_in - pos_out, f64::INFINITY)
            } else {
                (f64::INFINITY, pos_out - pos_in)
            }
        } else {
            let left = (pos_in - pos_out).rem_euclid(theta);
            (left, theta - left)
        };
        // Compare with T and decide whether the joint is exempt.
        match JointGeometry::new(tri, network, joint) {
            Ok(geometry) => {
                for (mine, side) in [(left, Side::Left), (right, Side::Right)] {
                    let theirs = geometry.angle(side);
                    if mine.is_finite() && theirs.is_finite() {
                        cert.max_wedge_mismatch =
                            cert.max_wedge_mismatch.max((mine - theirs).abs());
                    } else if mine.is_finite() != theirs.is_finite() {
                        cert.fault(format!(
                            "path {p}, joint {i} at vertex {v}: the {side:?} wedge crosses the \
                             boundary on only one of T and T0"
                        ));
                    }
                }
                let (status, _) = classify(&geometry, network, options.input_tolerance);
                if matches!(status, JointStatus::Pinned | JointStatus::BlockedByEnd) {
                    cert.n_exempt += 1;
                } else {
                    cert.n_joints += 1;
                    cert.min_joint_angle = cert.min_joint_angle.min(left.min(right));
                }
            }
            Err(e) => cert.fault(format!("path {p}, joint {i}: {e}")),
        }
    }
    cert
}

// ----- the network -----

impl<T: Triangulation + Sync> GeodesicNetwork<T> {
    /// The paths on the input mesh, one per path, segments in path order (module
    /// documentation). A contracted path has an empty polyline.
    ///
    /// # Errors
    /// As [`extract_with_report`](Self::extract_with_report).
    pub fn extract(&self, options: &ExtractOptions) -> Result<Vec<PathOnInput>, FlipError> {
        self.extract_with_report(options).map(|x| x.paths)
    }

    /// [`extract`](Self::extract) with the counts and, if `options.certify_input`, the
    /// certificate on `T⁰`.
    ///
    /// The distinct edges carrying segments are traced in parallel (rayon's current pool) and the
    /// results assembled in order, so the result is bitwise the same at any number of threads
    /// (F24).
    ///
    /// # Errors
    /// [`FlipError::InvalidOptions`]; [`FlipError::Extraction`] for the first segment (in path
    /// order) whose edge could not be traced, with the representation's error as the cause (a
    /// failed trace of IT N8 for signposts).
    pub fn extract_with_report(&self, options: &ExtractOptions) -> Result<Extraction, FlipError> {
        options.validate()?;
        let tri = &self.tri;
        let paths = self.network.paths();
        let mut edges: Vec<usize> = paths
            .iter()
            .flat_map(|p| p.segments.iter().map(|s| s.edge as usize))
            .collect();
        edges.sort_unstable();
        edges.dedup();
        let traced: Vec<Result<EdgePieces, FlipError>> = edges
            .par_iter()
            .map(|&e| pieces_of_edge(tri, e, &options.trace))
            .collect();
        let slot = |e: usize| edges.binary_search(&e).expect("every segment's edge was traced");
        // The first failure in path order.
        for (p, path) in paths.iter().enumerate() {
            for (index, s) in path.segments.iter().enumerate() {
                if let Err(cause) = &traced[slot(s.edge as usize)] {
                    return Err(FlipError::Extraction {
                        path: p,
                        index,
                        edge: s.edge as usize,
                        cause: Box::new(cause.clone()),
                    });
                }
            }
        }
        let traced: Vec<EdgePieces> = traced.into_iter().map(|r| r.expect("checked")).collect();

        let mut polylines = Vec::with_capacity(paths.len());
        let mut segment_starts = Vec::with_capacity(paths.len());
        let mut report = ExtractReport {
            edges: traced.iter().map(|t| t.info).collect(),
            n_pieces: 0,
            n_crossings: 0,
            max_relative_mismatch: 0.0,
            max_clamp: 0.0,
        };
        for info in &report.edges {
            report.max_relative_mismatch = report
                .max_relative_mismatch
                .max(info.mismatch / info.length);
            report.max_clamp = report.max_clamp.max(info.clamp);
        }
        for path in &paths {
            let mut points: Vec<SurfacePoint> = Vec::new();
            let mut faces: Vec<u32> = Vec::new();
            let mut starts = Vec::with_capacity(path.segments.len());
            let mut mismatch = 0.0f64;
            for s in &path.segments {
                let t = &traced[slot(s.edge as usize)];
                starts.push(points.len().saturating_sub(1));
                // The vertex shared with the previous segment is written once.
                let skip = usize::from(!points.is_empty());
                if s.forward {
                    points.extend(t.points[skip..].iter().copied());
                    faces.extend(t.faces.iter().copied());
                } else {
                    points.extend(t.points.iter().rev().skip(skip).copied());
                    faces.extend(t.faces.iter().rev().copied());
                }
                mismatch = mismatch.max(t.info.mismatch);
            }
            report.n_pieces += faces.len();
            report.n_crossings += faces.len().saturating_sub(path.segments.len());
            polylines.push(PathOnInput {
                points,
                faces,
                mismatch,
            });
            segment_starts.push(starts);
        }

        let certificate = options.certify_input.then(|| {
            let input = Input {
                mesh: tri.input().mesh(),
                lengths: tri.input().geometry().edge_lengths(),
                corner_angles: tri.input().geometry().corner_angles(),
            };
            let parts: Vec<PathCertificate> = (0..paths.len())
                .into_par_iter()
                .map(|p| {
                    certify_path(
                        tri,
                        &self.network,
                        &input,
                        p,
                        &polylines[p],
                        &segment_starts[p],
                        paths[p].length(tri),
                        options,
                    )
                })
                .collect();
            let mut c = InputCertificate {
                passed: true,
                tolerances: [
                    options.crossing_tolerance,
                    options.input_tolerance,
                    options.length_tolerance,
                ],
                n_crossings: 0,
                max_crossing_deviation: 0.0,
                n_joints: 0,
                n_exempt: 0,
                min_joint_angle: f64::INFINITY,
                max_wedge_mismatch: 0.0,
                max_length_error: 0.0,
                length_errors: Vec::with_capacity(parts.len()),
                n_faults: 0,
                faults: Vec::new(),
            };
            for part in parts {
                c.n_crossings += part.n_crossings;
                c.max_crossing_deviation = c.max_crossing_deviation.max(part.max_crossing_deviation);
                c.n_joints += part.n_joints;
                c.n_exempt += part.n_exempt;
                c.min_joint_angle = c.min_joint_angle.min(part.min_joint_angle);
                c.max_wedge_mismatch = c.max_wedge_mismatch.max(part.max_wedge_mismatch);
                c.max_length_error = c.max_length_error.max(part.length_error);
                c.length_errors.push(part.length_error);
                c.n_faults += part.n_faults;
                for f in part.faults {
                    if c.faults.len() < MAX_FAULTS {
                        c.faults.push(f);
                    }
                }
            }
            c.passed = c.n_faults == 0
                && c.max_crossing_deviation <= options.crossing_tolerance
                && c.min_joint_angle >= PI - options.input_tolerance
                && c.max_length_error <= options.length_tolerance;
            c
        });
        Ok(Extraction {
            paths: polylines,
            segment_starts,
            report,
            certificate,
        })
    }

    /// The polylines in space (for Python and Blender): the paths of
    /// [`extract`](Self::extract) with every point interpolated barycentrically from the
    /// vertex `positions` of the input mesh ([`polyline_3d`]).
    ///
    /// # Errors
    /// As [`extract`](Self::extract); [`FlipError::IndexOutOfRange`] if `positions` is not one
    /// per input vertex.
    pub fn polylines_3d(
        &self,
        positions: &[[f64; 3]],
        options: &ExtractOptions,
    ) -> Result<Vec<Vec<[f64; 3]>>, FlipError> {
        let paths = self.extract(options)?;
        let mesh = self.tri.input().mesh();
        paths
            .iter()
            .map(|p| polyline_3d(mesh, positions, p))
            .collect()
    }
}

/// The points of a path on the input mesh in space: a vertex is its position, a point of an
/// edge is `(1 − t) p_lo + t p_hi` (the edge's vertices in increasing order), a point of a face
/// is the barycentric combination of its corners.
///
/// # Errors
/// [`FlipError::IndexOutOfRange`] if `positions` does not have one entry per vertex of the mesh
/// or a point names an element that does not exist.
pub fn polyline_3d(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    path: &PathOnInput,
) -> Result<Vec<[f64; 3]>, FlipError> {
    if positions.len() != mesh.n_vertices() {
        return Err(FlipError::IndexOutOfRange {
            what: "position count",
            index: positions.len(),
            bound: mesh.n_vertices() + 1,
        });
    }
    let combine = |terms: &[(f64, usize)]| {
        let mut out = [0.0; 3];
        for &(w, v) in terms {
            for (o, x) in out.iter_mut().zip(positions[v]) {
                *o += w * x;
            }
        }
        out
    };
    path.points
        .iter()
        .map(|point| match *point {
            SurfacePoint::Vertex(v) => {
                if v as usize >= mesh.n_vertices() {
                    return Err(FlipError::IndexOutOfRange {
                        what: "vertex",
                        index: v as usize,
                        bound: mesh.n_vertices(),
                    });
                }
                Ok(positions[v as usize])
            }
            SurfacePoint::Edge { edge, t } => {
                if edge as usize >= mesh.n_edges() {
                    return Err(FlipError::IndexOutOfRange {
                        what: "input edge",
                        index: edge as usize,
                        bound: mesh.n_edges(),
                    });
                }
                let [lo, hi] = mesh.edge_vertices(edge as usize);
                Ok(combine(&[(1.0 - t, lo), (t, hi)]))
            }
            SurfacePoint::Face { face, bary } => {
                if face as usize >= mesh.n_faces() {
                    return Err(FlipError::IndexOutOfRange {
                        what: "input face",
                        index: face as usize,
                        bound: mesh.n_faces(),
                    });
                }
                let c = mesh.faces()[face as usize];
                Ok(combine(&[
                    (bary[0], c[0] as usize),
                    (bary[1], c[1] as usize),
                    (bary[2], c[2] as usize),
                ]))
            }
        })
        .collect()
}

/// For every point of a path on the input mesh that is a crossing of an input edge, the sine of
/// the angle between the path and that edge (`None` for the other points, and where the pieces
/// around the point are inconsistent). A position along an input edge is conditioned by this
/// sine: moving the crossing by `δ` along the edge moves the path by `δ sin θ` across it (SPEC
/// §7 layer 2 compares positions as `|Δ| sin θ`).
pub fn crossing_sines<T: Triangulation>(tri: &T, path: &PathOnInput) -> Vec<Option<f64>> {
    let input = Input {
        mesh: tri.input().mesh(),
        lengths: tri.input().geometry().edge_lengths(),
        corner_angles: tri.input().geometry().corner_angles(),
    };
    let n = path.points.len();
    let mut out = vec![None; n];
    for i in 1..n.saturating_sub(1) {
        let SurfacePoint::Edge { edge, .. } = path.points[i] else {
            continue;
        };
        let (fa, fb) = (path.faces[i - 1] as usize, path.faces[i] as usize);
        let (Ok(la), Ok(_)) = (input.layout(fa), input.layout(fb)) else {
            continue;
        };
        let Some(lb) = input.unfold(fa, &la, fb, edge as usize) else {
            continue;
        };
        let [lo, hi] = input.mesh.edge_vertices(edge as usize);
        let (Some(ka), Some(kb)) = (input.corner_of(fa, lo), input.corner_of(fa, hi)) else {
            continue;
        };
        let (Some(at), Some(after)) = (
            input.locate(fa, &la, &path.points[i]),
            input.locate(fb, &lb, &path.points[i + 1]),
        ) else {
            continue;
        };
        let (d, v) = (sub(la[kb], la[ka]), sub(after, at));
        out[i] = Some(cross(d, v).abs() / (norm(d) * norm(v)));
    }
    out
}
