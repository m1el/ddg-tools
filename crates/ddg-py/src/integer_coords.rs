//! `IntrinsicTriangulation`: the intrinsic triangulation of a mesh as a session object, on
//! integer coordinates (Gillespie, Sharp, Crane 2021, crate `ddg-integer-coords`) or on
//! signposts (Sharp, Soliman, Crane 2019, crate `ddg-intrinsic`). Integer Coordinates contract,
//! WP11.
//!
//! Plain Rust, as the rest of this crate: the PyO3 glue (`python/integer_coords.rs`) decodes
//! the buffers and calls [`Session`]. A session is built from positions and triangles
//! ([`Session::new`], optional mollification), flipped to Delaunay, refined, and then asked for
//!
//! * the **common subdivision** as arrays ([`Session::subdivision`]): its vertices with their
//!   places on the input mesh and on the triangulation, its convex polygons, and for every
//!   polygon its input face and its face of the triangulation (Fig. 12 of the integer paper);
//! * the **intrinsic edges as polylines on the input mesh** ([`Session::edge_polylines`],
//!   `ExtractIntrinsicEdge` for integers, tracing for signposts), for drawing; the same
//!   function serves FlipOut's "Show triangulation" later;
//! * **transfer** of vertex functions between the triangulation and the input
//!   ([`Session::to_input`], [`Session::to_intrinsic`], integer coordinates only);
//! * the **certificate** ([`Session::certify`], integer coordinates only).
//!
//! The default representation of the existing bindings is signposts; the representation of a
//! session is chosen by the caller (`"integer"` or `"signposts"`).

use ddg_integer_coords::{InitOptions, IntegerRefineOptions, IntegerTriangulation, Transfer};
use ddg_intrinsic::{
    DelaunayOptions, DelaunayStatus, EdgeTraceOptions, InsertionLimit, IntrinsicTriangulation,
    RefineOptions, RefineRules, RefineStatus, SurfacePoint, Triangulation,
};
use ddg_intrinsic_apps::representation::Representation;

use crate::error::{InputError, UnwrapError};
use crate::intrinsic_apps::{Mesh, MeshInput, bad_option, build, parse_representation, place};

enum State {
    Signposts(Box<IntrinsicTriangulation>),
    Integer(Box<IntegerTriangulation>),
}

/// An intrinsic triangulation of a mesh with the operations a front end needs.
pub struct Session {
    mesh: Mesh,
    state: State,
    transfer: Option<Transfer>,
    /// The relative mollification margin that was applied (integer coordinates only).
    pub mollified: Option<f64>,
}

/// Counts of a session.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Counts {
    /// Vertices of the triangulation (the input's first).
    pub n_vertices: usize,
    /// Edges.
    pub n_edges: usize,
    /// Faces.
    pub n_faces: usize,
    /// Vertices of the input mesh.
    pub n_input_vertices: usize,
    /// Faces of the input mesh.
    pub n_input_faces: usize,
    /// The smallest corner angle of the triangulation, radians.
    pub min_angle: f64,
    /// The largest corner angle, radians.
    pub max_angle: f64,
    /// Crossings of the intrinsic edges with input edges (integer coordinates; 0 for signposts,
    /// which do not store them).
    pub n_crossings: u64,
}

/// What [`Session::flip_to_delaunay`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct FlipOutput {
    /// `"delaunay"`, `"flip_limit"` or `"not_flippable"`.
    pub status: &'static str,
    /// Flips made.
    pub n_flips: usize,
    /// Edges that fail the criterion and cannot be flipped.
    pub n_not_delaunay: usize,
}

/// Options of [`Session::refine`].
#[derive(Debug, Clone, PartialEq)]
pub struct RefineArguments {
    /// The angle bound, degrees (0 to 35).
    pub min_angle_degrees: f64,
    /// `"chew"` (the integer paper's rules: boundary edges split at their midpoint, nearby
    /// inserted vertices removed) or `"ruppert"` (the signpost paper's rules).
    pub rules: String,
    /// At most this many insertions (`None`: ten per vertex, at least 1000).
    pub max_insertions: Option<usize>,
    /// A bound on the circumradius (`None`: none).
    pub max_circumradius: Option<f64>,
}

/// What [`Session::refine`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct RefineOutput {
    /// `"done"` or `"insertion_limit"`.
    pub status: &'static str,
    /// Vertices inserted.
    pub n_inserted: usize,
    /// Of them by splitting a face, and by splitting an edge.
    pub n_face_insertions: usize,
    /// Edge splits.
    pub n_edge_splits: usize,
    /// Flips.
    pub n_flips: usize,
    /// Inserted vertices removed again (`"chew"`), and removals refused.
    pub n_removals: usize,
    /// Removals refused.
    pub n_removals_refused: usize,
    /// The smallest angle over the faces refinement is responsible for (narrow-cone faces
    /// left out), radians.
    pub min_angle_valid: f64,
    /// Faces that stay below the bound.
    pub n_below_bound: usize,
    /// The certificate's counts after refinement (integer coordinates).
    pub certified: bool,
}

/// Points of the input mesh or of the triangulation, as arrays.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Places {
    /// 0 vertex, 1 edge, 2 face.
    pub kind: Vec<i32>,
    /// The vertex, edge or face id.
    pub id: Vec<u32>,
    /// Three numbers per point: `(0,0,0)` for a vertex, `(t,0,0)` for an edge (from the
    /// lower-numbered end of the input edge or the tail of the triangulation's edge), the
    /// barycentric coordinates for a face.
    pub params: Vec<f64>,
}

impl Places {
    pub(crate) fn push(&mut self, p: SurfacePoint) {
        let (kind, id, params) = match p {
            SurfacePoint::Vertex(v) => (0, v, [0.0; 3]),
            SurfacePoint::Edge { edge, t } => (1, edge, [t, 0.0, 0.0]),
            SurfacePoint::Face { face, bary } => (2, face, bary),
        };
        self.kind.push(kind);
        self.id.push(id);
        self.params.extend_from_slice(&params);
    }
}

/// The common subdivision of the input mesh and the triangulation.
#[derive(Debug, Clone, PartialEq)]
pub struct SubdivisionOutput {
    /// Positions in space of the vertices, 3 reals each (on the input surface).
    pub positions: Vec<f64>,
    /// The places of the vertices on the input mesh.
    pub on_input: Places,
    /// The places of the vertices on the triangulation.
    pub on_triangulation: Places,
    /// Where every polygon starts in `corners`; one more entry at the end.
    pub polygon_starts: Vec<u32>,
    /// The corners of all polygons, counter-clockwise, as vertex indices.
    pub corners: Vec<u32>,
    /// For every polygon, the input face it lies in.
    pub input_face: Vec<u32>,
    /// For every polygon, the face of the triangulation it lies in.
    pub intrinsic_face: Vec<u32>,
    /// Vertices of the triangulation come first in `positions`; this many.
    pub n_intrinsic_vertices: usize,
}

/// The intrinsic edges as polylines on the input mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct PolylinesOutput {
    /// The intrinsic edge ids drawn, in order.
    pub edges: Vec<u32>,
    /// The points of all polylines, 3 reals each (in space, on the input surface).
    pub points: Vec<f64>,
    /// Where each polyline starts in `points` (in points), one more entry at the end.
    pub offsets: Vec<u32>,
    /// For every polyline: the input edge it lies along (`-1` if none).
    pub along_input_edge: Vec<i64>,
    /// For every polyline: its number of crossings with input edges.
    pub n_crossings: Vec<u32>,
    /// The largest mismatch of an edge trace relative to the edge's length.
    pub max_relative_mismatch: f64,
}

/// The certificate of SPEC §2.4.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifyOutput {
    /// Faces checked.
    pub n_faces: usize,
    /// Edges that lie along an input edge.
    pub n_along: usize,
    /// Crossings.
    pub n_crossings: u64,
    /// Input edges extracted.
    pub n_input_edges: usize,
    /// Inserted vertices on input edges.
    pub n_inserted_on_edges: usize,
    /// The largest number of crossings on one input edge.
    pub max_crossings_per_input_edge: u64,
}

fn solve(e: ddg_intrinsic::IntrinsicError) -> UnwrapError {
    UnwrapError::Solve(e.to_string())
}

fn integer_error(e: ddg_integer_coords::IntegerError) -> UnwrapError {
    UnwrapError::Solve(e.to_string())
}

fn unsupported(operation: &str, representation: &str) -> UnwrapError {
    InputError::Rejected {
        name: "representation",
        reason: format!(
            "{operation} needs integer coordinates, this session uses {representation}"
        ),
    }
    .into()
}

impl Session {
    /// The triangulation `T¹ = T⁰` of the mesh.
    ///
    /// `mollify` is the relative margin `δ / h̄` of intrinsic mollification (integer
    /// coordinates only; the paper uses `1e-4`).
    ///
    /// # Errors
    /// [`UnwrapError::Input`] for malformed arrays, an unknown representation or a
    /// mollification with signposts, [`UnwrapError::Mesh`] for a mesh that is not a manifold
    /// triangle mesh.
    pub fn new(
        input: &MeshInput,
        representation: &str,
        mollify: Option<f64>,
    ) -> Result<Session, UnwrapError> {
        let kind = parse_representation(representation)?;
        if let Some(m) = mollify {
            if !(m.is_finite() && m >= 0.0) {
                return Err(bad_option(
                    "mollify",
                    m.to_string(),
                    "must be finite and >= 0",
                ));
            }
            if kind == Representation::Signposts {
                return Err(bad_option(
                    "mollify",
                    m.to_string(),
                    "only for representation=\"integer\"",
                ));
            }
        }
        let mesh = build(input)?;
        let mesh_error = |e: String| UnwrapError::Mesh(e);
        let state = match kind {
            Representation::Signposts => State::Signposts(Box::new(
                IntrinsicTriangulation::new(&mesh.mesh, &mesh.geometry)
                    .map_err(|e| mesh_error(e.to_string()))?,
            )),
            Representation::Integer => State::Integer(Box::new(match mollify {
                None => IntegerTriangulation::new(&mesh.mesh, &mesh.geometry)
                    .map_err(|e| mesh_error(e.to_string()))?,
                Some(m) => {
                    IntegerTriangulation::with_options(
                        &mesh.mesh,
                        mesh.geometry.edge_lengths(),
                        &InitOptions { mollify: Some(m) },
                    )
                    .map_err(|e| mesh_error(e.to_string()))?
                    .0
                }
            })),
        };
        Ok(Session {
            mesh,
            state,
            transfer: None,
            mollified: mollify,
        })
    }

    /// `"integer"` or `"signposts"`.
    pub fn representation(&self) -> &'static str {
        match self.state {
            State::Signposts(_) => "signposts",
            State::Integer(_) => "integer",
        }
    }

    fn with<R>(&self, f: impl FnOnce(&dyn Dyn) -> R) -> R {
        match &self.state {
            State::Signposts(t) => f(t.as_ref()),
            State::Integer(t) => f(t.as_ref()),
        }
    }

    /// Counts and angle range.
    pub fn counts(&self) -> Counts {
        let (n_crossings, tri): (u64, &dyn Dyn) = match &self.state {
            State::Signposts(t) => (0, t.as_ref()),
            State::Integer(t) => (
                u64::try_from(t.total_crossings()).unwrap_or(u64::MAX),
                t.as_ref(),
            ),
        };
        let c = tri.complex_ref();
        let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
        for h in 0..c.n_halfedges() {
            let a = tri.angle(h);
            lo = lo.min(a);
            hi = hi.max(a);
        }
        Counts {
            n_vertices: c.n_vertices(),
            n_edges: c.n_edges(),
            n_faces: c.n_faces(),
            n_input_vertices: self.mesh.mesh.n_vertices(),
            n_input_faces: self.mesh.mesh.n_faces(),
            min_angle: lo,
            max_angle: hi,
            n_crossings,
        }
    }

    /// Flips edges until every edge is Delaunay (paper Sec. 4.1).
    ///
    /// # Errors
    /// [`UnwrapError::Solve`] if a flip fails.
    pub fn flip_to_delaunay(&mut self) -> Result<FlipOutput, UnwrapError> {
        self.transfer = None;
        let options = DelaunayOptions::default();
        let report = match &mut self.state {
            State::Signposts(t) => t.flip_to_delaunay(&options),
            State::Integer(t) => t.flip_to_delaunay(&options),
        }
        .map_err(solve)?;
        Ok(FlipOutput {
            status: match report.status {
                DelaunayStatus::Delaunay => "delaunay",
                DelaunayStatus::FlipLimit => "flip_limit",
                _ => "not_flippable",
            },
            n_flips: report.n_flips,
            n_not_delaunay: report.n_not_delaunay,
        })
    }

    /// Delaunay refinement to an angle bound (paper Sec. 4.2; integer paper Supp. Alg. 10).
    ///
    /// # Errors
    /// [`UnwrapError::Input`] for an angle outside 0..35 degrees, unknown rules or a bad
    /// bound, [`UnwrapError::Solve`] if refinement fails.
    pub fn refine(&mut self, args: &RefineArguments) -> Result<RefineOutput, UnwrapError> {
        let rules = match args.rules.as_str() {
            "chew" => RefineRules::Chew,
            "ruppert" => RefineRules::Ruppert,
            other => {
                return Err(bad_option(
                    "rules",
                    format!("{other:?}"),
                    "expected \"chew\" or \"ruppert\"",
                ));
            }
        };
        if !(args.min_angle_degrees.is_finite() && (0.0..=35.0).contains(&args.min_angle_degrees)) {
            return Err(bad_option(
                "min_angle_degrees",
                args.min_angle_degrees.to_string(),
                "must be between 0 and 35",
            ));
        }
        if let Some(r) = args.max_circumradius
            && !(r.is_finite() && r > 0.0)
        {
            return Err(bad_option(
                "max_circumradius",
                r.to_string(),
                "must be finite and > 0",
            ));
        }
        self.transfer = None;
        let options = RefineOptions {
            min_angle: args.min_angle_degrees.to_radians(),
            max_circumradius: args.max_circumradius,
            max_insertions: args.max_insertions.map_or(
                RefineOptions::default().max_insertions,
                InsertionLimit::Vertices,
            ),
            ..RefineOptions::default()
        };
        let (report, certified) = match &mut self.state {
            State::Signposts(t) => (t.refine_with_rules(&options, rules).map_err(solve)?, false),
            State::Integer(t) => {
                let r = t
                    .refine(&IntegerRefineOptions {
                        refine: options,
                        rules,
                        ..IntegerRefineOptions::default()
                    })
                    .map_err(integer_error)?;
                (r.refine, r.certificate.is_some())
            }
        };
        Ok(RefineOutput {
            status: match report.status {
                RefineStatus::Done => "done",
                RefineStatus::InsertionLimit => "insertion_limit",
            },
            n_inserted: report.n_inserted,
            n_face_insertions: report.n_face_insertions,
            n_edge_splits: report.n_edge_splits,
            n_flips: report.n_flips,
            n_removals: report.n_removals,
            n_removals_refused: report.n_removals_refused,
            min_angle_valid: report.min_angle_valid,
            n_below_bound: report.below_bound.len(),
            certified,
        })
    }

    /// The faces of the triangulation (3 vertex ids each) and its edge lengths.
    pub fn triangles(&self) -> (Vec<u32>, Vec<f64>) {
        self.with(|t| {
            (
                t.complex_ref().faces().iter().flatten().copied().collect(),
                t.lengths().to_vec(),
            )
        })
    }

    /// The places on the input mesh of the vertices of the triangulation, with their positions
    /// in space.
    pub fn vertex_places(&self) -> (Places, Vec<f64>) {
        let mut places = Places::default();
        let mut xyz = Vec::new();
        self.with(|t| {
            for v in 0..t.complex_ref().n_vertices() {
                let p = t.location(v);
                places.push(p);
                xyz.extend_from_slice(&place(&self.mesh, p));
            }
        });
        (places, xyz)
    }

    /// The common subdivision of the input mesh and the triangulation (integer paper Sec. 4.5;
    /// signposts: the overlay of the traced edges).
    ///
    /// # Errors
    /// [`UnwrapError::Solve`] if the subdivision cannot be built.
    pub fn subdivision(&self) -> Result<SubdivisionOutput, UnwrapError> {
        let overlay = match &self.state {
            State::Integer(t) => t.common_subdivision().map_err(integer_error)?.overlay,
            State::Signposts(t) => {
                let traces = t.trace_edges(&EdgeTraceOptions::default()).map_err(solve)?;
                t.overlay(&traces).map_err(solve)?
            }
        };
        let mut on_input = Places::default();
        let mut on_triangulation = Places::default();
        let mut positions = Vec::with_capacity(3 * overlay.vertices.len());
        for v in &overlay.vertices {
            on_input.push(v.on_input);
            on_triangulation.push(v.on_triangulation);
            positions.extend_from_slice(&place(&self.mesh, v.on_input));
        }
        Ok(SubdivisionOutput {
            positions,
            on_input,
            on_triangulation,
            polygon_starts: overlay.face_starts,
            corners: overlay.corners,
            input_face: overlay.face_input,
            intrinsic_face: overlay.face_triangulation,
            n_intrinsic_vertices: self.counts().n_vertices,
        })
    }

    /// The intrinsic edges `edges` (all if `None`) as polylines on the input mesh: from the
    /// place of one end through the crossings with input edges to the other end.
    ///
    /// # Errors
    /// [`UnwrapError::Input`] for an edge id out of range, [`UnwrapError::Solve`] if a trace
    /// fails.
    pub fn edge_polylines(&self, edges: Option<&[u32]>) -> Result<PolylinesOutput, UnwrapError> {
        let n_edges = self.counts().n_edges;
        let ids: Vec<u32> = match edges {
            Some(e) => {
                if let Some(position) = e.iter().position(|&x| x as usize >= n_edges) {
                    return Err(InputError::IndexOutOfRange {
                        name: "edges",
                        position,
                        value: e[position],
                        count: n_edges,
                        of: "edges",
                    }
                    .into());
                }
                e.to_vec()
            }
            None => (0..n_edges as u32).collect(),
        };
        let mut out = PolylinesOutput {
            edges: ids.clone(),
            points: Vec::new(),
            offsets: Vec::with_capacity(ids.len() + 1),
            along_input_edge: Vec::with_capacity(ids.len()),
            n_crossings: Vec::with_capacity(ids.len()),
            max_relative_mismatch: 0.0,
        };
        let options = EdgeTraceOptions::default();
        let mut n_points = 0u32;
        for &e in &ids {
            let (trace, length) = match &self.state {
                State::Signposts(t) => (
                    Triangulation::trace_edge_on_input(t.as_ref(), e as usize, &options),
                    t.edge_lengths()[e as usize],
                ),
                State::Integer(t) => (
                    Triangulation::trace_edge_on_input(t.as_ref(), e as usize, &options),
                    t.edge_lengths()[e as usize],
                ),
            };
            let trace = trace.map_err(solve)?;
            out.offsets.push(n_points);
            let mut push = |p: [f64; 3]| {
                out.points.extend_from_slice(&p);
                n_points += 1;
            };
            let [a, b] = trace
                .vertices
                .map(|v| self.with(|t| t.location(v as usize)));
            push(place(&self.mesh, a));
            for c in &trace.crossings {
                push(place(
                    &self.mesh,
                    SurfacePoint::Edge {
                        edge: c.edge,
                        t: c.t,
                    },
                ));
            }
            push(place(&self.mesh, b));
            out.along_input_edge
                .push(trace.along_input_edge.map_or(-1, i64::from));
            out.n_crossings.push(trace.crossings.len() as u32);
            out.max_relative_mismatch = out.max_relative_mismatch.max(trace.mismatch / length);
        }
        out.offsets.push(n_points);
        Ok(out)
    }

    /// The certificate of the integer coordinates (SPEC §2.4).
    ///
    /// # Errors
    /// [`UnwrapError::Input`] for signposts, [`UnwrapError::Solve`] naming the item that fails.
    pub fn certify(&self) -> Result<CertifyOutput, UnwrapError> {
        match &self.state {
            State::Signposts(_) => Err(unsupported("certify", "signposts")),
            State::Integer(t) => {
                let r = t.certify().map_err(integer_error)?;
                Ok(CertifyOutput {
                    n_faces: r.n_faces,
                    n_along: r.n_along,
                    n_crossings: r.n_crossings,
                    n_input_edges: r.n_input_edges,
                    n_inserted_on_edges: r.n_inserted_on_edges,
                    max_crossings_per_input_edge: r.max_crossings_per_input_edge,
                })
            }
        }
    }

    fn transfer(&mut self) -> Result<&Transfer, UnwrapError> {
        let State::Integer(t) = &self.state else {
            return Err(unsupported("transfer", "signposts"));
        };
        if !self.transfer.as_ref().is_some_and(|x| x.is_current_for(t)) {
            self.transfer = Some(Transfer::from_state(t).map_err(integer_error)?);
        }
        Ok(self.transfer.as_ref().expect("just built"))
    }

    /// Transfers a function on the vertices of the triangulation to the input vertices
    /// (`L²` projection, Eq. 13-14).
    ///
    /// # Errors
    /// [`UnwrapError::Input`] for signposts, [`UnwrapError::Solve`] for a wrong length.
    pub fn to_input(&mut self, values: &[f64]) -> Result<Vec<f64>, UnwrapError> {
        self.transfer()?.to_input(values).map_err(integer_error)
    }

    /// Transfers a function on the input vertices to the vertices of the triangulation.
    ///
    /// # Errors
    /// As [`to_input`](Self::to_input).
    pub fn to_intrinsic(&mut self, values: &[f64]) -> Result<Vec<f64>, UnwrapError> {
        self.transfer()?.to_intrinsic(values).map_err(integer_error)
    }
}

/// What [`Session`] needs of both representations, object safe.
trait Dyn {
    fn complex_ref(&self) -> &ddg_intrinsic::Complex;
    fn lengths(&self) -> &[f64];
    fn location(&self, v: usize) -> SurfacePoint;
    fn angle(&self, h: usize) -> f64;
}

impl<T: Triangulation> Dyn for T {
    fn complex_ref(&self) -> &ddg_intrinsic::Complex {
        self.complex()
    }
    fn lengths(&self) -> &[f64] {
        self.edge_lengths()
    }
    fn location(&self, v: usize) -> SurfacePoint {
        self.vertex_location(v)
    }
    fn angle(&self, h: usize) -> f64 {
        self.corner_angle(h)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators::{grid_disk, icosphere};

    fn input(g: ddg_core::generators::GeneratedMesh) -> MeshInput {
        MeshInput {
            positions: g.positions.iter().flatten().copied().collect(),
            triangles: g.faces.iter().flatten().copied().collect(),
        }
    }

    fn refine_args(rules: &str) -> RefineArguments {
        RefineArguments {
            min_angle_degrees: 25.0,
            rules: rules.into(),
            max_insertions: None,
            max_circumradius: None,
        }
    }

    fn prepared(g: ddg_core::generators::GeneratedMesh, rep: &str, rules: &str) -> Session {
        let mut s = Session::new(&input(g), rep, None).unwrap();
        s.flip_to_delaunay().unwrap();
        s.refine(&refine_args(rules)).unwrap();
        s
    }

    #[test]
    fn both_representations_agree_on_the_subdivision() {
        for rules in ["ruppert", "chew"] {
            let a = prepared(icosphere(2).unwrap(), "signposts", rules);
            let b = prepared(icosphere(2).unwrap(), "integer", rules);
            let (ca, cb) = (a.counts(), b.counts());
            assert_eq!(
                (ca.n_vertices, ca.n_faces),
                (cb.n_vertices, cb.n_faces),
                "{rules}"
            );
            assert!((ca.min_angle - cb.min_angle).abs() < 1e-9);
            assert!(
                cb.min_angle.to_degrees() > 24.9,
                "{rules}: {}",
                cb.min_angle.to_degrees()
            );
            let (sa, sb) = (a.subdivision().unwrap(), b.subdivision().unwrap());
            assert_eq!(sa.polygon_starts, sb.polygon_starts, "{rules}");
            assert_eq!(sa.input_face, sb.input_face);
            assert_eq!(sa.intrinsic_face, sb.intrinsic_face);
            assert_eq!(sb.positions.len(), 3 * sb.on_input.kind.len());
            assert_eq!(
                *sb.polygon_starts.last().unwrap() as usize,
                sb.corners.len()
            );
            assert!(sb.intrinsic_face.iter().all(|&f| (f as usize) < cb.n_faces));
            assert!(
                sb.input_face
                    .iter()
                    .all(|&f| (f as usize) < cb.n_input_faces)
            );
            // Every face of the triangulation has pieces; every input face too.
            for n in [cb.n_faces, cb.n_input_faces] {
                let which = if n == cb.n_faces {
                    &sb.intrinsic_face
                } else {
                    &sb.input_face
                };
                let mut seen = vec![false; n];
                which.iter().for_each(|&f| seen[f as usize] = true);
                assert!(seen.iter().all(|&x| x));
            }
            for (x, y) in sa.positions.iter().zip(&sb.positions) {
                assert!((x - y).abs() < 1e-8, "{x} {y}");
            }
        }
    }

    #[test]
    fn polylines_start_and_end_at_the_vertices() {
        let s = prepared(icosphere(2).unwrap(), "integer", "chew");
        let p = s.edge_polylines(None).unwrap();
        let c = s.counts();
        assert_eq!(p.edges.len(), c.n_edges);
        assert_eq!(p.offsets.len(), c.n_edges + 1);
        assert_eq!(*p.offsets.last().unwrap() as usize * 3, p.points.len());
        assert!(p.max_relative_mismatch < 1e-6);
        let (_, xyz) = s.vertex_places();
        assert_eq!(xyz.len(), 3 * c.n_vertices);
        let total: u32 = p.n_crossings.iter().sum();
        assert_eq!(u64::from(total), c.n_crossings);
        let few = s.edge_polylines(Some(&[3, 5])).unwrap();
        assert_eq!(few.edges, vec![3, 5]);
        assert!(matches!(
            s.edge_polylines(Some(&[u32::MAX])).unwrap_err(),
            UnwrapError::Input(_)
        ));
    }

    #[test]
    fn certify_and_transfer() {
        let mut s = prepared(icosphere(2).unwrap(), "integer", "chew");
        let cert = s.certify().unwrap();
        assert_eq!(cert.n_faces, s.counts().n_faces);
        let c = s.counts();
        let ones = vec![1.0; c.n_vertices];
        let down = s.to_input(&ones).unwrap();
        assert_eq!(down.len(), c.n_input_vertices);
        assert!(down.iter().all(|x| (x - 1.0).abs() < 1e-10));
        let up = s.to_intrinsic(&down).unwrap();
        assert!(up.iter().all(|x| (x - 1.0).abs() < 1e-10));
        assert!(s.to_input(&[1.0]).is_err());
        // After a change the transfer is rebuilt.
        s.refine(&refine_args("chew")).unwrap();
        let _ = s.to_input(&vec![2.0; s.counts().n_vertices]).unwrap();
        let mut sp = prepared(icosphere(1).unwrap(), "signposts", "ruppert");
        assert!(matches!(sp.certify().unwrap_err(), UnwrapError::Input(_)));
        assert!(matches!(
            sp.to_input(&[]).unwrap_err(),
            UnwrapError::Input(_)
        ));
    }

    #[test]
    fn mollification_and_bad_options_name_the_parameter() {
        let g = || input(grid_disk(4, 4).unwrap());
        let s = Session::new(&g(), "integer", Some(1e-4)).unwrap();
        assert_eq!(s.mollified, Some(1e-4));
        for (rep, mollify, name) in [
            ("signposts", Some(1e-4), "mollify"),
            ("integer", Some(-1.0), "mollify"),
            ("nope", None, "representation"),
        ] {
            let e = Session::new(&g(), rep, mollify).err().unwrap().to_string();
            assert!(e.starts_with(name), "{e}");
        }
        let mut s = Session::new(&g(), "integer", None).unwrap();
        for (args, name) in [
            (
                RefineArguments {
                    rules: "x".into(),
                    ..refine_args("chew")
                },
                "rules",
            ),
            (
                RefineArguments {
                    min_angle_degrees: 40.0,
                    ..refine_args("chew")
                },
                "min_angle_degrees",
            ),
            (
                RefineArguments {
                    max_circumradius: Some(-1.0),
                    ..refine_args("chew")
                },
                "max_circumradius",
            ),
        ] {
            let e = s.refine(&args).unwrap_err().to_string();
            assert!(e.starts_with(name), "{e}");
        }
    }

    #[test]
    fn a_flat_disk_with_boundary_refines_on_both() {
        for rules in ["chew", "ruppert"] {
            let a = prepared(grid_disk(6, 6).unwrap(), "signposts", rules);
            let b = prepared(grid_disk(6, 6).unwrap(), "integer", rules);
            assert_eq!(a.counts().n_faces, b.counts().n_faces);
            b.certify().unwrap();
        }
    }
}
