//! The intrinsic triangulation: connectivity, lengths, signposts, and where every vertex is on
//! the input mesh.
//!
//! # Indexing conventions
//!
//! They are those of [`ddg_core::topology`], so that a triangulation fresh from
//! [`IntrinsicTriangulation::new`] has the vertex, edge, face and halfedge ids of the input
//! mesh:
//!
//! * **Halfedge** `h = 3 f + k` runs from `faces[f][k]` to `faces[f][(k + 1) % 3]`; `next` and
//!   `prev` are arithmetic. A corner and the halfedge leaving it share one id. Boundary edges
//!   have one halfedge only.
//! * **Edges** have ids of their own (they cannot be vertex pairs: two edges may join the same
//!   two vertices). An edge keeps its id when it is flipped. [`IntrinsicTriangulation::edge_halfedge`]
//!   gives every edge a direction, from the tail to the head of that halfedge; edge parameters
//!   and the frame of an edge point refer to it.
//! * A flip moves halfedges between the two faces of the edge, so **halfedge ids are not
//!   stable** under flips; vertex and edge ids are, and the two faces keep their ids.
//!
//! # What a triangulation may contain (paper Sec. 2.1, Fig. 6)
//!
//! Two edges between the same two vertices, an edge from a vertex to itself, a triangle that
//! uses a vertex two or three times or is glued to itself along an edge, vertices of degree
//! one. It is always an oriented manifold, every triangle satisfies the strict triangle
//! inequalities, and a boundary edge of the input mesh stays a boundary edge.

use std::f64::consts::{PI, TAU};
use std::sync::Arc;

use ddg_core::geometry::{corner_angle, triangle_area};
use ddg_core::topology::NONE;
use ddg_core::{Geometry, TriMesh};

use crate::error::IntrinsicError;
use crate::formulas::signpost_update;
use crate::input::InputSurface;
use crate::point::{SurfacePoint, point_defect};
use crate::surface::{Surface, next, prev};
use crate::trace::{Trace, TraceOptions, trace_query};

/// An intrinsic triangulation of a triangle mesh in the signpost representation of Sharp,
/// Soliman and Crane (2019, Sec. 3.1).
///
/// It stores a halfedge mesh (see the [module documentation](self) for the conventions), one
/// length per edge, one angle per halfedge (its *signpost*: the direction in which the
/// halfedge leaves its tail, in the frame of that vertex), the angle sum `Θ` of every vertex,
/// and the place of every vertex on the input mesh. The input mesh itself is held as a shared
/// immutable copy ([`InputSurface`]).
#[derive(Debug, Clone, PartialEq)]
pub struct IntrinsicTriangulation {
    pub(crate) input: Arc<InputSurface>,
    pub(crate) faces: Vec<[u32; 3]>,
    /// Per halfedge: the opposite halfedge, or `NONE`.
    pub(crate) twin: Vec<u32>,
    /// Per halfedge: its edge.
    pub(crate) halfedge_edge: Vec<u32>,
    /// Per edge: the halfedge that gives it its direction (the only one on the boundary).
    pub(crate) edge_halfedge: Vec<u32>,
    /// Per vertex: an outgoing halfedge; the one without twin at a boundary vertex.
    pub(crate) vertex_halfedge: Vec<u32>,
    pub(crate) edge_lengths: Vec<f64>,
    /// Per halfedge: normalized direction at its tail.
    pub(crate) signposts: Vec<f64>,
    /// Per vertex: the sum of its corner angles, fixed when the vertex is created.
    pub(crate) angle_sums: Vec<f64>,
    /// Per vertex: its place on the input mesh.
    pub(crate) locations: Vec<SurfacePoint>,
    /// Per edge: the input edge it lies along, or `NONE`.
    pub(crate) edge_on_input: Vec<u32>,
}

/// The halfedges leaving a vertex, counter-clockwise. See
/// [`IntrinsicTriangulation::outgoing_halfedges`].
#[derive(Debug, Clone)]
pub struct OutgoingHalfedges<'a> {
    tri: &'a IntrinsicTriangulation,
    start: u32,
    current: u32,
}

impl Iterator for OutgoingHalfedges<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        if self.current == NONE {
            return None;
        }
        let h = self.current as usize;
        let t = self.tri.twin[prev(h)];
        self.current = if t == self.start { NONE } else { t };
        Some(h)
    }
}

/// The ways in which a triangulation fails to be a simplicial complex. See
/// [`IntrinsicTriangulation::irregularities`].
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Irregularities {
    /// Edges from a vertex to itself, ascending.
    pub self_edges: Vec<u32>,
    /// Edges that are not self-edges and share both end points with another edge, ascending.
    pub multi_edges: Vec<u32>,
    /// Faces that use a vertex more than once, ascending.
    pub faces_with_repeated_vertex: Vec<u32>,
    /// Interior vertices with a single incident edge end (their one triangle is glued to
    /// itself along that edge), ascending.
    pub degree_one_vertices: Vec<u32>,
}

impl Irregularities {
    /// Whether there are none.
    pub fn is_empty(&self) -> bool {
        self.self_edges.is_empty()
            && self.multi_edges.is_empty()
            && self.faces_with_repeated_vertex.is_empty()
            && self.degree_one_vertices.is_empty()
    }
}

/// How well a triangulation satisfies the identities that hold exactly in exact arithmetic.
/// Returned by [`IntrinsicTriangulation::check`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConsistencyReport {
    /// Largest `|Σ corner angles at i − Θ_i|` over the vertices, in radians.
    pub max_angle_sum_error: f64,
    /// Largest difference, over consecutive outgoing halfedges of a vertex, between the stored
    /// signpost and the previous signpost plus the normalized corner angle, in normalized
    /// radians. For a boundary vertex the last corner is checked against `π`.
    pub max_signpost_error: f64,
    /// `(total area − total area of the input mesh) / total area of the input mesh`.
    pub relative_area_error: f64,
    /// `Σ_I (2π − Θ_i) + Σ_B (π − Θ_i) − 2π χ` with the stored `Θ` and `χ = V − E + F`.
    pub gauss_bonnet_residual: f64,
    /// Smallest value over the faces of `(b + c − a) / a` with `a` the longest side: how far
    /// the thinnest triangle is from violating the triangle inequality.
    pub min_triangle_slack: f64,
}

impl IntrinsicTriangulation {
    /// The intrinsic triangulation that is the input mesh itself (paper Alg. 4): same
    /// vertices, edges, faces and halfedges, with the same ids; signposts set; every vertex
    /// located at the input vertex of the same id.
    ///
    /// The input is copied (see [`InputSurface`]); the caller's mesh is not referenced
    /// afterwards.
    ///
    /// # Errors
    /// [`IntrinsicError::Geometry`] if the geometry was not built for the mesh.
    pub fn new(
        mesh: &TriMesh,
        geometry: &Geometry,
    ) -> Result<IntrinsicTriangulation, IntrinsicError> {
        let input = InputSurface::new(mesh, geometry)?;
        let n_halfedges = mesh.n_halfedges();
        let twin = (0..n_halfedges)
            .map(|h| mesh.twin(h).map_or(NONE, |t| t as u32))
            .collect();
        let halfedge_edge = (0..n_halfedges).map(|h| mesh.edge(h) as u32).collect();
        let edge_halfedge = (0..mesh.n_edges())
            .map(|e| mesh.edge_halfedge(e) as u32)
            .collect();
        let vertex_halfedge = (0..mesh.n_vertices())
            .map(|v| mesh.outgoing_halfedge(v) as u32)
            .collect();
        Ok(IntrinsicTriangulation {
            faces: mesh.faces().to_vec(),
            twin,
            halfedge_edge,
            edge_halfedge,
            vertex_halfedge,
            edge_lengths: geometry.edge_lengths().to_vec(),
            signposts: input.signposts().to_vec(),
            angle_sums: input.vertex_angle_sums().to_vec(),
            locations: (0..mesh.n_vertices() as u32)
                .map(SurfacePoint::Vertex)
                .collect(),
            edge_on_input: (0..mesh.n_edges() as u32).collect(),
            input: Arc::new(input),
        })
    }

    // ----- the input -----

    /// The input mesh this is a triangulation of.
    #[inline]
    pub fn input(&self) -> &InputSurface {
        &self.input
    }

    /// Number of vertices of the input mesh. The vertices `0..n_input_vertices()` of the
    /// triangulation are the input vertices, with the same ids.
    #[inline]
    pub fn n_input_vertices(&self) -> usize {
        self.input.mesh().n_vertices()
    }

    // ----- counts -----

    /// Number of vertices.
    #[inline]
    pub fn n_vertices(&self) -> usize {
        self.vertex_halfedge.len()
    }

    /// Number of edges.
    #[inline]
    pub fn n_edges(&self) -> usize {
        self.edge_halfedge.len()
    }

    /// Number of faces.
    #[inline]
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of halfedges, which is also the number of corners: `3 * n_faces`.
    #[inline]
    pub fn n_halfedges(&self) -> usize {
        self.twin.len()
    }

    /// Euler characteristic `χ = V − E + F`. No operation changes it.
    pub fn euler_characteristic(&self) -> i64 {
        self.n_vertices() as i64 - self.n_edges() as i64 + self.n_faces() as i64
    }

    // ----- halfedges -----

    /// The faces as vertex triples. A face may list a vertex more than once.
    #[inline]
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// The next halfedge in the same face.
    #[inline]
    pub fn next(&self, h: usize) -> usize {
        next(h)
    }

    /// The previous halfedge in the same face.
    #[inline]
    pub fn prev(&self, h: usize) -> usize {
        prev(h)
    }

    /// The oppositely directed halfedge of the same edge, or `None` on the boundary. It can
    /// lie in the same face.
    #[inline]
    pub fn twin(&self, h: usize) -> Option<usize> {
        let t = self.twin[h];
        (t != NONE).then_some(t as usize)
    }

    /// The face of halfedge (or corner) `h`: `h / 3`.
    #[inline]
    pub fn face(&self, h: usize) -> usize {
        h / 3
    }

    /// The vertex halfedge `h` starts at; equivalently the vertex of corner `h`.
    #[inline]
    pub fn from_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][h % 3] as usize
    }

    /// The vertex halfedge `h` points to.
    #[inline]
    pub fn to_vertex(&self, h: usize) -> usize {
        self.faces[h / 3][(h + 1) % 3] as usize
    }

    /// The edge halfedge `h` belongs to.
    #[inline]
    pub fn edge(&self, h: usize) -> usize {
        self.halfedge_edge[h] as usize
    }

    // ----- edges -----

    /// The halfedge that gives edge `e` its direction. Its twin, if any, is the other halfedge
    /// of the edge. After a flip it is the halfedge from the vertex opposite the edge on its
    /// right to the vertex opposite on its left (the edge is turned counter-clockwise).
    #[inline]
    pub fn edge_halfedge(&self, e: usize) -> usize {
        self.edge_halfedge[e] as usize
    }

    /// The end points of edge `e`: tail and head of [`edge_halfedge`](Self::edge_halfedge).
    /// They are equal for a self-edge.
    #[inline]
    pub fn edge_vertices(&self, e: usize) -> [usize; 2] {
        let h = self.edge_halfedge[e] as usize;
        [self.from_vertex(h), self.to_vertex(h)]
    }

    /// Whether edge `e` has only one incident triangle side.
    #[inline]
    pub fn is_boundary_edge(&self, e: usize) -> bool {
        self.twin[self.edge_halfedge[e] as usize] == NONE
    }

    /// The length of edge `e`.
    #[inline]
    pub fn edge_length(&self, e: usize) -> f64 {
        self.edge_lengths[e]
    }

    /// The lengths of all edges, indexed by edge id.
    #[inline]
    pub fn edge_lengths(&self) -> &[f64] {
        &self.edge_lengths
    }

    /// The input edge that edge `e` lies along, if it is known to: an edge that has not been
    /// flipped since the triangulation was created. Such an edge has the id, the end points
    /// and the length of the input edge. An edge that was flipped away and back is the same
    /// segment again but is not recognized here.
    #[inline]
    pub fn edge_on_input(&self, e: usize) -> Option<usize> {
        let i = self.edge_on_input[e];
        (i != NONE).then_some(i as usize)
    }

    // ----- vertices -----

    /// An outgoing halfedge of vertex `v`: the first one of
    /// [`outgoing_halfedges`](Self::outgoing_halfedges). At a boundary vertex it is the
    /// boundary halfedge leaving `v`, whose signpost is zero.
    #[inline]
    pub fn outgoing_halfedge(&self, v: usize) -> usize {
        self.vertex_halfedge[v] as usize
    }

    /// Whether vertex `v` lies on the boundary.
    #[inline]
    pub fn is_boundary_vertex(&self, v: usize) -> bool {
        self.twin[self.vertex_halfedge[v] as usize] == NONE
    }

    /// The halfedges leaving `v`, one per incident corner, counter-clockwise. At a boundary
    /// vertex the iteration starts with the boundary halfedge leaving `v`; the last boundary
    /// neighbour is then reached only by an incoming halfedge. A self-edge at `v` appears with
    /// each of its halfedges.
    pub fn outgoing_halfedges(&self, v: usize) -> OutgoingHalfedges<'_> {
        OutgoingHalfedges {
            tri: self,
            start: self.vertex_halfedge[v],
            current: self.vertex_halfedge[v],
        }
    }

    /// The degree of `v`: the number of edge ends at `v` (a self-edge counts twice).
    pub fn vertex_degree(&self, v: usize) -> usize {
        self.outgoing_halfedges(v).count() + usize::from(self.is_boundary_vertex(v))
    }

    /// The sum `Θ` of the corner angles at `v`: that of the input mesh for an input vertex,
    /// `2π` for a vertex inserted in the interior and `π` for one inserted on the boundary. It
    /// is stored when the vertex is created and never recomputed.
    #[inline]
    pub fn vertex_angle_sum(&self, v: usize) -> f64 {
        self.angle_sums[v]
    }

    /// The place of vertex `v` on the input mesh.
    #[inline]
    pub fn vertex_location(&self, v: usize) -> SurfacePoint {
        self.locations[v]
    }

    /// The places of all vertices on the input mesh.
    #[inline]
    pub fn vertex_locations(&self) -> &[SurfacePoint] {
        &self.locations
    }

    // ----- signposts and angles -----

    /// The signpost of halfedge `h`: the direction in which it leaves its tail, as a
    /// normalized angle in the frame of that vertex (paper Sec. 2.3, 3.1).
    ///
    /// *Normalized* means that the angles measured in the triangles around the vertex are
    /// multiplied by `2π / Θ` at an interior vertex, so that the directions fill `[0, 2π)`, and
    /// by `π / Θ` at a boundary vertex, where they fill `[0, π]` from the boundary halfedge
    /// that leaves the vertex (signpost `0`) to the boundary edge that arrives at it
    /// (direction `π`; it has no halfedge leaving the vertex). The frame of an input vertex is
    /// that of the input mesh ([`InputSurface`]), so a signpost can be handed to
    /// [`InputSurface::trace`] as it is.
    #[inline]
    pub fn signpost(&self, h: usize) -> f64 {
        self.signposts[h]
    }

    /// The signposts of all halfedges.
    #[inline]
    pub fn signposts(&self) -> &[f64] {
        &self.signposts
    }

    /// The lengths of the three halfedges of face `f`, by slot.
    #[inline]
    pub fn face_lengths(&self, f: usize) -> [f64; 3] {
        [
            self.edge_lengths[self.halfedge_edge[3 * f] as usize],
            self.edge_lengths[self.halfedge_edge[3 * f + 1] as usize],
            self.edge_lengths[self.halfedge_edge[3 * f + 2] as usize],
        ]
    }

    /// The interior angle of corner `h` (the corner of face `h / 3` at the tail of halfedge
    /// `h`), from the edge lengths.
    ///
    /// # Panics
    /// Never for a triangulation built by this crate: every operation keeps the strict
    /// triangle inequalities.
    pub fn corner_angle(&self, h: usize) -> f64 {
        let l = self.face_lengths(h / 3);
        let k = h % 3;
        corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3])
            .expect("every face keeps the triangle inequalities")
    }

    /// The area of face `f`, from the edge lengths.
    ///
    /// # Panics
    /// As [`corner_angle`](Self::corner_angle).
    pub fn face_area(&self, f: usize) -> f64 {
        let [a, b, c] = self.face_lengths(f);
        triangle_area(a, b, c).expect("every face keeps the triangle inequalities")
    }

    /// The sum of all face areas. It is the area of the input mesh up to rounding.
    pub fn total_area(&self) -> f64 {
        (0..self.n_faces()).map(|f| self.face_area(f)).sum()
    }

    // ----- tracing over the triangulation -----

    /// Follows the straight path that leaves `start` (a point of this triangulation) in
    /// `direction` over the intrinsic triangulation for `length` (paper Sec. 3.2.2). See
    /// [`Trace`] for what is returned and [`SurfacePoint`] for the frame `direction` refers
    /// to; from a vertex it is a signpost direction.
    ///
    /// # Errors
    /// As [`InputSurface::trace`].
    pub fn trace(
        &self,
        start: SurfacePoint,
        direction: f64,
        length: f64,
        options: &TraceOptions,
    ) -> Result<Trace, IntrinsicError> {
        trace_query(self, start, direction, length, options)
    }

    // ----- simplicial complexes -----

    /// Lists what keeps the triangulation from being a simplicial complex.
    pub fn irregularities(&self) -> Irregularities {
        let mut out = Irregularities::default();
        let mut pairs: Vec<(u32, u32, u32)> = Vec::with_capacity(self.n_edges());
        for e in 0..self.n_edges() {
            let [a, b] = self.edge_vertices(e);
            if a == b {
                out.self_edges.push(e as u32);
            } else {
                pairs.push((a.min(b) as u32, a.max(b) as u32, e as u32));
            }
        }
        pairs.sort_unstable();
        let mut i = 0;
        while i < pairs.len() {
            let mut j = i + 1;
            while j < pairs.len() && pairs[j].0 == pairs[i].0 && pairs[j].1 == pairs[i].1 {
                j += 1;
            }
            if j - i > 1 {
                out.multi_edges.extend(pairs[i..j].iter().map(|p| p.2));
            }
            i = j;
        }
        out.multi_edges.sort_unstable();
        for (f, &[a, b, c]) in self.faces.iter().enumerate() {
            if a == b || b == c || a == c {
                out.faces_with_repeated_vertex.push(f as u32);
            }
        }
        for v in 0..self.n_vertices() {
            let h = self.vertex_halfedge[v] as usize;
            if self.twin[prev(h)] == h as u32 {
                out.degree_one_vertices.push(v as u32);
            }
        }
        out
    }

    /// The triangulation as a `TriMesh` with a `Geometry`, if it is a simplicial complex.
    /// Vertices, faces, corners and halfedges keep their ids; edges are renumbered by the
    /// `TriMesh` (find them with `TriMesh::find_edge`).
    ///
    /// # Errors
    /// [`IntrinsicError::NotSimplicial`] naming the first self-edge, double edge or repeated
    /// vertex, or the reason `TriMesh::new` gives (two faces with the same three vertices).
    pub fn to_tri_mesh(&self) -> Result<(TriMesh, Geometry), IntrinsicError> {
        let irregular = self.irregularities();
        let not_simplicial = |reason: String| IntrinsicError::NotSimplicial {
            reason,
            mesh_error: None,
        };
        if let Some(&e) = irregular.self_edges.first() {
            let [a, _] = self.edge_vertices(e as usize);
            return Err(not_simplicial(format!(
                "edge {e} joins vertex {a} to itself ({} self-edges in total)",
                irregular.self_edges.len()
            )));
        }
        if let Some(&e) = irregular.multi_edges.first() {
            let [a, b] = self.edge_vertices(e as usize);
            return Err(not_simplicial(format!(
                "edge {e} = ({a}, {b}) is one of several edges between its end points ({} such \
                 edges in total)",
                irregular.multi_edges.len()
            )));
        }
        if let Some(&f) = irregular.faces_with_repeated_vertex.first() {
            let [a, b, c] = self.faces[f as usize];
            return Err(not_simplicial(format!(
                "face {f} = ({a}, {b}, {c}) uses a vertex more than once"
            )));
        }
        let mesh = TriMesh::new(self.n_vertices(), self.faces.clone()).map_err(|e| {
            IntrinsicError::NotSimplicial {
                reason: "its faces are not a valid simplicial mesh".to_string(),
                mesh_error: Some(e),
            }
        })?;
        let mut lengths = vec![0.0; mesh.n_edges()];
        for h in 0..self.n_halfedges() {
            lengths[mesh.edge(h)] = self.edge_lengths[self.halfedge_edge[h] as usize];
        }
        let geometry = Geometry::from_edge_lengths(&mesh, lengths)?;
        Ok((mesh, geometry))
    }

    // ----- consistency -----

    /// Checks the data structure and measures how well the triangulation satisfies the
    /// identities of the surface.
    ///
    /// Structure (any violation is an error): twins are mutual and run in opposite directions;
    /// every edge has one or two halfedges and [`edge_halfedge`](Self::edge_halfedge) is one
    /// of them; the corners of every vertex form one fan that
    /// [`outgoing_halfedge`](Self::outgoing_halfedge) starts, at a boundary vertex with the
    /// boundary halfedge; every face satisfies the strict triangle inequalities; lengths,
    /// signposts and angle sums are finite, signposts are in `[0, 2π]` (`[0, π]` at a boundary
    /// vertex, `0` for its boundary halfedge); every vertex location is a point of the input
    /// mesh.
    ///
    /// Identities (measured, returned in the [`ConsistencyReport`]): angle sums, signposts
    /// against corner angles, total area, Gauss–Bonnet.
    ///
    /// # Errors
    /// [`IntrinsicError::Inconsistent`] with the first violation found, and
    /// [`IntrinsicError::DegenerateFace`].
    pub fn check(&self) -> Result<ConsistencyReport, IntrinsicError> {
        let bad = |what: String| Err(IntrinsicError::Inconsistent { what });
        let (nv, ne, nf, nh) = (
            self.n_vertices(),
            self.n_edges(),
            self.n_faces(),
            self.n_halfedges(),
        );
        if nh != 3 * nf
            || self.halfedge_edge.len() != nh
            || self.signposts.len() != nh
            || self.edge_lengths.len() != ne
            || self.edge_on_input.len() != ne
            || self.angle_sums.len() != nv
            || self.locations.len() != nv
        {
            return bad("array lengths do not match the element counts".to_string());
        }
        // Faces, twins, edges.
        for (f, face) in self.faces.iter().enumerate() {
            if let Some(&v) = face.iter().find(|&&v| v as usize >= nv) {
                return bad(format!("face {f} references vertex {v} of {nv}"));
            }
        }
        let mut halfedges_of_edge = vec![0u8; ne];
        for h in 0..nh {
            let e = self.halfedge_edge[h];
            if e as usize >= ne {
                return bad(format!("halfedge {h} references edge {e} of {ne}"));
            }
            halfedges_of_edge[e as usize] += 1;
            let t = self.twin[h];
            if t == NONE {
                continue;
            }
            let t = t as usize;
            if t >= nh || t == h || self.twin[t] != h as u32 {
                return bad(format!("the twin of halfedge {h} is not mutual"));
            }
            if self.halfedge_edge[t] != e {
                return bad(format!(
                    "halfedge {h} and its twin {t} have different edges"
                ));
            }
            if self.from_vertex(t) != self.to_vertex(h) || self.to_vertex(t) != self.from_vertex(h)
            {
                return bad(format!(
                    "halfedge {h} and its twin {t} do not run in opposite directions"
                ));
            }
        }
        for (e, &n_halfedges) in halfedges_of_edge.iter().enumerate() {
            let h = self.edge_halfedge[e] as usize;
            if h >= nh || self.halfedge_edge[h] != e as u32 {
                return bad(format!("edge {e} does not own its halfedge {h}"));
            }
            let expected = if self.twin[h] == NONE { 1 } else { 2 };
            if n_halfedges != expected {
                return bad(format!(
                    "edge {e} has {} halfedges, expected {expected}",
                    n_halfedges
                ));
            }
            let l = self.edge_lengths[e];
            if !(l.is_finite() && l > 0.0) {
                return bad(format!("edge {e} has length {l}"));
            }
            let on = self.edge_on_input[e];
            if on != NONE && on as usize >= self.input.mesh().n_edges() {
                return bad(format!(
                    "edge {e} lies along input edge {on}, which does not exist"
                ));
            }
        }
        // Triangles.
        let mut min_slack = f64::INFINITY;
        let mut area = 0.0;
        for f in 0..nf {
            let l = self.face_lengths(f);
            match triangle_area(l[0], l[1], l[2]) {
                Some(a) => area += a,
                None => {
                    return Err(IntrinsicError::DegenerateFace {
                        surface: Self::NAME,
                        face: f as u32,
                        lengths: l,
                    });
                }
            }
            let longest = l[0].max(l[1]).max(l[2]);
            let shortest = l[0].min(l[1]).min(l[2]);
            let middle = l[0] + l[1] + l[2] - longest - shortest;
            min_slack = min_slack.min((shortest - (longest - middle)) / longest);
        }
        // Fans, signposts, angle sums.
        let mut corners_at = vec![0u32; nv];
        for face in &self.faces {
            for &v in face {
                corners_at[v as usize] += 1;
            }
        }
        let mut max_angle_sum_error = 0.0f64;
        let mut max_signpost_error = 0.0f64;
        let mut curvature = 0.0;
        for (v, &n_corners) in corners_at.iter().enumerate() {
            let h0 = self.vertex_halfedge[v] as usize;
            if h0 >= nh || self.from_vertex(h0) != v {
                return bad(format!("vertex {v} does not own its halfedge {h0}"));
            }
            let theta = self.angle_sums[v];
            if !(theta.is_finite() && theta > 0.0) {
                return bad(format!("vertex {v} has angle sum {theta}"));
            }
            let boundary = self.twin[h0] == NONE;
            let full = if boundary { PI } else { TAU };
            let mut count = 0u32;
            let mut sum = 0.0;
            let mut h = h0;
            loop {
                count += 1;
                if count > n_corners {
                    return bad(format!("the corners of vertex {v} do not form a fan"));
                }
                if self.from_vertex(h) != v {
                    return bad(format!(
                        "halfedge {h} in the fan of vertex {v} leaves another vertex"
                    ));
                }
                let phi = self.signposts[h];
                if !(phi.is_finite() && (0.0..=full).contains(&phi)) {
                    return bad(format!(
                        "halfedge {h} at vertex {v} has signpost {phi}, outside [0, {full}]"
                    ));
                }
                if h != h0 && self.twin[h] == NONE {
                    return bad(format!(
                        "boundary halfedge {h} is not the first halfedge of vertex {v}"
                    ));
                }
                let corner = self.corner_angle(h);
                sum += corner;
                let expected = signpost_update(phi, corner, theta, boundary);
                let t = self.twin[prev(h)];
                if t == NONE {
                    if !boundary {
                        return bad(format!("the fan of interior vertex {v} is open"));
                    }
                    max_signpost_error = max_signpost_error.max((expected - PI).abs());
                    break;
                }
                let d = (self.signposts[t as usize] - expected).abs();
                max_signpost_error = max_signpost_error.max(d.min((d - TAU).abs()));
                h = t as usize;
                if h == h0 {
                    if boundary {
                        return bad(format!("the fan of boundary vertex {v} is closed"));
                    }
                    break;
                }
            }
            if count != n_corners {
                return bad(format!(
                    "vertex {v} has {} corners but its fan has {count}",
                    n_corners
                ));
            }
            if boundary && self.signposts[h0] != 0.0 {
                return bad(format!(
                    "the boundary halfedge {h0} of vertex {v} has signpost {}",
                    self.signposts[h0]
                ));
            }
            max_angle_sum_error = max_angle_sum_error.max((sum - theta).abs());
            curvature += full - theta;
            let input = self.input.mesh();
            if let Some(reason) = point_defect(
                &self.locations[v],
                input.n_vertices(),
                input.n_edges(),
                input.n_faces(),
            ) {
                return bad(format!(
                    "vertex {v} is located at {}, but {reason}",
                    self.locations[v]
                ));
            }
        }
        let input_area = self.input.geometry().total_area();
        Ok(ConsistencyReport {
            max_angle_sum_error,
            max_signpost_error,
            relative_area_error: (area - input_area) / input_area,
            gauss_bonnet_residual: curvature - TAU * self.euler_characteristic() as f64,
            min_triangle_slack: min_slack,
        })
    }

    pub(crate) fn check_edge(&self, e: usize) -> Result<(), IntrinsicError> {
        if e < self.n_edges() {
            Ok(())
        } else {
            Err(IntrinsicError::IndexOutOfRange {
                what: "edge",
                index: e,
                bound: self.n_edges(),
            })
        }
    }
}

impl Surface for IntrinsicTriangulation {
    const NAME: &'static str = "intrinsic triangulation";

    #[inline]
    fn n_vertices(&self) -> usize {
        self.vertex_halfedge.len()
    }
    #[inline]
    fn n_edges(&self) -> usize {
        self.edge_halfedge.len()
    }
    #[inline]
    fn n_faces(&self) -> usize {
        self.faces.len()
    }
    #[inline]
    fn twin(&self, h: usize) -> u32 {
        self.twin[h]
    }
    #[inline]
    fn tail(&self, h: usize) -> u32 {
        self.faces[h / 3][h % 3]
    }
    #[inline]
    fn edge(&self, h: usize) -> u32 {
        self.halfedge_edge[h]
    }
    #[inline]
    fn length(&self, e: usize) -> f64 {
        self.edge_lengths[e]
    }
    #[inline]
    fn outgoing(&self, v: usize) -> u32 {
        self.vertex_halfedge[v]
    }
    #[inline]
    fn signpost(&self, h: usize) -> f64 {
        self.signposts[h]
    }
    #[inline]
    fn angle_sum(&self, v: usize) -> f64 {
        self.angle_sums[v]
    }
    #[inline]
    fn is_forward(&self, h: usize) -> bool {
        self.edge_halfedge[self.halfedge_edge[h] as usize] == h as u32
    }
    #[inline]
    fn edge_halfedge(&self, e: usize) -> u32 {
        self.edge_halfedge[e]
    }
}
