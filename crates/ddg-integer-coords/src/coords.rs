//! The integer representation: connectivity, lengths, normal coordinates, roundabouts and the
//! places of inserted vertices (SPEC §2), with its initialization (§4.1) and the derived
//! quantities of every face (§2.3).

use std::f64::consts::{PI, TAU};
use std::sync::Arc;

use ddg_core::mollify::{Mollified, mollify_edge_lengths};
use ddg_core::topology::NONE;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{Complex, InputSurface, SurfacePoint};

use crate::error::{IntegerError, Place, ValidityItem};
use crate::formulas::{FaceCoordinates, n_minus, n_plus};

/// The value of a roundabout for a halfedge whose tail is an inserted vertex (`V*`), which has
/// no input halfedges around it.
pub const NO_ROUNDABOUT: u32 = NONE;

/// The relative margin of intrinsic mollification the paper uses in its evaluation (Sec. 3.9:
/// `ε = 10⁻⁵ h̄`, `h̄` the mean edge length), for [`InitOptions::mollify`].
pub const PAPER_MOLLIFY_EPSILON: f64 = 1e-5;

/// The local indices of the input mesh (Supp. B; SPEC §2.2, I3): around every input vertex `a`
/// the input edges at `a` are numbered counter-clockwise `0 … deg₀(a) − 1`, starting with the
/// halfedge [`TriMesh::outgoing_halfedge`] (at a boundary vertex the boundary halfedge that
/// leaves it). At a boundary vertex the last index, `deg₀(a) − 1`, is the boundary edge that
/// arrives at `a` (it has no halfedge leaving `a`), and the numbering does not wrap.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalIndices {
    /// Per input vertex: where its directions start in `edges`; one more entry at the end.
    offsets: Vec<u32>,
    /// Per direction: the input edge.
    edges: Vec<u32>,
    /// Per input halfedge: its local index at its tail.
    of_halfedge: Vec<u32>,
    /// Per input edge: its local index at `edge_vertices(e)[0]` and `[1]`.
    of_edge: Vec<[u32; 2]>,
    /// Per input vertex: whether it lies on the boundary.
    boundary: Vec<bool>,
}

impl LocalIndices {
    /// The local indices of `mesh`.
    pub fn new(mesh: &TriMesh) -> LocalIndices {
        let nv = mesh.n_vertices();
        let mut offsets = Vec::with_capacity(nv + 1);
        let mut edges = Vec::with_capacity(2 * mesh.n_edges());
        let mut of_halfedge = vec![NONE; mesh.n_halfedges()];
        let mut of_edge = vec![[NONE; 2]; mesh.n_edges()];
        let mut boundary = Vec::with_capacity(nv);
        let set_edge = |e: usize, v: usize, k: u32, of_edge: &mut Vec<[u32; 2]>| {
            let side = usize::from(mesh.edge_vertices(e)[1] == v);
            of_edge[e][side] = k;
        };
        for v in 0..nv {
            offsets.push(edges.len() as u32);
            let mut k = 0u32;
            let mut last = NONE as usize;
            for h in mesh.outgoing_halfedges(v) {
                let e = mesh.edge(h);
                of_halfedge[h] = k;
                set_edge(e, v, k, &mut of_edge);
                edges.push(e as u32);
                k += 1;
                last = h;
            }
            let is_boundary = mesh.is_boundary_vertex(v);
            if is_boundary {
                // The boundary edge arriving at v: the edge of the halfedge before the last
                // outgoing one.
                let e = mesh.edge(mesh.prev(last));
                set_edge(e, v, k, &mut of_edge);
                edges.push(e as u32);
            }
            boundary.push(is_boundary);
        }
        offsets.push(edges.len() as u32);
        LocalIndices {
            offsets,
            edges,
            of_halfedge,
            of_edge,
            boundary,
        }
    }

    /// `deg₀(a)`: the number of input edges at input vertex `a`.
    #[inline]
    pub fn degree(&self, a: usize) -> u32 {
        self.offsets[a + 1] - self.offsets[a]
    }

    /// Whether input vertex `a` lies on the boundary (its numbering does not wrap).
    #[inline]
    pub fn is_boundary(&self, a: usize) -> bool {
        self.boundary[a]
    }

    /// The modulus of the local indices at `a`: `deg₀(a)` at an interior vertex, `0` (no
    /// reduction) at a boundary vertex, as the oracle's `DEG` argument.
    #[inline]
    pub fn modulus(&self, a: usize) -> u32 {
        if self.boundary[a] { 0 } else { self.degree(a) }
    }

    /// The input edge in direction `k` at `a` (`k < deg₀(a)`).
    #[inline]
    pub fn edge(&self, a: usize, k: u32) -> usize {
        self.edges[(self.offsets[a] + k) as usize] as usize
    }

    /// The local index of input halfedge `h` at its tail.
    #[inline]
    pub fn of_halfedge(&self, h: usize) -> u32 {
        self.of_halfedge[h]
    }

    /// The local index of input edge `e` at its end point `a`, or `None` if `a` is not one.
    pub fn of_edge_at(&self, mesh: &TriMesh, e: usize, a: usize) -> Option<u32> {
        let [x, y] = mesh.edge_vertices(e);
        if x == a {
            Some(self.of_edge[e][0])
        } else if y == a {
            Some(self.of_edge[e][1])
        } else {
            None
        }
    }
}

/// Options of [`IntegerTriangulation::with_options`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InitOptions {
    /// Intrinsic mollification of the input lengths before anything else (SPEC §7, I22): the
    /// relative margin `δ / h̄` of `ddg_core::mollify::mollify_edge_lengths`. `None` (the
    /// default) leaves the lengths as they are; the paper uses [`PAPER_MOLLIFY_EPSILON`].
    pub mollify: Option<f64>,
}

/// What [`IntegerTriangulation::with_options`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct InitReport {
    /// The record of mollification, if it was asked for; with `epsilon > 0` every edge length
    /// was changed by `epsilon` (the input surface is the mollified one).
    pub mollified: Option<Mollified>,
}

/// An intrinsic triangulation in the integer representation of Gillespie, Sharp and Crane
/// (2021): the connectivity ([`Complex`], shared with the signpost representation), edge
/// lengths, one normal coordinate per edge (`i64`, SPEC I2), one roundabout per halfedge whose
/// tail is an input vertex, and the place `q⁰` of every inserted vertex on the input mesh.
///
/// Vertices `0 … n_input_vertices() − 1` are the input vertices `V⁰` with their ids; the others
/// are inserted (`V*`). Read-only so far: the operations (flips, splits, removal, extraction of
/// the geometry) are later work packages; see the crate's README.
#[derive(Debug, Clone, PartialEq)]
pub struct IntegerTriangulation {
    pub(crate) input: Arc<InputSurface>,
    pub(crate) local: Arc<LocalIndices>,
    pub(crate) complex: Complex,
    pub(crate) lengths: Vec<f64>,
    /// Per edge: `n` (SPEC §2.2).
    pub(crate) normal: Vec<i64>,
    /// Per halfedge: `r`, or [`NO_ROUNDABOUT`] for a tail in `V*`.
    pub(crate) roundabouts: Vec<u32>,
    /// Per vertex: its place on the input mesh (`Vertex(v)` for `v ∈ V⁰`, `q⁰` otherwise).
    pub(crate) locations: Vec<SurfacePoint>,
}

/// The stored data of an [`IntegerTriangulation`], for building one from data that was stored
/// or computed elsewhere ([`IntegerTriangulation::from_parts`]) and for taking one apart.
#[derive(Debug, Clone, PartialEq)]
pub struct Parts {
    /// The input mesh.
    pub input: Arc<InputSurface>,
    /// The connectivity of the intrinsic triangulation.
    pub complex: Complex,
    /// One length per edge.
    pub lengths: Vec<f64>,
    /// One normal coordinate per edge.
    pub normal: Vec<i64>,
    /// One roundabout per halfedge ([`NO_ROUNDABOUT`] for a tail in `V*`).
    pub roundabouts: Vec<u32>,
    /// One place on the input mesh per vertex.
    pub locations: Vec<SurfacePoint>,
}

impl IntegerTriangulation {
    /// The intrinsic triangulation that is the input mesh itself (SPEC §2.2, §4.1): same
    /// connectivity and lengths with the same ids, `n ≡ −1`, every roundabout the local index of
    /// its halfedge, no inserted vertex.
    ///
    /// # Errors
    /// [`IntegerError::Intrinsic`] if the geometry was not built for the mesh.
    pub fn new(mesh: &TriMesh, geometry: &Geometry) -> Result<IntegerTriangulation, IntegerError> {
        let input = InputSurface::new(mesh, geometry)?;
        Ok(Self::of_input(Arc::new(input)))
    }

    /// [`new`](Self::new) from the input mesh and one length per edge, with options: optional
    /// intrinsic mollification of the lengths first (SPEC §7; the input surface is then the
    /// mollified one, so tracing and the triangulation agree, I22).
    ///
    /// # Errors
    /// [`IntegerError::Mollify`], [`IntegerError::Geometry`] for lengths that do not form
    /// triangles (after mollification, if asked for).
    pub fn with_options(
        mesh: &TriMesh,
        edge_lengths: &[f64],
        options: &InitOptions,
    ) -> Result<(IntegerTriangulation, InitReport), IntegerError> {
        let (lengths, mollified) = match options.mollify {
            Some(delta) => {
                let m = mollify_edge_lengths(mesh, edge_lengths, delta)?;
                (m.edge_lengths.clone(), Some(m))
            }
            None => (edge_lengths.to_vec(), None),
        };
        let geometry = Geometry::from_edge_lengths(mesh, lengths)?;
        let input = InputSurface::new(mesh, &geometry)?;
        Ok((Self::of_input(Arc::new(input)), InitReport { mollified }))
    }

    /// `T¹ = T⁰` over a given input surface.
    pub(crate) fn of_input(input: Arc<InputSurface>) -> IntegerTriangulation {
        let mesh = input.mesh();
        let local = LocalIndices::new(mesh);
        let roundabouts = (0..mesh.n_halfedges())
            .map(|h| local.of_halfedge(h))
            .collect();
        IntegerTriangulation {
            complex: Complex::from_tri_mesh(mesh),
            lengths: input.geometry().edge_lengths().to_vec(),
            normal: vec![-1; mesh.n_edges()],
            roundabouts,
            locations: (0..mesh.n_vertices() as u32)
                .map(SurfacePoint::Vertex)
                .collect(),
            local: Arc::new(local),
            input,
        }
    }

    /// A triangulation from its stored data. Only the sizes are checked here (one entry per
    /// element of `complex`, the input vertices first); whether the integers are valid is what
    /// [`certify`](Self::certify) decides, so data that fails it can still be built and
    /// examined.
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] (item [`ValidityItem::Ranges`]) naming the array
    /// whose length is wrong.
    pub fn from_parts(parts: Parts) -> Result<IntegerTriangulation, IntegerError> {
        let c = &parts.complex;
        let bad = |what: &str, got: usize, want: usize| {
            IntegerError::invalid(
                ValidityItem::Ranges,
                Place::Vertex(0),
                format!("{what} has {got} entries, expected {want}"),
            )
        };
        let checks = [
            ("lengths", parts.lengths.len(), c.n_edges()),
            ("normal", parts.normal.len(), c.n_edges()),
            ("roundabouts", parts.roundabouts.len(), c.n_halfedges()),
            ("locations", parts.locations.len(), c.n_vertices()),
        ];
        for (what, got, want) in checks {
            if got != want {
                return Err(bad(what, got, want));
            }
        }
        let n_input = parts.input.mesh().n_vertices();
        if c.n_vertices() < n_input {
            return Err(bad("the complex's vertices", c.n_vertices(), n_input));
        }
        c.check()?;
        let local = LocalIndices::new(parts.input.mesh());
        Ok(IntegerTriangulation {
            local: Arc::new(local),
            input: parts.input,
            complex: parts.complex,
            lengths: parts.lengths,
            normal: parts.normal,
            roundabouts: parts.roundabouts,
            locations: parts.locations,
        })
    }

    /// The stored data, for [`from_parts`](Self::from_parts).
    pub fn into_parts(self) -> Parts {
        Parts {
            input: self.input,
            complex: self.complex,
            lengths: self.lengths,
            normal: self.normal,
            roundabouts: self.roundabouts,
            locations: self.locations,
        }
    }

    // ----- accessors -----

    /// The connectivity (SPEC §2.1).
    #[inline]
    pub fn complex(&self) -> &Complex {
        &self.complex
    }

    /// The input mesh `T⁰`.
    #[inline]
    pub fn input(&self) -> &InputSurface {
        &self.input
    }

    /// The input mesh, shared.
    #[inline]
    pub fn input_arc(&self) -> &Arc<InputSurface> {
        &self.input
    }

    /// The local indices of the input mesh (SPEC I3).
    #[inline]
    pub fn local_indices(&self) -> &LocalIndices {
        &self.local
    }

    /// Number of input vertices; they are the vertices `0 … n_input_vertices() − 1`.
    #[inline]
    pub fn n_input_vertices(&self) -> usize {
        self.input.mesh().n_vertices()
    }

    /// Whether vertex `v` is an input vertex (`V⁰`) rather than an inserted one (`V*`).
    #[inline]
    pub fn is_input_vertex(&self, v: usize) -> bool {
        v < self.n_input_vertices()
    }

    /// The edge lengths, indexed by edge id.
    #[inline]
    pub fn edge_lengths(&self) -> &[f64] {
        &self.lengths
    }

    /// The normal coordinates `n`, indexed by edge id.
    #[inline]
    pub fn normal_coordinates(&self) -> &[i64] {
        &self.normal
    }

    /// The normal coordinate of edge `e`: the number of transversal crossings with input edges,
    /// or `−1` if it runs along an input edge.
    #[inline]
    pub fn normal_coordinate(&self, e: usize) -> i64 {
        self.normal[e]
    }

    /// The normal coordinate of the edge of halfedge `h`.
    #[inline]
    pub fn halfedge_normal(&self, h: usize) -> i64 {
        self.normal[self.complex.edge(h)]
    }

    /// The roundabouts, indexed by halfedge ([`NO_ROUNDABOUT`] for a tail in `V*`).
    #[inline]
    pub fn roundabouts(&self) -> &[u32] {
        &self.roundabouts
    }

    /// The roundabout of halfedge `h`: the local index at its tail of the first input halfedge
    /// at or after it counter-clockwise; `None` for a tail in `V*`.
    #[inline]
    pub fn roundabout(&self, h: usize) -> Option<u32> {
        let r = self.roundabouts[h];
        (r != NO_ROUNDABOUT).then_some(r)
    }

    /// The place of vertex `v` on the input mesh (`q⁰` for an inserted vertex).
    #[inline]
    pub fn vertex_location(&self, v: usize) -> SurfacePoint {
        self.locations[v]
    }

    /// The places of all vertices.
    #[inline]
    pub fn vertex_locations(&self) -> &[SurfacePoint] {
        &self.locations
    }

    /// The angle sum `Θ` of vertex `v`: that of the input mesh at an input vertex, `2π` at an
    /// inserted interior vertex, `π` at one inserted on the boundary.
    pub fn vertex_angle_sum(&self, v: usize) -> f64 {
        if self.is_input_vertex(v) {
            self.input.vertex_angle_sum(v)
        } else if self.complex.is_boundary_vertex(v) {
            PI
        } else {
            TAU
        }
    }

    /// The lengths of the three sides of face `f`, by slot.
    pub fn face_lengths(&self, f: usize) -> [f64; 3] {
        [0, 1, 2].map(|k| self.lengths[self.complex.edge(3 * f + k)])
    }

    // ----- derived quantities (SPEC §2.3) -----

    /// The normal coordinates of the three sides of face `f`, by slot (side `k` is halfedge
    /// `3f + k`).
    #[inline]
    pub fn face_normal(&self, f: usize) -> [i64; 3] {
        [0, 1, 2].map(|k| self.normal[self.complex.edge(3 * f + k)])
    }

    /// The corner coordinates and emanating counts of face `f` (Eq. 1, 2).
    ///
    /// # Errors
    /// [`IntegerError::InvalidCoordinates`] (item [`ValidityItem::FaceParity`]) naming the face
    /// and the corner whose numerator of Eq. 1 is odd.
    ///
    /// [formula: ic-corner-coordinates]
    pub fn face_coordinates(&self, f: usize) -> Result<FaceCoordinates, IntegerError> {
        face_coordinates_of(self.face_normal(f), f)
    }

    /// `c` at corner `h` (the corner of face `h / 3` at the tail of halfedge `h`): the number of
    /// arcs cutting it (Eq. 1).
    ///
    /// # Errors
    /// As [`face_coordinates`](Self::face_coordinates).
    ///
    /// [formula: ic-corner-coordinates]
    pub fn corner_coordinate(&self, h: usize) -> Result<i64, IntegerError> {
        Ok(self.face_coordinates(h / 3)?.c[h % 3])
    }

    /// `e` at corner `h`: the number of arcs from the tail of `h` to the opposite side of its
    /// face (Eq. 2).
    ///
    /// [formula: ic-emanating]
    pub fn emanating_count(&self, h: usize) -> i64 {
        crate::formulas::face_emanating(self.face_normal(h / 3), h % 3)
    }

    /// `n⁺` of the edge of halfedge `h`.
    ///
    /// [formula: ic-normal-sign]
    #[inline]
    pub fn n_plus(&self, h: usize) -> i64 {
        n_plus(self.halfedge_normal(h))
    }

    /// `n⁻` of the edge of halfedge `h`.
    ///
    /// [formula: ic-normal-sign]
    #[inline]
    pub fn n_minus(&self, h: usize) -> i64 {
        n_minus(self.halfedge_normal(h))
    }

    /// `Σ n⁺` over the edges: the number of crossings of the input edges with the edges of the
    /// triangulation (`|V^S| = |V¹| + Σ n⁺`).
    pub fn total_crossings(&self) -> u128 {
        self.normal.iter().map(|&n| n_plus(n) as u128).sum()
    }
}

/// [`FaceCoordinates::of`] with the error of this crate naming face `f`.
pub(crate) fn face_coordinates_of(n: [i64; 3], f: usize) -> Result<FaceCoordinates, IntegerError> {
    FaceCoordinates::of(n).map_err(|k| {
        IntegerError::invalid(
            ValidityItem::FaceParity,
            Place::Face(f as u32),
            format!(
                "the numerator of Eq. 1 at corner {k} is odd (side counts {} {} {})",
                n[0], n[1], n[2]
            ),
        )
    })
}
