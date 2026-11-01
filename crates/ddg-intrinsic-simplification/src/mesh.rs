//! The coarse mesh: a Δ-complex with edge lengths, angular coordinates and tracked points
//! (SPEC §2), its construction (§3.1), the derived quantities (§2.3) and its certificate (§2.4).

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::{corner_angle, corner_cotan, triangle_area};
use ddg_core::topology::NONE;
use ddg_core::{Geometry, SparseMatrix, TriMesh};
use ddg_intrinsic::formulas::delaunay_angle_excess;
use ddg_intrinsic::{Complex, DelaunayCriterion, FlipLimit, FlipShapeTest};

use crate::error::IemError;
use crate::formulas::wrap_angle;
use crate::undo::{Undo, UndoLog};

/// The tolerance of the certificate's angle checks (SPEC §2.4 item 3), in radians.
pub const ANGLE_CHECK_TOLERANCE: f64 = 1e-9;

/// Options fixed at construction and used by every operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoarseMeshOptions {
    /// The shape test of every flip (flips to degree three or two, Springborn's fallback,
    /// Delaunay flips) on top of strict convexity (SPEC §3.2). Default: IT-SPEC N5's.
    pub flip_test: FlipShapeTest,
    /// When an edge is not Delaunay (SPEC §3.3; IT-SPEC's default `AngleSum { 1e-10 }`).
    pub delaunay: DelaunayCriterion,
    /// The limit on the number of flips of one Delaunay pass.
    pub max_delaunay_flips: FlipLimit,
    /// The slack of the inside test of a tracked point (SPEC §12.1, S16), relative to the
    /// mean edge length of the region the point is re-expressed in: a signed distance to a
    /// side down to `−point_slack · L` is clamped to the side and counted; further out is an
    /// internal error.
    pub point_slack: f64,
    /// The flattening's stopping tolerance, relative to the target angle sum (SPEC §3.4).
    pub flatten_tolerance: f64,
    /// The flattening's iteration cap (Newton and bisection steps per attempt).
    pub flatten_max_iterations: usize,
    /// Whether a committed removal flips the edges among the neighbours to Delaunay
    /// (SPEC §3.3, S15; paper Alg. 1 line 12).
    pub delaunay_after_removal: bool,
}

impl Default for CoarseMeshOptions {
    fn default() -> Self {
        CoarseMeshOptions {
            flip_test: FlipShapeTest::N5,
            delaunay: DelaunayCriterion::default(),
            max_delaunay_flips: FlipLimit::default(),
            point_slack: 1e-12,
            flatten_tolerance: 1e-12,
            flatten_max_iterations: 100,
            delaunay_after_removal: true,
        }
    }
}

/// A tracked point on the input surface, given at construction: a face of the input mesh and
/// barycentric coordinates w.r.t. its corners (`faces[face][k]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InputPoint {
    /// The input face.
    pub face: u32,
    /// Barycentric coordinates, non-negative, summing to one (to `1e-12`).
    pub bary: [f64; 3],
}

/// Where a tracked point is (internal form; faces are complex face ids).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Location {
    /// The point is the fine vertex with this fine id, still a vertex of the coarse mesh.
    Vertex(u32),
    /// In face `face`, with barycentric coordinates w.r.t. its slots.
    Face { face: u32, bary: [f64; 3] },
}

/// Where a tracked point is on the coarse mesh.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PointLocation {
    /// It is (the fine vertex with) this fine id, a vertex of the coarse mesh.
    Vertex(u32),
    /// It lies in coarse face `face` (an id of [`CoarseMesh::faces`]) with barycentric
    /// coordinates `bary` w.r.t. the corners `corners` (fine ids, in the face's order; a
    /// corner may repeat in a Δ-complex).
    Face {
        /// The coarse face.
        face: u32,
        /// Its corners (fine ids).
        corners: [u32; 3],
        /// Barycentric coordinates, `≥ 0`, summing to one.
        bary: [f64; 3],
    },
}

/// What [`CoarseMesh::certify`] measured.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Certificate {
    /// The largest error of an angular coordinate difference against `s ×` the corner angle,
    /// in radians (item 3).
    pub max_angle_error: f64,
    /// `|Σ K + Σ κ − 2πχ|` (item 4).
    pub gauss_bonnet_error: f64,
    /// The bound item 4 holds it to, `1e-9 (1 + Σ|K| + Σ|κ|)`.
    pub gauss_bonnet_bound: f64,
    /// The largest `|Σ b − 1|` of a tracked point (item 5).
    pub max_bary_sum_error: f64,
    /// The largest opposite-angle sum minus `π` over interior edges (item 7, measured always;
    /// checked only when asked); `−∞` without interior edges.
    pub max_delaunay_excess: f64,
    /// Interior edges that fail the Delaunay test (item 7, counted always).
    pub n_not_delaunay: usize,
}

/// A coarse intrinsic triangulation of an input mesh, with the fine vertices and other points
/// tracked on it (SPEC §2). See the [crate documentation](crate).
#[derive(Debug, Clone)]
pub struct CoarseMesh {
    pub(crate) complex: Complex,
    /// Per edge (complex id): its length.
    pub(crate) length: Vec<f64>,
    /// Per halfedge: the normalized angle at its tail (S4).
    pub(crate) phi: Vec<f64>,
    /// Per complex vertex: its reference halfedge (S4).
    pub(crate) reference: Vec<u32>,
    /// Per complex vertex: the sum of its corner angles, as of the last change of its lengths.
    pub(crate) angle_sum: Vec<f64>,
    /// Per complex vertex: its fine id.
    pub(crate) fine_of: Vec<u32>,
    /// Per fine id: its complex vertex, or `NONE` once removed.
    pub(crate) complex_of: Vec<u32>,
    /// Per fine id: whether it is fixed.
    pub(crate) fixed: Vec<bool>,
    /// The number of input vertices (points `0..n_input` are the fine vertices).
    pub(crate) n_input: u32,
    /// Per point: where it is.
    pub(crate) points: Vec<Location>,
    /// Per point: its index in its face's list (`NONE` for a vertex).
    pub(crate) point_slot: Vec<u32>,
    /// Per complex face: the points in it.
    pub(crate) face_points: Vec<Vec<u32>>,
    pub(crate) options: CoarseMeshOptions,
    /// Points clamped onto a side so far (S16).
    pub(crate) n_clamped: u64,
    pub(crate) undo: UndoLog,
    pub(crate) pending: Option<crate::removal::Pending>,
    /// Scratch marks (epoch stamps) per vertex and per edge.
    pub(crate) vertex_mark: Vec<u32>,
    pub(crate) edge_mark: Vec<u32>,
    pub(crate) epoch: u32,
    /// Scratch: a fan with its corner angles.
    pub(crate) scratch_fan: Vec<(u32, f64)>,
}

/// The corner angle at the tail of halfedge `h` of a complex with these lengths, `None` if its
/// face is not a valid triangle.
#[inline]
pub(crate) fn corner_at(complex: &Complex, length: &[f64], h: usize) -> Option<f64> {
    let l = |g: usize| length[complex.edge(g)];
    corner_angle(l(complex.next(h)), l(h), l(complex.prev(h)))
}

/// The lengths of face `f`, `l[k]` from slot `k` to slot `k + 1`.
#[inline]
pub(crate) fn face_lengths(complex: &Complex, length: &[f64], f: usize) -> [f64; 3] {
    [0, 1, 2].map(|k| length[complex.edge(3 * f + k)])
}

/// The next outgoing halfedge counter-clockwise around its tail, `None` past the boundary.
#[inline]
pub(crate) fn ccw(complex: &Complex, h: usize) -> Option<usize> {
    complex.twin(complex.prev(h))
}

impl CoarseMesh {
    /// The coarse mesh of an input mesh with its lengths, before any simplification: the
    /// complex of `mesh` (same vertex, edge and face ids), the lengths of `geometry`, angular
    /// coordinates from the reference halfedges [`Complex::outgoing_halfedge`] (S4), every
    /// fine vertex a vertex, and the extra `points` located in their input faces (point id
    /// `n_vertices + index`). No vertex is fixed. Nothing is flipped (Delaunay flipping is
    /// [`flip_to_delaunay`](Self::flip_to_delaunay)).
    ///
    /// `mesh` is a validated manifold oriented mesh ([`TriMesh::new`]) and `geometry` validated
    /// lengths for it ([`Geometry::from_positions`], [`Geometry::from_edge_lengths`]: every
    /// face satisfies the strict triangle inequality); [`from_positions`](Self::from_positions)
    /// and [`from_lengths`](Self::from_lengths) do both.
    ///
    /// # Errors
    /// [`IemError::Geometry`] if `geometry` is not one of `mesh`; [`IemError::InvalidPoint`]
    /// for a point whose face is out of range or whose coordinates are not finite,
    /// non-negative and summing to one within `1e-12`.
    pub fn from_mesh(
        mesh: &TriMesh,
        geometry: &Geometry,
        points: &[InputPoint],
        options: CoarseMeshOptions,
    ) -> Result<CoarseMesh, IemError> {
        geometry.check_mesh(mesh)?;
        let complex = Complex::from_tri_mesh(mesh);
        let n = complex.n_vertices();
        let length = geometry.edge_lengths().to_vec();
        let mut face_points: Vec<Vec<u32>> = vec![Vec::new(); complex.n_faces()];
        let mut locations: Vec<Location> = (0..n as u32).map(Location::Vertex).collect();
        let mut point_slot = vec![NONE; n];
        for (index, p) in points.iter().enumerate() {
            let bad = |what: String| Err(IemError::InvalidPoint { index, what });
            if p.face as usize >= complex.n_faces() {
                return bad(format!(
                    "face {} is out of range (the mesh has {} faces)",
                    p.face,
                    complex.n_faces()
                ));
            }
            if !p.bary.iter().all(|b| b.is_finite() && *b >= 0.0) {
                return bad(format!(
                    "barycentric coordinates {:?} are not finite and non-negative",
                    p.bary
                ));
            }
            let sum = p.bary[0] + p.bary[1] + p.bary[2];
            if (sum - 1.0).abs() > 1e-12 {
                return bad(format!(
                    "barycentric coordinates {:?} sum to {sum}, not one",
                    p.bary
                ));
            }
            let id = (n + index) as u32;
            point_slot.push(face_points[p.face as usize].len() as u32);
            face_points[p.face as usize].push(id);
            locations.push(Location::Face {
                face: p.face,
                bary: p.bary,
            });
        }
        let mut mesh = CoarseMesh {
            reference: (0..n)
                .map(|v| complex.outgoing_halfedge(v) as u32)
                .collect(),
            phi: vec![0.0; complex.n_halfedges()],
            angle_sum: vec![0.0; n],
            fine_of: (0..n as u32).collect(),
            complex_of: (0..n as u32).collect(),
            fixed: vec![false; n],
            n_input: n as u32,
            points: locations,
            point_slot,
            face_points,
            options,
            n_clamped: 0,
            undo: UndoLog::default(),
            pending: None,
            vertex_mark: vec![0; n],
            edge_mark: vec![0; complex.n_edges()],
            epoch: 0,
            scratch_fan: Vec::new(),
            complex,
            length,
        };
        for v in 0..n {
            mesh.refresh_vertex(v);
        }
        Ok(mesh)
    }

    /// [`from_mesh`](Self::from_mesh) from a face list and vertex positions.
    ///
    /// # Errors
    /// [`IemError::Mesh`] (not a manifold oriented mesh), [`IemError::Geometry`] (a face
    /// without positive area, a non-finite position), and those of `from_mesh`.
    pub fn from_positions(
        positions: &[[f64; 3]],
        faces: &[[u32; 3]],
        points: &[InputPoint],
        options: CoarseMeshOptions,
    ) -> Result<CoarseMesh, IemError> {
        let mesh = TriMesh::new(positions.len(), faces.to_vec())?;
        let geometry = Geometry::from_positions(&mesh, positions)?;
        CoarseMesh::from_mesh(&mesh, &geometry, points, options)
    }

    /// [`from_mesh`](Self::from_mesh) from a mesh and one length per edge of it
    /// ([`TriMesh::edges`] order).
    ///
    /// # Errors
    /// [`IemError::Geometry`] (a face violating the strict triangle inequality, names the
    /// face), and those of `from_mesh`.
    pub fn from_lengths(
        mesh: &TriMesh,
        lengths: Vec<f64>,
        points: &[InputPoint],
        options: CoarseMeshOptions,
    ) -> Result<CoarseMesh, IemError> {
        let geometry = Geometry::from_edge_lengths(mesh, lengths)?;
        CoarseMesh::from_mesh(mesh, &geometry, points, options)
    }

    // ----- accessors -----

    /// The options.
    pub fn options(&self) -> &CoarseMeshOptions {
        &self.options
    }

    /// The connectivity (internal ids; vertex `v` of the complex is the fine vertex
    /// [`fine_vertex(v)`](Self::fine_vertex)).
    pub fn complex(&self) -> &Complex {
        &self.complex
    }

    /// The length of every edge of [`complex`](Self::complex).
    pub fn edge_lengths(&self) -> &[f64] {
        &self.length
    }

    /// The angular coordinate of every halfedge of [`complex`](Self::complex) at its tail
    /// (S4).
    pub fn angular_coordinates(&self) -> &[f64] {
        &self.phi
    }

    /// The reference halfedge of complex vertex `v` (S4).
    pub fn reference_halfedge(&self, v: usize) -> usize {
        self.reference[v] as usize
    }

    /// The number of fine (input) vertices.
    pub fn n_fine_vertices(&self) -> usize {
        self.n_input as usize
    }

    /// The number of vertices of the coarse mesh.
    pub fn n_vertices(&self) -> usize {
        self.complex.n_vertices()
    }

    /// The number of faces of the coarse mesh.
    pub fn n_faces(&self) -> usize {
        self.complex.n_faces()
    }

    /// The number of edges of the coarse mesh.
    pub fn n_edges(&self) -> usize {
        self.complex.n_edges()
    }

    /// The fine id of complex vertex `v`.
    pub fn fine_vertex(&self, v: usize) -> u32 {
        self.fine_of[v]
    }

    /// The complex vertex of fine vertex `fine`, `None` if it was removed (or is out of
    /// range).
    pub fn complex_vertex(&self, fine: u32) -> Option<usize> {
        match self.complex_of.get(fine as usize) {
            Some(&v) if v != NONE => Some(v as usize),
            _ => None,
        }
    }

    /// Whether fine vertex `fine` is a vertex of the coarse mesh.
    pub fn is_vertex(&self, fine: u32) -> bool {
        self.complex_vertex(fine).is_some()
    }

    /// The vertices of the coarse mesh as fine ids, in the order of the complex.
    pub fn vertices(&self) -> &[u32] {
        &self.fine_of
    }

    /// The faces of the coarse mesh with fine vertex ids (in complex face order; a Δ-face may
    /// repeat a vertex).
    pub fn faces(&self) -> Vec<[u32; 3]> {
        self.complex
            .faces()
            .iter()
            .map(|f| f.map(|v| self.fine_of[v as usize]))
            .collect()
    }

    /// Fixes or frees fine vertex `fine`: a fixed vertex is never removed (paper Sec. 8.6.3).
    ///
    /// # Errors
    /// [`IemError::IndexOutOfRange`].
    pub fn set_fixed(&mut self, fine: u32, fixed: bool) -> Result<(), IemError> {
        let bound = self.fixed.len();
        match self.fixed.get_mut(fine as usize) {
            Some(slot) => {
                *slot = fixed;
                Ok(())
            }
            None => Err(IemError::IndexOutOfRange {
                what: "fine vertex",
                index: fine as usize,
                bound,
            }),
        }
    }

    /// Whether fine vertex `fine` is fixed.
    pub fn is_fixed(&self, fine: u32) -> bool {
        self.fixed.get(fine as usize).copied().unwrap_or(false)
    }

    /// The number of tracked points: the fine vertices, then the extra points.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }

    /// Where tracked point `p` is.
    ///
    /// # Errors
    /// [`IemError::IndexOutOfRange`].
    pub fn point_location(&self, p: usize) -> Result<PointLocation, IemError> {
        match self.points.get(p) {
            None => Err(IemError::IndexOutOfRange {
                what: "point",
                index: p,
                bound: self.points.len(),
            }),
            Some(&Location::Vertex(v)) => Ok(PointLocation::Vertex(v)),
            Some(&Location::Face { face, bary }) => Ok(PointLocation::Face {
                face,
                corners: self.complex.faces()[face as usize].map(|v| self.fine_of[v as usize]),
                bary,
            }),
        }
    }

    /// The points in face `f` (internal order).
    pub fn points_in_face(&self, f: usize) -> &[u32] {
        &self.face_points[f]
    }

    /// The number of points clamped onto a side so far (S16).
    pub fn n_clamped_points(&self) -> u64 {
        self.n_clamped
    }

    // ----- derived quantities (SPEC §2.3) -----

    /// The corner angle at the tail of halfedge `h` (law of cosines, Kahan's form), `NaN` if
    /// its face is not a valid triangle (never the case for a certified mesh).
    pub fn corner_angle(&self, h: usize) -> f64 {
        corner_at(&self.complex, &self.length, h).unwrap_or(f64::NAN)
    }

    /// The area of face `f` (Heron, Kahan's form), `NaN` if it is not a valid triangle.
    pub fn face_area(&self, f: usize) -> f64 {
        let l = face_lengths(&self.complex, &self.length, f);
        triangle_area(l[0], l[1], l[2]).unwrap_or(f64::NAN)
    }

    /// The total area.
    pub fn total_area(&self) -> f64 {
        (0..self.n_faces()).map(|f| self.face_area(f)).sum()
    }

    /// The angle sum `Θ` of complex vertex `v`, recomputed from the lengths (counter-clockwise
    /// from its reference halfedge).
    pub fn angle_sum(&self, v: usize) -> f64 {
        self.fan(v)
            .map(|h| corner_at(&self.complex, &self.length, h).unwrap_or(f64::NAN))
            .sum()
    }

    /// The discrete curvature of complex vertex `v`: `K = 2π − Θ` inside (paper Eq. 1),
    /// `κ = π − Θ` on the boundary (Eq. 2), from the cached angle sum (bitwise what the
    /// operations use).
    pub fn curvature(&self, v: usize) -> f64 {
        let target = if self.complex.is_boundary_vertex(v) {
            PI
        } else {
            TAU
        };
        target - self.angle_sum[v]
    }

    /// The cotan Laplacian (positive semidefinite, `CLAUDE.md` conventions) over the complex
    /// vertices: `L_ij = −½ Σ cot` of the angles opposite the halfedges from `i` to `j`
    /// (self-edges contribute nothing), `L_ii = −Σ_{j≠i} L_ij`.
    ///
    /// [formula: cotan-laplacian]
    pub fn cotan_laplacian(&self) -> SparseMatrix {
        let c = &self.complex;
        let n = c.n_vertices();
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(6 * c.n_faces() + n);
        let mut diagonal = vec![0.0; n];
        for h in 0..c.n_halfedges() {
            let (i, j) = (c.from_vertex(h), c.to_vertex(h));
            if i == j {
                continue;
            }
            let l = face_lengths(c, &self.length, h / 3);
            let k = h % 3;
            let cot = corner_cotan(l[k], l[(k + 1) % 3], l[(k + 2) % 3]).unwrap_or(f64::NAN);
            let w = 0.5 * cot;
            triplets.push((i as u32, j as u32, -w));
            triplets.push((j as u32, i as u32, -w));
            diagonal[i] += w;
            diagonal[j] += w;
        }
        for (v, d) in diagonal.into_iter().enumerate() {
            triplets.push((v as u32, v as u32, d));
        }
        SparseMatrix::from_triplets(n, n, &triplets)
            .expect("vertex indices of the complex are in range")
    }

    /// The lumped masses `⅓ Σ A_f` over the complex vertices (a face that uses a vertex twice
    /// counts twice for it).
    ///
    /// [formula: lumped-mass]
    pub fn lumped_mass(&self) -> Vec<f64> {
        let mut mass = vec![0.0; self.n_vertices()];
        for (f, face) in self.complex.faces().iter().enumerate() {
            let third = self.face_area(f) / 3.0;
            for &v in face {
                mass[v as usize] += third;
            }
        }
        mass
    }

    // ----- fans and bookkeeping -----

    /// The outgoing halfedges of complex vertex `v`, counter-clockwise from its reference
    /// halfedge.
    pub(crate) fn fan(&self, v: usize) -> Fan<'_> {
        let start = self.reference[v] as usize;
        Fan {
            complex: &self.complex,
            start,
            current: Some(start),
        }
    }

    /// Recomputes the angle sum of complex vertex `v` and its angular coordinates from its
    /// reference halfedge (S4): the reference keeps its coordinate (`0` on the boundary), the
    /// others follow counter-clockwise by `s ×` the corner angles. Writes are recorded.
    pub(crate) fn refresh_vertex(&mut self, v: usize) {
        let boundary = self.complex.is_boundary_vertex(v);
        // One pass for the corner angles, summed in fan order (bitwise `angle_sum(v)`).
        let mut fan = std::mem::take(&mut self.scratch_fan);
        fan.clear();
        let mut theta = 0.0;
        for h in self.fan(v) {
            let corner = corner_at(&self.complex, &self.length, h).unwrap_or(f64::NAN);
            theta += corner;
            fan.push((h as u32, corner));
        }
        self.set_angle_sum(v, theta);
        let s = ddg_intrinsic::formulas::angle_scale(theta, boundary);
        let start = fan[0].0 as usize;
        if boundary && self.phi[start] != 0.0 {
            self.set_phi(start, 0.0);
        }
        let mut phi = self.phi[start];
        for k in 1..fan.len() {
            phi = wrap_angle(phi + s * fan[k - 1].1, boundary);
            self.set_phi(fan[k].0 as usize, phi);
        }
        self.scratch_fan = fan;
    }

    /// A new epoch for the scratch marks.
    pub(crate) fn next_epoch(&mut self) -> u32 {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.vertex_mark.iter_mut().for_each(|m| *m = 0);
            self.edge_mark.iter_mut().for_each(|m| *m = 0);
            self.epoch = 1;
        }
        if self.vertex_mark.len() < self.complex.n_vertices() {
            self.vertex_mark.resize(self.complex.n_vertices(), 0);
        }
        if self.edge_mark.len() < self.complex.n_edges() {
            self.edge_mark.resize(self.complex.n_edges(), 0);
        }
        self.epoch
    }

    // ----- recorded writes -----

    #[inline]
    pub(crate) fn set_length(&mut self, e: usize, value: f64) {
        self.undo.push(Undo::Length(e as u32, self.length[e]));
        self.length[e] = value;
    }

    #[inline]
    pub(crate) fn set_phi(&mut self, h: usize, value: f64) {
        self.undo.push(Undo::Phi(h as u32, self.phi[h]));
        self.phi[h] = value;
    }

    #[inline]
    pub(crate) fn set_reference(&mut self, v: usize, h: u32) {
        if self.reference[v] != h {
            self.undo.push(Undo::Reference(v as u32, self.reference[v]));
            self.reference[v] = h;
        }
    }

    #[inline]
    pub(crate) fn set_angle_sum(&mut self, v: usize, value: f64) {
        self.undo.push(Undo::AngleSum(v as u32, self.angle_sum[v]));
        self.angle_sum[v] = value;
    }

    /// Sets the coordinates of point `p`, which stays in its face.
    #[inline]
    pub(crate) fn set_bary(&mut self, p: u32, bary: [f64; 3]) {
        if let Location::Face { face, bary: old } = self.points[p as usize] {
            self.undo.push(Undo::Bary(p, old));
            self.points[p as usize] = Location::Face { face, bary };
        }
    }

    /// Moves point `p` (in a face) to the end of face `to`'s list with coordinates `bary`.
    pub(crate) fn move_point(&mut self, p: u32, to: u32, bary: [f64; 3]) {
        let Location::Face { face, bary: old } = self.points[p as usize] else {
            return;
        };
        let slot = self.point_slot[p as usize];
        self.undo.push(Undo::Moved {
            point: p,
            face,
            slot,
            bary: old,
        });
        let list = &mut self.face_points[face as usize];
        list.swap_remove(slot as usize);
        if let Some(&q) = list.get(slot as usize) {
            self.point_slot[q as usize] = slot;
        }
        self.point_slot[p as usize] = self.face_points[to as usize].len() as u32;
        self.face_points[to as usize].push(p);
        self.points[p as usize] = Location::Face { face: to, bary };
    }

    /// Replays the undo log backwards (SPEC S5): every recorded array is restored bitwise.
    pub(crate) fn replay_undo(&mut self) {
        let entries = std::mem::take(&mut self.undo.entries);
        for entry in entries.iter().rev() {
            match *entry {
                Undo::Length(e, v) => self.length[e as usize] = v,
                Undo::Phi(h, v) => self.phi[h as usize] = v,
                Undo::Reference(v, h) => self.reference[v as usize] = h,
                Undo::AngleSum(v, x) => self.angle_sum[v as usize] = x,
                Undo::Flip(i) => {
                    let u = self.undo.flips[i as usize];
                    self.complex.undo_flip(&u);
                }
                Undo::Bary(p, bary) => {
                    if let Location::Face { face, .. } = self.points[p as usize] {
                        self.points[p as usize] = Location::Face { face, bary };
                    }
                }
                Undo::Moved {
                    point,
                    face,
                    slot,
                    bary,
                } => {
                    // The point is the last of its current list (later moves were undone).
                    if let Location::Face { face: now, .. } = self.points[point as usize] {
                        self.face_points[now as usize].pop();
                    }
                    let list = &mut self.face_points[face as usize];
                    if slot as usize == list.len() {
                        list.push(point);
                    } else {
                        let q = list[slot as usize];
                        list.push(q);
                        list[slot as usize] = point;
                        self.point_slot[q as usize] = (list.len() - 1) as u32;
                    }
                    self.point_slot[point as usize] = slot;
                    self.points[point as usize] = Location::Face { face, bary };
                }
            }
        }
        self.undo.entries = entries;
        self.undo.entries.clear();
        self.undo.flips.clear();
    }

    // ----- certificate (SPEC §2.4) -----

    /// Checks the invariants of SPEC §2.4: (1) the complex and the fine/complex maps;
    /// (2) the strict triangle inequality of every face; (3) the angular coordinates against
    /// the corner angles (`1e-9` rad), the boundary halfedge at `0`; (4) Gauss–Bonnet;
    /// (5) every tracked point in exactly one place with valid coordinates; (6) all values
    /// finite; (7) only if `delaunay`: the Delaunay test of §3.3 at every interior edge.
    /// `O(|F| + #points)`.
    ///
    /// # Errors
    /// [`IemError::Certify`] with the first violation found.
    pub fn certify(&self, delaunay: bool) -> Result<Certificate, IemError> {
        let fail = |item: u8, what: String| Err(IemError::Certify { item, what });
        let c = &self.complex;
        let mut cert = Certificate {
            max_delaunay_excess: f64::NEG_INFINITY,
            ..Certificate::default()
        };
        // 1.
        if let Err(e) = c.check() {
            return fail(1, e.to_string());
        }
        let n = c.n_vertices();
        if self.fine_of.len() != n
            || self.reference.len() != n
            || self.angle_sum.len() != n
            || self.length.len() != c.n_edges()
            || self.phi.len() != c.n_halfedges()
            || self.face_points.len() != c.n_faces()
        {
            return fail(1, "array lengths do not match the complex".into());
        }
        for (v, &fine) in self.fine_of.iter().enumerate() {
            if self.complex_of.get(fine as usize) != Some(&(v as u32)) {
                return fail(
                    1,
                    format!("complex vertex {v} has fine id {fine}, which maps elsewhere"),
                );
            }
        }
        let live = self.complex_of.iter().filter(|&&v| v != NONE).count();
        if live != n {
            return fail(
                1,
                format!("{live} fine ids map to vertices, but the complex has {n}"),
            );
        }
        // 2. and 6.
        for (e, &l) in self.length.iter().enumerate() {
            if !(l.is_finite() && l > 0.0) {
                return fail(6, format!("edge {e} has length {l}"));
            }
        }
        for f in 0..c.n_faces() {
            let l = face_lengths(c, &self.length, f);
            if triangle_area(l[0], l[1], l[2]).is_none() {
                return fail(
                    2,
                    format!(
                        "face {f} {:?} has lengths {l:?}, which violate the strict triangle \
                         inequality",
                        self.faces_fine(f)
                    ),
                );
            }
        }
        // 3.
        for v in 0..n {
            let boundary = c.is_boundary_vertex(v);
            let r = self.reference[v] as usize;
            if r >= c.n_halfedges() || c.from_vertex(r) != v {
                return fail(
                    3,
                    format!("vertex {v}: reference halfedge {r} does not leave it"),
                );
            }
            if boundary && c.twin(r).is_some() {
                return fail(
                    3,
                    format!("boundary vertex {v}: reference {r} is not the boundary halfedge"),
                );
            }
            if boundary && self.phi[r] != 0.0 {
                return fail(
                    3,
                    format!("boundary vertex {v}: boundary halfedge at {}", self.phi[r]),
                );
            }
            let theta = self.angle_sum(v);
            let s = ddg_intrinsic::formulas::angle_scale(theta, boundary);
            let fan: Vec<usize> = self.fan(v).collect();
            for (k, &h) in fan.iter().enumerate() {
                if !self.phi[h].is_finite() {
                    return fail(6, format!("halfedge {h} has angle {}", self.phi[h]));
                }
                let corner = corner_at(c, &self.length, h).unwrap_or(f64::NAN);
                let expected = if k + 1 < fan.len() {
                    self.phi[fan[k + 1]]
                } else if boundary {
                    PI
                } else {
                    self.phi[r]
                };
                let got = self.phi[h] + s * corner;
                let err = if boundary {
                    (got - expected).abs()
                } else {
                    let d = (got - expected).rem_euclid(TAU);
                    d.min(TAU - d)
                };
                cert.max_angle_error = cert.max_angle_error.max(err);
                if err.is_nan() || err > ANGLE_CHECK_TOLERANCE {
                    return fail(
                        3,
                        format!(
                            "vertex {v} (fine {}): halfedge {h} at {} plus s × corner {corner} \
                             misses the next coordinate {expected} by {err:e}",
                            self.fine_of[v], self.phi[h]
                        ),
                    );
                }
                if boundary && !(0.0..=PI + ANGLE_CHECK_TOLERANCE).contains(&self.phi[h])
                    || !boundary && !(0.0..TAU).contains(&self.phi[h])
                {
                    return fail(
                        3,
                        format!("vertex {v}: halfedge {h} at {} out of range", self.phi[h]),
                    );
                }
            }
        }
        // 4.
        let mut sum = 0.0;
        let mut abs = 0.0;
        for v in 0..n {
            let target = if c.is_boundary_vertex(v) { PI } else { TAU };
            let k = target - self.angle_sum(v);
            sum += k;
            abs += k.abs();
        }
        let chi = c.euler_characteristic() as f64;
        cert.gauss_bonnet_error = (sum - TAU * chi).abs();
        cert.gauss_bonnet_bound = 1e-9 * (1.0 + abs);
        if cert.gauss_bonnet_error.is_nan() || cert.gauss_bonnet_error > cert.gauss_bonnet_bound {
            return fail(
                4,
                format!(
                    "the total curvature {sum} differs from 2πχ = {} by {:e}",
                    TAU * chi,
                    cert.gauss_bonnet_error
                ),
            );
        }
        // 5.
        let mut seen = vec![0u8; self.points.len()];
        for (f, list) in self.face_points.iter().enumerate() {
            for (slot, &p) in list.iter().enumerate() {
                let p = p as usize;
                if p >= self.points.len() {
                    return fail(5, format!("face {f} lists point {p}, out of range"));
                }
                seen[p] += 1;
                if self.point_slot[p] as usize != slot {
                    return fail(
                        5,
                        format!(
                            "point {p} is at slot {slot} of face {f}, recorded {}",
                            self.point_slot[p]
                        ),
                    );
                }
                match self.points[p] {
                    Location::Face { face, bary } if face as usize == f => {
                        let s = bary[0] + bary[1] + bary[2];
                        cert.max_bary_sum_error = cert.max_bary_sum_error.max((s - 1.0).abs());
                        if !bary.iter().all(|b| b.is_finite() && *b >= 0.0)
                            || (s - 1.0).abs() > 1e-12
                        {
                            return fail(
                                5,
                                format!("point {p} in face {f} has coordinates {bary:?}"),
                            );
                        }
                    }
                    other => {
                        return fail(
                            5,
                            format!("face {f} lists point {p}, which is at {other:?}"),
                        );
                    }
                }
            }
        }
        for (p, loc) in self.points.iter().enumerate() {
            match *loc {
                Location::Vertex(v) => {
                    if seen[p] != 0 {
                        return fail(5, format!("point {p} is a vertex and in a face list"));
                    }
                    if self.complex_vertex(v).is_none() {
                        return fail(5, format!("point {p} is vertex {v}, which was removed"));
                    }
                }
                Location::Face { .. } => {
                    if seen[p] != 1 {
                        return fail(5, format!("point {p} is in {} face lists", seen[p]));
                    }
                }
            }
        }
        // 7.
        for e in 0..c.n_edges() {
            if c.is_boundary_edge(e) {
                continue;
            }
            let Some(excess) = self.delaunay_excess(e) else {
                return fail(7, format!("edge {e}: its faces are not triangles"));
            };
            cert.max_delaunay_excess = cert.max_delaunay_excess.max(excess);
            if !self.is_delaunay_with(e, excess) {
                cert.n_not_delaunay += 1;
            }
        }
        if delaunay && cert.n_not_delaunay > 0 {
            return fail(
                7,
                format!(
                    "{} interior edges are not Delaunay (largest excess {:e} rad)",
                    cert.n_not_delaunay, cert.max_delaunay_excess
                ),
            );
        }
        Ok(cert)
    }

    /// Whether two coarse meshes hold bitwise the same state: connectivity, lengths, angular
    /// coordinates, references, angle sums, vertex maps, fixed mask, points (coordinates
    /// compared by bits) and per-face point lists in order, options and counters. Scratch
    /// marks are ignored. A testing aid for the exactness of reverts (SPEC S5).
    pub fn bitwise_eq(&self, other: &CoarseMesh) -> bool {
        let bits = |a: &[f64], b: &[f64]| {
            a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
        };
        let loc_eq = |a: &Location, b: &Location| match (a, b) {
            (Location::Vertex(x), Location::Vertex(y)) => x == y,
            (Location::Face { face: f, bary: p }, Location::Face { face: g, bary: q }) => {
                f == g && bits(p, q)
            }
            _ => false,
        };
        self.complex == other.complex
            && bits(&self.length, &other.length)
            && bits(&self.phi, &other.phi)
            && self.reference == other.reference
            && bits(&self.angle_sum, &other.angle_sum)
            && self.fine_of == other.fine_of
            && self.complex_of == other.complex_of
            && self.fixed == other.fixed
            && self.n_input == other.n_input
            && self.points.len() == other.points.len()
            && self
                .points
                .iter()
                .zip(&other.points)
                .all(|(a, b)| loc_eq(a, b))
            && self.point_slot == other.point_slot
            && self.face_points == other.face_points
            && self.options == other.options
            && self.n_clamped == other.n_clamped
            && self.pending.is_none() == other.pending.is_none()
    }

    /// The corners of face `f` as fine ids.
    pub(crate) fn faces_fine(&self, f: usize) -> [u32; 3] {
        self.complex.faces()[f].map(|v| self.fine_of[v as usize])
    }

    /// `θ_k + θ_l − π` of interior edge `e` ([`delaunay_angle_excess`]).
    pub(crate) fn delaunay_excess(&self, e: usize) -> Option<f64> {
        let c = &self.complex;
        let h0 = c.edge_halfedge(e);
        let h1 = c.twin(h0)?;
        let l = |h: usize| self.length[c.edge(h)];
        delaunay_angle_excess(
            l(h0),
            l(c.next(h0)),
            l(c.prev(h0)),
            l(c.next(h1)),
            l(c.prev(h1)),
        )
    }

    /// Whether an interior edge with this angle excess passes the Delaunay criterion.
    pub(crate) fn is_delaunay_with(&self, e: usize, excess: f64) -> bool {
        match self.options.delaunay {
            DelaunayCriterion::AngleSum { tolerance } => excess <= tolerance,
            DelaunayCriterion::CotanWeight { tolerance } => {
                let c = &self.complex;
                let h0 = c.edge_halfedge(e);
                let Some(h1) = c.twin(h0) else { return true };
                let l = |h: usize| self.length[c.edge(h)];
                match ddg_intrinsic::formulas::opposite_cotan_sum(
                    l(h0),
                    l(c.next(h0)),
                    l(c.prev(h0)),
                    l(c.next(h1)),
                    l(c.prev(h1)),
                ) {
                    Some(sum) => 0.5 * sum >= -tolerance,
                    None => false,
                }
            }
        }
    }
}

/// The outgoing halfedges of a vertex counter-clockwise from a start halfedge (see
/// [`CoarseMesh::fan`]).
pub(crate) struct Fan<'a> {
    complex: &'a Complex,
    start: usize,
    current: Option<usize>,
}

impl Iterator for Fan<'_> {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        let h = self.current?;
        self.current = match ccw(self.complex, h) {
            Some(n) if n != self.start => Some(n),
            _ => None,
        };
        Some(h)
    }
}
