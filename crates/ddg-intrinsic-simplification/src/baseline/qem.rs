//! Quadric error metrics simplification: the extrinsic baseline of the paper (Sec. 8.4,
//! Fig. 18, 19, 21, 22; SPEC §10.4, S23).
//!
//! Michael Garland and Paul S. Heckbert, *Surface Simplification Using Quadric Error Metrics*,
//! SIGGRAPH 1997, <https://doi.org/10.1145/258734.258849> (`papers/garland1997_qem.pdf`).
//! Written for this crate from the paper; the reference implementation of the intrinsic paper
//! has no QEM comparison (R43).
//!
//! # The algorithm
//!
//! * **Quadrics.** Every face plane `n·x + d = 0` (unit `n`) gives the quadric
//!   `w (n·x + d)²`, the squared distance to the plane times a weight; a vertex gets the sum
//!   over its faces ([`plane_quadric`], [`Quadric::error`]). The weight is the face area by
//!   default ([`QemOptions::area_weighted`]; the original paper has weight 1): then the error is
//!   an area-weighted integral of the squared distance and does not depend on how the same
//!   surface happens to be triangulated. A boundary edge adds the plane through the edge
//!   perpendicular to its face, with weight [`QemOptions::boundary_weight`] (times the squared
//!   edge length when area weighted), so that moving the boundary off itself is expensive.
//! * **Candidates** are the mesh edges only (no pairs of non-adjacent vertices). The collapse of
//!   `(u, v)` has the quadric `Q_u + Q_v` and goes to the point `p` that minimizes it
//!   ([`optimal_point`], the solution of a 3×3 system) when the system is well conditioned and
//!   `p` is within [`QemOptions::placement_radius`] edge lengths of the edge midpoint; otherwise
//!   to the best of `u`, `v` and the midpoint. Ties take the first of: the optimum, `u`, `v`,
//!   the midpoint.
//! * **Rejected collapses.** The link condition (the common neighbours of `u` and `v` are the
//!   opposite vertices of the edge, and no face at `u` or at `v` joins two common neighbours),
//!   two boundary vertices joined by an interior edge (the collapse would pinch the surface),
//!   and a collapse that turns the normal of a surviving face by more than
//!   [`QemOptions::min_normal_cos`] allows (or makes a face degenerate). These keep the
//!   result a manifold with the topology of the input.
//! * **Order.** A priority queue on the cost `Q(p)`, lazy deletion by per-vertex version stamps;
//!   after a collapse the edges at the new vertex and at its neighbours are queued again with
//!   new costs. Ties are broken by the vertex ids of the edge, then by insertion order: the
//!   result is deterministic. The run stops at the target number of vertices or when no
//!   collapse is allowed ([`QemReport::reached_target`]).
//!
//! # The maps
//!
//! * [`QemResult::vertex_map`]: the coarse vertex that every fine vertex collapsed into.
//! * [`QemResult::projection`]: the closest point of the coarse mesh to every fine vertex (a
//!   face and barycentric coordinates), "simple, documented as ours" (SPEC §10.4); found with a
//!   bounding-volume hierarchy over the coarse faces. It is the map the extrinsic result is
//!   given for prolongation-based comparisons.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use std::fmt;

use ddg_core::TriMesh;

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// A quadric `pᵀ A p + 2 bᵀ p + c` (Garland and Heckbert's `Q` as the matrix
/// `[[A, b], [bᵀ, c]]`, `A` symmetric stored as `a00 a01 a02 a11 a12 a22`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quadric {
    /// The upper triangle of `A`: `a00, a01, a02, a11, a12, a22`.
    pub a: [f64; 6],
    /// The vector `b`.
    pub b: V3,
    /// The constant `c`.
    pub c: f64,
}

impl Quadric {
    /// The zero quadric.
    pub const ZERO: Quadric = Quadric {
        a: [0.0; 6],
        b: [0.0; 3],
        c: 0.0,
    };

    /// The sum of two quadrics (the error of a point is the sum of the errors).
    pub fn plus(&self, o: &Quadric) -> Quadric {
        let mut a = self.a;
        for (x, y) in a.iter_mut().zip(o.a) {
            *x += y;
        }
        Quadric {
            a,
            b: add(self.b, o.b),
            c: self.c + o.c,
        }
    }

    /// The value `pᵀ A p + 2 bᵀ p + c` at `p`.
    ///
    /// [formula: iem-qem-quadric]
    pub fn error(&self, p: V3) -> f64 {
        let [a00, a01, a02, a11, a12, a22] = self.a;
        let ap = [
            a00 * p[0] + a01 * p[1] + a02 * p[2],
            a01 * p[0] + a11 * p[1] + a12 * p[2],
            a02 * p[0] + a12 * p[1] + a22 * p[2],
        ];
        dot(p, ap) + 2.0 * dot(self.b, p) + self.c
    }
}

/// The quadric of the plane through `point` with unit normal `normal`, times `weight`:
/// `weight · (n·x + d)²` with `d = −n·point`, i.e. `A = w n nᵀ`, `b = w d n`, `c = w d²`. Its
/// value at `x` is `weight` times the squared distance of `x` from the plane.
///
/// [formula: iem-qem-quadric]
pub fn plane_quadric(point: V3, normal: V3, weight: f64) -> Quadric {
    let n = normal;
    let d = -dot(n, point);
    Quadric {
        a: [
            weight * n[0] * n[0],
            weight * n[0] * n[1],
            weight * n[0] * n[2],
            weight * n[1] * n[1],
            weight * n[1] * n[2],
            weight * n[2] * n[2],
        ],
        b: scale(n, weight * d),
        c: weight * d * d,
    }
}

/// The point that minimizes the quadric: the solution of `A p = −b` by the adjugate, or `None`
/// when `A` is singular for the purposes of placement: `det A ≤ tolerance · (tr A / 3)³`
/// (for the positive semidefinite `A` of a sum of plane quadrics `(tr A / 3)³ ≥ det A`, the
/// ratio is 1 for three equal orthogonal planes and 0 when the planes do not determine a point)
/// or the result is not finite.
///
/// [formula: iem-qem-quadric]
pub fn optimal_point(q: &Quadric, tolerance: f64) -> Option<V3> {
    let [a, b, c, d, e, f] = q.a; // [[a b c] [b d e] [c e f]]
    let trace = a + d + f;
    if !(trace > 0.0 && trace.is_finite()) {
        return None;
    }
    let c00 = d * f - e * e;
    let c01 = c * e - b * f;
    let c02 = b * e - c * d;
    let det = a * c00 + b * c01 + c * c02;
    let mean = trace / 3.0;
    if !(det > tolerance * mean * mean * mean) {
        return None;
    }
    let c11 = a * f - c * c;
    let c12 = b * c - a * e;
    let c22 = a * d - b * b;
    // inverse = adjugate / det; the adjugate of a symmetric matrix is symmetric.
    let r = [
        c00 * q.b[0] + c01 * q.b[1] + c02 * q.b[2],
        c01 * q.b[0] + c11 * q.b[1] + c12 * q.b[2],
        c02 * q.b[0] + c12 * q.b[1] + c22 * q.b[2],
    ];
    let p = scale(r, -1.0 / det);
    p.iter().all(|x| x.is_finite()).then_some(p)
}

/// Options of [`qem_simplify`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QemOptions {
    /// The number of vertices to stop at (the run stops earlier if no collapse is allowed).
    pub target_vertices: usize,
    /// Weight the plane of a face by its area (default `true`), else by 1 (Garland–Heckbert).
    pub area_weighted: bool,
    /// The weight of the boundary constraint planes (default `1000`); times the squared edge
    /// length when [`area_weighted`](Self::area_weighted).
    pub boundary_weight: f64,
    /// A collapse is rejected when the cosine of the angle between the old and the new normal of
    /// a surviving face is below this (default `0.0`: a turn of 90° or more).
    pub min_normal_cos: f64,
    /// The conditioning test of [`optimal_point`] (default `1e-8`).
    pub singular_tolerance: f64,
    /// The optimal point is used only within this many edge lengths of the midpoint of the edge
    /// (default `4`); a nearly singular system otherwise throws vertices far away.
    pub placement_radius: f64,
}

impl QemOptions {
    /// The defaults, for a target of `target_vertices`.
    pub fn new(target_vertices: usize) -> QemOptions {
        QemOptions {
            target_vertices,
            area_weighted: true,
            boundary_weight: 1000.0,
            min_normal_cos: 0.0,
            singular_tolerance: 1e-8,
            placement_radius: 4.0,
        }
    }
}

/// The closest point of the coarse mesh to a fine vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectedPoint {
    /// The coarse face.
    pub face: u32,
    /// Barycentric coordinates in the corners `faces[face][0..3]`.
    pub bary: [f64; 3],
    /// The distance from the fine vertex.
    pub distance: f64,
}

/// What a run did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QemReport {
    /// Vertices and faces of the input.
    pub n_input_vertices: usize,
    /// Faces of the input.
    pub n_input_faces: usize,
    /// Vertices of the output.
    pub n_output_vertices: usize,
    /// Faces of the output.
    pub n_output_faces: usize,
    /// Whether the target was reached (else no collapse was allowed any more).
    pub reached_target: bool,
    /// Collapses done.
    pub n_collapses: usize,
    /// Queue entries dropped as stale (lazy deletion).
    pub n_stale: usize,
    /// Collapses rejected by the link condition (including the pinch of two boundary vertices).
    pub n_rejected_link: usize,
    /// Collapses rejected because a face normal turned too far or the face degenerated.
    pub n_rejected_flip: usize,
    /// Collapses whose edge was not shared by one or two faces, or whose placement was not
    /// finite.
    pub n_rejected_other: usize,
    /// Collapses that went to the optimal point, to an endpoint, to the midpoint.
    pub n_placed: [usize; 3],
    /// The sum of the costs `Q(p)` of the collapses done.
    pub total_cost: f64,
    /// Total area of the input and of the output.
    pub fine_area: f64,
    /// Total area of the output.
    pub coarse_area: f64,
    /// `100 |coarse − fine| / fine` (Fig. 19's "% error").
    pub area_error_percent: f64,
    /// The largest and the mean distance of a fine vertex from its closest point.
    pub max_projection_distance: f64,
    /// The mean distance.
    pub mean_projection_distance: f64,
}

/// The output of [`qem_simplify`].
#[derive(Debug, Clone, PartialEq)]
pub struct QemResult {
    /// The positions of the coarse vertices.
    pub positions: Vec<V3>,
    /// The coarse faces (orientation of the input kept).
    pub faces: Vec<[u32; 3]>,
    /// For every fine vertex, the coarse vertex it collapsed into (a surviving vertex maps to
    /// itself).
    pub vertex_map: Vec<u32>,
    /// For every fine vertex, its closest point on the coarse mesh.
    pub projection: Vec<ProjectedPoint>,
    /// Counts and measurements.
    pub report: QemReport,
}

/// Why [`qem_simplify`] failed.
#[derive(Debug, Clone, PartialEq)]
pub enum QemError {
    /// The faces do not form an oriented manifold mesh (the message says where).
    InvalidMesh(String),
    /// A position is NaN or infinite.
    NonFinitePosition {
        /// The vertex.
        vertex: usize,
    },
    /// An option is out of range.
    InvalidOptions(String),
}

impl fmt::Display for QemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QemError::InvalidMesh(m) => write!(f, "invalid mesh: {m}"),
            QemError::NonFinitePosition { vertex } => {
                write!(f, "vertex {vertex} has a non-finite position")
            }
            QemError::InvalidOptions(m) => write!(f, "invalid options: {m}"),
        }
    }
}

impl std::error::Error for QemError {}

/// A queue entry; the smallest cost first, then the smallest edge, then the oldest.
#[derive(Debug, Clone, Copy)]
struct Entry {
    cost: f64,
    edge: (u32, u32),
    seq: u64,
    versions: (u32, u32),
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}
impl Eq for Entry {}
impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.edge.cmp(&self.edge))
            .then_with(|| other.seq.cmp(&self.seq))
    }
}

enum Reject {
    Link,
    Flip,
    Other,
}

struct Work {
    pos: Vec<V3>,
    quad: Vec<Quadric>,
    alive: Vec<bool>,
    version: Vec<u32>,
    faces: Vec<[u32; 3]>,
    face_alive: Vec<bool>,
    vfaces: Vec<BTreeSet<u32>>,
    nbrs: Vec<BTreeSet<u32>>,
    parent: Vec<u32>,
    n_alive: usize,
}

impl Work {
    fn edge_faces(&self, a: u32, b: u32) -> Vec<u32> {
        self.vfaces[a as usize]
            .intersection(&self.vfaces[b as usize])
            .copied()
            .collect()
    }

    fn is_boundary_vertex(&self, a: u32) -> bool {
        self.nbrs[a as usize]
            .iter()
            .any(|&b| self.edge_faces(a, b).len() == 1)
    }

    fn recompute_nbrs(&mut self, w: u32) {
        let mut set = BTreeSet::new();
        for &f in &self.vfaces[w as usize] {
            for &x in &self.faces[f as usize] {
                if x != w {
                    set.insert(x);
                }
            }
        }
        self.nbrs[w as usize] = set;
    }

    /// The placement of the collapse of `(u, v)`: the point, its cost and its kind (0 optimum,
    /// 1 endpoint, 2 midpoint).
    fn placement(&self, u: u32, v: u32, o: &QemOptions) -> (V3, f64, usize) {
        let (pu, pv) = (self.pos[u as usize], self.pos[v as usize]);
        let q = self.quad[u as usize].plus(&self.quad[v as usize]);
        let mid = scale(add(pu, pv), 0.5);
        let len = norm(sub(pu, pv));
        let mut best: Option<(V3, f64, usize)> = None;
        let mut consider = |p: V3, kind: usize| {
            let e = q.error(p).max(0.0);
            if e.is_finite() && best.is_none_or(|(_, b, _)| e < b) {
                best = Some((p, e, kind));
            }
        };
        if let Some(p) = optimal_point(&q, o.singular_tolerance)
            && norm(sub(p, mid)) <= o.placement_radius * len
        {
            consider(p, 0);
        }
        consider(pu, 1);
        consider(pv, 1);
        consider(mid, 2);
        best.unwrap_or((mid, f64::INFINITY, 2))
    }

    /// The checks that keep the result a manifold of the input's topology and its normals sane.
    fn validate(&self, u: u32, v: u32, p: V3, o: &QemOptions) -> Result<Vec<u32>, Reject> {
        let ef = self.edge_faces(u, v);
        if ef.is_empty() || ef.len() > 2 || !p.iter().all(|x| x.is_finite()) {
            return Err(Reject::Other);
        }
        let opposite: BTreeSet<u32> = ef
            .iter()
            .filter_map(|&f| {
                self.faces[f as usize]
                    .iter()
                    .copied()
                    .find(|&x| x != u && x != v)
            })
            .collect();
        let common: BTreeSet<u32> = self.nbrs[u as usize]
            .intersection(&self.nbrs[v as usize])
            .copied()
            .collect();
        if common != opposite {
            return Err(Reject::Link);
        }
        // No face at u (or v) between two common neighbours: the edge of the link would be in
        // both links (the tetrahedron, or a triangle that the collapse would double).
        for &w in &[u, v] {
            for &f in &self.vfaces[w as usize] {
                let [a, b, c] = self.faces[f as usize];
                let others: Vec<u32> = [a, b, c].into_iter().filter(|&x| x != w).collect();
                if others.len() == 2 && others.iter().all(|x| common.contains(x)) {
                    return Err(Reject::Link);
                }
            }
        }
        if ef.len() == 2 && self.is_boundary_vertex(u) && self.is_boundary_vertex(v) {
            return Err(Reject::Link);
        }
        // A boundary edge whose triangle has its other two edges on the boundary too (a
        // triangle that is a whole boundary component's last face): the link condition with
        // the vertex at infinity.
        if let (1, Some(&c)) = (ef.len(), opposite.iter().next())
            && self.edge_faces(u, c).len() == 1
            && self.edge_faces(v, c).len() == 1
        {
            return Err(Reject::Link);
        }
        // Normals of the surviving faces.
        for &w in &[u, v] {
            for &f in &self.vfaces[w as usize] {
                if ef.contains(&f) {
                    continue;
                }
                let t = self.faces[f as usize];
                let before = t.map(|x| self.pos[x as usize]);
                let after = t.map(|x| if x == w { p } else { self.pos[x as usize] });
                let n0 = cross(sub(before[1], before[0]), sub(before[2], before[0]));
                let n1 = cross(sub(after[1], after[0]), sub(after[2], after[0]));
                let (l0, l1) = (norm(n0), norm(n1));
                if l0 > 0.0 && (l1 == 0.0 || dot(n0, n1) / (l0 * l1) < o.min_normal_cos) {
                    return Err(Reject::Flip);
                }
            }
        }
        Ok(ef)
    }
}

fn face_normal_area(p: [V3; 3]) -> (V3, f64) {
    let c = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let l = norm(c);
    if l > 0.0 && l.is_finite() {
        (scale(c, 1.0 / l), 0.5 * l)
    } else {
        ([0.0; 3], 0.0)
    }
}

fn total_area(positions: &[V3], faces: &[[u32; 3]]) -> f64 {
    faces
        .iter()
        .map(|f| face_normal_area(f.map(|x| positions[x as usize])).1)
        .sum()
}

/// Simplifies the triangle mesh `(positions, faces)` by quadric edge collapses to
/// `options.target_vertices` vertices. See the [module documentation](self).
///
/// # Errors
/// [`QemError`]: the faces are not an oriented manifold without unreferenced vertices, a
/// position is not finite, an option is out of range.
pub fn qem_simplify(
    positions: &[V3],
    faces: &[[u32; 3]],
    options: &QemOptions,
) -> Result<QemResult, QemError> {
    if !(options.boundary_weight.is_finite() && options.boundary_weight >= 0.0) {
        return Err(QemError::InvalidOptions(format!(
            "boundary_weight must be a finite number >= 0, not {}",
            options.boundary_weight
        )));
    }
    if !(options.min_normal_cos.is_finite()
        && options.singular_tolerance.is_finite()
        && options.singular_tolerance >= 0.0
        && options.placement_radius.is_finite()
        && options.placement_radius >= 0.0)
    {
        return Err(QemError::InvalidOptions(
            "min_normal_cos, singular_tolerance (>= 0) and placement_radius (>= 0) must be finite"
                .to_string(),
        ));
    }
    if let Some(vertex) = positions
        .iter()
        .position(|p| p.iter().any(|x| !x.is_finite()))
    {
        return Err(QemError::NonFinitePosition { vertex });
    }
    let n = positions.len();
    TriMesh::new(n, faces.to_vec()).map_err(|e| QemError::InvalidMesh(e.to_string()))?;
    // Quadrics.
    let mut quad = vec![Quadric::ZERO; n];
    let mut edge_faces: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
    for (fi, f) in faces.iter().enumerate() {
        let p = f.map(|x| positions[x as usize]);
        let (nrm, area) = face_normal_area(p);
        if area > 0.0 {
            let w = if options.area_weighted { area } else { 1.0 };
            let q = plane_quadric(p[0], nrm, w);
            for &x in f {
                quad[x as usize] = quad[x as usize].plus(&q);
            }
        }
        for k in 0..3 {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            edge_faces
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push(fi as u32);
        }
    }
    for (&(a, b), fs) in &edge_faces {
        if fs.len() != 1 {
            continue;
        }
        let (pa, pb) = (positions[a as usize], positions[b as usize]);
        let f = faces[fs[0] as usize];
        let (nrm, _) = face_normal_area(f.map(|x| positions[x as usize]));
        let e = sub(pb, pa);
        let len = norm(e);
        let nb = cross(e, nrm);
        let lb = norm(nb);
        if len > 0.0 && lb > 0.0 {
            let w = options.boundary_weight
                * if options.area_weighted {
                    len * len
                } else {
                    1.0
                };
            let q = plane_quadric(pa, scale(nb, 1.0 / lb), w);
            quad[a as usize] = quad[a as usize].plus(&q);
            quad[b as usize] = quad[b as usize].plus(&q);
        }
    }

    let mut vfaces = vec![BTreeSet::new(); n];
    for (fi, f) in faces.iter().enumerate() {
        for &x in f {
            vfaces[x as usize].insert(fi as u32);
        }
    }
    let mut w = Work {
        pos: positions.to_vec(),
        quad,
        alive: vec![true; n],
        version: vec![0; n],
        faces: faces.to_vec(),
        face_alive: vec![true; faces.len()],
        vfaces,
        nbrs: vec![BTreeSet::new(); n],
        parent: (0..n as u32).collect(),
        n_alive: n,
    };
    for v in 0..n as u32 {
        w.recompute_nbrs(v);
    }

    let mut report = QemReport {
        n_input_vertices: n,
        n_input_faces: faces.len(),
        n_output_vertices: n,
        n_output_faces: faces.len(),
        reached_target: n <= options.target_vertices,
        n_collapses: 0,
        n_stale: 0,
        n_rejected_link: 0,
        n_rejected_flip: 0,
        n_rejected_other: 0,
        n_placed: [0; 3],
        total_cost: 0.0,
        fine_area: total_area(positions, faces),
        coarse_area: 0.0,
        area_error_percent: 0.0,
        max_projection_distance: 0.0,
        mean_projection_distance: 0.0,
    };

    let mut heap = BinaryHeap::new();
    let mut seq = 0u64;
    let push = |w: &Work, heap: &mut BinaryHeap<Entry>, seq: &mut u64, a: u32, b: u32| {
        let (u, v) = (a.min(b), a.max(b));
        let (_, cost, _) = w.placement(u, v, options);
        *seq += 1;
        heap.push(Entry {
            cost,
            edge: (u, v),
            seq: *seq,
            versions: (w.version[u as usize], w.version[v as usize]),
        });
    };
    for v in 0..n as u32 {
        for &x in &w.nbrs[v as usize] {
            if v < x {
                push(&w, &mut heap, &mut seq, v, x);
            }
        }
    }

    while w.n_alive > options.target_vertices {
        let Some(entry) = heap.pop() else { break };
        let (u, v) = entry.edge;
        let (iu, iv) = (u as usize, v as usize);
        if !(w.alive[iu]
            && w.alive[iv]
            && w.version[iu] == entry.versions.0
            && w.version[iv] == entry.versions.1
            && w.nbrs[iu].contains(&v))
        {
            report.n_stale += 1;
            continue;
        }
        let (p, cost, kind) = w.placement(u, v, options);
        let ef = match w.validate(u, v, p, options) {
            Ok(ef) => ef,
            Err(Reject::Link) => {
                report.n_rejected_link += 1;
                continue;
            }
            Err(Reject::Flip) => {
                report.n_rejected_flip += 1;
                continue;
            }
            Err(Reject::Other) => {
                report.n_rejected_other += 1;
                continue;
            }
        };
        // Collapse v into u.
        let mut affected: BTreeSet<u32> = w.nbrs[iu].union(&w.nbrs[iv]).copied().collect();
        affected.remove(&v);
        affected.insert(u);
        for &f in &ef {
            w.face_alive[f as usize] = false;
            for k in 0..3 {
                let x = w.faces[f as usize][k] as usize;
                w.vfaces[x].remove(&f);
            }
        }
        let moved: Vec<u32> = w.vfaces[iv].iter().copied().collect();
        for f in moved {
            for x in w.faces[f as usize].iter_mut() {
                if *x == v {
                    *x = u;
                }
            }
            w.vfaces[iu].insert(f);
        }
        w.vfaces[iv].clear();
        w.nbrs[iv].clear();
        w.alive[iv] = false;
        w.parent[iv] = u;
        w.n_alive -= 1;
        w.quad[iu] = w.quad[iu].plus(&w.quad[iv]);
        w.pos[iu] = p;
        w.version[iu] += 1;
        for &a in &affected {
            w.recompute_nbrs(a);
        }
        report.n_collapses += 1;
        report.n_placed[kind] += 1;
        report.total_cost += cost;
        // Re-queue the edges at the new vertex and at its neighbours.
        let ring: Vec<u32> = w.nbrs[iu].iter().copied().collect();
        for &r in &ring {
            w.version[r as usize] += 1;
        }
        let mut edges = BTreeSet::new();
        for &r in &ring {
            edges.insert((u.min(r), u.max(r)));
            for &x in &w.nbrs[r as usize] {
                edges.insert((r.min(x), r.max(x)));
            }
        }
        for (a, b) in edges {
            push(&w, &mut heap, &mut seq, a, b);
        }
    }
    report.reached_target = w.n_alive <= options.target_vertices;

    // Compact.
    let mut new_id = vec![u32::MAX; n];
    let mut out_pos = Vec::with_capacity(w.n_alive);
    for v in 0..n {
        if w.alive[v] {
            new_id[v] = out_pos.len() as u32;
            out_pos.push(w.pos[v]);
        }
    }
    let out_faces: Vec<[u32; 3]> = w
        .faces
        .iter()
        .zip(&w.face_alive)
        .filter(|&(_, &a)| a)
        .map(|(f, _)| f.map(|x| new_id[x as usize]))
        .collect();
    let vertex_map: Vec<u32> = (0..n)
        .map(|v| {
            let mut r = v as u32;
            while w.parent[r as usize] != r {
                r = w.parent[r as usize];
            }
            new_id[r as usize]
        })
        .collect();

    let projection = project_onto(&out_pos, &out_faces, positions);
    report.n_output_vertices = out_pos.len();
    report.n_output_faces = out_faces.len();
    report.coarse_area = total_area(&out_pos, &out_faces);
    report.area_error_percent = if report.fine_area > 0.0 {
        100.0 * (report.coarse_area - report.fine_area).abs() / report.fine_area
    } else {
        0.0
    };
    report.max_projection_distance = projection.iter().fold(0.0, |m, p| m.max(p.distance));
    report.mean_projection_distance =
        projection.iter().map(|p| p.distance).sum::<f64>() / projection.len().max(1) as f64;
    Ok(QemResult {
        positions: out_pos,
        faces: out_faces,
        vertex_map,
        projection,
        report,
    })
}

// ----- closest point projection -----

/// The closest point of segment `ab` to `p`: the parameter along it.
fn segment_param(p: V3, a: V3, b: V3) -> f64 {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    if l2 > 0.0 {
        (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

/// The barycentric coordinates of the closest point of triangle `abc` to `p`, and the squared
/// distance (Ericson, *Real-Time Collision Detection*, 5.1.5; degenerate triangles fall back to
/// their three edges).
fn closest_on_triangle(p: V3, t: [V3; 3]) -> ([f64; 3], f64) {
    let [a, b, c] = t;
    let from_bary = |l: [f64; 3]| {
        let q = add(add(scale(a, l[0]), scale(b, l[1])), scale(c, l[2]));
        let d = sub(p, q);
        (l, dot(d, d))
    };
    let edges = || {
        let mut best: Option<([f64; 3], f64)> = None;
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            let s = segment_param(p, t[i], t[j]);
            let mut l = [0.0; 3];
            l[i] = 1.0 - s;
            l[j] = s;
            let r = from_bary(l);
            if best.is_none_or(|(_, d)| r.1 < d) {
                best = Some(r);
            }
        }
        best.expect("three edges")
    };
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return from_bary([1.0, 0.0, 0.0]);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return from_bary([0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let den = d1 - d3;
        if den > 0.0 {
            let v = d1 / den;
            return from_bary([1.0 - v, v, 0.0]);
        }
        return edges();
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return from_bary([0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let den = d2 - d6;
        if den > 0.0 {
            let w = d2 / den;
            return from_bary([1.0 - w, 0.0, w]);
        }
        return edges();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let den = (d4 - d3) + (d5 - d6);
        if den > 0.0 {
            let w = (d4 - d3) / den;
            return from_bary([0.0, 1.0 - w, w]);
        }
        return edges();
    }
    let denom = va + vb + vc;
    if !(denom > 0.0 && denom.is_finite()) {
        return edges();
    }
    let v = vb / denom;
    let w = vc / denom;
    from_bary([1.0 - v - w, v, w])
}

struct Node {
    lo: V3,
    hi: V3,
    /// Children, or `u32::MAX` for a leaf holding `order[start..start + count]`.
    left: u32,
    right: u32,
    start: u32,
    count: u32,
}

fn box_dist2(p: V3, lo: V3, hi: V3) -> f64 {
    let mut d2 = 0.0;
    for k in 0..3 {
        let d = (lo[k] - p[k]).max(0.0).max(p[k] - hi[k]);
        d2 += d * d;
    }
    d2
}

/// A bounding-volume hierarchy over triangles (median split of the centroids).
struct Bvh {
    nodes: Vec<Node>,
    order: Vec<u32>,
}

impl Bvh {
    fn build(positions: &[V3], faces: &[[u32; 3]]) -> Bvh {
        let boxes: Vec<(V3, V3, V3)> = faces
            .iter()
            .map(|f| {
                let p = f.map(|x| positions[x as usize]);
                let mut lo = p[0];
                let mut hi = p[0];
                for q in &p[1..] {
                    for k in 0..3 {
                        lo[k] = lo[k].min(q[k]);
                        hi[k] = hi[k].max(q[k]);
                    }
                }
                let c = scale(add(add(p[0], p[1]), p[2]), 1.0 / 3.0);
                (lo, hi, c)
            })
            .collect();
        let mut bvh = Bvh {
            nodes: Vec::new(),
            order: (0..faces.len() as u32).collect(),
        };
        if !faces.is_empty() {
            bvh.build_range(&boxes, 0, faces.len());
        }
        bvh
    }

    fn build_range(&mut self, boxes: &[(V3, V3, V3)], start: usize, end: usize) -> u32 {
        let mut lo = boxes[self.order[start] as usize].0;
        let mut hi = boxes[self.order[start] as usize].1;
        let mut clo = boxes[self.order[start] as usize].2;
        let mut chi = clo;
        for &i in &self.order[start..end] {
            let b = &boxes[i as usize];
            for k in 0..3 {
                lo[k] = lo[k].min(b.0[k]);
                hi[k] = hi[k].max(b.1[k]);
                clo[k] = clo[k].min(b.2[k]);
                chi[k] = chi[k].max(b.2[k]);
            }
        }
        let id = self.nodes.len() as u32;
        self.nodes.push(Node {
            lo,
            hi,
            left: u32::MAX,
            right: u32::MAX,
            start: start as u32,
            count: (end - start) as u32,
        });
        if end - start > 4 {
            let axis = (0..3)
                .max_by(|&a, &b| (chi[a] - clo[a]).total_cmp(&(chi[b] - clo[b])))
                .expect("three axes");
            self.order[start..end].sort_by(|&a, &b| {
                boxes[a as usize].2[axis]
                    .total_cmp(&boxes[b as usize].2[axis])
                    .then(a.cmp(&b))
            });
            let mid = (start + end) / 2;
            let l = self.build_range(boxes, start, mid);
            let r = self.build_range(boxes, mid, end);
            self.nodes[id as usize].left = l;
            self.nodes[id as usize].right = r;
        }
        id
    }

    /// The closest triangle to `p`, ties to the lowest face id: `(face, bary, d²)`.
    fn closest(&self, positions: &[V3], faces: &[[u32; 3]], p: V3) -> (u32, [f64; 3], f64) {
        let mut best = (u32::MAX, [0.0; 3], f64::INFINITY);
        let mut stack = if self.nodes.is_empty() {
            Vec::new()
        } else {
            vec![0u32]
        };
        while let Some(id) = stack.pop() {
            let node = &self.nodes[id as usize];
            if box_dist2(p, node.lo, node.hi) > best.2 {
                continue;
            }
            if node.left == u32::MAX {
                for &fi in &self.order[node.start as usize..(node.start + node.count) as usize] {
                    let t = faces[fi as usize].map(|x| positions[x as usize]);
                    let (bary, d2) = closest_on_triangle(p, t);
                    if d2 < best.2 || (d2 == best.2 && fi < best.0) {
                        best = (fi, bary, d2);
                    }
                }
            } else {
                stack.push(node.right);
                stack.push(node.left);
            }
        }
        best
    }
}

/// The closest point of the mesh `(positions, faces)` to every point of `points`.
fn project_onto(positions: &[V3], faces: &[[u32; 3]], points: &[V3]) -> Vec<ProjectedPoint> {
    let bvh = Bvh::build(positions, faces);
    points
        .iter()
        .map(|&p| {
            let (face, bary, d2) = bvh.closest(positions, faces, p);
            ProjectedPoint {
                face,
                bary,
                distance: d2.sqrt(),
            }
        })
        .collect()
}
