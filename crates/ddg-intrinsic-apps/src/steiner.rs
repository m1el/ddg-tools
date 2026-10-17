//! A short cut through all vertices: an approximate Steiner tree over the intrinsic Delaunay
//! triangulation (paper Sec. 5.1, Fig. 15).
//!
//! "Many applications require a short cut passing through all vertices of a polyhedron, so
//! that it can be isometrically unfolded into the plane" (Sec. 5.1). Cutting along a tree
//! that contains every vertex of a closed surface of genus zero leaves a disk whose interior
//! has no vertex, hence no curvature: it develops into the plane without distortion (with
//! possible overlaps). The common approach cuts along input edges, whose shortest such tree is
//! the minimal spanning tree (MST) of the edge graph (Sheffer and Hart 2002). Here, as in the
//! paper, the tree is made of intrinsic edges: the MST of the intrinsic Delaunay triangulation
//! (in the plane the Delaunay triangulation contains the Euclidean MST, Toussaint 1980), made
//! shorter by Steiner points, the Fermat point of a triangle inserted as an intrinsic vertex
//! with its three edges to the corners.
//!
//! # The algorithm (SPEC N30 to N32)
//!
//! The heuristic of Smith, Lee and Liebman (1981) as the paper summarizes it ("enqueue all
//! edges and per-triangle Steiner trees; greedily add the shortest element from the queue as
//! long as it does not form a cycle"), with the priority of the original, which the summary
//! leaves out (Zachariasen and Winter 1997, Sec. 1 and 6.2, describe it): a triangle's tree is
//! ranked by the ratio of its length to that of the MST of its three corners, the edges come
//! after all triangles in order of length. (By raw length a triangle's tree, which is longer
//! than each of its sides, would never be taken: its corners are joined by its two shorter
//! sides first.)
//!
//! 1. Flip to the intrinsic Delaunay triangulation (option).
//! 2. The MST of the edges (Kruskal; self-edges skipped). On a surface with boundary the
//!    vertices of every boundary loop count as joined from the start: the boundary is cut
//!    already (N31).
//! 3. Candidate triangles: three different corners, every angle below `120°` (otherwise the
//!    Fermat point is the obtuse corner and the tree is two sides), and with
//!    [`TriangleCandidates::MstAdjacent`] (default, as Smith, Lee and Liebman) two of their
//!    sides in the MST. The Fermat point has the barycentric coordinates `a_k csc(θ_k + π/3)`
//!    normalized (`a_k` the side opposite corner `k`); the paper prints `csc(θ_k + π/3)`,
//!    which are its trilinear coordinates (N30).
//! 4. Kruskal over the candidates by ratio, then over all edges by length.
//! 5. The Fermat point of every chosen triangle is inserted (`insert_vertex`); its three new
//!    edges join the tree. An insertion that is refused (a point too close to a side) drops
//!    the triangle and steps 4 and 5 are repeated without it (N32).
//!
//! With `MstAdjacent` each chosen triangle replaces its two MST sides by its shorter tree and
//! the remaining MST edges complete the tree, so the result is never longer than the MST.
//!
//! [`cut_to_disk`] cuts a triangulation along a set of edges and reports what the cut leaves:
//! the Euler characteristic, whether it is one disk, the length of its boundary, the angle
//! defect of every vertex not on the cut (zero if the disk is flat inside), and a layout of
//! every triangle in the plane found by unfolding across the edges that are not cut, with the
//! largest disagreement between the places of one vertex of the cut surface.
//!
//! # References
//!
//! A. Sheffer, J. C. Hart. *Seamster: inconspicuous low-distortion texture seam layout.* IEEE
//! Visualization 2002, 291–298. <https://doi.org/10.1109/VISUAL.2002.1183787>
//!
//! J. M. Smith, D. T. Lee, J. S. Liebman. *An O(n log n) heuristic for Steiner minimal tree
//! problems on the Euclidean metric.* Networks 11(1), 23–39, 1981.
//! <https://doi.org/10.1002/net.3230110104>
//!
//! M. Zachariasen, P. Winter. *Concatenation-based greedy heuristics for the Euclidean Steiner
//! tree problem.* Technical report DIKU-TR-97/20, University of Copenhagen, 1997; Algorithmica
//! 25, 418–437, 1999. <https://doi.org/10.1007/PL00009287>
//!
//! G. T. Toussaint. *The relative neighbourhood graph of a finite planar set.* Pattern
//! Recognition 12(4), 261–268, 1980. <https://doi.org/10.1016/0031-3203(80)90066-7>

use std::f64::consts::{FRAC_PI_3, PI, TAU};
use std::fmt;

use ddg_core::geometry::corner_angle;
use ddg_intrinsic::formulas::{corner_distances, triangle_layout, unfold_apex};
use ddg_intrinsic::{
    DelaunayOptions, DelaunayReport, EdgeTrace, EdgeTraceOptions, IntrinsicError,
    IntrinsicTriangulation,
};

/// Which triangles may contribute a Steiner point. See the [module documentation](self).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TriangleCandidates {
    /// Triangles with two sides in the minimal spanning tree (Smith, Lee and Liebman).
    #[default]
    MstAdjacent,
    /// Every triangle with three different corners and all angles below `120°`.
    All,
}

/// Options of [`steiner_tree`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteinerOptions {
    /// Flip to the intrinsic Delaunay triangulation first. `None` keeps the triangulation as it
    /// is: on a fresh one that is the tree over the input edges and input faces (the extrinsic
    /// version). Default: flip with the default options.
    pub delaunay: Option<DelaunayOptions>,
    /// Insert Fermat points. `false` returns the minimal spanning tree. Default `true`.
    pub steiner_points: bool,
    /// Which triangles are candidates. Default [`TriangleCandidates::MstAdjacent`].
    pub candidates: TriangleCandidates,
    /// Trace the tree's edges over the input mesh. Default: the default options.
    pub trace: Option<EdgeTraceOptions>,
}

impl Default for SteinerOptions {
    fn default() -> Self {
        SteinerOptions {
            delaunay: Some(DelaunayOptions::default()),
            steiner_points: true,
            candidates: TriangleCandidates::MstAdjacent,
            trace: Some(EdgeTraceOptions::default()),
        }
    }
}

/// What [`steiner_tree`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct SteinerReport {
    /// The report of the Delaunay flips, if they were asked for.
    pub delaunay: Option<DelaunayReport>,
    /// Length of the minimal spanning tree of the triangulation's edges (after the flips).
    pub mst_length: f64,
    /// Length of the result.
    pub length: f64,
    /// Number of candidate triangles.
    pub n_candidates: usize,
    /// Number of Fermat points inserted.
    pub n_steiner_points: usize,
    /// Number of triangles dropped because the insertion of their Fermat point was refused.
    pub n_refused: usize,
    /// Number of rounds of steps 4 and 5 (one unless an insertion was refused).
    pub n_rounds: usize,
    /// Number of boundary loops, whose vertices count as joined.
    pub n_boundary_loops: usize,
    /// The largest mismatch of an edge trace relative to its edge (zero without traces).
    pub max_trace_mismatch: f64,
}

/// A tree through all vertices. See [`steiner_tree`].
#[derive(Debug, Clone)]
pub struct SteinerTree {
    /// The triangulation the tree lives on: the input triangulation, flipped if asked, with the
    /// Fermat points inserted.
    pub triangulation: IntrinsicTriangulation,
    /// The edges of the tree (ids of `triangulation`), sorted.
    pub edges: Vec<u32>,
    /// The inserted Fermat points (vertex ids of `triangulation`).
    pub steiner_vertices: Vec<u32>,
    /// The total length of the edges.
    pub length: f64,
    /// The path of every tree edge over the input mesh, in the order of `edges` (empty when
    /// [`SteinerOptions::trace`] is `None`).
    pub traces: Vec<EdgeTrace>,
    /// What was done.
    pub report: SteinerReport,
}

/// Why [`steiner_tree`] or [`cut_to_disk`] failed.
#[derive(Debug, Clone, PartialEq)]
pub enum SteinerError {
    /// An operation of the triangulation failed (flips, insertion other than a refusal,
    /// tracing).
    Intrinsic(IntrinsicError),
    /// An edge id given to [`cut_to_disk`] does not exist.
    EdgeOutOfRange {
        /// The edge.
        edge: usize,
        /// The number of edges.
        n_edges: usize,
    },
}

impl fmt::Display for SteinerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SteinerError::Intrinsic(e) => write!(f, "{e}"),
            SteinerError::EdgeOutOfRange { edge, n_edges } => {
                write!(f, "edge {edge} does not exist (the triangulation has {n_edges} edges)")
            }
        }
    }
}

impl std::error::Error for SteinerError {}

impl From<IntrinsicError> for SteinerError {
    fn from(e: IntrinsicError) -> Self {
        SteinerError::Intrinsic(e)
    }
}

struct UnionFind(Vec<u32>);

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind((0..n as u32).collect())
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] as usize != x {
            self.0[x] = self.0[self.0[x] as usize];
            x = self.0[x] as usize;
        }
        x
    }
    /// Joins the sets of `a` and `b`; false if they were joined already.
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        self.0[ra.max(rb)] = ra.min(rb) as u32;
        true
    }
}

/// The vertices joined from the start: those of every boundary loop. Returns the union-find
/// and the number of loops.
fn joined_boundary(tri: &IntrinsicTriangulation) -> (UnionFind, usize) {
    let mut uf = UnionFind::new(tri.n_vertices());
    let mut loops = UnionFind::new(tri.n_vertices());
    let mut on_boundary = vec![false; tri.n_vertices()];
    for e in 0..tri.n_edges() {
        if tri.is_boundary_edge(e) {
            let [a, b] = tri.edge_vertices(e);
            uf.union(a, b);
            loops.union(a, b);
            on_boundary[a] = true;
            on_boundary[b] = true;
        }
    }
    let n_loops = (0..tri.n_vertices())
        .filter(|&v| on_boundary[v] && loops.find(v) == v)
        .count();
    (uf, n_loops)
}

/// The edges sorted by length (ties by id), self-edges and boundary edges left out.
fn edges_by_length(tri: &IntrinsicTriangulation) -> Vec<usize> {
    let mut edges: Vec<usize> = (0..tri.n_edges())
        .filter(|&e| {
            let [a, b] = tri.edge_vertices(e);
            a != b && !tri.is_boundary_edge(e)
        })
        .collect();
    edges.sort_by(|&x, &y| {
        tri.edge_length(x)
            .total_cmp(&tri.edge_length(y))
            .then(x.cmp(&y))
    });
    edges
}

/// The minimal spanning tree (forest, one tree per component; boundary loops joined from the
/// start) of the triangulation's edges by Kruskal's algorithm: its edges, sorted, and its
/// length.
pub fn minimum_spanning_tree(tri: &IntrinsicTriangulation) -> (Vec<u32>, f64) {
    let (mut uf, _) = joined_boundary(tri);
    let mut tree = Vec::new();
    let mut length = 0.0;
    for e in edges_by_length(tri) {
        let [a, b] = tri.edge_vertices(e);
        if uf.union(a, b) {
            tree.push(e as u32);
            length += tri.edge_length(e);
        }
    }
    tree.sort_unstable();
    (tree, length)
}

/// The Fermat point of a triangle with side lengths `l` (`l[k]` from corner `k` to `k + 1`)
/// whose angles are all below `120°`: its barycentric coordinates
/// `b_k ∝ a_k / sin(θ_k + π/3)` (`a_k = l[k+1]` the side opposite corner `k`, `θ_k` its angle),
/// its distances to the corners, and their sum, the length of the Steiner tree of the three
/// corners. `None` if an angle is `120°` or more or the lengths are not a triangle.
///
/// [formula: app-steiner-fermat-point]
pub fn fermat_point(l: [f64; 3]) -> Option<([f64; 3], [f64; 3], f64)> {
    let theta = [0, 1, 2].map(|k| corner_angle(l[(k + 1) % 3], l[k], l[(k + 2) % 3]));
    let mut b = [0.0; 3];
    for k in 0..3 {
        let t = theta[k]?;
        if t >= 2.0 * FRAC_PI_3 {
            return None;
        }
        b[k] = l[(k + 1) % 3] / (t + FRAC_PI_3).sin();
    }
    let sum = b[0] + b[1] + b[2];
    let bary = b.map(|x| x / sum);
    let d = corner_distances(l, bary)?;
    Some((bary, d, d[0] + d[1] + d[2]))
}

struct Candidate {
    face: usize,
    ratio: f64,
    bary: [f64; 3],
}

/// A short tree through all vertices of the triangulation (paper Sec. 5.1); see the
/// [module documentation](self). The input triangulation is not changed.
///
/// # Errors
/// [`SteinerError::Intrinsic`] when the flips, an insertion (other than a refusal, which drops
/// the triangle) or an edge trace fail.
pub fn steiner_tree(
    tri: &IntrinsicTriangulation,
    options: &SteinerOptions,
) -> Result<SteinerTree, SteinerError> {
    let mut base = tri.clone();
    let delaunay = match &options.delaunay {
        Some(d) => Some(base.flip_to_delaunay(d)?),
        None => None,
    };
    let (mst, mst_length) = minimum_spanning_tree(&base);
    let (_, n_boundary_loops) = joined_boundary(&base);
    let mut in_mst = vec![false; base.n_edges()];
    for &e in &mst {
        in_mst[e as usize] = true;
    }
    let mut candidates: Vec<Candidate> = Vec::new();
    if options.steiner_points {
        for (f, face) in base.faces().iter().enumerate() {
            if face[0] == face[1] || face[1] == face[2] || face[2] == face[0] {
                continue;
            }
            let sides = [0, 1, 2].map(|k| base.edge(3 * f + k));
            let n_mst = sides.iter().filter(|&&e| in_mst[e]).count();
            if options.candidates == TriangleCandidates::MstAdjacent && n_mst < 2 {
                continue;
            }
            let l = base.face_lengths(f);
            let Some((bary, _, s)) = fermat_point(l) else {
                continue;
            };
            let longest = l[0].max(l[1]).max(l[2]);
            let mst3 = l[0] + l[1] + l[2] - longest;
            candidates.push(Candidate {
                face: f,
                ratio: s / mst3,
                bary,
            });
        }
        candidates.sort_by(|x, y| x.ratio.total_cmp(&y.ratio).then(x.face.cmp(&y.face)));
    }
    let n_candidates = candidates.len();
    let edges = edges_by_length(&base);
    let mut excluded = vec![false; base.n_faces()];
    let mut n_refused = 0;
    let mut n_rounds = 0;
    loop {
        n_rounds += 1;
        let (mut uf, _) = joined_boundary(&base);
        let mut chosen: Vec<&Candidate> = Vec::new();
        for c in &candidates {
            if excluded[c.face] {
                continue;
            }
            let [a, b, d] = base.faces()[c.face].map(|v| v as usize);
            let (ra, rb, rd) = (uf.find(a), uf.find(b), uf.find(d));
            if ra != rb && rb != rd && ra != rd {
                uf.union(a, b);
                uf.union(b, d);
                chosen.push(c);
            }
        }
        let mut tree: Vec<u32> = Vec::new();
        for &e in &edges {
            let [a, b] = base.edge_vertices(e);
            if uf.union(a, b) {
                tree.push(e as u32);
            }
        }
        let mut result = base.clone();
        let mut steiner_vertices = Vec::new();
        let mut refused = false;
        for c in &chosen {
            match result.insert_vertex(c.face, c.bary) {
                Ok(insertion) => {
                    steiner_vertices.push(insertion.vertex);
                    tree.extend_from_slice(&insertion.new_edges);
                }
                Err(IntrinsicError::InsertionRefused { .. }) => {
                    excluded[c.face] = true;
                    n_refused += 1;
                    refused = true;
                }
                Err(e) => return Err(e.into()),
            }
        }
        if refused {
            continue;
        }
        tree.sort_unstable();
        let length: f64 = tree.iter().map(|&e| result.edge_length(e as usize)).sum();
        let mut traces = Vec::new();
        let mut max_trace_mismatch = 0.0f64;
        if let Some(trace_options) = &options.trace {
            for &e in &tree {
                let t = result.trace_edge(e as usize, trace_options)?;
                max_trace_mismatch =
                    max_trace_mismatch.max(t.mismatch / result.edge_length(e as usize));
                traces.push(t);
            }
        }
        let report = SteinerReport {
            delaunay,
            mst_length,
            length,
            n_candidates,
            n_steiner_points: steiner_vertices.len(),
            n_refused,
            n_rounds,
            n_boundary_loops,
            max_trace_mismatch,
        };
        return Ok(SteinerTree {
            triangulation: result,
            edges: tree,
            steiner_vertices,
            length,
            traces,
            report,
        });
    }
}

/// What cutting a triangulation along a set of edges leaves. See [`cut_to_disk`].
#[derive(Debug, Clone, PartialEq)]
pub struct CutReport {
    /// Vertices of the cut surface (a vertex on the cut splits into one per wedge between cut
    /// edges).
    pub n_vertices: usize,
    /// Edges of the cut surface (a cut interior edge counts twice).
    pub n_edges: usize,
    /// Faces.
    pub n_faces: usize,
    /// `n_vertices − n_edges + n_faces`.
    pub euler_characteristic: i64,
    /// Whether the faces are connected across edges that are not cut.
    pub connected: bool,
    /// Connected with Euler characteristic one: a disk (the surface is orientable).
    pub is_disk: bool,
    /// The length of the boundary of the cut surface: twice the cut interior edges, once the
    /// boundary edges.
    pub boundary_length: f64,
    /// The largest `|2π − angle sum|` of a vertex that is not on the cut or the boundary: zero
    /// when every vertex with curvature lies on the cut, so that the disk is flat inside and
    /// unfolds.
    pub max_interior_defect: f64,
    /// The sum of `π − wedge angle` over the vertices of the cut surface on its boundary: `2π`
    /// for a disk that is flat inside (Gauss–Bonnet).
    pub boundary_turning: f64,
    /// The corners of every face in the plane, from unfolding across the edges that are not
    /// cut (one breadth-first unfolding per connected piece).
    pub layout: Vec<[[f64; 2]; 3]>,
    /// The largest distance between two places of one vertex of the cut surface in the layout,
    /// relative to the mean edge length: zero up to rounding when the cut surface is a flat
    /// disk.
    pub max_layout_mismatch: f64,
}

/// Cuts the triangulation along `edges` (and its boundary) and reports what is left: see
/// [`CutReport`]. Paper Sec. 5.1: a tree through all vertices of a closed surface of genus
/// zero leaves a disk that unfolds isometrically into the plane.
///
/// # Errors
/// [`SteinerError::EdgeOutOfRange`].
pub fn cut_to_disk(tri: &IntrinsicTriangulation, edges: &[u32]) -> Result<CutReport, SteinerError> {
    let n_edges = tri.n_edges();
    let mut cut = vec![false; n_edges];
    for &e in edges {
        let e = e as usize;
        if e >= n_edges {
            return Err(SteinerError::EdgeOutOfRange { edge: e, n_edges });
        }
        cut[e] = true;
    }
    for (e, c) in cut.iter_mut().enumerate() {
        if tri.is_boundary_edge(e) {
            *c = true;
        }
    }
    // Vertices of the cut surface: one per wedge of corners between cut edges.
    let n_h = tri.n_halfedges();
    let mut corner_vertex = vec![u32::MAX; n_h];
    let mut wedge_angle: Vec<f64> = Vec::new();
    let mut wedge_on_boundary: Vec<bool> = Vec::new();
    for v in 0..tri.n_vertices() {
        let around: Vec<usize> = tri.outgoing_halfedges(v).collect();
        let d = around.len();
        // Between corner i and corner i + 1 lies edge(prev(h_i)) = edge(h_{i+1}).
        let boundary = tri.is_boundary_vertex(v);
        let start = if boundary {
            0
        } else {
            (0..d).find(|&i| cut[tri.edge(around[i])]).unwrap_or(0)
        };
        let any_cut = boundary || around.iter().any(|&h| cut[tri.edge(h)]);
        let mut current = wedge_angle.len();
        wedge_angle.push(0.0);
        wedge_on_boundary.push(any_cut);
        for step in 0..d {
            let i = (start + step) % d;
            let h = around[i];
            if step > 0 && cut[tri.edge(h)] {
                current = wedge_angle.len();
                wedge_angle.push(0.0);
                wedge_on_boundary.push(true);
            }
            corner_vertex[h] = current as u32;
            wedge_angle[current] += tri.corner_angle(h);
        }
    }
    let n_vertices = wedge_angle.len();
    let n_cut_interior = (0..n_edges)
        .filter(|&e| cut[e] && !tri.is_boundary_edge(e))
        .count();
    let n_cut_edges = n_edges + n_cut_interior;
    let n_faces = tri.n_faces();
    let euler_characteristic = n_vertices as i64 - n_cut_edges as i64 + n_faces as i64;
    let boundary_length: f64 = (0..n_edges)
        .filter(|&e| cut[e])
        .map(|e| {
            let factor = if tri.is_boundary_edge(e) { 1.0 } else { 2.0 };
            factor * tri.edge_length(e)
        })
        .sum();
    let mut max_interior_defect = 0.0f64;
    let mut boundary_turning = 0.0;
    for w in 0..n_vertices {
        if wedge_on_boundary[w] {
            boundary_turning += PI - wedge_angle[w];
        } else {
            max_interior_defect = max_interior_defect.max((TAU - wedge_angle[w]).abs());
        }
    }
    // Unfold across the edges that are not cut.
    let mut layout = vec![[[0.0; 2]; 3]; n_faces];
    let mut placed = vec![false; n_faces];
    let mut n_pieces = 0;
    for root in 0..n_faces {
        if placed[root] {
            continue;
        }
        n_pieces += 1;
        layout[root] = triangle_layout(tri.face_lengths(root)).expect("a triangle");
        placed[root] = true;
        let mut stack = vec![root];
        while let Some(f) = stack.pop() {
            for k in 0..3 {
                let h = 3 * f + k;
                if cut[tri.edge(h)] {
                    continue;
                }
                let Some(t) = tri.twin(h) else { continue };
                let g = t / 3;
                if placed[g] {
                    continue;
                }
                // The twin runs from the head of h to its tail; its face lies on its left.
                let (a, b) = (layout[f][(k + 1) % 3], layout[f][k]);
                let kt = t % 3;
                let l = tri.face_lengths(g);
                let c = unfold_apex(a, b, l[kt], l[(kt + 1) % 3], l[(kt + 2) % 3])
                    .expect("a triangle");
                let mut corners = [[0.0; 2]; 3];
                corners[kt] = a;
                corners[(kt + 1) % 3] = b;
                corners[(kt + 2) % 3] = c;
                layout[g] = corners;
                placed[g] = true;
                stack.push(g);
            }
        }
    }
    // Disagreement between the places of one vertex of the cut surface.
    let mut first: Vec<Option<[f64; 2]>> = vec![None; n_vertices];
    let mut mismatch = 0.0f64;
    for h in 0..n_h {
        let w = corner_vertex[h] as usize;
        let p = layout[h / 3][h % 3];
        match first[w] {
            None => first[w] = Some(p),
            Some(q) => mismatch = mismatch.max((p[0] - q[0]).hypot(p[1] - q[1])),
        }
    }
    let mean = tri.edge_lengths().iter().sum::<f64>() / n_edges.max(1) as f64;
    let connected = n_pieces == 1;
    Ok(CutReport {
        n_vertices,
        n_edges: n_cut_edges,
        n_faces,
        euler_characteristic,
        connected,
        is_disk: connected && euler_characteristic == 1,
        boundary_length,
        max_interior_defect,
        boundary_turning,
        layout,
        max_layout_mismatch: mismatch / mean,
    })
}
