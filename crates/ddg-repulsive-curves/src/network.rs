//! Curve networks: validated topology, components, edge quantities and the geometric checks of
//! `SPEC.md` §2.

use crate::error::RcError;
use crate::kernel::v3::{dot, norm, scale, sub};

/// Relative tolerance of [`RcError::EdgesIntersect`]: two disjoint edges closer than
/// `INTERSECTION_TOLERANCE * (mean edge length)` are rejected.
pub const INTERSECTION_TOLERANCE: f64 = 1e-12;

/// A coordinate axis; `Planar { axis }` mode of `SPEC.md` §2.2 (D18) holds one coordinate fixed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// The x coordinate.
    X,
    /// The y coordinate.
    Y,
    /// The z coordinate.
    Z,
}

impl Axis {
    /// The coordinate index 0, 1 or 2.
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }
}

/// Checks that every vertex has the same value of the coordinate `axis` as vertex 0 (planar
/// mode, `SPEC.md` §2.2). Exact comparison: the input must be exactly planar.
///
/// # Errors
/// [`RcError::NotPlanar`] with the first vertex that differs.
pub fn check_planar(positions: &[[f64; 3]], axis: Axis) -> Result<(), RcError> {
    let k = axis.index();
    if let Some(first) = positions.first() {
        for (v, p) in positions.iter().enumerate() {
            if p[k] != first[k] {
                return Err(RcError::NotPlanar {
                    vertex: v as u32,
                    axis: k,
                    value: p[k],
                });
            }
        }
    }
    Ok(())
}

/// A connected component of a curve network.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component {
    /// Vertex indices, increasing.
    pub vertices: Vec<u32>,
    /// Edge indices, increasing.
    pub edges: Vec<u32>,
    /// `true` if every vertex of the component has degree 2 (a closed loop).
    pub closed: bool,
}

/// Per-edge quantities of a geometry (`SPEC.md` §2.1).
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeGeometry {
    /// `ℓ_I`.
    pub length: Vec<f64>,
    /// `T_I`, unit tangent from the first to the second endpoint.
    pub tangent: Vec<[f64; 3]>,
    /// `x_I`, the midpoint.
    pub midpoint: Vec<[f64; 3]>,
}

/// A validated curve network: the graph `G = (V, E)` of `SPEC.md` §2.1 without positions.
///
/// Guarantees: every edge joins two distinct vertices in range, no two edges join the same
/// vertex pair, every vertex has degree `>= 1`, and at least one pair of edges is disjoint
/// (`𝒫 ≠ ∅`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CurveNetwork {
    n_vertices: usize,
    edges: Vec<[u32; 2]>,
    /// CSR incidence: edges of vertex `v` are `incident[offsets[v]..offsets[v + 1]]`, increasing.
    offsets: Vec<usize>,
    incident: Vec<u32>,
    components: Vec<Component>,
    component_of_vertex: Vec<u32>,
}

impl CurveNetwork {
    /// Builds and validates the topology (`SPEC.md` §2.3, the part that does not need positions).
    ///
    /// # Errors
    /// [`RcError::VertexOutOfRange`], [`RcError::SelfLoop`], [`RcError::DuplicateEdge`],
    /// [`RcError::IsolatedVertex`], [`RcError::NoInteractingPairs`], in this order of checking
    /// (the first offending element in index order is reported).
    pub fn new(n_vertices: usize, edges: Vec<[u32; 2]>) -> Result<Self, RcError> {
        for (i, &[a, b]) in edges.iter().enumerate() {
            for v in [a, b] {
                if v as usize >= n_vertices {
                    return Err(RcError::VertexOutOfRange {
                        edge: i,
                        vertex: v,
                        n_vertices,
                    });
                }
            }
            if a == b {
                return Err(RcError::SelfLoop { edge: i, vertex: a });
            }
        }
        // Duplicates: sort the (min, max, edge) triples.
        let mut keyed: Vec<(u32, u32, usize)> = edges
            .iter()
            .enumerate()
            .map(|(i, &[a, b])| (a.min(b), a.max(b), i))
            .collect();
        keyed.sort_unstable();
        for w in keyed.windows(2) {
            if w[0].0 == w[1].0 && w[0].1 == w[1].1 {
                return Err(RcError::DuplicateEdge {
                    a: w[0].0,
                    b: w[0].1,
                    edges: [w[0].2, w[1].2],
                });
            }
        }
        // CSR incidence.
        let mut offsets = vec![0usize; n_vertices + 1];
        for &[a, b] in &edges {
            offsets[a as usize + 1] += 1;
            offsets[b as usize + 1] += 1;
        }
        for v in 0..n_vertices {
            offsets[v + 1] += offsets[v];
        }
        let mut fill = offsets.clone();
        let mut incident = vec![0u32; 2 * edges.len()];
        for (i, &[a, b]) in edges.iter().enumerate() {
            for v in [a, b] {
                incident[fill[v as usize]] = i as u32;
                fill[v as usize] += 1;
            }
        }
        for v in 0..n_vertices {
            if offsets[v] == offsets[v + 1] {
                return Err(RcError::IsolatedVertex { vertex: v as u32 });
            }
        }
        // 𝒫 is empty iff every pair of distinct edges is adjacent (no duplicates, so two edges
        // share at most one vertex).
        let m = edges.len();
        let adjacent_pairs: usize = (0..n_vertices)
            .map(|v| {
                let d = offsets[v + 1] - offsets[v];
                d * (d - 1) / 2
            })
            .sum();
        if adjacent_pairs == m * m.saturating_sub(1) / 2 {
            return Err(RcError::NoInteractingPairs { n_edges: m });
        }
        // Components by flood fill, in increasing order of the smallest vertex.
        let mut component_of_vertex = vec![u32::MAX; n_vertices];
        let mut components = Vec::new();
        let mut stack = Vec::new();
        for s in 0..n_vertices {
            if component_of_vertex[s] != u32::MAX {
                continue;
            }
            let c = components.len() as u32;
            let mut vertices = Vec::new();
            component_of_vertex[s] = c;
            stack.push(s as u32);
            while let Some(v) = stack.pop() {
                vertices.push(v);
                for &e in &incident[offsets[v as usize]..offsets[v as usize + 1]] {
                    let [a, b] = edges[e as usize];
                    let w = if a == v { b } else { a };
                    if component_of_vertex[w as usize] == u32::MAX {
                        component_of_vertex[w as usize] = c;
                        stack.push(w);
                    }
                }
            }
            vertices.sort_unstable();
            let mut comp_edges: Vec<u32> = vertices
                .iter()
                .flat_map(|&v| {
                    incident[offsets[v as usize]..offsets[v as usize + 1]]
                        .iter()
                        .copied()
                })
                .collect();
            comp_edges.sort_unstable();
            comp_edges.dedup();
            let closed = vertices
                .iter()
                .all(|&v| offsets[v as usize + 1] - offsets[v as usize] == 2);
            components.push(Component {
                vertices,
                edges: comp_edges,
                closed,
            });
        }
        Ok(CurveNetwork {
            n_vertices,
            edges,
            offsets,
            incident,
            components,
            component_of_vertex,
        })
    }

    /// Number of vertices `|V|`.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// Number of edges `|E|`.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }

    /// The edges `(i₁, i₂)` with the orientation the caller gave.
    pub fn edges(&self) -> &[[u32; 2]] {
        &self.edges
    }

    /// The connected components, ordered by their smallest vertex.
    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// Index of the component that contains vertex `v`.
    pub fn component_of(&self, v: u32) -> usize {
        self.component_of_vertex[v as usize] as usize
    }

    /// Degree of vertex `v`.
    pub fn degree(&self, v: u32) -> usize {
        self.offsets[v as usize + 1] - self.offsets[v as usize]
    }

    /// The edges that contain vertex `v`, in increasing edge index.
    pub fn incident_edges(&self, v: u32) -> &[u32] {
        &self.incident[self.offsets[v as usize]..self.offsets[v as usize + 1]]
    }

    /// `true` if edges `i` and `j` share a vertex (or are the same edge): the pair is not in `𝒫`.
    #[inline]
    pub fn edges_adjacent(&self, i: usize, j: usize) -> bool {
        let [a, b] = self.edges[i];
        let [c, d] = self.edges[j];
        a == c || a == d || b == c || b == d
    }

    /// `ℓ_I`, `T_I`, `x_I` for every edge.
    ///
    /// [formula: rc-edge-quantities]
    ///
    /// The tangent of a zero-length edge is NaN; [`CurveNetwork::validate_positions`] excludes
    /// that case.
    pub fn edge_geometry(&self, positions: &[[f64; 3]]) -> EdgeGeometry {
        let m = self.edges.len();
        let mut length = Vec::with_capacity(m);
        let mut tangent = Vec::with_capacity(m);
        let mut midpoint = Vec::with_capacity(m);
        for &[a, b] in &self.edges {
            let (pa, pb) = (positions[a as usize], positions[b as usize]);
            let d = sub(pb, pa);
            let l = norm(d);
            length.push(l);
            tangent.push(scale(d, 1.0 / l));
            midpoint.push(scale([pa[0] + pb[0], pa[1] + pb[1], pa[2] + pb[2]], 0.5));
        }
        EdgeGeometry {
            length,
            tangent,
            midpoint,
        }
    }

    /// Dual lengths `m_v = ½ Σ_{I ∋ v} ℓ_I`; `Σ_v m_v = Σ_I ℓ_I`.
    ///
    /// [formula: rc-dual-length]
    pub fn dual_lengths(&self, positions: &[[f64; 3]]) -> Vec<f64> {
        (0..self.n_vertices)
            .map(|v| {
                0.5 * self
                    .incident_edges(v as u32)
                    .iter()
                    .map(|&e| {
                        let [a, b] = self.edges[e as usize];
                        norm(sub(positions[a as usize], positions[b as usize]))
                    })
                    .sum::<f64>()
            })
            .collect()
    }

    /// Mean edge length `ℓ̄₀` (the length unit of the relative tolerances, `SPEC.md` §2.1).
    pub fn mean_edge_length(&self, positions: &[[f64; 3]]) -> f64 {
        let total: f64 = self
            .edges
            .iter()
            .map(|&[a, b]| norm(sub(positions[a as usize], positions[b as usize])))
            .sum();
        total / self.edges.len() as f64
    }

    /// The geometric checks of `SPEC.md` §2.3 (D27), in this order: array length, finite
    /// coordinates, zero-length edges, coincident vertices, intersecting disjoint edges.
    ///
    /// # Errors
    /// [`RcError::PositionCountMismatch`], [`RcError::NonFinitePosition`],
    /// [`RcError::ZeroLengthEdge`], [`RcError::CoincidentVertices`], [`RcError::EdgesIntersect`].
    pub fn validate_positions(&self, positions: &[[f64; 3]]) -> Result<(), RcError> {
        if positions.len() != self.n_vertices {
            return Err(RcError::PositionCountMismatch {
                expected: self.n_vertices,
                found: positions.len(),
            });
        }
        for (v, p) in positions.iter().enumerate() {
            if !p.iter().all(|c| c.is_finite()) {
                return Err(RcError::NonFinitePosition { vertex: v as u32 });
            }
        }
        for (i, &[a, b]) in self.edges.iter().enumerate() {
            if positions[a as usize] == positions[b as usize] {
                return Err(RcError::ZeroLengthEdge { edge: i, a, b });
            }
        }
        // Coincident vertices: sort lexicographically (total order on finite values).
        let mut order: Vec<u32> = (0..self.n_vertices as u32).collect();
        order.sort_unstable_by(|&i, &j| {
            let (p, q) = (positions[i as usize], positions[j as usize]);
            p[0].total_cmp(&q[0])
                .then(p[1].total_cmp(&q[1]))
                .then(p[2].total_cmp(&q[2]))
                .then(i.cmp(&j))
        });
        let mut worst: Option<(u32, u32)> = None;
        for w in order.windows(2) {
            if positions[w[0] as usize] == positions[w[1] as usize] {
                // Report the pair with the smallest first index, deterministic.
                let pair = (w[0].min(w[1]), w[0].max(w[1]));
                if worst.is_none_or(|x| pair < x) {
                    worst = Some(pair);
                }
            }
        }
        if let Some((a, b)) = worst {
            return Err(RcError::CoincidentVertices { a, b });
        }
        let tol = INTERSECTION_TOLERANCE * self.mean_edge_length(positions);
        if let Some((a, b, distance)) = self.first_close_pair(positions, tol) {
            return Err(RcError::EdgesIntersect { a, b, distance });
        }
        Ok(())
    }

    /// The first pair `(I, J)`, `I < J`, of disjoint edges at distance `<= tol`
    /// (lexicographically smallest), by a bounding-box tree broad phase.
    fn first_close_pair(&self, positions: &[[f64; 3]], tol: f64) -> Option<(usize, usize, f64)> {
        let m = self.edges.len();
        let boxes: Vec<([f64; 3], [f64; 3])> = self
            .edges
            .iter()
            .map(|&[a, b]| {
                let (p, q) = (positions[a as usize], positions[b as usize]);
                (
                    [
                        p[0].min(q[0]) - tol,
                        p[1].min(q[1]) - tol,
                        p[2].min(q[2]) - tol,
                    ],
                    [
                        p[0].max(q[0]) + tol,
                        p[1].max(q[1]) + tol,
                        p[2].max(q[2]) + tol,
                    ],
                )
            })
            .collect();
        let tree = BoxTree::build(&boxes);
        let mut hits = Vec::new();
        for i in 0..m {
            let mut best: Option<(usize, f64)> = None;
            hits.clear();
            tree.query(&boxes[i], &mut hits);
            for &j in &hits {
                let j = j as usize;
                if j <= i || self.edges_adjacent(i, j) {
                    continue;
                }
                let [a, b] = self.edges[i];
                let [c, d] = self.edges[j];
                let dist = segment_distance(
                    positions[a as usize],
                    positions[b as usize],
                    positions[c as usize],
                    positions[d as usize],
                );
                if dist <= tol && best.is_none_or(|(bj, _)| j < bj) {
                    best = Some((j, dist));
                }
            }
            if let Some((j, dist)) = best {
                return Some((i, j, dist));
            }
        }
        None
    }
}

/// A median-split bounding-box tree over a list of boxes (broad phase of the validation only;
/// the BVH of `SPEC.md` §8.1 is `spatial::bvh`).
struct BoxTree {
    nodes: Vec<BoxNode>,
    perm: Vec<u32>,
}

struct BoxNode {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: range `start..end` of `perm`; inner: children `start`, `end` are node indices.
    start: u32,
    end: u32,
    leaf: bool,
}

impl BoxTree {
    const LEAF: usize = 4;

    fn build(boxes: &[([f64; 3], [f64; 3])]) -> BoxTree {
        let mut perm: Vec<u32> = (0..boxes.len() as u32).collect();
        let mut nodes = Vec::new();
        Self::build_rec(boxes, &mut perm, 0, boxes.len(), &mut nodes);
        BoxTree { nodes, perm }
    }

    fn build_rec(
        boxes: &[([f64; 3], [f64; 3])],
        perm: &mut [u32],
        start: usize,
        end: usize,
        nodes: &mut Vec<BoxNode>,
    ) -> u32 {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &i in &perm[start..end] {
            let (l, h) = boxes[i as usize];
            for k in 0..3 {
                lo[k] = lo[k].min(l[k]);
                hi[k] = hi[k].max(h[k]);
            }
        }
        let id = nodes.len() as u32;
        nodes.push(BoxNode {
            lo,
            hi,
            start: start as u32,
            end: end as u32,
            leaf: true,
        });
        if end - start > Self::LEAF {
            let ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
            let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                0
            } else if ext[1] >= ext[2] {
                1
            } else {
                2
            };
            let centre = |i: u32| boxes[i as usize].0[axis] + boxes[i as usize].1[axis];
            perm[start..end]
                .sort_unstable_by(|&a, &b| centre(a).total_cmp(&centre(b)).then(a.cmp(&b)));
            let mid = (start + end) / 2;
            let l = Self::build_rec(boxes, perm, start, mid, nodes);
            let r = Self::build_rec(boxes, perm, mid, end, nodes);
            let n = &mut nodes[id as usize];
            n.leaf = false;
            n.start = l;
            n.end = r;
        }
        id
    }

    fn query(&self, b: &([f64; 3], [f64; 3]), out: &mut Vec<u32>) {
        let mut stack = vec![0u32];
        while let Some(id) = stack.pop() {
            let n = &self.nodes[id as usize];
            if (0..3).any(|k| n.hi[k] < b.0[k] || n.lo[k] > b.1[k]) {
                continue;
            }
            if n.leaf {
                out.extend_from_slice(&self.perm[n.start as usize..n.end as usize]);
            } else {
                stack.push(n.start);
                stack.push(n.end);
            }
        }
    }
}

/// Distance between the segments `p1 q1` and `p2 q2` (closest points of two segments,
/// Ericson, *Real-Time Collision Detection*, §5.1.9).
pub(crate) fn segment_distance(p1: [f64; 3], q1: [f64; 3], p2: [f64; 3], q2: [f64; 3]) -> f64 {
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let r = sub(p1, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);
    let (s, t);
    if a <= 0.0 && e <= 0.0 {
        s = 0.0;
        t = 0.0;
    } else if a <= 0.0 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= 0.0 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            let mut s0 = if denom > 0.0 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t0 = (b * s0 + f) / e;
            if t0 < 0.0 {
                t0 = 0.0;
                s0 = (-c / a).clamp(0.0, 1.0);
            } else if t0 > 1.0 {
                t0 = 1.0;
                s0 = ((b - c) / a).clamp(0.0, 1.0);
            }
            s = s0;
            t = t0;
        }
    }
    let c1 = [p1[0] + s * d1[0], p1[1] + s * d1[1], p1[2] + s * d1[2]];
    let c2 = [p2[0] + t * d2[0], p2[1] + t * d2[1], p2[2] + t * d2[2]];
    norm(sub(c1, c2))
}
