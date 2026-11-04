//! A bounding volume hierarchy over triangles with closest-hit ray casting (`truth.rs` has no
//! ray queries; its BVH-free closest-point code is not reusable for rays).
//!
//! Median split on the longest centroid axis, leaves of at most 4 triangles, slab test for boxes,
//! Möller-Trumbore for triangles in `f64`. Rays that pass exactly through an edge shared by two
//! triangles are accepted by both (barycentric tolerance `1e-12`), so closed meshes have no
//! cracks; `t` is the smallest over all candidates.

use crate::PatError;

/// A ray hit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// Ray parameter (`origin + t·dir`).
    pub t: f64,
    /// Face index.
    pub face: u32,
}

#[derive(Debug, Clone)]
struct Node {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: first index into `order` and count; inner: `start` = left child, `count` = 0 and the
    /// right child is `left + 1` (children are allocated in pairs).
    start: u32,
    count: u32,
}

/// The hierarchy (immutable after construction; queries are `&self`).
#[derive(Debug, Clone)]
pub struct TriBvh {
    verts: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    order: Vec<u32>,
    nodes: Vec<Node>,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl TriBvh {
    /// Build over a triangle mesh (errors name an out-of-range index or a non-finite vertex).
    pub fn new(verts: &[[f64; 3]], faces: &[[u32; 3]]) -> Result<TriBvh, PatError> {
        if let Some(i) = verts.iter().position(|v| !v.iter().all(|c| c.is_finite())) {
            return Err(PatError::NonFinite { index: i, what: "mesh vertex" });
        }
        for (f, t) in faces.iter().enumerate() {
            if let Some(&v) = t.iter().find(|&&v| v as usize >= verts.len()) {
                return Err(PatError::Format { reason: format!("face {f} has vertex index {v}, but the mesh has {} vertices", verts.len()) });
            }
        }
        if faces.is_empty() {
            return Err(PatError::Format { reason: "mesh has no faces".into() });
        }
        let cent: Vec<[f64; 3]> = faces
            .iter()
            .map(|t| {
                let (a, b, c) = (verts[t[0] as usize], verts[t[1] as usize], verts[t[2] as usize]);
                [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0]
            })
            .collect();
        let mut order: Vec<u32> = (0..faces.len() as u32).collect();
        let mut bvh = TriBvh { verts: verts.to_vec(), faces: faces.to_vec(), order: Vec::new(), nodes: Vec::new() };
        bvh.nodes.push(Node { lo: [0.0; 3], hi: [0.0; 3], start: 0, count: 0 });
        bvh.build(0, 0, faces.len(), &mut order, &cent);
        bvh.order = order;
        Ok(bvh)
    }

    fn bounds(&self, order: &[u32], s: usize, e: usize) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &f in &order[s..e] {
            for &v in &self.faces[f as usize] {
                for a in 0..3 {
                    lo[a] = lo[a].min(self.verts[v as usize][a]);
                    hi[a] = hi[a].max(self.verts[v as usize][a]);
                }
            }
        }
        (lo, hi)
    }

    fn build(&mut self, node: usize, s: usize, e: usize, order: &mut [u32], cent: &[[f64; 3]]) {
        let (lo, hi) = self.bounds(order, s, e);
        self.nodes[node].lo = lo;
        self.nodes[node].hi = hi;
        if e - s <= 4 {
            self.nodes[node].start = s as u32;
            self.nodes[node].count = (e - s) as u32;
            return;
        }
        let mut clo = [f64::INFINITY; 3];
        let mut chi = [f64::NEG_INFINITY; 3];
        for &f in &order[s..e] {
            for a in 0..3 {
                clo[a] = clo[a].min(cent[f as usize][a]);
                chi[a] = chi[a].max(cent[f as usize][a]);
            }
        }
        let axis = (0..3).max_by(|&a, &b| (chi[a] - clo[a]).total_cmp(&(chi[b] - clo[b]))).unwrap_or(0);
        let mid = (s + e) / 2;
        order[s..e].sort_by(|&a, &b| cent[a as usize][axis].total_cmp(&cent[b as usize][axis]).then(a.cmp(&b)));
        let left = self.nodes.len();
        self.nodes.push(Node { lo: [0.0; 3], hi: [0.0; 3], start: 0, count: 0 });
        self.nodes.push(Node { lo: [0.0; 3], hi: [0.0; 3], start: 0, count: 0 });
        self.nodes[node].start = left as u32;
        self.nodes[node].count = 0;
        self.build(left, s, mid, order, cent);
        self.build(left + 1, mid, e, order, cent);
    }

    /// The vertices.
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.verts
    }

    /// The faces.
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// Unit normal of a face (right-hand rule; zero for a degenerate face).
    pub fn face_normal(&self, f: u32) -> [f64; 3] {
        let t = self.faces[f as usize];
        let (a, b, c) = (self.verts[t[0] as usize], self.verts[t[1] as usize], self.verts[t[2] as usize]);
        let n = cross(sub(b, a), sub(c, a));
        let l = dot(n, n).sqrt();
        if l > 0.0 { [n[0] / l, n[1] / l, n[2] / l] } else { [0.0; 3] }
    }

    fn hit_triangle(&self, f: u32, o: [f64; 3], d: [f64; 3]) -> Option<f64> {
        let t = self.faces[f as usize];
        let (a, b, c) = (self.verts[t[0] as usize], self.verts[t[1] as usize], self.verts[t[2] as usize]);
        let (e1, e2) = (sub(b, a), sub(c, a));
        let p = cross(d, e2);
        let det = dot(e1, p);
        if det.abs() < 1e-300 {
            return None;
        }
        let inv = 1.0 / det;
        let s = sub(o, a);
        let u = dot(s, p) * inv;
        const TOL: f64 = 1e-12;
        if !(-TOL..=1.0 + TOL).contains(&u) {
            return None;
        }
        let q = cross(s, e1);
        let v = dot(d, q) * inv;
        if v < -TOL || u + v > 1.0 + TOL {
            return None;
        }
        let t = dot(e2, q) * inv;
        if t > 0.0 { Some(t) } else { None }
    }

    /// The closest hit with `t ∈ (0, ∞)` along `origin + t·dir` (`dir` need not be unit; `t` is in
    /// units of `dir`), or `None`.
    pub fn intersect(&self, o: [f64; 3], d: [f64; 3]) -> Option<Hit> {
        let inv = [1.0 / d[0], 1.0 / d[1], 1.0 / d[2]];
        let mut best: Option<Hit> = None;
        let mut stack = [0u32; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let n = &self.nodes[stack[sp] as usize];
            // slab test
            let (mut t0, mut t1) = (0.0f64, best.map_or(f64::INFINITY, |h| h.t));
            let mut miss = false;
            for a in 0..3 {
                let (mut ta, mut tb) = ((n.lo[a] - o[a]) * inv[a], (n.hi[a] - o[a]) * inv[a]);
                if ta > tb {
                    std::mem::swap(&mut ta, &mut tb);
                }
                if ta.is_nan() || tb.is_nan() {
                    // ray parallel to the slab and on its boundary: treat as inside
                    if o[a] < n.lo[a] || o[a] > n.hi[a] {
                        miss = true;
                    }
                    continue;
                }
                t0 = t0.max(ta);
                t1 = t1.min(tb);
            }
            if miss || t0 > t1 * (1.0 + 1e-12) {
                continue;
            }
            if n.count > 0 {
                for &f in &self.order[n.start as usize..(n.start + n.count) as usize] {
                    if let Some(t) = self.hit_triangle(f, o, d)
                        && best.is_none_or(|h| t < h.t || (t == h.t && f < h.face))
                    {
                        best = Some(Hit { t, face: f });
                    }
                }
            } else if sp + 2 <= stack.len() {
                stack[sp] = n.start;
                stack[sp + 1] = n.start + 1;
                sp += 2;
            }
        }
        best
    }
}
