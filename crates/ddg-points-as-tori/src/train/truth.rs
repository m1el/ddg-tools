//! Ground-truth signed distance for training queries (SPEC.md §10.3): the exact distance to the
//! training mesh through a small BVH, signed by the generalized winding number summed exactly over
//! all faces (negative inside), as the reference (`signed_distance_to_mesh`,
//! `training/mesh_sampling_jax.py`: brute-force distance, `gwn > π` ⇒ inside).
//!
//! The reference computes in f32; [`GwnPrecision::ReferenceF32`] repeats its sign test in f32 for
//! parity, the default is f64.

/// Precision of the winding-number sign test.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GwnPrecision {
    /// f64 (default).
    #[default]
    F64,
    /// The reference's f32 arithmetic (`2·atan2` per face summed in f32, `> π`).
    ReferenceF32,
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn segment_d2(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let e = sub(b, a);
    let l = dot(e, e);
    let t = if l > 0.0 { (dot(sub(p, a), e) / l).clamp(0.0, 1.0) } else { 0.0 };
    let c = [a[0] + t * e[0], a[1] + t * e[1], a[2] + t * e[2]];
    let d = sub(p, c);
    dot(d, d)
}

/// Squared distance from `p` to the triangle `(a, b, c)` (Ericson, *Real-Time Collision
/// Detection*, §5.1.5; degenerate triangles fall back to their edges).
pub fn point_triangle_d2(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let n = cross(ab, ac);
    let nn = dot(n, n);
    let scale = dot(ab, ab).max(dot(ac, ac));
    if nn.is_nan() || nn <= 1e-24 * scale * scale {
        return segment_d2(p, a, b).min(segment_d2(p, b, c)).min(segment_d2(p, c, a));
    }
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    let q = if d1 <= 0.0 && d2 <= 0.0 {
        a
    } else {
        let bp = sub(p, b);
        let d3 = dot(ab, bp);
        let d4 = dot(ac, bp);
        if d3 >= 0.0 && d4 <= d3 {
            b
        } else {
            let vc = d1 * d4 - d3 * d2;
            if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
                let v = d1 / (d1 - d3);
                [a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]]
            } else {
                let cp = sub(p, c);
                let d5 = dot(ab, cp);
                let d6 = dot(ac, cp);
                if d6 >= 0.0 && d5 <= d6 {
                    c
                } else {
                    let vb = d5 * d2 - d1 * d6;
                    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
                        let w = d2 / (d2 - d6);
                        [a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]]
                    } else {
                        let va = d3 * d6 - d5 * d4;
                        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
                            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                            [b[0] + w * (c[0] - b[0]), b[1] + w * (c[1] - b[1]), b[2] + w * (c[2] - b[2])]
                        } else {
                            let denom = 1.0 / (va + vb + vc);
                            let v = vb * denom;
                            let w = vc * denom;
                            [a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w]
                        }
                    }
                }
            }
        }
    };
    let d = sub(p, q);
    dot(d, d)
}

/// Solid angle of the triangle seen from `q` (Van Oosterom–Strackee), as the reference's
/// `gwn_single_face`: `2·atan2(det[a,b,c], abc + (a·b)c + (b·c)a + (c·a)b)` with `a = v0 − q`, ….
pub fn solid_angle(q: [f64; 3], v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let a = sub(v0, q);
    let b = sub(v1, q);
    let c = sub(v2, q);
    let (la, lb, lc) = (dot(a, a).sqrt(), dot(b, b).sqrt(), dot(c, c).sqrt());
    let num = dot(a, cross(b, c));
    let den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    2.0 * num.atan2(den)
}

fn solid_angle_f32(q: [f64; 3], v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f32 {
    let f = |v: [f64; 3]| [v[0] as f32, v[1] as f32, v[2] as f32];
    let (q, v0, v1, v2) = (f(q), f(v0), f(v1), f(v2));
    let s = |x: [f32; 3], y: [f32; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let d = |x: [f32; 3], y: [f32; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let (a, b, c) = (s(v0, q), s(v1, q), s(v2, q));
    let (la, lb, lc) = (d(a, a).sqrt(), d(b, b).sqrt(), d(c, c).sqrt());
    // the reference's expanded determinant
    let num = -a[2] * b[1] * c[0] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1] - a[0] * b[2] * c[1] - a[1] * b[0] * c[2] + a[0] * b[1] * c[2];
    let den = la * lb * lc + d(a, b) * lc + d(b, c) * la + d(c, a) * lb;
    2.0 * num.atan2(den)
}

#[derive(Debug, Clone, Copy)]
struct Node {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: `start..start+count` into `order`; inner: children `left`, `left + 1`… stored as
    /// `start` = left child, `count` = 0.
    start: u32,
    count: u32,
    right: u32,
}

/// A triangle mesh with a BVH for closest-point queries.
#[derive(Debug, Clone)]
pub struct MeshTruth {
    vertices: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    order: Vec<u32>,
    nodes: Vec<Node>,
}

const LEAF: usize = 4;

impl MeshTruth {
    /// Build from vertices and faces (indices must be in range; [`super::data::load_file`]
    /// checks them).
    pub fn new(vertices: &[[f64; 3]], faces: &[[u32; 3]]) -> MeshTruth {
        let faces: Vec<[u32; 3]> = faces.iter().copied().filter(|f| f.iter().all(|&j| (j as usize) < vertices.len())).collect();
        let mut t = MeshTruth { vertices: vertices.to_vec(), faces, order: Vec::new(), nodes: Vec::new() };
        t.order = (0..t.faces.len() as u32).collect();
        if !t.faces.is_empty() {
            let cent: Vec<[f64; 3]> = t
                .faces
                .iter()
                .map(|f| {
                    let [a, b, c] = f.map(|j| t.vertices[j as usize]);
                    [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0]
                })
                .collect();
            let n = t.order.len();
            let mut order = std::mem::take(&mut t.order);
            t.build(&mut order, &cent, 0, n);
            t.order = order;
        }
        t
    }

    /// Number of faces.
    pub fn len(&self) -> usize {
        self.faces.len()
    }

    /// True without faces.
    pub fn is_empty(&self) -> bool {
        self.faces.is_empty()
    }

    fn bounds(&self, order: &[u32]) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &f in order {
            for &j in &self.faces[f as usize] {
                let v = self.vertices[j as usize];
                for a in 0..3 {
                    lo[a] = lo[a].min(v[a]);
                    hi[a] = hi[a].max(v[a]);
                }
            }
        }
        (lo, hi)
    }

    fn build(&mut self, order: &mut [u32], cent: &[[f64; 3]], start: usize, end: usize) -> u32 {
        let (lo, hi) = self.bounds(&order[start..end]);
        let id = self.nodes.len() as u32;
        self.nodes.push(Node { lo, hi, start: start as u32, count: (end - start) as u32, right: 0 });
        if end - start <= LEAF {
            return id;
        }
        let ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
        let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] { 0 } else if ext[1] >= ext[2] { 1 } else { 2 };
        let mid = (start + end) / 2;
        order[start..end].sort_unstable_by(|&x, &y| {
            cent[x as usize][axis].partial_cmp(&cent[y as usize][axis]).unwrap_or(std::cmp::Ordering::Equal).then(x.cmp(&y))
        });
        let left = self.build(order, cent, start, mid);
        let right = self.build(order, cent, mid, end);
        debug_assert_eq!(left, id + 1);
        let node = &mut self.nodes[id as usize];
        node.count = 0;
        node.start = left;
        node.right = right;
        id
    }

    fn box_d2(n: &Node, p: [f64; 3]) -> f64 {
        (0..3).map(|a| (n.lo[a] - p[a]).max(0.0).max(p[a] - n.hi[a])).map(|e| e * e).sum()
    }

    /// Unsigned distance to the mesh (`+∞` without faces).
    pub fn distance(&self, p: [f64; 3]) -> f64 {
        if self.nodes.is_empty() {
            return f64::INFINITY;
        }
        let mut best = f64::INFINITY;
        let mut stack: Vec<(u32, f64)> = Vec::with_capacity(64);
        stack.push((0, Self::box_d2(&self.nodes[0], p)));
        while let Some((id, bd)) = stack.pop() {
            if bd >= best {
                continue;
            }
            let n = self.nodes[id as usize];
            if n.count > 0 {
                for &f in &self.order[n.start as usize..(n.start + n.count) as usize] {
                    let [a, b, c] = self.faces[f as usize].map(|j| self.vertices[j as usize]);
                    best = best.min(point_triangle_d2(p, a, b, c));
                }
            } else {
                let (l, r) = (n.start, n.right);
                let (dl, dr) = (Self::box_d2(&self.nodes[l as usize], p), Self::box_d2(&self.nodes[r as usize], p));
                // nearer child last (popped first)
                if dl <= dr {
                    stack.push((r, dr));
                    stack.push((l, dl));
                } else {
                    stack.push((l, dl));
                    stack.push((r, dr));
                }
            }
        }
        best.sqrt()
    }

    /// Sum of the faces' solid angles at `p` (4π inside a closed outward-oriented mesh).
    pub fn solid_angle_sum(&self, p: [f64; 3], precision: GwnPrecision) -> f64 {
        match precision {
            GwnPrecision::F64 => self
                .faces
                .iter()
                .map(|f| {
                    let [a, b, c] = f.map(|j| self.vertices[j as usize]);
                    solid_angle(p, a, b, c)
                })
                .sum(),
            GwnPrecision::ReferenceF32 => {
                let mut s = 0f32;
                for f in &self.faces {
                    let [a, b, c] = f.map(|j| self.vertices[j as usize]);
                    s += solid_angle_f32(p, a, b, c);
                }
                s as f64
            }
        }
    }

    /// Signed distance: negative where the solid-angle sum exceeds `π` (winding number > ¼), as
    /// the reference.
    pub fn signed_distance(&self, p: [f64; 3], precision: GwnPrecision) -> f64 {
        let d = self.distance(p);
        let inside = match precision {
            GwnPrecision::F64 => self.solid_angle_sum(p, precision) > std::f64::consts::PI,
            GwnPrecision::ReferenceF32 => (self.solid_angle_sum(p, precision) as f32) > std::f32::consts::PI,
        };
        if inside { -d } else { d }
    }
}
