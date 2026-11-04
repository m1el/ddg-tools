//! A kd-tree over f64 points with exact, deterministically ordered kNN and radius queries
//! (SPEC.md §3.1).
//!
//! - Each node splits its points at the median of the widest axis of their bounding box (ties
//!   between axes go to the lower axis), with the median found by `select_nth_unstable_by` on
//!   `(coordinate, index)`. The tree's shape depends only on `n` and the bucket size; the
//!   permutation only on the points. Construction is parallel over subtrees and deterministic.
//! - Every node stores the tight bounding box of its points; queries prune a node when the
//!   squared distance to its box is *strictly* larger than the current k-th candidate, so that a
//!   point at equal distance with a smaller index is never missed.
//! - Results are ordered by `(distance², index)`.

use crate::PatError;

/// A neighbor found by a query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Neighbor {
    /// Squared Euclidean distance to the query.
    pub d2: f64,
    /// Index of the point in the input order.
    pub index: u32,
}

#[inline]
fn less(a: (f64, u32), b: (f64, u32)) -> bool {
    a.0 < b.0 || (a.0 == b.0 && a.1 < b.1)
}

#[derive(Debug, Clone, Copy)]
struct Node {
    lo: [f64; 3],
    hi: [f64; 3],
    start: u32,
    end: u32,
    /// Index of the right child; the left child is the next node. `u32::MAX` for a leaf.
    right: u32,
}

/// Options of [`KdTree::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KdTreeOptions {
    /// Maximum number of points in a leaf (default 16).
    pub bucket: usize,
}

impl Default for KdTreeOptions {
    fn default() -> Self {
        KdTreeOptions { bucket: 16 }
    }
}

/// The kd-tree. Immutable after construction; queries are `&self` and thread safe.
#[derive(Debug, Clone)]
pub struct KdTree {
    /// Points in tree order.
    pts: Vec<[f64; 3]>,
    /// Tree order → input index.
    idx: Vec<u32>,
    nodes: Vec<Node>,
}

fn node_count(m: usize, bucket: usize) -> usize {
    if m <= bucket { 1 } else { 1 + node_count(m / 2, bucket) + node_count(m - m / 2, bucket) }
}

const PAR_MIN: usize = 4096;

fn build_rec(points: &[[f64; 3]], idx: &mut [u32], nodes: &mut [Node], start: usize, first_node: usize, bucket: usize) {
    let m = idx.len();
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &i in idx.iter() {
        let p = points[i as usize];
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    let (head, rest) = nodes.split_first_mut().expect("node slice sized by node_count");
    *head = Node { lo, hi, start: start as u32, end: (start + m) as u32, right: u32::MAX };
    if m <= bucket {
        return;
    }
    let mut axis = 0;
    for a in 1..3 {
        if hi[a] - lo[a] > hi[axis] - lo[axis] {
            axis = a;
        }
    }
    let half = m / 2;
    idx.select_nth_unstable_by(half, |&a, &b| {
        let (ca, cb) = (points[a as usize][axis], points[b as usize][axis]);
        ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b))
    });
    let nl = node_count(half, bucket);
    head.right = (first_node + 1 + nl) as u32;
    let (il, ir) = idx.split_at_mut(half);
    let (nlft, nrgt) = rest.split_at_mut(nl);
    if m >= PAR_MIN {
        rayon::join(
            || build_rec(points, il, nlft, start, first_node + 1, bucket),
            || build_rec(points, ir, nrgt, start + half, first_node + 1 + nl, bucket),
        );
    } else {
        build_rec(points, il, nlft, start, first_node + 1, bucket);
        build_rec(points, ir, nrgt, start + half, first_node + 1 + nl, bucket);
    }
}

#[inline]
fn box_d2(n: &Node, q: &[f64; 3]) -> f64 {
    let mut s = 0.0;
    for a in 0..3 {
        let d = if q[a] < n.lo[a] {
            n.lo[a] - q[a]
        } else if q[a] > n.hi[a] {
            q[a] - n.hi[a]
        } else {
            0.0
        };
        s += d * d;
    }
    s
}

#[inline]
fn dist2(p: &[f64; 3], q: &[f64; 3]) -> f64 {
    let dx = p[0] - q[0];
    let dy = p[1] - q[1];
    let dz = p[2] - q[2];
    dx * dx + dy * dy + dz * dz
}

impl KdTree {
    /// Build the tree. Runs on the current rayon pool (call inside `pool.install` to choose it);
    /// the result does not depend on the number of threads.
    pub fn new(points: &[[f64; 3]], opts: KdTreeOptions) -> Result<Self, PatError> {
        if opts.bucket == 0 {
            return Err(PatError::InvalidOption { what: "KdTreeOptions::bucket = 0".into() });
        }
        if points.is_empty() {
            return Err(PatError::TooFewPoints { n: 0, need: 1 });
        }
        if points.len() > u32::MAX as usize - 1 {
            return Err(PatError::InvalidOption { what: format!("{} points exceed u32 indices", points.len()) });
        }
        for (i, p) in points.iter().enumerate() {
            if !p.iter().all(|c| c.is_finite()) {
                return Err(PatError::NonFinite { index: i, what: "point" });
            }
        }
        let n = points.len();
        let mut idx: Vec<u32> = (0..n as u32).collect();
        let nn = node_count(n, opts.bucket);
        let empty = Node { lo: [0.0; 3], hi: [0.0; 3], start: 0, end: 0, right: u32::MAX };
        let mut nodes = vec![empty; nn];
        build_rec(points, &mut idx, &mut nodes, 0, 0, opts.bucket);
        let pts = idx.iter().map(|&i| points[i as usize]).collect();
        Ok(KdTree { pts, idx, nodes })
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.pts.len()
    }

    /// Always false (a tree has at least one point).
    pub fn is_empty(&self) -> bool {
        self.pts.is_empty()
    }

    /// The `k` nearest points of `q`, ordered by `(d², index)`, written to `out` (cleared first).
    /// `k ≥ n` returns all points. Exact.
    pub fn knn(&self, q: [f64; 3], k: usize, out: &mut Vec<Neighbor>) {
        out.clear();
        let k = k.min(self.pts.len());
        if k == 0 {
            return;
        }
        let mut worst = f64::INFINITY;
        let mut stack: [(u32, f64); 64] = [(0, 0.0); 64];
        let mut sp = 1usize;
        stack[0] = (0, box_d2(&self.nodes[0], &q));
        while sp > 0 {
            sp -= 1;
            let (ni, bd) = stack[sp];
            if out.len() == k && bd > worst {
                continue;
            }
            let node = &self.nodes[ni as usize];
            if node.right == u32::MAX {
                for t in node.start as usize..node.end as usize {
                    let d2 = dist2(&self.pts[t], &q);
                    let id = self.idx[t];
                    if out.len() < k {
                        let pos = out.partition_point(|nb| less((nb.d2, nb.index), (d2, id)));
                        out.insert(pos, Neighbor { d2, index: id });
                    } else if less((d2, id), (worst, out[k - 1].index)) {
                        let pos = out.partition_point(|nb| less((nb.d2, nb.index), (d2, id)));
                        out.insert(pos, Neighbor { d2, index: id });
                        out.pop();
                    } else {
                        continue;
                    }
                    if out.len() == k {
                        worst = out[k - 1].d2;
                    }
                }
            } else {
                let l = ni + 1;
                let r = node.right;
                let dl = box_d2(&self.nodes[l as usize], &q);
                let dr = box_d2(&self.nodes[r as usize], &q);
                // push the farther child first so the nearer one is visited first
                let (a, da, b, db) = if dl <= dr { (r, dr, l, dl) } else { (l, dl, r, dr) };
                if !(out.len() == k && da > worst) {
                    stack[sp] = (a, da);
                    sp += 1;
                }
                if !(out.len() == k && db > worst) {
                    stack[sp] = (b, db);
                    sp += 1;
                }
            }
        }
    }

    /// All points with `d² ≤ r²`, ordered by `(d², index)`, written to `out` (cleared first).
    pub fn radius(&self, q: [f64; 3], r: f64, out: &mut Vec<Neighbor>) {
        out.clear();
        let r2 = r * r;
        let mut stack: [u32; 64] = [0; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let ni = stack[sp];
            let node = &self.nodes[ni as usize];
            if box_d2(node, &q) > r2 {
                continue;
            }
            if node.right == u32::MAX {
                for t in node.start as usize..node.end as usize {
                    let d2 = dist2(&self.pts[t], &q);
                    if d2 <= r2 {
                        out.push(Neighbor { d2, index: self.idx[t] });
                    }
                }
            } else {
                stack[sp] = node.right;
                stack[sp + 1] = ni + 1; // the left child follows its parent in pre-order
                sp += 2;
            }
        }
        out.sort_unstable_by(|a, b| {
            a.d2.partial_cmp(&b.d2).unwrap_or(std::cmp::Ordering::Equal).then(a.index.cmp(&b.index))
        });
    }
}
