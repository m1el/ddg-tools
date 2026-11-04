//! Fast winding numbers (Barill, Dickson, Schmidt, Levin, Jacobson 2018, Sec. 3.2–3.3,
//! Algorithm 1; SPEC.md §9). Written from the paper, not from libigl.
//!
//! The elements are area-weighted dipoles: points `(p_i, A_i n_i)` or triangles. A tree splits
//! them into clusters; each node stores its expansion center `p̃` (the area-weighted mean
//! position), its radius `r̃` (the largest distance from `p̃` to an element: a point or a triangle
//! corner) and the moments of its dipoles about `p̃`:
//!
//! `ñ = Σ m`, `M = Σ d ⊗ m`, `T = Σ d ⊗ d ⊗ m`, with `d = p − p̃`, `m = A n̂`
//!
//! (Eq. 13–15; for triangles the integrals of Eq. 18–20, computed exactly by the edge-midpoint
//! rule, which integrates quadratics exactly: a triangle is three dipoles `m/3` at its edge
//! midpoints). With `Y = x − p̃` and `G(y) = 1/(4π|y|)` the cluster's winding number is
//!
//! `w̃(x) = ñ·∇G(Y) − M : ∇²G(Y) + ½ T ⋮ ∇³G(Y)`
//!
//! (the order-0, 1, 2 Taylor terms of `m·∇G(x − p)` in `d`). A query uses `w̃` of a node when
//! `‖x − p̃‖ > β r̃` and otherwise recurses; leaves sum exactly (dipoles, or triangle solid
//! angles by van Oosterom–Strackee). Barill's `β`: 2 for triangles, 2.3 for points; the
//! reference (libigl) and SPEC §9 use 2 for both.
//!
//! **Trees.** Points: an octree on the bounding cube of the cloud, split while a cell holds more
//! than `leaf_size` points (Barill: 2), children in octant order, depth at most 48 (exact
//! duplicates). Triangles: an AABB tree, split at the median centroid of the widest centroid
//! axis (ties by triangle index), leaves of at most `leaf_size` (8) triangles. Construction is
//! serial and deterministic; queries run on the explicit pool, each in a fixed traversal order.

use super::{CHUNK, check_len, check_queries, cross, dot, norm, sub};
use crate::{PatError, PointCloud};
use rayon::prelude::*;
use std::f64::consts::PI;

/// Options of [`FastWinding`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FwnOptions {
    /// Accuracy parameter: a node is expanded when `‖x − p̃‖ > β r̃` (SPEC: 2).
    pub beta: f64,
    /// Taylor order: 0, 1 or 2 (SPEC: 2).
    pub order: u8,
    /// Elements per leaf (`None`: 2 for points as Barill, 8 for triangles).
    pub leaf_size: Option<usize>,
    /// Threads of batch evaluation (`None`: all cores).
    pub threads: Option<usize>,
}

impl Default for FwnOptions {
    fn default() -> Self {
        FwnOptions { beta: 2.0, order: 2, leaf_size: None, threads: None }
    }
}

/// What a batch evaluation did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FwnReport {
    /// Queries.
    pub queries: usize,
    /// Node expansions used.
    pub far_nodes: u64,
    /// Elements summed exactly.
    pub near_elements: u64,
    /// Dipole terms with `p = x` left out.
    pub skipped: u64,
    /// Wall time.
    pub seconds: f64,
}

#[derive(Debug, Clone)]
enum Elements {
    /// Positions and dipole vectors `A n`.
    Points(Vec<[f64; 3]>, Vec<[f64; 3]>),
    /// Triangle corners.
    Triangles(Vec<[[f64; 3]; 3]>),
}

#[derive(Debug, Clone, Default)]
struct Node {
    center: [f64; 3],
    radius: f64,
    n0: [f64; 3],
    m1: [[f64; 3]; 3],
    t2: [[[f64; 3]; 3]; 3],
    first_child: u32,
    n_children: u32,
    start: u32,
    end: u32,
}

/// Barill's hierarchical winding number of a point cloud or a triangle soup.
#[derive(Debug)]
pub struct FastWinding {
    elems: Elements,
    /// Original index of each element in tree order.
    original: Vec<u32>,
    nodes: Vec<Node>,
    children: Vec<u32>,
    opts: FwnOptions,
    pool: rayon::ThreadPool,
}

fn check_opts(o: &FwnOptions) -> Result<(), PatError> {
    if !(o.beta.is_finite() && o.beta > 0.0) {
        return Err(PatError::InvalidOption { what: format!("FwnOptions::beta = {} (must be > 0)", o.beta) });
    }
    if o.order > 2 {
        return Err(PatError::InvalidOption { what: format!("FwnOptions::order = {} (0, 1 or 2)", o.order) });
    }
    if o.leaf_size == Some(0) {
        return Err(PatError::InvalidOption { what: "FwnOptions::leaf_size = Some(0)".into() });
    }
    Ok(())
}

/// The dipole term `m·(p − x)/(4π‖p − x‖³)` (0 and `true` if `p = x`).
///
/// [formula: pat-winding-dipole]
#[inline]
pub fn dipole(p: [f64; 3], m: [f64; 3], x: [f64; 3]) -> (f64, bool) {
    let d = sub(p, x);
    let r2 = dot(d, d);
    if r2 == 0.0 {
        return (0.0, true);
    }
    (dot(m, d) / (4.0 * PI * r2 * r2.sqrt()), false)
}

/// The solid angle of triangle `abc` seen from `x`, divided by `4π` (van Oosterom–Strackee:
/// `Ω = 2 atan2(det[a,b,c], |a||b||c| + (a·b)|c| + (b·c)|a| + (c·a)|b|)`, corners relative to
/// `x`). Positive when `x` is on the side the normal `(b−a)×(c−a)` points away from.
///
/// [formula: pat-solid-angle]
#[inline]
pub fn triangle_winding(t: &[[f64; 3]; 3], x: [f64; 3]) -> f64 {
    let a = sub(t[0], x);
    let b = sub(t[1], x);
    let c = sub(t[2], x);
    let (la, lb, lc) = (norm(a), norm(b), norm(c));
    let det = dot(a, cross(b, c));
    let den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
    2.0 * det.atan2(den) / (4.0 * PI)
}

/// The far-field value `ñ·∇G(Y) − M:∇²G(Y) + ½ T⋮∇³G(Y)` of a node up to `order`, `Y = x − p̃ ≠ 0`.
///
/// [formula: pat-fwn-expansion]
fn expansion(nd: &Node, x: [f64; 3], order: u8) -> f64 {
    let y = sub(x, nd.center);
    let r2 = dot(y, y);
    let r = r2.sqrt();
    let c = 1.0 / (4.0 * PI * r2 * r);
    // order 0: ñ·∇G, ∇G = −Y/(4πr³)
    let mut w = -c * dot(nd.n0, y);
    if order >= 1 {
        // −M_jk ∂_j∂_k G, ∂_j∂_k G = (3Y_jY_k/r² − δ_jk)/(4πr³)
        let mut s = 0.0;
        let mut tr = 0.0;
        for j in 0..3 {
            tr += nd.m1[j][j];
            for k in 0..3 {
                s += nd.m1[j][k] * y[j] * y[k];
            }
        }
        w -= c * (3.0 * s / r2 - tr);
    }
    if order >= 2 {
        // ½ T_jlk ∂_j∂_l∂_k G, ∂_j∂_l∂_k G = (−15 Y_jY_lY_k/r⁴ + 3(δ_jl Y_k + δ_jk Y_l + δ_lk Y_j)/r²)/(4πr³)
        let mut cubic = 0.0;
        let mut lin = 0.0;
        for j in 0..3 {
            for l in 0..3 {
                for k in 0..3 {
                    cubic += nd.t2[j][l][k] * y[j] * y[l] * y[k];
                }
            }
            for (k, &yk) in y.iter().enumerate() {
                lin += (nd.t2[j][j][k] + nd.t2[j][k][j] + nd.t2[k][j][j]) * yk;
            }
        }
        w += 0.5 * c * (-15.0 * cubic / (r2 * r2) + 3.0 * lin / r2);
    }
    w
}

/// Sub-dipoles of an element for the moments: (position, dipole vector, weight for the center).
fn for_each_dipole(elems: &Elements, e: usize, mut f: impl FnMut([f64; 3], [f64; 3], f64)) {
    match elems {
        Elements::Points(p, m) => f(p[e], m[e], norm(m[e])),
        Elements::Triangles(t) => {
            let [a, b, c] = t[e];
            let m = cross(sub(b, a), sub(c, a));
            let m3 = [m[0] / 6.0, m[1] / 6.0, m[2] / 6.0];
            let w = norm(m) / 6.0;
            for (u, v) in [(a, b), (b, c), (c, a)] {
                f([0.5 * (u[0] + v[0]), 0.5 * (u[1] + v[1]), 0.5 * (u[2] + v[2])], m3, w);
            }
        }
    }
}

fn node_moments(elems: &Elements, start: usize, end: usize) -> Node {
    let mut nd = Node { start: start as u32, end: end as u32, ..Default::default() };
    let mut wsum = 0.0;
    let mut c = [0.0; 3];
    let mut mean = [0.0; 3];
    let mut cnt = 0.0;
    for e in start..end {
        for_each_dipole(elems, e, |p, _, w| {
            for t in 0..3 {
                c[t] += w * p[t];
                mean[t] += p[t];
            }
            wsum += w;
            cnt += 1.0;
        });
    }
    nd.center = if wsum > 0.0 { [c[0] / wsum, c[1] / wsum, c[2] / wsum] } else { [mean[0] / cnt, mean[1] / cnt, mean[2] / cnt] };
    let mut rad: f64 = 0.0;
    for e in start..end {
        match elems {
            Elements::Points(p, _) => rad = rad.max(norm(sub(p[e], nd.center))),
            Elements::Triangles(t) => {
                for v in t[e] {
                    rad = rad.max(norm(sub(v, nd.center)));
                }
            }
        }
        for_each_dipole(elems, e, |p, m, _| {
            let d = sub(p, nd.center);
            for j in 0..3 {
                nd.n0[j] += m[j];
                for (k, &mk) in m.iter().enumerate() {
                    nd.m1[j][k] += d[j] * mk;
                    for l in 0..3 {
                        nd.t2[j][l][k] += d[j] * d[l] * mk;
                    }
                }
            }
        });
    }
    nd.radius = rad;
    nd
}

struct Builder<'a> {
    nodes: Vec<Node>,
    children: Vec<u32>,
    leaf: usize,
    keys: &'a [[f64; 3]],
}

impl Builder<'_> {
    /// Octree over `perm[start..end]` in the cube `(center, half)`; returns the node index.
    fn octree(&mut self, perm: &mut [u32], start: usize, end: usize, center: [f64; 3], half: f64, depth: u32) -> u32 {
        let id = self.nodes.len();
        self.nodes.push(Node::default());
        let mut kids = Vec::new();
        if end - start > self.leaf && depth < 48 && half > 0.0 {
            let code = |i: u32| {
                let p = self.keys[i as usize];
                (p[0] >= center[0]) as usize | ((p[1] >= center[1]) as usize) << 1 | ((p[2] >= center[2]) as usize) << 2
            };
            let mut buckets: [Vec<u32>; 8] = Default::default();
            for &i in &perm[start..end] {
                buckets[code(i)].push(i);
            }
            let mut s = start;
            let mut ranges = Vec::new();
            for (o, b) in buckets.iter().enumerate() {
                perm[s..s + b.len()].copy_from_slice(b);
                if !b.is_empty() {
                    ranges.push((o, s, s + b.len()));
                }
                s += b.len();
            }
            let h = 0.5 * half;
            for (o, s, e) in ranges {
                let c = [
                    center[0] + if o & 1 != 0 { h } else { -h },
                    center[1] + if o & 2 != 0 { h } else { -h },
                    center[2] + if o & 4 != 0 { h } else { -h },
                ];
                kids.push(self.octree(perm, s, e, c, h, depth + 1));
            }
        }
        self.finish(id, kids, start, end)
    }

    /// Median-split AABB tree over `perm[start..end]` (keys: centroids).
    fn median(&mut self, perm: &mut [u32], start: usize, end: usize) -> u32 {
        let id = self.nodes.len();
        self.nodes.push(Node::default());
        let mut kids = Vec::new();
        if end - start > self.leaf {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for &i in &perm[start..end] {
                let k = self.keys[i as usize];
                for t in 0..3 {
                    lo[t] = lo[t].min(k[t]);
                    hi[t] = hi[t].max(k[t]);
                }
            }
            let ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
            let ax = if ext[0] >= ext[1] && ext[0] >= ext[2] { 0 } else if ext[1] >= ext[2] { 1 } else { 2 };
            let mid = start + (end - start) / 2;
            let keys = self.keys;
            perm[start..end].select_nth_unstable_by(mid - start, |&a, &b| {
                keys[a as usize][ax].total_cmp(&keys[b as usize][ax]).then(a.cmp(&b))
            });
            kids.push(self.median(perm, start, mid));
            kids.push(self.median(perm, mid, end));
        }
        self.finish(id, kids, start, end)
    }

    /// Record the ranges and children (the moments are computed on the reordered elements).
    fn finish(&mut self, id: usize, kids: Vec<u32>, start: usize, end: usize) -> u32 {
        let nd = &mut self.nodes[id];
        nd.start = start as u32;
        nd.end = end as u32;
        nd.first_child = self.children.len() as u32;
        nd.n_children = kids.len() as u32;
        self.children.extend(kids);
        id as u32
    }
}

impl FastWinding {
    /// The tree of an oriented point cloud with point areas `areas` (e.g. [`super::point_areas`]).
    pub fn from_cloud(cloud: &PointCloud, areas: &[f64], opts: &FwnOptions) -> Result<Self, PatError> {
        check_opts(opts)?;
        let n = cloud.len();
        if areas.len() != n {
            return Err(PatError::SizeMismatch { what: format!("{} areas for {n} points", areas.len()) });
        }
        if let Some(i) = areas.iter().position(|a| !a.is_finite()) {
            return Err(PatError::NonFinite { index: i, what: "point area" });
        }
        if n == 0 {
            return Err(PatError::TooFewPoints { n, need: 1 });
        }
        let p = cloud.positions();
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for q in p {
            for t in 0..3 {
                lo[t] = lo[t].min(q[t]);
                hi[t] = hi[t].max(q[t]);
            }
        }
        let center = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
        let half = 0.5 * (hi[0] - lo[0]).max(hi[1] - lo[1]).max(hi[2] - lo[2]);
        let mut perm: Vec<u32> = (0..n as u32).collect();
        let mut b = Builder { nodes: Vec::new(), children: Vec::new(), leaf: opts.leaf_size.unwrap_or(2), keys: p };
        b.octree(&mut perm, 0, n, center, half, 0);
        let pp: Vec<[f64; 3]> = perm.iter().map(|&i| p[i as usize]).collect();
        let mm: Vec<[f64; 3]> = perm
            .iter()
            .map(|&i| {
                let (a, nn) = (areas[i as usize], cloud.normals()[i as usize]);
                [a * nn[0], a * nn[1], a * nn[2]]
            })
            .collect();
        let elems = Elements::Points(pp, mm);
        Self::assemble(b.nodes, b.children, elems, perm, opts)
    }

    /// The tree of a triangle soup (`faces` index `vertices`; degenerate triangles are allowed and
    /// contribute nothing).
    pub fn from_triangles(vertices: &[[f64; 3]], faces: &[[u32; 3]], opts: &FwnOptions) -> Result<Self, PatError> {
        check_opts(opts)?;
        if let Some(i) = vertices.iter().position(|v| !v.iter().all(|c| c.is_finite())) {
            return Err(PatError::NonFinite { index: i, what: "mesh vertex" });
        }
        if faces.is_empty() {
            return Err(PatError::Format { reason: "mesh has no faces".into() });
        }
        for (f, t) in faces.iter().enumerate() {
            if let Some(&v) = t.iter().find(|&&v| v as usize >= vertices.len()) {
                return Err(PatError::Format {
                    reason: format!("face {f} has vertex index {v}, but the mesh has {} vertices", vertices.len()),
                });
            }
        }
        let tris: Vec<[[f64; 3]; 3]> =
            faces.iter().map(|t| [vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]]).collect();
        let keys: Vec<[f64; 3]> = tris
            .iter()
            .map(|t| {
                [(t[0][0] + t[1][0] + t[2][0]) / 3.0, (t[0][1] + t[1][1] + t[2][1]) / 3.0, (t[0][2] + t[1][2] + t[2][2]) / 3.0]
            })
            .collect();
        let n = tris.len();
        let mut perm: Vec<u32> = (0..n as u32).collect();
        let mut b = Builder { nodes: Vec::new(), children: Vec::new(), leaf: opts.leaf_size.unwrap_or(8), keys: &keys };
        b.median(&mut perm, 0, n);
        let elems = Elements::Triangles(perm.iter().map(|&i| tris[i as usize]).collect());
        Self::assemble(b.nodes, b.children, elems, perm, opts)
    }

    fn assemble(mut nodes: Vec<Node>, children: Vec<u32>, elems: Elements, perm: Vec<u32>, opts: &FwnOptions) -> Result<Self, PatError> {
        // moments on the reordered elements (the builder only fixed the ranges)
        for nd in nodes.iter_mut() {
            let m = node_moments(&elems, nd.start as usize, nd.end as usize);
            *nd = Node { first_child: nd.first_child, n_children: nd.n_children, ..m };
        }
        let pool = crate::pool::make_pool(opts.threads)?;
        Ok(FastWinding { elems, original: perm, nodes, children, opts: *opts, pool })
    }

    /// Nodes of the tree.
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Original element index of each element in tree order.
    pub fn tree_order(&self) -> &[u32] {
        &self.original
    }

    /// `(w, far nodes, near elements, skipped)` at `x` (finite).
    fn query(&self, x: [f64; 3], stack: &mut Vec<u32>) -> (f64, u32, u32, u32) {
        let mut w = 0.0;
        let (mut far, mut near, mut skipped) = (0, 0, 0);
        stack.clear();
        stack.push(0);
        while let Some(id) = stack.pop() {
            let nd = &self.nodes[id as usize];
            let dist = norm(sub(x, nd.center));
            if dist > self.opts.beta * nd.radius {
                w += expansion(nd, x, self.opts.order);
                far += 1;
            } else if nd.n_children == 0 {
                for e in nd.start as usize..nd.end as usize {
                    match &self.elems {
                        Elements::Points(p, m) => {
                            let (v, s) = dipole(p[e], m[e], x);
                            w += v;
                            skipped += s as u32;
                        }
                        Elements::Triangles(t) => w += triangle_winding(&t[e], x),
                    }
                    near += 1;
                }
            } else {
                let c = nd.first_child as usize;
                for k in (0..nd.n_children as usize).rev() {
                    stack.push(self.children[c + k]);
                }
            }
        }
        (w, far, near, skipped)
    }

    /// The approximate winding number at `x`.
    pub fn winding_number(&self, x: [f64; 3]) -> Result<f64, PatError> {
        check_queries(&[x])?;
        Ok(self.query(x, &mut Vec::new()).0)
    }

    /// The approximate winding number at every query, on the tree's pool.
    pub fn evaluate(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<FwnReport, PatError> {
        let t0 = std::time::Instant::now();
        check_queries(xs)?;
        check_len("winding", out.len(), xs.len())?;
        let parts: Vec<(u64, u64, u64)> = self.pool.install(|| {
            out.par_chunks_mut(CHUNK)
                .zip(xs.par_chunks(CHUNK))
                .map(|(o, q)| {
                    let mut stack = Vec::new();
                    let mut acc = (0u64, 0u64, 0u64);
                    for (o, &x) in o.iter_mut().zip(q) {
                        let (w, f, n, s) = self.query(x, &mut stack);
                        *o = w;
                        acc.0 += f as u64;
                        acc.1 += n as u64;
                        acc.2 += s as u64;
                    }
                    acc
                })
                .collect()
        });
        let mut rep = FwnReport { queries: xs.len(), ..Default::default() };
        for (f, n, s) in parts {
            rep.far_nodes += f;
            rep.near_elements += n;
            rep.skipped += s;
        }
        rep.seconds = t0.elapsed().as_secs_f64();
        Ok(rep)
    }
}
