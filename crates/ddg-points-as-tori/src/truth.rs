//! Ground truth for meshes (SPEC.md §9, §16.2): the exact signed distance to a triangle mesh.
//!
//! * **Distance.** The exact distance to the nearest triangle, by an AABB tree over the triangles
//!   (median split of the centroids on the widest axis, leaves of 4) and the closest point of a
//!   triangle by the Voronoi-region test of Ericson, *Real-Time Collision Detection* (2004),
//!   §5.1.5. Degenerate (zero-area) triangles are measured by their three edges. Ties between
//!   triangles go to the smaller face index, so the closest face is deterministic.
//! * **Sign.** The generalized winding number (Jacobson, Kavan, Sorkine-Hornung 2013): the sum of
//!   the triangles' solid angles over `4π` (van Oosterom–Strackee), exact by direct summation up
//!   to [`AUTO_EXACT_FACES`] faces and by Barill's fast winding number
//!   ([`crate::baselines::FastWinding`]) beyond (SPEC §9). Inside iff `w > threshold`; the
//!   threshold is `½` (smooth theory) or [`JAX_THRESHOLD`] `= ¼` (the reference's training code:
//!   `Σ Ω > π`).
//! * `φ(x) = −d(x)` inside, `+d(x)` outside (0 on the surface).
//!
//! The paper computes the distance with FCPW and the sign with the GWN (Sec. 5.2); the result is
//! the same function. Every query is independent and runs in a fixed traversal order, so batch
//! results are bitwise independent of the thread count.
//!
//! [`check_mesh`] reports boundary, non-manifold and inconsistently oriented edges, degenerate
//! faces and the signed volume. With [`TruthOptions::require_closed`] (default) a mesh that is not
//! closed, manifold, consistently oriented and of positive volume is rejected; reconstructed
//! meshes (SPSR) may be accepted as soups.

use crate::PatError;
use crate::baselines::fwn::triangle_winding;
use crate::baselines::{FastWinding, FwnOptions};
use rayon::prelude::*;

/// Faces up to which [`SignMode::Auto`] sums the solid angles exactly (SPEC §9: 20k).
pub const AUTO_EXACT_FACES: usize = 20_000;

/// The reference's training threshold `Σ Ω > π`, i.e. `w > ¼`.
pub const JAX_THRESHOLD: f64 = 0.25;

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn axpy(a: [f64; 3], s: f64, d: [f64; 3]) -> [f64; 3] {
    [a[0] + s * d[0], a[1] + s * d[1], a[2] + s * d[2]]
}
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// The closest point of the segment `ab` to `p`.
pub fn closest_point_on_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    if l2 == 0.0 {
        return a;
    }
    let t = (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0);
    axpy(a, t, ab)
}

/// The closest point of the triangle `abc` to `p` (Ericson 2004, §5.1.5: the Voronoi regions of
/// the corners, then of the edges, then the interior). A degenerate triangle is measured by its
/// edges.
///
/// [formula: pat-truth-closest-point]
pub fn closest_point_on_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 && d1 - d3 > 0.0 {
        return axpy(a, d1 / (d1 - d3), ab);
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 && d2 - d6 > 0.0 {
        return axpy(a, d2 / (d2 - d6), ac);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 && (d4 - d3) + (d5 - d6) > 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return axpy(b, w, sub(c, b));
    }
    let den = va + vb + vc;
    if den > 0.0 && den.is_finite() {
        let v = vb / den;
        let w = vc / den;
        return axpy(axpy(a, v, ab), w, ac);
    }
    // degenerate: the nearest of the three edges
    let mut best = closest_point_on_segment(p, a, b);
    let mut bd = dot(sub(p, best), sub(p, best));
    for (u, v) in [(b, c), (c, a)] {
        let q = closest_point_on_segment(p, u, v);
        let d = dot(sub(p, q), sub(p, q));
        if d < bd {
            bd = d;
            best = q;
        }
    }
    best
}

/// Topology and orientation of a triangle mesh.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MeshReport {
    /// Vertices.
    pub vertices: usize,
    /// Faces.
    pub faces: usize,
    /// Faces of zero area.
    pub degenerate_faces: usize,
    /// Edges with one face.
    pub boundary_edges: usize,
    /// Edges with more than two faces.
    pub nonmanifold_edges: usize,
    /// Edges with two faces that traverse it in the same direction.
    pub inconsistent_edges: usize,
    /// `Σ det[a, b, c]/6`: positive for a closed outward-oriented mesh.
    pub signed_volume: f64,
    /// The first offending edge `(u, v)` with its faces, for messages.
    pub first_problem: Option<((u32, u32), Vec<u32>)>,
}

impl MeshReport {
    /// Closed, manifold and consistently oriented.
    pub fn is_closed(&self) -> bool {
        self.boundary_edges == 0 && self.nonmanifold_edges == 0 && self.inconsistent_edges == 0
    }
}

/// Validate indices and coordinates and report the topology (errors name the element).
pub fn check_mesh(vertices: &[[f64; 3]], faces: &[[u32; 3]]) -> Result<MeshReport, PatError> {
    if let Some(i) = vertices.iter().position(|v| !v.iter().all(|c| c.is_finite())) {
        return Err(PatError::NonFinite { index: i, what: "mesh vertex" });
    }
    if faces.is_empty() {
        return Err(PatError::Format { reason: "mesh has no faces".into() });
    }
    let nv = vertices.len();
    for (f, t) in faces.iter().enumerate() {
        if let Some(&v) = t.iter().find(|&&v| v as usize >= nv) {
            return Err(PatError::Format { reason: format!("face {f} has vertex index {v}, but the mesh has {nv} vertices") });
        }
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            return Err(PatError::Format { reason: format!("face {f} repeats a vertex: {:?}", t) });
        }
    }
    let mut rep = MeshReport { vertices: nv, faces: faces.len(), ..Default::default() };
    // directed edges (min, max, forward?, face)
    let mut e: Vec<(u32, u32, bool, u32)> = Vec::with_capacity(3 * faces.len());
    for (f, t) in faces.iter().enumerate() {
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            e.push((u.min(v), u.max(v), u < v, f as u32));
        }
        let (a, b, c) = (vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]);
        let cr = cross(sub(b, a), sub(c, a));
        if dot(cr, cr) == 0.0 {
            rep.degenerate_faces += 1;
        }
        rep.signed_volume += dot(a, cross(b, c)) / 6.0;
    }
    e.sort_unstable();
    let mut i = 0;
    while i < e.len() {
        let mut j = i + 1;
        while j < e.len() && e[j].0 == e[i].0 && e[j].1 == e[i].1 {
            j += 1;
        }
        let problem = match j - i {
            1 => {
                rep.boundary_edges += 1;
                true
            }
            2 if e[i].2 == e[i + 1].2 => {
                rep.inconsistent_edges += 1;
                true
            }
            2 => false,
            _ => {
                rep.nonmanifold_edges += 1;
                true
            }
        };
        if problem && rep.first_problem.is_none() {
            rep.first_problem = Some(((e[i].0, e[i].1), e[i..j].iter().map(|x| x.3).collect()));
        }
        i = j;
    }
    Ok(rep)
}

/// How the sign is computed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SignMode {
    /// Direct sum of all solid angles.
    Exact,
    /// Barill's fast winding number with accuracy `beta`.
    Fast {
        /// Accuracy parameter (2).
        beta: f64,
    },
    /// Exact up to [`AUTO_EXACT_FACES`] faces, otherwise `Fast { beta: 2 }`.
    Auto,
}

/// Options of [`MeshTruth`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TruthOptions {
    /// Sign computation.
    pub sign: SignMode,
    /// Inside iff `w > threshold` (½; [`JAX_THRESHOLD`] for training parity).
    pub threshold: f64,
    /// Reject meshes that are not closed, manifold, consistently oriented and of positive volume.
    pub require_closed: bool,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
}

impl Default for TruthOptions {
    fn default() -> Self {
        TruthOptions { sign: SignMode::Auto, threshold: 0.5, require_closed: true, threads: None }
    }
}

/// What a batch evaluation did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TruthReport {
    /// Queries.
    pub queries: usize,
    /// Queries classified inside.
    pub inside: usize,
    /// Queries on the surface (distance 0).
    pub on_surface: usize,
    /// Queries whose winding number is within 0.1 of the threshold (ambiguous sign).
    pub ambiguous: usize,
    /// Wall time.
    pub seconds: f64,
}

#[derive(Debug, Clone, Copy)]
struct BNode {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: `start..end` into the ordered faces; inner: children `left`, `left + 1`… stored as
    /// `left` and `right`.
    left: u32,
    right: u32,
    start: u32,
    end: u32,
}

/// The exact signed distance to a triangle mesh.
#[derive(Debug)]
pub struct MeshTruth {
    /// Triangles in BVH order.
    tris: Vec<[[f64; 3]; 3]>,
    /// Original face index of each.
    face: Vec<u32>,
    nodes: Vec<BNode>,
    fwn: Option<FastWinding>,
    opts: TruthOptions,
    pool: rayon::ThreadPool,
}

const LEAF: usize = 4;
const NONE: u32 = u32::MAX;

fn box_d2(n: &BNode, p: [f64; 3]) -> f64 {
    let mut s = 0.0;
    for (k, &pk) in p.iter().enumerate() {
        let d = (n.lo[k] - pk).max(0.0).max(pk - n.hi[k]);
        s += d * d;
    }
    s
}

impl MeshTruth {
    /// Build the BVH (and the FWN tree if the sign mode needs it).
    pub fn new(vertices: &[[f64; 3]], faces: &[[u32; 3]], opts: &TruthOptions) -> Result<(Self, MeshReport), PatError> {
        if !opts.threshold.is_finite() {
            return Err(PatError::InvalidOption { what: format!("TruthOptions::threshold = {}", opts.threshold) });
        }
        let rep = check_mesh(vertices, faces)?;
        if opts.require_closed {
            if let Some(((u, v), fs)) = &rep.first_problem {
                return Err(PatError::Format {
                    reason: format!(
                        "mesh is not closed, manifold and consistently oriented: edge ({u}, {v}) has faces {fs:?} \
                         ({} boundary, {} non-manifold, {} inconsistent edges)",
                        rep.boundary_edges, rep.nonmanifold_edges, rep.inconsistent_edges
                    ),
                });
            }
            if rep.signed_volume <= 0.0 {
                return Err(PatError::Format {
                    reason: format!("mesh encloses signed volume {} ≤ 0 (inward orientation?)", rep.signed_volume),
                });
            }
        }
        let tri_all: Vec<[[f64; 3]; 3]> =
            faces.iter().map(|t| [vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]]).collect();
        let cent: Vec<[f64; 3]> = tri_all
            .iter()
            .map(|t| [(t[0][0] + t[1][0] + t[2][0]) / 3.0, (t[0][1] + t[1][1] + t[2][1]) / 3.0, (t[0][2] + t[1][2] + t[2][2]) / 3.0])
            .collect();
        let mut perm: Vec<u32> = (0..faces.len() as u32).collect();
        let mut nodes = Vec::new();
        build(&mut nodes, &mut perm, 0, faces.len(), &cent, &tri_all);
        let tris = perm.iter().map(|&f| tri_all[f as usize]).collect();
        let beta = match opts.sign {
            SignMode::Exact => None,
            SignMode::Fast { beta } => Some(beta),
            SignMode::Auto => (faces.len() > AUTO_EXACT_FACES).then_some(2.0),
        };
        let fwn = match beta {
            Some(beta) => Some(FastWinding::from_triangles(
                vertices,
                faces,
                &FwnOptions { beta, threads: Some(1), ..Default::default() },
            )?),
            None => None,
        };
        let pool = crate::pool::make_pool(opts.threads)?;
        Ok((MeshTruth { tris, face: perm, nodes, fwn, opts: *opts, pool }, rep))
    }

    /// `(distance, closest point, face)` of `x`.
    pub fn closest(&self, x: [f64; 3]) -> (f64, [f64; 3], u32) {
        let mut best = (f64::INFINITY, [0.0; 3], NONE);
        let mut stack = vec![0u32];
        while let Some(id) = stack.pop() {
            let n = &self.nodes[id as usize];
            if box_d2(n, x) > best.0 {
                continue;
            }
            if n.left == NONE {
                for t in n.start as usize..n.end as usize {
                    let [a, b, c] = self.tris[t];
                    let q = closest_point_on_triangle(x, a, b, c);
                    let d = sub(x, q);
                    let d2 = dot(d, d);
                    let f = self.face[t];
                    if d2 < best.0 || (d2 == best.0 && f < best.2) {
                        best = (d2, q, f);
                    }
                }
            } else {
                let (l, r) = (n.left, n.right);
                let (dl, dr) = (box_d2(&self.nodes[l as usize], x), box_d2(&self.nodes[r as usize], x));
                // visit the nearer child first (pushed last)
                if dl <= dr {
                    stack.push(r);
                    stack.push(l);
                } else {
                    stack.push(l);
                    stack.push(r);
                }
            }
        }
        (best.0.sqrt(), best.1, best.2)
    }

    /// The generalized winding number of the mesh at `x` (exact or fast per the options).
    pub fn winding(&self, x: [f64; 3]) -> f64 {
        match &self.fwn {
            Some(f) => f.winding_number(x).unwrap_or(f64::NAN),
            None => {
                let mut w = 0.0;
                for t in &self.tris {
                    w += triangle_winding(t, x);
                }
                w
            }
        }
    }

    /// `(φ, w)` at a finite `x`.
    fn eval(&self, x: [f64; 3]) -> (f64, f64) {
        let (d, _, _) = self.closest(x);
        if d == 0.0 {
            return (0.0, f64::NAN);
        }
        let w = self.winding(x);
        (if w > self.opts.threshold { -d } else { d }, w)
    }

    /// The signed distance at `x`.
    pub fn signed_distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        crate::baselines::check_queries(&[x])?;
        Ok(self.eval(x).0)
    }

    /// The signed distance (and optionally the winding number; NaN on the surface) at every query.
    pub fn evaluate(&self, xs: &[[f64; 3]], sdf: &mut [f64], mut winding: Option<&mut [f64]>) -> Result<TruthReport, PatError> {
        let t0 = std::time::Instant::now();
        crate::baselines::check_queries(xs)?;
        crate::baselines::check_len("sdf", sdf.len(), xs.len())?;
        if let Some(w) = winding.as_deref() {
            crate::baselines::check_len("winding", w.len(), xs.len())?;
        }
        let vals: Vec<(f64, f64)> = self.pool.install(|| xs.par_iter().with_min_len(64).map(|&x| self.eval(x)).collect());
        let mut rep = TruthReport { queries: xs.len(), ..Default::default() };
        for (i, &(phi, w)) in vals.iter().enumerate() {
            sdf[i] = phi;
            if let Some(o) = winding.as_deref_mut() {
                o[i] = w;
            }
            if phi == 0.0 {
                rep.on_surface += 1;
            } else {
                rep.inside += (phi < 0.0) as usize;
                rep.ambiguous += ((w - self.opts.threshold).abs() < 0.1) as usize;
            }
        }
        rep.seconds = t0.elapsed().as_secs_f64();
        Ok(rep)
    }
}

fn build(nodes: &mut Vec<BNode>, perm: &mut [u32], start: usize, end: usize, cent: &[[f64; 3]], tris: &[[[f64; 3]; 3]]) -> u32 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let mut clo = [f64::INFINITY; 3];
    let mut chi = [f64::NEG_INFINITY; 3];
    for &f in &perm[start..end] {
        for v in tris[f as usize] {
            for k in 0..3 {
                lo[k] = lo[k].min(v[k]);
                hi[k] = hi[k].max(v[k]);
            }
        }
        let c = cent[f as usize];
        for k in 0..3 {
            clo[k] = clo[k].min(c[k]);
            chi[k] = chi[k].max(c[k]);
        }
    }
    let id = nodes.len();
    nodes.push(BNode { lo, hi, left: NONE, right: NONE, start: start as u32, end: end as u32 });
    if end - start > LEAF {
        let ext = [chi[0] - clo[0], chi[1] - clo[1], chi[2] - clo[2]];
        let ax = if ext[0] >= ext[1] && ext[0] >= ext[2] { 0 } else if ext[1] >= ext[2] { 1 } else { 2 };
        let mid = start + (end - start) / 2;
        perm[start..end]
            .select_nth_unstable_by(mid - start, |&a, &b| cent[a as usize][ax].total_cmp(&cent[b as usize][ax]).then(a.cmp(&b)));
        let l = build(nodes, perm, start, mid, cent, tris);
        let r = build(nodes, perm, mid, end, cent, tris);
        nodes[id].left = l;
        nodes[id].right = r;
    }
    id as u32
}
