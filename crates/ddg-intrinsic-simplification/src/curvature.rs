//! Principal curvatures and directions at the vertices of a triangle mesh, at a chosen scale
//! (SPEC §8, S22), and the conversion of 3D directions to angles in the normalized frames of
//! the vertices.
//!
//! Davide Panozzo, Enrico Puppo, Luigi Rocca, *Efficient Multi-scale Curvature and Crease
//! Estimation*, Workshop on Computer Graphics, Computer Vision and Mathematics (GraVisMa)
//! 2010 (`papers/GraVisMa10-PanozzoPuppoRocca.pdf`), Sec. 3 and 4. Implemented from the
//! paper; the crease extraction of its Sec. 5 is not needed here. Libigl's
//! `principal_curvature` was not read or copied.
//!
//! # The estimator (paper Sec. 4)
//!
//! For every vertex `P`:
//!
//! 1. **Neighbourhood** at scale `r = k · ē` (`ē` the mean edge length, `k` =
//!    [`CurvatureOptions::scale`], the paper's `2e … 16e`). As in the paper: a breadth-first
//!    traversal of the mesh graph starting at `P` that collects the vertices within the
//!    Euclidean ball of radius `r` around `P` and expands only from collected vertices. (It is
//!    the paper's ball, not a `k`-ring: the number of rings adapts to the edge lengths. On a
//!    boundary the traversal stops at the boundary, so the neighbourhood is one-sided.)
//! 2. **Fitting frame**: the `w` axis is the mean of the pseudo-normals of the collected
//!    vertices (pseudo-normal: unit area-weighted average of the normals of the incident
//!    faces, computed once). Collected vertices whose pseudo-normal has a negative scalar
//!    product with `w` are discarded (the paper's rule against portions of surface that bend
//!    away); `P` itself is always kept. The axes `u`, `v` are an arbitrary orthonormal basis of
//!    the plane orthogonal to `w` (the result does not depend on the choice).
//! 3. **Monte-Carlo** ([`CurvatureOptions::max_samples`], paper default 50): if more vertices
//!    remain, a uniformly random subset of that size is fitted (deterministic: seeded by
//!    [`CurvatureOptions::seed`] and the vertex id). `None` fits all.
//! 4. **Quadric fit** (paper Eq. 4): `w = a u² + b v² + c u v + d u + e v` through `P`
//!    (no constant term), five coefficients by linear least squares (Householder QR on
//!    coordinates scaled by `r`). Fewer than 5 vertices or a rank-deficient system:
//!    [`CurvatureStatus::Underdetermined`].
//! 5. **Shape operator** at the origin of the graph `(u, v, f(u, v))` from the first and second
//!    fundamental forms (`E = 1 + d²`, `F = d e`, `G = 1 + e²`, `L = 2a/s`, `M = c/s`,
//!    `N = 2b/s`, `s = √(1 + d² + e²)`; paper Sec. 3): `S = I⁻¹ II`, a 2×2 matrix whose
//!    eigenvalues are the principal curvatures and whose eigenvectors (mapped to 3D with
//!    `t = ξ X_u + η X_v`) the principal directions.
//!
//! # Conventions
//!
//! * The normal is the outward one (counter-clockwise faces); **convex is positive**: a sphere
//!   of radius `ρ` has `k_max = k_min = 1/ρ`. (The paper's `S = −dN` with the outward normal
//!   gives the opposite sign; [`CurvatureOptions::convex_positive`] = `false` selects that.)
//! * `k_max ≥ k_min`; `d_max`, `d_min` are unit tangent vectors, orthogonal to each other and
//!   to [`PrincipalCurvatures::normal`] with `d_min = n × d_max`. They define lines, not
//!   vectors: the sign is arbitrary. At an umbilic point (`k_max ≈ k_min`) the directions are
//!   arbitrary (an orthonormal pair of the tangent plane).
//! * **Boundary vertices** are processed like the others on their one-sided neighbourhood;
//!   the values are less accurate (the fit is an extrapolation), and nothing marks them but
//!   [`PrincipalCurvatures::is_boundary`].
//! * **Degenerate input**: a vertex without a valid normal (all incident faces have zero area, or
//!   the normals cancel, or it has no face) gets [`CurvatureStatus::DegenerateNormal`]; one whose
//!   fit is rank deficient (a flat or collinear neighbourhood in the projection, fewer than 5
//!   points) gets [`CurvatureStatus::Underdetermined`]. In both cases `k_max = k_min = 0`,
//!   the normal is the pseudo-normal (or `+z`) and the directions an orthonormal pair of its
//!   tangent plane; nothing is NaN. A flat region is *not* degenerate: it fits `a = … = e = 0`.
//!
//! # Directions to angles
//!
//! [`direction_angles`] converts 3D directions at the vertices into the normalized intrinsic
//! frame of the vertex (signposts of `ddg_intrinsic::IntrinsicTriangulation` on the input
//! mesh; angles in `[0, 2π)` at an interior vertex, `[0, π]` at a boundary one). The rule is
//! documented there.

#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

use std::f64::consts::{PI, TAU};

use ddg_core::{GeometryError, TriMesh};
use ddg_core::generators::SplitMix64;
use ddg_core::Geometry;
use ddg_intrinsic::{IntrinsicError, IntrinsicTriangulation};

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
/// The unit vector of `a`, or `None` if it is not finite and non-zero.
fn unit(a: V3) -> Option<V3> {
    let n = norm(a);
    (n.is_finite() && n > 0.0).then(|| scale(a, 1.0 / n))
}
/// A unit vector orthogonal to the unit vector `n`.
fn any_orthogonal(n: V3) -> V3 {
    let axis = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    unit(cross(n, axis)).unwrap_or([1.0, 0.0, 0.0])
}

/// Errors of the curvature estimation and of the direction conversion.
#[derive(Debug, Clone, PartialEq)]
pub enum CurvatureError {
    /// Not one position per vertex.
    PositionCount {
        /// Vertices of the mesh.
        expected: usize,
        /// Positions given.
        got: usize,
    },
    /// A position has a non-finite coordinate.
    NonFinitePosition {
        /// The vertex.
        vertex: usize,
    },
    /// An option is out of range.
    BadOption {
        /// What is wrong.
        reason: String,
    },
    /// Not one direction per vertex, or a vertex id out of range.
    BadInput {
        /// What is wrong.
        reason: String,
    },
    /// A direction is zero or not finite.
    BadDirection {
        /// The vertex.
        vertex: usize,
    },
    /// The intrinsic triangulation of the input could not be built.
    Intrinsic(String),
}

impl std::fmt::Display for CurvatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CurvatureError::PositionCount { expected, got } => {
                write!(f, "{got} positions for a mesh with {expected} vertices")
            }
            CurvatureError::NonFinitePosition { vertex } => {
                write!(f, "vertex {vertex} has a non-finite coordinate")
            }
            CurvatureError::BadOption { reason } => write!(f, "bad curvature option: {reason}"),
            CurvatureError::BadInput { reason } => write!(f, "bad input: {reason}"),
            CurvatureError::BadDirection { vertex } => {
                write!(f, "the direction at vertex {vertex} is zero or not finite")
            }
            CurvatureError::Intrinsic(s) => write!(f, "intrinsic triangulation: {s}"),
        }
    }
}

impl std::error::Error for CurvatureError {}

impl From<IntrinsicError> for CurvatureError {
    fn from(e: IntrinsicError) -> Self {
        CurvatureError::Intrinsic(e.to_string())
    }
}
impl From<GeometryError> for CurvatureError {
    fn from(e: GeometryError) -> Self {
        CurvatureError::Intrinsic(e.to_string())
    }
}

/// Options of [`principal_curvatures`].
#[derive(Debug, Clone, PartialEq)]
pub struct CurvatureOptions {
    /// The scale `k`: the neighbourhood is the ball of radius `k ×` the mean edge length.
    /// Default 3 (the paper uses 2 to 16).
    pub scale: f64,
    /// Monte-Carlo cap on the number of fitted vertices (paper: 50). `None` fits all.
    /// Must be at least 5 when given. Default `Some(50)`.
    pub max_samples: Option<usize>,
    /// Seed of the (deterministic) sampling. Default 0.
    pub seed: u64,
    /// Convex is positive (see the module documentation). Default `true`.
    pub convex_positive: bool,
}

impl Default for CurvatureOptions {
    fn default() -> Self {
        CurvatureOptions {
            scale: 3.0,
            max_samples: Some(50),
            seed: 0,
            convex_positive: true,
        }
    }
}

/// How the estimate at one vertex came about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvatureStatus {
    /// A quadric was fitted.
    Fitted,
    /// Fewer than 5 vertices or a rank-deficient fit: curvatures set to 0.
    Underdetermined,
    /// No valid normal at the vertex: curvatures set to 0.
    DegenerateNormal,
}

/// Principal curvatures and directions at every vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct PrincipalCurvatures {
    /// Larger principal curvature.
    pub k_max: Vec<f64>,
    /// Smaller principal curvature.
    pub k_min: Vec<f64>,
    /// Unit direction of `k_max` (a line: the sign is arbitrary).
    pub d_max: Vec<[f64; 3]>,
    /// Unit direction of `k_min`: `normal × d_max`.
    pub d_min: Vec<[f64; 3]>,
    /// The unit fitting normal (the average pseudo-normal of the neighbourhood).
    pub normal: Vec<[f64; 3]>,
    /// Per-vertex status.
    pub status: Vec<CurvatureStatus>,
    /// Number of vertices fitted at each vertex (after the discarding and the sampling).
    pub n_fitted: Vec<u32>,
    /// Whether the vertex is on the boundary of the mesh.
    pub is_boundary: Vec<bool>,
    /// The radius `r = scale · ē` that was used.
    pub radius: f64,
    /// The mean edge length `ē`.
    pub mean_edge_length: f64,
}

impl PrincipalCurvatures {
    /// Number of vertices whose status is not [`CurvatureStatus::Fitted`].
    pub fn n_failed(&self) -> usize {
        self.status
            .iter()
            .filter(|&&s| s != CurvatureStatus::Fitted)
            .count()
    }
}

/// Least-squares solution of `A x = y` (`A` is `m × 5`, row-major) by Householder QR. `None` if
/// the system is rank deficient (relative diagonal of `R` below `1e-10`) or not finite.
fn lstsq5(a: &mut [[f64; 5]], y: &mut [f64]) -> Option<[f64; 5]> {
    let m = a.len();
    if m < 5 {
        return None;
    }
    let mut rdiag = [0.0f64; 5];
    for k in 0..5 {
        let mut nrm = 0.0;
        for i in k..m {
            nrm += a[i][k] * a[i][k];
        }
        let nrm = nrm.sqrt();
        rdiag[k] = nrm;
        if !(nrm > 0.0) {
            return None;
        }
        let alpha = if a[k][k] > 0.0 { -nrm } else { nrm };
        // v = x - alpha e_k, stored in column k below the diagonal (v_k separately).
        let mut v: Vec<f64> = (k..m).map(|i| a[i][k]).collect();
        v[0] -= alpha;
        let vn: f64 = v.iter().map(|t| t * t).sum();
        if vn > 0.0 {
            for j in k..5 {
                let s: f64 = (k..m).map(|i| v[i - k] * a[i][j]).sum::<f64>() * 2.0 / vn;
                for i in k..m {
                    a[i][j] -= s * v[i - k];
                }
            }
            let s: f64 = (k..m).map(|i| v[i - k] * y[i]).sum::<f64>() * 2.0 / vn;
            for i in k..m {
                y[i] -= s * v[i - k];
            }
        }
    }
    // Rank check: R is upper triangular now.
    let big = (0..5).map(|k| a[k][k].abs()).fold(0.0, f64::max);
    if !(big > 0.0) || (0..5).any(|k| a[k][k].abs() < 1e-10 * big) {
        return None;
    }
    let mut x = [0.0; 5];
    for k in (0..5).rev() {
        let mut s = y[k];
        for j in k + 1..5 {
            s -= a[k][j] * x[j];
        }
        x[k] = s / a[k][k];
    }
    x.iter().all(|t| t.is_finite()).then_some(x)
}

/// Unit area-weighted vertex normals, `None` where the sum vanishes.
fn pseudo_normals(mesh: &TriMesh, p: &[V3]) -> Vec<Option<V3>> {
    let mut acc = vec![[0.0; 3]; mesh.n_vertices()];
    for &[a, b, c] in mesh.faces() {
        let pa = p[a as usize];
        let n = cross(sub(p[b as usize], pa), sub(p[c as usize], pa));
        for v in [a, b, c] {
            acc[v as usize] = add(acc[v as usize], n);
        }
    }
    acc.into_iter().map(unit).collect()
}

fn mean_edge_length(mesh: &TriMesh, p: &[V3]) -> f64 {
    if mesh.n_edges() == 0 {
        return 0.0;
    }
    let s: f64 = mesh
        .edges()
        .iter()
        .map(|&[a, b]| norm(sub(p[a as usize], p[b as usize])))
        .sum();
    s / mesh.n_edges() as f64
}

/// Principal curvatures and directions at every vertex, at scale
/// `options.scale ×` the mean edge length, by the method of Panozzo, Puppo and Rocca (see the
/// module documentation).
///
/// # Errors
/// [`CurvatureError`]: wrong number of positions, non-finite positions, a scale that is not
/// positive and finite, `max_samples < 5`. Degenerate geometry is not an error (see the
/// per-vertex [`CurvatureStatus`]).
pub fn principal_curvatures(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    options: &CurvatureOptions,
) -> Result<PrincipalCurvatures, CurvatureError> {
    let n = mesh.n_vertices();
    if positions.len() != n {
        return Err(CurvatureError::PositionCount {
            expected: n,
            got: positions.len(),
        });
    }
    if let Some(vertex) = positions.iter().position(|q| q.iter().any(|t| !t.is_finite())) {
        return Err(CurvatureError::NonFinitePosition { vertex });
    }
    if !(options.scale > 0.0 && options.scale.is_finite()) {
        return Err(CurvatureError::BadOption {
            reason: format!("scale must be positive and finite, got {}", options.scale),
        });
    }
    if let Some(m) = options.max_samples
        && m < 5
    {
        return Err(CurvatureError::BadOption {
            reason: format!("max_samples must be at least 5 (the fit has 5 unknowns), got {m}"),
        });
    }
    let e_mean = mean_edge_length(mesh, positions);
    let radius = options.scale * e_mean;
    let normals = pseudo_normals(mesh, positions);

    let mut out = PrincipalCurvatures {
        k_max: vec![0.0; n],
        k_min: vec![0.0; n],
        d_max: vec![[1.0, 0.0, 0.0]; n],
        d_min: vec![[0.0, 1.0, 0.0]; n],
        normal: vec![[0.0, 0.0, 1.0]; n],
        status: vec![CurvatureStatus::Fitted; n],
        n_fitted: vec![0; n],
        is_boundary: (0..n).map(|v| mesh.is_boundary_vertex(v)).collect(),
        radius,
        mean_edge_length: e_mean,
    };

    // Reused buffers.
    let mut stamp = vec![u32::MAX; n];
    let mut queue: Vec<u32> = Vec::new();
    for v in 0..n {
        let pv = positions[v];
        let Some(nv) = normals[v] else {
            out.status[v] = CurvatureStatus::DegenerateNormal;
            continue;
        };
        out.normal[v] = nv;
        let tangent = |out: &mut PrincipalCurvatures, n: V3| {
            let t = any_orthogonal(n);
            out.d_max[v] = t;
            out.d_min[v] = cross(n, t);
        };
        tangent(&mut out, nv);

        // 1. Breadth-first traversal inside the ball.
        queue.clear();
        queue.push(v as u32);
        stamp[v] = v as u32;
        let mut head = 0;
        while head < queue.len() {
            let x = queue[head] as usize;
            head += 1;
            for y in mesh.vertex_neighbors(x) {
                if stamp[y] == v as u32 {
                    continue;
                }
                if norm(sub(positions[y], pv)) <= radius {
                    stamp[y] = v as u32;
                    queue.push(y as u32);
                }
            }
        }
        // 2. Fitting frame and the discarding of vertices that bend away.
        let mut wsum = [0.0; 3];
        let mut count = 0usize;
        for &q in &queue {
            if let Some(nq) = normals[q as usize] {
                wsum = add(wsum, nq);
                count += 1;
            }
        }
        let Some(w) = (count > 0).then(|| unit(wsum)).flatten() else {
            out.status[v] = CurvatureStatus::DegenerateNormal;
            continue;
        };
        let mut set: Vec<u32> = queue
            .iter()
            .copied()
            .filter(|&q| {
                q as usize != v && normals[q as usize].is_some_and(|nq| dot(nq, w) >= 0.0)
            })
            .collect();
        // 3. Monte-Carlo.
        if let Some(m) = options.max_samples
            && set.len() > m
        {
            let mut rng = SplitMix64::new(options.seed ^ (v as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            for i in 0..m {
                let j = i + (rng.next_u64() % (set.len() - i) as u64) as usize;
                set.swap(i, j);
            }
            set.truncate(m);
            set.sort_unstable();
        }
        out.n_fitted[v] = set.len() as u32 + 1;
        out.normal[v] = w;
        tangent(&mut out, w);
        let eu = any_orthogonal(w);
        let ev = cross(w, eu);
        if set.len() < 5 || !(radius > 0.0) {
            out.status[v] = CurvatureStatus::Underdetermined;
            continue;
        }
        // 4. Fit in coordinates scaled by r.
        let mut a = Vec::with_capacity(set.len());
        let mut y = Vec::with_capacity(set.len());
        for &q in &set {
            let d = sub(positions[q as usize], pv);
            let (u, vv, ww) = (dot(d, eu) / radius, dot(d, ev) / radius, dot(d, w) / radius);
            a.push([u * u, vv * vv, u * vv, u, vv]);
            y.push(ww);
        }
        let Some(c) = lstsq5(&mut a, &mut y) else {
            out.status[v] = CurvatureStatus::Underdetermined;
            continue;
        };
        // Back to unscaled coefficients: second order divided by r, first order unchanged.
        let (ca, cb, cc) = (c[0] / radius, c[1] / radius, c[2] / radius);
        let (cd, ce) = (c[3], c[4]);
        // 5. Shape operator.
        let s = (1.0 + cd * cd + ce * ce).sqrt();
        let (e1, f1, g1) = (1.0 + cd * cd, cd * ce, 1.0 + ce * ce);
        let (l2, m2, n2) = (2.0 * ca / s, cc / s, 2.0 * cb / s);
        let det = e1 * g1 - f1 * f1;
        let s11 = (g1 * l2 - f1 * m2) / det;
        let s12 = (g1 * m2 - f1 * n2) / det;
        let s21 = (e1 * m2 - f1 * l2) / det;
        let s22 = (e1 * n2 - f1 * m2) / det;
        let tr = s11 + s22;
        let dt = s11 * s22 - s12 * s21;
        let disc = (tr * tr / 4.0 - dt).max(0.0).sqrt();
        let (mut k1, mut k2) = (tr / 2.0 + disc, tr / 2.0 - disc);
        // Eigenvector of S for eigenvalue k (coordinates in the basis X_u, X_v).
        let eigvec = |k: f64| -> Option<(f64, f64)> {
            let r1 = (-s12, s11 - k);
            let r2 = (s22 - k, -s21);
            let (n1, n2) = (r1.0.hypot(r1.1), r2.0.hypot(r2.1));
            let r = if n1 >= n2 { r1 } else { r2 };
            (r.0.hypot(r.1) > 1e-12 * (1.0 + tr.abs())).then_some(r)
        };
        let to3 = |(xi, eta): (f64, f64)| -> Option<V3> {
            let xu = add(eu, scale(w, cd));
            let xv = add(ev, scale(w, ce));
            unit(add(scale(xu, xi), scale(xv, eta)))
        };
        let nfit = unit(sub(w, add(scale(eu, cd), scale(ev, ce)))).unwrap_or(w);
        // k1 >= k2 in the paper's convention (S = -dN, outward N); flip for convex positive.
        let mut dir1 = eigvec(k1).and_then(to3);
        if options.convex_positive {
            // Negating the eigenvalues swaps the roles.
            let (a1, a2) = (-k2, -k1);
            k1 = a1;
            k2 = a2;
            dir1 = eigvec(-k1).and_then(to3);
        }
        if let Some(t1) = dir1 {
            // Orthogonalize against the fitted normal, then d_min = n x d_max.
            let t1 = unit(sub(t1, scale(nfit, dot(t1, nfit)))).unwrap_or(t1);
            out.d_max[v] = t1;
            out.d_min[v] = cross(nfit, t1);
        } else {
            // Umbilic: any orthonormal pair of the tangent plane of the graph.
            let t1 = unit(sub(eu, scale(nfit, dot(eu, nfit)))).unwrap_or(eu);
            out.d_max[v] = t1;
            out.d_min[v] = cross(nfit, t1);
        }
        out.normal[v] = nfit;
        out.k_max[v] = k1;
        out.k_min[v] = k2;
    }
    Ok(out)
}

// ------------------------------------------------------------------------------------------
// Directions to angles
// ------------------------------------------------------------------------------------------

/// Converts 3D directions at the vertices to angles in the normalized intrinsic frames of the
/// input mesh: `directions[v]` is a direction at vertex `v` (any non-zero vector; it is
/// projected, so a tangent vector from [`principal_curvatures`] or any 3D vector works). The
/// output is the direction as a signpost-style angle: `[0, 2π)` at an interior vertex,
/// `[0, π]` at a boundary vertex, measured counter-clockwise from the direction of
/// `TriMesh::outgoing_halfedge(v)` and scaled by `2π/Θ` (interior) or `π/Θ` (boundary), the
/// frame of `IntrinsicTriangulation::signpost` (paper Sec. 3.1; SPEC S4).
///
/// **The rule** ([`direction_angle`]): the corners of the faces around `v` are wedges between
/// consecutive outgoing edges. The direction is projected onto the plane of each incident face;
/// a face *contains* it if the projection lies in the wedge (non-negative coefficients on the two
/// edge vectors). Among the containing faces the one with the largest projected length (the face
/// the direction is most parallel to) wins, ties to the lowest halfedge id. The angle is the
/// wedge's start signpost plus the fraction `φ/θ` of its width, with `φ` the planar angle in
/// the face from the first edge and `θ` the corner angle. If no face contains the direction
/// (it points out of a boundary, or into a gap at a saddle-like vertex where the wedges overlap
/// in projection) the nearest wedge edge in planar angle is taken, so the result is always
/// defined. A direction is a line for curvature: add `π` (Euclidean; `π·2π/Θ` normalized at an
/// interior vertex) for the opposite direction.
///
/// # Errors
/// [`CurvatureError`]: wrong number of positions or directions, non-finite positions, a zero or
/// non-finite direction, a mesh whose faces have no area (no intrinsic triangulation).
pub fn direction_angles(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    directions: &[[f64; 3]],
) -> Result<Vec<f64>, CurvatureError> {
    let n = mesh.n_vertices();
    if directions.len() != n {
        return Err(CurvatureError::BadInput {
            reason: format!("{} directions for {n} vertices", directions.len()),
        });
    }
    let geometry = Geometry::from_positions(mesh, positions)?;
    let tri = IntrinsicTriangulation::new(mesh, &geometry)?;
    (0..n)
        .map(|v| direction_angle(mesh, positions, &tri, v, directions[v]))
        .collect()
}

/// The angle of one direction at vertex `v`; see [`direction_angles`] for the rule. `tri` must
/// be `IntrinsicTriangulation::new(mesh, …)` of the same mesh (its halfedge ids are the mesh's).
///
/// # Errors
/// As [`direction_angles`], and a vertex out of range.
pub fn direction_angle(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    tri: &IntrinsicTriangulation,
    v: usize,
    d: [f64; 3],
) -> Result<f64, CurvatureError> {
    if v >= mesh.n_vertices() || positions.len() != mesh.n_vertices() {
        return Err(CurvatureError::BadInput {
            reason: format!("vertex {v} out of range or positions of the wrong length"),
        });
    }
    if tri.n_halfedges() != mesh.n_halfedges() {
        return Err(CurvatureError::BadInput {
            reason: "the intrinsic triangulation does not belong to the mesh".into(),
        });
    }
    let dn = norm(d);
    if !(dn.is_finite() && dn > 0.0) {
        return Err(CurvatureError::BadDirection { vertex: v });
    }
    let d = scale(d, 1.0 / dn);
    let pv = positions[v];
    let boundary = mesh.is_boundary_vertex(v);

    struct Wedge {
        start: f64,
        width: f64,
        theta: f64,
        ea: V3,
        eb: V3,
        nf: V3,
    }
    let mut wedges = Vec::new();
    let mut total = 0.0;
    for h in mesh.outgoing_halfedges(v) {
        let ea = sub(positions[mesh.to_vertex(h)], pv);
        let eb = sub(positions[mesh.from_vertex(mesh.prev(h))], pv);
        let cr = cross(ea, eb);
        let theta = norm(cr).atan2(dot(ea, eb));
        let nf = unit(cr).unwrap_or([0.0, 0.0, 1.0]);
        total += theta;
        wedges.push(Wedge {
            start: tri.signpost(h),
            width: 0.0,
            theta,
            ea,
            eb,
            nf,
        });
    }
    if wedges.is_empty() || !(total > 0.0) {
        return Err(CurvatureError::BadInput {
            reason: format!("vertex {v} has no non-degenerate face"),
        });
    }
    let s = (if boundary { PI } else { TAU }) / total;
    for w in &mut wedges {
        w.width = w.theta * s;
    }
    let wrap = |a: f64| if boundary { a } else { a.rem_euclid(TAU) };

    let mut best: Option<(f64, f64)> = None; // (score, angle)
    for w in &wedges {
        let p = sub(d, scale(w.nf, dot(d, w.nf)));
        let pl = norm(p);
        if pl < 1e-12 || w.theta <= 0.0 {
            continue;
        }
        // p = alpha ea + beta eb
        let (g11, g12, g22) = (dot(w.ea, w.ea), dot(w.ea, w.eb), dot(w.eb, w.eb));
        let (r1, r2) = (dot(p, w.ea), dot(p, w.eb));
        let det = g11 * g22 - g12 * g12;
        if det <= 0.0 {
            continue;
        }
        let alpha = (r1 * g22 - r2 * g12) / det;
        let beta = (r2 * g11 - r1 * g12) / det;
        if alpha >= -1e-12 && beta >= -1e-12 {
            let phi = norm(cross(w.ea, p)).atan2(dot(w.ea, p));
            let angle = w.start + (phi / w.theta).min(1.0) * w.width;
            if best.is_none_or(|(sc, _)| pl > sc + 1e-12) {
                best = Some((pl, angle));
            }
        }
    }
    if let Some((_, a)) = best {
        return Ok(wrap(a));
    }
    // Fallback: the nearest wedge edge in planar angle.
    let mut near: Option<(f64, f64)> = None; // (planar angle distance, angle)
    for w in &wedges {
        let p = sub(d, scale(w.nf, dot(d, w.nf)));
        if norm(p) < 1e-12 {
            continue;
        }
        let da = norm(cross(w.ea, p)).atan2(dot(w.ea, p));
        let db = norm(cross(w.eb, p)).atan2(dot(w.eb, p));
        for (dist, angle) in [(da, w.start), (db, w.start + w.width)] {
            if near.is_none_or(|(b, _)| dist < b) {
                near = Some((dist, angle));
            }
        }
    }
    // A direction orthogonal to every face plane (cannot happen on a valid mesh): take 0.
    Ok(wrap(near.map_or(0.0, |(_, a)| a)))
}
