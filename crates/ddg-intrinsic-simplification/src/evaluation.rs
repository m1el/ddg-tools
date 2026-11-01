//! Evaluation of a bijective map between a fine and a coarse mesh (paper Sec. 8.3 and 8.4;
//! SPEC §10.3, §10.4): the parametric distortion of the map, and the extrinsic baseline wrapped
//! as a [`CoarseApplications`] so that every comparison of Sec. 8.5 runs on it unchanged.
//!
//! # Distortion (SPEC §10.3, `[formula: iem-distortion]`)
//!
//! For a map `φ` from the fine mesh `M` to the coarse mesh `M̄` the paper approximates the
//! length of the image of a fine edge `ij` by the distance `ℓ̂_ij` between the image points
//! `φ(i)`, `φ(j)` *on `M̄`* ("since all fine edges are minimal geodesics ... we assume this
//! property is preserved"): in the layout of the coarse face if both are in one face, else the
//! geodesic distance. With `ℓ̂` on the three edges of a fine face, the triangle of those lengths
//! is compared to the input triangle: the affine map `J` from the input triangle to it has the
//! metric tensor `JᵀJ` with eigenvalues `γ_max ≥ γ_min` (`γ = σ²`), and, as Khodakovsky,
//! Litke, Schröder 2003 (Sec. 2, "Distortion Measures") define them:
//!
//! * anisotropic distortion `Γ_a = γ_max / γ_min` (scale invariant, `1` for a similarity),
//! * area distortion `Γ_d = γ_max γ_min = det g`; we report `log Γ_d / 2 = log(σ₁ σ₂)` (`0` for
//!   an area preserving map; the log is symmetric between shrinking and growing),
//! * and `σ₁/σ₂ = √Γ_a` for the convenience of the reader.
//!
//! Mean (weighted by the area of the input face) and maximum over the fine faces.
//!
//! ## The distance between two image points (our choice)
//!
//! The paper uses the exact geodesic algorithm of Mitchell et al. for every pair that is not in
//! one coarse face. One exact run per fine vertex is far too slow, and the exact runs of
//! `ddg-intrinsic-apps` give distances to *vertices*, not between points inside faces. We use a
//! **local unfolding search** instead, which is the exact geodesic distance whenever the shortest
//! path crosses at most [`DistortionOptions::max_depth`] faces and bends at no more than one
//! vertex, and an upper bound otherwise:
//!
//! 1. From the source image `s` the faces are unfolded across their edges into the plane,
//!    carrying the visibility cone from `s` (the window), exactly as the exact algorithm does
//!    for one source without merging. If the target image `t` lies in an unfolded face inside
//!    the cone, `|s − t|` in the plane is the length of a straight path on `M̄`; the minimum
//!    over all such unfoldings is the shortest path that touches no vertex. Branches whose edge
//!    is farther from `s` than the best length are pruned, so a pair in one face with both
//!    points far from the sides costs three distance tests. (Pairs in one face are searched
//!    too: on a sharp cone a path around the apex can be shorter than the straight segment.)
//! 2. Every coarse vertex seen inside a cone on the way gets its straight distance `d_s(c)` from
//!    `s`; the same search from `t` gives `d_t(c)`. `min_c d_s(c) + d_t(c)` is the shortest path
//!    with one bend at a vertex (a shortest path bends only at vertices of angle sum above
//!    `2π`).
//! 3. If both fail (the caps were hit), the vertex route `|s − a| + d(a, b) + |b − t|` through
//!    corners `a`, `b` of the two faces with the exact vertex-to-vertex distances of
//!    [`CoarseApplications::exact_distance`]. These pairs are counted as `n_fallback`.
//!
//! Sources on a side or a corner are moved by `1e-9` of the face towards its centre, so the
//! cones are never degenerate (a relative error `≤ 1e-9`). The test of the module checks the
//! search against the exact geodesic distance on a developable cone and a flat torus.
//!
//! Fine faces whose `ℓ̂` violate the triangle inequality (the image of the edges is not a
//! triangle) are counted in [`DistortionReport::n_invalid_faces`] with their area share and left
//! out of the means; none of the numbers hides them.
//!
//! # The extrinsic maps as applications
//!
//! [`extrinsic_applications`] turns a [`QemResult`] and one of its two maps (the closest-point
//! projection or the successive self-parameterization of Liu et al. 2021) into a
//! [`CoarseApplications`]: the coarse mesh with the lengths of its *extrinsic* edges, the fine
//! vertices tracked at the images. Poisson, geodesics, all pairs and Voronoi of Sec. 8.5 then run
//! unchanged, as the paper does for QEM in Fig. 21 and 22. The angular coordinates `phi` of
//! the output are those of the extrinsic lengths with the reference halfedge of each vertex
//! (they are not used by the scalar applications).

use std::collections::HashMap;
use std::fmt;

use ddg_core::sparse::ThreadPool;
use ddg_intrinsic_apps::exact::{ExactOptions, ExactSource};

use crate::applications::{AppError, CoarseApplications};
use crate::baseline::qem::{ProjectedPoint, QemResult};
use crate::mesh::{CoarseMesh, CoarseMeshOptions};
use crate::output::{CoarseOutput, CoarsePoint};

/// Why an evaluation failed.
#[derive(Debug, Clone, PartialEq)]
pub enum EvalError {
    /// The input is inconsistent (the message says what and where).
    Invalid(String),
    /// An application of the coarse mesh failed.
    App(String),
}

impl fmt::Display for EvalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EvalError::Invalid(m) => write!(f, "invalid evaluation input: {m}"),
            EvalError::App(m) => write!(f, "coarse mesh application failed: {m}"),
        }
    }
}

impl std::error::Error for EvalError {}

impl From<AppError> for EvalError {
    fn from(e: AppError) -> Self {
        EvalError::App(e.to_string())
    }
}

type P2 = [f64; 2];

fn sub2(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn cross2(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn len2(a: P2) -> f64 {
    (a[0] * a[0] + a[1] * a[1]).sqrt()
}

/// The planar layout of a face from its side lengths (`l[k]` from slot `k` to `k + 1`): slot 0
/// at the origin, slot 1 on the x axis, slot 2 above it. `None` if the lengths are not a
/// triangle.
fn layout(l: [f64; 3]) -> Option<[P2; 3]> {
    if !(l[0] > 0.0 && l[0].is_finite()) {
        return None;
    }
    let x = (l[0] * l[0] + l[2] * l[2] - l[1] * l[1]) / (2.0 * l[0]);
    let y2 = l[2] * l[2] - x * x;
    if !(y2 > 0.0) || !y2.is_finite() {
        return None;
    }
    Some([[0.0, 0.0], [l[0], 0.0], [x, y2.sqrt()]])
}

/// Options of [`distortion`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistortionOptions {
    /// The most faces an unfolding may cross (default 16).
    pub max_depth: usize,
    /// The most unfolded faces visited per search (default 4000).
    pub max_nodes: usize,
    /// The nudge of points on a side or corner into the face (default `1e-9`).
    pub nudge: f64,
}

impl Default for DistortionOptions {
    fn default() -> Self {
        DistortionOptions {
            max_depth: 16,
            max_nodes: 4000,
            nudge: 1e-9,
        }
    }
}

/// The distortion of a map (SPEC §10.3).
#[derive(Debug, Clone, PartialEq)]
pub struct DistortionReport {
    /// Fine faces, and those whose `ℓ̂` are not a triangle (left out of the statistics).
    pub n_faces: usize,
    /// Faces whose `ℓ̂` violate the (strict) triangle inequality or are not finite.
    pub n_invalid_faces: usize,
    /// The share of the input area they cover.
    pub invalid_area_fraction: f64,
    /// Area-weighted mean of the anisotropic distortion `γ_max/γ_min`.
    pub aniso_mean: f64,
    /// Maximum of `γ_max/γ_min`.
    pub aniso_max: f64,
    /// Area-weighted mean of `σ₁/σ₂`.
    pub sigma_ratio_mean: f64,
    /// Maximum of `σ₁/σ₂`.
    pub sigma_ratio_max: f64,
    /// Area-weighted mean of the absolute area distortion `|log(σ₁σ₂)|`.
    pub area_mean_abs: f64,
    /// Maximum of `|log(σ₁σ₂)|`.
    pub area_max_abs: f64,
    /// Area-weighted mean of the signed `log(σ₁σ₂)` (negative: the image is smaller).
    pub area_mean_signed: f64,
    /// Per fine face `[γ_max/γ_min, log(σ₁σ₂), σ₁/σ₂]`, NaN for invalid faces.
    pub per_face: Vec<[f64; 3]>,
    /// Unique fine edges.
    pub n_edges: usize,
    /// Edges whose two images lie in one coarse face and the straight segment was the answer.
    pub n_same_face: usize,
    /// Edges answered by an unfolding that touches no vertex.
    pub n_unfolded: usize,
    /// Edges answered by a path with one bend at a vertex.
    pub n_via_vertex: usize,
    /// Edges answered by the vertex-route fallback (an upper bound).
    pub n_fallback: usize,
    /// Edges for which nothing was found (their faces are invalid).
    pub n_unresolved: usize,
}

/// The distance between two points of a coarse mesh by the local unfolding search of the
/// [module documentation](self).
pub struct PointDistance<'a> {
    out: &'a CoarseOutput,
    lay: Vec<[P2; 3]>,
    opts: DistortionOptions,
}

impl<'a> PointDistance<'a> {
    /// Prepares the layouts of the faces of `out`.
    ///
    /// # Errors
    /// [`EvalError::Invalid`] naming a face whose lengths are not a triangle.
    pub fn new(
        out: &'a CoarseOutput,
        options: &DistortionOptions,
    ) -> Result<PointDistance<'a>, EvalError> {
        let lay = out
            .face_lengths
            .iter()
            .enumerate()
            .map(|(f, l)| {
                layout(*l).ok_or_else(|| {
                    EvalError::Invalid(format!("coarse face {f} has side lengths {l:?}"))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(PointDistance {
            out,
            lay,
            opts: *options,
        })
    }

    /// The distance between the points `(face, barycentric coordinates)` `a` and `b`; `None` if
    /// the search found no path (the caps were hit).
    pub fn distance(&self, a: (u32, [f64; 3]), b: (u32, [f64; 3])) -> Option<f64> {
        self.distance_kind(a, b).map(|(d, _)| d)
    }
}

#[derive(Clone, Copy)]
struct Frame {
    face: u32,
    q: [P2; 3],
    entry: i8,
    cone: Option<(P2, P2)>,
    depth: u32,
}

fn in_cone(cone: &Option<(P2, P2)>, s: P2, x: P2) -> bool {
    let Some((dr, dl)) = cone else { return true };
    let d = sub2(x, s);
    let tol = 1e-12;
    cross2(*dr, d) >= -tol * len2(*dr) * len2(d) && cross2(d, *dl) >= -tol * len2(*dl) * len2(d)
}

fn clip(old: &Option<(P2, P2)>, ra: P2, rb: P2) -> Option<(P2, P2)> {
    if cross2(ra, rb) <= 0.0 {
        return None;
    }
    let (dr, dl) = match old {
        None => (ra, rb),
        Some((dr0, dl0)) => (
            if cross2(*dr0, ra) >= 0.0 { ra } else { *dr0 },
            if cross2(rb, *dl0) >= 0.0 { rb } else { *dl0 },
        ),
    };
    (cross2(dr, dl) > 1e-14 * len2(dr) * len2(dl)).then_some((dr, dl))
}

fn seg_dist(s: P2, a: P2, b: P2) -> f64 {
    let ab = sub2(b, a);
    let l2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if l2 > 0.0 {
        ((sub2(s, a)[0] * ab[0] + sub2(s, a)[1] * ab[1]) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    len2(sub2(s, [a[0] + t * ab[0], a[1] + t * ab[1]]))
}

type FacePoint = (u32, [f64; 3]);

impl PointDistance<'_> {
    fn pos(&self, p: FacePoint) -> P2 {
        let l = &self.lay[p.0 as usize];
        [
            p.1[0] * l[0][0] + p.1[1] * l[1][0] + p.1[2] * l[2][0],
            p.1[0] * l[0][1] + p.1[1] * l[1][1] + p.1[2] * l[2][1],
        ]
    }

    fn nudge(&self, p: FacePoint) -> FacePoint {
        let e = self.opts.nudge;
        (p.0, p.1.map(|b| (1.0 - e) * b + e / 3.0))
    }

    /// The unfolding search from `src`. Updates `best` with the straight distances to `tgt`
    /// (if any) and appends `(vertex, distance)` of every vertex seen inside a cone. Returns
    /// whether a cap was hit.
    fn search(
        &self,
        src: FacePoint,
        tgt: Option<FacePoint>,
        best: &mut f64,
        seen: &mut Vec<(u32, f64)>,
    ) -> bool {
        let out = self.out;
        let s = self.pos(src);
        let t = tgt.map(|t| (t.0, self.pos(t)));
        let mut stack = vec![Frame {
            face: src.0,
            q: self.lay[src.0 as usize],
            entry: -1,
            cone: None,
            depth: 0,
        }];
        let mut nodes = 0usize;
        let mut capped = false;
        while let Some(fr) = stack.pop() {
            nodes += 1;
            if nodes > self.opts.max_nodes {
                return true;
            }
            let f = fr.face as usize;
            for c in 0..3 {
                if in_cone(&fr.cone, s, fr.q[c]) {
                    seen.push((out.faces[f][c], len2(sub2(fr.q[c], s))));
                }
            }
            if let Some((tf, tp)) = t
                && tf == fr.face
            {
                // The target's planar position is in the layout of its own face; in an unfolded
                // copy it is the same affine combination of the unfolded corners.
                let b = tgt.expect("target").1;
                let x = [
                    b[0] * fr.q[0][0] + b[1] * fr.q[1][0] + b[2] * fr.q[2][0],
                    b[0] * fr.q[0][1] + b[1] * fr.q[1][1] + b[2] * fr.q[2][1],
                ];
                let _ = tp;
                if in_cone(&fr.cone, s, x) {
                    *best = best.min(len2(sub2(x, s)));
                }
            }
            if fr.depth as usize >= self.opts.max_depth {
                if fr.depth as usize == self.opts.max_depth {
                    capped = true;
                }
                continue;
            }
            for k in 0..3 {
                if k as i8 == fr.entry {
                    continue;
                }
                let h = 3 * f + k;
                let Some(tw) = out.twin[h] else { continue };
                let (a, b) = (fr.q[k], fr.q[(k + 1) % 3]);
                let cone = match clip(&fr.cone, sub2(a, s), sub2(b, s)) {
                    Some(c) => Some(c),
                    None => continue,
                };
                if seg_dist(s, a, b) >= *best {
                    continue;
                }
                let g = (tw / 3) as usize;
                let kk = (tw % 3) as usize;
                let lg = out.face_lengths[g];
                // g's slot kk is at b, slot kk+1 at a; the third corner is on the far side.
                let d = len2(sub2(a, b));
                if !(d > 0.0) {
                    continue;
                }
                let u = [(a[0] - b[0]) / d, (a[1] - b[1]) / d];
                let n = [-u[1], u[0]];
                let l_bc = lg[(kk + 2) % 3]; // side from slot kk+2 to slot kk
                let l_ac = lg[(kk + 1) % 3]; // side from slot kk+1 to slot kk+2
                let x = (l_bc * l_bc + d * d - l_ac * l_ac) / (2.0 * d);
                let y2 = l_bc * l_bc - x * x;
                if !y2.is_finite() {
                    continue;
                }
                let y = y2.max(0.0).sqrt();
                let c = [b[0] + x * u[0] + y * n[0], b[1] + x * u[1] + y * n[1]];
                let mut q = [[0.0; 2]; 3];
                q[kk] = b;
                q[(kk + 1) % 3] = a;
                q[(kk + 2) % 3] = c;
                stack.push(Frame {
                    face: g as u32,
                    q,
                    entry: kk as i8,
                    cone,
                    depth: fr.depth + 1,
                });
            }
        }
        capped
    }

    /// The distance between two points of the coarse mesh, and how it was found.
    fn distance_kind(&self, s: FacePoint, t: FacePoint) -> Option<(f64, Kind)> {
        let (s, t) = (self.nudge(s), self.nudge(t));
        let mut best = f64::INFINITY;
        let same = s.0 == t.0;
        if same {
            best = len2(sub2(self.pos(s), self.pos(t)));
        }
        let planar = best;
        let mut seen_s = Vec::new();
        self.search(s, Some(t), &mut best, &mut seen_s);
        let direct = best;
        let mut seen_t = Vec::new();
        let mut bound = best;
        self.search(t, None, &mut bound, &mut seen_t);
        seen_s.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        seen_s.dedup_by_key(|x| x.0);
        seen_t.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        seen_t.dedup_by_key(|x| x.0);
        let mut via = f64::INFINITY;
        let (mut i, mut j) = (0, 0);
        while i < seen_s.len() && j < seen_t.len() {
            match seen_s[i].0.cmp(&seen_t[j].0) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    via = via.min(seen_s[i].1 + seen_t[j].1);
                    i += 1;
                    j += 1;
                }
            }
        }
        let total = direct.min(via);
        if !total.is_finite() {
            return None;
        }
        let kind = if via < direct {
            Kind::Via
        } else if same && direct >= planar - 1e-15 * planar.max(1.0) && direct == planar {
            Kind::Same
        } else {
            Kind::Unfolded
        };
        Some((total, kind))
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Same,
    Unfolded,
    Via,
}

/// The coarse location of every fine vertex as a face and barycentric coordinates.
fn image_points(out: &CoarseOutput) -> Result<Vec<FacePoint>, EvalError> {
    let first = out.face_of_points();
    (0..out.n_fine_vertices)
        .map(|i| match out.points.get(i) {
            Some(&CoarsePoint::Face { face, bary }) => Ok((face, bary)),
            Some(&CoarsePoint::Vertex(v)) => {
                let f = first[i].ok_or_else(|| {
                    EvalError::Invalid(format!("coarse vertex {v} (image of {i}) has no face"))
                })?;
                let k = out.faces[f as usize]
                    .iter()
                    .position(|&x| x == v)
                    .expect("vertex of its face");
                let mut b = [0.0; 3];
                b[k] = 1.0;
                Ok((f, b))
            }
            None => Err(EvalError::Invalid(format!(
                "fine vertex {i} is not tracked"
            ))),
        })
        .collect()
}

/// The singular values of the affine map from the planar triangle `a` to `b`, `σ₁ ≥ σ₂`.
fn singular_values(a: [P2; 3], b: [P2; 3]) -> Option<(f64, f64)> {
    let (a1, a2) = (sub2(a[1], a[0]), sub2(a[2], a[0]));
    let (b1, b2) = (sub2(b[1], b[0]), sub2(b[2], b[0]));
    let det = cross2(a1, a2);
    if !(det.abs() > 0.0) {
        return None;
    }
    // J = [b1 b2] [a1 a2]^{-1}
    let inv = [[a2[1] / det, -a2[0] / det], [-a1[1] / det, a1[0] / det]];
    let j = [
        [
            b1[0] * inv[0][0] + b2[0] * inv[1][0],
            b1[0] * inv[0][1] + b2[0] * inv[1][1],
        ],
        [
            b1[1] * inv[0][0] + b2[1] * inv[1][0],
            b1[1] * inv[0][1] + b2[1] * inv[1][1],
        ],
    ];
    let (p, q) = (j[0][0] + j[1][1], j[0][0] - j[1][1]);
    let (r, s) = (j[1][0] + j[0][1], j[1][0] - j[0][1]);
    let big_r = (p * p + s * s).sqrt();
    let big_s = (q * q + r * r).sqrt();
    Some((0.5 * (big_r + big_s), 0.5 * (big_r - big_s).abs()))
}

/// The distortion of the map recorded in `app` (the images of the fine vertices
/// `app.output().points[..n_fine]`) from the fine mesh `(positions, faces)` to the coarse mesh
/// (SPEC §10.3; see the [module documentation](self) for the definitions and the distance).
///
/// [formula: iem-distortion]
///
/// # Errors
/// [`EvalError::Invalid`] if the fine mesh does not match the tracked points, a fine face has a
/// non-positive area, or the coarse mesh has a degenerate face.
pub fn distortion(
    positions: &[[f64; 3]],
    faces: &[[u32; 3]],
    app: &CoarseApplications,
    pool: &ThreadPool,
    options: &DistortionOptions,
) -> Result<DistortionReport, EvalError> {
    let out = app.output();
    let n = out.n_fine_vertices;
    if positions.len() != n {
        return Err(EvalError::Invalid(format!(
            "{} fine positions for {n} tracked fine vertices",
            positions.len()
        )));
    }
    if let Some((f, c)) = faces
        .iter()
        .enumerate()
        .find(|(_, c)| c.iter().any(|&v| v as usize >= n))
    {
        return Err(EvalError::Invalid(format!(
            "fine face {f} = {c:?} has a vertex out of range ({n} vertices)"
        )));
    }
    let geo = PointDistance::new(out, options)?;
    let imgs = image_points(out)?;
    let mut keys: Vec<u64> = Vec::with_capacity(3 * faces.len());
    for f in faces {
        for k in 0..3 {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            keys.push((a.min(b) as u64) << 32 | a.max(b) as u64);
        }
    }
    keys.sort_unstable();
    keys.dedup();
    let results: Vec<Option<(f64, Kind)>> = pool.map(keys.len(), |e| {
        let (a, b) = ((keys[e] >> 32) as usize, (keys[e] & 0xffff_ffff) as usize);
        geo.distance_kind(imgs[a], imgs[b])
    });
    let mut report = DistortionReport {
        n_faces: faces.len(),
        n_invalid_faces: 0,
        invalid_area_fraction: 0.0,
        aniso_mean: 0.0,
        aniso_max: 0.0,
        sigma_ratio_mean: 0.0,
        sigma_ratio_max: 0.0,
        area_mean_abs: 0.0,
        area_max_abs: 0.0,
        area_mean_signed: 0.0,
        per_face: vec![[f64::NAN; 3]; faces.len()],
        n_edges: keys.len(),
        n_same_face: 0,
        n_unfolded: 0,
        n_via_vertex: 0,
        n_fallback: 0,
        n_unresolved: 0,
    };
    // Fallback for unresolved edges: the vertex route with exact vertex distances.
    let mut lhat: Vec<f64> = vec![f64::NAN; keys.len()];
    let mut unresolved: Vec<usize> = Vec::new();
    for (e, r) in results.iter().enumerate() {
        match r {
            Some((d, kind)) => {
                lhat[e] = *d;
                match kind {
                    Kind::Same => report.n_same_face += 1,
                    Kind::Unfolded => report.n_unfolded += 1,
                    Kind::Via => report.n_via_vertex += 1,
                }
            }
            None => unresolved.push(e),
        }
    }
    if !unresolved.is_empty() {
        let mut need: Vec<u32> = Vec::new();
        for &e in &unresolved {
            let (a, b) = ((keys[e] >> 32) as usize, (keys[e] & 0xffff_ffff) as usize);
            for p in [imgs[a], imgs[b]] {
                need.extend(out.faces[p.0 as usize]);
            }
        }
        need.sort_unstable();
        need.dedup();
        let rows: Vec<Result<Vec<f64>, AppError>> = pool.map(need.len(), |k| {
            app.exact_distance(
                &[ExactSource::Vertex(need[k] as usize)],
                &ExactOptions::default(),
            )
        });
        let mut row_of: HashMap<u32, Vec<f64>> = HashMap::new();
        for (k, r) in rows.into_iter().enumerate() {
            row_of.insert(need[k], r?);
        }
        for &e in &unresolved {
            let (a, b) = ((keys[e] >> 32) as usize, (keys[e] & 0xffff_ffff) as usize);
            let (pa, pb) = (imgs[a], imgs[b]);
            let (xa, xb) = (geo.pos(pa), geo.pos(pb));
            let mut best = f64::INFINITY;
            for ca in 0..3 {
                for cb in 0..3 {
                    let va = out.faces[pa.0 as usize][ca];
                    let vb = out.faces[pb.0 as usize][cb];
                    let la = geo.lay[pa.0 as usize][ca];
                    let lb = geo.lay[pb.0 as usize][cb];
                    let d = len2(sub2(xa, la)) + row_of[&va][vb as usize] + len2(sub2(xb, lb));
                    best = best.min(d);
                }
            }
            if best.is_finite() {
                lhat[e] = best;
                report.n_fallback += 1;
            } else {
                report.n_unresolved += 1;
            }
        }
    }
    // Per fine face.
    let mut area_total = 0.0;
    let mut valid_area = 0.0;
    let (mut a_sum, mut r_sum, mut d_sum, mut s_sum) = (0.0, 0.0, 0.0, 0.0);
    for (fi, f) in faces.iter().enumerate() {
        let p = f.map(|v| positions[v as usize]);
        let side = |k: usize| {
            let (a, b) = (p[k], p[(k + 1) % 3]);
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        };
        let l = [side(0), side(1), side(2)];
        let Some(la) = layout(l) else {
            return Err(EvalError::Invalid(format!(
                "fine face {fi} = {f:?} has side lengths {l:?}, not a triangle"
            )));
        };
        let area = 0.5 * la[1][0] * la[2][1];
        area_total += area;
        let edge = |k: usize| {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            let key = (a.min(b) as u64) << 32 | a.max(b) as u64;
            lhat[keys.binary_search(&key).expect("edge key")]
        };
        let lh = [edge(0), edge(1), edge(2)];
        let ok = lh.iter().all(|x| x.is_finite())
            && lh[0] < lh[1] + lh[2]
            && lh[1] < lh[0] + lh[2]
            && lh[2] < lh[0] + lh[1];
        let sv = ok
            .then(|| layout(lh))
            .flatten()
            .and_then(|lb| singular_values(la, lb));
        match sv {
            Some((s1, s2)) if s2 > 0.0 && s1.is_finite() => {
                let gamma = (s1 / s2).powi(2);
                let lg = (s1 * s2).ln();
                report.per_face[fi] = [gamma, lg, s1 / s2];
                valid_area += area;
                a_sum += area * gamma;
                r_sum += area * s1 / s2;
                d_sum += area * lg.abs();
                s_sum += area * lg;
                report.aniso_max = report.aniso_max.max(gamma);
                report.sigma_ratio_max = report.sigma_ratio_max.max(s1 / s2);
                report.area_max_abs = report.area_max_abs.max(lg.abs());
            }
            _ => {
                report.n_invalid_faces += 1;
                report.invalid_area_fraction += area;
            }
        }
    }
    if valid_area > 0.0 {
        report.aniso_mean = a_sum / valid_area;
        report.sigma_ratio_mean = r_sum / valid_area;
        report.area_mean_abs = d_sum / valid_area;
        report.area_mean_signed = s_sum / valid_area;
    }
    if area_total > 0.0 {
        report.invalid_area_fraction /= area_total;
    }
    Ok(report)
}

/// Which map of the extrinsic baseline to use (SPEC §10.4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtrinsicMap {
    /// The closest-point projection of the fine vertices onto the coarse mesh (ours).
    ClosestPoint,
    /// The successive self-parameterization of Liu et al. 2021 (needs
    /// `QemOptions::successive_map`).
    Successive,
}

/// The QEM result as a [`CoarseApplications`] on the extrinsic coarse mesh, the fine vertices
/// tracked at their images under `map`. See the [module documentation](self).
///
/// # Errors
/// [`EvalError::Invalid`] if the successive map was not built or the coarse mesh has a
/// degenerate face; [`EvalError::App`].
pub fn extrinsic_applications(
    qem: &QemResult,
    n_fine: usize,
    map: ExtrinsicMap,
) -> Result<CoarseApplications, EvalError> {
    let images: &[ProjectedPoint] = match map {
        ExtrinsicMap::ClosestPoint => &qem.projection,
        ExtrinsicMap::Successive => qem.successive.as_deref().ok_or_else(|| {
            EvalError::Invalid("the QEM run did not build the successive map".to_string())
        })?,
    };
    if images.len() != n_fine {
        return Err(EvalError::Invalid(format!(
            "{} image points for {n_fine} fine vertices",
            images.len()
        )));
    }
    let mesh = CoarseMesh::from_positions(
        &qem.positions,
        &qem.faces,
        &[],
        CoarseMeshOptions::default(),
    )
    .map_err(|e| EvalError::Invalid(e.to_string()))?;
    let mut out = mesh
        .compact()
        .map_err(|e| EvalError::Invalid(e.to_string()))?;
    if out.faces != qem.faces {
        return Err(EvalError::Invalid(
            "the coarse complex renumbered the faces of the QEM mesh".to_string(),
        ));
    }
    out.n_fine_vertices = n_fine;
    out.coarse_of_fine = vec![None; n_fine];
    out.points = images
        .iter()
        .map(|p| {
            let corner = qem.faces[p.face as usize];
            match p.bary.iter().position(|&b| b >= 1.0 - 1e-12) {
                Some(k) => CoarsePoint::Vertex(corner[k]),
                None => CoarsePoint::Face {
                    face: p.face,
                    bary: p.bary,
                },
            }
        })
        .collect();
    Ok(CoarseApplications::from_output(
        out,
        mesh.complex().clone(),
        mesh.edge_lengths().to_vec(),
    )?)
}

/// The total area of the coarse mesh (from its lengths): Fig. 19's quantity.
///
/// # Errors
/// [`EvalError::App`] for a degenerate face.
pub fn coarse_area(app: &CoarseApplications) -> Result<f64, EvalError> {
    Ok(app.output().face_areas()?.iter().sum())
}
