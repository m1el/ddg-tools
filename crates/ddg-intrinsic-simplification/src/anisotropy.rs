//! Anisotropic coarsening (paper Sec. 8.6.2, Fig. 27; SPEC §8).
//!
//! Given a direction `u_i` (a line in the tangent plane) at every vertex, every edge length is
//! scaled by
//!
//! ```text
//!   ℓ'_ij = ℓ_ij · ( (1 − τ) + (τ/2) ((u_i · ê_ij)² + (u_j · ê_ji)²) ),      τ ∈ [0, 1],
//! ```
//!
//! where `ê_ij ∈ T_i M` is the unit tangent vector along edge `ij` at `i`. The factor lies in
//! `[1 − τ, 1]` (`[formula: iem-anisotropic-scale]`): an edge parallel to `u` at both ends keeps
//! its length, an edge orthogonal to `u` at both ends shrinks to `(1 − τ) ℓ`. The coarse mesh is
//! built for the scaled metric ([`anisotropic_coarse_mesh`] = [`CoarseMesh::from_lengths`]
//! with the scaled lengths), so its triangles are about isotropic **in the scaled metric** and
//! therefore **elongated orthogonally to `u`** in the input metric: with `u` the maximum
//! principal direction they are long along the minimum principal direction (along a cylinder's
//! axis), with `u` the minimum direction they are long around it. (SPEC's measurement:
//! `tests/anisotropy.rs`.)
//!
//! # Conventions
//!
//! * **Directions are angles in the normalized frames** of the input mesh, as returned by
//!   [`crate::curvature::direction_angles`]: `[0, 2π)` at an interior vertex, `[0, π]` at a
//!   boundary vertex, measured from the direction of `TriMesh::outgoing_halfedge(v)` and scaled
//!   by `2π/Θ_v` (interior) or `π/Θ_v` (boundary). The angle of the edge in the same frame is
//!   rebuilt from the corner angles of the *input* lengths, summed counter-clockwise from the
//!   reference halfedge (this is `IntrinsicTriangulation::signpost / s`; a test compares them).
//! * The scalar product `u · ê` is taken **in the unnormalized tangent plane**, i.e. with the
//!   real angles: the normalized difference `δ̃` is multiplied by `Θ/2π` (interior) or `Θ/π`
//!   (boundary) to get the angle `δ` between `u` and the edge, `δ ∈ [0, Θ)`. At a cone vertex
//!   (`Θ ≠ 2π`) the two ways around differ; the **smaller angle** `min(δ, Θ − δ)` is used. For a
//!   flat vertex this is exactly `cos²` of the angle between two lines. On the boundary the
//!   angle is `|δ̃| Θ/π` (there is only one way around).
//! * `u` is a **line**: `(u · ê)²` is invariant under `u → −u`. In the normalized frame the
//!   opposite direction is the angle `+ π · 2π/Θ` (interior), as [`crate::curvature::direction_angles`]
//!   documents; at a cone vertex the smaller-angle rule above makes `cos²` of it differ from
//!   `cos²` of `u` only through the cone's two-sided ambiguity.
//! * A boundary edge has an endpoint-direction at each end; the one that is not an outgoing
//!   halfedge of its vertex (the closing edge of the fan) sits at `Θ` (normalized `π`).
//!
//! # Validity and the strength limit
//!
//! A scaled metric can violate the triangle inequality (Campen et al. 2013 Sec. 4.1, noted by
//! the paper), which limits `τ`. The scaled length of edge `e` is affine in `τ`,
//! `ℓ_e (1 − τ (1 − a_e))` with the *alignment* `a_e = (cos²_i + cos²_j)/2 ∈ [0, 1]`, so every
//! inequality `ℓ'_x < ℓ'_y + ℓ'_z` is affine in `τ` and holds on an interval `[0, τ*)`; the
//! valid set of the whole mesh is `[0, τ_max)` and [`max_valid_tau`] computes `τ_max` **exactly**
//! (the minimum of the roots), not by bisection. With
//! [`AnisotropyOptions::clamp_tau`] `false` (default) a `τ` above the limit is an error
//! naming the faces; with `true` the largest valid `τ' = min(τ, τ_max (1 − clamp_safety))` is
//! used and reported ([`AnisotropyReport`]); the clamp is never silent.

#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)]

use std::f64::consts::{PI, TAU};

use ddg_core::geometry::{corner_angle, triangle_area};
use ddg_core::{MeshError, TriMesh};

use crate::curvature::{
    CurvatureError, CurvatureOptions, PrincipalCurvatures, direction_angles, principal_curvatures,
};
use crate::{CoarseMesh, CoarseMeshOptions, IemError, InputPoint};

/// Errors of this module.
#[derive(Debug, Clone, PartialEq)]
pub enum AnisotropyError {
    /// `τ` is not a finite number in `[0, 1]`, or an option is out of range.
    BadOption {
        /// What is wrong.
        reason: String,
    },
    /// Wrong number of lengths or angles, or a non-finite / non-positive value.
    BadInput {
        /// What is wrong.
        reason: String,
    },
    /// The scaled lengths violate the strict triangle inequality at the requested `τ`
    /// (and `clamp_tau` is off).
    TriangleInequality {
        /// The requested `τ`.
        tau: f64,
        /// The largest `τ` for which every face is valid (the faces are valid on `[0, τ_max)`).
        tau_max: f64,
        /// The faces that violate it (all of them, ascending; at most 20 are named in the
        /// message).
        faces: Vec<u32>,
    },
    /// The input mesh is not a manifold oriented triangle mesh.
    Mesh(MeshError),
    /// The curvature estimation or the conversion of directions failed.
    Curvature(CurvatureError),
    /// The coarse mesh could not be built.
    Coarse(IemError),
}

impl std::fmt::Display for AnisotropyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AnisotropyError::BadOption { reason } => write!(f, "bad anisotropy option: {reason}"),
            AnisotropyError::BadInput { reason } => write!(f, "bad anisotropy input: {reason}"),
            AnisotropyError::TriangleInequality {
                tau,
                tau_max,
                faces,
            } => {
                let shown: Vec<String> = faces.iter().take(20).map(|f| f.to_string()).collect();
                write!(
                    f,
                    "tau = {tau} violates the triangle inequality of {} face(s): {}{}; every face \
                     is valid for tau < {tau_max} (use `clamp_tau` to lower it)",
                    faces.len(),
                    shown.join(", "),
                    if faces.len() > 20 { ", ..." } else { "" }
                )
            }
            AnisotropyError::Mesh(e) => write!(f, "{e}"),
            AnisotropyError::Curvature(e) => write!(f, "{e}"),
            AnisotropyError::Coarse(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for AnisotropyError {}

impl From<CurvatureError> for AnisotropyError {
    fn from(e: CurvatureError) -> Self {
        AnisotropyError::Curvature(e)
    }
}
impl From<IemError> for AnisotropyError {
    fn from(e: IemError) -> Self {
        AnisotropyError::Coarse(e)
    }
}
impl From<MeshError> for AnisotropyError {
    fn from(e: MeshError) -> Self {
        AnisotropyError::Mesh(e)
    }
}

/// Options of [`anisotropic_lengths`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnisotropyOptions {
    /// The strength `τ ∈ [0, 1]`. `0` returns the input lengths bitwise.
    pub tau: f64,
    /// If `τ` violates the triangle inequality: `false` (default) is an error, `true` lowers
    /// `τ` to the largest valid value (see `clamp_safety`) and reports it.
    pub clamp_tau: bool,
    /// The clamped `τ` is `τ_max (1 − clamp_safety)`, so that the worst face keeps a relative
    /// margin (its triangle inequality holds with slack `≈ clamp_safety` of its length scale).
    /// In `[0, 1)`. Default `1e-6`.
    pub clamp_safety: f64,
}

impl Default for AnisotropyOptions {
    fn default() -> Self {
        AnisotropyOptions {
            tau: 0.0,
            clamp_tau: false,
            clamp_safety: 1e-6,
        }
    }
}

impl AnisotropyOptions {
    /// Options with strength `tau`, no clamping.
    pub fn with_tau(tau: f64) -> Self {
        AnisotropyOptions {
            tau,
            ..Default::default()
        }
    }
}

/// What [`anisotropic_lengths`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct AnisotropyReport {
    /// The requested `τ`.
    pub tau_requested: f64,
    /// The `τ` used (equal to the requested one unless it was clamped).
    pub tau_used: f64,
    /// Whether `τ` was lowered.
    pub clamped: bool,
    /// The largest `τ'` such that every face is valid on `[0, τ')` (`+∞` if none limits it,
    /// i.e. no root; it can exceed 1).
    pub tau_max: f64,
    /// Number of faces violating the triangle inequality at the requested `τ`.
    pub n_violating_faces_requested: usize,
    /// Smallest factor `ℓ'/ℓ` over the edges at `τ_used`.
    pub min_factor: f64,
    /// Largest factor `ℓ'/ℓ`.
    pub max_factor: f64,
    /// Mean factor.
    pub mean_factor: f64,
}

// ----- the frame of the input mesh -----

/// The direction of every edge at both of its endpoints, as the unnormalized angle from the
/// reference halfedge of the vertex, plus the angle sum.
struct EdgeFrames {
    /// `angle[e][k]`: the unnormalized angle of edge `e` at its endpoint `edge_vertices(e)[k]`.
    angle: Vec<[f64; 2]>,
    /// `Θ_v`.
    angle_sum: Vec<f64>,
    boundary: Vec<bool>,
}

fn edge_frames(mesh: &TriMesh, lengths: &[f64]) -> Result<EdgeFrames, AnisotropyError> {
    let nv = mesh.n_vertices();
    let mut angle = vec![[f64::NAN; 2]; mesh.n_edges()];
    let mut angle_sum = vec![0.0; nv];
    let mut boundary = vec![false; nv];
    let mut put = |e: usize, v: usize, a: f64| {
        let [p, _] = mesh.edge_vertices(e);
        angle[e][usize::from(p != v)] = a;
    };
    for v in 0..nv {
        boundary[v] = mesh.is_boundary_vertex(v);
        let mut theta = 0.0;
        let mut last = None;
        for h in mesh.outgoing_halfedges(v) {
            put(mesh.edge(h), v, theta);
            let (l_h, l_opp, l_prev) = (
                lengths[mesh.edge(h)],
                lengths[mesh.edge(mesh.next(h))],
                lengths[mesh.edge(mesh.prev(h))],
            );
            let c = corner_angle(l_opp, l_h, l_prev).ok_or_else(|| AnisotropyError::BadInput {
                reason: format!(
                    "face {} (lengths {l_h}, {l_opp}, {l_prev}) violates the triangle \
                         inequality",
                    mesh.face(h)
                ),
            })?;
            theta += c;
            last = Some(h);
        }
        if boundary[v] {
            // The closing edge of the fan: the incoming boundary halfedge.
            if let Some(h) = last {
                put(mesh.edge(mesh.prev(h)), v, theta);
            }
        }
        angle_sum[v] = theta;
    }
    Ok(EdgeFrames {
        angle,
        angle_sum,
        boundary,
    })
}

/// `cos²` of the angle between the direction `alpha` (normalized frame angle) and the edge at
/// unnormalized angle `theta_edge`, at a vertex with angle sum `theta_sum` (module docs).
fn cos2_between(alpha: f64, theta_edge: f64, theta_sum: f64, boundary: bool) -> f64 {
    let (s, period) = if boundary { (PI, 0.0) } else { (TAU, TAU) };
    // The unnormalized angle of u: alpha * Theta / s.
    let a_u = alpha * theta_sum / s;
    let mut d = a_u - theta_edge;
    if !boundary {
        d = d.rem_euclid(theta_sum);
        d = d.min(theta_sum - d);
    } else {
        d = d.abs();
    }
    let _ = period;
    let c = d.cos();
    c * c
}

/// The alignment `a_e = (cos²_i + cos²_j)/2 ∈ [0, 1]` of every edge with the directions
/// (module docs); the length factor at strength `τ` is `(1 − τ) + τ a_e`
/// (`[formula: iem-anisotropic-scale]`).
///
/// `lengths` are the input edge lengths (`TriMesh::edges` order), `angles` one normalized frame
/// angle per vertex ([`crate::curvature::direction_angles`], [`principal_direction_field`]).
///
/// # Errors
/// [`AnisotropyError::BadInput`]: counts, non-finite values, a face violating the triangle
/// inequality.
pub fn edge_alignment(
    mesh: &TriMesh,
    lengths: &[f64],
    angles: &[f64],
) -> Result<Vec<f64>, AnisotropyError> {
    if lengths.len() != mesh.n_edges() {
        return Err(AnisotropyError::BadInput {
            reason: format!(
                "{} lengths for a mesh with {} edges",
                lengths.len(),
                mesh.n_edges()
            ),
        });
    }
    if angles.len() != mesh.n_vertices() {
        return Err(AnisotropyError::BadInput {
            reason: format!(
                "{} direction angles for a mesh with {} vertices",
                angles.len(),
                mesh.n_vertices()
            ),
        });
    }
    if let Some(e) = lengths.iter().position(|l| !(l.is_finite() && *l > 0.0)) {
        return Err(AnisotropyError::BadInput {
            reason: format!(
                "edge {e} {:?} has the non-positive or non-finite length {}",
                mesh.edge_vertices(e),
                lengths[e]
            ),
        });
    }
    if let Some(v) = angles.iter().position(|a| !a.is_finite()) {
        return Err(AnisotropyError::BadInput {
            reason: format!("the direction angle of vertex {v} is not finite"),
        });
    }
    let fr = edge_frames(mesh, lengths)?;
    Ok((0..mesh.n_edges())
        .map(|e| {
            let [a, b] = mesh.edge_vertices(e);
            let ci = cos2_between(angles[a], fr.angle[e][0], fr.angle_sum[a], fr.boundary[a]);
            let cj = cos2_between(angles[b], fr.angle[e][1], fr.angle_sum[b], fr.boundary[b]);
            (0.5 * (ci + cj)).clamp(0.0, 1.0)
        })
        .collect())
}

/// The length factor `(1 − τ) + τ a` of an edge of alignment `a`
/// (`[formula: iem-anisotropic-scale]`). It is exactly `1` for `τ = 0`.
#[inline]
pub fn length_factor(tau: f64, alignment: f64) -> f64 {
    (1.0 - tau) + tau * alignment
}

/// The faces that violate the strict triangle inequality of the lengths `ℓ_e (1−τ) + τ a_e ℓ_e`
/// (ascending ids).
pub fn violating_faces(mesh: &TriMesh, lengths: &[f64], alignment: &[f64], tau: f64) -> Vec<u32> {
    (0..mesh.n_faces())
        .filter(|&f| {
            let l = |k: usize| {
                let e = mesh.edge(3 * f + k);
                lengths[e] * length_factor(tau, alignment[e])
            };
            triangle_area(l(0), l(1), l(2)).is_none()
        })
        .map(|f| f as u32)
        .collect()
}

/// The supremum `τ_max` of the strengths for which every face satisfies the strict triangle
/// inequality: all faces are valid on `[0, τ_max)` (the inequalities are affine in `τ` and hold
/// at `0`), exactly, as the minimum over faces and sides of
/// `(ℓ_y + ℓ_z − ℓ_x) / ((d_y + d_z) − d_x)`, `d_k = ℓ_k (1 − a_k)`. `+∞` if no side limits it.
/// Not capped at `1`.
pub fn max_valid_tau(mesh: &TriMesh, lengths: &[f64], alignment: &[f64]) -> f64 {
    let mut tau_max = f64::INFINITY;
    for f in 0..mesh.n_faces() {
        let e = [0, 1, 2].map(|k| mesh.edge(3 * f + k));
        let l = e.map(|e| lengths[e]);
        let d = e.map(|e| lengths[e] * (1.0 - alignment[e]));
        for x in 0..3 {
            let (y, z) = ((x + 1) % 3, (x + 2) % 3);
            let g0 = l[y] + l[z] - l[x];
            let slope = d[y] + d[z] - d[x];
            if slope > 0.0 {
                tau_max = tau_max.min(g0 / slope);
            }
        }
    }
    tau_max
}

/// Scales the edge lengths anisotropically (paper Sec. 8.6.2; module docs).
///
/// Returns the scaled lengths (`TriMesh::edges` order) and a report. `τ = 0` returns the input
/// lengths bitwise. If the scaled lengths violate the triangle inequality: an error listing the
/// faces, or with [`AnisotropyOptions::clamp_tau`] the largest valid `τ` ([`max_valid_tau`]
/// times `1 − clamp_safety`), reported in the [`AnisotropyReport`].
///
/// # Errors
/// [`AnisotropyError::BadOption`] (`τ` not in `[0, 1]`, `clamp_safety` not in `[0, 1)`),
/// [`AnisotropyError::BadInput`] (see [`edge_alignment`]),
/// [`AnisotropyError::TriangleInequality`].
///
/// [formula: iem-anisotropic-scale]
pub fn anisotropic_lengths(
    mesh: &TriMesh,
    lengths: &[f64],
    angles: &[f64],
    options: &AnisotropyOptions,
) -> Result<(Vec<f64>, AnisotropyReport), AnisotropyError> {
    let tau = options.tau;
    if !(tau.is_finite() && (0.0..=1.0).contains(&tau)) {
        return Err(AnisotropyError::BadOption {
            reason: format!("tau = {tau} is not in [0, 1]"),
        });
    }
    if !(options.clamp_safety.is_finite() && (0.0..1.0).contains(&options.clamp_safety)) {
        return Err(AnisotropyError::BadOption {
            reason: format!("clamp_safety = {} is not in [0, 1)", options.clamp_safety),
        });
    }
    let alignment = edge_alignment(mesh, lengths, angles)?;
    let tau_max = max_valid_tau(mesh, lengths, &alignment);
    let bad = violating_faces(mesh, lengths, &alignment, tau);
    let n_violating = bad.len();
    let mut used = tau;
    if !bad.is_empty() {
        if !options.clamp_tau {
            return Err(AnisotropyError::TriangleInequality {
                tau,
                tau_max,
                faces: bad,
            });
        }
        used = tau.min(tau_max * (1.0 - options.clamp_safety));
        // The exact root can be off by rounding: back off until valid.
        let mut tries = 0;
        while !violating_faces(mesh, lengths, &alignment, used).is_empty() {
            tries += 1;
            if tries > 60 {
                return Err(AnisotropyError::BadInput {
                    reason: format!("no valid tau found below {used} (tau_max = {tau_max})"),
                });
            }
            used *= 1.0 - 1e-9 * f64::from(1u32 << tries.min(30));
        }
    }
    let factors: Vec<f64> = alignment.iter().map(|&a| length_factor(used, a)).collect();
    let scaled: Vec<f64> = lengths.iter().zip(&factors).map(|(l, f)| l * f).collect();
    let (mut lo, mut hi, mut sum) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
    for &f in &factors {
        lo = lo.min(f);
        hi = hi.max(f);
        sum += f;
    }
    let n = factors.len().max(1) as f64;
    Ok((
        scaled,
        AnisotropyReport {
            tau_requested: tau,
            tau_used: used,
            clamped: used < tau,
            tau_max,
            n_violating_faces_requested: n_violating,
            min_factor: lo,
            max_factor: hi,
            mean_factor: sum / n,
        },
    ))
}

/// A coarsenable mesh with the scaled metric, and what was done to `τ`.
///
/// The `CoarseMesh` is [`CoarseMesh::from_lengths`] of the **scaled** lengths with the tracked
/// points given on the **input** faces: a point is a face of the input mesh and barycentric
/// coordinates w.r.t. its corners, which do not depend on the metric, so the tracked fine
/// vertices and the extra points keep referring to the input surface while the metric being
/// coarsened is the scaled one (SPEC §8). Everything that follows (flattening, removal,
/// Delaunay flips, [`CoarseMesh::compact`], the prolongation `P`) works on the scaled metric;
/// `P` is the barycentric map between the scaled coarse and the (scaled) fine mesh, to be
/// applied to functions on the input vertices.
///
/// # Errors
/// As [`anisotropic_lengths`], and [`AnisotropyError::Coarse`] for a bad extra point.
pub fn anisotropic_coarse_mesh(
    mesh: &TriMesh,
    lengths: &[f64],
    angles: &[f64],
    points: &[InputPoint],
    mesh_options: CoarseMeshOptions,
    options: &AnisotropyOptions,
) -> Result<(CoarseMesh, AnisotropyReport), AnisotropyError> {
    let (scaled, report) = anisotropic_lengths(mesh, lengths, angles, options)?;
    let coarse = CoarseMesh::from_lengths(mesh, scaled, points, mesh_options)?;
    Ok((coarse, report))
}

// ----- principal directions -----

/// Which principal direction to emphasize.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrincipalDirection {
    /// The direction of the larger principal curvature `k_max`.
    Max,
    /// The direction of the smaller principal curvature `k_min`.
    Min,
}

/// A direction field from the principal curvature directions.
#[derive(Debug, Clone, PartialEq)]
pub struct PrincipalDirectionField {
    /// One normalized frame angle per vertex ([`direction_angles`]).
    pub angles: Vec<f64>,
    /// The 3D unit directions (lines) that were converted.
    pub directions: Vec<[f64; 3]>,
    /// The curvatures they come from.
    pub curvatures: PrincipalCurvatures,
}

/// Runs the principal curvature estimator ([`principal_curvatures`]) and converts the chosen
/// direction to angles in the normalized frames ([`direction_angles`]).
///
/// # Errors
/// [`AnisotropyError::Mesh`], [`AnisotropyError::Curvature`].
pub fn principal_direction_field(
    positions: &[[f64; 3]],
    faces: &[[u32; 3]],
    which: PrincipalDirection,
    options: &CurvatureOptions,
) -> Result<PrincipalDirectionField, AnisotropyError> {
    let mesh = TriMesh::new(positions.len(), faces.to_vec())?;
    let curvatures = principal_curvatures(&mesh, positions, options)?;
    let directions = match which {
        PrincipalDirection::Max => curvatures.d_max.clone(),
        PrincipalDirection::Min => curvatures.d_min.clone(),
    };
    let angles = direction_angles(&mesh, positions, &directions)?;
    Ok(PrincipalDirectionField {
        angles,
        directions,
        curvatures,
    })
}

// ----- a stand-in removal order and the elongation measurement -----

/// Coarsens `mesh` to at most `target` vertices by repeatedly removing the vertices of smallest
/// lumped area (`⅓ Σ A_f` in the current metric), in rounds: each round removes the
/// `max(1, n/6)` lightest vertices (at most `n − target`), skipping a vertex whose neighbour
/// was removed earlier in the round (its mass is stale). It is a simple, memoryless stand-in
/// for the ICE metric's order (`simplify`, WP3), independent of it: with the scaled lengths
/// it yields triangles of roughly equal *scaled* area. Vertices refused by the removal
/// (SPEC §3.5) are skipped. Returns the number of vertices left.
///
/// # Errors
/// An [`IemError`] of the removal other than a refusal.
pub fn coarsen_by_area(mesh: &mut CoarseMesh, target: usize) -> Result<usize, IemError> {
    loop {
        let n = mesh.n_vertices();
        if n <= target {
            return Ok(n);
        }
        let mass = mesh.lumped_mass();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| mass[a].total_cmp(&mass[b]).then(a.cmp(&b)));
        let budget = (n - target).min((n / 6).max(1));
        let ids: Vec<u32> = order.iter().map(|&i| mesh.vertices()[i]).collect();
        let mut stale = std::collections::HashSet::new();
        let mut removed = 0;
        for v in ids {
            if removed >= budget {
                break;
            }
            if stale.contains(&v) || !mesh.is_vertex(v) || mesh.is_fixed(v) {
                continue;
            }
            match mesh.remove_vertex(v) {
                Ok(r) => {
                    removed += 1;
                    stale.extend(r.prepared.neighbours.iter().map(|c| c.vertex));
                }
                Err(IemError::Removal(_)) => {}
                Err(e) => return Err(e),
            }
        }
        if removed == 0 {
            return Ok(mesh.n_vertices());
        }
    }
}

/// How well the coarse faces of a [`crate::CoarseOutput`] are elongated along (or across)
/// given directions (Fig. 27); see [`coarse_face_alignment`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElongationStats {
    /// Coarse faces that were measured (at least `min_samples` samples).
    pub n_faces: usize,
    /// Mean over those faces of the angle in degrees between the long axis of the face and
    /// the direction field (mean over the face's samples of the angle between the axis,
    /// projected on the tangent plane, and the line `d_v`), in `[0, 90]`. A field-independent
    /// face averages about 45; elongation **across** the field gives 90, **along** it 0.
    pub mean_angle_deg: f64,
    /// Mean over the faces of `√(λ₁/λ₂)`, the ratio of the standard deviations of the samples
    /// of the face along its long axis and along the orthogonal tangent direction.
    pub mean_aspect: f64,
}

/// Measures the elongation of the coarse faces in the input geometry. The fine vertices
/// (tracked points `0..n_fine`) are assigned to coarse faces by the tracking map
/// ([`crate::CoarseOutput::face_of_points`]); the samples of a face are those fine vertices plus
/// its three corners. The **long axis** of a face is the dominant eigenvector of the 2x2
/// covariance of its samples' positions projected on the tangent plane at their mean normal;
/// `λ₁/λ₂` are its eigenvalues. Faces with fewer than `min_samples` samples are skipped.
///
/// `positions`, `directions` (unit lines) and `normals` are per fine vertex.
pub fn coarse_face_alignment(
    out: &crate::CoarseOutput,
    positions: &[[f64; 3]],
    directions: &[[f64; 3]],
    normals: &[[f64; 3]],
    min_samples: usize,
) -> ElongationStats {
    type V = [f64; 3];
    let dot = |a: V, b: V| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cr = |a: V, b: V| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let n_fine = out
        .n_fine_vertices
        .min(positions.len())
        .min(directions.len())
        .min(normals.len());
    let fine_faces = out.faces_fine();
    let mut samples: Vec<Vec<u32>> = fine_faces.iter().map(|f| f.to_vec()).collect();
    for (p, f) in out.face_of_points().into_iter().enumerate().take(n_fine) {
        if let Some(f) = f {
            samples[f as usize].push(p as u32);
        }
    }
    let (mut sum_angle, mut sum_aspect, mut count) = (0.0, 0.0, 0usize);
    for s in &mut samples {
        s.retain(|&v| (v as usize) < n_fine);
        s.sort_unstable();
        s.dedup();
        if s.len() < min_samples.max(3) {
            continue;
        }
        let k = s.len() as f64;
        let mean: V =
            [0, 1, 2].map(|c| s.iter().map(|&v| positions[v as usize][c]).sum::<f64>() / k);
        let mut nrm: V = [0.0; 3];
        for &v in s.iter() {
            for c in 0..3 {
                nrm[c] += normals[v as usize][c];
            }
        }
        let nl = dot(nrm, nrm).sqrt();
        if !(nl > 1e-12) {
            continue;
        }
        let nrm = nrm.map(|x| x / nl);
        let helper: V = if nrm[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let t1 = cr(nrm, helper);
        let l1 = dot(t1, t1).sqrt();
        let t1 = t1.map(|x| x / l1);
        let t2 = cr(nrm, t1);
        let (mut cxx, mut cxy, mut cyy) = (0.0, 0.0, 0.0);
        for &v in s.iter() {
            let d: V = [0, 1, 2].map(|c| positions[v as usize][c] - mean[c]);
            let (x, y) = (dot(d, t1), dot(d, t2));
            cxx += x * x;
            cxy += x * y;
            cyy += y * y;
        }
        let tr = cxx + cyy;
        let disc = (0.25 * (cxx - cyy).powi(2) + cxy * cxy).sqrt();
        let (lam1, lam2) = (0.5 * tr + disc, (0.5 * tr - disc).max(0.0));
        if !(lam1 > 0.0) {
            continue;
        }
        let ang = 0.5 * (2.0 * cxy).atan2(cxx - cyy);
        let axis: V = [0, 1, 2].map(|c| ang.cos() * t1[c] + ang.sin() * t2[c]);
        let mut a_sum = 0.0;
        let mut a_n = 0.0;
        for &v in s.iter() {
            let n = normals[v as usize];
            let mut a: V = [0, 1, 2].map(|c| axis[c] - dot(axis, n) * n[c]);
            let al = dot(a, a).sqrt();
            if !(al > 1e-12) {
                continue;
            }
            a = a.map(|x| x / al);
            a_sum += dot(a, directions[v as usize]).abs().min(1.0).acos();
            a_n += 1.0;
        }
        if a_n == 0.0 {
            continue;
        }
        sum_angle += (a_sum / a_n).to_degrees();
        sum_aspect += (lam1 / lam2.max(1e-300)).sqrt();
        count += 1;
    }
    let c = count.max(1) as f64;
    ElongationStats {
        n_faces: count,
        mean_angle_deg: sum_angle / c,
        mean_aspect: sum_aspect / c,
    }
}
