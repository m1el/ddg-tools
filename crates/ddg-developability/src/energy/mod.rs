//! The six developability energies and their (sub)gradients (SPEC §3), and their evaluation on
//! a mesh (SPEC §3.9).
//!
//! `E(f) = Σ_{v ∈ D} e_v(f)` where `e_v` depends only on the star of `v`. Every `e_v` is a
//! function of the unit normals `N_σ`, the corner angles `θ` at `v` and (intrinsic energies) the
//! vector area `S_v` of the star; it is computed by a *star function* that also returns the
//! sensitivities `y_c = ∂e_v/∂N_σ`, `s_c = ∂e_v/∂θ_c` per corner of the star and
//! `z_v = ∂e_v/∂S_v`. [`crate::geometry::backprop_face`] then turns them into positional
//! gradients face by face.

mod combinatorial;
mod covariance;
mod intrinsic;

pub use combinatorial::{bipartitions, combinatorial, combinatorial_max};
pub use covariance::{
    circumcenter_value_gradient, covariance, covariance_max, spherical_circumcenter,
};
pub use intrinsic::{
    LOG_MAP_SINGULAR_ANGLE, LogMap, REASON_ANTIPODAL_NORMAL, REASON_ZERO_VECTOR_AREA, intrinsic,
    intrinsic_max, log_map, log_map_adjoint, log_map_g, log_map_h, tangent_basis, width_direction,
};

use ddg_core::TriMesh;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, Ordering};

use crate::error::DevError;
use crate::geometry::{FaceFrame, backprop_face, corner_angles, face_frame, v3};

/// Which per-vertex energy (SPEC §3.2–§3.7).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Energy {
    /// `E_P` (Eq. 1–2): minimum over bipartitions of the summed normal variance of the parts.
    Combinatorial,
    /// `π_max` (Eq. 5): minimum over bipartitions of the largest normal distance within a part.
    CombinatorialMax,
    /// `E_λ` (Eq. 3–4): smallest eigenvalue of the angle-weighted normal covariance.
    Covariance,
    /// `λ_max` (Eq. 6): `min_{|u|=1} max_σ ⟨u, N_σ⟩²`.
    CovarianceMax,
    /// Intrinsic width (App. B.5): `E_λ` of the log-mapped normals at the vertex normal.
    Intrinsic,
    /// Branching-free intrinsic width (App. B.5.1): squared minimum width of the log-mapped normals.
    IntrinsicMax,
}

impl Energy {
    /// All six energies, in declaration order.
    pub const ALL: [Energy; 6] = [
        Energy::Combinatorial,
        Energy::CombinatorialMax,
        Energy::Covariance,
        Energy::CovarianceMax,
        Energy::Intrinsic,
        Energy::IntrinsicMax,
    ];

    /// Whether the energy uses the vertex normal (and can be undefined, D12).
    pub fn is_intrinsic(self) -> bool {
        matches!(self, Energy::Intrinsic | Energy::IntrinsicMax)
    }
}

/// Options of an energy evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyOptions {
    /// The per-vertex energy.
    pub energy: Energy,
    /// Count boundary vertices (default `true`: the paper sums over all vertices, D8).
    pub include_boundary: bool,
    /// Leave out interior vertices with exactly 3 faces (default `true`; Sec. 4.1.4, D8).
    pub omit_valence_three: bool,
    /// Number of worker threads: `Some(n)` runs in a local pool of `n` threads (unless the call
    /// already runs inside a rayon pool of exactly `n` threads, which is then used); `None` runs
    /// on the current rayon pool (the caller's, or rayon's default pool). Results are bitwise
    /// identical for every thread count.
    pub threads: Option<usize>,
}

impl Default for EnergyOptions {
    fn default() -> Self {
        EnergyOptions {
            energy: Energy::Covariance,
            include_boundary: true,
            omit_valence_three: true,
            threads: None,
        }
    }
}

impl EnergyOptions {
    /// Options with the given energy and the other fields at their defaults.
    pub fn new(energy: Energy) -> Self {
        EnergyOptions {
            energy,
            ..Default::default()
        }
    }

    /// Check the options.
    ///
    /// # Errors
    /// [`DevError::InvalidOption`] if `threads == Some(0)`.
    pub fn validate(&self) -> Result<(), DevError> {
        if self.threads == Some(0) {
            return Err(DevError::InvalidOption {
                name: "threads",
                value: 0.0,
                requirement: "threads >= 1 (or None for the current pool)",
            });
        }
        Ok(())
    }
}

/// Result of [`evaluate`].
#[derive(Debug, Clone, PartialEq)]
pub struct Evaluation {
    /// `E = Σ_{v ∈ D} e_v`, summed in vertex order.
    pub energy: f64,
    /// `e_v` for every vertex, also for vertices outside the domain (`NaN` where the formula is
    /// undefined at such a vertex).
    pub vertex_energy: Vec<f64>,
    /// Whether the vertex is in the domain `D` (counts toward `energy` and `gradient`).
    pub counted: Vec<bool>,
    /// `∇E` per vertex.
    pub gradient: Vec<[f64; 3]>,
    /// Margin of the active piece per vertex (SPEC §3.8): second-best minus best key value, or the
    /// eigenvalue gap over the trace; `+∞` when there is no alternative; `NaN` where the energy is
    /// undefined.
    pub vertex_margin: Vec<f64>,
    /// Number of counted vertices with `margin ≤ 1e-12 · max(1, e_v)`.
    pub near_ties: usize,
}

/// Near-tie threshold of SPEC §3.8 (relative to `max(1, e_v)`).
pub const NEAR_TIE: f64 = 1e-12;

/// The star of one vertex, in the order of `TriMesh::outgoing_halfedges`.
#[derive(Debug, Clone, Copy)]
pub struct Star<'a> {
    /// Unit normals `N_σ` of the faces of the star.
    pub normals: &'a [[f64; 3]],
    /// Corner angles `θ_v^σ` at the vertex.
    pub angles: &'a [f64],
    /// `true` for an interior vertex (the faces form a cycle), `false` on the boundary (a path).
    pub closed: bool,
    /// Vector area `S_v = ½ Σ n_σ` (used by the intrinsic energies).
    pub vector_area: [f64; 3],
}

/// Value of a star function.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StarOut {
    /// `e_v`.
    pub value: f64,
    /// Margin of the active piece (SPEC §3.8).
    pub margin: f64,
    /// `z_v = ∂e_v/∂S_v` (zero for the extrinsic energies).
    pub z: [f64; 3],
}

/// Reusable buffers of the star functions.
#[derive(Default)]
pub struct Scratch {
    f1: Vec<f64>,
    v1: Vec<[f64; 3]>,
    logs: Vec<LogMap>,
    arcs: Vec<combinatorial::ArcMaxSlot>,
}

/// Evaluate the energy `energy` on one star: returns `e_v`, its margin and `z_v`, and writes the
/// sensitivities `y[m] = ∂e_v/∂N_σm`, `s[m] = ∂e_v/∂θ_m` (`y`, `s` have the star's length).
/// `Err(reason)` where the energy is undefined (intrinsic energies, D12).
///
/// The margin of a boundary star with at most two faces is `+∞`: every energy vanishes
/// identically there, so there is no tie to report (for `E_λ` with one face the two zero
/// eigenvalues would otherwise count as one).
pub fn star_energy(
    energy: Energy,
    star: &Star,
    y: &mut [[f64; 3]],
    s: &mut [f64],
    scratch: &mut Scratch,
) -> Result<StarOut, &'static str> {
    debug_assert!(star.normals.len() == star.angles.len());
    debug_assert!(y.len() == star.normals.len() && s.len() == star.normals.len());
    let mut out = match energy {
        Energy::Combinatorial => combinatorial(star, y, s),
        Energy::CombinatorialMax => combinatorial_max(star, y, s, scratch),
        Energy::Covariance => covariance(star, y, s)?,
        Energy::CovarianceMax => covariance_max(star, y, s),
        Energy::Intrinsic => intrinsic(star, y, s, scratch)?,
        Energy::IntrinsicMax => intrinsic_max(star, y, s, scratch)?,
    };
    if !star.closed && star.normals.len() <= 2 {
        // The energy vanishes identically on such stars (SPEC §3.1): no nonsmoothness.
        out.margin = f64::INFINITY;
    }
    Ok(out)
}

/// The domain `D` (SPEC §3.1, D8): `counted[v]` is false for boundary vertices when
/// `!include_boundary`, and for interior vertices with exactly three faces when
/// `omit_valence_three` (a valence-3 hinge is flat, Prop. A.3, so these "triple points" are left
/// free, Sec. 4.1.4).
///
/// [formula: dev-valence-three-flat]
pub fn domain(mesh: &TriMesh, include_boundary: bool, omit_valence_three: bool) -> Vec<bool> {
    (0..mesh.n_vertices())
        .into_par_iter()
        .with_min_len(VERTEX_CHUNK)
        .map(|v| {
            let closed = !mesh.is_boundary_vertex(v);
            let k = if closed && omit_valence_three {
                mesh.vertex_face_count(v)
            } else {
                0
            };
            in_domain(closed, k, include_boundary, omit_valence_three)
        })
        .collect()
}

/// The rule of [`domain`] for one vertex (`k` = number of faces).
#[inline]
fn in_domain(closed: bool, k: usize, include_boundary: bool, omit_valence_three: bool) -> bool {
    if closed {
        !(omit_valence_three && k == 3)
    } else {
        include_boundary
    }
}

/// Vertices per parallel work item. Fixed (not derived from the thread count) so that the
/// results do not depend on it.
const VERTEX_CHUNK: usize = 256;

/// Star layout: the corners of each vertex in `outgoing_halfedges` order (CSR) and the inverse map.
struct Stars {
    /// `offsets[v] .. offsets[v + 1]` are the slots of vertex `v`.
    offsets: Vec<usize>,
    /// Corner id (`3 f + slot`) per star slot.
    corners: Vec<u32>,
    /// Star slot per corner id (each written once, by the chunk owning the corner's vertex).
    slot: Vec<AtomicU32>,
    /// Interior vertex (closed fan).
    closed: Vec<bool>,
}

/// Split `data` into consecutive pieces of the given lengths.
fn split_lengths<T>(mut data: &mut [T], lengths: impl Iterator<Item = usize>) -> Vec<&mut [T]> {
    let mut out = Vec::new();
    for len in lengths {
        let (a, b) = std::mem::take(&mut data).split_at_mut(len);
        out.push(a);
        data = b;
    }
    out
}

impl Stars {
    /// One walk around every vertex (parallel over fixed vertex chunks), then the offsets, then
    /// the inverse map (parallel scatter; every corner belongs to exactly one star).
    fn new(mesh: &TriMesh) -> Stars {
        let nv = mesh.n_vertices();
        let nc = 3 * mesh.n_faces();
        let chunks: Vec<(Vec<u32>, Vec<u32>, Vec<bool>)> = (0..nv.div_ceil(VERTEX_CHUNK))
            .into_par_iter()
            .map(|ci| {
                let vs = ci * VERTEX_CHUNK..((ci + 1) * VERTEX_CHUNK).min(nv);
                let mut corners = Vec::with_capacity(7 * vs.len());
                let mut valence = Vec::with_capacity(vs.len());
                let mut closed = Vec::with_capacity(vs.len());
                for v in vs {
                    let before = corners.len();
                    corners.extend(mesh.outgoing_halfedges(v).map(|h| h as u32));
                    valence.push((corners.len() - before) as u32);
                    closed.push(!mesh.is_boundary_vertex(v));
                }
                (corners, valence, closed)
            })
            .collect();
        let mut offsets = Vec::with_capacity(nv + 1);
        offsets.push(0usize);
        let mut corners = Vec::with_capacity(nc);
        let mut closed = Vec::with_capacity(nv);
        for (c, val, cl) in &chunks {
            corners.extend_from_slice(c);
            closed.extend_from_slice(cl);
            for &k in val {
                offsets.push(offsets.last().copied().unwrap_or(0) + k as usize);
            }
        }
        let slot: Vec<AtomicU32> = (0..nc).map(|_| AtomicU32::new(0)).collect();
        (0..nv.div_ceil(VERTEX_CHUNK))
            .into_par_iter()
            .for_each(|ci| {
                let lo = offsets[ci * VERTEX_CHUNK];
                let hi = offsets[((ci + 1) * VERTEX_CHUNK).min(nv)];
                for (i, &c) in corners[lo..hi].iter().enumerate() {
                    slot[c as usize].store((lo + i) as u32, Ordering::Relaxed);
                }
            });
        Stars {
            offsets,
            corners,
            slot,
            closed,
        }
    }

    #[inline]
    fn slot(&self, corner: usize) -> usize {
        self.slot[corner].load(Ordering::Relaxed) as usize
    }
}

/// Number of star slots in each vertex chunk.
fn chunk_slot_lengths(offsets: &[usize], nv: usize) -> Vec<usize> {
    (0..nv.div_ceil(VERTEX_CHUNK))
        .map(|ci| offsets[((ci + 1) * VERTEX_CHUNK).min(nv)] - offsets[ci * VERTEX_CHUNK])
        .collect()
}

/// Per-vertex result of the star pass.
#[derive(Clone, Copy)]
struct VertexOut {
    value: f64,
    margin: f64,
    z: [f64; 3],
}

/// Run `f` with the requested number of threads (see [`EnergyOptions::threads`]).
pub(crate) fn with_threads<R: Send>(
    threads: Option<usize>,
    f: impl FnOnce() -> R + Send,
) -> Result<R, DevError> {
    match threads {
        None => Ok(f()),
        Some(n) => {
            if rayon::current_thread_index().is_some() && rayon::current_num_threads() == n {
                return Ok(f());
            }
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| {
                    DevError::Internal(format!("could not build a thread pool of {n} threads: {e}"))
                })?;
            Ok(pool.install(f))
        }
    }
}

/// Evaluate the energy, the per-vertex energies, the gradient and the margins (SPEC §3.9).
///
/// Steps (each parallel, results bitwise independent of the thread count, D25):
/// 1. per face: `n`, `|n|`, `N`, corner angles (a zero or non-finite `n` is
///    [`DevError::DegenerateFace`]);
/// 2. per star: `e_v` and the sensitivities `y_c`, `s_c`, `z_v` (uncounted vertices contribute
///    zero sensitivities);
/// 3. per face: [`backprop_face`] into three per-corner contributions;
/// 4. per vertex: `g_v = Σ` of its corners' contributions in star order; `E = Σ_{v∈D} e_v` in
///    vertex order.
///
/// All energies depend only on unit normals and angles, so they are invariant under rigid motions
/// and uniform scaling; hence `Σ_v ∇_v E = 0`, `Σ_v f_v × ∇_v E = 0` and `Σ_v f_v · ∇_v E = 0`.
///
/// # Errors
/// [`DevError::InvalidOption`], [`DevError::LengthMismatch`] (positions),
/// [`DevError::NonFinitePosition`], [`DevError::DegenerateFace`] (lowest face id),
/// [`DevError::UndefinedEnergy`] (intrinsic energies at a counted vertex; lowest vertex id).
///
/// [formula: dev-gradient-invariance]
pub fn evaluate(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &EnergyOptions,
) -> Result<Evaluation, DevError> {
    opts.validate()?;
    with_threads(opts.threads, || evaluate_impl(mesh, positions, opts, true))?
}

/// The energy and the per-vertex energies only (no gradient; cheaper). Same errors as
/// [`evaluate`].
pub fn evaluate_energy(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &EnergyOptions,
) -> Result<(f64, Vec<f64>), DevError> {
    opts.validate()?;
    with_threads(opts.threads, || evaluate_impl(mesh, positions, opts, false))?
        .map(|e| (e.energy, e.vertex_energy))
}

fn evaluate_impl(
    mesh: &TriMesh,
    positions: &[[f64; 3]],
    opts: &EnergyOptions,
    want_gradient: bool,
) -> Result<Evaluation, DevError> {
    let nv = mesh.n_vertices();
    let nf = mesh.n_faces();
    if positions.len() != nv {
        return Err(DevError::LengthMismatch {
            what: "positions",
            expected: nv,
            got: positions.len(),
        });
    }
    if let Some(v) = positions.par_iter().position_first(|p| !v3::is_finite(*p)) {
        return Err(DevError::NonFinitePosition { vertex: v as u32 });
    }
    let faces = mesh.faces();
    let corner_pos = |f: usize| {
        let t = faces[f];
        [
            positions[t[0] as usize],
            positions[t[1] as usize],
            positions[t[2] as usize],
        ]
    };

    // 1. Faces.
    let bad = FaceFrame {
        n: [f64::NAN; 3],
        n_len: f64::NAN,
        normal: [f64::NAN; 3],
    };
    let (frames, angles): (Vec<FaceFrame>, Vec<[f64; 3]>) = (0..nf)
        .into_par_iter()
        .with_min_len(1024)
        .map(|f| {
            let p = corner_pos(f);
            match face_frame(p) {
                Some(fr) => (fr, corner_angles(p)),
                None => (bad, [f64::NAN; 3]),
            }
        })
        .unzip();
    if let Some(f) = frames.par_iter().position_first(|fr| fr.n_len.is_nan()) {
        return Err(DevError::DegenerateFace { face: f as u32 });
    }

    // 2. Stars.
    let stars = Stars::new(mesh);
    let counted: Vec<bool> = (0..nv)
        .into_par_iter()
        .with_min_len(VERTEX_CHUNK)
        .map(|v| {
            let k = stars.offsets[v + 1] - stars.offsets[v];
            in_domain(
                stars.closed[v],
                k,
                opts.include_boundary,
                opts.omit_valence_three,
            )
        })
        .collect();
    let n_slots = stars.corners.len();
    let (mut ys, mut ss) = if want_gradient {
        (vec![[0.0; 3]; n_slots], vec![0.0; n_slots])
    } else {
        (Vec::new(), Vec::new())
    };
    let mut vout = vec![
        VertexOut {
            value: 0.0,
            margin: 0.0,
            z: [0.0; 3]
        };
        nv
    ];
    // Per-chunk slot ranges of the sensitivity arrays (empty without a gradient).
    let lens: Vec<usize> = chunk_slot_lengths(&stars.offsets, nv)
        .into_iter()
        .map(|l| if want_gradient { l } else { 0 })
        .collect();
    let y_chunks = split_lengths(&mut ys, lens.iter().copied());
    let s_chunks = split_lengths(&mut ss, lens.iter().copied());
    let energy = opts.energy;
    let errors: Vec<Option<DevError>> = vout
        .par_chunks_mut(VERTEX_CHUNK)
        .zip(y_chunks)
        .zip(s_chunks)
        .enumerate()
        .map(|(ci, ((vo, yc), sc))| {
            let mut scratch = Scratch::default();
            let mut normals: Vec<[f64; 3]> = Vec::new();
            let mut ang: Vec<f64> = Vec::new();
            let mut ybuf: Vec<[f64; 3]> = Vec::new();
            let mut sbuf: Vec<f64> = Vec::new();
            let v0 = ci * VERTEX_CHUNK;
            let base = stars.offsets[v0];
            for (i, out) in vo.iter_mut().enumerate() {
                let v = v0 + i;
                let (lo, hi) = (stars.offsets[v], stars.offsets[v + 1]);
                normals.clear();
                ang.clear();
                let mut sv = [0.0; 3];
                for &c in &stars.corners[lo..hi] {
                    let c = c as usize;
                    let fr = &frames[c / 3];
                    normals.push(fr.normal);
                    ang.push(angles[c / 3][c % 3]);
                    sv = v3::add(sv, fr.n);
                }
                let star = Star {
                    normals: &normals,
                    angles: &ang,
                    closed: stars.closed[v],
                    vector_area: v3::scale(sv, 0.5),
                };
                let k = hi - lo;
                let (y, s): (&mut [[f64; 3]], &mut [f64]) = if want_gradient {
                    (&mut yc[lo - base..hi - base], &mut sc[lo - base..hi - base])
                } else {
                    ybuf.resize(k, [0.0; 3]);
                    sbuf.resize(k, 0.0);
                    (&mut ybuf[..k], &mut sbuf[..k])
                };
                match star_energy(energy, &star, y, s, &mut scratch) {
                    Ok(r) => {
                        *out = VertexOut {
                            value: r.value,
                            margin: r.margin,
                            z: r.z,
                        };
                        if !counted[v] {
                            y.fill([0.0; 3]);
                            s.fill(0.0);
                            out.z = [0.0; 3];
                        }
                    }
                    Err(reason) => {
                        if counted[v] {
                            return Some(DevError::UndefinedEnergy {
                                vertex: v as u32,
                                reason,
                            });
                        }
                        *out = VertexOut {
                            value: f64::NAN,
                            margin: f64::NAN,
                            z: [0.0; 3],
                        };
                        y.fill([0.0; 3]);
                        s.fill(0.0);
                    }
                }
            }
            None
        })
        .collect();
    if let Some(e) = errors.into_iter().flatten().next() {
        return Err(e);
    }

    let mut total = 0.0;
    let mut near_ties = 0;
    for (o, &c) in vout.iter().zip(&counted) {
        if c {
            total += o.value;
            if o.margin <= NEAR_TIE * o.value.max(1.0) {
                near_ties += 1;
            }
        }
    }
    let vertex_energy: Vec<f64> = vout.iter().map(|o| o.value).collect();
    let vertex_margin: Vec<f64> = vout.iter().map(|o| o.margin).collect();

    let gradient = if want_gradient {
        // 3. Faces: back-propagation.
        let q: Vec<[[f64; 3]; 3]> = (0..nf)
            .into_par_iter()
            .with_min_len(1024)
            .map(|f| {
                let t = faces[f];
                let mut y = [0.0; 3];
                let mut z = [0.0; 3];
                let mut s = [0.0; 3];
                for m in 0..3 {
                    let sl = stars.slot(3 * f + m);
                    y = v3::add(y, ys[sl]);
                    s[m] = ss[sl];
                    z = v3::add(z, vout[t[m] as usize].z);
                }
                backprop_face(corner_pos(f), &frames[f], y, z, s)
            })
            .collect();
        // 4. Vertices: gather in star order.
        (0..nv)
            .into_par_iter()
            .with_min_len(VERTEX_CHUNK)
            .map(|v| {
                let mut g = [0.0; 3];
                for &c in &stars.corners[stars.offsets[v]..stars.offsets[v + 1]] {
                    let c = c as usize;
                    g = v3::add(g, q[c / 3][c % 3]);
                }
                g
            })
            .collect()
    } else {
        Vec::new()
    };

    Ok(Evaluation {
        energy: total,
        vertex_energy,
        counted,
        gradient,
        vertex_margin,
        near_ties,
    })
}
