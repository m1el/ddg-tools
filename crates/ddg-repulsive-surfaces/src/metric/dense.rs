//! Dense assembly of `L_τ`, `B`, `B⁰`, `A` and the L²/H¹/H² metrics (`SPEC.md` §4, Alg. 1).
//!
//! The three operators of Eq. 11–13 are double sums over ordered face pairs `S ≠ T`
//! (`ρ = |X_S − X_T|`):
//!
//! ```text
//! uᵀ L_τ v = Σ_{S≠T} (ū_S − ū_T)(v̄_S − v̄_T) a_S a_T ρ^{−(2τ+2)}
//! uᵀ B   v = Σ_{S≠T} ⟨D_f u(S) − D_f u(T), D_f v(S) − D_f v(T)⟩ a_S a_T ρ^{−2s}
//! uᵀ B⁰  v = Σ_{S≠T} (ū_S − ū_T)(v̄_S − v̄_T) K₂(S, T) a_S a_T ρ^{−2s},  K₂ = |N_S·(X_S − X_T)|² ρ^{−4}
//! ```
//!
//! Two assemblies are provided and agree to round-off (tests):
//!
//! * [`kernel_matrix_form`] (used by [`dense_metric`]): App. A's identity (§4.3)
//!   `2 Xᵀ (diag(a)⁻¹ diag(H a) − H) X` with the symmetric kernel matrix `H` (`B⁰` with the
//!   symmetrized `½(K₂(S,T) + K₂(T,S))`, `rs-low-order-symmetrisation`), `X = U` or `V`. It is
//!   `O(|F|² + |F||V|)` per matrix and is exactly what the hierarchical product of WP6 compresses,
//!   so WP6 checks its products against [`kernel_matrix_apply`] with the same `H`.
//! * [`pair_sum_matrix`]: Alg. 1 literally (every ordered pair, the at most six vertices of the
//!   pair), with any (not necessarily symmetric) kernel; the independent reference for tests.
//!
//! All parallel loops write disjoint rows in a fixed order: results are bitwise independent of
//! the size of the pool.

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use ddg_core::SparseMatrix;
use rayon::prelude::*;

use super::{Channel, DenseMatrix, FaceDiagonal, FaceOperators, Metric, L2Kind, SobolevWeights};
use crate::error::RsError;
use crate::geometry::v3::{dot, sub};
use crate::geometry::{FaceGeometry, cotan_stiffness, vertex_areas};
use crate::kernel::Exponents;
use crate::surface::Surface;

fn zero_matrix(n: usize) -> DenseMatrix {
    DenseMatrix::zeros(n)
}

/// `ρ² = |X_S − X_T|²`, or [`RsError::NonFiniteEnergy`] for coincident barycenters.
#[inline]
fn rho2(ops: &FaceOperators, s: usize, t: usize) -> Result<f64, RsError> {
    let d = sub(ops.barycenters[s], ops.barycenters[t]);
    let r2 = dot(d, d);
    if r2 > 0.0 && r2.is_finite() {
        Ok(r2)
    } else {
        Err(RsError::NonFiniteEnergy {
            faces: [s as u32, t as u32],
        })
    }
}

/// The `|F| × |F|` kernel matrix `H_ST = h(S, T)` for `S ≠ T`, zero diagonal, computed for
/// `S < T` and mirrored (so it is exactly symmetric), parallel over rows.
fn symmetric_kernel(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    h: impl Fn(usize, usize, f64) -> f64 + Sync,
) -> Result<DenseMatrix, RsError> {
    let nf = ops.n_faces();
    let mut m = zero_matrix(nf);
    // Check coincident barycenters first (deterministic: the first pair in row-major order).
    for s in 0..nf {
        for t in s + 1..nf {
            rho2(ops, s, t)?;
        }
    }
    pool.install(|| {
        m.as_mut_slice()
            .par_chunks_mut(nf.max(1))
            .enumerate()
            .for_each(|(s, row)| {
                for (t, x) in row.iter_mut().enumerate() {
                    if t == s {
                        continue;
                    }
                    let (a, b) = if s < t { (s, t) } else { (t, s) };
                    let d = sub(ops.barycenters[a], ops.barycenters[b]);
                    *x = h(a, b, dot(d, d));
                }
            });
    });
    Ok(m)
}

/// The high-order kernel matrix `H_ST = ρ^{−2s}` of `B` (§4.2).
///
/// # Errors
/// [`RsError::NonFiniteEnergy`] for two faces with the same barycenter.
pub fn high_order_kernel(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    exponents: &Exponents,
) -> Result<DenseMatrix, RsError> {
    let e = -exponents.s();
    symmetric_kernel(pool, ops, |_, _, r2| r2.powf(e))
}

/// The low-order kernel matrix `H_ST = ½(K₂(S,T) + K₂(T,S)) ρ^{−2s}` of `B⁰`,
/// `K₂(S,T) = |N_S·(X_S − X_T)|² ρ^{−4}` (fixed `α = 2`, `β = 4` whatever the energy's
/// exponents, §4.2). The ordered sum with `K₂(S,T)` equals the sum with this symmetric kernel.
///
/// [formula: rs-low-order-symmetrisation]
///
/// # Errors
/// [`RsError::NonFiniteEnergy`] for two faces with the same barycenter.
pub fn low_order_kernel(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    exponents: &Exponents,
) -> Result<DenseMatrix, RsError> {
    let e = -exponents.s();
    symmetric_kernel(pool, ops, |s, t, r2| {
        let d = sub(ops.barycenters[s], ops.barycenters[t]);
        let ns = dot(ops.normals[s], d);
        let nt = dot(ops.normals[t], d);
        0.5 * (ns * ns + nt * nt) / (r2 * r2) * r2.powf(e)
    })
}

/// The kernel matrix `H_ST = ρ^{−(2τ+2)}` of the fractional operator `L_τ` (order `2τ`, §4.2;
/// `τ = 2 − s` is the preconditioner's `L_{2−s}`, §6.3).
///
/// # Errors
/// [`RsError::NonFiniteEnergy`] for two faces with the same barycenter;
/// [`RsError::InvalidOption`] for a non-finite `τ`.
pub fn fractional_kernel(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    tau: f64,
) -> Result<DenseMatrix, RsError> {
    if !tau.is_finite() {
        return Err(RsError::InvalidOption {
            name: "tau",
            value: tau.to_string(),
            requirement: "a finite order",
        });
    }
    let e = -(tau + 1.0);
    symmetric_kernel(pool, ops, |_, _, r2| r2.powf(e))
}

/// Vertex → incident faces, ascending.
fn incident_faces(ops: &FaceOperators) -> Vec<Vec<u32>> {
    let mut out = vec![Vec::new(); ops.n_vertices];
    for (s, t) in ops.faces.iter().enumerate() {
        for &v in t {
            out[v as usize].push(s as u32);
        }
    }
    out
}

fn check_kernel(ops: &FaceOperators, h: &DenseMatrix) -> Result<(), RsError> {
    let nf = ops.n_faces();
    if h.n() != nf {
        return Err(RsError::InvalidOption {
            name: "kernel matrix",
            value: format!("{0} x {0}", h.n()),
            requirement: "|F| x |F|, matching the face operators",
        });
    }
    let scale = h.max_abs();
    for s in 0..nf {
        if h.get(s, s) != 0.0 {
            return Err(RsError::InvalidOption {
                name: "kernel matrix",
                value: format!("diagonal entry {s} = {:e}", h.get(s, s)),
                requirement: "a zero diagonal (S = T is excluded)",
            });
        }
    }
    let defect = h.symmetry_defect();
    if !(defect <= 1e-12 * scale) {
        return Err(RsError::InvalidOption {
            name: "kernel matrix",
            value: format!("symmetry defect {defect:e}"),
            requirement: "a symmetric kernel (symmetrize first, rs-low-order-symmetrisation)",
        });
    }
    Ok(())
}

/// `H a` (row sums weighted by the face areas).
fn kernel_times_areas(pool: &rayon::ThreadPool, ops: &FaceOperators, h: &DenseMatrix) -> Vec<f64> {
    let nf = ops.n_faces();
    pool.install(|| {
        (0..nf)
            .into_par_iter()
            .map(|s| h.row(s).iter().zip(&ops.areas).map(|(x, a)| x * a).sum())
            .collect()
    })
}

/// The dense `|V| × |V|` matrix of App. A's kernel-matrix form (`SPEC.md` §4.3)
///
/// ```text
/// 2 Xᵀ (diag(a)⁻¹ diag(H a) − H) X,     X = U (Average) or V (Gradient, summed over 3 channels)
/// ```
///
/// for a symmetric kernel matrix `H` (`|F| × |F|`) with zero diagonal. Its quadratic form is the
/// ordered double sum `Σ_{S,T} (ū_S − ū_T)(v̄_S − v̄_T) a_S H_ST a_T` (resp. with `D_f`).
/// Parallel over output rows; bitwise independent of the pool size.
///
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// [`RsError::InvalidOption`] if `h` is not `|F| × |F|`, has a nonzero diagonal entry or is
/// not symmetric (to `1e-12` relative).
pub fn kernel_matrix_form(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    h: &DenseMatrix,
    channel: Channel,
) -> Result<DenseMatrix, RsError> {
    check_kernel(ops, h)?;
    let nf = ops.n_faces();
    let nv = ops.n_vertices;
    let w = channel.width();
    let d = kernel_times_areas(pool, ops, h);
    // Y = H X, |F| × (|V| w), row-major; X row S has 3 vertices × w channels.
    let rows: Vec<[(usize, [f64; 3]); 3]> = (0..nf).map(|s| ops.row(s, channel)).collect();
    let mut y = vec![0.0; nf * nv * w];
    pool.install(|| {
        y.par_chunks_mut((nv * w).max(1))
            .enumerate()
            .for_each(|(s, ys)| {
                let hs = h.row(s);
                for (t, xt) in rows.iter().enumerate() {
                    let hst = hs[t];
                    if hst == 0.0 {
                        continue;
                    }
                    for &(v, c) in xt {
                        for k in 0..w {
                            ys[v * w + k] += hst * c[k];
                        }
                    }
                }
            });
    });
    let inc = incident_faces(ops);
    let mut out = zero_matrix(nv);
    pool.install(|| {
        out.as_mut_slice()
            .par_chunks_mut(nv.max(1))
            .enumerate()
            .for_each(|(p, row)| {
                for &s in &inc[p] {
                    let s = s as usize;
                    let xs = &rows[s];
                    let cp = xs.iter().find(|e| e.0 == p).map(|e| e.1).unwrap_or([0.0; 3]);
                    let ys = &y[s * nv * w..(s + 1) * nv * w];
                    for q in 0..nv {
                        let mut acc = 0.0;
                        for k in 0..w {
                            acc += cp[k] * ys[q * w + k];
                        }
                        row[q] -= 2.0 * acc;
                    }
                    let ds = d[s] / ops.areas[s];
                    for &(q, cq) in xs {
                        let mut acc = 0.0;
                        for k in 0..w {
                            acc += cp[k] * cq[k];
                        }
                        row[q] += 2.0 * ds * acc;
                    }
                }
            });
    });
    Ok(out)
}

/// The product of the kernel-matrix form with three columns (`|V| × 3`, one vertex per entry):
/// `2 Xᵀ (diag(a)⁻¹ diag(H a) − H) X u` for each column of `u` (`SPEC.md` §4.3, §8.3). This is
/// the operation WP6's hierarchical product approximates; with the same `H` it is exact.
///
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// As [`kernel_matrix_form`]; [`RsError::PositionCountMismatch`] if `u.len() != |V|`.
pub fn kernel_matrix_apply(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    h: &DenseMatrix,
    channel: Channel,
    u: &[[f64; 3]],
) -> Result<Vec<[f64; 3]>, RsError> {
    check_kernel(ops, h)?;
    if u.len() != ops.n_vertices {
        return Err(RsError::PositionCountMismatch {
            expected: ops.n_vertices,
            found: u.len(),
        });
    }
    let nf = ops.n_faces();
    let w = channel.width();
    let d = kernel_times_areas(pool, ops, h);
    // x = X u: per face, w channels × 3 columns.
    let rows: Vec<[(usize, [f64; 3]); 3]> = (0..nf).map(|s| ops.row(s, channel)).collect();
    let x: Vec<[[f64; 3]; 3]> = rows
        .iter()
        .map(|r| {
            let mut o = [[0.0; 3]; 3]; // o[k][col]
            for &(v, c) in r {
                for k in 0..w {
                    for col in 0..3 {
                        o[k][col] += c[k] * u[v][col];
                    }
                }
            }
            o
        })
        .collect();
    // z = (diag(d/a) − H) x.
    let z: Vec<[[f64; 3]; 3]> = pool.install(|| {
        (0..nf)
            .into_par_iter()
            .map(|s| {
                let hs = h.row(s);
                let mut o = [[0.0; 3]; 3];
                for (t, xt) in x.iter().enumerate() {
                    let hst = hs[t];
                    if hst == 0.0 {
                        continue;
                    }
                    for k in 0..w {
                        for col in 0..3 {
                            o[k][col] += hst * xt[k][col];
                        }
                    }
                }
                let ds = d[s] / ops.areas[s];
                for k in 0..w {
                    for col in 0..3 {
                        o[k][col] = ds * x[s][k][col] - o[k][col];
                    }
                }
                o
            })
            .collect()
    });
    // out = 2 Xᵀ z, faces in order.
    let mut out = vec![[0.0; 3]; ops.n_vertices];
    for (s, r) in rows.iter().enumerate() {
        for &(v, c) in r {
            for k in 0..w {
                for col in 0..3 {
                    out[v][col] += 2.0 * c[k] * z[s][k][col];
                }
            }
        }
    }
    Ok(out)
}

/// Alg. 1 literally: the `|V| × |V|` matrix of the ordered double sum
/// `Σ_{S≠T} ⟨X̂u(S) − X̂u(T), X̂v(S) − X̂v(T)⟩ a_S a_T k(S, T)` with `X̂ = ū` (Average) or
/// `D_f` (Gradient), for any kernel `k(S, T)` (not necessarily symmetric), looping over the
/// ordered pairs and the at most six vertices of each pair (`SPEC.md` §4.2). `O(36 |F|²)`;
/// the independent reference for [`kernel_matrix_form`] in tests.
///
/// [formula: rs-sobolev-forms]
///
/// # Errors
/// [`RsError::NonFiniteEnergy`] for two faces with the same barycenter.
pub fn pair_sum_matrix(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    channel: Channel,
    kernel: &(dyn Fn(usize, usize) -> f64 + Sync),
) -> Result<DenseMatrix, RsError> {
    let nf = ops.n_faces();
    let nv = ops.n_vertices;
    for s in 0..nf {
        for t in s + 1..nf {
            rho2(ops, s, t)?;
        }
    }
    let w = channel.width();
    let coef = |s: usize| -> [(usize, [f64; 3]); 3] {
        let t = ops.faces[s];
        let mut out = [(0usize, [0.0; 3]); 3];
        for c in 0..3 {
            out[c] = (
                t[c] as usize,
                match channel {
                    Channel::Average => [1.0 / 3.0, 0.0, 0.0],
                    Channel::Gradient => ops.df[s][c],
                },
            );
        }
        out
    };
    // Fixed chunks of S, each with its own accumulator, summed in chunk order.
    let n_chunks = nf.clamp(1, 16);
    let chunk = nf.div_ceil(n_chunks);
    let partial: Vec<Vec<f64>> = pool.install(|| {
        (0..n_chunks)
            .into_par_iter()
            .map(|ci| {
                let mut acc = vec![0.0; nv * nv];
                let mut es: Vec<(usize, [f64; 3])> = Vec::with_capacity(6);
                for s in ci * chunk..((ci + 1) * chunk).min(nf) {
                    let cs = coef(s);
                    for t in 0..nf {
                        if t == s {
                            continue;
                        }
                        let wgt = ops.areas[s] * ops.areas[t] * kernel(s, t);
                        es.clear();
                        for (sign, cc) in [(1.0, cs), (-1.0, coef(t))] {
                            for (v, c) in cc {
                                let c = [sign * c[0], sign * c[1], sign * c[2]];
                                match es.iter_mut().find(|e| e.0 == v) {
                                    Some(e) => {
                                        for k in 0..3 {
                                            e.1[k] += c[k];
                                        }
                                    }
                                    None => es.push((v, c)),
                                }
                            }
                        }
                        for &(p, cp) in &es {
                            for &(q, cq) in &es {
                                let mut x = 0.0;
                                for k in 0..w {
                                    x += cp[k] * cq[k];
                                }
                                acc[p * nv + q] += wgt * x;
                            }
                        }
                    }
                }
                acc
            })
            .collect()
    });
    let mut out = zero_matrix(nv);
    for p in partial {
        for (o, x) in out.as_mut_slice().iter_mut().zip(p) {
            *o += x;
        }
    }
    Ok(out)
}

/// The fractional operator `L_τ` (Eq. 11 with order `2τ`, kernel `ρ^{−(2τ+2)}`), dense.
///
/// [formula: rs-sobolev-forms]
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// As [`fractional_kernel`].
pub fn fractional_operator(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    tau: f64,
) -> Result<DenseMatrix, RsError> {
    let h = fractional_kernel(pool, ops, tau)?;
    kernel_matrix_form(pool, ops, &h, Channel::Average)
}

/// The fractional operator over **vertex** pairs (`FractionalPairs::Vertices`, §6.3 option):
/// `uᵀ L u = Σ_{v≠w} (u_v − u_w)² a_v a_w |f_v − f_w|^{−(2τ+2)}`, i.e. the kernel-matrix form
/// `2 (diag(a ⊙ H a) − diag(a) H diag(a))` with `H_vw = |f_v − f_w|^{−(2τ+2)}` and the vertex
/// areas `a`. Parallel over rows; bitwise independent of the pool size.
///
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// [`RsError::InvalidOption`] for a non-finite `τ`, mismatched lengths or two vertices at the
/// same position.
pub fn vertex_fractional_operator(
    pool: &rayon::ThreadPool,
    positions: &[[f64; 3]],
    mass: &[f64],
    tau: f64,
) -> Result<DenseMatrix, RsError> {
    let n = positions.len();
    if !tau.is_finite() || mass.len() != n {
        return Err(RsError::InvalidOption {
            name: "vertex fractional operator",
            value: format!("tau {tau}, {} masses for {n} vertices", mass.len()),
            requirement: "a finite order and one mass per vertex",
        });
    }
    let e = -(tau + 1.0);
    for v in 0..n {
        for w in v + 1..n {
            let d = sub(positions[v], positions[w]);
            if !(dot(d, d) > 0.0) {
                return Err(RsError::InvalidOption {
                    name: "positions",
                    value: format!("vertices {v} and {w} coincide"),
                    requirement: "distinct vertex positions for the vertex-pair operator",
                });
            }
        }
    }
    let mut out = zero_matrix(n);
    pool.install(|| {
        out.as_mut_slice()
            .par_chunks_mut(n.max(1))
            .enumerate()
            .for_each(|(v, row)| {
                let mut diag = 0.0;
                for w in 0..n {
                    if w == v {
                        continue;
                    }
                    let (a, b) = if v < w { (v, w) } else { (w, v) };
                    let d = sub(positions[a], positions[b]);
                    let h = 2.0 * mass[v] * mass[w] * dot(d, d).powf(e);
                    row[w] = -h;
                    diag += h;
                }
                row[v] = diag;
            });
    });
    Ok(out)
}

/// The high-order operator `B` (Eq. 12, kernel `ρ^{−2s}`, through `D_f`), dense.
///
/// [formula: rs-sobolev-forms]
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// As [`high_order_kernel`].
pub fn high_order_operator(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    exponents: &Exponents,
) -> Result<DenseMatrix, RsError> {
    let h = high_order_kernel(pool, ops, exponents)?;
    kernel_matrix_form(pool, ops, &h, Channel::Gradient)
}

/// The low-order operator `B⁰` (Eq. 13, fixed `K₂`, symmetrized), dense.
///
/// [formula: rs-sobolev-forms]
/// [formula: rs-low-order-symmetrisation]
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// As [`low_order_kernel`].
pub fn low_order_operator(
    pool: &rayon::ThreadPool,
    ops: &FaceOperators,
    exponents: &Exponents,
) -> Result<DenseMatrix, RsError> {
    let h = low_order_kernel(pool, ops, exponents)?;
    kernel_matrix_form(pool, ops, &h, Channel::Average)
}

/// `W M⁻¹ W` (sparse, symmetric) from the stiffness `W` and the vertex areas `M = diag(m)`.
///
/// # Errors
/// [`RsError::InvalidOption`] for a non-positive mass or mismatched sizes.
pub fn willmore_matrix(w: &SparseMatrix, mass: &[f64]) -> Result<SparseMatrix, RsError> {
    let n = w.ncols();
    if w.nrows() != n || mass.len() != n {
        return Err(RsError::InvalidOption {
            name: "willmore_matrix",
            value: format!("{} x {} and {} masses", w.nrows(), n, mass.len()),
            requirement: "a square stiffness and one mass per vertex",
        });
    }
    if let Some(v) = mass.iter().position(|m| !(*m > 0.0 && m.is_finite())) {
        return Err(RsError::InvalidOption {
            name: "mass",
            value: format!("vertex {v}: {:e}", mass[v]),
            requirement: "positive finite vertex areas",
        });
    }
    // (W M⁻¹ W)_{ij} = Σ_k W_ik W_kj / m_k; W symmetric: column j holds W_kj, column k holds W_ik.
    let (cp, ri, val) = (w.col_ptr(), w.row_idx(), w.values());
    let mut trip: Vec<(u32, u32, f64)> = Vec::new();
    for j in 0..n {
        for p in cp[j]..cp[j + 1] {
            let k = ri[p] as usize;
            let wkj = val[p] / mass[k];
            for q in cp[k]..cp[k + 1] {
                trip.push((ri[q], j as u32, val[q] * wkj));
            }
        }
    }
    SparseMatrix::from_triplets(n, n, &trip).map_err(|e| RsError::InvalidOption {
        name: "willmore_matrix",
        value: e.to_string(),
        requirement: "the sparse product must fit",
    })
}

/// A sparse matrix as a dense one.
pub fn to_dense(m: &SparseMatrix) -> DenseMatrix {
    let n = m.nrows();
    let mut out = zero_matrix(n);
    let (cp, ri, val) = (m.col_ptr(), m.row_idx(), m.values());
    let data = out.as_mut_slice();
    for j in 0..m.ncols().min(n) {
        for p in cp[j]..cp[j + 1] {
            data[ri[p] as usize * n + j] += val[p];
        }
    }
    out
}

/// The dense metric of one configuration (`SPEC.md` §4.4, §9.1 `dense_metric`).
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMetric {
    /// `B` (Eq. 12).
    pub b: DenseMatrix,
    /// `B⁰` (Eq. 13).
    pub b0: DenseMatrix,
    /// The cotan stiffness `W`.
    pub w: SparseMatrix,
    /// The vertex areas `a_v` (diagonal of `M`).
    pub mass: Vec<f64>,
    /// `W M⁻¹ W` (only when `willmore_weight > 0`).
    pub wmw: Option<SparseMatrix>,
    /// `A = w_TP (B + c_lo B⁰) + w_W W M⁻¹ W + O_diag` (`|V| × |V|`; the metric on positions is
    /// `A ⊗ I₃`).
    pub a: DenseMatrix,
    /// The weights used.
    pub weights: SobolevWeights,
}

/// Assembles the dense Sobolev metric `A = w_TP (B + c_lo B⁰) + w_W W M⁻¹ W + O_diag` (§4.4).
/// `obstacle` is the hook for WP6/WP7's mesh-obstacle diagonal ([`FaceDiagonal`]). `B⁰` is
/// assembled (and reported) also when `low_order` is off.
///
/// [formula: rs-sobolev-forms]
/// [formula: rs-kernel-matrix-form]
///
/// # Errors
/// The errors of [`Surface::validate_positions`]; [`RsError::NonFiniteEnergy`] for coincident
/// barycenters; [`RsError::InvalidOption`] for invalid weights or an obstacle diagonal of the
/// wrong length.
pub fn dense_metric(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    exponents: &Exponents,
    weights: &SobolevWeights,
    obstacle: Option<&FaceDiagonal>,
) -> Result<DenseMetric, RsError> {
    weights.validate()?;
    let geom = FaceGeometry::new(surface, positions)?;
    let mass = vertex_areas(surface, &geom.areas);
    let ops = FaceOperators::from_geometry(surface, positions, geom);
    let b = high_order_operator(pool, &ops, exponents)?;
    let b0 = low_order_operator(pool, &ops, exponents)?;
    let w = cotan_stiffness(surface, positions)?;
    let wmw = if weights.willmore_weight > 0.0 {
        Some(willmore_matrix(&w, &mass)?)
    } else {
        None
    };
    let nv = surface.n_vertices();
    let mut a = zero_matrix(nv);
    {
        let c_lo = weights.c_lo();
        let wt = weights.tp_weight;
        let (bs, b0s) = (b.as_slice(), b0.as_slice());
        for (i, x) in a.as_mut_slice().iter_mut().enumerate() {
            *x = wt * (bs[i] + c_lo * b0s[i]);
        }
    }
    if let Some(m) = &wmw {
        let ww = weights.willmore_weight;
        let (cp, ri, val) = (m.col_ptr(), m.row_idx(), m.values());
        let data = a.as_mut_slice();
        for j in 0..nv {
            for p in cp[j]..cp[j + 1] {
                data[ri[p] as usize * nv + j] += ww * val[p];
            }
        }
    }
    if let Some(o) = obstacle {
        add_face_diagonal(&ops, o, &mut a)?;
    }
    Ok(DenseMetric {
        b,
        b0,
        w,
        mass,
        wmw,
        a,
        weights: *weights,
    })
}

/// Adds `2 Vᵀ diag(hi/a) V + 2 Uᵀ diag(lo/a) U` to `a` (the obstacle hook, §4.4).
///
/// # Errors
/// [`RsError::InvalidOption`] if a non-empty `hi`/`lo` does not have one entry per face.
pub fn add_face_diagonal(
    ops: &FaceOperators,
    o: &FaceDiagonal,
    a: &mut DenseMatrix,
) -> Result<(), RsError> {
    let nf = ops.n_faces();
    let nv = ops.n_vertices;
    for (name, d, channel) in [
        ("obstacle diagonal hi", &o.hi, Channel::Gradient),
        ("obstacle diagonal lo", &o.lo, Channel::Average),
    ] {
        if d.is_empty() {
            continue;
        }
        if d.len() != nf {
            return Err(RsError::InvalidOption {
                name,
                value: format!("{} entries", d.len()),
                requirement: "one entry per face",
            });
        }
        let data = a.as_mut_slice();
        for s in 0..nf {
            let r = ops.row(s, channel);
            let f = 2.0 * d[s] / ops.areas[s];
            for &(p, cp) in &r {
                for &(q, cq) in &r {
                    let mut x = 0.0;
                    for k in 0..channel.width() {
                        x += cp[k] * cq[k];
                    }
                    data[p * nv + q] += f * x;
                }
            }
        }
    }
    Ok(())
}

/// The sparse matrix of a comparison metric (`SPEC.md` §4.5, D25; `|V| × |V|`, the metric on
/// positions is `⊗ I₃`): `L2(Mass)` = `M`, `L2(Identity)` = `I`, `H1` = `W`, `H2` = `W M⁻¹ W`.
///
/// # Errors
/// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] for
/// [`Metric::Sobolev`] (use [`dense_metric`] or the hierarchical back end).
pub fn comparison_metric(
    surface: &Surface,
    positions: &[[f64; 3]],
    metric: Metric,
) -> Result<SparseMatrix, RsError> {
    let geom = FaceGeometry::new(surface, positions)?;
    let n = surface.n_vertices();
    let diag = |d: &[f64]| {
        SparseMatrix::from_diagonal(d).map_err(|e| RsError::InvalidOption {
            name: "comparison metric",
            value: e.to_string(),
            requirement: "a representable diagonal",
        })
    };
    match metric {
        Metric::Sobolev { .. } => Err(RsError::InvalidOption {
            name: "metric",
            value: "Sobolev".into(),
            requirement: "a comparison metric (L2, H1, H2); the Sobolev metric is dense_metric",
        }),
        Metric::L2(L2Kind::Mass) => diag(&vertex_areas(surface, &geom.areas)),
        Metric::L2(L2Kind::Identity) => diag(&vec![1.0; n]),
        Metric::H1 => cotan_stiffness(surface, positions),
        Metric::H2 => {
            let w = cotan_stiffness(surface, positions)?;
            willmore_matrix(&w, &vertex_areas(surface, &geom.areas))
        }
    }
}
