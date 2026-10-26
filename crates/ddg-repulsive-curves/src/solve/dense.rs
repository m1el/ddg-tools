//! Dense saddle-point solver (`SPEC.md` §6.1–6.2, §8.5, §17l).
//!
//! The saddle matrix of one step,
//!
//! ```text
//! K = [ Ā  Cᵀ ]      Ā = A ⊗ I_d (interleaved, D29),
//!     [ C  0  ]
//! ```
//!
//! is factored once and reused for Eq. 22 and every chord iteration of the projection (Sec. 7.4),
//! by one of two methods ([`DenseMethod`]):
//!
//! * **Schur** (the default, `SPEC.md` §17l, the idea of §17b): `Ā = A ⊗ I_d` needs only the
//!   `|V| × |V|` matrix `A`. `A` is positive semidefinite with `ker A = span{𝟙}` (§4.5), so
//!   `M = A + ρ u uᵀ` with `u = 𝟙/√|V|`, `ρ = max_v A_vv`, is positive definite and
//!   `M u = ρ u`. With `U = u ⊗ I_d` and `μ = −ρ Uᵀ x`, `K [x; λ] = [a; b]` is
//!
//!   ```text
//!   M̄ x + Cᵀ λ + U μ = a,   C x = b,   Uᵀ x + μ/ρ = 0      (M̄ = M ⊗ I_d),
//!   ```
//!
//!   and eliminating `x = M̄⁻¹ (a − Cᵀ λ − U μ)` (with `M̄⁻¹ U = U/ρ`) leaves
//!
//!   ```text
//!   [ S   E ] [ λ ]   [ C M̄⁻¹ a − b ]
//!   [ Eᵀ  0 ] [ μ ] = [ Uᵀ a / ρ     ],    S = C M̄⁻¹ Cᵀ,  E = C U / ρ.
//!   ```
//!
//!   `S` (`k × k`) is positive definite when `C` has full row rank, and the `d × d` matrix
//!   `W = Eᵀ S⁻¹ E` is positive definite when the constraints fix the translations (D5), which is
//!   exactly when `K` is nonsingular for a positive definite `Ā` on `ker C`. Work: Cholesky of `M`
//!   (`n³/3`), its explicit inverse (`2n³/3`), `S` by sparse products (`O(n · nnz C)`), Cholesky
//!   of `S` (`k³/3`); about `1.3 n³` flops for the default constraints (`k ≈ n`) against
//!   `(2/3)(3n + k)³ ≈ 43 n³` for the LU. One step of iterative refinement on the full system
//!   follows every solve.
//! * **LU**: one LU with partial pivoting of the whole `(d|V| + k)` matrix, scaled to
//!   `[Ā/s Cᵀ; C 0]`, `s = max_v A_vv` (the scaling only balances the blocks for pivoting; `λ`
//!   is scaled back). `O((d|V| + k)³)`.
//!
//! The Schur method falls back to the LU when one of its factorizations fails or is close to
//! singular ([`SCHUR_PIVOT_TOLERANCE`]), so singular and nearly singular saddles are decided by
//! the LU as before. Both run on an explicit pool (`faer`); results are bitwise independent of
//! its size (D25).

use std::sync::Arc;
use std::time::Instant;

use ddg_core::sparse::SparseMatrix;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::llt::factor::{
    LltRegularization, cholesky_in_place, cholesky_in_place_scratch,
};
use faer::linalg::cholesky::llt::inverse::{
    inverse as llt_inverse, inverse_scratch as llt_inverse_scratch,
};
use faer::linalg::cholesky::llt::solve::solve_in_place as llt_solve;
use faer::linalg::lu::partial_pivoting::{factor as lu_factor, solve as lu_solve};
use faer::perm::PermRef;
use faer::{Mat, Par};
use rayon::prelude::*;

use super::{LinearReport, SaddleSolution, SaddleSolver, check_rhs};
use crate::error::RcError;
use crate::metric::{DenseMatrix, DofLayout};

/// Pivots with `|u_ii| ≤ PIVOT_TOLERANCE · max |u_jj|` mean the saddle matrix is singular.
pub const PIVOT_TOLERANCE: f64 = 1e-14;

/// The factorization of a [`DenseSaddle`] (`SPEC.md` §17l).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DenseMethod {
    /// Cholesky of `A + ρ u uᵀ` and the `k × k` Schur complement (the default), falling back to
    /// [`DenseMethod::Lu`] when it fails.
    Schur,
    /// LU with partial pivoting of the whole saddle matrix.
    Lu,
}

/// The Schur method falls back to the LU when a Cholesky factor of `M`, `S` or `W` has
/// `min L_ii² / max L_ii²` at or below this (nearly singular; the LU then decides, and reports a
/// singular saddle by [`PIVOT_TOLERANCE`] as before).
pub const SCHUR_PIVOT_TOLERANCE: f64 = 1e-13;

#[derive(Debug)]
struct LuFactor {
    scale: f64,
    lu: Mat<f64>,
    perm: Vec<usize>,
    perm_inv: Vec<usize>,
}

#[derive(Debug)]
struct SchurFactor {
    /// `ρ`.
    rho: f64,
    /// `M⁻¹` (`n × n`, symmetric, both triangles).
    minv: Mat<f64>,
    /// Cholesky factor of `S = C M̄⁻¹ Cᵀ` (lower).
    ls: Mat<f64>,
    /// `E = C U / ρ` (`k × d`).
    e: Mat<f64>,
    /// `Z = S⁻¹ E` (`k × d`).
    z: Mat<f64>,
    /// Cholesky factor of `W = Eᵀ Z` (`d × d`, lower).
    lw: Mat<f64>,
}

#[derive(Debug)]
enum Factor {
    Lu(LuFactor),
    Schur(SchurFactor),
}

/// The factored dense saddle matrix of one step.
#[derive(Debug)]
pub struct DenseSaddle {
    pool: Arc<rayon::ThreadPool>,
    layout: DofLayout,
    n_vertices: usize,
    k: usize,
    factor: Factor,
    /// Copies for the residual check.
    a: DenseMatrix,
    c: SparseMatrix,
    ct: SparseMatrix,
    pivot_ratio: f64,
    factor_seconds: f64,
}

fn par_of(pool: &rayon::ThreadPool) -> Par {
    let n = pool.current_num_threads();
    if n > 1 { Par::rayon(n) } else { Par::Seq }
}

/// `min L_ii² / max L_ii²` of a Cholesky factor (1 for an empty one).
fn cholesky_ratio(l: &Mat<f64>) -> f64 {
    let n = l.nrows();
    if n == 0 {
        return 1.0;
    }
    let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
    for i in 0..n {
        let v = l[(i, i)] * l[(i, i)];
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if hi > 0.0 && lo.is_finite() {
        lo / hi
    } else {
        0.0
    }
}

/// Cholesky in place on `pool`; the strict upper triangle is zeroed. `None` on a non-positive
/// pivot.
fn cholesky(pool: &rayon::ThreadPool, m: &mut Mat<f64>) -> Option<()> {
    let n = m.nrows();
    let par = par_of(pool);
    let ok = pool.install(|| {
        let mut mem = MemBuffer::new(cholesky_in_place_scratch::<f64>(n, par, Default::default()));
        cholesky_in_place(
            m.as_mut(),
            LltRegularization::default(),
            par,
            MemStack::new(&mut mem),
            Default::default(),
        )
        .is_ok()
    });
    if !ok {
        return None;
    }
    for j in 1..n {
        for i in 0..j {
            m[(i, j)] = 0.0;
        }
    }
    Some(())
}

impl DenseSaddle {
    /// Builds and factors the saddle matrix of the metric `a` (`|V| × |V|`, symmetric) and the
    /// constraint Jacobian `c` (`k × d|V|`) on `pool`, by the Schur method with the LU as its
    /// fallback ([`DenseMethod::Schur`]).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if the dimensions do not match; [`RcError::SingularMetric`]
    /// (`components: 0`) if `a` has no positive diagonal entry; [`RcError::SingularSaddle`] with
    /// the index of the smallest pivot if a pivot of the LU is below [`PIVOT_TOLERANCE`]
    /// relative to the largest (or not finite).
    pub fn new(
        pool: Arc<rayon::ThreadPool>,
        a: &DenseMatrix,
        c: &SparseMatrix,
        layout: DofLayout,
    ) -> Result<DenseSaddle, RcError> {
        Self::with_method(pool, a, c, layout, DenseMethod::Schur)
    }

    /// [`DenseSaddle::new`] with an explicit method ([`DenseMethod::Lu`]: the LU only).
    ///
    /// # Errors
    /// As [`DenseSaddle::new`].
    pub fn with_method(
        pool: Arc<rayon::ThreadPool>,
        a: &DenseMatrix,
        c: &SparseMatrix,
        layout: DofLayout,
        method: DenseMethod,
    ) -> Result<DenseSaddle, RcError> {
        let t0 = Instant::now();
        let n = a.n();
        let d = layout.dim();
        let nd = n * d;
        if c.ncols() != nd {
            return Err(RcError::InvalidOption {
                name: "constraint Jacobian",
                value: format!("{} columns", c.ncols()),
                requirement: "d|V| columns, matching the metric",
            });
        }
        let k = c.nrows();
        let scale = (0..n).map(|v| a.get(v, v)).fold(0.0f64, f64::max);
        if !(scale > 0.0 && scale.is_finite()) {
            return Err(RcError::SingularMetric { components: 0 });
        }
        let ct = c.transpose();
        let schur = match method {
            DenseMethod::Schur => factor_schur(&pool, a, c, &ct, layout, scale),
            DenseMethod::Lu => None,
        };
        let (factor, pivot_ratio) = match schur {
            Some((f, r)) => (Factor::Schur(f), r),
            None => {
                let (f, r) = factor_lu(&pool, a, c, &ct, layout, scale)?;
                (Factor::Lu(f), r)
            }
        };
        Ok(DenseSaddle {
            pool,
            layout,
            n_vertices: n,
            k,
            factor,
            a: a.clone(),
            c: c.clone(),
            ct,
            pivot_ratio,
            factor_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// The method in use (`Lu` also when the Schur method fell back to it).
    pub fn method(&self) -> DenseMethod {
        match self.factor {
            Factor::Lu(_) => DenseMethod::Lu,
            Factor::Schur(_) => DenseMethod::Schur,
        }
    }

    /// Seconds spent assembling and factoring.
    pub fn factor_seconds(&self) -> f64 {
        self.factor_seconds
    }

    /// LU: `min |u_ii| / max |u_ii|` of the LU factor; Schur: the smallest
    /// `min L_ii² / max L_ii²` of its Cholesky factors (of `M`, `S`, `W`).
    pub fn pivot_ratio(&self) -> f64 {
        self.pivot_ratio
    }

    /// The metric `A` the solver was built with.
    pub fn metric(&self) -> &DenseMatrix {
        &self.a
    }

    /// The constraint Jacobian `C` the solver was built with.
    pub fn jacobian(&self) -> &SparseMatrix {
        &self.c
    }

    /// The layout of the unknowns.
    pub fn layout(&self) -> DofLayout {
        self.layout
    }

    /// `Ā x` (`A ⊗ I_d` times an interleaved vector), parallel over vertices, fixed order.
    pub fn apply_metric(&self, x: &[f64]) -> Vec<f64> {
        let d = self.layout.dim();
        let n = self.n_vertices;
        let mut y = vec![0.0; n * d];
        self.pool.install(|| {
            y.par_chunks_mut(d).enumerate().for_each(|(v, yv)| {
                let row = self.a.row(v);
                for (w, &aw) in row.iter().enumerate() {
                    for s in 0..d {
                        yv[s] += aw * x[w * d + s];
                    }
                }
            });
        });
        y
    }

    /// `(M⁻¹ ⊗ I_d) x`, parallel over vertices, fixed order.
    fn apply_minv(&self, f: &SchurFactor, x: &[f64]) -> Vec<f64> {
        let d = self.layout.dim();
        let n = self.n_vertices;
        let mut y = vec![0.0; n * d];
        self.pool.install(|| {
            y.par_chunks_mut(d).enumerate().for_each(|(v, yv)| {
                // M⁻¹ is symmetric: row v = column v (contiguous).
                let col = f.minv.col_as_slice(v);
                for (w, &mw) in col.iter().enumerate() {
                    for s in 0..d {
                        yv[s] += mw * x[w * d + s];
                    }
                }
            });
        });
        y
    }

    /// One solve of the Schur system (module docs), without refinement.
    fn solve_schur(
        &self,
        f: &SchurFactor,
        a: &[f64],
        b: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), RcError> {
        let d = self.layout.dim();
        let n = self.n_vertices;
        let k = self.k;
        let su = 1.0 / (f.rho * (n as f64).sqrt());
        // r_C = C M̄⁻¹ a − b, r_U = Uᵀ a / ρ.
        let ga = self.apply_minv(f, a);
        let mut rc = vec![0.0; k];
        self.c
            .mul_vec(&ga, &mut rc)
            .map_err(|e| super::sparse_error(e, "C"))?;
        for (r, bi) in rc.iter_mut().zip(b) {
            *r -= bi;
        }
        let mut ru = vec![0.0; d];
        for v in 0..n {
            for s in 0..d {
                ru[s] += a[v * d + s];
            }
        }
        for r in &mut ru {
            *r *= su;
        }
        // μ = W⁻¹ (Zᵀ r_C − r_U).
        let mut mu = Mat::<f64>::from_fn(d, 1, |s, _| {
            let mut t = -ru[s];
            for i in 0..k {
                t += f.z[(i, s)] * rc[i];
            }
            t
        });
        llt_solve(f.lw.as_ref(), mu.as_mut(), Par::Seq, MemStack::new(&mut []));
        // λ = S⁻¹ (r_C − E μ).
        let mut lam = Mat::<f64>::from_fn(k, 1, |i, _| {
            let mut t = rc[i];
            for s in 0..d {
                t -= f.e[(i, s)] * mu[(s, 0)];
            }
            t
        });
        // Sequential triangular solves: deterministic (D25), memory bound.
        llt_solve(
            f.ls.as_ref(),
            lam.as_mut(),
            Par::Seq,
            MemStack::new(&mut []),
        );
        let lambda: Vec<f64> = (0..k).map(|i| lam[(i, 0)]).collect();
        // x = M̄⁻¹ (a − Cᵀ λ) − U μ / ρ.
        let mut ctl = vec![0.0; n * d];
        self.ct
            .mul_vec(&lambda, &mut ctl)
            .map_err(|e| super::sparse_error(e, "C^T"))?;
        let t: Vec<f64> = a.iter().zip(&ctl).map(|(p, q)| p - q).collect();
        let mut x = self.apply_minv(f, &t);
        for v in 0..n {
            for s in 0..d {
                x[v * d + s] -= mu[(s, 0)] * su;
            }
        }
        Ok((x, lambda))
    }

    /// The residual `(a − Ā x − Cᵀ λ, b − C x)` of the unscaled system.
    fn residual(
        &self,
        a: &[f64],
        b: &[f64],
        x: &[f64],
        lambda: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), RcError> {
        let nd = self.n_dofs();
        let ax = self.apply_metric(x);
        let mut ctl = vec![0.0; nd];
        self.ct
            .mul_vec(lambda, &mut ctl)
            .map_err(|e| super::sparse_error(e, "C^T"))?;
        let mut cx = vec![0.0; self.k];
        self.c
            .mul_vec(x, &mut cx)
            .map_err(|e| super::sparse_error(e, "C"))?;
        let r1 = (0..nd).map(|i| a[i] - ax[i] - ctl[i]).collect();
        let r2 = (0..self.k).map(|i| b[i] - cx[i]).collect();
        Ok((r1, r2))
    }
}

/// The LU of the scaled saddle matrix.
fn factor_lu(
    pool: &Arc<rayon::ThreadPool>,
    a: &DenseMatrix,
    c: &SparseMatrix,
    ct: &SparseMatrix,
    layout: DofLayout,
    scale: f64,
) -> Result<(LuFactor, f64), RcError> {
    let n = a.n();
    let d = layout.dim();
    let nd = n * d;
    let k = c.nrows();
    let dim = nd + k;
    let inv = 1.0 / scale;
    let mut lu = Mat::<f64>::zeros(dim, dim);
    let par = par_of(pool);
    let (ccp, cri, cval) = (c.col_ptr(), c.row_idx(), c.values());
    let (tcp, tri, tval) = (ct.col_ptr(), ct.row_idx(), ct.values());
    pool.install(|| {
        lu.par_col_iter_mut().enumerate().for_each(|(j, mut col)| {
            if j < nd {
                let (v, s) = (j / d, j % d);
                let row = a.row(v); // A is symmetric: column v = row v
                for w in 0..n {
                    col[w * d + s] = row[w] * inv;
                }
                for p in ccp[j]..ccp[j + 1] {
                    col[nd + cri[p] as usize] = cval[p];
                }
            } else {
                let r = j - nd;
                for p in tcp[r]..tcp[r + 1] {
                    col[tri[p] as usize] = tval[p];
                }
            }
        });
    });
    let mut perm = vec![0usize; dim];
    let mut perm_inv = vec![0usize; dim];
    pool.install(|| {
        let mut mem = MemBuffer::new(lu_factor::lu_in_place_scratch::<usize, f64>(
            dim,
            dim,
            par,
            Default::default(),
        ));
        lu_factor::lu_in_place(
            lu.as_mut(),
            &mut perm,
            &mut perm_inv,
            par,
            MemStack::new(&mut mem),
            Default::default(),
        );
    });
    let mut max_p: f64 = 0.0;
    let mut min_p = f64::INFINITY;
    let mut argmin = 0;
    for i in 0..dim {
        let p = lu[(i, i)].abs();
        if !p.is_finite() {
            return Err(RcError::SingularSaddle { pivot: i });
        }
        max_p = max_p.max(p);
        if p < min_p {
            min_p = p;
            argmin = i;
        }
    }
    let pivot_ratio = if dim == 0 { 1.0 } else { min_p / max_p };
    if pivot_ratio.is_nan() || pivot_ratio <= PIVOT_TOLERANCE {
        return Err(RcError::SingularSaddle { pivot: argmin });
    }
    Ok((
        LuFactor {
            scale,
            lu,
            perm,
            perm_inv,
        },
        pivot_ratio,
    ))
}

/// The Schur factorization (module docs); `None` when a factor fails or is nearly singular.
fn factor_schur(
    pool: &Arc<rayon::ThreadPool>,
    a: &DenseMatrix,
    c: &SparseMatrix,
    ct: &SparseMatrix,
    layout: DofLayout,
    rho: f64,
) -> Option<(SchurFactor, f64)> {
    let n = a.n();
    let d = layout.dim();
    let k = c.nrows();
    if n == 0 || k == 0 {
        // No constraint cannot fix the translations: the saddle is singular (the LU says so).
        return None;
    }
    let par = par_of(pool);
    // M = A + (ρ/n) 𝟙𝟙ᵀ.
    let shift = rho / n as f64;
    let mut m = Mat::<f64>::zeros(n, n);
    pool.install(|| {
        m.par_col_iter_mut().enumerate().for_each(|(j, mut col)| {
            let row = a.row(j);
            for i in 0..n {
                col[i] = row[i] + shift;
            }
        });
    });
    cholesky(pool, &mut m)?;
    let ratio_m = cholesky_ratio(&m);
    if ratio_m.is_nan() || ratio_m <= SCHUR_PIVOT_TOLERANCE {
        return None;
    }
    // M⁻¹ (lower triangle), then mirrored.
    let mut minv = Mat::<f64>::zeros(n, n);
    pool.install(|| {
        let mut mem = MemBuffer::new(llt_inverse_scratch::<f64>(n, par));
        llt_inverse(minv.as_mut(), m.as_ref(), par, MemStack::new(&mut mem));
    });
    drop(m);
    for j in 0..n {
        for i in 0..j {
            minv[(i, j)] = minv[(j, i)];
        }
    }
    if (0..n).any(|i| !minv[(i, i)].is_finite()) {
        return None;
    }
    // S = C M̄⁻¹ Cᵀ: column i of S is C (M̄⁻¹ c_iᵀ), c_i the i-th row of C (= column i of Cᵀ).
    let (tcp, tri, tval) = (ct.col_ptr(), ct.row_idx(), ct.values());
    let mut sm = Mat::<f64>::zeros(k, k);
    pool.install(|| {
        sm.par_col_iter_mut().enumerate().for_each_init(
            || (vec![0.0; n * d], vec![0.0; k]),
            |(g, out), (i, mut col)| {
                g.iter_mut().for_each(|x| *x = 0.0);
                for p in tcp[i]..tcp[i + 1] {
                    let (j, val) = (tri[p] as usize, tval[p]);
                    let (v, s) = (j / d, j % d);
                    let mcol = minv.col_as_slice(v);
                    for (w, &mw) in mcol.iter().enumerate() {
                        g[w * d + s] += val * mw;
                    }
                }
                // out = C g, row by row through Cᵀ's columns (fixed order).
                for (r, o) in out.iter_mut().enumerate() {
                    let mut t = 0.0;
                    for p in tcp[r]..tcp[r + 1] {
                        t += tval[p] * g[tri[p] as usize];
                    }
                    *o = t;
                }
                for r in 0..k {
                    col[r] = out[r];
                }
            },
        );
    });
    cholesky(pool, &mut sm)?;
    let ratio_s = cholesky_ratio(&sm);
    if ratio_s.is_nan() || ratio_s <= SCHUR_PIVOT_TOLERANCE {
        return None;
    }
    // E = C U / ρ: E[i, s] = Σ_v C[i, v d + s] / (ρ √n).
    let su = 1.0 / (rho * (n as f64).sqrt());
    let mut e = Mat::<f64>::zeros(k, d);
    for i in 0..k {
        for p in tcp[i]..tcp[i + 1] {
            let j = tri[p] as usize;
            e[(i, j % d)] += tval[p] * su;
        }
    }
    let mut z = e.clone();
    llt_solve(sm.as_ref(), z.as_mut(), Par::Seq, MemStack::new(&mut []));
    let mut lw = Mat::<f64>::from_fn(d, d, |s, t| {
        let mut acc = 0.0;
        for i in 0..k {
            acc += e[(i, s)] * z[(i, t)];
        }
        acc
    });
    // Symmetrize W exactly before its Cholesky.
    for s in 0..d {
        for t in 0..s {
            let v = 0.5 * (lw[(s, t)] + lw[(t, s)]);
            lw[(s, t)] = v;
            lw[(t, s)] = v;
        }
    }
    let seq = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .ok()?;
    cholesky(&seq, &mut lw)?;
    let ratio_w = cholesky_ratio(&lw);
    if ratio_w.is_nan() || ratio_w <= SCHUR_PIVOT_TOLERANCE {
        return None;
    }
    let ratio = ratio_m.min(ratio_s).min(ratio_w);
    Some((
        SchurFactor {
            rho,
            minv,
            ls: sm,
            e,
            z,
            lw,
        },
        ratio,
    ))
}

/// Steps of iterative refinement after every Schur solve (on the full unscaled system).
pub const SCHUR_REFINEMENT_STEPS: usize = 1;

impl DenseSaddle {
    fn solve_lu(&self, f: &LuFactor, a: &[f64], b: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let nd = self.n_dofs();
        let dim = nd + self.k;
        let inv = 1.0 / f.scale;
        let mut rhs =
            Mat::<f64>::from_fn(dim, 1, |i, _| if i < nd { a[i] * inv } else { b[i - nd] });
        // The triangular solves of one right-hand side run sequentially: they are memory bound
        // (no gain from threads) and faer's parallel triangular solve sums in an order that
        // depends on the thread count, which would break D25.
        let par = Par::Seq;
        let perm = PermRef::new_checked(&f.perm, &f.perm_inv, dim);
        let mut mem = MemBuffer::new(lu_solve::solve_in_place_scratch::<usize, f64>(dim, 1, par));
        lu_solve::solve_in_place(
            f.lu.as_ref(),
            f.lu.as_ref(),
            perm,
            rhs.as_mut(),
            par,
            MemStack::new(&mut mem),
        );
        let x: Vec<f64> = (0..nd).map(|i| rhs[(i, 0)]).collect();
        let lambda: Vec<f64> = (0..self.k).map(|r| rhs[(nd + r, 0)] * f.scale).collect();
        (x, lambda)
    }
}

impl SaddleSolver for DenseSaddle {
    fn n_dofs(&self) -> usize {
        self.n_vertices * self.layout.dim()
    }

    fn n_constraints(&self) -> usize {
        self.k
    }

    /// [formula: rc-saddle-gradient]
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        check_rhs(self, a, b)?;
        let t0 = Instant::now();
        let nd = self.n_dofs();
        let dim = nd + self.k;
        let (x, lambda, residual, name, iterations) = match &self.factor {
            Factor::Lu(f) => {
                let (x, lambda) = self.solve_lu(f, a, b);
                let r = self.residual(a, b, &x, &lambda)?;
                (x, lambda, r, "dense LU", 0)
            }
            Factor::Schur(f) => {
                let (mut x, mut lambda) = self.solve_schur(f, a, b)?;
                let mut r = self.residual(a, b, &x, &lambda)?;
                for _ in 0..SCHUR_REFINEMENT_STEPS {
                    let (dx, dl) = self.solve_schur(f, &r.0, &r.1)?;
                    for (p, q) in x.iter_mut().zip(&dx) {
                        *p += q;
                    }
                    for (p, q) in lambda.iter_mut().zip(&dl) {
                        *p += q;
                    }
                    r = self.residual(a, b, &x, &lambda)?;
                }
                (x, lambda, r, "dense Schur", 0)
            }
        };
        if let Some(i) = x.iter().chain(&lambda).position(|v| !v.is_finite()) {
            return Err(RcError::SingularSaddle { pivot: i });
        }
        let r2: f64 = residual.0.iter().chain(&residual.1).map(|r| r * r).sum();
        let rhs2: f64 = a.iter().chain(b).map(|v| v * v).sum();
        let relative_residual = if rhs2 > 0.0 {
            (r2 / rhs2).sqrt()
        } else {
            r2.sqrt()
        };
        Ok(SaddleSolution {
            x,
            lambda,
            report: LinearReport {
                solver: name,
                dimension: dim,
                iterations,
                relative_residual,
                fallback_used: false,
                levels: 0,
                pivot_ratio: Some(self.pivot_ratio),
                seconds: t0.elapsed().as_secs_f64(),
            },
        })
    }
}
