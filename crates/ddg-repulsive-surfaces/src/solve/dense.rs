//! Dense nested saddle solver and the dense preconditioner (`SPEC.md` §6.3, `Backend::Dense`).
//!
//! [`DenseSaddle`] factors `A′ = [[A ⊗ I₃, C_Qᵀ], [C_Q, 0]]` by LU with partial pivoting of the
//! equilibrated matrix (primal unknowns scaled by `1/√A_vv`, constraint rows to unit norm, as the
//! oracle), with one step of iterative refinement per solve. The positional rows of §5.1 act on
//! the three axes alike (`C_Q = C₁ ⊗ I₃` in groups of three rows); then `A′` is three copies of
//! the scalar saddle `[[A, C₁ᵀ], [C₁, 0]]` and one LU of size `|V| + k_Q/3` serves all three
//! axes (27× fewer flops). Any other `C_Q` is factored as a whole. The pivot ratio of the
//! equilibrated LU below [`PIVOT_TOLERANCE`] is [`RsError::SingularSaddle`].
//!
//! [`DensePreconditioner`] is the dense `P⁻¹ = L′⁻¹ blockdiag(L_{2−s} ⊗ I₃, 0) L′⁻¹` of §6.3 for
//! tests and the oracle (WP6 builds the sparse one); [`DenseSaddleOperator`] multiplies by `A′`
//! for GMRES on dense problems.

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::time::Instant;

use ddg_core::SparseMatrix;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::lu::partial_pivoting::{factor as lu_factor, solve as lu_solve};
use faer::perm::PermRef;
use faer::{Mat, Par};

use super::gmres::{LinearOperator, Preconditioner};
use super::{LinearReport, SaddleSolution, SaddleSolver, SolveRecord, check_rhs, norm2};
use crate::error::RsError;
use crate::geometry::{FaceGeometry, cotan_stiffness, vertex_areas};
use crate::kernel::Exponents;
use crate::metric::dense::{fractional_operator, to_dense, vertex_fractional_operator};
use crate::metric::{DenseMatrix, FaceOperators, FractionalPairs};
use crate::surface::Surface;

/// The equilibrated LU is singular when `min |u_ii| / max |u_ii|` is at most this (the oracle's
/// `SingularSaddle` threshold).
pub const PIVOT_TOLERANCE: f64 = 1e-13;

fn par_of(pool: &rayon::ThreadPool) -> Par {
    let n = pool.current_num_threads();
    if n > 1 { Par::rayon(n) } else { Par::Seq }
}

/// An LU factorization with partial pivoting of an equilibrated square matrix `D K D`.
#[derive(Debug)]
struct Lu {
    lu: Mat<f64>,
    perm: Vec<usize>,
    perm_inv: Vec<usize>,
    /// The equilibration `D`.
    d: Vec<f64>,
}

impl Lu {
    /// Factors `D K D` (`k` dense, column-major, consumed); returns the factor and the pivot ratio.
    fn new(pool: &rayon::ThreadPool, mut k: Mat<f64>, d: Vec<f64>) -> Result<(Lu, f64), RsError> {
        let n = k.nrows();
        for j in 0..n {
            for i in 0..n {
                k[(i, j)] *= d[i] * d[j];
            }
        }
        let mut perm = vec![0usize; n];
        let mut perm_inv = vec![0usize; n];
        let par = par_of(pool);
        pool.install(|| {
            let mut mem = MemBuffer::new(lu_factor::lu_in_place_scratch::<usize, f64>(
                n,
                n,
                par,
                Default::default(),
            ));
            lu_factor::lu_in_place(
                k.as_mut(),
                &mut perm,
                &mut perm_inv,
                par,
                MemStack::new(&mut mem),
                Default::default(),
            );
        });
        let (mut lo, mut hi, mut arg) = (f64::INFINITY, 0.0f64, 0usize);
        for i in 0..n {
            let p = k[(i, i)].abs();
            if !p.is_finite() {
                return Err(RsError::SingularSaddle { pivot: i });
            }
            hi = hi.max(p);
            if p < lo {
                lo = p;
                arg = i;
            }
        }
        let ratio = if n == 0 { 1.0 } else { lo / hi };
        if !(ratio > PIVOT_TOLERANCE) {
            return Err(RsError::SingularSaddle { pivot: arg });
        }
        Ok((
            Lu {
                lu: k,
                perm,
                perm_inv,
                d,
            },
            ratio,
        ))
    }

    /// Solves `K X = R` for the columns of `r` in place (sequential: bitwise reproducible).
    fn solve(&self, r: &mut Mat<f64>) {
        let n = self.lu.nrows();
        let m = r.ncols();
        for j in 0..m {
            for i in 0..n {
                r[(i, j)] *= self.d[i];
            }
        }
        let par = Par::Seq;
        let perm = PermRef::new_checked(&self.perm, &self.perm_inv, n);
        let mut mem = MemBuffer::new(lu_solve::solve_in_place_scratch::<usize, f64>(n, m, par));
        lu_solve::solve_in_place(
            self.lu.as_ref(),
            self.lu.as_ref(),
            perm,
            r.as_mut(),
            par,
            MemStack::new(&mut mem),
        );
        for j in 0..m {
            for i in 0..n {
                r[(i, j)] *= self.d[i];
            }
        }
    }
}

/// `C_Q = C₁ ⊗ I₃` in groups of three rows? Returns `C₁` (`k_Q/3 × |V|`) if so.
fn per_axis_rows(c: &SparseMatrix, nv: usize) -> Option<SparseMatrix> {
    let k = c.nrows();
    if !k.is_multiple_of(3) || c.ncols() != 3 * nv {
        return None;
    }
    let ct = c.transpose(); // column r = row r of C
    let (cp, ri, val) = (ct.col_ptr(), ct.row_idx(), ct.values());
    let mut trip = Vec::new();
    for g in 0..k / 3 {
        let base = cp[3 * g]..cp[3 * g + 1];
        for a in 0..3 {
            let r = 3 * g + a;
            let range = cp[r]..cp[r + 1];
            if range.len() != base.len() {
                return None;
            }
            for (p, q) in range.zip(base.clone()) {
                let (col, colb) = (ri[p] as usize, ri[q] as usize);
                if col % 3 != a || colb % 3 != 0 || col / 3 != colb / 3 || val[p] != val[q] {
                    return None;
                }
            }
        }
        for q in base {
            trip.push((g as u32, ri[q] / 3, val[q]));
        }
    }
    SparseMatrix::from_triplets(k / 3, nv, &trip).ok()
}

/// `y = A′ [x; λ]` for `A′ = [[A ⊗ I₃, C_Qᵀ], [C_Q, 0]]` (dense `A`, sparse `C_Q`).
#[derive(Debug, Clone)]
pub struct DenseSaddleOperator {
    a: DenseMatrix,
    c: SparseMatrix,
    ct: SparseMatrix,
}

impl DenseSaddleOperator {
    /// The operator of the scalar metric `a` (`|V| × |V|`) and the positional Jacobian `c`
    /// (`k_Q × 3|V|`).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] if the sizes do not match.
    pub fn new(a: &DenseMatrix, c: &SparseMatrix) -> Result<DenseSaddleOperator, RsError> {
        if c.ncols() != 3 * a.n() {
            return Err(RsError::InvalidOption {
                name: "positional Jacobian",
                value: format!("{} columns", c.ncols()),
                requirement: "3|V| columns, matching the metric",
            });
        }
        Ok(DenseSaddleOperator {
            a: a.clone(),
            c: c.clone(),
            ct: c.transpose(),
        })
    }

    /// `3|V|`.
    pub fn n_dofs(&self) -> usize {
        3 * self.a.n()
    }

    /// `k_Q`.
    pub fn n_constraints(&self) -> usize {
        self.c.nrows()
    }

    /// The scalar metric `A`.
    pub fn metric(&self) -> &DenseMatrix {
        &self.a
    }

    /// The positional Jacobian `C_Q`.
    pub fn jacobian(&self) -> &SparseMatrix {
        &self.c
    }

    fn apply_metric(&self, x: &[f64], y: &mut [f64]) {
        let nv = self.a.n();
        // Sequential: `O(|V|²)` per product, small next to the dense factorization.
        let body = |(v, yv): (usize, &mut [f64])| {
            let row = self.a.row(v);
            let mut s = [0.0; 3];
            for (w, &aw) in row.iter().enumerate() {
                s[0] += aw * x[3 * w];
                s[1] += aw * x[3 * w + 1];
                s[2] += aw * x[3 * w + 2];
            }
            yv.copy_from_slice(&s);
        };
        y[..3 * nv].chunks_mut(3).enumerate().for_each(body);
    }
}

impl LinearOperator for DenseSaddleOperator {
    fn dim(&self) -> usize {
        self.n_dofs() + self.n_constraints()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        let nd = self.n_dofs();
        let k = self.n_constraints();
        if x.len() != nd + k || y.len() != nd + k {
            return Err(RsError::InvalidOption {
                name: "saddle operator vector",
                value: format!("({}, {})", x.len(), y.len()),
                requirement: "3|V| + k_Q entries",
            });
        }
        self.apply_metric(x, y);
        let mut ctl = vec![0.0; nd];
        self.ct
            .mul_vec(&x[nd..], &mut ctl)
            .map_err(|e| sparse_error(e, "C^T"))?;
        for i in 0..nd {
            y[i] += ctl[i];
        }
        self.c
            .mul_vec(&x[..nd], &mut y[nd..])
            .map_err(|e| sparse_error(e, "C"))?;
        Ok(())
    }
}

fn sparse_error(e: ddg_core::sparse::SparseError, what: &'static str) -> RsError {
    RsError::InvalidOption {
        name: what,
        value: e.to_string(),
        requirement: "matching sparse dimensions",
    }
}

#[derive(Debug)]
enum Factor {
    /// Three copies of `[[A, C₁ᵀ], [C₁, 0]]`.
    PerAxis { lu: Lu, k1: usize },
    /// The whole `A′`.
    Full { lu: Lu },
}

/// The factored dense `A′` of one step (module docs).
#[derive(Debug)]
pub struct DenseSaddle {
    op: DenseSaddleOperator,
    factor: Factor,
    pivot_ratio: f64,
    factor_seconds: f64,
}

/// The equilibration of the oracle: `1/√K_jj` for primal unknowns (1 if `K_jj ≤ 0`), constraint
/// rows (after the primal scaling) to unit norm.
fn equilibration(k: &Mat<f64>, n_primal: usize) -> Vec<f64> {
    let n = k.nrows();
    let mut d = vec![1.0; n];
    for j in 0..n_primal {
        let x = k[(j, j)];
        if x > 0.0 && x.is_finite() {
            d[j] = 1.0 / x.sqrt();
        }
    }
    for r in n_primal..n {
        let s: f64 = (0..n_primal).map(|j| (k[(r, j)] * d[j]).powi(2)).sum();
        if s > 0.0 {
            d[r] = 1.0 / s.sqrt();
        }
    }
    d
}

impl DenseSaddle {
    /// Factors `A′` for the scalar metric `a` (`|V| × |V|`, symmetric) and the positional
    /// Jacobian `c` (`k_Q × 3|V|`) on `pool`.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] if the sizes do not match; [`RsError::SingularSaddle`] if the
    /// equilibrated pivot ratio is at most [`PIVOT_TOLERANCE`] (a component whose translations
    /// are not fixed, redundant rows, a singular metric on `ker C_Q`).
    pub fn new(
        pool: &rayon::ThreadPool,
        a: &DenseMatrix,
        c: &SparseMatrix,
    ) -> Result<DenseSaddle, RsError> {
        let t0 = Instant::now();
        let op = DenseSaddleOperator::new(a, c)?;
        let nv = a.n();
        let (factor, pivot_ratio) = match per_axis_rows(c, nv) {
            Some(c1) => {
                let k1 = c1.nrows();
                let lu = Self::factor_scalar(pool, a, &c1)?;
                (Factor::PerAxis { lu: lu.0, k1 }, lu.1)
            }
            None => {
                let lu = Self::factor_full(pool, a, c)?;
                (Factor::Full { lu: lu.0 }, lu.1)
            }
        };
        Ok(DenseSaddle {
            op,
            factor,
            pivot_ratio,
            factor_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    fn factor_scalar(
        pool: &rayon::ThreadPool,
        a: &DenseMatrix,
        c1: &SparseMatrix,
    ) -> Result<(Lu, f64), RsError> {
        let nv = a.n();
        let k1 = c1.nrows();
        let n = nv + k1;
        let mut m = Mat::<f64>::zeros(n, n);
        for j in 0..nv {
            let row = a.row(j); // symmetric: column j = row j
            for i in 0..nv {
                m[(i, j)] = row[i];
            }
        }
        let (cp, ri, val) = (c1.col_ptr(), c1.row_idx(), c1.values());
        for v in 0..nv {
            for p in cp[v]..cp[v + 1] {
                let r = nv + ri[p] as usize;
                m[(r, v)] = val[p];
                m[(v, r)] = val[p];
            }
        }
        let d = equilibration(&m, nv);
        Lu::new(pool, m, d)
    }

    fn factor_full(
        pool: &rayon::ThreadPool,
        a: &DenseMatrix,
        c: &SparseMatrix,
    ) -> Result<(Lu, f64), RsError> {
        let nv = a.n();
        let nd = 3 * nv;
        let k = c.nrows();
        let n = nd + k;
        let mut m = Mat::<f64>::zeros(n, n);
        for w in 0..nv {
            let row = a.row(w);
            for v in 0..nv {
                for s in 0..3 {
                    m[(3 * v + s, 3 * w + s)] = row[v];
                }
            }
        }
        let (cp, ri, val) = (c.col_ptr(), c.row_idx(), c.values());
        for j in 0..nd {
            for p in cp[j]..cp[j + 1] {
                let r = nd + ri[p] as usize;
                m[(r, j)] = val[p];
                m[(j, r)] = val[p];
            }
        }
        let d = equilibration(&m, nd);
        Lu::new(pool, m, d)
    }

    /// `min |u_ii| / max |u_ii|` of the equilibrated LU.
    pub fn pivot_ratio(&self) -> f64 {
        self.pivot_ratio
    }

    /// Whether the per-axis factorization is used.
    pub fn is_per_axis(&self) -> bool {
        matches!(self.factor, Factor::PerAxis { .. })
    }

    /// The operator `A′` (for residuals, GMRES tests).
    pub fn operator(&self) -> &DenseSaddleOperator {
        &self.op
    }

    /// One direct solve of `A′ [x; λ] = rhs` (no refinement).
    fn direct(&self, rhs: &[f64]) -> Vec<f64> {
        let nv = self.op.a.n();
        let nd = 3 * nv;
        match &self.factor {
            Factor::PerAxis { lu, k1 } => {
                let n = nv + k1;
                let mut r = Mat::<f64>::from_fn(n, 3, |i, axis| {
                    if i < nv {
                        rhs[3 * i + axis]
                    } else {
                        rhs[nd + 3 * (i - nv) + axis]
                    }
                });
                lu.solve(&mut r);
                let mut out = vec![0.0; nd + 3 * k1];
                for axis in 0..3 {
                    for i in 0..nv {
                        out[3 * i + axis] = r[(i, axis)];
                    }
                    for g in 0..*k1 {
                        out[nd + 3 * g + axis] = r[(nv + g, axis)];
                    }
                }
                out
            }
            Factor::Full { lu } => {
                let n = rhs.len();
                let mut r = Mat::<f64>::from_fn(n, 1, |i, _| rhs[i]);
                lu.solve(&mut r);
                (0..n).map(|i| r[(i, 0)]).collect()
            }
        }
    }

    /// `rhs − A′ y`.
    fn residual(&self, rhs: &[f64], y: &[f64]) -> Result<Vec<f64>, RsError> {
        let mut ay = vec![0.0; rhs.len()];
        self.op.apply(y, &mut ay)?;
        Ok(rhs.iter().zip(&ay).map(|(r, a)| r - a).collect())
    }

    /// Solves `A′ y = rhs` for a full right-hand side (`3|V| + k_Q`), with one step of
    /// iterative refinement; returns `y` and the relative residual.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a wrong length.
    pub fn solve_full(&self, rhs: &[f64]) -> Result<(Vec<f64>, f64), RsError> {
        let n = self.op.dim();
        if rhs.len() != n {
            return Err(RsError::InvalidOption {
                name: "saddle right-hand side",
                value: rhs.len().to_string(),
                requirement: "3|V| + k_Q entries",
            });
        }
        let mut y = self.direct(rhs);
        let r = self.residual(rhs, &y)?;
        let dy = self.direct(&r);
        for (a, b) in y.iter_mut().zip(&dy) {
            *a += b;
        }
        let r = self.residual(rhs, &y)?;
        let nb = norm2(rhs);
        let rel = if nb > 0.0 { norm2(&r) / nb } else { norm2(&r) };
        Ok((y, rel))
    }
}

impl SaddleSolver for DenseSaddle {
    fn n_dofs(&self) -> usize {
        self.op.n_dofs()
    }

    fn n_constraints(&self) -> usize {
        self.op.n_constraints()
    }

    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RsError> {
        check_rhs(self, a, b)?;
        let nd = self.n_dofs();
        let rhs: Vec<f64> = a.iter().chain(b).copied().collect();
        let (y, rel) = self.solve_full(&rhs)?;
        Ok(SaddleSolution {
            x: y[..nd].to_vec(),
            lambda: y[nd..].to_vec(),
            record: SolveRecord {
                iterations: 0,
                relative_residual: rel,
                converged: true,
                gmres: None,
            },
        })
    }

    fn base_report(&self) -> LinearReport {
        LinearReport {
            solver: if self.is_per_axis() {
                "dense LU (per axis)"
            } else {
                "dense LU"
            },
            dimension: self.op.dim(),
            pivot_ratio: Some(self.pivot_ratio),
            schur_pivot: 1.0,
            factor_seconds: self.factor_seconds,
            ..Default::default()
        }
    }
}

/// `ε = ε_rel · tr(W) / tr(M)` of the preconditioner's Laplacian (`SPEC.md` §6.3, D14;
/// scale-invariant).
pub fn preconditioner_epsilon(w: &SparseMatrix, mass: &[f64], eps_rel: f64) -> f64 {
    let trw: f64 = w.diagonal().iter().sum();
    let trm: f64 = mass.iter().sum();
    eps_rel * trw / trm
}

/// The dense fractional preconditioner of §6.3:
///
/// ```text
/// P⁻¹ r = L′⁻¹ · blockdiag(L_{2−s} ⊗ I₃, 0) · L′⁻¹ r,    L′ = [[(W + ε M) ⊗ I₃, C_Qᵀ], [C_Q, 0]]
/// ```
///
/// (the multiplier part is zeroed between the two solves; `ε` from
/// [`preconditioner_epsilon`]). With stiffness matrices the mass factors cancel:
/// `W⁻¹ (M Δ^{2−s}) W⁻¹ = Δ^{−s} M⁻¹` for `Δ = M⁻¹ W`, so `P⁻¹` has the order `−2s` of `A⁻¹`
/// (`rs-preconditioner-order`; the approximation of `B` by it is cited,
/// `rs-preconditioner-approx`). Dense: small meshes, tests, the oracle.
#[derive(Debug)]
pub struct DensePreconditioner {
    lprime: DenseSaddle,
    lfr: DenseMatrix,
    eps: f64,
}

impl DensePreconditioner {
    /// Builds `L′` and `L_{2−s}` at `positions` for the positional Jacobian `c`
    /// (`k_Q × 3|V|`), `L_{2−s}` over face pairs (the SPEC) or vertex pairs
    /// ([`FractionalPairs`]).
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`]; [`RsError::NonFiniteEnergy`] for
    /// coincident barycenters; [`RsError::InvalidOption`] for a non-finite or negative
    /// `eps_rel` or a Jacobian of the wrong width; [`RsError::SingularSaddle`] for a singular
    /// `L′`.
    pub fn new(
        pool: &rayon::ThreadPool,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        c: &SparseMatrix,
        eps_rel: f64,
        pairs: FractionalPairs,
    ) -> Result<DensePreconditioner, RsError> {
        if !(eps_rel.is_finite() && eps_rel >= 0.0) {
            return Err(RsError::InvalidOption {
                name: "eps_rel",
                value: eps_rel.to_string(),
                requirement: "a finite non-negative regularization",
            });
        }
        let geom = FaceGeometry::new(surface, positions)?;
        let mass = vertex_areas(surface, &geom.areas);
        let ops = FaceOperators::from_geometry(surface, positions, geom);
        let w = cotan_stiffness(surface, positions)?;
        let eps = preconditioner_epsilon(&w, &mass, eps_rel);
        let mut lap = to_dense(&w);
        let nv = surface.n_vertices();
        for v in 0..nv {
            lap.as_mut_slice()[v * nv + v] += eps * mass[v];
        }
        let tau = 2.0 - exponents.s();
        let lfr = match pairs {
            FractionalPairs::Faces => fractional_operator(pool, &ops, tau)?,
            FractionalPairs::Vertices => vertex_fractional_operator(pool, positions, &mass, tau)?,
        };
        let lprime = DenseSaddle::new(pool, &lap, c)?;
        Ok(DensePreconditioner { lprime, lfr, eps })
    }

    /// `ε` used in `L′`.
    pub fn epsilon(&self) -> f64 {
        self.eps
    }

    /// `L_{2−s}` (`|V| × |V|`).
    pub fn fractional(&self) -> &DenseMatrix {
        &self.lfr
    }

    /// The factored `L′`.
    pub fn laplacian_saddle(&self) -> &DenseSaddle {
        &self.lprime
    }

    /// `P⁻¹ r` for a full vector (`3|V| + k_Q`).
    ///
    /// [formula: rs-preconditioner-order]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a wrong length.
    pub fn apply_full(&self, r: &[f64]) -> Result<Vec<f64>, RsError> {
        let nv = self.lfr.n();
        let (y, _) = self.lprime.solve_full(r)?;
        let mut v = vec![0.0; r.len()];
        for i in 0..nv {
            let row = self.lfr.row(i);
            let mut s = [0.0; 3];
            for (j, &l) in row.iter().enumerate() {
                for k in 0..3 {
                    s[k] += l * y[3 * j + k];
                }
            }
            v[3 * i..3 * i + 3].copy_from_slice(&s);
        }
        let (z, _) = self.lprime.solve_full(&v)?;
        Ok(z)
    }
}

impl Preconditioner for DensePreconditioner {
    fn dim(&self) -> usize {
        self.lprime.operator().dim()
    }

    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError> {
        let out = self.apply_full(r)?;
        if z.len() != out.len() {
            return Err(RsError::InvalidOption {
                name: "preconditioner output",
                value: z.len().to_string(),
                requirement: "3|V| + k_Q entries",
            });
        }
        z.copy_from_slice(&out);
        Ok(())
    }
}
