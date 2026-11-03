//! The fractional preconditioner `P⁻¹ = L′⁻¹ · blockdiag(L_{2−s} ⊗ I₃, 0) · L′⁻¹` with a sparse
//! `L′` (`SPEC.md` §6.3, D14, D32; Sec. 5.2, 7.2).
//!
//! ```text
//! L′ = [ (W + ε M + E) ⊗ I₃   C_Qᵀ ]        ε = ε_rel · tr W / tr M
//!      [ C_Q                  0    ]
//! ```
//!
//! with the cotan stiffness `W`, the lumped mass `M`, the positional rows `C_Q` (barycenters,
//! pins) and an optional extra sparse term `E` (the hook for the weighted Willmore matrix
//! `w_W W M⁻¹ W`, D7). [`LaplacianSaddle`] factors `L′` once per step by a sparse `LDLᵀ`
//! (`ddg_core::sparse`, no pivoting) with the elimination order "fill-reducing order of the
//! Laplacian block, then the multipliers": the Laplacian block is positive definite (`ε > 0`), so
//! every pivot of that order is non-zero (positive, then negative for the multipliers). When
//! `C_Q = C₁ ⊗ I₃` (always for barycenters and pins) the scalar saddle `[[W + εM + E, C₁ᵀ],
//! [C₁, 0]]` of size `|V| + k_Q/3` is factored once and serves the three axes; any other `C_Q` is
//! factored as a whole. Every solve does one step of iterative refinement against the factored
//! matrix (the unpivoted factorization of a matrix with a near-null space of size `ε` loses
//! digits; the refinement restores them, measured in `tests/hierarchical.rs`).
//!
//! **`ε` (D14, §5.1).** `ε_rel` (default `1e-10`) when every component's translations are fixed
//! by positional rows that touch only that component (per-component barycenters, pins); otherwise
//! (a global barycenter on several components, `BarycenterMode::Global`) the per-component
//! constants of `W` would stay a near-null space of `L′`, and `ε_multi_rel` (default `1e-4`; the
//! reference uses an absolute `1e-2`, measured worse, `NOTES-WP6.md`) is used. [`SparsePreconditioner::epsilon`]
//! reports which.
//!
//! **The fractional factor (D32).** [`FractionalForm::VertexPairs`] (the default of
//! [`PreconditionerOptions`]) uses the vertex-pair `L^V_{2−s}` ([`VertexFractionalTree`]);
//! [`FractionalForm::FacePairs`] the face form of the paper and the reference
//! ([`FaceFractionalTree`]). Both hierarchical with the options' `χ`.
//!
//! **The form (proposed D33, WP6).** As printed (`PreconditionerForm::Congruence`), `P⁻¹` is
//! symmetric positive semidefinite and its output satisfies `C_Q z = 0`
//! (`rs-fractional-preconditioner`), but it is singular on `A′`: GMRES cannot produce the
//! multipliers (pin reaction forces, the net force of a non-translation-invariant gradient). The
//! default [`PreconditionerForm::Augmented`] adds the coupling of `L′⁻¹` back (see there): same
//! cost, nonsingular, `C_Q P⁻¹ [r; b] = b`. Its position block is the same `E L E`, of the order of
//! `A⁻¹` (`rs-preconditioner-order`).

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

use ddg_core::SparseMatrix;
use ddg_core::sparse::{
    EliminationOrder, Factorization, Method, Parallelism, SparseError, Symbolic,
    fill_reducing_ordering,
};

use super::coarse::CoarseSpace;
use super::dense::preconditioner_epsilon;
use super::gmres::Preconditioner;
use crate::error::RsError;
use crate::geometry::{FaceGeometry, cotan_stiffness, vertex_areas};
use crate::kernel::Exponents;
use crate::metric::hierarchical::{
    ClusterTreeOptions, FaceFractionalTree, FaceMetricTree, VertexFractionalTree,
};
use crate::metric::{FaceOperators, FractionalForm};
use crate::spatial::FaceTree;
use crate::surface::Surface;

fn sparse_error(e: SparseError, what: &'static str) -> RsError {
    match e {
        SparseError::ZeroPivot { step } => RsError::SingularSaddle { pivot: step },
        SparseError::NonFiniteSolution { index } => RsError::SingularSaddle { pivot: index },
        e => RsError::InvalidOption {
            name: what,
            value: e.to_string(),
            requirement: "a factorizable preconditioner Laplacian",
        },
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

/// Whether every component has a positional row supported on that component only (its
/// translations are then fixed per component; §5.1, D8).
pub fn translations_fixed_per_component(surface: &Surface, c: &SparseMatrix) -> bool {
    let vc = surface.vertex_components();
    let nc = surface.n_components();
    let mut fixed = vec![false; nc];
    let ct = c.transpose();
    let (cp, ri) = (ct.col_ptr(), ct.row_idx());
    for r in 0..c.nrows() {
        let mut comp: Option<u32> = None;
        let mut single = true;
        for p in cp[r]..cp[r + 1] {
            let v = ri[p] as usize / 3;
            let cv = vc[v];
            match comp {
                None => comp = Some(cv),
                Some(x) if x != cv => single = false,
                _ => {}
            }
        }
        if single && let Some(x) = comp {
            fixed[x as usize] = true;
        }
    }
    fixed.iter().all(|&f| f)
}

/// How `L′` is factored and solved (`NOTES-GPU2.md`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TriangularSolve {
    /// [`super::ldlt::SimplicialLdlt`]: faer's factor read out, solves in the fixed chunked order
    /// of [`super::ldlt`] (default; the GPU reproduces it bit for bit, GPU-2).
    #[default]
    Chunked,
    /// `ddg_core::sparse` (faer's supernodal solves; GPU-1 and before).
    Faer,
}

#[derive(Debug)]
enum LFactor {
    Faer(Factorization),
    Chunked(super::ldlt::SimplicialLdlt),
}

/// The factored `L′ = [[K ⊗ I₃, C_Qᵀ], [C_Q, 0]]` of one step (module docs), `K` the sparse
/// symmetric positive definite Laplacian block (`|V| × |V|`).
#[derive(Debug)]
pub struct LaplacianSaddle {
    matrix: SparseMatrix,
    factor: LFactor,
    /// `Some(k₁)` for the per-axis factorization of size `|V| + k₁`.
    per_axis: Option<usize>,
    n_vertices: usize,
    n_constraints: usize,
    /// Columns solved together (per-axis case: the three axes of one vector).
    group: usize,
    factor_seconds: f64,
}

impl LaplacianSaddle {
    /// Assembles and factors `L′` for the Laplacian block `k` and the positional rows `c`
    /// (`k_Q × 3|V|`).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for mismatched sizes; [`RsError::SingularSaddle`] if the
    /// factorization breaks down (a zero pivot: `K` not positive definite or dependent rows).
    pub fn new(k: &SparseMatrix, c: &SparseMatrix) -> Result<LaplacianSaddle, RsError> {
        Self::new_with(k, c, TriangularSolve::default())
    }

    /// [`LaplacianSaddle::new`] with the factorization `kind`.
    ///
    /// # Errors
    /// As [`LaplacianSaddle::new`].
    pub fn new_with(k: &SparseMatrix, c: &SparseMatrix, kind: TriangularSolve) -> Result<LaplacianSaddle, RsError> {
        let t0 = Instant::now();
        let nv = k.nrows();
        if k.ncols() != nv || c.ncols() != 3 * nv {
            return Err(RsError::InvalidOption {
                name: "preconditioner Laplacian",
                value: format!("{} x {} Laplacian, {} x {} rows", nv, k.ncols(), c.nrows(), c.ncols()),
                requirement: "a square |V| x |V| Laplacian and k_Q x 3|V| positional rows",
            });
        }
        let amd = fill_reducing_ordering(k).map_err(|e| sparse_error(e, "Laplacian ordering"))?;
        let ktrip: Vec<(u32, u32, f64)> = k.triplets().collect();
        let (matrix, order, per_axis) = match per_axis_rows(c, nv) {
            Some(c1) => {
                let k1 = c1.nrows();
                let mut trip = ktrip;
                for (r, v, x) in c1.triplets() {
                    let row = nv as u32 + r;
                    trip.push((row, v, x));
                    trip.push((v, row, x));
                }
                let m = SparseMatrix::from_triplets(nv + k1, nv + k1, &trip)
                    .map_err(|e| sparse_error(e, "preconditioner saddle"))?;
                let order: Vec<u32> = amd.iter().copied().chain(nv as u32..(nv + k1) as u32).collect();
                (m, order, Some(k1))
            }
            None => {
                let nd = 3 * nv;
                let kq = c.nrows();
                let mut trip = Vec::with_capacity(3 * ktrip.len() + 2 * c.nnz());
                for &(i, j, x) in &ktrip {
                    for a in 0..3 {
                        trip.push((3 * i + a, 3 * j + a, x));
                    }
                }
                for (r, col, x) in c.triplets() {
                    let row = nd as u32 + r;
                    trip.push((row, col, x));
                    trip.push((col, row, x));
                }
                let m = SparseMatrix::from_triplets(nd + kq, nd + kq, &trip)
                    .map_err(|e| sparse_error(e, "preconditioner saddle"))?;
                let order: Vec<u32> = amd
                    .iter()
                    .flat_map(|&v| [3 * v, 3 * v + 1, 3 * v + 2])
                    .chain(nd as u32..(nd + kq) as u32)
                    .collect();
                (m, order, None)
            }
        };
        let factor = match kind {
            TriangularSolve::Chunked => LFactor::Chunked(super::ldlt::SimplicialLdlt::factor(&matrix, &order)?),
            TriangularSolve::Faer => {
                let symbolic = Symbolic::analyze_with_order(&matrix, Method::Ldlt, EliminationOrder::Custom(&order))
                    .map_err(|e| sparse_error(e, "preconditioner analysis"))?;
                LFactor::Faer(
                    symbolic
                        .factor(&matrix, Parallelism::Sequential)
                        .map_err(|e| sparse_error(e, "preconditioner factorization"))?,
                )
            }
        };
        Ok(LaplacianSaddle {
            matrix,
            factor,
            per_axis,
            n_vertices: nv,
            n_constraints: c.nrows(),
            group: 3,
            factor_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// `3|V| + k_Q`.
    pub fn dim(&self) -> usize {
        3 * self.n_vertices + self.n_constraints
    }

    /// Whether the per-axis factorization is used.
    pub fn is_per_axis(&self) -> bool {
        self.per_axis.is_some()
    }

    /// Seconds of the assembly, analysis and factorization.
    pub fn factor_seconds(&self) -> f64 {
        self.factor_seconds
    }

    /// The factored matrix (the scalar saddle in the per-axis case).
    pub fn matrix(&self) -> &SparseMatrix {
        &self.matrix
    }

    /// The chunked-order factor ([`TriangularSolve::Chunked`]), if that is the one in use.
    pub fn chunked_factor(&self) -> Option<&super::ldlt::SimplicialLdlt> {
        match &self.factor {
            LFactor::Chunked(f) => Some(f),
            LFactor::Faer(_) => None,
        }
    }

    /// `|V|`.
    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    /// `k_Q`.
    pub fn n_constraints(&self) -> usize {
        self.n_constraints
    }

    /// Solves `M X = B` in place (`ncols` columns of `M`'s dimension, stacked), by the factor.
    fn factor_solve(&self, x: &mut [f64], ncols: usize) -> Result<(), RsError> {
        match &self.factor {
            LFactor::Faer(f) => f
                .solve_many_in_place(x, ncols, Parallelism::Sequential)
                .map_err(|e| sparse_error(e, "preconditioner solve")),
            LFactor::Chunked(f) => {
                let n = f.dim();
                let perm = f.perm();
                let mut xr = vec![0.0; n * ncols];
                for (new, &old) in perm.iter().enumerate() {
                    for c in 0..ncols {
                        xr[new * ncols + c] = x[c * n + old as usize];
                    }
                }
                f.solve_rows(&mut xr, ncols)?;
                for (new, &old) in perm.iter().enumerate() {
                    for c in 0..ncols {
                        x[c * n + old as usize] = xr[new * ncols + c];
                    }
                }
                Ok(())
            }
        }
    }

    /// Solves `M X = B` for `ncols` stacked columns at once (one pass over the factor for all
    /// of them) with one refinement step.
    fn solve_block(&self, b: &[f64], ncols: usize, refine: bool) -> Result<Vec<f64>, RsError> {
        let n = self.matrix.nrows();
        let mut x = b.to_vec();
        self.factor_solve(&mut x, ncols)?;
        if !refine {
            return Ok(x);
        }
        let mut r = vec![0.0; b.len()];
        for c in 0..ncols {
            let cols = c * n..(c + 1) * n;
            self.matrix
                .mul_vec(&x[cols.clone()], &mut r[cols.clone()])
                .map_err(|e| sparse_error(e, "preconditioner residual"))?;
        }
        for (ri, bi) in r.iter_mut().zip(b) {
            *ri = bi - *ri;
        }
        self.factor_solve(&mut r, ncols)?;
        for (xi, di) in x.iter_mut().zip(&r) {
            *xi += di;
        }
        Ok(x)
    }

    /// Solves `M X = B` for `nrhs` stacked columns in groups of `group` consecutive columns
    /// ([`LaplacianSaddle::solve_block`]; the groups in parallel on `pool` if given). The result
    /// depends on the grouping only, not on the pool.
    fn solve_columns(
        &self,
        pool: Option<&rayon::ThreadPool>,
        b: &[f64],
        nrhs: usize,
        group: usize,
        refine: bool,
    ) -> Result<Vec<f64>, RsError> {
        use rayon::prelude::*;
        let n = self.matrix.nrows();
        let group = group.max(1);
        let ng = nrhs.div_ceil(group);
        let solve = |g: usize| {
            let c0 = g * group;
            let c1 = (c0 + group).min(nrhs);
            self.solve_block(&b[c0 * n..c1 * n], c1 - c0, refine)
        };
        let cols: Vec<Result<Vec<f64>, RsError>> = match pool {
            Some(pool) if ng > 1 => pool.install(|| (0..ng).into_par_iter().map(solve).collect()),
            _ => (0..ng).map(solve).collect(),
        };
        let mut x = Vec::with_capacity(n * nrhs);
        for c in cols {
            x.extend(c?);
        }
        Ok(x)
    }

    /// Solves `L′ y = r` for a full vector (`3|V| + k_Q`, interleaved positions then
    /// multipliers).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a wrong length; [`RsError::SingularSaddle`] for a
    /// non-finite solution.
    pub fn solve(&self, r: &[f64]) -> Result<Vec<f64>, RsError> {
        self.solve_on(None, r)
    }

    /// As [`LaplacianSaddle::solve`], the three axes in parallel on `pool` (bitwise the same
    /// result).
    ///
    /// # Errors
    /// As [`LaplacianSaddle::solve`].
    pub fn solve_on(&self, pool: Option<&rayon::ThreadPool>, r: &[f64]) -> Result<Vec<f64>, RsError> {
        let mut out = self.solve_many_on(pool, &[r], true)?;
        Ok(out.pop().expect("one solution"))
    }

    /// [`LaplacianSaddle::solve`] for several vectors: the three axis columns of a vector (per-axis
    /// case) in one pass over the factor, the vectors in parallel on `pool`; with `refine` (as
    /// [`LaplacianSaddle::solve`]) one refinement step, without it the plain back substitution.
    /// A vector's result does not depend on the pool or on the other vectors.
    ///
    /// # Errors
    /// As [`LaplacianSaddle::solve`].
    pub fn solve_many_on(
        &self,
        pool: Option<&rayon::ThreadPool>,
        rs: &[&[f64]],
        refine: bool,
    ) -> Result<Vec<Vec<f64>>, RsError> {
        let nv = self.n_vertices;
        let nd = 3 * nv;
        for r in rs {
            if r.len() != self.dim() {
                return Err(RsError::InvalidOption {
                    name: "preconditioner right-hand side",
                    value: r.len().to_string(),
                    requirement: "3|V| + k_Q entries",
                });
            }
        }
        match self.per_axis {
            Some(k1) => {
                let n = nv + k1;
                let mut b = vec![0.0; 3 * n * rs.len()];
                for (q, r) in rs.iter().enumerate() {
                    for axis in 0..3 {
                        let col = &mut b[(3 * q + axis) * n..(3 * q + axis + 1) * n];
                        for v in 0..nv {
                            col[v] = r[3 * v + axis];
                        }
                        for g in 0..k1 {
                            col[nv + g] = r[nd + 3 * g + axis];
                        }
                    }
                }
                let x = self.solve_columns(pool, &b, 3 * rs.len(), self.group, refine)?;
                Ok((0..rs.len())
                    .map(|q| {
                        let mut out = vec![0.0; self.dim()];
                        for axis in 0..3 {
                            let col = &x[(3 * q + axis) * n..(3 * q + axis + 1) * n];
                            for v in 0..nv {
                                out[3 * v + axis] = col[v];
                            }
                            for g in 0..k1 {
                                out[nd + 3 * g + axis] = col[nv + g];
                            }
                        }
                        out
                    })
                    .collect())
            }
            None => {
                let n = self.dim();
                let b: Vec<f64> = rs.iter().flat_map(|r| r.iter().copied()).collect();
                let x = self.solve_columns(pool, &b, rs.len(), 1, refine)?;
                Ok(x.chunks(n).map(<[f64]>::to_vec).collect())
            }
        }
    }
}

/// Options of [`SparsePreconditioner`] (`SPEC.md` §6.3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreconditionerOptions {
    /// `ε_rel` of `ε = ε_rel tr W / tr M` (default `1e-10`).
    pub eps_rel: f64,
    /// `ε_multi_rel`, used when some component's translations are not fixed on their own
    /// (default `1e-4`; `1e-2` costs 1.5–6× the iterations, `NOTES-WP6.md`).
    pub eps_multi_rel: f64,
    /// The discretization of `L_{2−s}` (default [`FractionalForm::VertexPairs`], D32).
    pub fractional: FractionalForm,
    /// The block cluster tree of `L_{2−s}` (`χ`, leaf size).
    pub cluster: ClusterTreeOptions,
    /// How `L′⁻¹` and `L_{2−s}` are combined (default [`PreconditionerForm::Augmented`]).
    pub form: PreconditionerForm,
    /// The scale of `L_{2−s}` against the metric in the augmented form (default
    /// [`FractionalScale::Symbol`]; no effect on the congruence form, where it is a global factor).
    pub scale: FractionalScale,
    /// Add the coarse correction on the translations of components without positional rows
    /// of their own ([`super::coarse::CoarseSpace`]; default `true`; augmented form only).
    pub coarse: bool,
    /// One step of iterative refinement in every `L′` solve of an application (default
    /// `false`: the plain back substitution; the refinement doubles the cost of the solves and
    /// changed neither the GMRES iteration counts nor the solutions beyond the GMRES tolerance,
    /// `NOTES-WP17.md`). [`LaplacianSaddle::solve`] always refines.
    pub refine: bool,
    /// How `L′` is factored and solved (default [`TriangularSolve::Chunked`], GPU-2).
    pub triangular: TriangularSolve,
    /// With a GPU: run the whole fine preconditioner there (the `L′` solves too; default `true`;
    /// needs [`TriangularSolve::Chunked`], vertex pairs and `refine = false`, else only the
    /// fractional product runs on the GPU as in GPU-1). Bitwise the CPU's at determinism level 2.
    pub gpu_solves: bool,
    /// The smallest `|V|` for which [`PreconditionerOptions::gpu_solves`] applies (default
    /// [`GPU_SOLVES_MIN_VERTICES`]): below it the latency of the GPU solves (one launch per level of
    /// the elimination tree) exceeds the CPU's solves (`NOTES-GPU2.md`, Fig. 13 suite).
    pub gpu_solves_min_vertices: usize,
}

/// The default of [`PreconditionerOptions::gpu_solves_min_vertices`]: measured crossover of the
/// Fig. 13 suite between 10k and 16k faces (`NOTES-GPU2.md`).
pub const GPU_SOLVES_MIN_VERTICES: usize = 7000;

/// The scale of the fractional factor in the augmented preconditioner (`SPEC.md` §6.3,
/// `NOTES-SOLVERFIX.md`).
///
/// `L′⁻¹ L_{2−s} L′⁻¹` approximates `A⁻¹` only up to a constant: on a plane the symbols of
/// `B` (`Σ_{S≠T} |D u(S) − D u(T)|² a_S a_T ρ^{−2s}`), `L_{2−s}` and `W` are
/// `2 I(σ) |k|^{2s}`, `2 I(2−s) |k|^{2(2−s)}` and `|k|²` (`σ = s − 1`,
/// `I(t) = ∫_{R²} (1 − cos z₁) |z|^{−2−2t} dz = π Γ(1−t) / (t 4^t Γ(1+t))`), so
/// `W⁻¹ L_{2−s} W⁻¹ B ≈ 4 I(σ) I(2−s)` (≈ 200 at `s = 5/3`). The congruence form is blind to
/// it (a global factor), but the augmented form adds the exact constraint part (eigenvalue 1 on
/// the `2 k_Q` constraint directions) to the fractional part (eigenvalues ≈ that constant):
/// with many positional rows (pins) GMRES must resolve two clusters two orders of magnitude
/// apart, and restarted GMRES stagnates (the H of Fig. 22 with 24 pins: relative residual 0.96
/// after 200 iterations).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FractionalScale {
    /// `L_{2−s} / (w_TP · 4 I(σ) I(2−s))` (default).
    #[default]
    Symbol,
    /// `L_{2−s}` unscaled (WP6's augmented form).
    Unit,
}

/// `Γ(x)` for `x > 0` (Lanczos, g = 7, 9 terms; relative error below `1e-13` on `(0, 3]`).
fn gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * gamma(1.0 - x));
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    (2.0 * std::f64::consts::PI).sqrt() * t.powf(x + 0.5) * (-t).exp() * a
}

/// `I(t) = ∫_{R²} (1 − cos z₁) |z|^{−2−2t} dz = π Γ(1−t) / (t 4^t Γ(1+t))`, `t ∈ (0, 1)`.
fn gagliardo_constant(t: f64) -> f64 {
    std::f64::consts::PI * gamma(1.0 - t) / (t * 4f64.powf(t) * gamma(1.0 + t))
}

/// The planar symbol ratio `4 I(σ) I(2 − s)` of `W⁻¹ L_{2−s} W⁻¹ B` ([`FractionalScale`]).
///
/// [formula: rs-preconditioner-scale]
pub fn fractional_symbol_ratio(exponents: &Exponents) -> f64 {
    let s = exponents.s();
    4.0 * gagliardo_constant(s - 1.0) * gagliardo_constant(2.0 - s)
}

impl Default for PreconditionerOptions {
    fn default() -> Self {
        PreconditionerOptions {
            eps_rel: 1e-10,
            eps_multi_rel: 1e-4,
            fractional: FractionalForm::VertexPairs,
            cluster: ClusterTreeOptions::default(),
            form: PreconditionerForm::Augmented,
            scale: FractionalScale::Symbol,
            coarse: true,
            refine: false,
            triangular: TriangularSolve::Chunked,
            gpu_solves: true,
            gpu_solves_min_vertices: GPU_SOLVES_MIN_VERTICES,
        }
    }
}

/// How the preconditioner combines `L′⁻¹` and `L_{2−s}` (`SPEC.md` §6.3; module docs).
///
/// With `L′⁻¹ = [[E, F], [Fᵀ, G]]`:
///
/// * `Congruence` (§6.3 as printed, the paper, the reference): `P⁻¹ = L′⁻¹ D L′⁻¹`,
///   `D = blockdiag(L_{2−s} ⊗ I₃, 0)`. Symmetric PSD, `C_Q P⁻¹ r = 0`, but **singular on
///   `A′`**: its range has no multiplier component (`D` zeroes it, and `L_{2−s} 1 = 0` removes
///   the translations), so right-preconditioned GMRES cannot reach a solution with nonzero
///   multipliers (pins: the reaction forces; any force with `Σ g ≠ 0`) and stalls at the
///   residual of that component; left GMRES stops early with that true residual (measured:
///   0.08–0.7).
/// * `Augmented` (default): `P⁻¹ = L′⁻¹ + L′⁻¹ blockdiag((L_{2−s} − K) ⊗ I₃, 0) L′⁻¹
///   = [[E L E, F + E L F], [Fᵀ + Fᵀ L E, 2G + Fᵀ L F]]` (`L′⁻¹ blockdiag(K, 0) L′⁻¹ =
///   [[E, 0], [0, −G]]`, so the `E` block is replaced by `E L E` and the coupling `F` of the
///   positions to the multipliers is kept). It has the structure of `A′⁻¹` (the constraint
///   preconditioner of `Ã = K L⁻¹ K`, up to rank `k_Q`), is nonsingular and symmetric, and
///   reproduces the constraint rows: `C_Q P⁻¹ [r; b] = b`. Same cost (two `L′` solves, one
///   `L_{2−s}` product, one sparse product with `K`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PreconditionerForm {
    /// `L′⁻¹ + L′⁻¹ blockdiag(L_{2−s} − K, 0) L′⁻¹` (proposed, D33).
    #[default]
    Augmented,
    /// `L′⁻¹ blockdiag(L_{2−s}, 0) L′⁻¹` (the SPEC/paper/reference form).
    Congruence,
}

impl PreconditionerOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a non-finite or non-positive `eps_rel`/`eps_multi_rel`
    /// and invalid cluster options.
    pub fn validate(&self) -> Result<(), RsError> {
        for (name, v) in [("eps_rel", self.eps_rel), ("eps_multi_rel", self.eps_multi_rel)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(RsError::InvalidOption {
                    name,
                    value: v.to_string(),
                    requirement: "a finite positive regularization (the Laplacian block must be definite)",
                });
            }
        }
        self.cluster.validate()
    }
}

/// The fractional factor of the preconditioner.
#[derive(Debug, Clone)]
pub enum FractionalFactor {
    /// `L^V_{2−s}` over vertex pairs (D32).
    Vertices(VertexFractionalTree),
    /// `L_{2−s}` over face pairs (the paper, the reference).
    Faces(FaceFractionalTree),
}

/// The sparse fractional preconditioner of `SPEC.md` §6.3 (module docs), a [`Preconditioner`]
/// on `3|V| + k_Q` vectors.
///
/// [formula: rs-fractional-preconditioner]
#[derive(Debug)]
pub struct SparsePreconditioner {
    pool: Arc<rayon::ThreadPool>,
    lprime: LaplacianSaddle,
    /// The Laplacian block `K = W + εM + E` (for the augmented form).
    laplacian: SparseMatrix,
    form: PreconditionerForm,
    refine: bool,
    fractional: FractionalFactor,
    /// The factor applied to `L_{2−s}` in the augmented form ([`FractionalScale`]).
    fractional_scale: f64,
    /// Whether [`SparsePreconditioner::attach_coarse`] may build a coarse space.
    coarse_enabled: bool,
    coarse: Option<CoarseSpace>,
    eps: f64,
    multi: bool,
    fractional_seconds: f64,
    /// The GPU copy of the vertex-pair tree ([`SparsePreconditioner::with_gpu`]).
    gpu_fractional: Option<crate::gpu::GpuVertexFractional>,
    /// Whether [`PreconditionerOptions::gpu_solves`] asked for the solves on the GPU (and the
    /// mesh is large enough, [`PreconditionerOptions::gpu_solves_min_vertices`]).
    gpu_solves: bool,
    /// The fine preconditioner on the GPU (GPU-2).
    gpu_prec: Option<crate::gpu::GpuPreconditioner>,
    upload_seconds: f64,
    applications: AtomicUsize,
    nanos: AtomicU64,
}

impl SparsePreconditioner {
    /// Builds `L′` (sparse `LDLᵀ`) and the hierarchical `L_{2−s}` at `positions` for the
    /// positional rows `c` (`k_Q × 3|V|`). `extra_laplacian` (`|V| × |V|`, symmetric positive
    /// semidefinite, already weighted) is added to the Laplacian block (the Willmore hook, D7).
    /// `faces`: a face tree and the face operators of the same positions to reuse for
    /// [`FractionalForm::FacePairs`] (built here otherwise; unused for vertex pairs).
    ///
    /// # Errors
    /// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] for invalid
    /// options or sizes; [`RsError::SingularSaddle`] if `L′` cannot be factored;
    /// [`RsError::NonFiniteEnergy`] / [`RsError::InvalidOption`] for coincident barycenters /
    /// vertices.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        pool: Arc<rayon::ThreadPool>,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        c: &SparseMatrix,
        options: &PreconditionerOptions,
        extra_laplacian: Option<&SparseMatrix>,
        faces: Option<(&FaceTree, &FaceOperators)>,
    ) -> Result<SparsePreconditioner, RsError> {
        Self::new_on(pool, surface, positions, exponents, c, options, extra_laplacian, faces, None)
    }

    /// [`SparsePreconditioner::new`], with the vertex-pair tree's near field evaluated on the GPU
    /// `gpu_build` (level 1, `NOTES-GPU.md`; its products then run there).
    ///
    /// # Errors
    /// As [`SparsePreconditioner::new`]; [`RsError::Gpu`].
    #[allow(clippy::too_many_arguments)]
    pub fn new_on(
        pool: Arc<rayon::ThreadPool>,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        c: &SparseMatrix,
        options: &PreconditionerOptions,
        extra_laplacian: Option<&SparseMatrix>,
        faces: Option<(&FaceTree, &FaceOperators)>,
        gpu_build: Option<&Arc<crate::gpu::GpuContext>>,
    ) -> Result<SparsePreconditioner, RsError> {
        options.validate()?;
        let geom = FaceGeometry::new(surface, positions)?;
        let nv = surface.n_vertices();
        let mass = vertex_areas(surface, &geom.areas);
        let w = cotan_stiffness(surface, positions)?;
        let multi = !translations_fixed_per_component(surface, c);
        let eps_rel = if multi {
            options.eps_multi_rel
        } else {
            options.eps_rel
        };
        let eps = preconditioner_epsilon(&w, &mass, eps_rel);
        let mut trip: Vec<(u32, u32, f64)> = w.triplets().collect();
        for (v, &m) in mass.iter().enumerate() {
            trip.push((v as u32, v as u32, eps * m));
        }
        if let Some(e) = extra_laplacian {
            if e.nrows() != nv || e.ncols() != nv {
                return Err(RsError::InvalidOption {
                    name: "extra preconditioner term",
                    value: format!("{} x {}", e.nrows(), e.ncols()),
                    requirement: "|V| x |V|",
                });
            }
            trip.extend(e.triplets());
        }
        let k = SparseMatrix::from_triplets(nv, nv, &trip)
            .map_err(|e| sparse_error(e, "preconditioner Laplacian"))?;
        let lprime = LaplacianSaddle::new_with(&k, c, options.triangular)?;
        let t0 = Instant::now();
        let tau = 2.0 - exponents.s();
        let mut gpu_fractional = None;
        let fractional = match (options.fractional, gpu_build) {
            (FractionalForm::VertexPairs, Some(g)) => {
                let (t, dev) = VertexFractionalTree::build_on_gpu(&pool, positions, &mass, tau, &options.cluster, g)?;
                gpu_fractional = Some(dev);
                FractionalFactor::Vertices(t)
            }
            (FractionalForm::VertexPairs, None) => FractionalFactor::Vertices(VertexFractionalTree::build(
                &pool,
                positions,
                &mass,
                tau,
                &options.cluster,
            )?),
            (FractionalForm::FacePairs, _) => {
                let tree = match faces {
                    Some((t, ops)) => FaceFractionalTree::build(&pool, t, ops, tau, &options.cluster)?,
                    None => {
                        let t = FaceTree::build(surface, positions, &geom, options.cluster.leaf_size)?;
                        let ops = FaceOperators::from_geometry(surface, positions, geom);
                        FaceFractionalTree::build(&pool, &t, &ops, tau, &options.cluster)?
                    }
                };
                FractionalFactor::Faces(tree)
            }
        };
        let fractional_scale = match options.scale {
            FractionalScale::Symbol => 1.0 / fractional_symbol_ratio(exponents),
            FractionalScale::Unit => 1.0,
        };
        let mut sp = SparsePreconditioner {
            pool,
            lprime,
            laplacian: k,
            form: options.form,
            refine: options.refine,
            fractional,
            fractional_scale,
            coarse_enabled: options.coarse && options.form == PreconditionerForm::Augmented,
            coarse: None,
            eps,
            multi,
            fractional_seconds: t0.elapsed().as_secs_f64(),
            gpu_fractional,
            gpu_solves: options.gpu_solves && nv >= options.gpu_solves_min_vertices,
            gpu_prec: None,
            upload_seconds: 0.0,
            applications: AtomicUsize::new(0),
            nanos: AtomicU64::new(0),
        };
        if let Some(g) = gpu_build {
            sp.attach_gpu_solves(g)?;
        }
        Ok(sp)
    }

    /// Uploads `L′` with its schedule and `K` for [`PreconditionerOptions::gpu_solves`], if the
    /// configuration allows it (module docs of [`crate::gpu`], `NOTES-GPU2.md`).
    fn attach_gpu_solves(&mut self, gpu: &Arc<crate::gpu::GpuContext>) -> Result<(), RsError> {
        self.gpu_prec = None;
        if !self.gpu_solves || self.refine || self.gpu_fractional.is_none() {
            return Ok(());
        }
        let Some(l) = self.lprime.chunked_factor() else {
            return Ok(());
        };
        let t0 = Instant::now();
        self.gpu_prec = Some(crate::gpu::GpuPreconditioner::new(
            gpu,
            l,
            self.lprime.is_per_axis(),
            self.lprime.n_vertices,
            self.lprime.n_constraints,
            &self.laplacian,
            self.form == PreconditionerForm::Augmented,
            self.fractional_scale,
        )?);
        self.upload_seconds += t0.elapsed().as_secs_f64();
        Ok(())
    }

    /// The fine preconditioner on the GPU, if in use.
    pub fn gpu_preconditioner(&self) -> Option<&crate::gpu::GpuPreconditioner> {
        self.gpu_prec.as_ref()
    }

    /// Whether the `L′` solves run on the GPU.
    pub fn solves_on_gpu(&self) -> bool {
        self.gpu_prec.is_some()
    }

    /// As [`SparsePreconditioner::new`], reusing the face tree of a [`FaceMetricTree`] for
    /// [`FractionalForm::FacePairs`].
    ///
    /// # Errors
    /// As [`SparsePreconditioner::new`].
    #[allow(clippy::too_many_arguments)]
    pub fn with_metric_tree(
        pool: Arc<rayon::ThreadPool>,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        c: &SparseMatrix,
        options: &PreconditionerOptions,
        extra_laplacian: Option<&SparseMatrix>,
        tree: &FaceMetricTree,
    ) -> Result<SparsePreconditioner, RsError> {
        Self::new(
            pool,
            surface,
            positions,
            exponents,
            c,
            options,
            extra_laplacian,
            Some((tree.face_tree(), tree.operators())),
        )
    }

    /// [`SparsePreconditioner::with_metric_tree`] with the vertex-pair near field on the GPU
    /// `gpu_build` ([`SparsePreconditioner::new_on`]).
    ///
    /// # Errors
    /// As [`SparsePreconditioner::new_on`].
    #[allow(clippy::too_many_arguments)]
    pub fn with_metric_tree_on(
        pool: Arc<rayon::ThreadPool>,
        surface: &Surface,
        positions: &[[f64; 3]],
        exponents: &Exponents,
        c: &SparseMatrix,
        options: &PreconditionerOptions,
        extra_laplacian: Option<&SparseMatrix>,
        tree: &FaceMetricTree,
        gpu_build: Option<&Arc<crate::gpu::GpuContext>>,
    ) -> Result<SparsePreconditioner, RsError> {
        Self::new_on(
            pool,
            surface,
            positions,
            exponents,
            c,
            options,
            extra_laplacian,
            Some((tree.face_tree(), tree.operators())),
            gpu_build,
        )
    }

    /// Divides the fractional factor by the metric's weight `w_TP` too (`A = w_TP (B + …)`;
    /// [`FractionalScale::Symbol`] only).
    #[must_use]
    pub fn with_metric_weight(mut self, tp_weight: f64) -> Self {
        if tp_weight.is_finite() && tp_weight > 0.0 && self.fractional_scale != 1.0 {
            self.fractional_scale /= tp_weight;
            if let Some(g) = &mut self.gpu_prec {
                g.set_scale(self.fractional_scale);
            }
        }
        self
    }

    /// Runs the vertex-pair fractional products on the GPU of `gpu` (uploads the tree;
    /// `NOTES-GPU.md`). The face-pair form ([`FractionalForm::FacePairs`], the reference mode)
    /// stays on the CPU; so do the sparse `L′` solves.
    ///
    /// # Errors
    /// [`RsError::Gpu`].
    pub fn with_gpu(mut self, gpu: &Arc<crate::gpu::GpuContext>) -> Result<SparsePreconditioner, RsError> {
        if let FractionalFactor::Vertices(t) = &self.fractional {
            let t0 = Instant::now();
            self.gpu_fractional = Some(t.upload(gpu)?);
            self.upload_seconds = t0.elapsed().as_secs_f64();
        }
        self.attach_gpu_solves(gpu)?;
        Ok(self)
    }

    /// Seconds of the GPU upload (0 on the CPU).
    pub fn upload_seconds(&self) -> f64 {
        self.upload_seconds
    }

    /// The factor applied to `L_{2−s}` in the augmented form.
    pub fn fractional_scale(&self) -> f64 {
        self.fractional_scale
    }

    /// Builds the coarse correction on the translations of the components whose
    /// translations are not fixed by rows of their own ([`CoarseSpace`]), for the saddle
    /// operator `op` with positional rows `c`; nothing when disabled, when every component is
    /// fixed, or above [`super::coarse::COARSE_MAX_ENTRIES`]. Returns whether a coarse space
    /// is in use.
    ///
    /// # Errors
    /// The errors of [`CoarseSpace::build`].
    pub fn attach_coarse(
        &mut self,
        op: &dyn super::gmres::LinearOperator,
        surface: &Surface,
        c: &SparseMatrix,
    ) -> Result<bool, RsError> {
        self.coarse = None;
        if !self.coarse_enabled {
            return Ok(false);
        }
        let (free, rows) = CoarseSpace::free_components(surface, c);
        self.coarse = CoarseSpace::build(op, surface, c, &free, &rows)?;
        Ok(self.coarse.is_some())
    }

    /// The coarse space, if one is attached.
    pub fn coarse(&self) -> Option<&CoarseSpace> {
        self.coarse.as_ref()
    }

    /// `(ε, multi)`: the `ε` used in `L′` and whether it is the several-components value.
    pub fn epsilon(&self) -> (f64, bool) {
        (self.eps, self.multi)
    }

    /// The factored `L′`.
    pub fn laplacian_saddle(&self) -> &LaplacianSaddle {
        &self.lprime
    }

    /// The fractional factor.
    pub fn fractional(&self) -> &FractionalFactor {
        &self.fractional
    }

    /// Seconds of the `L_{2−s}` tree build.
    pub fn fractional_seconds(&self) -> f64 {
        self.fractional_seconds
    }

    /// Seconds of the `L′` assembly and factorization.
    pub fn factor_seconds(&self) -> f64 {
        self.lprime.factor_seconds()
    }

    /// `(applications, seconds)` of [`Preconditioner::apply`] so far.
    pub fn usage(&self) -> (usize, f64) {
        (
            self.applications.load(Ordering::Relaxed),
            self.nanos.load(Ordering::Relaxed) as f64 * 1e-9,
        )
    }

    /// `P⁻¹ r` for a full vector (`3|V| + k_Q`).
    ///
    /// [formula: rs-fractional-preconditioner]
    /// [formula: rs-preconditioner-order]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for a wrong length; [`RsError::SingularSaddle`].
    pub fn apply_full(&self, r: &[f64]) -> Result<Vec<f64>, RsError> {
        let mut out = self.apply_many(&[r])?;
        Ok(out.pop().expect("one result"))
    }

    /// [`SparsePreconditioner::apply_full`] for several vectors: the `L′` solves of all of them
    /// in parallel, one pass of the fractional tree; bitwise the results of `apply_full`.
    ///
    /// # Errors
    /// As [`SparsePreconditioner::apply_full`].
    pub fn apply_many(&self, rs: &[&[f64]]) -> Result<Vec<Vec<f64>>, RsError> {
        let t0 = Instant::now();
        let zs = match &self.coarse {
            Some(cs) => cs.apply_many(rs, |xs| self.apply_fine(xs))?,
            None => self.apply_fine(rs)?,
        };
        self.applications.fetch_add(rs.len(), Ordering::Relaxed);
        self.nanos
            .fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        Ok(zs)
    }

    /// `P⁻¹ r` for several vectors without the coarse correction.
    fn apply_fine(&self, rs: &[&[f64]]) -> Result<Vec<Vec<f64>>, RsError> {
        if let (Some(g), Some(f)) = (&self.gpu_prec, &self.gpu_fractional) {
            return g.apply_many(f, rs);
        }
        let nd = 3 * self.lprime.n_vertices;
        let ys = self.lprime.solve_many_on(Some(&self.pool), rs, self.refine)?;
        let mut vs: Vec<Vec<f64>> = rs.iter().map(|r| vec![0.0; r.len()]).collect();
        {
            let yp: Vec<&[f64]> = ys.iter().map(|y| &y[..nd]).collect();
            let mut vp: Vec<&mut [f64]> = vs.iter_mut().map(|v| &mut v[..nd]).collect();
            match (&self.fractional, &self.gpu_fractional) {
                (_, Some(g)) => g.apply_batch(&yp, &mut vp)?,
                (FractionalFactor::Vertices(t), None) => t.apply_batch(&self.pool, &yp, &mut vp)?,
                (FractionalFactor::Faces(t), None) => {
                    for (y, v) in yp.iter().zip(vp.iter_mut()) {
                        t.apply(&self.pool, y, v)?;
                    }
                }
            }
        }
        if self.form == PreconditionerForm::Augmented {
            // v = s L y − (K ⊗ I₃) y (s the fractional scale; K symmetric: column j = row j).
            let k = &self.laplacian;
            let (cp, ri, val) = (k.col_ptr(), k.row_idx(), k.values());
            for (y, v) in ys.iter().zip(vs.iter_mut()) {
                if self.fractional_scale != 1.0 {
                    v[..nd].iter_mut().for_each(|x| *x *= self.fractional_scale);
                }
                for j in 0..k.ncols() {
                    let mut acc = [0.0; 3];
                    for p in cp[j]..cp[j + 1] {
                        let i = 3 * ri[p] as usize;
                        for a in 0..3 {
                            acc[a] += val[p] * y[i + a];
                        }
                    }
                    for a in 0..3 {
                        v[3 * j + a] -= acc[a];
                    }
                }
            }
        }
        let vr: Vec<&[f64]> = vs.iter().map(Vec::as_slice).collect();
        let mut zs = self.lprime.solve_many_on(Some(&self.pool), &vr, self.refine)?;
        if self.form == PreconditionerForm::Augmented {
            for (z, y) in zs.iter_mut().zip(&ys) {
                for (zi, yi) in z.iter_mut().zip(y) {
                    *zi += yi;
                }
            }
        }
        Ok(zs)
    }
}

impl Preconditioner for SparsePreconditioner {
    fn dim(&self) -> usize {
        self.lprime.dim()
    }

    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError> {
        if z.len() != self.lprime.dim() {
            return Err(RsError::InvalidOption {
                name: "preconditioner output",
                value: z.len().to_string(),
                requirement: "3|V| + k_Q entries",
            });
        }
        let out = self.apply_full(r)?;
        z.copy_from_slice(&out);
        Ok(())
    }

    fn apply_batch(&self, rs: &[&[f64]], zs: &mut [&mut [f64]]) -> Result<(), RsError> {
        if zs.iter().any(|z| z.len() != self.lprime.dim()) || rs.len() != zs.len() {
            return Err(RsError::InvalidOption {
                name: "preconditioner output",
                value: format!("{} outputs", zs.len()),
                requirement: "3|V| + k_Q entries, one per input",
            });
        }
        let out = self.apply_many(rs)?;
        for (z, o) in zs.iter_mut().zip(&out) {
            z.copy_from_slice(o);
        }
        Ok(())
    }
}
