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
//! constants of `W` would stay a near-null space of `L′`, and `ε_multi_rel` (default `1e-2`, the
//! reference's absolute `1e-2` made scale-invariant) is used. [`SparsePreconditioner::epsilon`]
//! reports which.
//!
//! **The fractional factor (D32).** [`FractionalPairs::Vertices`] (the default of
//! [`PreconditionerOptions`]) uses the vertex-pair `L^V_{2−s}` ([`VertexFractionalTree`]);
//! [`FractionalPairs::Faces`] the face form of the paper and the reference
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

use super::dense::preconditioner_epsilon;
use super::gmres::Preconditioner;
use crate::error::RsError;
use crate::geometry::{FaceGeometry, cotan_stiffness, vertex_areas};
use crate::kernel::Exponents;
use crate::metric::hierarchical::{
    ClusterTreeOptions, FaceFractionalTree, FaceMetricTree, VertexFractionalTree,
};
use crate::metric::{FaceOperators, FractionalPairs};
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

/// The factored `L′ = [[K ⊗ I₃, C_Qᵀ], [C_Q, 0]]` of one step (module docs), `K` the sparse
/// symmetric positive definite Laplacian block (`|V| × |V|`).
#[derive(Debug)]
pub struct LaplacianSaddle {
    matrix: SparseMatrix,
    factor: Factorization,
    /// `Some(k₁)` for the per-axis factorization of size `|V| + k₁`.
    per_axis: Option<usize>,
    n_vertices: usize,
    n_constraints: usize,
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
        let symbolic = Symbolic::analyze_with_order(&matrix, Method::Ldlt, EliminationOrder::Custom(&order))
            .map_err(|e| sparse_error(e, "preconditioner analysis"))?;
        let factor = symbolic
            .factor(&matrix, Parallelism::Sequential)
            .map_err(|e| sparse_error(e, "preconditioner factorization"))?;
        Ok(LaplacianSaddle {
            matrix,
            factor,
            per_axis,
            n_vertices: nv,
            n_constraints: c.nrows(),
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

    /// Solves `M X = B` for `nrhs` stacked columns with one refinement step.
    fn solve_columns(&self, b: &[f64], nrhs: usize) -> Result<Vec<f64>, RsError> {
        let n = self.matrix.nrows();
        let mut x = b.to_vec();
        self.factor
            .solve_many_in_place(&mut x, nrhs, Parallelism::Sequential)
            .map_err(|e| sparse_error(e, "preconditioner solve"))?;
        let mut r = vec![0.0; n * nrhs];
        for col in 0..nrhs {
            let (xs, rs) = (&x[col * n..(col + 1) * n], &mut r[col * n..(col + 1) * n]);
            self.matrix
                .mul_vec(xs, rs)
                .map_err(|e| sparse_error(e, "preconditioner residual"))?;
            for (ri, bi) in rs.iter_mut().zip(&b[col * n..(col + 1) * n]) {
                *ri = bi - *ri;
            }
        }
        self.factor
            .solve_many_in_place(&mut r, nrhs, Parallelism::Sequential)
            .map_err(|e| sparse_error(e, "preconditioner solve"))?;
        for (xi, di) in x.iter_mut().zip(&r) {
            *xi += di;
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
        let nv = self.n_vertices;
        let nd = 3 * nv;
        if r.len() != self.dim() {
            return Err(RsError::InvalidOption {
                name: "preconditioner right-hand side",
                value: r.len().to_string(),
                requirement: "3|V| + k_Q entries",
            });
        }
        match self.per_axis {
            Some(k1) => {
                let n = nv + k1;
                let mut b = vec![0.0; 3 * n];
                for axis in 0..3 {
                    for v in 0..nv {
                        b[axis * n + v] = r[3 * v + axis];
                    }
                    for g in 0..k1 {
                        b[axis * n + nv + g] = r[nd + 3 * g + axis];
                    }
                }
                let x = self.solve_columns(&b, 3)?;
                let mut out = vec![0.0; r.len()];
                for axis in 0..3 {
                    for v in 0..nv {
                        out[3 * v + axis] = x[axis * n + v];
                    }
                    for g in 0..k1 {
                        out[nd + 3 * g + axis] = x[axis * n + nv + g];
                    }
                }
                Ok(out)
            }
            None => self.solve_columns(r, 1),
        }
    }
}

/// Options of [`SparsePreconditioner`] (`SPEC.md` §6.3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreconditionerOptions {
    /// `ε_rel` of `ε = ε_rel tr W / tr M` (default `1e-10`).
    pub eps_rel: f64,
    /// `ε_multi_rel`, used when some component's translations are not fixed on their own
    /// (default `1e-2`).
    pub eps_multi_rel: f64,
    /// The discretization of `L_{2−s}` (default [`FractionalPairs::Vertices`], D32).
    pub fractional: FractionalPairs,
    /// The block cluster tree of `L_{2−s}` (`χ`, leaf size).
    pub cluster: ClusterTreeOptions,
    /// How `L′⁻¹` and `L_{2−s}` are combined (default [`PreconditionerForm::Augmented`]).
    pub form: PreconditionerForm,
}

impl Default for PreconditionerOptions {
    fn default() -> Self {
        PreconditionerOptions {
            eps_rel: 1e-10,
            eps_multi_rel: 1e-2,
            fractional: FractionalPairs::Vertices,
            cluster: ClusterTreeOptions::default(),
            form: PreconditionerForm::Augmented,
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
    fractional: FractionalFactor,
    eps: f64,
    multi: bool,
    fractional_seconds: f64,
    applications: AtomicUsize,
    nanos: AtomicU64,
}

impl SparsePreconditioner {
    /// Builds `L′` (sparse `LDLᵀ`) and the hierarchical `L_{2−s}` at `positions` for the
    /// positional rows `c` (`k_Q × 3|V|`). `extra_laplacian` (`|V| × |V|`, symmetric positive
    /// semidefinite, already weighted) is added to the Laplacian block (the Willmore hook, D7).
    /// `faces`: a face tree and the face operators of the same positions to reuse for
    /// [`FractionalPairs::Faces`] (built here otherwise; unused for vertex pairs).
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
        let lprime = LaplacianSaddle::new(&k, c)?;
        let t0 = Instant::now();
        let tau = 2.0 - exponents.s();
        let fractional = match options.fractional {
            FractionalPairs::Vertices => FractionalFactor::Vertices(VertexFractionalTree::build(
                &pool,
                positions,
                &mass,
                tau,
                &options.cluster,
            )?),
            FractionalPairs::Faces => {
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
        Ok(SparsePreconditioner {
            pool,
            lprime,
            laplacian: k,
            form: options.form,
            fractional,
            eps,
            multi,
            fractional_seconds: t0.elapsed().as_secs_f64(),
            applications: AtomicUsize::new(0),
            nanos: AtomicU64::new(0),
        })
    }

    /// As [`SparsePreconditioner::new`], reusing the face tree of a [`FaceMetricTree`] for
    /// [`FractionalPairs::Faces`].
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
        let t0 = Instant::now();
        let nd = 3 * self.lprime.n_vertices;
        let y = self.lprime.solve(r)?;
        let mut v = vec![0.0; r.len()];
        match &self.fractional {
            FractionalFactor::Vertices(t) => t.apply(&self.pool, &y[..nd], &mut v[..nd])?,
            FractionalFactor::Faces(t) => t.apply(&self.pool, &y[..nd], &mut v[..nd])?,
        }
        if self.form == PreconditionerForm::Augmented {
            // v −= (K ⊗ I₃) y (K symmetric: column j = row j).
            let k = &self.laplacian;
            let (cp, ri, val) = (k.col_ptr(), k.row_idx(), k.values());
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
        let mut z = self.lprime.solve(&v)?;
        if self.form == PreconditionerForm::Augmented {
            for (zi, yi) in z.iter_mut().zip(&y) {
                *zi += yi;
            }
        }
        self.applications.fetch_add(1, Ordering::Relaxed);
        self.nanos
            .fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        Ok(z)
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
}
