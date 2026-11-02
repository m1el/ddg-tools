//! Iterative solves with the hierarchical metric (`SPEC.md` §6.3, §8.4; Sec. 5, 7.2).
//!
//! * [`HierarchicalOperator`]: the matrix-free saddle operator
//!
//!   ```text
//!   A′ [x; λ] = [ (w_TP (B̃ + c_lo B̃⁰) + w_W W M⁻¹ W + O_diag) ⊗ I₃ x + C_Qᵀ λ ;  C_Q x ]
//!   ```
//!
//!   on interleaved `3|V| + k_Q` vectors, with `B̃`, `B̃⁰` the block-cluster-tree products of
//!   [`FaceMetricTree`] (the three columns `x`, `y`, `z` and both kernels in one tree pass) and
//!   the two hooks of [`MetricExtras`]: the Willmore matrix (sparse, weighted by
//!   `SobolevWeights::willmore_weight`) and the mesh-obstacle diagonal [`FaceDiagonal`] (§4.4,
//!   D7, D23).
//! * [`hierarchical_saddle`]: GMRES ([`GmresSaddle`]) on that operator with the sparse
//!   preconditioner of [`super::precond`]; WP2's Schur layer ([`super::SchurComplement`]) runs on
//!   it unchanged.
//! * [`Backend`] and [`saddle_solver`]: `Dense`, `Hierarchical`, or `Auto` (dense up to
//!   `dense_max_faces` faces).
//!
//! What a step needs (after remeshing too, since the vertex set changes): the positions, the
//! exponents, the weights, the positional Jacobian `C_Q` of the step's constraint set, the extras;
//! [`saddle_solver`] builds the face tree, the block cluster tree, the factored `L′` and the
//! `L_{2−s}` tree from scratch. A Barnes–Hut tree of the same positions can donate its face tree
//! ([`hierarchical_saddle`]'s `face_tree`).

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::Instant;

use ddg_core::SparseMatrix;

use super::dense::DenseSaddle;
use super::gmres::{GmresOptions, GmresSaddle, LinearOperator};
use super::precond::{PreconditionerOptions, SparsePreconditioner};
use super::SaddleSolver;
use crate::error::RsError;
use crate::geometry::{FaceGeometry, cotan_stiffness, vertex_areas};
use crate::kernel::Exponents;
use crate::metric::dense::{dense_metric, willmore_matrix};
use crate::metric::hierarchical::{ClusterTreeOptions, FaceMetricTree};
use crate::metric::{FaceDiagonal, FaceOperators, FractionalForm, SobolevWeights};
use crate::spatial::FaceTree;
use crate::surface::Surface;

/// The optional terms of the metric that other work packages provide (`SPEC.md` §4.4).
#[derive(Debug, Clone, Default)]
pub struct MetricExtras {
    /// The (unweighted) Willmore matrix `W M⁻¹ W` (`|V| × |V|`, symmetric), used with
    /// `SobolevWeights::willmore_weight > 0`; `None` computes it from the cotan stiffness and
    /// the dual areas ([`willmore_matrix`]). WP7's penalty (interior rows of `W f`, §5.2) passes
    /// its own.
    pub willmore: Option<SparseMatrix>,
    /// The mesh-obstacle diagonal `O_diag` (D23), added unweighted.
    pub obstacle: Option<FaceDiagonal>,
}

impl MetricExtras {
    /// The weighted Willmore term `w_W W M⁻¹ W` (`None` for `w_W = 0`).
    ///
    /// # Errors
    /// The errors of [`cotan_stiffness`] and [`willmore_matrix`]; [`RsError::InvalidOption`] for
    /// a supplied matrix of the wrong size.
    pub fn weighted_willmore(
        &self,
        surface: &Surface,
        positions: &[[f64; 3]],
        weights: &SobolevWeights,
    ) -> Result<Option<SparseMatrix>, RsError> {
        let ww = weights.willmore_weight;
        if ww == 0.0 {
            return Ok(None);
        }
        let nv = surface.n_vertices();
        let mut m = match &self.willmore {
            Some(m) => {
                if m.nrows() != nv || m.ncols() != nv {
                    return Err(RsError::InvalidOption {
                        name: "Willmore matrix",
                        value: format!("{} x {}", m.nrows(), m.ncols()),
                        requirement: "|V| x |V|",
                    });
                }
                m.clone()
            }
            None => {
                let geom = FaceGeometry::new(surface, positions)?;
                let w = cotan_stiffness(surface, positions)?;
                willmore_matrix(&w, &vertex_areas(surface, &geom.areas))?
            }
        };
        m.values_mut().iter_mut().for_each(|v| *v *= ww);
        Ok(Some(m))
    }
}

/// `y += M ⊗ I₃ x` for a sparse symmetric `M` and interleaved vectors (columns in order).
fn add_kron3(m: &SparseMatrix, x: &[f64], y: &mut [f64]) {
    let (cp, ri, val) = (m.col_ptr(), m.row_idx(), m.values());
    for j in 0..m.ncols() {
        let xj = [x[3 * j], x[3 * j + 1], x[3 * j + 2]];
        for p in cp[j]..cp[j + 1] {
            let i = 3 * ri[p] as usize;
            let v = val[p];
            y[i] += v * xj[0];
            y[i + 1] += v * xj[1];
            y[i + 2] += v * xj[2];
        }
    }
}

fn sparse_err(e: ddg_core::SparseError, what: &'static str) -> RsError {
    RsError::InvalidOption {
        name: what,
        value: e.to_string(),
        requirement: "matching sparse dimensions",
    }
}

/// The hierarchical saddle operator `A′` (module docs), a [`LinearOperator`] on `3|V| + k_Q`
/// vectors.
#[derive(Debug)]
pub struct HierarchicalOperator {
    pool: Arc<rayon::ThreadPool>,
    tree: FaceMetricTree,
    weights: SobolevWeights,
    /// `w_W W M⁻¹ W` (weighted).
    willmore: Option<SparseMatrix>,
    obstacle: Option<FaceDiagonal>,
    c: SparseMatrix,
    ct: SparseMatrix,
    applications: AtomicUsize,
    nanos: AtomicU64,
}

impl HierarchicalOperator {
    /// The operator of the metric tree `tree`, the weights, the weighted Willmore term
    /// ([`MetricExtras::weighted_willmore`]), the obstacle diagonal and the positional rows `c`
    /// (`k_Q × 3|V|`).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for invalid weights or mismatched sizes.
    pub fn new(
        pool: Arc<rayon::ThreadPool>,
        tree: FaceMetricTree,
        weights: &SobolevWeights,
        willmore: Option<SparseMatrix>,
        obstacle: Option<FaceDiagonal>,
        c: &SparseMatrix,
    ) -> Result<HierarchicalOperator, RsError> {
        weights.validate()?;
        let nv = tree.n_vertices();
        if c.ncols() != 3 * nv {
            return Err(RsError::InvalidOption {
                name: "positional Jacobian",
                value: format!("{} columns", c.ncols()),
                requirement: "3|V| columns, matching the metric",
            });
        }
        if let Some(m) = &willmore
            && (m.nrows() != nv || m.ncols() != nv)
        {
            return Err(RsError::InvalidOption {
                name: "Willmore matrix",
                value: format!("{} x {}", m.nrows(), m.ncols()),
                requirement: "|V| x |V|",
            });
        }
        let nf = tree.operators().n_faces();
        if let Some(o) = &obstacle {
            for (name, d) in [("obstacle diagonal hi", &o.hi), ("obstacle diagonal lo", &o.lo)] {
                if !d.is_empty() && d.len() != nf {
                    return Err(RsError::InvalidOption {
                        name,
                        value: format!("{} entries", d.len()),
                        requirement: "one entry per face",
                    });
                }
            }
        }
        Ok(HierarchicalOperator {
            pool,
            tree,
            weights: *weights,
            willmore,
            obstacle,
            c: c.clone(),
            ct: c.transpose(),
            applications: AtomicUsize::new(0),
            nanos: AtomicU64::new(0),
        })
    }

    /// The metric tree.
    pub fn tree(&self) -> &FaceMetricTree {
        &self.tree
    }

    /// `3|V|`.
    pub fn n_dofs(&self) -> usize {
        3 * self.tree.n_vertices()
    }

    /// `k_Q`.
    pub fn n_constraints(&self) -> usize {
        self.c.nrows()
    }

    /// `(applications, seconds)` of [`LinearOperator::apply`] so far.
    pub fn usage(&self) -> (usize, f64) {
        (
            self.applications.load(Ordering::Relaxed),
            self.nanos.load(Ordering::Relaxed) as f64 * 1e-9,
        )
    }

    /// `y = A₃ x` (the metric block only) for interleaved `3|V|` vectors.
    ///
    /// [formula: rs-kernel-matrix-form]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths.
    pub fn apply_metric(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        let w = self.weights.tp_weight;
        self.tree.apply(
            &self.pool,
            w,
            w * self.weights.c_lo(),
            self.obstacle.as_ref(),
            x,
            y,
        )?;
        if let Some(m) = &self.willmore {
            add_kron3(m, x, y);
        }
        Ok(())
    }
}

impl LinearOperator for HierarchicalOperator {
    fn dim(&self) -> usize {
        self.n_dofs() + self.n_constraints()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        self.apply_batch(&[x], &mut [y])
    }

    /// One pass of the block cluster tree for all vectors (App. B); every output is bitwise
    /// that of [`LinearOperator::apply`].
    fn apply_batch(&self, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        let t0 = Instant::now();
        let nd = self.n_dofs();
        let n = self.dim();
        if xs.len() != ys.len() || xs.iter().zip(ys.iter()).any(|(x, y)| x.len() != n || y.len() != n) {
            return Err(RsError::InvalidOption {
                name: "saddle operator vector",
                value: format!("{} inputs, {} outputs", xs.len(), ys.len()),
                requirement: "3|V| + k_Q entries, one output per input",
            });
        }
        {
            let xp: Vec<&[f64]> = xs.iter().map(|x| &x[..nd]).collect();
            let mut yp: Vec<&mut [f64]> = ys.iter_mut().map(|y| &mut y[..nd]).collect();
            let w = self.weights.tp_weight;
            self.tree.apply_batch(
                &self.pool,
                w,
                w * self.weights.c_lo(),
                self.obstacle.as_ref(),
                &xp,
                &mut yp,
            )?;
        }
        for (x, y) in xs.iter().zip(ys.iter_mut()) {
            let (xp, xl) = x.split_at(nd);
            let (yp, yl) = y.split_at_mut(nd);
            if let Some(m) = &self.willmore {
                add_kron3(m, xp, yp);
            }
            if !xl.is_empty() {
                let mut ctl = vec![0.0; nd];
                self.ct
                    .mul_vec(xl, &mut ctl)
                    .map_err(|e| sparse_err(e, "C^T"))?;
                for (a, b) in yp.iter_mut().zip(&ctl) {
                    *a += b;
                }
                self.c.mul_vec(xp, yl).map_err(|e| sparse_err(e, "C"))?;
            }
        }
        self.applications.fetch_add(xs.len(), Ordering::Relaxed);
        self.nanos
            .fetch_add(t0.elapsed().as_nanos() as u64, Ordering::Relaxed);
        Ok(())
    }
}

/// Options of the hierarchical back end (`SPEC.md` §6.3, §8.4).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HierarchicalOptions {
    /// The block cluster trees of `B`, `B⁰` and `L_{2−s}` (`χ`, leaf size).
    pub cluster: ClusterTreeOptions,
    /// GMRES (right, restart 30, tolerance `1e-4`, 200 iterations).
    pub gmres: GmresOptions,
    /// `ε_rel`, `ε_multi_rel`, the fractional form (vertex pairs) — the cluster options of the
    /// preconditioner are replaced by [`HierarchicalOptions::cluster`].
    pub preconditioner: PreconditionerOptions,
    /// Add the weighted Willmore term to the preconditioner's Laplacian block too (D7; default
    /// `false`: measured harmful, 6 → 154 iterations at `w_W = 0.01` and divergence at `w_W ≥ 1` on
    /// a 10k torus, while without it the iterations stay 5–20, `NOTES-WP6.md`).
    pub willmore_in_preconditioner: bool,
    /// Add the mesh-obstacle term `O_diag` (as the sparse vertex matrix
    /// [`obstacle_vertex_matrix`]) to the preconditioner's Laplacian block (default `false`,
    /// `NOTES-WP17.md`).
    pub obstacle_in_preconditioner: bool,
}

impl HierarchicalOptions {
    /// The reference's preconditioner: face pairs, the congruence form, left GMRES.
    pub fn reference() -> Self {
        let mut o = Self::default();
        o.preconditioner.fractional = FractionalForm::FacePairs;
        o.preconditioner.form = super::precond::PreconditionerForm::Congruence;
        o.gmres.side = super::gmres::GmresSide::Left;
        o
    }
}

/// The solver of `A′` (`SPEC.md` §8.4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Backend {
    /// Dense assembly and LU (`O(|F|²)` memory; small meshes, the oracle).
    Dense,
    /// GMRES on the hierarchical operator with the sparse preconditioner.
    Hierarchical(HierarchicalOptions),
    /// Dense up to `dense_max_faces` faces, hierarchical above.
    Auto {
        /// The largest face count solved densely.
        dense_max_faces: usize,
        /// The options of the hierarchical back end.
        hierarchical: HierarchicalOptions,
    },
}

/// The default of [`Backend::Auto`]'s `dense_max_faces` (`SPEC.md` §6.3; measured by WP6,
/// `NOTES-WP6.md`).
pub const DEFAULT_DENSE_MAX_FACES: usize = 1000;

impl Default for Backend {
    fn default() -> Self {
        Backend::Auto {
            dense_max_faces: DEFAULT_DENSE_MAX_FACES,
            hierarchical: HierarchicalOptions::default(),
        }
    }
}

impl Backend {
    /// The hierarchical options this back end uses for a surface of `n_faces` faces, `None`
    /// for the dense path.
    pub fn resolve(&self, n_faces: usize) -> Option<HierarchicalOptions> {
        match self {
            Backend::Dense => None,
            Backend::Hierarchical(o) => Some(*o),
            Backend::Auto {
                dense_max_faces,
                hierarchical,
            } => (n_faces > *dense_max_faces).then_some(*hierarchical),
        }
    }
}

/// GMRES with the hierarchical operator and the sparse preconditioner.
pub type HierarchicalSaddle = GmresSaddle<HierarchicalOperator, SparsePreconditioner>;

/// Builds the hierarchical solver of one step: the face tree (or `face_tree`, which must be of
/// the same positions, e.g. the Barnes–Hut tree's), the block cluster tree of `B`, `B⁰`, the
/// operator `A′` with the extras, `L′` (sparse `LDLᵀ`) and the `L_{2−s}` tree.
///
/// # Errors
/// The errors of [`FaceMetricTree::build`], [`HierarchicalOperator::new`],
/// [`SparsePreconditioner::new`] and [`GmresSaddle::new`].
#[allow(clippy::too_many_arguments)]
pub fn hierarchical_saddle(
    pool: Arc<rayon::ThreadPool>,
    surface: &Surface,
    positions: &[[f64; 3]],
    exponents: &Exponents,
    weights: &SobolevWeights,
    extras: &MetricExtras,
    c_q: &SparseMatrix,
    options: &HierarchicalOptions,
    face_tree: Option<FaceTree>,
) -> Result<HierarchicalSaddle, RsError> {
    weights.validate()?;
    options.cluster.validate()?;
    let tree = match face_tree {
        Some(t) => {
            let geom = FaceGeometry::new(surface, positions)?;
            let ops = FaceOperators::from_geometry(surface, positions, geom);
            FaceMetricTree::from_parts(&pool, t, ops, exponents, &options.cluster)?
        }
        None => FaceMetricTree::build(&pool, surface, positions, exponents, &options.cluster)?,
    };
    let willmore = extras.weighted_willmore(surface, positions, weights)?;
    let mut popts = options.preconditioner;
    popts.cluster = options.cluster;
    let mut extra = if options.willmore_in_preconditioner {
        willmore.clone()
    } else {
        None
    };
    if options.obstacle_in_preconditioner
        && let Some(o) = &extras.obstacle
    {
        let om = obstacle_vertex_matrix(tree.operators(), o)?;
        extra = Some(match extra {
            None => om,
            Some(w) => {
                let mut t: Vec<(u32, u32, f64)> = w.triplets().collect();
                t.extend(om.triplets());
                SparseMatrix::from_triplets(w.nrows(), w.ncols(), &t)
                    .map_err(|e| sparse_err(e, "preconditioner extra term"))?
            }
        });
    }
    let prec = SparsePreconditioner::with_metric_tree(
        pool.clone(),
        surface,
        positions,
        exponents,
        c_q,
        &popts,
        extra.as_ref(),
        &tree,
    )?;
    let op = HierarchicalOperator::new(pool, tree, weights, willmore, extras.obstacle.clone(), c_q)?;
    let nd = op.n_dofs();
    Ok(GmresSaddle::new(op, prec, nd, options.gmres)?.with_label("gmres (hierarchical)"))
}

/// The obstacle term `O_diag = 2 Vᵀ diag(hi/a) V + 2 Uᵀ diag(lo/a) U` of `SPEC.md` §4.4 as a
/// scalar sparse vertex matrix `O` (`O_diag = O ⊗ I₃`):
/// `O_vw = 2 Σ_S a_S (hi_S g_c(v)(S) · g_c(w)(S) + lo_S / 9)` over the faces `S ∋ v, w`.
///
/// # Errors
/// [`RsError::InvalidOption`] for a diagonal of the wrong length.
pub fn obstacle_vertex_matrix(ops: &FaceOperators, d: &FaceDiagonal) -> Result<SparseMatrix, RsError> {
    let nf = ops.n_faces();
    for (name, v) in [("obstacle diagonal hi", &d.hi), ("obstacle diagonal lo", &d.lo)] {
        if !v.is_empty() && v.len() != nf {
            return Err(RsError::InvalidOption {
                name,
                value: format!("{} entries", v.len()),
                requirement: "one entry per face",
            });
        }
    }
    let mut trip = Vec::with_capacity(9 * nf);
    for s in 0..nf {
        let hi = d.hi.get(s).copied().unwrap_or(0.0);
        let lo = d.lo.get(s).copied().unwrap_or(0.0);
        let a = ops.areas[s];
        let g = &ops.df[s];
        let t = ops.faces[s];
        for i in 0..3 {
            for j in 0..3 {
                let gg = g[i][0] * g[j][0] + g[i][1] * g[j][1] + g[i][2] * g[j][2];
                trip.push((t[i], t[j], 2.0 * a * (hi * gg + lo / 9.0)));
            }
        }
    }
    SparseMatrix::from_triplets(ops.n_vertices, ops.n_vertices, &trip)
        .map_err(|e| sparse_err(e, "obstacle vertex matrix"))
}

/// The dense metric `A` with the extras (the supplied Willmore matrix replaces the computed one).
///
/// # Errors
/// The errors of [`dense_metric`] and [`MetricExtras::weighted_willmore`].
pub fn dense_metric_with_extras(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    exponents: &Exponents,
    weights: &SobolevWeights,
    extras: &MetricExtras,
) -> Result<crate::metric::DenseMatrix, RsError> {
    let mut w0 = *weights;
    w0.willmore_weight = 0.0;
    let mut a = dense_metric(pool, surface, positions, exponents, &w0, extras.obstacle.as_ref())?.a;
    if let Some(m) = extras.weighted_willmore(surface, positions, weights)? {
        let nv = a.n();
        let data = a.as_mut_slice();
        for (i, j, x) in m.triplets() {
            data[i as usize * nv + j as usize] += x;
        }
    }
    Ok(a)
}

/// The solver of `A′` for one step, by [`Backend`]: [`DenseSaddle`] or [`HierarchicalSaddle`]
/// (`Auto`: dense up to `dense_max_faces` faces).
///
/// # Errors
/// The errors of [`dense_metric_with_extras`] and [`DenseSaddle::new`], or of
/// [`hierarchical_saddle`].
#[allow(clippy::too_many_arguments)]
pub fn saddle_solver(
    pool: Arc<rayon::ThreadPool>,
    surface: &Surface,
    positions: &[[f64; 3]],
    exponents: &Exponents,
    weights: &SobolevWeights,
    extras: &MetricExtras,
    c_q: &SparseMatrix,
    backend: &Backend,
) -> Result<Box<dyn SaddleSolver + Send>, RsError> {
    match backend.resolve(surface.n_faces()) {
        None => {
            let a = dense_metric_with_extras(&pool, surface, positions, exponents, weights, extras)?;
            Ok(Box::new(DenseSaddle::new(&pool, &a, c_q)?))
        }
        Some(o) => Ok(Box::new(hierarchical_saddle(
            pool, surface, positions, exponents, weights, extras, c_q, &o, None,
        )?)),
    }
}
