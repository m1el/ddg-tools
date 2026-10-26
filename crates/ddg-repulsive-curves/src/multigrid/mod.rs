//! Multigrid solver of the saddle systems of Eq. 22 and 23 (`SPEC.md` §8.4, App. B.3; D17,
//! D25).
//!
//! * [`coarsen`]: the black/white coarsening of curve networks and the prolongation `J`
//!   [formula: rc-coarsen] [formula: rc-prolongation].
//! * [`vcycle`]: the generic V-cycle with a CG smoother, the full-multigrid start and the
//!   flexible CG fallback over a [`vcycle::Hierarchy`] [formula: rc-multigrid].
//! * [`CurveHierarchy`]: the levels of one step. Level 0 is the curve of the step; level `ℓ+1`
//!   is the coarsening of level `ℓ` at the fine positions of its black vertices (D17), with its
//!   own constraint rows `C_ℓ` (§8.4 "Constraints per level"; a rank-deficient coarse `C_ℓ` drops
//!   rows there, reported), its projector `P_ℓ = I − C_ℓ†C_ℓ` and its rediscretized metric `Ā_ℓ`
//!   (not Galerkin): a [`HierarchicalMetric`] on every level but the coarsest, which assembles the
//!   dense `Ā_L` and factors `[Ā_L C_Lᵀ; C_L 0]` ([`DenseSaddle`]). The operators are
//!   `H_ℓ = P_ℓ Ā_ℓ P_ℓ`, `Π_ℓ = P_ℓ (J_ℓ ⊗ I_d) P_{ℓ+1}` and `Π_ℓᵀ` (the reference's
//!   `MatrixProjectorOperator`); the coarsest solve gives `x ∈ ker C_L` with `H_L x = P_L b`.
//! * [`MultigridSaddle`]: a [`SaddleSolver`] (Braess–Sarazin as in B.3.3): Eq. 30 by multigrid
//!   for `b = 0` (the projected gradient, Eq. 22), Eq. 33 for `b ≠ 0` (each chord iteration of the
//!   constraint projection, Eq. 23; no factorization to reuse, Sec. 7.4), through
//!   [`Projector::solve_with_kernel_solver`]; then `λ = (C Cᵀ)⁻¹ C (a − Ā x)`.
//! * [`MultigridSolvers`]: the [`SolverFactory`] that plugs the solver into a [`crate::Session`]
//!   as [`crate::Solver::Custom`]. The classic metrics (L², H¹, H²) go to the sparse direct
//!   solver of [`DirectSolvers`].
//!
//! Every reduction runs in a fixed order and the products of [`HierarchicalMetric`] and
//! [`DenseSaddle`] are deterministic, so results are bitwise independent of the thread count
//! (D25).
//!
//! The level structure follows the reference's `ConstraintProjectorDomain` and
//! `MatrixProjectorOperator` (`include/multigrid/`, MIT, Copyright (c) 2019 Christopher Yu) and
//! libgmultigrid (MIT, Copyright (c) 2020 Christopher Yu, `LICENSES/libgmultigrid-MIT.txt`).

use std::sync::{Arc, Mutex};
use std::time::Instant;

use ddg_core::sparse::SparseMatrix;

use crate::constraints::{Constraint, ConstraintSet};
use crate::energy::exact::make_pool;
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::metric::hierarchical::HierarchicalMetric;
use crate::metric::{DofLayout, GramQuadrature, Metric, sobolev_gram};
use crate::network::CurveNetwork;
use crate::optimize::mean_free_rows;
use crate::solve::{
    DenseSaddle, LinearReport, Projector, SaddleSolution, SaddleSolver, check_rhs, sparse_error,
};
use crate::spatial::bct::BctOptions;
use crate::step::{DEFAULT_DENSE_THREADS, DirectSolvers, SolverFactory, SolverRequest};

pub mod coarsen;
pub mod vcycle;

pub use coarsen::{Coarsening, ProlongationRow, coarsen, forced_black, mark_black};
pub use vcycle::{CycleOptions, CycleReport, Hierarchy};

/// Options of the multigrid solver (`SPEC.md` §8.4, §9.1, §10).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MgOptions {
    /// The block cluster tree of every level but the coarsest.
    pub bct: BctOptions,
    /// Stop coarsening once a level has at most this many vertices.
    pub coarsest_vertices: usize,
    /// Maximal number of levels (the finest included).
    pub max_levels: usize,
    /// CG smoothing iterations before and after the coarse correction (`ν`).
    pub smoothing_steps: usize,
    /// CG iterations per level in the full-multigrid start.
    pub initial_smoothing_steps: usize,
    /// Relative residual `‖P(b − Āx)‖ / ‖Pb‖` at which the iteration stops.
    pub tolerance: f64,
    /// Maximal number of V-cycles per solve (FCG iterations included).
    pub max_cycles: usize,
    /// Flexible CG with the V-cycle preconditioner from the start (the default, measured), or
    /// (`false`) stationary V-cycles with the FCG fallback on stagnation as `SPEC.md` §8.4
    /// describes.
    pub flexible_cg: bool,
}

impl Default for MgOptions {
    /// Coarsest 512 vertices, 20 levels, 12 + 12 CG smoothing (24 per level in the start),
    /// tolerance `1e-3`, 20 cycles, flexible CG from the start.
    ///
    /// `SPEC.md` §10 lists 3 + 3 smoothing and stationary V-cycles; measured on the
    /// Freedman–He–Wang unknot (592 to 9472 vertices, `runs_tmp/rc-wp6-notes.md`), 3 + 3 needs 4
    /// to 15 cycles, growing with the depth, while 12 + 12 with FCG needs 2 to 3 at every size
    /// and is as fast or faster. The reference: 12 smoothing iterations, 24 in the start,
    /// tolerance `1e-2` in the ∞-norm, `max(1, log₂|V| − 4)` levels, stationary cycles.
    fn default() -> Self {
        MgOptions {
            bct: BctOptions::default(),
            coarsest_vertices: 512,
            max_levels: 20,
            smoothing_steps: 12,
            initial_smoothing_steps: 24,
            tolerance: 1e-3,
            max_cycles: 20,
            flexible_cg: true,
        }
    }
}

impl MgOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] naming the option.
    pub fn validate(&self) -> Result<(), RcError> {
        self.bct.validate()?;
        if !(self.tolerance.is_finite() && self.tolerance > 0.0) {
            return Err(RcError::InvalidOption {
                name: "multigrid.tolerance",
                value: self.tolerance.to_string(),
                requirement: "a finite positive tolerance",
            });
        }
        if self.max_levels == 0 {
            return Err(RcError::InvalidOption {
                name: "multigrid.max_levels",
                value: "0".into(),
                requirement: "at least 1 level",
            });
        }
        if self.coarsest_vertices == 0 {
            return Err(RcError::InvalidOption {
                name: "multigrid.coarsest_vertices",
                value: "0".into(),
                requirement: "at least 1 vertex",
            });
        }
        Ok(())
    }

    fn cycle(&self) -> CycleOptions {
        CycleOptions {
            smoothing_steps: self.smoothing_steps,
            initial_smoothing_steps: self.initial_smoothing_steps,
            tolerance: self.tolerance,
            max_cycles: self.max_cycles,
            flexible_cg: self.flexible_cg,
        }
    }
}

/// The metric of one level.
#[derive(Debug)]
enum LevelMetric {
    Hierarchical(HierarchicalMetric),
    /// The coarsest level: dense `Ā` and the factored saddle matrix.
    Dense(Mutex<DenseSaddle>),
}

/// One level of a [`CurveHierarchy`].
#[derive(Debug)]
pub struct Level {
    network: CurveNetwork,
    positions: Vec<[f64; 3]>,
    constraints: Option<ConstraintSet>,
    jacobian: SparseMatrix,
    projector: Projector,
    dropped_rows: Vec<usize>,
    metric: Option<LevelMetric>,
    /// To the next (coarser) level.
    coarsening: Option<Coarsening>,
}

impl Level {
    /// The network of the level.
    pub fn network(&self) -> &CurveNetwork {
        &self.network
    }

    /// The positions of the level (fine positions of its vertices, D17).
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// The constraint set of the level (`None` on soft problems, whose rows are the mean-free
    /// rows `M_c` of §6.4, and on level 0, whose rows are the caller's Jacobian).
    pub fn constraints(&self) -> Option<&ConstraintSet> {
        self.constraints.as_ref()
    }

    /// The constraint rows `C_ℓ` used on this level (after the rank drops).
    pub fn jacobian(&self) -> &SparseMatrix {
        &self.jacobian
    }

    /// The projector `P_ℓ` onto `ker C_ℓ`.
    pub fn projector(&self) -> &Projector {
        &self.projector
    }

    /// Rows of the level's constraint Jacobian dropped because the coarse `C_ℓ` was rank
    /// deficient (indices into the rows before dropping).
    pub fn dropped_rows(&self) -> &[usize] {
        &self.dropped_rows
    }

    /// The coarsening to the next level (`None` on the coarsest).
    pub fn coarsening(&self) -> Option<&Coarsening> {
        self.coarsening.as_ref()
    }

    /// Whether the level uses the dense metric (the coarsest level).
    pub fn is_dense(&self) -> bool {
        matches!(self.metric, Some(LevelMetric::Dense(_)))
    }

    /// The hierarchical metric (all levels but the coarsest).
    pub fn hierarchical_metric(&self) -> Option<&HierarchicalMetric> {
        match &self.metric {
            Some(LevelMetric::Hierarchical(h)) => Some(h),
            _ => None,
        }
    }

    /// `Ā_ℓ x`.
    fn apply_metric(&self, x: &[f64]) -> Vec<f64> {
        match &self.metric {
            Some(LevelMetric::Hierarchical(h)) => h.apply_unchecked(x),
            Some(LevelMetric::Dense(d)) => {
                d.lock().unwrap_or_else(|e| e.into_inner()).apply_metric(x)
            }
            None => vec![0.0; x.len()],
        }
    }
}

/// The multigrid levels of one step (`SPEC.md` §8.4).
#[derive(Debug)]
pub struct CurveHierarchy {
    levels: Vec<Level>,
    layout: DofLayout,
    build_seconds: f64,
}

/// What [`CurveHierarchy::build`] gets.
#[derive(Clone, Copy)]
pub struct HierarchyInput<'a> {
    /// The network of the step.
    pub network: &'a CurveNetwork,
    /// The base curve of the step.
    pub positions: &'a [[f64; 3]],
    /// The exponents of the metric.
    pub exponents: Exponents,
    /// The quadrature of the metric (D2).
    pub quadrature: GramQuadrature,
    /// The layout of the unknowns.
    pub layout: DofLayout,
    /// The constraint rows `C` of level 0 (`k × d|V|`).
    pub jacobian: &'a SparseMatrix,
    /// The hard constraints the rows come from (`None`: the rows are the mean-free rows of a
    /// soft problem, §6.4, and coarse levels use the mean-free rows of the coarse curve).
    pub constraints: Option<&'a ConstraintSet>,
}

/// The coarse constraint set; a tangent constraint whose fine direction is infeasible on the
/// coarse curve (`s·X ≤ 0`) takes the coarse tangent as its direction (the coarse rows only
/// precondition).
fn coarse_constraint_set(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    mut list: Vec<Constraint>,
    layout: DofLayout,
) -> Result<ConstraintSet, RcError> {
    for _ in 0..=list.len() {
        match ConstraintSet::new(network, positions, &list, layout.planar()) {
            Err(RcError::InfeasibleConstraint { constraint, .. }) => {
                if let Some(Constraint::Tangent { direction, .. }) = list.get_mut(constraint)
                    && direction.is_some()
                {
                    *direction = None;
                    continue;
                }
                return ConstraintSet::new(network, positions, &list, layout.planar());
            }
            r => return r,
        }
    }
    ConstraintSet::new(network, positions, &list, layout.planar())
}

/// The projector of `c`; on a coarse level (`drop`) a rank-deficient `c` loses the largest row
/// of each reported dependency until it has full rank.
fn projector_with_drops(
    c: &SparseMatrix,
    drop: bool,
) -> Result<(SparseMatrix, Projector, Vec<usize>), RcError> {
    let k = c.nrows();
    let mut keep: Vec<u32> = (0..k as u32).collect();
    let cols: Vec<u32> = (0..c.ncols() as u32).collect();
    let mut dropped = Vec::new();
    loop {
        let cc = if keep.len() == k {
            c.clone()
        } else {
            c.submatrix(&keep, &cols)
                .map_err(|e| sparse_error(e, "coarse constraint Jacobian"))?
        };
        match Projector::new(&cc, None) {
            Ok(p) => return Ok((cc, p, dropped)),
            Err(RcError::RedundantConstraints { rows }) if drop && !rows.is_empty() => {
                let r = *rows.iter().max().unwrap_or(&0);
                dropped.push(keep[r] as usize);
                keep.remove(r);
            }
            Err(e) => return Err(e),
        }
    }
}

impl CurveHierarchy {
    /// Builds the levels (`SPEC.md` §8.4): coarsens until a level has at most
    /// `coarsest_vertices` vertices, no vertex can be removed, the coarse network would be
    /// degenerate (no interacting pair), or `max_levels` levels exist; then builds a
    /// [`HierarchicalMetric`] on `pool` for every level but the coarsest and factors the dense
    /// saddle matrix of the coarsest on `dense_pool`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for bad options; the errors of [`Projector::new`] on level 0
    /// ([`RcError::RedundantConstraints`]); [`RcError::SingularMetric`] if the coarsest saddle
    /// matrix is singular; the errors of the metric builds.
    pub fn build(
        pool: &Arc<rayon::ThreadPool>,
        dense_pool: &Arc<rayon::ThreadPool>,
        input: &HierarchyInput<'_>,
        options: &MgOptions,
    ) -> Result<CurveHierarchy, RcError> {
        let t0 = Instant::now();
        options.validate()?;
        let layout = input.layout;
        let n0 = input.network.n_vertices();
        if input.positions.len() != n0 {
            return Err(RcError::PositionCountMismatch {
                expected: n0,
                found: input.positions.len(),
            });
        }
        if input.jacobian.ncols() != layout.n_dofs(n0) {
            return Err(RcError::InvalidOption {
                name: "constraint Jacobian",
                value: format!("{} columns", input.jacobian.ncols()),
                requirement: "d|V| columns",
            });
        }
        let projector = Projector::new(input.jacobian, None)?;
        let mut levels = vec![Level {
            network: input.network.clone(),
            positions: input.positions.to_vec(),
            constraints: input.constraints.cloned(),
            jacobian: input.jacobian.clone(),
            projector,
            dropped_rows: Vec::new(),
            metric: None,
            coarsening: None,
        }];
        loop {
            let last = levels.last().expect("at least one level");
            if last.network.n_vertices() <= options.coarsest_vertices
                || levels.len() >= options.max_levels
            {
                break;
            }
            let forced = forced_black(last.network.n_vertices(), last.constraints.as_ref());
            let coarsening = match coarsen(&last.network, &forced) {
                Ok(Some(c)) => c,
                // Nothing to remove, or the coarse network is degenerate: stop here.
                Ok(None) | Err(_) => break,
            };
            let net = coarsening.coarse().clone();
            let pos = coarsening.coarse_positions(&last.positions)?;
            let (constraints, c) = match &last.constraints {
                Some(set) => {
                    let list = coarsening.coarse_constraints(set);
                    let cs = coarse_constraint_set(&net, &pos, list, layout)?;
                    let c = cs.jacobian(&net, &pos)?;
                    (Some(cs), c)
                }
                None => (None, mean_free_rows(&net, &pos, layout)?),
            };
            let (c, projector, dropped_rows) = projector_with_drops(&c, true)?;
            levels.last_mut().expect("at least one level").coarsening = Some(coarsening);
            levels.push(Level {
                network: net,
                positions: pos,
                constraints,
                jacobian: c,
                projector,
                dropped_rows,
                metric: None,
                coarsening: None,
            });
        }
        let n_levels = levels.len();
        for (l, level) in levels.iter_mut().enumerate() {
            let metric = if l + 1 == n_levels {
                let gram = sobolev_gram(
                    dense_pool,
                    &level.network,
                    &level.positions,
                    input.exponents,
                    input.quadrature,
                );
                let s = DenseSaddle::new(dense_pool.clone(), &gram.a, &level.jacobian, layout)
                    .map_err(|e| match e {
                        RcError::SingularSaddle { .. } => RcError::SingularMetric {
                            components: level.network.components().len(),
                        },
                        e => e,
                    })?;
                LevelMetric::Dense(Mutex::new(s))
            } else {
                LevelMetric::Hierarchical(HierarchicalMetric::build(
                    pool.clone(),
                    &level.network,
                    &level.positions,
                    input.exponents,
                    input.quadrature,
                    layout,
                    options.bct,
                )?)
            };
            level.metric = Some(metric);
        }
        Ok(CurveHierarchy {
            levels,
            layout,
            build_seconds: t0.elapsed().as_secs_f64(),
        })
    }

    /// The levels, finest first.
    pub fn levels(&self) -> &[Level] {
        &self.levels
    }

    /// The layout of the unknowns.
    pub fn layout(&self) -> DofLayout {
        self.layout
    }

    /// Wall time of [`CurveHierarchy::build`].
    pub fn build_seconds(&self) -> f64 {
        self.build_seconds
    }

    /// `Ā_ℓ x` (no projection) on `level`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a bad level or length.
    pub fn apply_metric(&self, level: usize, x: &[f64]) -> Result<Vec<f64>, RcError> {
        self.check(level, x.len())?;
        Ok(self.levels[level].apply_metric(x))
    }

    fn check(&self, level: usize, len: usize) -> Result<(), RcError> {
        if level >= self.levels.len() || len != self.len(level) {
            return Err(RcError::InvalidOption {
                name: "multigrid level vector",
                value: format!("level {level}, {len} entries"),
                requirement: "an existing level and d|V_ℓ| entries",
            });
        }
        Ok(())
    }
}

impl Hierarchy for CurveHierarchy {
    fn n_levels(&self) -> usize {
        self.levels.len()
    }

    fn len(&self, level: usize) -> usize {
        self.layout.n_dofs(self.levels[level].network.n_vertices())
    }

    /// `H_ℓ x = P_ℓ Ā_ℓ P_ℓ x`.
    ///
    /// [formula: rc-projected-system]
    fn apply(&self, level: usize, x: &[f64]) -> Result<Vec<f64>, RcError> {
        self.check(level, x.len())?;
        let lv = &self.levels[level];
        lv.projector.projected_operator(|w| lv.apply_metric(w), x)
    }

    /// `Π_ℓ x = P_ℓ (J_ℓ ⊗ I_d) P_{ℓ+1} x`.
    ///
    /// [formula: rc-prolongation]
    fn prolong(&self, level: usize, coarse: &[f64]) -> Result<Vec<f64>, RcError> {
        self.check(level + 1, coarse.len())?;
        let (fine, next) = (&self.levels[level], &self.levels[level + 1]);
        let cg = fine.coarsening.as_ref().ok_or(RcError::InvalidOption {
            name: "multigrid level",
            value: level.to_string(),
            requirement: "a level with a coarser level below it",
        })?;
        let pc = next.projector.project(coarse)?;
        let y = cg.prolong_unchecked(&pc, self.layout.dim());
        fine.projector.project(&y)
    }

    /// `Π_ℓᵀ y = P_{ℓ+1} (J_ℓᵀ ⊗ I_d) P_ℓ y`.
    ///
    /// [formula: rc-prolongation]
    fn restrict(&self, level: usize, fine: &[f64]) -> Result<Vec<f64>, RcError> {
        self.check(level, fine.len())?;
        let (lv, next) = (&self.levels[level], self.levels.get(level + 1));
        let (Some(cg), Some(next)) = (lv.coarsening.as_ref(), next) else {
            return Err(RcError::InvalidOption {
                name: "multigrid level",
                value: level.to_string(),
                requirement: "a level with a coarser level below it",
            });
        };
        let pf = lv.projector.project(fine)?;
        let y = cg.restrict_unchecked(&pf, self.layout.dim());
        next.projector.project(&y)
    }

    /// `[Ā_L C_Lᵀ; C_L 0] [x; µ] = [b; 0]` by the dense LU: `x ∈ ker C_L`, `H_L x = P_L b`.
    fn coarse_solve(&self, b: &[f64]) -> Result<Vec<f64>, RcError> {
        let l = self.levels.len() - 1;
        self.check(l, b.len())?;
        let lv = &self.levels[l];
        match &lv.metric {
            Some(LevelMetric::Dense(d)) => {
                let mut d = d.lock().unwrap_or_else(|e| e.into_inner());
                let zero = vec![0.0; d.n_constraints()];
                Ok(d.solve(b, &zero)?.x)
            }
            _ => Err(RcError::InvalidOption {
                name: "multigrid coarsest level",
                value: l.to_string(),
                requirement: "a dense coarsest level",
            }),
        }
    }
}

/// The multigrid [`SaddleSolver`] of one step (`SPEC.md` §8.4, B.3.3).
#[derive(Debug)]
pub struct MultigridSaddle {
    hierarchy: CurveHierarchy,
    options: MgOptions,
    c: SparseMatrix,
    ct: SparseMatrix,
    last: Vec<CycleReport>,
}

impl MultigridSaddle {
    /// Builds the hierarchy of `input` ([`CurveHierarchy::build`]).
    ///
    /// # Errors
    /// The errors of [`CurveHierarchy::build`].
    pub fn new(
        pool: &Arc<rayon::ThreadPool>,
        dense_pool: &Arc<rayon::ThreadPool>,
        input: &HierarchyInput<'_>,
        options: MgOptions,
    ) -> Result<MultigridSaddle, RcError> {
        let hierarchy = CurveHierarchy::build(pool, dense_pool, input, &options)?;
        Ok(MultigridSaddle {
            hierarchy,
            options,
            c: input.jacobian.clone(),
            ct: input.jacobian.transpose(),
            last: Vec::new(),
        })
    }

    /// The levels.
    pub fn hierarchy(&self) -> &CurveHierarchy {
        &self.hierarchy
    }

    /// The options.
    pub fn options(&self) -> &MgOptions {
        &self.options
    }

    /// The cycle reports of the last [`SaddleSolver::solve`] (one per kernel solve).
    pub fn last_reports(&self) -> &[CycleReport] {
        &self.last
    }

    /// `Ā x` with the finest level's metric.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a wrong length.
    pub fn apply_metric(&self, x: &[f64]) -> Result<Vec<f64>, RcError> {
        self.hierarchy.apply_metric(0, x)
    }
}

impl SaddleSolver for MultigridSaddle {
    fn n_dofs(&self) -> usize {
        self.c.ncols()
    }

    fn n_constraints(&self) -> usize {
        self.c.nrows()
    }

    /// Eq. 30 by multigrid (`P Ā P y = P r`, `x = P y`) inside Eq. 33 (`z = C† b`,
    /// `x = z + y(a − Ā z)`), then `λ = (C Cᵀ)⁻¹ C (a − Ā x)`.
    ///
    /// [formula: rc-projected-system]
    /// [formula: rc-projection-transform]
    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RcError> {
        check_rhs(&*self, a, b)?;
        let t0 = Instant::now();
        let h = &self.hierarchy;
        let level0 = &h.levels[0];
        let p0 = &level0.projector;
        let cycle = self.options.cycle();
        let mut reports = Vec::new();
        let x = p0.solve_with_kernel_solver(
            |v| level0.apply_metric(v),
            |r| {
                let pr = p0.project(r)?;
                let (y, rep) = vcycle::solve(h, &pr, &cycle)?;
                reports.push(rep);
                p0.project(&y)
            },
            a,
            b,
        )?;
        let ax = level0.apply_metric(&x);
        let resid: Vec<f64> = a.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let mut cr = vec![0.0; self.c.nrows()];
        self.c
            .mul_vec(&resid, &mut cr)
            .map_err(|e| sparse_error(e, "C"))?;
        let lambda = p0.solve_gram(&cr)?;
        let mut ctl = vec![0.0; a.len()];
        self.ct
            .mul_vec(&lambda, &mut ctl)
            .map_err(|e| sparse_error(e, "C^T"))?;
        let mut cx = vec![0.0; b.len()];
        self.c
            .mul_vec(&x, &mut cx)
            .map_err(|e| sparse_error(e, "C"))?;
        let top: f64 = resid.iter().zip(&ctl).map(|(r, c)| (c - r) * (c - r)).sum();
        let bot: f64 = cx.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum();
        let rhs: f64 = a.iter().chain(b).map(|v| v * v).sum();
        let res = (top + bot).sqrt();
        let relative_residual = if rhs > 0.0 { res / rhs.sqrt() } else { res };
        let report = LinearReport {
            solver: "multigrid",
            dimension: a.len() + b.len(),
            iterations: reports.iter().map(CycleReport::total_cycles).sum(),
            relative_residual,
            fallback_used: reports.iter().any(|r| r.fallback_used),
            levels: h.levels.len(),
            pivot_ratio: None,
            seconds: t0.elapsed().as_secs_f64(),
        };
        self.last = reports;
        Ok(SaddleSolution { x, lambda, report })
    }
}

/// The multigrid back end as a [`SolverFactory`] for [`crate::Solver::Custom`]: one
/// [`MultigridSaddle`] per step for the Sobolev metric; the classic metrics use the sparse
/// direct solver of [`DirectSolvers`].
#[derive(Debug, Clone)]
pub struct MultigridSolvers {
    pool: Arc<rayon::ThreadPool>,
    dense_pool: Arc<rayon::ThreadPool>,
    options: MgOptions,
}

impl MultigridSolvers {
    /// The factory on `pool` (block cluster trees and products) and `dense_pool` (the coarsest
    /// level's dense assembly and LU).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for bad options.
    pub fn new(
        pool: Arc<rayon::ThreadPool>,
        dense_pool: Arc<rayon::ThreadPool>,
        options: MgOptions,
    ) -> Result<MultigridSolvers, RcError> {
        options.validate()?;
        Ok(MultigridSolvers {
            pool,
            dense_pool,
            options,
        })
    }

    /// The factory with its own pools: `threads` (`None`: all cores) and `dense_threads`
    /// (`None`: `min(DEFAULT_DENSE_THREADS, threads)`).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for bad options or a zero thread count.
    pub fn with_threads(
        threads: Option<usize>,
        dense_threads: Option<usize>,
        options: MgOptions,
    ) -> Result<MultigridSolvers, RcError> {
        let pool = Arc::new(make_pool(threads)?);
        let n = pool.current_num_threads();
        let dense_pool = Arc::new(make_pool(Some(
            dense_threads.unwrap_or(DEFAULT_DENSE_THREADS.min(n)),
        ))?);
        MultigridSolvers::new(pool, dense_pool, options)
    }

    /// The options.
    pub fn options(&self) -> &MgOptions {
        &self.options
    }
}

impl SolverFactory for MultigridSolvers {
    fn name(&self) -> &'static str {
        "multigrid"
    }

    fn build(&self, r: &SolverRequest<'_>) -> Result<Box<dyn SaddleSolver>, RcError> {
        match r.metric {
            Metric::Sobolev { quadrature } => {
                let input = HierarchyInput {
                    network: r.network,
                    positions: r.positions,
                    exponents: r.exponents,
                    quadrature,
                    layout: r.layout,
                    jacobian: r.jacobian,
                    constraints: r.constraints,
                };
                Ok(Box::new(MultigridSaddle::new(
                    &self.pool,
                    &self.dense_pool,
                    &input,
                    self.options,
                )?))
            }
            Metric::L2 | Metric::H1 | Metric::H2 => {
                DirectSolvers::new(self.pool.clone(), self.dense_pool.clone()).build(r)
            }
        }
    }
}
