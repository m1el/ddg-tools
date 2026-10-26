//! V-cycle, full-multigrid start and the flexible CG fallback (`SPEC.md` §8.4, App. B.3;
//! D17, D25).
//!
//! Generic over a [`Hierarchy`] of levels `ℓ = 0 … L−1` (0 the finest): an operator `H_ℓ` per
//! level, a prolongation `Π_ℓ` from level `ℓ+1` to `ℓ` with the restriction `Π_ℓᵀ`, and a direct
//! solve on the coarsest level. There is no curve type in any signature (`SPEC.md` §16); the
//! curve hierarchy (`H_ℓ = P_ℓ Ā_ℓ P_ℓ`, `Π_ℓ = P_ℓ (J_ℓ ⊗ I₃) P_{ℓ+1}`) is
//! [`super::CurveHierarchy`].
//!
//! * **Smoother:** `ν` conjugate-gradient iterations on `H_ℓ x = b_ℓ` from the current guess
//!   (paper: "a standard conjugate gradient smoother"; the reference uses Eigen's CG with 12
//!   iterations).
//! * **V-cycle:** pre-smooth, restrict the residual, recurse (zero start), prolong and correct,
//!   post-smooth; the coarsest level is solved directly.
//! * **Start (full multigrid, B.3.1 "Initialization"):** restrict the right side down to the
//!   coarsest level, solve there, then prolong level by level with `ν_init` smoothing iterations
//!   (the reference's `VCycleInitGuess`, 24 iterations).
//! * **Iteration:** V-cycles until `‖b − H₀x‖ ≤ tol ‖b‖` (Euclidean) or `max_cycles`. If
//!   [`STAGNATION_CYCLES`] cycles in a row each reduce the residual by less than 10 %, the
//!   remaining budget runs flexible CG (Notay 2000) with one V-cycle (zero start) as the variable
//!   preconditioner (ours, D17). With [`CycleOptions::flexible_cg`] the iteration is flexible CG
//!   from the start (the measured default of the curve solver).
//!
//! Every vector operation is sequential in a fixed order, so the result is bitwise independent
//! of the thread count when the hierarchy's operations are (D25).
//!
//! The cycle structure follows libgmultigrid (`multigrid_hierarchy.h`, `VCycleSolve`,
//! `VCycleInitGuess`; MIT, Copyright (c) 2020 Christopher Yu, notice in
//! `LICENSES/libgmultigrid-MIT.txt`); the stopping rule, the stagnation test and the FCG fallback
//! are ours.
//!
//! Convergence of the V-cycle and of FCG for this operator is cited, not proved
//! [formula: rc-multigrid].

use crate::error::RcError;

/// Cycles in a row with less than 10 % residual reduction that switch to flexible CG.
pub const STAGNATION_CYCLES: usize = 3;
/// A cycle "stagnates" when the residual is above this fraction of the previous one.
pub const STAGNATION_RATIO: f64 = 0.9;

/// A multigrid hierarchy: the operations the V-cycle needs on each level.
///
/// Level 0 is the finest, `n_levels() − 1` the coarsest. All operations must be deterministic
/// (bitwise independent of the thread count) for the solve to be (D25).
pub trait Hierarchy {
    /// Number of levels `L ≥ 1`.
    fn n_levels(&self) -> usize;
    /// Length of the vectors of `level`.
    fn len(&self, level: usize) -> usize;
    /// `H_ℓ x` on `level` (symmetric positive semidefinite, definite on the level's space).
    ///
    /// # Errors
    /// Back-end failures.
    fn apply(&self, level: usize, x: &[f64]) -> Result<Vec<f64>, RcError>;
    /// `Π_ℓ x`: a vector of level `level + 1` prolonged to `level`.
    ///
    /// # Errors
    /// Back-end failures.
    fn prolong(&self, level: usize, coarse: &[f64]) -> Result<Vec<f64>, RcError>;
    /// `Π_ℓᵀ y`: a vector of `level` restricted to `level + 1`.
    ///
    /// # Errors
    /// Back-end failures.
    fn restrict(&self, level: usize, fine: &[f64]) -> Result<Vec<f64>, RcError>;
    /// Solves `H x = b` on the coarsest level (directly).
    ///
    /// # Errors
    /// Back-end failures.
    fn coarse_solve(&self, b: &[f64]) -> Result<Vec<f64>, RcError>;
    /// The preconditioner of the smoother on `level` applied to `r` (symmetric positive
    /// definite on the level's space), or `None` for plain CG smoothing (the default).
    ///
    /// # Errors
    /// Back-end failures.
    fn precondition(&self, level: usize, r: &[f64]) -> Result<Option<Vec<f64>>, RcError> {
        let _ = (level, r);
        Ok(None)
    }
}

/// Options of [`solve`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CycleOptions {
    /// CG iterations before and after the coarse correction on each level (`ν`).
    pub smoothing_steps: usize,
    /// CG iterations per level when the full-multigrid start prolongs upward (`ν_init`).
    pub initial_smoothing_steps: usize,
    /// Stop when `‖b − H₀x‖ ≤ tolerance · ‖b‖`.
    pub tolerance: f64,
    /// Maximal number of V-cycles (FCG iterations included).
    pub max_cycles: usize,
    /// Run flexible CG with the V-cycle preconditioner from the start instead of stationary
    /// V-cycles (the fallback then never triggers).
    pub flexible_cg: bool,
}

/// What [`solve`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CycleReport {
    /// V-cycles run as stationary iterations.
    pub cycles: usize,
    /// Flexible CG iterations (one V-cycle each) after stagnation.
    pub fcg_iterations: usize,
    /// Whether the FCG fallback was used (stagnating stationary cycles).
    pub fallback_used: bool,
    /// `‖b − H₀x‖ / ‖b‖` after the start and after every cycle or FCG iteration.
    pub residuals: Vec<f64>,
    /// Whether the tolerance was reached.
    pub converged: bool,
}

impl CycleReport {
    /// V-cycles in total (stationary + FCG); the full-multigrid start is not counted.
    pub fn total_cycles(&self) -> usize {
        self.cycles + self.fcg_iterations
    }

    /// The last relative residual.
    pub fn final_residual(&self) -> f64 {
        self.residuals.last().copied().unwrap_or(0.0)
    }
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[inline]
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn axpy(y: &mut [f64], a: f64, x: &[f64]) {
    for (p, q) in y.iter_mut().zip(x) {
        *p += a * q;
    }
}

fn sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(p, q)| p - q).collect()
}

/// `steps` CG iterations on `H_level x = b` from `x`. Stops early on a zero residual or a
/// non-positive curvature `pᵀHp` (a smoother does not fail).
///
/// # Errors
/// The hierarchy's errors.
pub fn cg_smooth<H: Hierarchy + ?Sized>(
    h: &H,
    level: usize,
    b: &[f64],
    x: Vec<f64>,
    steps: usize,
) -> Result<Vec<f64>, RcError> {
    if steps == 0 {
        return Ok(x);
    }
    let r = sub(b, &h.apply(level, &x)?);
    Ok(cg_steps(h, level, x, r, steps)?.0)
}

/// CG iterations from `x` with its residual `r = b − Hx` given; returns `x` and the recursively
/// updated residual (so no extra product is needed for it).
fn cg_steps<H: Hierarchy + ?Sized>(
    h: &H,
    level: usize,
    mut x: Vec<f64>,
    mut r: Vec<f64>,
    steps: usize,
) -> Result<(Vec<f64>, Vec<f64>), RcError> {
    if steps == 0 {
        return Ok((x, r));
    }
    // Preconditioned CG when the level has a smoother preconditioner, plain CG otherwise.
    let pre = |r: &Vec<f64>| -> Result<Vec<f64>, RcError> {
        Ok(h.precondition(level, r)?.unwrap_or_else(|| r.clone()))
    };
    let mut z = pre(&r)?;
    let mut rz = dot(&r, &z);
    let mut p = z;
    for _ in 0..steps {
        if rz == 0.0 || !rz.is_finite() || rz < 0.0 {
            break;
        }
        let q = h.apply(level, &p)?;
        let pq = dot(&p, &q);
        if pq.is_nan() || pq <= 0.0 {
            break;
        }
        let alpha = rz / pq;
        axpy(&mut x, alpha, &p);
        axpy(&mut r, -alpha, &q);
        z = pre(&r)?;
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for (pk, zk) in p.iter_mut().zip(&z) {
            *pk = zk + beta * *pk;
        }
    }
    Ok((x, r))
}

/// One V-cycle on `level` for `H x = b` from the guess `x` (`None`: zero). The residuals are
/// carried through the smoothers (one product for the correction, none for a zero start).
///
/// # Errors
/// The hierarchy's errors.
pub fn v_cycle<H: Hierarchy + ?Sized>(
    h: &H,
    level: usize,
    b: &[f64],
    x: Option<Vec<f64>>,
    smoothing_steps: usize,
) -> Result<Vec<f64>, RcError> {
    if level + 1 >= h.n_levels() {
        return h.coarse_solve(b);
    }
    let (x, r) = match x {
        Some(x) => {
            let r = sub(b, &h.apply(level, &x)?);
            (x, r)
        }
        None => (vec![0.0; b.len()], b.to_vec()),
    };
    let (mut x, r) = cg_steps(h, level, x, r, smoothing_steps)?;
    let rc = h.restrict(level, &r)?;
    let xc = v_cycle(h, level + 1, &rc, None, smoothing_steps)?;
    let corr = h.prolong(level, &xc)?;
    let hc = h.apply(level, &corr)?;
    axpy(&mut x, 1.0, &corr);
    let r: Vec<f64> = r.iter().zip(&hc).map(|(p, q)| p - q).collect();
    Ok(cg_steps(h, level, x, r, smoothing_steps)?.0)
}

/// The full-multigrid start: restrict `b` to the coarsest level, solve there, then prolong
/// upward with `initial_smoothing_steps` CG iterations per level.
///
/// # Errors
/// The hierarchy's errors.
pub fn full_multigrid_start<H: Hierarchy + ?Sized>(
    h: &H,
    b: &[f64],
    initial_smoothing_steps: usize,
) -> Result<Vec<f64>, RcError> {
    let levels = h.n_levels();
    let mut bs = vec![b.to_vec()];
    for l in 0..levels - 1 {
        let next = h.restrict(l, &bs[l])?;
        bs.push(next);
    }
    let mut x = h.coarse_solve(&bs[levels - 1])?;
    for l in (0..levels - 1).rev() {
        x = h.prolong(l, &x)?;
        x = cg_smooth(h, l, &bs[l], x, initial_smoothing_steps)?;
    }
    Ok(x)
}

fn check_options(o: &CycleOptions) -> Result<(), RcError> {
    if !(o.tolerance.is_finite() && o.tolerance > 0.0) {
        return Err(RcError::InvalidOption {
            name: "multigrid tolerance",
            value: o.tolerance.to_string(),
            requirement: "a finite positive tolerance",
        });
    }
    Ok(())
}

/// Solves `H₀ x = b` (`SPEC.md` §8.4): full-multigrid start, V-cycles until the relative
/// residual is at most `tolerance` or `max_cycles` cycles ran, with the flexible CG fallback on
/// stagnation. Not reaching the tolerance is reported (`converged = false`), not an error.
///
/// [formula: rc-multigrid]
///
/// # Errors
/// [`RcError::InvalidOption`] for a bad tolerance or a right side of the wrong length;
/// [`RcError::MultigridDiverged`] when a residual becomes non-finite or the final residual
/// exceeds `‖b‖`; the hierarchy's errors.
pub fn solve<H: Hierarchy + ?Sized>(
    h: &H,
    b: &[f64],
    options: &CycleOptions,
) -> Result<(Vec<f64>, CycleReport), RcError> {
    solve_with_energy_check(h, b, options, false)
}

/// [`solve`]; with `energy_check` (flexible CG only) the iteration also requires the
/// energy-norm estimate `rᵀ z ≤ tolerance² · bᵀ x` (`z` the V-cycle of the residual, so
/// `rᵀ z ≈ ‖x* − x‖²_H` and `bᵀ x ≈ ‖x*‖²_H`) before it stops. The Euclidean residual alone can
/// be met while the smooth part of the error is still large when a few stiff couplings
/// (near-contacts) dominate `‖b‖` (`SPEC.md` §17g). The reported residuals stay Euclidean.
///
/// [formula: rc-multigrid]
///
/// # Errors
/// As [`solve`].
pub fn solve_with_energy_check<H: Hierarchy + ?Sized>(
    h: &H,
    b: &[f64],
    options: &CycleOptions,
    energy_check: bool,
) -> Result<(Vec<f64>, CycleReport), RcError> {
    check_options(options)?;
    if h.n_levels() == 0 || b.len() != h.len(0) {
        return Err(RcError::InvalidOption {
            name: "multigrid right-hand side",
            value: b.len().to_string(),
            requirement: "one entry per unknown of the finest level (and at least one level)",
        });
    }
    let mut report = CycleReport::default();
    let bn = norm(b);
    if bn == 0.0 {
        report.converged = true;
        report.residuals.push(0.0);
        return Ok((vec![0.0; b.len()], report));
    }
    let diverged = |report: &CycleReport| RcError::MultigridDiverged {
        residuals: report.residuals.clone(),
    };
    let mut x = full_multigrid_start(h, b, options.initial_smoothing_steps)?;
    let mut r = sub(b, &h.apply(0, &x)?);
    let mut rel = norm(&r) / bn;
    report.residuals.push(rel);
    if !rel.is_finite() {
        return Err(diverged(&report));
    }
    let nu = options.smoothing_steps;
    let mut stagnant = 0;
    while !options.flexible_cg
        && rel > options.tolerance
        && report.total_cycles() < options.max_cycles
    {
        if stagnant >= STAGNATION_CYCLES {
            report.fallback_used = true;
            break;
        }
        x = v_cycle(h, 0, b, Some(x), nu)?;
        r = sub(b, &h.apply(0, &x)?);
        let new = norm(&r) / bn;
        report.cycles += 1;
        report.residuals.push(new);
        if !new.is_finite() {
            return Err(diverged(&report));
        }
        stagnant = if new > STAGNATION_RATIO * rel {
            stagnant + 1
        } else {
            0
        };
        rel = new;
    }
    // The energy check runs in the flexible CG iteration only.
    let mut energy_ok = true;
    if report.fallback_used || options.flexible_cg {
        // Flexible CG (Notay 2000) with full orthogonalization against the previous search
        // directions (at most max_cycles of them), the V-cycle from zero as preconditioner.
        let mut dirs: Vec<(Vec<f64>, Vec<f64>, f64)> = Vec::new();
        let tol2 = options.tolerance * options.tolerance;
        energy_ok = !energy_check;
        while (rel > options.tolerance || !energy_ok) && report.total_cycles() < options.max_cycles
        {
            let z = v_cycle(h, 0, &r, None, nu)?;
            if energy_check {
                energy_ok = dot(&r, &z) <= tol2 * dot(b, &x).max(0.0);
                if energy_ok && rel <= options.tolerance {
                    break;
                }
            }
            let mut d = z.clone();
            for (dj, qj, dqj) in &dirs {
                let c = dot(&z, qj) / dqj;
                axpy(&mut d, -c, dj);
            }
            let q = h.apply(0, &d)?;
            let dq = dot(&d, &q);
            report.fcg_iterations += 1;
            if dq.is_nan() || dq <= 0.0 {
                report.residuals.push(rel);
                break;
            }
            let alpha = dot(&d, &r) / dq;
            axpy(&mut x, alpha, &d);
            axpy(&mut r, -alpha, &q);
            rel = norm(&r) / bn;
            report.residuals.push(rel);
            if !rel.is_finite() {
                return Err(diverged(&report));
            }
            dirs.push((d, q, dq));
        }
        // The true residual (the recursive one drifts by rounding).
        r = sub(b, &h.apply(0, &x)?);
        rel = norm(&r) / bn;
        if let Some(last) = report.residuals.last_mut() {
            *last = rel;
        }
    }
    if !rel.is_finite() || rel > 1.0 {
        return Err(diverged(&report));
    }
    report.converged = rel <= options.tolerance && energy_ok;
    Ok((x, report))
}
