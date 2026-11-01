//! Restarted GMRES, right or left preconditioned (`SPEC.md` §6.3, D13).
//!
//! Generic over a matrix-free [`LinearOperator`] and a [`Preconditioner`]. Initial guess
//! `x₀ = 0`, restart length `m` (default 30), Arnoldi with modified Gram–Schmidt, Givens
//! rotations for the small least-squares problem.
//!
//! * [`GmresSide::Right`] (default): GMRES on `A P⁻¹ u = b`, `x = P⁻¹ u`. It minimizes the true
//!   residual `‖b − A x‖` over the Krylov space, and stops when `‖b − A x‖ ≤ tol ‖b‖`.
//! * [`GmresSide::Left`] (the reference's Eigen GMRES): GMRES on `P⁻¹ A x = P⁻¹ b`; it stops on
//!   the preconditioned residual, `‖P⁻¹(b − A x)‖ ≤ tol ‖P⁻¹ b‖`. The true residual is reported
//!   as well.
//!
//! Inside a cycle the Arnoldi estimate decides; at the start of every cycle (and at the end) the
//! residual is recomputed explicitly, and only the explicit residual counts as converged. The
//! vectors are the interleaved `3|V| (+ k_Q)` unknowns, so the `x`, `y`, `z` columns go through
//! one Krylov space and every operator application sees all three columns at once (App. B; the
//! hierarchical product of WP6 multiplies them as one `|V| × 3` block). All reductions are
//! sequential in a fixed order: deterministic.
//!
//! Not converging within `max_iterations` is a reported condition ([`GmresReport::converged`]),
//! not an error, unless the true relative residual is above [`DIVERGENCE_THRESHOLD`]
//! ([`RsError::SolverDiverged`]).

use std::time::Instant;

use super::{LinearReport, SaddleSolution, SaddleSolver, SolveRecord, check_rhs, norm2};
use crate::error::RsError;

/// A square linear operator `y = A x`.
pub trait LinearOperator {
    /// The dimension `n`.
    fn dim(&self) -> usize;
    /// `y = A x` (`x.len() = y.len() = n`; `y` is overwritten).
    ///
    /// # Errors
    /// Back-end failures.
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError>;
}

/// A preconditioner `z = P⁻¹ r`.
pub trait Preconditioner {
    /// The dimension `n`.
    fn dim(&self) -> usize;
    /// `z = P⁻¹ r` (`z` is overwritten).
    ///
    /// # Errors
    /// Back-end failures.
    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError>;
}

/// `P = I`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IdentityPreconditioner {
    /// The dimension.
    pub n: usize,
}

impl Preconditioner for IdentityPreconditioner {
    fn dim(&self) -> usize {
        self.n
    }

    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError> {
        z.copy_from_slice(r);
        Ok(())
    }
}

/// A [`LinearOperator`] (and [`Preconditioner`]) from a closure.
pub struct FnOperator<F> {
    n: usize,
    f: F,
}

impl<F: Fn(&[f64], &mut [f64]) -> Result<(), RsError>> FnOperator<F> {
    /// The operator `y = f(x)` of dimension `n`.
    pub fn new(n: usize, f: F) -> Self {
        FnOperator { n, f }
    }
}

impl<F: Fn(&[f64], &mut [f64]) -> Result<(), RsError>> LinearOperator for FnOperator<F> {
    fn dim(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        (self.f)(x, y)
    }
}

impl<F: Fn(&[f64], &mut [f64]) -> Result<(), RsError>> Preconditioner for FnOperator<F> {
    fn dim(&self) -> usize {
        self.n
    }

    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError> {
        (self.f)(r, z)
    }
}

/// Which side the preconditioner is applied on (D13).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GmresSide {
    /// `A P⁻¹`; stops on the true residual (default).
    #[default]
    Right,
    /// `P⁻¹ A`; stops on the preconditioned residual (the reference).
    Left,
}

/// Options of [`gmres`] (`SPEC.md` §6.3, §9.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GmresOptions {
    /// Restart length `m` (default 30, as Eigen).
    pub restart: usize,
    /// Relative tolerance (default `1e-4`, the reference): of the true residual (right) or of
    /// the preconditioned residual (left).
    pub tolerance: f64,
    /// Maximal number of iterations (matrix–vector products in the Arnoldi process) over all
    /// cycles (default 200).
    pub max_iterations: usize,
    /// The side (default right).
    pub side: GmresSide,
}

impl Default for GmresOptions {
    fn default() -> Self {
        GmresOptions {
            restart: 30,
            tolerance: 1e-4,
            max_iterations: 200,
            side: GmresSide::Right,
        }
    }
}

impl GmresOptions {
    /// Checks the options.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for `restart = 0`, `max_iterations = 0` or a tolerance outside
    /// `(0, 1)`.
    pub fn validate(&self) -> Result<(), RsError> {
        if self.restart == 0 {
            return Err(RsError::InvalidOption {
                name: "gmres.restart",
                value: "0".into(),
                requirement: "at least 1",
            });
        }
        if self.max_iterations == 0 {
            return Err(RsError::InvalidOption {
                name: "gmres.max_iterations",
                value: "0".into(),
                requirement: "at least 1",
            });
        }
        if !(self.tolerance > 0.0 && self.tolerance < 1.0) {
            return Err(RsError::InvalidOption {
                name: "gmres.tolerance",
                value: self.tolerance.to_string(),
                requirement: "a relative tolerance in (0, 1)",
            });
        }
        Ok(())
    }
}

/// A true relative residual above this after `max_iterations` is [`RsError::SolverDiverged`]
/// (`SPEC.md` §6.3).
pub const DIVERGENCE_THRESHOLD: f64 = 1e-1;

/// How one GMRES solve went.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GmresReport {
    /// The side.
    pub side: GmresSide,
    /// Iterations (Arnoldi steps) over all cycles.
    pub iterations: usize,
    /// Number of cycles started (1 + restarts).
    pub cycles: usize,
    /// Whether the stopping measure reached the tolerance.
    pub converged: bool,
    /// The true relative residual `‖b − A x‖ / ‖b‖` at the end.
    pub relative_residual: f64,
    /// Left: `‖P⁻¹(b − A x)‖ / ‖P⁻¹ b‖` at the end (the stopping measure); right: `None`.
    pub preconditioned_relative_residual: Option<f64>,
    /// The stopping measure (relative) at the start of every cycle and at the end, explicit.
    pub cycle_residuals: Vec<f64>,
    /// The Arnoldi estimate of the stopping measure (relative) after every iteration.
    pub history: Vec<f64>,
    /// Products with the operator.
    pub operator_applications: usize,
    /// Applications of the preconditioner.
    pub preconditioner_applications: usize,
    /// Seconds.
    pub seconds: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solves `A x = b` by restarted GMRES (module docs); returns `x` and the report.
///
/// # Errors
/// [`RsError::InvalidOption`] for invalid options or mismatched dimensions;
/// [`RsError::SolverDiverged`] if the true relative residual is above
/// [`DIVERGENCE_THRESHOLD`] after `max_iterations` (or not finite); operator and preconditioner
/// errors.
pub fn gmres(
    op: &dyn LinearOperator,
    prec: &dyn Preconditioner,
    b: &[f64],
    options: &GmresOptions,
) -> Result<(Vec<f64>, GmresReport), RsError> {
    options.validate()?;
    let t0 = Instant::now();
    let n = op.dim();
    if prec.dim() != n || b.len() != n {
        return Err(RsError::InvalidOption {
            name: "gmres dimensions",
            value: format!("operator {n}, preconditioner {}, rhs {}", prec.dim(), b.len()),
            requirement: "equal dimensions",
        });
    }
    let mut rep = GmresReport {
        side: options.side,
        ..Default::default()
    };
    let mut x = vec![0.0; n];
    let bnorm = norm2(b);
    if bnorm == 0.0 {
        rep.converged = true;
        rep.preconditioned_relative_residual =
            (options.side == GmresSide::Left).then_some(0.0);
        rep.seconds = t0.elapsed().as_secs_f64();
        return Ok((x, rep));
    }
    let left = options.side == GmresSide::Left;
    // Reference norm of the stopping measure.
    let mut tmp = vec![0.0; n];
    let ref_norm = if left {
        prec.apply(b, &mut tmp)?;
        rep.preconditioner_applications += 1;
        norm2(&tmp)
    } else {
        bnorm
    };
    if !(ref_norm > 0.0 && ref_norm.is_finite()) {
        return Err(RsError::SolverDiverged {
            residual: f64::NAN,
        });
    }
    let m = options.restart.min(n.max(1));
    let tol = options.tolerance;
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
    // Hessenberg columns (j+2 entries each), Givens rotations, rhs of the LS problem.
    let mut h: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut cs = vec![0.0; m];
    let mut sn = vec![0.0; m];
    let mut g = vec![0.0; m + 1];
    let mut w = vec![0.0; n];
    let mut true_rel;
    loop {
        // Explicit residual r = b − A x (x = 0 in the first cycle).
        let mut r = b.to_vec();
        if rep.cycles > 0 {
            op.apply(&x, &mut tmp)?;
            rep.operator_applications += 1;
            for (ri, ai) in r.iter_mut().zip(&tmp) {
                *ri -= ai;
            }
        }
        true_rel = norm2(&r) / bnorm;
        let measure_vec = if left {
            prec.apply(&r, &mut tmp)?;
            rep.preconditioner_applications += 1;
            tmp.clone()
        } else {
            r
        };
        let beta = norm2(&measure_vec);
        let rel = beta / ref_norm;
        rep.cycle_residuals.push(rel);
        if !rel.is_finite() {
            break;
        }
        if rel <= tol {
            rep.converged = true;
            break;
        }
        if rep.iterations >= options.max_iterations {
            break;
        }
        rep.cycles += 1;
        basis.clear();
        h.clear();
        basis.push(measure_vec.iter().map(|v| v / beta).collect());
        g.iter_mut().for_each(|v| *v = 0.0);
        g[0] = beta;
        let mut k = 0;
        for j in 0..m {
            // w = A P⁻¹ v_j (right) or P⁻¹ A v_j (left).
            if left {
                op.apply(&basis[j], &mut tmp)?;
                prec.apply(&tmp, &mut w)?;
            } else {
                prec.apply(&basis[j], &mut tmp)?;
                op.apply(&tmp, &mut w)?;
            }
            rep.operator_applications += 1;
            rep.preconditioner_applications += 1;
            // Modified Gram–Schmidt.
            let mut col = vec![0.0; j + 2];
            for i in 0..=j {
                let hij = dot(&w, &basis[i]);
                col[i] = hij;
                for (wk, vk) in w.iter_mut().zip(&basis[i]) {
                    *wk -= hij * vk;
                }
            }
            let hn = norm2(&w);
            col[j + 1] = hn;
            // Previous rotations, then a new one to zero col[j+1].
            for i in 0..j {
                let t = cs[i] * col[i] + sn[i] * col[i + 1];
                col[i + 1] = -sn[i] * col[i] + cs[i] * col[i + 1];
                col[i] = t;
            }
            let denom = col[j].hypot(col[j + 1]);
            if denom == 0.0 {
                cs[j] = 1.0;
                sn[j] = 0.0;
            } else {
                cs[j] = col[j] / denom;
                sn[j] = col[j + 1] / denom;
            }
            col[j] = denom;
            col[j + 1] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            h.push(col);
            rep.iterations += 1;
            k = j + 1;
            let est = g[j + 1].abs() / ref_norm;
            rep.history.push(est);
            if hn == 0.0 || !hn.is_finite() || est <= tol || rep.iterations >= options.max_iterations
            {
                break;
            }
            basis.push(w.iter().map(|v| v / hn).collect());
        }
        // Back substitution H y = g (k × k upper triangular).
        let mut y = vec![0.0; k];
        for i in (0..k).rev() {
            let mut s = g[i];
            for l in i + 1..k {
                s -= h[l][i] * y[l];
            }
            y[i] = if h[i][i] != 0.0 { s / h[i][i] } else { 0.0 };
        }
        let mut u = vec![0.0; n];
        for (l, yl) in y.iter().enumerate() {
            for (ui, vi) in u.iter_mut().zip(&basis[l]) {
                *ui += yl * vi;
            }
        }
        if left {
            for (xi, ui) in x.iter_mut().zip(&u) {
                *xi += ui;
            }
        } else {
            prec.apply(&u, &mut tmp)?;
            rep.preconditioner_applications += 1;
            for (xi, ti) in x.iter_mut().zip(&tmp) {
                *xi += ti;
            }
        }
    }
    rep.relative_residual = true_rel;
    if left {
        rep.preconditioned_relative_residual = rep.cycle_residuals.last().copied();
    }
    rep.seconds = t0.elapsed().as_secs_f64();
    if !rep.converged && !(true_rel <= DIVERGENCE_THRESHOLD) {
        return Err(RsError::SolverDiverged { residual: true_rel });
    }
    Ok((x, rep))
}

/// A [`SaddleSolver`] that applies `A′⁻¹` by [`gmres`] with a matrix-free operator `A′`
/// (dimension `3|V| + k_Q`) and a preconditioner, so the Schur layer of
/// [`super::schur`] runs unchanged on top of an iterative back end (WP6: the hierarchical `A′`
/// and the sparse `L′` preconditioner).
pub struct GmresSaddle<O, P> {
    op: O,
    prec: P,
    n_dofs: usize,
    options: GmresOptions,
    label: &'static str,
}

impl<O: LinearOperator, P: Preconditioner> GmresSaddle<O, P> {
    /// The solver of `A′ = op` (with `n_dofs = 3|V|` position unknowns, the rest multipliers).
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for invalid options or `n_dofs > op.dim()` or a preconditioner
    /// of another dimension.
    pub fn new(op: O, prec: P, n_dofs: usize, options: GmresOptions) -> Result<Self, RsError> {
        options.validate()?;
        if n_dofs > op.dim() || prec.dim() != op.dim() {
            return Err(RsError::InvalidOption {
                name: "gmres saddle dimensions",
                value: format!("operator {}, preconditioner {}, dofs {n_dofs}", op.dim(), prec.dim()),
                requirement: "equal operator and preconditioner dimensions, at least 3|V|",
            });
        }
        Ok(GmresSaddle {
            op,
            prec,
            n_dofs,
            options,
            label: "gmres",
        })
    }

    /// Sets the solver name of the report.
    pub fn with_label(mut self, label: &'static str) -> Self {
        self.label = label;
        self
    }

    /// The operator.
    pub fn operator(&self) -> &O {
        &self.op
    }

    /// The preconditioner.
    pub fn preconditioner(&self) -> &P {
        &self.prec
    }
}

impl<O: LinearOperator, P: Preconditioner> SaddleSolver for GmresSaddle<O, P> {
    fn n_dofs(&self) -> usize {
        self.n_dofs
    }

    fn n_constraints(&self) -> usize {
        self.op.dim() - self.n_dofs
    }

    fn solve(&mut self, a: &[f64], b: &[f64]) -> Result<SaddleSolution, RsError> {
        check_rhs(self, a, b)?;
        let rhs: Vec<f64> = a.iter().chain(b).copied().collect();
        let (y, rep) = gmres(&self.op, &self.prec, &rhs, &self.options)?;
        Ok(SaddleSolution {
            x: y[..self.n_dofs].to_vec(),
            lambda: y[self.n_dofs..].to_vec(),
            record: SolveRecord {
                iterations: rep.iterations,
                relative_residual: rep.relative_residual,
                converged: rep.converged,
                gmres: Some(rep),
            },
        })
    }

    fn base_report(&self) -> LinearReport {
        LinearReport {
            solver: self.label,
            dimension: self.op.dim(),
            schur_pivot: 1.0,
            ..Default::default()
        }
    }
}

impl<T: LinearOperator + ?Sized> LinearOperator for &T {
    fn dim(&self) -> usize {
        (**self).dim()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), RsError> {
        (**self).apply(x, y)
    }
}

impl<T: Preconditioner + ?Sized> Preconditioner for &T {
    fn dim(&self) -> usize {
        (**self).dim()
    }

    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError> {
        (**self).apply(r, z)
    }
}
