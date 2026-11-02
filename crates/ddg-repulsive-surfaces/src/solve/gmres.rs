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

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

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
    /// `y_k = A x_k` for several vectors (default: one [`LinearOperator::apply`] each).
    /// Implementations that batch must give every `y_k` exactly as `apply` would.
    ///
    /// # Errors
    /// Back-end failures.
    fn apply_batch(&self, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        for (x, y) in xs.iter().zip(ys.iter_mut()) {
            self.apply(x, y)?;
        }
        Ok(())
    }
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
    /// `z_k = P⁻¹ r_k` for several vectors (default: one [`Preconditioner::apply`] each).
    ///
    /// # Errors
    /// Back-end failures.
    fn apply_batch(&self, rs: &[&[f64]], zs: &mut [&mut [f64]]) -> Result<(), RsError> {
        for (r, z) in rs.iter().zip(zs.iter_mut()) {
            self.apply(r, z)?;
        }
        Ok(())
    }
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

/// What a [`GmresRun`] needs next.
enum Request {
    /// `A v`.
    Op(Vec<f64>),
    /// `P⁻¹ v`.
    Prec(Vec<f64>),
    /// Finished.
    Done,
}

/// Where a [`GmresRun`] waits for the product it requested.
#[derive(Clone, Copy)]
enum Pending {
    /// Left: `P⁻¹ b`, the reference norm.
    RefNorm,
    /// `A x` for the explicit residual of a cycle.
    Residual,
    /// Left: `P⁻¹ r`, the measure of a cycle.
    Measure,
    /// The first factor of an Arnoldi step (`P⁻¹ v_j` right, `A v_j` left).
    ArnoldiFirst,
    /// The second factor (`A P⁻¹ v_j` right, `P⁻¹ A v_j` left).
    ArnoldiSecond,
    /// Right: `P⁻¹ u`, the update of `x` at the end of a cycle.
    Update,
    /// Finished.
    Done,
}

/// One restarted GMRES solve as a resumable computation (module docs): it asks for one operator
/// or preconditioner product at a time, so that several solves can share batched products
/// ([`gmres_batch`]). The arithmetic of a run does not depend on how its products are batched.
struct GmresRun<'a> {
    b: &'a [f64],
    options: GmresOptions,
    left: bool,
    n: usize,
    m: usize,
    t0: Instant,
    rep: GmresReport,
    x: Vec<f64>,
    bnorm: f64,
    ref_norm: f64,
    true_rel: f64,
    basis: Vec<Vec<f64>>,
    h: Vec<Vec<f64>>,
    cs: Vec<f64>,
    sn: Vec<f64>,
    g: Vec<f64>,
    j: usize,
    pending: Pending,
    error: Option<RsError>,
}

impl<'a> GmresRun<'a> {
    /// Starts a solve of dimension `n`; returns the run and its first request.
    fn start(n: usize, b: &'a [f64], options: &GmresOptions) -> Result<(GmresRun<'a>, Request), RsError> {
        options.validate()?;
        let t0 = Instant::now();
        let mut run = GmresRun {
            b,
            options: *options,
            left: options.side == GmresSide::Left,
            n,
            m: options.restart.min(n.max(1)),
            t0,
            rep: GmresReport {
                side: options.side,
                ..Default::default()
            },
            x: vec![0.0; n],
            bnorm: norm2(b),
            ref_norm: 0.0,
            true_rel: 0.0,
            basis: Vec::new(),
            h: Vec::new(),
            cs: Vec::new(),
            sn: Vec::new(),
            g: Vec::new(),
            j: 0,
            pending: Pending::Done,
            error: None,
        };
        if run.bnorm == 0.0 {
            run.rep.converged = true;
            run.rep.preconditioned_relative_residual = run.left.then_some(0.0);
            return Ok((run, Request::Done));
        }
        let m = run.m;
        run.basis = Vec::with_capacity(m + 1);
        run.h = Vec::with_capacity(m);
        run.cs = vec![0.0; m];
        run.sn = vec![0.0; m];
        run.g = vec![0.0; m + 1];
        if run.left {
            run.rep.preconditioner_applications += 1;
            run.pending = Pending::RefNorm;
            return Ok((run, Request::Prec(b.to_vec())));
        }
        run.ref_norm = run.bnorm;
        let req = run.after_ref_norm()?;
        Ok((run, req))
    }

    fn after_ref_norm(&mut self) -> Result<Request, RsError> {
        if !(self.ref_norm > 0.0 && self.ref_norm.is_finite()) {
            return Err(RsError::SolverDiverged {
                residual: f64::NAN,
            });
        }
        Ok(self.cycle_start())
    }

    /// The explicit residual `r = b − A x` of a cycle (`x = 0` in the first).
    fn cycle_start(&mut self) -> Request {
        if self.rep.cycles > 0 {
            self.rep.operator_applications += 1;
            self.pending = Pending::Residual;
            return Request::Op(self.x.clone());
        }
        let r = self.b.to_vec();
        self.after_residual(r)
    }

    fn after_residual(&mut self, r: Vec<f64>) -> Request {
        self.true_rel = norm2(&r) / self.bnorm;
        if self.left {
            self.rep.preconditioner_applications += 1;
            self.pending = Pending::Measure;
            return Request::Prec(r);
        }
        self.after_measure(r)
    }

    fn after_measure(&mut self, measure: Vec<f64>) -> Request {
        let beta = norm2(&measure);
        let rel = beta / self.ref_norm;
        self.rep.cycle_residuals.push(rel);
        if !rel.is_finite() {
            return self.done();
        }
        if rel <= self.options.tolerance {
            self.rep.converged = true;
            return self.done();
        }
        if self.rep.iterations >= self.options.max_iterations {
            return self.done();
        }
        self.rep.cycles += 1;
        self.basis.clear();
        self.h.clear();
        self.basis.push(measure.iter().map(|v| v / beta).collect());
        self.g.iter_mut().for_each(|v| *v = 0.0);
        self.g[0] = beta;
        self.j = 0;
        self.arnoldi_first()
    }

    /// `w = A P⁻¹ v_j` (right) or `P⁻¹ A v_j` (left), first factor.
    fn arnoldi_first(&mut self) -> Request {
        self.rep.operator_applications += 1;
        self.rep.preconditioner_applications += 1;
        self.pending = Pending::ArnoldiFirst;
        let v = self.basis[self.j].clone();
        if self.left { Request::Op(v) } else { Request::Prec(v) }
    }

    /// Modified Gram–Schmidt, the Givens rotations, the stop test of one Arnoldi step.
    fn after_arnoldi(&mut self, mut w: Vec<f64>) -> Request {
        let j = self.j;
        let mut col = vec![0.0; j + 2];
        for i in 0..=j {
            let hij = dot(&w, &self.basis[i]);
            col[i] = hij;
            for (wk, vk) in w.iter_mut().zip(&self.basis[i]) {
                *wk -= hij * vk;
            }
        }
        let hn = norm2(&w);
        col[j + 1] = hn;
        let (cs, sn) = (&mut self.cs, &mut self.sn);
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
        self.g[j + 1] = -sn[j] * self.g[j];
        self.g[j] *= cs[j];
        self.h.push(col);
        self.rep.iterations += 1;
        let est = self.g[j + 1].abs() / self.ref_norm;
        self.rep.history.push(est);
        if hn == 0.0
            || !hn.is_finite()
            || est <= self.options.tolerance
            || self.rep.iterations >= self.options.max_iterations
            || j + 1 >= self.m
        {
            return self.cycle_end(j + 1);
        }
        self.basis.push(w.iter().map(|v| v / hn).collect());
        self.j = j + 1;
        self.arnoldi_first()
    }

    /// Back substitution `H y = g` (`k × k`), `u = V y`, and the update of `x`.
    fn cycle_end(&mut self, k: usize) -> Request {
        let mut y = vec![0.0; k];
        for i in (0..k).rev() {
            let mut s = self.g[i];
            for l in i + 1..k {
                s -= self.h[l][i] * y[l];
            }
            y[i] = if self.h[i][i] != 0.0 { s / self.h[i][i] } else { 0.0 };
        }
        let mut u = vec![0.0; self.n];
        for (l, yl) in y.iter().enumerate() {
            for (ui, vi) in u.iter_mut().zip(&self.basis[l]) {
                *ui += yl * vi;
            }
        }
        if self.left {
            for (xi, ui) in self.x.iter_mut().zip(&u) {
                *xi += ui;
            }
            return self.cycle_start();
        }
        self.rep.preconditioner_applications += 1;
        self.pending = Pending::Update;
        Request::Prec(u)
    }

    fn done(&mut self) -> Request {
        self.pending = Pending::Done;
        Request::Done
    }

    /// Continues with the product that was requested.
    fn advance(&mut self, out: Vec<f64>) -> Result<Request, RsError> {
        Ok(match self.pending {
            Pending::RefNorm => {
                self.ref_norm = norm2(&out);
                return self.after_ref_norm();
            }
            Pending::Residual => {
                let mut r = self.b.to_vec();
                for (ri, ai) in r.iter_mut().zip(&out) {
                    *ri -= ai;
                }
                self.after_residual(r)
            }
            Pending::Measure => self.after_measure(out),
            Pending::ArnoldiFirst => {
                self.pending = Pending::ArnoldiSecond;
                if self.left { Request::Prec(out) } else { Request::Op(out) }
            }
            Pending::ArnoldiSecond => self.after_arnoldi(out),
            Pending::Update => {
                for (xi, ti) in self.x.iter_mut().zip(&out) {
                    *xi += ti;
                }
                self.cycle_start()
            }
            Pending::Done => Request::Done,
        })
    }

    fn finish(mut self) -> Result<(Vec<f64>, GmresReport), RsError> {
        if let Some(e) = self.error {
            return Err(e);
        }
        if self.bnorm == 0.0 {
            self.rep.seconds = self.t0.elapsed().as_secs_f64();
            return Ok((self.x, self.rep));
        }
        self.rep.relative_residual = self.true_rel;
        if self.left {
            self.rep.preconditioned_relative_residual = self.rep.cycle_residuals.last().copied();
        }
        self.rep.seconds = self.t0.elapsed().as_secs_f64();
        if !self.rep.converged && !(self.true_rel <= DIVERGENCE_THRESHOLD) {
            return Err(RsError::SolverDiverged {
                residual: self.true_rel,
            });
        }
        Ok((self.x, self.rep))
    }
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
    let mut out = gmres_batch(op, prec, &[b], options)?;
    out.pop().expect("one solve")
}

/// [`gmres`] for several right-hand sides at once: independent solves (each exactly as
/// [`gmres`] would do it) that advance in lockstep, so that the products they need at the same
/// time go through [`LinearOperator::apply_batch`] and [`Preconditioner::apply_batch`] together
/// (one pass over the hierarchical matrix for all of them, App. B). One result per right-hand
/// side.
///
/// # Errors
/// [`RsError::InvalidOption`] for invalid options or mismatched dimensions (the whole batch);
/// operator and preconditioner errors (the whole batch). A solve that diverges is an `Err` in its
/// own slot.
#[allow(clippy::type_complexity)]
pub fn gmres_batch(
    op: &dyn LinearOperator,
    prec: &dyn Preconditioner,
    bs: &[&[f64]],
    options: &GmresOptions,
) -> Result<Vec<Result<(Vec<f64>, GmresReport), RsError>>, RsError> {
    options.validate()?;
    let n = op.dim();
    for b in bs {
        if prec.dim() != n || b.len() != n {
            return Err(RsError::InvalidOption {
                name: "gmres dimensions",
                value: format!("operator {n}, preconditioner {}, rhs {}", prec.dim(), b.len()),
                requirement: "equal dimensions",
            });
        }
    }
    let mut runs = Vec::with_capacity(bs.len());
    let mut reqs = Vec::with_capacity(bs.len());
    for b in bs {
        let (r, q) = GmresRun::start(n, b, options)?;
        runs.push(r);
        reqs.push(q);
    }
    let mut ys: Vec<Vec<f64>> = Vec::new();
    loop {
        let mut progressed = false;
        for want_op in [false, true] {
            let idx: Vec<usize> = (0..runs.len())
                .filter(|&i| match &reqs[i] {
                    Request::Op(_) => want_op,
                    Request::Prec(_) => !want_op,
                    Request::Done => false,
                })
                .collect();
            if idx.is_empty() {
                continue;
            }
            progressed = true;
            ys.resize_with(idx.len(), || vec![0.0; n]);
            {
                let xs: Vec<&[f64]> = idx
                    .iter()
                    .map(|&i| match &reqs[i] {
                        Request::Op(v) | Request::Prec(v) => v.as_slice(),
                        Request::Done => unreachable!(),
                    })
                    .collect();
                let mut outs: Vec<&mut [f64]> = ys.iter_mut().map(|y| y.as_mut_slice()).collect();
                if want_op {
                    op.apply_batch(&xs, &mut outs)?;
                } else {
                    prec.apply_batch(&xs, &mut outs)?;
                }
            }
            for (k, &i) in idx.iter().enumerate() {
                let out = std::mem::replace(&mut ys[k], vec![0.0; n]);
                reqs[i] = match runs[i].advance(out) {
                    Ok(q) => q,
                    Err(e) => {
                        runs[i].error = Some(e);
                        Request::Done
                    }
                };
            }
        }
        if !progressed {
            break;
        }
    }
    Ok(runs.into_iter().map(GmresRun::finish).collect())
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

    fn solve_many(&mut self, rhs: &[(&[f64], &[f64])]) -> Result<Vec<SaddleSolution>, RsError> {
        let mut full = Vec::with_capacity(rhs.len());
        for (a, b) in rhs {
            check_rhs(self, a, b)?;
            full.push(a.iter().chain(b.iter()).copied().collect::<Vec<f64>>());
        }
        let refs: Vec<&[f64]> = full.iter().map(Vec::as_slice).collect();
        let outs = gmres_batch(&self.op, &self.prec, &refs, &self.options)?;
        outs.into_iter()
            .map(|r| {
                let (y, rep) = r?;
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
            })
            .collect()
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

    fn apply_batch(&self, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        (**self).apply_batch(xs, ys)
    }
}

impl<T: Preconditioner + ?Sized> Preconditioner for &T {
    fn dim(&self) -> usize {
        (**self).dim()
    }

    fn apply(&self, r: &[f64], z: &mut [f64]) -> Result<(), RsError> {
        (**self).apply(r, z)
    }

    fn apply_batch(&self, rs: &[&[f64]], zs: &mut [&mut [f64]]) -> Result<(), RsError> {
        (**self).apply_batch(rs, zs)
    }
}
