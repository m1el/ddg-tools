//! The generic optimizer: L-BFGS / gradient descent / fixed step on an [`Objective`].

use std::collections::VecDeque;
use std::fmt;
use std::time::Instant;

use super::options::{
    Control, FlatProgress, InvalidOption, LineSearch, Method, OptimizeOptions, OptimizeReport,
    OptimizeStatus,
};

/// Failure of one objective evaluation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvalError<E> {
    /// The energy is undefined at this point (SPEC D12): `+∞` inside a line search, an error at
    /// the starting point.
    Undefined(E),
    /// The point is outside the admissible set (D31/D32): handled like `Undefined` but counted in
    /// `inadmissible_trials`.
    Inadmissible(E),
    /// A real failure: aborts the run.
    Fatal(E),
}

/// A function to minimize, with its (sub)gradient.
pub trait Objective {
    /// Error type carried by [`EvalError`].
    type Error;
    /// Writes the gradient at `x` into `gradient` (same length as `x`) and returns the energy.
    /// A non-finite energy or gradient is treated like [`EvalError::Undefined`].
    fn evaluate(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64, EvalError<Self::Error>>;
    /// Marks the groups (vertices) that must not move from `x`: `blocked[i]` is set for group `i`.
    /// Their gradient entries are treated as zero and their direction entries are zeroed, as if
    /// they were pinned, until the next accepted step (D31: the vertices of a face at the
    /// admissibility floor, whose gradient would push the face through it). Default: none.
    fn blocked(&mut self, _x: &[f64], _blocked: &mut [bool]) {}
    /// Number of near ties of the last evaluation (reporting only).
    fn near_ties(&self) -> usize {
        0
    }
}

/// Error type of [`ClosureObjective`]: the closure returned `None`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Undefined;

/// Adapter: a closure `(x, gradient) -> Option<energy>` (`None` = undefined).
pub struct ClosureObjective<F>(pub F);

impl<F: FnMut(&[f64], &mut [f64]) -> Option<f64>> Objective for ClosureObjective<F> {
    type Error = Undefined;
    fn evaluate(&mut self, x: &[f64], g: &mut [f64]) -> Result<f64, EvalError<Undefined>> {
        (self.0)(x, g).ok_or(EvalError::Undefined(Undefined))
    }
}

/// Units and layout of the variables.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Problem {
    /// Length unit `ℓ̄ > 0` (D17).
    pub length_scale: f64,
    /// Coordinates per vertex (3 for meshes); norms are per-vertex Euclidean maxima.
    pub group: usize,
}

/// Errors of [`minimize`].
#[derive(Debug, Clone, PartialEq)]
pub enum OptimizeError<E> {
    /// Invalid option.
    Options(InvalidOption),
    /// `x0.len()` is not a positive multiple of `group`, or `group == 0`.
    BadLayout {
        /// Length of `x0`.
        len: usize,
        /// Group size.
        group: usize,
    },
    /// `length_scale` is not finite and positive.
    BadLengthScale(f64),
    /// A pinned id is `>= n_groups`.
    PinnedOutOfRange {
        /// Offending id.
        id: u32,
        /// Number of groups.
        n_groups: usize,
    },
    /// The starting point has a non-finite energy or gradient.
    NonFiniteStart,
    /// The energy is undefined at the starting point.
    UndefinedAtStart(E),
    /// Fatal objective error.
    Objective(E),
}

impl<E: fmt::Display> fmt::Display for OptimizeError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OptimizeError::Options(e) => e.fmt(f),
            OptimizeError::BadLayout { len, group } => {
                write!(f, "{len} variables is not a positive multiple of the group size {group}")
            }
            OptimizeError::BadLengthScale(l) => write!(f, "length scale {l} must be finite and > 0"),
            OptimizeError::PinnedOutOfRange { id, n_groups } => {
                write!(f, "pinned id {id} out of range ({n_groups} vertices)")
            }
            OptimizeError::NonFiniteStart => {
                write!(f, "energy or gradient is not finite at the starting point")
            }
            OptimizeError::UndefinedAtStart(e) => write!(f, "energy undefined at the start: {e}"),
            OptimizeError::Objective(e) => e.fmt(f),
        }
    }
}

impl<E: fmt::Debug + fmt::Display> std::error::Error for OptimizeError<E> {}

/// Result of [`minimize`].
#[derive(Debug, Clone, PartialEq)]
pub struct Minimized {
    /// Last accepted variables.
    pub x: Vec<f64>,
    /// Report.
    pub report: OptimizeReport,
}

type Res<T, O> = Result<T, OptimizeError<<O as Objective>::Error>>;

struct Memory {
    cap: usize,
    pairs: VecDeque<(Vec<f64>, Vec<f64>, f64)>, // s, y, rho
}

impl Memory {
    fn new(cap: usize) -> Self {
        Memory { cap, pairs: VecDeque::new() }
    }
    fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }
    fn reset(&mut self) {
        self.pairs.clear();
    }
    fn push(&mut self, s: Vec<f64>, y: Vec<f64>, sy: f64) {
        if self.pairs.len() == self.cap {
            self.pairs.pop_front();
        }
        self.pairs.push_back((s, y, 1.0 / sy));
    }

    /// Two-loop recursion, `d = -H g` with `H_0 = (sᵀy / yᵀy) I` of the newest pair.
    ///
    /// [formula: dev-lbfgs-two-loop] (N&W Alg. 7.4)
    fn direction(&self, g: &[f64], d: &mut [f64]) {
        let k = self.pairs.len();
        let mut q = g.to_vec();
        let mut alpha = vec![0.0; k];
        for i in (0..k).rev() {
            let (s, y, rho) = &self.pairs[i];
            alpha[i] = rho * dot(s, &q);
            for (qj, yj) in q.iter_mut().zip(y) {
                *qj -= alpha[i] * yj;
            }
        }
        let (s, y, _) = &self.pairs[k - 1];
        let gamma = dot(s, y) / dot(y, y);
        for qj in q.iter_mut() {
            *qj *= gamma;
        }
        for ((s, y, rho), a) in self.pairs.iter().zip(&alpha) {
            let beta = rho * dot(y, &q);
            for (qj, sj) in q.iter_mut().zip(s) {
                *qj += (a - beta) * sj;
            }
        }
        for (dj, qj) in d.iter_mut().zip(&q) {
            *dj = -qj;
        }
    }
}

/// Zeroes the entries of the blocked groups.
fn mask_blocked(v: &mut [f64], blocked: &[bool], group: usize) {
    for (c, &b) in v.chunks_exact_mut(group).zip(blocked) {
        if b {
            c.fill(0.0);
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Max over groups of the Euclidean norm.
fn norm_inf(v: &[f64], group: usize) -> f64 {
    v.chunks_exact(group)
        .map(|c| c.iter().map(|x| x * x).sum::<f64>().sqrt())
        .fold(0.0, f64::max)
}

enum Ls {
    Accepted { t: f64, f: f64, capped: bool, weak: bool },
    Failed,
}

struct Run<'a, O: Objective> {
    obj: &'a mut O,
    group: usize,
    free: Vec<bool>,
    // scratch
    tx: Vec<f64>,
    tg: Vec<f64>, // gradient of the last trial
    ag: Vec<f64>, // gradient at the best Armijo-satisfying step so far
    rep: OptimizeReport,
}

impl<O: Objective> Run<'_, O> {
    /// Evaluates `x + t d` into `tx`, `tg`. `None` if undefined (counted).
    fn trial(&mut self, x: &[f64], d: &[f64], t: f64) -> Res<Option<f64>, O> {
        for ((o, xi), di) in self.tx.iter_mut().zip(x).zip(d) {
            *o = xi + t * di;
        }
        self.rep.line_search_trials += 1;
        let tx = std::mem::take(&mut self.tx);
        let r = self.eval_at(&tx);
        self.tx = tx;
        r
    }

    fn eval_at(&mut self, x: &[f64]) -> Res<Option<f64>, O> {
        self.rep.evaluations += 1;
        let r = self.obj.evaluate(x, &mut self.tg);
        match r {
            Ok(f) => {
                for (g, &fr) in self.tg.iter_mut().zip(&self.free) {
                    if !fr {
                        *g = 0.0;
                    }
                }
                if f.is_finite() && self.tg.iter().all(|v| v.is_finite()) {
                    Ok(Some(f))
                } else {
                    self.rep.undefined_trials += 1;
                    Ok(None)
                }
            }
            Err(EvalError::Undefined(_)) => {
                self.rep.undefined_trials += 1;
                Ok(None)
            }
            Err(EvalError::Inadmissible(_)) => {
                self.rep.inadmissible_trials += 1;
                Ok(None)
            }
            Err(EvalError::Fatal(e)) => Err(OptimizeError::Objective(e)),
        }
    }

    /// Lewis–Overton Algorithm 2.6. On acceptance the gradient is in `self.tg`.
    ///
    /// [formula: dev-lewis-overton-line-search]; existence/termination L&O Thm. 2.5, 2.7
    /// [formula: dev-lewis-overton-termination]; the curvature test gives `sᵀy > 0`
    /// [formula: dev-weak-wolfe-curvature].
    #[allow(clippy::too_many_arguments)]
    fn weak_wolfe(
        &mut self,
        x: &[f64],
        f0: f64,
        d: &[f64],
        s: f64,
        t_max: f64,
        t_start: f64,
        (c1, c2, max_doublings, max_bisections): (f64, f64, u32, u32),
    ) -> Res<Ls, O> {
        let mut alpha = 0.0_f64;
        let mut f_alpha = f0;
        let mut beta = f64::INFINITY;
        let mut t = t_start.min(t_max);
        let (mut nd, mut nb) = (0u32, 0u32);
        loop {
            let ft = self.trial(x, d, t)?;
            let s_ok = matches!(ft, Some(f) if f - f0 < c1 * s * t);
            if !s_ok {
                beta = t;
            } else {
                let f = ft.unwrap_or(f0);
                if dot(&self.tg, d) > c2 * s {
                    return Ok(Ls::Accepted { t, f, capped: false, weak: false });
                }
                alpha = t;
                f_alpha = f;
                self.ag.copy_from_slice(&self.tg);
            }
            if beta.is_finite() {
                if nb >= max_bisections {
                    break;
                }
                nb += 1;
                t = 0.5 * (alpha + beta);
            } else if alpha >= t_max {
                std::mem::swap(&mut self.tg, &mut self.ag);
                return Ok(Ls::Accepted { t: alpha, f: f_alpha, capped: true, weak: false });
            } else {
                if nd >= max_doublings {
                    break;
                }
                nd += 1;
                t = (2.0 * alpha).min(t_max);
            }
        }
        if alpha > 0.0 {
            std::mem::swap(&mut self.tg, &mut self.ag);
            Ok(Ls::Accepted { t: alpha, f: f_alpha, capped: false, weak: true })
        } else {
            Ok(Ls::Failed)
        }
    }

    /// Armijo backtracking (N&W Alg. 3.1). [formula: dev-armijo-backtracking]
    #[allow(clippy::too_many_arguments)]
    fn backtracking(
        &mut self,
        x: &[f64],
        f0: f64,
        d: &[f64],
        s: f64,
        t_max: f64,
        t_start: f64,
        (c1, shrink, max_shrinks): (f64, f64, u32),
    ) -> Res<Ls, O> {
        let mut t = t_start.min(t_max);
        for _ in 0..=max_shrinks {
            if let Some(f) = self.trial(x, d, t)?
                && f - f0 <= c1 * s * t
            {
                return Ok(Ls::Accepted { t, f, capped: t >= t_max, weak: false });
            }
            t *= shrink;
        }
        Ok(Ls::Failed)
    }
}

/// Minimizes `obj` from `x0`.
///
/// `opts.pinned` are group (vertex) ids. Nothing is rescaled; step sizes and tolerances are in
/// units of `problem.length_scale` (SPEC §4.1). The returned `x` is the last accepted iterate.
///
/// # Errors
/// Invalid options or layout; an undefined/non-finite start; a fatal objective error.
pub fn minimize<O: Objective>(
    obj: &mut O,
    x0: &[f64],
    problem: &Problem,
    opts: &OptimizeOptions,
    mut progress: Option<&mut dyn FnMut(&FlatProgress) -> Control>,
) -> Res<Minimized, O> {
    let t_begin = Instant::now();
    opts.validate().map_err(OptimizeError::Options)?;
    let n = x0.len();
    let group = problem.group;
    if group == 0 || n == 0 || !n.is_multiple_of(group) {
        return Err(OptimizeError::BadLayout { len: n, group });
    }
    let ell = problem.length_scale;
    if !(ell.is_finite() && ell > 0.0) {
        return Err(OptimizeError::BadLengthScale(ell));
    }
    let mut free = vec![true; n];
    for &p in &opts.pinned {
        if p as usize >= n / group {
            return Err(OptimizeError::PinnedOutOfRange { id: p, n_groups: n / group });
        }
        free[p as usize * group..(p as usize + 1) * group].fill(false);
    }

    let mut run = Run {
        obj,
        group,
        free,
        tx: vec![0.0; n],
        tg: vec![0.0; n],
        ag: vec![0.0; n],
        rep: OptimizeReport {
            status: OptimizeStatus::MaxIterations,
            iterations: 0,
            evaluations: 0,
            line_search_trials: 0,
            undefined_trials: 0,
            inadmissible_trials: 0,
            capped_steps: 0,
            weak_steps: 0,
            skipped_pairs: 0,
            memory_resets: 0,
            energy_initial: 0.0,
            energy_final: 0.0,
            energy_history: Vec::new(),
            gradient_norm_history: Vec::new(),
            mean_edge_length: ell,
            max_vertex_displacement: 0.0,
            energy_increases: 0,
            near_ties: 0,
            seconds: 0.0,
        },
    };

    let mut x = x0.to_vec();
    // Initial evaluation: undefined is an error here.
    run.rep.evaluations += 1;
    let mut f = match run.obj.evaluate(&x, &mut run.tg) {
        Ok(f) => f,
        Err(EvalError::Undefined(e) | EvalError::Inadmissible(e)) => return Err(OptimizeError::UndefinedAtStart(e)),
        Err(EvalError::Fatal(e)) => return Err(OptimizeError::Objective(e)),
    };
    for (g, &fr) in run.tg.iter_mut().zip(&run.free) {
        if !fr {
            *g = 0.0;
        }
    }
    if !f.is_finite() || run.tg.iter().any(|v| !v.is_finite()) {
        return Err(OptimizeError::NonFiniteStart);
    }
    let mut g = run.tg.clone();
    let mut blocked = vec![false; n / group];
    run.obj.blocked(&x, &mut blocked);
    let mut ge = g.clone();
    mask_blocked(&mut ge, &blocked, group);
    let mut gnorm = ell * norm_inf(&ge, group);
    run.rep.energy_initial = f;
    run.rep.energy_history.push(f);
    run.rep.gradient_norm_history.push(gnorm);

    let (lbfgs, mem_cap) = match opts.method {
        Method::Lbfgs { memory } => (true, memory),
        _ => (false, 1),
    };
    let mut mem = Memory::new(mem_cap);
    let mut d = vec![0.0; n];
    let mut t_prev = 1.0_f64;
    let mut stalled = 0usize;
    let mut retry_pending = false;

    let status = loop {
        if gnorm <= opts.gradient_tolerance {
            break OptimizeStatus::GradientTolerance;
        }
        if run.rep.iterations >= opts.max_iterations {
            break OptimizeStatus::MaxIterations;
        }
        let gmax = norm_inf(&ge, group);
        let f_old = f;

        // ---- fixed step (§4.5)
        if let Method::FixedStep { step } = opts.method {
            let c = step * ell * ell;
            let x_new: Vec<f64> = x.iter().zip(&ge).map(|(xi, gi)| xi - c * gi).collect();
            let r = run.eval_at(&x_new)?;
            run.rep.line_search_trials += 1;
            let Some(f_new) = r else {
                break OptimizeStatus::LineSearchFailed;
            };
            if f_new > f_old {
                run.rep.energy_increases += 1;
            }
            x = x_new;
            f = f_new;
            g.copy_from_slice(&run.tg);
            run.obj.blocked(&x, &mut blocked);
            ge.copy_from_slice(&g);
            mask_blocked(&mut ge, &blocked, group);
        } else {
            // ---- direction
            let mut steepest = true;
            let scaled = |g: &[f64], d: &mut [f64]| {
                let c = opts.initial_step * ell / gmax;
                for (di, gi) in d.iter_mut().zip(g) {
                    *di = -c * gi;
                }
            };
            if lbfgs && !mem.is_empty() {
                mem.direction(&ge, &mut d);
                mask_blocked(&mut d, &blocked, group);
                steepest = false;
                let s = dot(&ge, &d);
                if !(s < 0.0 && s.is_finite()) {
                    mem.reset();
                    run.rep.memory_resets += 1;
                    steepest = true;
                }
            }
            if steepest {
                scaled(&ge, &mut d);
            }
            let s = dot(&ge, &d);
            if !(s < 0.0 && s.is_finite()) {
                break OptimizeStatus::NonDescent;
            }
            let dmax = norm_inf(&d, group);
            let t_max = opts.max_displacement * ell / dmax;
            let t_start = if lbfgs || matches!(opts.line_search, LineSearch::WeakWolfe { .. }) {
                1.0
            } else {
                2.0 * t_prev
            };
            let ls = match opts.line_search {
                LineSearch::WeakWolfe { c1, c2, max_doublings, max_bisections } => {
                    run.weak_wolfe(&x, f, &d, s, t_max, t_start, (c1, c2, max_doublings, max_bisections))?
                }
                LineSearch::Backtracking { c1, shrink, max_shrinks } => {
                    run.backtracking(&x, f, &d, s, t_max, t_start, (c1, shrink, max_shrinks))?
                }
            };
            match ls {
                Ls::Failed => {
                    if !steepest && !retry_pending {
                        mem.reset();
                        run.rep.memory_resets += 1;
                        retry_pending = true;
                        continue;
                    }
                    break OptimizeStatus::LineSearchFailed;
                }
                Ls::Accepted { t, f: f_new, capped, weak } => {
                    retry_pending = false;
                    run.rep.capped_steps += usize::from(capped);
                    run.rep.weak_steps += usize::from(weak);
                    t_prev = t;
                    let mut s_vec = vec![0.0; n];
                    let mut y_vec = vec![0.0; n];
                    for i in 0..n {
                        s_vec[i] = t * d[i];
                        y_vec[i] = run.tg[i] - g[i];
                        x[i] += s_vec[i];
                    }
                    mask_blocked(&mut y_vec, &blocked, group);
                    f = f_new;
                    g.copy_from_slice(&run.tg);
                    run.obj.blocked(&x, &mut blocked);
                    ge.copy_from_slice(&g);
                    mask_blocked(&mut ge, &blocked, group);
                    if lbfgs {
                        let sy = dot(&s_vec, &y_vec);
                        if sy > 0.0 && sy.is_finite() {
                            mem.push(s_vec, y_vec, sy);
                        } else {
                            run.rep.skipped_pairs += 1;
                        }
                    }
                }
            }
        }

        // ---- accepted iterate
        run.rep.iterations += 1;
        gnorm = ell * norm_inf(&ge, group);
        run.rep.energy_history.push(f);
        run.rep.gradient_norm_history.push(gnorm);
        if gnorm <= opts.gradient_tolerance {
            break OptimizeStatus::GradientTolerance;
        }
        if f_old - f <= opts.energy_tolerance * f_old.max(1e-300) {
            stalled += 1;
        } else {
            stalled = 0;
        }
        if stalled >= opts.stall_iterations {
            break OptimizeStatus::EnergyStalled;
        }
        if let Some(cb) = progress.as_deref_mut() {
            let info = FlatProgress { iteration: run.rep.iterations, energy: f, gradient_norm: gnorm, x: &x };
            if cb(&info) == Control::Break {
                break OptimizeStatus::Cancelled;
            }
        }
    };

    let mut rep = run.rep;
    rep.status = status;
    rep.energy_final = f;
    rep.near_ties = run.obj.near_ties();
    rep.max_vertex_displacement = x
        .chunks_exact(group)
        .zip(x0.chunks_exact(group))
        .map(|(a, b)| a.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum::<f64>().sqrt())
        .fold(0.0, f64::max)
        / ell;
    rep.seconds = t_begin.elapsed().as_secs_f64();
    let _ = run.group;
    Ok(Minimized { x, report: rep })
}
