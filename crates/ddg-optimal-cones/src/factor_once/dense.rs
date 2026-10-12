#![allow(clippy::needless_range_loop)]
//! The dense active-set method for the reduced problem on a working set of cone sites.
//!
//! ```text
//! minimize   ½ µᵀ H µ + qᵀ µ + Σ_j t_j |µ_j|      subject to   lo_j ≤ µ_j ≤ hi_j,   [Σ_j µ_j = c]
//! ```
//!
//! with `H` positive semidefinite and `t_j ≥ 0`. With the sum constraint, variable 0 is the
//! *pin*: its row and column of `H` and its `q` are zero (it does not influence `u`), and it
//! absorbs the sum.
//!
//! This is the classical method for ℓ1 problems (Osborne, Presnell, Turlach 2000) extended to
//! box constraints. Every variable is in one of three states: *zero*, *free* (strictly between
//! its breakpoints, with a fixed sign label) or *clamped* at a nonzero bound. The free
//! variables `B` are determined by the stationarity system with a Cholesky factor of `H_BB`
//! that is extended and shrunk in `O(a²)`. The iteration walks from the current point towards
//! the solution of that system and stops at the first breakpoint it meets (a free variable
//! reaching zero becomes zero, one reaching a bound becomes clamped); when the solution is
//! reached, the optimality conditions of all other variables are checked and violators enter
//! (zero ones become free, clamped ones are released). The objective decreases at every step.

/// Counters of the dense iteration.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct DenseStats {
    /// Iterations (one linear system each).
    pub iterations: usize,
    /// Variables that entered the free set.
    pub adds: usize,
    /// Variables that left it (to zero or to a bound).
    pub drops: usize,
    /// Entries that were refused because `H_BB` would become singular.
    pub failed_adds: usize,
}

const ZERO: u8 = 0;
const FREE: u8 = 1;
const CLAMPED: u8 = 2;

pub(super) struct Dense {
    sum: Option<f64>,
    /// Starting value of the pin.
    pin_start: f64,
    dim: usize,
    cap: usize,
    h: Vec<f64>,
    q: Vec<f64>,
    pub mu: Vec<f64>,
    sign: Vec<f64>,
    state: Vec<u8>,
    t: Vec<f64>,
    lo: Vec<f64>,
    hi: Vec<f64>,
    blocked: Vec<bool>,
    /// Free variables except the pin, in the order of the Cholesky factor.
    act: Vec<usize>,
    clamped: Vec<usize>,
    chol: Vec<f64>,
    chol_ld: usize,
    /// The multiplier of the sum constraint at the last solution (0 without the constraint).
    pub nu: f64,
    /// Up to this fraction of the nonzero variables may enter per iteration (at least one).
    pub multi_add: f64,
    pub stats: DenseStats,
}

/// Whether a variable has a kink or a bound at zero: then its sign label restricts it to one
/// side of zero. A variable with `t = 0` and `lo < 0 < hi` passes through zero freely.
#[inline]
fn kinked(t: f64, lo: f64, hi: f64) -> bool {
    t > 0.0 || lo == 0.0 || hi == 0.0
}

/// How much the optimality condition of a variable that is at zero is violated, and the sign it
/// would take: with `v = ∂(distortion)/∂µ_j + ν`, positive values need `v < −t` (and
/// `hi > 0`), negative ones `v > t` (and `lo < 0`).
///
/// [formula: reduced-optimality-conditions]
pub(super) fn violation_at_zero(t: f64, lo: f64, hi: f64, v: f64) -> Option<(f64, f64)> {
    let pos = if hi > 0.0 { -v - t } else { f64::NEG_INFINITY };
    let neg = if lo < 0.0 { v - t } else { f64::NEG_INFINITY };
    if pos >= neg {
        pos.is_finite().then_some((pos, 1.0))
    } else {
        neg.is_finite().then_some((neg, -1.0))
    }
}

impl Dense {
    /// `sum`: `Some(c)` for the constraint `Σ µ = c`, then the first variable added by
    /// [`Dense::extend`] is the pin, which starts free with the value `pin_start`.
    pub fn new(sum: Option<f64>, pin_start: f64) -> Dense {
        Dense {
            sum,
            pin_start,
            dim: 0,
            cap: 0,
            h: vec![],
            q: vec![],
            mu: vec![],
            sign: vec![],
            state: vec![],
            t: vec![],
            lo: vec![],
            hi: vec![],
            blocked: vec![],
            act: vec![],
            clamped: vec![],
            chol: vec![],
            chol_ld: 0,
            nu: 0.0,
            multi_add: 0.1,
            stats: DenseStats::default(),
        }
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Bytes held by the dense matrices at the capacity that [`Dense::extend`] will allocate
    /// for `dim` variables (the matrix `H` and the Cholesky factor).
    pub fn bytes_for(dim: usize) -> usize {
        let cap = Self::capacity_for(dim);
        8 * 2 * cap * cap
    }

    fn capacity_for(dim: usize) -> usize {
        (dim * 3 / 2).max(64)
    }

    fn pin_free(&self) -> bool {
        self.sum.is_some() && self.state[0] == FREE
    }

    /// Variables whose value is not zero (free, clamped, and the pin if it is free).
    pub fn nonzeros(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.act.iter().chain(&self.clamped).copied().collect();
        if self.pin_free() {
            v.push(0);
        }
        v
    }

    pub fn n_nonzero(&self) -> usize {
        self.act.len() + self.clamped.len() + usize::from(self.pin_free())
    }

    /// Appends `k` variables. `cross` is `(dim + k) × k` column major: column `s` holds
    /// `H[·, dim + s]` for all rows of the enlarged matrix (rows `dim..` are the new block).
    /// `t`, `lo`, `hi` describe the new variables.
    pub fn extend(
        &mut self,
        cross: &[f64],
        k: usize,
        q_new: &[f64],
        t: &[f64],
        lo: &[f64],
        hi: &[f64],
    ) {
        let old = self.dim;
        let new = old + k;
        if new > self.cap {
            let cap = Self::capacity_for(new);
            let mut h = vec![0.0; cap * cap];
            for j in 0..old {
                h[j * cap..j * cap + old]
                    .copy_from_slice(&self.h[j * self.cap..j * self.cap + old]);
            }
            self.h = h;
            self.cap = cap;
        }
        for s in 0..k {
            for i in 0..new {
                let v = if cross.is_empty() {
                    0.0
                } else {
                    cross[s * new + i]
                };
                self.h[(old + s) * self.cap + i] = v;
                self.h[i * self.cap + old + s] = v;
            }
        }
        self.q.extend_from_slice(q_new);
        self.mu.extend(std::iter::repeat_n(0.0, k));
        self.sign.extend(std::iter::repeat_n(0.0, k));
        self.state.extend(std::iter::repeat_n(ZERO, k));
        self.t.extend_from_slice(t);
        self.lo.extend_from_slice(lo);
        self.hi.extend_from_slice(hi);
        self.blocked.extend(std::iter::repeat_n(false, k));
        self.dim = new;
        if let (0, Some(_)) = (old, self.sum) {
            // The pin starts free with the value the caller chose (within its bounds).
            let c = self.pin_start;
            self.mu[0] = c;
            self.sign[0] = if c > 0.0 || (c == 0.0 && self.hi[0] > 0.0) {
                1.0
            } else {
                -1.0
            };
            self.state[0] = FREE;
        }
    }

    /// Starts variable `j` clamped at the nonzero bound `bound` (a feasible starting point for a
    /// sum that no single variable can carry).
    pub fn clamp(&mut self, j: usize, bound: f64) {
        self.state[j] = CLAMPED;
        self.mu[j] = bound;
        self.clamped.push(j);
    }

    #[inline]
    fn hh(&self, i: usize, j: usize) -> f64 {
        self.h[j * self.cap + i]
    }

    /// The interval a free variable may move in with its current sign label.
    #[inline]
    fn interval(&self, j: usize) -> (f64, f64) {
        if kinked(self.t[j], self.lo[j], self.hi[j]) {
            if self.sign[j] > 0.0 {
                (0.0, self.hi[j])
            } else {
                (self.lo[j], 0.0)
            }
        } else {
            (self.lo[j], self.hi[j])
        }
    }

    fn chol_append(&mut self, j: usize) -> bool {
        let a = self.act.len();
        if a + 1 > self.chol_ld {
            let ld = ((a + 1) * 3 / 2).max(64);
            let mut c = vec![0.0; ld * ld];
            for i in 0..a {
                c[i * ld..i * ld + i + 1]
                    .copy_from_slice(&self.chol[i * self.chol_ld..i * self.chol_ld + i + 1]);
            }
            self.chol = c;
            self.chol_ld = ld;
        }
        let ld = self.chol_ld;
        let mut y = vec![0.0; a];
        for i in 0..a {
            let mut s = self.hh(self.act[i], j);
            for t in 0..i {
                s -= self.chol[i * ld + t] * y[t];
            }
            y[i] = s / self.chol[i * ld + i];
        }
        let d2 = self.hh(j, j) - y.iter().map(|v| v * v).sum::<f64>();
        if d2.is_nan() || d2 <= 1e-13 * self.hh(j, j).abs() || !d2.is_finite() {
            return false;
        }
        self.chol[a * ld..a * ld + a].copy_from_slice(&y);
        self.chol[a * ld + a] = d2.sqrt();
        self.act.push(j);
        true
    }

    fn chol_delete(&mut self, k: usize) {
        let a = self.act.len();
        let ld = self.chol_ld;
        for i in k..a - 1 {
            for t in 0..=i + 1 {
                self.chol[i * ld + t] = self.chol[(i + 1) * ld + t];
            }
        }
        for i in k..a - 1 {
            let x = self.chol[i * ld + i];
            let y = self.chol[i * ld + i + 1];
            let r = x.hypot(y);
            let (c, s) = (x / r, y / r);
            for t in i..a - 1 {
                let xt = self.chol[t * ld + i];
                let yt = self.chol[t * ld + i + 1];
                self.chol[t * ld + i] = c * xt + s * yt;
                self.chol[t * ld + i + 1] = -s * xt + c * yt;
            }
            self.chol[i * ld + i + 1] = 0.0;
        }
        for t in 0..a {
            self.chol[(a - 1) * ld + t] = 0.0;
        }
        for t in 0..a - 1 {
            self.chol[t * ld + a - 1] = 0.0;
        }
        self.act.remove(k);
    }

    fn chol_solve(&self, b: &[f64]) -> Vec<f64> {
        let a = self.act.len();
        let ld = self.chol_ld;
        let mut x = b.to_vec();
        for i in 0..a {
            let mut s = x[i];
            for t in 0..i {
                s -= self.chol[i * ld + t] * x[t];
            }
            x[i] = s / self.chol[i * ld + i];
        }
        for i in (0..a).rev() {
            let mut s = x[i];
            for t in i + 1..a {
                s -= self.chol[t * ld + i] * x[t];
            }
            x[i] = s / self.chol[i * ld + i];
        }
        x
    }

    /// Takes a free variable out of the free set, to zero or clamped at `bound`.
    fn leave(&mut self, j: usize, to: Option<f64>) {
        if self.sum.is_some() && j == 0 {
            // the pin has no row in the factor
        } else {
            let pos = self
                .act
                .iter()
                .position(|&x| x == j)
                .expect("free variable");
            self.chol_delete(pos);
        }
        match to {
            None => {
                self.state[j] = ZERO;
                self.sign[j] = 0.0;
                self.mu[j] = 0.0;
            }
            Some(b) => {
                self.state[j] = CLAMPED;
                self.mu[j] = b;
                self.clamped.push(j);
            }
        }
    }

    /// Re-anchors the linear term so that the dense gradient `q + H µ` equals the exact
    /// gradient `exact` (one entry per variable) at the current point. The Gram matrix of a
    /// working set with clustered cones is badly conditioned, and the dense solution inherits an
    /// error of the order of `cond(H) · ε` in its stationarity conditions; the exact gradient,
    /// from solves with the Laplacian factor, does not. Solving again with the anchored `q` is
    /// iterative refinement with the Cholesky factor of `H_BB` as the approximate inverse.
    pub fn anchor_gradient(&mut self, exact: &[f64]) {
        let mut g = self.q.clone();
        for i in self.act.iter().chain(&self.clamped) {
            let m = self.mu[*i];
            if m != 0.0 {
                for (j, gj) in g.iter_mut().enumerate() {
                    *gj += self.hh(j, *i) * m;
                }
            }
        }
        for j in 0..self.dim {
            if self.sum.is_some() && j == 0 {
                continue;
            }
            self.q[j] += exact[j] - g[j];
        }
    }

    /// How much the optimality condition of variable `j` is violated, given the exact value
    /// `v = ∂(distortion)/∂µ_j + ν` (free variables satisfy theirs by construction): at zero as
    /// [`violation_at_zero`]; clamped at `hi > 0` the slope `v + t` must be `≤ 0`, clamped at
    /// `lo < 0` the slope `v − t` must be `≥ 0`.
    ///
    /// [formula: reduced-optimality-conditions]
    pub fn exact_violation(&self, j: usize, v: f64) -> Option<(f64, f64)> {
        match self.state[j] {
            FREE => None,
            ZERO => violation_at_zero(self.t[j], self.lo[j], self.hi[j], v),
            _ if self.mu[j] > 0.0 => Some((v + self.t[j], 1.0)),
            _ => Some((self.t[j] - v, -1.0)),
        }
    }

    /// Runs the active-set iteration until no variable violates its optimality condition by
    /// more than `tol`. Returns false if `max_iter` iterations did not suffice.
    pub fn solve(&mut self, tol: f64, max_iter: usize) -> bool {
        let mut unblocks = 0;
        for _ in 0..max_iter {
            self.stats.iterations += 1;
            let b: Vec<usize> = self.act.clone();
            let mut rhs: Vec<f64> = b
                .iter()
                .map(|&j| -self.q[j] - self.t[j] * self.sign[j])
                .collect();
            for &c in &self.clamped {
                let m = self.mu[c];
                for (r, &j) in rhs.iter_mut().zip(&b) {
                    *r -= self.hh(j, c) * m;
                }
            }
            let sum_clamped: f64 = self.clamped.iter().map(|&c| self.mu[c]).sum();
            let pin_free = self.pin_free();
            let mut star = vec![0.0; self.dim];
            let nu;
            match self.sum {
                None => {
                    let x = self.chol_solve(&rhs);
                    for (s, &j) in b.iter().enumerate() {
                        star[j] = x[s];
                    }
                    nu = 0.0;
                }
                Some(c) if pin_free => {
                    nu = -self.t[0] * self.sign[0];
                    let r2: Vec<f64> = rhs.iter().map(|v| v - nu).collect();
                    let x = self.chol_solve(&r2);
                    let mut sum = 0.0;
                    for (s, &j) in b.iter().enumerate() {
                        star[j] = x[s];
                        sum += x[s];
                    }
                    star[0] = c - sum - sum_clamped;
                }
                Some(c) => {
                    // the pin is at zero or clamped: the sum constraint binds the free variables
                    let pin_value = if self.state[0] == CLAMPED {
                        self.mu[0]
                    } else {
                        0.0
                    };
                    if b.is_empty() {
                        nu = 0.0;
                    } else {
                        let x = self.chol_solve(&rhs);
                        let y = self.chol_solve(&vec![1.0; b.len()]);
                        let sx: f64 = x.iter().sum();
                        let sy: f64 = y.iter().sum();
                        // the pin is clamped: its value is part of `sum_clamped` already
                        let target = c
                            - sum_clamped
                            - if self.state[0] == CLAMPED {
                                0.0
                            } else {
                                pin_value
                            };
                        nu = (sx - target) / sy;
                        for (s, &j) in b.iter().enumerate() {
                            star[j] = x[s] - nu * y[s];
                        }
                    }
                }
            }
            // line search over the free variables (and the free pin)
            let mut actl = b.clone();
            if pin_free {
                actl.push(0);
            }
            let mut t_step = 1.0f64;
            let mut blocking: Option<(usize, f64, bool)> = None; // variable, bound, is_zero
            for &i in &actl {
                let (a, bnd) = self.interval(i);
                let (x, m) = (star[i], self.mu[i]);
                if x > bnd {
                    let ti = ((bnd - m) / (x - m)).max(0.0);
                    if ti < t_step {
                        t_step = ti;
                        blocking = Some((
                            i,
                            bnd,
                            bnd == 0.0 && kinked(self.t[i], self.lo[i], self.hi[i]),
                        ));
                    }
                } else if x < a {
                    let ti = ((a - m) / (x - m)).max(0.0);
                    if ti < t_step {
                        t_step = ti;
                        blocking =
                            Some((i, a, a == 0.0 && kinked(self.t[i], self.lo[i], self.hi[i])));
                    }
                }
            }
            if let Some((d, bound, is_zero)) = blocking {
                for &i in &actl {
                    self.mu[i] += t_step * (star[i] - self.mu[i]);
                }
                self.leave(d, if is_zero { None } else { Some(bound) });
                self.stats.drops += 1;
                if t_step < 1e-12 {
                    self.blocked[d] = true;
                } else {
                    self.blocked.iter_mut().for_each(|x| *x = false);
                    unblocks = 0;
                }
                continue;
            }
            for &i in &actl {
                self.mu[i] = star[i];
            }
            self.nu = nu;
            // optimality conditions of everything else
            let mut g = self.q.clone();
            for &i in actl.iter().chain(&self.clamped) {
                let m = self.mu[i];
                if m != 0.0 {
                    for (j, gj) in g.iter_mut().enumerate() {
                        *gj += self.hh(j, i) * m;
                    }
                }
            }
            let mut viol: Vec<(usize, f64, f64)> = Vec::new();
            let mut viol_blocked: Vec<(usize, f64, f64)> = Vec::new();
            for j in 0..self.dim {
                if self.state[j] == FREE {
                    continue;
                }
                let v = g[j] + nu;
                let found = if self.state[j] == ZERO {
                    violation_at_zero(self.t[j], self.lo[j], self.hi[j], v)
                } else if self.mu[j] > 0.0 {
                    // clamped at the upper bound: it may move down if the slope is positive
                    Some((v + self.t[j], 1.0))
                } else {
                    Some((self.t[j] - v, -1.0))
                };
                if let Some((ex, s)) = found.filter(|f| f.0 > tol) {
                    if self.blocked[j] {
                        viol_blocked.push((j, ex, s));
                    } else {
                        viol.push((j, ex, s));
                    }
                }
            }
            if viol.is_empty() && !viol_blocked.is_empty() && unblocks < 50 {
                // Everything else is settled: a variable that was blocked after a zero-length
                // step may enter now (a limited number of times in a row, against cycling).
                unblocks += 1;
                self.blocked.iter_mut().for_each(|x| *x = false);
                viol = viol_blocked;
            }
            if viol.is_empty() {
                return true;
            }
            viol.sort_by(|a, b| b.1.total_cmp(&a.1));
            let m = ((self.multi_add * self.n_nonzero() as f64) as usize).max(1);
            for &(j, _, s) in viol.iter().take(m) {
                let was_clamped = self.state[j] == CLAMPED;
                if self.sum.is_some() && j == 0 {
                    if was_clamped {
                        self.clamped.retain(|&x| x != 0);
                    }
                    self.state[0] = FREE;
                    self.sign[0] = s;
                    self.stats.adds += 1;
                } else if self.chol_append(j) {
                    if was_clamped {
                        self.clamped.retain(|&x| x != j);
                    }
                    self.state[j] = FREE;
                    self.sign[j] = s;
                    self.stats.adds += 1;
                } else {
                    self.blocked[j] = true;
                    self.stats.failed_adds += 1;
                }
            }
        }
        false
    }
}
