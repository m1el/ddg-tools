//! The two Newton iterations of the paper (Fig. 12): the regularized problems of Stage I and
//! the exact problem of Stage II, with the safeguards and the closed-mesh treatment the paper
//! does not have.
//!
//! Both iterations work on the unknowns `(u, φ, µ)` of the optimality system
//!
//! ```text
//! A u = µ − Ω,      A φ = −W u,      µ_i = clamp(P_{t_i}(φ_i + µ_i), lo_i, hi_i)
//! ```
//!
//! and every step solves the linear system of one choice of pieces (see [`crate::linsys`]).

use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use ddg_core::sparse::Parallelism;

use crate::linsys::{LinearFailure, NewtonLinearSolver, Row};
use crate::problem::Discrete;
use crate::site::{Piece, Site};

/// The unknowns of the optimality system, on the unknown vertices.
#[derive(Debug, Clone)]
pub(crate) struct Iterate {
    /// Log scale factor `u`.
    pub u: Vec<f64>,
    /// Adjoint state `φ`.
    pub phi: Vec<f64>,
    /// Cone angles `µ`.
    pub mu: Vec<f64>,
}

impl Iterate {
    pub fn zeros(n: usize) -> Iterate {
        Iterate {
            u: vec![0.0; n],
            phi: vec![0.0; n],
            mu: vec![0.0; n],
        }
    }
}

/// Why a Newton iteration stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageStop {
    /// A full Newton step left every vertex on its piece (the active sets did not change), so
    /// the new iterate solves the stage's optimality system exactly (paper Sec. 5.3.3).
    ActiveSetsUnchanged,
    /// The residual of the stage's optimality system is at rounding level (regularized
    /// problems: a step changed pieces only within the rounding error of `φ`) or below the
    /// requested tolerance (exact problem).
    Tolerance,
    /// The line search of a regularized problem found no descent along the Newton direction
    /// although the residual is not small: the Newton system is too ill-conditioned. This
    /// happens when the minimizer is (nearly) not unique, e.g. with zero distortion weights.
    Stalled,
    /// The iteration limit was reached.
    IterationLimit,
    /// The exact iteration returned to a set of pieces it had already tried.
    Cycle,
    /// A Newton system could not be solved (singular or too ill-conditioned).
    LinearSolverFailure,
}

impl StageStop {
    /// Whether the stage solved its problem.
    pub fn is_success(self) -> bool {
        matches!(self, StageStop::ActiveSetsUnchanged | StageStop::Tolerance)
    }
}

/// What one Newton iteration did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NewtonOutcome {
    /// Newton steps, each one factorization of the Newton matrix.
    pub iterations: usize,
    /// Steps shortened by the line search (Stage I only).
    pub damped_steps: usize,
    /// Shifts of the adjoint constant that restored Gauss–Bonnet (closed meshes only).
    pub constant_shifts: usize,
    /// Why it stopped.
    pub stop: StageStop,
    /// Stage II on a closed mesh: the interval of constants that can be added to the returned
    /// `φ` without violating the optimality system. `(0, 0)` when `φ` is unique.
    pub adjoint_shift: (f64, f64),
}

/// Relative rounding error assumed for an adjoint value when deciding whether a vertex with a
/// fixed cone angle has really left its piece (64 machine epsilons).
const PIECE_NOISE: f64 = 64.0 * f64::EPSILON;
/// A regularized iteration whose line search finds no descent counts as solved if the residual
/// of the state equation is at most this, relative to the size of its terms.
const STALL_TOLERANCE: f64 = 1e-9;

fn max_abs(x: &[f64]) -> f64 {
    x.iter().fold(0.0, |m, v| m.max(v.abs()))
}

fn classify(sites: &[Site], k: &[f64], args: &[f64], out: &mut Vec<Piece>) {
    out.clear();
    out.extend(
        sites
            .iter()
            .zip(k)
            .zip(args)
            .map(|((site, &k), &x)| site.piece(k, x)),
    );
}

/// The constant `c` of smallest magnitude with
/// `Σ_i clamp(k_i · P_{t_i}(x_i + c), lo_i, hi_i) = target`.
///
/// On a closed mesh the cone angles must sum to `Σ Ω = 2πχ`. The left side is the sum of the
/// cone angles that the adjoint `φ + c` asks for; it is continuous, nondecreasing and piecewise
/// linear in `c`, so the equation is solved exactly on its breakpoints. For a regularized
/// problem this is the exact minimization of the dual objective along the constants, the null
/// direction of the Laplacian.
///
/// Returns `0` if the sum already matches (up to rounding) and the nearest attainable value if
/// the target is out of reach.
///
/// [formula: gauss-bonnet-adjoint-shift]
pub(crate) fn gauss_bonnet_shift(sites: &[Site], k: &[f64], x: &[f64], target: f64) -> f64 {
    let sum_at = |c: f64| -> f64 {
        sites
            .iter()
            .zip(k)
            .zip(x)
            .map(|((site, &k), &x)| site.measure(k, x + c))
            .sum()
    };
    let scale: f64 = sites
        .iter()
        .zip(k)
        .zip(x)
        .map(|((site, &k), &x)| site.measure(k, x).abs())
        .sum::<f64>()
        + target.abs();
    let h0 = sum_at(0.0);
    let tolerance = 64.0 * f64::EPSILON * scale;
    if (h0 - target).abs() <= tolerance {
        return 0.0;
    }
    // The sum must grow (direction +1) or shrink (direction −1).
    let direction = if h0 < target { 1.0 } else { -1.0 };
    let mut breakpoints: Vec<f64> = Vec::new();
    for ((site, &k), &x) in sites.iter().zip(k).zip(x) {
        for kink in site.kinks(k) {
            let c = kink - x;
            if c.is_finite() && direction * c > 0.0 {
                breakpoints.push(direction * c);
            }
        }
    }
    breakpoints.sort_unstable_by(f64::total_cmp);
    // `reached(c)`: moving by `c` in the chosen direction gets the sum to the target or past it.
    let reached = |c: f64| direction * (sum_at(direction * c) - target) >= 0.0;
    let (mut lo, mut hi) = (0, breakpoints.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if reached(breakpoints[mid]) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    if lo == breakpoints.len() {
        // Beyond the last breakpoint the sum is linear: constant if every vertex is clamped
        // there (the target is out of reach), sloped if some vertex has no bound on that side.
        let c_last = breakpoints.last().copied().unwrap_or(0.0);
        let h_last = sum_at(direction * c_last);
        let per_unit = sum_at(direction * (c_last + 1.0)) - h_last;
        return if direction * per_unit > 0.0 {
            direction * (c_last + (target - h_last) / per_unit)
        } else {
            direction * c_last
        };
    }
    let c1 = breakpoints[lo];
    let h1 = sum_at(direction * c1);
    let (c0, h0) = if lo == 0 {
        (0.0, h0)
    } else {
        (breakpoints[lo - 1], sum_at(direction * breakpoints[lo - 1]))
    };
    // Linear between two consecutive breakpoints.
    let c = if h1 == h0 {
        c1
    } else {
        c0 + (c1 - c0) * (target - h0) / (h1 - h0)
    };
    direction * c.clamp(c0, c1)
}

/// The derivative of the regularized dual objective along a Newton direction.
///
/// For a regularized problem with stiffness `k_i` and centre `µ̄_i` the dual objective is
///
/// ```text
/// Φ(φ) = ½ Σ w_i u_i² − Σ Ω_i φ_i + Σ e_i(φ_i),    W u = −A φ,
/// e_i'(φ) = clamp(k_i · P_{t_i}(φ + µ̄_i / k_i), lo_i, hi_i)
/// ```
///
/// which is convex and continuously differentiable; its minimizers are the adjoint states of
/// the regularized problem, and `∇Φ = −(A u − µ + Ω)`. Along a direction `(du, dφ)` with
/// `W du = −A dφ` its derivative at step length `α` is
///
/// ```text
/// g(α) = Σ w_i du_i (u_i + α du_i) − Σ Ω_i dφ_i + Σ e_i'(φ_i + α dφ_i) dφ_i
/// ```
///
/// which is nondecreasing and piecewise linear in `α`.
///
/// [formula: regularized-dual-objective]
struct DualSlope<'a> {
    sites: &'a [Site],
    k: &'a [f64],
    /// `φ_i + µ̄_i / k_i`.
    args: &'a [f64],
    dphi: &'a [f64],
    /// `Σ w_i du_i u_i − Σ Ω_i dφ_i`.
    constant: f64,
    /// `Σ w_i du_i²`.
    quadratic: f64,
}

impl DualSlope<'_> {
    fn at(&self, alpha: f64) -> f64 {
        let cones: f64 = self
            .sites
            .iter()
            .zip(self.k)
            .zip(self.args)
            .zip(self.dphi)
            .map(|(((site, &k), &x), &d)| site.measure(k, x + alpha * d) * d)
            .sum();
        self.constant + alpha * self.quadratic + cones
    }

    /// The step length in `[0, 1]` that minimizes the dual objective along the direction: `1`
    /// if the objective still decreases there, otherwise the root of `g`.
    fn minimizing_step(&self) -> f64 {
        let g1 = self.at(1.0);
        if g1 <= 0.0 {
            return 1.0;
        }
        let g0 = self.at(0.0);
        if g0 >= 0.0 {
            return 0.0;
        }
        let mut breakpoints: Vec<f64> = Vec::new();
        for (((site, &k), &x), &d) in self.sites.iter().zip(self.k).zip(self.args).zip(self.dphi) {
            if d == 0.0 {
                continue;
            }
            for kink in site.kinks(k) {
                let alpha = (kink - x) / d;
                if alpha > 0.0 && alpha < 1.0 {
                    breakpoints.push(alpha);
                }
            }
        }
        breakpoints.sort_unstable_by(f64::total_cmp);
        // First breakpoint at which the slope is no longer negative.
        let (mut lo, mut hi) = (0, breakpoints.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.at(breakpoints[mid]) >= 0.0 {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        let (a1, s1) = if lo == breakpoints.len() {
            (1.0, g1)
        } else {
            (breakpoints[lo], self.at(breakpoints[lo]))
        };
        let (a0, s0) = if lo == 0 {
            (0.0, g0)
        } else {
            (breakpoints[lo - 1], self.at(breakpoints[lo - 1]))
        };
        if s1 == s0 {
            return a0;
        }
        (a0 + (a1 - a0) * (-s0) / (s1 - s0)).clamp(a0, a1)
    }
}

/// Solves one regularized problem (paper Sec. 4.6, `ComputeRegularizedMeasure` of Fig. 12):
///
/// ```text
/// minimize  ½ Σ w_i u_i² + Σ g_i(µ_i) + Σ (µ_i − µ̄_i)² / (2 k_i)     subject to  A u = µ − Ω
/// ```
///
/// with stiffness `k_i > 0` and centre `µ̄` (zero in Stage I). Its optimality condition is
/// `µ_i = clamp(k_i · P_{t_i}(φ_i + µ̄_i/k_i), lo_i, hi_i)`; a fixed point `µ = µ̄` solves the
/// exact problem.
///
/// Each step solves the Newton system of the current pieces (paper Eqs. 22 and 25, with `µ`
/// eliminated). The step is the paper's whenever it leaves the pieces unchanged or the dual
/// objective is still decreasing at the full step; otherwise it is cut at the minimizer of the
/// dual objective along the step, which makes the iteration monotone.
///
/// On entry and exit the adjoint equation `W u + A φ = 0` holds for `x`.
///
/// [formula: newton-stopping-rule]
pub(crate) fn solve_regularized(
    d: &Discrete,
    solver: &mut dyn NewtonLinearSolver,
    k: &[f64],
    center: &[f64],
    x: &mut Iterate,
    max_iterations: usize,
    par: Parallelism<'_>,
) -> Result<NewtonOutcome, LinearFailure> {
    let n = d.n;
    let sites = &d.sites;
    let mut outcome = NewtonOutcome {
        iterations: 0,
        damped_steps: 0,
        constant_shifts: 0,
        stop: StageStop::IterationLimit,
        adjoint_shift: (0.0, 0.0),
    };
    let args_of = |phi: &[f64], out: &mut Vec<f64>| {
        out.clear();
        out.extend((0..n).map(|i| phi[i] + center[i] / k[i]));
    };
    let mut args = Vec::with_capacity(n);
    let mut pieces = Vec::with_capacity(n);
    let mut new_args = Vec::with_capacity(n);
    let mut rows = Vec::with_capacity(n);
    let mut au = vec![0.0; n];
    args_of(&x.phi, &mut args);
    classify(sites, k, &args, &mut pieces);
    while outcome.iterations < max_iterations {
        let mut singular = d.closed && !pieces.iter().any(|p| p.is_active());
        if singular {
            // No vertex pins the constant of φ: choose it to satisfy Gauss–Bonnet.
            let c = gauss_bonnet_shift(sites, k, &args, d.total_curvature);
            if c != 0.0 {
                for phi in &mut x.phi {
                    *phi += c;
                }
                args_of(&x.phi, &mut args);
                classify(sites, k, &args, &mut pieces);
                outcome.constant_shifts += 1;
                singular = !pieces.iter().any(|p| p.is_active());
            }
        }
        rows.clear();
        rows.extend((0..n).map(|i| match pieces[i] {
            Piece::Pos => Row::Sloped {
                k: k[i],
                offset: sites[i].t,
                base: center[i],
            },
            Piece::Neg => Row::Sloped {
                k: k[i],
                offset: -sites[i].t,
                base: center[i],
            },
            piece => Row::Flat {
                mu: sites[i].value_on(piece, k[i], args[i]),
            },
        }));
        // Without an active vertex on a closed mesh the system determines φ only up to a
        // constant: keep φ_0 and drop the (redundant) state equation of vertex 0.
        let kept = rows[0];
        if singular {
            rows[0] = Row::Pinned { phi: x.phi[0] };
        }
        let (u, phi) = solver.solve(&rows, par)?;
        outcome.iterations += 1;
        if singular {
            rows[0] = kept;
        }
        d.a.mul_vec(&u, &mut au).expect("dimensions match");
        let mu: Vec<f64> = (0..n)
            .map(|i| match rows[i] {
                Row::Flat { mu } => mu,
                _ => au[i] + d.omega[i],
            })
            .collect();
        args_of(&phi, &mut new_args);
        let new_pieces: Vec<Piece> = (0..n)
            .map(|i| sites[i].next_piece(k[i], pieces[i], new_args[i], mu[i]))
            .collect();
        if new_pieces == pieces {
            *x = Iterate { u, phi, mu };
            outcome.stop = StageStop::ActiveSetsUnchanged;
            return Ok(outcome);
        }
        // Pieces that changed only within rounding do not count. An active vertex is judged
        // by its cone angle, which the step determines accurately, so its change is real. A
        // vertex with a fixed cone angle is judged by `φ_i`, and a rounding error `ε` of `φ_i`
        // looks like a change `k_i ε` of the cone angle its piece asks for.
        let within_rounding = (0..n).all(|i| {
            new_pieces[i] == pieces[i]
                || (!pieces[i].is_active()
                    && (sites[i].measure(k[i], new_args[i]) - mu[i]).abs()
                        <= PIECE_NOISE * k[i] * (new_args[i].abs() + sites[i].t))
        });
        if within_rounding {
            *x = Iterate { u, phi, mu };
            outcome.stop = StageStop::Tolerance;
            return Ok(outcome);
        }
        let du: Vec<f64> = u.iter().zip(&x.u).map(|(a, b)| a - b).collect();
        let dphi: Vec<f64> = phi.iter().zip(&x.phi).map(|(a, b)| a - b).collect();
        let slope = DualSlope {
            sites,
            k,
            args: &args,
            dphi: &dphi,
            constant: (0..n)
                .map(|i| d.w[i] * du[i] * x.u[i] - d.omega[i] * dphi[i])
                .sum(),
            quadratic: (0..n).map(|i| d.w[i] * du[i] * du[i]).sum(),
        };
        let alpha = slope.minimizing_step();
        if alpha >= 1.0 {
            *x = Iterate { u, phi, mu };
            pieces = new_pieces;
            std::mem::swap(&mut args, &mut new_args);
        } else if alpha <= 0.0 {
            // The dual objective cannot be decreased along the Newton direction. At rounding
            // level that means the iterate is optimal; otherwise the direction is not accurate
            // enough to be a descent direction (an ill-conditioned system).
            d.a.mul_vec(&x.u, &mut au).expect("dimensions match");
            let residual = (0..n)
                .map(|i| (au[i] - x.mu[i] + d.omega[i]).abs())
                .fold(0.0, f64::max);
            let scale = 1.0 + max_abs(&d.omega) + max_abs(&x.mu);
            outcome.stop = if residual <= STALL_TOLERANCE * scale {
                StageStop::Tolerance
            } else {
                StageStop::Stalled
            };
            return Ok(outcome);
        } else {
            outcome.damped_steps += 1;
            for i in 0..n {
                x.u[i] += alpha * du[i];
                x.phi[i] += alpha * dphi[i];
            }
            args_of(&x.phi, &mut args);
            classify(sites, k, &args, &mut pieces);
            for i in 0..n {
                x.mu[i] = sites[i].value_on(pieces[i], k[i], args[i]);
            }
        }
    }
    Ok(outcome)
}

fn hash_pieces(pieces: &[Piece]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    pieces.hash(&mut hasher);
    hasher.finish()
}

/// The constants `c` for which `φ + c` satisfies the inclusion `φ_i + c ∈ ∂g_i(µ_i)` at every
/// vertex: the interval `[max_i (a_i − φ_i), min_i (b_i − φ_i)]` with `∂g_i(µ_i) = [a_i, b_i]`.
/// It is empty (first entry larger than the second) if no constant works.
///
/// On a closed mesh the adjoint equation determines `φ` only up to such a constant; with an
/// active vertex (`φ_i = ±t_i`) the interval is a point.
///
/// [formula: adjoint-constant-interval]
pub(crate) fn adjoint_constant_interval(sites: &[Site], phi: &[f64], mu: &[f64]) -> (f64, f64) {
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for ((site, &phi), &mu) in sites.iter().zip(phi).zip(mu) {
        let (a, b) = site.subdifferential(mu);
        lo = lo.max(a - phi);
        hi = hi.min(b - phi);
    }
    (lo, hi)
}

/// Solves the exact problem (`ComputeOptimalMeasure` of Fig. 12): the semismooth Newton method
/// on `F(u, φ, µ) = (A u − µ + Ω, A φ + W u, µ − clamp(P_t(φ + µ)))` (paper Eqs. 20 and 26).
///
/// The pieces are read off `φ + µ`; a step solves the linear system of these pieces exactly
/// (an active vertex gets `φ_i = ±t_i` and a free `µ_i`, every other vertex a fixed `µ_i`).
/// It stops when a step leaves the pieces unchanged, which means the iterate solves the
/// optimality system, or when the cone residual `max_i |µ_i − clamp(P_t(φ_i + µ_i))|` of a
/// solved step is at most `tolerance` (ties between pieces at rounding level).
///
/// Between two steps a vertex moves to a neighbouring piece only (lower bound – negative cone
/// – no cone – positive cone – upper bound). This differs from the paper's update only when a
/// vertex would skip a piece, e.g. turn from a positive into a negative cone in one step, which
/// is what makes the plain iteration diverge from a start with clustered cones. A step that
/// changes no piece is still exactly a step after which every vertex lies on its piece.
///
/// There is no convergence guarantee: the iteration is abandoned at the iteration limit or as
/// soon as a set of pieces repeats, and `x` is then left unchanged.
///
/// [formula: newton-stopping-rule]
pub(crate) fn solve_exact(
    d: &Discrete,
    solver: &mut dyn NewtonLinearSolver,
    x: &mut Iterate,
    max_iterations: usize,
    tolerance: f64,
    par: Parallelism<'_>,
) -> Result<NewtonOutcome, LinearFailure> {
    let n = d.n;
    let sites = &d.sites;
    let ones = vec![1.0; n];
    let mut outcome = NewtonOutcome {
        iterations: 0,
        damped_steps: 0,
        constant_shifts: 0,
        stop: StageStop::IterationLimit,
        adjoint_shift: (0.0, 0.0),
    };
    let mut y = x.clone();
    let mut seen: HashSet<u64> = HashSet::new();
    let mut z: Vec<f64> = y.phi.iter().zip(&y.mu).map(|(p, m)| p + m).collect();
    let mut pieces: Vec<Piece> = Vec::with_capacity(n);
    classify(sites, &ones, &z, &mut pieces);
    let mut rows = Vec::with_capacity(n);
    let mut au = vec![0.0; n];
    let stop = loop {
        if outcome.iterations >= max_iterations {
            break StageStop::IterationLimit;
        }
        let mut singular = d.closed && !pieces.iter().any(|p| p.is_active());
        if singular {
            // No vertex pins the constant of φ: choose it to satisfy Gauss–Bonnet.
            let c = gauss_bonnet_shift(sites, &ones, &z, d.total_curvature);
            if c != 0.0 {
                for i in 0..n {
                    y.phi[i] += c;
                    z[i] += c;
                }
                classify(sites, &ones, &z, &mut pieces);
                outcome.constant_shifts += 1;
                singular = !pieces.iter().any(|p| p.is_active());
            }
        }
        if !seen.insert(hash_pieces(&pieces)) {
            break StageStop::Cycle;
        }
        rows.clear();
        rows.extend((0..n).map(|i| match pieces[i] {
            Piece::Pos => Row::Pinned { phi: sites[i].t },
            Piece::Neg => Row::Pinned { phi: -sites[i].t },
            piece => Row::Flat {
                mu: sites[i].value_on(piece, 1.0, z[i]),
            },
        }));
        // Without an active vertex on a closed mesh the system determines φ only up to a
        // constant: keep φ_0 and drop the (redundant) state equation of vertex 0.
        let kept = rows[0];
        if singular {
            rows[0] = Row::Pinned { phi: y.phi[0] };
        }
        let (u, mut phi) = solver.solve(&rows, par)?;
        outcome.iterations += 1;
        if singular {
            rows[0] = kept;
        }
        d.a.mul_vec(&u, &mut au).expect("dimensions match");
        let mu: Vec<f64> = (0..n)
            .map(|i| match rows[i] {
                Row::Flat { mu } => mu,
                _ => au[i] + d.omega[i],
            })
            .collect();
        outcome.adjoint_shift = (0.0, 0.0);
        if singular {
            // φ is determined up to a constant. Take the admissible one whose mass-weighted
            // mean is closest to zero; if none is admissible, the one that violates the
            // inclusion least, so that the next pieces activate a vertex.
            let (lo, hi) = adjoint_constant_interval(sites, &phi, &mu);
            let mean: f64 = phi.iter().zip(&d.mass).map(|(p, m)| p * m).sum::<f64>()
                / d.mass.iter().sum::<f64>();
            let c = if lo <= hi {
                (-mean).clamp(lo, hi)
            } else {
                0.5 * (lo + hi)
            };
            for p in &mut phi {
                *p += c;
            }
            if lo <= hi {
                outcome.adjoint_shift = (lo - c, hi - c);
            }
        }
        for i in 0..n {
            z[i] = phi[i] + mu[i];
        }
        // The piece each vertex now lies on, approached one piece at a time.
        let new_pieces: Vec<Piece> = (0..n)
            .map(|i| {
                let target = sites[i].next_piece(1.0, pieces[i], z[i], mu[i]);
                sites[i].toward(pieces[i], target)
            })
            .collect();
        let residual = (0..n)
            .map(|i| (mu[i] - sites[i].prox(z[i])).abs())
            .fold(0.0, f64::max);
        y = Iterate { u, phi, mu };
        if new_pieces == pieces {
            break StageStop::ActiveSetsUnchanged;
        }
        if residual <= tolerance {
            break StageStop::Tolerance;
        }
        pieces = new_pieces;
    };
    outcome.stop = stop;
    if stop.is_success() {
        *x = y;
    }
    Ok(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;

    const INF: f64 = f64::INFINITY;

    #[test]
    fn shift_solves_the_sum_equation() {
        let sites = vec![
            Site::penalized(0.5, -INF, INF),
            Site::penalized(0.2, -1.0, 0.7),
            Site::penalized(0.0, -INF, INF),
            Site::fixed(0.3),
            Site::penalized(1.0, 0.0, 2.0),
        ];
        let k = vec![1.0, 3.0, 0.5, 1.0, 2.0];
        let x = vec![0.1, -0.4, 0.2, 5.0, 0.3];
        let sum = |c: f64| -> f64 {
            sites
                .iter()
                .zip(&k)
                .zip(&x)
                .map(|((s, &k), &x)| s.measure(k, x + c))
                .sum()
        };
        for target in [-3.0, -0.9, 0.0, 0.3, 0.4, 1.0, 2.5, 7.0] {
            let c = gauss_bonnet_shift(&sites, &k, &x, target);
            assert!(
                (sum(c) - target).abs() < 1e-12,
                "target {target}: c = {c}, sum = {}",
                sum(c)
            );
            // No constant of smaller magnitude does it.
            for f in [0.0, 0.5, 0.99] {
                if c != 0.0 {
                    assert!((sum(f * c) - target).abs() > 1e-12 || f * c == c);
                }
            }
        }
    }

    #[test]
    fn shift_reports_the_nearest_attainable_sum() {
        // Sum is confined to [0.3 − 1, 0.3 + 2].
        let sites = vec![Site::penalized(0.5, -1.0, 2.0), Site::fixed(0.3)];
        let k = vec![1.0, 1.0];
        let x = vec![0.0, 0.0];
        let c = gauss_bonnet_shift(&sites, &k, &x, 10.0);
        assert_eq!(sites[0].measure(1.0, c), 2.0);
        // All fixed: nothing to shift.
        let fixed = vec![Site::fixed(1.0), Site::fixed(2.0)];
        assert_eq!(gauss_bonnet_shift(&fixed, &k, &x, 5.0), 0.0);
    }

    #[test]
    fn constant_interval() {
        let sites = vec![
            Site::penalized(0.5, -INF, INF),
            Site::penalized(0.5, -INF, 1.0),
            Site::fixed(2.0),
        ];
        // Vertex 0 has no cone: φ + c ∈ [−0.5, 0.5]. Vertex 1 sits at its upper bound:
        // φ + c ≥ 0.5. Vertex 2 is fixed: no condition.
        let (lo, hi) = adjoint_constant_interval(&sites, &[0.1, 0.9, 7.0], &[0.0, 1.0, 2.0]);
        assert!((lo - (-0.4)).abs() < 1e-15 && (hi - 0.4).abs() < 1e-15);
        // An active vertex pins the constant.
        let (lo, hi) = adjoint_constant_interval(&sites, &[0.5, 0.9, 7.0], &[0.3, 1.0, 2.0]);
        assert_eq!((lo, hi), (0.0, 0.0));
    }
}
