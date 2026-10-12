//! The factor-once method: one Cholesky factorization of the Laplacian, a working set of cone
//! sites, and a dense active-set method on it.
//!
//! Eliminating `u` from the canonical problem gives a problem in the cone angles alone (see
//! `reduced.rs`): with `Z` the columns `z_j = A⁻¹ e_j` of the free sites and `u = Z µ − w0`,
//!
//! ```text
//! minimize   ½ µᵀ G µ + qᵀ µ + Σ_j t_j |µ_j|     subject to  lo_j ≤ µ_j ≤ hi_j   [and Σ µ = c]
//! G = Zᵀ W Z,    q = −Zᵀ W w0
//! ```
//!
//! whose solution is that of the canonical problem. The method never forms `G`: it grows a
//! *working set* of sites, computes their columns with the fixed factor (independent solves,
//! on the thread pool), keeps the Gram matrix of the working set, and solves the problem on it
//! with a dense active-set method. After each dense solve, one solve for `u` and one for the
//! gradient `y = A⁻¹ W u` give the optimality conditions of *all* sites; violators join the
//! working set. When none remains, `(µ, u, y)` satisfies the optimality system of the full
//! problem, up to the tolerance. There is no convergence proof for this loop: it ends in a
//! result that the certificate checks, or it gives up (see [`FallbackReason`]).

mod dense;
mod reduced;

use std::time::Instant;

use ddg_core::sparse::{Parallelism, ThreadPool};

use crate::certificate::Certificate;
use crate::newton::Iterate;
use crate::problem::Discrete;
pub use dense::DenseStats;
use dense::{Dense, violation_at_zero};
use reduced::Reduced;

/// Which method [`solve`](crate::solve) uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SolveMethod {
    /// The factor-once method whenever the problem is within what it supports, and the
    /// reference method otherwise, or when the factor-once method gives up (see
    /// [`FallbackReason`]) or its result fails the certificate. The report says which method
    /// produced the result, and why. The default.
    #[default]
    Auto,
    /// The paper's regularized stages and semismooth Newton method (a saddle-point
    /// factorization per step), always.
    Reference,
    /// The factor-once method, always: no fallback to the reference method. A problem outside
    /// what it supports is an error ([`ConeError::InvalidOptions`](crate::ConeError)); if it gives up or its
    /// result fails the certificate, the best point it has is returned with
    /// [`Status::NotConverged`](crate::Status) and the reason in the report.
    FactorOnce,
}

/// Which method produced a result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethodUsed {
    /// The reference method.
    Reference,
    /// The factor-once method.
    FactorOnce,
}

/// Why a problem is outside what the factor-once method supports.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Unsupported {
    /// `λ = 0`: every site is free, the method has nothing to be sparse about.
    ZeroLambda,
    /// A distortion weight `wE_i m_i` is zero at an unknown (a vertex of the unknown set, with
    /// the vertex number): the Gram matrix `Zᵀ W Z` is singular and the minimizer may not be
    /// unique.
    ZeroDistortionWeight {
        /// The vertex.
        vertex: u32,
    },
    /// A mesh whose Laplacian block is singular (closed, or an optimized boundary) where the
    /// bounds of the free sites do not add up to the required angle sum (only possible within the
    /// feasibility tolerance of the problem), so there is no feasible starting point.
    NoFeasiblePin,
    /// Nothing is unknown.
    NoUnknowns,
}

impl std::fmt::Display for Unsupported {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unsupported::ZeroLambda => write!(f, "lambda is 0"),
            Unsupported::ZeroDistortionWeight { vertex } => write!(
                f,
                "the distortion weight of vertex {vertex} is 0 (the Gram matrix would be singular)"
            ),
            Unsupported::NoFeasiblePin => write!(
                f,
                "the bounds of the free cone sites do not add up to the total curvature"
            ),
            Unsupported::NoUnknowns => write!(f, "the problem has no unknowns"),
        }
    }
}

/// Why the reference method produced the result although the factor-once method was allowed
/// (or why the factor-once method stopped without a certified result).
#[derive(Debug, Clone, PartialEq)]
pub enum FallbackReason {
    /// The problem is outside what the factor-once method supports.
    Unsupported(Unsupported),
    /// The working set would have grown past [`FactorOnceOptions::max_working_set`]: beyond it the
    /// dense part and the stored columns cost more than the reference method.
    WorkingSetLimit {
        /// The size the working set would have reached.
        working_set: usize,
        /// The limit.
        limit: usize,
    },
    /// The stored columns and dense matrices would have exceeded
    /// [`FactorOnceOptions::memory_limit_bytes`].
    MemoryLimit {
        /// Bytes the next step would have needed (estimate).
        required_bytes: usize,
        /// The limit.
        limit_bytes: usize,
    },
    /// The outer loop reached [`FactorOnceOptions::max_rounds`], or the dense active-set
    /// method repeatedly reached its iteration cap.
    IterationLimit {
        /// Rounds done.
        rounds: usize,
    },
    /// The factor-once method ended, but the certificate of its result does not hold with
    /// [`SolveOptions::certificate_tolerance`](crate::SolveOptions::certificate_tolerance).
    CertificateFailed {
        /// The certificate of the discarded result.
        certificate: Box<Certificate>,
    },
    /// The factorization or a solve failed.
    FactorizationFailed {
        /// What failed.
        message: String,
    },
}

impl From<Unsupported> for FallbackReason {
    fn from(u: Unsupported) -> FallbackReason {
        FallbackReason::Unsupported(u)
    }
}

impl std::fmt::Display for FallbackReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FallbackReason::Unsupported(u) => write!(f, "problem not supported: {u}"),
            FallbackReason::WorkingSetLimit { working_set, limit } => write!(
                f,
                "working set of {working_set} sites would exceed the limit of {limit}"
            ),
            FallbackReason::MemoryLimit {
                required_bytes,
                limit_bytes,
            } => write!(
                f,
                "stored columns would need {required_bytes} bytes, the limit is {limit_bytes}"
            ),
            FallbackReason::IterationLimit { rounds } => {
                write!(f, "iteration limit reached after {rounds} rounds")
            }
            FallbackReason::CertificateFailed { certificate } => write!(
                f,
                "certificate does not hold (state {:.1e}, adjoint {:.1e}, cone {:.1e}, bound {:.1e}, dual {:.1e}, gap {:.1e})",
                certificate.state_residual,
                certificate.adjoint_residual,
                certificate.cone_residual,
                certificate.bound_violation,
                certificate.dual_infeasibility,
                certificate.duality_gap
            ),
            FallbackReason::FactorizationFailed { message } => {
                write!(f, "linear algebra failed: {message}")
            }
        }
    }
}

/// Options of the factor-once method ([`SolveOptions::factor_once`](crate::SolveOptions)).
#[derive(Debug, Clone, PartialEq)]
pub struct FactorOnceOptions {
    /// A site outside the working set enters it if its optimality condition is violated by more
    /// than this times `λ`. Looser values (`1e-8`) leave the certificate at `1e-9` failing on
    /// meshes with many near-zero cones (measured on `camelshotnoise`, see the README); the
    /// default costs nothing. Default `1e-12`.
    pub tolerance: f64,
    /// Give up when the working set (the sites whose columns are stored) would exceed this
    /// size. `None` (the default) uses the crossover with the reference method, measured on
    /// meshes of 11 000 to 153 000 vertices on one thread (README, "The factor-once method"):
    /// about `1300 · (n / 15000)^0.2` sites for `n` free sites, and up to twice that when the
    /// pool has 8 or more threads (`1 + (threads − 1)/7`, at most 2): the dense active-set part
    /// and the Gram product dominate beyond it, and they do not parallelize as well as the
    /// column solves. A result that needs more is left to the reference method.
    pub max_working_set: Option<usize>,
    /// Give up when the stored columns (`8 n` bytes each) and dense matrices would exceed this
    /// many bytes. Default 4 GiB.
    pub memory_limit_bytes: usize,
    /// Give up after this many rounds (a round is one dense solve, one gradient and one batch of
    /// new columns). Default 200.
    pub max_rounds: usize,
}

impl Default for FactorOnceOptions {
    fn default() -> Self {
        FactorOnceOptions {
            tolerance: 1e-12,
            max_working_set: None,
            memory_limit_bytes: 4 << 30,
            max_rounds: 200,
        }
    }
}

/// Wall times of the factor-once method, in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FactorOnceTimings {
    /// Factorization and `w0`.
    pub setup: f64,
    /// Column solves.
    pub columns: f64,
    /// Gram products.
    pub gram: f64,
    /// Dense active-set solves.
    pub dense: f64,
    /// Gradient (`u` and `y`).
    pub gradient: f64,
    /// Final accurate `u` and adjoint.
    pub finish: f64,
    /// Everything.
    pub total: f64,
}

/// What the factor-once method did.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorOnceReport {
    /// Rounds (dense solve, gradient, new columns).
    pub rounds: usize,
    /// Solves with the Laplacian factor for columns.
    pub column_solves: usize,
    /// All other solves (gradients, `w0`, the final accurate `u` and adjoint, refinement).
    pub other_solves: usize,
    /// Sites whose columns were computed.
    pub working_set: usize,
    /// Nonzero cone angles at the end (including clamped ones).
    pub final_nonzero: usize,
    /// Largest number of nonzero cone angles during the run.
    pub peak_nonzero: usize,
    /// Free sites: the number of variables of the reduced problem.
    pub sites: usize,
    /// Largest violation (relative to `λ`) of an optimality condition of a site in the working
    /// set, as seen by the exact gradient, when the loop ended.
    pub max_working_set_excess: f64,
    /// Estimated peak memory of the stored columns and dense matrices, in bytes.
    pub peak_bytes: usize,
    /// Counters of the dense iteration.
    pub dense: DenseStats,
    /// Times.
    pub seconds: FactorOnceTimings,
}

/// The outcome of one attempt of the method.
pub(crate) struct Attempt {
    /// The iterate (when it was asked for, or the attempt succeeded).
    pub x: Option<Iterate>,
    pub report: FactorOnceReport,
    pub gave_up: Option<FallbackReason>,
}

fn pool_of<'a>(par: Parallelism<'a>) -> Option<&'a ThreadPool> {
    match par {
        Parallelism::Pool(p) => Some(p),
        Parallelism::Sequential => None,
    }
}

/// Checks whether the factor-once method supports a problem (a discrete form of it).
pub(crate) fn check_support(d: &Discrete, lambda: f64) -> Result<(), Unsupported> {
    Reduced::check(d, lambda).map(|_| ())
}

/// Runs the method. `want_iterate`: also produce the iterate if the method gave up.
///
/// [formula: reduced-lasso-problem]
/// [formula: reduced-optimality-conditions]
/// [formula: reduced-sum-constraint-pin]
pub(crate) fn attempt(
    d: &Discrete,
    lambda: f64,
    options: &FactorOnceOptions,
    par: Parallelism<'_>,
    want_iterate: bool,
) -> Attempt {
    let t_start = Instant::now();
    let mut report = FactorOnceReport {
        rounds: 0,
        column_solves: 0,
        other_solves: 0,
        working_set: 0,
        final_nonzero: 0,
        peak_nonzero: 0,
        sites: 0,
        max_working_set_excess: 0.0,
        peak_bytes: 0,
        dense: DenseStats::default(),
        seconds: FactorOnceTimings::default(),
    };
    let red = match Reduced::new(d, lambda) {
        Ok(r) => r,
        Err(reason) => {
            report.seconds.total = t_start.elapsed().as_secs_f64();
            return Attempt {
                x: None,
                report,
                gave_up: Some(reason),
            };
        }
    };
    report.seconds.setup = t_start.elapsed().as_secs_f64();
    let mut run = Run::new(&red, lambda, options, par, report);
    let gave_up = run.rounds().err();
    let x = if gave_up.is_none() || want_iterate {
        match run.finish() {
            Ok(x) => Some(x),
            Err(reason) => {
                return run.into_attempt(None, Some(gave_up.unwrap_or(reason)), t_start);
            }
        }
    } else {
        None
    };
    run.into_attempt(x, gave_up, t_start)
}

/// The working-set size at which the factor-once method stops being faster than the reference
/// method: `1300 · (n / 15000)^0.2`, times `1 + (threads − 1)/7` (at most 2).
///
/// Measured crossovers (README): one thread, `cow1` (15 049 vertices) 1 200 sites, `rhino` and
/// `bunny` 1 300 to 1 400, `brain` (152 893) about 1 900; with 8 threads `cow1` 1 900 and `brain`
/// beyond 3 000. The constant is the round figure between them.
pub(crate) fn automatic_working_set_limit(n: usize, threads: usize) -> usize {
    let thread_factor = 1.0 + (threads.clamp(1, 8) - 1) as f64 / 7.0;
    (1300.0 * (n as f64 / 15000.0).powf(0.2) * thread_factor).round() as usize
}

struct Run<'a> {
    red: &'a Reduced<'a>,
    lambda: f64,
    options: &'a FactorOnceOptions,
    pool: Option<&'a ThreadPool>,
    threads: usize,
    chunk: usize,
    report: FactorOnceReport,
    dense: Dense,
    /// Dense index -> slot.
    w_slots: Vec<usize>,
    /// Slot -> dense index.
    in_w: Vec<usize>,
    /// Stored columns: `(data, columns)`; dense index of the first is `offset`.
    blocks: Vec<(Vec<f64>, usize)>,
    offset: usize,
    limit: usize,
    nu: f64,
    y: Vec<f64>,
}

impl<'a> Run<'a> {
    fn new(
        red: &'a Reduced<'a>,
        lambda: f64,
        options: &'a FactorOnceOptions,
        par: Parallelism<'a>,
        mut report: FactorOnceReport,
    ) -> Run<'a> {
        let pool = pool_of(par);
        let threads = pool.map_or(1, |p| p.n_threads());
        let chunk = if threads <= 1 {
            16
        } else if threads <= 4 {
            8
        } else {
            4
        };
        let ns = red.n_slots();
        report.sites = ns;
        let limit = options
            .max_working_set
            .unwrap_or_else(|| automatic_working_set_limit(ns, threads));
        let mut dense = Dense::new(red.sum, red.plan.as_ref().map_or(0.0, |p| p.pin_value));
        let mut w_slots = Vec::new();
        let mut in_w = vec![usize::MAX; ns];
        if red.pin.is_some() {
            let s = &red.d.sites[red.slots[0] as usize];
            dense.extend(&[0.0], 1, &[0.0], &[s.t], &[s.lo], &[s.hi]);
            w_slots.push(0);
            in_w[0] = 0;
        }
        let offset = w_slots.len();
        Run {
            red,
            lambda,
            options,
            pool,
            threads,
            chunk,
            report,
            dense,
            w_slots,
            in_w,
            blocks: Vec::new(),
            offset,
            limit,
            nu: 0.0,
            y: Vec::new(),
        }
    }

    fn site(&self, slot: usize) -> &crate::site::Site {
        &self.red.d.sites[self.red.slots[slot] as usize]
    }

    fn fail(e: ddg_core::sparse::SparseError) -> FallbackReason {
        FallbackReason::FactorizationFailed {
            message: e.to_string(),
        }
    }

    fn rounds(&mut self) -> Result<(), FallbackReason> {
        let red = self.red;
        let ns = red.n_slots();
        let tol = self.options.tolerance * self.lambda;
        let mut cap_hits = 0;
        let mut refinements = 0;
        // No single site can carry the sum: start with the sites of largest curvature clamped.
        if let Some(plan) = red.plan.as_ref().filter(|p| !p.clamped.is_empty()) {
            let slots: Vec<usize> = plan
                .clamped
                .iter()
                .map(|&(i, _)| red.slot_of[i] as usize)
                .collect();
            self.add_columns(&slots)?;
            for &(i, bound) in &plan.clamped {
                let j = self.in_w[red.slot_of[i] as usize];
                self.dense.clamp(j, bound);
            }
        }
        loop {
            if self.report.rounds >= self.options.max_rounds {
                return Err(FallbackReason::IterationLimit {
                    rounds: self.report.rounds,
                });
            }
            self.report.rounds += 1;
            let t = Instant::now();
            let converged = self.dense.solve(tol * 0.1, 20 * (self.dense.dim() + 100));
            self.report.seconds.dense += t.elapsed().as_secs_f64();
            if converged {
                cap_hits = 0;
            } else {
                cap_hits += 1;
                if cap_hits >= 3 {
                    return Err(FallbackReason::IterationLimit {
                        rounds: self.report.rounds,
                    });
                }
            }
            self.nu = self.dense.nu;
            let nz = self.dense.nonzeros();
            self.report.peak_nonzero = self.report.peak_nonzero.max(nz.len());
            // The exact gradient.
            let t = Instant::now();
            let use_gemv = (nz.len() as f64)
                <= red.factor_per_row * if self.pool.is_some() { 4.0 } else { 1.0 };
            let u = if use_gemv {
                let mut starts = Vec::with_capacity(self.blocks.len());
                let mut acc = self.offset;
                for b in &self.blocks {
                    starts.push(acc);
                    acc += b.1;
                }
                let terms: Vec<(usize, usize, f64)> = nz
                    .iter()
                    .filter(|&&j| j >= self.offset)
                    .map(|&j| {
                        let b = starts.partition_point(|&s| s <= j) - 1;
                        (b, j - starts[b], self.dense.mu[j])
                    })
                    .collect();
                red.u_from_columns(&self.blocks, &terms, self.pool)
            } else {
                let mut mu = red.mu_fixed.clone();
                for &j in &nz {
                    mu[red.slots[self.w_slots[j]] as usize] = self.dense.mu[j];
                }
                // No refinement inside the loop; the answer is refined at the end.
                let rhs: Vec<f64> = mu.iter().zip(&red.d.omega).map(|(m, o)| m - o).collect();
                let mut x = red.gather(&rhs);
                red.solve_rows(&mut x, 1).map_err(Self::fail)?;
                let mut u = red.scatter(&x);
                red.project(&mut u);
                u
            };
            // `u_from_columns` ignores the fixed sites' constant; it is inside w0.
            let y = red.gradient(&u, false).map_err(Self::fail)?;
            self.report.seconds.gradient += t.elapsed().as_secs_f64();
            // Violators among the sites.
            let mut cand: Vec<(usize, f64)> = Vec::new();
            let mut max_in_w = 0.0f64;
            for s in 0..ns {
                let v = red.grad_of_slot(&y, s) + self.nu;
                let site = self.site(s);
                let (t_s, lo, hi) = (site.t, site.lo, site.hi);
                let found = if self.in_w[s] == usize::MAX {
                    violation_at_zero(t_s, lo, hi, v)
                } else {
                    self.dense.exact_violation(self.in_w[s], v)
                };
                let Some((ex, _)) = found else { continue };
                if ex > tol {
                    if self.in_w[s] == usize::MAX {
                        cand.push((s, ex));
                    } else {
                        max_in_w = max_in_w.max(ex);
                    }
                }
            }
            self.report.max_working_set_excess = max_in_w / self.lambda;
            // Re-anchor the dense gradient to the exact one (iterative refinement).
            let exact: Vec<f64> = (0..self.dense.dim())
                .map(|j| red.grad_of_slot(&y, self.w_slots[j]))
                .collect();
            self.dense.anchor_gradient(&exact);
            if cand.is_empty() {
                if max_in_w <= tol {
                    return Ok(());
                }
                refinements += 1;
                if refinements > 8 {
                    // The conditions of the working set do not improve; the certificate decides.
                    return Ok(());
                }
                continue;
            }
            refinements = 0;
            // A batch of new columns, spread out over the mesh.
            cand.sort_by(|a, b| b.1.total_cmp(&a.1));
            let k = ((0.5 * self.dense.n_nonzero() as f64) as usize).max(16);
            let mut chosen: Vec<usize> = Vec::new();
            let mut blocked = vec![false; red.n];
            for &(s, _) in &cand {
                if chosen.len() >= k {
                    break;
                }
                let i = red.slots[s] as usize;
                if blocked[i] {
                    continue;
                }
                chosen.push(s);
                let mut frontier = vec![i];
                blocked[i] = true;
                for _ in 0..3 {
                    let mut next = Vec::new();
                    for &f in &frontier {
                        let cp = red.d.a.col_ptr();
                        for &w in &red.d.a.row_idx()[cp[f]..cp[f + 1]] {
                            if !blocked[w as usize] {
                                blocked[w as usize] = true;
                                next.push(w as usize);
                            }
                        }
                    }
                    frontier = next;
                }
            }
            self.add_columns(&chosen)?;
        }
    }

    /// Computes the columns of the given slots and adds them to the working set. Gives up if
    /// the working set or the memory would exceed the limits.
    fn add_columns(&mut self, chosen: &[usize]) -> Result<(), FallbackReason> {
        let red = self.red;
        let new_size = self.w_slots.len() + chosen.len();
        if new_size > self.limit {
            return Err(FallbackReason::WorkingSetLimit {
                working_set: new_size,
                limit: self.limit,
            });
        }
        let bytes = 8 * red.n * (new_size - self.offset) + Dense::bytes_for(new_size);
        self.report.peak_bytes = self.report.peak_bytes.max(bytes);
        if bytes > self.options.memory_limit_bytes {
            return Err(FallbackReason::MemoryLimit {
                required_bytes: bytes,
                limit_bytes: self.options.memory_limit_bytes,
            });
        }
        let t = Instant::now();
        let zn = red
            .columns(chosen, self.pool, self.chunk)
            .map_err(Self::fail)?;
        self.report.seconds.columns += t.elapsed().as_secs_f64();
        self.report.column_solves += chosen.len();
        let t = Instant::now();
        let (cross, q) = red.gram_cross(
            &self.blocks,
            self.offset,
            &zn,
            chosen.len(),
            self.pool,
            self.threads,
        );
        self.report.seconds.gram += t.elapsed().as_secs_f64();
        let ts: Vec<f64> = chosen.iter().map(|&s| self.site(s).t).collect();
        let los: Vec<f64> = chosen.iter().map(|&s| self.site(s).lo).collect();
        let his: Vec<f64> = chosen.iter().map(|&s| self.site(s).hi).collect();
        self.dense.extend(&cross, chosen.len(), &q, &ts, &los, &his);
        for &s in chosen {
            self.in_w[s] = self.w_slots.len();
            self.w_slots.push(s);
        }
        self.blocks.push((zn, chosen.len()));
        Ok(())
    }

    /// The iterate (accurate `u`, adjoint, cone angles over the unknowns) of the current point.
    fn finish(&mut self) -> Result<Iterate, FallbackReason> {
        let t = Instant::now();
        let red = self.red;
        let mut mu = red.mu_fixed.clone();
        for j in self.dense.nonzeros() {
            mu[red.slots[self.w_slots[j]] as usize] = self.dense.mu[j];
        }
        let u = red.u_of_mu(&mu).map_err(Self::fail)?;
        let yr = red.gradient(&u, true).map_err(Self::fail)?;
        let mut phi = vec![0.0; red.n];
        for (r, &i) in red.rows.iter().enumerate() {
            phi[i as usize] = -(yr[r] + self.nu);
        }
        if let Some(p) = red.pin {
            phi[p] = -self.nu;
        }
        // At a cone strictly inside its bounds the optimality condition is φ = ±t exactly (as
        // in the reference method); the adjoint computed from the gradient is that up to
        // rounding (64 ulp at most: a larger change would show in `L φ`, whose entries can be
        // large on needle-like triangles). Anything further off is left for the certificate.
        for (i, site) in red.d.sites.iter().enumerate() {
            let m = mu[i];
            if !site.fixed && site.t > 0.0 && m != 0.0 && site.lo < m && m < site.hi {
                let target = site.t.copysign(m);
                if (phi[i] - target).abs() <= 64.0 * f64::EPSILON * site.t {
                    phi[i] = target;
                }
            }
        }
        self.report.seconds.finish = t.elapsed().as_secs_f64();
        Ok(Iterate { u, phi, mu })
    }

    fn into_attempt(
        mut self,
        x: Option<Iterate>,
        gave_up: Option<FallbackReason>,
        t_start: Instant,
    ) -> Attempt {
        let red = self.red;
        self.report.working_set = self.w_slots.len() - self.offset;
        self.report.final_nonzero = self.dense.n_nonzero();
        self.report.dense = self.dense.stats;
        let solves = red.solves.load(std::sync::atomic::Ordering::Relaxed);
        self.report.other_solves = solves.saturating_sub(self.report.column_solves);
        self.report.seconds.total = t_start.elapsed().as_secs_f64();
        let _ = &self.y;
        Attempt {
            x,
            report: self.report,
            gave_up,
        }
    }
}
