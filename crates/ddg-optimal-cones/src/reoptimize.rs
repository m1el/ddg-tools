//! Re-optimization of the angles of given cones: the canonical problem with the given vertices
//! as free cone sites and every other interior vertex forbidden.
//!
//! With cone costs `0` at the given vertices and `+∞` elsewhere, the `ℓ¹` term of the canonical
//! problem is gone: the problem is the quadratic program
//!
//! ```text
//! minimize     ½ Σ_i wE_i m_i u_i²
//! subject to   (L u)_i = µ_i − Ω_i (i unknown),  µ_i = a_i (i forbidden or prescribed),
//!              lo_i ≤ µ_i ≤ hi_i (i a given cone)
//! ```
//!
//! (and the boundary sites of the treatment). If no bound is active at the optimum, its
//! optimality system is linear: `W u + A φ = 0`, `(A u)_i = a_i − Ω_i` at the fixed vertices,
//! `φ_i = 0` at the free ones (a free site's subdifferential is `{0}`); this is the saddle-point
//! system of one Newton step of the reference method with every piece known in advance
//! (`saddle-point-newton-system`), and a solution of it is a global minimizer
//! (`optimality-system-sufficient`). One sparse factorization and one solve give the answer.
//! The angles of the given cones are then `µ_i = (A u)_i + Ω_i`. If an angle comes out beyond
//! its bound, a bound is active and the problem is a quadratic program with box constraints;
//! then [`solve_with_boundary`](crate::solve_with_boundary) solves it (the factor-once method,
//! whose dense active-set part handles the bounds exactly, or the reference method).

use std::time::Instant;

use ddg_core::io::Cone;
use ddg_core::sparse::Parallelism;

use crate::certificate::{Certificate, certify_instance};
use crate::error::{ConeError, PrescribedConeIssue};
use crate::factor_once;
use crate::linsys::{LinearSolverReport, NewtonLinearSolver, Row, SaddleSolver};
use crate::problem::{BoundaryTreatment, ConeProblem, Discrete, Instance, split_pieces};
use crate::solve::{SolveOptions, SolveReport, Status, solve_part};

/// How [`reoptimize_cones`] solves the problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReoptimizeMethod {
    /// The faster of the two for the problem at hand: the generic solver (in practice its
    /// factor-once method, i.e. the reduced system of the free sites) when a component has at
    /// most [`REDUCED_SITES_LIMIT`] free sites and the factor-once method supports it; the direct
    /// solve otherwise (and the generic solver where that does not apply). The default.
    #[default]
    Auto,
    /// The direct solve of the linear optimality system wherever it applies (no bound active,
    /// no cone cost left, not singular); the generic solver otherwise.
    DirectKkt,
    /// Always the generic solver ([`solve_with_boundary`](crate::solve_with_boundary) on the
    /// problem of [`reoptimization_problem`]).
    Solve,
}

/// [`ReoptimizeMethod::Auto`] uses the generic solver on a component with at most this many free
/// sites (given cones, plus the boundary vertices with an optimized boundary). There the
/// factor-once method solves the reduced system: one Cholesky factorization of the `n × n`
/// Laplacian block, one solve per free site, a dense system of that size; the direct solve
/// factorizes the `2n × 2n` saddle-point matrix instead, which costs about as much as 30 to 45
/// extra solves. Measured on the dataset (one thread, `examples/reoptimize.rs`): the generic solver
/// is 1.3 to 2.6 times faster up to 29 cones, the two are even at 36 (`rhino`) and 45 (`brain`),
/// and the direct solve is 6.6 and 34 times faster at 439 and 13 836 cones.
pub const REDUCED_SITES_LIMIT: usize = 32;

/// Which method produced the result of a component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReoptimizeMethodUsed {
    /// One sparse factorization and solve of the linear optimality system.
    DirectKkt,
    /// The generic solver (with the options in [`ReoptimizeOptions::solve`]).
    Solve,
}

/// Options of [`reoptimize_cones`].
#[derive(Debug, Clone, Default)]
pub struct ReoptimizeOptions<'a> {
    /// What happens at the boundary (as in [`solve_with_boundary`](crate::solve_with_boundary)).
    /// Default [`BoundaryTreatment::Fixed`]: `u = 0` there.
    pub boundary: BoundaryTreatment,
    /// The method. Default [`ReoptimizeMethod::Auto`].
    pub method: ReoptimizeMethod,
    /// Options of the generic solver where it is used; also the certificate tolerance
    /// (`certificate_tolerance`), the linear solver of the direct solve (`linear_solver`) and the
    /// thread pool (`parallelism`; with several components they are solved concurrently on it).
    pub solve: SolveOptions<'a>,
}

/// What happened on one connected component.
#[derive(Debug, Clone, PartialEq)]
pub struct ReoptimizeComponentReport {
    /// The number of the component (`TriMesh::face_components`; `0` for a connected mesh).
    pub component: u32,
    /// Its vertices.
    pub n_vertices: usize,
    /// Its unknowns (interior vertices, or all with an optimized boundary or on a closed
    /// component).
    pub n_unknowns: usize,
    /// The given cones on it.
    pub n_cones: usize,
    /// Whether every vertex is an unknown (closed, or an optimized boundary).
    pub closed: bool,
    /// The method that produced its result.
    pub method: ReoptimizeMethodUsed,
    /// Why the direct solve was not used (or not accepted), if it was not.
    pub direct_not_used: Option<String>,
    /// The status of its part of the problem (its certificate decides).
    pub status: Status,
    /// The certificate of its part of the problem (its share of the objective).
    pub certificate: Certificate,
    /// Its part of the distortion `½ Σ wE m u²` (masses divided by the area of the whole mesh).
    pub distortion: f64,
    /// The report of the generic solver, if it was used.
    pub solve: Option<SolveReport>,
}

/// What [`reoptimize_cones`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct ReoptimizeReport {
    /// [`Status::Converged`] iff the certificate of the whole problem holds with
    /// `options.solve.certificate_tolerance`.
    pub status: Status,
    /// [`ReoptimizeMethodUsed::DirectKkt`] if every component with unknowns was solved directly.
    pub method: ReoptimizeMethodUsed,
    /// The certificate of the whole problem of [`reoptimization_problem`] for the returned
    /// `(log_scale, adjoint, measure)`: `certify_with_boundary(&reoptimization_problem(..)?,
    /// boundary, ..)` reproduces it.
    pub certificate: Certificate,
    /// The factorizations and solves of the direct solves (summed over the components).
    pub linear_solver: LinearSolverReport,
    /// One entry per connected component (one for a connected mesh).
    pub components: Vec<ReoptimizeComponentReport>,
    /// Wall time in seconds.
    pub seconds: f64,
    /// What did not go to plan, if anything.
    pub message: Option<String>,
}

/// The result of [`reoptimize_cones`].
#[derive(Debug, Clone, PartialEq)]
pub struct Reoptimization {
    /// The given cones with their optimal angles, in the order given (an angle may be `0`).
    pub cones: Vec<Cone>,
    /// `µ` per vertex: the angle of a given cone, the prescribed angle, `0` elsewhere; with an
    /// optimized boundary `ν` at the boundary vertices.
    pub measure: Vec<f64>,
    /// The log scale factor `u` per vertex.
    pub log_scale: Vec<f64>,
    /// The adjoint `φ` per vertex (zero at the given cones, whose angle is free).
    pub adjoint: Vec<f64>,
    /// `ν` at the boundary vertices (as [`Solution::boundary_curvature`](crate::Solution)).
    pub boundary_curvature: Vec<f64>,
    /// The distortion `½ Σ wE m u²` (masses divided by the total area).
    pub distortion: f64,
    /// The objective of the problem: the distortion, plus the boundary cost with
    /// [`BoundaryTreatment::Penalized`].
    pub objective: f64,
    /// `Σ |angle|` of the given cones.
    pub total_cone_angle: f64,
    /// What was done, and the certificate.
    pub report: ReoptimizeReport,
}

/// The canonical problem that [`reoptimize_cones`] solves: `problem` with cone cost `0` at the
/// given vertices and `+∞` at every other interior vertex. Everything else is kept: the
/// distortion weights, the bounds (they apply to the given cones), the prescribed cones, and at
/// boundary vertices the cone costs of `problem` (they matter only for
/// [`BoundaryTreatment::Penalized`]). The cone costs of `problem` at interior vertices are not
/// used.
///
/// # Errors
/// [`ConeError::InvalidCone`] (with angle `0`) for a vertex that is out of range, listed twice,
/// on the boundary (no cone can be placed there; the treatment governs the boundary), or carries
/// a prescribed cone; [`ConeError::WrongLength`] / [`ConeError::InvalidConeCost`] for the cone
/// costs of `problem`.
pub fn reoptimization_problem<'a>(
    problem: &ConeProblem<'a>,
    cone_vertices: &[u32],
) -> Result<ConeProblem<'a>, ConeError> {
    let mesh = problem.mesh;
    let nv = mesh.n_vertices();
    let mut is_cone = vec![false; nv];
    let mut prescribed = vec![false; nv];
    for cone in &problem.prescribed_cones {
        if let Some(p) = prescribed.get_mut(cone.vertex as usize) {
            *p = true;
        }
    }
    for (index, &v) in cone_vertices.iter().enumerate() {
        let issue = |reason| ConeError::InvalidCone {
            index,
            vertex: v,
            angle: 0.0,
            reason,
        };
        let i = v as usize;
        if i >= nv {
            return Err(issue(PrescribedConeIssue::VertexOutOfRange));
        }
        if mesh.is_boundary_vertex(i) {
            return Err(issue(PrescribedConeIssue::BoundaryVertex));
        }
        if is_cone[i] {
            return Err(issue(PrescribedConeIssue::DuplicateVertex));
        }
        if prescribed[i] {
            return Err(issue(PrescribedConeIssue::PrescribedVertex));
        }
        is_cone[i] = true;
    }
    if let Some(costs) = &problem.cone_costs
        && costs.len() != nv
    {
        return Err(ConeError::WrongLength {
            what: "cone_costs",
            expected: nv,
            got: costs.len(),
        });
    }
    let cone_costs: Vec<f64> = (0..nv)
        .map(|v| {
            if mesh.is_boundary_vertex(v) {
                problem.cone_costs.as_ref().map_or(1.0, |c| c[v])
            } else if is_cone[v] {
                0.0
            } else {
                f64::INFINITY
            }
        })
        .collect();
    Ok(ConeProblem {
        cone_costs: Some(cone_costs),
        ..problem.clone()
    })
}

/// The result of one component, per vertex of the component.
struct PartOutcome {
    measure: Vec<f64>,
    log_scale: Vec<f64>,
    adjoint: Vec<f64>,
    linear: LinearSolverReport,
    report: ReoptimizeComponentReport,
}

/// The direct solve of the linear optimality system of a problem without cone costs, or why it
/// does not apply. `Ok(Err(reason))`: not applicable or not accepted.
///
/// [formula: saddle-point-newton-system]
/// [formula: optimality-system-sufficient]
#[allow(clippy::type_complexity)]
fn direct(
    d: &Discrete,
    options: &ReoptimizeOptions<'_>,
) -> Result<Result<(Vec<f64>, Vec<f64>, Vec<f64>, LinearSolverReport), String>, ConeError> {
    if d.n == 0 {
        return Ok(Ok((
            Vec::new(),
            Vec::new(),
            Vec::new(),
            LinearSolverReport::default(),
        )));
    }
    if let Some(i) = d.sites.iter().position(|s| !s.fixed && s.t > 0.0) {
        return Ok(Err(format!(
            "vertex {} has a cone cost (a penalized boundary): the problem keeps an l1 term",
            d.labels[i]
        )));
    }
    if d.closed && d.sites.iter().all(|s| s.fixed) {
        return Ok(Err(
            "every angle is fixed on a closed component: the linear system is singular".to_string(),
        ));
    }
    let rows: Vec<Row> = d
        .sites
        .iter()
        .map(|s| {
            if s.fixed {
                Row::Flat { mu: s.lo }
            } else {
                Row::Pinned { phi: 0.0 }
            }
        })
        .collect();
    let mut solver = SaddleSolver::new(d, options.solve.linear_solver)?;
    let (u, phi) = match solver.solve(&rows, Parallelism::Sequential) {
        Ok(x) => x,
        Err(failure) => {
            return Ok(Err(format!(
                "the linear system could not be solved: {}",
                failure.message
            )));
        }
    };
    let mut au = vec![0.0; d.n];
    d.a.mul_vec(&u, &mut au)?;
    let mut mu = vec![0.0; d.n];
    for i in 0..d.n {
        let s = &d.sites[i];
        mu[i] = if s.fixed { s.lo } else { au[i] + d.omega[i] };
        // Written so that a NaN counts as outside.
        let inside = s.lo <= mu[i] && mu[i] <= s.hi;
        if !s.fixed && !inside {
            return Ok(Err(format!(
                "the unconstrained angle {} at vertex {} is outside its bounds [{}, {}]: a bound is \
                 active",
                mu[i], d.labels[i], s.lo, s.hi
            )));
        }
    }
    Ok(Ok((u, phi, mu, solver.report())))
}

/// Re-optimizes one connected problem (the whole problem, or a component).
fn reoptimize_part(
    instance: &Instance<'_>,
    lambda: f64,
    n_cones: usize,
    options: &ReoptimizeOptions<'_>,
) -> Result<PartOutcome, ConeError> {
    let d = Discrete::from_instance(instance)?;
    let tolerance = options.solve.certificate_tolerance;
    let why;
    let n_free = d.sites.iter().filter(|s| !s.fixed).count();
    let reduced_is_faster = options.method == ReoptimizeMethod::Auto
        && n_free <= REDUCED_SITES_LIMIT
        && factor_once::check_support(&d, lambda).is_ok();
    if reduced_is_faster {
        why = Some(format!(
            "{n_free} free sites (at most {REDUCED_SITES_LIMIT}): the reduced system of the \
             factor-once method is faster"
        ));
    } else if options.method != ReoptimizeMethod::Solve {
        match direct(&d, options)? {
            Ok((u, phi, mu, linear)) => {
                let (u, phi, mu) = (d.to_vertices(&u), d.to_vertices(&phi), d.to_vertices(&mu));
                let certificate = certify_instance(instance, &u, &phi, &mu);
                if certificate.holds(tolerance) {
                    let distortion = 0.5
                        * instance
                            .data
                            .w_e
                            .iter()
                            .zip(&instance.mass)
                            .zip(&u)
                            .map(|((w, m), u)| w * m * u * u)
                            .sum::<f64>();
                    return Ok(PartOutcome {
                        report: ReoptimizeComponentReport {
                            component: instance.component.unwrap_or(0),
                            n_vertices: instance.mesh.n_vertices(),
                            n_unknowns: d.n,
                            n_cones,
                            closed: d.closed,
                            method: ReoptimizeMethodUsed::DirectKkt,
                            direct_not_used: None,
                            status: Status::Converged,
                            certificate,
                            distortion,
                            solve: None,
                        },
                        measure: mu,
                        log_scale: u,
                        adjoint: phi,
                        linear,
                    });
                }
                why = Some(format!(
                    "the certificate of the direct solve does not hold (state residual {:.1e}, \
                     gap {:.1e})",
                    certificate.state_residual, certificate.duality_gap
                ));
            }
            Err(reason) => why = Some(reason),
        }
    } else {
        why = Some("ReoptimizeMethod::Solve was requested".to_string());
    }
    let solve_options = SolveOptions {
        parallelism: Parallelism::Sequential,
        ..options.solve.clone()
    };
    let part = solve_part(instance, &d, lambda, &solve_options)?;
    let report = part.report.report;
    Ok(PartOutcome {
        report: ReoptimizeComponentReport {
            component: instance.component.unwrap_or(0),
            n_vertices: instance.mesh.n_vertices(),
            n_unknowns: d.n,
            n_cones,
            closed: d.closed,
            method: ReoptimizeMethodUsed::Solve,
            direct_not_used: why,
            status: report.status,
            certificate: report.certificate,
            distortion: part.report.distortion_optimal,
            solve: Some(report),
        },
        measure: part.measure,
        log_scale: part.log_scale_optimal,
        adjoint: part.adjoint,
        linear: LinearSolverReport::default(),
    })
}

/// Finds the optimal cone **angles** for given cone **positions**: the canonical problem with
/// the given vertices as free cone sites (cost `0`), cones forbidden at every other interior
/// vertex (cost `+∞`), and the angle bounds, distortion weights and prescribed cones of
/// `problem` respected ([`reoptimization_problem`]). `problem.lambda` matters only for
/// [`BoundaryTreatment::Penalized`] (the cost of the boundary curvature); otherwise nothing is
/// penalized and the result does not depend on `λ`.
///
/// **Method.** Without an `ℓ¹` term the problem is a quadratic program. Its optimality system is
/// linear when no bound is active: one sparse `LDLᵀ` factorization of the `2n × 2n`
/// saddle-point matrix `[[W, A], [A, 0]]` (the rows of the given cones replaced by `φ_i = 0`) and
/// one solve with iterative refinement give `u` and `φ`, and `µ_i = (A u)_i + Ω_i` at the given
/// cones ([`ReoptimizeMethodUsed::DirectKkt`]). The cost does not depend on the number of cones.
/// With few free sites the reduced system is cheaper: the generic solver's factor-once method
/// factorizes the `n × n` Laplacian block once, solves once per free site, and solves the dense
/// system of the free sites exactly (bounds included, by its active-set method); with
/// [`ReoptimizeMethod::Auto`] it is used up to [`REDUCED_SITES_LIMIT`] free sites per component
/// (the crossover measured on the dataset). If the direct solution violates a bound (a bound
/// is active), or a boundary site has a cost (`Penalized`), or the system is singular (a closed
/// component whose angles are all fixed), the generic solver
/// ([`solve_with_boundary`](crate::solve_with_boundary) on [`reoptimization_problem`], with
/// `options.solve`) solves the component instead ([`ReoptimizeMethodUsed::Solve`]); the report
/// says why. Every result carries the certificate of the problem.
///
/// The mesh may have several components (each is solved on its own as part of the whole
/// problem, as in [`solve`](crate::solve)). On a closed component the angles must sum to its
/// `2πχ`: with at least one given cone (unbounded) that is automatic; with none, or with bounds
/// that do not allow it, the error is [`ConeError::Infeasible`] naming the component.
///
/// # Errors
/// As [`reoptimization_problem`] for the cone vertices (an empty list is allowed: no cones),
/// the errors of the problem's validation and of the options, [`ConeError::Infeasible`].
///
/// [formula: saddle-point-newton-system]
/// [formula: component-separation]
pub fn reoptimize_cones(
    problem: &ConeProblem<'_>,
    cone_vertices: &[u32],
    options: &ReoptimizeOptions<'_>,
) -> Result<Reoptimization, ConeError> {
    let start = Instant::now();
    options.solve.validate()?;
    let reopt = reoptimization_problem(problem, cone_vertices)?;
    let mesh = problem.mesh;
    let nv = mesh.n_vertices();
    let whole = Instance::new(&reopt, options.boundary)?;
    let pieces = if mesh.n_components() > 1 {
        split_pieces(mesh, problem.geometry)?
    } else {
        Vec::new()
    };
    let component_of = mesh.vertex_components();
    let mut cones_per = vec![0usize; mesh.n_components()];
    for &v in cone_vertices {
        cones_per[component_of[v as usize] as usize] += 1;
    }

    let outcomes: Vec<PartOutcome> = if pieces.is_empty() {
        vec![reoptimize_part(
            &whole,
            reopt.lambda,
            cone_vertices.len(),
            options,
        )?]
    } else {
        let parts: Vec<Instance<'_>> = pieces
            .iter()
            .enumerate()
            .map(|(c, piece)| whole.part(piece, c as u32))
            .collect();
        // Validate every component (feasibility) before solving any.
        for part in &parts {
            Discrete::from_instance(part)?;
        }
        let one = |c: usize| reoptimize_part(&parts[c], reopt.lambda, cones_per[c], options);
        let results: Vec<Result<PartOutcome, ConeError>> = match options.solve.parallelism {
            Parallelism::Pool(pool) => pool.map(parts.len(), one),
            Parallelism::Sequential => (0..parts.len()).map(one).collect(),
        };
        results.into_iter().collect::<Result<_, _>>()?
    };

    let mut measure = vec![0.0; nv];
    let mut log_scale = vec![0.0; nv];
    let mut adjoint = vec![0.0; nv];
    let mut linear_solver = LinearSolverReport::default();
    let mut distortion = 0.0;
    let mut components = Vec::with_capacity(outcomes.len());
    let mut messages = Vec::new();
    for (c, outcome) in outcomes.into_iter().enumerate() {
        match pieces.get(c) {
            Some(piece) => {
                for (l, &v) in piece.sub.vertex_to_parent.iter().enumerate() {
                    let v = v as usize;
                    measure[v] = outcome.measure[l];
                    log_scale[v] = outcome.log_scale[l];
                    adjoint[v] = outcome.adjoint[l];
                }
            }
            None => {
                measure = outcome.measure;
                log_scale = outcome.log_scale;
                adjoint = outcome.adjoint;
            }
        }
        linear_solver.accumulate(&outcome.linear);
        distortion += outcome.report.distortion;
        if let Some(message) = outcome
            .report
            .solve
            .as_ref()
            .and_then(|r| r.message.as_ref())
        {
            messages.push(format!("component {c}: {message}"));
        }
        components.push(outcome.report);
    }
    let certificate = certify_instance(&whole, &log_scale, &adjoint, &measure);
    let status = if certificate.holds(options.solve.certificate_tolerance) {
        Status::Converged
    } else if components
        .iter()
        .any(|c| c.status == Status::LinearSolverFailed)
    {
        Status::LinearSolverFailed
    } else {
        Status::NotConverged
    };
    let method = if components
        .iter()
        .all(|c| c.n_unknowns == 0 || c.method == ReoptimizeMethodUsed::DirectKkt)
    {
        ReoptimizeMethodUsed::DirectKkt
    } else {
        ReoptimizeMethodUsed::Solve
    };
    // The cost of the boundary curvature (Penalized only; every other site costs nothing).
    let cone_cost: f64 = whole
        .data
        .sites
        .iter()
        .zip(&whole.data.unknown)
        .zip(&measure)
        .filter(|((_, unknown), _)| **unknown)
        .map(|((site, _), &mu)| site.cost(mu))
        .sum();
    let cones: Vec<Cone> = cone_vertices
        .iter()
        .map(|&v| Cone {
            vertex: v,
            angle: measure[v as usize],
        })
        .collect();
    let total_cone_angle = cones.iter().map(|c| c.angle.abs()).sum();
    let boundary_curvature =
        crate::solve::boundary_curvature_of(&reopt, options.boundary, &log_scale, &measure)?;
    Ok(Reoptimization {
        cones,
        measure,
        log_scale,
        adjoint,
        boundary_curvature,
        distortion,
        objective: distortion + cone_cost,
        total_cone_angle,
        report: ReoptimizeReport {
            status,
            method,
            certificate,
            linear_solver,
            components,
            seconds: start.elapsed().as_secs_f64(),
            message: (!messages.is_empty()).then(|| messages.join("; ")),
        },
    })
}
