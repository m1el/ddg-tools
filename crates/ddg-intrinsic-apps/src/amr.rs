//! Adaptive mesh refinement on intrinsic triangulations (paper Sec. 5.4).
//!
//! The paper's idea: use a local a posteriori error estimate "as an additional refinement
//! criterion in iDR" (Sec. 4.2), with the estimators of Morin, Nochetto, Siebert (2002) and
//! Mekchay, Nochetto (2005) as its supplement (Sec. 2) describes them. Here:
//!
//! * [`solve`]: piecewise linear finite elements on the triangulation (cotan stiffness,
//!   Galerkin mass) for the model problems of the paper ([`Problem`]): the harmonic Green's
//!   function, one backward Euler step of the heat equation from a point source (the
//!   short-time heat kernel), and a harmonic map of a disk onto the unit disk (Fig. 20).
//! * [`estimate`]: the residual indicator `η_T` and the oscillation `osc_T` of every face
//!   (Mekchay–Nochetto Eqs. 1.3.5, 1.3.8; supplement Sec. 2).
//! * [`mark`]: the bulk criterion of the supplement (Dörfler marking for `η` and `osc`).
//! * [`adaptive_solve`]: solve, estimate, mark, refine with `ddg-intrinsic`'s Delaunay
//!   refinement and the marked faces as an extra criterion
//!   ([`ddg_intrinsic::refine::refine`] with [`RefineRules::Ruppert`] (`IntrinsicTriangulation::refine_with`)), until a tolerance, a vertex budget or an
//!   iteration limit; every iteration is reported.
//!
//! Conventions: `Δ` is the positive semidefinite Laplace–Beltrami operator and `L` its cotan
//! matrix (the repository's `CLAUDE.md`); the PDE is `(c + Δ) u = f`. Decisions: SPEC N46 to
//! N49.

use std::fmt;
use std::time::Instant;

use crate::representation::InputPoints;
use ddg_core::{SparseError, SparseMatrix};
use ddg_intrinsic::{IntrinsicError, RefineOptions, RefineReport, RefineRules, Triangulation};

use crate::fem::galerkin_mass_matrix;
pub use crate::vector::Factorized;
use crate::vector::linsolve;

/// Errors of this module.
#[derive(Debug, Clone, PartialEq)]
pub enum AmrError {
    /// An option or a problem parameter is out of range.
    InvalidOptions {
        /// What is wrong.
        reason: String,
    },
    /// The triangulation does not fit the problem (a disk map on a surface that is not a disk,
    /// a source that is not a vertex, a disconnected surface).
    InvalidInput {
        /// What is wrong.
        reason: String,
    },
    /// Refinement failed.
    Intrinsic(IntrinsicError),
    /// A sparse factorization or solve failed.
    Sparse(SparseError),
}

impl fmt::Display for AmrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AmrError::InvalidOptions { reason } => write!(f, "invalid options: {reason}"),
            AmrError::InvalidInput { reason } => write!(f, "invalid input: {reason}"),
            AmrError::Intrinsic(e) => write!(f, "refinement: {e}"),
            AmrError::Sparse(e) => write!(f, "sparse solver: {e}"),
        }
    }
}

impl std::error::Error for AmrError {}

impl From<SparseError> for AmrError {
    fn from(e: SparseError) -> Self {
        AmrError::Sparse(e)
    }
}

impl From<IntrinsicError> for AmrError {
    fn from(e: IntrinsicError) -> Self {
        AmrError::Intrinsic(e)
    }
}

/// A model problem `(c + Δ) u = f`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Problem {
    /// The harmonic Green's function of an input vertex: `Δu = δ_s − 1/|M|` with `∫u = 0`
    /// (natural boundary conditions on a boundary). `c = 0`, `f = −1/|M|` inside the faces.
    GreenFunction {
        /// The source, an input vertex (its id is the same in every triangulation).
        source: u32,
    },
    /// The short-time heat kernel: one backward Euler step `(1 + tΔ) u = δ_s` of the heat
    /// equation from a unit point source (natural boundary conditions). `c = 1/t`, `f = 0`
    /// inside the faces.
    HeatKernel {
        /// The source, an input vertex.
        source: u32,
        /// The time step `t > 0`.
        time: f64,
    },
    /// A harmonic map of a disk onto the unit disk (paper Fig. 20): both coordinates harmonic,
    /// the boundary mapped to the unit circle by arc length (from the input boundary vertex
    /// with the smallest id, counter-clockwise). `c = 0`, `f = 0`, Dirichlet conditions.
    HarmonicDisk,
}

impl Problem {
    fn source(&self) -> Option<u32> {
        match *self {
            Problem::GreenFunction { source } | Problem::HeatKernel { source, .. } => Some(source),
            Problem::HarmonicDisk => None,
        }
    }

    fn dirichlet(&self) -> bool {
        matches!(self, Problem::HarmonicDisk)
    }
}

/// A finite element solution: one or two components (the disk map has two) per vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// `components[k][v]`: component `k` at vertex `v` of the triangulation.
    pub components: Vec<Vec<f64>>,
    /// The factorization used.
    pub factorization: Factorized,
    /// The largest relative residual of the linear systems.
    pub residual: f64,
}

fn check_connected<T: Triangulation + ?Sized>(tri: &T) -> Result<(), AmrError> {
    let n = tri.complex().n_vertices();
    let mut seen = vec![false; n];
    let mut stack = vec![0usize];
    seen[0] = true;
    let mut count = 1;
    while let Some(v) = stack.pop() {
        let mut visit = |w: usize| {
            if !seen[w] {
                seen[w] = true;
                count += 1;
                stack.push(w);
            }
        };
        for h in tri.complex().outgoing_halfedges(v) {
            visit(tri.complex().to_vertex(h));
            visit(tri.complex().to_vertex(tri.complex().next(h)));
        }
    }
    if count != n {
        return Err(AmrError::InvalidInput {
            reason: format!(
                "the surface is not connected: {count} of {n} vertices reachable from vertex 0"
            ),
        });
    }
    Ok(())
}

/// The boundary vertices of a disk in order, with their arc-length parameter in `[0, 1)`.
///
/// # Errors
/// [`AmrError::InvalidInput`] unless the triangulation has one boundary loop and Euler
/// characteristic 1.
pub fn disk_boundary<T: Triangulation + ?Sized>(tri: &T) -> Result<Vec<(u32, f64)>, AmrError> {
    let boundary: Vec<usize> = (0..tri.complex().n_vertices())
        .filter(|&v| tri.complex().is_boundary_vertex(v))
        .collect();
    let Some(&start) = boundary.first() else {
        return Err(AmrError::InvalidInput {
            reason: "the surface has no boundary; a disk map needs a disk".to_string(),
        });
    };
    if tri.complex().euler_characteristic() != 1 {
        return Err(AmrError::InvalidInput {
            reason: format!(
                "the surface has Euler characteristic {}, a disk has 1",
                tri.complex().euler_characteristic()
            ),
        });
    }
    let mut out = Vec::with_capacity(boundary.len());
    let mut v = start;
    let mut s = 0.0;
    loop {
        out.push((v as u32, s));
        let h = tri.complex().outgoing_halfedge(v);
        s += tri.edge_length(tri.complex().edge(h));
        v = tri.complex().to_vertex(h);
        if v == start || out.len() > boundary.len() {
            break;
        }
    }
    if out.len() != boundary.len() {
        return Err(AmrError::InvalidInput {
            reason: format!(
                "the boundary has more than one loop: the loop from vertex {start} has {} of {} boundary vertices",
                out.len(),
                boundary.len()
            ),
        });
    }
    for x in &mut out {
        x.1 /= s;
    }
    Ok(out)
}

/// Solves `problem` with piecewise linear finite elements on `tri`.
///
/// # Errors
/// [`AmrError::InvalidOptions`] for a time step that is not positive and finite;
/// [`AmrError::InvalidInput`] for a source that is not an input vertex, a disconnected surface,
/// or a disk map on a surface that is not a disk; [`AmrError::Sparse`].
pub fn solve<T: Triangulation + ?Sized>(tri: &T, problem: &Problem) -> Result<Solution, AmrError> {
    let n = tri.complex().n_vertices();
    if let Some(s) = problem.source()
        && s as usize >= tri.n_input_vertices()
    {
        return Err(AmrError::InvalidInput {
            reason: format!(
                "the source {s} is not an input vertex (the input mesh has {})",
                tri.n_input_vertices()
            ),
        });
    }
    check_connected(tri)?;
    let l = ddg_intrinsic::generic_ops::cotan_laplacian(tri);
    match *problem {
        Problem::GreenFunction { source } => {
            // L u = e_s − m/|M|, pinned at one vertex other than the source, then the mean
            // removed. The right-hand side sums to zero, so pinning changes nothing else.
            let mass = ddg_intrinsic::generic_ops::lumped_mass(tri);
            let total: f64 = mass.iter().sum();
            let mut b: Vec<f64> = mass.iter().map(|m| -m / total).collect();
            b[source as usize] += 1.0;
            let pin = if source == 0 { 1 } else { 0 };
            let keep: Vec<u32> = (0..n as u32).filter(|&v| v as usize != pin).collect();
            let l_red = l.submatrix(&keep, &keep)?;
            let mut x: Vec<f64> = keep.iter().map(|&v| b[v as usize]).collect();
            let (factor, factorization) = linsolve::factor(&l_red)?;
            let rhs = x.clone();
            linsolve::solve(&factor, &mut x)?;
            let residual = relative_residual(&l_red, &x, &rhs);
            let mut u = vec![0.0; n];
            for (k, &v) in keep.iter().enumerate() {
                u[v as usize] = x[k];
            }
            let mean = u.iter().zip(&mass).map(|(a, m)| a * m).sum::<f64>() / total;
            for a in &mut u {
                *a -= mean;
            }
            Ok(Solution {
                components: vec![u],
                factorization,
                residual,
            })
        }
        Problem::HeatKernel { source, time } => {
            if !(time.is_finite() && time > 0.0) {
                return Err(AmrError::InvalidOptions {
                    reason: format!("the time step must be positive and finite, got {time}"),
                });
            }
            let m = galerkin_mass_matrix(tri);
            let mut t: Vec<(u32, u32, f64)> = m.triplets().collect();
            t.extend(l.triplets().map(|(i, j, v)| (i, j, time * v)));
            let a = SparseMatrix::from_triplets(n, n, &t)?;
            let mut x = vec![0.0; n];
            x[source as usize] = 1.0;
            let rhs = x.clone();
            let (factor, factorization) = linsolve::factor(&a)?;
            linsolve::solve(&factor, &mut x)?;
            let residual = relative_residual(&a, &x, &rhs);
            Ok(Solution {
                components: vec![x],
                factorization,
                residual,
            })
        }
        Problem::HarmonicDisk => {
            let boundary = disk_boundary(tri)?;
            let mut is_b = vec![false; n];
            let mut values = [vec![0.0; n], vec![0.0; n]];
            for &(v, s) in &boundary {
                is_b[v as usize] = true;
                let (sin, cos) = (std::f64::consts::TAU * s).sin_cos();
                values[0][v as usize] = cos;
                values[1][v as usize] = sin;
            }
            let interior: Vec<u32> = (0..n as u32).filter(|&v| !is_b[v as usize]).collect();
            let bnd: Vec<u32> = boundary.iter().map(|x| x.0).collect();
            if interior.is_empty() {
                return Ok(Solution {
                    components: values.to_vec(),
                    factorization: Factorized::Cholesky,
                    residual: 0.0,
                });
            }
            let l_ii = l.submatrix(&interior, &interior)?;
            let l_ib = l.submatrix(&interior, &bnd)?;
            let (factor, factorization) = linsolve::factor(&l_ii)?;
            let mut residual = 0.0f64;
            for comp in &mut values {
                let xb: Vec<f64> = bnd.iter().map(|&v| comp[v as usize]).collect();
                let rhs: Vec<f64> = linsolve::mul(&l_ib, &xb).iter().map(|x| -x).collect();
                let mut x = rhs.clone();
                linsolve::solve(&factor, &mut x)?;
                residual = residual.max(relative_residual(&l_ii, &x, &rhs));
                for (k, &v) in interior.iter().enumerate() {
                    comp[v as usize] = x[k];
                }
            }
            Ok(Solution {
                components: values.to_vec(),
                factorization,
                residual,
            })
        }
    }
}

fn relative_residual(a: &SparseMatrix, x: &[f64], b: &[f64]) -> f64 {
    let ax = linsolve::mul(a, x);
    let r: Vec<f64> = ax.iter().zip(b).map(|(p, q)| p - q).collect();
    linsolve::norm(&r) / linsolve::norm(b).max(f64::MIN_POSITIVE)
}

/// The local error indicators of a solution.
#[derive(Debug, Clone, PartialEq)]
pub struct Estimate {
    /// Per face, `η_T` (not squared).
    pub indicator: Vec<f64>,
    /// Per face, `osc_T` (not squared).
    pub oscillation: Vec<f64>,
    /// `(Σ η_T²)^{1/2}`, the global estimator of Mekchay–Nochetto.
    pub total: f64,
    /// `(Σ osc_T²)^{1/2}`.
    pub total_oscillation: f64,
}

/// The flux `ℓ_e ∇u_T · ν` of a linear function through side `k` of face `f` (the side of
/// halfedge `3f + k`, outward normal `ν`): `cot θ_{k+1} (u_k − u_o) + cot θ_k (u_{k+1} − u_o)`
/// with `o` the opposite corner.
///
/// [formula: app-amr-edge-flux]
pub fn side_flux<T: Triangulation + ?Sized>(tri: &T, f: usize, k: usize, u: &[f64]) -> f64 {
    let v = tri.complex().faces()[f];
    let (i, j, o) = (
        v[k] as usize,
        v[(k + 1) % 3] as usize,
        v[(k + 2) % 3] as usize,
    );
    let l = tri.face_lengths(f);
    let cot = |c: usize| {
        ddg_core::geometry::corner_cotan(l[(c + 1) % 3], l[c], l[(c + 2) % 3])
            .expect("faces keep the triangle inequalities")
    };
    cot((k + 1) % 3) * (u[i] - u[o]) + cot(k) * (u[j] - u[o])
}

/// The residual indicator and the oscillation of every face for the solution `u` of `problem`
/// (supplement Sec. 2; Mekchay–Nochetto Eqs. 1.3.3 to 1.3.5 and 1.3.8):
///
/// ```text
/// η_T² = h_T² ‖R‖²_{L²(T)} + Σ_{e ⊂ ∂T} ℓ_e ‖J_e‖²_{L²(e)},   osc_T = h_T ‖R − R̄_T‖_{L²(T)}
/// ```
///
/// with `R = f − c u` (`Δu` vanishes inside a linear element), `J_e` the jump of the normal
/// derivative across `e` (on the boundary, for the natural conditions of the Green's function
/// and the heat kernel, the normal derivative itself; nothing for Dirichlet conditions), `R̄_T`
/// the mean of `R` over `T`, and `h_T` the diameter of `T` (its longest side). For several
/// components the squares add.
///
/// [formula: app-amr-residual-estimator]
/// [formula: app-amr-element-norms]
///
/// # Errors
/// [`AmrError::InvalidInput`] if the solution does not have one value per vertex.
pub fn estimate<T: Triangulation + ?Sized>(
    tri: &T,
    problem: &Problem,
    solution: &Solution,
) -> Result<Estimate, AmrError> {
    for c in &solution.components {
        if c.len() != tri.complex().n_vertices() {
            return Err(AmrError::InvalidInput {
                reason: format!(
                    "a solution with {} values on a triangulation with {} vertices",
                    c.len(),
                    tri.complex().n_vertices()
                ),
            });
        }
    }
    let nf = tri.complex().n_faces();
    let mut eta2 = vec![0.0; nf];
    let mut osc2 = vec![0.0; nf];
    let total_area = tri.total_area();
    for u in &solution.components {
        // Fluxes per halfedge.
        let flux: Vec<f64> = (0..tri.complex().n_halfedges())
            .map(|h| side_flux(tri, h / 3, h % 3, u))
            .collect();
        for f in 0..nf {
            let l = tri.face_lengths(f);
            let area = tri.face_area(f);
            let diameter = l[0].max(l[1]).max(l[2]);
            // The element residual at the corners (linear in the face).
            let v = tri.complex().faces()[f];
            let r: [f64; 3] = match *problem {
                Problem::GreenFunction { .. } => [-1.0 / total_area; 3],
                Problem::HeatKernel { time, .. } => v.map(|x| -u[x as usize] / time),
                Problem::HarmonicDisk => [0.0; 3],
            };
            let norm2 = area / 6.0
                * (r[0] * r[0]
                    + r[1] * r[1]
                    + r[2] * r[2]
                    + r[0] * r[1]
                    + r[1] * r[2]
                    + r[2] * r[0]);
            let mean = (r[0] + r[1] + r[2]) / 3.0;
            let osc = (norm2 - area * mean * mean).max(0.0);
            let mut jumps = 0.0;
            for k in 0..3 {
                let h = 3 * f + k;
                let lj = match tri.complex().twin(h) {
                    Some(t) => flux[h] + flux[t],
                    None if problem.dirichlet() => continue,
                    None => flux[h],
                };
                // ℓ_e ‖J‖² over the side = ℓ_e · ℓ_e J² = (ℓ_e J)².
                jumps += lj * lj;
            }
            eta2[f] += diameter * diameter * norm2 + jumps;
            osc2[f] += diameter * diameter * osc;
        }
    }
    let total = eta2.iter().sum::<f64>().sqrt();
    let total_oscillation = osc2.iter().sum::<f64>().sqrt();
    Ok(Estimate {
        indicator: eta2.into_iter().map(f64::sqrt).collect(),
        oscillation: osc2.into_iter().map(f64::sqrt).collect(),
        total,
        total_oscillation,
    })
}

/// Which sums the bulk criterion compares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BulkSum {
    /// `Σ_{marked} η_T > θ Σ η_T` (and the same for `osc`): the supplement's criterion as
    /// printed.
    Indicators,
    /// `Σ_{marked} η_T² ≥ θ² Σ η_T²`: Dörfler's marking as Mekchay–Nochetto (1.3.9) and Morin
    /// et al. state it.
    Squares,
}

/// The faces to refine: the fewest faces of largest indicator whose indicators reach the
/// fraction `bulk` of the total, then, if needed, faces of largest oscillation until the
/// oscillations reach it too (supplement Sec. 2; Mekchay–Nochetto, marking strategies E and
/// O). Ascending face ids.
///
/// [formula: app-amr-bulk-marking]
pub fn mark(estimate: &Estimate, bulk: f64, sum: BulkSum) -> Vec<u32> {
    let nf = estimate.indicator.len();
    let mut marked = vec![false; nf];
    let value = |x: f64| match sum {
        BulkSum::Indicators => x,
        BulkSum::Squares => x * x,
    };
    let fraction = match sum {
        BulkSum::Indicators => bulk,
        BulkSum::Squares => bulk * bulk,
    };
    for data in [&estimate.indicator, &estimate.oscillation] {
        let total: f64 = data.iter().map(|&x| value(x)).sum();
        if total <= 0.0 {
            continue;
        }
        let mut have: f64 = (0..nf).filter(|&f| marked[f]).map(|f| value(data[f])).sum();
        let mut order: Vec<usize> = (0..nf).collect();
        order.sort_by(|&a, &b| data[b].total_cmp(&data[a]).then(a.cmp(&b)));
        for f in order {
            let done = match sum {
                BulkSum::Indicators => have > fraction * total,
                BulkSum::Squares => have >= fraction * total,
            };
            if done {
                break;
            }
            if !marked[f] {
                marked[f] = true;
                have += value(data[f]);
            }
        }
    }
    (0..nf as u32).filter(|&f| marked[f as usize]).collect()
}

/// What identifies a face of a triangulation as unchanged: its id and its three edges.
fn face_key<T: Triangulation + ?Sized>(tri: &T, f: usize) -> (u32, [u32; 3]) {
    (
        f as u32,
        [0, 1, 2].map(|k| tri.complex().edge(3 * f + k) as u32),
    )
}

/// Refines `tri` (Delaunay refinement with `options`) with the faces `marked` as an extra
/// criterion: a marked face is refined while it is unchanged (same id, same three edges), so
/// each marked face gets a vertex or is removed by another insertion
/// ([`ddg_intrinsic::refine::refine`] with [`RefineRules::Ruppert`] (`IntrinsicTriangulation::refine_with`)).
///
/// # Errors
/// [`AmrError::InvalidInput`] for a face id out of range; the errors of refinement.
pub fn refine_marked<T: Triangulation>(
    tri: &mut T,
    options: &RefineOptions,
    marked: &[u32],
) -> Result<RefineReport, AmrError> {
    if let Some(&f) = marked
        .iter()
        .find(|&&f| f as usize >= tri.complex().n_faces())
    {
        return Err(AmrError::InvalidInput {
            reason: format!(
                "face {f} does not exist (the triangulation has {})",
                tri.complex().n_faces()
            ),
        });
    }
    let keys: std::collections::HashSet<(u32, [u32; 3])> =
        marked.iter().map(|&f| face_key(tri, f as usize)).collect();
    let criterion = |t: &T, f: usize| keys.contains(&face_key(t, f));
    Ok(ddg_intrinsic::refine::refine(
        tri,
        options,
        RefineRules::Ruppert,
        &criterion,
    )?)
}

/// Options of [`adaptive_solve`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AmrOptions {
    /// The refinement every step uses (angle bound, criterion, insertion limit per step). Its
    /// `max_circumradius` applies too.
    pub refine: RefineOptions,
    /// The bulk fraction `θ` (supplement: `1/2`).
    pub bulk: f64,
    /// Which sums the bulk criterion compares.
    pub bulk_sum: BulkSum,
    /// Stop when the global estimator is at most this (`0`: never).
    pub tolerance: f64,
    /// Stop before a step that would start with more vertices than this.
    pub max_vertices: usize,
    /// The limit on solve–estimate–mark–refine steps.
    pub max_iterations: usize,
}

impl Default for AmrOptions {
    /// The default refinement (25°), the supplement's criterion with `θ = 1/2`, no tolerance,
    /// at most 100 000 vertices and 60 steps.
    fn default() -> Self {
        AmrOptions {
            refine: RefineOptions::default(),
            bulk: 0.5,
            bulk_sum: BulkSum::Indicators,
            tolerance: 0.0,
            max_vertices: 100_000,
            max_iterations: 60,
        }
    }
}

/// How [`adaptive_solve`] ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AmrStatus {
    /// The estimator reached the tolerance.
    Tolerance,
    /// The next step would exceed the vertex budget.
    VertexLimit,
    /// The iteration limit was reached.
    IterationLimit,
    /// A refinement step inserted nothing (every marked face was refused).
    Stalled,
}

/// One solve–estimate–mark–refine step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AmrIteration {
    /// Vertices of the triangulation that was solved on.
    pub n_vertices: usize,
    /// Its faces.
    pub n_faces: usize,
    /// The global estimator.
    pub estimate: f64,
    /// The global oscillation.
    pub oscillation: f64,
    /// Faces marked.
    pub n_marked: usize,
    /// Vertices the refinement inserted (0 in the last step, which only solves).
    pub n_inserted: usize,
    /// Seconds spent in the solve.
    pub seconds_solve: f64,
    /// Seconds spent in the estimate and the marking.
    pub seconds_estimate: f64,
    /// Seconds spent in the refinement.
    pub seconds_refine: f64,
}

/// What [`adaptive_solve`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct AmrReport {
    /// How it ended.
    pub status: AmrStatus,
    /// The initial refinement to the angle bound (before the first solve).
    pub initial: RefineReport,
    /// Seconds of the initial refinement.
    pub seconds_initial: f64,
    /// Every step; the last one's solution is the result.
    pub iterations: Vec<AmrIteration>,
    /// Seconds in all, including the initial refinement.
    pub seconds_total: f64,
}

/// Adaptive refinement (paper Sec. 5.4): refine `tri` to the angle bound of
/// `options.refine`, then repeat: solve, estimate, stop if done, mark (the bulk criterion),
/// refine the marked faces together with every face below the angle bound
/// ([`ddg_intrinsic::refine::refine`] with [`RefineRules::Ruppert`] (`IntrinsicTriangulation::refine_with`), where a face counts as marked while it is
/// unchanged). Returns the last solution, on the final `tri`.
///
/// # Errors
/// [`AmrError::InvalidOptions`] for `bulk` outside `(0, 1]` or a negative tolerance; the errors
/// of [`solve`] and of refinement.
pub fn adaptive_solve<T: Triangulation>(
    tri: &mut T,
    problem: &Problem,
    options: &AmrOptions,
) -> Result<(Solution, AmrReport), AmrError> {
    if !(options.bulk > 0.0 && options.bulk <= 1.0) {
        return Err(AmrError::InvalidOptions {
            reason: format!("the bulk fraction must be in (0, 1], got {}", options.bulk),
        });
    }
    if options.tolerance.is_nan() || options.tolerance < 0.0 {
        return Err(AmrError::InvalidOptions {
            reason: format!(
                "the tolerance must be non-negative, got {}",
                options.tolerance
            ),
        });
    }
    let start = Instant::now();
    let initial =
        ddg_intrinsic::refine::refine(tri, &options.refine, RefineRules::Ruppert, &|_: &T, _| {
            false
        })?;
    let seconds_initial = start.elapsed().as_secs_f64();
    let mut iterations = Vec::new();
    loop {
        let t0 = Instant::now();
        let solution = solve(tri, problem)?;
        let t1 = Instant::now();
        let est = estimate(tri, problem, &solution)?;
        let mut it = AmrIteration {
            n_vertices: tri.complex().n_vertices(),
            n_faces: tri.complex().n_faces(),
            estimate: est.total,
            oscillation: est.total_oscillation,
            n_marked: 0,
            n_inserted: 0,
            seconds_solve: (t1 - t0).as_secs_f64(),
            seconds_estimate: 0.0,
            seconds_refine: 0.0,
        };
        let stop = if est.total <= options.tolerance {
            Some(AmrStatus::Tolerance)
        } else if iterations.len() + 1 >= options.max_iterations {
            Some(AmrStatus::IterationLimit)
        } else {
            None
        };
        if let Some(status) = stop {
            it.seconds_estimate = t1.elapsed().as_secs_f64();
            iterations.push(it);
            return Ok((
                solution,
                AmrReport {
                    status,
                    initial,
                    seconds_initial,
                    iterations,
                    seconds_total: start.elapsed().as_secs_f64(),
                },
            ));
        }
        let marked = mark(&est, options.bulk, options.bulk_sum);
        it.n_marked = marked.len();
        let t2 = Instant::now();
        it.seconds_estimate = (t2 - t1).as_secs_f64();
        let before = tri.complex().n_vertices();
        // Each marked face gets at least one vertex: stop if that alone exceeds the budget.
        if before + marked.len() > options.max_vertices {
            iterations.push(it);
            return Ok((
                solution,
                AmrReport {
                    status: AmrStatus::VertexLimit,
                    initial,
                    seconds_initial,
                    iterations,
                    seconds_total: start.elapsed().as_secs_f64(),
                },
            ));
        }
        refine_marked(tri, &options.refine, &marked)?;
        it.n_inserted = tri.complex().n_vertices() - before;
        it.seconds_refine = t2.elapsed().as_secs_f64();
        iterations.push(it);
        if it.n_inserted == 0 {
            let solution = solve(tri, problem)?;
            return Ok((
                solution,
                AmrReport {
                    status: AmrStatus::Stalled,
                    initial,
                    seconds_initial,
                    iterations,
                    seconds_total: start.elapsed().as_secs_f64(),
                },
            ));
        }
    }
}

/// The relative `L²` difference of two solutions at the input vertices (which every
/// triangulation of the same input shares, with the same ids), weighted by the lumped masses
/// of the input mesh, leaving out the vertices in `exclude` (a point source, where the Green's
/// function is singular): `(Σ m_v |u_v − r_v|²)^{1/2} / (Σ m_v |r_v|²)^{1/2}` over the
/// components.
///
/// # Errors
/// [`AmrError::InvalidInput`] if a solution has fewer values than the input has vertices, or
/// the two have different numbers of components.
pub fn input_vertex_error<T: Triangulation + ?Sized>(
    tri_input: &T,
    solution: &Solution,
    reference: &Solution,
    exclude: &[u32],
) -> Result<f64, AmrError> {
    let input = tri_input.input();
    let n = input.mesh().n_vertices();
    if solution.components.len() != reference.components.len()
        || solution
            .components
            .iter()
            .chain(&reference.components)
            .any(|c| c.len() < n)
    {
        return Err(AmrError::InvalidInput {
            reason: "the solutions do not cover the input vertices with the same components"
                .to_string(),
        });
    }
    let mass = ddg_core::operators::lumped_mass(input.mesh(), input.geometry()).map_err(|e| {
        AmrError::InvalidInput {
            reason: e.to_string(),
        }
    })?;
    let mut num = 0.0;
    let mut den = 0.0;
    for (u, r) in solution.components.iter().zip(&reference.components) {
        for v in 0..n {
            if exclude.contains(&(v as u32)) {
                continue;
            }
            num += mass[v] * (u[v] - r[v]).powi(2);
            den += mass[v] * r[v] * r[v];
        }
    }
    Ok((num / den.max(f64::MIN_POSITIVE)).sqrt())
}

/// The value of a piecewise linear function of `tri` at a point of `tri`.
fn interpolate<T: Triangulation + ?Sized>(
    tri: &T,
    u: &[f64],
    point: ddg_intrinsic::SurfacePoint,
) -> f64 {
    use ddg_intrinsic::SurfacePoint;
    match point {
        SurfacePoint::Vertex(v) => u[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = tri.complex().edge_vertices(edge as usize);
            (1.0 - t) * u[a] + t * u[b]
        }
        SurfacePoint::Face { face, bary } => {
            let v = tri.complex().faces()[face as usize];
            (0..3).map(|k| bary[k] * u[v[k] as usize]).sum()
        }
    }
}

/// The relative `L²` error of `solution` (on `tri`) against `reference` (on `reference_tri`,
/// a finer triangulation of the same input mesh): `‖u − r‖ / ‖r‖`, both integrals by the lumped
/// quadrature of the reference (its vertex areas), with `u` evaluated at the places of the
/// reference's vertices (point queries through the input mesh). For several components the
/// squares add.
///
/// # Errors
/// [`AmrError::InvalidInput`] if the components do not match; the errors of a point query.
pub fn l2_error<R: Triangulation + ?Sized, T: InputPoints + ?Sized>(
    reference_tri: &R,
    reference: &Solution,
    tri: &T,
    solution: &Solution,
) -> Result<f64, AmrError> {
    if reference.components.len() != solution.components.len() {
        return Err(AmrError::InvalidInput {
            reason: "the solutions have different numbers of components".to_string(),
        });
    }
    let mass = ddg_intrinsic::generic_ops::lumped_mass(reference_tri);
    let mut num = 0.0;
    let mut den = 0.0;
    for v in 0..reference_tri.complex().n_vertices() {
        let place = reference_tri.vertex_location(v);
        let here = tri.locate_on_triangulation(place)?;
        for (u, r) in solution.components.iter().zip(&reference.components) {
            let x = interpolate(tri, u, here);
            num += mass[v] * (x - r[v]).powi(2);
            den += mass[v] * r[v] * r[v];
        }
    }
    Ok((num / den.max(f64::MIN_POSITIVE)).sqrt())
}

/// Injectivity and distortion of a map of the triangulation to the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MapQuality {
    /// Faces whose image has non-positive signed area.
    pub n_flipped: usize,
    /// The smallest ratio of image area to face area.
    pub min_area_ratio: f64,
    /// The area-weighted mean of the conformal distortion `σ₁/σ₂` over unflipped faces.
    pub mean_conformal: f64,
    /// Its largest value.
    pub max_conformal: f64,
    /// The area-weighted quantile at 99 % of the conformal distortion.
    pub conformal_q99: f64,
}

/// The quality of a planar map given per vertex (`uv[0][v], uv[1][v]`): flipped faces and the
/// conformal distortion `σ₁/σ₂` of the affine map of every face (singular values of its
/// Jacobian).
///
/// [formula: app-amr-map-distortion]
pub fn map_quality<T: Triangulation + ?Sized>(tri: &T, uv: &Solution) -> MapQuality {
    let mut n_flipped = 0;
    let mut min_area_ratio = f64::INFINITY;
    let mut samples: Vec<(f64, f64)> = Vec::with_capacity(tri.complex().n_faces());
    for f in 0..tri.complex().n_faces() {
        let p = ddg_intrinsic::formulas::triangle_layout(tri.face_lengths(f))
            .expect("faces keep the triangle inequalities");
        let v = tri.complex().faces()[f];
        let q = v.map(|x| [uv.components[0][x as usize], uv.components[1][x as usize]]);
        let (e1, e2) = (
            [p[1][0] - p[0][0], p[1][1] - p[0][1]],
            [p[2][0] - p[0][0], p[2][1] - p[0][1]],
        );
        let (g1, g2) = (
            [q[1][0] - q[0][0], q[1][1] - q[0][1]],
            [q[2][0] - q[0][0], q[2][1] - q[0][1]],
        );
        let det_p = e1[0] * e2[1] - e1[1] * e2[0];
        let det_q = g1[0] * g2[1] - g1[1] * g2[0];
        let area = 0.5 * det_p;
        min_area_ratio = min_area_ratio.min(det_q / det_p);
        if det_q <= 0.0 {
            n_flipped += 1;
            continue;
        }
        // J = G P⁻¹ with P = [e1 e2], G = [g1 g2].
        let inv = [
            [e2[1] / det_p, -e2[0] / det_p],
            [-e1[1] / det_p, e1[0] / det_p],
        ];
        let j = [
            [
                g1[0] * inv[0][0] + g2[0] * inv[1][0],
                g1[0] * inv[0][1] + g2[0] * inv[1][1],
            ],
            [
                g1[1] * inv[0][0] + g2[1] * inv[1][0],
                g1[1] * inv[0][1] + g2[1] * inv[1][1],
            ],
        ];
        // J = conformal part + anticonformal part: σ₁ = |a| + |b|, σ₂ = |a| − |b| (no
        // cancellation for a nearly conformal map, unlike the closed form via the trace).
        let a = (0.5 * (j[0][0] + j[1][1])).hypot(0.5 * (j[1][0] - j[0][1]));
        let b = (0.5 * (j[0][0] - j[1][1])).hypot(0.5 * (j[1][0] + j[0][1]));
        samples.push(((a + b) / (a - b), area));
    }
    let total: f64 = samples.iter().map(|s| s.1).sum();
    let mean_conformal =
        samples.iter().map(|s| s.0 * s.1).sum::<f64>() / total.max(f64::MIN_POSITIVE);
    let max_conformal = samples.iter().fold(0.0f64, |m, s| m.max(s.0));
    samples.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut acc = 0.0;
    let mut conformal_q99 = max_conformal;
    for s in &samples {
        acc += s.1;
        if acc >= 0.99 * total {
            conformal_q99 = s.0;
            break;
        }
    }
    MapQuality {
        n_flipped,
        min_area_ratio,
        mean_conformal,
        max_conformal,
        conformal_q99,
    }
}
