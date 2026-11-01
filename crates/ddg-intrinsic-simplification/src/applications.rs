//! Geometry algorithms on the coarse mesh, prolongated to the fine mesh (paper Sec. 8.5;
//! SPEC §10.1, S20, S21).
//!
//! Everything here takes a [`CoarseMesh`] (any coarse mesh: the output of `simplify`, or one made
//! by [`greedy_coarsening`]), builds [`CoarseApplications`] from it once, and works on the
//! compacted coarse mesh ([`CoarseOutput`]) and the scalar prolongation `P` (`|V| × |V̄|`,
//! SPEC §7.2). The fine mesh is never touched: a quantity of the fine mesh is `P` applied to the
//! same quantity of the coarse mesh.
//!
//! # Operators (item 1)
//!
//! [`CoarseOutput::cotan_laplacian`] and [`CoarseOutput::lumped_mass`] are the cotan Laplacian
//! (symmetric positive semidefinite, the workspace convention) and the lumped mass of the Δ-complex
//! from its lengths, per halfedge as `ddg_intrinsic::generic_ops` does: a halfedge `i → j` with
//! `i ≠ j` adds `−½ cot` of its opposite angle to `L_ij` and `L_ji`, a self-edge adds nothing, a
//! face that uses a vertex twice counts twice for it in the mass.
//! [formula: cotan-laplacian] [formula: lumped-mass]
//!
//! # Poisson (item 2, Fig. 21)
//!
//! [`PoissonSolver`] solves `L x = b` with Dirichlet values on a set of vertices, or, with none,
//! on a closed surface with Neumann conditions: the right-hand side is made compatible (its
//! `M`-weighted mean removed, `rhs_shift` in the report) and the solution is the one with
//! `M`-weighted mean zero. It works for any `L`, `M`, so the same code solves the fine problem.
//! [`CoarseApplications::poisson`] takes data `f̄` on the coarse mesh (`L̄ x̄ = M̄ f̄`) and returns
//! `x = P x̄`. **Fine data** enter by one of two restrictions, the choice being the caller's:
//! [`CoarseApplications::sample`] (`f̄_i = f` at the fine vertex of coarse vertex `i`; the pointwise
//! restriction of a function, then `b̄ = M̄ f̄`), or [`CoarseApplications::restrict_load`]
//! (`b̄ = Pᵀ M f`, the Galerkin restriction of the load; it conserves `Σ b = Σ M f` since the rows
//! of `P` sum to one). The tests and the example use the sample (the paper solves `L̄ x̄ = M̄ f̄`
//! for a function `f`).
//!
//! # Geodesic distance (item 3, Fig. 22, 23)
//!
//! A tracked point (usually a fine vertex) is mapped to its coarse location by
//! [`CoarseApplications::source_of`]: a coarse vertex or a point of a coarse face.
//!
//! * Exactly: [`CoarseApplications::exact_distance`] runs
//!   `ddg_intrinsic_apps::exact::exact_distance_complex` (window propagation, generalized to
//!   Δ-complexes and face-point sources, S20) on the coarse complex.
//! * By heat: [`CoarseHeat`] **reimplements the three steps of the heat method** (Crane, Weischedel,
//!   Wardetzky 2013) for a Δ-complex with lengths, because `ddg-intrinsic-apps`' heat solver needs
//!   an `IntrinsicTriangulation` of a mesh with positions, which the coarse mesh does not have. The
//!   operators of item 1 are used; the gradient and divergence are per face from the planar layout
//!   of its lengths. `(M̄ + t L̄) u = δ`, `X = −∇u/|∇u|`, `L̄ φ = −∇·X` (our `L` is the positive
//!   semidefinite one, so `Δ = −L`), `t = (time factor) × (mean edge length)²`, shifted so that the
//!   mean of `φ` at the sources is zero. Both factorizations are computed once per [`CoarseHeat`].
//!
//! The coarse distances are prolongated, `d̃ = P d̄`.
//!
//! # All pairs (item 4, Fig. 24), Voronoi (item 5, Fig. 25)
//!
//! [`CoarseApplications::all_pairs`] runs the exact distance from every coarse vertex, in parallel
//! over the sources on an explicit [`ThreadPool`], each result stored at its own row (so the matrix
//! is bitwise independent of the thread count). [`CoarseApplications::prolongate_rows`] evaluates
//! rows of `D̃ = P D̄ Pᵀ` (the dense matrix is `|V|²` numbers: [`CoarseApplications::prolongate_dense`]
//! for small meshes). [`CoarseApplications::voronoi`] labels each fine vertex by the nearest of the
//! given sources, from the prolongated distances to each source (ties to the lower source index);
//! the labels of the coarse vertices come with it.
//!
//! # Exponential map (item 6, Sec. 8.5.4, S21)
//!
//! [`CoarseApplications::exp_map_face`] and [`CoarseApplications::exp_map_vertex`] trace the
//! straightest geodesic of a given length on the coarse mesh by unfolding the faces one after the
//! other along the ray. A face has the planar frame of [`triangle_layout`] (corner 0 at the origin,
//! corner 1 on the positive x-axis, counter-clockwise); a tangent vector at a face point is given in
//! this frame. A direction at a vertex is an angle `θ` in `[0, Θ)` (`[0, Θ]` on the boundary), the
//! *unnormalized* angle of the angular coordinates of SPEC §2.2 (`θ = φ/s`, `s = 2π/Θ` inside,
//! `π/Θ` on the boundary), counter-clockwise from the vertex's `φ = 0` direction. When the ray
//! hits a vertex (within a relative `1e-10` of a side), it continues "straight" by the convention of
//! Polthier and Schmies (1998): the outgoing direction divides the vertex's total angle `Θ` equally
//! on both sides of the incoming ray, i.e. `θ_out = θ_back + Θ/2`. The ray stops at the boundary.
//! The result is a coarse point (the map is fine → coarse, so no fine point is produced, S21), the
//! traced length, the faces crossed and the final direction in the frame of the final face
//! (reversible: tracing the opposite direction for the traced length returns to the start).
//!
//! # A stand-in coarsening
//!
//! [`greedy_coarsening`] removes vertices in the order of the smallest memoryless cost `m_i Σ_j α_ij ℓ_ij`
//! (SPEC §4.4, `α_ij = (K̃_j − K_j)/K_i` clamped to `[0, ∞)`, the curvature `K_i` redistributed by the
//! flattening, with the current lumped mass) — a cheap order that coarsens flat regions first, to make
//! coarse meshes for the tests and the example before `simplify` (ICE) is merged. **Numbers measured
//! with it are not those of the ICE order.** Only the types `CoarseMesh`/`CoarseOutput` are used by
//! the applications, so `simplify`'s output plugs in unchanged.
//!
//! # Relative errors
//!
//! [`DistanceError`] reports for an approximation `d̃` of `d` over the fine vertices: the relative L1
//! error `Σ|d̃ − d| / Σ d`, the relative L2 error, and the maximum error over the maximum `d`.

use std::collections::BinaryHeap;
use std::f64::consts::{PI, TAU};
use std::fmt;

use ddg_core::geometry::corner_cotan;
use ddg_core::sparse::{
    Factorization, Method, Parallelism, SparseError, SparseMatrix, Symbolic, ThreadPool,
};
use ddg_intrinsic::Complex;
use ddg_intrinsic::formulas::triangle_layout;
use ddg_intrinsic_apps::exact::{ExactError, ExactOptions, ExactSource, exact_distance_complex};

use crate::error::IemError;
use crate::mesh::CoarseMesh;
use crate::output::{CoarseOutput, CoarsePoint};
use crate::prolongation::ProlongationError;

/// What goes wrong in the applications.
#[derive(Debug, Clone, PartialEq)]
pub enum AppError {
    /// The coarse mesh could not be compacted (a removal is pending).
    Mesh(IemError),
    /// The prolongation could not be built.
    Prolongation(ProlongationError),
    /// A sparse operation failed (message of the underlying error, with what was being done).
    Sparse(String),
    /// The exact distance failed.
    Exact(ExactError),
    /// An argument is invalid (message names it).
    Invalid(String),
    /// A face of the coarse mesh has no planar layout (lengths violating the triangle inequality).
    Degenerate(String),
    /// The exponential map did not terminate within the step limit.
    StepLimit {
        /// The limit.
        limit: usize,
    },
}

impl fmt::Display for AppError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AppError::Mesh(e) => write!(f, "{e}"),
            AppError::Prolongation(e) => write!(f, "{e}"),
            AppError::Sparse(m) => write!(f, "sparse linear algebra failed: {m}"),
            AppError::Exact(e) => write!(f, "exact geodesic distance failed: {e}"),
            AppError::Invalid(m) => write!(f, "invalid argument: {m}"),
            AppError::Degenerate(m) => write!(f, "degenerate coarse mesh: {m}"),
            AppError::StepLimit { limit } => write!(
                f,
                "the exponential map crossed more than {limit} faces without finishing"
            ),
        }
    }
}

impl std::error::Error for AppError {}

impl From<ExactError> for AppError {
    fn from(e: ExactError) -> Self {
        AppError::Exact(e)
    }
}

fn sparse_err(what: &str, e: SparseError) -> AppError {
    AppError::Sparse(format!("{what}: {e}"))
}

// ---------------------------------------------------------------------------------------------
// Operators on the compacted coarse mesh
// ---------------------------------------------------------------------------------------------

/// The cotangent of the corner at slot `j` of a face with side lengths `l` (`l[k]` from slot `k`
/// to `k + 1`), `None` for a degenerate face.
fn slot_cotan(l: [f64; 3], j: usize) -> Option<f64> {
    corner_cotan(l[(j + 1) % 3], l[(j + 2) % 3], l[j])
}

impl CoarseOutput {
    /// The cotan Laplacian of the coarse Δ-complex from its lengths: `n̄ × n̄`, symmetric positive
    /// semidefinite with zero row sums (see the [module documentation](self)).
    ///
    /// [formula: cotan-laplacian]
    ///
    /// # Errors
    /// [`AppError::Degenerate`] naming the face if its lengths do not form a triangle.
    pub fn cotan_laplacian(&self) -> Result<SparseMatrix, AppError> {
        let n = self.vertices.len();
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(6 * self.faces.len() + n);
        let mut diagonal = vec![0.0; n];
        for (f, face) in self.faces.iter().enumerate() {
            let l = self.face_lengths[f];
            for k in 0..3 {
                let (i, j) = (face[k], face[(k + 1) % 3]);
                if i == j {
                    continue;
                }
                // The halfedge `k` is opposite the corner at slot `k + 2`.
                let cot = slot_cotan(l, (k + 2) % 3).ok_or_else(|| {
                    AppError::Degenerate(format!("face {f} has side lengths {l:?}"))
                })?;
                let w = 0.5 * cot;
                triplets.push((i, j, -w));
                triplets.push((j, i, -w));
                diagonal[i as usize] += w;
                diagonal[j as usize] += w;
            }
        }
        for (v, d) in diagonal.into_iter().enumerate() {
            triplets.push((v as u32, v as u32, d));
        }
        SparseMatrix::from_triplets(n, n, &triplets).map_err(|e| sparse_err("Laplacian", e))
    }

    /// The area of every coarse face (Heron, stable form).
    ///
    /// # Errors
    /// [`AppError::Degenerate`] naming the face.
    pub fn face_areas(&self) -> Result<Vec<f64>, AppError> {
        self.face_lengths
            .iter()
            .enumerate()
            .map(|(f, l)| {
                ddg_core::geometry::triangle_area(l[0], l[1], l[2])
                    .ok_or_else(|| AppError::Degenerate(format!("face {f} has side lengths {l:?}")))
            })
            .collect()
    }

    /// The lumped masses `⅓ Σ A_f` of the coarse vertices (a face that uses a vertex twice counts
    /// twice for it).
    ///
    /// [formula: lumped-mass]
    ///
    /// # Errors
    /// [`AppError::Degenerate`] naming the face.
    pub fn lumped_mass(&self) -> Result<Vec<f64>, AppError> {
        let areas = self.face_areas()?;
        let mut mass = vec![0.0; self.vertices.len()];
        for (f, face) in self.faces.iter().enumerate() {
            for &v in face {
                mass[v as usize] += areas[f] / 3.0;
            }
        }
        Ok(mass)
    }
}

// ---------------------------------------------------------------------------------------------
// Poisson
// ---------------------------------------------------------------------------------------------

/// What a Poisson solve did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PoissonReport {
    /// Unknowns solved for.
    pub n_free: usize,
    /// Fixed (Dirichlet) vertices; 0 means the Neumann/mean-zero problem.
    pub n_fixed: usize,
    /// `Σ b` before the compatibility projection (Neumann only, else 0): `b` was changed by
    /// `−(Σ b) m_i / Σ m`.
    pub rhs_shift: f64,
    /// `‖L x − b‖_∞` over the free rows (Dirichlet) or all rows (Neumann, after the projection).
    pub residual_max: f64,
}

/// A Poisson solver `L x = b` with a Cholesky factorization computed once (any `L`, `M`; used for
/// the coarse and the fine problem alike). See the [module documentation](self).
pub struct PoissonSolver {
    n: usize,
    free: Vec<u32>,
    fixed: Vec<u32>,
    neumann: bool,
    l: SparseMatrix,
    l_free_fixed: SparseMatrix,
    factor: Factorization,
    mass: Vec<f64>,
}

impl fmt::Debug for PoissonSolver {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PoissonSolver")
            .field("n", &self.n)
            .field("n_fixed", &self.fixed.len())
            .field("neumann", &self.neumann)
            .finish()
    }
}

impl PoissonSolver {
    /// Factors `L` with the vertices `fixed` held (Dirichlet), or, if `fixed` is empty, with
    /// vertex 0 pinned internally for the mean-zero problem.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a bad size, an out-of-range or repeated fixed vertex;
    /// [`AppError::Sparse`] if the free block is not positive definite (a component without a
    /// fixed vertex).
    pub fn new(l: &SparseMatrix, mass: &[f64], fixed: &[u32]) -> Result<PoissonSolver, AppError> {
        let n = l.nrows();
        if l.ncols() != n || mass.len() != n {
            return Err(AppError::Invalid(format!(
                "Laplacian is {} x {} and the mass has {} entries",
                l.nrows(),
                l.ncols(),
                mass.len()
            )));
        }
        if n == 0 {
            return Err(AppError::Invalid("the mesh has no vertices".into()));
        }
        let neumann = fixed.is_empty();
        let held: Vec<u32> = if neumann { vec![0] } else { fixed.to_vec() };
        let mut is_held = vec![false; n];
        for &v in &held {
            if v as usize >= n {
                return Err(AppError::Invalid(format!(
                    "fixed vertex {v} is out of range (the mesh has {n} vertices)"
                )));
            }
            if is_held[v as usize] {
                return Err(AppError::Invalid(format!(
                    "fixed vertex {v} is listed twice"
                )));
            }
            is_held[v as usize] = true;
        }
        let free: Vec<u32> = (0..n as u32).filter(|&v| !is_held[v as usize]).collect();
        let a = l
            .submatrix(&free, &free)
            .map_err(|e| sparse_err("free block", e))?;
        let l_free_fixed = l
            .submatrix(&free, &held)
            .map_err(|e| sparse_err("coupling block", e))?;
        let factor = if free.is_empty() {
            // Nothing to solve for: a dummy 1 x 1 factorization is never used.
            Symbolic::analyze(
                &SparseMatrix::from_diagonal(&[1.0]).map_err(|e| sparse_err("dummy", e))?,
                Method::Llt,
            )
            .and_then(|s| {
                s.factor(
                    &SparseMatrix::from_diagonal(&[1.0]).expect("1 x 1"),
                    Parallelism::Sequential,
                )
            })
            .map_err(|e| sparse_err("dummy factorization", e))?
        } else {
            Symbolic::analyze(&a, Method::Llt)
                .and_then(|s| s.factor(&a, Parallelism::Sequential))
                .map_err(|e| {
                    sparse_err(
                        "Cholesky of the free block (is every component constrained?)",
                        e,
                    )
                })?
        };
        Ok(PoissonSolver {
            n,
            free,
            fixed: held,
            neumann,
            l: l.clone(),
            l_free_fixed,
            factor,
            mass: mass.to_vec(),
        })
    }

    /// The number of unknowns of the full problem.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Solves `L x = M f` (data `f`, load `b = M f`); see [`solve_load`](Self::solve_load).
    ///
    /// # Errors
    /// As [`solve_load`](Self::solve_load).
    pub fn solve(&self, f: &[f64], values: &[f64]) -> Result<(Vec<f64>, PoissonReport), AppError> {
        if f.len() != self.n {
            return Err(AppError::Invalid(format!(
                "data has {} entries for {} vertices",
                f.len(),
                self.n
            )));
        }
        let b: Vec<f64> = f.iter().zip(&self.mass).map(|(f, m)| f * m).collect();
        self.solve_load(&b, values)
    }

    /// Solves `L x = b`. `values` are the Dirichlet values at the `fixed` vertices, in their order
    /// (ignored, and may be empty, in the Neumann case). In the Neumann case `b` is projected to
    /// sum zero ([`PoissonReport::rhs_shift`]) and `x` has `M`-weighted mean zero.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a wrong length or a non-finite entry; [`AppError::Sparse`].
    pub fn solve_load(
        &self,
        b: &[f64],
        values: &[f64],
    ) -> Result<(Vec<f64>, PoissonReport), AppError> {
        if b.len() != self.n {
            return Err(AppError::Invalid(format!(
                "load has {} entries for {} vertices",
                b.len(),
                self.n
            )));
        }
        if !self.neumann && values.len() != self.fixed.len() {
            return Err(AppError::Invalid(format!(
                "{} Dirichlet values for {} fixed vertices",
                values.len(),
                self.fixed.len()
            )));
        }
        if let Some(i) = b.iter().position(|x| !x.is_finite()) {
            return Err(AppError::Invalid(format!("load entry {i} is not finite")));
        }
        let mut b = b.to_vec();
        let mut rhs_shift = 0.0;
        let mut xf = vec![0.0; self.fixed.len()];
        if self.neumann {
            let total_mass: f64 = self.mass.iter().sum();
            rhs_shift = b.iter().sum::<f64>();
            for (bi, m) in b.iter_mut().zip(&self.mass) {
                *bi -= rhs_shift * m / total_mass;
            }
        } else {
            if let Some(i) = values.iter().position(|x| !x.is_finite()) {
                return Err(AppError::Invalid(format!(
                    "Dirichlet value {i} is not finite"
                )));
            }
            xf.copy_from_slice(values);
        }
        let mut rhs: Vec<f64> = self.free.iter().map(|&v| b[v as usize]).collect();
        if !self.free.is_empty() {
            let mut coupling = vec![0.0; self.free.len()];
            self.l_free_fixed
                .mul_vec(&xf, &mut coupling)
                .map_err(|e| sparse_err("coupling", e))?;
            for (r, c) in rhs.iter_mut().zip(&coupling) {
                *r -= c;
            }
            self.factor
                .solve_in_place(&mut rhs, Parallelism::Sequential)
                .map_err(|e| sparse_err("solve", e))?;
        }
        let mut x = vec![0.0; self.n];
        for (k, &v) in self.free.iter().enumerate() {
            x[v as usize] = rhs[k];
        }
        for (k, &v) in self.fixed.iter().enumerate() {
            x[v as usize] = xf[k];
        }
        if self.neumann {
            let total_mass: f64 = self.mass.iter().sum();
            let mean = x.iter().zip(&self.mass).map(|(x, m)| x * m).sum::<f64>() / total_mass;
            for xi in &mut x {
                *xi -= mean;
            }
        }
        let mut lx = vec![0.0; self.n];
        self.l
            .mul_vec(&x, &mut lx)
            .map_err(|e| sparse_err("residual", e))?;
        let mut is_fixed = vec![false; self.n];
        if !self.neumann {
            for &v in &self.fixed {
                is_fixed[v as usize] = true;
            }
        }
        let residual_max = (0..self.n)
            .filter(|&v| !is_fixed[v])
            .map(|v| (lx[v] - b[v]).abs())
            .fold(0.0, f64::max);
        let report = PoissonReport {
            n_free: self.free.len(),
            n_fixed: if self.neumann { 0 } else { self.fixed.len() },
            rhs_shift,
            residual_max,
        };
        Ok((x, report))
    }
}

// ---------------------------------------------------------------------------------------------
// The applications object
// ---------------------------------------------------------------------------------------------

/// A dense row-major matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    /// Number of rows.
    pub n_rows: usize,
    /// Number of columns.
    pub n_cols: usize,
    /// Row-major entries, `data[r * n_cols + c]`.
    pub data: Vec<f64>,
}

impl DenseMatrix {
    /// Entry `(r, c)`.
    pub fn get(&self, r: usize, c: usize) -> f64 {
        self.data[r * self.n_cols + c]
    }

    /// Row `r`.
    pub fn row(&self, r: usize) -> &[f64] {
        &self.data[r * self.n_cols..(r + 1) * self.n_cols]
    }
}

/// Approximation errors of a fine-vertex quantity against a reference (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceError {
    /// `Σ|d̃ − d| / Σ|d|`.
    pub rel_l1: f64,
    /// `‖d̃ − d‖₂ / ‖d‖₂`.
    pub rel_l2: f64,
    /// `max|d̃ − d| / max|d|`.
    pub max_over_max: f64,
}

impl DistanceError {
    /// The errors of `approx` against `reference`.
    ///
    /// # Errors
    /// [`AppError::Invalid`] if the lengths differ or a reference value is not finite.
    pub fn of(approx: &[f64], reference: &[f64]) -> Result<DistanceError, AppError> {
        if approx.len() != reference.len() {
            return Err(AppError::Invalid(format!(
                "{} values against {} reference values",
                approx.len(),
                reference.len()
            )));
        }
        if reference.iter().any(|x| !x.is_finite()) {
            return Err(AppError::Invalid(
                "a reference value is not finite (unreached vertex)".into(),
            ));
        }
        let (mut e1, mut e2, mut r1, mut r2, mut emax, mut rmax) =
            (0.0, 0.0, 0.0, 0.0, 0.0f64, 0.0f64);
        for (a, r) in approx.iter().zip(reference) {
            let e = (a - r).abs();
            e1 += e;
            e2 += e * e;
            r1 += r.abs();
            r2 += r * r;
            emax = emax.max(e);
            rmax = rmax.max(r.abs());
        }
        let ratio = |a: f64, b: f64| if b > 0.0 { a / b } else { 0.0 };
        Ok(DistanceError {
            rel_l1: ratio(e1, r1),
            rel_l2: ratio(e2.sqrt(), r2.sqrt()),
            max_over_max: ratio(emax, rmax),
        })
    }
}

/// The result of [`CoarseApplications::voronoi`].
#[derive(Debug, Clone, PartialEq)]
pub struct Voronoi {
    /// Per fine vertex: the index (into the sources) of the nearest source by `P d̄`.
    pub fine_labels: Vec<u32>,
    /// Per coarse vertex: the nearest source by `d̄`.
    pub coarse_labels: Vec<u32>,
    /// Per fine vertex: the prolongated distance to that nearest source.
    pub fine_distance: Vec<f64>,
}

/// The coarse mesh with everything the applications need: the compacted output, the Δ-complex and
/// its lengths (for the exact distance), the prolongation, the operators of item 1.
#[derive(Debug, Clone)]
pub struct CoarseApplications {
    output: CoarseOutput,
    complex: Complex,
    lengths: Vec<f64>,
    prolongation: SparseMatrix,
    laplacian: SparseMatrix,
    mass: Vec<f64>,
    layouts: Vec<[[f64; 2]; 3]>,
    /// Per fine vertex: its (coarse vertex, weight) entries of `P`.
    p_rows: Vec<Vec<(u32, f64)>>,
    /// Per coarse vertex: the (face, slot) corners.
    corner_ptr: Vec<usize>,
    corners: Vec<(u32, u8)>,
}

impl CoarseApplications {
    /// Builds the applications on a coarse mesh (compacts it; builds `P`, `L̄`, `M̄`).
    ///
    /// # Errors
    /// [`AppError::Mesh`] while a removal is pending; [`AppError::Degenerate`];
    /// [`AppError::Prolongation`].
    pub fn new(mesh: &CoarseMesh) -> Result<CoarseApplications, AppError> {
        let output = mesh.compact().map_err(AppError::Mesh)?;
        let prolongation = output
            .scalar_prolongation()
            .map_err(AppError::Prolongation)?;
        let laplacian = output.cotan_laplacian()?;
        let mass = output.lumped_mass()?;
        let layouts = output
            .face_lengths
            .iter()
            .enumerate()
            .map(|(f, l)| {
                triangle_layout(*l)
                    .ok_or_else(|| AppError::Degenerate(format!("face {f} has side lengths {l:?}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let p_rows = output
            .points
            .iter()
            .take(output.n_fine_vertices)
            .map(|p| match *p {
                CoarsePoint::Vertex(v) => vec![(v, 1.0)],
                CoarsePoint::Face { face, bary } => (0..3)
                    .map(|k| (output.faces[face as usize][k], bary[k]))
                    .collect(),
            })
            .collect();
        let (corner_ptr, corners) = output.vertex_slots();
        Ok(CoarseApplications {
            lengths: mesh.edge_lengths().to_vec(),
            complex: mesh.complex().clone(),
            output,
            prolongation,
            laplacian,
            mass,
            layouts,
            p_rows,
            corner_ptr,
            corners,
        })
    }

    /// The compacted coarse mesh.
    pub fn output(&self) -> &CoarseOutput {
        &self.output
    }

    /// The coarse Δ-complex and its edge lengths.
    pub fn complex(&self) -> (&Complex, &[f64]) {
        (&self.complex, &self.lengths)
    }

    /// The scalar prolongation `P` (`|V| × |V̄|`).
    pub fn prolongation(&self) -> &SparseMatrix {
        &self.prolongation
    }

    /// The coarse cotan Laplacian `L̄`.
    pub fn laplacian(&self) -> &SparseMatrix {
        &self.laplacian
    }

    /// The coarse lumped mass `M̄` (diagonal).
    pub fn mass(&self) -> &[f64] {
        &self.mass
    }

    /// `x = P x̄`: coarse values to the fine vertices.
    ///
    /// # Errors
    /// [`AppError::Invalid`] if `coarse.len() != |V̄|`.
    pub fn prolong(&self, coarse: &[f64]) -> Result<Vec<f64>, AppError> {
        if coarse.len() != self.output.n_vertices() {
            return Err(AppError::Invalid(format!(
                "{} coarse values for {} coarse vertices",
                coarse.len(),
                self.output.n_vertices()
            )));
        }
        Ok(self
            .p_rows
            .iter()
            .map(|r| r.iter().map(|&(c, w)| w * coarse[c as usize]).sum())
            .collect())
    }

    /// Pointwise restriction: the value of the fine data `f` at the fine vertex of every coarse
    /// vertex.
    ///
    /// # Errors
    /// [`AppError::Invalid`] if `f.len() != |V|`.
    pub fn sample(&self, f: &[f64]) -> Result<Vec<f64>, AppError> {
        if f.len() != self.output.n_fine_vertices {
            return Err(AppError::Invalid(format!(
                "{} fine values for {} fine vertices",
                f.len(),
                self.output.n_fine_vertices
            )));
        }
        Ok(self
            .output
            .vertices
            .iter()
            .map(|&v| f[v as usize])
            .collect())
    }

    /// Galerkin restriction of a fine load: `b̄ = Pᵀ b` (with `b = M f` for data `f`; sums are
    /// conserved).
    ///
    /// # Errors
    /// [`AppError::Invalid`] if `b.len() != |V|`.
    pub fn restrict_load(&self, b: &[f64]) -> Result<Vec<f64>, AppError> {
        if b.len() != self.output.n_fine_vertices {
            return Err(AppError::Invalid(format!(
                "{} fine values for {} fine vertices",
                b.len(),
                self.output.n_fine_vertices
            )));
        }
        let mut out = vec![0.0; self.output.n_vertices()];
        for (row, bi) in self.p_rows.iter().zip(b) {
            for &(c, w) in row {
                out[c as usize] += w * bi;
            }
        }
        Ok(out)
    }

    /// A Poisson solver on the coarse mesh with the coarse vertices `fixed` held.
    ///
    /// # Errors
    /// As [`PoissonSolver::new`].
    pub fn poisson_solver(&self, fixed: &[u32]) -> Result<PoissonSolver, AppError> {
        PoissonSolver::new(&self.laplacian, &self.mass, fixed)
    }

    /// Poisson (Fig. 21): solves `L̄ x̄ = M̄ f̄` for coarse data `f̄` with Dirichlet values at the
    /// coarse vertices `fixed` (none: the mean-zero problem), and returns `(x̄, x = P x̄, report)`.
    ///
    /// # Errors
    /// As [`PoissonSolver::new`] and [`PoissonSolver::solve`].
    pub fn poisson(
        &self,
        fbar: &[f64],
        fixed: &[u32],
        values: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>, PoissonReport), AppError> {
        let solver = self.poisson_solver(fixed)?;
        let (xbar, report) = solver.solve(fbar, values)?;
        let x = self.prolong(&xbar)?;
        Ok((xbar, x, report))
    }

    // ----- geodesic distance -------------------------------------------------------------

    /// The coarse location of a tracked point as a source of the exact distance.
    ///
    /// # Errors
    /// [`AppError::Invalid`] if `p` is not a tracked point.
    pub fn source_of(&self, p: usize) -> Result<ExactSource, AppError> {
        match self.output.points.get(p) {
            None => Err(AppError::Invalid(format!(
                "point {p} is not tracked (there are {})",
                self.output.points.len()
            ))),
            Some(&CoarsePoint::Vertex(v)) => Ok(ExactSource::Vertex(v as usize)),
            Some(&CoarsePoint::Face { face, bary }) => Ok(ExactSource::Point {
                face: face as usize,
                bary,
            }),
        }
    }

    /// The exact geodesic distance on the coarse Δ-complex from the nearest of `sources` to every
    /// coarse vertex (S20).
    ///
    /// # Errors
    /// [`AppError::Exact`].
    pub fn exact_distance(
        &self,
        sources: &[ExactSource],
        options: &ExactOptions,
    ) -> Result<Vec<f64>, AppError> {
        Ok(exact_distance_complex(&self.complex, &self.lengths, sources, options)?.distances)
    }

    /// Single-source geodesic distance by the exact algorithm: the tracked point `source` (a fine
    /// vertex, usually) mapped to the coarse mesh, `d̄` there, `d̃ = P d̄` on the fine vertices.
    /// Returns `(d̄, d̃)`.
    ///
    /// # Errors
    /// As [`source_of`](Self::source_of) and [`exact_distance`](Self::exact_distance).
    pub fn geodesic_exact(
        &self,
        source: usize,
        options: &ExactOptions,
    ) -> Result<(Vec<f64>, Vec<f64>), AppError> {
        let s = self.source_of(source)?;
        let dbar = self.exact_distance(&[s], options)?;
        let d = self.prolong(&dbar)?;
        Ok((dbar, d))
    }

    /// Single-source geodesic distance by the heat method on the coarse mesh (see
    /// [`CoarseHeat`]; the factorizations are made for this one call: build a [`CoarseHeat`] to
    /// reuse them). Returns `(d̄, d̃)`.
    ///
    /// # Errors
    /// As [`CoarseHeat::new`] and [`CoarseHeat::distance`].
    pub fn geodesic_heat(
        &self,
        source: usize,
        time_factor: f64,
    ) -> Result<(Vec<f64>, Vec<f64>), AppError> {
        let s = self.source_of(source)?;
        let heat = CoarseHeat::new(self, time_factor)?;
        let dbar = heat.distance(&[s])?;
        let d = self.prolong(&dbar)?;
        Ok((dbar, d))
    }

    /// All-pairs coarse distances `D̄` (`|V̄| × |V̄|`, row `i` = the exact distances from coarse
    /// vertex `i`), the sources distributed over the `pool`'s threads (deterministic: row `i`
    /// depends on `i` alone).
    ///
    /// # Errors
    /// [`AppError::Exact`] (the first failing row).
    pub fn all_pairs(
        &self,
        pool: &ThreadPool,
        options: &ExactOptions,
    ) -> Result<DenseMatrix, AppError> {
        let n = self.output.n_vertices();
        let rows = pool.map(n, |i| {
            exact_distance_complex(
                &self.complex,
                &self.lengths,
                &[ExactSource::Vertex(i)],
                options,
            )
            .map(|r| r.distances)
        });
        let mut data = Vec::with_capacity(n * n);
        for r in rows {
            data.extend(r?);
        }
        Ok(DenseMatrix {
            n_rows: n,
            n_cols: n,
            data,
        })
    }

    /// Rows `rows` (fine vertices) of `D̃ = P D̄ Pᵀ`, as a `rows.len() × |V|` matrix.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a row out of range or a `D̄` of the wrong size.
    pub fn prolongate_rows(
        &self,
        dbar: &DenseMatrix,
        rows: &[usize],
    ) -> Result<DenseMatrix, AppError> {
        let nc = self.output.n_vertices();
        if dbar.n_rows != nc || dbar.n_cols != nc {
            return Err(AppError::Invalid(format!(
                "D̄ is {} x {} for {} coarse vertices",
                dbar.n_rows, dbar.n_cols, nc
            )));
        }
        let nf = self.output.n_fine_vertices;
        let mut data = Vec::with_capacity(rows.len() * nf);
        let mut y = vec![0.0; nc];
        for &r in rows {
            if r >= nf {
                return Err(AppError::Invalid(format!(
                    "row {r} is out of range ({nf} fine vertices)"
                )));
            }
            y.iter_mut().for_each(|v| *v = 0.0);
            for &(c, w) in &self.p_rows[r] {
                for (yk, d) in y.iter_mut().zip(dbar.row(c as usize)) {
                    *yk += w * d;
                }
            }
            for prow in &self.p_rows {
                data.push(prow.iter().map(|&(c, w)| w * y[c as usize]).sum());
            }
        }
        Ok(DenseMatrix {
            n_rows: rows.len(),
            n_cols: nf,
            data,
        })
    }

    /// The dense `D̃ = P D̄ Pᵀ` (`|V| × |V|`: memory `8 |V|²` bytes), rows in parallel on `pool`.
    ///
    /// # Errors
    /// As [`prolongate_rows`](Self::prolongate_rows).
    pub fn prolongate_dense(
        &self,
        dbar: &DenseMatrix,
        pool: &ThreadPool,
    ) -> Result<DenseMatrix, AppError> {
        let nf = self.output.n_fine_vertices;
        let chunks = pool.map(nf, |r| self.prolongate_rows(dbar, &[r]));
        let mut data = Vec::with_capacity(nf * nf);
        for c in chunks {
            data.extend(c?.data);
        }
        Ok(DenseMatrix {
            n_rows: nf,
            n_cols: nf,
            data,
        })
    }

    /// Geodesic Voronoi labels (Fig. 25): the exact coarse distance from each of the tracked
    /// points `sources` (one run per source, in parallel on `pool`), prolonged, and the nearest
    /// source of every fine vertex (ties to the lower index).
    ///
    /// # Errors
    /// [`AppError::Invalid`] for no sources; [`AppError::Exact`].
    pub fn voronoi(
        &self,
        sources: &[usize],
        pool: &ThreadPool,
        options: &ExactOptions,
    ) -> Result<Voronoi, AppError> {
        if sources.is_empty() {
            return Err(AppError::Invalid("no Voronoi sources".into()));
        }
        let exact: Vec<ExactSource> = sources
            .iter()
            .map(|&s| self.source_of(s))
            .collect::<Result<_, _>>()?;
        let runs = pool.map(exact.len(), |k| self.exact_distance(&exact[k..=k], options));
        let mut coarse_d = Vec::with_capacity(sources.len());
        for r in runs {
            coarse_d.push(r?);
        }
        let nc = self.output.n_vertices();
        let mut coarse_labels = vec![0u32; nc];
        for v in 0..nc {
            let mut best = 0;
            for k in 1..coarse_d.len() {
                if coarse_d[k][v] < coarse_d[best][v] {
                    best = k;
                }
            }
            coarse_labels[v] = best as u32;
        }
        let fine_d: Vec<Vec<f64>> = coarse_d
            .iter()
            .map(|d| self.prolong(d))
            .collect::<Result<_, _>>()?;
        let nf = self.output.n_fine_vertices;
        let mut fine_labels = vec![0u32; nf];
        let mut fine_distance = vec![0.0; nf];
        for v in 0..nf {
            let mut best = 0;
            for k in 1..fine_d.len() {
                if fine_d[k][v] < fine_d[best][v] {
                    best = k;
                }
            }
            fine_labels[v] = best as u32;
            fine_distance[v] = fine_d[best][v];
        }
        Ok(Voronoi {
            fine_labels,
            coarse_labels,
            fine_distance,
        })
    }
}

// ---------------------------------------------------------------------------------------------
// Heat method on the Δ-complex
// ---------------------------------------------------------------------------------------------

/// The heat method (Crane, Weischedel, Wardetzky 2013) on the coarse Δ-complex with lengths: the
/// factorizations of `M̄ + t L̄` and of `L̄` are made once, each call to
/// [`distance`](Self::distance) is two solves and a pass over the faces. See the
/// [module documentation](self) for the discretization.
pub struct CoarseHeat<'a> {
    app: &'a CoarseApplications,
    heat: Factorization,
    poisson: PoissonSolver,
    /// The diffusion time `t`.
    pub time: f64,
}

impl fmt::Debug for CoarseHeat<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CoarseHeat")
            .field("time", &self.time)
            .finish()
    }
}

impl<'a> CoarseHeat<'a> {
    /// Factors the two systems with `t = time_factor × h²`, `h` the mean edge length (the
    /// paper's recommendation is a factor of 1).
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a non-positive factor; [`AppError::Sparse`].
    pub fn new(app: &'a CoarseApplications, time_factor: f64) -> Result<CoarseHeat<'a>, AppError> {
        if !(time_factor.is_finite() && time_factor > 0.0) {
            return Err(AppError::Invalid(format!(
                "the heat time factor must be positive, not {time_factor}"
            )));
        }
        let h = app.output.edge_lengths.iter().sum::<f64>()
            / app.output.edge_lengths.len().max(1) as f64;
        let time = time_factor * h * h;
        let mut a = app.laplacian.clone();
        for v in a.values_mut() {
            *v *= time;
        }
        a.add_to_diagonal(&app.mass)
            .map_err(|e| sparse_err("M + tL", e))?;
        let heat = Symbolic::analyze(&a, Method::Llt)
            .and_then(|s| s.factor(&a, Parallelism::Sequential))
            .map_err(|e| sparse_err("Cholesky of M + tL", e))?;
        let poisson = app.poisson_solver(&[])?;
        Ok(CoarseHeat {
            app,
            heat,
            poisson,
            time,
        })
    }

    /// The heat-method distance from the nearest of `sources` to every coarse vertex.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for no sources or a bad source; [`AppError::Sparse`].
    pub fn distance(&self, sources: &[ExactSource]) -> Result<Vec<f64>, AppError> {
        let out = &self.app.output;
        let n = out.n_vertices();
        if sources.is_empty() {
            return Err(AppError::Invalid("no heat source".into()));
        }
        // Step 1: heat flow from the sources.
        let mut u = vec![0.0; n];
        for s in sources {
            for (v, w) in self.source_weights(s)? {
                u[v as usize] += w;
            }
        }
        self.heat
            .solve_in_place(&mut u, Parallelism::Sequential)
            .map_err(|e| sparse_err("heat solve", e))?;
        // Step 2: normalized negative gradient per face; step 3: its divergence at the vertices.
        let mut div = vec![0.0; n];
        for (f, face) in out.faces.iter().enumerate() {
            let p = self.app.layouts[f];
            let x = unit_negative_gradient(
                &p,
                [
                    u[face[0] as usize],
                    u[face[1] as usize],
                    u[face[2] as usize],
                ],
            );
            let l = out.face_lengths[f];
            for k in 0..3 {
                let ea = sub(p[(k + 1) % 3], p[k]);
                let eb = sub(p[(k + 2) % 3], p[k]);
                let cot_a = slot_cotan(l, (k + 2) % 3).unwrap_or(0.0);
                let cot_b = slot_cotan(l, (k + 1) % 3).unwrap_or(0.0);
                div[face[k] as usize] += 0.5 * (cot_a * dot(ea, x) + cot_b * dot(eb, x));
            }
        }
        // `L` is positive semidefinite (`Δ = −L`), so `Δ φ = ∇·X` reads `L φ = −∇·X`.
        let b: Vec<f64> = div.iter().map(|d| -d).collect();
        let (mut phi, _) = self.poisson.solve_load(&b, &[])?;
        // Shift: the mean value at the sources is zero.
        let mut shift = 0.0;
        let mut count = 0.0;
        for s in sources {
            for (v, w) in self.source_weights(s)? {
                shift += w * phi[v as usize];
            }
            count += 1.0;
        }
        shift /= count;
        for p in &mut phi {
            *p -= shift;
        }
        Ok(phi)
    }

    /// The vertices and weights (summing to one) that represent a source.
    fn source_weights(&self, s: &ExactSource) -> Result<Vec<(u32, f64)>, AppError> {
        let out = &self.app.output;
        match *s {
            ExactSource::Vertex(v) => {
                if v >= out.n_vertices() {
                    return Err(AppError::Invalid(format!(
                        "source vertex {v} is out of range ({} coarse vertices)",
                        out.n_vertices()
                    )));
                }
                Ok(vec![(v as u32, 1.0)])
            }
            ExactSource::Point { face, bary } => {
                if face >= out.n_faces() {
                    return Err(AppError::Invalid(format!(
                        "source face {face} is out of range ({} coarse faces)",
                        out.n_faces()
                    )));
                }
                Ok((0..3).map(|k| (out.faces[face][k], bary[k])).collect())
            }
        }
    }
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn rot90(a: [f64; 2]) -> [f64; 2] {
    [-a[1], a[0]]
}
fn norm(a: [f64; 2]) -> f64 {
    dot(a, a).sqrt()
}

/// `−∇u/|∇u|` of the linear function with values `u` at the corners of the planar triangle `p`
/// (zero where the gradient vanishes).
fn unit_negative_gradient(p: &[[f64; 2]; 3], u: [f64; 3]) -> [f64; 2] {
    let two_area = cross(sub(p[1], p[0]), sub(p[2], p[0]));
    let mut g = [0.0, 0.0];
    for k in 0..3 {
        let e = rot90(sub(p[(k + 2) % 3], p[(k + 1) % 3]));
        g[0] += u[k] * e[0];
        g[1] += u[k] * e[1];
    }
    g[0] /= two_area;
    g[1] /= two_area;
    let n = norm(g);
    if n > 0.0 && n.is_finite() {
        [-g[0] / n, -g[1] / n]
    } else {
        [0.0, 0.0]
    }
}

// ---------------------------------------------------------------------------------------------
// Exponential map
// ---------------------------------------------------------------------------------------------

/// Why an exponential map stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExpStop {
    /// The requested length was traced.
    Completed,
    /// The ray reached the boundary of the mesh.
    Boundary,
}

/// The result of an exponential map on the coarse mesh (S21).
#[derive(Debug, Clone, PartialEq)]
pub struct ExpMapResult {
    /// Where the geodesic ends.
    pub point: CoarsePoint,
    /// The length traced (the requested one unless it stopped at the boundary).
    pub traced: f64,
    /// Whether it completed.
    pub stop: ExpStop,
    /// The faces visited, in order (the start face first).
    pub faces: Vec<u32>,
    /// Vertices passed through (hit by the ray), in order.
    pub vertices_hit: Vec<u32>,
    /// The unit direction at the end, in the planar frame of the last face of `faces`.
    pub direction: [f64; 2],
}

const VERTEX_TOLERANCE: f64 = 1e-10;

impl CoarseApplications {
    /// The planar frame of a face: corners 0, 1, 2 (counter-clockwise).
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a face out of range.
    pub fn face_layout(&self, face: u32) -> Result<[[f64; 2]; 3], AppError> {
        self.layouts.get(face as usize).copied().ok_or_else(|| {
            AppError::Invalid(format!(
                "face {face} is out of range ({} coarse faces)",
                self.layouts.len()
            ))
        })
    }

    /// The position in the frame of `face` of the point with barycentric coordinates `bary`.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a face out of range.
    pub fn face_position(&self, face: u32, bary: [f64; 3]) -> Result<[f64; 2], AppError> {
        let p = self.face_layout(face)?;
        Ok([
            bary[0] * p[0][0] + bary[1] * p[1][0] + bary[2] * p[2][0],
            bary[0] * p[0][1] + bary[1] * p[1][1] + bary[2] * p[2][1],
        ])
    }

    fn bary_of(&self, face: usize, q: [f64; 2]) -> [f64; 3] {
        let p = &self.layouts[face];
        let a = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        let b0 = cross(sub(p[1], q), sub(p[2], q)) / a;
        let b1 = cross(sub(p[2], q), sub(p[0], q)) / a;
        let b2 = 1.0 - b0 - b1;
        [b0, b1, b2]
    }

    /// Exponential map from the point `bary` of coarse face `face` with the tangent vector
    /// `vector` in the frame of the face (its length is the geodesic length). Zero length returns
    /// the start.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a bad face, non-finite input or barycentric coordinates that are
    /// negative by more than `1e-9` or do not sum to one; [`AppError::StepLimit`].
    pub fn exp_map_face(
        &self,
        face: u32,
        bary: [f64; 3],
        vector: [f64; 2],
    ) -> Result<ExpMapResult, AppError> {
        self.face_layout(face)?;
        let bary_ok = bary.iter().all(|b| b.is_finite() && *b >= -1e-9)
            && (bary.iter().sum::<f64>() - 1.0).abs() <= 1e-9;
        if !bary_ok || !vector.iter().all(|x| x.is_finite()) {
            return Err(AppError::Invalid(format!(
                "barycentric coordinates {bary:?} and vector {vector:?} must be finite, the \
                 coordinates non-negative and summing to one"
            )));
        }
        let length = norm(vector);
        let start = self.face_position(face, bary)?;
        let dir = if length > 0.0 {
            [vector[0] / length, vector[1] / length]
        } else {
            [1.0, 0.0]
        };
        self.trace(face as usize, start, dir, length, Vec::new())
    }

    /// Exponential map from coarse vertex `vertex` in the direction `angle` (unnormalized, from
    /// the vertex's `φ = 0` direction, in `[0, Θ)`, or `[0, Θ]` on the boundary; interior angles
    /// are taken modulo `Θ`) for `length`.
    ///
    /// # Errors
    /// [`AppError::Invalid`] for a bad vertex, a negative or non-finite length or angle, or an
    /// angle outside `[0, Θ]` at a boundary vertex; [`AppError::StepLimit`].
    pub fn exp_map_vertex(
        &self,
        vertex: u32,
        angle: f64,
        length: f64,
    ) -> Result<ExpMapResult, AppError> {
        if vertex as usize >= self.output.n_vertices() {
            return Err(AppError::Invalid(format!(
                "vertex {vertex} is out of range ({} coarse vertices)",
                self.output.n_vertices()
            )));
        }
        if !(length.is_finite() && length >= 0.0 && angle.is_finite()) {
            return Err(AppError::Invalid(format!(
                "length {length} must be finite and non-negative and the angle {angle} finite"
            )));
        }
        let (face, slot, dir) = self.corner_for_angle(vertex, angle)?;
        let start = self.layouts[face][slot];
        let mut r = self.trace(face, start, dir, length, vec![vertex])?;
        // The start vertex is not a "hit".
        r.vertices_hit.remove(0);
        Ok(r)
    }

    /// The unnormalized angle (module docs) of the unit direction `d` (face frame) leaving the
    /// corner at `slot` of `face`, measured like the angular coordinates.
    fn vertex_angle_of(&self, face: usize, slot: usize, d: [f64; 2]) -> f64 {
        let p = &self.layouts[face];
        let v = self.output.faces[face][slot] as usize;
        let u = sub(p[(slot + 1) % 3], p[slot]);
        let a = cross(u, d).atan2(dot(u, d));
        let scale = self.vertex_scale(v);
        self.output.phi[3 * face + slot] / scale + a
    }

    fn vertex_scale(&self, v: usize) -> f64 {
        let theta = self.output.angle_sum[v];
        if self.output.boundary[v] {
            PI / theta
        } else {
            TAU / theta
        }
    }

    /// The corner (face, slot) at `vertex` whose wedge contains the direction `angle`, and the
    /// unit direction in that face's frame.
    fn corner_for_angle(
        &self,
        vertex: u32,
        angle: f64,
    ) -> Result<(usize, usize, [f64; 2]), AppError> {
        let v = vertex as usize;
        let theta = self.output.angle_sum[v];
        let boundary = self.output.boundary[v];
        let scale = self.vertex_scale(v);
        let range = self.corner_ptr[v]..self.corner_ptr[v + 1];
        if range.is_empty() {
            return Err(AppError::Invalid(format!("vertex {vertex} has no face")));
        }
        if boundary && !(-1e-12..=theta + 1e-12).contains(&angle) {
            return Err(AppError::Invalid(format!(
                "angle {angle} is outside the boundary vertex {vertex}'s range [0, {theta}]"
            )));
        }
        let mut fallback = None;
        for &(f, k) in &self.corners[range] {
            let (f, k) = (f as usize, k as usize);
            let start = self.output.phi[3 * f + k] / scale;
            let width = self.output.corner_angles[f][k];
            let a = if boundary {
                angle - start
            } else {
                (angle - start).rem_euclid(theta)
            };
            let dist_to_end = (a - width).abs();
            let hit = a >= -1e-12 && a < width;
            let hit_end = boundary && a >= -1e-12 && dist_to_end <= 1e-12;
            if hit || hit_end {
                let p = &self.layouts[f];
                let u = sub(p[(k + 1) % 3], p[k]);
                let un = norm(u);
                let (c, s) = (a.cos(), a.sin());
                let d = [(u[0] * c - u[1] * s) / un, (u[0] * s + u[1] * c) / un];
                return Ok((f, k, d));
            }
            if fallback.is_none() {
                fallback = Some((f, k, a));
            }
        }
        // Rounding put the angle between two wedges: use the first one.
        let (f, k, _) = fallback.expect("range is non-empty");
        let p = &self.layouts[f];
        let u = sub(p[(k + 1) % 3], p[k]);
        let un = norm(u);
        Ok((f, k, [u[0] / un, u[1] / un]))
    }

    /// The tracing loop. `start` is a point of `face` (in its frame), `dir` a unit direction
    /// pointing into the face (or along it).
    fn trace(
        &self,
        face: usize,
        start: [f64; 2],
        dir: [f64; 2],
        length: f64,
        mut vertices_hit: Vec<u32>,
    ) -> Result<ExpMapResult, AppError> {
        let limit = 100 * self.output.n_faces() + 1000;
        let mut f = face;
        let mut p = start;
        let mut d = dir;
        let mut remaining = length;
        let mut faces = vec![f as u32];
        let mut traced = 0.0;
        let mut steps = 0usize;
        let finish = |app: &Self,
                      f: usize,
                      p: [f64; 2],
                      d: [f64; 2],
                      traced: f64,
                      stop: ExpStop,
                      faces: Vec<u32>,
                      vertices_hit: Vec<u32>| {
            let mut bary = app.bary_of(f, p);
            for b in &mut bary {
                if *b < 0.0 {
                    *b = 0.0;
                }
            }
            let s: f64 = bary.iter().sum();
            for b in &mut bary {
                *b /= s;
            }
            let point = match bary.iter().position(|b| *b >= 1.0 - 1e-12) {
                Some(k) => CoarsePoint::Vertex(app.output.faces[f][k]),
                None => CoarsePoint::Face {
                    face: f as u32,
                    bary,
                },
            };
            ExpMapResult {
                point,
                traced,
                stop,
                faces,
                vertices_hit,
                direction: d,
            }
        };
        loop {
            steps += 1;
            if steps > limit {
                return Err(AppError::StepLimit { limit });
            }
            if remaining <= 0.0 {
                return Ok(finish(
                    self,
                    f,
                    p,
                    d,
                    traced,
                    ExpStop::Completed,
                    faces,
                    vertices_hit,
                ));
            }
            let tri = &self.layouts[f];
            // The exit side: smallest parameter among the sides the ray crosses outward.
            let mut best: Option<(f64, usize)> = None;
            for k in 0..3 {
                let a = tri[k];
                let e = sub(tri[(k + 1) % 3], a);
                let slope = cross(e, d);
                if slope < 0.0 {
                    let c0 = cross(e, sub(p, a)).max(0.0);
                    let s = c0 / -slope;
                    if best.is_none_or(|(bs, _)| s < bs) {
                        best = Some((s, k));
                    }
                }
            }
            let Some((s, k)) = best else {
                // A direction that never leaves the face cannot happen in a triangle.
                return Err(AppError::Degenerate(format!(
                    "the ray in face {f} does not leave it"
                )));
            };
            if s >= remaining {
                let q = [p[0] + remaining * d[0], p[1] + remaining * d[1]];
                traced += remaining;
                return Ok(finish(
                    self,
                    f,
                    q,
                    d,
                    traced,
                    ExpStop::Completed,
                    faces,
                    vertices_hit,
                ));
            }
            let q = [p[0] + s * d[0], p[1] + s * d[1]];
            traced += s;
            remaining -= s;
            // Did the ray hit an end point of the side (a vertex)?
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let side = norm(sub(b, a));
            let t = dot(sub(q, a), sub(b, a)) / (side * side);
            let corner = if t * side <= VERTEX_TOLERANCE * side {
                Some(k)
            } else if (1.0 - t) * side <= VERTEX_TOLERANCE * side {
                Some((k + 1) % 3)
            } else {
                None
            };
            if let Some(slot) = corner {
                let v = self.output.faces[f][slot];
                vertices_hit.push(v);
                if self.output.boundary[v as usize] {
                    let pos = tri[slot];
                    return Ok(finish(
                        self,
                        f,
                        pos,
                        d,
                        traced,
                        ExpStop::Boundary,
                        faces,
                        vertices_hit,
                    ));
                }
                // Straight on: the direction back along the ray plus half the total angle.
                let back = [-d[0], -d[1]];
                let theta_back = self.vertex_angle_of(f, slot, back);
                let theta = self.output.angle_sum[v as usize];
                let out = theta_back + 0.5 * theta;
                let (g, gk, gd) = self.corner_for_angle(v, out)?;
                f = g;
                p = self.layouts[g][gk];
                d = gd;
                faces.push(f as u32);
                continue;
            }
            // Cross the side into the neighbouring face.
            let h = 3 * f + k;
            let Some(tw) = self.output.twin[h] else {
                return Ok(finish(
                    self,
                    f,
                    q,
                    d,
                    traced,
                    ExpStop::Boundary,
                    faces,
                    vertices_hit,
                ));
            };
            let (g, j) = ((tw / 3) as usize, (tw % 3) as usize);
            // f's side runs a -> b; g's side runs from a' = b to b' = a.
            let e = unit(sub(a, b));
            let n = rot90(e);
            let a2 = self.layouts[g][j];
            let b2 = self.layouts[g][(j + 1) % 3];
            let e2 = unit(sub(b2, a2));
            let n2 = rot90(e2);
            let rel = sub(q, b);
            let (x, y) = (dot(rel, e), dot(rel, n));
            p = [a2[0] + x * e2[0] + y * n2[0], a2[1] + x * e2[1] + y * n2[1]];
            let (dx, dy) = (dot(d, e), dot(d, n));
            d = [dx * e2[0] + dy * n2[0], dx * e2[1] + dy * n2[1]];
            f = g;
            faces.push(f as u32);
        }
    }
}

fn unit(a: [f64; 2]) -> [f64; 2] {
    let n = norm(a);
    [a[0] / n, a[1] / n]
}

// ---------------------------------------------------------------------------------------------
// A stand-in coarsening order
// ---------------------------------------------------------------------------------------------

/// A heap entry: the cost (min first), the vertex, and the version of its cost.
#[derive(Debug, Clone, Copy)]
struct Entry {
    cost: f64,
    vertex: u32,
    version: u32,
}

impl PartialEq for Entry {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == std::cmp::Ordering::Equal
    }
}
impl Eq for Entry {}
impl PartialOrd for Entry {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Entry {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        // Reversed: BinaryHeap is a max-heap. Ties to the lower vertex id (deterministic).
        o.cost
            .total_cmp(&self.cost)
            .then(o.vertex.cmp(&self.vertex))
    }
}

/// The memoryless cost `m_i Σ_j α_ij ℓ_ij` of removing fine vertex `v` now (tentatively; the mesh
/// is unchanged afterwards), `None` if it cannot be removed.
fn memoryless_cost(mesh: &mut CoarseMesh, v: u32, mass: f64) -> Option<f64> {
    let prepared = mesh.tentative_removal(v).ok()?;
    let k = prepared.curvature;
    if k.abs() < 1e-12 {
        // A flat vertex: no curvature to move, the cost is the area times the shortest edge scale.
        let l = prepared
            .neighbours
            .iter()
            .map(|n| n.length)
            .fold(f64::INFINITY, f64::min);
        return Some(1e-9 * mass * l);
    }
    let mut cost = 0.0;
    for n in &prepared.neighbours {
        let alpha = ((n.curvature_after - n.curvature_before) / k).max(0.0);
        cost += alpha * n.length;
    }
    Some(mass * cost)
}

/// Removes vertices from `mesh`, cheapest memoryless cost first (lazy queue; the costs of the
/// neighbours of a removed vertex and their neighbours are recomputed), until `target` vertices
/// are left or none can be removed. Returns the number of vertices removed. See the
/// [module documentation](self): a stand-in for `simplify`.
///
/// # Errors
/// [`IemError`] from the mesh (never for a valid mesh without pending removal).
pub fn greedy_coarsening(mesh: &mut CoarseMesh, target: usize) -> Result<usize, IemError> {
    let nf = mesh.n_fine_vertices() as u32;
    let mut version = vec![0u32; nf as usize];
    let mut heap: BinaryHeap<Entry> = BinaryHeap::new();
    let mass_of = |mesh: &CoarseMesh, v: u32| -> f64 {
        mesh.complex_vertex(v)
            .map(|c| {
                let cx = mesh.complex();
                cx.outgoing_halfedges(c)
                    .map(|h| mesh.face_area(h / 3) / 3.0)
                    .sum::<f64>()
            })
            .unwrap_or(0.0)
    };
    for v in 0..nf {
        let m = mass_of(mesh, v);
        if let Some(cost) = memoryless_cost(mesh, v, m) {
            heap.push(Entry {
                cost,
                vertex: v,
                version: 0,
            });
        }
    }
    let mut removed = 0;
    while mesh.n_vertices() > target {
        let Some(e) = heap.pop() else { break };
        if !mesh.is_vertex(e.vertex) || e.version != version[e.vertex as usize] {
            continue;
        }
        let report = match mesh.remove_vertex(e.vertex) {
            Ok(r) => r,
            Err(IemError::Removal(_)) => continue,
            Err(err) => return Err(err),
        };
        removed += 1;
        // Recompute the neighbours and their neighbours.
        let mut touched: Vec<u32> = report
            .prepared
            .neighbours
            .iter()
            .map(|n| n.vertex)
            .collect();
        let first = touched.clone();
        for &u in &first {
            if let Some(c) = mesh.complex_vertex(u) {
                for h in mesh.complex().outgoing_halfedges(c) {
                    touched.push(mesh.fine_vertex(mesh.complex().to_vertex(h)));
                }
            }
        }
        touched.sort_unstable();
        touched.dedup();
        for u in touched {
            if !mesh.is_vertex(u) {
                continue;
            }
            version[u as usize] += 1;
            let m = mass_of(mesh, u);
            if let Some(cost) = memoryless_cost(mesh, u, m) {
                heap.push(Entry {
                    cost,
                    vertex: u,
                    version: version[u as usize],
                });
            }
        }
    }
    Ok(removed)
}
