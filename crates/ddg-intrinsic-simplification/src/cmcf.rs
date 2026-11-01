//! Conformalized mean curvature flow (Kazhdan, Solomon, Ben-Chen 2012, Sec. 3.2 and 5)
//! accelerated by the intrinsic multigrid hierarchy (paper Sec. 8.7, Fig. 30; SPEC §9.3).
//!
//! > Michael Kazhdan, Jake Solomon, Mirela Ben-Chen. *Can Mean-Curvature Flow be Modified to be
//! > Non-singular?* Computer Graphics Forum 31(5), 2012.
//!
//! The semi-implicit step of the conformalized flow keeps the stiffness matrix of the input and
//! updates only the mass matrix: `(M_t − h L̃₀) X_{t+1} = M_t X_t`, with `L̃₀` the (negative
//! semidefinite) stiffness matrix of the input. With this crate's positive semidefinite cotan
//! Laplacian `L₀` (`CLAUDE.md`) the system is
//!
//! ```text
//!     (M_t + h L₀) X_{t+1} = M_t X_t            (symmetric positive definite)
//! ```
//!
//! `L₀` is the cotan Laplacian of the input positions (fixed), `M_t` the diagonal lumped
//! (barycentric) mass matrix of the current positions, **not** normalized, so `h` has the units
//! of an area (length²) of the input; a sphere of radius `r` shrinks as `r² = r₀² − 4t` in the
//! continuum limit. The three coordinates share the matrix.
//!
//! Two solvers, same state: [`CmcfFlow::step_direct`] (sparse Cholesky; the symbolic analysis is
//! computed once, the numeric factorization every step) and [`CmcfFlow::step_multigrid`]
//! (conjugate gradients preconditioned by one V-cycle of the Galerkin multigrid of
//! [`crate::multigrid`] on the **input's** intrinsic [`Hierarchy`], built once; the prolongations
//! `P_k` are reused and the Galerkin products `P_kᵀ A P_k` are rebuilt every step because `M_t`
//! changes). The previous positions are the initial guess. The reports separate assembly, setup
//! (Galerkin products and coarsest factorization, or the numeric factorization) and solve.

use std::time::Instant;

use ddg_core::operators::{cotan_laplacian, lumped_mass};
use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use ddg_core::{Geometry, SparseError, SparseMatrix, TriMesh};

use crate::hierarchy::Hierarchy;
use crate::multigrid::{MultigridError, MultigridOptions};

/// Errors of the flow.
#[derive(Debug, Clone, PartialEq)]
pub enum CmcfError {
    /// The positions or the step are not valid (wrong count, non-finite, a degenerate face, a
    /// non-positive step).
    BadInput {
        /// What is wrong.
        reason: String,
    },
    /// The hierarchy is not one of this mesh (its finest level has other vertices).
    HierarchyMismatch {
        /// What is wrong.
        reason: String,
    },
    /// The multigrid setup or solve failed.
    Multigrid(MultigridError),
    /// The sparse backend failed (the direct solve, the matrix assembly).
    Sparse(SparseError),
}

impl std::fmt::Display for CmcfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CmcfError::BadInput { reason } => write!(f, "cMCF: {reason}"),
            CmcfError::HierarchyMismatch { reason } => write!(f, "cMCF: hierarchy: {reason}"),
            CmcfError::Multigrid(e) => write!(f, "cMCF: {e}"),
            CmcfError::Sparse(e) => write!(f, "cMCF: sparse backend: {e}"),
        }
    }
}

impl std::error::Error for CmcfError {}

impl From<MultigridError> for CmcfError {
    fn from(e: MultigridError) -> Self {
        CmcfError::Multigrid(e)
    }
}

impl From<SparseError> for CmcfError {
    fn from(e: SparseError) -> Self {
        CmcfError::Sparse(e)
    }
}

/// Options of [`CmcfFlow::step_multigrid`].
#[derive(Debug, Clone, PartialEq)]
pub struct CmcfMultigridOptions {
    /// Smoothing sweeps and stopping tolerance (relative residual) of the PCG.
    pub multigrid: MultigridOptions,
    /// Cap on the PCG iterations per coordinate (default 100).
    pub max_iterations: usize,
}

impl Default for CmcfMultigridOptions {
    fn default() -> Self {
        CmcfMultigridOptions {
            multigrid: MultigridOptions {
                tolerance: 1e-10,
                ..MultigridOptions::default()
            },
            max_iterations: 100,
        }
    }
}

/// What one step did.
#[derive(Debug, Clone, PartialEq)]
pub struct CmcfStepReport {
    /// `"direct"` or `"multigrid"`.
    pub solver: &'static str,
    /// Seconds assembling `A = M_t + h L₀` and the right-hand sides.
    pub assemble_seconds: f64,
    /// Seconds in the setup: the numeric factorization, or the Galerkin hierarchy and the
    /// coarsest factorization.
    pub setup_seconds: f64,
    /// Seconds in the solve (substitutions, or the PCG iterations).
    pub solve_seconds: f64,
    /// PCG iterations per coordinate (empty for the direct solver).
    pub iterations: Vec<usize>,
    /// Whether every PCG reached its tolerance (`true` for the direct solver).
    pub converged: bool,
    /// Largest relative residual `‖b − A x‖ / ‖b‖` over the coordinates, computed after the
    /// solve for both solvers.
    pub relative_residual: f64,
}

impl CmcfStepReport {
    /// Total seconds of the step.
    pub fn total_seconds(&self) -> f64 {
        self.assemble_seconds + self.setup_seconds + self.solve_seconds
    }
}

/// The flow's state.
#[derive(Clone)]
pub struct CmcfFlow {
    mesh: TriMesh,
    positions: Vec<[f64; 3]>,
    l0: SparseMatrix,
    step: f64,
    symbolic: Symbolic,
    steps_done: usize,
}

fn masses(mesh: &TriMesh, pos: &[[f64; 3]]) -> Result<Vec<f64>, CmcfError> {
    let g = Geometry::from_positions(mesh, pos).map_err(|e| CmcfError::BadInput {
        reason: format!("the current positions: {e}"),
    })?;
    lumped_mass(mesh, &g).map_err(|e| CmcfError::BadInput {
        reason: e.to_string(),
    })
}

impl CmcfFlow {
    /// Starts a flow from `positions` with the time step `h` (an area; see the module docs).
    ///
    /// # Errors
    /// [`CmcfError::BadInput`] for a wrong count, non-finite positions, a degenerate face or
    /// `h ≤ 0`; the sparse analysis failing.
    pub fn new(mesh: &TriMesh, positions: &[[f64; 3]], step: f64) -> Result<CmcfFlow, CmcfError> {
        if !(step > 0.0 && step.is_finite()) {
            return Err(CmcfError::BadInput {
                reason: format!("the time step must be positive and finite, got {step}"),
            });
        }
        let g = Geometry::from_positions(mesh, positions).map_err(|e| CmcfError::BadInput {
            reason: format!("the input positions: {e}"),
        })?;
        let l0 = cotan_laplacian(mesh, &g).map_err(|e| CmcfError::BadInput {
            reason: e.to_string(),
        })?;
        let a = Self::system_matrix(&l0, step, &masses(mesh, positions)?)?;
        let symbolic = Symbolic::analyze(&a, Method::Llt)?;
        Ok(CmcfFlow {
            mesh: mesh.clone(),
            positions: positions.to_vec(),
            l0,
            step,
            symbolic,
            steps_done: 0,
        })
    }

    fn system_matrix(l0: &SparseMatrix, h: f64, m: &[f64]) -> Result<SparseMatrix, CmcfError> {
        let mut a = l0.clone();
        for v in a.values_mut() {
            *v *= h;
        }
        a.add_to_diagonal(m)?;
        Ok(a)
    }

    /// The current positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// The mesh.
    pub fn mesh(&self) -> &TriMesh {
        &self.mesh
    }

    /// The fixed stiffness matrix `L₀` (positive semidefinite) of the input.
    pub fn laplacian(&self) -> &SparseMatrix {
        &self.l0
    }

    /// The time step `h`.
    pub fn step_size(&self) -> f64 {
        self.step
    }

    /// Steps taken.
    pub fn steps_done(&self) -> usize {
        self.steps_done
    }

    /// The system of the next step: `A` and the three right-hand sides `M_t X_t`, column after
    /// column.
    fn assemble(&self) -> Result<(SparseMatrix, Vec<f64>), CmcfError> {
        let m = masses(&self.mesh, &self.positions)?;
        let a = Self::system_matrix(&self.l0, self.step, &m)?;
        let n = m.len();
        let mut b = vec![0.0; 3 * n];
        for (i, p) in self.positions.iter().enumerate() {
            for c in 0..3 {
                b[i + n * c] = m[i] * p[c];
            }
        }
        Ok((a, b))
    }

    fn store(&mut self, x: &[f64]) {
        let n = self.positions.len();
        for (i, p) in self.positions.iter_mut().enumerate() {
            *p = [x[i], x[i + n], x[i + 2 * n]];
        }
        self.steps_done += 1;
    }

    fn residual(a: &SparseMatrix, x: &[f64], b: &[f64]) -> Result<f64, CmcfError> {
        let n = a.nrows();
        let mut worst = 0.0f64;
        let mut ax = vec![0.0; n];
        for c in 0..3 {
            a.mul_vec(&x[c * n..(c + 1) * n], &mut ax)?;
            let (mut r, mut s) = (0.0, 0.0);
            for i in 0..n {
                r += (b[c * n + i] - ax[i]).powi(2);
                s += b[c * n + i].powi(2);
            }
            worst = worst.max((r / s.max(f64::MIN_POSITIVE)).sqrt());
        }
        Ok(worst)
    }

    /// One step with a direct sparse Cholesky solve.
    ///
    /// # Errors
    /// [`CmcfError`]: a degenerate current mesh, a failed factorization.
    pub fn step_direct(&mut self) -> Result<CmcfStepReport, CmcfError> {
        let t0 = Instant::now();
        let (a, b) = self.assemble()?;
        let assemble_seconds = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let f: Factorization = self.symbolic.factor(&a, Parallelism::Sequential)?;
        let setup_seconds = t1.elapsed().as_secs_f64();
        let t2 = Instant::now();
        let mut x = b.clone();
        f.solve_many_in_place(&mut x, 3, Parallelism::Sequential)?;
        let solve_seconds = t2.elapsed().as_secs_f64();
        let relative_residual = Self::residual(&a, &x, &b)?;
        self.store(&x);
        Ok(CmcfStepReport {
            solver: "direct",
            assemble_seconds,
            setup_seconds,
            solve_seconds,
            iterations: Vec::new(),
            converged: true,
            relative_residual,
        })
    }

    /// One step by PCG with a V-cycle of the Galerkin multigrid on `hierarchy` (the intrinsic
    /// hierarchy of the input, built once; its `P_k` are reused, the Galerkin products are
    /// rebuilt for this step's `M_t`).
    ///
    /// # Errors
    /// [`CmcfError::HierarchyMismatch`] if the finest level of `hierarchy` is not this mesh's
    /// vertex set; multigrid setup failures; a degenerate current mesh. Not reaching the
    /// tolerance is reported (`converged = false`), not an error.
    pub fn step_multigrid(
        &mut self,
        hierarchy: &Hierarchy,
        options: &CmcfMultigridOptions,
    ) -> Result<CmcfStepReport, CmcfError> {
        let n = self.positions.len();
        let level0 = &hierarchy.levels()[0];
        if level0.n_vertices() != n
            || level0
                .vertices
                .iter()
                .enumerate()
                .any(|(i, &v)| v as usize != i)
        {
            return Err(CmcfError::HierarchyMismatch {
                reason: format!(
                    "the finest level has {} vertices (ids 0..{n} expected), the mesh {n}",
                    level0.n_vertices(),
                ),
            });
        }
        let t0 = Instant::now();
        let (a, b) = self.assemble()?;
        let assemble_seconds = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let mg = hierarchy.multigrid(&a, &[], options.multigrid.clone())?;
        let setup_seconds = t1.elapsed().as_secs_f64();
        let t2 = Instant::now();
        let mut x = vec![0.0; 3 * n];
        let mut iterations = Vec::with_capacity(3);
        let mut converged = true;
        for c in 0..3 {
            let guess: Vec<f64> = self.positions.iter().map(|p| p[c]).collect();
            let (xc, rep) = mg.pcg(
                &b[c * n..(c + 1) * n],
                &[],
                Some(&guess),
                options.max_iterations,
            )?;
            iterations.push(rep.cycles);
            converged &= rep.converged;
            x[c * n..(c + 1) * n].copy_from_slice(&xc);
        }
        let solve_seconds = t2.elapsed().as_secs_f64();
        let relative_residual = Self::residual(&a, &x, &b)?;
        self.store(&x);
        Ok(CmcfStepReport {
            solver: "multigrid",
            assemble_seconds,
            setup_seconds,
            solve_seconds,
            iterations,
            converged,
            relative_residual,
        })
    }
}
