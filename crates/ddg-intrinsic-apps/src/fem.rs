//! Finite elements on an intrinsic triangulation (paper Sec. 5.2, Fig. 16, 17).
//!
//! Piecewise linear elements on any triangulation of the surface are assembled from edge
//! lengths alone, "exactly as with an ordinary triangle mesh" (paper Sec. 5.2): the stiffness
//! matrix is the cotan Laplacian (MacNeal 1949), the mass matrix the standard Galerkin mass
//! matrix or its lumped (diagonal) form. An intrinsic triangulation decouples the elements from
//! the mesh that describes the shape, so Delaunay flipping, refinement and the optimal
//! Delaunay triangulation of `ddg-intrinsic` give better elements on the same surface.
//!
//! What is here:
//!
//! * [`stiffness_matrix`], [`galerkin_mass_matrix`], [`lumped_mass`] on an
//!   [`IntrinsicTriangulation`] (signs and conventions of the workspace: `L` positive
//!   semidefinite with zero row sums).
//! * [`regular_subdivision`]: the regular 4-1 subdivision of the input mesh, the extrinsic
//!   refinement of Fig. 17 (every triangle split into four similar ones; the surface and the
//!   angles do not change).
//! * [`condition_number`]: the largest and the smallest non-zero eigenvalue of `L x = λ M x`
//!   for `M` the identity (the stiffness matrix alone), the lumped or the Galerkin mass matrix,
//!   by the Lanczos method with full reorthogonalization: in the `M` inner product on
//!   `M⁻¹ L` for the largest, and on `L⁺ M` (shift-and-invert at zero, the constants of every
//!   component projected out) for the smallest.
//! * [`PinnedLaplacian`]: solves `L x = b` exactly for a right-hand side with zero sum on every
//!   component, by a Cholesky factorization of `L` with one vertex per component removed. It
//!   is shared with [`heat`](crate::heat).

use std::fmt;

use ddg_core::generators::SplitMix64;
use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use ddg_core::{Geometry, GeometryError, MeshError, SparseError, SparseMatrix, TriMesh};
use ddg_intrinsic::IntrinsicTriangulation;

/// Why a function of this module failed.
#[derive(Debug, Clone, PartialEq)]
pub enum FemError {
    /// A sparse factorization or solve failed.
    Sparse(SparseError),
    /// The subdivided mesh could not be built.
    Mesh(MeshError),
    /// The lengths of the subdivided mesh were refused.
    Geometry(GeometryError),
    /// The matrices do not fit together or an option is out of range.
    Invalid(String),
}

impl fmt::Display for FemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FemError::Sparse(e) => write!(f, "{e}"),
            FemError::Mesh(e) => write!(f, "{e}"),
            FemError::Geometry(e) => write!(f, "{e}"),
            FemError::Invalid(m) => write!(f, "{m}"),
        }
    }
}

impl std::error::Error for FemError {}

impl From<SparseError> for FemError {
    fn from(e: SparseError) -> Self {
        FemError::Sparse(e)
    }
}
impl From<MeshError> for FemError {
    fn from(e: MeshError) -> Self {
        FemError::Mesh(e)
    }
}
impl From<GeometryError> for FemError {
    fn from(e: GeometryError) -> Self {
        FemError::Geometry(e)
    }
}

/// The stiffness matrix of piecewise linear elements on the triangulation: the cotan
/// Laplacian, `L_ij = −½ Σ cot θ` over the corners opposite the edges between `i` and `j`,
/// `L_ii = −Σ_{j≠i} L_ij` (positive semidefinite). It is
/// [`IntrinsicTriangulation::cotan_laplacian`].
///
/// [formula: cotan-fem-stiffness]
pub fn stiffness_matrix(tri: &IntrinsicTriangulation) -> SparseMatrix {
    tri.cotan_laplacian()
}

/// The lumped mass: a third of the area of every triangle at each of its corners
/// ([`IntrinsicTriangulation::lumped_mass`]). Not normalized: the masses sum to the area.
pub fn lumped_mass(tri: &IntrinsicTriangulation) -> Vec<f64> {
    tri.lumped_mass()
}

/// The Galerkin (consistent) mass matrix `M_ij = ∫ φ_i φ_j` of piecewise linear elements:
/// per triangle of area `A`, `A/6` for a corner with itself and `A/12` for two different
/// corners, summed over the corners at the vertices. A triangle that uses a vertex twice adds
/// both corners to it, which is the integral of its hat function there. Row `i` sums to the
/// lumped mass of `i`; the matrix is symmetric positive definite.
///
/// [formula: app-fem-galerkin-mass]
pub fn galerkin_mass_matrix(tri: &IntrinsicTriangulation) -> SparseMatrix {
    let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(9 * tri.n_faces());
    for (f, face) in tri.faces().iter().enumerate() {
        let area = tri.face_area(f);
        for a in 0..3 {
            for b in 0..3 {
                let w = if a == b { area / 6.0 } else { area / 12.0 };
                triplets.push((face[a], face[b], w));
            }
        }
    }
    let n = tri.n_vertices();
    SparseMatrix::from_triplets(n, n, &triplets)
        .expect("vertex indices of a valid triangulation are in range")
}

/// The regular 4-1 subdivision, applied `levels` times: a vertex at the midpoint of every
/// edge, and every triangle replaced by the three triangles at its corners and the one in the
/// middle. Computed from edge lengths alone: the halves of an edge have half its length and an
/// edge between two midpoints has half the length of the side it is parallel to. The surface,
/// the angles and the angle sums do not change; vertices `0..n` are those of the input, vertex
/// `n + e` is the midpoint of edge `e` (at each level).
///
/// [formula: app-fem-midpoint-subdivision]
///
/// # Errors
/// [`FemError::Geometry`] if the geometry is not that of the mesh.
pub fn regular_subdivision(
    mesh: &TriMesh,
    geometry: &Geometry,
    levels: usize,
) -> Result<(TriMesh, Geometry), FemError> {
    geometry.check_mesh(mesh)?;
    let mut mesh = mesh.clone();
    let mut geometry = geometry.clone();
    for _ in 0..levels {
        let n = mesh.n_vertices() as u32;
        let lengths = geometry.edge_lengths();
        let mut faces = Vec::with_capacity(4 * mesh.n_faces());
        // (vertex pair, length) of every new edge.
        let mut new_lengths: Vec<((u32, u32), f64)> = Vec::with_capacity(9 * mesh.n_faces());
        for (f, &[a, b, c]) in mesh.faces().iter().enumerate() {
            let e = [0, 1, 2].map(|k| mesh.edge(3 * f + k));
            let [mab, mbc, mca] = e.map(|e| n + e as u32);
            let [lab, lbc, lca] = e.map(|e| lengths[e]);
            faces.push([a, mab, mca]);
            faces.push([mab, b, mbc]);
            faces.push([mca, mbc, c]);
            faces.push([mab, mbc, mca]);
            for (p, q, l) in [
                (a, mab, lab / 2.0),
                (mab, b, lab / 2.0),
                (b, mbc, lbc / 2.0),
                (mbc, c, lbc / 2.0),
                (c, mca, lca / 2.0),
                (mca, a, lca / 2.0),
                (mab, mbc, lca / 2.0),
                (mbc, mca, lab / 2.0),
                (mca, mab, lbc / 2.0),
            ] {
                new_lengths.push(((p.min(q), p.max(q)), l));
            }
        }
        let fine = TriMesh::new(n as usize + mesh.n_edges(), faces)?;
        let mut l = vec![0.0; fine.n_edges()];
        for ((p, q), length) in new_lengths {
            let e = fine
                .find_edge(p as usize, q as usize)
                .expect("every listed pair is an edge of the subdivided mesh");
            l[e] = length;
        }
        geometry = Geometry::from_edge_lengths(&fine, l)?;
        mesh = fine;
    }
    Ok((mesh, geometry))
}

/// The connected components of the graph of a square matrix's pattern: a component id per
/// row, numbered in the order of their smallest row.
fn pattern_components(a: &SparseMatrix) -> (Vec<u32>, usize) {
    let n = a.nrows();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for (i, j, _) in a.triplets() {
        let (ri, rj) = (find(&mut parent, i as usize), find(&mut parent, j as usize));
        if ri != rj {
            parent[ri.max(rj)] = ri.min(rj);
        }
    }
    let mut id = vec![u32::MAX; n];
    let mut component = vec![0u32; n];
    let mut count = 0;
    for v in 0..n {
        let r = find(&mut parent, v);
        if id[r] == u32::MAX {
            id[r] = count as u32;
            count += 1;
        }
        component[v] = id[r];
    }
    (component, count)
}

/// `L x = b` for a symmetric positive semidefinite `L` whose kernel is the constants on every
/// component of its pattern (the cotan Laplacian of a triangulation with non-degenerate
/// triangles), solved with `x = 0` at the first vertex of every component by a Cholesky
/// factorization of `L` without those rows and columns. For a right-hand side that sums to zero
/// on every component this is an exact solution of the whole system (the removed rows hold
/// because all rows sum to zero); otherwise the removed rows are off by the sums, which
/// [`PinnedLaplacian::inconsistency`] measures.
pub struct PinnedLaplacian {
    n: usize,
    component: Vec<u32>,
    pinned: Vec<u32>,
    /// Row of the reduced system for every vertex (`u32::MAX` for a pinned one).
    reduced: Vec<u32>,
    factor: Option<Factorization>,
}

impl PinnedLaplacian {
    /// Factors `L` with one vertex per component removed.
    ///
    /// # Errors
    /// [`FemError::Invalid`] for a matrix that is not square; [`FemError::Sparse`] when the
    /// reduced matrix is not positive definite (a degenerate triangle, or a matrix that is not
    /// a Laplacian).
    pub fn new(l: &SparseMatrix) -> Result<PinnedLaplacian, FemError> {
        if l.nrows() != l.ncols() {
            return Err(FemError::Invalid(format!(
                "the Laplacian must be square, not {} x {}",
                l.nrows(),
                l.ncols()
            )));
        }
        let n = l.nrows();
        let (component, count) = pattern_components(l);
        let mut pinned = vec![u32::MAX; count];
        for v in 0..n {
            let c = component[v] as usize;
            if pinned[c] == u32::MAX {
                pinned[c] = v as u32;
            }
        }
        let mut reduced = vec![u32::MAX; n];
        let mut keep: Vec<u32> = Vec::with_capacity(n - count);
        for v in 0..n {
            if pinned[component[v] as usize] != v as u32 {
                reduced[v] = keep.len() as u32;
                keep.push(v as u32);
            }
        }
        let factor = if keep.is_empty() {
            None
        } else {
            let a = l.submatrix(&keep, &keep)?;
            let symbolic = Symbolic::analyze(&a, Method::Llt)?;
            Some(symbolic.factor(&a, Parallelism::Sequential)?)
        };
        Ok(PinnedLaplacian {
            n,
            component,
            pinned,
            reduced,
            factor,
        })
    }

    /// The number of components (and of pinned vertices).
    pub fn n_components(&self) -> usize {
        self.pinned.len()
    }

    /// The component of every vertex.
    pub fn components(&self) -> &[u32] {
        &self.component
    }

    /// The largest absolute sum of `b` over a component: zero (to rounding) when `L x = b` is
    /// solvable.
    pub fn inconsistency(&self, b: &[f64]) -> f64 {
        let mut sums = vec![0.0; self.pinned.len()];
        for (v, &x) in b.iter().enumerate() {
            sums[self.component[v] as usize] += x;
        }
        sums.iter().fold(0.0f64, |m, s| m.max(s.abs()))
    }

    /// Overwrites `b` with the solution that is zero at the pinned vertices.
    ///
    /// # Errors
    /// [`FemError::Invalid`] for a vector of the wrong length, [`FemError::Sparse`] if the
    /// solve fails.
    pub fn solve_in_place(&self, b: &mut [f64]) -> Result<(), FemError> {
        if b.len() != self.n {
            return Err(FemError::Invalid(format!(
                "the right-hand side has {} entries, the Laplacian {} rows",
                b.len(),
                self.n
            )));
        }
        let mut rhs: Vec<f64> = (0..self.n)
            .filter(|&v| self.reduced[v] != u32::MAX)
            .map(|v| b[v])
            .collect();
        if let Some(factor) = &self.factor {
            factor.solve_in_place(&mut rhs, Parallelism::Sequential)?;
        }
        for v in 0..self.n {
            let r = self.reduced[v];
            b[v] = if r == u32::MAX { 0.0 } else { rhs[r as usize] };
        }
        Ok(())
    }
}

/// The mass matrix of the eigenvalue problem `L x = λ M x` in [`condition_number`].
#[derive(Debug, Clone, PartialEq)]
pub enum MassMatrix {
    /// `M = I`: the eigenvalues of the stiffness matrix itself.
    Identity,
    /// A diagonal mass matrix (e.g. [`lumped_mass`]); every entry positive.
    Lumped(Vec<f64>),
    /// A symmetric positive definite mass matrix (e.g. [`galerkin_mass_matrix`]) with the
    /// pattern of the stiffness matrix or a part of it.
    Galerkin(SparseMatrix),
}

/// Options of [`condition_number`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EigenOptions {
    /// The largest number of Lanczos steps for each of the two eigenvalues. Default 300.
    pub max_steps: usize,
    /// Converged when the residual bound `β_k |s_k|` of the Ritz value is below this fraction
    /// of it. Default `1e-10`.
    pub tolerance: f64,
    /// Seed of the random start vector. Default 1.
    pub seed: u64,
}

impl Default for EigenOptions {
    fn default() -> Self {
        EigenOptions {
            max_steps: 300,
            tolerance: 1e-10,
            seed: 1,
        }
    }
}

/// One extreme eigenvalue found by the Lanczos method.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Eigenvalue {
    /// The eigenvalue (the Ritz value).
    pub value: f64,
    /// Lanczos steps taken.
    pub steps: usize,
    /// The residual bound relative to the value: `|L x − λ M x|_{M⁻¹} / |λ|` for the Ritz
    /// vector `x` (`M`-normalized).
    pub relative_residual: f64,
    /// Whether the residual met the tolerance (or the Krylov space became invariant).
    pub converged: bool,
}

/// The extreme eigenvalues of `L x = λ M x` and their ratio. See [`condition_number`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConditionReport {
    /// The largest eigenvalue.
    pub lambda_max: Eigenvalue,
    /// The smallest eigenvalue on the complement of the kernel (the constants of every
    /// component).
    pub lambda_min: Eigenvalue,
    /// `lambda_max / lambda_min`.
    pub condition_number: f64,
    /// The number of components of the pattern (the dimension of the kernel).
    pub n_components: usize,
}

/// The largest and the smallest non-zero eigenvalue of `L x = λ M x` and the condition number
/// `λ_max / λ_min` (paper Fig. 17), for the stiffness matrix `stiffness` (positive
/// semidefinite, kernel the constants per component) and the mass matrix `mass`.
///
/// Lanczos with full reorthogonalization in the `M` inner product: on `M⁻¹ L` for `λ_max`
/// (a solve with `M`: none, a division, or a Cholesky factorization), and on `P L⁺ M` for
/// `λ_min`, where `L⁺` is applied with [`PinnedLaplacian`] and `P` projects out the constants
/// of every component `M`-orthogonally, so that its largest eigenvalue is `1 / λ_min`. Both start
/// from a random vector (`options.seed`). Not converging within `max_steps` is reported in
/// [`Eigenvalue::converged`], not an error.
///
/// # Errors
/// [`FemError::Invalid`] for matrices of different sizes or a non-positive lumped mass,
/// [`FemError::Sparse`] for failed factorizations.
pub fn condition_number(
    stiffness: &SparseMatrix,
    mass: &MassMatrix,
    options: &EigenOptions,
) -> Result<ConditionReport, FemError> {
    let n = stiffness.nrows();
    if stiffness.ncols() != n || n == 0 {
        return Err(FemError::Invalid(format!(
            "the stiffness matrix must be square and not empty, not {} x {}",
            stiffness.nrows(),
            stiffness.ncols()
        )));
    }
    if options.max_steps == 0 || !(options.tolerance > 0.0) {
        return Err(FemError::Invalid(
            "EigenOptions: max_steps must be positive and tolerance a positive number".into(),
        ));
    }
    let mass_factor = match mass {
        MassMatrix::Identity => None,
        MassMatrix::Lumped(m) => {
            if m.len() != n {
                return Err(FemError::Invalid(format!(
                    "{} lumped masses for {n} vertices",
                    m.len()
                )));
            }
            if let Some((i, x)) = m.iter().enumerate().find(|(_, x)| !(**x > 0.0)) {
                return Err(FemError::Invalid(format!(
                    "the lumped mass of vertex {i} is {x}, not positive"
                )));
            }
            None
        }
        MassMatrix::Galerkin(m) => {
            if m.nrows() != n || m.ncols() != n {
                return Err(FemError::Invalid(format!(
                    "the mass matrix is {} x {}, the stiffness matrix {n} x {n}",
                    m.nrows(),
                    m.ncols()
                )));
            }
            let symbolic = Symbolic::analyze(m, Method::Llt)?;
            Some(symbolic.factor(m, Parallelism::Sequential)?)
        }
    };
    let mul_mass = |x: &[f64], y: &mut [f64]| -> Result<(), FemError> {
        match mass {
            MassMatrix::Identity => y.copy_from_slice(x),
            MassMatrix::Lumped(m) => {
                for i in 0..n {
                    y[i] = m[i] * x[i];
                }
            }
            MassMatrix::Galerkin(m) => m.mul_vec(x, y)?,
        }
        Ok(())
    };
    let mut scratch = vec![0.0; n];
    let mut m_dot = |x: &[f64], y: &[f64]| -> Result<f64, FemError> {
        mul_mass(y, &mut scratch)?;
        Ok(x.iter().zip(&scratch).map(|(a, b)| a * b).sum())
    };

    // Largest: M⁻¹ L.
    let mut tmp = vec![0.0; n];
    let lambda_max = lanczos_largest(
        n,
        |x, y| {
            stiffness.mul_vec(x, y)?;
            match mass {
                MassMatrix::Identity => {}
                MassMatrix::Lumped(m) => {
                    for i in 0..n {
                        y[i] /= m[i];
                    }
                }
                MassMatrix::Galerkin(_) => {
                    let f = mass_factor.as_ref().expect("factored above");
                    f.solve_in_place(y, Parallelism::Sequential)?;
                }
            }
            Ok(())
        },
        &mut m_dot,
        |_| Ok(()),
        options,
    )?;

    // Smallest non-zero: P L⁺ M, constants projected out in the M inner product.
    let pinned = PinnedLaplacian::new(stiffness)?;
    let comp = pinned.components().to_vec();
    let n_comp = pinned.n_components();
    // M 1_c for every component, and 1_cᵀ M 1_c.
    let mut ones_mass = vec![0.0; n];
    mul_mass(&vec![1.0; n], &mut ones_mass)?;
    let mut comp_mass = vec![0.0; n_comp];
    for v in 0..n {
        comp_mass[comp[v] as usize] += ones_mass[v];
    }
    let project = |x: &mut [f64]| -> Result<(), FemError> {
        let mut dots = vec![0.0; n_comp];
        for v in 0..n {
            dots[comp[v] as usize] += ones_mass[v] * x[v];
        }
        for v in 0..n {
            let c = comp[v] as usize;
            x[v] -= dots[c] / comp_mass[c];
        }
        Ok(())
    };
    let inverse = lanczos_largest(
        n,
        |x, y| {
            mul_mass(x, &mut tmp)?;
            y.copy_from_slice(&tmp);
            pinned.solve_in_place(y)?;
            project(y)
        },
        &mut m_dot,
        project,
        options,
    )?;
    let lambda_min = Eigenvalue {
        value: 1.0 / inverse.value,
        ..inverse
    };
    Ok(ConditionReport {
        lambda_max,
        lambda_min,
        condition_number: lambda_max.value / lambda_min.value,
        n_components: n_comp,
    })
}

/// The largest eigenvalue of an operator that is self-adjoint in the inner product `dot`, by
/// Lanczos with full reorthogonalization from a random start (projected by `project`).
fn lanczos_largest(
    n: usize,
    mut apply: impl FnMut(&[f64], &mut [f64]) -> Result<(), FemError>,
    dot: &mut impl FnMut(&[f64], &[f64]) -> Result<f64, FemError>,
    project: impl Fn(&mut [f64]) -> Result<(), FemError>,
    options: &EigenOptions,
) -> Result<Eigenvalue, FemError> {
    let mut rng = SplitMix64::new(options.seed);
    let mut q: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    project(&mut q)?;
    let norm = dot(&q, &q)?.sqrt();
    if !(norm > 0.0) {
        return Err(FemError::Invalid(
            "the start vector vanishes after projection (no eigenvalue outside the kernel)"
                .into(),
        ));
    }
    q.iter_mut().for_each(|x| *x /= norm);
    let mut basis: Vec<Vec<f64>> = vec![q];
    let mut alpha: Vec<f64> = Vec::new();
    let mut beta: Vec<f64> = Vec::new();
    let mut w = vec![0.0; n];
    let mut best = Eigenvalue {
        value: 0.0,
        steps: 0,
        relative_residual: f64::INFINITY,
        converged: false,
    };
    for k in 0..options.max_steps.min(n) {
        apply(&basis[k], &mut w)?;
        let a = dot(&w, &basis[k])?;
        alpha.push(a);
        // Full reorthogonalization, twice.
        for _ in 0..2 {
            for v in &basis {
                let c = dot(&w, v)?;
                w.iter_mut().zip(v).for_each(|(x, y)| *x -= c * y);
            }
        }
        let b = dot(&w, &w)?.max(0.0).sqrt();
        let steps = k + 1;
        let check = steps % 5 == 0 || steps == options.max_steps.min(n) || b == 0.0;
        if check {
            let (theta, last) = tridiagonal_largest(&alpha, &beta);
            let residual = (b * last).abs() / theta.abs();
            best = Eigenvalue {
                value: theta,
                steps,
                relative_residual: residual,
                converged: residual <= options.tolerance,
            };
            if best.converged {
                break;
            }
        }
        if b <= 1e-14 * alpha.iter().fold(0.0f64, |m, x| m.max(x.abs())) {
            // An invariant subspace: the Ritz values are eigenvalues.
            let (theta, _) = tridiagonal_largest(&alpha, &beta);
            best = Eigenvalue {
                value: theta,
                steps,
                relative_residual: 0.0,
                converged: true,
            };
            break;
        }
        beta.push(b);
        basis.push(w.iter().map(|x| x / b).collect());
    }
    Ok(best)
}

/// The largest eigenvalue of the symmetric tridiagonal matrix with diagonal `alpha` and
/// off-diagonal `beta` (one shorter), and the last component of its unit eigenvector (cyclic
/// Jacobi on the dense matrix; the size is that of the Krylov space).
fn tridiagonal_largest(alpha: &[f64], beta: &[f64]) -> (f64, f64) {
    let k = alpha.len();
    let mut a = vec![0.0; k * k];
    let mut v = vec![0.0; k * k];
    for i in 0..k {
        a[i * k + i] = alpha[i];
        v[i * k + i] = 1.0;
        if i + 1 < k {
            a[i * k + i + 1] = beta[i];
            a[(i + 1) * k + i] = beta[i];
        }
    }
    for _sweep in 0..100 {
        let off: f64 = (0..k)
            .flat_map(|i| (0..k).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i * k + j] * a[i * k + j])
            .sum();
        let diag: f64 = (0..k).map(|i| a[i * k + i] * a[i * k + i]).sum();
        if off <= 1e-30 * diag {
            break;
        }
        for p in 0..k {
            for r in p + 1..k {
                let apr = a[p * k + r];
                if apr == 0.0 {
                    continue;
                }
                let (app, arr) = (a[p * k + p], a[r * k + r]);
                let theta = (arr - app) / (2.0 * apr);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for i in 0..k {
                    let (aip, air) = (a[i * k + p], a[i * k + r]);
                    a[i * k + p] = c * aip - s * air;
                    a[i * k + r] = s * aip + c * air;
                }
                for i in 0..k {
                    let (api, ari) = (a[p * k + i], a[r * k + i]);
                    a[p * k + i] = c * api - s * ari;
                    a[r * k + i] = s * api + c * ari;
                }
                for i in 0..k {
                    let (vip, vir) = (v[i * k + p], v[i * k + r]);
                    v[i * k + p] = c * vip - s * vir;
                    v[i * k + r] = s * vip + c * vir;
                }
            }
        }
    }
    let mut best = 0;
    for i in 1..k {
        if a[i * k + i] > a[best * k + best] {
            best = i;
        }
    }
    (a[best * k + best], v[(k - 1) * k + best])
}
