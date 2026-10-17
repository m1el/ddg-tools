//! Tangent vector field processing on intrinsic triangulations (paper Sec. 5.5).
//!
//! The signposts of an intrinsic triangulation are tangent frames at its vertices, and the
//! vertex update keeps them consistent through flips, insertions and moves. That is exactly
//! the data of a discrete connection: the angle by which a vector turns when it is carried
//! along an edge ([`transport_angle`], Sharp, Soliman, Crane 2019, "The Vector Heat Method"
//! (VHM), Eq. 10). The frame of an input vertex is that of the input mesh, so a field computed
//! on any triangulation is read off at the input vertices as it is.
//!
//! * [`connection_matrices`]: the connection Laplacian and its mass matrix, by finite
//!   differences (VHM Sec. 5.3: cotan weights times transport, lumped mass) or by finite
//!   elements (Knöppel, Crane, Pinkall, Schröder 2013, "Globally Optimal Direction Fields"
//!   (GODF), Sec. 6.1.1: piecewise linear basis sections on triangles of constant curvature).
//! * [`harmonic_field`]: vectors prescribed on the boundary, `L∇ u = 0` inside (paper,
//!   "flip-free vector fields"); [`check_maximum_principle`] tests the vector maximum
//!   principle of VHM Sec. 5.4 (every vector in the convex cone of its transported
//!   neighbours), which holds on an intrinsic Delaunay triangulation.
//! * [`smoothest_field`]: the smoothest unit-norm `n`-direction field, the eigenvector of the
//!   smallest eigenvalue of `A u = λ M u` (GODF Sec. 6.1), by inverse iteration with a Cholesky
//!   factorization of the real form of `A + εM`; [`singularities`] its singular faces (GODF
//!   Eq. 19).
//! * [`TangentField::input_vectors_3d`]: the field at the input vertices as 3D vectors.
//!
//! # Conventions (SPEC N40 to N44)
//!
//! An `n`-direction field is stored by its `n`-th power: `u_v ∈ C`, with angle measured in
//! the tangent frame of `v`. The frame of an interior vertex is its signpost frame (angles
//! normalized to `2π`); at a boundary vertex the angles are **not** normalized (`s = 1`, GODF
//! Eq. 11), so that the connection of a flat domain is trivial: angle 0 is the boundary
//! halfedge leaving the vertex, `Θ` the boundary edge arriving at it. Complex matrices are
//! stored in real form (`vector/complex.rs`): entry `a + ib` becomes `[[a, −b], [b, a]]`.

use std::f64::consts::{PI, TAU};
use std::fmt;

use ddg_core::{SparseError, SparseMatrix};
use ddg_intrinsic::IntrinsicTriangulation;

mod complex;
mod connection;
mod eigen;
mod knoeppel;
pub(crate) mod linsolve;

pub use complex::C64;
pub use connection::{face_holonomy, transport_angle};
pub use linsolve::Factorized;

use complex::{from_real, to_real};
use connection::{fd_laplacian, fe_matrices, frame_factor, lumped_mass, tangent_scale, transport};

/// Errors of this module.
#[derive(Debug, Clone, PartialEq)]
pub enum VectorError {
    /// An option is out of range.
    InvalidOptions {
        /// What is wrong.
        reason: String,
    },
    /// Input data that does not fit the triangulation.
    InvalidInput {
        /// What is wrong.
        reason: String,
    },
    /// A sparse factorization or solve failed.
    Sparse(SparseError),
    /// The energy matrix is not positive definite, so the smallest eigenvector cannot be found
    /// by inverse iteration with a Cholesky factorization (a finite-difference connection
    /// Laplacian with negative cotan weights).
    NotPositiveDefinite {
        /// Which matrix.
        what: String,
    },
}

impl fmt::Display for VectorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VectorError::InvalidOptions { reason } => write!(f, "invalid options: {reason}"),
            VectorError::InvalidInput { reason } => write!(f, "invalid input: {reason}"),
            VectorError::Sparse(e) => write!(f, "sparse solver: {e}"),
            VectorError::NotPositiveDefinite { what } => {
                write!(f, "{what} is not positive definite")
            }
        }
    }
}

impl std::error::Error for VectorError {}

impl From<SparseError> for VectorError {
    fn from(e: SparseError) -> Self {
        VectorError::Sparse(e)
    }
}

/// Which discretization of the connection Laplacian.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConnectionLaplacian {
    /// Cotan weights times the transport of the edge, lumped vertex areas as the mass
    /// (VHM Sec. 5.3). Positive semidefinite when every cotan weight is non-negative (an
    /// intrinsic Delaunay triangulation), possibly indefinite otherwise.
    FiniteDifference,
    /// Piecewise linear basis sections on triangles of constant curvature, Galerkin mass
    /// (GODF Sec. 6.1.1). Always positive semidefinite (a Gram matrix).
    FiniteElement,
}

/// The connection Laplacian of a triangulation and its mass matrix, in real form
/// (`2 n_vertices` rows; see the [module documentation](self)).
#[derive(Debug, Clone, PartialEq)]
pub struct ConnectionMatrices {
    /// The symmetry: 1 for vectors, `n` for `n`-direction fields.
    pub n_sym: u32,
    /// The discretization.
    pub kind: ConnectionLaplacian,
    /// The stiffness (energy) matrix `A`: `u* A u` is the Dirichlet energy.
    pub stiffness: SparseMatrix,
    /// The mass matrix `M`: `u* M u` is the squared `L²` norm.
    pub mass: SparseMatrix,
}

impl ConnectionMatrices {
    /// `u* A u` for a field given by its `n`-th powers.
    pub fn energy(&self, values: &[C64]) -> f64 {
        let x = to_real(values);
        linsolve::dot(&x, &linsolve::mul(&self.stiffness, &x))
    }

    /// `u* M u`.
    pub fn norm_squared(&self, values: &[C64]) -> f64 {
        let x = to_real(values);
        linsolve::dot(&x, &linsolve::mul(&self.mass, &x))
    }
}

/// The connection Laplacian and mass matrix of `tri` for `n_sym`-vectors.
///
/// # Errors
/// [`VectorError::InvalidOptions`] for `n_sym = 0`; [`VectorError::Sparse`] if the matrix is
/// too large.
pub fn connection_matrices(
    tri: &IntrinsicTriangulation,
    n_sym: u32,
    kind: ConnectionLaplacian,
) -> Result<ConnectionMatrices, VectorError> {
    if n_sym == 0 {
        return Err(VectorError::InvalidOptions {
            reason: "the symmetry n must be at least 1".to_string(),
        });
    }
    let (stiffness, mass) = match kind {
        ConnectionLaplacian::FiniteDifference => {
            connection::build_pair(&fd_laplacian(tri, n_sym), &lumped_mass(tri))?
        }
        ConnectionLaplacian::FiniteElement => {
            let (m, a) = fe_matrices(tri, n_sym);
            connection::build_pair(&a, &m)?
        }
    };
    Ok(ConnectionMatrices {
        n_sym,
        kind,
        stiffness,
        mass,
    })
}

/// A field of `n`-vectors at the vertices of a triangulation.
#[derive(Debug, Clone, PartialEq)]
pub struct TangentField {
    /// The symmetry `n`.
    pub n_sym: u32,
    /// Per vertex of the triangulation, the `n`-th power of the vector in the tangent frame of
    /// the vertex (see the [module documentation](self)). The first
    /// `n_input_vertices` entries are the input vertices, in the frames of the input mesh.
    pub values: Vec<C64>,
}

impl TangentField {
    /// One of the `n` vectors of the field at every vertex: angle `arg(u)/n`, length `|u|`,
    /// in the tangent frame of the vertex.
    pub fn representative(&self) -> Vec<C64> {
        let n = self.n_sym as f64;
        self.values
            .iter()
            .map(|u| C64::from_angle(u.arg() / n).scale(u.abs()))
            .collect()
    }

    /// The field at the input vertices as 3D vectors: one of its `n` vectors, of length `|u|`,
    /// in the plane of the input face that the direction points into. `positions` are the
    /// positions of the input mesh's vertices.
    ///
    /// # Errors
    /// [`VectorError::InvalidInput`] if `positions` does not have one entry per input vertex
    /// or the field does not have one value per vertex of `tri`.
    pub fn input_vectors_3d(
        &self,
        tri: &IntrinsicTriangulation,
        positions: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, VectorError> {
        let input = tri.input();
        let mesh = input.mesh();
        if positions.len() != mesh.n_vertices() {
            return Err(VectorError::InvalidInput {
                reason: format!(
                    "{} positions for an input mesh with {} vertices",
                    positions.len(),
                    mesh.n_vertices()
                ),
            });
        }
        if self.values.len() != tri.n_vertices() {
            return Err(VectorError::InvalidInput {
                reason: format!(
                    "a field with {} values on a triangulation with {} vertices",
                    self.values.len(),
                    tri.n_vertices()
                ),
            });
        }
        let reps = self.representative();
        let mut out = Vec::with_capacity(mesh.n_vertices());
        for (v, &u) in reps.iter().enumerate().take(mesh.n_vertices()) {
            let boundary = mesh.is_boundary_vertex(v);
            let theta = input.vertex_angle_sum(v);
            let scale = ddg_intrinsic::formulas::angle_scale(theta, boundary);
            // The direction as a normalized angle (the signpost convention).
            let mut phi = u.arg() / frame_factor_input(theta, boundary);
            if !boundary {
                phi = phi.rem_euclid(TAU);
            }
            let fan: Vec<usize> = mesh.outgoing_halfedges(v).collect();
            // The wedge: the last halfedge whose signpost is not above phi (the first if none).
            let mut h = fan[0];
            for &g in &fan {
                if input.signpost(g) <= phi {
                    h = g;
                }
            }
            let offset = (phi - input.signpost(h)) / scale;
            let [a, b] = [mesh.from_vertex(h), mesh.to_vertex(h)];
            let c = mesh.to_vertex(mesh.next(h));
            let e = sub3(positions[b], positions[a]);
            let normal = cross3(e, sub3(positions[c], positions[a]));
            let (e, normal) = (unit3(e), unit3(normal));
            let perp = cross3(normal, e);
            let (s, co) = offset.sin_cos();
            let len = u.abs();
            out.push([0, 1, 2].map(|i| len * (co * e[i] + s * perp[i])));
        }
        Ok(out)
    }
}

fn frame_factor_input(theta: f64, boundary: bool) -> f64 {
    if boundary { theta / PI } else { 1.0 }
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit3(a: [f64; 3]) -> [f64; 3] {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if n > 0.0 { a.map(|x| x / n) } else { a }
}

// ----- the smoothest direction field -----

/// Options of [`smoothest_field`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectionFieldOptions {
    /// The symmetry `n` (1: vectors, 2: line fields, 4: cross fields).
    pub n_sym: u32,
    /// The discretization of the energy.
    pub laplacian: ConnectionLaplacian,
    /// The limit on inverse iterations.
    pub max_iterations: usize,
    /// Stop when the relative eigen-residual `‖A u − λ M u‖ / (‖A u‖ + (|λ| + ε) ‖M u‖)` is
    /// below this.
    pub tolerance: f64,
    /// The shift `ε` of GODF Alg. 1 (`A + εM`), relative: `ε = shift · tr A / tr M`. It
    /// makes the matrix definite when the field can be parallel (a flat torus) and does not
    /// change the eigenvectors.
    pub shift: f64,
    /// The number of complex start vectors of the block inverse iteration (each with its
    /// multiple by `i`); 1 is GODF Alg. 2. More converge faster when the smallest eigenvalues
    /// nearly coincide.
    pub block: usize,
    /// Seed of the start vectors (deterministic).
    pub seed: u64,
}

impl Default for DirectionFieldOptions {
    /// Vectors, the finite element energy (the paper's choice for Fig. 21), a block of 4, at
    /// most 1000 iterations to a residual of `1e-10`, shift `1e-10`, seed 1.
    fn default() -> Self {
        DirectionFieldOptions {
            n_sym: 1,
            laplacian: ConnectionLaplacian::FiniteElement,
            max_iterations: 1000,
            tolerance: 1e-10,
            shift: 1e-10,
            block: 4,
            seed: 1,
        }
    }
}

/// How the inverse iteration ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EigenStatus {
    /// The residual is below the tolerance.
    Converged,
    /// The iteration limit was reached first.
    IterationLimit,
}

/// What [`smoothest_field`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectionFieldReport {
    /// How it ended.
    pub status: EigenStatus,
    /// Block inverse iterations made.
    pub iterations: usize,
    /// The Rayleigh quotient `λ = u* A u / u* M u` of the result: the Dirichlet energy of the
    /// field normalized to `‖u‖_{L²} = 1` (unshifted).
    pub eigenvalue: f64,
    /// `λ` times the total area: the same energy for a field of mean square length one,
    /// independent of the scale of the surface.
    pub energy: f64,
    /// The final relative residual.
    pub residual: f64,
    /// The absolute shift `ε` used.
    pub shift: f64,
}

/// The smoothest `n`-direction field (GODF Sec. 6.1, Alg. 2): the eigenvector of
/// `A u = λ M u` of the smallest eigenvalue, normalized to `u* M u = 1`.
///
/// # Errors
/// [`VectorError::InvalidOptions`]; [`VectorError::NotPositiveDefinite`] if `A + εM` is not
/// (the finite-difference energy on a triangulation with negative cotan weights);
/// [`VectorError::Sparse`].
pub fn smoothest_field(
    tri: &IntrinsicTriangulation,
    options: &DirectionFieldOptions,
) -> Result<(TangentField, DirectionFieldReport), VectorError> {
    if !(options.tolerance > 0.0 && options.shift >= 0.0 && options.shift.is_finite()) {
        return Err(VectorError::InvalidOptions {
            reason: format!(
                "the tolerance must be positive and the shift non-negative, got {} and {}",
                options.tolerance, options.shift
            ),
        });
    }
    let matrices = connection_matrices(tri, options.n_sym, options.laplacian)?;
    let (a, m) = (&matrices.stiffness, &matrices.mass);
    let trace_ratio = a.diagonal().iter().sum::<f64>() / m.diagonal().iter().sum::<f64>();
    let eps = options.shift * trace_ratio;
    let shifted = {
        let mut triplets: Vec<(u32, u32, f64)> = a.triplets().collect();
        triplets.extend(m.triplets().map(|(i, j, v)| (i, j, eps * v)));
        SparseMatrix::from_triplets(a.nrows(), a.ncols(), &triplets)?
    };
    let factor = match linsolve::factor_spd(&shifted) {
        Ok(f) => f,
        Err(SparseError::NotPositiveDefinite { .. }) => {
            return Err(VectorError::NotPositiveDefinite {
                what: format!(
                    "the shifted energy matrix A + εM ({:?} discretization)",
                    options.laplacian
                ),
            });
        }
        Err(e) => return Err(e.into()),
    };
    let pair = eigen::smallest(
        a,
        m,
        &factor,
        options.block,
        options.max_iterations,
        options.tolerance,
        eps,
        options.seed,
    )?;
    let (u, lambda, residual, iterations) = (pair.vector, pair.value, pair.residual, pair.iterations);
    let status = if pair.converged {
        EigenStatus::Converged
    } else {
        EigenStatus::IterationLimit
    };
    let field = TangentField {
        n_sym: options.n_sym,
        values: from_real(&u),
    };
    Ok((
        field,
        DirectionFieldReport {
            status,
            iterations,
            eigenvalue: lambda,
            energy: lambda * tri.total_area(),
            residual,
            shift: eps,
        },
    ))
}

// ----- singularities -----

/// A face around which the field turns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Singularity {
    /// The face.
    pub face: u32,
    /// The index in units of `1/n` (GODF Eq. 19): `+1` is a singularity of index `1/n`.
    pub index: i32,
}

/// The singularities of a field and their sum.
#[derive(Debug, Clone, PartialEq)]
pub struct Singularities {
    /// Every face with a non-zero index, ascending.
    pub faces: Vec<Singularity>,
    /// The sum of all indices (units of `1/n`); on a closed surface `n χ` (Poincaré–Hopf).
    pub index_sum: i64,
    /// Faces left out because the field vanishes at one of their corners.
    pub n_skipped: usize,
    /// The largest distance of `(ω_ij + ω_jk + ω_ki + nΩ)/2π` from an integer: rounding only.
    pub max_rounding: f64,
}

/// The index of every face (GODF Eq. 19, with the face holonomy of [`face_holonomy`]):
/// `(ω_ij + ω_jk + ω_ki + n Ω_ijk) / 2π`, where `ω_ij ∈ (−π, π]` is the angle with
/// `u_j/|u_j| = e^{iω_ij} r_ij u_i/|u_i|`.
///
/// [formula: app-vector-singularity-index]
///
/// # Errors
/// [`VectorError::InvalidInput`] if the field does not have one value per vertex.
pub fn singularities(
    tri: &IntrinsicTriangulation,
    field: &TangentField,
) -> Result<Singularities, VectorError> {
    if field.values.len() != tri.n_vertices() {
        return Err(VectorError::InvalidInput {
            reason: format!(
                "a field with {} values on a triangulation with {} vertices",
                field.values.len(),
                tri.n_vertices()
            ),
        });
    }
    let n = field.n_sym;
    let mut faces = Vec::new();
    let mut index_sum = 0i64;
    let mut n_skipped = 0;
    let mut max_rounding = 0.0f64;
    for f in 0..tri.n_faces() {
        let verts = tri.faces()[f];
        if verts.iter().any(|&v| field.values[v as usize].abs() == 0.0) {
            n_skipped += 1;
            continue;
        }
        let mut total = n as f64 * face_holonomy(tri, f);
        for k in 0..3 {
            let h = 3 * f + k;
            let (i, j) = (tri.from_vertex(h), tri.to_vertex(h));
            let ui = field.values[i];
            let uj = field.values[j];
            total += (uj / (transport(tri, h, n) * ui)).arg();
        }
        let x = total / TAU;
        let index = x.round();
        max_rounding = max_rounding.max((x - index).abs());
        if index != 0.0 {
            faces.push(Singularity {
                face: f as u32,
                index: index as i32,
            });
            index_sum += index as i64;
        }
    }
    Ok(Singularities {
        faces,
        index_sum,
        n_skipped,
        max_rounding,
    })
}

// ----- harmonic fields with boundary values -----

/// The vectors prescribed at the boundary for [`harmonic_field`].
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryValues {
    /// The unit tangent of the boundary (the direction in which the surface is on the left),
    /// to the `n`-th power.
    Tangent,
    /// The unit inward normal of the boundary, to the `n`-th power.
    InwardNormal,
    /// Values at the boundary vertices of the input mesh: `(input vertex, n-th power in the
    /// frame of that vertex)`, every input boundary vertex exactly once. A vertex inserted on
    /// the boundary gets the linear interpolation, by length along the boundary, of the values
    /// at the input vertices before and after it, both transported to it along the boundary.
    AtInputVertices(Vec<(u32, C64)>),
}

/// Options of [`harmonic_field`].
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicFieldOptions {
    /// The symmetry `n`.
    pub n_sym: u32,
    /// The discretization (the paper uses finite differences here).
    pub laplacian: ConnectionLaplacian,
    /// The boundary values.
    pub boundary: BoundaryValues,
}

impl Default for HarmonicFieldOptions {
    /// Vectors, finite differences, the inward normal on the boundary.
    fn default() -> Self {
        HarmonicFieldOptions {
            n_sym: 1,
            laplacian: ConnectionLaplacian::FiniteDifference,
            boundary: BoundaryValues::InwardNormal,
        }
    }
}

/// What [`harmonic_field`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicFieldReport {
    /// Interior vertices (unknowns).
    pub n_interior: usize,
    /// Boundary vertices (prescribed).
    pub n_boundary: usize,
    /// The factorization: Cholesky, or LU when the interior block was not positive definite.
    pub factorization: Factorized,
    /// `‖L_II u_I + L_IB u_B‖ / ‖L_IB u_B‖`.
    pub residual: f64,
}

/// The boundary value of every boundary vertex of `tri` (zero elsewhere).
fn boundary_values(
    tri: &IntrinsicTriangulation,
    n_sym: u32,
    values: &BoundaryValues,
) -> Result<Vec<C64>, VectorError> {
    let n = n_sym as f64;
    let mut out = vec![C64::default(); tri.n_vertices()];
    let boundary: Vec<usize> = (0..tri.n_vertices())
        .filter(|&v| tri.is_boundary_vertex(v))
        .collect();
    match values {
        BoundaryValues::Tangent | BoundaryValues::InwardNormal => {
            for &v in &boundary {
                let theta = tri.vertex_angle_sum(v);
                let alpha = match values {
                    BoundaryValues::Tangent => 0.5 * (theta - PI),
                    _ => 0.5 * theta,
                };
                out[v] = C64::from_angle(n * alpha);
            }
        }
        BoundaryValues::AtInputVertices(given) => {
            let n_input = tri.n_input_vertices();
            let mut set = vec![false; n_input];
            for &(v, u) in given {
                let v = v as usize;
                if v >= n_input || !tri.is_boundary_vertex(v) {
                    return Err(VectorError::InvalidInput {
                        reason: format!("vertex {v} is not a boundary vertex of the input mesh"),
                    });
                }
                if set[v] {
                    return Err(VectorError::InvalidInput {
                        reason: format!("boundary vertex {v} is given twice"),
                    });
                }
                if !(u.re.is_finite() && u.im.is_finite()) {
                    return Err(VectorError::InvalidInput {
                        reason: format!("the value at boundary vertex {v} is not finite"),
                    });
                }
                set[v] = true;
                out[v] = u;
            }
            if let Some(v) = (0..n_input).find(|&v| tri.is_boundary_vertex(v) && !set[v]) {
                return Err(VectorError::InvalidInput {
                    reason: format!("no value given for input boundary vertex {v}"),
                });
            }
            // Inserted boundary vertices: walk each run between two input vertices.
            for &start in boundary.iter().filter(|&&v| v < n_input) {
                // From an input vertex forward to the next input vertex.
                let mut run: Vec<(usize, f64, C64)> = Vec::new(); // vertex, length, transport
                let mut v = start;
                let mut length = 0.0;
                let mut carry = C64::new(1.0, 0.0);
                loop {
                    let h = tri.outgoing_halfedge(v);
                    length += tri.edge_length(tri.edge(h));
                    carry = transport(tri, h, n_sym) * carry;
                    v = tri.to_vertex(h);
                    if v < n_input {
                        break;
                    }
                    run.push((v, length, carry));
                    if run.len() > tri.n_vertices() {
                        return Err(VectorError::InvalidInput {
                            reason: "a boundary loop of inserted vertices only".to_string(),
                        });
                    }
                }
                let end = v;
                let total = length;
                let to_end = carry; // transport from start to end
                for &(w, s, from_start) in &run {
                    // Transport from `end` back to `w` is the inverse of w → end.
                    let w_to_end = to_end / from_start;
                    let a = from_start * out[start];
                    let b = out[end] / w_to_end;
                    let t = s / total;
                    out[w] = a.scale(1.0 - t) + b.scale(t);
                }
            }
        }
    }
    Ok(out)
}

/// The harmonic `n`-vector field with the given boundary values: `(L∇ u)_i = 0` at every
/// interior vertex, `u` prescribed on the boundary (paper Sec. 5.5, "flip-free vector
/// fields").
///
/// # Errors
/// [`VectorError::InvalidInput`] for a closed surface, a surface without interior vertices,
/// or boundary values that do not fit; [`VectorError::Sparse`].
pub fn harmonic_field(
    tri: &IntrinsicTriangulation,
    options: &HarmonicFieldOptions,
) -> Result<(TangentField, HarmonicFieldReport), VectorError> {
    let matrices = connection_matrices(tri, options.n_sym, options.laplacian)?;
    let ub = boundary_values(tri, options.n_sym, &options.boundary)?;
    let interior: Vec<usize> = (0..tri.n_vertices())
        .filter(|&v| !tri.is_boundary_vertex(v))
        .collect();
    let bnd: Vec<usize> = (0..tri.n_vertices())
        .filter(|&v| tri.is_boundary_vertex(v))
        .collect();
    if bnd.is_empty() {
        return Err(VectorError::InvalidInput {
            reason: "the surface has no boundary to prescribe vectors on".to_string(),
        });
    }
    if interior.is_empty() {
        return Err(VectorError::InvalidInput {
            reason: "the surface has no interior vertex".to_string(),
        });
    }
    let real_index = |vs: &[usize]| -> Vec<u32> {
        vs.iter()
            .flat_map(|&v| [2 * v as u32, 2 * v as u32 + 1])
            .collect()
    };
    let (ri, rb) = (real_index(&interior), real_index(&bnd));
    let l_ii = matrices.stiffness.submatrix(&ri, &ri)?;
    let l_ib = matrices.stiffness.submatrix(&ri, &rb)?;
    let xb: Vec<f64> = bnd.iter().flat_map(|&v| [ub[v].re, ub[v].im]).collect();
    let lb = linsolve::mul(&l_ib, &xb);
    let mut rhs: Vec<f64> = lb.iter().map(|x| -x).collect();
    let (factor, factorization) = linsolve::factor(&l_ii)?;
    linsolve::solve(&factor, &mut rhs)?;
    let check = linsolve::mul(&l_ii, &rhs);
    let res: Vec<f64> = check.iter().zip(&lb).map(|(x, y)| x + y).collect();
    let residual = linsolve::norm(&res) / linsolve::norm(&lb).max(f64::MIN_POSITIVE);
    let mut values = ub;
    for (k, &v) in interior.iter().enumerate() {
        values[v] = C64::new(rhs[2 * k], rhs[2 * k + 1]);
    }
    Ok((
        TangentField {
            n_sym: options.n_sym,
            values,
        },
        HarmonicFieldReport {
            n_interior: interior.len(),
            n_boundary: bnd.len(),
            factorization,
            residual,
        },
    ))
}

// ----- the vector maximum principle -----

/// The result of [`check_maximum_principle`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaximumPrincipleReport {
    /// Interior vertices checked (those with a non-zero value).
    pub n_checked: usize,
    /// Interior vertices whose value is not in the closed convex cone of the transported values
    /// of their neighbours (beyond `tolerance` radians).
    pub n_outside_cone: usize,
    /// The smallest margin over the checked vertices: how far, in radians, the value is inside
    /// the cone (negative: outside). `π` for a vertex whose neighbours span more than a half
    /// plane in every direction (the cone is the plane); `+∞` if nothing was checked.
    pub min_margin: f64,
    /// The vertex with the smallest margin.
    pub worst_vertex: Option<u32>,
}

/// Checks the vector maximum principle (VHM Sec. 5.4) at every interior vertex `i`: `u_i` lies
/// in the convex cone of `{r_ji u_j}` over the neighbours `j` (the `n`-th powers, transported
/// to `i`). It holds for a harmonic field ([`harmonic_field`]) when every cotan weight at `i` is
/// non-negative, as on an intrinsic Delaunay triangulation, since then `u_i` is the weighted
/// mean `Σ w_ij r_ji u_j / Σ w_ij`.
///
/// The margin of `u_i` is `π − (a⁺ − a⁻)`, where `a⁺ ≥ 0` is the smallest and `a⁻ ≤ 0` the
/// largest angle of a neighbour's value relative to `u_i` (in `(−π, π]`); `u_i` is in the cone
/// exactly when some neighbour points its way or the margin is positive. A vertex counts as
/// outside when its margin is below `−tolerance`.
///
/// [formula: app-vector-maximum-principle]
pub fn check_maximum_principle(
    tri: &IntrinsicTriangulation,
    field: &TangentField,
    tolerance: f64,
) -> MaximumPrincipleReport {
    let mut report = MaximumPrincipleReport {
        n_checked: 0,
        n_outside_cone: 0,
        min_margin: f64::INFINITY,
        worst_vertex: None,
    };
    for i in 0..tri.n_vertices() {
        if tri.is_boundary_vertex(i) {
            continue;
        }
        let ui = field.values[i];
        if ui.abs() == 0.0 {
            continue;
        }
        report.n_checked += 1;
        let mut a_plus = f64::INFINITY;
        let mut a_minus = f64::NEG_INFINITY;
        for h in tri.outgoing_halfedges(i) {
            let j = tri.to_vertex(h);
            let vj = transport(tri, h, field.n_sym).conj() * field.values[j];
            if vj.abs() == 0.0 {
                continue;
            }
            let a = (vj / ui).arg();
            if a >= 0.0 {
                a_plus = a_plus.min(a);
            }
            if a <= 0.0 {
                a_minus = a_minus.max(a);
            }
        }
        let margin = if a_plus == 0.0 || a_minus == 0.0 {
            PI
        } else if a_plus.is_finite() && a_minus.is_finite() {
            PI - (a_plus - a_minus)
        } else {
            // Every neighbour on one side: outside by the angle to the nearest.
            -(a_plus.min(-a_minus))
        };
        if margin < -tolerance {
            report.n_outside_cone += 1;
        }
        if margin < report.min_margin {
            report.min_margin = margin;
            report.worst_vertex = Some(i as u32);
        }
    }
    report
}

/// The tangent-space scale of every vertex (`2π/Θ` inside, `1` on the boundary), for callers
/// that convert between corner angles and tangent angles.
pub fn tangent_scales(tri: &IntrinsicTriangulation) -> Vec<f64> {
    (0..tri.n_vertices()).map(|v| tangent_scale(tri, v)).collect()
}

/// The direction of every halfedge at its tail in the tangent frame of the tail (the signpost,
/// unnormalized at boundary vertices).
pub fn halfedge_directions(tri: &IntrinsicTriangulation) -> Vec<f64> {
    (0..tri.n_halfedges())
        .map(|h| tri.signpost(h) * frame_factor(tri, tri.from_vertex(h)))
        .collect()
}
