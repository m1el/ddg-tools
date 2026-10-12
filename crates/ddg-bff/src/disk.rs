//! Boundary First Flattening of a topological disk: the three steps of the paper (Sec. 5) and
//! the boundary data for automatic flattening, prescribed scale factors, prescribed angles and
//! the unit disk (Sec. 6.1, 6.2, 6.5).
//!
//! Ported from `BFF::flatten`, `BFF::flattenToDisk` and their helpers in `Bff.cpp` of the
//! reference implementation (MIT, © 2017 Rohan Sawhney and Keenan Crane).
//!
//! # Conventions
//!
//! * `L` is the positive semidefinite cotan Laplacian (`A` in the paper).
//! * Boundary data is indexed by position along the boundary loop: entry `i` belongs to
//!   [`DiskBff::boundary`]`[i]`, and boundary edge `i` runs from that vertex to the next one.
//!   The loop is walked with the surface on its left.
//! * `u` is the log scale factor: lengths grow by `e^u`. `k` is the exterior angle of the
//!   boundary (integrated geodesic curvature), `Ω` the angle defect.
//! * The Cherrier formula in integrated form (paper Eq. 11 and 12, with the sign of `L` used
//!   here): `(L u)_i = −Ω_i` at interior vertices and `k̃_i = k_i − h_i` at boundary vertices,
//!   where `h` is the Neumann data of `u`.
//!
//! # Where the layout sits
//!
//! Nothing is recentred, rotated or rescaled after the fact. The first boundary vertex
//! ([`DiskBff::boundary`]`[0]`) is at the origin and the first boundary edge points along
//! `+x`. With [`Extension::Conjugate`] the second coordinate is determined only up to a
//! constant, which is chosen such that it agrees on average, over the boundary vertices, with
//! the boundary curve.

use std::f64::consts::PI;

use ddg_core::geometry::Geometry;
use ddg_core::topology::TriMesh;

use crate::chart::{
    ChartReport, PoissonSystem, Surface, check_finite_layout, halfedge_lengths, layout_quality,
};
use crate::curve::{
    CurveOptions, LoopData, best_fit_curves, dual_boundary_lengths, target_boundary_lengths,
};
use crate::error::{BffError, check_values};

/// The diagonal shift the reference implementation adds to its Laplacian (`1e-8`).
pub const REFERENCE_DIAGONAL_SHIFT: f64 = 1e-8;

/// How the boundary curve is extended to the interior (paper Sec. 4.4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Extension {
    /// Both coordinates are harmonic functions that interpolate the boundary curve. The
    /// boundary of the layout is exactly the constructed curve, so the prescribed exterior
    /// angles are reproduced exactly.
    Harmonic,
    /// The first coordinate interpolates the boundary curve and the second is its harmonic
    /// conjugate (discrete Hilbert transform). The map is as conformal as a piecewise linear
    /// map can be, but its boundary only approximates the constructed curve.
    Conjugate,
}

/// Options of [`DiskBff`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiskOptions {
    /// A constant added to the diagonal of the Laplacian before factorization.
    ///
    /// `0.0` (the default) solves the Poisson problems of the method exactly; the pure Neumann
    /// problem, which is singular, is then solved with one vertex pinned. The reference
    /// implementation instead adds [`REFERENCE_DIAGONAL_SHIFT`] to make every matrix positive
    /// definite, which perturbs each solve by about `shift · n_vertices / 20` relative; set
    /// this to reproduce its numbers.
    pub diagonal_shift: f64,
    /// Prescribed exterior angles must sum to `2π` within this absolute tolerance.
    pub angle_sum_tolerance: f64,
    /// Largest fraction by which a boundary length may be changed to close the boundary
    /// curve. If closing needs more, the curve is closed as far as this bound allows and the
    /// remaining gap is spread along it; see [`crate::curve::ClosureReport`].
    pub max_length_adjustment: f64,
    /// The boundary vertex at which the boundary loop starts: it becomes
    /// [`DiskBff::boundary`]`[0]` and is placed at the origin. `None` (the default) starts
    /// where `TriMesh::boundary_loops` does. For exact solves the layout depends on this
    /// choice only by a rigid motion.
    pub boundary_start: Option<u32>,
}

impl Default for DiskOptions {
    fn default() -> Self {
        DiskOptions {
            diagonal_shift: 0.0,
            angle_sum_tolerance: 1e-8,
            max_length_adjustment: 0.5,
            boundary_start: None,
        }
    }
}

/// A flattening of a disk.
#[derive(Debug, Clone, PartialEq)]
pub struct DiskFlattening {
    /// The position of every vertex in the plane.
    pub uv: Vec<[f64; 2]>,
    /// The log scale factors at the boundary vertices that the layout was built from
    /// (prescribed, or computed from the prescribed angles), by boundary position.
    pub boundary_scale_factors: Vec<f64>,
    /// The exterior angles of the boundary curve (prescribed, or computed from the prescribed
    /// scale factors), by boundary position. They sum to `2π`.
    pub exterior_angles: Vec<f64>,
    /// How the interior was filled in.
    pub extension: Extension,
    /// Largest relative residual of the linear solves that produced the layout.
    pub solve_residual: f64,
    /// Closure of the boundary curve, distortion and overlap.
    pub report: ChartReport,
}

/// The pure Neumann problem `L a = b`, made non-singular either by a diagonal shift or by
/// pinning one vertex to zero.
#[derive(Debug)]
struct NeumannSystem {
    system: PoissonSystem,
    pinned: Option<usize>,
}

/// Boundary First Flattening of one disk: set up once, then flatten with any boundary data.
///
/// The factorization of the Dirichlet Laplacian is computed on first use and reused; so is
/// the factorization of the Neumann Laplacian, which only prescribed angles and
/// [`Extension::Conjugate`] need.
#[derive(Debug)]
pub struct DiskBff<'a> {
    mesh: &'a TriMesh,
    geometry: &'a Geometry,
    options: DiskOptions,
    halfedge_lengths: Vec<f64>,
    /// Boundary vertices in loop order.
    boundary: Vec<u32>,
    is_boundary: Vec<bool>,
    /// Length of boundary edge `i`.
    boundary_lengths: Vec<f64>,
    /// Angle defect per vertex (0 on the boundary).
    omega: Vec<f64>,
    /// Exterior angle per boundary position.
    kappa: Vec<f64>,
    mass: Vec<f64>,
    dirichlet: Option<PoissonSystem>,
    neumann: Option<NeumannSystem>,
}

impl<'a> DiskBff<'a> {
    /// Prepares the flattening of `mesh` with the metric `geometry`.
    ///
    /// # Errors
    /// [`BffError::Geometry`] if the geometry does not belong to the mesh;
    /// [`BffError::NotADisk`] unless the mesh is connected with Euler characteristic 1 and one
    /// boundary loop; [`BffError::StartNotOnBoundary`] for a bad `boundary_start`.
    pub fn new(
        mesh: &'a TriMesh,
        geometry: &'a Geometry,
        options: DiskOptions,
    ) -> Result<Self, BffError> {
        geometry.check_mesh(mesh)?;
        if mesh.n_components() != 1
            || mesh.euler_characteristic() != 1
            || mesh.boundary_loops().len() != 1
        {
            return Err(BffError::NotADisk {
                components: mesh.n_components(),
                euler_characteristic: mesh.euler_characteristic(),
                boundary_loops: mesh.boundary_loops().len(),
            });
        }
        if !options.diagonal_shift.is_finite() || options.diagonal_shift < 0.0 {
            return Err(BffError::NonFiniteInput {
                what: "diagonal_shift (must be finite and non-negative)",
                index: 0,
                value: options.diagonal_shift,
            });
        }
        let mut boundary = mesh.boundary_loops()[0].clone();
        if let Some(start) = options.boundary_start {
            let position = boundary
                .iter()
                .position(|&v| v == start)
                .ok_or(BffError::StartNotOnBoundary { vertex: start })?;
            boundary.rotate_left(position);
        }
        let mut is_boundary = vec![false; mesh.n_vertices()];
        for &v in &boundary {
            is_boundary[v as usize] = true;
        }
        let halfedge_lengths = halfedge_lengths(mesh, geometry);
        let boundary_lengths: Vec<f64> = boundary
            .iter()
            .map(|&v| halfedge_lengths[mesh.outgoing_halfedge(v as usize)])
            .collect();
        let mut disk = DiskBff {
            mesh,
            geometry,
            options,
            halfedge_lengths,
            boundary,
            is_boundary,
            boundary_lengths,
            omega: Vec::new(),
            kappa: Vec::new(),
            mass: Vec::new(),
            dirichlet: None,
            neumann: None,
        };
        let surface = disk.surface();
        let angle_sums = surface.angle_sums();
        let mass = surface.masses();
        disk.omega = curvatures(&angle_sums, &disk.is_boundary);
        disk.kappa = disk
            .boundary
            .iter()
            .map(|&v| PI - angle_sums[v as usize])
            .collect();
        disk.mass = mass;
        Ok(disk)
    }

    fn surface(&self) -> Surface<'_> {
        Surface {
            mesh: self.mesh,
            cotans: self.geometry.corner_cotans(),
            angles: self.geometry.corner_angles(),
            areas: self.geometry.face_areas(),
            halfedge_lengths: &self.halfedge_lengths,
        }
    }

    /// The boundary vertices in loop order, with the surface on the left. All boundary data
    /// of this type is indexed by position in this list.
    pub fn boundary(&self) -> &[u32] {
        &self.boundary
    }

    /// Lengths of the boundary edges: entry `i` is the edge from boundary vertex `i` to
    /// `i + 1`.
    pub fn boundary_edge_lengths(&self) -> &[f64] {
        &self.boundary_lengths
    }

    /// The exterior angles `k_i = π − Σ θ` of the boundary of the input surface, by boundary
    /// position.
    ///
    /// [formula: boundary-geodesic-curvature]
    pub fn geodesic_curvatures(&self) -> &[f64] {
        &self.kappa
    }

    /// The angle defects `Ω_i = 2π − Σ θ` per vertex (zero at boundary vertices).
    ///
    /// [formula: angle-defect]
    pub fn angle_defects(&self) -> &[f64] {
        &self.omega
    }

    fn dirichlet(&mut self) -> Result<(), BffError> {
        if self.dirichlet.is_none() {
            self.dirichlet = Some(PoissonSystem::new(
                "the Dirichlet Laplacian of the disk",
                &self.surface(),
                &self.is_boundary,
                self.options.diagonal_shift,
            )?);
        }
        Ok(())
    }

    fn neumann(&mut self) -> Result<(), BffError> {
        if self.neumann.is_none() {
            let mut fixed = vec![false; self.mesh.n_vertices()];
            let pinned = if self.options.diagonal_shift > 0.0 {
                None
            } else {
                let p = *self.boundary.last().expect("a loop has vertices") as usize;
                fixed[p] = true;
                Some(p)
            };
            let system = PoissonSystem::new(
                "the Neumann Laplacian of the disk",
                &self.surface(),
                &fixed,
                self.options.diagonal_shift,
            )?;
            self.neumann = Some(NeumannSystem { system, pinned });
        }
        Ok(())
    }

    /// The Dirichlet-to-Neumann map of the Poisson problem (paper Eq. 14, Alg. 3): given a
    /// source `φ` (one value per vertex) and boundary values `g`, solves
    /// `L_II a_I = φ_I − L_IB g` and returns the Neumann data
    /// `h = φ_B − L_IBᵀ a_I − L_BB g` together with the solution `a` at all vertices
    /// (`a = g` on the boundary).
    ///
    /// # Errors
    /// [`BffError::LengthMismatch`], [`BffError::NonFiniteInput`], [`BffError::Solve`].
    ///
    /// [formula: bff-dirichlet-to-neumann]
    pub fn dirichlet_to_neumann(
        &mut self,
        source: &[f64],
        boundary_values: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), BffError> {
        check_values("source", source, self.mesh.n_vertices())?;
        check_values("boundary values", boundary_values, self.boundary.len())?;
        self.dirichlet()?;
        let mut a = source.to_vec();
        for (&v, &g) in self.boundary.iter().zip(boundary_values) {
            a[v as usize] = g;
        }
        let surface = self.surface();
        self.dirichlet
            .as_ref()
            .expect("factorized above")
            .solve(&surface, &mut a, 1)?;
        let mut product = vec![0.0; a.len()];
        surface.apply_laplacian(&a, &mut product);
        let h = self
            .boundary
            .iter()
            .map(|&v| source[v as usize] - product[v as usize])
            .collect();
        Ok((h, a))
    }

    /// The Neumann-to-Dirichlet map (paper Eq. 5, Alg. 4): solves `L a = φ − [0; h]` for a
    /// source `φ` (per vertex) and Neumann data `h` (per boundary position) and returns `a` at
    /// all vertices. The problem has a solution only if `Σ φ = Σ h`, and then only up to an
    /// additive constant; the returned solution has area-weighted mean zero
    /// (`Σ m_i a_i = 0` with the lumped masses `m`).
    ///
    /// # Errors
    /// [`BffError::LengthMismatch`], [`BffError::NonFiniteInput`], [`BffError::Solve`].
    ///
    /// [formula: bff-neumann-to-dirichlet]
    /// [formula: closed-mesh-optimal-constant]
    pub fn neumann_to_dirichlet(
        &mut self,
        source: &[f64],
        neumann: &[f64],
    ) -> Result<Vec<f64>, BffError> {
        check_values("source", source, self.mesh.n_vertices())?;
        check_values("Neumann data", neumann, self.boundary.len())?;
        let mut a = source.to_vec();
        for (&v, &h) in self.boundary.iter().zip(neumann) {
            a[v as usize] -= h;
        }
        self.solve_neumann(&mut a)?;
        let total: f64 = self.mass.iter().sum();
        let mean: f64 = a.iter().zip(&self.mass).map(|(x, m)| x * m).sum::<f64>() / total;
        for x in &mut a {
            *x -= mean;
        }
        Ok(a)
    }

    /// Solves `L a = rhs` in place (pinned or shifted). Returns the relative residual.
    fn solve_neumann(&mut self, values: &mut [f64]) -> Result<f64, BffError> {
        self.neumann()?;
        let surface = self.surface();
        let neumann = self.neumann.as_ref().expect("factorized above");
        if let Some(p) = neumann.pinned {
            values[p] = 0.0;
        }
        neumann.system.solve(&surface, values, 1)
    }

    /// The exterior angles that are compatible with the given boundary scale factors:
    /// `k̃ = k − h`, where `h` is the Neumann data of the solution of `L u = −Ω` with `u` on
    /// the boundary as given (Step I of the paper, discrete Cherrier formula, Eq. 12). They
    /// sum to `2π` up to rounding (paper Proposition 1). Also returns `u` at all vertices.
    ///
    /// [formula: bff-cherrier-boundary-curvature]
    pub fn compatible_exterior_angles(
        &mut self,
        boundary_scale_factors: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), BffError> {
        let source: Vec<f64> = self.omega.iter().map(|&w| -w).collect();
        let (h, u) = self.dirichlet_to_neumann(&source, boundary_scale_factors)?;
        let k_tilde = self.kappa.iter().zip(&h).map(|(k, h)| k - h).collect();
        Ok((k_tilde, u))
    }

    /// The scale factors (at all vertices) that are compatible with the given exterior
    /// angles: the solution of `L u = −Ω` with Neumann data `h = k − k̃` (Step I of the paper).
    /// The additive constant, which only sets the size of the layout, is chosen such that the
    /// area-weighted mean of `u` is zero; this is the choice with the least area distortion
    /// `Σ m_i u_i²`.
    ///
    /// # Errors
    /// [`BffError::ExteriorAnglesDoNotClose`] if the angles do not sum to `2π`.
    ///
    /// [formula: bff-cherrier-boundary-curvature]
    pub fn compatible_scale_factors(
        &mut self,
        exterior_angles: &[f64],
    ) -> Result<Vec<f64>, BffError> {
        check_values("exterior angles", exterior_angles, self.boundary.len())?;
        let sum: f64 = exterior_angles.iter().sum();
        if (sum - 2.0 * PI).abs() > self.options.angle_sum_tolerance {
            return Err(BffError::ExteriorAnglesDoNotClose {
                sum,
                tolerance: self.options.angle_sum_tolerance,
            });
        }
        let source: Vec<f64> = self.omega.iter().map(|&w| -w).collect();
        let h: Vec<f64> = self
            .kappa
            .iter()
            .zip(exterior_angles)
            .map(|(k, kt)| k - kt)
            .collect();
        self.neumann_to_dirichlet(&source, &h)
    }

    /// Flattens with prescribed log scale factors at the boundary vertices (by boundary
    /// position). The reference implementation uses [`Extension::Conjugate`] for this case.
    ///
    /// # Errors
    /// [`BffError::LengthMismatch`], [`BffError::NonFiniteInput`], [`BffError::Solve`],
    /// [`BffError::NonFiniteResult`].
    pub fn flatten_with_scale_factors(
        &mut self,
        boundary_scale_factors: &[f64],
        extension: Extension,
    ) -> Result<DiskFlattening, BffError> {
        let (k_tilde, _) = self.compatible_exterior_angles(boundary_scale_factors)?;
        self.flatten_with_boundary_data(boundary_scale_factors, &k_tilde, extension)
    }

    /// Flattens with prescribed exterior angles at the boundary vertices (by boundary
    /// position), which must sum to `2π`. The size of the layout is fixed by the convention of
    /// [`DiskBff::compatible_scale_factors`]. The reference implementation uses
    /// [`Extension::Harmonic`] for this case, which reproduces the angles exactly.
    ///
    /// # Errors
    /// As for [`DiskBff::flatten_with_scale_factors`], and
    /// [`BffError::ExteriorAnglesDoNotClose`].
    pub fn flatten_with_exterior_angles(
        &mut self,
        exterior_angles: &[f64],
        extension: Extension,
    ) -> Result<DiskFlattening, BffError> {
        let u = self.compatible_scale_factors(exterior_angles)?;
        let boundary_u: Vec<f64> = self.boundary.iter().map(|&v| u[v as usize]).collect();
        self.flatten_with_boundary_data(&boundary_u, exterior_angles, extension)
    }

    /// The automatic flattening with the least area distortion: zero scale factors on the
    /// boundary (paper Sec. 6.1), extended holomorphically ([`Extension::Conjugate`]), as the
    /// reference implementation does by default.
    pub fn flatten_automatic(&mut self) -> Result<DiskFlattening, BffError> {
        let zero = vec![0.0; self.boundary.len()];
        self.flatten_with_scale_factors(&zero, Extension::Conjugate)
    }

    /// Flattens onto a disk that is as round as possible (paper Sec. 6.5, `BFF::flattenToDisk`):
    /// starting from `u = 0`, the exterior angles are repeatedly set proportional to the
    /// current dual boundary lengths, `k̃_i = 2π ℓ̃_i / Σ ℓ̃`, and the compatible scale factors
    /// recomputed. The reference implementation uses 10 iterations.
    ///
    /// [formula: bff-uniformization-angles]
    pub fn flatten_to_disk(&mut self, iterations: usize) -> Result<DiskFlattening, BffError> {
        let n = self.boundary.len();
        let mut u = vec![0.0; n];
        let mut k_tilde = vec![2.0 * PI / n as f64; n];
        for _ in 0..iterations.max(1) {
            let target = target_boundary_lengths(&u, &self.boundary_lengths);
            let dual = dual_boundary_lengths(&target);
            let total: f64 = dual.iter().sum();
            for (k, &l) in k_tilde.iter_mut().zip(&dual) {
                *k = 2.0 * PI * l / total;
            }
            let full = self.compatible_scale_factors(&k_tilde)?;
            for (ui, &v) in u.iter_mut().zip(&self.boundary) {
                *ui = full[v as usize];
            }
        }
        self.flatten_with_boundary_data(&u, &k_tilde, Extension::Harmonic)
    }

    /// Steps II and III of the paper for boundary data that is already complete: builds the
    /// closed curve with exterior angles `exterior_angles` and edge lengths close to
    /// `e^{(u_i + u_j)/2} ℓ_ij`, then extends it over the disk. The two arrays should be
    /// compatible (see [`DiskBff::compatible_exterior_angles`]); this is not checked.
    ///
    /// # Errors
    /// [`BffError::LengthMismatch`], [`BffError::NonFiniteInput`], [`BffError::Solve`],
    /// [`BffError::NonFiniteResult`].
    pub fn flatten_with_boundary_data(
        &mut self,
        boundary_scale_factors: &[f64],
        exterior_angles: &[f64],
        extension: Extension,
    ) -> Result<DiskFlattening, BffError> {
        let n_boundary = self.boundary.len();
        check_values("boundary scale factors", boundary_scale_factors, n_boundary)?;
        check_values("exterior angles", exterior_angles, n_boundary)?;
        let target = target_boundary_lengths(boundary_scale_factors, &self.boundary_lengths);
        let curve = best_fit_curves(
            &[LoopData {
                dof: (0..n_boundary as u32).collect(),
                length: self.boundary_lengths.clone(),
                target,
                exterior: exterior_angles.to_vec(),
            }],
            &CurveOptions {
                max_length_adjustment: self.options.max_length_adjustment,
                max_coupled_loops: 1,
            },
        )
        .pop()
        .expect("one loop in, one curve out");

        let n = self.mesh.n_vertices();
        self.dirichlet()?;
        let mut solve_residual;
        let uv: Vec<[f64; 2]> = match extension {
            Extension::Harmonic => {
                // Columns: real part, imaginary part. Source zero, boundary values the curve.
                let mut values = vec![0.0; 2 * n];
                for (&v, p) in self.boundary.iter().zip(&curve.positions) {
                    values[v as usize] = p[0];
                    values[n + v as usize] = p[1];
                }
                let surface = self.surface();
                solve_residual = self
                    .dirichlet
                    .as_ref()
                    .expect("factorized above")
                    .solve(&surface, &mut values, 2)?;
                (0..n).map(|v| [values[v], values[n + v]]).collect()
            }
            Extension::Conjugate => {
                let mut a = vec![0.0; n];
                for (&v, p) in self.boundary.iter().zip(&curve.positions) {
                    a[v as usize] = p[0];
                }
                let surface = self.surface();
                solve_residual = self
                    .dirichlet
                    .as_ref()
                    .expect("factorized above")
                    .solve(&surface, &mut a, 1)?;
                let mut b = self.hilbert_transform(&a);
                solve_residual = solve_residual.max(self.solve_neumann(&mut b)?);
                // The conjugate is determined up to a constant: match the curve on average.
                let offset = self
                    .boundary
                    .iter()
                    .zip(&curve.positions)
                    .map(|(&v, p)| p[1] - b[v as usize])
                    .sum::<f64>()
                    / n_boundary as f64;
                (0..n).map(|v| [a[v], b[v] + offset]).collect()
            }
        };
        check_finite_layout(self.mesh, &uv)?;
        let surface = self.surface();
        let (mut quality, _) = layout_quality(&surface, &uv);
        let (distortion, overlap) = quality.pop().expect("one component");
        Ok(DiskFlattening {
            uv,
            boundary_scale_factors: boundary_scale_factors.to_vec(),
            exterior_angles: exterior_angles.to_vec(),
            extension,
            solve_residual,
            report: ChartReport {
                first_face: 0,
                n_vertices: n,
                n_boundary_vertices: n_boundary,
                closure: curve.report,
                distortion,
                overlap,
            },
        })
    }

    /// The right-hand side of the Neumann problem for the harmonic conjugate `b` of `a`
    /// (paper Eq. 17): minimizing the conformal energy `E_C(a, b) = ½ (aᵀLa + bᵀLb) − A(a, b)`,
    /// with `A` the signed area of the image of the boundary, over `b` gives
    /// `(L b)_j = ½ (a_i − a_k)` at a boundary vertex `j` with predecessor `i` and successor
    /// `k` along the boundary (surface on the left), and `(L b)_j = 0` in the interior.
    ///
    /// [formula: bff-hilbert-transform]
    fn hilbert_transform(&self, a: &[f64]) -> Vec<f64> {
        let n_boundary = self.boundary.len();
        let mut rhs = vec![0.0; a.len()];
        for j in 0..n_boundary {
            let previous = self.boundary[(j + n_boundary - 1) % n_boundary] as usize;
            let next = self.boundary[(j + 1) % n_boundary] as usize;
            rhs[self.boundary[j] as usize] = 0.5 * (a[previous] - a[next]);
        }
        rhs
    }
}

/// Angle defects from angle sums: `2π − Σ θ` at interior vertices, zero on the boundary.
///
/// [formula: angle-defect]
pub(crate) fn curvatures(angle_sums: &[f64], is_boundary: &[bool]) -> Vec<f64> {
    angle_sums
        .iter()
        .zip(is_boundary)
        .map(|(&sum, &boundary)| if boundary { 0.0 } else { 2.0 * PI - sum })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quality::fit_similarity;
    use ddg_core::generators;

    fn setup(g: &generators::GeneratedMesh) -> (TriMesh, Geometry) {
        let mesh = g.topology().unwrap();
        let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
        (mesh, geometry)
    }

    #[test]
    fn a_flat_disk_maps_to_itself() {
        let mut g = generators::grid_disk(7, 5).unwrap();
        generators::perturb(&mut g.positions, [0.03, 0.03, 0.0], 5);
        let (mesh, geometry) = setup(&g);
        let flat: Vec<[f64; 2]> = g.positions.iter().map(|p| [p[0], p[1]]).collect();
        for shift in [0.0, REFERENCE_DIAGONAL_SHIFT] {
            let mut bff = DiskBff::new(
                &mesh,
                &geometry,
                DiskOptions {
                    diagonal_shift: shift,
                    ..DiskOptions::default()
                },
            )
            .unwrap();
            let k: Vec<f64> = bff.geodesic_curvatures().to_vec();
            let zero = vec![0.0; k.len()];
            let tolerance = if shift == 0.0 { 1e-11 } else { 1e-6 };
            for extension in [Extension::Harmonic, Extension::Conjugate] {
                let by_scale = bff.flatten_with_scale_factors(&zero, extension).unwrap();
                let by_angle = bff.flatten_with_exterior_angles(&k, extension).unwrap();
                for result in [&by_scale, &by_angle] {
                    let fit = fit_similarity(&result.uv, &flat).unwrap();
                    assert!(fit.max < tolerance, "{extension:?}: {fit:?}");
                    assert!((fit.scale - 1.0).abs() < tolerance);
                    assert_eq!(result.report.distortion.flipped_triangles, 0);
                    assert!(result.report.distortion.quasi_conformal.max < 1.0 + tolerance);
                    assert!(result.report.is_injective());
                    assert!(result.report.closure.gap_before < tolerance);
                }
                // Compatible data of a flat disk: angles k for u = 0 and u = 0 for angles k.
                for (a, b) in by_scale.exterior_angles.iter().zip(&k) {
                    assert!((a - b).abs() < tolerance);
                }
                for u in &by_angle.boundary_scale_factors {
                    assert!(u.abs() < tolerance);
                }
            }
        }
    }

    #[test]
    fn compatible_angles_sum_to_two_pi_and_round_trip() {
        let g = generators::bumpy_grid_disk(9, 8, 0.15).unwrap();
        let (mesh, geometry) = setup(&g);
        let mut bff = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap();
        let n = bff.boundary().len();
        let u: Vec<f64> = (0..n).map(|i| 0.2 * (i as f64 * 0.5).sin()).collect();
        let (k_tilde, u_full) = bff.compatible_exterior_angles(&u).unwrap();
        assert!((k_tilde.iter().sum::<f64>() - 2.0 * PI).abs() < 1e-11);
        // Neumann-to-Dirichlet inverts Dirichlet-to-Neumann up to the additive constant.
        let back = bff.compatible_scale_factors(&k_tilde).unwrap();
        let offset = back[0] - u_full[0];
        for (a, b) in back.iter().zip(&u_full) {
            assert!((a - b - offset).abs() < 1e-10);
        }
        let total: f64 = bff.mass.iter().sum();
        let mean: f64 = back.iter().zip(&bff.mass).map(|(x, m)| x * m).sum::<f64>() / total;
        assert!(mean.abs() < 1e-13);
    }

    #[test]
    fn harmonic_extension_reproduces_the_prescribed_angles() {
        let g = generators::bumpy_grid_disk(8, 8, 0.2).unwrap();
        let (mesh, geometry) = setup(&g);
        let mut bff = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap();
        let n = bff.boundary().len();
        // A rectangle: four right angles at the corners of the grid, straight elsewhere.
        let corners = [0u32, 8, 80, 72];
        let k: Vec<f64> = bff
            .boundary()
            .iter()
            .map(|v| if corners.contains(v) { PI / 2.0 } else { 0.0 })
            .collect();
        let result = bff.flatten_with_exterior_angles(&k, Extension::Harmonic).unwrap();
        let boundary = bff.boundary().to_vec();
        for i in 0..n {
            let (a, b, c) = (
                result.uv[boundary[(i + n - 1) % n] as usize],
                result.uv[boundary[i] as usize],
                result.uv[boundary[(i + 1) % n] as usize],
            );
            let turn = ((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]))
                .atan2((b[0] - a[0]) * (c[0] - b[0]) + (b[1] - a[1]) * (c[1] - b[1]));
            assert!((turn - k[i]).abs() < 1e-9, "vertex {i}: {turn} vs {}", k[i]);
        }
        assert_eq!(result.report.distortion.flipped_triangles, 0);
        assert!(result.report.closure.gap_after < 1e-12);
    }

    #[test]
    fn automatic_flattening_of_a_bumpy_disk_is_nearly_conformal() {
        let g = generators::bumpy_grid_disk(24, 24, 0.08).unwrap();
        let (mesh, geometry) = setup(&g);
        let mut bff = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap();
        let result = bff.flatten_automatic().unwrap();
        assert_eq!(result.extension, Extension::Conjugate);
        assert_eq!(result.report.distortion.flipped_triangles, 0);
        assert!(result.report.distortion.quasi_conformal.mean < 1.02);
        assert!(result.report.is_injective());
        assert!(result.solve_residual < 1e-10);
        // u = 0 on the boundary: boundary lengths are preserved up to the closing adjustment.
        let boundary = bff.boundary().to_vec();
        let n = boundary.len();
        for i in 0..n {
            let (a, b) = (
                result.uv[boundary[i] as usize],
                result.uv[boundary[(i + 1) % n] as usize],
            );
            let ratio = (b[0] - a[0]).hypot(b[1] - a[1]) / bff.boundary_edge_lengths()[i];
            assert!((ratio - 1.0).abs() < 0.05, "edge {i}: {ratio}");
        }
    }

    #[test]
    fn flatten_to_disk_gives_a_round_boundary() {
        let g = generators::bumpy_grid_disk(16, 16, 0.1).unwrap();
        let (mesh, geometry) = setup(&g);
        let mut bff = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap();
        let result = bff.flatten_to_disk(10).unwrap();
        let boundary: Vec<[f64; 2]> = bff.boundary().iter().map(|&v| result.uv[v as usize]).collect();
        let n = boundary.len() as f64;
        let centre = boundary
            .iter()
            .fold([0.0, 0.0], |c, p| [c[0] + p[0] / n, c[1] + p[1] / n]);
        let radii: Vec<f64> = boundary
            .iter()
            .map(|p| (p[0] - centre[0]).hypot(p[1] - centre[1]))
            .collect();
        let mean = radii.iter().sum::<f64>() / n;
        for r in &radii {
            assert!((r / mean - 1.0).abs() < 2e-2, "radius {r} vs mean {mean}");
        }
        assert_eq!(result.report.distortion.flipped_triangles, 0);
    }

    #[test]
    fn bad_input_is_rejected() {
        let g = generators::annulus(8, 2, 1.0, 2.0).unwrap();
        let (mesh, geometry) = setup(&g);
        let err = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap_err();
        assert!(matches!(err, BffError::NotADisk { boundary_loops: 2, .. }), "{err}");

        let g = generators::grid_disk(3, 3).unwrap();
        let (mesh, geometry) = setup(&g);
        let mut bff = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap();
        let n = bff.boundary().len();
        let err = bff
            .flatten_with_scale_factors(&vec![0.0; n + 1], Extension::Harmonic)
            .unwrap_err();
        assert!(matches!(err, BffError::LengthMismatch { .. }), "{err}");
        let mut u = vec![0.0; n];
        u[2] = f64::NAN;
        let err = bff
            .flatten_with_scale_factors(&u, Extension::Harmonic)
            .unwrap_err();
        assert!(matches!(err, BffError::NonFiniteInput { index: 2, .. }), "{err}");
        let err = bff
            .flatten_with_exterior_angles(&vec![0.5; n], Extension::Harmonic)
            .unwrap_err();
        assert!(matches!(err, BffError::ExteriorAnglesDoNotClose { .. }), "{err}");
    }
}
