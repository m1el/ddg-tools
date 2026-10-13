//! The static evaluation of a cut (SPEC §3.1–§3.7, §6): from a level-set function to the scale
//! factor, the gradient density on the cut, the speed of the cut, and the energy.

use ddg_core::operators::vertex_angle_sums;
use ddg_core::sparse::{Factorization, Parallelism, Symbolic};
use ddg_core::{Geometry, TriMesh};

use crate::curvature::{face_curvatures_from, feature_size_weight};
use crate::cut::Cut;
use crate::cut_mesh::{BaseGeometry, CutMesh, CutTriangles, EdgeKind};
use crate::error::CutError;
use crate::problem::{
    CurvatureSource, CutProblem, Energy, EvalOptions, PointwiseEnergy, ScaleWeight, UnitAreaScale,
};
use crate::recovery::FluxSolver;
use crate::system::LaplaceSystem;

/// The terms of the discrete energy (SPEC §6), on the unit-area surface:
///
/// ```text
/// E = Σ_r [ w_r J_r + β_r Area(S_r) ] + Σ_segments ½ (w_p + w_q) |pq|,     w = α_L + η
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyTerms {
    /// The distortion `J_r` of every region, without its weight `w_r`. A closed component that
    /// the cut does not touch contributes nothing (SPEC D9c).
    pub distortion: Vec<f64>,
    /// The area `Area(S_r)` of every region, without the closed components that the cut does
    /// not touch (SPEC D9c; their areas are in [`EvalReport::excluded_area`]). The two together
    /// sum to 1.
    pub area: Vec<f64>,
    /// The length `Σ |pq|` of the cut, counted once.
    pub length: f64,
    /// The penalty `Σ ½ (η_p + η_q) |pq|`, with `η` interpolated linearly to the cut points.
    pub penalty: f64,
    /// The energy `E = Σ_r (w_r J_r + β_r Area_r) + α_L · length + penalty`.
    pub total: f64,
}

/// Counters and residuals of one evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct EvalReport {
    /// Number of cut points (crossed edges).
    pub n_cut_points: usize,
    /// Number of cut segments (crossed faces).
    pub n_segments: usize,
    /// Number of crossings that were moved away from a vertex (SPEC §3.2, D9b).
    pub n_clamped_crossings: usize,
    /// The smallest distance of a crossing from the nearer end of its edge, as a fraction of
    /// the edge, before any was moved; `0.5` if there is no crossing.
    pub smallest_crossing_offset: f64,
    /// Number of cut points that are sliver vertices (SPEC §3.7).
    pub n_sliver_points: usize,
    /// The largest aspect ratio among the triangles of cut faces (1 if no face is cut).
    pub max_aspect_ratio: f64,
    /// Number of cut faces in which the two diagonals of the quad agree to a relative `10⁻⁹`,
    /// so that rounding decides which one splits it (SPEC §3.3).
    pub n_diagonal_near_ties: usize,
    /// The closed connected components of the input mesh that the cut does not touch, by
    /// component id ([`TriMesh::face_components`]). The Yamabe equation has no solution there;
    /// `u = p = 0` is reported and they are left out of the energy: of the distortion and of
    /// the area term (SPEC §3.4, D9c).
    pub uncut_closed_components: Vec<u32>,
    /// For every region, the area of its part in [`EvalReport::uncut_closed_components`].
    pub excluded_area: Vec<f64>,
    /// The residual `max_i |(L u)_i + g_i| / max_i |g_i|` of the Yamabe system over the
    /// interior vertices (`0` if `g = 0`).
    pub yamabe_residual: f64,
    /// The same for the adjoint system `L p = r`; `None` for the Dirichlet energy, which needs
    /// no adjoint solve.
    pub adjoint_residual: Option<f64>,
}

/// Everything one evaluation of a cut produces.
///
/// **Units.** All numbers are on the surface scaled to unit area ([`Evaluation::scale`]
/// converts). `log_scale`, `flux` and `adjoint_flux` do not depend on the unit.
///
/// **Indexing.** Per-vertex arrays are indexed by vertex of the *cut mesh*
/// ([`CutTriangles`]): the vertices of the input mesh keep their index, and every cut point
/// has one copy per side, so each entry belongs to exactly one region
/// ([`CutTriangles::vertex_region`]). The *boundary vertices* of a region are the copies of cut
/// points and the vertices on the boundary of the input mesh.
#[derive(Debug, Clone)]
pub struct Evaluation {
    /// The scaling to unit area that all numbers refer to.
    pub scale: UnitAreaScale,
    /// The cut: points on edges, segments across faces, the region of every input vertex.
    pub cut: Cut,
    /// The cut mesh in raw form. [`Evaluation::cut_mesh`] turns it into a mesh with geometry.
    pub triangles: CutTriangles,
    /// `u`: the log conformal scale factor (SPEC §3.4). Zero at boundary vertices.
    pub log_scale: Vec<f64>,
    /// `v`: the integrated outward normal derivative of `u` at boundary vertices (SPEC §3.5);
    /// zero elsewhere.
    pub flux: Vec<f64>,
    /// `ℓ`: half the length of the boundary edges of the region at a boundary vertex
    /// (SPEC §3.5); zero elsewhere.
    pub dual_length: Vec<f64>,
    /// `∂u/∂n`: the outward normal derivative of `u` at boundary vertices, per unit length of
    /// the unit-area surface; zero elsewhere. At the copies of the cut points it comes from the
    /// [`FluxRecovery`](crate::FluxRecovery) of the options (SPEC §3.5); at the vertices on the
    /// boundary of the input mesh it is `flux / dual_length`.
    pub normal_derivative: Vec<f64>,
    /// `p`: the adjoint state (SPEC §3.6). `None` for the Dirichlet energy, where `p = −2u`.
    pub adjoint: Option<Vec<f64>>,
    /// `q`: the integrated outward normal derivative of `p` at boundary vertices. `None` for
    /// the Dirichlet energy, where `q = −2v`.
    pub adjoint_flux: Option<Vec<f64>>,
    /// `∂p/∂n`, obtained from `q` exactly as `normal_derivative` from `flux`. `None` for the
    /// Dirichlet energy, where it is `−2 ∂u/∂n`.
    pub adjoint_normal_derivative: Option<Vec<f64>>,
    /// `G_r = w_r · G + β_r`: the gradient density of the region's energy at the copies of the
    /// cut points (SPEC §1, §3.6); zero at every other vertex. Moving the boundary of the
    /// region outward with speed `σ` changes its energy at the rate `∫ G_r σ ds`.
    pub gradient: Vec<f64>,
    /// `s` at every cut point before sliver smoothing (SPEC §3.7): the gradient density of the
    /// region with the lower id minus that of the region with the higher id. With two regions
    /// this is `G_0 − G_1`, the rate of change of the level-set function on the cut.
    pub speed_raw: Vec<f64>,
    /// `s` at every cut point after sliver smoothing.
    pub speed: Vec<f64>,
    /// The energy.
    pub energy: EnergyTerms,
    /// Counters and residuals.
    pub report: EvalReport,
}

impl Evaluation {
    /// The cut mesh as a validated mesh with geometry (in the caller's length unit) and maps to
    /// the input mesh (SPEC §7).
    ///
    /// # Errors
    /// See [`CutTriangles::to_mesh`].
    pub fn cut_mesh(&self) -> Result<CutMesh, CutError> {
        self.triangles.to_mesh(&self.scale)
    }
}

/// Maps level-set functions on one mesh to their [`Evaluation`]. It owns the symbolic
/// factorization of the mesh (one analysis serves every cut, SPEC §3.4) and keeps the numeric
/// factorization of the last evaluation for further solves with the same matrix.
#[derive(Debug)]
pub struct Evaluator {
    mesh: TriMesh,
    scale: UnitAreaScale,
    base: BaseGeometry,
    /// Gaussian curvature density per face, on the unit-area surface.
    gaussian: Vec<f64>,
    /// Lumped vertex areas of the unit-area surface.
    vertex_area: Vec<f64>,
    /// Sum of the corner angles at every vertex.
    angle_sums: Vec<f64>,
    vertex_component: Vec<u32>,
    component_closed: Vec<bool>,
    system: LaplaceSystem,
    /// The curvature-based weight computed last, with its `(r_eps, smoothing)`.
    feature_weight: Option<(f64, f64, Vec<f64>)>,
}

/// The pointwise energy density in use.
enum Density<'a> {
    Hencky,
    User(&'a PointwiseEnergy),
}

impl Density<'_> {
    #[inline]
    fn value(&self, u: f64) -> f64 {
        match self {
            Density::Hencky => u * u,
            Density::User(b) => b.value(u),
        }
    }

    #[inline]
    fn derivative(&self, u: f64) -> f64 {
        match self {
            Density::Hencky => 2.0 * u,
            Density::User(b) => b.derivative(u),
        }
    }
}

/// `max |residual| / max |rhs|` (0 if the right-hand side vanishes).
fn relative(residual: f64, rhs: f64) -> f64 {
    if rhs > 0.0 { residual / rhs } else { residual }
}

impl Evaluator {
    /// Prepares the evaluation of cuts on `mesh` with the intrinsic `geometry`: scales the
    /// surface to unit area, computes the curvature per face, and analyses the sparsity pattern
    /// of the linear systems. The mesh may have any number of components and any topology.
    ///
    /// # Errors
    /// [`CutError::Geometry`] if the geometry does not belong to the mesh,
    /// [`CutError::InvalidTotalArea`], or a resource error of the symbolic analysis.
    pub fn new(mesh: &TriMesh, geometry: &Geometry) -> Result<Evaluator, CutError> {
        geometry.check_mesh(mesh)?;
        let scale = UnitAreaScale::from_total_area(geometry.total_area())?;
        let base = BaseGeometry::new(mesh, geometry, &scale);
        let angle_sums = vertex_angle_sums(mesh, geometry)?;
        let gaussian = face_curvatures_from(mesh, geometry.corner_angles(), &angle_sums)
            .iter()
            .zip(&base.area)
            .map(|(omega, area)| omega / area)
            .collect();
        let mut vertex_area = vec![0.0; mesh.n_vertices()];
        for (face, area) in mesh.faces().iter().zip(&base.area) {
            for &v in face {
                vertex_area[v as usize] += area / 3.0;
            }
        }
        let vertex_component = mesh.vertex_components();
        let mut component_closed = vec![true; mesh.n_components()];
        for (v, &c) in vertex_component.iter().enumerate() {
            if mesh.is_boundary_vertex(v) {
                component_closed[c as usize] = false;
            }
        }
        let system = LaplaceSystem::new(mesh, |v| !mesh.is_boundary_vertex(v))?;
        Ok(Evaluator {
            mesh: mesh.clone(),
            scale,
            base,
            gaussian,
            vertex_area,
            angle_sums,
            vertex_component,
            component_closed,
            system,
            feature_weight: None,
        })
    }

    /// The mesh this evaluator was built for.
    #[inline]
    pub fn mesh(&self) -> &TriMesh {
        &self.mesh
    }

    /// The scaling of the input surface to unit area.
    #[inline]
    pub fn scale(&self) -> UnitAreaScale {
        self.scale
    }

    /// The Gaussian curvature density `k_f = Ω_f / A_f` of every face on the unit-area surface
    /// (SPEC §3.1): what [`CurvatureSource::Gaussian`] uses.
    #[inline]
    pub fn gaussian_curvature(&self) -> &[f64] {
        &self.gaussian
    }

    /// The lumped vertex areas `A_i` of the unit-area surface; they sum to 1.
    #[inline]
    pub fn vertex_areas(&self) -> &[f64] {
        &self.vertex_area
    }

    /// For every vertex of the input mesh, its index among the unknowns of the linear system
    /// (the interior vertices), or [`NONE`](ddg_core::topology::NONE).
    #[inline]
    pub fn unknown_of_vertex(&self) -> &[u32] {
        self.system.unknown_of_vertex()
    }

    /// For every unknown of the linear system, its vertex of the input mesh.
    #[inline]
    pub fn vertex_of_unknown(&self) -> &[u32] {
        self.system.vertex_of_unknown()
    }

    /// The symbolic analysis shared by all evaluations; `None` if the mesh has no interior
    /// vertex.
    #[inline]
    pub fn symbolic(&self) -> Option<&Symbolic> {
        self.system.symbolic()
    }

    /// The numeric factorization of the last successful [`Evaluator::evaluate`]: the cotan
    /// Laplacian of that cut mesh on the interior vertices of the input mesh (indexed by
    /// [`Evaluator::unknown_of_vertex`]), with identity rows for the vertices of
    /// [`EvalReport::uncut_closed_components`]. On a closed mesh this is also the matrix of the
    /// harmonic extension (SPEC §3.8). `None` before the first evaluation, after a failed one,
    /// or if the mesh has no interior vertex.
    #[inline]
    pub fn factorization(&self) -> Option<&Factorization> {
        self.system.factorization()
    }

    /// The linear system on the interior vertices, for the flow step.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn system(&self) -> &LaplaceSystem {
        &self.system
    }

    /// The curvature-based weight `a = r³` of the scaled Dirichlet energy for every vertex
    /// (paper Sec. 3.5.1, 5.6; SPEC §3.6), on the unit-area surface. The result is cached: a
    /// second call with the same parameters costs nothing.
    ///
    /// # Errors
    /// [`CutError::InvalidFeatureSize`], [`CutError::NonPositiveFeatureSize`], or a failure of
    /// the sparse solve.
    ///
    /// [formula: cut-feature-size]
    pub fn feature_size_weight(
        &mut self,
        r_eps: f64,
        smoothing: f64,
        par: Parallelism<'_>,
    ) -> Result<&[f64], CutError> {
        let cached = matches!(&self.feature_weight, Some((r, c, _)) if *r == r_eps && *c == smoothing);
        if !cached {
            let weight = feature_size_weight(
                &self.mesh,
                &self.base.cot,
                &self.vertex_area,
                &self.angle_sums,
                r_eps,
                smoothing,
                par,
            )?;
            self.feature_weight = Some((r_eps, smoothing, weight));
        }
        Ok(&self
            .feature_weight
            .as_ref()
            .expect("the weight was stored above")
            .2)
    }

    /// Evaluates the cut of the level-set function `phi` (one value per vertex; region 0 is
    /// where `φ ≤ 0`) for the energy described by `problem`: SPEC §3.2–§3.7 and the energy of
    /// §6.
    ///
    /// Only the signs and ratios of `phi` matter, so its length unit is irrelevant here.
    ///
    /// # Errors
    /// Every [`CutError`] that concerns the level-set function, the problem or the options;
    /// [`CutError::DegenerateCutTriangle`] if a crossing sits numerically on a vertex (only
    /// with a `min_crossing_offset` near zero); [`CutError::Sparse`] if the Laplacian of the
    /// cut mesh is numerically not positive definite.
    pub fn evaluate(
        &mut self,
        phi: &[f64],
        problem: &CutProblem,
        options: &EvalOptions<'_>,
    ) -> Result<Evaluation, CutError> {
        options.validate()?;
        let nv = self.mesh.n_vertices();
        problem.validate(nv, self.mesh.n_faces())?;
        let par = options.parallelism;
        let n_regions = problem.regions.len();

        // §3.2, §3.3: regions, crossings, cut mesh.
        let cut = Cut::from_level_set(&self.mesh, phi, options.min_crossing_offset)?;
        let tri = CutTriangles::build(&self.mesh, &self.base, &cut, options.sliver_aspect_ratio)?;
        let n = tri.n_vertices();

        // D9c: closed components without a cut are taken out of the system.
        let mut touched = vec![false; self.component_closed.len()];
        for point in &cut.points {
            let [a, _] = self.mesh.edge_vertices(point.edge as usize);
            touched[self.vertex_component[a] as usize] = true;
        }
        let excluded: Vec<bool> = self
            .component_closed
            .iter()
            .zip(&touched)
            .map(|(&closed, &touched)| closed && !touched)
            .collect();
        let pinned: Vec<bool> = self
            .vertex_component
            .iter()
            .map(|&c| excluded[c as usize])
            .collect();

        // §3.4: the source term g_i = Σ k_t A_t / 3 and the lumped vertex areas of the regions.
        let mut source = vec![0.0; n];
        let mut area_of_vertex = vec![0.0; n];
        let mut area = vec![0.0; n_regions];
        let mut excluded_area = vec![0.0; n_regions];
        for (t, face) in tri.faces.iter().enumerate() {
            let face_area = tri.area[t];
            let parent = tri.parent[t] as usize;
            if pinned[self.mesh.faces()[parent][0] as usize] {
                excluded_area[tri.region[t] as usize] += face_area;
                continue;
            }
            area[tri.region[t] as usize] += face_area;
            let density = match &problem.curvature {
                CurvatureSource::Gaussian => self.gaussian[parent],
                CurvatureSource::Uniform(c) => *c,
                CurvatureSource::PerFace(k) => k[parent],
            };
            for &v in face {
                source[v as usize] += density * face_area / 3.0;
                area_of_vertex[v as usize] += face_area / 3.0;
            }
        }

        // §3.4: L_II u_I = −g_I, u_B = 0.
        self.system.assemble(&self.mesh, &tri, &pinned);
        self.system.factor(par)?;
        let mut log_scale = vec![0.0; n];
        self.system.solve(|v| -source[v], &mut log_scale, par)?;

        // §3.5: v_i = (L u)_i + g_i on the boundary vertices; the same expression is the
        // residual of the system at the interior vertices.
        let mut product = vec![0.0; n];
        tri.apply_laplacian(None, &log_scale, &mut product);
        let mut flux = vec![0.0; n];
        let (mut worst, mut largest) = (0.0f64, 0.0f64);
        for v in 0..n {
            if self.system.is_unknown(v) {
                if !pinned[v] {
                    worst = worst.max((product[v] + source[v]).abs());
                    largest = largest.max(source[v].abs());
                }
            } else {
                flux[v] = normal_flux(product[v], source[v]);
            }
        }
        let yamabe_residual = relative(worst, largest);
        let dual_length = dual_lengths(&tri);
        let recovery =
            FluxSolver::new(&tri, &self.base, &dual_length, options.flux_recovery, par)?;
        let normal_derivative = recovery.recover(&flux, &dual_length, par)?;

        // §3.6: energy, adjoint, gradient density G at the copies of the cut points.
        let boundary: Vec<bool> = (0..n).map(|v| !self.system.is_unknown(v)).collect();
        let mut distortion = vec![0.0; n_regions];
        let mut density_on_cut = vec![0.0; n];
        let mut adjoint = None;
        let mut adjoint_flux = None;
        let mut adjoint_normal_derivative = None;
        let mut adjoint_residual = None;
        match &problem.energy {
            Energy::Dirichlet => {
                for v in 0..nv {
                    distortion[tri.vertex_region[v] as usize] +=
                        dirichlet_energy_term(log_scale[v], source[v]);
                }
                for v in nv..n {
                    density_on_cut[v] = gradient_dirichlet(normal_derivative[v]);
                }
            }
            Energy::ScaledDirichlet(weight) => {
                let weight: Vec<f64> = match weight {
                    ScaleWeight::PerVertex(a) => a.clone(),
                    ScaleWeight::FeatureSize { r_eps, smoothing } => {
                        self.feature_size_weight(*r_eps, *smoothing, par)?.to_vec()
                    }
                };
                let weight = tri.extend(&self.mesh, &cut, &weight);
                let mut weighted = vec![0.0; n];
                tri.apply_laplacian(Some(&weight), &log_scale, &mut weighted);
                for v in 0..nv {
                    distortion[tri.vertex_region[v] as usize] += log_scale[v] * weighted[v];
                }
                let rhs = |v: usize| adjoint_source_scaled(weighted[v], 0.0, 0.0);
                let mut p = vec![0.0; n];
                self.system.solve(rhs, &mut p, par)?;
                tri.apply_laplacian(None, &p, &mut product);
                let mut q = vec![0.0; n];
                let (mut worst, mut largest) = (0.0f64, 0.0f64);
                for v in 0..n {
                    if boundary[v] {
                        q[v] = product[v] - adjoint_source_scaled(weighted[v], weight[v], flux[v]);
                    } else if !pinned[v] {
                        worst = worst.max((product[v] - rhs(v)).abs());
                        largest = largest.max(rhs(v).abs());
                    }
                }
                let dp_dn = recovery.recover(&q, &dual_length, par)?;
                for v in nv..n {
                    density_on_cut[v] = gradient_scaled(weight[v], normal_derivative[v], dp_dn[v]);
                }
                adjoint_residual = Some(relative(worst, largest));
                adjoint = Some(p);
                adjoint_flux = Some(q);
                adjoint_normal_derivative = Some(dp_dn);
            }
            Energy::Hencky | Energy::Pointwise(_) => {
                let density = match &problem.energy {
                    Energy::Pointwise(b) => Density::User(b),
                    _ => Density::Hencky,
                };
                let mut rhs = vec![0.0; n];
                for v in 0..n {
                    let weight = area_of_vertex[v];
                    if weight == 0.0 {
                        continue;
                    }
                    let (b, db) = (
                        density.value(log_scale[v]),
                        density.derivative(log_scale[v]),
                    );
                    for (what, value) in [("b", b), ("b'", db)] {
                        if !value.is_finite() {
                            return Err(CutError::NonFiniteEnergyDensity {
                                vertex: v as u32,
                                log_scale: log_scale[v],
                                what,
                                value,
                            });
                        }
                    }
                    distortion[tri.vertex_region[v] as usize] += weight * b;
                    rhs[v] = adjoint_source_pointwise(weight, db);
                }
                let mut p = vec![0.0; n];
                self.system.solve(|v| rhs[v], &mut p, par)?;
                tri.apply_laplacian(None, &p, &mut product);
                let mut q = vec![0.0; n];
                let (mut worst, mut largest) = (0.0f64, 0.0f64);
                for v in 0..n {
                    if boundary[v] {
                        q[v] = product[v] - rhs[v];
                    } else if !pinned[v] {
                        worst = worst.max((product[v] - rhs[v]).abs());
                        largest = largest.max(rhs[v].abs());
                    }
                }
                let dp_dn = recovery.recover(&q, &dual_length, par)?;
                let at_zero = density.value(0.0);
                for v in nv..n {
                    density_on_cut[v] = gradient_pointwise(at_zero, normal_derivative[v], dp_dn[v]);
                }
                adjoint_residual = Some(relative(worst, largest));
                adjoint = Some(p);
                adjoint_flux = Some(q);
                adjoint_normal_derivative = Some(dp_dn);
            }
        }

        // Region totals G_r = w_r G + β_r at the copies of the cut points, and §3.7.
        let mut gradient = vec![0.0; n];
        for v in nv..n {
            let weights = &problem.regions[tri.vertex_region[v] as usize];
            gradient[v] = weights.distortion * density_on_cut[v] + weights.area;
        }
        let speed_raw: Vec<f64> = (0..cut.points.len())
            .map(|k| {
                let [lower, higher] = tri.point_vertices_by_region(k);
                cut_speed(gradient[lower], gradient[higher])
            })
            .collect();
        let speed = smooth_slivers(
            &cut,
            &tri.segment_length,
            &tri.sliver_point,
            &speed_raw,
            options.smoothing_passes,
        );

        // §6: the energy.
        let length: f64 = tri.segment_length.iter().sum();
        let penalty = match &problem.penalty {
            Some(eta) => {
                let eta = cut.interpolate(&self.mesh, eta);
                cut.segments
                    .iter()
                    .zip(&tri.segment_length)
                    .map(|(segment, length)| {
                        let [p, q] = segment.points;
                        0.5 * (eta[p as usize] + eta[q as usize]) * length
                    })
                    .sum()
            }
            None => 0.0,
        };
        let mut total = problem.length_weight * length + penalty;
        for (r, weights) in problem.regions.iter().enumerate() {
            total += weights.distortion * distortion[r] + weights.area * area[r];
        }

        let report = EvalReport {
            n_cut_points: cut.points.len(),
            n_segments: cut.segments.len(),
            n_clamped_crossings: cut.n_clamped,
            smallest_crossing_offset: cut.smallest_offset,
            n_sliver_points: tri.sliver_point.iter().filter(|&&s| s).count(),
            max_aspect_ratio: tri.max_aspect_ratio,
            n_diagonal_near_ties: tri.n_diagonal_near_ties,
            uncut_closed_components: (0..excluded.len() as u32)
                .filter(|&c| excluded[c as usize])
                .collect(),
            excluded_area,
            yamabe_residual,
            adjoint_residual,
        };
        Ok(Evaluation {
            scale: self.scale,
            cut,
            triangles: tri,
            log_scale,
            flux,
            dual_length,
            normal_derivative,
            adjoint,
            adjoint_flux,
            adjoint_normal_derivative,
            gradient,
            speed_raw,
            speed,
            energy: EnergyTerms {
                distortion,
                area,
                length,
                penalty,
                total,
            },
            report,
        })
    }
}

/// The integrated outward normal derivative of `u` at a boundary vertex `i` of a region
/// (SPEC §3.5, D6): `v_i = (L u)_i + g_i`, from the weak form of `Δu = −K` tested with the hat
/// function of `i`. Summed over a stretch of the boundary it approximates `∫ ∂u/∂n ds` over
/// that stretch; [`FluxRecovery`](crate::FluxRecovery) turns it into `∂u/∂n` per vertex.
///
/// [formula: cut-normal-derivative]
#[inline]
pub fn normal_flux(laplacian_u: f64, source: f64) -> f64 {
    laplacian_u + source
}

/// The dual length `ℓ_i` of every vertex of the cut mesh: half the sum of the lengths of the
/// boundary edges of its region at the vertex (cut segments and boundary edges of the input
/// mesh, SPEC §3.5, D9d); zero at interior vertices.
///
/// [formula: cut-normal-derivative]
pub fn dual_lengths(tri: &CutTriangles) -> Vec<f64> {
    let mut dual = vec![0.0; tri.n_vertices()];
    for ((face, len2), kinds) in tri.faces.iter().zip(&tri.len2).zip(&tri.edge_kind) {
        for k in 0..3 {
            if kinds[k] != EdgeKind::Interior {
                let half = 0.5 * len2[k].sqrt();
                dual[face[k] as usize] += half;
                dual[face[(k + 1) % 3] as usize] += half;
            }
        }
    }
    dual
}

/// The contribution `−u_i g_i` of an interior vertex to the Dirichlet energy
/// `J = uᵀ L u = −u_I · g_I` (SPEC §3.6).
///
/// [formula: cut-dirichlet-energy]
#[inline]
pub fn dirichlet_energy_term(log_scale: f64, source: f64) -> f64 {
    -log_scale * source
}

/// The gradient density of the Dirichlet energy, `G = (∂u/∂n)²` (paper Eq. 13).
///
/// [formula: cut-gradient-dirichlet]
#[inline]
pub fn gradient_dirichlet(du_dn: f64) -> f64 {
    du_dn * du_dn
}

/// The source `r_i` of the adjoint system of the scaled Dirichlet energy (SPEC §3.6, D6):
/// `r_i = −2 [ (L_a u)_i − a_i v_i ]`. At an interior vertex pass `a_i v_i = 0`; there
/// `L_II p_I = r_I`. At a boundary vertex `q_i = (L p)_i − r_i`.
///
/// [formula: cut-adjoint-scaled]
#[inline]
pub fn adjoint_source_scaled(weighted_laplacian_u: f64, weight: f64, flux: f64) -> f64 {
    -2.0 * (weighted_laplacian_u - weight * flux)
}

/// The gradient density of the scaled Dirichlet energy,
/// `G = −a (∂u/∂n)² − (∂u/∂n)(∂p/∂n)` (paper Eq. 18, Table 1).
///
/// [formula: cut-gradient-scaled]
#[inline]
pub fn gradient_scaled(weight: f64, du_dn: f64, dp_dn: f64) -> f64 {
    -weight * du_dn * du_dn - du_dn * dp_dn
}

/// The source `r_i = −A_i b′(u_i)` of the adjoint system of a pointwise energy (SPEC §3.6),
/// with `A_i` the lumped area of the vertex within its region. At an interior vertex
/// `L_II p_I = r_I`; at a boundary vertex `q_i = (L p)_i − r_i`.
///
/// [formula: cut-adjoint-pointwise]
#[inline]
pub fn adjoint_source_pointwise(vertex_area: f64, db: f64) -> f64 {
    -vertex_area * db
}

/// The gradient density of a pointwise energy on the boundary, where `u = 0`:
/// `G = b(0) − (∂u/∂n)(∂p/∂n)` (paper Table 1).
///
/// [formula: cut-gradient-pointwise]
#[inline]
pub fn gradient_pointwise(b_at_zero: f64, du_dn: f64, dp_dn: f64) -> f64 {
    b_at_zero - du_dn * dp_dn
}

/// The speed of the cut at a cut point (SPEC §3.7, D2): the gradient density of the region with
/// the lower id minus that of the region with the higher id, `s = G_0 − G_1` for two regions.
/// The level-set function of the lower region changes at this rate, so a region whose own
/// gradient density is the larger one shrinks.
///
/// [formula: cut-speed]
#[inline]
pub fn cut_speed(gradient_lower: f64, gradient_higher: f64) -> f64 {
    gradient_lower - gradient_higher
}

/// Sliver smoothing (SPEC §3.7, D7): in each of `passes` rounds every sliver point that has two
/// neighbours along the cut takes, simultaneously with all others, the arclength-weighted mean
/// `(ℓ₂ s₁ + ℓ₁ s₂)/(ℓ₁ + ℓ₂)` of their values, `ℓ_k` being the length of the segment to
/// neighbour `k`. A cut point on the boundary of the mesh (one neighbour) keeps its value.
///
/// [formula: cut-sliver-smoothing]
pub fn smooth_slivers(
    cut: &Cut,
    segment_length: &[f64],
    sliver_point: &[bool],
    speed: &[f64],
    passes: usize,
) -> Vec<f64> {
    let mut out = speed.to_vec();
    if passes == 0 || !sliver_point.iter().any(|&s| s) {
        return out;
    }
    // The (at most two) neighbours of every point along the cut, with the segment lengths.
    let mut neighbours: Vec<[(u32, f64); 2]> = vec![[(u32::MAX, 0.0); 2]; cut.points.len()];
    let mut count = vec![0u8; cut.points.len()];
    for (segment, &length) in cut.segments.iter().zip(segment_length) {
        let [p, q] = segment.points;
        for (from, to) in [(p, q), (q, p)] {
            let slot = &mut count[from as usize];
            if *slot < 2 {
                neighbours[from as usize][*slot as usize] = (to, length);
            }
            *slot += 1;
        }
    }
    let mut previous = out.clone();
    for _ in 0..passes {
        previous.copy_from_slice(&out);
        for k in 0..out.len() {
            if sliver_point[k] && count[k] == 2 {
                let [(n1, l1), (n2, l2)] = neighbours[k];
                out[k] = (l2 * previous[n1 as usize] + l1 * previous[n2 as usize]) / (l1 + l2);
            }
        }
    }
    out
}
