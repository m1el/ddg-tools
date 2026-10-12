//! Flattening with cone singularities (paper Sec. 6.4, `BFF::flattenWithCones`).
//!
//! Three steps:
//!
//! 1. Solve the Cherrier problem for the log scale factor `u` on the uncut surface:
//!    `(L u)_i = µ_i − Ω_i` at interior vertices, `u = 0` on the boundary, where `µ` is zero
//!    except at the cones ([`solve_scale_factors`]). The metric `e^{2u} g` is then flat away
//!    from the cones. A precomputed `u` can be passed in instead.
//! 2. Cut the surface into disks through all cones ([`crate::plan`], [`crate::cut`]).
//! 3. For every disk, read off the boundary data of its flattening from `u` — target lengths
//!    `ℓ* = e^{(u_i + u_j)/2} ℓ` and exterior angles `k̃ = k − h`, where the Neumann data `h`
//!    is evaluated with the Laplacian of the *cut* surface — build the closed boundary curve
//!    with one length per pair of cut sides, and extend it harmonically.
//!
//! Because the two sides of a cut edge get the same length, and because the exterior angles
//! at the copies of a vertex add up to what its cone angle demands (paper Proposition 2), the
//! layouts on the two sides of every cut differ by a rigid motion: the map is seamless.
//!
//! Ported from `Bff.cpp` of the reference implementation (MIT, © 2017 Rohan Sawhney and
//! Keenan Crane); see the crate documentation for the deviations.

use std::f64::consts::PI;
use std::time::Instant;

use ddg_core::geometry::Geometry;
pub use ddg_core::io::Cone;
use ddg_core::topology::{NONE, TriMesh};

use crate::chart::{
    ChartReport, PoissonSystem, Surface, check_finite_layout, halfedge_lengths, layout_quality,
};
use crate::curve::{CurveOptions, LoopData, best_fit_curves};
use crate::cut::{ComponentTopology, CutMesh, component_topology, cut_mesh};
use crate::disk::curvatures;
use crate::error::{BffError, check_values};
use crate::plan::{CutPlan, CutPlanInput, check_disk_charts, plan_and_cut};

/// What to do when the cone angles on a closed connected component do not sum to `2πχ`, the
/// condition under which `L u = µ − Ω` has a solution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GaussBonnetPolicy {
    /// Return [`BffError::GaussBonnetViolated`] if the discrepancy exceeds the tolerance.
    /// The default.
    #[default]
    Reject,
    /// Change the cone angles: each cone of the component receives a share of the discrepancy
    /// in proportion to the absolute value of its angle, so that the sum becomes `2πχ` and
    /// the curvature stays concentrated at the cones.
    DistributeOverCones,
    /// Keep the cone angles and spread the discrepancy uniformly over all vertices of the
    /// component (each gets the same small curvature). This is the least-squares solution of
    /// the inconsistent system, and what the reference implementation computes when it is
    /// handed such angles, because of the diagonal shift in its Laplacian.
    ///
    /// With a non-zero [`ScaleFactorOptions::diagonal_shift`] the shifted system is solved as
    /// it stands, exactly like the reference: the discrepancy `d` then also shows up as a
    /// constant `d / (shift · n_vertices)` in `u`, which scales the whole layout.
    SpreadOverVertices,
}

/// Options of [`solve_scale_factors`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScaleFactorOptions {
    /// What to do with cone angles that violate Gauss–Bonnet on a closed component.
    pub gauss_bonnet: GaussBonnetPolicy,
    /// Absolute tolerance on `|Σµ − 2πχ|` for [`GaussBonnetPolicy::Reject`]. A discrepancy
    /// within the tolerance is spread uniformly over the vertices of the component.
    pub gauss_bonnet_tolerance: f64,
    /// A constant added to the diagonal of the Laplacian, see
    /// [`crate::disk::DiskOptions::diagonal_shift`]. With `0.0` (the default) the system is
    /// solved exactly: on a closed component the right-hand side is made compatible, one
    /// vertex is pinned, and the solution is shifted to an area-weighted mean of zero. With a
    /// positive shift the regularized system `(L + shift · I) u = µ − Ω` is solved as it
    /// stands, which is what the reference implementation does; `u` is then not normalized.
    pub diagonal_shift: f64,
}

impl Default for ScaleFactorOptions {
    fn default() -> Self {
        ScaleFactorOptions {
            gauss_bonnet: GaussBonnetPolicy::Reject,
            gauss_bonnet_tolerance: 1e-8,
            diagonal_shift: 0.0,
        }
    }
}

/// How the cone angles of one closed connected component relate to Gauss–Bonnet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosedComponentReport {
    /// Topology of the component.
    pub topology: ComponentTopology,
    /// Number of cones on it.
    pub n_cones: usize,
    /// Sum of the cone angles as given.
    pub cone_sum: f64,
    /// `2πχ`.
    pub expected: f64,
    /// `cone_sum − expected`.
    pub discrepancy: f64,
    /// What was done about the discrepancy.
    pub applied: GaussBonnetPolicy,
    /// The constant that was added to the pinned solution on this component to give `u` an
    /// area-weighted mean of zero. It only sets the overall size of the layout. Zero when a
    /// diagonal shift is used (nothing is added then).
    pub additive_constant: f64,
}

/// The solution of the Cherrier problem.
#[derive(Debug, Clone, PartialEq)]
pub struct ScaleFactors {
    /// The log scale factor `u` per vertex: zero on the boundary; on a closed component its
    /// area-weighted mean is zero (unless a diagonal shift is used).
    pub u: Vec<f64>,
    /// The cone angle `µ` per vertex that was used (zero away from the cones). Differs from
    /// the input only under [`GaussBonnetPolicy::DistributeOverCones`].
    pub cone_angles: Vec<f64>,
    /// One entry per closed connected component.
    pub closed_components: Vec<ClosedComponentReport>,
    /// Number of unknowns of the linear system.
    pub n_unknowns: usize,
    /// `max_i |(L u)_i + shift · u_i − (µ_i − Ω_i − s_i)|` over the interior vertices, where
    /// `s_i` is the uniformly spread discrepancy of the component (zero with a diagonal shift
    /// and for compatible cone angles).
    pub max_residual: f64,
    /// Relative residual of the linear solve.
    pub relative_residual: f64,
}

/// Validates a cone list and returns the cone angle per vertex and the cone flags.
fn cone_angles_per_vertex(
    mesh: &TriMesh,
    cones: &[Cone],
) -> Result<(Vec<f64>, Vec<u32>), BffError> {
    let n = mesh.n_vertices();
    let mut angle = vec![0.0; n];
    let mut position_of = vec![NONE; n];
    for (position, cone) in cones.iter().enumerate() {
        let v = cone.vertex as usize;
        if v >= n {
            return Err(BffError::ConeVertexOutOfRange {
                position,
                vertex: cone.vertex,
                n_vertices: n,
            });
        }
        if position_of[v] != NONE {
            return Err(BffError::DuplicateCone {
                vertex: cone.vertex,
                first: position_of[v] as usize,
                second: position,
            });
        }
        if mesh.is_boundary_vertex(v) {
            return Err(BffError::ConeOnBoundary {
                position,
                vertex: cone.vertex,
            });
        }
        if !cone.angle.is_finite() || cone.angle >= 2.0 * PI {
            return Err(BffError::InvalidConeAngle {
                position,
                vertex: cone.vertex,
                angle: cone.angle,
            });
        }
        position_of[v] = position as u32;
        angle[v] = cone.angle;
    }
    Ok((angle, position_of))
}

/// Solves `(L u)_i = µ_i − Ω_i` at the interior vertices with `u = 0` on the boundary, for
/// the cone angles `µ` given by `cones` (zero elsewhere).
///
/// On a connected component with boundary the solution is unique. On a closed component a
/// solution exists only if `Σ µ = 2πχ` (Gauss–Bonnet), which is checked and handled as
/// `options.gauss_bonnet` says; the solution is then unique up to a constant, which is fixed
/// by `Σ m_i u_i = 0` with the lumped masses `m` (the constant of least area distortion).
///
/// # Errors
/// [`BffError::ConeVertexOutOfRange`], [`BffError::DuplicateCone`],
/// [`BffError::ConeOnBoundary`], [`BffError::InvalidConeAngle`],
/// [`BffError::GaussBonnetViolated`], [`BffError::NoConesToAdjust`], [`BffError::Geometry`],
/// [`BffError::Solve`].
///
/// [formula: cone-scale-factor-equation]
/// [formula: closed-mesh-cone-sum]
/// [formula: closed-mesh-optimal-constant]
pub fn solve_scale_factors(
    mesh: &TriMesh,
    geometry: &Geometry,
    cones: &[Cone],
    options: &ScaleFactorOptions,
) -> Result<ScaleFactors, BffError> {
    geometry.check_mesh(mesh)?;
    let (mut mu, cone_position) = cone_angles_per_vertex(mesh, cones)?;
    let lengths = halfedge_lengths(mesh, geometry);
    let surface = Surface {
        mesh,
        cotans: geometry.corner_cotans(),
        angles: geometry.corner_angles(),
        areas: geometry.face_areas(),
        halfedge_lengths: &lengths,
    };
    let n = mesh.n_vertices();
    let is_boundary: Vec<bool> = (0..n).map(|v| mesh.is_boundary_vertex(v)).collect();
    let omega = curvatures(&surface.angle_sums(), &is_boundary);
    let component = mesh.vertex_components();
    let topology = component_topology(mesh);

    // Gauss-Bonnet on the closed components.
    let mut closed_index = vec![NONE; topology.len()];
    let mut closed: Vec<ClosedComponentReport> = Vec::new();
    for (c, t) in topology.iter().enumerate().filter(|(_, t)| t.is_closed()) {
        closed_index[c] = closed.len() as u32;
        let expected = 2.0 * PI * t.euler_characteristic as f64;
        closed.push(ClosedComponentReport {
            topology: *t,
            n_cones: 0,
            cone_sum: 0.0,
            expected,
            discrepancy: 0.0,
            applied: options.gauss_bonnet,
            additive_constant: 0.0,
        });
    }
    let mut absolute_sum = vec![0.0; closed.len()];
    for v in (0..n).filter(|&v| cone_position[v] != NONE) {
        let k = closed_index[component[v] as usize];
        if k != NONE {
            closed[k as usize].n_cones += 1;
            closed[k as usize].cone_sum += mu[v];
            absolute_sum[k as usize] += mu[v].abs();
        }
    }
    for (k, report) in closed.iter_mut().enumerate() {
        report.discrepancy = report.cone_sum - report.expected;
        let violated = report.discrepancy.abs() > options.gauss_bonnet_tolerance;
        match options.gauss_bonnet {
            GaussBonnetPolicy::Reject if violated => {
                return Err(BffError::GaussBonnetViolated {
                    component: closed_index
                        .iter()
                        .position(|&x| x == k as u32)
                        .unwrap_or(0) as u32,
                    first_face: report.topology.first_face,
                    euler_characteristic: report.topology.euler_characteristic,
                    n_cones: report.n_cones,
                    cone_sum: report.cone_sum,
                    expected: report.expected,
                    tolerance: options.gauss_bonnet_tolerance,
                });
            }
            GaussBonnetPolicy::DistributeOverCones if absolute_sum[k] == 0.0 && violated => {
                return Err(BffError::NoConesToAdjust {
                    component: closed_index
                        .iter()
                        .position(|&x| x == k as u32)
                        .unwrap_or(0) as u32,
                    discrepancy: report.discrepancy,
                });
            }
            _ => {}
        }
    }
    if options.gauss_bonnet == GaussBonnetPolicy::DistributeOverCones {
        for v in (0..n).filter(|&v| cone_position[v] != NONE) {
            let k = closed_index[component[v] as usize];
            if k != NONE && absolute_sum[k as usize] > 0.0 {
                mu[v] -= closed[k as usize].discrepancy * mu[v].abs() / absolute_sum[k as usize];
            }
        }
    }

    // Right-hand side µ − Ω; on closed components made exactly compatible by removing its
    // mean (the uniformly spread remainder).
    let mut rhs: Vec<f64> = (0..n)
        .map(|v| {
            if is_boundary[v] {
                0.0
            } else {
                mu[v] - omega[v]
            }
        })
        .collect();
    let exact = options.diagonal_shift == 0.0;
    if exact {
        let mut sums = vec![0.0; closed.len()];
        for v in 0..n {
            let k = closed_index[component[v] as usize];
            if k != NONE {
                sums[k as usize] += rhs[v];
            }
        }
        for v in 0..n {
            let k = closed_index[component[v] as usize];
            if k != NONE {
                rhs[v] -= sums[k as usize] / closed[k as usize].topology.n_vertices as f64;
            }
        }
    }

    // Unknowns: interior vertices; without a diagonal shift one vertex per closed component
    // is pinned to make the matrix definite.
    let mut fixed = is_boundary.clone();
    if exact {
        let mut pinned = vec![false; closed.len()];
        for v in 0..n {
            let k = closed_index[component[v] as usize];
            if k != NONE && !pinned[k as usize] {
                pinned[k as usize] = true;
                fixed[v] = true;
            }
        }
    }
    let system = PoissonSystem::new(
        "the Laplacian of the surface (scale factors)",
        &surface,
        &fixed,
        options.diagonal_shift,
    )?;
    let mut u: Vec<f64> = (0..n)
        .map(|v| if fixed[v] { 0.0 } else { rhs[v] })
        .collect();
    let relative_residual = system.solve(&surface, &mut u, 1)?;

    // Area-weighted mean zero on closed components.
    let mass = surface.masses();
    let mut weighted = vec![0.0; closed.len()];
    let mut total = vec![0.0; closed.len()];
    for v in 0..n {
        let k = closed_index[component[v] as usize];
        if k != NONE {
            weighted[k as usize] += mass[v] * u[v];
            total[k as usize] += mass[v];
        }
    }
    for (k, report) in closed.iter_mut().enumerate() {
        if exact {
            report.additive_constant = -weighted[k] / total[k];
        }
    }
    for v in 0..n {
        let k = closed_index[component[v] as usize];
        if k != NONE {
            u[v] += closed[k as usize].additive_constant;
        }
    }

    let mut product = vec![0.0; n];
    surface.apply_laplacian(&u, &mut product);
    let max_residual = (0..n)
        .filter(|&v| !is_boundary[v])
        .map(|v| (product[v] + options.diagonal_shift * u[v] - rhs[v]).abs())
        .fold(0.0, f64::max);
    if !max_residual.is_finite() {
        return Err(BffError::NonFiniteResult {
            what: "the scale factor",
            chart: None,
        });
    }
    Ok(ScaleFactors {
        u,
        cone_angles: mu,
        closed_components: closed,
        n_unknowns: system.n_free(),
        max_residual,
        relative_residual,
    })
}

/// Where the cuts come from.
#[derive(Debug, Clone, Copy)]
pub enum CutSource<'a> {
    /// Let [`crate::plan::plan_cut`] choose them.
    Plan {
        /// Per edge: whether it may be cut. `None` allows every edge.
        may_cut: Option<&'a [bool]>,
        /// Edges that must be cut (seams). They may split the surface into several charts.
        must_cut: &'a [u32],
    },
    /// Cut exactly these edges. The result must consist of disks with every cone on a
    /// boundary; otherwise an error is returned.
    Edges(&'a [u32]),
}

/// Input of [`flatten_with_cones`].
#[derive(Debug, Clone, Copy)]
pub struct ConeFlattenInput<'a> {
    /// The surface.
    pub mesh: &'a TriMesh,
    /// Its metric.
    pub geometry: &'a Geometry,
    /// The cones: interior vertices with their cone angles `µ = 2π − (angle sum in the
    /// plane)`. With `scale_factors` given, only the vertices are used (to route the cut and
    /// for the report); the angles are then implied by the scale factors.
    pub cones: &'a [Cone],
    /// A precomputed log scale factor per vertex, or `None` to solve for it from the cones.
    pub scale_factors: Option<&'a [f64]>,
    /// The cuts.
    pub cut: CutSource<'a>,
    /// Optional anchors, as corner ids `3 f + k` of the input mesh, at most one per chart: the
    /// boundary curve of the chart that contains face `f` starts at the copy of vertex
    /// `faces[f][k]` that this corner belongs to, which must lie on the cut or on the
    /// boundary. That vertex is placed at the origin and the boundary edge leaving it (chart
    /// on the left) points along `+x`. A chart without anchor starts at the first vertex of
    /// its boundary loop in the cut mesh. For exact solves the choice changes the layout of a
    /// chart only by a rigid motion.
    pub anchors: &'a [u32],
}

/// Options of [`flatten_with_cones`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConeFlattenOptions {
    /// Options of the scale factor solve (ignored if scale factors are passed in), including
    /// the diagonal shift, which is also applied to the Laplacian of the cut surface.
    pub scale_factors: ScaleFactorOptions,
    /// Scale factors passed in by the caller must be flat away from the cut: the curvature
    /// `Ω_i + (L u)_i` they imply at a vertex that is interior to its chart may not exceed
    /// this in absolute value.
    pub flatness_tolerance: f64,
    /// Largest fraction by which a boundary length may be changed to close a boundary curve;
    /// see [`crate::disk::DiskOptions::max_length_adjustment`].
    pub max_length_adjustment: f64,
    /// A cut edge between two different charts has one length shared by both, which couples
    /// the closing problems of the two charts. Groups of more than this many coupled charts
    /// are closed chart by chart instead; the two sides of such edges then differ in length
    /// by the closing adjustment (reported in [`SeamReport`]).
    pub max_coupled_charts: usize,
}

impl Default for ConeFlattenOptions {
    fn default() -> Self {
        ConeFlattenOptions {
            scale_factors: ScaleFactorOptions::default(),
            flatness_tolerance: 1e-6,
            max_length_adjustment: 0.5,
            max_coupled_charts: 512,
        }
    }
}

/// What happened at one cone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConeReport {
    /// The cone vertex.
    pub vertex: u32,
    /// The cone angle that was asked for (after any Gauss–Bonnet adjustment).
    pub angle: f64,
    /// The cone angle implied by the scale factors, `Ω + (L u)` at the vertex. Equals `angle`
    /// up to the solver residual when the scale factors were solved for here.
    pub implied_angle: f64,
    /// The cone angle of the layout: `2π` minus the sum of the signed corner angles of the
    /// layout at all copies of the vertex.
    pub realized_angle: f64,
    /// Number of copies of the vertex in the cut surface: the number of cut edges at it.
    pub n_wedges: usize,
    /// The largest angle that the layout spans around one copy. More than `2π` means that
    /// the triangles around that copy overlap, which a cone of negative angle at the end of
    /// a cut always produces.
    pub max_wedge_angle: f64,
}

/// Agreement of the layout across the cut edges.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeamReport {
    /// Number of cut edges.
    pub n_cut_edges: usize,
    /// Number of cut edges whose two sides lie in different charts.
    pub n_between_charts: usize,
    /// `max |ℓ_a − ℓ_b| / max(ℓ_a, ℓ_b)` over the cut edges, with `ℓ_a`, `ℓ_b` the lengths of
    /// the two sides in the layout. Rounding level when lengths are shared.
    pub max_length_mismatch: f64,
    /// The cut edge with the largest mismatch.
    pub worst_edge: Option<u32>,
}

/// Wall-clock seconds spent in the stages of [`flatten_with_cones`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Timings {
    /// Solving for the scale factors (zero if they were passed in).
    pub scale_factors: f64,
    /// Planning the cut and cutting.
    pub cut: f64,
    /// Boundary data and boundary curves.
    pub curves: f64,
    /// Factorizing the Laplacian of the cut surface and extending the curves.
    pub extension: f64,
    /// Computing the reports.
    pub report: f64,
    /// Everything.
    pub total: f64,
}

/// A flattening with cones.
#[derive(Debug, Clone, PartialEq)]
pub struct ConeFlattening {
    /// The position in the plane of every triangle corner: entry `3 f + k` belongs to the
    /// corner of face `f` at vertex `faces[f][k]` of the input mesh. A vertex on the cut has
    /// different positions in different corners.
    ///
    /// Each chart is placed with the first vertex of its boundary loop (see
    /// `TriMesh::boundary_loops` of `cut.mesh`, or its anchor) at the origin and its first
    /// boundary edge along `+x`. Charts are not packed, recentred or rescaled: lengths in the plane are
    /// about `e^u` times the lengths on the surface in every chart alike, so charts lie on
    /// top of each other and their relative size is meaningful.
    pub corner_uvs: Vec<[f64; 2]>,
    /// The position of every vertex of the cut mesh (`corner_uvs` by vertex).
    pub vertex_uvs: Vec<[f64; 2]>,
    /// The cut surface with its maps to the input mesh; [`CutMesh::cut_edges`] lists the cut
    /// edges.
    pub cut: CutMesh,
    /// The chart of every face. Charts are the connected components of the cut surface,
    /// numbered by their lowest face.
    pub face_chart: Vec<u32>,
    /// One report per chart.
    pub charts: Vec<ChartReport>,
    /// The planned cut, if the cut was planned here.
    pub plan: Option<CutPlan>,
    /// The log scale factor per vertex of the input mesh that the layout was built from.
    pub scale_factors: Vec<f64>,
    /// The scale factor solve, if it happened here.
    pub scale_factor_solve: Option<ScaleFactors>,
    /// The largest curvature `|Ω_i + (L u)_i|` implied by the scale factors at a vertex that
    /// is interior to its chart; zero for an exactly flat cone metric.
    pub max_interior_curvature: f64,
    /// Relative residual of the harmonic extension.
    pub extension_residual: f64,
    /// One report per cone, in the order of the input.
    pub cones: Vec<ConeReport>,
    /// Agreement across the cut.
    pub seams: SeamReport,
    /// Where the time went.
    pub timings: Timings,
}

impl ConeFlattening {
    /// The edges of the input mesh that were cut, ascending.
    pub fn cut_edges(&self) -> &[u32] {
        &self.cut.cut_edges
    }

    /// Number of charts.
    pub fn n_charts(&self) -> usize {
        self.charts.len()
    }
}

/// Flattens a surface with cone singularities. See the [module documentation](self).
///
/// # Errors
/// Those of [`solve_scale_factors`], [`crate::plan::plan_cut`] and [`crate::cut::cut_mesh`];
/// [`BffError::ChartNotDisk`] and [`BffError::ConeNotOnCut`] for explicit cuts that do not
/// work; [`BffError::AnchorOutOfRange`], [`BffError::AnchorNotOnBoundary`] and
/// [`BffError::DuplicateAnchor`] for bad anchors; [`BffError::ScaleFactorsNotFlat`] for scale factors that are not flat away from the
/// cut; [`BffError::Solve`]; [`BffError::NonFiniteResult`] instead of returning NaN.
pub fn flatten_with_cones(
    input: &ConeFlattenInput<'_>,
    options: &ConeFlattenOptions,
) -> Result<ConeFlattening, BffError> {
    let start = Instant::now();
    let mut timings = Timings::default();
    let (mesh, geometry) = (input.mesh, input.geometry);
    geometry.check_mesh(mesh)?;
    let n = mesh.n_vertices();
    let shift = options.scale_factors.diagonal_shift;
    if !shift.is_finite() || shift < 0.0 {
        return Err(BffError::NonFiniteInput {
            what: "diagonal_shift (must be finite and non-negative)",
            index: 0,
            value: shift,
        });
    }

    // Step 1: scale factors.
    let (cone_angle, cone_position) = cone_angles_per_vertex(mesh, input.cones)?;
    let stage = Instant::now();
    let (u, solve) = match input.scale_factors {
        Some(u) => {
            check_values("scale factors", u, n)?;
            (u.to_vec(), None)
        }
        None => {
            let solution =
                solve_scale_factors(mesh, geometry, input.cones, &options.scale_factors)?;
            (solution.u.clone(), Some(solution))
        }
    };
    timings.scale_factors = stage.elapsed().as_secs_f64();
    let cone_angle = solve.as_ref().map_or(cone_angle, |s| s.cone_angles.clone());

    // Step 2: cut.
    let stage = Instant::now();
    let is_cone: Vec<bool> = cone_position.iter().map(|&p| p != NONE).collect();
    let (plan, cut) = match input.cut {
        CutSource::Plan { may_cut, must_cut } => {
            let cone_vertices: Vec<u32> = input.cones.iter().map(|c| c.vertex).collect();
            let (plan, cut) = plan_and_cut(&CutPlanInput {
                mesh,
                edge_lengths: geometry.edge_lengths(),
                cone_vertices: &cone_vertices,
                may_cut,
                must_cut,
            })?;
            (Some(plan), cut)
        }
        CutSource::Edges(edges) => {
            let cut = cut_mesh(mesh, edges)?;
            check_disk_charts(&cut, &is_cone)?;
            (None, cut)
        }
    };
    timings.cut = stage.elapsed().as_secs_f64();

    // Step 3a: boundary data from u, with the Laplacian of the cut surface.
    let stage = Instant::now();
    let lengths = halfedge_lengths(mesh, geometry);
    let surface = Surface {
        mesh: &cut.mesh,
        cotans: geometry.corner_cotans(),
        angles: geometry.corner_angles(),
        areas: geometry.face_areas(),
        halfedge_lengths: &lengths,
    };
    let n_cut = cut.mesh.n_vertices();
    let u_cut = cut.lift_vertex_values(&u);
    let mut laplacian_u = vec![0.0; n_cut];
    surface.apply_laplacian(&u_cut, &mut laplacian_u);
    let angle_sums = surface.angle_sums();
    let on_boundary: Vec<bool> = (0..n_cut).map(|v| cut.mesh.is_boundary_vertex(v)).collect();

    // The cone metric must be flat inside the charts.
    let mut max_interior_curvature = 0.0f64;
    let mut worst_vertex = 0usize;
    for v in (0..n_cut).filter(|&v| !on_boundary[v]) {
        let curvature = (2.0 * PI - angle_sums[v] + laplacian_u[v]).abs();
        if curvature > max_interior_curvature || curvature.is_nan() {
            max_interior_curvature = curvature;
            worst_vertex = v;
        }
    }
    let not_flat =
        max_interior_curvature.is_nan() || max_interior_curvature > options.flatness_tolerance;
    if input.scale_factors.is_some() && not_flat {
        return Err(BffError::ScaleFactorsNotFlat {
            vertex: cut.vertex_to_original[worst_vertex],
            curvature: 2.0 * PI - angle_sums[worst_vertex] + laplacian_u[worst_vertex],
            tolerance: options.flatness_tolerance,
        });
    }

    let face_chart = cut.mesh.face_components().to_vec();
    let n_charts = cut.mesh.n_components();
    let boundary_loops = anchored_loops(mesh, &cut, &face_chart, input.anchors)?;
    let mut loop_of_chart = vec![NONE; n_charts];
    let mut loops: Vec<LoopData> = Vec::with_capacity(boundary_loops.len());
    for (index, boundary) in boundary_loops.iter().enumerate() {
        let chart = face_chart[cut.mesh.outgoing_halfedge(boundary[0] as usize) / 3];
        loop_of_chart[chart as usize] = index as u32;
        loops.push(cone_loop_data(
            mesh,
            &cut,
            boundary,
            &lengths,
            &u_cut,
            &angle_sums,
            &laplacian_u,
            shift,
        ));
    }
    let curves = best_fit_curves(
        &loops,
        &CurveOptions {
            max_length_adjustment: options.max_length_adjustment,
            max_coupled_loops: options.max_coupled_charts,
        },
    );
    timings.curves = stage.elapsed().as_secs_f64();

    // Step 3b: harmonic extension on all charts at once (their Laplacians do not interact).
    let stage = Instant::now();
    let mut values = vec![0.0; 2 * n_cut];
    for (boundary, curve) in boundary_loops.iter().zip(&curves) {
        for (&v, p) in boundary.iter().zip(&curve.positions) {
            values[v as usize] = p[0];
            values[n_cut + v as usize] = p[1];
        }
    }
    let system = PoissonSystem::new(
        "the Dirichlet Laplacian of the cut surface",
        &surface,
        &on_boundary,
        shift,
    )?;
    let extension_residual = system.solve(&surface, &mut values, 2)?;
    let vertex_uvs: Vec<[f64; 2]> = (0..n_cut).map(|v| [values[v], values[n_cut + v]]).collect();
    check_finite_layout(&cut.mesh, &vertex_uvs)?;
    let corner_uvs: Vec<[f64; 2]> = cut
        .mesh
        .faces()
        .iter()
        .flat_map(|face| face.map(|v| vertex_uvs[v as usize]))
        .collect();
    timings.extension = stage.elapsed().as_secs_f64();

    // Reports.
    let stage = Instant::now();
    let (quality, uv_angle_sums) = layout_quality(&surface, &vertex_uvs);
    let topology = component_topology(&cut.mesh);
    let charts: Vec<ChartReport> = quality
        .into_iter()
        .enumerate()
        .map(|(chart, (distortion, overlap))| {
            let loop_index = loop_of_chart[chart] as usize;
            ChartReport {
                first_face: topology[chart].first_face,
                n_vertices: topology[chart].n_vertices,
                n_boundary_vertices: boundary_loops[loop_index].len(),
                closure: curves[loop_index].report,
                distortion,
                overlap,
            }
        })
        .collect();

    // Cones: sums over the copies of each cone vertex.
    let mut wedges = vec![0usize; input.cones.len()];
    let mut uv_sum = vec![0.0; input.cones.len()];
    let mut max_wedge = vec![0.0f64; input.cones.len()];
    let mut implied = vec![2.0 * PI; input.cones.len()];
    for (v, &original) in cut.vertex_to_original.iter().enumerate() {
        let position = cone_position[original as usize];
        if position != NONE {
            let p = position as usize;
            wedges[p] += 1;
            uv_sum[p] += uv_angle_sums[v];
            max_wedge[p] = max_wedge[p].max(uv_angle_sums[v]);
            implied[p] += laplacian_u[v] - angle_sums[v];
        }
    }
    let cones: Vec<ConeReport> = input
        .cones
        .iter()
        .enumerate()
        .map(|(p, cone)| ConeReport {
            vertex: cone.vertex,
            angle: cone_angle[cone.vertex as usize],
            implied_angle: implied[p],
            realized_angle: 2.0 * PI - uv_sum[p],
            n_wedges: wedges[p],
            max_wedge_angle: max_wedge[p],
        })
        .collect();

    let seams = seam_report(mesh, &cut, &corner_uvs, &face_chart);
    timings.report = stage.elapsed().as_secs_f64();
    timings.total = start.elapsed().as_secs_f64();

    Ok(ConeFlattening {
        corner_uvs,
        vertex_uvs,
        cut,
        face_chart,
        charts,
        plan,
        scale_factors: u,
        scale_factor_solve: solve,
        max_interior_curvature,
        extension_residual,
        cones,
        seams,
        timings,
    })
}

/// The boundary loops of the cut surface, each rotated to start at its chart's anchor.
fn anchored_loops(
    mesh: &TriMesh,
    cut: &CutMesh,
    face_chart: &[u32],
    anchors: &[u32],
) -> Result<Vec<Vec<u32>>, BffError> {
    let mut loops: Vec<Vec<u32>> = cut.mesh.boundary_loops().to_vec();
    if anchors.is_empty() {
        return Ok(loops);
    }
    let n_charts = cut.mesh.n_components();
    let mut anchor_of_chart = vec![NONE; n_charts];
    for (position, &corner) in anchors.iter().enumerate() {
        if corner as usize >= mesh.n_halfedges() {
            return Err(BffError::AnchorOutOfRange {
                position,
                corner,
                n_corners: mesh.n_halfedges(),
            });
        }
        let chart = face_chart[corner as usize / 3];
        let previous = anchor_of_chart[chart as usize];
        if previous != NONE {
            return Err(BffError::DuplicateAnchor {
                chart,
                first: previous,
                second: corner,
            });
        }
        anchor_of_chart[chart as usize] = corner;
    }
    for boundary in &mut loops {
        let chart = face_chart[cut.mesh.outgoing_halfedge(boundary[0] as usize) / 3];
        let corner = anchor_of_chart[chart as usize];
        if corner == NONE {
            continue;
        }
        let vertex = cut.mesh.from_vertex(corner as usize) as u32;
        let position =
            boundary
                .iter()
                .position(|&v| v == vertex)
                .ok_or(BffError::AnchorNotOnBoundary {
                    corner,
                    vertex: cut.vertex_to_original[vertex as usize],
                })?;
        boundary.rotate_left(position);
    }
    Ok(loops)
}

/// The boundary data of one loop of the cut surface.
///
/// Exterior angles: `k̃_j = k_j − h_j` with `k_j = π − Σ θ` over the corners of the copy `j`
/// and Neumann data `h_j = −(L_cut u)_j` (paper Eq. 12 and 13 with the scale factors of the
/// uncut surface on both sides of the cut). Summed over the copies of a vertex this gives
/// `Σ_p (π − k̃_j^p) = 2π − µ_j` (paper Proposition 2). Target lengths:
/// `ℓ* = e^{(u_i + u_j)/2} ℓ`; the two sides of a cut edge share the degree of freedom of the
/// original edge.
///
/// [formula: bff-cherrier-boundary-curvature]
/// [formula: bff-target-boundary-length]
/// [formula: bff-cone-angle-sum]
#[allow(clippy::too_many_arguments)]
fn cone_loop_data(
    mesh: &TriMesh,
    cut: &CutMesh,
    boundary: &[u32],
    halfedge_lengths: &[f64],
    u_cut: &[f64],
    angle_sums: &[f64],
    laplacian_u: &[f64],
    shift: f64,
) -> LoopData {
    let n = boundary.len();
    let mut data = LoopData {
        dof: Vec::with_capacity(n),
        length: Vec::with_capacity(n),
        target: Vec::with_capacity(n),
        exterior: Vec::with_capacity(n),
    };
    for (i, &v) in boundary.iter().enumerate() {
        let (v, next) = (v as usize, boundary[(i + 1) % n] as usize);
        // The boundary halfedge leaving v; halfedge ids are those of the input mesh.
        let h = cut.mesh.outgoing_halfedge(v);
        let length = halfedge_lengths[h];
        data.dof.push(mesh.edge(h) as u32);
        data.length.push(length);
        data.target
            .push((0.5 * (u_cut[v] + u_cut[next])).exp() * length);
        // The reference implementation evaluates h with its shifted matrix.
        data.exterior
            .push(PI - angle_sums[v] + laplacian_u[v] + shift * u_cut[v]);
    }
    data
}

/// Compares the lengths of the two sides of every cut edge in the layout.
fn seam_report(
    mesh: &TriMesh,
    cut: &CutMesh,
    corner_uvs: &[[f64; 2]],
    face_chart: &[u32],
) -> SeamReport {
    let side = |h: usize| {
        let (a, b) = (corner_uvs[h], corner_uvs[mesh.next(h)]);
        (b[0] - a[0]).hypot(b[1] - a[1])
    };
    let mut report = SeamReport {
        n_cut_edges: cut.cut_edges.len(),
        n_between_charts: 0,
        max_length_mismatch: 0.0,
        worst_edge: None,
    };
    for &e in &cut.cut_edges {
        let h = mesh.edge_halfedge(e as usize);
        let t = mesh.twin(h).expect("cut edges are interior edges");
        if face_chart[h / 3] != face_chart[t / 3] {
            report.n_between_charts += 1;
        }
        let (a, b) = (side(h), side(t));
        let mismatch = (a - b).abs() / a.max(b);
        if mismatch > report.max_length_mismatch {
            report.max_length_mismatch = mismatch;
            report.worst_edge = Some(e);
        }
    }
    report
}
