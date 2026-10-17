//! The canonical discrete problem: the public problem description, its validation, and the
//! internal form on the unknown vertices.

use ddg_core::io::Cone;
use ddg_core::operators::{
    angle_defects, boundary_geodesic_curvatures, cotan_laplacian, lumped_mass_normalized,
};
use ddg_core::topology::SubMesh;
use ddg_core::{Geometry, SparseMatrix, TriMesh};

use crate::error::{ConeError, PrescribedConeIssue};
use crate::site::Site;

/// Tolerance of the Gauss–Bonnet feasibility test on closed meshes: the bounds are accepted if
/// they allow the required sum `2πχ` up to this absolute slack (so that bounds that allow
/// exactly `2πχ`, like `0 ≤ µ ≤ π/2` on a cube, are not rejected by rounding).
pub(crate) const FEASIBILITY_TOLERANCE: f64 = 1e-9;

/// What the solver does with the boundary of a mesh (paper Sec. 8.3, Fig. 10). It matters only
/// on meshes with boundary; on a closed mesh all variants give the same problem.
///
/// With boundary vertices `B` and interior vertices `I`, the boundary geodesic curvature of the
/// mesh is `k_i = π − Σ(corner angles at i)` and `ν_i` is the geodesic curvature of the
/// flattened boundary at `i`. Every treatment but [`Fixed`](BoundaryTreatment::Fixed) makes
/// `u_i` an unknown at the boundary vertices as well, with the constraint
///
/// ```text
/// (L u)_i = ν_i − k_i                                  (i ∈ B)
/// ```
///
/// next to `(L u)_i = µ_i − Ω_i` at the interior vertices, and a *site* for `ν_i` exactly like
/// the site of a cone angle `µ_i`. The distortion `½ Σ wE_i m_i u_i²` then runs over all
/// vertices, and `u` is defined up to a constant that is optimized too, as on a closed mesh.
/// `Σ_I µ + Σ_B ν = 2πχ` holds (Gauss–Bonnet).
///
/// [formula: boundary-optimality-system]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BoundaryTreatment {
    /// `u_i = 0` at the boundary vertices (paper Fig. 10b): the canonical problem of the crate
    /// documentation, and the behaviour of [`solve`](crate::solve).
    #[default]
    Fixed,
    /// `ν_i` is free at no cost: optimized boundary values (Fig. 10c). The adjoint vanishes on
    /// the boundary (a natural boundary condition).
    Free,
    /// `ν_i` costs `λ·wR_i·|ν_i|`, with `wR_i = ` [`ConeProblem::cone_costs`]`[i]` (default 1;
    /// `+∞` makes the boundary vertex straight, `ν_i = 0`): a polygonal boundary (Fig. 10d).
    Penalized,
    /// `ν_i ≥ 0` at no cost: a convex boundary (Fig. 10e).
    Convex,
}

/// The canonical discrete problem (see the [crate documentation](crate)):
///
/// ```text
/// minimize     ½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |µ_i|
/// subject to   (L u)_i = µ_i − Ω_i  (i interior),   u_i = 0  (i boundary),   lo_i ≤ µ_i ≤ hi_i
/// ```
///
/// Every per-vertex array has one entry per mesh vertex. Entries of boundary vertices are
/// validated like the others but have no effect: `u = 0` and `µ = 0` there. (The
/// `*_with_boundary` functions take a [`BoundaryTreatment`] that changes this; this struct
/// does not carry it, so that it stays constructible by struct literal as before. Then
/// `distortion_weights` apply at the boundary vertices, and for
/// [`BoundaryTreatment::Penalized`] `cone_costs` too.)
///
/// `m` is the lumped vertex mass divided by the total area, so `λ` is dimensionless and the
/// problem does not depend on the length unit.
#[derive(Debug, Clone)]
pub struct ConeProblem<'a> {
    /// The mesh: with or without boundary, with any number of connected components (each
    /// component is solved as its part of the whole problem).
    pub mesh: &'a TriMesh,
    /// Its intrinsic geometry (edge lengths, angles, areas).
    pub geometry: &'a Geometry,
    /// The parameter `λ ≥ 0` that trades distortion against total cone angle. The paper's
    /// default is `0.5`, its useful range about `0.1` to `1` (Sec. 6.2).
    pub lambda: f64,
    /// Distortion weights `wE_i ≥ 0`; `None` means 1 everywhere. A zero excludes the vertex
    /// from the distortion measure (then the minimizer need not be unique there).
    pub distortion_weights: Option<Vec<f64>>,
    /// Cone costs `wR_i ∈ [0, +∞]`; `None` means 1 everywhere. `0.0` makes the vertex a free
    /// cone site, `f64::INFINITY` forbids a cone there.
    pub cone_costs: Option<Vec<f64>>,
    /// Lower bounds `lo_i ≤ 0` of the cone angles; `None` means unbounded
    /// (`f64::NEG_INFINITY`) everywhere.
    pub lower_bounds: Option<Vec<f64>>,
    /// Upper bounds `hi_i ≥ 0` of the cone angles; `None` means unbounded (`f64::INFINITY`)
    /// everywhere.
    pub upper_bounds: Option<Vec<f64>>,
    /// Cones whose vertex and angle are fixed, at no cost. A prescribed cone overrides the cone
    /// cost and the bounds of its vertex. Each must sit at an interior vertex.
    pub prescribed_cones: Vec<Cone>,
}

impl<'a> ConeProblem<'a> {
    /// The basic problem of the paper: unit weights, no bounds, no prescribed cones.
    pub fn new(mesh: &'a TriMesh, geometry: &'a Geometry, lambda: f64) -> ConeProblem<'a> {
        ConeProblem {
            mesh,
            geometry,
            lambda,
            distortion_weights: None,
            cone_costs: None,
            lower_bounds: None,
            upper_bounds: None,
            prescribed_cones: Vec::new(),
        }
    }
}

/// The validated per-vertex data of a problem, in vertex numbering.
pub(crate) struct VertexData {
    /// Whether the vertex is an unknown (interior vertex, or any vertex of a closed mesh).
    pub unknown: Vec<bool>,
    /// Distortion weight `wE` per vertex.
    pub w_e: Vec<f64>,
    /// Site per vertex (a forbidden site at boundary vertices).
    pub sites: Vec<Site>,
    /// Number of forbidden unknown vertices.
    pub n_forbidden: usize,
    /// Whether the vertex carries a prescribed cone.
    pub prescribed: Vec<bool>,
}

impl VertexData {
    /// The data of the vertices `vertices` (in this order), e.g. those of one component.
    fn restrict(&self, vertices: &[u32]) -> VertexData {
        let pick = |v: &u32| *v as usize;
        let unknown: Vec<bool> = vertices.iter().map(|v| self.unknown[pick(v)]).collect();
        let sites: Vec<Site> = vertices.iter().map(|v| self.sites[pick(v)]).collect();
        let prescribed: Vec<bool> = vertices.iter().map(|v| self.prescribed[pick(v)]).collect();
        let n_forbidden = (0..vertices.len())
            .filter(|&i| unknown[i] && sites[i].fixed && !prescribed[i])
            .count();
        VertexData {
            unknown,
            w_e: vertices.iter().map(|v| self.w_e[pick(v)]).collect(),
            sites,
            n_forbidden,
            prescribed,
        }
    }
}

fn check_length(what: &'static str, array: &Option<Vec<f64>>, n: usize) -> Result<(), ConeError> {
    match array {
        Some(a) if a.len() != n => Err(ConeError::WrongLength {
            what,
            expected: n,
            got: a.len(),
        }),
        _ => Ok(()),
    }
}

/// Validates everything about a problem that does not need the geometry's values.
pub(crate) fn validate_with(
    problem: &ConeProblem<'_>,
    boundary: BoundaryTreatment,
) -> Result<VertexData, ConeError> {
    let mesh = problem.mesh;
    problem.geometry.check_mesh(mesh)?;
    let lambda = problem.lambda;
    if !(lambda.is_finite() && lambda >= 0.0) {
        return Err(ConeError::InvalidLambda { lambda });
    }
    let nv = mesh.n_vertices();
    check_length("distortion_weights", &problem.distortion_weights, nv)?;
    check_length("cone_costs", &problem.cone_costs, nv)?;
    check_length("lower_bounds", &problem.lower_bounds, nv)?;
    check_length("upper_bounds", &problem.upper_bounds, nv)?;

    let fixed_boundary = boundary == BoundaryTreatment::Fixed;
    let unknown: Vec<bool> = (0..nv)
        .map(|v| !fixed_boundary || !mesh.is_boundary_vertex(v))
        .collect();
    let mut w_e = vec![1.0; nv];
    let mut sites = Vec::with_capacity(nv);
    let mut n_forbidden = 0;
    for v in 0..nv {
        if let Some(weights) = &problem.distortion_weights {
            let value = weights[v];
            if !(value.is_finite() && value >= 0.0) {
                return Err(ConeError::InvalidDistortionWeight {
                    vertex: v as u32,
                    value,
                });
            }
            w_e[v] = value;
        }
        let cost = problem.cone_costs.as_ref().map_or(1.0, |c| c[v]);
        // `+∞` is allowed; NaN is not.
        if cost.is_nan() || cost < 0.0 {
            return Err(ConeError::InvalidConeCost {
                vertex: v as u32,
                value: cost,
            });
        }
        let lower = problem
            .lower_bounds
            .as_ref()
            .map_or(f64::NEG_INFINITY, |b| b[v]);
        let upper = problem
            .upper_bounds
            .as_ref()
            .map_or(f64::INFINITY, |b| b[v]);
        if !(lower <= 0.0 && upper >= 0.0) {
            return Err(ConeError::InvalidBounds {
                vertex: v as u32,
                lower,
                upper,
            });
        }
        if !unknown[v] {
            sites.push(Site::fixed(0.0));
        } else if mesh.is_boundary_vertex(v) {
            // A boundary vertex that is an unknown: its site is that of the treatment.
            match boundary {
                BoundaryTreatment::Fixed => unreachable!("boundary vertices are not unknowns"),
                BoundaryTreatment::Free => {
                    sites.push(Site::penalized(0.0, f64::NEG_INFINITY, f64::INFINITY))
                }
                BoundaryTreatment::Penalized if cost == f64::INFINITY => {
                    n_forbidden += 1;
                    sites.push(Site::fixed(0.0));
                }
                BoundaryTreatment::Penalized => sites.push(Site::penalized(
                    lambda * cost,
                    f64::NEG_INFINITY,
                    f64::INFINITY,
                )),
                BoundaryTreatment::Convex => sites.push(Site::penalized(0.0, 0.0, f64::INFINITY)),
            }
        } else if cost == f64::INFINITY {
            n_forbidden += 1;
            sites.push(Site::fixed(0.0));
        } else {
            sites.push(Site::penalized(lambda * cost, lower, upper));
        }
    }

    let mut seen = vec![false; nv];
    let mut n_prescribed_forbidden = 0;
    for (index, cone) in problem.prescribed_cones.iter().enumerate() {
        let issue = |reason| ConeError::InvalidPrescribedCone {
            index,
            vertex: cone.vertex,
            angle: cone.angle,
            reason,
        };
        let v = cone.vertex as usize;
        if v >= nv {
            return Err(issue(PrescribedConeIssue::VertexOutOfRange));
        }
        if !cone.angle.is_finite() {
            return Err(issue(PrescribedConeIssue::NonFiniteAngle));
        }
        if mesh.is_boundary_vertex(v) {
            return Err(issue(PrescribedConeIssue::BoundaryVertex));
        }
        if seen[v] {
            return Err(issue(PrescribedConeIssue::DuplicateVertex));
        }
        seen[v] = true;
        if sites[v].fixed {
            // A prescribed cone overrides "forbidden".
            n_prescribed_forbidden += 1;
        }
        sites[v] = Site::fixed(cone.angle);
    }
    Ok(VertexData {
        unknown,
        w_e,
        sites,
        n_forbidden: n_forbidden - n_prescribed_forbidden,
        prescribed: seen,
    })
}

/// One connected component of a mesh, as a mesh of its own (faces keep their corner slots) with
/// the geometry of the parent's edge lengths.
pub(crate) struct Piece {
    /// The component, with the maps to the parent's vertices and faces.
    pub sub: SubMesh,
    /// Its geometry: the parent's edge lengths, so every corner angle, cotangent and face area
    /// is bit for bit the parent's.
    pub geometry: Geometry,
}

/// Splits a mesh into its connected components ([`TriMesh::split_components`], numbered as
/// [`TriMesh::face_components`]: in order of their lowest face).
pub(crate) fn split_pieces(mesh: &TriMesh, geometry: &Geometry) -> Result<Vec<Piece>, ConeError> {
    geometry.check_mesh(mesh)?;
    let lengths = geometry.edge_lengths();
    mesh.split_components()
        .into_iter()
        .map(|sub| {
            let mut piece_lengths = vec![0.0; sub.mesh.n_edges()];
            for (f, &parent) in sub.face_to_parent.iter().enumerate() {
                for k in 0..3 {
                    piece_lengths[sub.mesh.edge(3 * f + k)] =
                        lengths[mesh.edge(3 * parent as usize + k)];
                }
            }
            let geometry = Geometry::from_edge_lengths(&sub.mesh, piece_lengths)?;
            Ok(Piece { sub, geometry })
        })
        .collect()
}

/// A validated problem on one mesh, or on one connected component of the mesh of a problem: the
/// per-vertex data the solver, the certificate and the evaluation work from.
///
/// For a component, every array is that of the whole problem restricted to the component's
/// vertices; in particular the masses are those of the whole mesh, **divided by the total area
/// of the whole mesh**, not by the component's area (see [`Discrete::from_instance`]).
pub(crate) struct Instance<'m> {
    /// The mesh (the whole mesh, or a component of it).
    pub mesh: &'m TriMesh,
    /// Its geometry.
    pub geometry: &'m Geometry,
    /// The boundary treatment.
    pub boundary: BoundaryTreatment,
    /// The validated per-vertex data.
    pub data: VertexData,
    /// The right-hand side `Ω` per vertex (`Ω + k` with unknown boundary vertices).
    pub omega: Vec<f64>,
    /// The lumped masses divided by the total area of the whole mesh.
    pub mass: Vec<f64>,
    /// For a component: the vertex of the whole mesh of each of its vertices.
    pub to_parent: Option<&'m [u32]>,
    /// For a component: its number among the components of the whole mesh.
    pub component: Option<u32>,
    /// Number and sum of the prescribed cones on this mesh (for the infeasibility message).
    pub prescribed: (usize, f64),
}

impl<'m> Instance<'m> {
    /// Validates a problem and evaluates its per-vertex data. The mesh may have several
    /// components; then this is the data of the whole problem, and [`Instance::part`] gives
    /// that of one component.
    ///
    /// [formula: lambda-area-normalization]
    pub fn new(
        problem: &ConeProblem<'m>,
        boundary: BoundaryTreatment,
    ) -> Result<Instance<'m>, ConeError> {
        let data = validate_with(problem, boundary)?;
        let omega = curvatures(problem.mesh, problem.geometry, boundary)?;
        let mass = lumped_mass_normalized(problem.mesh, problem.geometry)?;
        Ok(Instance {
            mesh: problem.mesh,
            geometry: problem.geometry,
            boundary,
            data,
            omega,
            mass,
            to_parent: None,
            component: None,
            prescribed: (
                problem.prescribed_cones.len(),
                problem.prescribed_cones.iter().map(|c| c.angle).sum(),
            ),
        })
    }

    /// The data of component `component` of the whole problem `self`, on its own mesh `piece`.
    /// The masses stay those of the whole mesh (divided by the whole area).
    ///
    /// [formula: component-separation]
    pub fn part<'p>(&self, piece: &'p Piece, component: u32) -> Instance<'p> {
        let vertices = &piece.sub.vertex_to_parent;
        let pick =
            |values: &[f64]| -> Vec<f64> { vertices.iter().map(|&v| values[v as usize]).collect() };
        let data = self.data.restrict(vertices);
        let prescribed = (0..vertices.len())
            .filter(|&i| data.prescribed[i])
            .fold((0, 0.0), |(n, s), i| (n + 1, s + data.sites[i].lo));
        Instance {
            mesh: &piece.sub.mesh,
            geometry: &piece.geometry,
            boundary: self.boundary,
            omega: pick(&self.omega),
            mass: pick(&self.mass),
            data,
            to_parent: Some(vertices),
            component: Some(component),
            prescribed,
        }
    }

    /// The vertex of the caller's mesh that vertex `v` of this instance is.
    pub fn label(&self, v: u32) -> u32 {
        self.to_parent.map_or(v, |p| p[v as usize])
    }
}

/// The problem restricted to its unknowns: the interior vertices of a mesh with boundary
/// (`u = 0` on the boundary eliminated), or all vertices of a closed mesh. Always on a
/// connected mesh: a problem on a mesh with several components is one of these per component
/// ([`Discrete::from_instance`] of [`Instance::part`]).
#[derive(Clone)]
pub(crate) struct Discrete {
    /// Number of unknown vertices.
    pub n: usize,
    /// Number of mesh vertices.
    pub n_vertices: usize,
    /// Whether every vertex is an unknown (a closed mesh, or a boundary treatment other than
    /// `Fixed`); then `a` is the whole Laplacian, singular: it annihilates constants.
    pub closed: bool,
    /// Per unknown: whether its vertex is a boundary vertex (only with an unknown boundary).
    pub is_boundary: Vec<bool>,
    /// For each unknown, its vertex (ascending).
    pub unknowns: Vec<u32>,
    /// For each unknown, its vertex in the caller's mesh: `unknowns` itself, unless this is a
    /// component of a larger mesh. Used in messages.
    pub labels: Vec<u32>,
    /// The Laplacian block of the unknowns: `L_II`, or all of `L` on a closed mesh.
    pub a: SparseMatrix,
    /// Angle defects `Ω` of the unknowns (the boundary geodesic curvature `k` at an unknown
    /// boundary vertex).
    pub omega: Vec<f64>,
    /// Area-normalized lumped masses `m` of the unknowns (divided by the area of the whole mesh
    /// of the problem, also for a component).
    pub mass: Vec<f64>,
    /// Distortion weights times mass, `w_i = wE_i m_i`.
    pub w: Vec<f64>,
    /// The sites of the unknowns.
    pub sites: Vec<Site>,
    /// `Σ Ω_i` over the unknowns; on a closed mesh this is `2πχ` up to rounding, and the cone
    /// angles must sum to it.
    pub total_curvature: f64,
    /// A bound for the rounding error of `total_curvature`.
    pub curvature_rounding: f64,
}

impl Discrete {
    /// Assembles the discrete form of a validated instance on a connected mesh: the whole
    /// problem of a connected mesh, or one component of a disconnected one.
    ///
    /// **Why a problem on a mesh with several components is solved component by component.**
    /// The vertices of different components share no edge, so `L` (and its block on the
    /// unknowns) is block diagonal; the objective `½ Σ wE m u² + λ Σ wR |µ|` is a sum over the
    /// components, and every constraint involves the vertices of one component only. So `(u, µ)`
    /// is a minimizer iff its restriction to every component minimizes that component's part of
    /// the objective under that component's constraints (`isMinimizer_sum_iff`; the optimality
    /// systems separate likewise, `kkt_sum_iff`). On a closed component the additive constant of
    /// `u` is a free variable of its own, and the cone angles of every closed component sum to
    /// its own `2πχ`. The part of the objective that belongs to a component has the masses of
    /// the whole problem, `a_i / A` with `A` the total area of the mesh, and that is what
    /// [`Instance::part`] keeps: each component is solved with the given `λ` and these masses.
    /// The equivalent alternative, a component normalized by its own area `A_c`, needs
    /// `λ · A / A_c` for the same minimizers, and its adjoint is `A / A_c` times the adjoint of
    /// the whole problem (`isMinimizer_scaleMass_withLam_iff`, `kkt_scaleMass_withLam_iff`).
    ///
    /// With a treatment other than `Fixed` the boundary vertices are unknowns whose "angle
    /// defect" is the boundary geodesic curvature `k`, so that one formula covers `µ` and `ν`.
    ///
    /// [formula: component-separation]
    /// [formula: boundary-optimality-system]
    pub fn from_instance(instance: &Instance<'_>) -> Result<Discrete, ConeError> {
        let data = &instance.data;
        let mesh = instance.mesh;
        let geometry = instance.geometry;
        let nv = mesh.n_vertices();
        let unknowns: Vec<u32> = (0..nv as u32)
            .filter(|&v| data.unknown[v as usize])
            .collect();
        let n = unknowns.len();
        let closed = n == nv;

        let laplacian = cotan_laplacian(mesh, geometry)?;
        let a = if closed {
            laplacian
        } else {
            laplacian.submatrix(&unknowns, &unknowns)?
        };
        let omega_v = &instance.omega;
        let mass_v = &instance.mass;
        let pick =
            |values: &[f64]| -> Vec<f64> { unknowns.iter().map(|&v| values[v as usize]).collect() };
        let omega = pick(omega_v);
        let mass = pick(mass_v);
        let w: Vec<f64> = unknowns
            .iter()
            .map(|&v| data.w_e[v as usize] * mass_v[v as usize])
            .collect();
        let sites: Vec<Site> = unknowns.iter().map(|&v| data.sites[v as usize]).collect();
        let total_curvature: f64 = omega.iter().sum();
        // Each angle defect is 2π minus a sum of corner angles, each good to a few ulp.
        let curvature_rounding = 64.0
            * f64::EPSILON
            * omega
                .iter()
                .map(|o| o.abs() + 2.0 * std::f64::consts::PI)
                .sum::<f64>();

        if closed {
            if w.iter().all(|&x| x == 0.0) {
                return Err(ConeError::NoDistortionWeight {
                    component: instance.component,
                });
            }
            check_gauss_bonnet_feasible(instance, &sites, total_curvature)?;
        }
        let is_boundary: Vec<bool> = unknowns
            .iter()
            .map(|&v| mesh.is_boundary_vertex(v as usize))
            .collect();
        let labels = unknowns.iter().map(|&v| instance.label(v)).collect();
        Ok(Discrete {
            n,
            n_vertices: nv,
            closed,
            is_boundary,
            unknowns,
            labels,
            a,
            omega,
            mass,
            w,
            sites,
            total_curvature,
            curvature_rounding,
        })
    }

    /// The same problem with the cone angle of every interior unknown fixed to that of `cones`
    /// (zero where there is none), at no cost; the sites of the boundary vertices are kept.
    /// With unknown boundary values this is "optimize the boundary for given cones".
    pub fn with_fixed_interior(&self, cones: &[Cone]) -> Result<Discrete, ConeError> {
        let mut unknown_of = vec![u32::MAX; self.n_vertices];
        for (i, &v) in self.unknowns.iter().enumerate() {
            unknown_of[v as usize] = i as u32;
        }
        let mut angle = vec![0.0; self.n];
        for (index, cone) in cones.iter().enumerate() {
            let i = unknown_of
                .get(cone.vertex as usize)
                .copied()
                .unwrap_or(u32::MAX);
            if i == u32::MAX || self.is_boundary[i as usize] || !cone.angle.is_finite() {
                return Err(ConeError::InvalidCone {
                    index,
                    vertex: cone.vertex,
                    angle: cone.angle,
                    reason: PrescribedConeIssue::BoundaryVertex,
                });
            }
            angle[i as usize] = cone.angle;
        }
        let mut out = self.clone();
        for ((site, &boundary), &a) in out.sites.iter_mut().zip(&self.is_boundary).zip(&angle) {
            if !boundary {
                *site = Site::fixed(a);
            }
        }
        Ok(out)
    }

    /// Spreads a vector on the unknowns to all vertices (zero on the boundary).
    pub fn to_vertices(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.n_vertices];
        for (&v, &value) in self.unknowns.iter().zip(x) {
            out[v as usize] = value;
        }
        out
    }

    /// The distortion `½ Σ w_i u_i²`.
    pub fn distortion(&self, u: &[f64]) -> f64 {
        0.5 * self.w.iter().zip(u).map(|(&w, &u)| w * u * u).sum::<f64>()
    }

    /// The cone cost `Σ t_i |µ_i|` (fixed sites cost nothing).
    pub fn cone_cost(&self, mu: &[f64]) -> f64 {
        self.sites
            .iter()
            .zip(mu)
            .map(|(site, &mu)| site.cost(mu))
            .sum()
    }
}

/// The right-hand-side data of the state equation per vertex: the angle defect `Ω` at interior
/// vertices, and with unknown boundary vertices the geodesic curvature `k` at those (so that
/// the entries sum to `2πχ`). Zero at boundary vertices for [`BoundaryTreatment::Fixed`].
pub(crate) fn curvatures(
    mesh: &TriMesh,
    geometry: &Geometry,
    boundary: BoundaryTreatment,
) -> Result<Vec<f64>, ddg_core::GeometryError> {
    let mut out = angle_defects(mesh, geometry)?;
    if boundary != BoundaryTreatment::Fixed {
        let k = boundary_geodesic_curvatures(mesh, geometry)?;
        for (o, k) in out.iter_mut().zip(k) {
            *o += k;
        }
    }
    Ok(out)
}

/// On a closed mesh, summing the state equation gives `Σ µ_i = Σ Ω_i = 2πχ`; the constraints
/// must allow that sum. On a mesh with several components this holds for every closed
/// component on its own (its block of `L` annihilates the constants of that component).
///
/// [formula: closed-mesh-cone-sum]
fn check_gauss_bonnet_feasible(
    instance: &Instance<'_>,
    sites: &[Site],
    required: f64,
) -> Result<(), ConeError> {
    let min_sum: f64 = sites.iter().map(|s| s.lo).sum();
    let max_sum: f64 = sites.iter().map(|s| s.hi).sum();
    if min_sum - FEASIBILITY_TOLERANCE <= required && required <= max_sum + FEASIBILITY_TOLERANCE {
        return Ok(());
    }
    Err(ConeError::Infeasible {
        required,
        min_sum,
        max_sum,
        euler_characteristic: instance.mesh.euler_characteristic(),
        n_forbidden: instance.data.n_forbidden,
        n_prescribed: instance.prescribed.0,
        prescribed_sum: instance.prescribed.1,
        component: instance.component,
    })
}
