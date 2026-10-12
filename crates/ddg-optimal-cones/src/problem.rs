//! The canonical discrete problem: the public problem description, its validation, and the
//! internal form on the unknown vertices.

use ddg_core::io::Cone;
use ddg_core::operators::{angle_defects, cotan_laplacian, lumped_mass_normalized};
use ddg_core::{Geometry, SparseMatrix, TriMesh};

use crate::error::{ConeError, PrescribedConeIssue};
use crate::site::Site;

/// Tolerance of the Gauss–Bonnet feasibility test on closed meshes: the bounds are accepted if
/// they allow the required sum `2πχ` up to this absolute slack (so that bounds that allow
/// exactly `2πχ`, like `0 ≤ µ ≤ π/2` on a cube, are not rejected by rounding).
pub(crate) const FEASIBILITY_TOLERANCE: f64 = 1e-9;

/// The canonical discrete problem (see the [crate documentation](crate)):
///
/// ```text
/// minimize     ½ Σ_i wE_i m_i u_i²  +  λ Σ_i wR_i |µ_i|
/// subject to   (L u)_i = µ_i − Ω_i  (i interior),   u_i = 0  (i boundary),   lo_i ≤ µ_i ≤ hi_i
/// ```
///
/// Every per-vertex array has one entry per mesh vertex. Entries of boundary vertices are
/// validated like the others but have no effect: `u = 0` and `µ = 0` there.
///
/// `m` is the lumped vertex mass divided by the total area, so `λ` is dimensionless and the
/// problem does not depend on the length unit.
#[derive(Debug, Clone)]
pub struct ConeProblem<'a> {
    /// The mesh: connected, with or without boundary.
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
pub(crate) fn validate(problem: &ConeProblem<'_>) -> Result<VertexData, ConeError> {
    let mesh = problem.mesh;
    problem.geometry.check_mesh(mesh)?;
    if mesh.n_components() != 1 {
        return Err(ConeError::NotConnected {
            components: mesh.n_components(),
        });
    }
    let lambda = problem.lambda;
    if !(lambda.is_finite() && lambda >= 0.0) {
        return Err(ConeError::InvalidLambda { lambda });
    }
    let nv = mesh.n_vertices();
    check_length("distortion_weights", &problem.distortion_weights, nv)?;
    check_length("cone_costs", &problem.cone_costs, nv)?;
    check_length("lower_bounds", &problem.lower_bounds, nv)?;
    check_length("upper_bounds", &problem.upper_bounds, nv)?;

    let unknown: Vec<bool> = (0..nv).map(|v| !mesh.is_boundary_vertex(v)).collect();
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
        } else if cost == f64::INFINITY {
            n_forbidden += 1;
            sites.push(Site::fixed(0.0));
        } else {
            sites.push(Site::penalized(lambda * cost, lower, upper));
        }
    }

    let mut seen = vec![false; nv];
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
        if !unknown[v] {
            return Err(issue(PrescribedConeIssue::BoundaryVertex));
        }
        if seen[v] {
            return Err(issue(PrescribedConeIssue::DuplicateVertex));
        }
        seen[v] = true;
        if sites[v].fixed {
            // A prescribed cone overrides "forbidden".
            n_forbidden -= 1;
        }
        sites[v] = Site::fixed(cone.angle);
    }
    Ok(VertexData {
        unknown,
        w_e,
        sites,
        n_forbidden,
    })
}

/// The problem restricted to its unknowns: the interior vertices of a mesh with boundary
/// (`u = 0` on the boundary eliminated), or all vertices of a closed mesh.
pub(crate) struct Discrete {
    /// Number of unknown vertices.
    pub n: usize,
    /// Number of mesh vertices.
    pub n_vertices: usize,
    /// Whether the mesh has no boundary (then `a` is singular: it annihilates constants).
    pub closed: bool,
    /// For each unknown, its vertex (ascending).
    pub unknowns: Vec<u32>,
    /// The Laplacian block of the unknowns: `L_II`, or all of `L` on a closed mesh.
    pub a: SparseMatrix,
    /// Angle defects `Ω` of the unknowns.
    pub omega: Vec<f64>,
    /// Area-normalized lumped masses `m` of the unknowns.
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
    /// Validates the problem and assembles its discrete form.
    ///
    /// The masses are divided by the total area, which makes `λ` dimensionless.
    ///
    /// [formula: lambda-area-normalization]
    pub fn new(problem: &ConeProblem<'_>) -> Result<Discrete, ConeError> {
        let data = validate(problem)?;
        let mesh = problem.mesh;
        let geometry = problem.geometry;
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
        let omega_v = angle_defects(mesh, geometry)?;
        let mass_v = lumped_mass_normalized(mesh, geometry)?;
        let pick =
            |values: &[f64]| -> Vec<f64> { unknowns.iter().map(|&v| values[v as usize]).collect() };
        let omega = pick(&omega_v);
        let mass = pick(&mass_v);
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
                return Err(ConeError::NoDistortionWeight);
            }
            check_gauss_bonnet_feasible(problem, &data, &sites, total_curvature)?;
        }
        Ok(Discrete {
            n,
            n_vertices: nv,
            closed,
            unknowns,
            a,
            omega,
            mass,
            w,
            sites,
            total_curvature,
            curvature_rounding,
        })
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

/// On a closed mesh, summing the state equation gives `Σ µ_i = Σ Ω_i = 2πχ`; the constraints
/// must allow that sum.
///
/// [formula: closed-mesh-cone-sum]
fn check_gauss_bonnet_feasible(
    problem: &ConeProblem<'_>,
    data: &VertexData,
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
        euler_characteristic: problem.mesh.euler_characteristic(),
        n_forbidden: data.n_forbidden,
        n_prescribed: problem.prescribed_cones.len(),
        prescribed_sum: problem.prescribed_cones.iter().map(|c| c.angle).sum(),
    })
}
