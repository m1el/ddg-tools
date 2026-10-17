//! Conformal map of a closed surface of genus 0 to the unit sphere (paper Sec. 7.5,
//! "Topology": "For topological spheres one can remove a single vertex star and map to the
//! circular disk, then to the unit sphere via stereographic projection").
//!
//! Ported from `BFF::mapToSphere`, `BFF::projectStereographically` and `BFF::centerMobius` in
//! `Bff.cpp` of the reference implementation (MIT, © 2017 Rohan Sawhney and Keenan Crane). The
//! centering is Algorithm 1 of A. Baden, K. Crane, M. Kazhdan, *Möbius Registration*, Computer
//! Graphics Forum 37(5) (SGP 2018), which the reference cites.
//!
//! # The steps
//!
//! 1. **Pole.** One vertex, the *pole*, is chosen ([`PoleChoice`]) and its star (the pole and
//!    the faces around it) is removed. What is left is a disk whose boundary is the link of the
//!    pole.
//! 2. **Disk.** The disk is flattened to a round disk ([`DiskBff::flatten_to_disk`], paper
//!    Sec. 6.5), then translated and scaled so that its boundary vertices are on average at
//!    distance 1 from the origin.
//! 3. **Stereographic projection.** The disk, scaled by a radius `R`, is mapped to the sphere
//!    by the inverse stereographic projection from the south pole,
//!    `σ(X, Y) = (2X, 2Y, 1 − X² − Y²) / (1 + X² + Y²)`; the pole goes to the south pole
//!    `(0, 0, −1)`. `R` is found by a search that maximizes the area-weighted mean of the log
//!    area ratio, as in the reference ([`RadiusSearch`]).
//! 4. **Möbius centering.** Inversions of the sphere `x ↦ (1 − |c|²)(x + c)/|x + c|² + c` are
//!    applied until the center of mass of the surface's area, carried to the sphere, is at the
//!    origin ([`mobius_center`]). This makes the map canonical up to a rotation and removes the
//!    "low-frequency" part of the area distortion.
//!
//! Every step but the first and the last is conformal in the smooth setting (stereographic
//! projection and inversions are conformal), so the result is as conformal as the disk
//! flattening.
//!
//! # Where the result sits
//!
//! Positions are on the unit sphere up to rounding ([`SphereReport::max_radius_error`]); they
//! are not renormalized. With centering the map is determined up to a rotation of the sphere;
//! the rotation that comes out is whatever the steps above give (no alignment to the input).
//! Without centering the pole is at `(0, 0, −1)`.

use std::f64::consts::PI;
use std::time::Instant;

use ddg_core::geometry::Geometry;
use ddg_core::topology::TriMesh;

use crate::chart::ChartReport;
use crate::cones::check_option;
use crate::disk::{DiskBff, DiskOptions};
use crate::error::BffError;
use crate::quality::{DistortionAccumulator, DistortionReport};

/// Which vertex is removed (and becomes the pole of the stereographic projection).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PoleChoice {
    /// The vertex whose star is the most regular (the default); see [`most_regular_star`].
    /// Ties go to the lowest index. The choice depends on the edge lengths only, so it does not
    /// change under rigid motions and scaling of the input.
    #[default]
    MostRegularStar,
    /// This vertex. The reference implementation removes the first vertex of its mesh, i.e.
    /// `Vertex(0)` for a closed mesh.
    Vertex(u32),
}

/// The search for the radius of the disk before stereographic projection (reference
/// `BFF::mapToSphere`).
///
/// The disk, normalized to mean boundary radius 1, is scaled by `R` and projected; the objective
/// is the area-weighted mean over all faces of `ln(a_f / A_f)`, with `A_f` the area on the
/// surface and `a_f` the area of the flat triangle spanned by the image vertices. It is
/// maximized by ternary search over `ln R` in `[ln min, ln max]`, until the interval is shorter
/// than `tolerance`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiusSearch {
    /// Smallest radius tried (`> 0`). `1` puts the disk on the northern hemisphere.
    pub min: f64,
    /// Largest radius tried (`>= min`).
    pub max: f64,
    /// Length of the final interval of `ln R` (`> 0`).
    pub tolerance: f64,
}

impl Default for RadiusSearch {
    fn default() -> Self {
        RadiusSearch {
            min: 1.0,
            max: 1e4,
            tolerance: 1e-3,
        }
    }
}

/// Options of [`mobius_center`] (Baden, Crane, Kazhdan 2018, Algorithm 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenteringOptions {
    /// Largest number of inversions.
    pub max_iterations: usize,
    /// Stop when the center of mass `µ` (of area weights that sum to 1) has `|µ| <=` this.
    pub tolerance: f64,
    /// Largest `|c|` of one inversion, in `(0, 1)`. A Gauss–Newton step that is longer is
    /// shortened to this length (the paper suggests a shorter step "for larger or more
    /// challenging models"); [`CenteringReport::shortened_steps`] counts how often.
    pub max_step: f64,
}

impl Default for CenteringOptions {
    fn default() -> Self {
        CenteringOptions {
            max_iterations: 100,
            tolerance: 1e-10,
            max_step: 0.5,
        }
    }
}

/// Options of [`map_to_sphere`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereOptions {
    /// Which vertex becomes the pole.
    pub pole: PoleChoice,
    /// Iterations of the flattening to a round disk ([`DiskBff::flatten_to_disk`]); the
    /// reference uses 10.
    pub disk_iterations: usize,
    /// See [`DiskOptions::diagonal_shift`]: `0` solves exactly,
    /// [`crate::REFERENCE_DIAGONAL_SHIFT`] reproduces the reference.
    pub diagonal_shift: f64,
    /// The search for the radius of the disk; `None` projects the normalized disk as it is
    /// (`R = 1`).
    pub radius_search: Option<RadiusSearch>,
    /// The Möbius centering; `None` leaves the pole at `(0, 0, −1)`.
    pub centering: Option<CenteringOptions>,
}

impl Default for SphereOptions {
    fn default() -> Self {
        SphereOptions {
            pole: PoleChoice::MostRegularStar,
            disk_iterations: 10,
            diagonal_shift: 0.0,
            radius_search: Some(RadiusSearch::default()),
            centering: Some(CenteringOptions::default()),
        }
    }
}

/// What the radius search found.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiusReport {
    /// The radius `R` used (relative to the normalized disk).
    pub radius: f64,
    /// Number of steps of the ternary search (two evaluations each).
    pub iterations: usize,
    /// The objective at `R`: area-weighted mean of `ln(a_f / A_f)`.
    pub objective: f64,
}

/// Why the centering stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CenteringStatus {
    /// `|µ| <= tolerance`.
    Converged,
    /// `max_iterations` inversions were applied and `|µ|` is still above the tolerance.
    MaxIterations,
    /// No step along the Gauss–Newton direction, halved up to 40 times, decreased `|µ|`
    /// (seen on meshes of a handful of faces, where the face centers are a poor stand-in for
    /// the surface). The positions are those of the last accepted step.
    Stalled,
    /// The Jacobian could not be inverted (all face centers on one line through the origin,
    /// or not finite). The positions are those before the failed step.
    SingularJacobian,
}

/// What the Möbius centering did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenteringReport {
    /// Why it stopped.
    pub status: CenteringStatus,
    /// Number of inversions applied.
    pub iterations: usize,
    /// `|µ|` before the first inversion.
    pub initial_center: f64,
    /// `|µ|` at the end.
    pub final_center: f64,
    /// Number of steps that were shortened to `max_step`.
    pub shortened_steps: usize,
    /// Number of halvings of the step by the line search (a step is accepted only if it
    /// decreases `|µ|`).
    pub halvings: usize,
}

impl CenteringReport {
    /// Whether the center of mass reached the tolerance.
    pub fn converged(&self) -> bool {
        self.status == CenteringStatus::Converged
    }
}

/// Wall-clock seconds spent in the stages of [`map_to_sphere`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SphereTimings {
    /// Choosing the pole and building the punctured disk.
    pub puncture: f64,
    /// Flattening the disk to a round disk.
    pub disk: f64,
    /// Radius search and stereographic projection.
    pub projection: f64,
    /// Möbius centering.
    pub centering: f64,
    /// The distortion report.
    pub report: f64,
    /// Everything.
    pub total: f64,
}

/// The distortion of a map to the sphere, measured on the flat triangles spanned by the image
/// vertices (as the reference measures it).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereDistortion {
    /// Per face: quasi-conformal error `Q = σ_max / σ_min` of the affine map from the triangle
    /// on the surface to the image triangle, and the log scale factor
    /// `u = ½ ln((a_f / Σa) / (A_f / ΣA))`, i.e. relative to a uniform scaling of the total
    /// area `ΣA` to the total image area `Σa`, so that `u = 0` for no area distortion.
    /// `flipped_triangles` counts faces whose image is oriented clockwise seen from outside
    /// the sphere (`det(x_0, x_1, x_2) < 0`); `plane_area` is the signed total image area
    /// after that scaling (equal to `surface_area` if nothing is flipped).
    pub faces: DistortionReport,
    /// `Σa`, the total area of the image triangles (slightly below `4π`).
    pub image_area: f64,
    /// The sum of the signed areas of the geodesic triangles with the image vertices, divided
    /// by `4π`: the degree of the map. `1` for a map onto the sphere that is injective.
    pub covering_degree: f64,
    /// Largest `| |x| − 1 |` over the image vertices.
    pub max_radius_error: f64,
}

impl SphereDistortion {
    /// Whether the map is a bijection onto the sphere (with geodesic triangles): no flipped or
    /// degenerate face and covering degree 1 (within `1e-6`). An orientation preserving map of
    /// a sphere to the sphere of degree 1 covers every point exactly once.
    ///
    /// [formula: sphere-injectivity-criterion]
    pub fn is_injective(&self) -> bool {
        self.faces.flipped_triangles == 0
            && self.faces.degenerate_triangles == 0
            && (self.covering_degree - 1.0).abs() < 1e-6
    }
}

/// What [`map_to_sphere`] reports.
#[derive(Debug, Clone, PartialEq)]
pub struct SphereReport {
    /// The vertex that was removed and mapped to the pole.
    pub pole: u32,
    /// The flattening of the punctured surface to a round disk (indices of the punctured mesh).
    pub disk: ChartReport,
    /// Largest relative residual of the linear solves of the disk flattening.
    pub disk_solve_residual: f64,
    /// The radius search, if it ran.
    pub radius: Option<RadiusReport>,
    /// The Möbius centering, if it ran.
    pub centering: Option<CenteringReport>,
    /// Angle and area distortion, flipped triangles, degree.
    pub distortion: SphereDistortion,
    /// Where the time went.
    pub timings: SphereTimings,
}

impl SphereReport {
    /// See [`SphereDistortion::is_injective`].
    pub fn is_injective(&self) -> bool {
        self.distortion.is_injective()
    }
}

/// A map of a closed surface of genus 0 to the unit sphere.
#[derive(Debug, Clone, PartialEq)]
pub struct SphereMap {
    /// One point of the unit sphere per vertex of the input mesh (the pole included).
    pub positions: Vec<[f64; 3]>,
    /// What the computation says about itself.
    pub report: SphereReport,
}

/// Maps a closed, connected surface of genus 0 conformally to the unit sphere; see the
/// [module documentation](self).
///
/// # Errors
/// [`BffError::Geometry`] if the geometry does not belong to the mesh; [`BffError::NotASphere`]
/// unless the mesh is connected and closed with Euler characteristic 2 (a disconnected mesh is
/// refused: map its components one by one, `TriMesh::split_components`);
/// [`BffError::PoleOutOfRange`]; [`BffError::NonFiniteInput`] for options out of their domain;
/// whatever the disk flattening reports ([`BffError::Solve`], [`BffError::NonFiniteResult`],
/// [`BffError::DegenerateLayout`]).
pub fn map_to_sphere(
    mesh: &TriMesh,
    geometry: &Geometry,
    options: &SphereOptions,
) -> Result<SphereMap, BffError> {
    let start = Instant::now();
    geometry.check_mesh(mesh)?;
    check_sphere(mesh)?;
    check_sphere_options(options)?;

    // 1. The pole and the punctured disk.
    let pole = match options.pole {
        PoleChoice::Vertex(v) => {
            if v as usize >= mesh.n_vertices() {
                return Err(BffError::PoleOutOfRange {
                    vertex: v,
                    n_vertices: mesh.n_vertices(),
                });
            }
            v
        }
        PoleChoice::MostRegularStar => most_regular_star(mesh, geometry),
    };
    let kept: Vec<u32> = (0..mesh.n_faces() as u32)
        .filter(|&f| !mesh.faces()[f as usize].contains(&pole))
        .collect();
    let punctured = mesh.submesh(&kept).map_err(BffError::Mesh)?;
    let lengths: Vec<f64> = punctured
        .mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let (a, b) = (
                punctured.vertex_to_parent[a as usize] as usize,
                punctured.vertex_to_parent[b as usize] as usize,
            );
            mesh.find_edge(a, b)
                .map(|e| geometry.edge_lengths()[e])
                .ok_or(BffError::Internal {
                    what: "an edge of the punctured mesh is not an edge of the input",
                })
        })
        .collect::<Result<_, _>>()?;
    let punctured_geometry = Geometry::from_edge_lengths(&punctured.mesh, lengths)?;
    let t_puncture = start.elapsed().as_secs_f64();

    // 2. The round disk.
    let t = Instant::now();
    let mut bff = DiskBff::new(
        &punctured.mesh,
        &punctured_geometry,
        DiskOptions {
            diagonal_shift: options.diagonal_shift,
            overlap_samples_per_triangle: 0.0,
            ..DiskOptions::default()
        },
    )?;
    let disk = bff.flatten_to_disk(options.disk_iterations)?;
    let boundary = bff.boundary().to_vec();
    let t_disk = t.elapsed().as_secs_f64();

    // 3. Normalize, search the radius, project.
    let t = Instant::now();
    let n = mesh.n_vertices();
    let mut plane = vec![[0.0; 2]; n];
    {
        let nb = boundary.len() as f64;
        let center = boundary.iter().fold([0.0, 0.0], |c, &v| {
            let p = disk.uv[v as usize];
            [c[0] + p[0] / nb, c[1] + p[1] / nb]
        });
        let radius = boundary
            .iter()
            .map(|&v| {
                let p = disk.uv[v as usize];
                (p[0] - center[0]).hypot(p[1] - center[1])
            })
            .sum::<f64>()
            / nb;
        if !(radius > 0.0 && radius.is_finite()) {
            return Err(BffError::NonFiniteResult {
                what: "the radius of the flattened disk",
                chart: Some(0),
            });
        }
        for (v, p) in disk.uv.iter().enumerate() {
            plane[punctured.vertex_to_parent[v] as usize] =
                [(p[0] - center[0]) / radius, (p[1] - center[1]) / radius];
        }
    }
    let project = |scale: f64, out: &mut Vec<[f64; 3]>| {
        out.clear();
        out.extend((0..n).map(|v| {
            if v == pole as usize {
                [0.0, 0.0, -1.0]
            } else {
                let p = plane[v];
                inverse_stereographic([scale * p[0], scale * p[1]])
            }
        }));
    };
    let mut positions = Vec::with_capacity(n);
    let radius_report = match options.radius_search {
        None => {
            project(1.0, &mut positions);
            None
        }
        Some(search) => {
            let objective = |positions: &[[f64; 3]]| log_area_ratio_mean(mesh, geometry, positions);
            let (mut lo, mut hi) = (search.min.ln(), search.max.ln());
            let mut iterations = 0;
            while hi - lo > search.tolerance {
                let left = lo + (hi - lo) / 3.0;
                let right = hi - (hi - lo) / 3.0;
                project(left.exp(), &mut positions);
                let f_left = objective(&positions);
                project(right.exp(), &mut positions);
                let f_right = objective(&positions);
                iterations += 1;
                // Maximize; a NaN on either side ends the search where it is, as in the
                // reference.
                if f_left.is_nan() || f_right.is_nan() {
                    break;
                }
                if f_left < f_right {
                    lo = left;
                } else {
                    hi = right;
                }
            }
            let radius = (0.5 * (lo + hi)).exp();
            project(radius, &mut positions);
            Some(RadiusReport {
                radius,
                iterations,
                objective: objective(&positions),
            })
        }
    };
    let t_projection = t.elapsed().as_secs_f64();

    // 4. Möbius centering.
    let t = Instant::now();
    let centering = match &options.centering {
        Some(centering) => Some(mobius_center(mesh, geometry, &mut positions, centering)?),
        None => None,
    };
    let t_centering = t.elapsed().as_secs_f64();

    if positions.iter().any(|p| !p.iter().all(|x| x.is_finite())) {
        return Err(BffError::NonFiniteResult {
            what: "a position on the sphere",
            chart: None,
        });
    }
    let t = Instant::now();
    let distortion = sphere_distortion(mesh, geometry, &positions)?;
    let t_report = t.elapsed().as_secs_f64();

    Ok(SphereMap {
        positions,
        report: SphereReport {
            pole,
            disk: disk.report,
            disk_solve_residual: disk.solve_residual,
            radius: radius_report,
            centering,
            distortion,
            timings: SphereTimings {
                puncture: t_puncture,
                disk: t_disk,
                projection: t_projection,
                centering: t_centering,
                report: t_report,
                total: start.elapsed().as_secs_f64(),
            },
        },
    })
}

/// Refuses meshes that are not a connected closed surface of genus 0.
fn check_sphere(mesh: &TriMesh) -> Result<(), BffError> {
    if mesh.n_components() != 1
        || !mesh.boundary_loops().is_empty()
        || mesh.euler_characteristic() != 2
    {
        return Err(BffError::NotASphere {
            components: mesh.n_components(),
            euler_characteristic: mesh.euler_characteristic(),
            boundary_loops: mesh.boundary_loops().len(),
        });
    }
    Ok(())
}

fn check_sphere_options(options: &SphereOptions) -> Result<(), BffError> {
    check_option(
        "diagonal_shift (must be finite and non-negative)",
        options.diagonal_shift,
        true,
    )?;
    if let Some(search) = &options.radius_search {
        let bad = |what, value| BffError::NonFiniteInput {
            what,
            index: 0,
            value,
        };
        if !(search.min > 0.0 && search.min.is_finite()) {
            return Err(bad(
                "radius_search.min (must be positive and finite)",
                search.min,
            ));
        }
        if !(search.max >= search.min && search.max.is_finite()) {
            return Err(bad(
                "radius_search.max (must be finite and >= min)",
                search.max,
            ));
        }
        if !(search.tolerance > 0.0 && search.tolerance.is_finite()) {
            return Err(bad(
                "radius_search.tolerance (must be positive and finite)",
                search.tolerance,
            ));
        }
    }
    if let Some(centering) = &options.centering {
        check_centering_options(centering)?;
    }
    Ok(())
}

fn check_centering_options(centering: &CenteringOptions) -> Result<(), BffError> {
    check_option(
        "centering.tolerance (must be finite and non-negative)",
        centering.tolerance,
        true,
    )?;
    if !(centering.max_step > 0.0 && centering.max_step < 1.0) {
        return Err(BffError::NonFiniteInput {
            what: "centering.max_step (must be in (0, 1))",
            index: 0,
            value: centering.max_step,
        });
    }
    Ok(())
}

/// The vertex whose star is the most regular, the default pole.
///
/// For a vertex with `n` corners of angles `θ_j` and `n` edges of lengths `r_j` the score is
/// `|Ω| / 2π + max_j |n θ_j / Σθ − 1| + (max_j r_j / min_j r_j − 1)`, with `Ω = 2π − Σθ` the
/// angle defect: zero for the star of a regular polygon in the plane. The lowest score wins;
/// ties go to the lowest index. The removed star becomes a round cap of the sphere, so a flat
/// star with equal angles and equal edges is the one the round disk can take best. This is a
/// heuristic of this crate; the reference removes its first vertex.
pub fn most_regular_star(mesh: &TriMesh, geometry: &Geometry) -> u32 {
    let n = mesh.n_vertices();
    let mut angles: Vec<Vec<f64>> = vec![Vec::new(); n];
    for (f, face) in mesh.faces().iter().enumerate() {
        for (k, &v) in face.iter().enumerate() {
            angles[v as usize].push(geometry.corner_angles()[f][k]);
        }
    }
    let mut lengths: Vec<(f64, f64)> = vec![(f64::INFINITY, 0.0); n];
    for (e, &[a, b]) in mesh.edges().iter().enumerate() {
        let l = geometry.edge_lengths()[e];
        for v in [a, b] {
            let entry = &mut lengths[v as usize];
            entry.0 = entry.0.min(l);
            entry.1 = entry.1.max(l);
        }
    }
    let mut best = (f64::INFINITY, 0u32);
    for v in 0..n {
        let theta = &angles[v];
        let sum: f64 = theta.iter().sum();
        let count = theta.len() as f64;
        let spread = theta
            .iter()
            .map(|t| (count * t / sum - 1.0).abs())
            .fold(0.0, f64::max);
        let (lo, hi) = lengths[v];
        let score = (2.0 * PI - sum).abs() / (2.0 * PI) + spread + (hi / lo - 1.0);
        // `<` keeps the lowest index on ties; a NaN score never wins.
        if score < best.0 {
            best = (score, v as u32);
        }
    }
    best.1
}

/// The inverse stereographic projection from the south pole: the point `w = (X, Y)` of the
/// plane goes to `(2X, 2Y, 1 − |w|²) / (1 + |w|²)` on the unit sphere. The origin goes to the
/// north pole, infinity to the south pole `(0, 0, −1)`. The map is conformal and preserves
/// orientation (counter-clockwise in the plane is counter-clockwise seen from outside).
///
/// [formula: stereographic-projection]
/// [formula: stereographic-conformal]
pub fn inverse_stereographic(w: [f64; 2]) -> [f64; 3] {
    let s = w[0] * w[0] + w[1] * w[1];
    let d = 1.0 + s;
    [2.0 * w[0] / d, 2.0 * w[1] / d, (1.0 - s) / d]
}

/// The stereographic projection from the south pole, the inverse of
/// [`inverse_stereographic`]: `(x, y, z) ↦ (x, y) / (1 + z)` for a point of the unit sphere
/// other than `(0, 0, −1)`.
///
/// [formula: stereographic-projection]
pub fn stereographic(p: [f64; 3]) -> [f64; 2] {
    let d = 1.0 + p[2];
    [p[0] / d, p[1] / d]
}

/// The inversion of the unit sphere with center `c` (`|c| < 1`), Baden et al. 2018, Alg. 1
/// line 7: `x ↦ (1 − |c|²)(x + c)/|x + c|² + c`. It maps the unit sphere to itself, is the
/// identity for `c = 0`, and moves points towards `c`.
///
/// [formula: mobius-sphere-inversion]
pub fn sphere_inversion(c: [f64; 3], x: [f64; 3]) -> [f64; 3] {
    let y = [x[0] + c[0], x[1] + c[1], x[2] + c[2]];
    let y2 = dot(y, y);
    let s = (1.0 - dot(c, c)) / y2;
    [s * y[0] + c[0], s * y[1] + c[1], s * y[2] + c[2]]
}

/// Möbius centering (Baden, Crane, Kazhdan 2018, Algorithm 1): repeatedly computes the center
/// of mass `µ = Σ A_f C_f` of the faces, with `A_f` the areas **on the surface** (normalized to
/// sum 1) and `C_f` the centroid of the image triangle projected to the unit sphere, and applies
/// the inversion with center `c = −J⁻¹ µ`, `J = 2 Σ A_f (I − C_f C_fᵀ)` (the Jacobian of `µ`
/// with respect to `c` at `c = 0`; one Gauss–Newton step), until `|µ| <= tolerance`. A step
/// longer than `max_step` is shortened to it, and a step that does not decrease `|µ|` is halved
/// until it does (the line search the paper suggests for "more challenging models").
/// `positions` must be on the unit sphere; they are changed in place.
///
/// # Errors
/// [`BffError::Geometry`], [`BffError::LengthMismatch`] (one position per vertex),
/// [`BffError::NonFiniteInput`] for options out of their domain or a non-finite position.
///
/// [formula: mobius-centering-jacobian]
/// [formula: mobius-centering-existence]
/// [formula: mobius-sphere-inversion]
pub fn mobius_center(
    mesh: &TriMesh,
    geometry: &Geometry,
    positions: &mut [[f64; 3]],
    options: &CenteringOptions,
) -> Result<CenteringReport, BffError> {
    geometry.check_mesh(mesh)?;
    check_centering_options(options)?;
    check_positions(mesh, positions)?;
    let total: f64 = geometry.face_areas().iter().sum();
    let weights: Vec<f64> = geometry.face_areas().iter().map(|a| a / total).collect();
    let center_of_mass = |positions: &[[f64; 3]], jacobian: Option<&mut [[f64; 3]; 3]>| {
        let mut mu = [0.0; 3];
        let mut j = [[0.0; 3]; 3];
        let want_j = jacobian.is_some();
        for (face, &w) in mesh.faces().iter().zip(&weights) {
            let [a, b, c] = face.map(|v| positions[v as usize]);
            let mut centroid = [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]];
            let norm = dot(centroid, centroid).sqrt();
            if norm > 0.0 {
                centroid = centroid.map(|x| x / norm);
            }
            for i in 0..3 {
                mu[i] += w * centroid[i];
            }
            if want_j {
                for i in 0..3 {
                    for k in 0..3 {
                        let id = if i == k { 1.0 } else { 0.0 };
                        j[i][k] += 2.0 * w * (id - centroid[i] * centroid[k]);
                    }
                }
            }
        }
        if let Some(out) = jacobian {
            *out = j;
        }
        mu
    };
    let mut j = [[0.0; 3]; 3];
    let mut mu = center_of_mass(positions, Some(&mut j));
    let initial_center = dot(mu, mu).sqrt();
    let mut report = CenteringReport {
        status: CenteringStatus::MaxIterations,
        iterations: 0,
        initial_center,
        final_center: initial_center,
        shortened_steps: 0,
        halvings: 0,
    };
    let mut trial = positions.to_vec();
    loop {
        let norm = dot(mu, mu).sqrt();
        report.final_center = norm;
        if norm <= options.tolerance {
            report.status = CenteringStatus::Converged;
            break;
        }
        if report.iterations >= options.max_iterations {
            report.status = CenteringStatus::MaxIterations;
            break;
        }
        let Some(step) = solve3(&j, mu) else {
            report.status = CenteringStatus::SingularJacobian;
            break;
        };
        let mut c = step.map(|x| -x);
        let length = dot(c, c).sqrt();
        if length > options.max_step {
            c = c.map(|x| x * options.max_step / length);
            report.shortened_steps += 1;
        }
        // Backtracking: halve the step until |µ| decreases.
        let mut accepted = false;
        for _ in 0..=MAX_HALVINGS {
            for (t, p) in trial.iter_mut().zip(positions.iter()) {
                *t = sphere_inversion(c, *p);
            }
            let mut j_new = [[0.0; 3]; 3];
            let mu_new = center_of_mass(&trial, Some(&mut j_new));
            if dot(mu_new, mu_new).sqrt() < norm {
                positions.copy_from_slice(&trial);
                mu = mu_new;
                j = j_new;
                accepted = true;
                break;
            }
            c = c.map(|x| 0.5 * x);
            report.halvings += 1;
        }
        if !accepted {
            report.status = CenteringStatus::Stalled;
            break;
        }
        report.iterations += 1;
    }
    Ok(report)
}

/// Largest number of halvings of one centering step before the centering gives up.
const MAX_HALVINGS: usize = 40;

/// The distortion of a map of the mesh to the unit sphere; see [`SphereDistortion`].
///
/// # Errors
/// [`BffError::Geometry`], [`BffError::LengthMismatch`], [`BffError::NonFiniteInput`].
///
/// [formula: quasi-conformal-error]
/// [formula: face-log-scale-factor]
/// [formula: spherical-triangle-area]
pub fn sphere_distortion(
    mesh: &TriMesh,
    geometry: &Geometry,
    positions: &[[f64; 3]],
) -> Result<SphereDistortion, BffError> {
    geometry.check_mesh(mesh)?;
    check_positions(mesh, positions)?;
    let surface_total: f64 = geometry.face_areas().iter().sum();
    let mut image_area = 0.0;
    let mut solid_angle = 0.0;
    for face in mesh.faces() {
        let [a, b, c] = face.map(|v| positions[v as usize]);
        let n = cross(sub(b, a), sub(c, a));
        image_area += 0.5 * dot(n, n).sqrt();
        // Van Oosterom and Strackee (1983): tan(E/2) = det / (1 + a·b + b·c + c·a).
        let det = dot(a, cross(b, c));
        solid_angle += 2.0 * det.atan2(1.0 + dot(a, b) + dot(b, c) + dot(c, a));
    }
    // Lengths in the image are scaled by `scale` so that the total areas agree.
    let scale = (surface_total / image_area).sqrt();
    let mut accumulator = DistortionAccumulator::new();
    for (f, face) in mesh.faces().iter().enumerate() {
        let [a, b, c] = face.map(|v| positions[v as usize]);
        let lengths = [0, 1, 2].map(|k| geometry.edge_lengths()[mesh.edge(3 * f + k)]);
        accumulator.add(
            lengths,
            geometry.face_areas()[f],
            plane_frame(a, b, c, scale),
        );
    }
    let max_radius_error = positions
        .iter()
        .map(|p| (dot(*p, *p).sqrt() - 1.0).abs())
        .fold(0.0, f64::max);
    Ok(SphereDistortion {
        faces: accumulator.finish(),
        image_area,
        covering_degree: solid_angle / (4.0 * PI),
        max_radius_error,
    })
}

/// The image triangle `a, b, c` (in space), scaled by `scale`, laid out in the plane with `a` at
/// the origin and `b` on the positive `x` axis, counter-clockwise if `det(a, b, c) > 0`
/// (counter-clockwise seen from outside the sphere) and clockwise if it is negative. All zero if
/// `a = b`.
fn plane_frame(a: [f64; 3], b: [f64; 3], c: [f64; 3], scale: f64) -> [[f64; 2]; 3] {
    let e = sub(b, a);
    let d = sub(c, a);
    let l = dot(e, e).sqrt();
    if l == 0.0 {
        return [[0.0; 2]; 3];
    }
    let n = cross(e, d);
    let mut y = dot(n, n).sqrt() / l;
    if dot(a, cross(b, c)) < 0.0 {
        y = -y;
    }
    let x = dot(d, e) / l;
    [[0.0, 0.0], [scale * l, 0.0], [scale * x, scale * y]]
}

/// The area-weighted mean over all faces of `ln(a_f / A_f)` (the reference's
/// `Distortion::computeAreaScaling(...)[2]`, the objective of the radius search).
fn log_area_ratio_mean(mesh: &TriMesh, geometry: &Geometry, positions: &[[f64; 3]]) -> f64 {
    let mut total = 0.0;
    let mut weight = 0.0;
    for (f, face) in mesh.faces().iter().enumerate() {
        let [a, b, c] = face.map(|v| positions[v as usize]);
        let n = cross(sub(b, a), sub(c, a));
        let area = geometry.face_areas()[f];
        total += area * (0.5 * dot(n, n).sqrt() / area).ln();
        weight += area;
    }
    total / weight
}

fn check_positions(mesh: &TriMesh, positions: &[[f64; 3]]) -> Result<(), BffError> {
    if positions.len() != mesh.n_vertices() {
        return Err(BffError::LengthMismatch {
            what: "positions on the sphere",
            expected: mesh.n_vertices(),
            got: positions.len(),
        });
    }
    if let Some(index) = positions
        .iter()
        .position(|p| !p.iter().all(|x| x.is_finite()))
    {
        return Err(BffError::NonFiniteInput {
            what: "positions on the sphere",
            index,
            value: positions[index]
                .into_iter()
                .find(|x| !x.is_finite())
                .unwrap_or(f64::NAN),
        });
    }
    Ok(())
}

/// Solves the symmetric 3×3 system `J x = b` by Cramer's rule; `None` if `det J` is zero or not
/// finite.
fn solve3(j: &[[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = dot(j[0], cross(j[1], j[2]));
    if !(det.is_finite() && det != 0.0) {
        return None;
    }
    // Columns of J are its rows (symmetric); x_i = det(J with column i replaced by b) / det.
    let cols = [
        [j[0][0], j[1][0], j[2][0]],
        [j[0][1], j[1][1], j[2][1]],
        [j[0][2], j[1][2], j[2][2]],
    ];
    let x = [
        dot(b, cross(cols[1], cols[2])) / det,
        dot(cols[0], cross(b, cols[2])) / det,
        dot(cols[0], cross(cols[1], b)) / det,
    ];
    x.iter().all(|v| v.is_finite()).then_some(x)
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stereographic_round_trip_and_unit_norm() {
        for w in [[0.0, 0.0], [0.3, -2.0], [1e3, 5.0], [-0.7, 0.1]] {
            let p = inverse_stereographic(w);
            assert!((dot(p, p) - 1.0).abs() < 1e-15);
            let back = stereographic(p);
            // `1 + z = 2 / (1 + |w|²)` is computed by cancellation near the south pole, so the
            // round trip loses a factor `1 + |w|²` of the rounding error.
            let tol = 1e-15 * (1.0 + w[0] * w[0] + w[1] * w[1]);
            assert!((back[0] - w[0]).abs() < tol * (1.0 + w[0].abs()));
            assert!((back[1] - w[1]).abs() < tol * (1.0 + w[1].abs()));
        }
        assert_eq!(inverse_stereographic([0.0, 0.0]), [0.0, 0.0, 1.0]);
    }

    #[test]
    fn inversion_keeps_the_sphere_and_is_identity_at_zero() {
        let x = inverse_stereographic([0.4, -1.3]);
        assert_eq!(sphere_inversion([0.0; 3], x), x);
        for c in [[0.2, -0.1, 0.4], [0.0, 0.0, 0.9], [-0.5, 0.5, 0.1]] {
            let y = sphere_inversion(c, x);
            assert!((dot(y, y) - 1.0).abs() < 1e-14);
        }
    }

    #[test]
    fn solve3_inverts() {
        let j = [[2.0, 0.3, 0.1], [0.3, 1.5, -0.2], [0.1, -0.2, 1.0]];
        let x = solve3(&j, [1.0, 2.0, 3.0]).unwrap();
        for (row, b) in j.iter().zip([1.0, 2.0, 3.0]) {
            assert!((dot(*row, x) - b).abs() < 1e-14);
        }
        assert!(solve3(&[[0.0; 3]; 3], [1.0, 0.0, 0.0]).is_none());
    }
}
