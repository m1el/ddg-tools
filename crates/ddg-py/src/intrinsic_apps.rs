//! `geodesic_distance` and `steiner_tree`: two applications of intrinsic triangulations
//! (Sharp, Soliman, Crane 2019, Sec. 5.3 and 5.1), from the crate `ddg-intrinsic-apps`.
//!
//! The Python glue decodes the buffers into a [`MeshInput`]; the functions here validate it
//! (every message names the parameter), build the mesh, the intrinsic triangulation and the
//! result as plain vectors. No Python is involved, so `cargo test -p ddg-py` covers all of it.
//!
//! * [`geodesic_distance`]: the heat method on the input mesh, its intrinsic Delaunay
//!   triangulation or its Delaunay refinement (default), with the distance at every input
//!   vertex.
//! * [`steiner_tree`]: a short cut through all vertices; its edges as polylines on the surface
//!   (points in space, one polyline per edge of the tree), its length and the check that the
//!   cut leaves a disk.

use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, RefineOptions, RefineRules, SurfacePoint, Triangulation,
};
use ddg_intrinsic_apps::representation::{InputPoints, Representation};
use ddg_intrinsic_apps::heat::{HeatBoundary, HeatOptions, HeatSolver, HeatSource};
use ddg_intrinsic_apps::steiner::{SteinerOptions, cut_to_disk};

use crate::error::{InputError, UnwrapError};

/// A triangle mesh as decoded from the buffers.
#[derive(Debug, Clone, Default)]
pub struct MeshInput {
    /// 3 reals per vertex.
    pub positions: Vec<f64>,
    /// 3 vertex indices per triangle.
    pub triangles: Vec<u32>,
}

/// The validated mesh with its positions.
pub(crate) struct Mesh {
    /// The connectivity.

    pub(crate) mesh: TriMesh,
    /// The geometry.
    pub(crate) geometry: Geometry,
    /// The positions.
    pub(crate) positions: Vec<[f64; 3]>,
}

pub(crate) fn build(input: &MeshInput) -> Result<Mesh, UnwrapError> {
    if !input.positions.len().is_multiple_of(3) {
        return Err(InputError::LengthNotMultiple {
            name: "positions",
            len: input.positions.len(),
            width: 3,
            row: "vertex",
        }
        .into());
    }
    if let Some(position) = input.positions.iter().position(|x| !x.is_finite()) {
        return Err(InputError::NonFinite {
            name: "positions",
            position,
            value: input.positions[position],
        }
        .into());
    }
    let n_vertices = input.positions.len() / 3;
    if !input.triangles.len().is_multiple_of(3) {
        return Err(InputError::LengthNotMultiple {
            name: "triangles",
            len: input.triangles.len(),
            width: 3,
            row: "triangle",
        }
        .into());
    }
    if let Some(position) = input
        .triangles
        .iter()
        .position(|&v| v as usize >= n_vertices)
    {
        return Err(InputError::IndexOutOfRange {
            name: "triangles",
            position,
            value: input.triangles[position],
            count: n_vertices,
            of: "vertices",
        }
        .into());
    }
    let faces: Vec<[u32; 3]> = input
        .triangles
        .chunks_exact(3)
        .map(|t| [t[0], t[1], t[2]])
        .collect();
    let positions: Vec<[f64; 3]> = input
        .positions
        .chunks_exact(3)
        .map(|p| [p[0], p[1], p[2]])
        .collect();
    let mesh = TriMesh::new(n_vertices, faces)
        .map_err(|e| UnwrapError::Mesh(format!("triangles: not a manifold mesh: {e}")))?;
    let geometry = Geometry::from_positions(&mesh, &positions)
        .map_err(|e| UnwrapError::Mesh(format!("positions: {e}")))?;
    Ok(Mesh {
        mesh,
        geometry,
        positions,
    })
}

pub(crate) fn bad_option(name: &'static str, value: String, requirement: &'static str) -> UnwrapError {
    InputError::BadOption {
        name,
        value,
        requirement,
    }
    .into()
}

/// Which triangulation the heat method runs on: the input, its intrinsic Delaunay
/// triangulation or its Delaunay refinement (Ruppert's rules), on either representation.
fn retriangulate<T: Triangulation>(tri: &mut T, triangulation: &str) -> Result<(), UnwrapError> {
    let solve = |e: ddg_intrinsic::IntrinsicError| UnwrapError::Solve(e.to_string());
    match triangulation {
        "input" => {}
        "delaunay" => {
            tri.flip_to_delaunay(&DelaunayOptions::default())
                .map_err(solve)?;
        }
        "refined" => {
            tri.flip_to_delaunay(&DelaunayOptions::default())
                .map_err(solve)?;
            tri.refine_with_rules(&RefineOptions::default(), RefineRules::Ruppert)
                .map_err(solve)?;
        }
        other => {
            return Err(bad_option(
                "triangulation",
                format!("{other:?}"),
                "expected \"input\", \"delaunay\" or \"refined\"",
            ));
        }
    }
    Ok(())
}

/// Options of [`geodesic_distance`].
#[derive(Debug, Clone)]
pub struct HeatArguments {
    /// The source vertices (input ids).
    pub sources: Vec<u32>,
    /// `"input"`, `"delaunay"` or `"refined"`.
    pub triangulation: String,
    /// `m` in `t = m h²`; finite and positive.
    pub time_factor: f64,
    /// `"neumann"`, `"dirichlet"` or `"average"`.
    pub boundary: String,
    /// `"signposts"` (default) or `"integer"`: the representation of the intrinsic
    /// triangulation (Integer Coordinates, WP11).
    pub representation: String,
}

/// The heat-method distance and what was done.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatOutput {
    /// The distance at every input vertex.
    pub distances: Vec<f64>,
    /// Vertices of the triangulation the method ran on.
    pub n_triangulation_vertices: usize,
    /// The time step.
    pub time: f64,
    /// Seconds to build and factor.
    pub seconds_factor: f64,
    /// Seconds of the distance query.
    pub seconds_solve: f64,
}

/// The geodesic distance from `args.sources` by the heat method (paper Sec. 5.3).
///
/// # Errors
/// [`UnwrapError::Input`] for malformed arguments (empty or out-of-range sources, an unknown
/// option, a source on the boundary with Dirichlet conditions), [`UnwrapError::Mesh`] for a mesh
/// that is not a manifold triangle mesh, [`UnwrapError::Solve`] for a failure of the
/// retriangulation or the solves.
pub fn geodesic_distance(input: &MeshInput, args: &HeatArguments) -> Result<HeatOutput, UnwrapError> {
    let boundary = match args.boundary.as_str() {
        "neumann" => HeatBoundary::Neumann,
        "dirichlet" => HeatBoundary::Dirichlet,
        "average" => HeatBoundary::Average,
        other => {
            return Err(bad_option(
                "boundary",
                format!("{other:?}"),
                "expected \"neumann\", \"dirichlet\" or \"average\"",
            ));
        }
    };
    if !(args.time_factor.is_finite() && args.time_factor > 0.0) {
        return Err(bad_option(
            "time_factor",
            args.time_factor.to_string(),
            "must be finite and > 0",
        ));
    }
    let mesh = build(input)?;
    let n = mesh.mesh.n_vertices();
    if args.sources.is_empty() {
        return Err(bad_option("sources", "[]".into(), "at least one source vertex"));
    }
    if let Some(position) = args.sources.iter().position(|&v| v as usize >= n) {
        return Err(InputError::IndexOutOfRange {
            name: "sources",
            position,
            value: args.sources[position],
            count: n,
            of: "vertices",
        }
        .into());
    }
    match parse_representation(&args.representation)? {
        Representation::Signposts => {
            let mut tri = IntrinsicTriangulation::new(&mesh.mesh, &mesh.geometry)
                .map_err(|e| UnwrapError::Mesh(e.to_string()))?;
            retriangulate(&mut tri, &args.triangulation)?;
            heat_on(&tri, args, boundary, n)
        }
        Representation::Integer => {
            let mut tri = IntegerTriangulation::new(&mesh.mesh, &mesh.geometry)
                .map_err(|e| UnwrapError::Mesh(e.to_string()))?;
            retriangulate(&mut tri, &args.triangulation)?;
            heat_on(&tri, args, boundary, n)
        }
    }
}

/// The representation named by the option `representation`: `"signposts"` or `"integer"`.
pub(crate) fn parse_representation(name: &str) -> Result<Representation, UnwrapError> {
    match name {
        "signposts" => Ok(Representation::Signposts),
        "integer" => Ok(Representation::Integer),
        other => Err(bad_option(
            "representation",
            format!("{other:?}"),
            "expected \"signposts\" or \"integer\"",
        )),
    }
}

fn heat_on<T: InputPoints + Clone>(
    tri: &T,
    args: &HeatArguments,
    boundary: HeatBoundary,
    n: usize,
) -> Result<HeatOutput, UnwrapError> {
    let options = HeatOptions {
        time_factor: args.time_factor,
        time: None,
        boundary,
    };
    let solver =
        HeatSolver::new(tri, &options).map_err(|e| UnwrapError::Solve(e.to_string()))?;
    let sources: Vec<HeatSource> = args
        .sources
        .iter()
        .map(|&v| HeatSource::Vertex(v as usize))
        .collect();
    let result = solver.distance(&sources).map_err(|e| match e {
        ddg_intrinsic_apps::heat::HeatError::SourceOnDirichletBoundary { .. } => {
            UnwrapError::Input(InputError::Rejected {
                name: "sources",
                reason: e.to_string(),
            })
        }
        other => UnwrapError::Solve(other.to_string()),
    })?;
    Ok(HeatOutput {
        distances: result.distances[..n].to_vec(),
        n_triangulation_vertices: tri.complex().n_vertices(),
        time: solver.setup().time,
        seconds_factor: solver.setup().seconds,
        seconds_solve: result.report.seconds,
    })
}

/// Options of [`steiner_tree`].
#[derive(Debug, Clone)]
pub struct SteinerArguments {
    /// Over the intrinsic Delaunay triangulation (`true`) or the input edges.
    pub intrinsic: bool,
    /// Insert Fermat points (`false`: the minimal spanning tree).
    pub steiner_points: bool,
    /// `"signposts"` (default) or `"integer"`.
    pub representation: String,
}

/// A tree through all vertices, as polylines in space.
#[derive(Debug, Clone, PartialEq)]
pub struct SteinerOutput {
    /// The points of all polylines, 3 reals each.
    pub points: Vec<f64>,
    /// Where the polyline of every tree edge starts in `points` (in points), with a last entry
    /// equal to the number of points.
    pub offsets: Vec<u32>,
    /// The total length.
    pub length: f64,
    /// The length of the minimal spanning tree of the same triangulation's edges.
    pub mst_length: f64,
    /// Number of Fermat points inserted.
    pub n_steiner_points: usize,
    /// Whether cutting along the tree (and the boundary) leaves a disk.
    pub is_disk: bool,
    /// The largest angle defect of a vertex off the cut (zero when the disk is flat).
    pub max_interior_defect: f64,
}

pub(crate) fn place(mesh: &Mesh, point: SurfacePoint) -> [f64; 3] {
    let p = &mesh.positions;
    match point {
        SurfacePoint::Vertex(v) => p[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = mesh.mesh.edge_vertices(edge as usize);
            [0, 1, 2].map(|k| (1.0 - t) * p[a][k] + t * p[b][k])
        }
        SurfacePoint::Face { face, bary } => {
            let f = mesh.mesh.faces()[face as usize];
            [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * p[f[c] as usize][k]).sum())
        }
    }
}

/// A short cut through all vertices (paper Sec. 5.1).
///
/// # Errors
/// [`UnwrapError::Mesh`] for a mesh that is not a manifold triangle mesh, [`UnwrapError::Solve`]
/// when the triangulation or a trace fails.
pub fn steiner_tree(
    input: &MeshInput,
    args: &SteinerArguments,
) -> Result<SteinerOutput, UnwrapError> {
    let mesh = build(input)?;
    match parse_representation(&args.representation)? {
        Representation::Signposts => {
            let tri = IntrinsicTriangulation::new(&mesh.mesh, &mesh.geometry)
                .map_err(|e| UnwrapError::Mesh(e.to_string()))?;
            steiner_on(&mesh, &tri, args)
        }
        Representation::Integer => {
            let tri = IntegerTriangulation::new(&mesh.mesh, &mesh.geometry)
                .map_err(|e| UnwrapError::Mesh(e.to_string()))?;
            steiner_on(&mesh, &tri, args)
        }
    }
}

fn steiner_on<T: Triangulation + Clone>(
    mesh: &Mesh,
    tri: &T,
    args: &SteinerArguments,
) -> Result<SteinerOutput, UnwrapError> {
    let options = SteinerOptions {
        delaunay: args.intrinsic.then(DelaunayOptions::default),
        steiner_points: args.steiner_points,
        ..SteinerOptions::default()
    };
    let tree = ddg_intrinsic_apps::steiner::steiner_tree(tri, &options)
        .map_err(|e| UnwrapError::Solve(e.to_string()))?;
    let cut = cut_to_disk(&tree.triangulation, &tree.edges)
        .map_err(|e| UnwrapError::Solve(e.to_string()))?;
    let mut points = Vec::new();
    let mut offsets = Vec::with_capacity(tree.traces.len() + 1);
    let mut n_points = 0u32;
    for trace in &tree.traces {
        offsets.push(n_points);
        let [a, b] = trace.vertices.map(|v| tree.triangulation.vertex_location(v as usize));
        let mut push = |p: [f64; 3]| {
            points.extend_from_slice(&p);
            n_points += 1;
        };
        push(place(mesh, a));
        for c in &trace.crossings {
            push(place(mesh, SurfacePoint::Edge { edge: c.edge, t: c.t }));
        }
        push(place(mesh, b));
    }
    offsets.push(n_points);
    Ok(SteinerOutput {
        points,
        offsets,
        length: tree.length,
        mst_length: tree.report.mst_length,
        n_steiner_points: tree.report.n_steiner_points,
        is_disk: cut.is_disk,
        max_interior_defect: cut.max_interior_defect,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators::icosphere;

    fn sphere() -> MeshInput {
        let g = icosphere(2).unwrap();
        MeshInput {
            positions: g.positions.iter().flatten().copied().collect(),
            triangles: g.faces.iter().flatten().copied().collect(),
        }
    }

    #[test]
    fn heat_distance_on_a_sphere() {
        let input = sphere();
        for triangulation in ["input", "delaunay", "refined"] {
            let out = geodesic_distance(
                &input,
                &HeatArguments {
                    sources: vec![0],
                    triangulation: triangulation.into(),
                    time_factor: 1.0,
                    boundary: "neumann".into(),
                    representation: "signposts".into(),
                },
            )
            .unwrap();
            assert_eq!(out.distances.len(), 162);
            assert!(out.distances[0].abs() < 1e-12);
            // The unit sphere: the antipode is at distance π.
            let far = out.distances.iter().fold(0.0f64, |m, &x| m.max(x));
            assert!((far - std::f64::consts::PI).abs() < 0.1, "{triangulation}: {far}");
        }
    }

    #[test]
    fn the_integer_representation_gives_the_same_results() {
        let input = sphere();
        for triangulation in ["input", "delaunay", "refined"] {
            let run = |representation: &str| {
                geodesic_distance(
                    &input,
                    &HeatArguments {
                        sources: vec![0],
                        triangulation: triangulation.into(),
                        time_factor: 1.0,
                        boundary: "neumann".into(),
                        representation: representation.into(),
                    },
                )
                .unwrap()
            };
            let (a, b) = (run("signposts"), run("integer"));
            assert_eq!(a.n_triangulation_vertices, b.n_triangulation_vertices);
            for (x, y) in a.distances.iter().zip(&b.distances) {
                assert!((x - y).abs() < 1e-10, "{triangulation}: {x} {y}");
            }
        }
        let tree = |representation: &str| {
            steiner_tree(
                &input,
                &SteinerArguments {
                    intrinsic: true,
                    steiner_points: true,
                    representation: representation.into(),
                },
            )
            .unwrap()
        };
        let (a, b) = (tree("signposts"), tree("integer"));
        assert_eq!(a.n_steiner_points, b.n_steiner_points);
        assert!((a.length - b.length).abs() < 1e-10);
        assert_eq!(a.is_disk, b.is_disk);
    }

    #[test]
    fn steiner_tree_on_a_sphere_cuts_a_disk() {
        let out = steiner_tree(
            &sphere(),
            &SteinerArguments {
                intrinsic: true,
                steiner_points: true,
                representation: "signposts".into(),
            },
        )
        .unwrap();
        assert!(out.is_disk && out.max_interior_defect < 1e-9);
        assert!(out.length < out.mst_length);
        assert_eq!(*out.offsets.last().unwrap() as usize * 3, out.points.len());
        let along: f64 = out
            .offsets
            .windows(2)
            .map(|w| {
                (w[0] as usize..w[1] as usize - 1)
                    .map(|i| {
                        let (p, q) = (&out.points[3 * i..3 * i + 3], &out.points[3 * i + 3..3 * i + 6]);
                        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
                    })
                    .sum::<f64>()
            })
            .sum();
        assert!((along - out.length).abs() < 1e-9 * out.length);
    }

    #[test]
    fn bad_arguments_name_the_parameter() {
        let input = sphere();
        let args = HeatArguments {
            sources: vec![1000],
            triangulation: "refined".into(),
            time_factor: 1.0,
            boundary: "neumann".into(),
            representation: "signposts".into(),
        };
        let e = geodesic_distance(&input, &args).unwrap_err().to_string();
        assert!(e.starts_with("sources"), "{e}");
        for (args, name) in [
            (
                HeatArguments {
                    sources: vec![0],
                    triangulation: "nope".into(),
                    ..args.clone()
                },
                "triangulation",
            ),
            (
                HeatArguments {
                    sources: vec![0],
                    boundary: "nope".into(),
                    ..args.clone()
                },
                "boundary",
            ),
            (
                HeatArguments {
                    sources: vec![0],
                    time_factor: -1.0,
                    ..args.clone()
                },
                "time_factor",
            ),
            (
                HeatArguments {
                    sources: vec![],
                    ..args.clone()
                },
                "sources",
            ),
        ] {
            let e = geodesic_distance(&input, &args).unwrap_err().to_string();
            assert!(e.starts_with(name), "{e}");
        }
    }
}
