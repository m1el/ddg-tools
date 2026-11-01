//! `geodesic_path`, `geodesic_network` and `geodesic_bezier`: geodesic paths, loops, networks and
//! Bézier curves by edge flips on the surface (Sharp, Crane 2020, *You Can Find Geodesic Paths
//! in Triangle Meshes by Just Flipping Edges*, ACM TOG 39(6)), from the crate
//! `ddg-flip-geodesics`, on either representation of the intrinsic triangulation. Contract:
//! `crates/ddg-flip-geodesics/SPEC.md` §13.
//!
//! The Python glue decodes the buffers into a [`MeshInput`] and a [`GeodesicTask`]; [`run`]
//! validates them (every message names the parameter), builds the triangulation, the initial
//! paths, shortens or subdivides, extracts the paths on the input mesh and returns plain vectors.
//! No Python is involved, so `cargo test -p ddg-py` covers all of it.
//!
//! # The representation
//!
//! [`GeodesicArguments::representation`] is `None` for [`Representation::DEFAULT`] (integer
//! coordinates, since their evaluation passed), `Some("signposts")` or
//! `Some("integer")`. Nothing falls back silently: an extraction that fails on signposts is an
//! error that says to try `representation="integer"`.

use std::collections::HashMap;

use ddg_flip_geodesics::{
    BezierOptions, BezierReport, CertifyOptions, ExtractOptions, FlipError, GeodesicCertificate,
    GeodesicNetwork, InitialPaths, InputCertificate, ShortenOptions, ShortenReport, ShortenStatus,
    geodesic_bezier, network_from_edges, polyline_3d, waypoint_path,
};
use ddg_integer_coords::{InitOptions, IntegerTriangulation};
use ddg_intrinsic::{EdgeTraceOptions, IntrinsicTriangulation, SurfacePoint, Triangulation};
use ddg_intrinsic_apps::representation::Representation;

use crate::error::{InputError, UnwrapError};
use crate::integer_coords::{Places, PolylinesOutput};
use crate::intrinsic_apps::{Mesh, MeshInput, bad_option, build, parse_representation, place};

/// What to compute.
#[derive(Debug, Clone, PartialEq)]
pub enum GeodesicTask {
    /// A path (or loop) through vertices of the input mesh.
    Path {
        /// The vertices, in order.
        vertices: Vec<u32>,
        /// Close the path from the last vertex back to the first.
        closed: bool,
        /// Keep the path through every given vertex; otherwise they only seed it.
        pin_waypoints: bool,
    },
    /// A network of input edges, each given by its two vertices.
    Network {
        /// Pairs of vertices, 2 per edge.
        edges: Vec<u32>,
        /// Vertices that are pinned in addition to the junctions.
        pinned: Vec<u32>,
    },
    /// A geodesic Bézier curve.
    Bezier {
        /// The control vertices.
        control_vertices: Vec<u32>,
        /// Number of subdivision rounds.
        rounds: usize,
    },
}

/// The options common to the three tasks.
#[derive(Debug, Clone, PartialEq)]
pub struct GeodesicArguments {
    /// `None`: [`Representation::DEFAULT`]; or `"signposts"`, `"integer"`.
    pub representation: Option<String>,
    /// Relative intrinsic mollification (integer coordinates only).
    pub mollify: Option<f64>,
    /// A joint is straight when its smaller side exceeds `π − tolerance`.
    pub tolerance: f64,
    /// Iteration limit of the shortening.
    pub max_iterations: Option<usize>,
    /// Stop once the length is below this fraction of the initial length.
    pub stop_length_ratio: Option<f64>,
    /// Also compute the certificate on the input mesh.
    pub certify_input: bool,
    /// Also return the edges of the final triangulation as polylines on the input mesh.
    pub triangulation_edges: bool,
}

impl Default for GeodesicArguments {
    fn default() -> Self {
        GeodesicArguments {
            representation: None,
            mollify: None,
            tolerance: ShortenOptions::default().straight_tolerance,
            max_iterations: None,
            stop_length_ratio: None,
            certify_input: true,
            triangulation_edges: false,
        }
    }
}

/// The report of the shortening (or the sum over the shortenings of a Bézier curve).
#[derive(Debug, Clone, PartialEq)]
pub struct ShorteningReport {
    /// `"geodesic"`, `"stalled"`, `"iteration_limit"`, `"length_ratio"`, `"min_angle"`,
    /// `"sweep_limit"`; for a Bézier curve `"geodesic"` if every shortening was.
    pub status: &'static str,
    /// Joints shortened.
    pub iterations: usize,
    /// Edge flips made.
    pub flips: usize,
    /// Flips the triangulation refused.
    pub refused_flips: usize,
    /// Joints that stalled.
    pub stalled: usize,
    /// Length before.
    pub initial_length: f64,
    /// Length after.
    pub final_length: f64,
}

/// The counts of a Bézier subdivision.
#[derive(Debug, Clone, PartialEq)]
pub struct BezierSummary {
    /// The rounds asked for.
    pub rounds: usize,
    /// The control polygon of the result: vertices of the final triangulation, in path order.
    pub control_points: Vec<u32>,
    /// Their positions in space, 3 reals each.
    pub control_polygon: Vec<f64>,
    /// Edge splits made.
    pub n_splits: usize,
    /// Midpoints snapped to a vertex.
    pub n_snapped: usize,
    /// Splits redone after a Delaunay pass.
    pub n_delaunay_retries: usize,
    /// Midpoints nudged off an input edge.
    pub n_nudged: usize,
    /// Number of shortenings.
    pub n_shortenings: usize,
    /// Shortenings that did not end geodesic.
    pub n_not_geodesic: usize,
}

/// The certificate of the triangulation `T` (SPEC §7 layer 1a).
#[derive(Debug, Clone, PartialEq)]
pub struct TriangulationCertificate {
    /// Every joint side is straight or exempt.
    pub geodesic: bool,
    /// The tolerance used, radians.
    pub tolerance: f64,
    /// Joints.
    pub n_joints: usize,
    /// Pinned joints.
    pub n_pinned: usize,
    /// Straight joints.
    pub n_straight: usize,
    /// Joints blocked by a path end or another path.
    pub n_blocked: usize,
    /// The smallest side, radians.
    pub min_angle: f64,
    /// The failing joints: `(vertex, angle, flexible)`.
    pub failures: Vec<(usize, f64, bool)>,
}

/// The certificate on the input mesh (SPEC §5), as plain numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct InputCertificateOutput {
    /// Whether nothing failed.
    pub passed: bool,
    /// Crossings of input edges checked.
    pub n_crossings: usize,
    /// The largest deviation from `π` at a crossing, radians.
    pub max_crossing_deviation: f64,
    /// Joints checked.
    pub n_joints: usize,
    /// Joints exempt.
    pub n_exempt: usize,
    /// The smallest side over the joints checked (`+∞` if none), radians.
    pub min_joint_angle: f64,
    /// Largest difference between a wedge measured on the input and on `T`, radians.
    pub max_wedge_mismatch: f64,
    /// Largest relative length difference.
    pub max_length_error: f64,
    /// Number of inconsistencies.
    pub n_faults: usize,
    /// The first of them, in words.
    pub faults: Vec<String>,
}

/// The result of a task.
#[derive(Debug, Clone, PartialEq)]
pub struct GeodesicOutput {
    /// `"signposts"` or `"integer"`.
    pub representation: &'static str,
    /// Per path: its points in space, 3 reals each (a closed path repeats its first point).
    pub polylines: Vec<Vec<f64>>,
    /// Per path: the places of its points on the input mesh.
    pub surface_points: Vec<Places>,
    /// Per path: its length.
    pub lengths: Vec<f64>,
    /// Per path: whether it is closed.
    pub closed: Vec<bool>,
    /// What the shortening did.
    pub shortening: ShorteningReport,
    /// The Bézier counts, for a Bézier task.
    pub bezier: Option<BezierSummary>,
    /// The certificate of `T`.
    pub certificate: TriangulationCertificate,
    /// The certificate on the input mesh, if asked for.
    pub input_certificate: Option<InputCertificateOutput>,
    /// Number of input-edge crossings over all polylines.
    pub n_crossings: usize,
    /// The largest `mismatch / length` of the extraction.
    pub max_relative_mismatch: f64,
    /// Edges of the final triangulation, if asked for.
    pub triangulation_edges: Option<PolylinesOutput>,
    /// Vertices of the final triangulation (input vertices first).
    pub n_triangulation_vertices: usize,
}

const STATUSES: [(ShortenStatus, &str); 6] = [
    (ShortenStatus::Geodesic, "geodesic"),
    (ShortenStatus::Stalled, "stalled"),
    (ShortenStatus::IterationLimit, "iteration_limit"),
    (ShortenStatus::LengthRatio, "length_ratio"),
    (ShortenStatus::MinAngle, "min_angle"),
    (ShortenStatus::SweepLimit, "sweep_limit"),
];

fn status_name(s: ShortenStatus) -> &'static str {
    STATUSES
        .iter()
        .find(|(x, _)| *x == s)
        .map_or("unknown", |(_, n)| n)
}

/// Where a [`FlipError`] came from, for its message.
#[derive(Clone, Copy)]
enum Context {
    Path,
    Network,
    Bezier,
}

impl Context {
    fn parameter(self) -> &'static str {
        match self {
            Context::Path => "vertices",
            Context::Network => "edges",
            Context::Bezier => "control_vertices",
        }
    }
}

/// A user-facing error: arguments that cannot work are an input error, the rest a solve error
/// with the advice that fits.
fn flip_error(e: FlipError, context: Context, representation: Representation) -> UnwrapError {
    let text = e.to_string();
    let input = |name: &'static str| {
        UnwrapError::Input(InputError::Rejected {
            name,
            reason: text.clone(),
        })
    };
    match e {
        FlipError::IndexOutOfRange { .. }
        | FlipError::TooFewVertices { .. }
        | FlipError::SameVertex { .. }
        | FlipError::WaypointLeg { .. }
        | FlipError::MaskSize { .. }
        | FlipError::NotAdjacent { .. } => input(context.parameter()),
        FlipError::InvalidOptions { .. } => input("options"),
        FlipError::Crossing { .. } => UnwrapError::Solve(match context {
            Context::Bezier => format!(
                "{text}: the shortest paths between the control points touch or cross; \
                 move a control point"
            ),
            Context::Path => format!(
                "{text}: the shortest paths between the vertices touch or cross; \
                 move a vertex or use fewer vertices"
            ),
            Context::Network => format!("{text}: the edges of the network cross"),
        }),
        FlipError::Bezier { .. } if matches!(context, Context::Bezier) => {
            UnwrapError::Solve(format!(
                "{text} (a control point is repeated, or the curve passes through one twice; move a control point)"
            ))
        }
        FlipError::Extraction { .. } if representation == Representation::Signposts => {
            UnwrapError::Solve(format!(
                "{text}; the signpost representation could not trace this path on the input \
                 mesh, try representation=\"integer\""
            ))
        }
        _ => UnwrapError::Solve(text),
    }
}

/// The edge of the triangulation with ends `a` and `b`.
fn edge_lookup<T: Triangulation>(tri: &T) -> HashMap<(u32, u32), usize> {
    let c = tri.complex();
    (0..c.n_edges())
        .map(|e| {
            let [a, b] = c.edge_vertices(e);
            ((a.min(b) as u32, a.max(b) as u32), e)
        })
        .collect()
}

fn check_indices(name: &'static str, values: &[u32], count: usize) -> Result<(), UnwrapError> {
    match values.iter().position(|&v| v as usize >= count) {
        Some(position) => Err(InputError::IndexOutOfRange {
            name,
            position,
            value: values[position],
            count,
            of: "vertices",
        }
        .into()),
        None => Ok(()),
    }
}

/// The three tasks on the input mesh.
///
/// # Errors
/// [`UnwrapError::Input`] for malformed arguments (every message names the parameter),
/// [`UnwrapError::Mesh`] for a mesh that is not a manifold triangle mesh,
/// [`UnwrapError::Solve`] when the paths cross, the curve cannot be built or the extraction on
/// the input mesh fails.
pub fn run(
    input: &MeshInput,
    task: &GeodesicTask,
    args: &GeodesicArguments,
) -> Result<GeodesicOutput, UnwrapError> {
    let representation = match &args.representation {
        None => Representation::DEFAULT,
        Some(name) => parse_representation(name)?,
    };
    if let Some(m) = args.mollify {
        if !(m.is_finite() && m >= 0.0) {
            return Err(bad_option(
                "mollify",
                m.to_string(),
                "must be finite and >= 0",
            ));
        }
        if representation == Representation::Signposts {
            return Err(bad_option(
                "mollify",
                m.to_string(),
                "only for representation=\"integer\"",
            ));
        }
    }
    let shorten = ShortenOptions {
        straight_tolerance: args.tolerance,
        max_iterations: args.max_iterations,
        stop_length_ratio: args.stop_length_ratio,
        ..ShortenOptions::default()
    };
    if !(args.tolerance.is_finite() && args.tolerance >= 0.0) {
        return Err(bad_option(
            "tolerance",
            args.tolerance.to_string(),
            "must be finite and >= 0",
        ));
    }
    if let Some(r) = args.stop_length_ratio
        && !(r.is_finite() && r > 0.0 && r <= 1.0)
    {
        return Err(bad_option(
            "stop_length_ratio",
            r.to_string(),
            "must be in (0, 1]",
        ));
    }
    let mesh = build(input)?;
    let n = mesh.mesh.n_vertices();
    match task {
        GeodesicTask::Path { vertices, .. } => check_indices("vertices", vertices, n)?,
        GeodesicTask::Network { edges, pinned } => {
            check_indices("edges", edges, n)?;
            check_indices("pinned", pinned, n)?;
        }
        GeodesicTask::Bezier {
            control_vertices, ..
        } => check_indices("control_vertices", control_vertices, n)?,
    }
    let mesh_error = |e: String| UnwrapError::Mesh(e);
    match representation {
        Representation::Signposts => {
            let tri = IntrinsicTriangulation::new(&mesh.mesh, &mesh.geometry)
                .map_err(|e| mesh_error(e.to_string()))?;
            execute(tri, &mesh, task, args, &shorten, representation)
        }
        Representation::Integer => {
            let tri = match args.mollify {
                None => IntegerTriangulation::new(&mesh.mesh, &mesh.geometry)
                    .map_err(|e| mesh_error(e.to_string()))?,
                Some(m) => {
                    IntegerTriangulation::with_options(
                        &mesh.mesh,
                        mesh.geometry.edge_lengths(),
                        &InitOptions { mollify: Some(m) },
                    )
                    .map_err(|e| mesh_error(e.to_string()))?
                    .0
                }
            };
            execute(tri, &mesh, task, args, &shorten, representation)
        }
    }
}

fn execute<T: Triangulation + Sync>(
    tri: T,
    mesh: &Mesh,
    task: &GeodesicTask,
    args: &GeodesicArguments,
    shorten: &ShortenOptions,
    representation: Representation,
) -> Result<GeodesicOutput, UnwrapError> {
    let context = match task {
        GeodesicTask::Path { .. } => Context::Path,
        GeodesicTask::Network { .. } => Context::Network,
        GeodesicTask::Bezier { .. } => Context::Bezier,
    };
    let fail = |e: FlipError| flip_error(e, context, representation);
    let (network, shortening, bezier_report, rounds) = match task {
        GeodesicTask::Path {
            vertices,
            closed,
            pin_waypoints,
        } => {
            let points: Vec<usize> = vertices.iter().map(|&v| v as usize).collect();
            let initial = waypoint_path(&tri, &points, *closed, *pin_waypoints).map_err(fail)?;
            let (network, report) = shortened(tri, &initial, shorten).map_err(fail)?;
            (network, summary(&report), None, 0)
        }
        GeodesicTask::Network { edges, pinned } => {
            let lookup = edge_lookup(&tri);
            let mut mask = vec![false; tri.complex().n_edges()];
            for (i, pair) in edges.chunks_exact(2).enumerate() {
                let key = (pair[0].min(pair[1]), pair[0].max(pair[1]));
                let e = lookup.get(&key).ok_or_else(|| {
                    UnwrapError::Input(InputError::Rejected {
                        name: "edges",
                        reason: format!(
                            "edge {i} joins vertices {} and {}, which are not joined by an edge \
                             of the mesh",
                            pair[0], pair[1]
                        ),
                    })
                })?;
                mask[*e] = true;
            }
            if mask.iter().all(|&m| !m) {
                return Err(bad_option("edges", "[]".into(), "at least one edge"));
            }
            let pinned: Vec<usize> = pinned.iter().map(|&v| v as usize).collect();
            let initial = network_from_edges(&tri, &mask, &pinned).map_err(fail)?;
            let (network, report) = shortened(tri, &initial, shorten).map_err(fail)?;
            (network, summary(&report), None, 0)
        }
        GeodesicTask::Bezier {
            control_vertices,
            rounds,
        } => {
            let controls: Vec<usize> = control_vertices.iter().map(|&v| v as usize).collect();
            let options = BezierOptions {
                shorten: *shorten,
                ..BezierOptions::default()
            };
            let (network, report) =
                geodesic_bezier(tri, &controls, *rounds, &options).map_err(fail)?;
            let shortening = bezier_shortening(&report);
            (network, shortening, Some(report), *rounds)
        }
    };
    finish(
        network,
        mesh,
        args,
        shortening,
        bezier_report,
        rounds,
        context,
        representation,
    )
}

fn shortened<T: Triangulation>(
    tri: T,
    initial: &InitialPaths,
    options: &ShortenOptions,
) -> Result<(GeodesicNetwork<T>, ShortenReport), FlipError> {
    let mut network = GeodesicNetwork::new(tri, &initial.paths, &initial.pinned, None)?;
    let report = network.shorten(options)?;
    Ok((network, report))
}

fn summary(r: &ShortenReport) -> ShorteningReport {
    ShorteningReport {
        status: status_name(r.status),
        iterations: r.iterations,
        flips: r.flips,
        refused_flips: r.refused_flips,
        stalled: r.stalled,
        initial_length: r.initial_length,
        final_length: r.final_length,
    }
}

fn bezier_shortening(r: &BezierReport) -> ShorteningReport {
    let status = if r.n_not_geodesic == 0 {
        "geodesic"
    } else {
        r.other_statuses
            .iter()
            .copied()
            .find(|&s| s != ShortenStatus::Geodesic)
            .map_or("stalled", status_name)
    };
    ShorteningReport {
        status,
        iterations: r.iterations,
        flips: r.flips,
        refused_flips: r.refused_flips,
        stalled: r.stalled,
        initial_length: r.initial_length,
        final_length: r.final_length,
    }
}

#[allow(clippy::too_many_arguments)]
fn finish<T: Triangulation + Sync>(
    network: GeodesicNetwork<T>,
    mesh: &Mesh,
    args: &GeodesicArguments,
    shortening: ShorteningReport,
    bezier_report: Option<BezierReport>,
    rounds: usize,
    context: Context,
    representation: Representation,
) -> Result<GeodesicOutput, UnwrapError> {
    let fail = |e: FlipError| flip_error(e, context, representation);
    let certify_options = CertifyOptions {
        tolerance: CertifyOptions::default().tolerance.max(args.tolerance),
    };
    let certificate = network.certify(&certify_options).map_err(fail)?;
    let extract_options = ExtractOptions {
        certify_input: args.certify_input,
        ..ExtractOptions::default()
    };
    let extraction = network
        .extract_with_report(&extract_options)
        .map_err(fail)?;
    let mut polylines = Vec::new();
    let mut surface_points = Vec::new();
    for path in &extraction.paths {
        let line = polyline_3d(&mesh.mesh, &mesh.positions, path).map_err(fail)?;
        polylines.push(line.into_iter().flatten().collect());
        let mut places = Places::default();
        for &p in &path.points {
            places.push(p);
        }
        surface_points.push(places);
    }
    let n_paths = network.network().n_paths();
    let lengths = (0..n_paths).map(|p| network.path_length(p)).collect();
    let closed = (0..n_paths)
        .map(|p| network.network().is_closed(p))
        .collect();
    let tri = network.triangulation();
    let bezier = bezier_report.map(|r| {
        let control_points: Vec<u32> = r.control_points.iter().map(|&v| v as u32).collect();
        let control_polygon = r
            .control_points
            .iter()
            .flat_map(|&v| place(mesh, tri.vertex_location(v)))
            .collect();
        BezierSummary {
            rounds,
            control_points,
            control_polygon,
            n_splits: r.n_splits,
            n_snapped: r.n_snapped,
            n_delaunay_retries: r.n_delaunay_retries,
            n_nudged: r.n_nudged,
            n_shortenings: r.n_shortenings,
            n_not_geodesic: r.n_not_geodesic,
        }
    });
    let triangulation_edges = if args.triangulation_edges {
        Some(triangulation_edge_polylines(tri, mesh).map_err(fail_solve)?)
    } else {
        None
    };
    Ok(GeodesicOutput {
        representation: match representation {
            Representation::Signposts => "signposts",
            Representation::Integer => "integer",
        },
        polylines,
        surface_points,
        lengths,
        closed,
        shortening,
        bezier,
        certificate: TriangulationCertificate::from(&certificate),
        input_certificate: extraction
            .certificate
            .as_ref()
            .map(InputCertificateOutput::from),
        n_crossings: extraction.report.n_crossings,
        max_relative_mismatch: extraction.report.max_relative_mismatch,
        triangulation_edges,
        n_triangulation_vertices: tri.complex().n_vertices(),
    })
}

fn fail_solve(e: ddg_intrinsic::IntrinsicError) -> UnwrapError {
    UnwrapError::Solve(format!("tracing the edges of the triangulation: {e}"))
}

impl From<&GeodesicCertificate> for TriangulationCertificate {
    fn from(c: &GeodesicCertificate) -> Self {
        TriangulationCertificate {
            geodesic: c.geodesic,
            tolerance: c.tolerance,
            n_joints: c.n_joints,
            n_pinned: c.n_pinned,
            n_straight: c.n_straight,
            n_blocked: c.n_blocked,
            min_angle: c.min_angle,
            failures: c
                .failures
                .iter()
                .map(|f| (f.vertex, f.angle, f.flexible))
                .collect(),
        }
    }
}

impl From<&InputCertificate> for InputCertificateOutput {
    fn from(c: &InputCertificate) -> Self {
        InputCertificateOutput {
            passed: c.passed,
            n_crossings: c.n_crossings,
            max_crossing_deviation: c.max_crossing_deviation,
            n_joints: c.n_joints,
            n_exempt: c.n_exempt,
            min_joint_angle: c.min_joint_angle,
            max_wedge_mismatch: c.max_wedge_mismatch,
            max_length_error: c.max_length_error,
            n_faults: c.n_faults,
            faults: c.faults.clone(),
        }
    }
}

/// Every edge of `tri` as a polyline on the input mesh (the drawing of "Show triangulation").
fn triangulation_edge_polylines<T: Triangulation>(
    tri: &T,
    mesh: &Mesh,
) -> Result<PolylinesOutput, ddg_intrinsic::IntrinsicError> {
    let n_edges = tri.complex().n_edges();
    let options = EdgeTraceOptions::default();
    let mut out = PolylinesOutput {
        edges: (0..n_edges as u32).collect(),
        points: Vec::new(),
        offsets: Vec::with_capacity(n_edges + 1),
        along_input_edge: Vec::with_capacity(n_edges),
        n_crossings: Vec::with_capacity(n_edges),
        max_relative_mismatch: 0.0,
    };
    let mut n_points = 0u32;
    for e in 0..n_edges {
        let trace = tri.trace_edge_on_input(e, &options)?;
        out.offsets.push(n_points);
        let mut push = |p: [f64; 3]| {
            out.points.extend_from_slice(&p);
            n_points += 1;
        };
        let [a, b] = trace.vertices.map(|v| tri.vertex_location(v as usize));
        push(place(mesh, a));
        for c in &trace.crossings {
            push(place(
                mesh,
                SurfacePoint::Edge {
                    edge: c.edge,
                    t: c.t,
                },
            ));
        }
        push(place(mesh, b));
        out.along_input_edge
            .push(trace.along_input_edge.map_or(-1, i64::from));
        out.n_crossings.push(trace.crossings.len() as u32);
        let length = tri.edge_length(e);
        if length > 0.0 {
            out.max_relative_mismatch = out.max_relative_mismatch.max(trace.mismatch / length);
        }
    }
    out.offsets.push(n_points);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ddg_core::generators::{GeneratedMesh, grid_disk, icosphere};

    fn input(g: &GeneratedMesh) -> MeshInput {
        MeshInput {
            positions: g.positions.iter().flatten().copied().collect(),
            triangles: g.faces.iter().flatten().copied().collect(),
        }
    }

    fn distance(g: &GeneratedMesh, a: usize, b: usize) -> f64 {
        (0..3)
            .map(|k| (g.positions[a][k] - g.positions[b][k]).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    fn args(representation: &str) -> GeodesicArguments {
        GeodesicArguments {
            representation: Some(representation.into()),
            ..GeodesicArguments::default()
        }
    }

    const REPRESENTATIONS: [&str; 2] = ["signposts", "integer"];

    #[test]
    fn a_flat_path_is_the_straight_line() {
        let g = grid_disk(6, 6).unwrap();
        for rep in REPRESENTATIONS {
            let task = GeodesicTask::Path {
                vertices: vec![0, 48],
                closed: false,
                pin_waypoints: true,
            };
            let out = run(&input(&g), &task, &args(rep)).unwrap();
            assert_eq!(out.representation, rep);
            assert!((out.lengths[0] - distance(&g, 0, 48)).abs() < 1e-9, "{rep}");
            assert!(out.certificate.geodesic && out.input_certificate.as_ref().unwrap().passed);
            assert_eq!(out.polylines[0].len() % 3, 0);
            assert_eq!(out.polylines[0].len() / 3, out.surface_points[0].kind.len());
        }
    }

    #[test]
    fn the_default_representation_is_the_constant() {
        let g = grid_disk(3, 3).unwrap();
        let task = GeodesicTask::Path {
            vertices: vec![0, 15],
            closed: false,
            pin_waypoints: true,
        };
        let out = run(&input(&g), &task, &GeodesicArguments::default()).unwrap();
        assert_eq!(
            out.representation,
            match Representation::DEFAULT {
                Representation::Signposts => "signposts",
                Representation::Integer => "integer",
            }
        );
    }

    #[test]
    fn a_bezier_curve_on_a_plane_has_a_length_and_a_polygon() {
        let g = grid_disk(8, 8).unwrap();
        for rep in REPRESENTATIONS {
            let task = GeodesicTask::Bezier {
                control_vertices: vec![0, 4, 40, 80],
                rounds: 3,
            };
            let out = run(&input(&g), &task, &args(rep)).unwrap();
            let b = out.bezier.unwrap();
            assert_eq!(b.control_polygon.len(), 3 * b.control_points.len());
            assert!(b.control_points.len() > 2usize.pow(3));
            assert!(out.certificate.geodesic, "{rep}");
        }
    }

    #[test]
    fn a_loop_and_a_network() {
        let g = icosphere(2).unwrap();
        for rep in REPRESENTATIONS {
            let task = GeodesicTask::Path {
                vertices: vec![0, 5, 20],
                closed: true,
                pin_waypoints: true,
            };
            let out = run(&input(&g), &task, &args(rep)).unwrap();
            assert!(out.closed[0]);
            let p = &out.polylines[0];
            assert_eq!(p[..3], p[p.len() - 3..]);
        }
        let flat = grid_disk(6, 6).unwrap();
        for rep in REPRESENTATIONS {
            let task = GeodesicTask::Network {
                edges: vec![0, 1, 1, 2, 1, 8],
                pinned: vec![],
            };
            let out = run(&input(&flat), &task, &args(rep)).unwrap();
            assert_eq!(out.polylines.len(), 3, "{rep}");
            assert!((out.lengths.iter().sum::<f64>() - 3.0 * distance(&flat, 0, 1)).abs() < 1e-9);
        }
    }

    #[test]
    fn errors_name_the_parameter() {
        let g = grid_disk(3, 3).unwrap();
        let path = |vertices: Vec<u32>| GeodesicTask::Path {
            vertices,
            closed: false,
            pin_waypoints: true,
        };
        let e = run(&input(&g), &path(vec![0]), &args("signposts")).unwrap_err();
        assert!(e.to_string().contains("vertices"), "{e}");
        let e = run(&input(&g), &path(vec![0, 0]), &args("integer")).unwrap_err();
        assert!(e.to_string().contains("vertices"), "{e}");
        let e = run(&input(&g), &path(vec![0, 999]), &args("integer")).unwrap_err();
        assert!(e.to_string().contains("vertices"), "{e}");
        let e = run(
            &input(&g),
            &GeodesicTask::Network {
                edges: vec![0, 15],
                pinned: vec![],
            },
            &args("integer"),
        )
        .unwrap_err();
        assert!(e.to_string().contains("edges"), "{e}");
        let mut a = args("signposts");
        a.mollify = Some(1e-4);
        assert!(run(&input(&g), &path(vec![0, 15]), &a).is_err());
        let mut a = args("nope");
        a.mollify = None;
        assert!(run(&input(&g), &path(vec![0, 15]), &a).is_err());
    }

    #[test]
    fn the_triangulation_edges_are_drawn_on_the_input() {
        let g = grid_disk(4, 4).unwrap();
        let task = GeodesicTask::Path {
            vertices: vec![0, 24],
            closed: false,
            pin_waypoints: true,
        };
        for rep in REPRESENTATIONS {
            let mut a = args(rep);
            a.triangulation_edges = true;
            let out = run(&input(&g), &task, &a).unwrap();
            let t = out.triangulation_edges.unwrap();
            assert_eq!(t.offsets.len(), t.edges.len() + 1);
        }
    }
}
