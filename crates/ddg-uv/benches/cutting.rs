//! End-to-end benchmarks of the cut flow on polygon meshes: start, a fixed number of flow steps,
//! description of the cut, insertion, unwrap, and the whole chain.
//!
//! ```text
//! cargo bench -p ddg-uv --bench cutting                  # everything that is available
//! cargo bench -p ddg-uv --bench cutting -- cut/steps     # one stage
//! ```
//!
//! The generated meshes need no files. The model meshes are read from `DDG_MODELS` (default
//! `data/models/spot/spot`) and the dataset meshes from `DDG_DATASET` (default
//! `data/optimal-cones/optimal-cones-data/data`); a missing file is skipped with a notice.
//! Every run uses the default options, one thread, [`STEPS`] steps (the step count is fixed, not
//! a stop rule, so runs are comparable) and the same start on a mesh: the geodesic circle of
//! radius `0.2 · sqrt(area)` around the vertex at a third of the vertex list.

use std::path::PathBuf;
use std::time::Duration;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use ddg_core::io::read_obj;
use ddg_uv::cutting::{
    ConePlacement, CutOnPolygons, CutUnwrapOptions, DEFAULT_CROSSING_OFFSET, FlowKind, Session,
    Start, insert,
};
use ddg_uv::primitives::{grid_cube, torus, uv_sphere};
use ddg_uv::{MeshInput, PolygonMesh, unwrap};
use ddg_variational_cuts::{CutProblem, Energy, FlowOptions, StepSize, length_weight_for_ripples};
use ddg_variational_cuts::{PAPER_FEATURE_SMOOTHING, ScaleWeight};

/// Flow steps per run.
const STEPS: usize = 10;

struct Fixture {
    name: String,
    mesh: PolygonMesh,
    start: Start,
    problem: CutProblem,
}

fn area(mesh: &PolygonMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|c| mesh.positions[mesh.corner_vertices[c as usize] as usize]);
            let (u, v) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
            );
            let x = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
        })
        .sum()
}

fn fixture(name: &str, mesh: PolygonMesh) -> Fixture {
    let start = Start::GeodesicCircle {
        vertex: (mesh.positions.len() / 3) as u32,
        radius: 0.2 * area(&mesh).sqrt(),
    };
    let mut session = Session::new(&mesh.input(), 2, &start).expect("the benchmark meshes load");
    // The length weight from the roughening rate of the start (a ripple wavelength of a quarter
    // of the unit-area surface's size).
    let mut problem = CutProblem::new(Energy::Hencky, 0.0);
    let omega = session
        .roughening_rate(&problem, &Default::default())
        .expect("the start has a cut");
    problem.length_weight = length_weight_for_ripples(omega, 0.25);
    Fixture {
        name: name.to_string(),
        mesh,
        start,
        problem,
    }
}

fn polygons_of(path: PathBuf) -> Option<PolygonMesh> {
    let Ok(obj) = read_obj(&path) else {
        eprintln!("{} not found: skipped", path.display());
        return None;
    };
    let polygons: Vec<Vec<u32>> = (0..obj.n_faces())
        .map(|f| obj.corner_verts[obj.face_corners(f)].to_vec())
        .collect();
    Some(PolygonMesh::from_polygons(obj.positions.clone(), &polygons))
}

fn fixtures() -> Vec<Fixture> {
    let models = std::env::var_os("DDG_MODELS").map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/spot/spot"),
        PathBuf::from,
    );
    let dataset = std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    );
    let mut out = vec![
        fixture("uv_sphere_64x32", uv_sphere(64, 32, 1.0).expect("valid")),
        fixture("grid_cube_16", grid_cube(16, 2.0)),
        fixture("torus_96x24", torus(96, 24, 1.0, 0.25).expect("valid")),
    ];
    for (name, path) in [
        ("spot_quads", models.join("spot_quadrangulated.obj")),
        ("spot_triangles", models.join("spot_triangulated.obj")),
        ("camel", dataset.join("camel.obj")),
        ("bunny", dataset.join("bunny.obj")),
    ] {
        if let Some(mesh) = polygons_of(path) {
            out.push(fixture(name, mesh));
        }
    }
    out
}

fn session(f: &Fixture) -> Session {
    Session::new(&f.mesh.input(), 2, &f.start).expect("session")
}

/// A session after `STEPS` steps.
fn stepped(f: &Fixture, options: &FlowOptions<'_>) -> Session {
    let mut s = session(f);
    for _ in 0..STEPS {
        s.step(&f.problem, options).expect("step");
    }
    s
}

fn group<'a>(
    c: &'a mut Criterion,
    name: &str,
) -> criterion::BenchmarkGroup<'a, criterion::measurement::WallTime> {
    let mut g = c.benchmark_group(name);
    g.sample_size(10).warm_up_time(Duration::from_secs(1));
    g
}

fn benches(c: &mut Criterion) {
    let options = FlowOptions::new(StepSize::default());
    let fixtures = fixtures();
    let off = CutUnwrapOptions {
        cones: ConePlacement::Off,
        ..Default::default()
    };
    let optimal = CutUnwrapOptions::default();

    // Start: validation, the seam cut, the maps, the flow's symbolic analysis, the distance.
    let mut g = group(c, "cut/start");
    for f in &fixtures {
        g.bench_function(&f.name, |b| b.iter(|| session(f)));
    }
    g.finish();

    // A fixed number of steps from the start (the session is built outside the timing).
    let mut g = group(c, "cut/steps");
    for f in &fixtures {
        g.bench_function(format!("{}/{STEPS}", f.name), |b| {
            b.iter_batched(
                || session(f),
                |mut s| {
                    for _ in 0..STEPS {
                        s.step(&f.problem, &options).expect("step");
                    }
                    s
                },
                BatchSize::LargeInput,
            )
        });
    }
    g.finish();

    // The stages after the flow, on the state after STEPS steps.
    let states: Vec<(&Fixture, Session, CutOnPolygons)> = fixtures
        .iter()
        .map(|f| {
            let start_cut = session(f).describe_cut(DEFAULT_CROSSING_OFFSET).expect("cut");
            let s = stepped(f, &options);
            let cut = s.describe_cut(DEFAULT_CROSSING_OFFSET).expect("cut");
            eprintln!(
                "{}: {} polygons, {} triangles; cut at the start {} edge points, after {STEPS} steps {} edge points in {} paths",
                f.name,
                f.mesh.n_polygons(),
                f.mesh.triangles.len(),
                start_cut.edge_points.len(),
                cut.edge_points.len(),
                cut.paths.len()
            );
            (f, s, cut)
        })
        .collect();
    let mut g = group(c, "cut/preview");
    for (f, s, _) in &states {
        g.bench_function(&f.name, |b| {
            b.iter(|| s.preview(DEFAULT_CROSSING_OFFSET).expect("preview"))
        });
    }
    g.finish();
    let mut g = group(c, "cut/describe");
    for (f, s, _) in &states {
        g.bench_function(&f.name, |b| {
            b.iter(|| s.describe_cut(DEFAULT_CROSSING_OFFSET).expect("describe"))
        });
    }
    g.finish();
    let mut g = group(c, "cut/insert");
    for (f, _, cut) in &states {
        g.bench_function(&f.name, |b| {
            b.iter(|| insert(&f.mesh, cut).expect("insert"))
        });
    }
    g.finish();
    // `unwrap` of the inserted mesh with the cut as seams (the insertion itself is timed above).
    let mut g = group(c, "cut/unwrap");
    for (f, _, cut) in &states {
        let inserted = insert(&f.mesh, cut).expect("insert");
        let seams = inserted.seams_after(&[]);
        let all: Vec<u32> = (0..inserted.mesh.positions.len() as u32).collect();
        let allowed = MeshInput {
            seam_edges: &seams,
            ..inserted.mesh.input()
        };
        let forbidden = MeshInput {
            forbidden_vertices: &all,
            ..allowed
        };
        g.bench_function(format!("{}/cones_off", f.name), |b| {
            b.iter(|| unwrap(&forbidden, &off.unwrap).expect("unwrap"))
        });
        g.bench_function(format!("{}/cones_optimal", f.name), |b| {
            b.iter(|| unwrap(&allowed, &optimal.unwrap).expect("unwrap"))
        });
    }
    g.finish();

    // The whole chain.
    let mut g = group(c, "cut/end_to_end");
    for f in &fixtures {
        g.bench_function(format!("{}/{STEPS}_steps_cones_off", f.name), |b| {
            b.iter(|| {
                let s = stepped(f, &options);
                s.unwrap(DEFAULT_CROSSING_OFFSET, &off).expect("unwrap")
            })
        });
    }
    g.finish();
}

// ----- N regions: the atlas of the paper (Fig. 16) -----

/// Regions of the atlas runs: the normals of the mesh clustered into this many.
const ATLAS_REGIONS: usize = 6;

struct AtlasFixture {
    name: String,
    mesh: PolygonMesh,
    problem: CutProblem,
    options: FlowOptions<'static>,
}

fn atlas_session(mesh: &PolygonMesh) -> Session {
    Session::new(
        &mesh.input(),
        ATLAS_REGIONS,
        &Start::NormalClusters {
            max_iterations: 100,
        },
    )
    .expect("the benchmark meshes load")
}

fn atlas_fixture(name: &str, mesh: PolygonMesh) -> AtlasFixture {
    let mut s = atlas_session(&mesh);
    let options = s.recommended_options(StepSize::default());
    let mut problem = s.new_problem(
        Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
            r_eps: 1.0,
            smoothing: PAPER_FEATURE_SMOOTHING,
        }),
        0.0,
    );
    let omega = s
        .roughening_rate(&problem, &options.eval)
        .expect("a start with a cut");
    problem.length_weight = length_weight_for_ripples(omega, 0.5);
    AtlasFixture {
        name: name.to_string(),
        mesh,
        problem,
        options,
    }
}

fn atlas_stepped(f: &AtlasFixture) -> Session {
    let mut s = atlas_session(&f.mesh);
    for _ in 0..STEPS {
        s.step(&f.problem, &f.options).expect("step");
    }
    s
}

/// The same stages for `ATLAS_REGIONS` regions started from clusters of the vertex normals:
/// `FlowOptions::for_regions`, the rescaled Dirichlet energy, the length weight of the ripple
/// rule with wavelength 0.5, [`STEPS`] steps.
fn region_benches(c: &mut Criterion) {
    let models = std::env::var_os("DDG_MODELS").map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/spot/spot"),
        PathBuf::from,
    );
    let mut fixtures = vec![
        atlas_fixture("uv_sphere_64x32", uv_sphere(64, 32, 1.0).expect("valid")),
        atlas_fixture("grid_cube_16", grid_cube(16, 2.0)),
    ];
    if let Some(mesh) = polygons_of(models.join("spot_quadrangulated.obj")) {
        fixtures.push(atlas_fixture("spot_quads", mesh));
    }
    let off = CutUnwrapOptions {
        cones: ConePlacement::Off,
        ..Default::default()
    };

    let mut g = group(c, "regions/start");
    for f in &fixtures {
        g.bench_function(&f.name, |b| b.iter(|| atlas_session(&f.mesh)));
    }
    g.finish();

    let mut g = group(c, "regions/steps");
    for f in &fixtures {
        g.bench_function(format!("{}/{STEPS}", f.name), |b| {
            b.iter_batched(
                || atlas_session(&f.mesh),
                |mut s| {
                    for _ in 0..STEPS {
                        s.step(&f.problem, &f.options).expect("step");
                    }
                    s
                },
                BatchSize::LargeInput,
            )
        });
    }
    g.finish();

    let states: Vec<(&AtlasFixture, Session, CutOnPolygons)> = fixtures
        .iter()
        .map(|f| {
            let s = atlas_stepped(f);
            let cut = s.describe_cut(DEFAULT_CROSSING_OFFSET).expect("cut");
            let inserted = insert(&f.mesh, &cut).expect("insert");
            eprintln!(
                "{}: {} regions, {} polygons; after {STEPS} steps {} edge points, {} junctions, {} paths; {} polygons after insertion",
                f.name,
                ATLAS_REGIONS,
                f.mesh.n_polygons(),
                cut.edge_points.len(),
                cut.junctions.len(),
                cut.paths.len(),
                inserted.mesh.n_polygons()
            );
            (f, s, cut)
        })
        .collect();
    let mut g = group(c, "regions/preview");
    for (f, s, _) in &states {
        g.bench_function(&f.name, |b| {
            b.iter(|| s.preview(DEFAULT_CROSSING_OFFSET).expect("preview"))
        });
    }
    g.finish();
    let mut g = group(c, "regions/describe");
    for (f, s, _) in &states {
        g.bench_function(&f.name, |b| {
            b.iter(|| s.describe_cut(DEFAULT_CROSSING_OFFSET).expect("describe"))
        });
    }
    g.finish();
    let mut g = group(c, "regions/insert");
    for (f, _, cut) in &states {
        g.bench_function(&f.name, |b| {
            b.iter(|| insert(&f.mesh, cut).expect("insert"))
        });
    }
    g.finish();
    let mut g = group(c, "regions/unwrap");
    for (f, _, cut) in &states {
        let inserted = insert(&f.mesh, cut).expect("insert");
        let seams = inserted.seams_after(&[]);
        let all: Vec<u32> = (0..inserted.mesh.positions.len() as u32).collect();
        let allowed = MeshInput {
            seam_edges: &seams,
            ..inserted.mesh.input()
        };
        let forbidden = MeshInput {
            forbidden_vertices: &all,
            ..allowed
        };
        g.bench_function(format!("{}/cones_off", f.name), |b| {
            b.iter(|| unwrap(&forbidden, &off.unwrap).expect("unwrap"))
        });
        g.bench_function(format!("{}/cones_optimal", f.name), |b| {
            b.iter(|| unwrap(&allowed, &ddg_uv::UnwrapOptions::default()).expect("unwrap"))
        });
    }
    g.finish();

    let mut g = group(c, "regions/end_to_end");
    for f in &fixtures {
        g.bench_function(format!("{}/{STEPS}_steps_cones_off", f.name), |b| {
            b.iter(|| {
                let s = atlas_stepped(f);
                s.unwrap(DEFAULT_CROSSING_OFFSET, &off).expect("unwrap")
            })
        });
    }
    g.finish();
}

// ----- the triangulation the flow runs on -----

/// A sheared grid (the largest angle 165 degrees): every triangle is badly shaped.
fn sheared_grid(n: usize, sigma: f64) -> PolygonMesh {
    let h = 1.0 / n as f64;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            positions.push([(i as f64 + sigma * j as f64) * h, j as f64 * h, 0.0]);
        }
    }
    let at = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    PolygonMesh::from_triangles(positions, &faces)
}

/// The chain (start, 10 steps, description, insertion, unwrap with cones off) with the caller's
/// triangulation, with Delaunay flips only and with refinement to 25 degrees.
fn intrinsic_benches(c: &mut Criterion) {
    use ddg_uv::cutting::{IntrinsicOptions, Triangulation};
    let models = std::env::var_os("DDG_MODELS").map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/spot/spot"),
        PathBuf::from,
    );
    let mut meshes: Vec<(String, PolygonMesh)> = vec![
        (
            "uv_sphere_64x32".into(),
            uv_sphere(64, 32, 1.0).expect("valid"),
        ),
        (
            "sheared_165_64".into(),
            sheared_grid(64, 1.0 / 15f64.to_radians().tan()),
        ),
    ];
    if let Some(mesh) = polygons_of(models.join("spot_quadrangulated.obj")) {
        meshes.push(("spot_quads".into(), mesh));
    }
    let off = CutUnwrapOptions {
        cones: ConePlacement::Off,
        ..Default::default()
    };
    let triangulations = [
        ("caller", Triangulation::Caller),
        (
            "flips",
            Triangulation::Intrinsic(IntrinsicOptions::flips_only()),
        ),
        (
            "refined",
            Triangulation::Intrinsic(IntrinsicOptions::refined()),
        ),
    ];
    let mut g = group(c, "intrinsic");
    for (name, mesh) in &meshes {
        let start = Start::GeodesicCircle {
            vertex: (mesh.positions.len() / 3) as u32,
            radius: 0.2 * area(mesh).sqrt(),
        };
        for (label, tri) in &triangulations {
            let make = || {
                Session::with_triangulation(&mesh.input(), 2, &start, FlowKind::Auto, tri)
                    .expect("session")
            };
            let mut probe = make();
            let mut problem = CutProblem::new(Energy::Hencky, 0.0);
            let options = FlowOptions::new(StepSize::default());
            let omega = probe
                .roughening_rate(&problem, &options.eval)
                .expect("a start with a cut");
            problem.length_weight = length_weight_for_ripples(omega, 0.25);
            if let Some(r) = probe.intrinsic_report() {
                eprintln!(
                    "{name}/{label}: {} vertices ({} added), angles {:.1} to {:.1} degrees",
                    r.n_vertices,
                    r.n_added,
                    r.min_angle.to_degrees(),
                    r.max_angle.to_degrees()
                );
            }
            g.bench_function(format!("{name}/{label}/start"), |b| b.iter(&make));
            g.bench_function(format!("{name}/{label}/steps_{STEPS}"), |b| {
                b.iter_batched(
                    &make,
                    |mut s| {
                        for _ in 0..STEPS {
                            s.step(&problem, &options).expect("step");
                        }
                        s
                    },
                    BatchSize::LargeInput,
                )
            });
            let mut stepped = make();
            for _ in 0..STEPS {
                stepped.step(&problem, &options).expect("step");
            }
            g.bench_function(format!("{name}/{label}/describe"), |b| {
                b.iter(|| {
                    stepped
                        .describe_cut(DEFAULT_CROSSING_OFFSET)
                        .expect("describe")
                })
            });
            let cut = stepped
                .describe_cut(DEFAULT_CROSSING_OFFSET)
                .expect("describe");
            g.bench_function(format!("{name}/{label}/insert"), |b| {
                b.iter(|| insert(mesh, &cut).expect("insert"))
            });
            g.bench_function(format!("{name}/{label}/unwrap_cones_off"), |b| {
                b.iter(|| {
                    stepped
                        .unwrap(DEFAULT_CROSSING_OFFSET, &off)
                        .expect("unwrap")
                })
            });
            g.bench_function(format!("{name}/{label}/chain"), |b| {
                b.iter(|| {
                    let mut s = make();
                    for _ in 0..STEPS {
                        s.step(&problem, &options).expect("step");
                    }
                    s.unwrap(DEFAULT_CROSSING_OFFSET, &off).expect("unwrap")
                })
            });
        }
    }
    g.finish();
}

criterion_group!(cutting, benches, region_benches, intrinsic_benches);
criterion_main!(cutting);
