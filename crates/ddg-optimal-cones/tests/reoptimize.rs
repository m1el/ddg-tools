//! `reoptimize_cones`: the optimal angles for given cone positions, against the generic solve of
//! the same problem (`solve_with_boundary` on `reoptimization_problem`), against the generic
//! convex solver, under moving one cone, on disconnected meshes, at the edge cases, and on the
//! published dataset (the re-optimized angles never give more distortion than the published
//! ones at the published positions).

mod common;

use std::f64::consts::PI;
use std::path::PathBuf;

use common::{
    Features, Reference, TestMesh, max_abs_diff, random_mesh, reference_solve_with_boundary,
};
use ddg_core::generators::SplitMix64;
use ddg_core::io::{read_cones, read_obj};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{
    BoundaryTreatment, Cone, ConeError, ConeProblem, Parallelism, PrescribedConeIssue,
    Reoptimization, ReoptimizeMethod, ReoptimizeMethodUsed, ReoptimizeOptions, SolveMethod,
    SolveOptions, Status, certify_with_boundary, evaluate_cones, reoptimization_problem,
    reoptimize_cones, solve_with_boundary,
};

const TREATMENTS: [BoundaryTreatment; 4] = [
    BoundaryTreatment::Fixed,
    BoundaryTreatment::Free,
    BoundaryTreatment::Convex,
    BoundaryTreatment::Penalized,
];

fn options(boundary: BoundaryTreatment) -> ReoptimizeOptions<'static> {
    ReoptimizeOptions {
        boundary,
        ..ReoptimizeOptions::default()
    }
}

/// `k` distinct random interior vertices.
fn random_cones(mesh: &TestMesh, k: usize, rng: &mut SplitMix64) -> Vec<u32> {
    let mut interior = mesh.interior();
    let mut out = Vec::new();
    while out.len() < k && !interior.is_empty() {
        let i = (rng.next_u64() as usize) % interior.len();
        out.push(interior.swap_remove(i) as u32);
    }
    out
}

/// Random steering that keeps the minimizer unique: positive weights, bounds sometimes, a
/// prescribed cone sometimes (never at a given cone).
fn random_steering(mesh: &TestMesh, cones: &[u32], rng: &mut SplitMix64) -> Features {
    let nv = mesh.mesh.n_vertices();
    let mut f = Features::plain(10f64.powf(rng.uniform(-2.0, 0.0)));
    if rng.next_f64() < 0.5 {
        f.distortion_weights = Some((0..nv).map(|_| 2f64.powf(rng.uniform(-2.0, 2.0))).collect());
    }
    if rng.next_f64() < 0.3 {
        let cap = rng.uniform(0.3, 3.0);
        f.lower_bounds = Some(vec![-cap; nv]);
        f.upper_bounds = Some(vec![cap; nv]);
    } else if rng.next_f64() < 0.2 {
        f.lower_bounds = Some(vec![0.0; nv]);
    }
    if rng.next_f64() < 0.3 {
        let free: Vec<usize> = mesh
            .interior()
            .into_iter()
            .filter(|v| !cones.contains(&(*v as u32)))
            .collect();
        if let Some(&v) = free.first() {
            f.prescribed_cones.push(Cone {
                vertex: v as u32,
                angle: rng.uniform(-1.0, 1.0),
            });
        }
    }
    f
}

/// Re-optimizes and checks the result against the generic solve of the same problem, by both
/// methods. `None` if the problem is infeasible.
fn check_against_solve(
    mesh: &TestMesh,
    f: &Features,
    cones: &[u32],
    boundary: BoundaryTreatment,
) -> Option<Reoptimization> {
    let problem = f.problem(mesh);
    let direct = ReoptimizeOptions {
        method: ReoptimizeMethod::DirectKkt,
        ..options(boundary)
    };
    let r = match reoptimize_cones(&problem, cones, &direct) {
        Ok(r) => r,
        Err(ConeError::Infeasible { .. }) => return None,
        Err(e) => panic!("{}: {e}", mesh.name),
    };
    let tag = format!("{} {boundary:?} cones {cones:?}", mesh.name);
    assert_eq!(r.report.status, Status::Converged, "{tag}: {:?}", r.report);
    assert_eq!(r.cones.len(), cones.len());
    for (c, &v) in r.cones.iter().zip(cones) {
        assert_eq!(c.vertex, v);
        assert_eq!(c.angle, r.measure[v as usize]);
    }
    // The certificate is that of the problem with costs 0 and +∞.
    let reopt = reoptimization_problem(&problem, cones).unwrap();
    let cert =
        certify_with_boundary(&reopt, boundary, &r.log_scale, &r.adjoint, &r.measure).unwrap();
    assert_eq!(cert, r.report.certificate);
    // The same problem through the generic solver, both methods.
    for method in [SolveMethod::Auto, SolveMethod::Reference] {
        let s = solve_with_boundary(
            &reopt,
            boundary,
            &SolveOptions {
                method,
                ..SolveOptions::default()
            },
        )
        .unwrap();
        assert_eq!(s.report.status, Status::Converged, "{tag}");
        let dmu = max_abs_diff(&s.measure, &r.measure);
        let du = max_abs_diff(&s.log_scale_optimal, &r.log_scale);
        assert!(
            dmu < 1e-7 && du < 1e-7,
            "{tag} {method:?}: |Δµ| {dmu:e}, |Δu| {du:e}"
        );
        let rel = (s.objective - r.objective).abs() / s.objective.abs().max(1e-9);
        assert!(
            rel < 1e-9,
            "{tag}: objective {} against {}",
            r.objective,
            s.objective
        );
    }
    // The generic route and the automatic choice of reoptimize_cones give the same.
    for method in [ReoptimizeMethod::Solve, ReoptimizeMethod::Auto] {
        let g = reoptimize_cones(
            &problem,
            cones,
            &ReoptimizeOptions {
                method,
                ..options(boundary)
            },
        )
        .unwrap();
        assert_eq!(g.report.status, Status::Converged, "{tag}");
        if method == ReoptimizeMethod::Solve {
            assert_eq!(g.report.method, ReoptimizeMethodUsed::Solve);
        }
        assert!(max_abs_diff(&g.measure, &r.measure) < 1e-7, "{tag}");
    }
    Some(r)
}

#[test]
fn agrees_with_the_generic_solve_on_random_problems() {
    let mut rng = SplitMix64::new(5);
    let mut direct = 0;
    let mut generic = 0;
    for seed in 0..60u64 {
        let mesh = random_mesh(seed as usize, 1000 + seed);
        let k = 1 + (rng.next_u64() % 6) as usize;
        let cones = random_cones(&mesh, k, &mut rng);
        let f = random_steering(&mesh, &cones, &mut rng);
        for boundary in TREATMENTS {
            if mesh.is_closed() && boundary != BoundaryTreatment::Fixed {
                continue;
            }
            if let Some(r) = check_against_solve(&mesh, &f, &cones, boundary) {
                match r.report.method {
                    ReoptimizeMethodUsed::DirectKkt => direct += 1,
                    ReoptimizeMethodUsed::Solve => generic += 1,
                }
            }
        }
    }
    eprintln!("{direct} solved directly, {generic} by the generic solver (a bound active)");
    assert!(
        direct >= 60 && generic >= 5,
        "{direct} direct, {generic} generic"
    );
}

#[test]
fn agrees_with_the_generic_convex_solver() {
    let mut rng = SplitMix64::new(9);
    let mut compared = 0;
    for seed in 0..40u64 {
        let mesh = random_mesh(seed as usize, 2000 + seed);
        let cones = random_cones(&mesh, 1 + (seed % 4) as usize, &mut rng);
        let f = random_steering(&mesh, &cones, &mut rng);
        let boundary = if mesh.is_closed() {
            BoundaryTreatment::Fixed
        } else {
            TREATMENTS[(seed % 3) as usize] // Fixed, Free, Convex
        };
        let problem = f.problem(&mesh);
        let reopt = reoptimization_problem(&problem, &cones).unwrap();
        let g = Features {
            lambda: reopt.lambda,
            distortion_weights: reopt.distortion_weights.clone(),
            cone_costs: reopt.cone_costs.clone(),
            lower_bounds: reopt.lower_bounds.clone(),
            upper_bounds: reopt.upper_bounds.clone(),
            prescribed_cones: reopt.prescribed_cones.clone(),
        };
        match (
            reoptimize_cones(&problem, &cones, &options(boundary)),
            reference_solve_with_boundary(&mesh, &g, boundary),
        ) {
            (
                Ok(r),
                Reference::Solved {
                    objective,
                    log_scale,
                    measure,
                },
            ) => {
                let rel = (r.objective - objective).abs() / objective.abs().max(1e-6);
                assert!(
                    rel < 1e-6,
                    "{}: {} against {objective}",
                    mesh.name,
                    r.objective
                );
                assert!(max_abs_diff(&r.measure, &measure) < 1e-4, "{}", mesh.name);
                assert!(
                    max_abs_diff(&r.log_scale, &log_scale) < 1e-4,
                    "{}",
                    mesh.name
                );
                compared += 1;
            }
            (Err(ConeError::Infeasible { .. }), Reference::Infeasible) => {}
            (_, Reference::Unreliable(why)) => eprintln!("{}: Clarabel {why}", mesh.name),
            (ours, theirs) => panic!("{}: {ours:?} against {theirs:?}", mesh.name),
        }
    }
    eprintln!("{compared} of 40 compared with Clarabel");
    assert!(compared >= 25);
}

#[test]
fn moving_one_cone_changes_the_distortion_as_expected() {
    // With S the cones, S' = S with one cone moved to a neighbour w: adding a free site can
    // only lower the distortion and removing one can only raise it, so
    //   D(S ∪ {w}) ≤ min(D(S), D(S'))   and   max(D(S), D(S')) ≤ D(S \ {v}).
    // On a fine smooth mesh the move changes the distortion by a little only.
    let shape = ddg_core::generators::icosphere(3).unwrap();
    let mut positions = shape.positions.clone();
    for p in &mut positions {
        p[0] *= 1.6;
        p[2] *= 0.7;
    }
    let mesh = TestMesh::from_positions("ellipsoid", positions, shape.faces.clone());
    let problem = Features::plain(0.5).problem(&mesh);
    let o = options(BoundaryTreatment::Fixed);
    let d = |cones: &[u32]| reoptimize_cones(&problem, cones, &o).unwrap().distortion;
    let mut rng = SplitMix64::new(3);
    for _ in 0..10 {
        let cones = random_cones(&mesh, 6, &mut rng);
        let v = cones[0];
        let w = mesh
            .mesh
            .vertex_neighbors(v as usize)
            .map(|x| x as u32)
            .find(|x| !cones.contains(x))
            .unwrap();
        let moved: Vec<u32> = std::iter::once(w)
            .chain(cones[1..].iter().copied())
            .collect();
        let with_both: Vec<u32> = cones.iter().copied().chain([w]).collect();
        let (d_s, d_moved, d_both, d_without) =
            (d(&cones), d(&moved), d(&with_both), d(&cones[1..]));
        let slack = 1e-12;
        assert!(
            d_both <= d_s.min(d_moved) + slack,
            "{d_both} {d_s} {d_moved}"
        );
        assert!(
            d_s.max(d_moved) <= d_without + slack,
            "{d_s} {d_moved} {d_without}"
        );
        let change = (d_moved - d_s).abs() / d_s;
        assert!(
            change < 0.25,
            "a move by one edge changed the distortion by {change}"
        );
        assert!(d_moved != d_s);
    }
}

#[test]
fn edge_cases() {
    let disk = random_mesh(0, 7);
    let sphere = random_mesh(3, 8);
    let torus = random_mesh(4, 9);
    let o = options(BoundaryTreatment::Fixed);
    // No cones on a disk: µ = 0, u solves L u = −Ω.
    let p = Features::plain(0.5).problem(&disk);
    let r = reoptimize_cones(&p, &[], &o).unwrap();
    assert!(r.cones.is_empty() && r.measure.iter().all(|&m| m == 0.0));
    assert_eq!(r.report.status, Status::Converged);
    let e = evaluate_cones(&p, &[], Parallelism::Sequential).unwrap();
    assert!(max_abs_diff(&e.log_scale, &r.log_scale) < 1e-10);
    // No cones on a sphere: Gauss–Bonnet cannot hold.
    let p = Features::plain(0.5).problem(&sphere);
    assert!(matches!(
        reoptimize_cones(&p, &[], &o),
        Err(ConeError::Infeasible {
            component: None,
            ..
        })
    ));
    // No cones on a torus: χ = 0, feasible.
    let p = Features::plain(0.5).problem(&torus);
    assert_eq!(
        reoptimize_cones(&p, &[], &o).unwrap().report.status,
        Status::Converged
    );
    // One cone on a sphere takes all of 4π; three share it.
    let p = Features::plain(0.5).problem(&sphere);
    let r = reoptimize_cones(&p, &[0], &o).unwrap();
    assert!((r.cones[0].angle - 4.0 * PI).abs() < 1e-9);
    let r = reoptimize_cones(&p, &[0, 5, 9], &o).unwrap();
    let sum: f64 = r.cones.iter().map(|c| c.angle).sum();
    assert!((sum - 4.0 * PI).abs() < 1e-9, "{sum}");
    // Three free sites: the automatic choice is the reduced system of the generic solver.
    assert_eq!(r.report.method, ReoptimizeMethodUsed::Solve);
    assert!(
        r.report.components[0]
            .direct_not_used
            .as_ref()
            .unwrap()
            .contains("reduced")
    );
    // Bounds that cannot hold 4π on three cones.
    let capped = Features {
        upper_bounds: Some(vec![1.0; sphere.mesh.n_vertices()]),
        ..Features::plain(0.5)
    };
    assert!(matches!(
        reoptimize_cones(&capped.problem(&sphere), &[0, 5, 9], &o),
        Err(ConeError::Infeasible { .. })
    ));
    // Errors name the cone and the reason.
    let p = Features::plain(0.5).problem(&disk);
    let boundary_vertex = disk.mesh.boundary_vertices()[0];
    let interior = disk.interior()[0] as u32;
    let issue = |cones: &[u32]| match reoptimize_cones(&p, cones, &o) {
        Err(ConeError::InvalidCone { index, reason, .. }) => (index, reason),
        other => panic!("{other:?}"),
    };
    assert_eq!(
        issue(&[interior, boundary_vertex]),
        (1, PrescribedConeIssue::BoundaryVertex)
    );
    assert_eq!(
        issue(&[interior, interior]),
        (1, PrescribedConeIssue::DuplicateVertex)
    );
    assert_eq!(
        issue(&[1_000_000]),
        (0, PrescribedConeIssue::VertexOutOfRange)
    );
    let prescribed = Features {
        prescribed_cones: vec![Cone {
            vertex: interior,
            angle: 0.3,
        }],
        ..Features::plain(0.5)
    };
    assert!(matches!(
        reoptimize_cones(&prescribed.problem(&disk), &[interior], &o),
        Err(ConeError::InvalidCone {
            reason: PrescribedConeIssue::PrescribedVertex,
            ..
        })
    ));
    // A boundary vertex is an error with every treatment.
    for b in TREATMENTS {
        assert!(reoptimize_cones(&p, &[boundary_vertex], &options(b)).is_err());
    }
    // Penalized boundary: an l1 term remains, so the generic solver runs, and says why.
    let forced = |b| ReoptimizeOptions {
        method: ReoptimizeMethod::DirectKkt,
        ..options(b)
    };
    let r = reoptimize_cones(&p, &[interior], &forced(BoundaryTreatment::Penalized)).unwrap();
    assert_eq!(r.report.method, ReoptimizeMethodUsed::Solve);
    assert!(
        r.report.components[0]
            .direct_not_used
            .as_ref()
            .unwrap()
            .contains("l1")
    );
    // An active bound: the generic solver, which respects it.
    let tight = Features {
        upper_bounds: Some(vec![0.05; disk.mesh.n_vertices()]),
        lower_bounds: Some(vec![-0.05; disk.mesh.n_vertices()]),
        ..Features::plain(0.5)
    };
    let all = disk
        .interior()
        .iter()
        .map(|&v| v as u32)
        .collect::<Vec<_>>();
    let r = reoptimize_cones(
        &tight.problem(&disk),
        &all,
        &forced(BoundaryTreatment::Fixed),
    )
    .unwrap();
    assert_eq!(r.report.status, Status::Converged);
    if r.report.method == ReoptimizeMethodUsed::Solve {
        assert!(
            r.report.components[0]
                .direct_not_used
                .as_ref()
                .unwrap()
                .contains("bound")
        );
    }
    assert!(r.cones.iter().all(|c| c.angle.abs() <= 0.05));
}

#[test]
fn disconnected_meshes() {
    let (a, b, c) = (random_mesh(0, 41), random_mesh(3, 42), random_mesh(2, 43));
    let mut positions = Vec::new();
    let mut faces = Vec::new();
    let mut offsets = Vec::new();
    for part in [&a, &b, &c] {
        let offset = positions.len() as u32;
        offsets.push(offset);
        positions.extend_from_slice(&part.positions);
        faces.extend(part.mesh.faces().iter().map(|f| f.map(|v| v + offset)));
    }
    let whole = TestMesh::from_positions("union", positions, faces);
    let mut rng = SplitMix64::new(4);
    let mut cones = Vec::new();
    for (part, &offset) in [&a, &b, &c].iter().zip(&offsets) {
        cones.extend(
            random_cones(part, 3, &mut rng)
                .into_iter()
                .map(|v| v + offset),
        );
    }
    for boundary in [BoundaryTreatment::Fixed, BoundaryTreatment::Free] {
        let r = check_against_solve(&whole, &Features::plain(0.5), &cones, boundary).unwrap();
        assert_eq!(r.report.components.len(), 3);
        // Each part alone gives the same angles (the angles do not depend on the normalization).
        for (k, part) in [&a, &b, &c].iter().enumerate() {
            let local: Vec<u32> = cones[3 * k..3 * k + 3]
                .iter()
                .map(|v| v - offsets[k])
                .collect();
            let alone = reoptimize_cones(
                &Features::plain(0.5).problem(part),
                &local,
                &options(boundary),
            )
            .unwrap();
            for (x, y) in alone.cones.iter().zip(&r.cones[3 * k..3 * k + 3]) {
                assert!((x.angle - y.angle).abs() < 1e-9);
            }
        }
        // The sphere's cones sum to 4π on their own.
        let sum: f64 = r.cones[3..6].iter().map(|c| c.angle).sum();
        assert!((sum - 4.0 * PI).abs() < 1e-9);
    }
    // No cone on the sphere: infeasible, naming it.
    let without: Vec<u32> = cones
        .iter()
        .copied()
        .filter(|&v| !(offsets[1]..offsets[2]).contains(&v))
        .collect();
    assert!(matches!(
        reoptimize_cones(
            &Features::plain(0.5).problem(&whole),
            &without,
            &options(BoundaryTreatment::Fixed)
        ),
        Err(ConeError::Infeasible {
            component: Some(1),
            ..
        })
    ));
}

// ----- the published dataset -----

fn dataset_dir() -> Option<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/optimal-cones/optimal-cones-data/data");
    if dir.is_dir() && !cfg!(debug_assertions) {
        Some(dir)
    } else {
        eprintln!("SKIPPED: the dataset is missing, or this is a debug build");
        None
    }
}

struct Row {
    name: String,
    n_cones: usize,
    closed: bool,
    published: f64,
    reoptimized: f64,
    published_free: Option<f64>,
    reoptimized_free: Option<f64>,
    max_change: f64,
}

fn dataset_row(dir: &std::path::Path, name: &str) -> Row {
    let obj = read_obj(dir.join(format!("{name}.obj"))).unwrap();
    let mesh = TriMesh::new(obj.positions.len(), obj.triangles().unwrap()).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    let published = read_cones(dir.join(format!("{name}-cones.txt")), mesh.n_vertices()).unwrap();
    let problem = ConeProblem::new(&mesh, &geometry, 0.5);
    let vertices: Vec<u32> = published.iter().map(|c| c.vertex).collect();
    let closed = mesh.n_boundary_edges() == 0;

    // The distortion of the published angles at the published positions. On a closed mesh the
    // file's angles sum to 2πχ only to six digits; `evaluate_cones` spreads the defect over the
    // surface, so compare with the angles shifted to sum exactly (the smallest change).
    let mut angles = published.clone();
    if closed {
        let defect: f64 = angles.iter().map(|c| c.angle).sum::<f64>()
            - 2.0 * PI * mesh.euler_characteristic() as f64;
        for c in &mut angles {
            c.angle -= defect / published.len() as f64;
        }
    }
    let evaluated = evaluate_cones(&problem, &angles, Parallelism::Sequential).unwrap();

    let r = reoptimize_cones(&problem, &vertices, &ReoptimizeOptions::default()).unwrap();
    assert_eq!(r.report.status, Status::Converged, "{name}");
    let d = reoptimize_cones(
        &problem,
        &vertices,
        &ReoptimizeOptions {
            method: ReoptimizeMethod::DirectKkt,
            ..ReoptimizeOptions::default()
        },
    )
    .unwrap();
    assert_eq!(d.report.method, ReoptimizeMethodUsed::DirectKkt, "{name}");
    assert!(max_abs_diff(&d.measure, &r.measure) < 1e-7, "{name}");
    let reopt = reoptimization_problem(&problem, &vertices).unwrap();
    let s =
        solve_with_boundary(&reopt, BoundaryTreatment::Fixed, &SolveOptions::default()).unwrap();
    assert_eq!(s.report.status, Status::Converged, "{name}");
    assert!(
        (s.objective - r.objective).abs() <= 1e-9 * r.objective,
        "{name}: {} against {}",
        s.objective,
        r.objective
    );
    // Lower or equal: the published angles are a feasible point of the same problem.
    assert!(
        r.distortion <= evaluated.distortion * (1.0 + 1e-12),
        "{name}: {} against the published {}",
        r.distortion,
        evaluated.distortion
    );
    let max_change = r
        .cones
        .iter()
        .zip(&published)
        .map(|(x, y)| (x.angle - y.angle).abs())
        .fold(0.0, f64::max);

    // With optimized boundary values: the published angles with the boundary optimized for
    // them, against the re-optimized angles with the boundary optimized.
    let (published_free, reoptimized_free) = if closed {
        (None, None)
    } else {
        let fixed = ConeProblem {
            cone_costs: Some(vec![f64::INFINITY; mesh.n_vertices()]),
            prescribed_cones: published.clone(),
            ..problem.clone()
        };
        let p =
            solve_with_boundary(&fixed, BoundaryTreatment::Free, &SolveOptions::default()).unwrap();
        assert_eq!(p.report.status, Status::Converged, "{name}");
        let o = ReoptimizeOptions {
            boundary: BoundaryTreatment::Free,
            ..ReoptimizeOptions::default()
        };
        let q = reoptimize_cones(&problem, &vertices, &o).unwrap();
        assert_eq!(q.report.status, Status::Converged, "{name}");
        assert!(
            q.distortion <= p.distortion_optimal * (1.0 + 1e-12),
            "{name}"
        );
        (Some(p.distortion_optimal), Some(q.distortion))
    };
    Row {
        name: name.to_string(),
        n_cones: published.len(),
        closed,
        published: evaluated.distortion,
        reoptimized: r.distortion,
        published_free,
        reoptimized_free,
        max_change,
    }
}

/// The published cone files: at the published positions, the re-optimized angles give at most
/// the distortion of the published angles (with `u = 0` on the boundary, and with optimized
/// boundary values). Prints the distortion columns of the README table (one thread per mesh; the
/// times of the README come from `examples/reoptimize.rs`, run alone).
#[test]
fn published_cones_reoptimized() {
    let Some(dir) = dataset_dir() else { return };
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            e.file_name()
                .to_str()
                .and_then(|n| n.strip_suffix(".obj"))
                .map(str::to_string)
        })
        .collect();
    names.sort();
    let rows: Vec<Row> = std::thread::scope(|scope| {
        let handles: Vec<_> = names
            .iter()
            .map(|name| scope.spawn(|| dataset_row(&dir, name)))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    println!(
        "| file (cones) | D, published angles | D, re-optimized | ratio | largest angle change | \
         D, published, boundary optimized | D, re-optimized, boundary optimized | ratio |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|");
    for r in &rows {
        let free = match (r.published_free, r.reoptimized_free) {
            (Some(p), Some(q)) => format!("{p:.5} | {q:.5} | {:.4}", q / p),
            _ => "– | – | –".to_string(),
        };
        println!(
            "| {}{} ({}) | {:.5} | {:.5} | {:.4} | {:.1e} | {free} |",
            r.name,
            if r.closed { ", closed" } else { "" },
            r.n_cones,
            r.published,
            r.reoptimized,
            r.reoptimized / r.published,
            r.max_change,
        );
    }
}
