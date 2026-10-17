//! The measurements of paper Sec. 5.1 to 5.3 on one mesh: the tables of this crate's README.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-apps --example apps_study -- TASK MESH
//!   TASK   steiner | fem | heat
//!   MESH   an OBJ file, or a generated mesh: sheared_torus (60 x 16, shift 5),
//!          sheared_square (16 x 16 lattice sheared by 4 cells per row, flat, with boundary),
//!          icosphere_jitter (level 3, vertices moved by 35 % of an edge)
//! ```
//!
//! One `row |` line per result, fields separated by `|`:
//!
//! * `steiner`: vertices, Euler characteristic; lengths of the extrinsic MST (input edges), the
//!   extrinsic Steiner tree (input edges and Fermat points in input faces), the intrinsic MST
//!   (edges of the iDT) and the intrinsic Steiner tree; the last over the first; the number of
//!   Fermat points; whether the intrinsic cut is a flat disk; seconds of the intrinsic tree
//!   (flips, tree, insertions, traces).
//! * `fem`: per triangulation (input, iDT, the regular 4-1 subdivision at one and two levels,
//!   and at each of those vertex counts iDR to 25° with a circumradius bound and iODT with a
//!   target length, both found by bisection): vertices, smallest and largest angle, and the
//!   condition numbers `λ_max / λ_min` of `L x = λ M x` for the Galerkin mass, the lumped mass
//!   and `M = I`.
//! * `heat`: per triangulation (input, iDT, iDR to 25°, iDR to 25° with circumradius at most
//!   half the mean edge of the iDT): vertices, seconds to factor, seconds per query, and the
//!   errors against the exact polyhedral distance at the input vertices for three sources
//!   (vertices `0`, `n/3`, `2n/3`): the mean relative error `|d − d*| / d*` (vertices with
//!   `d* > 0`), the mean absolute error over the largest distance, the largest absolute error
//!   over the largest distance. Timings are single-thread; run one process at a time when they
//!   matter.

use std::f64::consts::TAU;
use std::time::Instant;

use ddg_core::generators::{SplitMix64, icosphere};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayOptions, InsertionLimit, IntrinsicTriangulation, OptimalDelaunayOptions, RefineOptions,
    TargetLength,
};
use ddg_intrinsic_apps::exact::{ExactOptions, exact_distance};
use ddg_intrinsic_apps::fem::{
    EigenOptions, MassMatrix, condition_number, galerkin_mass_matrix, lumped_mass,
    regular_subdivision, stiffness_matrix,
};
use ddg_intrinsic_apps::heat::{HeatOptions, HeatSolver, HeatSource};
use ddg_intrinsic_apps::steiner::{SteinerOptions, cut_to_disk, steiner_tree};

type Soup = (Vec<[f64; 3]>, Vec<[u32; 3]>);

fn sheared_torus(nu: usize, nv: usize, big: f64, small: f64, shift: usize) -> Soup {
    let id = |i: usize, j: usize| ((j % nv) * nu + (i % nu)) as u32;
    let mut positions = Vec::new();
    for j in 0..nv {
        let v = TAU * j as f64 / nv as f64;
        for i in 0..nu {
            let u = TAU * i as f64 / nu as f64;
            let r = big + small * v.cos();
            positions.push([r * u.cos(), r * u.sin(), small * v.sin()]);
        }
    }
    let mut faces = Vec::new();
    for j in 0..nv {
        for i in 0..nu {
            let (a, b) = (id(i, j), id(i + 1, j));
            let (c, d) = (id(i + 1 + shift, j + 1), id(i + shift, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (positions, faces)
}

fn sheared_square(n: usize, sigma: f64) -> Soup {
    let h = 1.0 / n as f64;
    let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            positions.push([(i as f64 + sigma * j as f64) * h, j as f64 * h, 0.0]);
        }
    }
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j + 1));
            faces.push([a, b, d]);
            faces.push([a, d, c]);
        }
    }
    (positions, faces)
}

fn load(name: &str) -> Result<(TriMesh, Geometry), String> {
    let (positions, faces) = match name {
        "sheared_torus" => sheared_torus(60, 16, 1.0, 0.45, 5),
        "sheared_square" => sheared_square(16, 4.0),
        "icosphere_jitter" => {
            let g = icosphere(3).map_err(|e| e.to_string())?;
            let mut rng = SplitMix64::new(3);
            let h = 0.35 * 0.08;
            let positions = g
                .positions
                .iter()
                .map(|p| p.map(|x| x + h * rng.uniform(-1.0, 1.0)))
                .collect();
            (positions, g.faces)
        }
        path => {
            let obj = read_obj(path).map_err(|e| e.to_string())?;
            let mut faces = Vec::new();
            for f in 0..obj.n_faces() {
                let c = &obj.corner_verts[obj.face_corners(f)];
                for k in 1..c.len().saturating_sub(1) {
                    faces.push([c[0], c[k], c[k + 1]]);
                }
            }
            (obj.positions, faces)
        }
    };
    let mesh = TriMesh::new(positions.len(), faces).map_err(|e| e.to_string())?;
    let geometry = Geometry::from_positions(&mesh, &positions).map_err(|e| e.to_string())?;
    Ok((mesh, geometry))
}

fn short_name(name: &str) -> String {
    std::path::Path::new(name)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(name)
        .to_string()
}

fn idt(base: &IntrinsicTriangulation) -> Result<IntrinsicTriangulation, String> {
    let mut tri = base.clone();
    tri.flip_to_delaunay(&DelaunayOptions::default())
        .map_err(|e| e.to_string())?;
    Ok(tri)
}

fn steiner(name: &str, mesh: &TriMesh, geometry: &Geometry) -> Result<(), String> {
    let tri = IntrinsicTriangulation::new(mesh, geometry).map_err(|e| e.to_string())?;
    let run = |delaunay: bool, points: bool| {
        steiner_tree(
            &tri,
            &SteinerOptions {
                delaunay: delaunay.then(DelaunayOptions::default),
                steiner_points: points,
                ..SteinerOptions::default()
            },
        )
        .map_err(|e| e.to_string())
    };
    let ext_mst = run(false, false)?;
    let ext_steiner = run(false, true)?;
    let int_mst = run(true, false)?;
    let start = Instant::now();
    let int_steiner = run(true, true)?;
    let seconds = start.elapsed().as_secs_f64();
    let cut = cut_to_disk(&int_steiner.triangulation, &int_steiner.edges).map_err(|e| e.to_string())?;
    let flat_disk = cut.is_disk && cut.max_interior_defect < 1e-9;
    println!(
        "row | steiner | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} | {:.3}",
        short_name(name),
        mesh.n_vertices(),
        mesh.euler_characteristic(),
        ext_mst.length,
        ext_steiner.length,
        int_mst.length,
        int_steiner.length,
        int_steiner.length / ext_mst.length,
        int_steiner.report.n_steiner_points,
        flat_disk,
        seconds
    );
    Ok(())
}

fn conditions(tri: &IntrinsicTriangulation) -> Result<[f64; 3], String> {
    let l = stiffness_matrix(tri);
    let options = EigenOptions::default();
    let mut out = [0.0; 3];
    for (k, mass) in [
        MassMatrix::Galerkin(galerkin_mass_matrix(tri)),
        MassMatrix::Lumped(lumped_mass(tri)),
        MassMatrix::Identity,
    ]
    .into_iter()
    .enumerate()
    {
        let r = condition_number(&l, &mass, &options).map_err(|e| e.to_string())?;
        if !(r.lambda_max.converged && r.lambda_min.converged) {
            eprintln!("note: not converged: {r:?}");
        }
        out[k] = r.condition_number;
    }
    Ok(out)
}

fn fem_row(name: &str, kind: &str, tri: &IntrinsicTriangulation) -> Result<(), String> {
    let q = tri.quality(0.0);
    let [g, l, i] = conditions(tri)?;
    println!(
        "row | fem | {} | {kind} | {} | {:.2} | {:.2} | {g:.3e} | {l:.3e} | {i:.3e}",
        short_name(name),
        tri.n_vertices(),
        q.min_angle.to_degrees(),
        q.max_angle.to_degrees()
    );
    Ok(())
}

/// Bisection on a positive parameter `x` so that `run(x)` has about `goal` vertices; the
/// vertex count is assumed to decrease with `x`.
fn matched(
    goal: usize,
    (mut lo, mut hi): (f64, f64),
    run: impl Fn(f64) -> Result<IntrinsicTriangulation, String>,
) -> Result<IntrinsicTriangulation, String> {
    let mut best: Option<(usize, IntrinsicTriangulation)> = None;
    for _ in 0..16 {
        let x = (lo * hi).sqrt();
        let tri = run(x)?;
        let count = tri.n_vertices();
        let miss = count.abs_diff(goal);
        if count > goal {
            lo = x;
        } else {
            hi = x;
        }
        let done = miss * 50 < goal;
        if best.as_ref().is_none_or(|b| miss < b.0) {
            best = Some((miss, tri));
        }
        if done {
            break;
        }
    }
    Ok(best.ok_or("no run")?.1)
}

fn fem(name: &str, mesh: &TriMesh, geometry: &Geometry) -> Result<(), String> {
    let input = IntrinsicTriangulation::new(mesh, geometry).map_err(|e| e.to_string())?;
    fem_row(name, "input", &input)?;
    let delaunay = idt(&input)?;
    fem_row(name, "idt", &delaunay)?;
    let mean = delaunay.edge_lengths().iter().sum::<f64>() / delaunay.n_edges() as f64;
    for level in 1..=2 {
        let (m, g) = regular_subdivision(mesh, geometry, level).map_err(|e| e.to_string())?;
        let sub = IntrinsicTriangulation::new(&m, &g).map_err(|e| e.to_string())?;
        fem_row(name, &format!("4-1 x{level}"), &sub)?;
        let goal = sub.n_vertices();
        let idr = matched(goal, (0.01 * mean, 10.0 * mean), |r| {
            let mut tri = delaunay.clone();
            tri.refine(&RefineOptions {
                max_circumradius: Some(r),
                max_insertions: InsertionLimit::Vertices(4 * goal),
                ..RefineOptions::default()
            })
            .map_err(|e| e.to_string())?;
            Ok(tri)
        })?;
        fem_row(name, &format!("iDR @4-1 x{level}"), &idr)?;
        let iodt = matched(goal, (0.01 * mean, 10.0 * mean), |h| {
            let mut tri = delaunay.clone();
            tri.optimal_delaunay(&OptimalDelaunayOptions {
                target_length: TargetLength::Absolute(h),
                max_insertions: InsertionLimit::Vertices(4 * goal),
                ..OptimalDelaunayOptions::default()
            })
            .map_err(|e| e.to_string())?;
            Ok(tri)
        })?;
        fem_row(name, &format!("iODT @4-1 x{level}"), &iodt)?;
    }
    Ok(())
}

fn heat(name: &str, mesh: &TriMesh, geometry: &Geometry) -> Result<(), String> {
    let n = mesh.n_vertices();
    let sources = [0, n / 3, 2 * n / 3];
    let mut exact = Vec::new();
    let start = Instant::now();
    for &s in &sources {
        let e = exact_distance(mesh, geometry, &[s], &ExactOptions::default())
            .map_err(|e| e.to_string())?;
        exact.push(e.distances);
    }
    eprintln!(
        "{}: exact distances in {:.2} s",
        short_name(name),
        start.elapsed().as_secs_f64()
    );
    let input = IntrinsicTriangulation::new(mesh, geometry).map_err(|e| e.to_string())?;
    let delaunay = idt(&input)?;
    let mut idr = delaunay.clone();
    idr.refine(&RefineOptions::default())
        .map_err(|e| e.to_string())?;
    let mean = delaunay.edge_lengths().iter().sum::<f64>() / delaunay.n_edges() as f64;
    let mut fine = delaunay.clone();
    fine.refine(&RefineOptions {
        max_circumradius: Some(0.5 * mean),
        ..RefineOptions::default()
    })
    .map_err(|e| e.to_string())?;
    for (kind, tri) in [
        ("input", &input),
        ("iDT", &delaunay),
        ("iDR", &idr),
        ("iDR r/2", &fine),
    ] {
        let solver = HeatSolver::new(tri, &HeatOptions::default()).map_err(|e| e.to_string())?;
        let (mut rel, mut abs, mut worst, mut count, mut query) = (0.0, 0.0, 0.0f64, 0usize, 0.0);
        let mut abs_count = 0usize;
        for (k, &s) in sources.iter().enumerate() {
            let d = solver
                .distance(&[HeatSource::Vertex(s)])
                .map_err(|e| e.to_string())?;
            query += d.report.seconds;
            let reference = &exact[k];
            let scale = reference.iter().fold(0.0f64, |m, &x| if x.is_finite() { m.max(x) } else { m });
            for (v, &e) in reference.iter().enumerate().take(n) {
                if !e.is_finite() {
                    continue;
                }
                let err = (d.distances[v] - e).abs();
                abs += err / scale;
                abs_count += 1;
                worst = worst.max(err / scale);
                if e > 0.0 {
                    rel += err / e;
                    count += 1;
                }
            }
        }
        println!(
            "row | heat | {} | {kind} | {} | {:.4} | {:.5} | {:.4} | {:.4} | {:.4}",
            short_name(name),
            tri.n_vertices(),
            solver.setup().seconds,
            query / sources.len() as f64,
            rel / count as f64,
            abs / abs_count as f64,
            worst
        );
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() != 2 {
        eprintln!("usage: apps_study steiner|fem|heat MESH");
        std::process::exit(2);
    }
    let (task, name) = (args[0].as_str(), args[1].as_str());
    let result = load(name).and_then(|(mesh, geometry)| match task {
        "steiner" => steiner(name, &mesh, &geometry),
        "fem" => fem(name, &mesh, &geometry),
        "heat" => heat(name, &mesh, &geometry),
        other => Err(format!("unknown task {other}")),
    });
    if let Err(e) = result {
        eprintln!("{name}: {e}");
        std::process::exit(1);
    }
}
