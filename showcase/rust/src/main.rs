//! Data for the showcase renders, as JSON.
//!
//! ```text
//! showcase-data cones <mesh.json> <out.json> <lambda>...
//! showcase-data intrinsic <mesh.json> <out.json> <flip-steps> <refine-steps> [odt-target]
//! showcase-data exact <mesh.json> <out.json> <source-vertex> [odt-target-factor]
//! ```
//!
//! `mesh.json` is `{"positions": [[x, y, z], ...], "triangles": [[a, b, c], ...]}`.
//!
//! `cones`: the raw optimum of the canonical problem of "Optimal Cone Singularities for
//! Conformal Flattening" (Soliman, Slepčev, Crane 2018) for every lambda, any sign, through
//! `ddg_optimal_cones::solve`.
//!
//! `intrinsic`: snapshots of an intrinsic triangulation ("Navigating Intrinsic Triangulations",
//! Sharp, Soliman, Crane 2019) traced over the input mesh: Delaunay flipping stopped after a
//! growing number of flips, then Delaunay refinement stopped after a growing number of
//! insertions, then the intrinsic optimal Delaunay triangulation.

use std::fs;

use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, FlipLimit, InsertionLimit, IntrinsicTriangulation,
    OptimalDelaunayOptions, RefineOptions, SurfacePoint, TargetLength,
};
use ddg_optimal_cones::{ConeProblem, SolveOptions, solve};
use serde_json::{Value, json};

type Res<T> = Result<T, Box<dyn std::error::Error>>;

fn read_mesh(path: &str) -> Res<(Vec<[f64; 3]>, Vec<[u32; 3]>)> {
    let v: Value = serde_json::from_str(&fs::read_to_string(path)?)?;
    let pos = v["positions"]
        .as_array()
        .ok_or("positions")?
        .iter()
        .map(|p| {
            let a = p.as_array().unwrap();
            [a[0].as_f64().unwrap(), a[1].as_f64().unwrap(), a[2].as_f64().unwrap()]
        })
        .collect();
    let tri = v["triangles"]
        .as_array()
        .ok_or("triangles")?
        .iter()
        .map(|t| {
            let a = t.as_array().unwrap();
            [a[0].as_u64().unwrap() as u32, a[1].as_u64().unwrap() as u32, a[2].as_u64().unwrap() as u32]
        })
        .collect();
    Ok((pos, tri))
}

fn cones(args: &[String]) -> Res<()> {
    let (pos, tri) = read_mesh(&args[0])?;
    let mesh = TriMesh::new(pos.len(), tri)?;
    let geometry = Geometry::from_positions(&mesh, &pos)?;
    let mut runs = Vec::new();
    for l in &args[2..] {
        let lambda: f64 = l.parse()?;
        let problem = ConeProblem::new(&mesh, &geometry, lambda);
        let t = std::time::Instant::now();
        let sol = solve(&problem, &SolveOptions::default())?;
        let (cv, ca): (Vec<u32>, Vec<f64>) = sol
            .measure
            .iter()
            .enumerate()
            .filter(|(_, m)| **m != 0.0)
            .map(|(i, m)| (i as u32, *m))
            .unzip();
        eprintln!(
            "lambda {lambda}: {} cones, distortion {:.4}, {:.2} s, {:?}",
            cv.len(),
            sol.distortion_optimal,
            t.elapsed().as_secs_f64(),
            sol.report.status
        );
        runs.push(json!({
            "lambda": lambda,
            "cone_vertices": cv,
            "cone_angles": ca,
            "u": sol.log_scale_optimal,
            "distortion": sol.distortion_optimal,
            "status": format!("{:?}", sol.report.status),
        }));
    }
    fs::write(&args[1], serde_json::to_string(&json!({ "runs": runs }))?)?;
    Ok(())
}

/// A point of the input mesh in space.
fn place(p: &SurfacePoint, pos: &[[f64; 3]], mesh: &TriMesh) -> [f64; 3] {
    let lerp = |a: [f64; 3], b: [f64; 3], t: f64| {
        [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])]
    };
    match *p {
        SurfacePoint::Vertex(v) => pos[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = mesh.edge_vertices(edge as usize);
            lerp(pos[a], pos[b], t)
        }
        SurfacePoint::Face { face, bary } => {
            let f = mesh.faces()[face as usize];
            let mut q = [0.0; 3];
            for k in 0..3 {
                for d in 0..3 {
                    q[d] += bary[k] * pos[f[k] as usize][d];
                }
            }
            q
        }
    }
}

/// The edges of the triangulation as polylines over the input mesh, and its vertices.
fn snapshot(tri: &IntrinsicTriangulation, pos: &[[f64; 3]], mesh: &TriMesh, label: &str) -> Res<Value> {
    let traces = tri.trace_edges(&EdgeTraceOptions::default())?;
    let mut lines = Vec::with_capacity(tri.n_edges());
    let mut is_original = Vec::with_capacity(tri.n_edges());
    for e in 0..tri.n_edges() {
        let h = tri.edge_halfedge(e);
        let (a, b) = (tri.from_vertex(h), tri.to_vertex(h));
        let mut pl = vec![place(&tri.vertex_location(a), pos, mesh)];
        let cr = traces.crossings(e);
        for c in cr {
            let [p, q] = mesh.edge_vertices(c.edge as usize);
            let t = c.t;
            pl.push([
                pos[p][0] + t * (pos[q][0] - pos[p][0]),
                pos[p][1] + t * (pos[q][1] - pos[p][1]),
                pos[p][2] + t * (pos[q][2] - pos[p][2]),
            ]);
        }
        pl.push(place(&tri.vertex_location(b), pos, mesh));
        is_original.push(cr.is_empty() && a < mesh.n_vertices() && b < mesh.n_vertices());
        lines.push(pl);
    }
    let verts: Vec<[f64; 3]> =
        (0..tri.n_vertices()).map(|v| place(&tri.vertex_location(v), pos, mesh)).collect();
    let q = tri.quality(1e-12);
    Ok(json!({
        "label": label,
        "lines": lines,
        "straight": is_original,
        "vertices": verts,
        "n_input_vertices": mesh.n_vertices(),
        "min_angle": q.min_angle.to_degrees(), "max_angle": q.max_angle.to_degrees(),
    }))
}

fn intrinsic(args: &[String]) -> Res<()> {
    let (pos, tri) = read_mesh(&args[0])?;
    let flip_steps: usize = args[2].parse()?;
    let refine_steps: usize = args[3].parse()?;
    let mesh = TriMesh::new(pos.len(), tri)?;
    let geometry = Geometry::from_positions(&mesh, &pos)?;
    let fresh = IntrinsicTriangulation::new(&mesh, &geometry)?;

    let mut out = Vec::new();
    // Total number of flips to Delaunay.
    let mut full = fresh.clone();
    let rep = full.flip_to_delaunay(&DelaunayOptions::default())?;
    let n_flips = rep.n_flips;
    eprintln!("delaunay: {n_flips} flips, {:?}", rep.status);
    for k in 0..=flip_steps {
        let limit = (n_flips * k + flip_steps - 1) / flip_steps.max(1);
        let mut t = fresh.clone();
        if limit > 0 {
            t.flip_to_delaunay(&DelaunayOptions { max_flips: FlipLimit::Flips(limit), ..DelaunayOptions::default() })?;
        }
        out.push(snapshot(&t, &pos, &mesh, &format!("flip {limit}/{n_flips}"))?);
    }
    // Refinement.
    let mut refined = full.clone();
    let rrep = refined.refine(&RefineOptions::default())?;
    let n_ins = refined.n_vertices() - mesh.n_vertices();
    eprintln!("refine: {n_ins} insertions, {:?}", rrep.status);
    for k in 1..=refine_steps {
        let limit = (n_ins * k + refine_steps - 1) / refine_steps;
        let mut t = full.clone();
        t.refine(&RefineOptions { max_insertions: InsertionLimit::Vertices(limit), ..RefineOptions::default() })?;
        out.push(snapshot(&t, &pos, &mesh, &format!("refine {limit}/{n_ins}"))?);
    }
    // Optimal Delaunay.
    let mut odt = full.clone();
    let odt_factor: f64 = args.get(4).map(|a| a.parse()).transpose()?.unwrap_or(1.0);
    let orep = odt.optimal_delaunay(&OptimalDelaunayOptions {
        target_length: TargetLength::MeanInputEdge(odt_factor),
        ..OptimalDelaunayOptions::default()
    })?;
    eprintln!("odt: {} inserted, {} iterations, {:?}", orep.n_inserted, orep.n_iterations, orep.status);
    out.push(snapshot(&odt, &pos, &mesh, "odt")?);
    fs::write(&args[1], serde_json::to_string(&json!({ "stages": out }))?)?;
    Ok(())
}

/// Exact polyhedral geodesic distance from one vertex (`ddg_intrinsic_apps::exact`).
fn exact(args: &[String]) -> Res<()> {
    let (pos, tri) = read_mesh(&args[0])?;
    let source: usize = args[2].parse()?;
    let mesh = TriMesh::new(pos.len(), tri)?;
    let geometry = Geometry::from_positions(&mesh, &pos)?;
    let d = ddg_intrinsic_apps::exact::exact_distance(
        &mesh,
        &geometry,
        &[source],
        &ddg_intrinsic_apps::exact::ExactOptions::default(),
    )?;
    fs::write(&args[1], serde_json::to_string(&json!({ "distances": d.distances }))?)?;
    Ok(())
}

fn main() -> Res<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("cones") => cones(&args[1..]),
        Some("intrinsic") => intrinsic(&args[1..]),
        Some("exact") => exact(&args[1..]),
        _ => Err("usage: showcase-data cones|intrinsic ...".into()),
    }
}
