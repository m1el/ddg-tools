//! Scaling of `simplify` with tracking and the prolongation `P` (paper Fig. 16; SPEC §10.2).
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example scaling -- MESH.obj \
//!     [--levels k] [--fraction f] [--sign paper|cetm]
//! ```
//!
//! Subdivides MESH `k` times (4-1 midpoint subdivision from edge lengths,
//! `ddg_intrinsic_apps::fem::regular_subdivision`: the surface does not change), builds the
//! coarse mesh with a tracked point per fine vertex, runs `simplify` to `f · |V|` vertices
//! (default `0.1`) and builds `P`. **One size per process, single thread**: run it pinned
//! (`taskset -c N`) and note the machine load (the example prints the load average from
//! `/proc/loadavg` before and after). Prints one TSV line:
//!
//! `mesh levels n_vertices n_faces fraction n_coarse t_build t_simplify t_p t_total loadavg_start loadavg_end`
//!
//! `t_simplify` is the report's total (initial Delaunay, masses, costs, removals, requeues with
//! the point tracking), `t_p` the compaction plus `P`.

use std::time::Instant;

use ddg_core::Geometry;
use ddg_core::TriMesh;
use ddg_core::io::read_obj;
use ddg_intrinsic_apps::fem::regular_subdivision;
use ddg_intrinsic_simplification::{
    CoarseMesh, CoarseMeshOptions, ProjectiveSign, SimplifyOptions, Target, simplify,
};

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn die(msg: String) -> ! {
    eprintln!("{msg}");
    std::process::exit(1);
}

fn loadavg() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .map(|s| s.split_whitespace().take(3).collect::<Vec<_>>().join("/"))
        .unwrap_or_else(|_| "?".to_string())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        die("usage: scaling MESH.obj [--levels k] [--fraction f] [--sign paper|cetm]".into());
    };
    let levels: usize = arg(&args, "--levels")
        .and_then(|x| x.parse().ok())
        .unwrap_or(0);
    let fraction: f64 = arg(&args, "--fraction")
        .and_then(|x| x.parse().ok())
        .unwrap_or(0.1);
    let sign = match arg(&args, "--sign").as_deref() {
        None | Some("paper") => ProjectiveSign::Paper,
        Some("cetm") => ProjectiveSign::Cetm,
        Some(x) => die(format!("--sign: paper or cetm, not {x}")),
    };
    let obj = read_obj(path).unwrap_or_else(|e| die(format!("{path}: {e}")));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let load0 = loadavg();
    let mesh0 = TriMesh::new(obj.positions.len(), faces).unwrap_or_else(|e| die(format!("{e}")));
    let geo0 = Geometry::from_positions(&mesh0, &obj.positions)
        .unwrap_or_else(|e| die(format!("{e}")));
    let (mesh, geo) =
        regular_subdivision(&mesh0, &geo0, levels).unwrap_or_else(|e| die(format!("{e}")));
    let n = mesh.n_vertices();
    let t = Instant::now();
    let mut coarse = CoarseMesh::from_lengths(
        &mesh,
        geo.edge_lengths().to_vec(),
        &[],
        CoarseMeshOptions {
            projective_sign: sign,
            ..CoarseMeshOptions::default()
        },
    )
    .unwrap_or_else(|e| die(format!("{e}")));
    let t_build = t.elapsed().as_secs_f64();
    let opts = SimplifyOptions {
        target: Target::Ratio(fraction),
        ..SimplifyOptions::default()
    };
    let t = Instant::now();
    let report = simplify(&mut coarse, &opts).unwrap_or_else(|e| die(format!("{e}")));
    let t_simplify = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let out = coarse.compact().unwrap_or_else(|e| die(format!("{e}")));
    let p = out
        .scalar_prolongation()
        .unwrap_or_else(|e| die(format!("{e}")));
    let t_p = t.elapsed().as_secs_f64();
    let _ = p;
    let name = std::path::Path::new(path)
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    println!(
        "{name}\t{levels}\t{n}\t{}\t{fraction}\t{}\t{t_build:.3}\t{t_simplify:.3}\t{t_p:.3}\t{:.3}\t{load0}\t{}\t(report total {:.3}, stop {:?})",
        mesh.n_faces(),
        report.n_vertices_after,
        t_build + t_simplify + t_p,
        loadavg(),
        report.timings.total,
        report.stop
    );
}
