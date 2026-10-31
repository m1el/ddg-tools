//! WP2: the cost of reading the geometry of the crossings, after Delaunay flipping.
//!
//! For each mesh the signpost triangulation is flipped to Delaunay and converted to integer
//! coordinates. Timed, single thread, best of five:
//!
//! * `extract_intrinsic_edge` of every edge (`ExtractIntrinsicEdge`, the path over the input mesh),
//!   each edge on its own (no state shared between edges);
//! * `extract_edge_geometry` of every input edge (Alg. 4 with the positions of Alg. 2), and of
//!   that the combinatorial `extract_edge` alone (WP1);
//! * the signposts' `trace_edges` (every edge traced over the input mesh), the yardstick;
//! * `certify()`, for scale.
//!
//! ```text
//! cargo run --release -p ddg-integer-coords --example extract_cost
//! ```

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::{IntegerTriangulation, StripReport};
use ddg_intrinsic::{DelaunayOptions, EdgeTraceOptions, IntrinsicTriangulation};

fn load(relative: &str) -> Option<(String, TriMesh, Geometry)> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(relative);
    let obj = read_obj(&path).ok()?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len() - 1 {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).ok()?;
    let geometry = Geometry::from_positions(&mesh, &obj.positions).ok()?;
    Some((path.file_stem()?.to_str()?.to_string(), mesh, geometry))
}

fn best_of<T>(runs: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut best = f64::INFINITY;
    let mut last = None;
    for _ in 0..runs {
        let t0 = Instant::now();
        let out = f();
        best = best.min(t0.elapsed().as_secs_f64());
        last = Some(out);
    }
    (best, last.unwrap())
}

fn main() {
    let mut meshes = Vec::new();
    meshes.extend(load("models/spot/spot/spot_triangulated.obj"));
    meshes.extend(load("optimal-cones/optimal-cones-data/data/camelbad.obj"));
    let options = EdgeTraceOptions::default();
    println!(
        "{:<20} {:>7} {:>8} {:>7} | {:>10} {:>9} {:>8} | {:>10} {:>9} {:>9} | {:>10} {:>9} | {:>9}",
        "mesh",
        "|E|",
        "Σn⁺",
        "along",
        "intr. ms",
        "ns/(E+Σ)",
        "walk",
        "input ms",
        "ns/(E⁰+Σ)",
        "comb. ms",
        "signp. ms",
        "ns/(E+Σ)",
        "certify"
    );
    for (name, mesh, geometry) in &meshes {
        let mut s = IntrinsicTriangulation::new(mesh, geometry).unwrap();
        s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        let (t, _) = IntegerTriangulation::from_intrinsic(&s, &options).unwrap();
        let ne = t.complex().n_edges();
        let sum = t.total_crossings() as f64;
        let (intrinsic, report) = best_of(5, || {
            let mut report = StripReport::default();
            let mut n_along = 0usize;
            for e in 0..ne {
                let (trace, r) = t.extract_intrinsic_edge_with_report(e, &options).unwrap();
                n_along += usize::from(trace.along_input_edge.is_some());
                report.merge(&r);
            }
            (report, n_along)
        });
        let (report, n_along) = report;
        let input_mesh = t.input().mesh();
        let (input, input_report) = best_of(5, || {
            let mut report = StripReport::default();
            for ie in 0..input_mesh.n_edges() {
                let g = t
                    .extract_edge_geometry(input_mesh.edge_halfedge(ie))
                    .unwrap();
                report.merge(&g.report);
            }
            report
        });
        let (combinatorial, _) = best_of(5, || {
            for ie in 0..input_mesh.n_edges() {
                t.extract_edge(input_mesh.edge_halfedge(ie)).unwrap();
            }
        });
        let (signposts, traces) = best_of(5, || s.trace_edges(&options).unwrap());
        assert_eq!(traces.report.n_crossings as f64, sum);
        let (certify, _) = best_of(5, || t.certify().unwrap());
        let per = |secs: f64, n: f64| secs * 1e9 / n;
        println!(
            "{:<20} {:>7} {:>8} {:>7} | {:>10.2} {:>9.0} {:>8} | {:>10.2} {:>9.0} {:>9.2} | \
             {:>10.2} {:>9.0} | {:>7.2}ms",
            name,
            ne,
            sum,
            n_along,
            intrinsic * 1e3,
            per(intrinsic, ne as f64 + sum),
            report.walk_steps,
            input * 1e3,
            per(input, input_mesh.n_edges() as f64 + sum),
            combinatorial * 1e3,
            signposts * 1e3,
            per(signposts, ne as f64 + sum),
            certify * 1e3,
        );
        println!(
            "{:<20} clamped {} + {}, largest relative length mismatch {:.1e} (intrinsic) {:.1e} \
             (input pieces)",
            "",
            report.clamped_crossings,
            input_report.clamped_crossings,
            report.max_relative_mismatch,
            input_report.max_relative_mismatch
        );
    }
}
