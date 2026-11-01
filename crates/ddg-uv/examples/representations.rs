//! The intrinsic triangulation of a cut session on signposts and on integer coordinates: how
//! long building it takes (flips, refinement, the optimal Delaunay iteration, the export) and
//! what it gives, per mesh and option (Integer Coordinates SPEC I29).
//!
//! ```text
//! cargo run --release -p ddg-uv --example representations [-- MESH.obj ...]
//! ```
//!
//! Without arguments: spot (quads), bob, blub and camelbad from the workspace's `data/`. Each
//! build runs alone, one thread; the best of three is printed.

use std::path::PathBuf;

use ddg_core::io::read_obj;
use ddg_uv::PolygonMesh;
use ddg_uv::cutting::{FlowKind, IntrinsicOptions, Representation, Session, Start, Triangulation};

fn load(path: &PathBuf) -> Option<PolygonMesh> {
    let obj = read_obj(path).ok()?;
    let polygons: Vec<Vec<u32>> = (0..obj.n_faces())
        .map(|f| obj.corner_verts[obj.face_corners(f)].to_vec())
        .collect();
    Some(PolygonMesh::from_polygons(obj.positions.clone(), &polygons))
}

fn main() {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data");
    let args: Vec<PathBuf> = std::env::args().skip(1).map(PathBuf::from).collect();
    let paths = if args.is_empty() {
        vec![
            data.join("models/spot/spot/spot_quadrangulated.obj"),
            data.join("models/bob/bob/bob_quad.obj"),
            data.join("models/blub/blub/blub_quadrangulated.obj"),
            data.join("optimal-cones/optimal-cones-data/data/camelbad.obj"),
        ]
    } else {
        args
    };
    println!(
        "{:<24} {:<9} {:<9} {:>7} {:>7} {:>6} {:>8} {:>10}",
        "mesh", "options", "repr", "verts", "added", "splits", "min°", "seconds"
    );
    for path in &paths {
        let Some(mesh) = load(path) else {
            println!("{}: not found, skipped", path.display());
            continue;
        };
        let name = path.file_stem().unwrap().to_string_lossy().to_string();
        for (label, options) in [
            ("flips", IntrinsicOptions::flips_only()),
            ("refined", IntrinsicOptions::refined()),
            ("optimized", IntrinsicOptions::optimized()),
        ] {
            for representation in [Representation::Signposts, Representation::Integer] {
                let tri = Triangulation::Intrinsic(IntrinsicOptions {
                    representation,
                    ..options
                });
                let mut best = f64::INFINITY;
                let mut report = None;
                for _ in 0..3 {
                    let session = Session::with_triangulation(
                        &mesh.input(),
                        2,
                        &Start::GeodesicCircle {
                            vertex: 0,
                            radius: 0.1,
                        },
                        FlowKind::Auto,
                        &tri,
                    );
                    match session {
                        Ok(s) => {
                            let r = s.intrinsic_report().cloned().expect("intrinsic");
                            best = best.min(r.seconds);
                            report = Some(r);
                        }
                        Err(e) => {
                            println!("{name} {label} {representation:?}: {e}");
                            break;
                        }
                    }
                }
                if let Some(r) = report {
                    println!(
                        "{:<24} {:<9} {:<9} {:>7} {:>7} {:>6} {:>8.3} {:>10.4}",
                        name,
                        label,
                        format!("{representation:?}"),
                        r.n_vertices,
                        r.n_added,
                        r.n_splits,
                        r.min_angle.to_degrees(),
                        best
                    );
                }
            }
        }
    }
}
