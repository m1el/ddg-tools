//! Timings of pre-refinement and how often post-refinement is refused (SPEC §6; WP5).
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example refine_study -- MESH.obj \
//!     [--points] [--signposts] [--fractions 0.5,0.25,0.1,0.05]
//! ```
//!
//! Prints the pre-refinement (build `T⁰`, refine, locate) times and vertex counts for the chosen
//! representation, then, for each fraction, coarsens the pre-refined mesh to that fraction of
//! its vertices in a shuffled order (no error metric: `simplify` is WP3) and tries the
//! post-refinement: its time and result, or the reason it was refused. Single thread.

use std::time::Instant;

use ddg_core::generators::SplitMix64;
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic_simplification::refine::{RefineError, RefinementOptions, Representation};
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions, InputPoint};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!("usage: refine_study MESH.obj [--points] [--signposts] [--fractions a,b,c]");
        std::process::exit(2);
    };
    let with_points = args.iter().any(|a| a == "--points");
    let representation = if args.iter().any(|a| a == "--signposts") {
        Representation::Signposts
    } else {
        Representation::DEFAULT
    };
    let fractions: Vec<f64> = args
        .iter()
        .position(|a| a == "--fractions")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or_else(|| vec![0.5, 0.25, 0.1, 0.05]);
    let obj = read_obj(path).unwrap_or_else(|e| {
        eprintln!("{e}");
        std::process::exit(1)
    });
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces.clone()).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1)
    });
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1)
    });
    let points: Vec<InputPoint> = if with_points {
        (0..faces.len() as u32)
            .map(|face| InputPoint {
                face,
                bary: [0.2, 0.3, 0.5],
            })
            .collect()
    } else {
        Vec::new()
    };
    let options = RefinementOptions {
        representation,
        make_simplicial_flips: args.iter().any(|a| a == "--repair").then_some(1_000_000),
        ..RefinementOptions::default()
    };
    let t0 = Instant::now();
    let (coarse, report) = CoarseMesh::from_refined(
        &mesh,
        &geometry,
        &points,
        CoarseMeshOptions::default(),
        &options,
    )
    .unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1)
    });
    println!(
        "{path} {representation:?}: {} vertices, {} faces, {} extra points",
        mesh.n_vertices(),
        mesh.n_faces(),
        points.len()
    );
    println!(
        "pre-refinement: +{} vertices ({} inserted, {} removed, {} flips), refine {:.3} s, \
         locate+build {:.3} s, total {:.3} s, min angle {:.2} deg (reachable {:.2}), \
         guarantee {}",
        report.n_added,
        report.refine.n_inserted,
        report.refine.n_removals,
        report.refine.n_flips,
        report.seconds_refine,
        report.seconds_locate,
        t0.elapsed().as_secs_f64(),
        report.min_angle.to_degrees(),
        report.refine.min_angle_reachable.to_degrees(),
        report.guarantee_applies
    );
    let n = coarse.n_vertices();
    for fraction in fractions {
        let mut c = coarse.clone();
        let target = ((n as f64 * fraction) as usize).max(4);
        let mut order = c.vertices().to_vec();
        let mut rng = SplitMix64::new(42);
        for i in (1..order.len()).rev() {
            order.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let t = Instant::now();
        for v in order {
            if c.n_vertices() <= target {
                break;
            }
            let _ = c.remove_vertex(v);
        }
        let t_remove = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let result = c.post_refine(&options);
        let dt = t.elapsed().as_secs_f64();
        match result {
            Ok((_, r)) => println!(
                "fraction {fraction}: {} vertices ({t_remove:.2} s) -> post-refinement +{} in \
                 {dt:.3} s (refine {:.3}, locate {:.3}), min angle {:.2} deg, guarantee {}, \
                 repair {:?}",
                c.n_vertices(),
                r.n_added,
                r.seconds_refine,
                r.seconds_locate,
                r.min_angle.to_degrees(),
                r.guarantee_applies,
                r.repair
            ),
            Err(RefineError::NotSimplicial { what }) => println!(
                "fraction {fraction}: {} vertices -> REFUSED: {what}",
                c.n_vertices()
            ),
            Err(e) => println!(
                "fraction {fraction}: {} vertices -> ERROR: {e}",
                c.n_vertices()
            ),
        }
    }
}
