//! Removes every removable vertex of a mesh in passes by ascending fine id (no error metric:
//! the order is the id), with the fine vertices tracked, and prints the time. Single thread.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-simplification --example remove_all -- MESH.obj \
//!     [--no-delaunay] [--centroids] [--target N] [--verbose]
//! ```
//!
//! `--target N` stops once `N` vertices are left.

use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_intrinsic_simplification::{
    CoarseMesh, CoarseMeshOptions, IemError, InputPoint, RemovalKind,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!("usage: remove_all MESH.obj [--no-delaunay] [--centroids]");
        std::process::exit(2);
    };
    let delaunay = !args.iter().any(|a| a == "--no-delaunay");
    let centroids = args.iter().any(|a| a == "--centroids");
    let verbose = args.iter().any(|a| a == "--verbose");
    let target: usize = args
        .iter()
        .position(|a| a == "--target")
        .and_then(|i| args.get(i + 1))
        .and_then(|x| x.parse().ok())
        .unwrap_or(0);
    let obj = match read_obj(path) {
        Ok(obj) => obj,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len().saturating_sub(1) {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let points: Vec<InputPoint> = if centroids {
        (0..faces.len() as u32)
            .map(|face| InputPoint {
                face,
                bary: [1.0 / 3.0, 1.0 / 3.0, 1.0 - 2.0 / 3.0],
            })
            .collect()
    } else {
        Vec::new()
    };
    let options = CoarseMeshOptions {
        delaunay_after_removal: delaunay,
        ..CoarseMeshOptions::default()
    };
    let t0 = Instant::now();
    let mut mesh = match CoarseMesh::from_positions(&obj.positions, &faces, &points, options) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("{path}: {e}");
            std::process::exit(1);
        }
    };
    let t_build = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let initial = mesh.flip_to_delaunay().map(|r| r.flips.len()).unwrap_or(0);
    let t_delaunay = t1.elapsed().as_secs_f64();
    let n0 = mesh.n_vertices();
    // The tentative removal of every vertex (what the error metric's initial costs take).
    let t_tent = Instant::now();
    let mut n_tentative = 0usize;
    for v in 0..mesh.n_fine_vertices() as u32 {
        if mesh.tentative_removal(v).is_ok() {
            n_tentative += 1;
        }
    }
    let t_tent = t_tent.elapsed().as_secs_f64();
    println!(
        "tentative removals of all {} vertices: {n_tentative} possible, {t_tent:.3} s ({:.2} \
         us each)",
        mesh.n_fine_vertices(),
        1e6 * t_tent / mesh.n_fine_vertices() as f64
    );
    let mut removed = [0usize; 3];
    let mut skipped = 0usize;
    let mut passes = 0;
    let (mut degree_flips, mut delaunay_flips, mut fallback) = (0usize, 0usize, 0usize);
    let t2 = Instant::now();
    loop {
        passes += 1;
        let mut progress = false;
        for v in 0..mesh.n_fine_vertices() as u32 {
            if !mesh.is_vertex(v) {
                continue;
            }
            if mesh.n_vertices() <= target {
                break;
            }
            match mesh.remove_vertex(v) {
                Ok(r) => {
                    progress = true;
                    removed[match r.prepared.kind {
                        RemovalKind::Interior => 0,
                        RemovalKind::Boundary => 1,
                        RemovalKind::Ear => 2,
                    }] += 1;
                    degree_flips += r.prepared.degree_flips.len();
                    fallback += r.prepared.fallback_flips.len();
                    delaunay_flips += r.delaunay.map_or(0, |d| d.flips.len());
                }
                Err(IemError::Removal(f)) => {
                    skipped += 1;
                    if verbose {
                        println!("{f} ({} vertices left)", mesh.n_vertices());
                    }
                }
                Err(e) => {
                    eprintln!("{path}: removal of {v}: {e}");
                    std::process::exit(1);
                }
            }
        }
        if !progress || mesh.n_vertices() <= target {
            break;
        }
    }
    let t_remove = t2.elapsed().as_secs_f64();
    let n_removed: usize = removed.iter().sum();
    let t3 = Instant::now();
    let cert = mesh.certify(false);
    let t_certify = t3.elapsed().as_secs_f64();
    println!(
        "{path}: {n0} -> {} vertices in {passes} passes; removed {n_removed} (interior {}, \
         boundary {}, ear {}), refusals {skipped}; flips: degree {degree_flips}, fallback \
         {fallback}, Delaunay {delaunay_flips}; initial Delaunay {initial} flips",
        mesh.n_vertices(),
        removed[0],
        removed[1],
        removed[2],
    );
    println!(
        "time: build {t_build:.3} s, initial Delaunay {t_delaunay:.3} s, removals {t_remove:.3} s \
         ({:.0} removals/s, {:.2} us/removal), certify {t_certify:.3} s ({})",
        n_removed as f64 / t_remove,
        1e6 * t_remove / n_removed as f64,
        if cert.is_ok() { "ok" } else { "FAILED" }
    );
    let per_face: Vec<usize> = (0..mesh.n_faces())
        .map(|f| mesh.points_in_face(f).len())
        .collect();
    println!(
        "points per face at the end: max {}, mean {:.1}",
        per_face.iter().max().unwrap_or(&0),
        per_face.iter().sum::<usize>() as f64 / per_face.len().max(1) as f64
    );
    if let Err(e) = cert {
        println!("{e}");
    }
}
