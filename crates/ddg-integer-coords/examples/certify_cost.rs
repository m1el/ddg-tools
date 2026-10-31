//! SPEC I4: the cost of the certificate, measured after many random flips.
//!
//! For each mesh and each number of random flips (as a multiple of `|E|`), the signpost
//! triangulation is flipped, converted to integer coordinates, and `certify()` is timed
//! (best of five, single thread). The yardstick is the time of one signpost flip on the same
//! state, so the last column says how many flips one certificate costs.
//!
//! ```text
//! cargo run --release -p ddg-integer-coords --example certify_cost
//! ```

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::generators::{self, SplitMix64};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{EdgeTraceOptions, IntrinsicTriangulation};

fn shape(
    name: &str,
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
) -> (String, TriMesh, Geometry) {
    let mesh = TriMesh::new(positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    (name.to_string(), mesh, geometry)
}

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
    let name = path.file_stem()?.to_str()?.to_string();
    Some(shape(&name, obj.positions, faces))
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
    let g = generators::icosphere(4).unwrap();
    meshes.push(shape("icosphere 4", g.positions, g.faces));
    let g = generators::torus(48, 16, 1.0, 0.4).unwrap();
    meshes.push(shape("torus 48x16", g.positions, g.faces));
    let g = generators::bumpy_grid_disk(40, 40, 0.25).unwrap();
    meshes.push(shape("bumpy grid 40x40", g.positions, g.faces));
    meshes.extend(load("models/spot/spot/spot_triangulated.obj"));
    meshes.extend(load("optimal-cones/optimal-cones-data/data/camelbad.obj"));

    println!(
        "{:<18} {:>6} {:>7} {:>9} {:>11} {:>10} {:>11} {:>9} {:>11} {:>9}",
        "mesh",
        "flips",
        "|E|",
        "Σn⁺",
        "convert ms",
        "certify ms",
        "ns/(E+Σn⁺)",
        "flip µs",
        "cert/flips",
        "trace ms"
    );
    for (name, mesh, geometry) in &meshes {
        for factor in [0usize, 10, 50] {
            let mut t = IntrinsicTriangulation::new(mesh, geometry).unwrap();
            let mut rng = SplitMix64::new(7);
            let ne = t.n_edges();
            let n_flips = factor * ne;
            let t0 = Instant::now();
            let mut done = 0usize;
            for _ in 0..n_flips {
                let e = (rng.next_u64() % ne as u64) as usize;
                if t.flip_edge(e).is_ok() {
                    done += 1;
                }
            }
            let flip_us = if n_flips > 0 {
                1e6 * t0.elapsed().as_secs_f64() / n_flips as f64
            } else {
                f64::NAN
            };
            let t1 = Instant::now();
            let converted = IntegerTriangulation::from_intrinsic(&t, &EdgeTraceOptions::default());
            let convert_ms = 1e3 * t1.elapsed().as_secs_f64();
            let (it, _) = match converted {
                Ok(x) => x,
                Err(e) => {
                    println!("{name:<18} {factor:>5}E conversion failed: {e}");
                    continue;
                }
            };
            let (secs, report) = best_of(5, || it.certify().unwrap());
            let size = ne as f64 + report.n_crossings as f64;
            // The signposts' own trace of every edge, for comparison.
            let (trace_secs, _) =
                best_of(1, || t.trace_edges(&EdgeTraceOptions::default()).unwrap());
            println!(
                "{:<18} {:>5}E {:>7} {:>9} {:>11.1} {:>10.3} {:>11.1} {:>9.3} {:>11.0} {:>9.1}",
                name,
                factor,
                ne,
                report.n_crossings,
                convert_ms,
                1e3 * secs,
                1e9 * secs / size,
                flip_us,
                if flip_us.is_finite() {
                    1e6 * secs / flip_us
                } else {
                    f64::NAN
                },
                1e3 * trace_secs,
            );
            let _ = done;
        }
    }
}
