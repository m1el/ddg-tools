//! The cost of a tangent field on integer coordinates against signposts (Integer Coordinates
//! I26): the same mesh, flipped to Delaunay (and refined to 25° by the signposts' rules) on
//! both representations, then the frames, the connection matrices, the smoothest field and its
//! singularities, timed on one thread, best of `ROUNDS` interleaved rounds.
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-apps --example vector_cost -- <mesh>...
//!   <mesh>  a path below the workspace's data/ directory
//! ```

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, RefineOptions, RefineRules, Triangulation,
};
use ddg_intrinsic_apps::vector::*;

const ROUNDS: usize = 5;

fn load(name: &str) -> (TriMesh, Geometry) {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(name);
    let obj = read_obj(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..c.len() - 1 {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    (mesh, geometry)
}

fn prepare<T: Triangulation>(tri: &mut T, refine: bool) {
    tri.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
    if refine {
        tri.refine_with_rules(&RefineOptions::default(), RefineRules::Ruppert)
            .unwrap();
    }
}

/// Seconds of: frames, connection matrices (FE, n = 1), smoothest field (FE, n = 1),
/// singularities; and the eigenvalue.
fn run<T: TangentFrames>(tri: &T) -> ([f64; 4], f64) {
    let t = Instant::now();
    let frames = tri.tangent_frames().unwrap();
    std::hint::black_box(&frames);
    let t_frames = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let m = connection_matrices(tri, 1, ConnectionLaplacian::FiniteElement).unwrap();
    std::hint::black_box(&m);
    let t_matrices = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let (field, report) = smoothest_field(tri, &DirectionFieldOptions::default()).unwrap();
    let t_field = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let s = singularities(tri, &field).unwrap();
    std::hint::black_box(&s);
    let t_sing = t.elapsed().as_secs_f64();
    ([t_frames, t_matrices, t_field, t_sing], report.eigenvalue)
}

fn main() {
    let names: Vec<String> = std::env::args().skip(1).collect();
    println!(
        "{:<44} {:>8} | {:>9} {:>9} {:>9} {:>9} | {:>9} {:>9} {:>9} {:>9} | {:>6} {:>6}",
        "mesh",
        "vertices",
        "frames s",
        "matrix s",
        "field s",
        "sing s",
        "frames i",
        "matrix i",
        "field i",
        "sing i",
        "field+",
        "total"
    );
    for name in &names {
        let (mesh, geometry) = load(name);
        for refine in [false, true] {
            let mut s = IntrinsicTriangulation::new(&mesh, &geometry).unwrap();
            let mut i = IntegerTriangulation::new(&mesh, &geometry).unwrap();
            prepare(&mut s, refine);
            prepare(&mut i, refine);
            let same = s.complex() == i.complex();
            let mut best_s = [f64::INFINITY; 4];
            let mut best_i = [f64::INFINITY; 4];
            let (mut ls, mut li) = (0.0, 0.0);
            for _ in 0..ROUNDS {
                let (ts, l) = run(&s);
                ls = l;
                let (ti, l) = run(&i);
                li = l;
                for k in 0..4 {
                    best_s[k] = best_s[k].min(ts[k]);
                    best_i[k] = best_i[k].min(ti[k]);
                }
            }
            // The field call computes the frames once more, singularities once more.
            let total_s: f64 = best_s[2] + best_s[3];
            let total_i: f64 = best_i[2] + best_i[3];
            println!(
                "{:<44} {:>8} | {:>9.4} {:>9.4} {:>9.4} {:>9.4} | {:>9.4} {:>9.4} {:>9.4} {:>9.4} | {:>5.2}x {:>5.2}x{}",
                format!("{name} {}", if refine { "refined" } else { "iDT" }),
                s.complex().n_vertices(),
                best_s[0],
                best_s[1],
                best_s[2],
                best_s[3],
                best_i[0],
                best_i[1],
                best_i[2],
                best_i[3],
                best_i[2] / best_s[2],
                total_i / total_s,
                if same {
                    format!("  (λ {ls:.6} / {li:.6})")
                } else {
                    "  (different triangulations)".to_string()
                }
            );
        }
    }
}
