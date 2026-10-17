//! Tangent vector fields on the input mesh, its intrinsic Delaunay triangulation and its
//! intrinsic Delaunay refinement (paper Sec. 5.5, Fig. 21 and 22).
//!
//! ```text
//! cargo run --release -p ddg-intrinsic-apps --example vector_study -- smoothest <mesh>...
//! cargo run --release -p ddg-intrinsic-apps --example vector_study -- harmonic <mesh>...
//!   <mesh>  a path below the workspace's data/ directory, or `sheared:<n>:<shear>:<jitter>:<seed>`
//! ```
//!
//! `smoothest`: the smoothest 1-, 2- and 4-direction fields (finite element energy) and their
//! energy `λ · area` (the Rayleigh quotient of the field normalized to mean square length one),
//! the number of singular faces, the index sum and the time, on the three triangulations.
//! `harmonic`: the harmonic vector field with the inward normal on the boundary (finite
//! differences), and how many interior vertices violate the vector maximum principle.

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::generators::{self, SplitMix64};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{IntrinsicTriangulation, RefineOptions};
use ddg_intrinsic_apps::vector::*;

fn load(name: &str) -> (TriMesh, Geometry) {
    let (positions, faces) = if let Some(spec) = name.strip_prefix("sheared:") {
        let p: Vec<f64> = spec.split(':').map(|x| x.parse().unwrap()).collect();
        let (n, shear, jitter, seed) = (p[0] as usize, p[1], p[2], p[3] as u64);
        let mut g = generators::grid_disk(n, n).unwrap();
        let mut rng = SplitMix64::new(seed);
        let h = 1.0 / n as f64;
        for q in &mut g.positions {
            let inside = q[0] > 1e-9 && q[0] < 1.0 - 1e-9 && q[1] > 1e-9 && q[1] < 1.0 - 1e-9;
            if inside {
                q[0] += rng.uniform(-jitter, jitter) * h;
                q[1] += rng.uniform(-jitter, jitter) * h;
            }
            q[0] += shear * q[1];
        }
        (g.positions, g.faces)
    } else {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data").join(name);
        let obj = read_obj(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let mut faces = Vec::new();
        for f in 0..obj.n_faces() {
            let c = &obj.corner_verts[obj.face_corners(f)];
            for k in 1..c.len() - 1 {
                faces.push([c[0], c[k], c[k + 1]]);
            }
        }
        (obj.positions, faces)
    };
    let mesh = TriMesh::new(positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
    (mesh, geometry)
}

fn triangulations(mesh: &TriMesh, geometry: &Geometry) -> Vec<(&'static str, IntrinsicTriangulation, f64)> {
    let input = IntrinsicTriangulation::new(mesh, geometry).unwrap();
    let t = Instant::now();
    let mut idt = input.clone();
    idt.flip_to_delaunay(&Default::default()).unwrap();
    let s_idt = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut idr = input.clone();
    idr.refine(&RefineOptions::default()).unwrap();
    let s_idr = t.elapsed().as_secs_f64();
    vec![("input", input, 0.0), ("iDT", idt, s_idt), ("iDR", idr, s_idr)]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("");
    for name in &args[2..] {
        let (mesh, geometry) = load(name);
        println!(
            "## {name}: {} vertices, {} faces, {} boundary loops",
            mesh.n_vertices(),
            mesh.n_faces(),
            mesh.boundary_loops().len()
        );
        for (what, tri, seconds) in triangulations(&mesh, &geometry) {
            let q = tri.quality(0.0);
            match mode {
                "smoothest" => {
                    for n_sym in [1u32, 2, 4] {
                        let t = Instant::now();
                        let options = DirectionFieldOptions {
                            n_sym,
                            ..Default::default()
                        };
                        let (field, report) = smoothest_field(&tri, &options).unwrap();
                        let s = singularities(&tri, &field).unwrap();
                        println!(
                            "| {name} | {what} | {} | {n_sym} | {:.2} | {} | {} | {} | {:?} | {:.3} |",
                            tri.n_vertices(),
                            report.energy,
                            s.faces.len(),
                            s.index_sum,
                            report.iterations,
                            report.status,
                            seconds + t.elapsed().as_secs_f64()
                        );
                    }
                }
                "harmonic" => {
                    for laplacian in [ConnectionLaplacian::FiniteDifference, ConnectionLaplacian::FiniteElement] {
                        let options = HarmonicFieldOptions {
                            laplacian,
                            ..Default::default()
                        };
                        let t = Instant::now();
                        let (field, report) = harmonic_field(&tri, &options).unwrap();
                        let s = t.elapsed().as_secs_f64();
                        let mp = check_maximum_principle(&tri, &field, 1e-9);
                        println!(
                            "| {name} | {what} | {laplacian:?} | {} | {} | {:.2e} | {} | {} | {:.3} | {:?} | {:.3} |",
                            tri.n_vertices(),
                            q.n_negative_weights,
                            q.min_cotan_weight,
                            mp.n_checked,
                            mp.n_outside_cone,
                            mp.min_margin,
                            report.factorization,
                            seconds + s
                        );
                    }
                }
                _ => panic!("mode: smoothest | harmonic"),
            }
        }
    }
}
