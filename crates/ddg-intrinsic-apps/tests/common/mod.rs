//! Meshes shared by the integration tests of the amr and vector modules.
#![allow(dead_code)]

use std::path::PathBuf;

use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::IntrinsicTriangulation;

/// A mesh with its positions.
pub struct Shape {
    pub name: String,
    pub mesh: TriMesh,
    pub geometry: Geometry,
    pub positions: Vec<[f64; 3]>,
}

impl Shape {
    pub fn new(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Shape {
        let mesh = TriMesh::new(positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry =
            Geometry::from_positions(&mesh, &positions).unwrap_or_else(|e| panic!("{name}: {e}"));
        Shape {
            name: name.to_string(),
            mesh,
            geometry,
            positions,
        }
    }

    pub fn generated(name: &str, g: GeneratedMesh) -> Shape {
        Shape::new(name, g.positions, g.faces)
    }

    pub fn triangulation(&self) -> IntrinsicTriangulation {
        IntrinsicTriangulation::new(&self.mesh, &self.geometry)
            .unwrap_or_else(|e| panic!("{}: {e}", self.name))
    }
}

/// A flat grid on `[0, 1] × [0, 1]`, sheared by `shear` (every triangle obtuse for a shear of
/// 1.5) and jittered.
pub fn sheared_grid(n: usize, shear: f64, jitter: f64, seed: u64) -> Shape {
    let mut g = generators::grid_disk(n, n).unwrap();
    let mut rng = SplitMix64::new(seed);
    let h = 1.0 / n as f64;
    for p in &mut g.positions {
        let inside = p[0] > 1e-9 && p[0] < 1.0 - 1e-9 && p[1] > 1e-9 && p[1] < 1.0 - 1e-9;
        if inside {
            p[0] += rng.uniform(-jitter, jitter) * h;
            p[1] += rng.uniform(-jitter, jitter) * h;
        }
        p[0] += shear * p[1];
    }
    Shape::generated(&format!("sheared_grid_{n}"), g)
}

/// A bumpy disk-like grid (curved), sheared.
pub fn bumpy_sheared_grid(n: usize, amplitude: f64, shear: f64) -> Shape {
    let mut g = generators::bumpy_grid_disk(n, n, amplitude).unwrap();
    for p in &mut g.positions {
        p[0] += shear * p[1];
    }
    Shape::generated(&format!("bumpy_sheared_grid_{n}"), g)
}

pub fn icosphere(level: u32) -> Shape {
    Shape::generated(&format!("icosphere_{level}"), generators::icosphere(level).unwrap())
}

/// An icosphere with vertices moved randomly along the sphere (bad triangles, same surface
/// up to the chordal error).
pub fn jittered_icosphere(level: u32, amplitude: f64, seed: u64) -> Shape {
    let mut g = generators::icosphere(level).unwrap();
    let mut rng = SplitMix64::new(seed);
    for p in &mut g.positions {
        for x in p.iter_mut() {
            *x += rng.uniform(-amplitude, amplitude);
        }
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        for x in p.iter_mut() {
            *x /= n;
        }
    }
    Shape::generated(&format!("jittered_icosphere_{level}"), g)
}

fn workspace_data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data")
}

/// Loads a mesh below the workspace's `data/` directory (polygons split into fans); `None`
/// with a notice when absent.
pub fn load_data(relative: &str) -> Option<Shape> {
    let path = workspace_data().join(relative);
    if !path.is_file() {
        eprintln!("SKIPPED: {} not found; nothing was checked for it", path.display());
        return None;
    }
    let obj = read_obj(&path).unwrap_or_else(|e| panic!("{e}"));
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let name = path.file_stem().unwrap().to_str().unwrap().to_string();
    Some(Shape::new(&name, obj.positions, faces))
}

/// A mesh of the optimal-cones data set.
pub fn dataset(name: &str) -> Option<Shape> {
    load_data(&format!("optimal-cones/optimal-cones-data/data/{name}.obj"))
}
