//! Meshes and helpers shared by the integration tests.
#![allow(dead_code)]

use std::path::PathBuf;

use ddg_core::generators::{GeneratedMesh, SplitMix64, grid_disk};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::{IntrinsicTriangulation, SurfacePoint};

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

pub fn dist(p: [f64; 3], q: [f64; 3]) -> f64 {
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// The position in space of a point of the input mesh.
pub fn position(shape: &Shape, point: SurfacePoint) -> [f64; 3] {
    let p = &shape.positions;
    match point {
        SurfacePoint::Vertex(v) => p[v as usize],
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = shape.mesh.edge_vertices(edge as usize);
            [0, 1, 2].map(|k| (1.0 - t) * p[a][k] + t * p[b][k])
        }
        SurfacePoint::Face { face, bary } => {
            let f = shape.mesh.faces()[face as usize];
            [0, 1, 2].map(|k| (0..3).map(|c| bary[c] * p[f[c] as usize][k]).sum())
        }
    }
}

/// The positions of all vertices of a triangulation of `shape`.
pub fn vertex_positions(shape: &Shape, tri: &IntrinsicTriangulation) -> Vec<[f64; 3]> {
    tri.vertex_locations()
        .iter()
        .map(|&p| position(shape, p))
        .collect()
}

/// The lattice `x = (i + σ j) h, y = j h`, `0 ≤ i, j ≤ n`: a flat parallelogram (convex) whose
/// triangles all have the largest angle `180° − atan(1/σ)`.
pub fn sheared_parallelogram(n: usize, sigma: f64) -> Shape {
    let h = 1.0 / n as f64;
    let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            positions.push([(i as f64 + sigma * j as f64) * h, j as f64 * h, 0.0]);
        }
    }
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j + 1));
            faces.push([a, b, d]);
            faces.push([a, d, c]);
        }
    }
    Shape::new(&format!("sheared parallelogram {n} sigma {sigma}"), positions, faces)
}

/// The flat unit square `n × n` with interior vertices moved by up to `amplitude` cells.
pub fn jittered_grid(n: usize, amplitude: f64, seed: u64) -> Shape {
    let g = grid_disk(n, n).unwrap();
    let h = 1.0 / n as f64;
    let mut rng = SplitMix64::new(seed);
    let mut positions = g.positions;
    for (i, p) in positions.iter_mut().enumerate() {
        let (x, y) = (i / (n + 1), i % (n + 1));
        if x > 0 && x < n && y > 0 && y < n {
            p[0] += amplitude * h * rng.uniform(-1.0, 1.0);
            p[1] += amplitude * h * rng.uniform(-1.0, 1.0);
        }
    }
    Shape::new(&format!("jittered grid {n}"), positions, g.faces)
}

/// The surface of the unit cube `[0, 1]³`, every face an `n × n` grid of squares split into two
/// triangles (diagonals alternating).
pub fn cube(n: usize) -> Shape {
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut index = std::collections::HashMap::new();
    let mut id = |p: [i64; 3], positions: &mut Vec<[f64; 3]>| -> u32 {
        *index.entry(p).or_insert_with(|| {
            positions.push(p.map(|c| c as f64 / n as f64));
            (positions.len() - 1) as u32
        })
    };
    let n = n as i64;
    let mut faces = Vec::new();
    // Each face: an origin and two axes (u, v) with u × v pointing outwards.
    let sides: [([i64; 3], [i64; 3], [i64; 3]); 6] = [
        ([0, 0, 0], [0, 1, 0], [1, 0, 0]), // z = 0, normal −z
        ([0, 0, n], [1, 0, 0], [0, 1, 0]), // z = 1, normal +z
        ([0, 0, 0], [1, 0, 0], [0, 0, 1]), // y = 0, normal −y
        ([0, n, 0], [0, 0, 1], [1, 0, 0]), // y = 1, normal +y
        ([0, 0, 0], [0, 0, 1], [0, 1, 0]), // x = 0, normal −x
        ([n, 0, 0], [0, 1, 0], [0, 0, 1]), // x = 1, normal +x
    ];
    for (o, u, v) in sides {
        let at = |i: i64, j: i64| [0, 1, 2].map(|k| o[k] + i * u[k] + j * v[k]);
        for i in 0..n {
            for j in 0..n {
                let a = id(at(i, j), &mut positions);
                let b = id(at(i + 1, j), &mut positions);
                let c = id(at(i + 1, j + 1), &mut positions);
                let d = id(at(i, j + 1), &mut positions);
                if (i + j) % 2 == 0 {
                    faces.push([a, b, c]);
                    faces.push([a, c, d]);
                } else {
                    faces.push([a, b, d]);
                    faces.push([b, c, d]);
                }
            }
        }
    }
    Shape::new(&format!("cube {n}"), positions, faces)
}

fn workspace_data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data")
}

/// Loads a mesh below the workspace's `data/` directory, splitting polygons into fans. Prints
/// a notice and returns `None` when the file is absent.
pub fn load_data(relative: &str) -> Option<Shape> {
    let path = workspace_data().join(relative);
    if !path.is_file() {
        eprintln!(
            "SKIPPED: {} not found; nothing was checked for it",
            path.display()
        );
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
