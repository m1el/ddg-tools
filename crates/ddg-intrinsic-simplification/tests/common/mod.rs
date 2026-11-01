//! Meshes and helpers shared by the integration tests.
#![allow(dead_code)]

use std::path::PathBuf;

use ddg_core::generators::{self, SplitMix64};
use ddg_core::io::read_obj;
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions, InputPoint, PointLocation};

/// A mesh with its positions.
pub struct Shape {
    pub name: String,
    pub positions: Vec<[f64; 3]>,
    pub faces: Vec<[u32; 3]>,
}

impl Shape {
    pub fn new(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Shape {
        Shape {
            name: name.to_string(),
            positions,
            faces,
        }
    }

    /// The face centroids as extra tracked points.
    pub fn centroids(&self) -> Vec<InputPoint> {
        (0..self.faces.len() as u32)
            .map(|f| InputPoint {
                face: f,
                bary: [1.0 / 3.0, 1.0 / 3.0, 1.0 - 2.0 / 3.0],
            })
            .collect()
    }

    /// Random points, one per face.
    pub fn random_points(&self, seed: u64) -> Vec<InputPoint> {
        let mut rng = SplitMix64::new(seed);
        (0..self.faces.len() as u32)
            .map(|f| {
                let a = rng.next_f64();
                let b = rng.next_f64();
                let (a, b) = if a + b > 1.0 {
                    (1.0 - a, 1.0 - b)
                } else {
                    (a, b)
                };
                InputPoint {
                    face: f,
                    bary: [1.0 - a - b, a, b],
                }
            })
            .collect()
    }

    pub fn coarse(&self, points: &[InputPoint]) -> CoarseMesh {
        CoarseMesh::from_positions(
            &self.positions,
            &self.faces,
            points,
            CoarseMeshOptions::default(),
        )
        .unwrap_or_else(|e| panic!("{}: {e}", self.name))
    }

    /// The position in space of an input point.
    pub fn input_position(&self, p: &InputPoint) -> [f64; 3] {
        let f = self.faces[p.face as usize];
        [0, 1, 2].map(|k| {
            (0..3)
                .map(|c| p.bary[c] * self.positions[f[c] as usize][k])
                .sum()
        })
    }

    /// The position of a tracked point reconstructed from its coarse location with the fine
    /// positions of the coarse vertices (meaningful when the metric is still the planar one).
    pub fn reconstruct(&self, loc: PointLocation) -> [f64; 3] {
        match loc {
            PointLocation::Vertex(v) => self.positions[v as usize],
            PointLocation::Face { corners, bary, .. } => [0, 1, 2].map(|k| {
                (0..3)
                    .map(|c| bary[c] * self.positions[corners[c] as usize][k])
                    .sum()
            }),
        }
    }
}

pub fn dist(p: [f64; 3], q: [f64; 3]) -> f64 {
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
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

/// Keenan Crane's models of `data/models/` that are present (spot, bob, blub, the crumpled
/// developable sheet).
pub fn crane_models() -> Vec<Shape> {
    [
        "models/spot/spot/spot_triangulated.obj",
        "models/bob/bob/bob_tri.obj",
        "models/blub/blub/blub_triangulated.obj",
        "models/crumpleddevelopable/CrumpledDevelopable.obj",
    ]
    .iter()
    .filter_map(|p| load_data(p))
    .collect()
}

/// The reference's meshes of `data/intrinsic-simplification/` that are present, up to
/// `max_vertices`.
pub fn reference_meshes(max_vertices: usize) -> Vec<Shape> {
    [
        "spot",
        "boxpart",
        "capsule",
        "spot_subdiv",
        "bs_rest",
        "dragon_fat",
    ]
    .iter()
    .filter_map(|n| load_data(&format!("intrinsic-simplification/{n}.obj")))
    .filter(|s| s.positions.len() <= max_vertices)
    .collect()
}

/// A jittered planar grid of `(n + 1)²` vertices in the unit square (disk topology).
pub fn planar_grid(n: usize, jitter: f64, seed: u64) -> Shape {
    let mut rng = SplitMix64::new(seed);
    let h = 1.0 / n as f64;
    let id = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let mut x = i as f64 * h;
            let mut y = j as f64 * h;
            if i > 0 && i < n {
                x += jitter * h * (rng.next_f64() - 0.5);
            }
            if j > 0 && j < n {
                y += jitter * h * (rng.next_f64() - 0.5);
            }
            positions.push([x, y, 0.0]);
        }
    }
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i, j + 1), id(i + 1, j + 1));
            if (i + j) % 2 == 0 {
                faces.push([a, b, d]);
                faces.push([a, d, c]);
            } else {
                faces.push([a, b, c]);
                faces.push([b, d, c]);
            }
        }
    }
    Shape::new(&format!("planar grid {n}"), positions, faces)
}

/// Small generated shapes, closed and with boundary, whose random flips produce Δ-complexes.
pub fn small_shapes() -> Vec<Shape> {
    let g = |name: &str, m: generators::GeneratedMesh| Shape::new(name, m.positions, m.faces);
    vec![
        g("icosphere 1", generators::icosphere(1).unwrap()),
        g("torus 5x3", generators::torus(5, 3, 1.0, 0.25).unwrap()),
        g("torus 12x6", generators::torus(12, 6, 1.0, 0.4).unwrap()),
        g("grid 6x4", generators::grid_disk(6, 4).unwrap()),
        g(
            "bumpy grid 8x8",
            generators::bumpy_grid_disk(8, 8, 0.25).unwrap(),
        ),
        g(
            "annulus 12x3",
            generators::annulus(12, 3, 0.3, 1.0).unwrap(),
        ),
        g(
            "open cylinder 9x4",
            generators::open_cylinder(9, 4, 0.5, 2.0).unwrap(),
        ),
        planar_grid(6, 0.3, 7),
    ]
}
