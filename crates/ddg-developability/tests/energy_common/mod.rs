//! Test meshes and helpers shared by the energy tests.
#![allow(dead_code)]

use std::f64::consts::PI;
use std::path::PathBuf;

use ddg_core::TriMesh;
use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};

pub type Mesh = (TriMesh, Vec<[f64; 3]>);

pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
pub fn normalize(a: [f64; 3]) -> [f64; 3] {
    scale(a, 1.0 / norm(a))
}

/// Interior star: vertex 0 at the origin (plus jitter), ring of `k` vertices at random angles and
/// heights; faces `(0, m, m+1)`. Ring vertices are boundary vertices with 2 faces.
pub fn interior_star(k: usize, seed: u64, height: f64) -> Mesh {
    let mut r = SplitMix64::new(seed);
    let mut pos = vec![[
        r.uniform(-0.05, 0.05),
        r.uniform(-0.05, 0.05),
        r.uniform(-0.05, 0.05),
    ]];
    for m in 0..k {
        let t = 2.0 * PI * (m as f64 + r.uniform(-0.3, 0.3)) / k as f64;
        let rad = r.uniform(0.7, 1.3);
        pos.push([rad * t.cos(), rad * t.sin(), r.uniform(-height, height)]);
    }
    let faces = (0..k)
        .map(|m| [0, 1 + m as u32, 1 + ((m + 1) % k) as u32])
        .collect();
    (TriMesh::new(k + 1, faces).unwrap(), pos)
}

/// Boundary fan: vertex 0 on the boundary with `k` faces spanning a random opening angle.
pub fn boundary_fan(k: usize, seed: u64, height: f64) -> Mesh {
    let mut r = SplitMix64::new(seed);
    let open = r.uniform(0.6, 1.8) * PI;
    let mut pos = vec![[0.0, 0.0, 0.0]];
    for m in 0..=k {
        let t =
            open * (m as f64 + r.uniform(-0.25, 0.25) * (m > 0 && m < k) as u8 as f64) / k as f64;
        let rad = r.uniform(0.7, 1.3);
        pos.push([rad * t.cos(), rad * t.sin(), r.uniform(-height, height)]);
    }
    let faces = (0..k).map(|m| [0, 1 + m as u32, 2 + m as u32]).collect();
    (TriMesh::new(k + 2, faces).unwrap(), pos)
}

/// Regular cone: apex `(0, 0, h)` (vertex 0), ring `(cos 2πm/k, sin 2πm/k, 0)`.
pub fn regular_cone(k: usize, h: f64) -> Mesh {
    let mut pos = vec![[0.0, 0.0, h]];
    for m in 0..k {
        let t = 2.0 * PI * m as f64 / k as f64;
        pos.push([t.cos(), t.sin(), 0.0]);
    }
    let faces = (0..k)
        .map(|m| [0, 1 + m as u32, 1 + ((m + 1) % k) as u32])
        .collect();
    (TriMesh::new(k + 1, faces).unwrap(), pos)
}

pub fn perturbed_icosphere(level: u32, amp: f64, seed: u64) -> Mesh {
    let g = icosphere(level).unwrap();
    let mut p = g.positions.clone();
    perturb(&mut p, [amp; 3], seed);
    (g.topology().unwrap(), p)
}

pub fn bumpy_grid(n: usize, amp: f64, seed: u64) -> Mesh {
    let g = bumpy_grid_disk(n, n, amp).unwrap();
    let mut p = g.positions.clone();
    perturb(&mut p, [0.2 / n as f64; 3], seed);
    (g.topology().unwrap(), p)
}

pub fn perturbed_torus(seed: u64) -> Mesh {
    let g = torus(12, 8, 1.0, 0.4).unwrap();
    let mut p = g.positions.clone();
    perturb(&mut p, [0.03; 3], seed);
    (g.topology().unwrap(), p)
}

pub fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// Load an OBJ from `data/models/<rel>`, fan-triangulating polygons. `None` if missing.
pub fn load_model(rel: &str) -> Option<Mesh> {
    let path = workspace_root().join("data/models").join(rel);
    if !path.exists() {
        eprintln!("skipping: {} not found", path.display());
        return None;
    }
    let obj = ddg_core::io::read_obj(&path).unwrap();
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let r = obj.face_corners(f);
        let c: Vec<u32> = r.map(|i| obj.corner_verts[i]).collect();
        for i in 1..c.len() - 1 {
            faces.push([c[0], c[i], c[i + 1]]);
        }
    }
    Some((
        TriMesh::new(obj.positions.len(), faces).unwrap(),
        obj.positions,
    ))
}

pub fn mean_edge_length(mesh: &TriMesh, pos: &[[f64; 3]]) -> f64 {
    let e = mesh.edges();
    e.iter()
        .map(|&[a, b]| norm(sub(pos[a as usize], pos[b as usize])))
        .sum::<f64>()
        / e.len() as f64
}

/// A random rotation matrix (from a random unit quaternion).
pub fn random_rotation(seed: u64) -> [[f64; 3]; 3] {
    let mut r = SplitMix64::new(seed);
    let q = [
        r.uniform(-1.0, 1.0),
        r.uniform(-1.0, 1.0),
        r.uniform(-1.0, 1.0),
        r.uniform(-1.0, 1.0),
    ];
    let l = (q.iter().map(|x| x * x).sum::<f64>()).sqrt();
    let [w, x, y, z] = q.map(|v| v / l);
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

pub fn rotate(r: &[[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
    [dot(r[0], p), dot(r[1], p), dot(r[2], p)]
}
