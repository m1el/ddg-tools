//! Test meshes, copied from `crates/ddg-intrinsic/tests/common/mod.rs` (same workspace).
#![allow(dead_code)]

use std::collections::HashMap;
use std::f64::consts::PI;
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
    /// Whether the surface is a region of a plane (every straight path is a straight segment).
    pub flat: bool,
}

impl Shape {
    pub fn new(name: &str, positions: Vec<[f64; 3]>, faces: Vec<[u32; 3]>, flat: bool) -> Shape {
        let mesh = TriMesh::new(positions.len(), faces).unwrap_or_else(|e| panic!("{name}: {e}"));
        let geometry =
            Geometry::from_positions(&mesh, &positions).unwrap_or_else(|e| panic!("{name}: {e}"));
        Shape {
            name: name.to_string(),
            mesh,
            geometry,
            positions,
            flat,
        }
    }

    pub fn generated(name: &str, g: GeneratedMesh, flat: bool) -> Shape {
        Shape::new(name, g.positions, g.faces, flat)
    }

    pub fn triangulation(&self) -> IntrinsicTriangulation {
        IntrinsicTriangulation::new(&self.mesh, &self.geometry)
            .unwrap_or_else(|e| panic!("{}: {e}", self.name))
    }

    pub fn mean_edge_length(&self) -> f64 {
        let l = self.geometry.edge_lengths();
        l.iter().sum::<f64>() / l.len() as f64
    }
}

pub fn dist(p: [f64; 3], q: [f64; 3]) -> f64 {
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

// ----- generators -----

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
    Shape::new(
        &format!("sheared parallelogram {n} sigma {sigma}"),
        positions,
        faces,
        true,
    )
}

/// Removes faces until no vertex is a pinch (its faces form one fan): at a vertex with several
/// fans the largest is kept. As in `ddg-variational-cuts/examples/bad_mesh_study.rs`.
fn prune_pinches(mut faces: Vec<[u32; 3]>) -> Vec<[u32; 3]> {
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    loop {
        let nf = faces.len();
        let mut parent: Vec<usize> = (0..3 * nf).collect();
        let mut directed: HashMap<(u32, u32), usize> = HashMap::new();
        for (f, face) in faces.iter().enumerate() {
            for k in 0..3 {
                directed.insert((face[k], face[(k + 1) % 3]), f);
            }
        }
        let corner = |f: usize, v: u32| 3 * f + faces[f].iter().position(|&x| x == v).unwrap();
        for (f, face) in faces.iter().enumerate() {
            for k in 0..3 {
                let (a, b) = (face[k], face[(k + 1) % 3]);
                if let Some(&g) = directed.get(&(b, a)) {
                    for v in [a, b] {
                        let (rx, ry) = (
                            find(&mut parent, corner(f, v)),
                            find(&mut parent, corner(g, v)),
                        );
                        parent[rx] = ry;
                    }
                }
            }
        }
        // The fans of every vertex, with their sizes.
        let mut fans: HashMap<u32, HashMap<usize, usize>> = HashMap::new();
        for (f, face) in faces.iter().enumerate() {
            for (k, &v) in face.iter().enumerate() {
                let root = find(&mut parent, 3 * f + k);
                *fans.entry(v).or_default().entry(root).or_insert(0) += 1;
            }
        }
        let mut drop = vec![false; nf];
        let mut any = false;
        for (&v, fan) in &fans {
            if fan.len() < 2 {
                continue;
            }
            any = true;
            let best = *fan
                .iter()
                .max_by_key(|(root, n)| (**n, std::cmp::Reverse(**root)))
                .unwrap()
                .0;
            for (f, face) in faces.iter().enumerate() {
                if let Some(k) = face.iter().position(|&x| x == v)
                    && find(&mut parent, 3 * f + k) != best
                {
                    drop[f] = true;
                }
            }
        }
        if !any {
            return faces;
        }
        faces = faces
            .into_iter()
            .zip(drop)
            .filter(|(_, d)| !d)
            .map(|(f, _)| f)
            .collect();
    }
}

/// The lattice `x = (i + σ j) h, y = j h` clipped to the unit square, then scaled to area 1:
/// every triangle has the largest angle `180° − atan(1/σ)` (153° for `σ = 2`), and the left
/// and right boundaries are ragged. The `sheared` family of
/// `ddg-variational-cuts/examples/bad_mesh_study.rs`.
pub fn sheared_square(n: usize, sigma: f64) -> Shape {
    let h = 1.0 / n as f64;
    let mut index: HashMap<(i64, i64), u32> = HashMap::new();
    let mut positions = Vec::new();
    let eps = 1e-9;
    for j in 0..=n as i64 {
        let lo = (-sigma * j as f64).floor() as i64 - 1;
        let hi = (n as f64 - sigma * j as f64).ceil() as i64 + 1;
        for i in lo..=hi {
            let x = (i as f64 + sigma * j as f64) * h;
            if x >= -eps && x <= 1.0 + eps {
                index.insert((i, j), positions.len() as u32);
                positions.push([x, j as f64 * h, 0.0]);
            }
        }
    }
    let mut faces = Vec::new();
    for (&(i, j), &a) in &index {
        let (b, c, d) = (
            index.get(&(i + 1, j)),
            index.get(&(i + 1, j + 1)),
            index.get(&(i, j + 1)),
        );
        if let (Some(&b), Some(&c), Some(&d)) = (b, c, d) {
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    faces.sort_unstable();
    let faces = prune_pinches(faces);
    // Only the vertices that are used, in the order of their first use.
    let mut map = vec![u32::MAX; positions.len()];
    let mut used = Vec::new();
    let faces: Vec<[u32; 3]> = faces
        .iter()
        .map(|f| {
            f.map(|v| {
                if map[v as usize] == u32::MAX {
                    map[v as usize] = used.len() as u32;
                    used.push(positions[v as usize]);
                }
                map[v as usize]
            })
        })
        .collect();
    let area: f64 = faces
        .iter()
        .map(|f| {
            let p = f.map(|v| used[v as usize]);
            0.5 * ((p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]))
        })
        .sum();
    let scale = 1.0 / area.sqrt();
    for p in &mut used {
        p[0] *= scale;
        p[1] *= scale;
    }
    Shape::new(
        &format!("sheared square {n} sigma {sigma}"),
        used,
        faces,
        true,
    )
}

/// Rows that alternate between thin strips (caps: one angle of 168°) and tall ones (needles:
/// apex 19°), shifted by half a cell; the flat unit square. As in
/// `ddg-variational-cuts/examples/bad_mesh_study.rs`.
pub fn needles_and_caps(n: usize) -> Shape {
    let d = 1.0 / n as f64;
    let (ha, hb) = (0.05 * d, 3.0 * d);
    let pairs = ((1.0 / (ha + hb)).round() as usize).max(1);
    let scale = 1.0 / (pairs as f64 * (ha + hb));
    let rows = 2 * pairs + 1;
    let mut y = vec![0.0f64];
    for j in 0..2 * pairs {
        let last = y[j];
        y.push(last + scale * if j % 2 == 0 { ha } else { hb });
    }
    *y.last_mut().unwrap() = 1.0;
    let mut positions = Vec::new();
    let mut row_x: Vec<Vec<f64>> = Vec::new();
    let mut start = Vec::new();
    for (j, &yj) in y.iter().enumerate().take(rows) {
        let xs: Vec<f64> = if j % 2 == 0 {
            (0..=n).map(|i| i as f64 * d).collect()
        } else {
            std::iter::once(0.0)
                .chain((0..n).map(|i| (i as f64 + 0.5) * d))
                .chain(std::iter::once(1.0))
                .collect()
        };
        start.push(positions.len() as u32);
        for &x in &xs {
            positions.push([x, yj, 0.0]);
        }
        row_x.push(xs);
    }
    let mut faces = Vec::new();
    for j in 0..rows - 1 {
        let (xa, xb) = (&row_x[j], &row_x[j + 1]);
        let (mut i, mut m) = (0usize, 0usize);
        while i < xa.len() - 1 || m < xb.len() - 1 {
            let bottom = if i == xa.len() - 1 {
                false
            } else if m == xb.len() - 1 {
                true
            } else {
                xa[i + 1] <= xb[m + 1]
            };
            let (a, b) = (start[j] + i as u32, start[j + 1] + m as u32);
            if bottom {
                faces.push([a, a + 1, b]);
                i += 1;
            } else {
                faces.push([a, b + 1, b]);
                m += 1;
            }
        }
    }
    Shape::new(&format!("needles and caps {n}"), positions, faces, true)
}

/// An open tube, long and thin, each ring turned by `twist` steps against the previous one
/// (`scripts/intrinsic-reference/gen_meshes.py`, `thin_cylinder`). Two boundary loops; it is
/// developable, so its intrinsic geometry is flat, but paths wind around it.
pub fn thin_cylinder(
    n_around: usize,
    n_along: usize,
    radius: f64,
    length: f64,
    twist: f64,
) -> Shape {
    tube(n_around, n_along, radius, length, twist, 0.0)
}

/// [`thin_cylinder`] with the vertices of every ring spread along the axis over `stagger`
/// ring distances, so that they lie on a helix. A straight path around the tube from a vertex
/// back to itself then passes no other vertex, and for a thin tube it is an edge of the
/// Delaunay triangulation: a self-edge.
pub fn tube(
    n_around: usize,
    n_along: usize,
    radius: f64,
    length: f64,
    twist: f64,
    stagger: f64,
) -> Shape {
    let mut positions = Vec::new();
    for i in 0..=n_along {
        for j in 0..n_around {
            let w = 2.0 * PI * (j as f64 + twist * i as f64) / n_around as f64;
            let along = i as f64 + stagger * j as f64 / n_around as f64;
            positions.push([
                radius * w.cos(),
                radius * w.sin(),
                length * along / n_along as f64,
            ]);
        }
    }
    let mut faces = Vec::new();
    for i in 0..n_along {
        for j in 0..n_around {
            let a = (i * n_around + j) as u32;
            let b = (i * n_around + (j + 1) % n_around) as u32;
            let (c, d) = (a + n_around as u32, b + n_around as u32);
            faces.push([a, b, d]);
            faces.push([a, d, c]);
        }
    }
    Shape::new(
        &format!("tube {n_around}x{n_along} r {radius} twist {twist} stagger {stagger}"),
        positions,
        faces,
        false,
    )
}

/// A coarse torus with a thin tube, each ring turned by `twist` steps
/// (`gen_meshes.py`, `long_torus`). Closed, genus 1.
pub fn long_torus(n_major: usize, n_minor: usize, big: f64, small: f64, twist: f64) -> Shape {
    let mut positions = Vec::new();
    for i in 0..n_major {
        let u = 2.0 * PI * i as f64 / n_major as f64;
        for j in 0..n_minor {
            let w = 2.0 * PI * (j as f64 + twist * i as f64) / n_minor as f64;
            let rr = big + small * w.cos();
            positions.push([rr * u.cos(), rr * u.sin(), small * w.sin()]);
        }
    }
    let mut faces = Vec::new();
    for i in 0..n_major {
        for j in 0..n_minor {
            let a = (i * n_minor + j) as u32;
            let b = (i * n_minor + (j + 1) % n_minor) as u32;
            let c = (((i + 1) % n_major) * n_minor + j) as u32;
            let d = (((i + 1) % n_major) * n_minor + (j + 1) % n_minor) as u32;
            faces.push([a, c, d]);
            faces.push([a, d, b]);
        }
    }
    Shape::new(
        &format!("long torus {n_major}x{n_minor} R {big} r {small} twist {twist}"),
        positions,
        faces,
        false,
    )
}

/// A tall tetrahedron whose base is tilted: the apex has a small cone angle and one base
/// vertex is much closer to it than the others. Closed, genus 0, four vertices.
pub fn needle_tetrahedron(height: f64, lift: f64) -> Shape {
    let positions = vec![
        [1.0, 0.0, lift],
        [-0.5, 0.75f64.sqrt(), 0.0],
        [-0.5, -(0.75f64.sqrt()), 0.0],
        [0.0, 0.0, height],
    ];
    let faces = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
    Shape::new(
        &format!("needle tetrahedron h {height} lift {lift}"),
        positions,
        faces,
        false,
    )
}

/// An icosphere whose vertices are moved along the sphere by up to `amplitude` mean edge
/// lengths (seeded).
pub fn jittered_icosphere(level: u32, amplitude: f64, seed: u64) -> Shape {
    let g = generators::icosphere(level).unwrap();
    let mut positions = g.positions.clone();
    let h = {
        let mut sum = 0.0;
        let mut count = 0usize;
        for f in &g.faces {
            for k in 0..3 {
                sum += dist(
                    g.positions[f[k] as usize],
                    g.positions[f[(k + 1) % 3] as usize],
                );
                count += 1;
            }
        }
        sum / count as f64
    };
    let mut rng = SplitMix64::new(seed);
    for p in &mut positions {
        for x in p.iter_mut() {
            *x += amplitude * h * rng.uniform(-1.0, 1.0);
        }
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        for x in p.iter_mut() {
            *x /= n;
        }
    }
    Shape::new(
        &format!("jittered icosphere {level} amp {amplitude}"),
        positions,
        g.faces,
        false,
    )
}

/// A flat grid with its interior vertices moved in the plane (seeded).
pub fn jittered_grid(n: usize, amplitude: f64, seed: u64) -> Shape {
    let g = generators::grid_disk(n, n).unwrap();
    let h = 1.0 / n as f64;
    let mut positions = g.positions.clone();
    let mut rng = SplitMix64::new(seed);
    for p in &mut positions {
        let inside = |x: f64| x > 0.5 * h && x < 1.0 - 0.5 * h;
        if inside(p[0]) && inside(p[1]) {
            p[0] += amplitude * h * rng.uniform(-1.0, 1.0);
            p[1] += amplitude * h * rng.uniform(-1.0, 1.0);
        }
    }
    let shape = Shape::new(
        &format!("jittered grid {n} amp {amplitude}"),
        positions,
        g.faces,
        true,
    );
    assert_not_folded(&shape);
    shape
}

/// Panics if a face of a planar shape is not counter-clockwise in the plane: the surface the
/// edge lengths describe would then be a folded sheet, not the planar region.
pub fn assert_not_folded(shape: &Shape) {
    let p = &shape.positions;
    for face in shape.mesh.faces() {
        let [a, b, c] = face.map(|v| p[v as usize]);
        let twice_area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        assert!(
            twice_area > 0.0,
            "{}: face {face:?} is folded over",
            shape.name
        );
    }
}

/// Closed surfaces and surfaces with boundary, well shaped and badly shaped.
pub fn standard_shapes() -> Vec<Shape> {
    vec![
        Shape::generated("icosphere 2", generators::icosphere(2).unwrap(), false),
        jittered_icosphere(2, 0.35, 5),
        Shape::generated(
            "torus 24x8",
            generators::torus(24, 8, 1.0, 0.4).unwrap(),
            false,
        ),
        Shape::generated(
            "torus 5x3",
            generators::torus(5, 3, 1.0, 0.25).unwrap(),
            false,
        ),
        long_torus(8, 3, 4.0, 0.2, 0.37),
        needle_tetrahedron(8.0, 4.0),
        Shape::generated("grid 6x4", generators::grid_disk(6, 4).unwrap(), true),
        jittered_grid(10, 0.3, 3),
        Shape::generated(
            "bumpy grid 12x12",
            generators::bumpy_grid_disk(12, 12, 0.25).unwrap(),
            false,
        ),
        Shape::generated(
            "annulus 16x3",
            generators::annulus(16, 3, 0.3, 1.0).unwrap(),
            true,
        ),
        Shape::generated(
            "open cylinder 9x4",
            generators::open_cylinder(9, 4, 0.5, 2.0).unwrap(),
            false,
        ),
        thin_cylinder(4, 12, 0.1, 10.0, 0.37),
        tube(3, 6, 0.02, 10.0, 0.37, 1.0),
        sheared_parallelogram(8, 1.0),
        sheared_parallelogram(8, 3.0),
        sheared_parallelogram(10, 6.0),
        needles_and_caps(6),
    ]
}

// ----- files -----

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
    Some(Shape::new(&name, obj.positions, faces, false))
}

/// The models of `data/models/` (spot, bob, blub) that are present.
pub fn models() -> Vec<Shape> {
    [
        "models/spot/spot/spot_triangulated.obj",
        "models/bob/bob/bob_tri.obj",
        "models/blub/blub/blub_triangulated.obj",
    ]
    .iter()
    .filter_map(|p| load_data(p))
    .collect()
}

/// The deliberately bad meshes of the optimal-cones data set that are present.
pub fn bad_meshes() -> Vec<Shape> {
    [
        "camelbad",
        "camelbadnonuniform",
        "camel_noise1",
        "camel_noise2",
        "camelshotnoise",
    ]
    .iter()
    .filter_map(|name| load_data(&format!("optimal-cones/optimal-cones-data/data/{name}.obj")))
    .collect()
}
