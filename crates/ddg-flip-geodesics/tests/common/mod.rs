//! Meshes and helpers shared by the integration tests.
#![allow(dead_code)]

use std::f64::consts::PI;
use std::path::PathBuf;

use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_flip_geodesics::{
    CertifyOptions, EdgePath, GeodesicCertificate, GeodesicNetwork, Segment, ShortenOptions,
    ShortenReport, ShortenStatus, Shortening, path_from_vertices,
};
use ddg_intrinsic::{DelaunayOptions, IntrinsicTriangulation, Triangulation};

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

    /// The intrinsic Delaunay triangulation.
    pub fn delaunay(&self) -> IntrinsicTriangulation {
        let mut tri = self.triangulation();
        tri.flip_to_delaunay(&DelaunayOptions::default())
            .unwrap_or_else(|e| panic!("{}: {e}", self.name));
        tri
    }
}

pub fn dist(p: [f64; 3], q: [f64; 3]) -> f64 {
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// The angle at `b` between `a` and `c` in space.
pub fn angle_at(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let u = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let v = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2])
        .sqrt()
        .atan2(dot)
}

/// The edge path through `vertices` (panics if two are not adjacent).
pub fn vpath<T: Triangulation>(tri: &T, vertices: &[usize], closed: bool) -> EdgePath {
    path_from_vertices(tri, vertices, closed).unwrap_or_else(|e| panic!("{vertices:?}: {e}"))
}

/// The segment from `a` to `b` (panics if they are not adjacent).
pub fn seg<T: Triangulation>(tri: &T, a: usize, b: usize) -> Segment {
    vpath(tri, &[a, b], false).segments[0]
}

// ----- generators -----

/// The flat grid of `ddg_core::generators::grid_disk`; vertex `(i, j)` is `i (ny + 1) + j`.
pub fn grid(nx: usize, ny: usize) -> Shape {
    Shape::generated(
        &format!("grid {nx}x{ny}"),
        generators::grid_disk(nx, ny).unwrap(),
    )
}

/// The vertex `(i, j)` of [`grid`].
pub fn gv(ny: usize, i: usize, j: usize) -> usize {
    i * (ny + 1) + j
}

/// A fan of `n` triangles around vertex 0 at the origin, ring vertices `1..=n` at the given
/// positions; closed fan (the centre is interior), the ring is the boundary.
pub fn fan(name: &str, ring: &[[f64; 3]]) -> Shape {
    let n = ring.len() as u32;
    let mut positions = vec![[0.0, 0.0, 0.0]];
    positions.extend_from_slice(ring);
    let faces = (0..n).map(|i| [0, 1 + i, 1 + (i + 1) % n]).collect();
    Shape::new(name, positions, faces)
}

/// A regular hexagon fan in the plane: flat, the centre has angle sum `2π`.
pub fn hexagon() -> Shape {
    let ring: Vec<[f64; 3]> = (0..6)
        .map(|i| {
            let w = PI / 3.0 * i as f64;
            [w.cos(), w.sin(), 0.0]
        })
        .collect();
    fan("hexagon", &ring)
}

/// A fan of six right isosceles triangles: the centre is a saddle with `Θ = 3π`.
pub fn saddle() -> Shape {
    let ring = vec![
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, -1.0],
    ];
    fan("saddle 3pi", &ring)
}

/// A closed square pyramid of height `h` over the square `[-1, 1]²` (apex is vertex 4, base
/// corners 0..4 counter-clockwise, base split along 0–2).
pub fn pyramid(h: f64) -> Shape {
    let positions = vec![
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 0.0, h],
    ];
    let faces = vec![
        [0, 1, 4],
        [1, 2, 4],
        [2, 3, 4],
        [3, 0, 4],
        [0, 2, 1],
        [0, 3, 2],
    ];
    Shape::new(&format!("pyramid h {h}"), positions, faces)
}

/// The lattice `x = (i + σ j) h, y = j h`: a flat parallelogram of obtuse triangles.
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
    )
}

/// A tall tetrahedron with a tilted base (as `ddg-intrinsic/tests/common`): its Delaunay
/// triangulation has a vertex of degree one and a self-edge.
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
    )
}

/// A coarse torus with a thin tube (as `ddg-intrinsic/tests/common`): its Delaunay
/// triangulation has self-edges around the tube.
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
    )
}

/// Generated shapes of every topology class, well and badly shaped.
pub fn generated_shapes() -> Vec<Shape> {
    let mut jittered = generators::icosphere(2).unwrap();
    generators::perturb(&mut jittered.positions, [0.05; 3], 7);
    vec![
        Shape::generated("icosphere 2", generators::icosphere(2).unwrap()),
        Shape::generated("jittered icosphere 2", jittered),
        Shape::generated("torus 12x6", generators::torus(12, 6, 1.0, 0.3).unwrap()),
        grid(6, 4),
        Shape::generated(
            "bumpy grid 8x8",
            generators::bumpy_grid_disk(8, 8, 0.2).unwrap(),
        ),
        Shape::generated(
            "annulus 12x3",
            generators::annulus(12, 3, 0.4, 1.0).unwrap(),
        ),
        Shape::generated(
            "open cylinder 9x4",
            generators::open_cylinder(9, 4, 0.5, 2.0).unwrap(),
        ),
        sheared_parallelogram(8, 3.0),
        sheared_parallelogram(10, 6.0),
        needle_tetrahedron(8.0, 4.0),
        long_torus(8, 3, 4.0, 0.2, 0.37),
        pyramid(0.5),
        saddle(),
    ]
}

// ----- files -----

fn workspace_data() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data")
}

/// Loads a mesh below the workspace's `data/`, splitting polygons into fans; `None` (with a
/// notice) when the file is absent.
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

/// Spot (`data/models`) and the bad camel of the optimal-cones data (IT §5), when present.
pub fn data_shapes() -> Vec<Shape> {
    [
        "models/spot/spot/spot_triangulated.obj",
        "optimal-cones/optimal-cones-data/data/camelbad.obj",
    ]
    .iter()
    .filter_map(|p| load_data(p))
    .collect()
}

/// Pairs of distinct vertices, seeded.
pub fn random_pairs(n_vertices: usize, n: usize, seed: u64) -> Vec<(usize, usize)> {
    let mut rng = SplitMix64::new(seed);
    let mut out = Vec::new();
    while out.len() < n {
        let a = (rng.next_u64() % n_vertices as u64) as usize;
        let b = (rng.next_u64() % n_vertices as u64) as usize;
        if a != b {
            out.push((a, b));
        }
    }
    out
}

/// Flips `n` random flippable edges, seeded.
pub fn random_flips(tri: &mut IntrinsicTriangulation, n: usize, seed: u64) -> usize {
    let mut rng = SplitMix64::new(seed);
    let mut done = 0;
    for _ in 0..n {
        let e = (rng.next_u64() % tri.n_edges() as u64) as usize;
        if tri.flip_edge(e).is_ok() {
            done += 1;
        }
    }
    done
}

// ----- generators for the shortening tests (F2) -----

/// The surface of the cube `[0, a]³`, every face an `n × n` grid of squares split into two
/// triangles (diagonals alternating); as `ddg-intrinsic-apps/tests/common`.
pub fn cube(n: usize, a: f64) -> Shape {
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut index = std::collections::HashMap::new();
    let mut id = |p: [i64; 3], positions: &mut Vec<[f64; 3]>| -> u32 {
        *index.entry(p).or_insert_with(|| {
            positions.push(p.map(|c| a * c as f64 / n as f64));
            (positions.len() - 1) as u32
        })
    };
    let n = n as i64;
    let mut faces = Vec::new();
    let sides: [([i64; 3], [i64; 3], [i64; 3]); 6] = [
        ([0, 0, 0], [0, 1, 0], [1, 0, 0]),
        ([0, 0, n], [1, 0, 0], [0, 1, 0]),
        ([0, 0, 0], [1, 0, 0], [0, 0, 1]),
        ([0, n, 0], [0, 0, 1], [1, 0, 0]),
        ([0, 0, 0], [0, 0, 1], [0, 1, 0]),
        ([n, 0, 0], [0, 1, 0], [0, 0, 1]),
    ];
    for (o, u, v) in sides {
        let at = |i: i64, j: i64| [0, 1, 2].map(|k| o[k] + i * u[k] + j * v[k]);
        for i in 0..n {
            for j in 0..n {
                let p = id(at(i, j), &mut positions);
                let q = id(at(i + 1, j), &mut positions);
                let r = id(at(i + 1, j + 1), &mut positions);
                let s = id(at(i, j + 1), &mut positions);
                if (i + j) % 2 == 0 {
                    faces.push([p, q, r]);
                    faces.push([p, r, s]);
                } else {
                    faces.push([p, q, s]);
                    faces.push([q, r, s]);
                }
            }
        }
    }
    Shape::new(&format!("cube {n} side {a}"), positions, faces)
}

/// The vertex of [`cube`] at a position (exact lattice point), by search.
pub fn vertex_at(shape: &Shape, p: [f64; 3]) -> usize {
    shape
        .positions
        .iter()
        .position(|&q| dist(p, q) < 1e-12)
        .unwrap_or_else(|| panic!("{}: no vertex at {p:?}", shape.name))
}

/// The lateral surface of the square pyramid of [`pyramid`] (no base): a polyhedral cone
/// with apex 4 and the base square as its boundary.
pub fn open_pyramid(h: f64) -> Shape {
    let positions = vec![
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 0.0, h],
    ];
    let faces = vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
    Shape::new(&format!("open pyramid h {h}"), positions, faces)
}

/// A planar L: the square `[0, 2]²` without the open quadrant `(1, 2]²`, `n` cells per unit,
/// diagonals alternating. The reflex corner `(1, 1)` is a boundary vertex with interior angle
/// `3π/2`.
pub fn l_shape(n: usize) -> Shape {
    let h = 1.0 / n as f64;
    let m = 2 * n;
    let inside = |i: usize, j: usize| !(i > n && j > n); // lattice points of the L
    let mut index = vec![u32::MAX; (m + 1) * (m + 1)];
    let mut positions = Vec::new();
    for i in 0..=m {
        for j in 0..=m {
            if inside(i, j) {
                index[i * (m + 1) + j] = positions.len() as u32;
                positions.push([i as f64 * h, j as f64 * h, 0.0]);
            }
        }
    }
    let at = |i: usize, j: usize| index[i * (m + 1) + j];
    let mut faces = Vec::new();
    for i in 0..m {
        for j in 0..m {
            if i >= n && j >= n {
                continue;
            }
            let (p, q, r, s) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            if (i + j) % 2 == 0 {
                faces.push([p, q, r]);
                faces.push([p, r, s]);
            } else {
                faces.push([p, q, s]);
                faces.push([q, r, s]);
            }
        }
    }
    Shape::new(&format!("L-shape {n}"), positions, faces)
}

/// The distance in the plane L of [`l_shape`] between two of its points: straight if the
/// segment stays in the L (it misses the open quadrant `(1, 2)²`), else through the reflex
/// corner `(1, 1)` (the visibility graph of a polygon with one reflex vertex).
pub fn l_distance(p: [f64; 3], q: [f64; 3]) -> f64 {
    // Clip the segment against the open square (1, 2)² (Liang-Barsky).
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    let d = [q[0] - p[0], q[1] - p[1]];
    let mut visible = false;
    for k in 0..2 {
        for (lo, hi) in [(1.0, 2.0)] {
            if d[k].abs() < 1e-300 {
                if p[k] <= lo || p[k] >= hi {
                    visible = true;
                }
            } else {
                let (a, b) = ((lo - p[k]) / d[k], (hi - p[k]) / d[k]);
                t0 = t0.max(a.min(b));
                t1 = t1.min(a.max(b));
            }
        }
    }
    if visible || t1 - t0 <= 1e-12 {
        dist(p, q)
    } else {
        let r = [1.0, 1.0, 0.0];
        dist(p, r) + dist(r, q)
    }
}

/// A tube of revolution around the `z` axis: ring `j` has radius `radii[j]` at height
/// `j · dz`, `n` vertices; vertex `(j, k)` is `j n + k`. Ring neighbours `(j, k)–(j, k + 1)`,
/// verticals `(j, k)–(j + 1, k)` and diagonals `(j, k)–(j + 1, k + 1)` are edges. With
/// constant radii it is a prism tube.
pub fn tube(n: usize, radii: &[f64], dz: f64) -> Shape {
    let mut positions = Vec::new();
    for (j, &r) in radii.iter().enumerate() {
        for k in 0..n {
            let w = 2.0 * PI * k as f64 / n as f64;
            positions.push([r * w.cos(), r * w.sin(), j as f64 * dz]);
        }
    }
    let v = |j: usize, k: usize| (j * n + k % n) as u32;
    let mut faces = Vec::new();
    for j in 0..radii.len() - 1 {
        for k in 0..n {
            let (a, b, c, d) = (v(j, k), v(j, k + 1), v(j + 1, k), v(j + 1, k + 1));
            faces.push([a, b, d]);
            faces.push([a, d, c]);
        }
    }
    Shape::new(&format!("tube {n}x{}", radii.len()), positions, faces)
}

/// A tube with a waist: `2 w + 1` rings, radius `1 + 0.15 (j − w)²`, so ring `w` is the
/// shortest closed geodesic around it.
pub fn waist_tube(n: usize, w: usize) -> Shape {
    let radii: Vec<f64> = (0..=2 * w)
        .map(|j| 1.0 + 0.15 * (j as f64 - w as f64).powi(2))
        .collect();
    tube(n, &radii, 0.3)
}

// ----- running the shortening (F2) -----

/// A shortened network with its report and certificate.
pub struct Outcome<T: Triangulation> {
    pub network: GeodesicNetwork<T>,
    pub report: ShortenReport,
    pub certificate: GeodesicCertificate,
    /// Iterations seen by the observer, by kind.
    pub replaced: usize,
    pub contracted: usize,
    pub loop_moves: usize,
    pub stalled: usize,
}

/// Shortens `paths` on `tri` with `options` (and `check_invariants`), checking after every
/// iteration (SPEC §7 tests 5, 6): a replacement strictly shortens, nothing but an App. B
/// move lengthens, the running length agrees with the length summed afresh; then certifies.
pub fn run<T: Triangulation + Clone>(
    tri: T,
    paths: &[EdgePath],
    pinned: &[usize],
    options: &ShortenOptions,
) -> Outcome<T> {
    let mut network =
        GeodesicNetwork::new(tri, paths, pinned, None).unwrap_or_else(|e| panic!("network: {e}"));
    let mut last = network.length();
    let (mut replaced, mut contracted, mut loop_moves, mut stalled) = (0, 0, 0, 0);
    let options = ShortenOptions {
        check_invariants: true,
        ..*options
    };
    let report = network
        .shorten_with(&options, |step, net| {
            match &step.shortening {
                Shortening::Replaced {
                    old_length,
                    new_length,
                    ..
                } => {
                    assert!(new_length < old_length, "step {}: {step:?}", step.iteration);
                    replaced += 1;
                }
                Shortening::Contracted { .. } => contracted += 1,
                Shortening::LoopMoved { .. } => loop_moves += 1,
                Shortening::Stalled { .. } => stalled += 1,
            }
            let fresh = net.recompute_length();
            assert!(
                (fresh - step.length).abs() <= 1e-9 * fresh.max(1.0),
                "step {}: running length {} vs {fresh}",
                step.iteration,
                step.length
            );
            if !matches!(step.shortening, Shortening::LoopMoved { .. }) {
                assert!(step.length <= last, "step {}: length grew", step.iteration);
            }
            last = step.length;
        })
        .unwrap_or_else(|e| panic!("shorten: {e}"));
    let certificate = network
        .certify(&CertifyOptions::default())
        .unwrap_or_else(|e| panic!("certify: {e}"));
    if report.status == ShortenStatus::Geodesic && options.straight_tolerance <= 1e-4 {
        assert!(
            certificate.geodesic,
            "status Geodesic but the certificate fails: {:?}",
            certificate.failures
        );
    }
    Outcome {
        network,
        report,
        certificate,
        replaced,
        contracted,
        loop_moves,
        stalled,
    }
}

/// The vertices of path `p` of a network.
pub fn path_vertices<T: Triangulation>(net: &GeodesicNetwork<T>, p: usize) -> Vec<usize> {
    net.network()
        .path(p)
        .vertices(net.triangulation().complex())
}
