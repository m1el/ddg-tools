//! The exact distance on Δ-complexes with edge lengths (N37b): equal to the `TriMesh` entry
//! point on simplicial meshes, against brute-force unfolding and an analytic lattice distance
//! on coarsened complexes (doubled edges, self-edges), and face-point sources against
//! splitting the face.

mod common;

use common::{Shape, cube, jittered_grid};
use ddg_core::TriMesh;
use ddg_core::generators::{annulus, bumpy_grid_disk, icosphere, torus};
use ddg_core::geometry::triangle_area;
use ddg_intrinsic::Complex;
use ddg_intrinsic_apps::exact::{
    ExactError, ExactOptions, ExactSource, exact_distance, exact_distance_complex,
};

type V2 = [f64; 2];

fn sub(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn cross(a: V2, b: V2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn norm(a: V2) -> f64 {
    a[0].hypot(a[1])
}

/// The apex above the x-axis of the triangle on `[0, l_ab]` with sides `l_bc`, `l_ca`.
fn apex(l_ab: f64, l_bc: f64, l_ca: f64) -> V2 {
    let x = (l_ab * l_ab + l_ca * l_ca - l_bc * l_bc) / (2.0 * l_ab);
    [x, (l_ca * l_ca - x * x).max(0.0).sqrt()]
}

/// A Δ-complex with its edge lengths, built from a mesh by flips and removals.
#[derive(Clone)]
struct Surface {
    c: Complex,
    l: Vec<f64>,
    /// The id the vertex had in the mesh we started from.
    orig: Vec<usize>,
}

impl Surface {
    fn from_mesh(mesh: &TriMesh, l: Vec<f64>) -> Surface {
        Surface {
            c: Complex::from_tri_mesh(mesh),
            l,
            orig: (0..mesh.n_vertices()).collect(),
        }
    }

    fn length(&self, h: usize) -> f64 {
        self.l[self.c.edge(h)]
    }

    fn valid(&self) -> bool {
        (0..self.c.n_faces()).all(|f| {
            triangle_area(
                self.length(3 * f),
                self.length(3 * f + 1),
                self.length(3 * f + 2),
            )
            .is_some_and(|a| a > 1e-9)
        })
    }

    /// Flips edge `e` when the unfolded quadrilateral is comfortably convex.
    fn flip(&mut self, e: usize) -> bool {
        if self.c.flip_refusal(e).is_some() {
            return false;
        }
        let h0 = self.c.edge_halfedge(e);
        let h1 = self.c.twin(h0).unwrap();
        let a = self.l[e];
        // i = (0, 0), j = (a, 0); k above, l below.
        let k = apex(
            a,
            self.length(self.c.next(h0)),
            self.length(self.c.prev(h0)),
        );
        let lp = apex(
            a,
            self.length(self.c.prev(h1)),
            self.length(self.c.next(h1)),
        );
        let l = [lp[0], -lp[1]];
        let t = k[1] / (k[1] - l[1]);
        let x = k[0] + t * (l[0] - k[0]);
        if !(x > 0.1 * a && x < 0.9 * a) {
            return false;
        }
        let m = norm(sub(k, l));
        let mut after = self.clone();
        after.c.flip_edge(e).unwrap();
        after.l[e] = m;
        if !after.valid() {
            return false;
        }
        *self = after;
        true
    }

    /// Removes vertex `v` if it has degree three and the outer triangle is a triangle.
    fn remove(&mut self, v: usize) -> bool {
        let mut after = self.clone();
        let Ok(r) = after.c.remove_vertex(v) else {
            return false;
        };
        r.apply_edge_data(&mut after.l);
        r.apply_vertex_data(&mut after.orig);
        if !after.valid() {
            return false;
        }
        *self = after;
        true
    }

    /// Flips down to degree three and removes, until no vertex can go or `target` is reached.
    fn coarsen(&mut self, target: usize) {
        let mut progress = true;
        while progress && self.c.n_vertices() > target {
            progress = false;
            for v in 0..self.c.n_vertices() {
                if self.c.is_boundary_vertex(v) {
                    continue;
                }
                let before = self.clone();
                let mut budget = 40;
                while self.c.vertex_degree(v) > 3 && budget > 0 {
                    budget -= 1;
                    let edges: Vec<usize> = self
                        .c
                        .outgoing_halfedges(v)
                        .map(|h| self.c.edge(h))
                        .collect();
                    if !edges.into_iter().any(|e| self.flip(e)) {
                        break;
                    }
                }
                if self.c.vertex_degree(v) == 3 && self.remove(v) {
                    progress = true;
                    break;
                }
                *self = before;
            }
        }
    }

    fn has_self_edge(&self) -> bool {
        (0..self.c.n_edges()).any(|e| {
            let [a, b] = self.c.edge_vertices(e);
            a == b
        })
    }

    fn has_double_edge(&self) -> bool {
        let mut seen = std::collections::HashSet::new();
        (0..self.c.n_edges()).any(|e| {
            let [a, b] = self.c.edge_vertices(e);
            !seen.insert((a.min(b), a.max(b)))
        })
    }

    fn has_repeated_corner(&self) -> bool {
        self.c
            .faces()
            .iter()
            .any(|f| f[0] == f[1] || f[1] == f[2] || f[0] == f[2])
    }
}

fn edge_lengths_of(shape: &Shape) -> Vec<f64> {
    shape.geometry.edge_lengths().to_vec()
}

fn dist_opts() -> ExactOptions {
    ExactOptions::default()
}

fn vertex_sources(s: &[usize]) -> Vec<ExactSource> {
    s.iter().map(|&v| ExactSource::Vertex(v)).collect()
}

// ----- brute-force unfolding -----

/// Where a strip unfolding starts.
#[derive(Clone, Copy)]
enum Start {
    Corner { face: usize, k: usize },
    Inside { face: usize, bary: [f64; 3] },
}

struct Brute<'a> {
    s: &'a Surface,
    p: V2,
    best: Vec<f64>,
    lo_hi_eps: f64,
}

impl Brute<'_> {
    fn offer(&mut self, v: usize, x: V2) {
        let d = norm(sub(x, self.p));
        if d < self.best[v] {
            self.best[v] = d;
        }
    }

    fn in_cone(&self, x: V2, lo: V2, hi: V2) -> bool {
        let q = sub(x, self.p);
        let n = norm(q);
        cross(lo, q) >= -self.lo_hi_eps * n * norm(lo)
            && cross(q, hi) >= -self.lo_hi_eps * n * norm(hi)
    }

    fn cross_edge(&mut self, h: usize, q1: V2, q2: V2, lo: V2, hi: V2, depth: usize) {
        let c = &self.s.c;
        let Some(t) = c.twin(h) else { return };
        if depth == 0 {
            return;
        }
        let l = self.s.length(h);
        let ex = [(q1[0] - q2[0]) / l, (q1[1] - q2[1]) / l];
        let ey = [-ex[1], ex[0]];
        let ap = apex(l, self.s.length(c.next(t)), self.s.length(c.prev(t)));
        let x = [
            q2[0] + ap[0] * ex[0] + ap[1] * ey[0],
            q2[1] + ap[0] * ex[1] + ap[1] * ey[1],
        ];
        let xv = c.to_vertex(c.next(t));
        if self.in_cone(x, lo, hi) {
            self.offer(xv, x);
        }
        let exits = [(c.next(t), q1, x), (c.prev(t), x, q2)];
        for (e, a, b) in exits {
            let (ra, rb) = (sub(a, self.p), sub(b, self.p));
            let (lo2, hi2) = if cross(ra, rb) > 0.0 {
                (ra, rb)
            } else {
                (rb, ra)
            };
            let lo3 = if cross(lo, lo2) >= 0.0 { lo2 } else { lo };
            let hi3 = if cross(hi, hi2) >= 0.0 { hi } else { hi2 };
            if cross(lo3, hi3) > 1e-12 * norm(lo3) * norm(hi3) {
                self.cross_edge(e, a, b, lo3, hi3, depth - 1);
            }
        }
    }
}

/// The straight-segment distances from a start to every vertex, through strips of at most
/// `depth` faces.
fn strips(s: &Surface, start: Start, depth: usize) -> Vec<f64> {
    let c = &s.c;
    let face = match start {
        Start::Corner { face, .. } | Start::Inside { face, .. } => face,
    };
    let l0 = s.length(3 * face);
    let pos = [
        [0.0, 0.0],
        [l0, 0.0],
        apex(l0, s.length(3 * face + 1), s.length(3 * face + 2)),
    ];
    let p = match start {
        Start::Corner { k, .. } => pos[k],
        Start::Inside { bary, .. } => [
            bary[0] * pos[0][0] + bary[1] * pos[1][0] + bary[2] * pos[2][0],
            bary[0] * pos[0][1] + bary[1] * pos[1][1] + bary[2] * pos[2][1],
        ],
    };
    let mut b = Brute {
        s,
        p,
        best: vec![f64::INFINITY; c.n_vertices()],
        lo_hi_eps: 1e-9,
    };
    for k in 0..3 {
        b.offer(c.faces()[face][k] as usize, pos[k]);
    }
    for k in 0..3 {
        let (q1, q2) = (pos[k], pos[(k + 1) % 3]);
        let (ra, rb) = (sub(q1, p), sub(q2, p));
        if cross(ra, rb) > 1e-12 * norm(ra) * norm(rb) {
            b.cross_edge(3 * face + k, q1, q2, ra, rb, depth);
        }
    }
    b.best
}

/// The geodesic distance to every vertex from the sources, by Dijkstra over the vertices with
/// brute-force strips as the edges of the graph.
fn brute_force(s: &Surface, sources: &[ExactSource], depth: usize) -> Vec<f64> {
    let n = s.c.n_vertices();
    let mut d = vec![f64::INFINITY; n];
    for src in sources {
        match *src {
            ExactSource::Vertex(v) => d[v] = 0.0,
            ExactSource::Point { face, bary } => {
                let st = strips(s, Start::Inside { face, bary }, depth);
                for v in 0..n {
                    d[v] = d[v].min(st[v]);
                }
            }
        }
    }
    let mut done = vec![false; n];
    loop {
        let Some(u) = (0..n)
            .filter(|&v| !done[v] && d[v].is_finite())
            .min_by(|&a, &b| d[a].total_cmp(&d[b]))
        else {
            break;
        };
        done[u] = true;
        for f in 0..s.c.n_faces() {
            for k in 0..3 {
                if s.c.faces()[f][k] as usize == u {
                    let st = strips(s, Start::Corner { face: f, k }, depth);
                    for v in 0..n {
                        d[v] = d[v].min(d[u] + st[v]);
                    }
                }
            }
        }
    }
    d
}

fn assert_close(a: &[f64], b: &[f64], tol: f64, what: &str) {
    assert_eq!(a.len(), b.len(), "{what}");
    for (v, (x, y)) in a.iter().zip(b).enumerate() {
        assert!(
            (x - y).abs() <= tol * (1.0 + y.abs()),
            "{what}: vertex {v}: {x} != {y}"
        );
    }
}

// ----- tests -----

fn simplicial_cases() -> Vec<(Shape, Vec<usize>)> {
    vec![
        (
            Shape::generated("icosphere2", icosphere(2).unwrap()),
            vec![0],
        ),
        (
            Shape::generated("torus", torus(12, 7, 1.0, 0.4).unwrap()),
            vec![5, 20],
        ),
        (
            Shape::generated("bumpy", bumpy_grid_disk(14, 9, 0.3).unwrap()),
            vec![0, 40],
        ),
        (
            Shape::generated("annulus", annulus(16, 4, 0.5, 1.0).unwrap()),
            vec![2],
        ),
        (cube(3), vec![0]),
        (jittered_grid(10, 0.3, 5), vec![17]),
    ]
}

#[test]
fn simplicial_complex_equals_trimesh_entry_point_bitwise() {
    for (shape, sources) in simplicial_cases() {
        let old = exact_distance(&shape.mesh, &shape.geometry, &sources, &dist_opts()).unwrap();
        let complex = Complex::from_tri_mesh(&shape.mesh);
        let new = exact_distance_complex(
            &complex,
            shape.geometry.edge_lengths(),
            &vertex_sources(&sources),
            &dist_opts(),
        )
        .unwrap();
        let bits = |d: &[f64]| d.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&old.distances), bits(&new.distances), "{}", shape.name);
        assert_eq!(old.report, new.report, "{}", shape.name);
    }
}

/// Splits `face` of a surface at `bary` and returns the split surface and the new vertex.
fn split_at(s: &Surface, face: usize, bary: [f64; 3]) -> (Surface, usize) {
    let l0 = s.length(3 * face);
    let pos = [
        [0.0, 0.0],
        [l0, 0.0],
        apex(l0, s.length(3 * face + 1), s.length(3 * face + 2)),
    ];
    let p = [0, 1].map(|i| bary[0] * pos[0][i] + bary[1] * pos[1][i] + bary[2] * pos[2][i]);
    let mut t = s.clone();
    let split = t.c.split_face(face).unwrap();
    t.l.resize(t.c.n_edges(), 0.0);
    for k in 0..3 {
        t.l[split.new_edges[k] as usize] = norm(sub(pos[k], p));
    }
    t.orig.push(usize::MAX);
    (t, split.vertex as usize)
}

fn check_point_sources(s: &Surface, tol: f64, what: &str) {
    for face in 0..s.c.n_faces().min(12) {
        for bary in [[0.2, 0.3, 0.5], [0.6, 0.1, 0.3], [0.05, 0.9, 0.05]] {
            let direct = exact_distance_complex(
                &s.c,
                &s.l,
                &[ExactSource::Point { face, bary }],
                &dist_opts(),
            )
            .unwrap();
            let (t, p) = split_at(s, face, bary);
            let via = exact_distance_complex(&t.c, &t.l, &[ExactSource::Vertex(p)], &dist_opts())
                .unwrap();
            let n = s.c.n_vertices();
            assert_close(
                &direct.distances,
                &via.distances[..n],
                tol,
                &format!("{what} face {face}"),
            );
        }
    }
}

#[test]
fn face_point_source_equals_splitting_the_face_simplicial() {
    for (shape, _) in simplicial_cases() {
        let s = Surface::from_mesh(&shape.mesh, edge_lengths_of(&shape));
        check_point_sources(&s, 1e-11, &shape.name);
    }
}

/// A flat periodic `n × m` grid of unit squares with diagonals, and its edge lengths.
fn flat_torus(n: usize, m: usize) -> Surface {
    let id = |i: usize, j: usize| ((i % n) * m + (j % m)) as u32;
    let mut faces = Vec::new();
    for i in 0..n {
        for j in 0..m {
            faces.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            faces.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    let mesh = TriMesh::new(n * m, faces).unwrap();
    let l = mesh
        .edges()
        .iter()
        .map(|&[a, b]| {
            let (ai, aj, bi, bj) = (
                a as usize / m,
                a as usize % m,
                b as usize / m,
                b as usize % m,
            );
            let di = (bi + n - ai) % n;
            let dj = (bj + m - aj) % m;
            let step = |d: usize, k: usize| {
                if d == 0 {
                    0
                } else if d == 1 {
                    1
                } else {
                    assert_eq!(d, k - 1);
                    -1
                }
            };
            let (si, sj) = (step(di, n), step(dj, m));
            match (si, sj) {
                (0, 0) => unreachable!(),
                (0, _) | (_, 0) => 1.0,
                (1, 1) | (-1, -1) => 2f64.sqrt(),
                _ => unreachable!("anti-diagonal"),
            }
        })
        .collect();
    Surface::from_mesh(&mesh, l)
}

/// The distance between two lattice points on the torus `R² / (nZ × mZ)`.
fn lattice_distance(n: usize, m: usize, a: usize, b: usize) -> f64 {
    let (ai, aj, bi, bj) = (
        (a / m) as f64,
        (a % m) as f64,
        (b / m) as f64,
        (b % m) as f64,
    );
    let mut best = f64::INFINITY;
    for si in -2..=2 {
        for sj in -2..=2 {
            let dx = bi - ai + si as f64 * n as f64;
            let dy = bj - aj + sj as f64 * m as f64;
            best = best.min(dx.hypot(dy));
        }
    }
    best
}

#[test]
fn flat_torus_coarsened_to_few_vertices_matches_the_lattice_distance() {
    let (n, m) = (5, 4);
    let mut found_features = false;
    for target in [8, 4, 3, 2] {
        let mut s = flat_torus(n, m);
        s.coarsen(target);
        assert!(
            s.c.n_vertices() <= target.max(2) + 2,
            "target {target}: {}",
            s.c.n_vertices()
        );
        s.c.check().unwrap();
        found_features |= s.has_self_edge() || s.has_double_edge() || s.has_repeated_corner();
        for src in 0..s.c.n_vertices() {
            let d = exact_distance_complex(&s.c, &s.l, &[ExactSource::Vertex(src)], &dist_opts())
                .unwrap();
            for v in 0..s.c.n_vertices() {
                let e = lattice_distance(n, m, s.orig[src], s.orig[v]);
                assert!(
                    (d.distances[v] - e).abs() < 1e-9,
                    "target {target}: {src} -> {v}: {} != {e}",
                    d.distances[v]
                );
            }
        }
        check_point_sources(&s, 1e-11, &format!("flat torus {target}"));
    }
    assert!(
        found_features,
        "no doubled edge, self-edge or repeated corner was produced"
    );
}

#[test]
fn curved_complexes_coarsened_match_brute_force_unfolding() {
    let mut features = 0;
    let cases: Vec<(Shape, usize)> = vec![
        (
            Shape::generated("torus", torus(6, 4, 1.0, 0.45).unwrap()),
            3,
        ),
        (Shape::generated("icosphere1", icosphere(1).unwrap()), 5),
        (
            Shape::generated("torus", torus(5, 5, 1.0, 0.35).unwrap()),
            2,
        ),
    ];
    for (shape, target) in cases {
        let mut s = Surface::from_mesh(&shape.mesh, edge_lengths_of(&shape));
        s.coarsen(target);
        s.c.check().unwrap();
        if s.has_self_edge() || s.has_double_edge() || s.has_repeated_corner() {
            features += 1;
        }
        let n = s.c.n_vertices();
        let mut srcs = vec![vec![ExactSource::Vertex(0)]];
        if n > 1 {
            srcs.push(vec![ExactSource::Vertex(n - 1)]);
        }
        srcs.push(vec![ExactSource::Point {
            face: s.c.n_faces() - 1,
            bary: [0.3, 0.3, 0.4],
        }]);
        for sources in srcs {
            let d = exact_distance_complex(&s.c, &s.l, &sources, &dist_opts()).unwrap();
            let b = brute_force(&s, &sources, 12);
            assert_close(
                &d.distances,
                &b,
                1e-7,
                &format!("{} coarsened to {n}", shape.name),
            );
        }
    }
    assert!(
        features > 0,
        "the coarsened complexes have no Δ-complex features"
    );
}

#[test]
fn bad_input_is_an_error() {
    let shape = Shape::generated("icosphere1", icosphere(1).unwrap());
    let c = Complex::from_tri_mesh(&shape.mesh);
    let l = shape.geometry.edge_lengths().to_vec();
    let go = |l: &[f64], s: &[ExactSource]| exact_distance_complex(&c, l, s, &dist_opts());
    assert_eq!(go(&l, &[]).unwrap_err(), ExactError::NoSource);
    assert!(matches!(
        go(&l[1..], &vertex_sources(&[0])),
        Err(ExactError::Lengths(_))
    ));
    let mut bad = l.clone();
    bad[3] = f64::NAN;
    assert!(matches!(
        go(&bad, &vertex_sources(&[0])),
        Err(ExactError::Lengths(_))
    ));
    bad[3] = 0.0;
    assert!(matches!(
        go(&bad, &vertex_sources(&[0])),
        Err(ExactError::Lengths(_))
    ));
    bad[3] = 100.0;
    assert!(
        matches!(go(&bad, &vertex_sources(&[0])), Err(ExactError::Lengths(m)) if m.contains("not a triangle"))
    );
    assert!(matches!(
        go(&l, &vertex_sources(&[999])),
        Err(ExactError::SourceOutOfRange { .. })
    ));
    for point in [
        ExactSource::Point {
            face: 999,
            bary: [0.2, 0.3, 0.5],
        },
        ExactSource::Point {
            face: 0,
            bary: [0.5, 0.5, 0.5],
        },
        ExactSource::Point {
            face: 0,
            bary: [f64::NAN, 0.5, 0.5],
        },
        ExactSource::Point {
            face: 0,
            bary: [-0.5, 1.0, 0.5],
        },
    ] {
        assert!(
            matches!(go(&l, &[point]), Err(ExactError::InvalidPoint(_))),
            "{point:?}"
        );
    }
}

#[test]
fn point_on_an_edge_and_at_a_corner_work() {
    let shape = Shape::generated("icosphere2", icosphere(2).unwrap());
    let c = Complex::from_tri_mesh(&shape.mesh);
    let l = shape.geometry.edge_lengths();
    let base = exact_distance_complex(
        &c,
        l,
        &[ExactSource::Vertex(c.faces()[0][0] as usize)],
        &dist_opts(),
    )
    .unwrap();
    // At the corner: the same distances up to the nudge.
    let corner = exact_distance_complex(
        &c,
        l,
        &[ExactSource::Point {
            face: 0,
            bary: [1.0, 0.0, 0.0],
        }],
        &dist_opts(),
    )
    .unwrap();
    assert_close(&corner.distances, &base.distances, 1e-9, "corner");
    // On an edge: against the split of the edge's faces is not available; against the face
    // point just inside, both within the nudge.
    let on_edge = exact_distance_complex(
        &c,
        l,
        &[ExactSource::Point {
            face: 0,
            bary: [0.5, 0.5, 0.0],
        }],
        &dist_opts(),
    )
    .unwrap();
    let inside = exact_distance_complex(
        &c,
        l,
        &[ExactSource::Point {
            face: 0,
            bary: [0.5, 0.5 - 1e-9, 1e-9],
        }],
        &dist_opts(),
    )
    .unwrap();
    assert_close(&on_edge.distances, &inside.distances, 1e-7, "edge");
}
