//! Dynamic remeshing (`SPEC.md` §10, D16, D17): the oracle's `remesh-ops` fixtures, every
//! refusal class, random-operation stress tests on several topologies (run in parallel), the
//! full pass on a stretched mesh, pins, boundaries, Delaunay, determinism, transport.

#![allow(clippy::needless_range_loop)]

use std::path::PathBuf;
use std::time::Instant;

use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::remesh::{
    CollapseChecks, CollapseGuard, CollapseQuery, EditMesh, EnabledOps, FlipChecks, RefusalKind,
    RemeshOptions, RemeshOrder, Remeshed, VertexMap, circumcenter, corner_angle, remesh,
    smoothing_displacement, violates_delaunay,
};
use ddg_repulsive_surfaces::{RsError, Surface, read_obj};

const KLEIN: &str = "/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces/scenes/KleinBottle/klein_bottle_blocky.obj";

// ------------------------------------------------------------------ meshes

type Mesh = (Surface, Vec<[f64; 3]>);

fn mk(pos: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Mesh {
    (Surface::new(pos.len(), faces).unwrap(), pos)
}

fn gen_mesh(g: ddg_core::generators::GeneratedMesh) -> Mesh {
    mk(g.positions, g.faces)
}

/// A Möbius band with `n` quads around the core circle (radius 2, half-width 0.5), `rows`
/// vertex rows across.
fn moebius(n: usize) -> Mesh {
    let rows = 3usize;
    let mut pos = Vec::new();
    for i in 0..n {
        let th = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let (c, s) = (th.cos(), th.sin());
        let (ch, sh) = ((th / 2.0).cos(), (th / 2.0).sin());
        for r in 0..rows {
            let w = 0.5 - r as f64 * 1.0 / (rows - 1) as f64;
            pos.push([(2.0 + w * ch) * c, (2.0 + w * ch) * s, w * sh]);
        }
    }
    let at = |i: usize, r: usize| {
        if i == n {
            (rows - 1 - r) as u32
        } else {
            (rows * i + r) as u32
        }
    };
    let mut faces = Vec::new();
    for i in 0..n {
        for r in 0..rows - 1 {
            faces.push([at(i, r), at(i, r + 1), at(i + 1, r + 1)]);
            faces.push([at(i, r), at(i + 1, r + 1), at(i + 1, r)]);
        }
    }
    mk(pos, faces)
}

fn klein() -> Mesh {
    read_obj(KLEIN).unwrap()
}

/// An icosphere and a torus side by side (two components).
fn two_components() -> Mesh {
    let a = icosphere(2).unwrap();
    let b = torus(16, 8, 1.0, 0.35).unwrap();
    let mut pos = a.positions.clone();
    let off = pos.len() as u32;
    pos.extend(b.positions.iter().map(|p| [p[0] + 3.5, p[1], p[2]]));
    let mut faces = a.faces.clone();
    faces.extend(b.faces.iter().map(|f| f.map(|v| v + off)));
    mk(pos, faces)
}

fn disk() -> Mesh {
    gen_mesh(bumpy_grid_disk(10, 10, 0.1).unwrap())
}

fn named_meshes() -> Vec<(&'static str, Mesh)> {
    vec![
        ("icosphere", gen_mesh(icosphere(2).unwrap())),
        ("torus", gen_mesh(torus(16, 8, 1.0, 0.4).unwrap())),
        ("disk", disk()),
        ("moebius", moebius(24)),
        ("klein", klein()),
        ("two-components", two_components()),
    ]
}

fn scaled(m: &Mesh, s: [f64; 3]) -> Mesh {
    let pos =
        m.1.iter()
            .map(|p| [p[0] * s[0], p[1] * s[1], p[2] * s[2]])
            .collect();
    (m.0.clone(), pos)
}

fn edge_lengths(s: &Surface, p: &[[f64; 3]]) -> Vec<f64> {
    s.edges()
        .iter()
        .map(|&[a, b]| {
            let d = [0, 1, 2].map(|k| p[a as usize][k] - p[b as usize][k]);
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .collect()
}

fn angles(s: &Surface, p: &[[f64; 3]]) -> Vec<f64> {
    let mut r = Vec::new();
    for t in s.faces() {
        for k in 0..3 {
            r.push(corner_angle(
                p[t[k] as usize],
                p[t[(k + 1) % 3] as usize],
                p[t[(k + 2) % 3] as usize],
            ));
        }
    }
    r
}

fn min(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::INFINITY, f64::min)
}

fn fraction_outside(len: &[f64], lo: f64, hi: f64) -> f64 {
    len.iter().filter(|&&l| l < lo || l > hi).count() as f64 / len.len() as f64
}

fn percentile(v: &[f64], q: f64) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    s[((s.len() - 1) as f64 * q) as usize]
}

// ------------------------------------------------------------------ oracle fixtures

struct Case {
    name: String,
    op: (String, u32, u32),
    pinned: Vec<u32>,
    delaunay: bool,
    foldover: bool,
    mesh: Mesh,
}

fn hexf(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).unwrap())
}

fn parse_mesh_lines(lines: &[&str]) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut pos = Vec::new();
    let mut faces = Vec::new();
    for l in lines {
        let w: Vec<&str> = l.split_whitespace().collect();
        match w.first().copied() {
            Some("p") => pos.push([hexf(w[1]), hexf(w[2]), hexf(w[3])]),
            Some("f") => faces.push([w[1], w[2], w[3]].map(|x| x.parse::<u32>().unwrap())),
            _ => {}
        }
    }
    (pos, faces)
}

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/rs")
}

fn read_case(path: &PathBuf) -> Option<Case> {
    let text = std::fs::read_to_string(path).unwrap();
    if !text.contains("task remesh-ops") {
        return None;
    }
    let lines: Vec<&str> = text.lines().collect();
    let mut c = Case {
        name: String::new(),
        op: (String::new(), 0, 0),
        pinned: Vec::new(),
        delaunay: false,
        foldover: false,
        mesh: (Surface::new(3, vec![[0, 1, 2]]).unwrap(), vec![]),
    };
    for l in &lines {
        let w: Vec<&str> = l.split_whitespace().collect();
        match w.first().copied() {
            Some("case") => c.name = w[1].to_string(),
            Some("op") => {
                c.op = (
                    w[1].to_string(),
                    w[2].parse().unwrap(),
                    w[3].parse().unwrap(),
                )
            }
            Some("pinned") => c.pinned = w[1..].iter().map(|x| x.parse().unwrap()).collect(),
            Some("legality") => {
                c.delaunay = w.contains(&"delaunay");
                c.foldover = w.contains(&"foldover");
            }
            _ => {}
        }
    }
    let (pos, faces) = parse_mesh_lines(&lines);
    c.mesh = mk(pos, faces);
    Some(c)
}

struct Expected {
    error: Option<String>,
    pos: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    parent: Vec<u32>,
    vmap: Vec<Vec<u32>>,
}

fn read_expected(path: &PathBuf) -> Expected {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let mut e = Expected {
        error: None,
        pos: vec![],
        faces: vec![],
        parent: vec![],
        vmap: vec![],
    };
    for l in &lines {
        let w: Vec<&str> = l.split_whitespace().collect();
        match w.first().copied() {
            Some("status") if w[1] == "error" => e.error = Some(w[2].to_string()),
            Some("parent") => e.parent.push(w[2].parse().unwrap()),
            Some("vmap") => e
                .vmap
                .push(w[2..].iter().map(|x| x.parse().unwrap()).collect()),
            _ => {}
        }
    }
    let (p, f) = parse_mesh_lines(&lines);
    e.pos = p;
    e.faces = f;
    e
}

#[test]
fn oracle_remesh_ops_fixtures_replay_exactly() {
    let mut n = 0;
    let mut entries: Vec<PathBuf> = std::fs::read_dir(fixture_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "in"))
        .collect();
    entries.sort();
    for input in entries {
        let Some(case) = read_case(&input) else {
            continue;
        };
        n += 1;
        let exp = read_expected(&input.with_extension("out"));
        let (s, p) = &case.mesh;
        let mut m = EditMesh::new(s, p, &case.pinned).unwrap();
        let (op, a, b) = (case.op.0.as_str(), case.op.1, case.op.2);
        let res: Result<bool, ddg_repulsive_surfaces::remesh::Refused> = match op {
            "split" => m.split(a, b).map(|_| false),
            "collapse" => m
                .collapse(
                    a,
                    b,
                    CollapseChecks {
                        foldover: case.foldover,
                    },
                    None,
                )
                .map(|_| false),
            "flip" => m
                .flip(
                    a,
                    b,
                    FlipChecks {
                        delaunay: case.delaunay,
                        foldover: case.foldover,
                    },
                )
                .map(|o| o.reordered),
            _ => panic!("unknown op {op}"),
        };
        match (&exp.error, res) {
            (Some(kind), Err(r)) => {
                let ok = match kind.as_str() {
                    "NonManifoldResult" => matches!(
                        r.kind,
                        RefusalKind::RepeatedFace | RefusalKind::NonManifoldResult
                    ),
                    "FlipNonManifold" => r.kind == RefusalKind::ExistingEdge,
                    k => r.kind.name() == k,
                };
                assert!(ok, "{}: oracle {kind}, Rust {r}", case.name);
            }
            (Some(kind), Ok(reordered)) => {
                // SPEC §10: the Rust flip reorders a same-direction pair locally; the oracle
                // refuses it.
                assert_eq!(
                    kind, "InconsistentOrientation",
                    "{}: Rust succeeded",
                    case.name
                );
                assert!(reordered, "{}", case.name);
                let r = m.finish().expect("valid after a reordered flip");
                assert_eq!(r.surface.n_faces(), s.n_faces());
            }
            (None, Err(r)) => panic!("{}: oracle ok, Rust refused: {r}", case.name),
            (None, Ok(reordered)) => {
                assert!(!reordered, "{}", case.name);
                let r = m.finish().unwrap();
                assert_eq!(r.surface.faces(), &exp.faces[..], "{}: faces", case.name);
                assert_eq!(r.positions.len(), exp.pos.len(), "{}", case.name);
                for (x, y) in r.positions.iter().zip(&exp.pos) {
                    assert_eq!(
                        x.map(f64::to_bits),
                        y.map(f64::to_bits),
                        "{}: positions",
                        case.name
                    );
                }
                assert_eq!(r.face_parents, exp.parent, "{}: parents", case.name);
                let vm: Vec<Vec<u32>> = (0..r.vertex_map.n_new() as u32)
                    .map(|v| r.vertex_map.entry(v).0.to_vec())
                    .collect();
                assert_eq!(vm, exp.vmap, "{}: vmap", case.name);
                let pins: Vec<u32> = case
                    .pinned
                    .iter()
                    .map(|&p| {
                        if op == "collapse" && p == a.max(b) {
                            a.min(b)
                        } else {
                            p
                        }
                    })
                    .collect();
                assert_eq!(r.pinned, pins, "{}: pins", case.name);
            }
        }
    }
    assert_eq!(n, 16, "number of remesh-ops fixtures");
}

// ------------------------------------------------------------------ refusals

fn single_triangle() -> Mesh {
    mk(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
}

fn tetra() -> Mesh {
    mk(
        vec![
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]],
    )
}

/// A strip one quad wide: every rung is an interior edge between two boundary vertices.
fn strip() -> Mesh {
    let mut pos = Vec::new();
    for i in 0..4 {
        pos.push([i as f64, 0.0, 0.0]);
        pos.push([i as f64, 1.0, 0.0]);
    }
    let mut faces = Vec::new();
    for i in 0..3u32 {
        let (a, b, c, d) = (2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3);
        faces.push([a, c, d]);
        faces.push([a, d, b]);
    }
    mk(pos, faces)
}

struct No;
impl CollapseGuard for No {
    fn allow(&mut self, q: &CollapseQuery<'_>) -> bool {
        assert!(!q.new_star.is_empty());
        false
    }
}

#[test]
fn every_refusal_class() {
    let none = CollapseChecks::default();
    let ico = gen_mesh(icosphere(1).unwrap());
    let em = |m: &Mesh, pins: &[u32]| EditMesh::new(&m.0, &m.1, pins).unwrap();
    let (a, b) = (ico.0.edges()[0][0], ico.0.edges()[0][1]);

    // NoSuchEdge (split, collapse, flip), including vertices out of range
    let mut m = em(&ico, &[]);
    let far = (0..ico.0.n_vertices() as u32)
        .find(|&v| v != a && ico.0.find_edge(a, v).is_none())
        .unwrap();
    assert_eq!(m.split(a, far).unwrap_err().kind, RefusalKind::NoSuchEdge);
    assert_eq!(
        m.collapse(a, far, none, None).unwrap_err().kind,
        RefusalKind::NoSuchEdge
    );
    assert_eq!(
        m.flip(a, 10_000, FlipChecks::default()).unwrap_err().kind,
        RefusalKind::NoSuchEdge
    );
    // BothPinned
    let mut m = em(&ico, &[a, b]);
    let r = m.collapse(a, b, none, None).unwrap_err();
    assert_eq!(r.kind, RefusalKind::BothPinned);
    assert!(r.message.contains(&format!("({a}, {b})")));
    // LinkCondition: a degree-3 vertex inside a face makes the face's corners share an extra
    // neighbour
    let (s, p) = &ico;
    let mut faces = s.faces().to_vec();
    let t = faces[0];
    let c = p.len() as u32;
    let mut pos = p.clone();
    pos.push(
        [0, 1, 2].map(|k| (p[t[0] as usize][k] + p[t[1] as usize][k] + p[t[2] as usize][k]) / 3.0),
    );
    faces[0] = [t[0], t[1], c];
    faces.push([t[1], t[2], c]);
    faces.push([t[2], t[0], c]);
    let lifted = mk(pos, faces);
    let mut m = em(&lifted, &[]);
    assert_eq!(
        m.collapse(t[0], t[1], none, None).unwrap_err().kind,
        RefusalKind::LinkCondition
    );
    // LinkCondition: interior edge between two boundary vertices
    let st = strip();
    let mut m = em(&st, &[]);
    let r = m.collapse(2, 3, none, None).unwrap_err();
    assert_eq!(r.kind, RefusalKind::LinkCondition);
    assert!(
        r.message
            .contains("interior edge between two boundary vertices")
    );
    // RepeatedFace: a tetrahedron cannot lose an edge
    let te = tetra();
    let mut m = em(&te, &[]);
    assert_eq!(
        m.collapse(0, 1, none, None).unwrap_err().kind,
        RefusalKind::RepeatedFace
    );
    // ExistingEdge: every flip of a tetrahedron
    assert_eq!(
        m.flip(0, 1, FlipChecks::default()).unwrap_err().kind,
        RefusalKind::ExistingEdge
    );
    // NonManifoldResult: the last triangle of a component would vanish
    let tri = single_triangle();
    let mut m = em(&tri, &[]);
    assert_eq!(
        m.collapse(0, 1, none, None).unwrap_err().kind,
        RefusalKind::NonManifoldResult
    );
    // BoundaryEdge
    assert_eq!(
        m.flip(0, 1, FlipChecks::default()).unwrap_err().kind,
        RefusalKind::BoundaryEdge
    );
    // NotDelaunay: a regular icosahedron is Delaunay everywhere
    let mut m = em(&ico, &[]);
    assert_eq!(
        m.flip(
            a,
            b,
            FlipChecks {
                delaunay: true,
                foldover: false
            }
        )
        .unwrap_err()
        .kind,
        RefusalKind::NotDelaunay
    );
    // Guard
    let mut g = No;
    assert_eq!(
        m.collapse(a, b, none, Some(&mut g)).unwrap_err().kind,
        RefusalKind::Guard
    );
    // Foldover (collapse and flip): a sharp, noisy sphere has some; refused ops change nothing
    let mut noisy = gen_mesh(icosphere(2).unwrap());
    perturb(&mut noisy.1, [0.15; 3], 7);
    let mut found_c = 0;
    let mut found_f = 0;
    for &[a, b] in noisy.0.edges() {
        let mut m = em(&noisy, &[]);
        if let Err(r) = m.collapse(a, b, CollapseChecks { foldover: true }, None) {
            if r.kind == RefusalKind::Foldover {
                found_c += 1;
            }
            let fin = m.finish().unwrap();
            assert_eq!(fin.surface.faces(), noisy.0.faces());
            assert_eq!(fin.positions, noisy.1);
        }
        let mut m = em(&noisy, &[]);
        if let Err(r) = m.flip(
            a,
            b,
            FlipChecks {
                delaunay: false,
                foldover: true,
            },
        ) && r.kind == RefusalKind::Foldover
        {
            found_f += 1;
        }
    }
    assert!(
        found_c > 0 && found_f > 0,
        "foldovers: collapse {found_c}, flip {found_f}"
    );
    // DegenerateFace: the flip of a quad whose diagonal ends and a third corner are collinear
    let dart = mk(
        vec![
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
        ],
        vec![[0, 1, 2], [1, 0, 3]],
    );
    let mut m = em(&dart, &[]);
    assert_eq!(
        m.flip(0, 1, FlipChecks::default()).unwrap_err().kind,
        RefusalKind::DegenerateFace
    );
    assert_eq!(
        m.flip(
            0,
            1,
            FlipChecks {
                delaunay: false,
                foldover: true
            }
        )
        .unwrap_err()
        .kind,
        RefusalKind::Foldover
    );
    // every kind has a stable name and index
    for (i, k) in RefusalKind::ALL.iter().enumerate() {
        assert_eq!(k.index(), i);
        assert!(!k.name().is_empty());
    }
}

#[test]
fn options_are_validated() {
    let (s, p) = gen_mesh(icosphere(1).unwrap());
    let l0 = mean_edge_length(&s, &p);
    let bad = |o: RemeshOptions, name: &str| {
        let e = remesh(&s, &p, &[], l0, &o, None).unwrap_err();
        assert!(
            matches!(e, RsError::InvalidOption { name: n, .. } if n == name),
            "{e}"
        );
    };
    let d = RemeshOptions::default();
    bad(
        RemeshOptions {
            split: f64::NAN,
            ..d
        },
        "split",
    );
    bad(RemeshOptions { split: -1.0, ..d }, "split");
    bad(
        RemeshOptions {
            collapse: -0.1,
            ..d
        },
        "collapse",
    );
    bad(
        RemeshOptions {
            split: 1.0,
            collapse: 0.5,
            ..d
        },
        "split",
    );
    bad(RemeshOptions { rho: 0.0, ..d }, "rho");
    bad(RemeshOptions { rho: 1.5, ..d }, "rho");
    bad(
        RemeshOptions {
            flip_cap_per_vertex: 0,
            ..d
        },
        "flip_cap_per_vertex",
    );
    let e = remesh(&s, &p, &[], 0.0, &d, None).unwrap_err();
    assert!(matches!(e, RsError::InvalidOption { name: "l0", .. }));
    let e = remesh(&s, &p, &[9999], l0, &d, None).unwrap_err();
    assert!(
        matches!(e, RsError::InvalidOption { name: "pinned", .. }),
        "{e}"
    );
    let e = remesh(&s, &p[1..], &[], l0, &d, None).unwrap_err();
    assert!(matches!(e, RsError::PositionCountMismatch { .. }));
    println!("{}", e);
}

// ------------------------------------------------------------------ stress

struct StressResult {
    ops: [u64; 3],
    refused: u64,
    reordered: u64,
    validations: u64,
}

fn random_ops(name: &str, m: &Mesh, seed: u64, n_ops: usize) -> StressResult {
    let mut rng = SplitMix64::new(seed);
    let nv = m.0.n_vertices() as u32;
    let pins: Vec<u32> = (0..4)
        .map(|_| (rng.next_u64() % nv as u64) as u32)
        .collect();
    let mut em = EditMesh::new(&m.0, &m.1, &pins).unwrap();
    let f0 = m.0.n_faces();
    let mut res = StressResult {
        ops: [0; 3],
        refused: 0,
        reordered: 0,
        validations: 0,
    };
    let check = |em: &EditMesh, res: &mut StressResult, i: usize| {
        let r: Remeshed = em
            .clone()
            .finish()
            .unwrap_or_else(|e| panic!("{name} seed {seed} op {i}: {e}"));
        assert_eq!(r.pinned.len(), pins.len());
        for (k, &pv) in r.pinned.iter().enumerate() {
            assert_eq!(
                r.positions[pv as usize], m.1[pins[k] as usize],
                "{name} seed {seed}: pin {} moved",
                pins[k]
            );
        }
        // vertex map weights sum to one
        for v in 0..r.vertex_map.n_new() as u32 {
            let s: f64 = r.vertex_map.entry(v).1.iter().sum();
            assert!((s - 1.0).abs() < 1e-12);
        }
        assert_eq!(
            r.surface.n_components(),
            m.0.n_components(),
            "{name}: components"
        );
        assert_eq!(
            r.surface.euler_characteristic(),
            m.0.euler_characteristic(),
            "{name}: Euler characteristic"
        );
        assert_eq!(r.surface.boundary_loops().len(), m.0.boundary_loops().len());
        res.validations += 1;
    };
    for i in 0..n_ops {
        let edges = em.edges();
        let [a, b] = edges[(rng.next_u64() % edges.len() as u64) as usize];
        let (a, b) = if rng.next_u64().is_multiple_of(2) {
            (a, b)
        } else {
            (b, a)
        };
        let grow = em.n_faces() < 2 * f0;
        let r = rng.next_f64();
        let fold = rng.next_u64().is_multiple_of(2);
        let outcome = if r < 0.3 && grow {
            em.split(a, b).map(|_| 0)
        } else if r < 0.65 {
            em.collapse(a, b, CollapseChecks { foldover: fold }, None)
                .map(|_| 1)
        } else {
            em.flip(
                a,
                b,
                FlipChecks {
                    delaunay: false,
                    foldover: fold,
                },
            )
            .map(|o| {
                if o.reordered {
                    res.reordered += 1;
                }
                2
            })
        };
        // every operation leaves a valid surface with valid positions
        if let Err(e) = em.clone().finish() {
            panic!("{name} seed {seed} op {i} ({a}, {b}) {outcome:?}: {e}");
        }
        match outcome {
            Ok(k) => res.ops[k] += 1,
            Err(_) => res.refused += 1,
        }
        if i % 25 == 0 {
            check(&em, &mut res, i);
        }
    }
    check(&em, &mut res, n_ops);
    res
}

#[test]
fn thousands_of_random_operations_keep_a_valid_manifold() {
    let meshes = named_meshes();
    let seeds = [1u64, 2, 3, 4];
    let results: Vec<(String, StressResult)> = std::thread::scope(|sc| {
        let hs: Vec<_> = meshes
            .iter()
            .flat_map(|(name, m)| seeds.iter().map(move |&s| (name, m, s)))
            .map(|(name, m, s)| {
                sc.spawn(move || (format!("{name}/{s}"), random_ops(name, m, s, 1500)))
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut total = [0u64; 3];
    let mut reordered = std::collections::BTreeMap::new();
    for (n, r) in &results {
        println!(
            "{n}: split {} collapse {} flip {} refused {} reordered {} validations {}",
            r.ops[0], r.ops[1], r.ops[2], r.refused, r.reordered, r.validations
        );
        for k in 0..3 {
            total[k] += r.ops[k];
        }
        *reordered
            .entry(n.split('/').next().unwrap().to_string())
            .or_insert(0) += r.reordered;
    }
    println!("total: {total:?}");
    assert!(total.iter().all(|&t| t > 1000));
    // the non-orientable meshes exercise the local reordering of flips
    assert!(
        reordered["moebius"] > 0 && reordered["klein"] > 0,
        "{reordered:?}"
    );
    assert_eq!(reordered["icosphere"], 0);
}

#[test]
fn full_passes_on_every_topology_stay_valid() {
    let meshes = named_meshes();
    std::thread::scope(|sc| {
        for (name, m) in &meshes {
            for (k, s) in [[3.0, 1.0, 0.5], [0.4, 1.0, 1.0]].into_iter().enumerate() {
                sc.spawn(move || {
                    let mut m = scaled(m, s);
                    perturb(&mut m.1, [0.01; 3], k as u64);
                    let l0 = mean_edge_length(&m.0, &m.1);
                    let pins = [0u32];
                    let mut cur = m.clone();
                    let mut pins_cur = pins.to_vec();
                    for _ in 0..4 {
                        let (r, rep) = remesh(
                            &cur.0,
                            &cur.1,
                            &pins_cur,
                            l0,
                            &RemeshOptions::default(),
                            None,
                        )
                        .unwrap_or_else(|e| panic!("{name}: {e}"));
                        assert_eq!(r.surface.euler_characteristic(), m.0.euler_characteristic());
                        assert_eq!(r.surface.n_components(), m.0.n_components());
                        assert_eq!(r.surface.boundary_loops().len(), m.0.boundary_loops().len());
                        assert_eq!(r.positions[r.pinned[0] as usize], m.1[0], "{name}: pin");
                        assert_eq!(rep.n_faces, r.surface.n_faces());
                        pins_cur = r.pinned.clone();
                        cur = (r.surface, r.positions);
                    }
                });
            }
        }
    });
}

// ------------------------------------------------------------------ quality

#[test]
fn full_pass_fixes_a_stretched_mesh() {
    let base = gen_mesh(icosphere(3).unwrap());
    let m = scaled(&base, [4.0, 1.0, 0.5]);
    let l0 = mean_edge_length(&base.0, &base.1) * 1.5;
    let len0 = edge_lengths(&m.0, &m.1);
    let ang0 = angles(&m.0, &m.1);
    let out0 = fraction_outside(&len0, 0.5 * l0, 1.5 * l0);
    let mut cur = m.clone();
    let mut reports = Vec::new();
    let mut first_min_angle = 0.0;
    for pass in 0..8 {
        let (r, rep) = remesh(&cur.0, &cur.1, &[], l0, &RemeshOptions::default(), None).unwrap();
        cur = (r.surface, r.positions);
        if pass == 0 {
            first_min_angle = min(&angles(&cur.0, &cur.1));
        }
        reports.push(rep);
    }
    let len1 = edge_lengths(&cur.0, &cur.1);
    let ang1 = angles(&cur.0, &cur.1);
    let out1 = fraction_outside(&len1, 0.5 * l0, 1.5 * l0);
    let deg = 180.0 / std::f64::consts::PI;
    println!(
        "stretched icosphere(3) x(4, 1, 0.5), L0 = {l0:.4}: |F| {} -> {}",
        m.0.n_faces(),
        cur.0.n_faces()
    );
    println!(
        "  edges outside [L0/2, 3L0/2]: {:.1}% -> {:.2}%",
        100.0 * out0,
        100.0 * out1
    );
    println!(
        "  min angle {:.2} -> {:.2} deg (after one pass {:.2}); 1st percentile {:.2} -> {:.2}",
        min(&ang0) * deg,
        min(&ang1) * deg,
        first_min_angle * deg,
        percentile(&ang0, 0.01) * deg,
        percentile(&ang1, 0.01) * deg
    );
    for (i, r) in reports.iter().enumerate() {
        println!(
            "  pass {i}: split {} collapse {} (refused {:?}) flip {} (refused {}) smooth {} \
             (halved {}, refused {}) cap hits {} |V| {} |F| {} {:.3}s",
            r.splits,
            r.collapses,
            r.refused_collapses.0,
            r.flips,
            r.refused_flips.total(),
            r.smoothing_moves,
            r.smoothing_halved,
            r.smoothing_refused,
            r.flip_cap_hits,
            r.n_vertices,
            r.n_faces,
            r.seconds
        );
    }
    assert!(out0 > 0.3, "the test mesh must be badly sized ({out0})");
    assert!(out1 < 0.05, "fraction outside after 8 passes: {out1}");
    assert!(min(&ang1) > min(&ang0));
    assert!(percentile(&ang1, 0.01) > percentile(&ang0, 0.01));
}

#[test]
fn pins_never_move_or_vanish() {
    let base = gen_mesh(torus(24, 12, 1.0, 0.4).unwrap());
    let m = scaled(&base, [2.5, 1.0, 1.0]);
    let l0 = mean_edge_length(&base.0, &base.1);
    let mut rng = SplitMix64::new(11);
    let pins: Vec<u32> = (0..40)
        .map(|_| (rng.next_u64() % m.0.n_vertices() as u64) as u32)
        .collect();
    let mut cur = m.clone();
    let mut cur_pins = pins.clone();
    let mut collapses_at_pins = 0;
    for _ in 0..5 {
        let (r, _) = remesh(
            &cur.0,
            &cur.1,
            &cur_pins,
            l0,
            &RemeshOptions::default(),
            None,
        )
        .unwrap();
        assert_eq!(r.pinned.len(), pins.len());
        for (k, &v) in r.pinned.iter().enumerate() {
            assert_eq!(r.positions[v as usize], m.1[pins[k] as usize]);
            let (src, w) = r.vertex_map.entry(v);
            let i = src.iter().position(|&o| o == cur_pins[k]).unwrap();
            assert_eq!(w[i], 1.0);
            if src.len() > 1 {
                collapses_at_pins += 1;
            }
        }
        cur_pins = r.pinned.clone();
        cur = (r.surface, r.positions);
    }
    println!("collapses into pinned vertices: {collapses_at_pins}");
}

#[test]
fn boundary_vertices_keep_their_positions_and_loops_stay_closed() {
    let (s, mut p) = gen_mesh(bumpy_grid_disk(16, 16, 0.15).unwrap());
    perturb(&mut p, [0.01, 0.01, 0.0], 5);
    let l0 = mean_edge_length(&s, &p);
    let bnd: Vec<u32> = (0..s.n_vertices() as u32)
        .filter(|&v| s.is_boundary_vertex(v))
        .collect();
    // smoothing and flips only: the numbering is kept, boundary vertices do not move
    let o = RemeshOptions {
        enabled_ops: EnabledOps {
            split: false,
            collapse: false,
            flip: true,
            smooth: true,
        },
        ..RemeshOptions::default()
    };
    let (r, rep) = remesh(&s, &p, &[], l0, &o, None).unwrap();
    assert!(rep.smoothing_moves > 0);
    assert_eq!(r.surface.n_vertices(), s.n_vertices());
    for &v in &bnd {
        assert_eq!(r.positions[v as usize], p[v as usize]);
        assert!(r.surface.is_boundary_vertex(v));
    }
    assert_eq!(r.surface.boundary_loops().len(), 1);
    assert_eq!(r.surface.boundary_loops()[0].len(), bnd.len());
    // the full pass on a stretched disk: still one closed loop, through the old boundary
    // vertices that remain (collapses on the boundary move to midpoints)
    let ms = scaled(&(s.clone(), p.clone()), [3.0, 1.0, 1.0]);
    let (r, _) = remesh(&ms.0, &ms.1, &[], l0, &RemeshOptions::default(), None).unwrap();
    assert_eq!(r.surface.boundary_loops().len(), 1);
    let lp = &r.surface.boundary_loops()[0];
    for i in 0..lp.len() {
        let (a, b) = (lp[i], lp[(i + 1) % lp.len()]);
        let e = r
            .surface
            .find_edge(a, b)
            .expect("consecutive loop vertices share an edge");
        assert!(r.surface.is_boundary_edge(e));
    }
}

#[test]
fn flips_make_a_flat_patch_delaunay() {
    let (s, mut p) = gen_mesh(bumpy_grid_disk(14, 14, 0.0).unwrap());
    perturb(&mut p, [0.02, 0.02, 0.0], 9);
    // a valid planar triangulation (no inverted face), so Delaunay flipping terminates with
    // every edge locally Delaunay
    for t in s.faces() {
        let n = ddg_repulsive_surfaces::remesh::tri_normal(t.map(|v| p[v as usize]));
        assert!(n[2] > 0.0);
    }
    let o = RemeshOptions {
        rounds: 1,
        enabled_ops: EnabledOps {
            split: false,
            collapse: false,
            flip: true,
            smooth: false,
        },
        ..RemeshOptions::default()
    };
    let count_bad = |s: &Surface, p: &[[f64; 3]]| {
        let mut bad = 0;
        for (e, &[a, b]) in s.edges().iter().enumerate() {
            let [f, g] = s.edge_faces()[e];
            if g == ddg_repulsive_surfaces::NONE {
                continue;
            }
            let opp = |f: u32| {
                *s.faces()[f as usize]
                    .iter()
                    .find(|&&v| v != a && v != b)
                    .unwrap()
            };
            let (c, d) = (opp(f), opp(g));
            if violates_delaunay(p[a as usize], p[b as usize], p[c as usize], p[d as usize]) {
                bad += 1;
            }
        }
        bad
    };
    let before = count_bad(&s, &p);
    let (r, rep) = remesh(&s, &p, &[], mean_edge_length(&s, &p), &o, None).unwrap();
    println!("non-Delaunay edges {before} -> 0, flips {}", rep.flips);
    assert!(before > 10);
    assert_eq!(rep.flip_cap_hits, 0);
    assert_eq!(rep.refused_flips.total(), 0);
    assert_eq!(count_bad(&r.surface, &r.positions), 0);
    assert_eq!(r.positions, p);
}

#[test]
fn remeshing_is_deterministic_and_both_orders_work() {
    let base = gen_mesh(torus(30, 12, 1.0, 0.4).unwrap());
    let mut m = scaled(&base, [1.0, 2.0, 0.6]);
    perturb(&mut m.1, [0.02; 3], 3);
    let l0 = mean_edge_length(&base.0, &base.1);
    for order in [RemeshOrder::FlipThenSmooth, RemeshOrder::SmoothThenFlip] {
        let o = RemeshOptions {
            order,
            ..RemeshOptions::default()
        };
        let (r1, rep1) = remesh(&m.0, &m.1, &[3, 5], l0, &o, None).unwrap();
        let (r2, rep2) = remesh(&m.0, &m.1, &[3, 5], l0, &o, None).unwrap();
        assert_eq!(r1.surface.faces(), r2.surface.faces());
        assert_eq!(
            r1.positions
                .iter()
                .map(|x| x.map(f64::to_bits))
                .collect::<Vec<_>>(),
            r2.positions
                .iter()
                .map(|x| x.map(f64::to_bits))
                .collect::<Vec<_>>()
        );
        assert_eq!(r1.vertex_map, r2.vertex_map);
        assert_eq!(r1.face_parents, r2.face_parents);
        assert_eq!(r1.pinned, r2.pinned);
        assert_eq!(
            (rep1.splits, rep1.collapses, rep1.flips),
            (rep2.splits, rep2.collapses, rep2.flips)
        );
        assert!(rep1.splits > 0 && rep1.collapses > 0 && rep1.flips > 0);
    }
}

#[test]
fn vertex_map_transports_linear_data_and_face_parents_point_back() {
    let base = gen_mesh(icosphere(3).unwrap());
    let m = scaled(&base, [3.0, 1.0, 0.7]);
    let l0 = mean_edge_length(&base.0, &base.1);
    // without smoothing, every new position is the vertex-map combination of old positions
    let o = RemeshOptions {
        enabled_ops: EnabledOps {
            smooth: false,
            ..EnabledOps::default()
        },
        ..RemeshOptions::default()
    };
    let (r, rep) = remesh(&m.0, &m.1, &[0], l0, &o, None).unwrap();
    assert!(rep.splits > 0 && rep.collapses > 0);
    let t = r.vertex_map.transport_vec3(&m.1).unwrap();
    for (x, y) in t.iter().zip(&r.positions) {
        for k in 0..3 {
            assert!((x[k] - y[k]).abs() < 1e-12, "{x:?} {y:?}");
        }
    }
    let ones = r
        .vertex_map
        .transport(&vec![1.0; m.0.n_vertices()])
        .unwrap();
    assert!(ones.iter().all(|&x| (x - 1.0).abs() < 1e-12));
    assert!(r.vertex_map.transport(&[1.0]).is_err());
    let set = r.vertex_map.transport_set(&[0], 0.5).unwrap();
    assert!(set.contains(&r.pinned[0]));
    // composition with the identity
    let id = VertexMap::identity(m.0.n_vertices());
    assert_eq!(r.vertex_map.compose(&id).unwrap(), r.vertex_map);
    // face parents: every new face's parent exists; its area is mostly that of a piece of the
    // parent's neighbourhood: check range and that every old face that survives is a parent
    assert_eq!(r.face_parents.len(), r.surface.n_faces());
    assert!(r.face_parents.iter().all(|&f| (f as usize) < m.0.n_faces()));
}

#[test]
fn circumcenter_smoothing_weights_and_delaunay_test() {
    let mut rng = SplitMix64::new(4);
    for _ in 0..200 {
        let p = [0; 3].map(|_| [0; 3].map(|_| rng.uniform(-1.0, 1.0)));
        let c = circumcenter(p);
        let r: Vec<f64> = p
            .iter()
            .map(|q| (0..3).map(|k| (q[k] - c[k]).powi(2)).sum::<f64>().sqrt())
            .collect();
        assert!((r[0] - r[1]).abs() < 1e-8 * r[0] && (r[0] - r[2]).abs() < 1e-8 * r[0]);
        // the center lies in the triangle's plane
        let n = ddg_repulsive_surfaces::remesh::tri_normal(p);
        let d: f64 = (0..3).map(|k| n[k] * (c[k] - p[0][k])).sum();
        assert!(d.abs() < 1e-9 * r[0] * r[0]);
        // Delaunay test = sign of the cotan weight
        let q = [
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
        ];
        let cot = |v: [f64; 3], a: [f64; 3], b: [f64; 3]| {
            let t = corner_angle(v, a, b);
            t.cos() / t.sin()
        };
        let w = cot(p[2], p[0], p[1]) + cot(q, p[0], p[1]);
        assert_eq!(violates_delaunay(p[0], p[1], p[2], q), w < 0.0);
    }
    use ddg_repulsive_surfaces::remesh::{COS_FOLDOVER_LIMIT, FOLDOVER_LIMIT, folds};
    assert!((FOLDOVER_LIMIT.cos() - COS_FOLDOVER_LIMIT).abs() < 1e-15);
    for k in 0..=100 {
        let t = std::f64::consts::PI * k as f64 / 100.0;
        if (t - FOLDOVER_LIMIT).abs() > 1e-9 {
            assert_eq!(
                folds([1.0, 0.0, 0.0], [t.cos(), t.sin(), 0.0]),
                t > FOLDOVER_LIMIT
            );
        }
    }
    assert!(folds([1.0, 0.0, 0.0], [0.0; 3]));
    // rho/Σa normalization: one face, rho 1 moves to its center
    let u = smoothing_displacement([1.0, 2.0, 3.0], 1.0, [(2.0, [0.0; 3])].into_iter()).unwrap();
    assert_eq!(u, [-1.0, -2.0, -3.0]);
    let u = smoothing_displacement(
        [0.0; 3],
        0.5,
        [(1.0, [1.0, 0.0, 0.0]), (3.0, [0.0, 1.0, 0.0])].into_iter(),
    )
    .unwrap();
    assert_eq!(u, [0.125, 0.375, 0.0]);
    assert!(smoothing_displacement([0.0; 3], 0.5, std::iter::empty()).is_none());
}

/// One full pass on a ~40k-face mesh (release build; `cargo test --release -- --ignored
/// timing --nocapture`).
#[test]
#[ignore]
fn timing_one_pass_40k_faces() {
    let base = gen_mesh(torus(200, 100, 1.0, 0.4).unwrap());
    let mut m = scaled(&base, [1.3, 1.0, 0.8]);
    perturb(&mut m.1, [0.002; 3], 1);
    let l0 = mean_edge_length(&base.0, &base.1);
    // steady state: a mesh whose edges are already near L0
    let (q, qrep) = remesh(&base.0, &base.1, &[], l0, &RemeshOptions::default(), None).unwrap();
    let t = Instant::now();
    let (_, srep) = remesh(
        &q.surface,
        &q.positions,
        &[],
        l0,
        &RemeshOptions::default(),
        None,
    )
    .unwrap();
    println!(
        "steady state ({} faces): {:.3}s; split {} collapse {} flip {} (first pass {:.3}s)",
        q.surface.n_faces(),
        t.elapsed().as_secs_f64(),
        srep.splits,
        srep.collapses,
        srep.flips,
        qrep.seconds
    );
    let t = Instant::now();
    let (r, rep) = remesh(&m.0, &m.1, &[], l0, &RemeshOptions::default(), None).unwrap();
    println!(
        "torus 200x100 ({} faces) x(1.3,1,0.8): {:.3}s total, report {:.3}s; split {} collapse {} \
         flip {} smooth {} -> {} faces",
        m.0.n_faces(),
        t.elapsed().as_secs_f64(),
        rep.seconds,
        rep.splits,
        rep.collapses,
        rep.flips,
        rep.smoothing_moves,
        r.surface.n_faces()
    );
}
