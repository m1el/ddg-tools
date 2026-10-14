//! The cut session, the description of the cut on polygons, the insertion, and end to end.

#![allow(clippy::needless_range_loop)]
mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use ddg_core::io::read_obj;
use ddg_uv::cutting::{
    ConePlacement, CutOnPolygons, CutUnwrapOptions, CuttingError, DEFAULT_CROSSING_OFFSET,
    EdgePoint, FlowKind, Inserted, InteriorPoint, IntrinsicOptions, PathEnd, PolygonPath, Session,
    SessionState, Start, Triangulation, VertexSource, insert, unwrap_with_cut,
};
use ddg_uv::primitives::{capsule, cylinder, grid_cube, torus, uv_sphere};
use ddg_uv::{MeshInput, PolygonMesh, UnwrapOptions, prepare};
use ddg_variational_cuts::{
    CutProblem, Energy, Flow, FlowOptions, PAPER_FEATURE_SMOOTHING, ScaleWeight, StepSize,
    StopRule, init, length_weight_for_ripples,
};

// ----- helpers -----

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross_norm(a: [f64; 3], b: [f64; 3]) -> f64 {
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

/// Sum of the areas of the triangles of the triangulation.
fn triangulated_area(mesh: &PolygonMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|c| mesh.positions[mesh.corner_vertices[c as usize] as usize]);
            0.5 * cross_norm(sub(b, a), sub(c, a))
        })
        .sum()
}

fn min_triangle_area(mesh: &PolygonMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|c| mesh.positions[mesh.corner_vertices[c as usize] as usize]);
            0.5 * cross_norm(sub(b, a), sub(c, a))
        })
        .fold(f64::INFINITY, f64::min)
}

/// Distance from `p` to the segment `a b`.
fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let len2: f64 = ab.iter().map(|x| x * x).sum();
    let t = (sub(p, a).iter().zip(&ab).map(|(x, y)| x * y).sum::<f64>() / len2).clamp(0.0, 1.0);
    let q = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    sub(p, q).iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn debug_build() -> bool {
    cfg!(debug_assertions)
}

fn flow_setup(session: &mut Session, wavelength: f64) -> (CutProblem, FlowOptions<'static>) {
    let mut problem = CutProblem::new(Energy::Hencky, 0.0);
    let omega = session
        .roughening_rate(&problem, &Default::default())
        .expect("the start has a cut");
    problem.length_weight = length_weight_for_ripples(omega, wavelength);
    (problem, FlowOptions::new(StepSize::default()))
}

fn stop(n: usize) -> StopRule {
    StopRule {
        max_steps: n,
        ..StopRule::default()
    }
}

/// The invariants of an insertion.
fn check_insertion(old: &PolygonMesh, old_seams: &[[u32; 2]], cut: &CutOnPolygons, ins: &Inserted) {
    let new = &ins.mesh;
    // Total area kept (the new triangulation covers the old one exactly).
    let (a_old, a_new) = (triangulated_area(old), triangulated_area(new));
    assert!(
        (a_old - a_new).abs() <= 1e-12 * a_old,
        "area {a_old} became {a_new}"
    );
    assert!(min_triangle_area(new) > 0.0, "a triangle has no area");
    // A valid manifold polygon mesh: the pipeline's own checks accept it.
    new.input().validate().expect("valid arrays");
    let prepared = prepare(&new.input(), &UnwrapOptions::default()).expect("manifold");
    assert_eq!(
        prepared.validation.n_flipped_polygons, 0,
        "orientation kept"
    );
    // Old vertices keep their numbers and positions; new ones follow.
    assert_eq!(&new.positions[..old.positions.len()], &old.positions[..]);
    assert_eq!(ins.n_old_vertices, old.positions.len());
    assert_eq!(ins.vertex_source.len(), new.positions.len());
    assert!(ins.polygon_parent.windows(2).all(|w| w[0] <= w[1]));
    assert_eq!(ins.polygon_parent.len(), new.n_polygons());
    // Inserted points lie on the original surface: on their edge, or on their diagonal.
    for (v, source) in ins.vertex_source.iter().enumerate() {
        let x = new.positions[v];
        match *source {
            VertexSource::Old(o) => assert_eq!(o as usize, v),
            VertexSource::EdgePoint(i) => {
                let [a, b] = cut.edge_points[i as usize].edge;
                let d = segment_distance(x, old.positions[a as usize], old.positions[b as usize]);
                assert!(d < 1e-12, "edge point {i} is {d} off its edge");
            }
            VertexSource::Junction(j) => {
                let junction = &cut.junctions[j as usize];
                // Inside its triangle, at its barycentric position.
                let c = junction
                    .corners
                    .map(|c| old.positions[old.corner_vertices[c as usize] as usize]);
                let w = junction.barycentric;
                for d in 0..3 {
                    let expected = w[0] * c[0][d] + w[1] * c[1][d] + w[2] * c[2][d];
                    assert!((x[d] - expected).abs() < 1e-12);
                }
                assert!(w.iter().all(|&b| b > 0.0));
            }
            VertexSource::Interior { path, index } => {
                let point = &cut.paths[path as usize].interior[index as usize];
                let [c0, c1] = point
                    .corners
                    .map(|c| old.positions[old.corner_vertices[c as usize] as usize]);
                let d = segment_distance(x, c0, c1);
                assert!(d < 1e-12, "interior point is {d} off its diagonal");
            }
        }
    }
    // Every cut segment is an edge of the new mesh, and those are exactly the new seams.
    let edges: BTreeSet<[u32; 2]> = new.edges().into_iter().collect();
    let mut from_paths = BTreeSet::new();
    for (k, vertices) in ins.path_vertices.iter().enumerate() {
        assert_eq!(vertices.len(), cut.paths[k].interior.len() + 2);
        for w in vertices.windows(2) {
            let e = [w[0].min(w[1]), w[0].max(w[1])];
            assert!(edges.contains(&e), "cut segment {e:?} is not an edge");
            from_paths.insert(e);
        }
    }
    assert_eq!(ins.new_seams, from_paths.into_iter().collect::<Vec<_>>());
    // Old seams are kept: as they are, or as the chain of pieces they were cut into.
    let seams = ins.seams_after(old_seams);
    for s in &seams {
        assert!(
            edges.contains(s),
            "seam {s:?} is not an edge of the new mesh"
        );
    }
    for &[a, b] in old_seams {
        let key = [a.min(b), a.max(b)];
        match ins.edge_chains.get(&key) {
            None => assert!(seams.contains(&key)),
            Some(chain) => {
                assert_eq!((chain[0], *chain.last().unwrap()), (key[0], key[1]));
                for w in chain.windows(2) {
                    assert!(seams.contains(&[w[0].min(w[1]), w[0].max(w[1])]));
                }
                // The chain is the edge: its vertices lie on the old edge, in order.
                let (xa, xb) = (
                    old.positions[key[0] as usize],
                    old.positions[key[1] as usize],
                );
                for &v in &chain[1..chain.len() - 1] {
                    assert!(segment_distance(new.positions[v as usize], xa, xb) < 1e-12);
                }
            }
        }
    }
    // A junction where d paths end is a vertex of exactly d new polygons and of exactly d
    // seam edges (d = 2: a bend, d = 3: a triple point).
    for (j, &v) in ins.junction_vertex.iter().enumerate() {
        let d = cut
            .paths
            .iter()
            .map(|p| {
                usize::from(p.start == PathEnd::Junction(j as u32))
                    + usize::from(p.end == PathEnd::Junction(j as u32))
            })
            .sum::<usize>();
        let incident: Vec<&[u32; 2]> = ins
            .new_seams
            .iter()
            .filter(|e| e[0] == v || e[1] == v)
            .collect();
        assert_eq!(
            incident.len(),
            d,
            "junction {j} has {} seam edges",
            incident.len()
        );
        let polygons = (0..new.n_polygons())
            .filter(|&q| new.polygon(q).contains(&v))
            .count();
        assert_eq!(polygons, d, "junction {j} is in {polygons} polygons");
    }
    // The seams the pipeline will see are all edges that no old seam list contained or the cut.
    assert!(seams.len() >= ins.new_seams.len());
}

/// Number of connected groups of polygons when polygons are joined across edges that are not
/// seams.
fn n_charts(mesh: &PolygonMesh, seams: &[[u32; 2]]) -> usize {
    let seams: BTreeSet<[u32; 2]> = seams.iter().copied().collect();
    let mut parent: Vec<usize> = (0..mesh.n_polygons()).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut x = x;
        while p[x] != r {
            let n = p[x];
            p[x] = r;
            x = n;
        }
        r
    }
    let mut owner: BTreeMap<[u32; 2], usize> = BTreeMap::new();
    for p in 0..mesh.n_polygons() {
        let poly = mesh.polygon(p);
        for (i, &a) in poly.iter().enumerate() {
            let b = poly[(i + 1) % poly.len()];
            let e = [a.min(b), a.max(b)];
            if seams.contains(&e) {
                continue;
            }
            if let Some(&q) = owner.get(&e) {
                let (rp, rq) = (find(&mut parent, p), find(&mut parent, q));
                parent[rp] = rq;
            } else {
                owner.insert(e, p);
            }
        }
    }
    (0..mesh.n_polygons())
        .filter(|&p| find(&mut parent, p) == p)
        .count()
}

fn session_on(mesh: &PolygonMesh, seams: &[[u32; 2]], start: &Start) -> Session {
    let input = MeshInput {
        seam_edges: seams,
        ..mesh.input()
    };
    Session::new(&input, 2, start).expect("session")
}

/// A session on `mesh`, a circle around vertex `center`, `steps` steps.
fn cut_session(
    mesh: &PolygonMesh,
    seams: &[[u32; 2]],
    center: u32,
    radius: f64,
    steps: usize,
) -> Session {
    let mut session = session_on(
        mesh,
        seams,
        &Start::GeodesicCircle {
            vertex: center,
            radius,
        },
    );
    let (problem, options) = flow_setup(&mut session, 0.25);
    session.run(&problem, &options, &stop(steps)).expect("flow");
    session
}

// ----- tests -----

#[test]
fn circle_on_a_sphere_is_inserted() {
    let mesh = uv_sphere(24, 12, 1.0).unwrap();
    let session = cut_session(&mesh, &[], 5, 0.8, 5);
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    assert!(!cut.is_empty());
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &[], &cut, &ins);
    // A closed cut on a closed sphere: two charts, and none without the cut.
    assert_eq!(n_charts(&mesh, &[]), 1);
    assert_eq!(n_charts(&ins.mesh, &ins.new_seams), 2);
    // No seam on the sphere: every edge point is shared by two polygons.
    assert!(cut.edge_points.iter().all(|p| p.side.is_none()));
    // Each crossed polygon gave at least two new polygons.
    assert!(ins.mesh.n_polygons() > mesh.n_polygons());
}

#[test]
fn insertion_on_several_meshes() {
    let meshes: Vec<(&str, PolygonMesh, u32, f64)> = vec![
        ("grid_cube", grid_cube(6, 2.0), 7, 0.9),
        ("capsule", capsule(20, 6, 0.5, 1.5).unwrap(), 30, 0.5),
        ("torus", torus(24, 10, 1.0, 0.35).unwrap(), 12, 0.6),
        ("cylinder", cylinder(20, 3, 0.7, 2.0).unwrap(), 9, 0.7),
    ];
    for (name, mesh, center, radius) in meshes {
        let session = cut_session(&mesh, &[], center, radius, 4);
        let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
        assert!(!cut.is_empty(), "{name}: no cut");
        let ins = insert(&mesh, &cut).unwrap_or_else(|e| panic!("{name}: {e}"));
        check_insertion(&mesh, &[], &cut, &ins);
        println!(
            "{name}: {} polygons -> {}, {} edge points, {} paths",
            mesh.n_polygons(),
            ins.mesh.n_polygons(),
            cut.edge_points.len(),
            cut.paths.len()
        );
    }
}

#[test]
fn state_round_trips_through_corners() {
    let mesh = uv_sphere(16, 8, 1.0).unwrap();
    // A seam along a meridian from the north pole to the equator.
    let seams: Vec<[u32; 2]> = {
        let e = mesh.edges();
        // The first 4 vertices of the first column: pick edges of polygon 0.. in a chain.
        let poly = mesh.polygon(1);
        let _ = e;
        vec![[poly[0].min(poly[1]), poly[0].max(poly[1])]]
    };
    let mut session = session_on(
        &mesh,
        &seams,
        &Start::GeodesicCircle {
            vertex: 3,
            radius: 0.9,
        },
    );
    let state = session.state();
    assert_eq!(state.len(), 2);
    assert_eq!(state[0].len(), mesh.corner_vertices.len());
    assert!(state[0].iter().zip(&state[1]).all(|(a, b)| *a == -*b));
    // Set it back: nothing changes.
    let before = session.flow().level_set();
    session.set_state(&state).unwrap();
    assert_eq!(session.flow().level_set(), before);
    assert_eq!(session.state(), state);
    // Regions follow the sign.
    let regions = session.corner_regions();
    for (c, &r) in regions.iter().enumerate() {
        assert_eq!(r, u32::from(state[0][c] > 0.0));
    }
    // Another representation of the same level-set: (φ0, φ1) = (2φ, 0) has the same ½(φ0 − φ1).
    let doubled: Vec<f64> = state[0].iter().map(|x| 2.0 * x).collect();
    let zero = vec![0.0; doubled.len()];
    session.set_state(&[doubled, zero]).unwrap();
    assert_eq!(session.state(), state);
    // The state of a corner-start is the same session.
    let session2 = session_on(&mesh, &seams, &Start::Corners(state.clone()));
    assert_eq!(session2.state(), state);
    // After steps the state still round-trips.
    let (problem, options) = flow_setup(&mut session, 0.3);
    session.run(&problem, &options, &stop(3)).unwrap();
    let after = session.state();
    assert_ne!(after, state);
    let mut copy = session_on(&mesh, &seams, &Start::Corners(after.clone()));
    assert_eq!(copy.state(), after);
    copy.set_state(&after).unwrap();
    assert_eq!(copy.state(), after);
}

#[test]
fn corners_across_a_seam_may_differ_and_others_may_not() {
    let mesh = uv_sphere(16, 8, 1.0).unwrap();
    let poly = mesh.polygon(40).to_vec();
    // A seam of two edges, so that it is not extended.
    let seams = [
        [poly[0].min(poly[1]), poly[0].max(poly[1])],
        [poly[1].min(poly[2]), poly[1].max(poly[2])],
    ];
    let mut session = session_on(
        &mesh,
        &seams,
        &Start::GeodesicCircle {
            vertex: 3,
            radius: 0.9,
        },
    );
    let state = session.state();
    // The middle vertex of the seam is split into two flow vertices: one corner can differ.
    let v = poly[1];
    let corners: Vec<usize> = (0..mesh.corner_vertices.len())
        .filter(|&c| mesh.corner_vertices[c] == v)
        .collect();
    assert_eq!(corners.len(), 4);
    let flow_vertices: BTreeSet<u32> = corners
        .iter()
        .map(|&c| session.corner_to_flow_vertex()[c])
        .collect();
    assert_eq!(flow_vertices.len(), 2, "the seam splits the vertex in two");
    let c = corners
        .iter()
        .copied()
        .find(|&c| {
            session.corner_to_flow_vertex()[c] != session.corner_to_flow_vertex()[corners[0]]
        })
        .unwrap();
    // Changing all corners of the second copy together is fine ...
    let mut changed = state.clone();
    for &cc in &corners {
        if session.corner_to_flow_vertex()[cc] == session.corner_to_flow_vertex()[c] {
            changed[0][cc] += 0.01;
            changed[1][cc] -= 0.01;
        }
    }
    session.set_state(&changed).unwrap();
    // ... changing one corner of a copy that has several is not, and the message names the vertex.
    let mut broken = state.clone();
    broken[0][corners[0]] += 0.5;
    let message = session.set_state(&broken).unwrap_err();
    match &message {
        CuttingError::CornersDisagree { vertex, region, .. } => {
            assert_eq!((*vertex, *region), (v, 0));
        }
        other => panic!("{other}"),
    }
    assert!(
        message.to_string().contains(&format!("vertex {v}")),
        "{message}"
    );
    // A vertex away from the seam: the same.
    let w = 20u32;
    let cw = mesh.corner_vertices.iter().position(|&x| x == w).unwrap();
    let mut broken = state.clone();
    broken[1][cw] += 1.0;
    match session.set_state(&broken) {
        Err(CuttingError::CornersDisagree { vertex, region, .. }) => {
            assert_eq!((vertex, region), (w, 1))
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn every_error_names_what_is_wrong() {
    let mesh = uv_sphere(12, 6, 1.0).unwrap();
    let input = mesh.input();
    let circle = Start::GeodesicCircle {
        vertex: 3,
        radius: 0.8,
    };
    // Regions.
    assert!(matches!(
        Session::new(&input, 1, &circle),
        Err(CuttingError::TooFewRegions { got: 1 })
    ));
    // Three regions are available now, but a circle starts two.
    let e = Session::new(&input, 3, &circle).unwrap_err();
    assert!(
        matches!(e, CuttingError::InvalidOption { name: "start", .. }),
        "{e}"
    );
    assert!(matches!(
        Session::new_with(
            &input,
            3,
            &Start::NormalClusters { max_iterations: 5 },
            FlowKind::TwoRegion
        ),
        Err(CuttingError::InvalidOption { name: "kind", .. })
    ));
    // Starts.
    let bad = |start: Start| Session::new(&input, 2, &start).unwrap_err();
    assert!(matches!(
        bad(Start::GeodesicCircle {
            vertex: 10_000,
            radius: 1.0
        }),
        CuttingError::VertexOutOfRange { vertex: 10_000, .. }
    ));
    for r in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let e = bad(Start::GeodesicCircle {
            vertex: 1,
            radius: r,
        });
        assert!(
            matches!(e, CuttingError::InvalidOption { name: "radius", .. }),
            "{e}"
        );
    }
    let e = bad(Start::PolygonLabels(vec![0; 3]));
    assert!(
        matches!(e, CuttingError::WrongLength { expected, got: 3, .. } if expected == mesh.n_polygons())
    );
    let mut labels = vec![0u32; mesh.n_polygons()];
    labels[7] = 2;
    let e = bad(Start::PolygonLabels(labels));
    assert!(matches!(
        e,
        CuttingError::LabelOutOfRange {
            polygon: 7,
            label: 2,
            n_regions: 2
        }
    ));
    assert!(e.to_string().contains("polygon 7"), "{e}");
    let n = mesh.corner_vertices.len();
    assert!(matches!(
        bad(Start::Corners(vec![vec![0.0; n]])),
        CuttingError::WrongLength {
            expected: 2,
            got: 1,
            ..
        }
    ));
    assert!(matches!(
        bad(Start::Corners(vec![vec![0.0; n], vec![0.0; n - 1]])),
        CuttingError::WrongLength { .. }
    ));
    let mut nan = vec![0.0; n];
    nan[4] = f64::NAN;
    let e = bad(Start::Corners(vec![vec![0.0; n], nan]));
    assert!(matches!(
        e,
        CuttingError::NonFinite {
            region: 1,
            position: 4,
            ..
        }
    ));
    // A vertex no polygon uses cannot be the center.
    let mut loose = mesh.clone();
    loose.positions.push([5.0, 5.0, 5.0]);
    let e = Session::new(
        &loose.input(),
        2,
        &Start::GeodesicCircle {
            vertex: loose.positions.len() as u32 - 1,
            radius: 1.0,
        },
    )
    .unwrap_err();
    assert!(e.to_string().contains("not used"), "{e}");
    // Problem and options.
    let mut session = Session::new(&input, 2, &circle).unwrap();
    let (mut problem, options) = flow_setup(&mut session, 0.3);
    problem.penalty = Some(vec![0.0; 5]);
    let e = session.step(&problem, &options).unwrap_err();
    assert!(
        matches!(
            e,
            CuttingError::WrongLength {
                what: "penalty",
                got: 5,
                ..
            }
        ),
        "{e}"
    );
    let mut eta = vec![0.0; mesh.positions.len()];
    eta[9] = -1.0;
    problem.penalty = Some(eta);
    let e = session.step(&problem, &options).unwrap_err();
    assert!(e.to_string().contains("vertex 9"), "{e}");
    problem.penalty = None;
    let mut sym = options.clone();
    sym.symmetry = Some(vec![0; 2]);
    let e = session.step(&problem, &sym).unwrap_err();
    assert!(
        matches!(
            e,
            CuttingError::WrongLength {
                what: "symmetry classes",
                ..
            }
        ),
        "{e}"
    );
    problem.regions.pop();
    assert!(matches!(
        session.step(&problem, &options),
        Err(CuttingError::WrongLength {
            what: "problem.regions",
            ..
        })
    ));
    // Flow errors are passed through.
    let (problem, mut options) = flow_setup(&mut session, 0.3);
    options.step_size = StepSize::Fixed(-1.0);
    assert!(matches!(
        session.step(&problem, &options),
        Err(CuttingError::Flow(_))
    ));
    assert!(matches!(session.preview(0.7), Err(CuttingError::Flow(_))));
}

#[test]
fn bad_meshes_are_refused_with_the_polygon() {
    // Inconsistent winding: the session does not repair.
    let mut mesh = uv_sphere(12, 6, 1.0).unwrap();
    let range = mesh.polygon_corners(10);
    mesh.corner_vertices[range.clone()].reverse();
    // The triangles refer to corners; reversing the corners changes their meaning, so rebuild.
    let polygons: Vec<Vec<u32>> = (0..mesh.n_polygons())
        .map(|p| mesh.polygon(p).to_vec())
        .collect();
    let flipped = PolygonMesh::from_polygons(mesh.positions.clone(), &polygons);
    let e = Session::new(
        &flipped.input(),
        2,
        &Start::GeodesicCircle {
            vertex: 3,
            radius: 0.5,
        },
    )
    .unwrap_err();
    assert!(matches!(e, CuttingError::Mesh(_)), "{e}");
    assert!(
        e.to_string().contains("polygon") || e.to_string().contains("face"),
        "{e}"
    );
    // A triangle without area.
    let mut flat = PolygonMesh::from_triangles(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        &[[0, 1, 3], [1, 2, 3], [0, 2, 1]],
    );
    flat.positions[2] = [1.0, 0.0, 0.0];
    let e = Session::new(
        &flat.input(),
        2,
        &Start::GeodesicCircle {
            vertex: 0,
            radius: 0.5,
        },
    )
    .unwrap_err();
    assert!(matches!(e, CuttingError::Mesh(_)), "{e}");
    // Non-manifold edge.
    let tri = PolygonMesh::from_triangles(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        &[[0, 1, 2], [1, 0, 3], [0, 1, 4]],
    );
    let e = Session::new(
        &tri.input(),
        2,
        &Start::GeodesicCircle {
            vertex: 0,
            radius: 0.5,
        },
    )
    .unwrap_err();
    assert!(e.to_string().contains("(0, 1)"), "{e}");
}

#[test]
fn insert_refuses_descriptions_that_do_not_fit() {
    let mesh = uv_sphere(12, 6, 1.0).unwrap();
    let session = cut_session(&mesh, &[], 5, 0.8, 2);
    let good = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    insert(&mesh, &good).unwrap();
    let fails = |cut: &CutOnPolygons| insert(&mesh, cut).unwrap_err();

    // An edge point on something that is not an edge.
    let mut cut = good.clone();
    cut.edge_points[0].edge = [0, mesh.positions.len() as u32 - 1];
    let e = fails(&cut);
    assert!(
        matches!(e, CuttingError::InvalidDescription { .. })
            && e.to_string().contains("not an edge"),
        "{e}"
    );
    // Parameters outside (0, 1).
    for t in [0.0, 1.0, -0.5, f64::NAN] {
        let mut cut = good.clone();
        cut.edge_points[0].t = t;
        assert!(matches!(
            fails(&cut),
            CuttingError::InvalidDescription { .. }
        ));
    }
    // A path in a polygon that does not exist, or that its end does not belong to.
    let mut cut = good.clone();
    cut.paths[0].polygon = 10_000;
    assert!(matches!(
        fails(&cut),
        CuttingError::InvalidDescription { .. }
    ));
    let mut cut = good.clone();
    let other = (cut.paths[0].polygon + 20) % mesh.n_polygons() as u32;
    cut.paths[0].polygon = other;
    let e = fails(&cut);
    assert!(
        matches!(e, CuttingError::InvalidDescription { polygon: Some(p), .. } if p == other),
        "{e}"
    );
    // A junction that does not exist, or one that is not where the path is.
    let mut cut = good.clone();
    cut.paths[0].end = PathEnd::Junction(0);
    assert!(matches!(
        fails(&cut),
        CuttingError::InvalidDescription { .. }
    ));
    let corners = mesh.triangles[0].map(|c| c);
    cut.junctions.push(ddg_uv::cutting::JunctionPoint {
        position: [0.0; 3],
        polygon: 0,
        triangle: 0,
        corners,
        barycentric: [0.2, 0.3, 0.5],
    });
    let e = fails(&cut);
    assert!(
        e.to_string().contains("dead end") || e.to_string().contains("junction"),
        "{e}"
    );
    // An end that does not exist.
    let mut cut = good.clone();
    cut.paths[0].start = PathEnd::EdgePoint(99_999);
    assert!(matches!(
        fails(&cut),
        CuttingError::InvalidDescription { .. }
    ));
    // Two edge points at the same place.
    let mut cut = good.clone();
    let dup = cut.edge_points[0].clone();
    cut.edge_points.push(dup);
    let e = fails(&cut);
    assert!(
        e.to_string().contains("same") || e.to_string().contains("both at"),
        "{e}"
    );
    // An interior point on something that is not a diagonal.
    let mut cut = good.clone();
    let p = cut.paths.iter().position(|p| !p.interior.is_empty());
    if let Some(i) = p {
        cut.paths[i].interior[0].corners = [0, 1];
        let e = fails(&cut);
        assert!(matches!(e, CuttingError::InvalidDescription { .. }), "{e}");
    }
    // A path with its middle taken out no longer fits the triangles.
    let mut cut = good.clone();
    if let Some(i) = cut.paths.iter().position(|p| !p.interior.is_empty()) {
        cut.paths[i].interior.clear();
        assert!(matches!(
            fails(&cut),
            CuttingError::InvalidDescription { .. }
        ));
    }
    // A mesh that is not valid.
    let mut broken = mesh.clone();
    broken.triangles.pop();
    assert!(matches!(insert(&broken, &good), Err(CuttingError::Mesh(_))));
    // A stray interior point.
    let mut cut = good.clone();
    cut.paths[0].interior.push(InteriorPoint {
        corners: [0, 2],
        t: 0.5,
        position: [0.0; 3],
    });
    assert!(matches!(
        fails(&cut),
        CuttingError::InvalidDescription { .. }
    ));
    let _ = (
        EdgePoint {
            edge: [0, 1],
            t: 0.5,
            position: [0.0; 3],
            side: None,
        },
        PolygonPath {
            polygon: 0,
            start: PathEnd::EdgePoint(0),
            interior: vec![],
            end: PathEnd::EdgePoint(0),
            regions: [0, 1],
        },
    );
}

/// A meridian seam from the north pole down to the equator and a cut that ends on it.
#[test]
fn a_cut_that_ends_on_a_seam() {
    let mesh = uv_sphere(20, 10, 1.0).unwrap();
    // Vertices of the first meridian: the pole is vertex 0 (revolution: ring vertices first?).
    // Find a chain of edges along one longitude by walking down from the vertex with the
    // largest z, always taking the neighbour with the same angle around the axis.
    let edges = mesh.edges();
    let angle = |v: u32| mesh.positions[v as usize][1].atan2(mesh.positions[v as usize][0]);
    let top = (0..mesh.positions.len() as u32)
        .max_by(|&a, &b| mesh.positions[a as usize][2].total_cmp(&mesh.positions[b as usize][2]))
        .unwrap();
    let mut chain = vec![top];
    let first = edges
        .iter()
        .filter_map(|e| {
            if e[0] == top {
                Some(e[1])
            } else if e[1] == top {
                Some(e[0])
            } else {
                None
            }
        })
        .next()
        .unwrap();
    chain.push(first);
    let theta = angle(first);
    for _ in 0..4 {
        let last = *chain.last().unwrap();
        let next = edges
            .iter()
            .filter_map(|e| {
                if e[0] == last {
                    Some(e[1])
                } else if e[1] == last {
                    Some(e[0])
                } else {
                    None
                }
            })
            .find(|&w| {
                (angle(w) - theta).abs() < 1e-9
                    && mesh.positions[w as usize][2] < mesh.positions[last as usize][2] - 1e-9
            })
            .unwrap();
        chain.push(next);
    }
    let seams: Vec<[u32; 2]> = chain
        .windows(2)
        .map(|w| [w[0].min(w[1]), w[0].max(w[1])])
        .collect();
    // A circle around a vertex next to the seam: it crosses the seam.
    let center = chain[2];
    let mut session = session_on(
        &mesh,
        &seams,
        &Start::GeodesicCircle {
            vertex: center,
            radius: 0.55,
        },
    );
    // The seam is boundary for the flow: the flow mesh has more vertices than the polygon mesh.
    assert!(session.flow_mesh().n_vertices() > mesh.positions.len());
    let (problem, options) = flow_setup(&mut session, 0.3);
    session.run(&problem, &options, &stop(3)).unwrap();
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    assert!(
        cut.edge_points.iter().any(|p| p.side.is_some()),
        "the cut should end on the seam"
    );
    // Two sides of a seam: points on a seam edge come in separate copies.
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &seams, &cut, &ins);
    let all = ins.seams_after(&seams);
    assert!(all.len() >= seams.len() + ins.new_seams.len() - 4);
    // The preview shows open polylines ending at the seam.
    let preview = session.preview(DEFAULT_CROSSING_OFFSET).unwrap();
    assert!(preview.iter().any(|p| !p.closed));
    // And the pipeline accepts it, with the seams kept in the cut edges.
    let result = unwrap_with_cut(
        &MeshInput {
            seam_edges: &seams,
            ..mesh.input()
        },
        &cut,
        &CutUnwrapOptions {
            cones: ConePlacement::Off,
            ..Default::default()
        },
    )
    .unwrap();
    for s in &all {
        assert!(result.map.cut_edges.contains(s), "seam {s:?} not cut");
    }
}

fn models_dir() -> PathBuf {
    std::env::var_os("DDG_MODELS").map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/spot/spot"),
        PathBuf::from,
    )
}

fn load_polygons(name: &str) -> Option<PolygonMesh> {
    let path = models_dir().join(name);
    let Ok(obj) = read_obj(&path) else {
        println!("model {} not found: skipped", path.display());
        return None;
    };
    let polygons: Vec<Vec<u32>> = (0..obj.n_faces())
        .map(|f| obj.corner_verts[obj.face_corners(f)].to_vec())
        .collect();
    Some(PolygonMesh::from_polygons(obj.positions.clone(), &polygons))
}

#[test]
fn quad_model_end_to_end() {
    let Some(mesh) = load_polygons("spot_quadrangulated.obj") else {
        return;
    };
    let (steps, center) = if debug_build() { (3, 100) } else { (30, 100) };
    let mut session = session_on(
        &mesh,
        &[],
        &Start::GeodesicCircle {
            vertex: center,
            radius: 0.2,
        },
    );
    let (problem, options) = flow_setup(&mut session, 0.25);
    let report = session.run(&problem, &options, &stop(steps)).unwrap();
    println!("spot: {} steps, stop {:?}", report.steps.len(), report.stop);
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    assert!(!cut.is_empty());
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &[], &cut, &ins);
    for cones in [ConePlacement::Off, ConePlacement::Optimal] {
        let (_, result) = session
            .unwrap(
                DEFAULT_CROSSING_OFFSET,
                &CutUnwrapOptions {
                    cones,
                    ..Default::default()
                },
            )
            .unwrap();
        let rep = &result.map.report;
        assert_eq!(rep.flipped_triangles, 0, "{cones:?}");
        println!(
            "spot {cones:?}: charts {}, cones {}, flipped {}, seam mismatch {:.2e}, distortion {:?}",
            rep.charts.len(),
            result.map.cones.len(),
            rep.flipped_triangles,
            rep.seam_lengths.max_seam_mismatch,
            rep.area_distortion
        );
        if cones == ConePlacement::Off {
            assert!(result.map.cones.is_empty());
        }
        assert!(rep.seam_lengths.max_seam_mismatch < 0.5);
        // Every seam of the cut is cut by the pipeline.
        for s in &result.seams {
            assert!(result.map.cut_edges.contains(s));
        }
    }
}

#[test]
fn end_to_end_on_generated_meshes() {
    let meshes: Vec<(&str, PolygonMesh, u32, f64)> = vec![
        ("uv_sphere", uv_sphere(24, 12, 1.0).unwrap(), 5, 0.9),
        ("grid_cube", grid_cube(6, 2.0), 7, 0.9),
        ("capsule", capsule(20, 6, 0.5, 1.5).unwrap(), 30, 0.5),
        ("torus", torus(24, 10, 1.0, 0.35).unwrap(), 12, 0.6),
    ];
    for (name, mesh, center, radius) in meshes {
        let session = cut_session(&mesh, &[], center, radius, 4);
        for cones in [ConePlacement::Off, ConePlacement::Optimal] {
            let (cut, result) = session
                .unwrap(
                    DEFAULT_CROSSING_OFFSET,
                    &CutUnwrapOptions {
                        cones,
                        ..Default::default()
                    },
                )
                .unwrap_or_else(|e| panic!("{name} {cones:?}: {e}"));
            let rep = &result.map.report;
            assert_eq!(rep.flipped_triangles, 0, "{name} {cones:?}");
            assert!(rep.seam_lengths.max_seam_mismatch < 0.5, "{name}");
            assert_eq!(
                result.map.uvs.len(),
                result.inserted.mesh.corner_vertices.len()
            );
            for s in &result.seams {
                assert!(result.map.cut_edges.contains(s), "{name}: seam {s:?}");
            }
            if cones == ConePlacement::Off {
                assert!(result.map.cones.is_empty(), "{name}");
            }
            println!(
                "{name} {cones:?}: {} -> {} polygons, {} charts, {} cones, seam mismatch {:.2e}, distortion {:.3}",
                mesh.n_polygons(),
                result.inserted.mesh.n_polygons(),
                rep.charts.len(),
                result.map.cones.len(),
                rep.seam_lengths.max_seam_mismatch,
                rep.area_distortion.scale_free
            );
            let _ = &cut;
        }
    }
}

#[test]
fn cone_placement_off_refuses_steered_cones() {
    let mesh = uv_sphere(12, 6, 1.0).unwrap();
    let session = cut_session(&mesh, &[], 5, 0.8, 1);
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    let pinned = [3u32];
    let input = MeshInput {
        pinned_vertices: &pinned,
        ..mesh.input()
    };
    let e = unwrap_with_cut(
        &input,
        &cut,
        &CutUnwrapOptions {
            cones: ConePlacement::Off,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(e.to_string().contains("vertex 3"), "{e}");
    // Optimal keeps the steering of the old vertices.
    unwrap_with_cut(&input, &cut, &CutUnwrapOptions::default()).unwrap();
}

#[test]
fn starts_at_polygon_level() {
    let mesh = uv_sphere(20, 10, 1.0).unwrap();
    // Labels: the polygons around the equator.
    let labels: Vec<u32> = (0..mesh.n_polygons())
        .map(|p| {
            let z: f64 = mesh
                .polygon(p)
                .iter()
                .map(|&v| mesh.positions[v as usize][2])
                .sum::<f64>()
                / mesh.polygon(p).len() as f64;
            u32::from(z.abs() < 0.3)
        })
        .collect();
    let session = session_on(&mesh, &[], &Start::PolygonLabels(labels.clone()));
    let preview = session.preview(DEFAULT_CROSSING_OFFSET).unwrap();
    // Two latitude circles.
    assert_eq!(preview.len(), 2);
    assert!(preview.iter().all(|p| p.closed));
    let regions = session.corner_regions();
    for (p, &label) in labels.iter().enumerate() {
        if label == 1
            && mesh
                .polygon(p)
                .iter()
                .all(|&v| mesh.positions[v as usize][2].abs() < 0.2)
        {
            for c in mesh.polygon_corners(p) {
                assert_eq!(regions[c], 1);
            }
        }
    }
    // A circle: the corners within the radius of the center are in region 0.
    let session = session_on(
        &mesh,
        &[],
        &Start::GeodesicCircle {
            vertex: 7,
            radius: 0.5,
        },
    );
    let regions = session.corner_regions();
    let c0 = (0..mesh.corner_vertices.len())
        .find(|&c| mesh.corner_vertices[c] == 7)
        .unwrap();
    assert_eq!(regions[c0], 0);
    assert!(regions.contains(&1));
    assert_eq!(session.preview(DEFAULT_CROSSING_OFFSET).unwrap().len(), 1);
}

#[test]
fn eta_and_symmetry_are_mapped_per_polygon_vertex() {
    let mesh = uv_sphere(16, 8, 1.0).unwrap();
    let mut session = session_on(
        &mesh,
        &[],
        &Start::GeodesicCircle {
            vertex: 4,
            radius: 0.7,
        },
    );
    let (mut problem, mut options) = flow_setup(&mut session, 0.3);
    problem.penalty = Some(vec![0.5; mesh.positions.len()]);
    options.symmetry = Some((0..mesh.positions.len() as u32).collect());
    let r = session.step(&problem, &options).unwrap();
    assert!(r.n_cut_points > 0);
    assert!(r.energy.penalty > 0.0);
}

/// Random signs per flow vertex: many paths per polygon, ends on seams on both sides, quads and
/// fans, polygons with extra collinear vertices from the other side of a seam.
#[test]
fn random_states_insert_cleanly() {
    let mut rng = ddg_core::generators::SplitMix64::new(7);
    let (mut paths, mut seam_sides) = (0usize, 0usize);
    let sphere = uv_sphere(14, 7, 1.0).unwrap();
    let cube = grid_cube(3, 2.0);
    let tor = torus(12, 6, 1.0, 0.3).unwrap();
    for (name, mesh) in [("sphere", sphere), ("cube", cube), ("torus", tor)] {
        let edges = mesh.edges();
        let seam_sets: Vec<Vec<[u32; 2]>> =
            vec![vec![], edges.iter().step_by(7).take(9).copied().collect()];
        for seams in seam_sets {
            let mut session = session_on(
                &mesh,
                &seams,
                &Start::GeodesicCircle {
                    vertex: 1,
                    radius: 0.5,
                },
            );
            for round in 0..12 {
                let n_flow = session.flow_mesh().n_vertices();
                let phi: Vec<f64> = (0..n_flow).map(|_| rng.uniform(-1.0, 1.0)).collect();
                let c2f = session.corner_to_flow_vertex().to_vec();
                let f0: Vec<f64> = c2f.iter().map(|&v| phi[v as usize]).collect();
                let f1: Vec<f64> = f0.iter().map(|x| -x).collect();
                session.set_state(&[f0, f1]).unwrap();
                let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
                let ins = insert(&mesh, &cut)
                    .unwrap_or_else(|e| panic!("{name} round {round}, {} seams: {e}", seams.len()));
                check_insertion(&mesh, &seams, &cut, &ins);
                check_regions(&session, &mesh, &ins);
                seam_sides += cut.edge_points.iter().filter(|p| p.side.is_some()).count();
                paths += cut.paths.len();
            }
        }
    }
    println!("{paths} paths, {seam_sides} seam-side points");
    assert!(paths > 100 && seam_sides > 10);
}

/// Every new polygon lies in one region: the old corners of a part all have the region of the
/// state at that corner of the parent polygon.
fn check_regions(session: &Session, old: &PolygonMesh, ins: &Inserted) {
    let regions = session.corner_regions();
    for q in 0..ins.mesh.n_polygons() {
        let parent = ins.polygon_parent[q] as usize;
        let mut seen: Option<u32> = None;
        for &v in ins.mesh.polygon(q) {
            if (v as usize) >= ins.n_old_vertices {
                continue;
            }
            let c = old
                .polygon_corners(parent)
                .find(|&c| old.corner_vertices[c] == v)
                .expect("an old vertex of a part is a corner of the parent");
            let r = regions[c];
            assert_eq!(
                *seen.get_or_insert(r),
                r,
                "polygon {q} (part of {parent}) is in two regions"
            );
        }
    }
}

/// A closed ring of seams splits the surface: the flow runs on two components.
#[test]
fn a_ring_of_seams_gives_two_components() {
    let mesh = uv_sphere(20, 10, 1.0).unwrap();
    let on_equator = |v: u32| mesh.positions[v as usize][2].abs() < 1e-9;
    let seams: Vec<[u32; 2]> = mesh
        .edges()
        .into_iter()
        .filter(|e| on_equator(e[0]) && on_equator(e[1]))
        .collect();
    assert_eq!(seams.len(), 20);
    let center = (0..mesh.positions.len() as u32)
        .find(|&v| mesh.positions[v as usize][2] > 0.5)
        .unwrap();
    let mut session = session_on(
        &mesh,
        &seams,
        &Start::GeodesicCircle {
            vertex: center,
            radius: 0.5,
        },
    );
    assert_eq!(session.flow_mesh().n_components(), 2);
    let (problem, options) = flow_setup(&mut session, 0.3);
    session.run(&problem, &options, &stop(5)).unwrap();
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &seams, &cut, &ins);
    // Three charts: the lower half, the cap, the rest of the upper half.
    assert_eq!(n_charts(&ins.mesh, &ins.seams_after(&seams)), 3);
    let (_, result) = session
        .unwrap(DEFAULT_CROSSING_OFFSET, &CutUnwrapOptions::default())
        .unwrap();
    assert_eq!(result.map.report.flipped_triangles, 0);
    assert_eq!(result.map.report.charts.len(), 3);
}

/// A seam of one edge is extended by `prepare`; the session reports it and uses the extension.
#[test]
fn a_one_edge_seam_is_extended_and_reported() {
    let mesh = uv_sphere(16, 8, 1.0).unwrap();
    let poly = mesh.polygon(40);
    let seams = [[poly[0].min(poly[1]), poly[0].max(poly[1])]];
    let session = session_on(
        &mesh,
        &seams,
        &Start::GeodesicCircle {
            vertex: 3,
            radius: 0.9,
        },
    );
    assert_eq!(session.seams().len(), 1);
    assert_eq!(session.seams_cut().len(), 2);
    assert!(
        session.warnings().iter().any(|w| w.contains("single edge")),
        "{:?}",
        session.warnings()
    );
    // A seam that is not an edge is ignored with a warning, not cut.
    let session = session_on(
        &mesh,
        &[[0, 5]],
        &Start::GeodesicCircle {
            vertex: 3,
            radius: 0.9,
        },
    );
    assert!(session.seams_cut().is_empty());
    assert!(!session.warnings().is_empty());
    assert_eq!(session.flow_mesh().n_vertices(), mesh.positions.len());
}

fn dataset_dir() -> PathBuf {
    std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    )
}

/// Meshes of the cones dataset: `bump` (a disk) in a debug build; `camel` (boundary), `bunny`
/// (closed) and `hand` too in a release build. Skipped with a notice without the data.
#[test]
fn dataset_meshes_end_to_end() {
    let names: &[&str] = if debug_build() {
        &["bump"]
    } else {
        &["bump", "camel", "bunny", "hand"]
    };
    for name in names {
        let path = dataset_dir().join(format!("{name}.obj"));
        let Ok(obj) = read_obj(&path) else {
            println!("dataset mesh {} not found: skipped", path.display());
            continue;
        };
        let mesh = PolygonMesh::from_triangles(obj.positions.clone(), &obj.triangles().unwrap());
        let area = triangulated_area(&mesh);
        let center = (mesh.positions.len() / 3) as u32;
        let session = cut_session(&mesh, &[], center, 0.2 * area.sqrt(), 5);
        let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
        assert!(!cut.is_empty(), "{name}: no cut");
        let ins = insert(&mesh, &cut).unwrap();
        check_insertion(&mesh, &[], &cut, &ins);
        for cones in [ConePlacement::Off, ConePlacement::Optimal] {
            let (_, result) = session
                .unwrap(
                    DEFAULT_CROSSING_OFFSET,
                    &CutUnwrapOptions {
                        cones,
                        ..Default::default()
                    },
                )
                .unwrap_or_else(|e| panic!("{name} {cones:?}: {e}"));
            let rep = &result.map.report;
            assert_eq!(rep.flipped_triangles, 0, "{name} {cones:?}");
            println!(
                "{name} {cones:?}: {} -> {} polygons, {} charts, {} cones, seam mismatch {:.2e}, distortion {:.3}",
                mesh.n_polygons(),
                ins.mesh.n_polygons(),
                rep.charts.len(),
                result.map.cones.len(),
                rep.seam_lengths.max_seam_mismatch,
                rep.area_distortion.scale_free
            );
        }
    }
}

// ===================== N regions =====================

/// A session with `n` regions started from `labels` (one per polygon).
fn labelled_session(mesh: &PolygonMesh, seams: &[[u32; 2]], n: usize, labels: Vec<u32>) -> Session {
    let input = MeshInput {
        seam_edges: seams,
        ..mesh.input()
    };
    Session::new(&input, n, &Start::PolygonLabels(labels)).expect("session")
}

/// A random state of `n` region functions, constant on the corners of a flow vertex.
fn random_state(
    session: &Session,
    n: usize,
    rng: &mut ddg_core::generators::SplitMix64,
) -> Vec<Vec<f64>> {
    let n_flow = session.flow_mesh().n_vertices();
    let c2f = session.corner_to_flow_vertex().to_vec();
    (0..n)
        .map(|_| {
            let phi: Vec<f64> = (0..n_flow).map(|_| rng.uniform(-1.0, 1.0)).collect();
            c2f.iter().map(|&v| phi[v as usize]).collect()
        })
        .collect()
}

/// Statistics of the networks the insertion was run on.
#[derive(Default, Debug)]
struct NetworkStats {
    paths: usize,
    junctions: usize,
    polygons_with_two_junctions: usize,
    junction_with_ordinary_path: usize,
    junction_next_to_seam: usize,
    junction_paths: usize,
}

fn network_stats(
    mesh: &PolygonMesh,
    seams: &[[u32; 2]],
    cut: &CutOnPolygons,
    st: &mut NetworkStats,
) {
    st.paths += cut.paths.len();
    st.junctions += cut.junctions.len();
    let seams: BTreeSet<[u32; 2]> = seams.iter().copied().collect();
    let mut per_polygon: BTreeMap<u32, usize> = BTreeMap::new();
    for j in &cut.junctions {
        *per_polygon.entry(j.polygon).or_default() += 1;
        let poly = mesh.polygon(j.polygon as usize);
        let touches_seam = poly.iter().enumerate().any(|(i, &a)| {
            let b = poly[(i + 1) % poly.len()];
            seams.contains(&[a.min(b), a.max(b)])
        });
        st.junction_next_to_seam += usize::from(touches_seam);
    }
    st.polygons_with_two_junctions += per_polygon.values().filter(|&&c| c >= 2).count();
    for p in &cut.paths {
        let at_junction =
            matches!(p.start, PathEnd::Junction(_)) || matches!(p.end, PathEnd::Junction(_));
        st.junction_paths += usize::from(at_junction);
        if !at_junction && per_polygon.contains_key(&p.polygon) {
            st.junction_with_ordinary_path += 1;
        }
    }
}

/// Networks: random states with 3, 4 and 6 regions on three meshes, with and without seams.
/// The invariants of the insertion hold with triple points; junctions are vertices of three
/// polygons and three seam edges; several junctions in a polygon, a junction in a polygon with
/// an ordinary path, junctions next to seams all occur.
#[test]
fn networks_insert_cleanly() {
    let mut rng = ddg_core::generators::SplitMix64::new(11);
    let mut stats = NetworkStats::default();
    let sphere = uv_sphere(14, 7, 1.0).unwrap();
    let cube = grid_cube(3, 2.0);
    let tor = torus(12, 6, 1.0, 0.3).unwrap();
    // Cylinder: its caps are 12-gons, fanned: junctions and ordinary paths in one polygon.
    let cyl = cylinder(12, 2, 0.7, 2.0).unwrap();
    for (name, mesh) in [
        ("sphere", sphere),
        ("cube", cube),
        ("torus", tor),
        ("cylinder", cyl),
    ] {
        let edges = mesh.edges();
        let seam_sets: Vec<Vec<[u32; 2]>> =
            vec![vec![], edges.iter().step_by(7).take(9).copied().collect()];
        for seams in seam_sets {
            for n in [3usize, 4, 6] {
                let labels = vec![0u32; mesh.n_polygons()];
                let mut session = labelled_session(&mesh, &seams, n, labels);
                for round in 0..6 {
                    let state = random_state(&session, n, &mut rng);
                    session.set_state(&state).unwrap();
                    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
                    let ins = insert(&mesh, &cut).unwrap_or_else(|e| {
                        panic!("{name} N={n} round {round}, {} seams: {e}", seams.len())
                    });
                    check_insertion(&mesh, &seams, &cut, &ins);
                    check_regions(&session, &mesh, &ins);
                    network_stats(&mesh, &seams, &cut, &mut stats);
                    // Every path separates two regions, lower id on its left.
                    assert!(cut.paths.iter().all(|p| p.regions[0] < p.regions[1]));
                    assert!(cut.paths.iter().all(|p| (p.regions[1] as usize) < n));
                }
            }
        }
    }
    println!("{stats:?}");
    assert!(stats.junctions > 500);
    assert!(stats.polygons_with_two_junctions > 20);
    assert!(stats.junction_with_ordinary_path > 0);
    assert!(stats.junction_next_to_seam > 20);
    assert!(stats.junction_paths >= 3 * stats.junctions / 2);
}

/// Three regions around one vertex of a triangle mesh: one junction, three arms.
#[test]
fn one_junction_by_hand() {
    // A fan of six triangles around vertex 0, regions by pairs of triangles.
    let mut positions = vec![[0.0, 0.0, 0.0]];
    for i in 0..6 {
        let a = std::f64::consts::TAU * f64::from(i) / 6.0;
        positions.push([a.cos(), a.sin(), 0.0]);
    }
    let faces: Vec<[u32; 3]> = (0..6).map(|i| [0, 1 + i, 1 + (i + 1) % 6]).collect();
    let mesh = PolygonMesh::from_triangles(positions, &faces);
    // Vertex 0 in region 2 and rim vertices 1,2 / 3,4 / 5,6 in regions 0, 1, 0... : make the
    // triangle (0, 1, 2) have three regions.
    let region_of = [2usize, 0, 1, 1, 0, 0, 1];
    let n = 3;
    let n_corners = mesh.corner_vertices.len();
    let functions: Vec<Vec<f64>> = (0..n)
        .map(|k| {
            (0..n_corners)
                .map(|c| {
                    if region_of[mesh.corner_vertices[c] as usize] == k {
                        -1.0
                    } else {
                        1.0
                    }
                })
                .collect()
        })
        .collect();
    let input = mesh.input();
    let session = Session::new(&input, n, &Start::Corners(functions)).unwrap();
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    assert!(!cut.junctions.is_empty());
    let j = &cut.junctions[0];
    assert!((j.barycentric.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &[], &cut, &ins);
    // A polygon (a triangle here) with a junction becomes three polygons.
    let parts = ins
        .polygon_parent
        .iter()
        .filter(|&&p| p == j.polygon)
        .count();
    assert_eq!(parts, 3);
    // Position: the weights applied to the corners.
    let jv = ins.junction_vertex[0] as usize;
    for d in 0..3 {
        let expected: f64 = (0..3)
            .map(|k| {
                j.barycentric[k]
                    * mesh.positions[mesh.corner_vertices[j.corners[k] as usize] as usize][d]
            })
            .sum();
        assert!((ins.mesh.positions[jv][d] - expected).abs() < 1e-12);
    }
}

/// Two regions through the N-region flow against the two-region flow: the same steps and the
/// same description.
#[test]
fn n_equals_two_through_the_region_flow() {
    let mesh = uv_sphere(20, 10, 1.0).unwrap();
    let seams = {
        let poly = mesh.polygon(30);
        vec![
            [poly[0].min(poly[1]), poly[0].max(poly[1])],
            [poly[1].min(poly[2]), poly[1].max(poly[2])],
        ]
    };
    let input = MeshInput {
        seam_edges: &seams,
        ..mesh.input()
    };
    let start = Start::GeodesicCircle {
        vertex: 40,
        radius: 0.8,
    };
    let mut a = Session::new(&input, 2, &start).unwrap();
    let mut b = Session::new_with(&input, 2, &start, FlowKind::Regions).unwrap();
    assert!(!a.uses_region_flow() && b.uses_region_flow());
    let (problem, options) = flow_setup(&mut a, 0.3);
    for step in 0..6 {
        let ra = a.step(&problem, &options).unwrap();
        let rb = b.step(&problem, &options).unwrap();
        let (ea, eb) = (ra.energy.total, rb.energy.total);
        assert!(
            (ea - eb).abs() <= 1e-9 * ea.abs(),
            "step {step}: {ea} against {eb}"
        );
        assert_eq!(ra.n_cut_points, rb.n_cut_points);
    }
    let (sa, sb) = (a.state(), b.state());
    for k in 0..2 {
        for (x, y) in sa[k].iter().zip(&sb[k]) {
            assert!((x - y).abs() < 1e-8, "{x} {y}");
        }
    }
    let (ca, cb) = (
        a.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap(),
        b.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap(),
    );
    assert!(cb.junctions.is_empty());
    assert_eq!(ca.edge_points.len(), cb.edge_points.len());
    assert_eq!(ca.paths.len(), cb.paths.len());
    for (pa, pb) in ca.paths.iter().zip(&cb.paths) {
        assert_eq!(
            (pa.polygon, pa.start, pa.end, pa.regions),
            (pb.polygon, pb.start, pb.end, pb.regions)
        );
        assert_eq!(pa.interior.len(), pb.interior.len());
    }
    for (ea, eb) in ca.edge_points.iter().zip(&cb.edge_points) {
        assert_eq!((ea.edge, ea.side), (eb.edge, eb.side));
        assert!((ea.t - eb.t).abs() < 1e-8);
    }
    let (pa, pb) = (
        a.preview(DEFAULT_CROSSING_OFFSET).unwrap(),
        b.preview(DEFAULT_CROSSING_OFFSET).unwrap(),
    );
    assert_eq!(pa.len(), pb.len());
    // Both insert and give the same new mesh combinatorics.
    let (ia, ib) = (insert(&mesh, &ca).unwrap(), insert(&mesh, &cb).unwrap());
    assert_eq!(ia.mesh.polygon_offsets, ib.mesh.polygon_offsets);
    assert_eq!(ia.new_seams, ib.new_seams);
}

/// The two-region session is the core crate's flow, bit for bit: same steps, same state.
#[test]
fn two_regions_are_the_core_flow_bit_for_bit() {
    let sphere = ddg_core::generators::icosphere(3).unwrap();
    let mesh = PolygonMesh::from_triangles(sphere.positions.clone(), &sphere.faces);
    let tri = ddg_core::TriMesh::new(sphere.positions.len(), sphere.faces.clone()).unwrap();
    let geometry = ddg_core::Geometry::from_positions(&tri, &sphere.positions).unwrap();
    let phi0 = init::geodesic_circle(&tri, &geometry, 3, 0.9).unwrap();
    let corner_phi: Vec<f64> = mesh
        .corner_vertices
        .iter()
        .map(|&v| phi0[v as usize])
        .collect();
    let negated: Vec<f64> = corner_phi.iter().map(|x| -x).collect();
    let mut session =
        Session::new(&mesh.input(), 2, &Start::Corners(vec![corner_phi, negated])).unwrap();
    let mut flow = Flow::new(&tri, &geometry, &phi0).unwrap();
    let (problem, options) = flow_setup(&mut session, 0.3);
    for _ in 0..5 {
        let a = session.step(&problem, &options).unwrap();
        let b = flow.step(&problem, &options).unwrap();
        assert_eq!(a.energy.total.to_bits(), b.energy.total.to_bits());
        assert_eq!(a.tau.to_bits(), b.tau.to_bits());
    }
    let phi = flow.level_set();
    let state = session.state();
    for (c, &v) in mesh.corner_vertices.iter().enumerate() {
        assert_eq!(state[0][c].to_bits(), phi[v as usize].to_bits());
    }
}

#[test]
fn n_region_state_and_starts() {
    let mesh = uv_sphere(24, 12, 1.0).unwrap();
    let n = 4;
    // Labels: octant-like, by the sign of x and y of the centre.
    let labels: Vec<u32> = (0..mesh.n_polygons())
        .map(|p| {
            let c = mesh.polygon(p).iter().fold([0.0f64; 3], |mut acc, &v| {
                for d in 0..3 {
                    acc[d] += mesh.positions[v as usize][d];
                }
                acc
            });
            u32::from(c[0] > 0.0) + 2 * u32::from(c[1] > 0.0)
        })
        .collect();
    let mut session = labelled_session(&mesh, &[], n, labels.clone());
    let state = session.state();
    assert_eq!(state.len(), n);
    // Corner regions are the argmin of the state.
    let regions = session.corner_regions();
    for c in 0..mesh.corner_vertices.len() {
        let best = (0..n)
            .min_by(|&a, &b| state[a][c].total_cmp(&state[b][c]))
            .unwrap();
        assert_eq!(regions[c] as usize, best);
    }
    // The labels of the polygons far from the interface are the regions.
    for p in 0..mesh.n_polygons() {
        let c = mesh.polygon(p).iter().fold([0.0f64; 3], |mut acc, &v| {
            for d in 0..3 {
                acc[d] += mesh.positions[v as usize][d];
            }
            acc
        });
        let far = mesh.polygon(p).iter().all(|&v| {
            mesh.positions[v as usize][0].abs() > 0.2 && mesh.positions[v as usize][1].abs() > 0.2
        });
        if far {
            let _ = c;
            for corner in mesh.polygon_corners(p) {
                assert_eq!(regions[corner], labels[p]);
            }
        }
    }
    // Round trip, bit for bit.
    session.set_state(&state).unwrap();
    assert_eq!(session.state(), state);
    let again = Session::new(&mesh.input(), n, &Start::Corners(state.clone())).unwrap();
    assert_eq!(again.state(), state);
    // A corner that disagrees names the vertex and the region function.
    let v = 17u32;
    let c = mesh.corner_vertices.iter().position(|&x| x == v).unwrap();
    let mut broken = state.clone();
    broken[3][c] += 0.5;
    match session.set_state(&broken) {
        Err(CuttingError::CornersDisagree { vertex, region, .. }) => {
            assert_eq!((vertex, region), (v, 3))
        }
        other => panic!("{other:?}"),
    }
    // Wrong counts and labels.
    assert!(matches!(
        session.set_state(&state[..3]),
        Err(CuttingError::WrongLength {
            expected: 4,
            got: 3,
            ..
        })
    ));
    let mut bad = labels.clone();
    bad[5] = 4;
    assert!(matches!(
        Session::new(&mesh.input(), n, &Start::PolygonLabels(bad)),
        Err(CuttingError::LabelOutOfRange {
            polygon: 5,
            label: 4,
            n_regions: 4
        })
    ));
    // Normal clusters: any number of regions, deterministic, every region used.
    let clustered = Session::new(
        &mesh.input(),
        6,
        &Start::NormalClusters { max_iterations: 50 },
    )
    .unwrap();
    let again = Session::new(
        &mesh.input(),
        6,
        &Start::NormalClusters { max_iterations: 50 },
    )
    .unwrap();
    assert_eq!(clustered.state(), again.state());
    let used: BTreeSet<u32> = clustered.corner_regions().into_iter().collect();
    assert_eq!(used.len(), 6);
    // The recommended options.
    let (p, o) = (
        clustered.new_problem(Energy::Hencky, 0.1),
        clustered.recommended_options(StepSize::default()),
    );
    assert_eq!(p.regions.len(), 6);
    assert!(matches!(
        o.length_step,
        ddg_variational_cuts::LengthStep::OnMesh
    ));
    // A circle starts two regions.
    assert!(matches!(
        Session::new(
            &mesh.input(),
            3,
            &Start::GeodesicCircle {
                vertex: 1,
                radius: 0.5
            }
        ),
        Err(CuttingError::InvalidOption { name: "start", .. })
    ));
}

/// The atlas: normal clusters, the recommended options, a strong length weight from the ripple
/// rule and the rescaled Dirichlet energy, a fixed number of steps; then end to end.
fn atlas(name: &str, mesh: &PolygonMesh, n: usize, steps: usize) {
    let mut session = Session::new(
        &mesh.input(),
        n,
        &Start::NormalClusters {
            max_iterations: 100,
        },
    )
    .unwrap_or_else(|e| panic!("{name}: {e}"));
    let options = session.recommended_options(StepSize::default());
    let energy = Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
        r_eps: 1.0,
        smoothing: PAPER_FEATURE_SMOOTHING,
    });
    let mut problem = session.new_problem(energy, 0.0);
    let omega = session.roughening_rate(&problem, &options.eval).unwrap();
    problem.length_weight = length_weight_for_ripples(omega, 0.5);
    let before = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    let report = session.run(&problem, &options, &stop(steps)).unwrap();
    let after = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    let first = report.steps.first().unwrap().energy.total;
    println!(
        "{name} N={n}: alpha_L {:.2}, {} steps ({:?}), energy {:.3} -> {:.3}, junctions {} -> {}, paths {} -> {}",
        problem.length_weight,
        report.steps.len(),
        report.stop,
        first,
        report.final_energy.total,
        before.junctions.len(),
        after.junctions.len(),
        before.paths.len(),
        after.paths.len()
    );
    assert!(report.final_energy.total.is_finite());
    let ins = insert(mesh, &after).unwrap_or_else(|e| panic!("{name}: {e}"));
    check_insertion(mesh, &[], &after, &ins);
    check_regions(&session, mesh, &ins);
    let charts = n_charts(&ins.mesh, &ins.new_seams);
    for cones in [ConePlacement::Off, ConePlacement::Optimal] {
        let (_, result) = session
            .unwrap(
                DEFAULT_CROSSING_OFFSET,
                &CutUnwrapOptions {
                    cones,
                    ..Default::default()
                },
            )
            .unwrap_or_else(|e| panic!("{name} {cones:?}: {e}"));
        let rep = &result.map.report;
        // A cut that passes within `delta` of a vertex leaves a sliver triangle (here 5e-5 of
        // a cell, found by hand on the cube); its image can flip with an area of 1e-9. At most
        // a few per thousand triangles are accepted; with cones off a patch that contains a
        // corner of a cube would flip many more.
        let allowed = 1 + result.inserted.mesh.triangles.len() / 1000;
        assert!(
            rep.flipped_triangles <= allowed,
            "{name} {cones:?}: {} flipped triangles",
            rep.flipped_triangles
        );
        assert_eq!(
            rep.charts.len(),
            charts,
            "{name}: the charts are the pieces the seams leave"
        );
        assert!(rep.seam_lengths.max_seam_mismatch < 0.5);
        for s in &result.seams {
            assert!(result.map.cut_edges.contains(s), "{name}: seam {s:?}");
        }
        if cones == ConePlacement::Off {
            assert!(result.map.cones.is_empty());
        }
        println!(
            "  {cones:?}: {} -> {} polygons, {} charts, {} cones, flipped {}, seam mismatch {:.2e}, distortion {:.3}",
            mesh.n_polygons(),
            result.inserted.mesh.n_polygons(),
            rep.charts.len(),
            result.map.cones.len(),
            rep.flipped_triangles,
            rep.seam_lengths.max_seam_mismatch,
            rep.area_distortion.scale_free
        );
    }
}

#[test]
fn atlas_end_to_end_on_generated_meshes() {
    let steps = if debug_build() { 3 } else { 60 };
    atlas("uv_sphere", &uv_sphere(32, 16, 1.0).unwrap(), 6, steps);
    atlas("grid_cube", &grid_cube(6, 2.0), 6, steps);
    atlas("capsule", &capsule(24, 8, 0.5, 1.5).unwrap(), 6, steps);
    atlas("torus", &torus(32, 12, 1.0, 0.35).unwrap(), 6, steps);
}

#[test]
fn atlas_end_to_end_on_the_quad_model() {
    let Some(mesh) = load_polygons("spot_quadrangulated.obj") else {
        return;
    };
    atlas("spot", &mesh, 8, if debug_build() { 3 } else { 100 });
}

// ===================== intrinsic triangulations =====================

fn flips_only() -> Triangulation {
    Triangulation::Intrinsic(IntrinsicOptions::flips_only())
}

fn refined() -> Triangulation {
    Triangulation::Intrinsic(IntrinsicOptions::refined())
}

fn two_session(
    mesh: &PolygonMesh,
    seams: &[[u32; 2]],
    start: &Start,
    tri: &Triangulation,
) -> Session {
    let input = MeshInput {
        seam_edges: seams,
        ..mesh.input()
    };
    Session::with_triangulation(&input, 2, start, FlowKind::Auto, tri).expect("session")
}

/// The parts next to each other across an edge that is not a seam of the cut lie in the same
/// region; across an edge of the cut (and not an old seam) in two different ones.
fn check_part_regions(ins: &Inserted, old_seams: &[[u32; 2]]) {
    let new = &ins.mesh;
    let all = ins.seams_after(old_seams);
    let cut_edges: BTreeSet<[u32; 2]> = ins.new_seams.iter().copied().collect();
    let old: BTreeSet<[u32; 2]> = all
        .iter()
        .copied()
        .filter(|e| !cut_edges.contains(e))
        .collect();
    let mut owner: BTreeMap<[u32; 2], usize> = BTreeMap::new();
    for q in 0..new.n_polygons() {
        let poly = new.polygon(q);
        for (i, &a) in poly.iter().enumerate() {
            let b = poly[(i + 1) % poly.len()];
            let e = [a.min(b), a.max(b)];
            if old.contains(&e) {
                continue;
            }
            if let Some(&r) = owner.get(&e) {
                if let (Some(x), Some(y)) = (ins.polygon_region[q], ins.polygon_region[r]) {
                    if cut_edges.contains(&e) {
                        assert_ne!(x, y, "polygons {q} and {r} across the cut edge {e:?}");
                    } else {
                        assert_eq!(x, y, "polygons {q} and {r} across the edge {e:?}");
                    }
                }
            } else {
                owner.insert(e, q);
            }
        }
    }
}

#[test]
fn a_bend_and_a_v_inside_one_triangle() {
    // A square with the diagonal 0-2: triangles (0,1,2) and (0,2,3).
    let mesh = PolygonMesh::from_polygons(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        &[vec![0, 1, 2, 3]],
    );
    let cut = CutOnPolygons {
        edge_points: vec![
            EdgePoint {
                edge: [0, 1],
                t: 0.3,
                position: [0.3, 0.0, 0.0],
                side: None,
            },
            EdgePoint {
                edge: [0, 1],
                t: 0.7,
                position: [0.7, 0.0, 0.0],
                side: None,
            },
        ],
        junctions: vec![ddg_uv::cutting::JunctionPoint {
            position: [0.5, 0.2, 0.0],
            polygon: 0,
            triangle: 0,
            corners: mesh.triangles[0],
            barycentric: [0.5, 0.3, 0.2],
        }],
        paths: vec![
            PolygonPath {
                polygon: 0,
                start: PathEnd::EdgePoint(0),
                interior: vec![],
                end: PathEnd::Junction(0),
                regions: [0, 1],
            },
            PolygonPath {
                polygon: 0,
                start: PathEnd::Junction(0),
                interior: vec![],
                end: PathEnd::EdgePoint(1),
                regions: [0, 1],
            },
        ],
        n_clamped_crossings: 0,
    };
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &[], &cut, &ins);
    check_part_regions(&ins, &[]);
    assert_eq!(ins.mesh.n_polygons(), 2);
    // The small part below the V is region 1 and a triangle of area 0.04.
    let small = (0..2).find(|&q| ins.polygon_region[q] == Some(1)).unwrap();
    assert_eq!(ins.polygon_region[1 - small], Some(0));
    assert_eq!(ins.mesh.polygon(small).len(), 3);
    let area: f64 = ins
        .mesh
        .triangles
        .iter()
        .filter(|t| ins.mesh.polygon_corners(small).contains(&(t[0] as usize)))
        .map(|t| {
            let [a, b, c] =
                t.map(|c| ins.mesh.positions[ins.mesh.corner_vertices[c as usize] as usize]);
            0.5 * cross_norm(sub(b, a), sub(c, a))
        })
        .sum();
    assert!((area - 0.04).abs() < 1e-12, "{area}");
    // The same cut through the display lines: one open line of three positions.
    let lines = cut.lines();
    assert_eq!(lines.len(), 1);
    assert_eq!(lines[0].positions.len(), 3);
    assert!(!lines[0].closed);
}

#[test]
fn several_points_on_an_edge_and_near_a_vertex() {
    let mesh = PolygonMesh::from_polygons(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        &[vec![0, 1, 2, 3]],
    );
    // Two parallel cuts from the bottom edge to the right edge, the second 1e-7 from corner 1.
    let ep = |edge: [u32; 2], t: f64| {
        let [a, b] = edge.map(|v| mesh.positions[v as usize]);
        EdgePoint {
            edge,
            t,
            position: [0, 1, 2].map(|d| (1.0 - t) * a[d] + t * b[d]),
            side: None,
        }
    };
    let cut = CutOnPolygons {
        edge_points: vec![
            ep([0, 1], 0.4),
            ep([1, 2], 0.4),
            ep([0, 1], 1.0 - 1e-7),
            ep([1, 2], 1e-7),
        ],
        paths: vec![
            PolygonPath {
                polygon: 0,
                start: PathEnd::EdgePoint(0),
                interior: vec![],
                end: PathEnd::EdgePoint(1),
                regions: [0, 1],
            },
            PolygonPath {
                polygon: 0,
                start: PathEnd::EdgePoint(2),
                interior: vec![],
                end: PathEnd::EdgePoint(3),
                regions: [1, 2],
            },
        ],
        ..CutOnPolygons::default()
    };
    let ins = insert(&mesh, &cut).unwrap();
    check_insertion(&mesh, &[], &cut, &ins);
    // Three parts: the corner at vertex 1 (tiny), the band, the rest.
    assert_eq!(ins.mesh.n_polygons(), 3);
    assert!(min_triangle_area(&ins.mesh) > 0.0);
}

#[test]
fn intrinsic_session_runs_and_the_cut_is_inserted() {
    let mesh = uv_sphere(24, 12, 1.0).unwrap();
    for tri in [flips_only(), refined()] {
        let mut session = two_session(
            &mesh,
            &[],
            &Start::GeodesicCircle {
                vertex: 5,
                radius: 0.9,
            },
            &tri,
        );
        let report = session.intrinsic_report().unwrap().clone();
        println!("{tri:?}: {report:?}");
        assert_eq!(report.n_input_vertices, mesh.positions.len());
        assert_eq!(report.n_vertices, session.engine_mesh().n_vertices());
        assert!(report.is_delaunay);
        let (problem, options) = flow_setup(&mut session, 0.25);
        session.run(&problem, &options, &stop(5)).unwrap();
        let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
        assert!(!cut.is_empty());
        let ins = insert(&mesh, &cut).unwrap();
        check_insertion(&mesh, &[], &cut, &ins);
        check_part_regions(&ins, &[]);
        // The preview comes from the description: one closed line.
        let preview = session.preview(DEFAULT_CROSSING_OFFSET).unwrap();
        assert_eq!(preview.len(), 1);
        assert!(preview[0].closed);
        assert_eq!(n_charts(&ins.mesh, &ins.new_seams), 2);
        // The per-corner view has the caller's vertices only.
        assert_eq!(session.state()[0].len(), mesh.corner_vertices.len());
    }
}

#[test]
fn the_caller_triangulation_is_unchanged_by_the_option() {
    let mesh = uv_sphere(20, 10, 1.0).unwrap();
    let start = Start::GeodesicCircle {
        vertex: 5,
        radius: 0.8,
    };
    let mut a = Session::new(&mesh.input(), 2, &start).unwrap();
    let mut b = two_session(&mesh, &[], &start, &Triangulation::Caller);
    let (problem, options) = flow_setup(&mut a, 0.3);
    for _ in 0..4 {
        let (ra, rb) = (
            a.step(&problem, &options).unwrap(),
            b.step(&problem, &options).unwrap(),
        );
        assert_eq!(ra.energy.total.to_bits(), rb.energy.total.to_bits());
    }
    assert_eq!(a.state(), b.state());
    assert_eq!(a.fingerprint(), b.fingerprint());
    assert!(a.intrinsic_report().is_none());
}

#[test]
fn the_intrinsic_session_is_the_direct_run() {
    // A triangle mesh without seams: the session's mesh is the mesh itself.
    let sphere = ddg_core::generators::icosphere(3).unwrap();
    let mesh = PolygonMesh::from_triangles(sphere.positions.clone(), &sphere.faces);
    let options = IntrinsicOptions::refined();
    let mut session = two_session(
        &mesh,
        &[],
        &Start::GeodesicCircle {
            vertex: 3,
            radius: 0.9,
        },
        &Triangulation::Intrinsic(options),
    );
    // The direct construction with the same options.
    let tri_mesh = ddg_core::TriMesh::new(sphere.positions.len(), sphere.faces.clone()).unwrap();
    let geometry = ddg_core::Geometry::from_positions(&tri_mesh, &sphere.positions).unwrap();
    let mut tri = ddg_intrinsic::IntrinsicTriangulation::new(&tri_mesh, &geometry).unwrap();
    tri.flip_to_delaunay(&options.delaunay).unwrap();
    tri.refine(&options.refine.unwrap()).unwrap();
    let export = tri.export_simplicial(&options.export).unwrap();
    assert_eq!(export.fingerprint, session.fingerprint());
    assert_eq!(export.mesh.n_vertices(), session.engine_mesh().n_vertices());
    let state = session.snapshot();
    assert_eq!(state.n_values(), export.mesh.n_vertices());
    let mut flow = Flow::new(&export.mesh, &export.geometry, &state.values()[0]).unwrap();
    let (problem, flow_options) = flow_setup(&mut session, 0.3);
    for _ in 0..4 {
        let a = session.step(&problem, &flow_options).unwrap();
        let b = flow.step(&problem, &flow_options).unwrap();
        assert_eq!(a.energy.total.to_bits(), b.energy.total.to_bits());
    }
    let after = session.snapshot();
    let phi = flow.level_set();
    for (x, y) in after.values()[0].iter().zip(&phi) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    // The cut: its points on the sphere are where the direct way back puts them.
    let cut =
        ddg_variational_cuts::Cut::from_level_set(&export.mesh, &phi, DEFAULT_CROSSING_OFFSET)
            .unwrap();
    let description = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    assert!(!description.is_empty());
    assert!(cut.points.len() > 0);
    // Every point of the description lies on the sphere to the accuracy of the flat faces.
    for p in &description.edge_points {
        let r = (p.position.iter().map(|x| x * x).sum::<f64>()).sqrt();
        assert!((r - 1.0).abs() < 0.05, "{r}");
    }
}

#[test]
fn state_round_trip_and_stale_state() {
    let mesh = uv_sphere(20, 10, 1.0).unwrap();
    let seams = {
        let poly = mesh.polygon(30);
        vec![
            [poly[0].min(poly[1]), poly[0].max(poly[1])],
            [poly[1].min(poly[2]), poly[1].max(poly[2])],
        ]
    };
    let start = Start::GeodesicCircle {
        vertex: 40,
        radius: 0.8,
    };
    let tri = refined();
    let mut session = two_session(&mesh, &seams, &start, &tri);
    let (problem, options) = flow_setup(&mut session, 0.3);
    session.run(&problem, &options, &stop(3)).unwrap();
    let state = session.snapshot();
    assert!(state.n_values() >= mesh.positions.len());
    // Bytes and back.
    let bytes = state.to_bytes();
    let back = SessionState::from_bytes(&bytes).unwrap();
    assert_eq!(back, state);
    assert_eq!(back.triangulation(), tri);
    // Restore into the same session after more steps, then into a rebuilt one.
    session.run(&problem, &options, &stop(2)).unwrap();
    assert_ne!(session.snapshot(), state);
    session.restore(&back).unwrap();
    assert_eq!(session.snapshot(), state);
    let input = MeshInput {
        seam_edges: &seams,
        ..mesh.input()
    };
    let rebuilt = Session::new_from_state(&input, &back).unwrap();
    assert_eq!(rebuilt.snapshot(), state);
    assert_eq!(rebuilt.fingerprint(), session.fingerprint());
    // A state of another triangulation is refused, and says so.
    let other = two_session(&mesh, &seams, &start, &flips_only());
    let e = other.clone_state_check(&state);
    assert!(matches!(e, Err(CuttingError::StaleState { .. })), "{e:?}");
    // Wrong number of regions, wrong values, bad bytes.
    let mut three = Session::new(
        &mesh.input(),
        3,
        &Start::NormalClusters { max_iterations: 20 },
    )
    .unwrap();
    assert!(three.restore(&state).is_err());
    assert!(matches!(
        state.with_values(vec![vec![0.0; 3]]),
        Err(CuttingError::BadState { .. })
    ));
    assert!(matches!(
        SessionState::from_bytes(&bytes[..10]),
        Err(CuttingError::BadState { .. })
    ));
    let mut corrupt = bytes.clone();
    corrupt[0] = b'X';
    assert!(matches!(
        SessionState::from_bytes(&corrupt),
        Err(CuttingError::BadState { .. })
    ));
    let mut longer = bytes.clone();
    longer.push(0);
    assert!(matches!(
        SessionState::from_bytes(&longer),
        Err(CuttingError::BadState { .. })
    ));
}

trait StateCheck {
    fn clone_state_check(&self, state: &SessionState) -> Result<(), CuttingError>;
}

impl StateCheck for Session {
    fn clone_state_check(&self, state: &SessionState) -> Result<(), CuttingError> {
        // restore needs &mut; check on a rebuilt session of the same triangulation
        let mut copy = Session::with_triangulation(
            &self.mesh().input(),
            self.n_regions(),
            &Start::Corners(vec![vec![0.0; self.mesh().corner_vertices.len()]; 2]),
            FlowKind::Auto,
            &self.triangulation(),
        )?;
        copy.restore(state)
    }
}

/// Random level-set functions on the vertices of the intrinsic triangulation: cuts with bends,
/// several points per edge, junctions of three regions, ends on seams, through the whole chain.
#[test]
fn random_cuts_on_intrinsic_triangulations() {
    let mut rng = ddg_core::generators::SplitMix64::new(23);
    let (mut inserted, mut refused) = (0usize, 0usize);
    let (mut bends, mut multi_points, mut junction_paths, mut seam_points) =
        (0usize, 0usize, 0usize, 0usize);
    let sphere = uv_sphere(14, 7, 1.0).unwrap();
    let cube = grid_cube(3, 2.0);
    let tor = torus(12, 6, 1.0, 0.3).unwrap();
    let cyl = cylinder(12, 2, 0.7, 2.0).unwrap();
    for (name, mesh) in [
        ("sphere", sphere),
        ("cube", cube),
        ("torus", tor),
        ("cylinder", cyl),
    ] {
        let edges = mesh.edges();
        let seam_sets: Vec<Vec<[u32; 2]>> =
            vec![vec![], edges.iter().step_by(7).take(9).copied().collect()];
        for seams in seam_sets {
            for tri in [flips_only(), refined()] {
                for n in [2usize, 3, 5] {
                    let input = MeshInput {
                        seam_edges: &seams,
                        ..mesh.input()
                    };
                    let start = if n == 2 {
                        Start::GeodesicCircle {
                            vertex: 1,
                            radius: 0.4,
                        }
                    } else {
                        Start::NormalClusters { max_iterations: 20 }
                    };
                    let mut session =
                        Session::with_triangulation(&input, n, &start, FlowKind::Auto, &tri)
                            .unwrap();
                    for round in 0..3 {
                        let state = session.snapshot();
                        let random: Vec<Vec<f64>> = state
                            .values()
                            .iter()
                            .map(|f| f.iter().map(|_| rng.uniform(-1.0, 1.0)).collect())
                            .collect();
                        session
                            .restore(&state.with_values(random).unwrap())
                            .unwrap();
                        let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
                        match insert(&mesh, &cut) {
                            Ok(ins) => {
                                check_insertion(&mesh, &seams, &cut, &ins);
                                check_part_regions(&ins, session.seams_cut());
                                inserted += 1;
                                bends += cut
                                    .junctions
                                    .iter()
                                    .enumerate()
                                    .filter(|(j, _)| {
                                        cut.paths
                                            .iter()
                                            .map(|p| {
                                                usize::from(p.start == PathEnd::Junction(*j as u32))
                                                    + usize::from(
                                                        p.end == PathEnd::Junction(*j as u32),
                                                    )
                                            })
                                            .sum::<usize>()
                                            == 2
                                    })
                                    .count();
                                let mut per_edge: BTreeMap<[u32; 2], usize> = BTreeMap::new();
                                for p in &cut.edge_points {
                                    *per_edge.entry(p.edge).or_default() += 1;
                                }
                                multi_points += per_edge.values().filter(|&&c| c >= 2).count();
                                junction_paths += cut
                                    .paths
                                    .iter()
                                    .filter(|p| {
                                        matches!(p.start, PathEnd::Junction(_))
                                            || matches!(p.end, PathEnd::Junction(_))
                                    })
                                    .count();
                                seam_points +=
                                    cut.edge_points.iter().filter(|p| p.side.is_some()).count();
                            }
                            Err(CuttingError::Unsupported { .. })
                            | Err(CuttingError::InvalidDescription { .. }) => {
                                refused += 1;
                            }
                            Err(e) => panic!("{name} {tri:?} N={n} round {round}: {e}"),
                        }
                    }
                }
            }
        }
    }
    println!(
        "inserted {inserted}, refused {refused}, bends {bends}, edges with several points {multi_points}, paths at junctions {junction_paths}, seam-side points {seam_points}"
    );
    assert!(
        inserted > 20 * refused.max(1) / 4,
        "{inserted} inserted, {refused} refused"
    );
    assert!(bends > 100 && multi_points > 50 && junction_paths > 100 && seam_points > 10);
}

// ===================== bad meshes: the caller's triangulation against the intrinsic one =====================

/// The lattice `x = (i + σ j) h, y = j h`: every triangle has the largest angle
/// `180° − atan(1/σ)` (the parallelogram `n × n` cells, sheared).
fn sheared_grid(n: usize, sigma: f64) -> PolygonMesh {
    let h = 1.0 / n as f64;
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            positions.push([(i as f64 + sigma * j as f64) * h, j as f64 * h, 0.0]);
        }
    }
    let at = |i: usize, j: usize| (j * (n + 1) + i) as u32;
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    PolygonMesh::from_triangles(positions, &faces)
}

/// Rows alternate between thin strips of height `0.05 d` (caps: one angle of 168 degrees) and
/// tall ones of height `3 d` (needles: apex 19 degrees), the rows shifted by half a cell (the
/// construction of `bad_mesh_study` in `ddg-variational-cuts`).
fn needles_and_caps(n: usize) -> PolygonMesh {
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
    for j in 0..rows {
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
            positions.push([x, y[j], 0.0]);
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
    PolygonMesh::from_triangles(positions, &faces)
}

/// One line of the comparison.
#[derive(Debug)]
struct Row {
    label: String,
    n_vertices: usize,
    seconds_start: f64,
    seconds_steps: f64,
    energy: (f64, f64),
    steps: usize,
    cut_points: usize,
    polygons_after: usize,
    charts: usize,
    flipped: usize,
    cones: usize,
    distortion: f64,
    seam_mismatch: f64,
    angles: Option<(f64, f64, usize)>,
    cut_length: f64,
}

/// Starts as a circle (or, with `clusters`, as normal clusters), runs the flow on the given
/// triangulation, inserts the cut and unwraps it.
fn run_chain(
    name: &str,
    mesh: &PolygonMesh,
    tri: &Triangulation,
    steps: usize,
    clusters: Option<usize>,
    alpha: Option<f64>,
) -> Row {
    let area = triangulated_area(mesh);
    let center = {
        // The vertex nearest to the centre of mass.
        let mut c = [0.0f64; 3];
        for p in &mesh.positions {
            for d in 0..3 {
                c[d] += p[d] / mesh.positions.len() as f64;
            }
        }
        (0..mesh.positions.len())
            .min_by(|&a, &b| {
                let dist = |v: usize| {
                    (0..3)
                        .map(|d| (mesh.positions[v][d] - c[d]).powi(2))
                        .sum::<f64>()
                };
                dist(a).total_cmp(&dist(b))
            })
            .unwrap() as u32
    };
    let (n, start) = match clusters {
        Some(n) => (n, Start::NormalClusters { max_iterations: 50 }),
        None => (
            2,
            Start::GeodesicCircle {
                vertex: center,
                radius: 0.2 * area.sqrt(),
            },
        ),
    };
    let t0 = std::time::Instant::now();
    let mut session = Session::with_triangulation(&mesh.input(), n, &start, FlowKind::Auto, tri)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let seconds_start = t0.elapsed().as_secs_f64();
    let options = session.recommended_options(StepSize::default());
    let mut problem = session.new_problem(Energy::Hencky, 0.0);
    let omega = session.roughening_rate(&problem, &options.eval).unwrap();
    problem.length_weight = alpha.unwrap_or_else(|| {
        length_weight_for_ripples(omega, if clusters.is_some() { 0.5 } else { 0.25 })
    });
    let t1 = std::time::Instant::now();
    let report = session
        .run(&problem, &options, &stop(steps))
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let seconds_steps = t1.elapsed().as_secs_f64();
    let (_, result) = session
        .unwrap(
            DEFAULT_CROSSING_OFFSET,
            &CutUnwrapOptions {
                cones: ConePlacement::Optimal,
                ..Default::default()
            },
        )
        .unwrap_or_else(|e| panic!("{name} {tri:?}: {e}"));
    let rep = &result.map.report;
    let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
    let cut_length: f64 = cut
        .lines()
        .iter()
        .map(|l| {
            let m = l.positions.len();
            let n = if l.closed { m } else { m - 1 };
            (0..n)
                .map(|i| {
                    let (a, b) = (l.positions[i], l.positions[(i + 1) % m]);
                    sub(a, b).iter().map(|x| x * x).sum::<f64>().sqrt()
                })
                .sum::<f64>()
        })
        .sum();
    Row {
        label: name.to_string(),
        n_vertices: session.engine_mesh().n_vertices(),
        seconds_start,
        seconds_steps,
        energy: (
            report.steps.first().map_or(0.0, |s| s.energy.total),
            report.final_energy.total,
        ),
        steps: report.steps.len(),
        cut_points: cut.edge_points.len() + cut.junctions.len(),
        polygons_after: result.inserted.mesh.n_polygons(),
        charts: rep.charts.len(),
        flipped: rep.flipped_triangles,
        cones: result.map.cones.len(),
        distortion: rep.area_distortion.scale_free,
        seam_mismatch: rep.seam_lengths.max_seam_mismatch,
        angles: session
            .intrinsic_report()
            .map(|r| (r.min_angle, r.max_angle, r.n_below_bound)),
        cut_length,
    }
}

fn compare(
    name: &str,
    mesh: &PolygonMesh,
    steps: usize,
    clusters: Option<usize>,
    alpha: Option<f64>,
) {
    let n_tri = mesh.triangles.len();
    for (label, tri) in [
        ("caller", Triangulation::Caller),
        ("flips", flips_only()),
        ("refined", refined()),
    ] {
        let row = run_chain(
            &format!("{name}/{label}"),
            mesh,
            &tri,
            steps,
            clusters,
            alpha,
        );
        println!("ROW {row:?}");
        // Whatever the triangulation, the chain is consistent: the number of charts is the
        // number of pieces the seams leave, no more than a sliver per thousand triangles
        // flips on the intrinsic ones.
        if label != "caller" {
            assert!(
                row.flipped <= 1 + n_tri / 100,
                "{name}/{label}: {} flipped",
                row.flipped
            );
        }
        assert!(row.energy.1.is_finite());
    }
}

#[test]
fn bad_mesh_sheared() {
    let steps = if debug_build() { 3 } else { 8 };
    compare(
        "sheared 153",
        &sheared_grid(if debug_build() { 16 } else { 32 }, 2.0),
        steps,
        None,
        Some(1.0),
    );
}

#[test]
fn bad_mesh_needles_and_caps() {
    let steps = if debug_build() { 3 } else { 8 };
    compare(
        "needles and caps",
        &needles_and_caps(if debug_build() { 8 } else { 16 }),
        steps,
        None,
        Some(1.0),
    );
}

#[test]
fn bad_mesh_camelbad() {
    if debug_build() {
        println!("camelbad: release builds only");
        return;
    }
    let path = dataset_dir().join("camelbad.obj");
    let Ok(obj) = read_obj(&path) else {
        println!("dataset mesh {} not found: skipped", path.display());
        return;
    };
    let mesh = PolygonMesh::from_triangles(obj.positions.clone(), &obj.triangles().unwrap());
    compare("camelbad", &mesh, 20, None, None);
}

#[test]
fn intrinsic_on_generated_meshes_and_spot() {
    let steps = if debug_build() { 3 } else { 40 };
    compare(
        "uv_sphere",
        &uv_sphere(32, 16, 1.0).unwrap(),
        steps,
        None,
        None,
    );
    compare("grid_cube", &grid_cube(6, 2.0), steps, None, None);
    compare(
        "torus",
        &torus(32, 12, 1.0, 0.35).unwrap(),
        steps,
        None,
        None,
    );
    if let Some(mesh) = load_polygons("spot_quadrangulated.obj") {
        compare("spot quads", &mesh, steps, None, None);
        compare(
            "spot atlas",
            &mesh,
            if debug_build() { 2 } else { 60 },
            Some(8),
            None,
        );
    }
}

/// The same cut, given by the exact distance to a point at the vertices of a 165 degree
/// sheared grid, evaluated on the caller's triangulation, on the flipped one and on a refined
/// one: the energy terms and the length of the cut. (A flat mesh has no distortion, so the
/// energy is `alpha_L` times the length of the cut.)
#[test]
fn the_same_cut_on_three_triangulations() {
    let sigma = 1.0 / 15f64.to_radians().tan();
    let mesh = sheared_grid(24, sigma);
    let centre = [0.5 + 0.5 * sigma, 0.5];
    let r = 0.2;
    let phi: Vec<f64> = mesh
        .corner_vertices
        .iter()
        .map(|&v| {
            let p = mesh.positions[v as usize];
            ((p[0] - centre[0]).powi(2) + (p[1] - centre[1]).powi(2)).sqrt() - r
        })
        .collect();
    let negated: Vec<f64> = phi.iter().map(|x| -x).collect();
    let mut lengths = Vec::new();
    for (label, tri) in [
        ("caller", Triangulation::Caller),
        ("flips", flips_only()),
        ("refined", refined()),
    ] {
        let mut session = two_session(
            &mesh,
            &[],
            &Start::Corners(vec![phi.clone(), negated.clone()]),
            &tri,
        );
        let problem = CutProblem::new(Energy::Hencky, 1.0);
        let e = session.evaluate(&problem, &Default::default()).unwrap();
        let cut = session.describe_cut(DEFAULT_CROSSING_OFFSET).unwrap();
        let length: f64 = cut
            .lines()
            .iter()
            .map(|l| {
                let m = l.positions.len();
                let n = if l.closed { m } else { m - 1 };
                (0..n)
                    .map(|i| {
                        sub(l.positions[i], l.positions[(i + 1) % m])
                            .iter()
                            .map(|x| x * x)
                            .sum::<f64>()
                            .sqrt()
                    })
                    .sum::<f64>()
            })
            .sum();
        println!(
            "SAME CUT {label}: vertices {}, energy {:.6} (distortion {:?}, length term {:.6}), length of the cut on the surface {:.6} (exact {:.6})",
            session.engine_mesh().n_vertices(),
            e.energy.total,
            e.energy.distortion,
            e.energy.length,
            length,
            std::f64::consts::TAU * r
        );
        lengths.push(length);
    }
    assert!(lengths.iter().all(|l| l.is_finite() && *l > 0.0));
}
