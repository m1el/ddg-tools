//! The simplicial export, the way back of segments and the optimal Delaunay iteration on
//! integer coordinates (SPEC §4.10, §4.11, I29), against the signpost representation.
//!
//! Run with `--nocapture` for the tables.

mod common;

use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{Shape, all_shapes, bad_meshes, inner_bary, models};
use ddg_core::generators::SplitMix64;
use ddg_integer_coords::{IntegerOdtOptions, IntegerRefineOptions, IntegerTriangulation};
use ddg_intrinsic::{
    DelaunayOptions, EdgeTraceOptions, ExportOptions, IntrinsicError, IntrinsicTriangulation,
    OptimalDelaunayOptions, PathsOnInput, RefineOptions, RefineRules, Segment, SimplicialExport,
    SurfacePoint, Triangulation,
};

fn par_map<S: Sync, R: Send>(specs: &[S], f: impl Fn(&S) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(32);
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Option<R>>> = Mutex::new((0..specs.len()).map(|_| None).collect());
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= specs.len() {
                        break;
                    }
                    let r = f(&specs[i]);
                    out.lock().unwrap()[i] = Some(r);
                }
            });
        }
    });
    out.into_inner()
        .unwrap()
        .into_iter()
        .map(|r| r.unwrap())
        .collect()
}

fn position(shape: &Shape, p: SurfacePoint) -> [f64; 3] {
    let x = |v: usize| shape.positions[v];
    match p {
        SurfacePoint::Vertex(v) => x(v as usize),
        SurfacePoint::Edge { edge, t } => {
            let [a, b] = shape.mesh.edge_vertices(edge as usize);
            let (pa, pb) = (x(a), x(b));
            [0, 1, 2].map(|c| (1.0 - t) * pa[c] + t * pb[c])
        }
        SurfacePoint::Face { face, bary } => {
            let f = shape.mesh.faces()[face as usize];
            let p = f.map(|v| x(v as usize));
            [0, 1, 2].map(|c| bary[0] * p[0][c] + bary[1] * p[1][c] + bary[2] * p[2][c])
        }
    }
}

fn distance(shape: &Shape, a: SurfacePoint, b: SurfacePoint) -> f64 {
    common::dist(position(shape, a), position(shape, b))
}

/// The largest distance of a point of path `b` from the polyline of path `a`, in space.
fn path_distance(shape: &Shape, a: &[SurfacePoint], b: &[SurfacePoint]) -> f64 {
    let pa: Vec<[f64; 3]> = a.iter().map(|&x| position(shape, x)).collect();
    b.iter()
        .map(|&x| {
            let q = position(shape, x);
            pa.windows(2)
                .map(|w| {
                    let d = [0, 1, 2].map(|c| w[1][c] - w[0][c]);
                    let r = [0, 1, 2].map(|c| q[c] - w[0][c]);
                    let dd: f64 = d.iter().map(|x| x * x).sum();
                    let t = if dd > 0.0 {
                        ((0..3).map(|c| r[c] * d[c]).sum::<f64>() / dd).clamp(0.0, 1.0)
                    } else {
                        0.0
                    };
                    common::dist(q, [0, 1, 2].map(|c| w[0][c] + t * d[c]))
                })
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max)
}

/// What is done before the export.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Before {
    /// Delaunay flips.
    Flips,
    /// Delaunay flips and refinement to 25° (the signposts' rules on both).
    Refined,
}

fn refine_options() -> RefineOptions {
    RefineOptions {
        max_insertions: ddg_intrinsic::InsertionLimit::Vertices(4000),
        ..RefineOptions::default()
    }
}

fn export_options(before: Before) -> ExportOptions {
    match before {
        Before::Flips => ExportOptions::default(),
        Before::Refined => ExportOptions {
            refine: Some(refine_options()),
            ..ExportOptions::default()
        },
    }
}

/// Both representations of `shape`, prepared and exported the same way.
#[allow(clippy::type_complexity)]
fn both(
    shape: &Shape,
    before: Before,
) -> Result<
    (
        (IntrinsicTriangulation, SimplicialExport),
        (IntegerTriangulation, SimplicialExport),
    ),
    String,
> {
    let mut s = shape.triangulation();
    let mut i =
        IntegerTriangulation::new(&shape.mesh, &shape.geometry).map_err(|e| e.to_string())?;
    s.flip_to_delaunay(&DelaunayOptions::default())
        .map_err(|e| e.to_string())?;
    Triangulation::flip_to_delaunay(&mut i, &DelaunayOptions::default())
        .map_err(|e| e.to_string())?;
    if before == Before::Refined {
        s.refine(&refine_options()).map_err(|e| e.to_string())?;
        i.refine(&IntegerRefineOptions {
            refine: refine_options(),
            rules: RefineRules::Ruppert,
            ..Default::default()
        })
        .map_err(|e| e.to_string())?;
    }
    let es = s
        .export_simplicial(&export_options(before))
        .map_err(|e| format!("signposts: {e}"))?;
    let ei = i
        .export_simplicial(&export_options(before))
        .map_err(|e| format!("integers: {e}"))?;
    Ok(((s, es), (i, ei)))
}

/// How the two exports compare.
#[derive(Debug, Default)]
struct ExportComparison {
    same_connectivity: bool,
    same_lengths: bool,
    same_counts: bool,
    max_location_difference: f64,
    n_vertices: usize,
    n_splits: usize,
    n_off_midpoint: usize,
}

fn compare_exports(shape: &Shape, a: &SimplicialExport, b: &SimplicialExport) -> ExportComparison {
    let same_connectivity = a.mesh.faces() == b.mesh.faces();
    let same_lengths = same_connectivity
        && a.geometry
            .edge_lengths()
            .iter()
            .zip(b.geometry.edge_lengths())
            .all(|(x, y)| x.to_bits() == y.to_bits());
    let same_counts = a.n_splits == b.n_splits
        && a.n_flips == b.n_flips
        && a.n_refine_rounds == b.n_refine_rounds
        && a.n_refined == b.n_refined
        && a.is_delaunay == b.is_delaunay;
    let max_location_difference = if a.locations.len() == b.locations.len() {
        a.locations
            .iter()
            .zip(&b.locations)
            .map(|(&p, &q)| distance(shape, p, q))
            .fold(0.0, f64::max)
    } else {
        f64::INFINITY
    };
    ExportComparison {
        same_connectivity,
        same_lengths,
        same_counts,
        max_location_difference,
        n_vertices: a.mesh.n_vertices(),
        n_splits: a.n_splits,
        n_off_midpoint: b.n_off_midpoint_splits,
    }
}

fn shapes() -> Vec<Shape> {
    let mut shapes = all_shapes();
    shapes.extend(models());
    shapes.extend(bad_meshes());
    shapes
}

#[test]
fn the_export_is_the_signposts_export_up_to_the_places() {
    let shapes = shapes();
    let specs: Vec<(usize, Before)> = (0..shapes.len())
        .flat_map(|i| [(i, Before::Flips), (i, Before::Refined)])
        .collect();
    let results = par_map(&specs, |&(i, before)| {
        let shape = &shapes[i];
        both(shape, before).map(|((_, es), (it, ei))| {
            let fingerprint_ok =
                ei.fingerprint == it.fingerprint() && ei.fingerprint != es.fingerprint;
            let certified = it.certify().is_ok();
            (compare_exports(shape, &es, &ei), fingerprint_ok, certified)
        })
    });
    println!(
        "{:<36} {:<8} {:>6} {:>6} {:>4} {:>5} {:>5} {:>5} {:>10}",
        "shape", "before", "verts", "splits", "off", "conn", "len", "count", "max |Δq|"
    );
    let (mut n_same, mut n_runs) = (0, 0);
    for (&(i, before), r) in specs.iter().zip(&results) {
        let shape = &shapes[i];
        match r {
            Ok((c, fingerprint_ok, certified)) => {
                n_runs += 1;
                println!(
                    "{:<36} {:<8} {:>6} {:>6} {:>4} {:>5} {:>5} {:>5} {:>10.2e}",
                    shape.name,
                    format!("{before:?}"),
                    c.n_vertices,
                    c.n_splits,
                    c.n_off_midpoint,
                    c.same_connectivity,
                    c.same_lengths,
                    c.same_counts,
                    c.max_location_difference / shape.mean_edge_length()
                );
                assert!(fingerprint_ok, "{}: fingerprint", shape.name);
                assert!(
                    certified,
                    "{}: the exported integers do not certify",
                    shape.name
                );
                if c.same_connectivity && c.same_lengths && c.same_counts {
                    n_same += 1;
                    assert!(
                        c.max_location_difference <= 1e-9 * shape.mean_edge_length(),
                        "{} {before:?}: places differ by {:e}",
                        shape.name,
                        c.max_location_difference
                    );
                }
                if before == Before::Flips && c.n_off_midpoint == 0 {
                    // Flips and midpoint splits depend on the lengths alone.
                    assert!(
                        c.same_connectivity && c.same_lengths && c.same_counts,
                        "{}: the flips-only export differs",
                        shape.name
                    );
                }
            }
            Err(e) => println!("{:<36} {:<8} error: {e}", shape.name, format!("{before:?}")),
        }
    }
    println!("{n_same} of {n_runs} exports identical in connectivity, lengths and counts");
    assert!(n_runs > 0);
}

/// Segments of every kind in some faces of an export: between points of two sides, from a
/// point inside to a side, from a corner to the opposite side or to a point inside, along a
/// side.
fn random_segments(
    export: &SimplicialExport,
    rng: &mut SplitMix64,
    max_faces: usize,
) -> (Vec<SurfacePoint>, Vec<Segment>) {
    let mesh = &export.mesh;
    let nf = mesh.n_faces();
    let mut nodes = Vec::new();
    let mut segments = Vec::new();
    let edge_node = |nodes: &mut Vec<SurfacePoint>, h: usize, t: f64| {
        let e = mesh.edge(h);
        nodes.push(export.edge_point(e, t));
        (nodes.len() - 1) as u32
    };
    let step = (nf / max_faces).max(1);
    for f in (0..nf).step_by(step) {
        let face = f as u32;
        let t = |rng: &mut SplitMix64| rng.uniform(0.05, 0.95);
        let a = edge_node(&mut nodes, 3 * f, t(rng));
        let b = edge_node(&mut nodes, 3 * f + 1, t(rng));
        let c = edge_node(&mut nodes, 3 * f + 2, t(rng));
        nodes.push(SurfacePoint::Face {
            face,
            bary: inner_bary(rng),
        });
        let m = (nodes.len() - 1) as u32;
        let corner = mesh.faces()[f][0];
        nodes.push(SurfacePoint::Vertex(corner));
        let v = (nodes.len() - 1) as u32;
        // Two points of the same side, for a segment along it.
        let (t1, t2) = (rng.uniform(0.05, 0.45), rng.uniform(0.55, 0.95));
        let s1 = edge_node(&mut nodes, 3 * f + 1, t1);
        let s2 = edge_node(&mut nodes, 3 * f + 1, t2);
        for pair in [
            [a, b],
            [b, c],
            [c, a],
            [m, a],
            [b, m],
            [v, b],
            [m, v],
            [s2, s1],
            [v, a],
        ] {
            segments.push(Segment { face, nodes: pair });
        }
    }
    (nodes, segments)
}

/// How the paths of the two representations compare.
#[derive(Debug, Default, Clone)]
struct PathComparison {
    n_segments: usize,
    /// Segments whose crossed input edges are the same, in the same order.
    n_same_edges: usize,
    /// Segments whose pieces' input faces are the same.
    n_same_faces: usize,
    /// Segments whose paths lie within `1e-9` mean edge lengths of each other in space (the
    /// others differ by crossings of an input edge that a node lies on for one representation
    /// and beside for the other).
    n_equivalent: usize,
    max_node_difference: f64,
    max_crossing_difference: f64,
    max_relative_mismatch: f64,
}

fn compare_paths(shape: &Shape, a: &PathsOnInput, b: &PathsOnInput) -> PathComparison {
    let mut c = PathComparison {
        n_segments: a.paths.len(),
        max_relative_mismatch: b.max_relative_mismatch,
        ..PathComparison::default()
    };
    for (&p, &q) in a.nodes.iter().zip(&b.nodes) {
        c.max_node_difference = c.max_node_difference.max(distance(shape, p, q));
    }
    for (pa, pb) in a.paths.iter().zip(&b.paths) {
        let near = 1e-9 * shape.mean_edge_length();
        let edges = |p: &ddg_intrinsic::PathOnInput| -> Vec<u32> {
            let (first, last) = (p.points[0], p.points[p.points.len() - 1]);
            p.points[1..p.points.len() - 1]
                .iter()
                .filter(|&&x| distance(shape, x, first) > near && distance(shape, x, last) > near)
                .map(|x| match x {
                    SurfacePoint::Edge { edge, .. } => *edge,
                    _ => u32::MAX,
                })
                .collect()
        };
        let comparable =
            |p: &ddg_intrinsic::PathOnInput| -> Vec<SurfacePoint> {
                let (first, last) = (p.points[0], p.points[p.points.len() - 1]);
                let mut out = vec![first];
                out.extend(p.points[1..p.points.len() - 1].iter().filter(|&&x| {
                    distance(shape, x, first) > near && distance(shape, x, last) > near
                }));
                out.push(last);
                out
            };
        if edges(pa) != edges(pb) && std::env::var("I29_DIFF").is_ok() {
            println!(
                "    {}: signposts {:?}\n    integers  {:?}",
                shape.name, pa.points, pb.points
            );
        }
        if edges(pa) == edges(pb) {
            c.n_same_edges += 1;
            for (&x, &y) in comparable(pa).iter().zip(&comparable(pb)) {
                c.max_crossing_difference = c.max_crossing_difference.max(distance(shape, x, y));
            }
        }
        if pa.faces == pb.faces {
            c.n_same_faces += 1;
        }
        if path_distance(shape, &pa.points, &pb.points)
            .max(path_distance(shape, &pb.points, &pa.points))
            <= near
        {
            c.n_equivalent += 1;
        }
    }
    c
}

#[test]
fn segments_cross_the_input_edges_the_signposts_cross() {
    let shapes = shapes();
    let specs: Vec<(usize, Before)> = (0..shapes.len())
        .flat_map(|i| [(i, Before::Flips), (i, Before::Refined)])
        .collect();
    let results = par_map(&specs, |&(i, before)| {
        let shape = &shapes[i];
        let ((s, es), (it, ei)) = both(shape, before)?;
        let mut rng = SplitMix64::new(29 + i as u64);
        let (nodes, segments) = random_segments(&es, &mut rng, 60);
        let options = EdgeTraceOptions::default();
        let integers = it
            .segments_on_input(&nodes, &segments, &options)
            .map_err(|e| format!("integers: {e}"))?;
        // The signposts, segment by segment (a trace that does not arrive is theirs).
        let mut a = PathsOnInput {
            nodes: Vec::new(),
            paths: Vec::new(),
            max_relative_mismatch: 0.0,
        };
        let mut b = a.clone();
        let mut n_signpost_failures = 0;
        for (k, segment) in segments.iter().enumerate() {
            match s.segments_on_input(&nodes, std::slice::from_ref(segment), &options) {
                Ok(p) => {
                    if a.nodes.is_empty() {
                        a.nodes = p.nodes.clone();
                        b.nodes = integers.nodes.clone();
                    }
                    a.paths.push(p.paths[0].clone());
                    b.paths.push(integers.paths[k].clone());
                }
                Err(_) => n_signpost_failures += 1,
            }
        }
        b.max_relative_mismatch = integers.max_relative_mismatch;
        let same_mesh = es.mesh.faces() == ei.mesh.faces() && ei.n_off_midpoint_splits == 0;
        Ok::<_, String>((compare_paths(shape, &a, &b), n_signpost_failures, same_mesh))
    });
    println!(
        "{:<36} {:<8} {:>5} {:>5} {:>5} {:>5} {:>4} {:>10} {:>10} {:>10}",
        "shape",
        "before",
        "segs",
        "edges",
        "faces",
        "equiv",
        "fail",
        "max |Δn|",
        "max |Δx|",
        "mismatch"
    );
    let (mut total, mut same, mut equivalent) = (0, 0, 0);
    for (&(i, before), r) in specs.iter().zip(&results) {
        let shape = &shapes[i];
        let h = shape.mean_edge_length();
        match r {
            Ok((c, failures, same_mesh)) => {
                println!(
                    "{:<36} {:<8} {:>5} {:>5} {:>5} {:>5} {:>4} {:>10.2e} {:>10.2e} {:>10.2e}",
                    shape.name,
                    format!("{before:?}"),
                    c.n_segments,
                    c.n_same_edges,
                    c.n_same_faces,
                    c.n_equivalent,
                    failures,
                    c.max_node_difference / h,
                    c.max_crossing_difference / h,
                    c.max_relative_mismatch
                );
                assert!(c.max_relative_mismatch < 1e-9, "{}: mismatch", shape.name);
                if *same_mesh {
                    total += c.n_segments;
                    same += c.n_same_edges;
                    equivalent += c.n_equivalent;
                    assert!(c.max_node_difference < 1e-8 * h, "{}: nodes", shape.name);
                    // Crossings at a small angle are ill-conditioned in both representations
                    // (camelshotnoise: 2.4e-7 mean edge lengths).
                    assert!(
                        c.max_crossing_difference < 1e-6 * h,
                        "{}: crossings",
                        shape.name
                    );
                }
            }
            Err(e) => panic!("{} {before:?}: {e}", shape.name),
        }
    }
    println!(
        "{same} of {total} segments cross the same input edges, {equivalent} lie within 1e-9 \
         mean edge lengths of the signposts' path"
    );
    assert!(total > 0);
    assert_eq!(equivalent, total);
}

#[test]
fn bad_segments_are_the_signposts_errors() {
    let shape = &all_shapes()[0];
    let ((s, es), (it, _)) = both(shape, Before::Flips).unwrap();
    let options = EdgeTraceOptions::default();
    let e0 = es.triangulation_edge(0);
    let f_other = (0..es.mesh.n_faces())
        .find(|&f| (0..3).all(|k| es.mesh.edge(3 * f + k) != 0))
        .unwrap() as u32;
    let cases: Vec<(Vec<SurfacePoint>, Vec<Segment>)> = vec![
        // The end of an edge given as a point of the edge.
        (
            vec![SurfacePoint::Edge {
                edge: e0 as u32,
                t: 0.0,
            }],
            vec![],
        ),
        // A corner given as a point of the face.
        (
            vec![SurfacePoint::Face {
                face: 0,
                bary: [1.0, 0.0, 0.0],
            }],
            vec![],
        ),
        // A vertex that does not exist.
        (vec![SurfacePoint::Vertex(u32::MAX)], vec![]),
        // A node outside the face of its segment.
        (
            vec![
                SurfacePoint::Edge {
                    edge: e0 as u32,
                    t: 0.5,
                },
                SurfacePoint::Face {
                    face: f_other,
                    bary: [0.3, 0.3, 0.4],
                },
            ],
            vec![Segment {
                face: f_other,
                nodes: [0, 1],
            }],
        ),
        // A node index out of range, a face out of range.
        (
            vec![SurfacePoint::Face {
                face: 0,
                bary: [0.3, 0.3, 0.4],
            }],
            vec![Segment {
                face: 0,
                nodes: [0, 7],
            }],
        ),
        (
            vec![SurfacePoint::Face {
                face: 0,
                bary: [0.3, 0.3, 0.4],
            }],
            vec![Segment {
                face: u32::MAX,
                nodes: [0, 0],
            }],
        ),
        // A segment of length zero.
        (
            vec![SurfacePoint::Face {
                face: 0,
                bary: [0.3, 0.3, 0.4],
            }],
            vec![Segment {
                face: 0,
                nodes: [0, 0],
            }],
        ),
    ];
    for (nodes, segments) in cases {
        let a = s
            .segments_on_input(&nodes, &segments, &options)
            .unwrap_err();
        let b = it
            .segments_on_input(&nodes, &segments, &options)
            .unwrap_err();
        let b = ddg_integer_coords::to_intrinsic_error(b, "segments_on_input");
        assert_eq!(
            std::mem::discriminant(&a),
            std::mem::discriminant(&b),
            "{a} / {b}"
        );
        assert_eq!(a.to_string(), b.to_string());
    }
    let bad = EdgeTraceOptions {
        tolerance: f64::NAN,
        ..options
    };
    let a = s.segments_on_input(&[], &[], &bad).unwrap_err();
    let b = it.segments_on_input(&[], &[], &bad).unwrap_err();
    assert_eq!(
        a.to_string(),
        ddg_integer_coords::to_intrinsic_error(b, "").to_string()
    );
    assert!(matches!(a, IntrinsicError::InvalidOptions { .. }));
}

#[test]
fn the_export_is_deterministic() {
    let shape = &all_shapes()[1];
    let (_, (_, a)) = both(shape, Before::Refined).unwrap();
    let (_, (_, b)) = both(shape, Before::Refined).unwrap();
    assert_eq!(a, b);
}

#[test]
fn optimal_delaunay_on_integers_follows_the_signposts() {
    let shapes = shapes();
    let options = OptimalDelaunayOptions {
        max_insertions: ddg_intrinsic::InsertionLimit::Vertices(3000),
        max_iterations: 20,
        ..OptimalDelaunayOptions::default()
    };
    let results = par_map(&(0..shapes.len()).collect::<Vec<_>>(), |&i| {
        let shape = &shapes[i];
        let mut s = shape.triangulation();
        let mut it = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
        s.flip_to_delaunay(&DelaunayOptions::default()).unwrap();
        Triangulation::flip_to_delaunay(&mut it, &DelaunayOptions::default()).unwrap();
        let start = std::time::Instant::now();
        let rs = s
            .optimal_delaunay(&options)
            .map_err(|e| format!("signposts: {e}"))?;
        let ts = start.elapsed().as_secs_f64();
        let start = std::time::Instant::now();
        let ri = it
            .optimal_delaunay(&IntegerOdtOptions {
                odt: options,
                certify: true,
            })
            .map_err(|e| format!("integers: {e}"))?;
        let ti = start.elapsed().as_secs_f64();
        let not_moved = |r: &ddg_intrinsic::OptimalDelaunayReport| -> usize {
            r.iterations.iter().map(|it| it.n_not_moved).sum()
        };
        let halved = |r: &ddg_intrinsic::OptimalDelaunayReport| -> usize {
            r.iterations.iter().map(|it| it.n_halved).sum()
        };
        // The exports of both, then the places of the vertices.
        let es = s
            .export_simplicial(&ExportOptions::default())
            .map_err(|e| e.to_string())?;
        let ei = it
            .export_simplicial(&ExportOptions::default())
            .map_err(|e| e.to_string())?;
        let same = es.mesh.n_vertices() == ei.mesh.n_vertices();
        let max_place = if same {
            es.locations
                .iter()
                .zip(&ei.locations)
                .map(|(&p, &q)| distance(shape, p, q))
                .fold(0.0, f64::max)
        } else {
            f64::NAN
        };
        Ok::<_, String>((
            rs.clone(),
            ri.odt.clone(),
            (not_moved(&rs), not_moved(&ri.odt)),
            (halved(&rs), halved(&ri.odt)),
            max_place,
            (ts, ti),
            ri.certificate.is_some(),
        ))
    });
    println!(
        "{:<36} {:>5} {:>5} {:>5} {:>5} {:>7} {:>7} {:>7} {:>7} {:>9} {:>7} {:>7}",
        "shape",
        "it s",
        "it i",
        "ins s",
        "ins i",
        "stay s",
        "stay i",
        "half s",
        "half i",
        "max |Δq|",
        "s (s)",
        "s (i)"
    );
    for (i, r) in results.iter().enumerate() {
        let shape = &shapes[i];
        match r {
            Ok((rs, ri, nm, hv, place, (ts, ti), certified)) => {
                println!(
                    "{:<36} {:>5} {:>5} {:>5} {:>5} {:>7} {:>7} {:>7} {:>7} {:>9.2e} {:>7.3} {:>7.3}",
                    shape.name,
                    rs.n_iterations,
                    ri.n_iterations,
                    rs.n_inserted,
                    ri.n_inserted,
                    nm.0,
                    nm.1,
                    hv.0,
                    hv.1,
                    place / shape.mean_edge_length(),
                    ts,
                    ti
                );
                assert!(certified);
                println!(
                    "{:<36} min angle {:.2} / {:.2} deg, area cv {:.3} / {:.3}",
                    "",
                    rs.statistics.angles.min.to_degrees(),
                    ri.statistics.angles.min.to_degrees(),
                    rs.statistics.area_cv,
                    ri.statistics.area_cv
                );
            }
            Err(e) => println!("{:<36} error: {e}", shape.name),
        }
    }
}

#[test]
fn a_refused_move_leaves_the_state_unchanged() {
    // Every inserted vertex of a refined state, moved by a step that ends inside the star:
    // either it moves and the result certifies, or the state is exactly as before.
    let shape = &all_shapes()[1];
    let mut it = IntegerTriangulation::new(&shape.mesh, &shape.geometry).unwrap();
    it.refine(&IntegerRefineOptions {
        refine: refine_options(),
        rules: RefineRules::Ruppert,
        ..Default::default()
    })
    .unwrap();
    let mut rng = SplitMix64::new(7);
    let (mut moved, mut refused) = (0, 0);
    for v in it.n_input_vertices()..it.complex().n_vertices() {
        if it.complex().is_boundary_vertex(v) {
            continue;
        }
        let Ok(goal) = it.odt_target(v, true) else {
            continue;
        };
        let scale = rng.uniform(0.1, 0.6);
        let before = it.clone();
        match ddg_intrinsic::OdtMoves::move_vertex_with_shape(
            &mut it,
            v,
            [scale * goal[0], scale * goal[1]],
            1e-2,
        ) {
            Ok(_) => {
                moved += 1;
                it.certify()
                    .unwrap_or_else(|e| panic!("vertex {v}: the moved star does not certify: {e}"));
            }
            Err(_) => {
                refused += 1;
                assert_eq!(it, before, "vertex {v}: a refused move changed the state");
            }
        }
    }
    it.certify().unwrap();
    println!("{moved} moved, {refused} refused");
    assert!(moved > 0);
}
