//! Input-aligned tangent frames on integer coordinates (SPEC §4.12, I26): on the same
//! triangulation they are the signposts' `signpost(h)` to rounding, after conversions of signpost
//! states, after the same operations run on both representations (Delaunay flipping, insertions,
//! edge splits, refinement, the optimal Delaunay iteration); a flip changes the directions of
//! the flipped edge only, bitwise; bad integers are errors.
//!
//! Run with `--nocapture` for the tables.

mod common;

use std::f64::consts::{PI, TAU};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use common::{SCENARIOS, SEEDS, Shape, all_shapes, convert, inner_bary, models, state};
use ddg_core::generators::SplitMix64;
use ddg_integer_coords::{IntegerOdtOptions, IntegerRefineOptions, IntegerTriangulation};
use ddg_intrinsic::{
    DelaunayOptions, IntrinsicTriangulation, OptimalDelaunayOptions, RefineOptions, RefineRules,
    SurfacePoint, Triangulation,
};

/// At most this many threads (the machine is shared).
const THREADS: usize = 16;

fn par_map<S: Sync, R: Send>(specs: &[S], f: impl Fn(&S) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(THREADS);
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

/// The difference of two directions of a vertex: modulo `2π` inside, as they are on the
/// boundary (where directions do not wrap).
fn direction_difference(a: f64, b: f64, boundary: bool) -> f64 {
    let d = (a - b).abs();
    if boundary { d } else { d.min(TAU - d) }
}

/// How the integers' frames compare with the signposts' on one triangulation.
#[derive(Debug, Clone, Copy, Default)]
struct Comparison {
    n_halfedges: usize,
    /// The largest difference, in radians of the normalized frame.
    max_difference: f64,
    /// Where: a halfedge and whether its tail is an input vertex.
    worst: Option<(usize, bool)>,
    /// The largest difference at input vertices and at inserted vertices.
    max_input: f64,
    max_inserted: f64,
    n_clamped: usize,
    max_mismatch: f64,
}

fn compare(s: &IntrinsicTriangulation, i: &IntegerTriangulation) -> Comparison {
    let (dirs, report) = i.frame_directions().unwrap();
    assert_eq!(dirs.len(), s.n_halfedges());
    let mut c = Comparison {
        n_halfedges: dirs.len(),
        n_clamped: report.n_clamped,
        max_mismatch: report.max_relative_mismatch,
        ..Default::default()
    };
    assert_eq!(report.n_along + report.n_from_segment, dirs.len());
    for (h, &d) in dirs.iter().enumerate() {
        let v = s.from_vertex(h);
        let boundary = s.is_boundary_vertex(v);
        assert!(
            if boundary {
                (0.0..=PI).contains(&d)
            } else {
                (0.0..TAU).contains(&d)
            },
            "halfedge {h}: direction {d} out of range"
        );
        let x = direction_difference(d, s.signpost(h), boundary);
        let input = v < s.n_input_vertices();
        if input {
            c.max_input = c.max_input.max(x);
        } else {
            c.max_inserted = c.max_inserted.max(x);
        }
        if x > c.max_difference {
            c.max_difference = x;
            c.worst = Some((h, input));
        }
    }
    c
}

/// Whether two triangulations are the same: connectivity and the kind and element of every
/// place (places themselves may differ by rounding).
fn same_triangulation(s: &IntrinsicTriangulation, i: &IntegerTriangulation) -> bool {
    if s.complex().faces() != i.complex().faces() || s.n_vertices() != i.complex().n_vertices() {
        return false;
    }
    (0..s.n_vertices()).all(|v| match (s.vertex_location(v), i.vertex_location(v)) {
        (SurfacePoint::Vertex(a), SurfacePoint::Vertex(b)) => a == b,
        (SurfacePoint::Edge { edge: a, .. }, SurfacePoint::Edge { edge: b, .. }) => a == b,
        (SurfacePoint::Face { face: a, .. }, SurfacePoint::Face { face: b, .. }) => a == b,
        _ => false,
    })
}

/// The largest difference of the places of the vertices, in their barycentric coordinates (the
/// triangulations are the same, [`same_triangulation`]).
fn max_place_difference(s: &IntrinsicTriangulation, i: &IntegerTriangulation) -> f64 {
    (0..s.n_vertices())
        .map(|v| match (s.vertex_location(v), i.vertex_location(v)) {
            (SurfacePoint::Edge { t: a, .. }, SurfacePoint::Edge { t: b, .. }) => (a - b).abs(),
            (SurfacePoint::Face { bary: a, .. }, SurfacePoint::Face { bary: b, .. }) => {
                (0..3).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max)
            }
            _ => 0.0,
        })
        .fold(0.0, f64::max)
}

const TOLERANCE: f64 = 1e-8;

#[test]
fn frames_equal_the_signposts_on_converted_states() {
    let shapes = all_shapes();
    let mut specs = Vec::new();
    for s in 0..shapes.len() {
        for scenario in SCENARIOS {
            for seed in SEEDS {
                specs.push((s, scenario, seed));
            }
        }
    }
    let results = par_map(&specs, |&(s, scenario, seed)| {
        let shape = &shapes[s];
        let signposts = state(shape, scenario, seed);
        let Ok((integers, _)) = convert(&signposts) else {
            return None;
        };
        Some(compare(&signposts, &integers))
    });
    let mut total = Comparison::default();
    let mut n_states = 0;
    let mut n_unconverted = 0;
    println!(
        "{:<34} {:<26} {:>6} {:>8} {:>10} {:>10} {:>5} {:>9}",
        "shape", "scenario", "seed", "halfedges", "max V0", "max V*", "clamp", "mismatch"
    );
    for (&(s, scenario, seed), r) in specs.iter().zip(&results) {
        let Some(c) = r else {
            n_unconverted += 1;
            continue;
        };
        n_states += 1;
        println!(
            "{:<34} {:<26} {:>6} {:>8} {:>10.2e} {:>10.2e} {:>5} {:>9.1e}",
            shapes[s].name,
            format!("{scenario:?}"),
            seed,
            c.n_halfedges,
            c.max_input,
            c.max_inserted,
            c.n_clamped,
            c.max_mismatch
        );
        assert!(
            c.max_difference < TOLERANCE,
            "{} {scenario:?} {seed}: directions differ by {} at {:?}",
            shapes[s].name,
            c.max_difference,
            c.worst
        );
        total.n_halfedges += c.n_halfedges;
        total.max_input = total.max_input.max(c.max_input);
        total.max_inserted = total.max_inserted.max(c.max_inserted);
        total.n_clamped += c.n_clamped;
        total.max_mismatch = total.max_mismatch.max(c.max_mismatch);
    }
    println!(
        "{n_states} states ({n_unconverted} not converted), {} halfedges: max difference {:.2e} at \
         input vertices, {:.2e} at inserted vertices; {} clamps; relative length mismatch {:.1e}",
        total.n_halfedges, total.max_input, total.max_inserted, total.n_clamped, total.max_mismatch
    );
    assert!(n_states > 0);
}

/// The same operations on both representations from the input mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Native {
    Delaunay,
    /// Random insertions in faces and edge splits, the same on both, then Delaunay flipping.
    Insertions,
    /// Delaunay refinement to 25° with the signposts' rules.
    Refined,
    /// The optimal Delaunay iteration.
    OptimalDelaunay,
}

fn native(
    shape: &Shape,
    op: Native,
) -> Result<(IntrinsicTriangulation, IntegerTriangulation), String> {
    let mut s = shape.triangulation();
    let mut i =
        IntegerTriangulation::new(&shape.mesh, &shape.geometry).map_err(|e| e.to_string())?;
    let delaunay = DelaunayOptions::default();
    s.flip_to_delaunay(&delaunay).map_err(|e| e.to_string())?;
    Triangulation::flip_to_delaunay(&mut i, &delaunay).map_err(|e| e.to_string())?;
    match op {
        Native::Delaunay => {}
        Native::Insertions => {
            let mut rng = SplitMix64::new(7);
            let n = (s.n_faces() / 4).clamp(3, 200);
            for k in 0..n {
                if k % 3 == 2 {
                    let e = (rng.next_u64() % s.n_edges() as u64) as usize;
                    let t = rng.uniform(0.2, 0.8);
                    let a = Triangulation::split_edge(&mut s, e, t);
                    let b = Triangulation::split_edge(&mut i, e, t);
                    if a.is_ok() != b.is_ok() {
                        return Err(format!("edge split {e} at {t}: {a:?} vs {b:?}"));
                    }
                } else {
                    let f = (rng.next_u64() % s.n_faces() as u64) as usize;
                    let bary = inner_bary(&mut rng);
                    let a = Triangulation::insert_vertex(&mut s, f, bary);
                    let b = Triangulation::insert_vertex(&mut i, f, bary);
                    if a.is_ok() != b.is_ok() {
                        return Err(format!("insertion in face {f}: {a:?} vs {b:?}"));
                    }
                }
            }
            s.flip_to_delaunay(&delaunay).map_err(|e| e.to_string())?;
            Triangulation::flip_to_delaunay(&mut i, &delaunay).map_err(|e| e.to_string())?;
        }
        Native::Refined => {
            let options = RefineOptions {
                max_insertions: ddg_intrinsic::InsertionLimit::Vertices(4000),
                ..RefineOptions::default()
            };
            s.refine(&options).map_err(|e| e.to_string())?;
            i.refine(&IntegerRefineOptions {
                refine: options,
                rules: RefineRules::Ruppert,
                ..Default::default()
            })
            .map_err(|e| e.to_string())?;
        }
        Native::OptimalDelaunay => {
            let options = OptimalDelaunayOptions {
                max_insertions: ddg_intrinsic::InsertionLimit::Vertices(2000),
                max_iterations: 10,
                ..OptimalDelaunayOptions::default()
            };
            s.optimal_delaunay(&options).map_err(|e| e.to_string())?;
            i.optimal_delaunay(&IntegerOdtOptions {
                odt: options,
                certify: true,
            })
            .map_err(|e| e.to_string())?;
        }
    }
    Ok((s, i))
}

#[test]
fn frames_equal_the_signposts_after_the_same_operations() {
    let mut shapes = all_shapes();
    shapes.extend(models());
    let ops = [
        Native::Delaunay,
        Native::Insertions,
        Native::Refined,
        Native::OptimalDelaunay,
    ];
    let specs: Vec<(usize, Native)> = (0..shapes.len())
        .flat_map(|s| ops.iter().map(move |&op| (s, op)))
        .collect();
    let results = par_map(&specs, |&(s, op)| {
        let (a, b) = native(&shapes[s], op)?;
        if !same_triangulation(&a, &b) {
            return Ok(None);
        }
        Ok::<_, String>(Some((compare(&a, &b), max_place_difference(&a, &b))))
    });
    let (mut n_compared, mut n_different, mut max_input, mut max_inserted) = (0, 0, 0.0f64, 0.0f64);
    let mut n_halfedges = 0;
    for (&(s, op), r) in specs.iter().zip(&results) {
        match r {
            Ok(Some((c, place))) => {
                println!(
                    "{:<34} {:<16} {:>8} halfedges: max {:.2e} at V0, {:.2e} at V*, {} clamps; \
                     places within {place:.1e}",
                    shapes[s].name,
                    format!("{op:?}"),
                    c.n_halfedges,
                    c.max_input,
                    c.max_inserted,
                    c.n_clamped
                );
                assert!(
                    c.max_difference < TOLERANCE,
                    "{} {op:?}: directions differ by {} at {:?}",
                    shapes[s].name,
                    c.max_difference,
                    c.worst
                );
                n_compared += 1;
                n_halfedges += c.n_halfedges;
                max_input = max_input.max(c.max_input);
                max_inserted = max_inserted.max(c.max_inserted);
            }
            Ok(None) => {
                println!(
                    "{:<34} {:<16} not the same triangulation",
                    shapes[s].name,
                    format!("{op:?}")
                );
                n_different += 1;
            }
            Err(e) => panic!("{} {op:?}: {e}", shapes[s].name),
        }
    }
    println!(
        "{n_compared} runs compared ({n_different} ended on different triangulations), \
         {n_halfedges} halfedges: max difference {max_input:.2e} at input vertices, \
         {max_inserted:.2e} at inserted vertices"
    );
    assert!(n_compared > specs.len() / 2);
}

/// The directions at both ends of every edge, by edge id (halfedge ids change under flips).
fn by_edge(i: &IntegerTriangulation) -> Vec<[Option<(usize, f64)>; 2]> {
    let (dirs, _) = i.frame_directions().unwrap();
    let c = i.complex();
    (0..c.n_edges())
        .map(|e| {
            let h = c.edge_halfedge(e);
            [
                Some((c.from_vertex(h), dirs[h])),
                c.twin(h).map(|t| (c.from_vertex(t), dirs[t])),
            ]
        })
        .collect()
}

#[test]
fn a_flip_changes_the_frames_of_the_flipped_edge_only() {
    let shapes = all_shapes();
    let specs: Vec<usize> = (0..shapes.len()).collect();
    let counts = par_map(&specs, |&s| {
        let shape = &shapes[s];
        let mut t = state(shape, common::Scenario::Mixed, 5);
        let Ok((mut i, _)) = convert(&t) else {
            return (0, 0);
        };
        let mut rng = SplitMix64::new(3);
        let mut before = by_edge(&i);
        let (mut flipped, mut compared) = (0, 0);
        for _ in 0..30 {
            let e = (rng.next_u64() % i.complex().n_edges() as u64) as usize;
            if Triangulation::flip_edge(&mut i, e).is_err() {
                continue;
            }
            t.flip_edge(e).unwrap();
            flipped += 1;
            let after = by_edge(&i);
            for (k, (x, y)) in before.iter().zip(&after).enumerate() {
                if k != e {
                    assert_eq!(x, y, "{}: edge {k} after flipping edge {e}", shape.name);
                    compared += 1;
                }
            }
            before = after;
        }
        // And the signposts still agree after the flips.
        let c = compare(&t, &i);
        assert!(c.max_difference < TOLERANCE, "{}: {c:?}", shape.name);
        (flipped, compared)
    });
    let (flips, edges) = counts.iter().fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
    println!("{flips} flips, {edges} edges unchanged bitwise");
    assert!(flips > 0);
}

#[test]
fn corrupt_integers_are_errors() {
    let shape = &all_shapes()[0];
    let t = state(shape, common::Scenario::Mixed, 11);
    let (i, _) = convert(&t).unwrap();
    let mut rng = SplitMix64::new(1);
    let mut n_errors = 0;
    for k in 0..200 {
        let mut parts = i.clone().into_parts();
        if k % 2 == 0 {
            let e = (rng.next_u64() % parts.normal.len() as u64) as usize;
            parts.normal[e] += if k % 4 == 0 { 1 } else { -1 };
        } else {
            let h = (rng.next_u64() % parts.roundabouts.len() as u64) as usize;
            if parts.roundabouts[h] == ddg_integer_coords::NO_ROUNDABOUT {
                continue;
            }
            parts.roundabouts[h] = parts.roundabouts[h].wrapping_add(1 + k as u32 % 3);
        }
        let Ok(bad) = IntegerTriangulation::from_parts(parts) else {
            continue;
        };
        // Never a panic; an error names what is wrong.
        if let Err(e) = bad.frame_directions() {
            assert!(!e.to_string().is_empty());
            n_errors += 1;
        }
    }
    println!("{n_errors} of 200 mutations are errors");
    assert!(n_errors > 0);
}
