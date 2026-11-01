//! The ICE metric and Alg. 1 (SPEC §4, §5): the transport against the planar centre of mass,
//! the key check of S5, targets, fixed vertices, snapshots, the report, bad options.

mod common;

use std::f64::consts::PI;

use common::{Shape, crane_models, planar_grid, small_shapes};
use ddg_intrinsic_simplification::{
    Channel, ChannelSource, ChannelState, ChannelWeights, CoarseMesh, Frames, IemError,
    MetricOptions, PointLocation, Requeue, SimplifyOptions, StopReason, Target, metric, simplify,
};

fn area_only() -> ChannelWeights {
    ChannelWeights {
        channels: vec![Channel {
            source: ChannelSource::Area,
            weight: 1.0,
        }],
    }
}

/// The planar direction of a halfedge of a coarse mesh drawn on a planar input (the coarse
/// metric stays planar: every edge is the straight segment between its fine end points).
fn planar_angle(shape: &Shape, mesh: &CoarseMesh, h: usize) -> f64 {
    let c = mesh.complex();
    let a = shape.positions[mesh.fine_vertex(c.from_vertex(h)) as usize];
    let b = shape.positions[mesh.fine_vertex(c.to_vertex(h)) as usize];
    (b[1] - a[1]).atan2(b[0] - a[0])
}

/// On a planar grid the normalized frames are rigid (`Θ = 2π` inside, `π` on the straight
/// sides, the scale is 1; the four fixed corners are not checked), so Eq. 7/8 is exact: the error vector of every vertex points to the centre of
/// mass of the fine area it has collected (`iem-center-of-mass-update`). The mass bookkeeping is
/// replayed in the plane, independently of frames, and compared after every removal.
#[test]
fn error_vectors_point_to_the_planar_centre_of_mass() {
    for (n, jitter, seed) in [(8usize, 0.3, 4u64), (16, 0.5, 5)] {
        let shape = planar_grid(n, jitter, seed);
        let mut mesh = shape.coarse(&[]);
        for corner in [0, n, n * (n + 1), (n + 1) * (n + 1) - 1] {
            mesh.set_fixed(corner as u32, true).unwrap();
        }
        mesh.flip_to_delaunay().unwrap();
        let mut state = ChannelState::new(&mesh, &area_only(), false).unwrap();
        let nv = shape.positions.len();
        let mut mass: Vec<f64> = (0..nv as u32).map(|v| state.mass(v, 0)).collect();
        let mut centre: Vec<[f64; 2]> = shape.positions.iter().map(|p| [p[0], p[1]]).collect();
        let options = MetricOptions::default();
        let mut worst = 0.0f64;
        let mut removed = [0usize; 2];
        loop {
            let mut progress = false;
            for v in 0..nv as u32 {
                let Ok(p) = mesh.begin_removal(v) else {
                    continue;
                };
                state.apply_removal(&p, &options, Frames::AfterFlattening);
                let share = 1.0 / p.neighbours.len() as f64;
                assert!(p.curvature.abs() <= 1e-12, "a planar vertex is flat");
                let (mi, ci) = (mass[v as usize], centre[v as usize]);
                for nb in &p.neighbours {
                    let j = nb.vertex as usize;
                    let m = mass[j] + share * mi;
                    centre[j] = [0, 1].map(|k| (share * mi * ci[k] + mass[j] * centre[j][k]) / m);
                    mass[j] = m;
                }
                mesh.finish_removal().unwrap();
                removed
                    [(p.kind != ddg_intrinsic_simplification::RemovalKind::Interior) as usize] += 1;
                progress = true;
                // Every coarse vertex: its vector in the plane against its centre of mass.
                for cv in 0..mesh.n_vertices() {
                    let f = mesh.fine_vertex(cv);
                    if mesh.is_fixed(f) {
                        // A corner (Θ = π/2): its normalized frame is not the planar one.
                        continue;
                    }
                    let r = mesh.reference_halfedge(cv);
                    let offset = planar_angle(&shape, &mesh, r) - mesh.angular_coordinates()[r];
                    let t = state.vector(f, 0);
                    let (s, c) = offset.sin_cos();
                    let planar = [c * t[0] - s * t[1], s * t[0] + c * t[1]];
                    let x = shape.positions[f as usize];
                    let want = [centre[f as usize][0] - x[0], centre[f as usize][1] - x[1]];
                    let err = (planar[0] - want[0]).hypot(planar[1] - want[1]);
                    worst = worst.max(err);
                    assert!(
                        (state.mass(f, 0) - mass[f as usize]).abs() <= 1e-14,
                        "mass of {f}"
                    );
                }
            }
            if !progress {
                break;
            }
        }
        eprintln!(
            "grid {n}: removed {removed:?} (interior, boundary), {} left, worst error-vector \
             error {worst:e}",
            mesh.n_vertices()
        );
        assert!(removed[1] > 0 && removed[0] > 0);
        assert!(worst < 1e-12, "{worst}");
        assert!((state.total_mass(0) - 1.0).abs() < 1e-12);
    }
}

/// Eq. 9 and the memoryless cost on a prepared removal: the weights are the normalized
/// curvature changes, and with `t = 0` the memory-based cost equals `Σ α m_i ℓ_ij`.
#[test]
fn memoryless_cost_is_the_first_memory_cost() {
    for shape in small_shapes() {
        let mut mesh = shape.coarse(&[]);
        let state = ChannelState::new(&mesh, &ChannelWeights::mixture(0.3), false).unwrap();
        let memory = MetricOptions::default();
        let memoryless = MetricOptions {
            memory: false,
            ..MetricOptions::default()
        };
        let mut checked = 0;
        for v in 0..mesh.n_fine_vertices() as u32 {
            let Ok(p) = mesh.tentative_removal(v) else {
                continue;
            };
            let a = state.removal_cost(&p, &memory, Frames::AfterFlattening);
            let b = state.removal_cost(&p, &memoryless, Frames::AfterFlattening);
            assert!(
                (a - b).abs() <= 1e-12 * a.abs().max(1e-300),
                "{}: {a} {b}",
                shape.name
            );
            let before: Vec<f64> = p.neighbours.iter().map(|n| n.curvature_before).collect();
            let after: Vec<f64> = p.neighbours.iter().map(|n| n.curvature_after).collect();
            let alpha = metric::redistribution_weights(p.curvature, &before, &after);
            assert!((alpha.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            let lengths: Vec<f64> = p.neighbours.iter().map(|n| n.length).collect();
            let direct: f64 = (0..state.n_channels())
                .map(|c| {
                    state.weight(c) * metric::memoryless_cost(state.mass(v, c), &alpha, &lengths)
                })
                .sum();
            assert!((direct - b).abs() <= 1e-12 * b.abs().max(1e-300));
            checked += 1;
        }
        assert!(checked > 0, "{}", shape.name);
    }
}

fn run(
    shape: &Shape,
    options: &SimplifyOptions,
) -> (CoarseMesh, ddg_intrinsic_simplification::SimplifyReport) {
    let mut mesh = shape.coarse(&shape.centroids());
    let report = simplify(&mut mesh, options).unwrap_or_else(|e| panic!("{}: {e}", shape.name));
    (mesh, report)
}

/// Alg. 1 on every test shape to the minimum: the certificate with Delaunay, the S5 key check
/// (a key computed after the last change of the vertex and its neighbours is the commit cost,
/// bitwise), mass conservation, counts that add up.
#[test]
fn simplify_to_the_minimum_keeps_the_invariants() {
    let mut shapes = small_shapes();
    shapes.extend(crane_models());
    for shape in &shapes {
        for options in [
            SimplifyOptions {
                record_sequence: true,
                ..SimplifyOptions::default()
            },
            SimplifyOptions {
                channels: ChannelWeights::mixture(0.5),
                requeue: Requeue::OneRing,
                record_sequence: true,
                ..SimplifyOptions::default()
            },
            SimplifyOptions {
                metric: MetricOptions {
                    memory: false,
                    ..MetricOptions::default()
                },
                ..SimplifyOptions::default()
            },
        ] {
            let mut mesh = shape.coarse(&shape.centroids());
            let state0 = ChannelState::new(&mesh, &options.channels, false).unwrap();
            let r = simplify(&mut mesh, &options).unwrap();
            let removed = r.n_interior + r.n_boundary + r.n_ear;
            assert_eq!(r.n_vertices_before - removed, r.n_vertices_after);
            assert_eq!(r.n_vertices_after, mesh.n_vertices());
            assert_eq!(r.stop, StopReason::NothingRemovable, "{}", shape.name);
            assert!(
                r.commit_failures.is_empty(),
                "{}: {:?}",
                shape.name,
                r.commit_failures
            );
            assert_eq!(r.keys.n_clean_mismatch, 0, "{}: {:?}", shape.name, r.keys);
            assert!(r.keys.n_clean > 0, "{}", shape.name);
            let cert = mesh
                .certify(true)
                .unwrap_or_else(|e| panic!("{}: {e}", shape.name));
            assert_eq!(cert.n_not_delaunay, 0);
            for c in 0..r.channels.n_channels() {
                let (a, b) = (state0.total_mass(c), r.channels.total_mass(c));
                assert!(
                    (a - b).abs() <= 1e-9 * (1.0 + a.abs()),
                    "{}: {a} {b}",
                    shape.name
                );
            }
            // Removed vertices hold no mass; everything is finite.
            for v in 0..mesh.n_fine_vertices() as u32 {
                for c in 0..r.channels.n_channels() {
                    let m = r.channels.mass(v, c);
                    assert!(m.is_finite() && m >= 0.0);
                    if !mesh.is_vertex(v) {
                        assert_eq!(m, 0.0);
                    }
                }
            }
            if options.record_sequence {
                assert_eq!(r.sequence.len(), removed);
                assert!(r.sequence.iter().all(|s| s.cost.is_finite()));
            }
            eprintln!(
                "{}: {} -> {} ({:?}), pops {} stale {}, evaluations {}, keys {:?}, skips {:?}",
                shape.name,
                r.n_vertices_before,
                r.n_vertices_after,
                r.stop,
                r.n_pops,
                r.n_stale_pops,
                r.n_evaluations,
                r.keys,
                r.skips
            );
        }
    }
}

/// The same input gives the same removals (deterministic ties).
#[test]
fn simplify_is_deterministic() {
    for shape in crane_models().iter().take(2) {
        let options = SimplifyOptions {
            target: Target::Ratio(0.2),
            record_sequence: true,
            ..SimplifyOptions::default()
        };
        let (a, ra) = run(shape, &options);
        let (b, rb) = run(shape, &options);
        assert_eq!(ra.sequence, rb.sequence);
        assert!(a.bitwise_eq(&b));
    }
}

/// The targets: a vertex count, a ratio, per component, a maximum cost; fixed vertices stay.
#[test]
fn targets_and_fixed_vertices() {
    let Some(shape) = crane_models().into_iter().next() else {
        return;
    };
    let n = shape.positions.len();
    let (mesh, r) = run(
        &shape,
        &SimplifyOptions {
            target: Target::Vertices(300),
            ..SimplifyOptions::default()
        },
    );
    assert_eq!(mesh.n_vertices(), 300);
    assert_eq!(r.stop, StopReason::TargetReached);
    let (mesh, _) = run(
        &shape,
        &SimplifyOptions {
            target: Target::Ratio(0.25),
            ..SimplifyOptions::default()
        },
    );
    assert_eq!(mesh.n_vertices(), (0.25 * n as f64).round() as usize);
    let (_, all) = run(
        &shape,
        &SimplifyOptions {
            target: Target::Vertices(0),
            record_sequence: true,
            ..SimplifyOptions::default()
        },
    );
    let median = all.sequence[all.sequence.len() / 2].key;
    let (_, r) = run(
        &shape,
        &SimplifyOptions {
            target: Target::MaxCost(median),
            record_sequence: true,
            ..SimplifyOptions::default()
        },
    );
    assert_eq!(r.stop, StopReason::MaxCost);
    assert!(r.sequence.iter().all(|s| s.key <= median));
    let fixed: Vec<u32> = (0..n as u32).step_by(7).collect();
    let (mesh, r) = run(
        &shape,
        &SimplifyOptions {
            fixed: fixed.clone(),
            ..SimplifyOptions::default()
        },
    );
    assert!(fixed.iter().all(|&v| mesh.is_vertex(v)));
    assert!(r.skips.get("fixed").copied().unwrap_or(0) >= fixed.len());
}

/// Two components: each is coarsened to its own ratio.
#[test]
fn per_component_target() {
    let shapes = small_shapes();
    let (a, b) = (&shapes[0], &shapes[2]);
    let offset = a.positions.len() as u32;
    let mut positions = a.positions.clone();
    positions.extend(b.positions.iter().map(|p| [p[0] + 5.0, p[1], p[2]]));
    let mut faces = a.faces.clone();
    faces.extend(b.faces.iter().map(|f| f.map(|v| v + offset)));
    let shape = Shape::new("two components", positions, faces);
    let (_, r) = run(
        &shape,
        &SimplifyOptions {
            target: Target::PerComponent(0.5),
            ..SimplifyOptions::default()
        },
    );
    assert_eq!(r.components.len(), 2);
    for c in &r.components {
        assert_eq!(
            c.after,
            (0.5 * c.before as f64).round() as usize,
            "{:?}",
            r.components
        );
    }
    assert_eq!(r.stop, StopReason::TargetReached);
}

/// Snapshots at given counts: compact meshes whose `previous` locate the previous level's
/// vertices (a vertex kept is a vertex; the others lie in faces with valid coordinates).
#[test]
fn snapshots_form_a_hierarchy() {
    let Some(shape) = crane_models().into_iter().next() else {
        return;
    };
    let counts = vec![2000, 800, 200];
    let (mesh, r) = run(
        &shape,
        &SimplifyOptions {
            target: Target::Vertices(100),
            snapshots: counts.clone(),
            ..SimplifyOptions::default()
        },
    );
    assert_eq!(r.snapshots.len(), 3);
    let mut previous: Vec<u32> = (0..shape.positions.len() as u32).collect();
    for (s, &n) in r.snapshots.iter().zip(&counts) {
        assert_eq!(s.vertices.len(), n);
        assert_eq!(s.previous_vertices, previous);
        assert_eq!(s.previous.len(), previous.len());
        let kept: std::collections::HashSet<u32> = s.vertices.iter().copied().collect();
        for (&v, loc) in previous.iter().zip(&s.previous) {
            match *loc {
                PointLocation::Vertex(w) => {
                    assert_eq!(w, v);
                    assert!(kept.contains(&v));
                }
                PointLocation::Face {
                    face,
                    corners,
                    bary,
                } => {
                    assert!(!kept.contains(&v));
                    let f = s.faces[face as usize];
                    assert_eq!(f.map(|k| s.vertices[k as usize]), corners);
                    assert!(bary.iter().all(|&b| b >= 0.0));
                    assert!((bary.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                }
            }
        }
        for (f, l) in s.faces.iter().zip(&s.lengths) {
            assert!(f.iter().all(|&k| (k as usize) < s.vertices.len()));
            assert!(l.iter().all(|&x| x > 0.0));
        }
        let mut next = s.vertices.clone();
        next.sort_unstable();
        previous = next;
    }
    assert_eq!(mesh.n_vertices(), 100);
}

/// The reference's option set runs, gives keys and commit costs that differ only by the
/// reference's frames (R16), and keeps the invariants.
#[test]
fn reference_compat_runs() {
    for shape in crane_models() {
        let options = SimplifyOptions {
            record_sequence: true,
            ..SimplifyOptions::reference_compat(0.0, Target::Ratio(0.1))
        };
        let (mesh, r) = run(&shape, &options);
        assert_eq!(r.n_initial_delaunay_flips, 0);
        assert!(r.initial_delaunay_status.is_none());
        mesh.certify(false).unwrap();
        assert_eq!(r.sequence.len(), r.n_interior + r.n_boundary + r.n_ear);
        // Delaunay is restored locally, as in the reference (2-ring seeding, tolerance 1e-5).
        assert_eq!(
            mesh.options().delaunay,
            ddg_intrinsic_simplification::CoarseMeshOptions::default().delaunay
        );
    }
}

/// Bad options are errors naming the problem; the mesh is unchanged.
#[test]
fn bad_options_are_errors() {
    let shape = planar_grid(4, 0.0, 1);
    let mut mesh = shape.coarse(&[]);
    let before = mesh.clone();
    let cases = [
        SimplifyOptions {
            target: Target::Ratio(1.5),
            ..SimplifyOptions::default()
        },
        SimplifyOptions {
            target: Target::MaxCost(f64::NAN),
            ..SimplifyOptions::default()
        },
        SimplifyOptions {
            snapshots: vec![10, 12],
            ..SimplifyOptions::default()
        },
        SimplifyOptions {
            fixed: vec![1000],
            ..SimplifyOptions::default()
        },
        SimplifyOptions {
            channels: ChannelWeights::mixture(1.5),
            ..SimplifyOptions::default()
        },
        SimplifyOptions {
            channels: ChannelWeights::curvature().with(ChannelSource::PerVertex(vec![1.0; 3]), 1.0),
            ..SimplifyOptions::default()
        },
    ];
    for o in &cases {
        match simplify(&mut mesh, o) {
            Err(IemError::InvalidOptions { what }) => eprintln!("{what}"),
            other => panic!("{other:?}"),
        }
        assert!(mesh.bitwise_eq(&before));
    }
    let _ = PI;
}

/// A signed per-vertex array becomes two channels; the curvature channels follow the sign.
#[test]
fn per_vertex_channels_split_by_sign() {
    let shape = planar_grid(4, 0.2, 3);
    let mesh = shape.coarse(&[]);
    let n = mesh.n_fine_vertices();
    let values: Vec<f64> = (0..n).map(|v| v as f64 - 10.0).collect();
    let w = ChannelWeights {
        channels: vec![Channel {
            source: ChannelSource::PerVertex(values.clone()),
            weight: 2.0,
        }],
    };
    let s = ChannelState::new(&mesh, &w, false).unwrap();
    assert_eq!(s.n_channels(), 2);
    for v in 0..n as u32 {
        assert_eq!(s.mass(v, 0), values[v as usize].max(0.0));
        assert_eq!(s.mass(v, 1), (-values[v as usize]).max(0.0));
    }
    assert_eq!(s.weight(1), 2.0);
    let folded = ChannelState::new(&mesh, &w, true).unwrap();
    assert_eq!(folded.weight(0), 1.0);
    assert_eq!(folded.mass(15, 0), 10.0);
}
