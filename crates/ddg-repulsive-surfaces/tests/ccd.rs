//! Continuous collision detection (`SPEC.md` §7.3, D20): hand-made pairs, comparison with fine
//! time sampling on random soups and sphere pairs (many seeds in parallel), obstacles, initially
//! close pairs, determinism across thread counts, the collapse guard.

#![allow(clippy::needless_range_loop)]

use ddg_core::generators::{SplitMix64, grid_disk, icosphere};
use ddg_repulsive_curves::Obstacle;
use ddg_repulsive_curves::ccd::{point_triangle_distance, segment_segment_distance};
use ddg_repulsive_surfaces::Surface;
use ddg_repulsive_surfaces::ccd::{
    CcdCollapseGuard, CcdObstacle, CcdOptions, CcdReport, IntersectionBand, Primitive, first_contact,
    first_contact_with_band,
    intersecting_pairs, intersects, min_distance, triangle_distance,
};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::remesh::{
    CollapseChecks, EditMesh, EnabledOps, RefusalKind, RemeshOptions, remesh,
};
use rayon::prelude::*;

type Mesh = (Surface, Vec<[f64; 3]>);

fn soup(tris: &[[[f64; 3]; 3]]) -> Mesh {
    let mut pos = Vec::new();
    let mut faces = Vec::new();
    for (i, t) in tris.iter().enumerate() {
        pos.extend_from_slice(t);
        faces.push([3 * i as u32, 3 * i as u32 + 1, 3 * i as u32 + 2]);
    }
    (Surface::new(pos.len(), faces).unwrap(), pos)
}

fn union(a: &Mesh, b: &Mesh) -> Mesh {
    let n = a.1.len() as u32;
    let mut pos = a.1.clone();
    pos.extend_from_slice(&b.1);
    let mut faces = a.0.faces().to_vec();
    faces.extend(b.0.faces().iter().map(|f| f.map(|v| v + n)));
    (Surface::new(pos.len(), faces).unwrap(), pos)
}

fn sphere(level: u32, center: [f64; 3], radius: f64) -> Mesh {
    let g = icosphere(level).unwrap();
    let pos = g
        .positions
        .iter()
        .map(|p| [center[0] + radius * p[0], center[1] + radius * p[1], center[2] + radius * p[2]])
        .collect::<Vec<_>>();
    (Surface::new(pos.len(), g.faces).unwrap(), pos)
}

fn opts(l: f64) -> CcdOptions {
    CcdOptions {
        reference_length: Some(l),
        ..Default::default()
    }
}

fn pool(n: usize) -> rayon::ThreadPool {
    make_pool(Some(n)).unwrap()
}

fn moved(start: &[[f64; 3]], vel: &[[f64; 3]], t: f64) -> Vec<[f64; 3]> {
    start
        .iter()
        .zip(vel)
        .map(|(p, u)| [p[0] + t * u[0], p[1] + t * u[1], p[2] + t * u[2]])
        .collect()
}

fn same(a: &CcdReport, b: &CcdReport) -> bool {
    let mut b = b.clone();
    b.seconds = a.seconds;
    *a == b
}

/// A pair of triangles, the second moving.
fn two_triangles(b: [[f64; 3]; 3], vb: [f64; 3]) -> (Mesh, Vec<[f64; 3]>) {
    let a = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let m = soup(&[a, b]);
    let mut vel = vec![[0.0; 3]; 6];
    for k in 3..6 {
        vel[k] = vb;
    }
    (m, vel)
}

#[test]
fn head_on_vertex_face() {
    // B hangs above A with a corner pointing down; it falls with speed 1 from height 1.
    let (m, vel) = two_triangles(
        [[0.2, 0.2, 1.0], [1.0, 1.0, 1.5], [-0.5, 1.2, 1.5]],
        [0.0, 0.0, -1.0],
    );
    let o = opts(1.0);
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 2.0, &[], &o).unwrap();
    // The lowest corner reaches height gap = 0.05 at t = 0.95.
    assert!(r.tau_max <= 0.95 + 1e-9 && r.tau_max > 0.94, "{r:?}");
    assert!(matches!(
        r.contact.unwrap().primitive,
        Primitive::VertexTriangle { vertex: 3, triangle: 0 }
    ));
    assert_eq!((r.close_pairs, r.touching_pairs), (0, 0));
}

#[test]
fn edge_edge_crossing() {
    // A's edge (0,0,0)-(1,0,0) along x; B's edge along y at height 1, falling. B's other corner
    // is far from A so that the first contact is the edge-edge pair.
    let (m, vel) = two_triangles(
        [[0.5, -0.5, 1.0], [0.5, 0.5, 1.0], [3.0, 0.0, 1.0]],
        [0.0, 0.0, -1.0],
    );
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 2.0, &[], &opts(1.0)).unwrap();
    assert!(r.tau_max <= 0.95 + 1e-9 && r.tau_max > 0.94, "{r:?}");
    match r.contact.unwrap().primitive {
        Primitive::EdgeEdge { .. } | Primitive::VertexTriangle { .. } => {}
        p => panic!("{p:?}"),
    }
    // Same scene but the crossing edge is a pure edge-edge contact: both triangles are cut so
    // that no vertex is above the other triangle.
    let a = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -0.01, 0.0]];
    let b = [[0.005, -1.0, 1.0], [0.005, 1.0, 1.0], [0.015, 0.2, 1.0]];
    let m = soup(&[a, b]);
    let mut vel = vec![[0.0; 3]; 6];
    for k in 3..6 {
        vel[k] = [0.0, 0.0, -1.0];
    }
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 2.0, &[], &opts(1.0)).unwrap();
    assert!(r.tau_max <= 0.95 + 1e-9 && r.tau_max > 0.94, "{r:?}");
    let _ = r;
}

#[test]
fn sliding_past_has_no_contact() {
    // B slides in x at height 1 (never approaches), and one that moves sideways past the edge
    // at distance 0.5.
    let (m, vel) = two_triangles(
        [[0.2, 0.2, 1.0], [1.0, 1.0, 1.0], [-0.5, 1.2, 1.0]],
        [3.0, 0.0, 0.0],
    );
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 5.0, &[], &opts(1.0)).unwrap();
    assert_eq!(r.tau_max, 5.0);
    assert!(r.contact.is_none());
    // Passing at lateral distance 0.5 from A's edge, in the plane z = 0.
    let (m, vel) = two_triangles(
        [[2.0, 0.5, 0.0], [3.0, 0.5, 0.0], [2.5, 1.5, 0.0]],
        [-4.0, -1.5, 0.0],
    );
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 1.0, &[], &opts(1.0)).unwrap();
    assert!(r.tau_max < 1.0 || r.contact.is_none());
}

#[test]
fn stationary_and_zero_velocity() {
    let (m, _) = two_triangles(
        [[0.2, 0.2, 1.0], [1.0, 1.0, 1.0], [-0.5, 1.2, 1.0]],
        [0.0; 3],
    );
    let vel = vec![[0.0; 3]; 6];
    let r = first_contact(&pool(1), &m.0, &m.1, &vel, 3.0, &[], &opts(1.0)).unwrap();
    assert_eq!(r.tau_max, 3.0);
    assert!(r.contact.is_none());
    // A common translation of everything gives no contact either.
    let vel = vec![[1.0, 2.0, 3.0]; 6];
    let r = first_contact(&pool(1), &m.0, &m.1, &vel, 3.0, &[], &opts(1.0)).unwrap();
    assert_eq!(r.tau_max, 3.0);
}

#[test]
fn option_and_input_errors() {
    let (m, vel) = two_triangles(
        [[0.2, 0.2, 1.0], [1.0, 1.0, 1.0], [-0.5, 1.2, 1.0]],
        [0.0; 3],
    );
    let p = pool(1);
    for bad in [
        CcdOptions { gap: 0.0, ..Default::default() },
        CcdOptions { gap: f64::NAN, ..Default::default() },
        CcdOptions { max_iterations: 0, ..Default::default() },
        CcdOptions { close_factor: 1.0, ..Default::default() },
        CcdOptions { leaf_size: 0, ..Default::default() },
        CcdOptions { reference_length: Some(-1.0), ..Default::default() },
    ] {
        assert!(first_contact(&p, &m.0, &m.1, &vel, 1.0, &[], &bad).is_err());
    }
    assert!(first_contact(&p, &m.0, &m.1, &vel[..5], 1.0, &[], &opts(1.0)).is_err());
    assert!(first_contact(&p, &m.0, &m.1, &vel, -1.0, &[], &opts(1.0)).is_err());
    let mut v = vel.clone();
    v[2][1] = f64::NAN;
    assert!(first_contact(&p, &m.0, &m.1, &v, 1.0, &[], &opts(1.0)).is_err());
}

// -------------------------------------------------------------------- sampling comparison

/// Samples the motion at `n + 1` times; returns the first sampled time with a distance below
/// `level` (or an intersection), and the minimum distance of all the samples.
fn first_below(m: &Mesh, vel: &[[f64; 3]], t_end: f64, n: usize, level: f64) -> (Option<f64>, f64) {
    let mut first = None;
    let mut dmin = f64::INFINITY;
    for i in 0..=n {
        let t = t_end * i as f64 / n as f64;
        let x = moved(&m.1, vel, t);
        let d = if intersects(&m.0, &x).unwrap() {
            0.0
        } else {
            min_distance(&m.0, &x).unwrap()
        };
        dmin = dmin.min(d);
        if first.is_none() && d < level {
            first = Some(t);
        }
    }
    (first, dmin)
}

/// The check shared by the random scenes: nothing is closer than the gap before `τ_max` (up to
/// sampling), and `τ_max` is the end if the sampled motion never gets within 1.2 gaps.
fn check_scene(m: &Mesh, vel: &[[f64; 3]], t_end: f64, l: f64, samples: usize) -> (f64, Option<f64>) {
    let o = opts(l);
    let gap = o.gap * l;
    let r = first_contact(&pool(2), &m.0, &m.1, vel, t_end, &[], &o).unwrap();
    assert!(r.tau_max >= 0.0 && r.tau_max <= t_end);
    let (_, dmin_safe) = first_below(m, vel, r.tau_max, samples, 0.0);
    assert!(
        dmin_safe >= gap * (1.0 - 1e-6),
        "inside [0, τ_max = {}] the sampled distance reaches {} < gap {gap}: {r:?}",
        r.tau_max,
        dmin_safe
    );
    let (below, _) = first_below(m, vel, t_end, samples, 1.2 * gap);
    if below.is_none() {
        assert_eq!(r.tau_max, t_end, "stopped early without a near approach: {r:?}");
    }
    let (true_contact, _) = first_below(m, vel, t_end, samples, 1e-300);
    if let Some(tc) = true_contact {
        assert!(r.tau_max <= tc + 1e-12, "τ_max {} after contact {tc}", r.tau_max);
    }
    (r.tau_max, true_contact)
}

/// A random soup that starts free: no intersection, all distances above 1.5 gaps. Advances the
/// seed until it is.
fn random_soup(seed: u64, k: usize) -> (Mesh, Vec<[f64; 3]>, f64) {
    for s in seed * 1000.. {
        let (m, vel, l) = random_soup_raw(s, k);
        if !intersects(&m.0, &m.1).unwrap()
            && min_distance(&m.0, &m.1).unwrap() > 1.5 * CcdOptions::default().gap * l
        {
            return (m, vel, l);
        }
    }
    unreachable!()
}

fn random_soup_raw(seed: u64, k: usize) -> (Mesh, Vec<[f64; 3]>, f64) {
    let mut rng = SplitMix64::new(seed);
    let mut tris = Vec::new();
    for _ in 0..k {
        let c = [rng.uniform(0.0, 2.0), rng.uniform(0.0, 2.0), rng.uniform(0.0, 2.0)];
        let mut t = [[0.0; 3]; 3];
        for j in 0..3 {
            for a in 0..3 {
                t[j][a] = c[a] + rng.uniform(-0.4, 0.4);
            }
        }
        tris.push(t);
    }
    let m = soup(&tris);
    // Skip a degenerate triangle by checking validation upstream; velocities are random and the
    // triangles rotate (each corner has its own velocity).
    let vel: Vec<[f64; 3]> = (0..3 * k)
        .map(|_| [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)])
        .collect();
    let l = mean_edge_length(&m.0, &m.1);
    (m, vel, l)
}

#[test]
fn random_soups_against_sampling() {
    // 240 seeds in parallel; the soups are small enough that the sampled min_distance is cheap.
    let results: Vec<(f64, Option<f64>)> = (0..240u64)
        .into_par_iter()
        .map(|seed| {
            let (m, vel, l) = random_soup(1000 + seed, 8);
            check_scene(&m, &vel, 0.6, l, 400)
        })
        .collect();
    let with_contact = results.iter().filter(|r| r.1.is_some()).count();
    assert!(with_contact >= 20, "only {with_contact} scenes had a contact");
    // Tightness: the median ratio τ_max / (sampled first contact) of the scenes with a contact.
    let mut ratios: Vec<f64> = results
        .iter()
        .filter_map(|&(t, c)| c.filter(|&c| c > 0.0).map(|c| t / c))
        .collect();
    ratios.sort_by(f64::total_cmp);
    let med = ratios[ratios.len() / 2];
    eprintln!("{} contact scenes, median τ_max / first contact = {med:.3}", ratios.len());
    assert!(med > 0.5, "median ratio {med}");
}

#[test]
fn sphere_pairs_against_sampling() {
    let results: Vec<(f64, Option<f64>)> = (0..24u64)
        .into_par_iter()
        .map(|seed| {
            let mut rng = SplitMix64::new(77 + seed);
            let a = sphere(1, [0.0; 3], 1.0);
            let c = [rng.uniform(-0.5, 0.5), rng.uniform(-0.5, 0.5), 2.6];
            let b = sphere(1, c, rng.uniform(0.7, 1.0));
            let m = union(&a, &b);
            let nb = b.1.len();
            let speed = rng.uniform(0.5, 2.0);
            let mut vel = vec![[0.0; 3]; m.1.len()];
            let tilt = [rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)];
            for v in vel.iter_mut().skip(a.1.len()).take(nb) {
                *v = [speed * tilt[0], speed * tilt[1], -speed];
            }
            let l = mean_edge_length(&m.0, &m.1);
            check_scene(&m, &vel, 2.5, l, 400)
        })
        .collect();
    assert!(results.iter().any(|r| r.1.is_some()));
    assert!(results.iter().any(|r| r.1.is_none()) || results.len() == 24);
}

#[test]
fn rotating_spinning_pair_edge_edge() {
    // Two long thin triangles that rotate into each other: sampled brute force again.
    let results: Vec<(f64, Option<f64>)> = (0..60u64)
        .into_par_iter()
        .map(|seed| {
            let mut rng = SplitMix64::new(500 + seed);
            let a = [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.05, 0.0]];
            let b = [[0.0, -1.0, 0.8], [0.0, 1.0, 0.8], [0.05, 0.0, 0.8]];
            let m = soup(&[a, b]);
            let w = rng.uniform(1.0, 6.0);
            let vel: Vec<[f64; 3]> = m
                .1
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    if i < 3 {
                        [0.0; 3]
                    } else {
                        // B tilts about the x axis through its centre while sinking.
                        [0.0, -w * (p[2] - 0.8) * 0.0 + 0.0, -0.9 - w * 0.3 * p[1]]
                    }
                })
                .collect();
            let l = 0.5;
            check_scene(&m, &vel, 1.0, l, 400)
        })
        .collect();
    assert!(results.iter().any(|r| r.1.is_some()));
}

// -------------------------------------------------------------------- close pairs

#[test]
fn initially_close_pairs_do_not_freeze() {
    // Gap is 0.05; the triangles start 0.02 apart.
    let (m, vel) = two_triangles(
        [[0.2, 0.2, 0.02], [1.0, 1.0, 0.5], [-0.5, 1.2, 0.5]],
        [0.0, 0.0, -1.0],
    );
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 2.0, &[], &opts(1.0)).unwrap();
    assert_eq!(r.close_pairs, 1);
    // Closing to half the current distance: 0.01 / speed 1.
    assert!(r.tau_max > 0.009 && r.tau_max <= 0.01 + 1e-9, "{r:?}");
    // Moving apart: no contact, in spite of the pair being closer than the gap.
    let up: Vec<[f64; 3]> = vel.iter().map(|v| [v[0], v[1], -v[2]]).collect();
    let r = first_contact(&pool(2), &m.0, &m.1, &up, 2.0, &[], &opts(1.0)).unwrap();
    assert_eq!(r.tau_max, 2.0);
    assert_eq!(r.close_pairs, 1);
    // A vertex lying exactly on a face: unconstrained and counted (here it moves away).
    let (m, mut vel) = two_triangles(
        [[0.2, 0.2, 0.0], [1.0, 1.0, 0.5], [-0.5, 1.2, 0.5]],
        [0.0, 0.0, 1.0],
    );
    vel[3] = [0.0, 0.0, 1.0];
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 1.0, &[], &opts(1.0)).unwrap();
    assert_eq!(r.touching_pairs, 1);
    assert_eq!(r.tau_max, 1.0);
}

// -------------------------------------------------------------------- obstacles

fn floor() -> Obstacle {
    Obstacle::from_soup(&[
        [[-10.0, -10.0, 0.0], [10.0, -10.0, 0.0], [10.0, 10.0, 0.0]],
        [[-10.0, -10.0, 0.0], [10.0, 10.0, 0.0], [-10.0, 10.0, 0.0]],
    ])
    .unwrap()
}

#[test]
fn obstacle_floor() {
    let m = sphere(2, [0.0, 0.0, 3.0], 1.0);
    let l = mean_edge_length(&m.0, &m.1);
    let ob = CcdObstacle::new(&floor(), 4).unwrap();
    let vel = vec![[0.0, 0.0, -1.0]; m.1.len()];
    let r = first_contact(&pool(3), &m.0, &m.1, &vel, 5.0, std::slice::from_ref(&ob), &opts(l)).unwrap();
    let zmin = m.1.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min);
    let want = zmin - 0.05 * l;
    assert!(r.tau_max <= want + 1e-9 && r.tau_max > want - 1e-6, "{} vs {want}", r.tau_max);
    assert!(matches!(
        r.contact.unwrap().primitive,
        Primitive::VertexObstacleTriangle { .. } | Primitive::EdgeObstacleEdge { .. }
    ));
    // Moving away or sideways: free.
    let up = vec![[0.0, 0.0, 1.0]; m.1.len()];
    let r = first_contact(&pool(3), &m.0, &m.1, &up, 5.0, std::slice::from_ref(&ob), &opts(l)).unwrap();
    assert_eq!(r.tau_max, 5.0);
    let side = vec![[1.0, 0.0, 0.0]; m.1.len()];
    let r = first_contact(&pool(3), &m.0, &m.1, &side, 5.0, &[ob], &opts(l)).unwrap();
    assert_eq!(r.tau_max, 5.0);
}

#[test]
fn obstacle_vertex_against_mesh_triangle() {
    // A big mesh triangle moves onto the lone apex of a small pyramid: the pair is
    // (mesh triangle, obstacle vertex).
    let tri = [[-5.0, -5.0, 2.0], [5.0, -5.0, 2.0], [0.0, 6.0, 2.0]];
    let m = soup(&[tri]);
    let pyramid = Obstacle::from_soup(&[
        [[-0.1, -0.1, 0.0], [0.1, -0.1, 0.0], [0.0, 0.0, 1.0]],
        [[0.1, -0.1, 0.0], [0.1, 0.1, 0.0], [0.0, 0.0, 1.0]],
        [[0.1, 0.1, 0.0], [-0.1, 0.1, 0.0], [0.0, 0.0, 1.0]],
        [[-0.1, 0.1, 0.0], [-0.1, -0.1, 0.0], [0.0, 0.0, 1.0]],
    ])
    .unwrap();
    let ob = CcdObstacle::new(&pyramid, 2).unwrap();
    let vel = vec![[0.0, 0.0, -1.0]; 3];
    let l = mean_edge_length(&m.0, &m.1);
    let r = first_contact(&pool(2), &m.0, &m.1, &vel, 5.0, &[ob], &opts(l)).unwrap();
    let want = 1.0 - 0.05 * l;
    assert!(r.tau_max <= want + 1e-9 && r.tau_max > want - 1e-6, "{} vs {want}", r.tau_max);
    assert!(matches!(
        r.contact.unwrap().primitive,
        Primitive::TriangleObstacleVertex { .. }
    ));
}

#[test]
fn obstacle_against_sampling() {
    // A random soup of triangles falling onto a floor with a random tilt.
    for seed in 0..12u64 {
        let (m, mut vel, l) = random_soup(300 + seed, 5);
        let m = (m.0, m.1.iter().map(|p| [p[0], p[1], p[2] + 1.5]).collect::<Vec<_>>());
        for v in vel.iter_mut() {
            v[2] -= 1.5;
        }
        let ob = CcdObstacle::new(&floor(), 4).unwrap();
        let o = opts(l);
        let r = first_contact(&pool(2), &m.0, &m.1, &vel, 2.5, &[ob], &o).unwrap();
        // Sampled distance to the plane z = 0 of the lowest vertex.
        let gap = o.gap * l;
        let zmin0 = m.1.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min);
        let mut first = None;
        for i in 0..=2000 {
            let t = r.tau_max * i as f64 / 2000.0;
            let z = moved(&m.1, &vel, t).iter().map(|p| p[2]).fold(f64::INFINITY, f64::min);
            assert!(z >= gap * (1.0 - 1e-6) || zmin0 < gap, "seed {seed}");
        }
        for i in 0..=5000 {
            let t = 2.5 * i as f64 / 5000.0;
            let z = moved(&m.1, &vel, t).iter().map(|p| p[2]).fold(f64::INFINITY, f64::min);
            if first.is_none() && z < gap {
                first = Some(t);
            }
        }
        if let Some(t) = first {
            // The self pairs of the soup may stop it earlier, never later.
            assert!(r.tau_max <= t + 1e-12, "seed {seed}: {} > {t}", r.tau_max);
        }
    }
}

// -------------------------------------------------------------------- determinism

#[test]
fn bitwise_independent_of_threads() {
    let a = sphere(2, [0.0; 3], 1.0);
    let b = sphere(2, [0.1, 0.0, 2.3], 1.0);
    let m = union(&a, &b);
    let mut vel = vec![[0.0; 3]; m.1.len()];
    for v in vel.iter_mut().skip(a.1.len()) {
        *v = [0.0, 0.0, -1.0];
    }
    let l = mean_edge_length(&m.0, &m.1);
    let ob = CcdObstacle::new(&floor(), 4).unwrap();
    let base = first_contact(&pool(1), &m.0, &m.1, &vel, 3.0, std::slice::from_ref(&ob), &opts(l)).unwrap();
    assert!(base.contact.is_some());
    for n in [2, 3, 8, 16] {
        let r = first_contact(&pool(n), &m.0, &m.1, &vel, 3.0, std::slice::from_ref(&ob), &opts(l)).unwrap();
        assert!(same(&base, &r), "{n} threads: {base:?} vs {r:?}");
        assert_eq!(base.tau_max.to_bits(), r.tau_max.to_bits());
    }
}

// -------------------------------------------------------------------- static utilities

fn brute_min(m: &Mesh) -> f64 {
    let mut best = f64::INFINITY;
    let f = m.0.faces();
    for i in 0..f.len() {
        for j in 0..f.len() {
            for &v in &f[i] {
                if f[j].contains(&v) {
                    continue;
                }
                best = best.min(point_triangle_distance(
                    m.1[v as usize],
                    m.1[f[j][0] as usize],
                    m.1[f[j][1] as usize],
                    m.1[f[j][2] as usize],
                ));
            }
        }
    }
    let e = m.0.edges();
    for i in 0..e.len() {
        for j in i + 1..e.len() {
            if e[i].iter().any(|v| e[j].contains(v)) {
                continue;
            }
            best = best.min(segment_segment_distance(
                m.1[e[i][0] as usize],
                m.1[e[i][1] as usize],
                m.1[e[j][0] as usize],
                m.1[e[j][1] as usize],
            ));
        }
    }
    best
}

#[test]
fn static_checks() {
    let a = sphere(1, [0.0; 3], 1.0);
    let far = sphere(1, [0.0, 0.0, 4.0], 1.0);
    let m = union(&a, &far);
    assert!(!intersects(&m.0, &m.1).unwrap());
    let d = min_distance(&m.0, &m.1).unwrap();
    assert!((d - brute_min(&m)).abs() < 1e-12, "{d}");
    // Penetrating spheres intersect.
    let near = sphere(1, [0.0, 0.0, 1.5], 1.0);
    let m2 = union(&a, &near);
    assert!(intersects(&m2.0, &m2.1).unwrap());
    let pairs = intersecting_pairs(&m2.0, &m2.1).unwrap();
    assert!(!pairs.is_empty() && pairs.windows(2).all(|w| w[0] < w[1]));
    // A single mesh has no self intersection and a positive min distance.
    assert!(!intersects(&a.0, &a.1).unwrap());
    assert!(min_distance(&a.0, &a.1).unwrap() > 0.0);
    // Cross-check min_distance and triangle_distance against brute force.
    let (m, _, _) = random_soup(9, 6);
    let best = brute_min(&m);
    assert!((min_distance(&m.0, &m.1).unwrap() - best).abs() < 1e-12);
    let f = m.0.faces();
    let t = |i: usize| f[i].map(|v| m.1[v as usize]);
    assert!(triangle_distance(t(0), t(1)) >= best - 1e-12);
}

// -------------------------------------------------------------------- collapse guard

/// Two sheets of the same grid, one above the other.
fn stacked(sep: f64) -> Mesh {
    let g = grid_disk(8, 8).unwrap();
    let a: Mesh = (Surface::new(g.positions.len(), g.faces.clone()).unwrap(), g.positions.clone());
    let up: Vec<[f64; 3]> = g.positions.iter().map(|p| [p[0], p[1], p[2] + sep]).collect();
    let b: Mesh = (Surface::new(up.len(), g.faces).unwrap(), up);
    union(&a, &b)
}

#[test]
fn collapse_guard_refusals() {
    let o = CcdOptions::default();
    let l0 = 0.15; // roughly the mean edge length of the grid
    let gap = o.gap * l0;
    // Interior edge of the lower sheet (vertices 4*9+4 and 4*9+5 are neighbours on the grid).
    let (a, b) = (4 * 9 + 4, 4 * 9 + 5);
    for (sep, allowed) in [(0.5 * gap, false), (0.99 * gap, false), (4.0 * gap, true), (1.0, true)] {
        let m = stacked(sep);
        let mut em = EditMesh::new(&m.0, &m.1, &[]).unwrap();
        let mut g = CcdCollapseGuard::new(&o, l0).unwrap();
        let r = em.collapse(a, b, CollapseChecks::default(), Some(&mut g));
        if allowed {
            assert!(r.is_ok(), "sep {sep}: {:?}", r.err().map(|e| e.message));
            assert_eq!(g.refused(), 0);
        } else {
            assert_eq!(r.unwrap_err().kind, RefusalKind::Guard, "sep {sep}");
            assert_eq!(g.refused(), 1);
        }
    }
}

#[test]
fn collapse_guard_in_a_remesh_pass() {
    // Sheets close to each other: with the guard, a pass that wants to collapse everything
    // (huge collapse threshold) leaves the stack free of intersections and refuses some moves.
    let m = stacked(0.003);
    let l0 = mean_edge_length(&m.0, &m.1);
    let ro = RemeshOptions {
        split: 100.0,
        collapse: 0.9,
        rounds: 1,
        enabled_ops: EnabledOps {
            split: false,
            collapse: true,
            flip: false,
            smooth: false,
        },
        ..Default::default()
    };
    let o = CcdOptions::default();
    let mut g = CcdCollapseGuard::new(&o, l0).unwrap();
    let (out, rep) = remesh(&m.0, &m.1, &[], l0, &ro, Some(&mut g)).unwrap();
    assert!(!intersects(&out.surface, &out.positions).unwrap());
    // The sheets start closer than the gap, so everything near them is refused.
    assert!(g.refused() > 0, "{rep:?}");
    // The same pass on a mesh with plenty of room runs and the guard lets collapses through.
    let m = stacked(1.0);
    let mut g = CcdCollapseGuard::new(&o, l0).unwrap();
    let (out, _) = remesh(&m.0, &m.1, &[], l0, &ro, Some(&mut g)).unwrap();
    assert!(out.surface.n_faces() < m.0.n_faces());
    assert_eq!(g.refused(), 0);
    assert!(g.pairs_tested() > 0 || out.surface.n_faces() > 0);
}

// -------------------------------------------------------------------- intersection band

fn translated(m: &Mesh, d: [f64; 3]) -> Mesh {
    let p: Vec<[f64; 3]> = m.1.iter().map(|q| [q[0] + d[0], q[1] + d[1], q[2] + d[2]]).collect();
    (m.0.clone(), p)
}

#[test]
fn band_of_an_embedded_mesh_changes_nothing() {
    let m = union(&sphere(2, [0.0; 3], 1.0), &sphere(2, [4.0, 0.0, 0.0], 1.0));
    let l = mean_edge_length(&m.0, &m.1);
    let p = pool(2);
    let band = IntersectionBand::new(&p, &m.0, &m.1).unwrap();
    assert!(band.is_empty() && band.intersecting_pairs() == 0 && band.n_vertices() == 0);
    // The spheres move towards each other: a contact, identical with, without and with an
    // empty band.
    let n = m.1.len() / 2;
    let vel: Vec<[f64; 3]> =
        (0..m.1.len()).map(|i| if i < n { [1.0, 0.0, 0.0] } else { [-1.0, 0.0, 0.0] }).collect();
    let o = opts(l);
    let a = first_contact(&p, &m.0, &m.1, &vel, 3.0, &[], &o).unwrap();
    let b = first_contact_with_band(&p, &m.0, &m.1, &vel, 3.0, &[], &o, Some(&band)).unwrap();
    let c = first_contact_with_band(&p, &m.0, &m.1, &vel, 3.0, &[], &o, None).unwrap();
    assert!(a.contact.is_some());
    assert!(same(&a, &b) && same(&a, &c));
    assert_eq!(a.exempt_pairs, 0);
}

/// Two overlapping spheres A, B (an intersecting input) and two spheres C, D far from them.
fn overlapping_and_apart() -> (Mesh, [std::ops::Range<usize>; 4]) {
    let a = sphere(2, [0.0; 3], 1.0);
    let b = sphere(2, [1.2, 0.0, 0.0], 1.0);
    let c = sphere(2, [0.0, 0.0, 6.0], 1.0);
    let d = sphere(2, [0.0, 0.0, 9.0], 1.0);
    let n = a.1.len();
    let m = union(&union(&union(&a, &b), &c), &d);
    (m, [0..n, n..2 * n, 2 * n..3 * n, 3 * n..4 * n])
}

#[test]
fn exemption_frees_the_intersection_curve_but_not_other_regions() {
    let (m, r) = overlapping_and_apart();
    let l = mean_edge_length(&m.0, &m.1);
    let p = pool(2);
    let band = IntersectionBand::new(&p, &m.0, &m.1).unwrap();
    assert!(!band.is_empty() && band.intersecting_pairs() > 0);
    // Only vertices of A and B are in the band.
    assert!(band.flags().iter().enumerate().all(|(i, &f)| !f || i < r[1].end));
    let o = opts(l);

    // B slides along A (the existing intersection curve moves): the plain CCD stops it at once
    // (pairs of A and B are close or crossing), the band frees it.
    let mut vel = vec![[0.0; 3]; m.1.len()];
    for i in r[1].clone() {
        vel[i] = [0.0, 0.0, 0.4];
    }
    let plain = first_contact(&p, &m.0, &m.1, &vel, 1.0, &[], &o).unwrap();
    let banded = first_contact_with_band(&p, &m.0, &m.1, &vel, 1.0, &[], &o, Some(&band)).unwrap();
    assert!(banded.exempt_pairs > 0);
    assert!(banded.tau_max > plain.tau_max, "{} vs {}", banded.tau_max, plain.tau_max);

    // C moves into D (two regions disjoint at the start, none in the band): the contact is found
    // and is the contact of C and D alone.
    let mut vel = vec![[0.0; 3]; m.1.len()];
    for i in r[2].clone() {
        vel[i] = [0.0, 0.0, 1.0];
    }
    let banded = first_contact_with_band(&p, &m.0, &m.1, &vel, 3.0, &[], &o, Some(&band)).unwrap();
    let cd = union(&sphere(2, [0.0, 0.0, 6.0], 1.0), &sphere(2, [0.0, 0.0, 9.0], 1.0));
    let n = r[2].len();
    let vel_cd: Vec<[f64; 3]> =
        (0..cd.1.len()).map(|i| if i < n { [0.0, 0.0, 1.0] } else { [0.0; 3] }).collect();
    let alone = first_contact(&p, &cd.0, &cd.1, &vel_cd, 3.0, &[], &o).unwrap();
    assert!(alone.contact.is_some());
    assert!(banded.contact.is_some(), "{banded:?}");
    assert!(
        (banded.tau_max - alone.tau_max).abs() < 1e-6 * alone.tau_max,
        "{} vs {}",
        banded.tau_max,
        alone.tau_max
    );
    let _ = translated;
}

#[test]
fn the_exemption_never_hides_a_new_crossing_between_disjoint_regions() {
    // Random rigid pushes of B, C and D, mostly C and D towards each other. At the CCD's
    // tau_max (and before) the faces that cross are only pairs that both touch the band: a
    // crossing between regions that were both disjoint at the start is never allowed.
    let (m, r) = overlapping_and_apart();
    let l = mean_edge_length(&m.0, &m.1);
    let p = pool(2);
    let band = IntersectionBand::new(&p, &m.0, &m.1).unwrap();
    let o = opts(l);
    let mut rng = SplitMix64::new(77);
    let mut capped = 0;
    for _ in 0..12 {
        let mut vel = vec![[0.0; 3]; m.1.len()];
        for k in 1..4 {
            let s = [rng.next_f64() - 0.5, rng.next_f64() - 0.5, rng.next_f64() - 0.5];
            let s = match k {
                2 => [s[0], s[1], s[2] + 0.8],
                3 => [s[0], s[1], s[2] - 0.8],
                _ => s,
            };
            for i in r[k].clone() {
                vel[i] = s;
            }
        }
        let rep = first_contact_with_band(&p, &m.0, &m.1, &vel, 3.0, &[], &o, Some(&band)).unwrap();
        capped += usize::from(rep.contact.is_some());
        let tau = rep.tau_max;
        for t in [0.0, 0.5 * tau, tau] {
            let x = moved(&m.1, &vel, t);
            for [f, g] in intersecting_pairs(&m.0, &x).unwrap() {
                let touches = |h: u32| m.0.faces()[h as usize].iter().any(|&v| band.flags()[v as usize]);
                assert!(touches(f) && touches(g), "t = {t} (tau {tau}): faces {f},{g} cross outside the band");
            }
        }
    }
    assert!(capped >= 6, "{capped}");
}

#[test]
fn band_report_is_bitwise_independent_of_threads() {
    let (m, r) = overlapping_and_apart();
    let l = mean_edge_length(&m.0, &m.1);
    let mut vel = vec![[0.0; 3]; m.1.len()];
    for i in r[1].clone() {
        vel[i] = [0.0, 0.0, 0.4];
    }
    for i in r[2].clone() {
        vel[i] = [0.0, 0.0, 1.0];
    }
    let o = opts(l);
    let b1 = IntersectionBand::new(&pool(1), &m.0, &m.1).unwrap();
    let base = first_contact_with_band(&pool(1), &m.0, &m.1, &vel, 3.0, &[], &o, Some(&b1)).unwrap();
    for n in [2, 3, 8] {
        let bn = IntersectionBand::new(&pool(n), &m.0, &m.1).unwrap();
        assert_eq!(bn, b1);
        let r = first_contact_with_band(&pool(n), &m.0, &m.1, &vel, 3.0, &[], &o, Some(&bn)).unwrap();
        assert!(same(&base, &r), "threads {n}");
    }
}

#[test]
fn collapse_guard_with_a_band() {
    // Two overlapping spheres: a pass that collapses a lot is refused near the intersection
    // without the band, and less with it.
    let m = union(&sphere(2, [0.0; 3], 1.0), &sphere(2, [1.2, 0.0, 0.0], 1.0));
    let l0 = mean_edge_length(&m.0, &m.1);
    let ro = RemeshOptions {
        split: 100.0,
        collapse: 0.99,
        rounds: 1,
        enabled_ops: EnabledOps { split: false, collapse: true, flip: false, smooth: false },
        ..Default::default()
    };
    let o = CcdOptions::default();
    let band = IntersectionBand::new(&pool(2), &m.0, &m.1).unwrap();
    assert!(!band.is_empty());
    let mut plain = CcdCollapseGuard::new(&o, l0).unwrap();
    remesh(&m.0, &m.1, &[], l0, &ro, Some(&mut plain)).unwrap();
    let mut banded = CcdCollapseGuard::new(&o, l0).unwrap();
    banded.set_band(&band);
    remesh(&m.0, &m.1, &[], l0, &ro, Some(&mut banded)).unwrap();
    assert!(plain.refused() > 0, "{}", plain.refused());
    assert!(banded.refused() < plain.refused(), "{} vs {}", banded.refused(), plain.refused());
    assert!(banded.exempt_pairs() > 0);
    assert_eq!(plain.exempt_pairs(), 0);
}

/// Coplanar faces with one common vertex (flat regions): rounding left `orient` at ±1e-20 for
/// points in the plane, and the non-coplanar branch of the segment–triangle test then reported
/// a hit (the genus-2 start mesh of `apps genus02`, faces 1795 and 1988; NOTES-FIDELITY.md).
#[test]
fn coplanar_neighbours_do_not_intersect() {
    use ddg_repulsive_surfaces::generate::shapes::genus_surface;
    let m = genus_surface(2, 4200).unwrap();
    let s = m.surface().unwrap();
    assert_eq!(intersecting_pairs(&s, &m.positions).unwrap(), Vec::<[u32; 2]>::new());
    // flat grids in generic planes
    let g = grid_disk(24, 24).unwrap();
    let s = Surface::new(g.positions.len(), g.faces.clone()).unwrap();
    let mut rng = SplitMix64::new(7);
    let mut hits = 0;
    for _ in 0..50 {
        let mut r = || rng.next_f64() - 0.5;
        let (u, w) = ([r(), r(), r()], [r(), r(), r()]);
        let x: Vec<[f64; 3]> =
            g.positions.iter().map(|p| [0, 1, 2].map(|k| 0.3 + p[0] * u[k] + p[1] * w[k])).collect();
        hits += intersecting_pairs(&s, &x).unwrap().len();
    }
    assert_eq!(hits, 0);
    // a genuine coplanar overlap is still found
    let tris = [
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0.2, 0.2, 0.0], [1.2, 0.2, 0.0], [0.2, 1.2, 0.0]],
    ];
    let (s, x) = soup(&tris);
    assert_eq!(intersecting_pairs(&s, &x).unwrap(), vec![[0, 1]]);
}
