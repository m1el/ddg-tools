//! Barnes–Hut energy and its exact gradient for frozen lists (`SPEC.md` §3.4, §8.2, §13.2,
//! §13.1 (4)).

#![allow(clippy::needless_range_loop)]

use std::time::Instant;

use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::energy::barnes_hut::BhTree;
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::kernel::kernel_of_v;
use ddg_repulsive_surfaces::spatial::{FaceTree, bh_admissible, box_distance_sq};
use ddg_repulsive_surfaces::{
    BhOptions, EnergyEval, Exponents, FaceGeometry, RsError, Surface, energy, read_obj,
};

const KLEIN: &str = "/var/data/discrete-geometry/RepulsiveSurfaces/repulsive-surfaces/scenes/KleinBottle/klein_bottle_blocky.obj";

/// A Möbius band with `n` quads (2 triangles each), one half twist.
fn moebius(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut pos = Vec::new();
    for i in 0..n {
        let th = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let (c, s) = (th.cos(), th.sin());
        let (ch, sh) = ((th / 2.0).cos(), (th / 2.0).sin());
        for w in [0.5, 0.0, -0.5] {
            pos.push([(2.0 + w * ch) * c, (2.0 + w * ch) * s, w * sh]);
        }
    }
    // Rows 0, 1, 2 across the band; the twist maps row r to row 2 − r.
    let at = |i: usize, r: usize| {
        if i == n {
            (2 - r) as u32
        } else {
            (3 * i + r) as u32
        }
    };
    let mut faces = Vec::new();
    for i in 0..n {
        for r in 0..2 {
            let (a, b, c, d) = (at(i, r), at(i, r + 1), at(i + 1, r + 1), at(i + 1, r));
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (pos, faces)
}

fn two_spheres(level: u32) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let m = icosphere(level).unwrap();
    let mut p = m.positions.clone();
    let mut f = m.faces.clone();
    let off = p.len() as u32;
    p.extend(
        m.positions
            .iter()
            .map(|q| [0.7 * q[0] + 2.2, 0.7 * q[1] + 0.3, 0.7 * q[2]]),
    );
    f.extend(m.faces.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    (p, f)
}

type Mesh = (String, Surface, Vec<[f64; 3]>);

fn mesh(name: &str, mut pos: Vec<[f64; 3]>, faces: Vec<[u32; 3]>, amp: f64, seed: u64) -> Mesh {
    let s = Surface::new(pos.len(), faces).unwrap();
    let l = mean_edge_length(&s, &pos);
    perturb(&mut pos, [amp * l; 3], seed);
    (name.to_string(), s, pos)
}

/// Small meshes (a few hundred to ~1300 faces) on which Barnes–Hut lumps nodes at θ = 0.5.
fn test_meshes() -> Vec<Mesh> {
    let mut out = Vec::new();
    let m = icosphere(3).unwrap();
    out.push(mesh("icosphere", m.positions, m.faces, 0.1, 1));
    let m = torus(24, 12, 1.5, 0.6).unwrap();
    out.push(mesh("torus", m.positions, m.faces, 0.1, 2));
    let m = bumpy_grid_disk(16, 16, 0.3).unwrap();
    out.push(mesh("disk", m.positions, m.faces, 0.1, 3));
    let (p, f) = moebius(60);
    out.push(mesh("moebius", p, f, 0.1, 4));
    let (p, f) = two_spheres(2);
    out.push(mesh("two components", p, f, 0.1, 5));
    if std::path::Path::new(KLEIN).exists() {
        let (s, p) = read_obj(KLEIN).unwrap();
        let f = s.faces().to_vec();
        out.push(mesh("klein bottle", p, f, 0.05, 6));
    } else {
        eprintln!("skipping the Klein bottle: {KLEIN} not found");
    }
    out
}

fn bh(theta: f64, leaf_size: usize) -> BhOptions {
    BhOptions { theta, leaf_size }
}

fn rel_diff(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let (mut e2, mut n2) = (0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        for k in 0..3 {
            e2 += (x[k] - y[k]).powi(2);
            n2 += y[k].powi(2);
        }
    }
    (e2 / n2).sqrt()
}

#[test]
fn theta_zero_is_the_exact_energy_and_gradient() {
    let mut worst = (0.0f64, 0.0f64);
    for (name, s, p) in test_meshes() {
        for e in [Exponents::default(), Exponents::new(4.5, 9.5).unwrap()] {
            let ex = energy(&s, &p, e, &EnergyEval::Exact, &[], Some(4)).unwrap();
            for leaf in [1, 3, 8] {
                let r = energy(
                    &s,
                    &p,
                    e,
                    &EnergyEval::BarnesHut(bh(0.0, leaf)),
                    &[],
                    Some(4),
                )
                .unwrap();
                let nf = s.n_faces() as u64;
                assert_eq!(r.nodes_lumped, 0, "{name}");
                assert_eq!(r.pairs_exact, nf * (nf - 1), "{name}");
                let de = (r.energy - ex.energy).abs() / ex.energy;
                let dg = rel_diff(&r.differential, &ex.differential);
                worst = (worst.0.max(de), worst.1.max(dg));
                assert!(
                    de < 1e-13 && dg < 1e-12,
                    "{name} leaf {leaf}: {de:e} {dg:e}"
                );
            }
        }
    }
    eprintln!(
        "theta = 0 vs exact: energy rel {:.2e}, gradient rel {:.2e}",
        worst.0, worst.1
    );
}

#[test]
fn gradient_matches_finite_differences_with_frozen_lists() {
    let cases: Vec<(Mesh, f64, Exponents)> = test_meshes()
        .into_iter()
        .flat_map(|m| {
            [
                (m.clone(), 0.5, Exponents::default()),
                (m.clone(), 1.0, Exponents::new(4.5, 9.5).unwrap()),
            ]
        })
        .collect();
    std::thread::scope(|scope| {
        let handles: Vec<_> = cases
            .iter()
            .map(|((name, s, pos), theta, e)| {
                scope.spawn(move || {
                    let pool = make_pool(Some(2)).unwrap();
                    let mut tree = BhTree::build(&pool, s, pos, bh(*theta, 8)).unwrap();
                    assert!(tree.lists().n_lumped() > 0, "{name}: nothing lumped");
                    let ev = tree.energy_and_differential(&pool, e).unwrap();
                    let g = ev.differential;
                    let h = 1e-6 * mean_edge_length(s, pos);
                    let mut err2 = 0.0;
                    let mut norm2 = 0.0;
                    let mut p = pos.clone();
                    // Every vertex in release builds; about 80 evenly spaced ones in debug.
                    let stride = if cfg!(debug_assertions) {
                        (p.len() / 80).max(1)
                    } else {
                        1
                    };
                    for v in (0..p.len()).step_by(stride) {
                        for k in 0..3 {
                            let x0 = p[v][k];
                            p[v][k] = x0 + h;
                            tree.refit(&p).unwrap();
                            let ep = tree.energy(&pool, e);
                            p[v][k] = x0 - h;
                            tree.refit(&p).unwrap();
                            let em = tree.energy(&pool, e);
                            p[v][k] = x0;
                            let fd = (ep - em) / (2.0 * h);
                            err2 += (fd - g[v][k]).powi(2);
                            norm2 += g[v][k].powi(2);
                        }
                    }
                    let rel = (err2 / norm2).sqrt();
                    eprintln!(
                        "FD {name} theta {theta} ({}, {}), {} faces, {} lumped, {} exact pairs: {rel:.2e}",
                        e.alpha(),
                        e.beta(),
                        s.n_faces(),
                        tree.lists().n_lumped(),
                        tree.lists().pairs_exact()
                    );
                    assert!(rel < 1e-6, "{name} theta {theta}: relative FD error {rel:e}");
                })
            })
            .collect();
        for h in handles {
            h.join().unwrap();
        }
    });
}

/// Lists: for every leaf, the lumped nodes and the exact leaves partition the faces, lumped
/// nodes satisfy the strict admissibility test and exact leaves fail it.
#[test]
fn lists_partition_the_faces_and_respect_admissibility() {
    let pool = make_pool(Some(4)).unwrap();
    for (name, s, p) in test_meshes() {
        for theta in [0.25, 0.5, 1.0, 2.0] {
            let tree = BhTree::build(&pool, &s, &p, bh(theta, 8)).unwrap();
            let ft = tree.face_tree();
            let cl = ft.clusters();
            let mut pairs = 0u64;
            for (k, &l) in ft.leaves().iter().enumerate() {
                let mut seen = vec![0u8; s.n_faces()];
                for &c in tree.lists().lumped_of(k) {
                    assert!(bh_admissible(&cl[l as usize], &cl[c as usize], theta));
                    assert!(box_distance_sq(&cl[l as usize].bbox, &cl[c as usize].bbox) > 0.0);
                    for &f in ft.faces_of(c as usize) {
                        seen[f as usize] += 1;
                    }
                }
                let ns = ft.faces_of(l as usize).len() as u64;
                for &c in tree.lists().exact_of(k) {
                    assert!(ft.bvh().nodes()[c as usize].is_leaf());
                    assert!(!bh_admissible(&cl[l as usize], &cl[c as usize], theta));
                    let nc = ft.faces_of(c as usize).len() as u64;
                    pairs += if c == l { ns * (ns - 1) } else { ns * nc };
                    for &f in ft.faces_of(c as usize) {
                        seen[f as usize] += 1;
                    }
                }
                assert!(
                    seen.iter().all(|&x| x == 1),
                    "{name} theta {theta} leaf {l}"
                );
            }
            assert_eq!(pairs, tree.lists().pairs_exact());
        }
    }
}

/// Node data of `SPEC.md` §3.4 against direct sums.
#[test]
fn cluster_data() {
    let (_, s, p) = test_meshes().remove(4);
    let g = FaceGeometry::new(&s, &p).unwrap();
    let ft = FaceTree::build(&s, &p, &g, 8).unwrap();
    assert!(
        ft.leaves()
            .iter()
            .all(|&l| ft.faces_of(l as usize).len() <= 8)
    );
    for (n, c) in ft.clusters().iter().enumerate() {
        let faces = ft.faces_of(n);
        let a: f64 = faces.iter().map(|&f| g.areas[f as usize]).sum();
        assert!((c.area - a).abs() <= 1e-14 * a);
        let mut x = [0.0; 3];
        let mut pp = [0.0; 6];
        for &f in faces {
            let (af, nf, xf) = (
                g.areas[f as usize],
                g.normals[f as usize],
                g.barycenters[f as usize],
            );
            for k in 0..3 {
                x[k] += af * xf[k] / a;
            }
            let outer = [
                nf[0] * nf[0],
                nf[0] * nf[1],
                nf[0] * nf[2],
                nf[1] * nf[1],
                nf[1] * nf[2],
                nf[2] * nf[2],
            ];
            for k in 0..6 {
                pp[k] += af * outer[k] / a;
            }
            for &v in &s.faces()[f as usize] {
                let q = p[v as usize];
                assert!(c.bbox.contains_point(q));
                let d = ((q[0] - c.centroid[0]).powi(2)
                    + (q[1] - c.centroid[1]).powi(2)
                    + (q[2] - c.centroid[2]).powi(2))
                .sqrt();
                assert!(d <= c.radius * (1.0 + 1e-14));
            }
        }
        for k in 0..3 {
            assert!((c.centroid[k] - x[k]).abs() < 1e-13);
        }
        for k in 0..6 {
            assert!((c.projector[k] - pp[k]).abs() < 1e-13);
        }
        // trace P̄ = 1; vᵀP̄v = area average of (N·v)².
        assert!((c.projector[0] + c.projector[3] + c.projector[5] - 1.0).abs() < 1e-13);
        let v = [0.3, -0.7, 0.2];
        let avg: f64 = faces
            .iter()
            .map(|&f| {
                let nf = g.normals[f as usize];
                g.areas[f as usize] * (nf[0] * v[0] + nf[1] * v[1] + nf[2] * v[2]).powi(2)
            })
            .sum::<f64>()
            / a;
        assert!((c.projector_quadratic_form(v) - avg).abs() < 1e-13);
    }
    // Box distance: a lower bound of point distances, zero on overlap.
    let cl = ft.clusters();
    let (a, b) = (
        &cl[ft.leaves()[0] as usize],
        &cl[*ft.leaves().last().unwrap() as usize],
    );
    let d2 = box_distance_sq(&a.bbox, &b.bbox);
    let corners = |bx: &ddg_repulsive_curves::spatial::Aabb<3>| {
        (0..8)
            .map(|m| {
                std::array::from_fn::<f64, 3, _>(
                    |k| if m >> k & 1 == 1 { bx.hi[k] } else { bx.lo[k] },
                )
            })
            .collect::<Vec<_>>()
    };
    for x in corners(&a.bbox) {
        for y in corners(&b.bbox) {
            let dd: f64 = (0..3).map(|k| (x[k] - y[k]).powi(2)).sum();
            assert!(d2 <= dd);
        }
    }
    assert_eq!(box_distance_sq(&a.bbox, &a.bbox), 0.0);
}

/// `Ẽ` with the lists of a tree, summed naively from the faces (aggregates recomputed from the
/// member faces at the current positions).
fn naive_bh_energy(tree: &BhTree, s: &Surface, p: &[[f64; 3]], e: &Exponents) -> f64 {
    let g = FaceGeometry::new(s, p).unwrap();
    let ft = tree.face_tree();
    let mut sum = 0.0;
    for (k, &l) in ft.leaves().iter().enumerate() {
        for &c in tree.lists().lumped_of(k) {
            let tf = ft.faces_of(c as usize);
            let a: f64 = tf.iter().map(|&t| g.areas[t as usize]).sum();
            let mut x = [0.0; 3];
            for &t in tf {
                for q in 0..3 {
                    x[q] += g.areas[t as usize] * g.barycenters[t as usize][q];
                }
            }
            let x = x.map(|v| v / a);
            for &sf in ft.faces_of(l as usize) {
                let sf = sf as usize;
                let v = std::array::from_fn(|q| g.barycenters[sf][q] - x[q]);
                sum += kernel_of_v(e, v, g.normals[sf]) * g.areas[sf] * a;
            }
        }
        for &c in tree.lists().exact_of(k) {
            for &sf in ft.faces_of(l as usize) {
                for &t in ft.faces_of(c as usize) {
                    if sf == t {
                        continue;
                    }
                    let (sf, t) = (sf as usize, t as usize);
                    let v = std::array::from_fn(|q| g.barycenters[sf][q] - g.barycenters[t][q]);
                    sum += kernel_of_v(e, v, g.normals[sf]) * g.areas[sf] * g.areas[t];
                }
            }
        }
    }
    sum
}

#[test]
fn refit_equals_a_fresh_evaluation_with_the_same_lists() {
    let pool = make_pool(Some(4)).unwrap();
    let e = Exponents::default();
    for (name, s, p0) in test_meshes() {
        let mut tree = BhTree::build(&pool, &s, &p0, bh(0.5, 8)).unwrap();
        let e0 = tree.energy(&pool, &e);
        let ev0 = tree.energy_and_differential(&pool, &e).unwrap();
        assert_eq!(e0.to_bits(), ev0.energy.to_bits(), "{name}");
        let n0 = naive_bh_energy(&tree, &s, &p0, &e);
        assert!((e0 - n0).abs() <= 1e-12 * n0, "{name}: {e0} vs {n0}");
        // A small smooth motion plus noise.
        let l = mean_edge_length(&s, &p0);
        let mut rng = SplitMix64::new(9);
        let p1: Vec<[f64; 3]> = p0
            .iter()
            .map(|q| {
                [
                    q[0] + 0.3 * l * (q[1]).sin() + 0.02 * l * rng.uniform(-1.0, 1.0),
                    q[1] + 0.3 * l * (q[2]).cos() + 0.02 * l * rng.uniform(-1.0, 1.0),
                    q[2] + 0.3 * l * (q[0]).sin() + 0.02 * l * rng.uniform(-1.0, 1.0),
                ]
            })
            .collect();
        tree.refit(&p1).unwrap();
        let e1 = tree.energy(&pool, &e);
        let n1 = naive_bh_energy(&tree, &s, &p1, &e);
        assert!((e1 - n1).abs() <= 1e-12 * n1, "{name}: {e1} vs {n1}");
        assert_ne!(e1, e0);
        let ev1 = tree.energy_and_differential(&pool, &e).unwrap();
        assert_eq!(e1.to_bits(), ev1.energy.to_bits());
        // The geometry, aggregates and positions follow the refit exactly.
        let g1 = FaceGeometry::new(&s, &p1).unwrap();
        assert_eq!(tree.geometry(), &g1);
        // A tree built at p1 has its own lists; a refit back to p0 restores everything bitwise.
        tree.refit(&p0).unwrap();
        let ev2 = tree.energy_and_differential(&pool, &e).unwrap();
        assert_eq!(ev2, ev0, "{name}");
    }
}

#[test]
fn bitwise_identical_across_thread_counts() {
    let m = torus(60, 30, 1.5, 0.6).unwrap();
    let (_, s, p) = mesh("torus", m.positions, m.faces, 0.1, 7);
    let e = Exponents::default();
    let run = |t: usize| {
        energy(
            &s,
            &p,
            e,
            &EnergyEval::BarnesHut(BhOptions::default()),
            &[],
            Some(t),
        )
        .unwrap()
    };
    let r1 = run(1);
    assert!(r1.nodes_lumped > 0);
    for t in [3, 8] {
        let r = run(t);
        assert_eq!(r.energy.to_bits(), r1.energy.to_bits(), "{t} threads");
        assert_eq!(r.differential, r1.differential, "{t} threads");
        assert_eq!(
            (r.pairs_exact, r.nodes_lumped),
            (r1.pairs_exact, r1.nodes_lumped)
        );
    }
    // The trial energy too.
    let vals: Vec<u64> = [1, 3, 8]
        .iter()
        .map(|&t| {
            let pool = make_pool(Some(t)).unwrap();
            BhTree::build(&pool, &s, &p, BhOptions::default())
                .unwrap()
                .energy(&pool, &e)
                .to_bits()
        })
        .collect();
    assert!(vals.iter().all(|&v| v == vals[0]));
    assert_eq!(vals[0], r1.energy.to_bits());
}

#[test]
fn errors_and_non_finite_trials() {
    let pool = make_pool(Some(2)).unwrap();
    // Two coincident triangles in two components: coincident barycenters.
    let s = Surface::new(6, vec![[0, 1, 2], [3, 4, 5]]).unwrap();
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let r = energy(
        &s,
        &p,
        Exponents::default(),
        &EnergyEval::BarnesHut(BhOptions::default()),
        &[],
        Some(2),
    );
    assert!(
        matches!(r, Err(RsError::NonFiniteEnergy { faces: [0, 1] })),
        "{r:?}"
    );
    for (theta, leaf) in [(-1.0, 8), (f64::NAN, 8), (0.5, 0)] {
        let r = energy(
            &s,
            &p,
            Exponents::default(),
            &EnergyEval::BarnesHut(bh(theta, leaf)),
            &[],
            Some(2),
        );
        assert!(matches!(r, Err(RsError::InvalidOption { .. })), "{r:?}");
    }
    // A degenerate trial is a non-finite energy, not an error; a bad count is an error.
    let (_, s, p) = test_meshes().remove(0);
    let mut tree = BhTree::build(&pool, &s, &p, BhOptions::default()).unwrap();
    let mut q = p.clone();
    let t = s.faces()[0];
    q[t[1] as usize] = q[t[0] as usize];
    tree.refit(&q).unwrap();
    assert!(!tree.energy(&pool, &Exponents::default()).is_finite());
    assert!(matches!(
        tree.energy_and_differential(&pool, &Exponents::default()),
        Err(RsError::NonFiniteEnergy { .. })
    ));
    assert!(matches!(
        tree.refit(&q[1..]),
        Err(RsError::PositionCountMismatch { .. })
    ));
    q[0][1] = f64::INFINITY;
    assert!(matches!(
        tree.refit(&q),
        Err(RsError::NonFinitePosition { vertex: 0 })
    ));
}

/// Error of `Ẽ` and of its gradient against the exact ones as `θ` decreases (paper Fig. 10
/// right: `O(θ²)`), on a 10k-face torus. Printed as a table; the decay rate of the energy error at
/// every halving of `θ` must be at least 1.5 (measured: 1.9 to 7.4, faster at small `θ`).
#[test]
fn error_decays_like_theta_squared() {
    let m = torus(100, 50, 1.5, 0.6).unwrap();
    let (_, s, p) = mesh("torus", m.positions, m.faces, 0.02, 8);
    let e = Exponents::default();
    let ex = energy(&s, &p, e, &EnergyEval::Exact, &[], Some(16)).unwrap();
    let thetas = [1.0, 0.5, 0.25, 0.125, 0.0625];
    let mut errs = Vec::new();
    eprintln!("{} faces; exact energy {:.10e}", s.n_faces(), ex.energy);
    eprintln!("theta  | energy rel err | gradient rel err | lumped   | exact pairs");
    for &theta in &thetas {
        let r = energy(
            &s,
            &p,
            e,
            &EnergyEval::BarnesHut(bh(theta, 8)),
            &[],
            Some(16),
        )
        .unwrap();
        let de = (r.energy - ex.energy).abs() / ex.energy;
        let dg = rel_diff(&r.differential, &ex.differential);
        eprintln!(
            "{theta:<6} | {de:14.3e} | {dg:16.3e} | {:8} | {}",
            r.nodes_lumped, r.pairs_exact
        );
        errs.push((de, dg));
    }
    let n = errs.len();
    for i in 1..n {
        let rate_e = (errs[i - 1].0 / errs[i].0).log2();
        let rate_g = (errs[i - 1].1 / errs[i].1).log2();
        eprintln!(
            "rate {} -> {}: energy {rate_e:.2}, gradient {rate_g:.2}",
            thetas[i - 1],
            thetas[i]
        );
    }
    for i in 1..n {
        let rate = (errs[i - 1].0 / errs[i].0).log2();
        assert!(
            rate > 1.5,
            "energy error rate {rate} at theta {}",
            thetas[i]
        );
    }
    // θ = 0.5 (the default): measured 2.3e-2 (energy), 9.2e-3 (gradient).
    assert!(errs[1].0 < 5e-2 && errs[1].1 < 5e-2);
}

/// Speed sanity (`SPEC.md` §14.5 is WP17's): a 40k-face torus at θ = 0.5, build + energy +
/// gradient at 1 and 8 threads. Run alone, in release:
/// `cargo test --release -p ddg-repulsive-surfaces --test barnes_hut -- --ignored --nocapture`.
#[test]
#[ignore]
fn timing_40k_faces() {
    let m = torus(200, 100, 1.5, 0.6).unwrap();
    let (_, s, p) = mesh("torus", m.positions, m.faces, 0.05, 10);
    let e = Exponents::default();
    for threads in [1, 8] {
        let pool = make_pool(Some(threads)).unwrap();
        let mut best = [f64::INFINITY; 3];
        for _ in 0..3 {
            let t0 = Instant::now();
            let mut tree = BhTree::build(&pool, &s, &p, BhOptions::default()).unwrap();
            let t1 = Instant::now();
            let ev = tree.energy_and_differential(&pool, &e).unwrap();
            let t2 = Instant::now();
            tree.refit(&p).unwrap();
            let en = tree.energy(&pool, &e);
            let t3 = Instant::now();
            assert_eq!(en.to_bits(), ev.energy.to_bits());
            best[0] = best[0].min((t1 - t0).as_secs_f64());
            best[1] = best[1].min((t2 - t1).as_secs_f64());
            best[2] = best[2].min((t3 - t2).as_secs_f64());
        }
        let tree = BhTree::build(&pool, &s, &p, BhOptions::default()).unwrap();
        eprintln!(
            "{} faces, {threads} threads: build + lists {:.3} s, energy + gradient {:.3} s, refit + energy {:.3} s ({} lumped, {} exact pairs)",
            s.n_faces(),
            best[0],
            best[1],
            best[2],
            tree.lists().n_lumped(),
            tree.lists().pairs_exact()
        );
    }
}
