//! Penalties, obstacles, implicit surfaces, Willmore and Coulomb (`SPEC.md` §5.2, §5.3, §3.5).

#![allow(clippy::needless_range_loop)]

use std::sync::Arc;

use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_curves::obstacle::Obstacle;
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::{FaceGeometry, mean_edge_length, vertex_areas};
use ddg_repulsive_surfaces::implicit::{Cylinder, FramedTorus, ImplicitSurface, Sphere};
use ddg_repulsive_surfaces::obstacle::{CrossEvaluator, ObstacleKind, StaticSet, obstacle_diagonal};
use ddg_repulsive_surfaces::penalties::{
    Penalty, Quadrature, coulomb_energy, evaluate_one, willmore_and_gradient, willmore_matrix,
};
use ddg_repulsive_surfaces::{BhOptions, EnergyEval, Exponents, RsError, Surface, energy};

type Mesh = (&'static str, Surface, Vec<[f64; 3]>);

fn moebius(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut pos = Vec::new();
    for i in 0..n {
        let th = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let (c, s) = (th.cos(), th.sin());
        let (ch, sh) = ((th / 2.0).cos(), (th / 2.0).sin());
        for w in [0.5, -0.5] {
            pos.push([(2.0 + w * ch) * c, (2.0 + w * ch) * s, w * sh]);
        }
    }
    let t = |i: usize| (2 * (i % n)) as u32;
    let b = |i: usize| (2 * (i % n) + 1) as u32;
    let mut faces = Vec::new();
    for i in 0..n {
        let (ti, bi) = (t(i), b(i));
        let (tj, bj) = if i + 1 == n { (b(0), t(0)) } else { (t(i + 1), b(i + 1)) };
        faces.push([ti, bi, bj]);
        faces.push([ti, bj, tj]);
    }
    (pos, faces)
}

fn make(
    name: &'static str,
    mut pos: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    seed: u64,
) -> Mesh {
    let s = Surface::new(pos.len(), faces).unwrap();
    let l = mean_edge_length(&s, &pos);
    perturb(&mut pos, [0.1 * l; 3], seed);
    (name, s, pos)
}

fn union(
    a: (Vec<[f64; 3]>, Vec<[u32; 3]>),
    b: (Vec<[f64; 3]>, Vec<[u32; 3]>),
    shift: [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let (mut p, mut f) = a;
    let off = p.len() as u32;
    p.extend(b.0.iter().map(|q| [q[0] + shift[0], q[1] + shift[1], q[2] + shift[2]]));
    f.extend(b.1.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    (p, f)
}

fn ico(level: u32) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let m = icosphere(level).unwrap();
    (m.positions, m.faces)
}

fn disk() -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let m = bumpy_grid_disk(4, 4, 0.3).unwrap();
    (m.positions, m.faces)
}

/// Closed meshes (volume penalties allowed).
fn closed_meshes() -> Vec<Mesh> {
    let t = torus(8, 5, 1.5, 0.6).unwrap();
    let (p, f) = union(ico(0), ico(0), [3.0, 0.2, 0.1]);
    vec![
        {
            let (p, f) = ico(1);
            make("icosphere", p, f, 1)
        },
        make("torus", t.positions, t.faces, 2),
        make("two spheres", p, f, 3),
    ]
}

/// Meshes with boundary.
fn open_meshes() -> Vec<Mesh> {
    let (p, f) = union(disk(), disk(), [0.3, 2.9, 0.5]);
    let (mp, mf) = moebius(10);
    vec![
        {
            let (p, f) = disk();
            make("disk", p, f, 4)
        },
        make("two disks", p, f, 5),
        make("moebius", mp, mf, 6),
    ]
}

fn all_meshes() -> Vec<Mesh> {
    let mut v = closed_meshes();
    v.extend(open_meshes());
    v
}

/// Central finite-difference check of the weighted value and gradient of one penalty;
/// returns the relative error `max|g − fd| / max|g|`.
fn fd_error(surface: &Surface, pos: &[[f64; 3]], p: &Penalty, eval: &EnergyEval) -> f64 {
    let pool = make_pool(Some(2)).unwrap();
    let e = Exponents::default();
    let (_, g) = evaluate_one(&pool, surface, pos, p, &e, eval).unwrap();
    let h = 1e-6 * mean_edge_length(surface, pos);
    let n = pos.len();
    let stride = (n / 40).max(1);
    let (mut max_err, mut max_g) = (0.0f64, 0.0f64);
    for v in (0..n).step_by(stride) {
        for c in 0..3 {
            let mut q = pos.to_vec();
            q[v][c] += h;
            let (fp, _) = evaluate_one(&pool, surface, &q, p, &e, eval).unwrap();
            q[v][c] -= 2.0 * h;
            let (fm, _) = evaluate_one(&pool, surface, &q, p, &e, eval).unwrap();
            let fd = (fp - fm) / (2.0 * h);
            max_err = max_err.max((fd - g[v][c]).abs());
        }
    }
    for v in (0..n).step_by(stride) {
        for c in 0..3 {
            max_g = max_g.max(g[v][c].abs());
        }
    }
    assert!(max_g > 0.0, "gradient vanishes");
    max_err / max_g
}

fn assert_fd(label: &str, surface: &Surface, pos: &[[f64; 3]], p: &Penalty) {
    let err = fd_error(surface, pos, p, &EnergyEval::Exact);
    println!("{label}: fd error {err:.2e}");
    assert!(err < 1e-6, "{label}: relative FD error {err:e}");
}

fn area_of(s: &Surface, p: &[[f64; 3]]) -> f64 {
    FaceGeometry::new(s, p).unwrap().total_area()
}

#[test]
fn geometric_penalties_match_finite_differences() {
    for (name, s, pos) in all_meshes() {
        let a0 = area_of(&s, &pos);
        assert_fd(
            &format!("{name} AreaDeviation"),
            &s,
            &pos,
            &Penalty::AreaDeviation { weight: 3.0, target: 0.8 * a0 },
        );
        for power in [1.0, 2.0, 1.5] {
            assert_fd(
                &format!("{name} TotalArea {power}"),
                &s,
                &pos,
                &Penalty::TotalArea { weight: 0.7, power },
            );
        }
        let mut rng = SplitMix64::new(9);
        let targets: Vec<[f64; 3]> = pos
            .iter()
            .map(|q| std::array::from_fn(|k| q[k] + rng.uniform(-0.2, 0.2)))
            .collect();
        assert_fd(
            &format!("{name} SquaredError"),
            &s,
            &pos,
            &Penalty::SquaredError { weight: 2.0, targets },
        );
        // Every vertex of the Moebius strip is a boundary vertex: no interior rows.
        if name != "moebius" {
            assert_fd(
                &format!("{name} Willmore"),
                &s,
                &pos,
                &Penalty::Willmore { weight: 1.3 },
            );
        }
    }
}

#[test]
fn volume_penalties_match_finite_differences() {
    for (name, s, pos) in closed_meshes() {
        let (v0, _) = ddg_repulsive_surfaces::penalties::simple::volume_and_gradient(&s, &pos);
        assert!(v0.abs() > 0.0);
        assert_fd(
            &format!("{name} VolumeDeviation"),
            &s,
            &pos,
            &Penalty::VolumeDeviation { weight: 2.0, target: 1.3 * v0 },
        );
        for power in [1.0, 2.0, 3.0] {
            assert_fd(
                &format!("{name} TotalVolume {power}"),
                &s,
                &pos,
                &Penalty::TotalVolume { weight: 0.5, power },
            );
        }
    }
    // Non-orientable component: error naming the component and an edge.
    let (mp, mf) = moebius(10);
    let s = Surface::new(mp.len(), mf).unwrap();
    let pool = make_pool(Some(1)).unwrap();
    let err = evaluate_one(
        &pool,
        &s,
        &mp,
        &Penalty::TotalVolume { weight: 1.0, power: 1.0 },
        &Exponents::default(),
        &EnergyEval::Exact,
    )
    .unwrap_err();
    assert!(matches!(err, RsError::VolumeNeedsOrientation { component: 0, .. }), "{err}");
}

#[test]
fn boundary_penalties_match_finite_differences() {
    for (name, s, pos) in open_meshes() {
        let (l, _) = ddg_repulsive_surfaces::penalties::simple::boundary_length_and_gradient(&s, &pos);
        assert_fd(
            &format!("{name} BoundaryLength"),
            &s,
            &pos,
            &Penalty::BoundaryLength { weight: 1.5, target: 0.7 * l },
        );
        for dual_length in [true, false] {
            assert_fd(
                &format!("{name} BoundaryCurvature dual={dual_length}"),
                &s,
                &pos,
                &Penalty::BoundaryCurvature { weight: 0.9, dual_length },
            );
        }
    }
}

#[test]
fn scaling_and_invariance_sanity() {
    let pool = make_pool(Some(2)).unwrap();
    let e = Exponents::default();
    let (_, s, pos) = closed_meshes().remove(1);
    let c = 2.7;
    let scaled: Vec<[f64; 3]> = pos.iter().map(|p| [c * p[0], c * p[1], c * p[2]]).collect();
    let a0 = area_of(&s, &pos);
    let val = |s: &Surface, p: &[[f64; 3]], pen: &Penalty| {
        evaluate_one(&pool, s, p, pen, &e, &EnergyEval::Exact).unwrap().0
    };
    // AreaDeviation is invariant when the target scales with the area.
    let dev = |a: f64, q: &[[f64; 3]]| val(&s, q, &Penalty::AreaDeviation { weight: 1.0, target: a });
    let (d1, d2) = (dev(0.8 * a0, &pos), dev(0.8 * a0 * c * c, &scaled));
    assert!((d1 - d2).abs() < 1e-12 * d1, "{d1} {d2}");
    // TotalArea scales as c^(2 power).
    let ta = |q: &[[f64; 3]]| val(&s, q, &Penalty::TotalArea { weight: 1.0, power: 2.0 });
    assert!((ta(&scaled) / ta(&pos) - c.powi(4)).abs() < 1e-12 * c.powi(4));
    // Willmore is scale invariant.
    let w = |q: &[[f64; 3]]| val(&s, q, &Penalty::Willmore { weight: 1.0 });
    assert!((w(&scaled) / w(&pos) - 1.0).abs() < 1e-11);
    // Willmore is rigid-motion invariant.
    let moved: Vec<[f64; 3]> = pos.iter().map(|p| [-p[1] + 5.0, p[0] - 1.0, p[2] + 0.3]).collect();
    assert!((w(&moved) / w(&pos) - 1.0).abs() < 1e-11);
    // Weights scale values linearly; zero weight gives zero; bad weights are errors.
    let w2 = val(&s, &pos, &Penalty::Willmore { weight: 2.0 });
    assert!((w2 / w(&pos) - 2.0).abs() < 1e-14);
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let r = evaluate_one(&pool, &s, &pos, &Penalty::Willmore { weight: bad }, &e, &EnergyEval::Exact);
        assert!(matches!(r, Err(RsError::InvalidOption { name: "weight", .. })));
    }
}

#[test]
fn willmore_of_icosphere_is_close_to_the_smooth_value() {
    // P = fᵀ W M⁻¹ W f ≈ ∫ |Δf|² dA = 4 ∫ H² dA = 16π R⁰ on a round sphere (H = 1/R mean
    // curvature, |Δf| = 2/R). The reference's own value would be ¼ of this (4π).
    let want = 16.0 * std::f64::consts::PI;
    let mut errs = Vec::new();
    for level in 2..=4 {
        let (p, f) = ico(level);
        let s = Surface::new(p.len(), f).unwrap();
        for r in [1.0, 3.0] {
            let q: Vec<[f64; 3]> = p.iter().map(|x| [r * x[0], r * x[1], r * x[2]]).collect();
            let (v, _) = willmore_and_gradient(&s, &q).unwrap();
            let err = (v / want - 1.0).abs();
            println!("level {level} R {r}: P = {v:.6}, 16π = {want:.6}, rel {err:.2e}");
            if r == 1.0 {
                errs.push(err);
            }
        }
    }
    assert!(errs[2] < 0.01, "{errs:?}");
    assert!(errs[0] > errs[1] && errs[1] > errs[2], "{errs:?}");
}

#[test]
fn willmore_matrix_reproduces_the_value() {
    let (_, s, pos) = open_meshes().remove(0);
    let (value, _) = willmore_and_gradient(&s, &pos).unwrap();
    let m = willmore_matrix(&s, &pos, true).unwrap();
    assert!(m.symmetry_defect().unwrap() < 1e-12 * m.values().iter().fold(0.0f64, |a, b| a.max(b.abs())));
    let n = pos.len();
    let mut q = 0.0;
    for c in 0..3 {
        let x: Vec<f64> = pos.iter().map(|p| p[c]).collect();
        let mut y = vec![0.0; n];
        m.mul_vec(&x, &mut y).unwrap();
        q += x.iter().zip(&y).map(|(a, b)| a * b).sum::<f64>();
    }
    assert!((q - value).abs() < 1e-10 * value, "{q} {value}");
    // On a closed mesh, both variants agree; the full matrix has the constants in its kernel.
    let (_, s, pos) = closed_meshes().remove(0);
    let a = willmore_matrix(&s, &pos, true).unwrap();
    let b = willmore_matrix(&s, &pos, false).unwrap();
    assert_eq!(a, b);
    let ones = vec![1.0; pos.len()];
    let mut y = vec![0.0; pos.len()];
    b.mul_vec(&ones, &mut y).unwrap();
    assert!(y.iter().all(|v| v.abs() < 1e-10));
}

// ---------------------------------------------------------------------------------------------
// implicit surfaces

#[test]
fn cylinder_and_framed_torus_are_signed_distances_with_exact_gradients() {
    let cyl = Cylinder::new([0.3, -0.2, 0.1], [1.0, 2.0, -0.5], 0.7).unwrap();
    let ft = FramedTorus::new([0.1, 0.2, -0.3], [0.3, 1.0, 0.5], 1.6, 0.4).unwrap();
    let mut rng = SplitMix64::new(77);
    for (name, s) in [("cylinder", &cyl as &dyn ImplicitSurface), ("torus", &ft)] {
        for _ in 0..200 {
            let x: [f64; 3] = std::array::from_fn(|_| rng.uniform(-3.0, 3.0));
            let g = s.gradient(x);
            let gn = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            assert!((gn - 1.0).abs() < 1e-12, "{name}: |grad| = {gn}");
            for c in 0..3 {
                let h = 1e-6;
                let (mut a, mut b) = (x, x);
                a[c] += h;
                b[c] -= h;
                let fd = (s.value(a) - s.value(b)) / (2.0 * h);
                assert!((fd - g[c]).abs() < 1e-7, "{name}: {fd} {}", g[c]);
            }
        }
    }
    // Known values.
    let c = Cylinder::new([0.0; 3], [0.0, 0.0, 1.0], 1.0).unwrap();
    assert_eq!(c.value([3.0, 0.0, 9.0]), 2.0);
    let t = FramedTorus::new([0.0; 3], [0.0, 1.0, 0.0], 2.0, 0.5).unwrap();
    assert!((t.value([2.0, 0.5, 0.0])).abs() < 1e-15);
    assert!((t.value([0.0, 0.0, 0.0]) - 1.5).abs() < 1e-15);
    // Reference torus agreement.
    let rt = ddg_repulsive_curves::surfaces::Torus::new([0.0; 3], 1.0, 0.25).unwrap();
    let ft2 = FramedTorus::new([0.0; 3], [0.0, 1.0, 0.0], 1.0, 0.25).unwrap();
    for _ in 0..50 {
        let x: [f64; 3] = std::array::from_fn(|_| rng.uniform(-2.0, 2.0));
        assert!((rt.value(x) - ft2.value(x)).abs() < 1e-14);
    }
    // Singular points: finite (zero) gradients.
    for p in [[0.0, 0.0, 5.0], [0.3, -0.2, 0.1]] {
        assert!(cyl.gradient(p).iter().all(|x| x.is_finite()));
    }
    for p in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [2.0, 0.0, 0.0]] {
        assert!(t.gradient(p).iter().all(|x| x.is_finite()), "{p:?}");
    }
    assert_eq!(t.gradient([2.0, 0.0, 0.0]), [0.0; 3]); // core circle
    // Errors.
    assert!(Cylinder::new([0.0; 3], [0.0; 3], 1.0).is_err());
    assert!(Cylinder::new([0.0; 3], [1.0, 0.0, 0.0], -1.0).is_err());
    assert!(FramedTorus::new([0.0; 3], [0.0, 1.0, 0.0], 1.0, f64::NAN).is_err());
}

fn far_sphere() -> Arc<dyn ImplicitSurface> {
    Arc::new(Sphere::new([6.0, 1.0, 0.5], 1.0).unwrap())
}

#[test]
fn implicit_potentials_match_finite_differences() {
    let shapes: Vec<(&str, Arc<dyn ImplicitSurface>)> = vec![
        ("sphere", far_sphere()),
        ("cylinder", Arc::new(Cylinder::new([6.0, 1.0, 0.0], [0.0, 0.3, 1.0], 1.0).unwrap())),
        ("torus", Arc::new(FramedTorus::new([6.5, 1.0, 0.0], [0.2, 1.0, 0.1], 1.0, 0.3).unwrap())),
    ];
    for (name, s, pos) in all_meshes() {
        for (sn, shape) in &shapes {
            for q in [Quadrature::DualArea, Quadrature::PerVertex] {
                assert_fd(
                    &format!("{name} obstacle {sn} {q:?}"),
                    &s,
                    &pos,
                    &Penalty::ImplicitObstacle {
                        weight: 2.0,
                        surface: shape.clone(),
                        power: 2.0,
                        quadrature: q,
                    },
                );
                let sub: Vec<u32> = (0..pos.len() as u32).filter(|v| v % 3 != 1).collect();
                assert_fd(
                    &format!("{name} attractor {sn} {q:?}"),
                    &s,
                    &pos,
                    &Penalty::ImplicitAttractor {
                        weight: 0.5,
                        surface: shape.clone(),
                        power: 3.0,
                        quadrature: q,
                        vertices: Some(sub),
                    },
                );
            }
        }
    }
}

#[test]
fn implicit_obstacle_infinity_and_quadrature() {
    let pool = make_pool(Some(2)).unwrap();
    let e = Exponents::default();
    let (p, f) = ico(2);
    let s = Surface::new(p.len(), f).unwrap();
    // A sphere of radius 0.5 inside the unit icosphere: d = 0.5 on every vertex (up to the
    // icosphere's exact radius 1).
    let inner: Arc<dyn ImplicitSurface> = Arc::new(Sphere::new([0.0; 3], 0.5).unwrap());
    let run = |q: Quadrature, shape: &Arc<dyn ImplicitSurface>, pos: &[[f64; 3]]| {
        evaluate_one(
            &pool,
            &s,
            pos,
            &Penalty::ImplicitObstacle { weight: 1.0, surface: shape.clone(), power: 2.0, quadrature: q },
            &e,
            &EnergyEval::Exact,
        )
        .unwrap()
    };
    let (v_pv, _) = run(Quadrature::PerVertex, &inner, &p);
    assert!((v_pv - p.len() as f64 * 4.0).abs() < 1e-9 * v_pv, "{v_pv}");
    let a = area_of(&s, &p);
    let (v_da, _) = run(Quadrature::DualArea, &inner, &p);
    assert!((v_da - a * 4.0).abs() < 1e-10 * v_da, "{v_da}");
    // Refinement stability: the dual-area value converges to 4π·4, the per-vertex value grows.
    let (p3, f3) = ico(3);
    let s3 = Surface::new(p3.len(), f3).unwrap();
    let r3 = evaluate_one(
        &pool,
        &s3,
        &p3,
        &Penalty::ImplicitObstacle { weight: 1.0, surface: inner.clone(), power: 2.0, quadrature: Quadrature::DualArea },
        &e,
        &EnergyEval::Exact,
    )
    .unwrap()
    .0;
    assert!((r3 / v_da - 1.0).abs() < 0.05);
    // A vertex inside the obstacle: +infinity and no gradient.
    let outer: Arc<dyn ImplicitSurface> = Arc::new(Sphere::new([0.0; 3], 1.2).unwrap());
    let (v, g) = run(Quadrature::DualArea, &outer, &p);
    assert_eq!(v, f64::INFINITY);
    assert!(g.iter().flatten().all(|x| *x == 0.0));
    let (v, _) = run(Quadrature::PerVertex, &outer, &p);
    assert_eq!(v, f64::INFINITY);
    // Through energy(): infinite penalty, finite differential.
    let rep = energy(&s, &p, e, &EnergyEval::Exact, &[Penalty::ImplicitObstacle {
        weight: 1.0, surface: outer, power: 2.0, quadrature: Quadrature::PerVertex }], Some(2)).unwrap();
    assert_eq!(rep.penalties, vec![f64::INFINITY]);
    assert!(rep.differential.iter().flatten().all(|x| x.is_finite()));
    // Attractor vertex validation.
    let r = evaluate_one(
        &pool, &s, &p,
        &Penalty::ImplicitAttractor { weight: 1.0, surface: inner.clone(), power: 2.0,
            quadrature: Quadrature::PerVertex, vertices: Some(vec![1, 1]) },
        &e, &EnergyEval::Exact);
    assert!(matches!(r, Err(RsError::InvalidOption { .. })));
    // Attractor on all vertices equals the subset of all.
    let all: Vec<u32> = (0..p.len() as u32).collect();
    let mk = |vertices| Penalty::ImplicitAttractor { weight: 1.0, surface: inner.clone(), power: 2.0,
        quadrature: Quadrature::DualArea, vertices };
    let a1 = evaluate_one(&pool, &s, &p, &mk(None), &e, &EnergyEval::Exact).unwrap();
    let a2 = evaluate_one(&pool, &s, &p, &mk(Some(all)), &e, &EnergyEval::Exact).unwrap();
    assert_eq!(a1.0, a2.0);
    assert!((a1.0 - a * 0.25).abs() < 1e-10);
}

// ---------------------------------------------------------------------------------------------
// obstacles

fn shifted_ico(level: u32, scale: f64, shift: [f64; 3]) -> Obstacle {
    let (p, f) = ico(level);
    let p: Vec<[f64; 3]> = p
        .iter()
        .map(|q| [scale * q[0] + shift[0], scale * q[1] + shift[1], scale * q[2] + shift[2]])
        .collect();
    Obstacle::new(p, f).unwrap()
}

fn bh(theta: f64, leaf: usize) -> BhOptions {
    BhOptions { theta, leaf_size: leaf }
}

/// Brute-force double sum of the cross energies.
fn brute(
    s: &Surface,
    pos: &[[f64; 3]],
    set: &StaticSet,
    normals: Option<&[[f64; 3]]>,
    kind: ObstacleKind,
    e: &Exponents,
) -> f64 {
    use ddg_repulsive_surfaces::kernel::kernel;
    let g = FaceGeometry::new(s, pos).unwrap();
    let mut tot = 0.0;
    for f in 0..s.n_faces() {
        for t in 0..set.len() {
            let (xs, xt) = (g.barycenters[f], set.points()[t]);
            let w = g.areas[f] * set.weights()[t];
            match kind {
                ObstacleKind::TangentPoint => {
                    tot += w * kernel(e, xs, g.normals[f], xt);
                    if let Some(n) = normals {
                        tot += w * kernel(e, xt, n[t], xs);
                    }
                }
                ObstacleKind::Coulomb { power } => {
                    let d = ((xs[0] - xt[0]).powi(2) + (xs[1] - xt[1]).powi(2) + (xs[2] - xt[2]).powi(2)).sqrt();
                    tot += w * d.powf(-power);
                }
            }
        }
    }
    tot
}

fn triangle_normals(o: &Obstacle) -> Vec<[f64; 3]> {
    let v = o.vertices();
    o.triangle_indices()
        .iter()
        .map(|t| ddg_repulsive_surfaces::geometry::face_quantities(t.map(|i| v[i as usize])).1)
        .collect()
}

fn cases() -> Vec<(String, Surface, Vec<[f64; 3]>)> {
    let mut v = Vec::new();
    for (n, s, p) in all_meshes() {
        v.push((n.to_string(), s, p));
    }
    v
}

#[test]
fn obstacle_energy_at_theta_zero_is_the_exact_double_sum() {
    let pool = make_pool(Some(3)).unwrap();
    let e = Exponents::default();
    let ob = shifted_ico(1, 0.8, [2.5, 0.4, 0.3]);
    let normals = triangle_normals(&ob);
    let tri = Arc::new(StaticSet::from_obstacle(&ob, 3).unwrap());
    let cloud_pts: Vec<[f64; 3]> = ob.vertices().to_vec();
    let cloud = Arc::new(StaticSet::from_points(&cloud_pts, 4).unwrap());
    for (name, s, pos) in cases() {
        for leaf in [1, 4, 8] {
            for (label, set, nrm, kind) in [
                ("tp", &tri, Some(normals.as_slice()), ObstacleKind::TangentPoint),
                ("coulomb", &tri, None, ObstacleKind::Coulomb { power: 3.0 }),
                ("cloud", &cloud, None, ObstacleKind::TangentPoint),
            ] {
                let ev = CrossEvaluator::build(&pool, &s, &pos, set.clone(), kind, &e, bh(0.0, leaf)).unwrap();
                let want = brute(&s, &pos, set, nrm, kind, &e);
                let got = ev.energy(&pool);
                let (g2, grad) = ev.energy_and_differential(&pool);
                assert!((got / want - 1.0).abs() < 1e-12, "{name} {label} leaf {leaf}: {got} {want}");
                assert!((g2 / want - 1.0).abs() < 1e-12);
                assert!(grad.iter().flatten().all(|x| x.is_finite()));
                assert_eq!(ev.n_lumped(), 0);
            }
        }
    }
}

#[test]
fn obstacle_theta_error_decreases() {
    let pool = make_pool(Some(3)).unwrap();
    let e = Exponents::default();
    let ob = shifted_ico(2, 1.0, [3.0, 0.0, 0.0]);
    let tri = Arc::new(StaticSet::from_obstacle(&ob, 8).unwrap());
    let (_, s, pos) = closed_meshes().remove(0);
    let exact = CrossEvaluator::build(&pool, &s, &pos, tri.clone(), ObstacleKind::TangentPoint, &e, bh(0.0, 8))
        .unwrap()
        .energy(&pool);
    let mut prev = f64::INFINITY;
    for theta in [1.0, 0.5, 0.25, 0.125] {
        let ev = CrossEvaluator::build(&pool, &s, &pos, tri.clone(), ObstacleKind::TangentPoint, &e, bh(theta, 8)).unwrap();
        let err = (ev.energy(&pool) / exact - 1.0).abs();
        println!("theta {theta}: rel err {err:.2e}, lumped {}", ev.n_lumped());
        assert!(err < prev * 1.01 || err < 1e-12);
        prev = err;
    }
    assert!(prev < 1e-3);
}

#[test]
fn obstacle_gradient_matches_finite_differences_on_frozen_lists() {
    let pool = make_pool(Some(3)).unwrap();
    let e = Exponents::default();
    let ob = shifted_ico(2, 0.9, [2.6, 0.3, 0.2]);
    let tri = Arc::new(StaticSet::from_obstacle(&ob, 4).unwrap());
    let cloud = Arc::new(StaticSet::from_points(ob.vertices(), 4).unwrap());
    let mut worst: f64 = 0.0;
    for (name, s, pos) in cases() {
        for theta in [0.5, 1.0] {
            for (label, set, kind) in [
                ("tp", &tri, ObstacleKind::TangentPoint),
                ("coulomb", &tri, ObstacleKind::Coulomb { power: 4.0 }),
                ("cloud", &cloud, ObstacleKind::TangentPoint),
            ] {
                let mut ev = CrossEvaluator::build(&pool, &s, &pos, set.clone(), kind, &e, bh(theta, 4)).unwrap();
                assert!(ev.n_lumped() > 0 || theta == 0.0, "{name} {label}: nothing lumped");
                let (en, g) = ev.energy_and_differential(&pool);
                assert!((ev.energy(&pool) / en - 1.0).abs() < 1e-13);
                let h = 1e-6 * mean_edge_length(&s, &pos);
                let stride = (pos.len() / 40).max(1);
                let (mut err, mut mg) = (0.0f64, 0.0f64);
                for v in (0..pos.len()).step_by(stride) {
                    for c in 0..3 {
                        let mut q = pos.clone();
                        q[v][c] += h;
                        ev.refit(&q).unwrap();
                        let fp = ev.energy(&pool);
                        q[v][c] -= 2.0 * h;
                        ev.refit(&q).unwrap();
                        let fm = ev.energy(&pool);
                        err = err.max(((fp - fm) / (2.0 * h) - g[v][c]).abs());
                        mg = mg.max(g[v][c].abs());
                    }
                }
                ev.refit(&pos).unwrap();
                let rel = err / mg;
                worst = worst.max(rel);
                println!("{name} {label} theta {theta}: fd rel {rel:.2e}");
                assert!(rel < 1e-6, "{name} {label} theta {theta}: {rel:e}");
            }
        }
    }
    println!("worst obstacle fd error {worst:.2e}");
}

#[test]
fn obstacle_penalties_through_energy() {
    let e = Exponents::default();
    let ob = shifted_ico(1, 0.8, [2.5, 0.4, 0.3]);
    let tri = Arc::new(StaticSet::from_obstacle(&ob, 4).unwrap());
    let cloud = Arc::new(StaticSet::from_points(ob.vertices(), 4).unwrap());
    let (_, s, pos) = closed_meshes().remove(0);
    for pen in [
        Penalty::MeshObstacle { weight: 0.7, obstacle: tri.clone(), kind: ObstacleKind::TangentPoint },
        Penalty::MeshObstacle { weight: 0.7, obstacle: tri.clone(), kind: ObstacleKind::Coulomb { power: 2.0 } },
        Penalty::PointCloudObstacle { weight: 0.7, points: cloud.clone() },
    ] {
        // theta = 0 (exact evaluation): the penalty gradient is the exact one.
        assert_fd(pen_name(&pen), &s, &pos, &pen);
    }
    // Wiring: energy() adds the weighted differential, reports the weighted value.
    let pen = Penalty::MeshObstacle { weight: 0.7, obstacle: tri.clone(), kind: ObstacleKind::TangentPoint };
    let base = energy(&s, &pos, e, &EnergyEval::Exact, &[], Some(2)).unwrap();
    let with = energy(&s, &pos, e, &EnergyEval::Exact, std::slice::from_ref(&pen), Some(2)).unwrap();
    let pool = make_pool(Some(2)).unwrap();
    let (v, g) = evaluate_one(&pool, &s, &pos, &pen, &e, &EnergyEval::Exact).unwrap();
    assert_eq!(with.energy, base.energy);
    assert_eq!(with.penalties, vec![v]);
    for i in 0..pos.len() {
        for c in 0..3 {
            assert!((with.differential[i][c] - base.differential[i][c] - g[i][c]).abs() < 1e-12 * (1.0 + g[i][c].abs()));
        }
    }
    // Errors.
    let r = evaluate_one(&pool, &s, &pos, &Penalty::MeshObstacle { weight: 1.0, obstacle: cloud.clone(), kind: ObstacleKind::TangentPoint }, &e, &EnergyEval::Exact);
    assert!(matches!(r, Err(RsError::ObstacleError { .. })));
    let r = evaluate_one(&pool, &s, &pos, &Penalty::PointCloudObstacle { weight: 1.0, points: tri.clone() }, &e, &EnergyEval::Exact);
    assert!(matches!(r, Err(RsError::ObstacleError { .. })));
    assert!(StaticSet::from_points(&[], 4).is_err());
}

fn pen_name(p: &Penalty) -> &'static str {
    p.name()
}

#[test]
fn obstacle_diagonal_matches_the_brute_force_formula() {
    let pool = make_pool(Some(3)).unwrap();
    let e = Exponents::default();
    let ob = shifted_ico(1, 0.8, [2.2, 0.4, 0.3]);
    let normals = triangle_normals(&ob);
    let set = Arc::new(StaticSet::from_obstacle(&ob, 3).unwrap());
    let (_, s, pos) = closed_meshes().remove(0);
    let g = FaceGeometry::new(&s, &pos).unwrap();
    let w = 1.7;
    let sx = e.s();
    let d = obstacle_diagonal(&pool, &s, &pos, set.clone(), &e, bh(0.0, 4), w).unwrap();
    for f in 0..s.n_faces() {
        let (mut hi, mut lo) = (0.0, 0.0);
        for t in 0..set.len() {
            let v = [
                g.barycenters[f][0] - set.points()[t][0],
                g.barycenters[f][1] - set.points()[t][1],
                g.barycenters[f][2] - set.points()[t][2],
            ];
            let r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            let dn = |n: [f64; 3]| (n[0] * v[0] + n[1] * v[1] + n[2] * v[2]).powi(2);
            let rs = r2.powf(-sx);
            hi += set.weights()[t] * rs;
            lo += set.weights()[t] * rs * 0.5 * (dn(g.normals[f]) + dn(normals[t])) / (r2 * r2);
        }
        let c = 2.0 * w / g.areas[f];
        assert!((d.high[f] / (c * hi) - 1.0).abs() < 1e-12);
        assert!((d.low[f] / (c * lo) - 1.0).abs() < 1e-12);
    }
    // Barnes-Hut version is close and positive.
    let b = obstacle_diagonal(&pool, &s, &pos, set.clone(), &e, bh(0.5, 4), w).unwrap();
    for f in 0..s.n_faces() {
        assert!(b.high[f] > 0.0 && b.low[f] >= 0.0);
        assert!((b.high[f] / d.high[f] - 1.0).abs() < 0.2);
    }
    // A point cloud has the metric term too (the reference builds the obstacle BCT for it): unit
    // weights and zero obstacle projectors, so `K₂(T, S) = 0`.
    let cloud = Arc::new(StaticSet::from_points(ob.vertices(), 4).unwrap());
    for theta in [0.0, 0.5] {
        let d = obstacle_diagonal(&pool, &s, &pos, cloud.clone(), &e, bh(theta, 4), w).unwrap();
        for f in 0..s.n_faces() {
            let (mut hi, mut lo) = (0.0, 0.0);
            for t in 0..cloud.len() {
                let p = cloud.points()[t];
                let v = [
                    g.barycenters[f][0] - p[0],
                    g.barycenters[f][1] - p[1],
                    g.barycenters[f][2] - p[2],
                ];
                let r2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                let nv = g.normals[f];
                let dn = (nv[0] * v[0] + nv[1] * v[1] + nv[2] * v[2]).powi(2);
                let rs = r2.powf(-sx);
                hi += rs;
                lo += rs * 0.5 * dn / (r2 * r2);
            }
            let c = 2.0 * w / g.areas[f];
            let tol = if theta == 0.0 { 1e-12 } else { 0.2 };
            assert!(d.high[f].is_finite() && d.low[f].is_finite());
            assert!((d.high[f] / (c * hi) - 1.0).abs() < tol, "{theta} {f}");
            assert!((d.low[f] / (c * lo) - 1.0).abs() < tol.max(0.3) || theta == 0.0, "{theta} {f}");
        }
    }
}

#[test]
fn coulomb_energy_matches_the_double_sum_and_finite_differences() {
    let pool = make_pool(Some(3)).unwrap();
    for (name, s, pos) in all_meshes() {
        for alpha in [1.0, 2.0, 3.5] {
            let (en, g) = coulomb_energy(&pool, &s, &pos, alpha).unwrap();
            let gm = FaceGeometry::new(&s, &pos).unwrap();
            let a = vertex_areas(&s, &gm.areas);
            let mut want = 0.0;
            for v in 0..pos.len() {
                for w in 0..pos.len() {
                    if v != w {
                        let d = (0..3).map(|c| (pos[v][c] - pos[w][c]).powi(2)).sum::<f64>().sqrt();
                        want += a[v] * a[w] * d.powf(-alpha);
                    }
                }
            }
            assert!((en / want - 1.0).abs() < 1e-12, "{name} {alpha}");
            let h = 1e-6 * mean_edge_length(&s, &pos);
            let stride = (pos.len() / 30).max(1);
            let (mut err, mut mg) = (0.0f64, 0.0f64);
            for v in (0..pos.len()).step_by(stride) {
                for c in 0..3 {
                    let mut q = pos.clone();
                    q[v][c] += h;
                    let fp = coulomb_energy(&pool, &s, &q, alpha).unwrap().0;
                    q[v][c] -= 2.0 * h;
                    let fm = coulomb_energy(&pool, &s, &q, alpha).unwrap().0;
                    err = err.max(((fp - fm) / (2.0 * h) - g[v][c]).abs());
                    mg = mg.max(g[v][c].abs());
                }
            }
            println!("{name} coulomb {alpha}: fd rel {:.2e}", err / mg);
            assert!(err / mg < 1e-6, "{name} {alpha}: {:e}", err / mg);
        }
    }
    let (_, s, pos) = closed_meshes().remove(0);
    for bad in [0.0, -1.0, f64::NAN] {
        assert!(matches!(coulomb_energy(&pool, &s, &pos, bad), Err(RsError::InvalidOption { .. })));
    }
}

#[test]
fn results_do_not_depend_on_the_thread_count() {
    let e = Exponents::default();
    let ob = shifted_ico(2, 0.9, [2.6, 0.3, 0.2]);
    let tri = Arc::new(StaticSet::from_obstacle(&ob, 4).unwrap());
    let (_, s, pos) = closed_meshes().remove(1);
    let run = |threads: usize| {
        let pool = make_pool(Some(threads)).unwrap();
        let ev = CrossEvaluator::build(&pool, &s, &pos, tri.clone(), ObstacleKind::TangentPoint, &e, bh(0.5, 4)).unwrap();
        let (en, g) = ev.energy_and_differential(&pool);
        let d = ev.metric_diagonal(&pool, e.s(), 1.0).unwrap();
        let c = coulomb_energy(&pool, &s, &pos, 2.0).unwrap();
        (en, g, d, c)
    };
    let r1 = run(1);
    for t in [3, 8] {
        assert_eq!(run(t), r1, "threads {t}");
    }
}
