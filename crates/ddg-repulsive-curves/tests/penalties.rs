//! Penalties, obstacles, implicit surfaces and vector fields (`SPEC.md` §5.3, §13.2, §13.5;
//! work package WP7): gradients against central differences, the surface potential with the BVH
//! against brute force, the barrier and the guarded flow, field alignment, length difference,
//! surfaces and the surface constraint in short flows, subdivision of penalties, OBJ obstacles
//! and their errors, determinism, and (ignored, release) the cost of the surface potential.

#![allow(clippy::needless_range_loop)]

use std::f64::consts::{PI, TAU};
use std::sync::Arc;

use ddg_repulsive_curves::ccd::{CcdOptions, PreparedObstacle, first_contact};
use ddg_repulsive_curves::obstacle::{Obstacle, closest_point_on_triangle};
use ddg_repulsive_curves::penalties::{
    FieldAlignment, LengthDifference, SurfacePotential, TotalLengthPenalty, calibrate_hencky_weight,
};
use ddg_repulsive_curves::step::stop::Status;
use ddg_repulsive_curves::step::{FlowOptions, PenaltyTerm, Session};
use ddg_repulsive_curves::subdivide::{SubdivisionMap, subdivide};
use ddg_repulsive_curves::surfaces::{
    CircularField, ConstantField, DoubleTorus, Metaballs, Plane, RoundedBox, SmoothUnion, Sphere,
    SphereTangentField, Torus, Union, VectorField,
};
use ddg_repulsive_curves::{
    CollisionObstacle, Constraint, CurveNetwork, ImplicitSurface, LengthSchedule, RcError,
};

type P = [f64; 3];

// ---------------------------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------------------------

fn sub3(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm3(a: P) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn loop_edges(n: usize) -> Vec<[u32; 2]> {
    (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect()
}

fn chain_edges(n_vertices: usize) -> Vec<[u32; 2]> {
    (0..n_vertices as u32 - 1).map(|i| [i, i + 1]).collect()
}

/// A wobbly, irregularly sampled closed curve in R³.
fn wobbly_loop(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|i| {
            let t = TAU * (i as f64 + 0.3 * (3.0 * i as f64).sin()) / n as f64;
            [
                t.cos() + 0.2 * (3.0 * t).cos(),
                t.sin() + 0.1 * (2.0 * t).sin(),
                0.3 * (2.0 * t).sin() + 0.1 * (5.0 * t).cos(),
            ]
        })
        .collect();
    (CurveNetwork::new(n, loop_edges(n)).unwrap(), pos)
}

/// A "Y": three arms meeting at a degree-3 vertex, one arm with a bend, positions generic.
fn y_network() -> (CurveNetwork, Vec<P>) {
    let edges = vec![[0, 1], [1, 2], [2, 3], [3, 4], [3, 5], [5, 6], [6, 7]];
    let pos = vec![
        [0.0, 0.0, 0.0],
        [0.5, 0.1, 0.2],
        [1.1, -0.1, 0.1],
        [1.6, 0.3, 0.0],
        [2.0, 0.9, 0.3],
        [2.2, -0.2, -0.1],
        [2.9, -0.4, 0.2],
        [3.3, -0.1, 0.5],
    ];
    (CurveNetwork::new(8, edges).unwrap(), pos)
}

/// Central differences of the value against the gradient (which must be *added* to the array).
fn check_gradient(pen: &dyn PenaltyTerm, net: &CurveNetwork, pos: &[P], tol: f64) {
    let mut g = vec![[0.5, -0.25, 0.125]; pos.len()];
    let v0 = pen.evaluate(net, pos, Some(&mut g)).unwrap();
    assert_eq!(v0, pen.evaluate(net, pos, None).unwrap());
    for gi in &mut g {
        *gi = sub3(*gi, [0.5, -0.25, 0.125]);
    }
    let scale = g.iter().map(|x| norm3(*x)).fold(0.0, f64::max).max(1.0);
    let mut worst = 0.0f64;
    for v in 0..pos.len() {
        for k in 0..3 {
            let h = 1e-6 * (1.0 + pos[v][k].abs());
            let mut p = pos.to_vec();
            p[v][k] += h;
            let fp = pen.evaluate(net, &p, None).unwrap();
            p[v][k] -= 2.0 * h;
            let fm = pen.evaluate(net, &p, None).unwrap();
            let fd = (fp - fm) / (2.0 * h);
            worst = worst.max((fd - g[v][k]).abs() / scale);
        }
    }
    assert!(
        worst < tol,
        "{}: gradient vs central differences: {worst:.2e} (tolerance {tol:.0e})",
        pen.name()
    );
}

/// A UV sphere as an obstacle (outward normals irrelevant).
fn uv_sphere(center: P, radius: f64, n_lat: usize, n_lon: usize) -> Obstacle {
    let mut v = vec![[center[0], center[1], center[2] + radius]];
    for i in 1..n_lat {
        let th = PI * i as f64 / n_lat as f64;
        for j in 0..n_lon {
            let ph = TAU * j as f64 / n_lon as f64;
            v.push([
                center[0] + radius * th.sin() * ph.cos(),
                center[1] + radius * th.sin() * ph.sin(),
                center[2] + radius * th.cos(),
            ]);
        }
    }
    v.push([center[0], center[1], center[2] - radius]);
    let south = (v.len() - 1) as u32;
    let ring = |i: usize, j: usize| (1 + (i - 1) * n_lon + j % n_lon) as u32;
    let mut t = Vec::new();
    for j in 0..n_lon {
        t.push([0, ring(1, j), ring(1, j + 1)]);
    }
    for i in 1..n_lat - 1 {
        for j in 0..n_lon {
            t.push([ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)]);
            t.push([ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)]);
        }
    }
    for j in 0..n_lon {
        t.push([ring(n_lat - 1, j), south, ring(n_lat - 1, j + 1)]);
    }
    Obstacle::new(v, t).unwrap()
}

/// The big horizontal quad `[-2, 3] × [-2, 2]` at height `z` as two triangles.
fn quad_at(z: f64) -> Obstacle {
    Obstacle::new(
        vec![
            [-2.0, -2.0, z],
            [3.0, -2.0, z],
            [3.0, 2.0, z],
            [-2.0, 2.0, z],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap()
}

/// A tiny deterministic generator (SplitMix64) for the sampled tests.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }
}

// ---------------------------------------------------------------------------------------------
// Gradients against central differences (every penalty)
// ---------------------------------------------------------------------------------------------

#[test]
fn total_length_gradient_and_value() {
    let (net, pos) = wobbly_loop(23);
    for w in [1.0, -0.7] {
        let pen = TotalLengthPenalty::new(w).unwrap();
        let geo = net.edge_geometry(&pos);
        let total: f64 = geo.length.iter().sum();
        let v = pen.evaluate(&net, &pos, None).unwrap();
        assert!((v - w * total).abs() < 1e-12);
        check_gradient(&pen, &net, &pos, 1e-8);
    }
    let (ny, py) = y_network();
    check_gradient(&TotalLengthPenalty::new(2.0).unwrap(), &ny, &py, 1e-8);
}

#[test]
fn length_difference_gradient_and_value() {
    let (net, pos) = wobbly_loop(23);
    let geo = net.edge_geometry(&pos);
    let pen = LengthDifference::new(3.0).unwrap();
    // On a loop every vertex has degree 2: w Σ_v (ℓ_I − ℓ_J)^2 over adjacent pairs.
    let expect: f64 = (0..23)
        .map(|i| (geo.length[i] - geo.length[(i + 1) % 23]).powi(2))
        .sum();
    assert!((pen.evaluate(&net, &pos, None).unwrap() - 3.0 * expect).abs() < 1e-12);
    check_gradient(&pen, &net, &pos, 1e-7);
    check_gradient(&LengthDifference::new(-2.0).unwrap(), &net, &pos, 1e-7);
    // Degree-1 and degree-3 vertices contribute nothing.
    let (ny, py) = y_network();
    check_gradient(&pen, &ny, &py, 1e-7);
    let gy = ny.edge_geometry(&py);
    let mut e = 0.0;
    for v in 0..ny.n_vertices() as u32 {
        let inc = ny.incident_edges(v);
        if inc.len() == 2 {
            e += (gy.length[inc[0] as usize] - gy.length[inc[1] as usize]).powi(2);
        }
    }
    assert!((pen.evaluate(&ny, &py, None).unwrap() - 3.0 * e).abs() < 1e-12);
}

#[test]
fn field_alignment_value_and_gradient() {
    let (net, pos) = wobbly_loop(19);
    let geo = net.edge_geometry(&pos);
    let fields: Vec<(&str, Arc<dyn VectorField>)> = vec![
        (
            "constant",
            Arc::new(ConstantField::new([0.3, -0.5, 0.8]).unwrap()),
        ),
        ("circular", Arc::new(CircularField)),
        ("sphere", Arc::new(SphereTangentField)),
    ];
    for (name, f) in &fields {
        let pen = FieldAlignment::new(f.clone(), 2.5).unwrap();
        // value = w Σ ℓ |T × X(x)|²
        let mut e = 0.0;
        for i in 0..net.n_edges() {
            let x = f.value(geo.midpoint[i]);
            let t = geo.tangent[i];
            let c = [
                t[1] * x[2] - t[2] * x[1],
                t[2] * x[0] - t[0] * x[2],
                t[0] * x[1] - t[1] * x[0],
            ];
            e += geo.length[i] * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        }
        let v = pen.evaluate(&net, &pos, None).unwrap();
        assert!((v - 2.5 * e).abs() < 1e-12 * (1.0 + e), "{name}");
        check_gradient(&pen, &net, &pos, 1e-6);
        let (ny, py) = y_network();
        check_gradient(&pen, &ny, &py, 1e-6);
    }
}

#[test]
fn vector_field_jacobians_match_finite_differences() {
    let fields: Vec<(&str, Box<dyn VectorField>)> = vec![
        (
            "constant",
            Box::new(ConstantField::unit([1.0, 2.0, -1.0]).unwrap()),
        ),
        ("circular", Box::new(CircularField)),
        ("sphere", Box::new(SphereTangentField)),
    ];
    let mut rng = Rng(7);
    for (name, f) in &fields {
        for _ in 0..20 {
            let x = [
                rng.range(-1.5, 1.5),
                rng.range(-1.5, 1.5),
                rng.range(-1.5, 1.5),
            ];
            let j = f.jacobian(x);
            for k in 0..3 {
                let h = 1e-6;
                let (mut a, mut b) = (x, x);
                a[k] += h;
                b[k] -= h;
                let (fa, fb) = (f.value(a), f.value(b));
                for i in 0..3 {
                    let fd = (fa[i] - fb[i]) / (2.0 * h);
                    assert!(
                        (fd - j[i][k]).abs() < 1e-6,
                        "{name}: J[{i}][{k}] {fd} vs {}",
                        j[i][k]
                    );
                }
            }
        }
    }
    // The tangent field is a unit field tangent to spheres about the origin.
    let x = [0.3, 0.4, -0.8];
    let v = SphereTangentField.value(x);
    assert!((norm3(v) - 1.0).abs() < 1e-12);
    assert!((v[0] * x[0] + v[1] * x[1] + v[2] * x[2]).abs() < 1e-12);
}

#[test]
fn surface_potential_gradient_with_and_without_lumping() {
    let obs = Arc::new(uv_sphere([0.0, 0.0, 0.0], 1.0, 12, 20));
    // a curve outside the sphere, at various distances
    let n = 17;
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let t = TAU * i as f64 / n as f64;
            let r = 1.4 + 0.5 * (3.0 * t).sin().abs();
            [r * t.cos(), r * t.sin(), 0.6 * (2.0 * t).sin()]
        })
        .collect();
    let net = CurveNetwork::new(n, loop_edges(n)).unwrap();
    for theta in [0.0, 0.25, 0.5] {
        let pen = SurfacePotential::new(obs.clone(), 0.7, 3.0, theta).unwrap();
        // the lumped function is differentiated exactly (the opening test does not switch for
        // perturbations of 1e-6 at these points)
        check_gradient(&pen, &net, &pos, 1e-6);
    }
    // inside the sphere as well (the soup does not need to be closed or oriented)
    let inner: Vec<P> = pos
        .iter()
        .map(|p| [0.4 * p[0], 0.4 * p[1], 0.4 * p[2]])
        .collect();
    check_gradient(
        &SurfacePotential::new(obs, 1.0, 2.5, 0.0).unwrap(),
        &net,
        &inner,
        1e-6,
    );
}

#[test]
fn hencky_gradient_and_calibration_rule() {
    use ddg_repulsive_curves::penalties::hencky::Hencky;
    let (net, pos) = wobbly_loop(15);
    let h = Hencky::at_rest(2.0, &net, &pos).unwrap();
    let moved: Vec<P> = pos.iter().map(|p| [p[0] * 1.1, p[1], p[2] * 0.9]).collect();
    check_gradient(&h, &net, &moved, 1e-7);
    // smallest 10^k whose median strain is <= 1e-2
    let w = calibrate_hencky_weight(-3..=6, 1e-2, |w| Ok(1.0 / w.sqrt())).unwrap();
    assert_eq!(w, Some(1e4));
    assert_eq!(
        calibrate_hencky_weight(0..=2, 1e-9, |_| Ok(1.0)).unwrap(),
        None
    );
    assert!(calibrate_hencky_weight(0..=2, -1.0, |_| Ok(1.0)).is_err());
}

#[test]
fn penalty_options_are_validated() {
    let obs = Arc::new(quad_at(0.0));
    assert!(TotalLengthPenalty::new(f64::NAN).is_err());
    assert!(LengthDifference::new(f64::INFINITY).is_err());
    assert!(FieldAlignment::new(Arc::new(CircularField), -1.0).is_err());
    assert!(SurfacePotential::new(obs.clone(), -1.0, 3.0, 0.25).is_err());
    assert!(SurfacePotential::new(obs.clone(), 1.0, 0.0, 0.25).is_err());
    assert!(SurfacePotential::new(obs.clone(), 1.0, 3.0, -0.1).is_err());
    assert!(SurfacePotential::new(obs.clone(), 1.0, 3.0, f64::NAN).is_err());
    // wrong array sizes are errors, not panics
    let (net, pos) = wobbly_loop(8);
    let pen = TotalLengthPenalty::new(1.0).unwrap();
    assert!(matches!(
        pen.evaluate(&net, &pos[..5], None),
        Err(RcError::PositionCountMismatch { .. })
    ));
    let mut g = vec![[0.0; 3]; 3];
    assert!(matches!(
        pen.evaluate(&net, &pos, Some(&mut g)),
        Err(RcError::PositionCountMismatch { .. })
    ));
    // a vertex on the obstacle: infinite value, gradient request is an error
    let on = vec![
        [-2.0, -2.0, 0.0],
        [0.5, 0.0, 0.3],
        [0.0, 0.5, 0.4],
        [0.5, 0.5, 0.5],
    ];
    let n4 = CurveNetwork::new(4, loop_edges(4)).unwrap();
    let sp = SurfacePotential::new(obs, 1.0, 3.0, 0.25).unwrap();
    assert_eq!(sp.evaluate(&n4, &on, None).unwrap(), f64::INFINITY);
    let mut g4 = vec![[0.0; 3]; 4];
    assert!(sp.evaluate(&n4, &on, Some(&mut g4)).is_err());
}

// ---------------------------------------------------------------------------------------------
// Closest points, distance, BVH against brute force
// ---------------------------------------------------------------------------------------------

#[test]
fn closest_point_on_triangle_is_the_minimizer() {
    let mut rng = Rng(11);
    for _ in 0..2000 {
        let mut pt = || {
            [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ]
        };
        let (a, b, c, p) = (pt(), pt(), pt(), pt());
        let q = closest_point_on_triangle(p, a, b, c);
        let d = norm3(sub3(p, q));
        // optimality: (p - q)·(r - q) <= 0 for r in the triangle (corners suffice by convexity),
        // and q is in the triangle (barycentric)
        for r in [a, b, c] {
            let rq = sub3(r, q);
            let pq = sub3(p, q);
            assert!(pq[0] * rq[0] + pq[1] * rq[1] + pq[2] * rq[2] <= 1e-12);
        }
        // sampled points of the triangle are never closer
        for i in 0..=10 {
            for j in 0..=10 - i {
                let (u, v) = (i as f64 / 10.0, j as f64 / 10.0);
                let s = [
                    a[0] + u * (b[0] - a[0]) + v * (c[0] - a[0]),
                    a[1] + u * (b[1] - a[1]) + v * (c[1] - a[1]),
                    a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2]),
                ];
                assert!(norm3(sub3(p, s)) >= d - 1e-12);
            }
        }
        // consistent with the CCD's distance
        let d2 = ddg_repulsive_curves::ccd::point_triangle_distance(p, a, b, c);
        assert!((d - d2).abs() < 1e-10, "{d} vs {d2}");
    }
}

#[test]
fn distance_gradient_is_the_unit_direction_to_the_closest_point() {
    // [formula: rc-distance-gradient]
    let obs = uv_sphere([0.0; 3], 1.0, 8, 12);
    let mut rng = Rng(3);
    for _ in 0..100 {
        let x = [
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
        ];
        let (d, t) = obs.distance_to(x);
        if d < 0.05 {
            continue;
        }
        let t = t as usize;
        let (dt, q) = obs.triangle_distance(x, t);
        assert_eq!(d, dt);
        let n = [(x[0] - q[0]) / d, (x[1] - q[1]) / d, (x[2] - q[2]) / d];
        for k in 0..3 {
            let h = 1e-6;
            let (mut a, mut b) = (x, x);
            a[k] += h;
            b[k] -= h;
            let fd = (obs.triangle_distance(a, t).0 - obs.triangle_distance(b, t).0) / (2.0 * h);
            assert!((fd - n[k]).abs() < 1e-6);
        }
    }
}

#[test]
fn distance_query_agrees_with_brute_force() {
    let obs = uv_sphere([0.2, -0.1, 0.3], 1.0, 30, 50);
    let mut rng = Rng(5);
    for _ in 0..300 {
        let x = [
            rng.range(-3.0, 3.0),
            rng.range(-3.0, 3.0),
            rng.range(-3.0, 3.0),
        ];
        let (d, t) = obs.distance_to(x);
        let b = obs.distance_to_brute(x);
        assert_eq!(d, b);
        assert_eq!(obs.triangle_distance(x, t as usize).0, b);
    }
}

#[test]
fn surface_potential_bvh_against_brute_force() {
    let obs = uv_sphere([0.0; 3], 1.0, 40, 60); // 4680 triangles
    assert!(obs.n_triangles() > 4000);
    let mut rng = Rng(9);
    let p = 3.0;
    let thetas = [0.0, 0.05, 0.1, 0.25, 0.5];
    let mut max_err = vec![0.0f64; thetas.len()];
    let mut max_gerr = vec![0.0f64; thetas.len()];
    for _ in 0..400 {
        // points outside, log-uniformly in the distance to the sphere
        let dist = 10f64.powf(rng.range(-2.0, 0.8));
        let dir = loop {
            let d = [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ];
            let n = norm3(d);
            if n > 0.1 && n < 1.0 {
                break [d[0] / n, d[1] / n, d[2] / n];
            }
        };
        let x = [
            dir[0] * (1.0 + dist),
            dir[1] * (1.0 + dist),
            dir[2] * (1.0 + dist),
        ];
        let exact = obs.potential_at_brute(x, p, true);
        for (i, &theta) in thetas.iter().enumerate() {
            let q = obs.potential_at(x, p, theta, true);
            let rel = (q.value - exact.value).abs() / exact.value;
            let grel = norm3(sub3(q.gradient, exact.gradient)) / norm3(exact.gradient);
            max_err[i] = max_err[i].max(rel);
            max_gerr[i] = max_gerr[i].max(grel);
            if theta == 0.0 {
                assert_eq!(q.far_nodes, 0);
                assert!(
                    rel < 1e-12 && grel < 1e-11,
                    "theta 0 is exact: {rel:e} {grel:e}"
                );
            } else {
                // Every lumped node has r_N <= θ|x − c_N| and all its triangles lie in its box,
                // so |δ − |x − c_N|| <= 2 r_N <= 2 θ |x − c_N|, and each lumped term is off by
                // a factor in [(1 + 2θ)^{−p}, (1 − 2θ)^{−p}]: a rigorous relative bound.
                if theta < 0.5 {
                    let bound = (1.0 - 2.0 * theta).powf(-p) - 1.0;
                    assert!(rel <= bound, "theta {theta}: {rel:e} > {bound:e}");
                }
            }
        }
    }
    eprintln!("potential error vs brute force (p = 3, 400 points, 4680 triangles):");
    for (i, th) in thetas.iter().enumerate() {
        eprintln!(
            "  theta {th:>5}: value {:.2e}, gradient {:.2e}",
            max_err[i], max_gerr[i]
        );
    }
    // the error decreases with theta; the default opening 0.1 is within 10 %
    assert!(max_err[1] < max_err[2] && max_err[2] < max_err[3] && max_err[3] < max_err[4]);
    assert!(max_err[2] < 0.1, "{}", max_err[2]);
    assert!(max_gerr[2] < 0.15, "{}", max_gerr[2]);
}

#[test]
fn surface_potential_near_field_is_exact() {
    // Within one triangle size of the soup the nearest triangles are exact, so the lumped
    // evaluation agrees with brute force to far better than the far-field bound.
    let obs = uv_sphere([0.0; 3], 1.0, 40, 60);
    let mut rng = Rng(21);
    for _ in 0..100 {
        let d = [
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        ];
        let n = norm3(d);
        if n < 0.2 {
            continue;
        }
        let x = [d[0] / n * 1.003, d[1] / n * 1.003, d[2] / n * 1.003];
        let a = obs.potential_at(x, 3.0, 0.25, false);
        let b = obs.potential_at_brute(x, 3.0, false);
        assert!(a.near_triangles > 0 && a.far_nodes > 0);
        assert!(
            (a.value - b.value).abs() / b.value < 2e-3,
            "{}",
            (a.value - b.value) / b.value
        );
    }
}

#[test]
fn surface_potential_is_deterministic_across_thread_counts() {
    let obs = Arc::new(uv_sphere([0.0; 3], 1.0, 20, 30));
    let (net, base) = wobbly_loop(200);
    let pos: Vec<P> = base
        .iter()
        .map(|p| [p[0] * 1.5, p[1] * 1.5, p[2] * 1.5 + 0.1])
        .collect();
    let pen = SurfacePotential::new(obs, 1.3, 3.0, 0.25).unwrap();
    let mut reference: Option<(u64, Vec<[u64; 3]>)> = None;
    for threads in [1usize, 7, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let (v, g) = pool.install(|| {
            let mut g = vec![[0.0; 3]; pos.len()];
            let v = pen.evaluate(&net, &pos, Some(&mut g)).unwrap();
            (v, g)
        });
        let bits = (
            v.to_bits(),
            g.iter().map(|x| x.map(f64::to_bits)).collect::<Vec<_>>(),
        );
        match &reference {
            None => reference = Some(bits),
            Some(r) => assert_eq!(r, &bits, "{threads} threads"),
        }
    }
}

// ---------------------------------------------------------------------------------------------
// The barrier and the guarded flow
// ---------------------------------------------------------------------------------------------

#[test]
fn potential_blows_up_at_the_interior_of_a_face() {
    // A coarse obstacle (two huge triangles): the barrier is on the face, whose vertices are far
    // away from the curve (D19; sampling the mesh vertices would give a bounded potential).
    let obs = Arc::new(quad_at(0.0));
    let net = CurveNetwork::new(4, loop_edges(4)).unwrap();
    let p = 3.0;
    let pen = SurfacePotential::new(obs.clone(), 1.0, p, 0.25).unwrap();
    let at = |h: f64| -> f64 {
        let pos = vec![[0.4, 0.3, h], [0.6, 0.3, h], [0.6, 0.5, h], [0.4, 0.5, h]];
        pen.evaluate(&net, &pos, None).unwrap()
    };
    let mut prev = at(1e-1);
    for k in 2..=5 {
        let v = at(10f64.powi(-k));
        // ~ 10^p per decade (the face is an almost infinite plane)
        let ratio = v / prev;
        assert!(
            ratio > 0.5 * 10f64.powf(p) && ratio < 2.0 * 10f64.powf(p),
            "{k}: {ratio}"
        );
        prev = v;
    }
    assert!(at(1e-5).is_finite() && at(1e-5) > 1e9);
    // symmetric from below: the obstacle is a soup without orientation
    let below = vec![
        [0.4, 0.3, -1e-3],
        [0.6, 0.3, -1e-3],
        [0.6, 0.5, -1e-3],
        [0.4, 0.5, -1e-3],
    ];
    assert!((pen.evaluate(&net, &below, None).unwrap() - at(1e-3)).abs() < 1e-9 * at(1e-3));
}

/// Gravity for the flow tests: `w Σ_I ℓ_I (z_a + z_b)/2` (the potential of a uniform field on a
/// chain of uniform line density), exact gradient.
#[derive(Debug)]
struct Gravity(f64);
impl PenaltyTerm for Gravity {
    fn name(&self) -> &str {
        "Gravity"
    }
    fn evaluate(
        &self,
        network: &CurveNetwork,
        positions: &[P],
        mut gradient: Option<&mut [P]>,
    ) -> Result<f64, RcError> {
        let mut total = 0.0;
        for &[a, b] in network.edges() {
            let (pa, pb) = (positions[a as usize], positions[b as usize]);
            let d = sub3(pb, pa);
            let l = norm3(d);
            let zm = 0.5 * (pa[2] + pb[2]);
            total += l * zm;
            if let Some(g) = gradient.as_deref_mut() {
                for k in 0..3 {
                    let t = d[k] / l * zm;
                    g[b as usize][k] += self.0 * t;
                    g[a as usize][k] -= self.0 * t;
                }
                g[a as usize][2] += self.0 * 0.5 * l;
                g[b as usize][2] += self.0 * 0.5 * l;
            }
        }
        Ok(self.0 * total)
    }
    fn subdivided(&self, _map: &SubdivisionMap) -> Result<Arc<dyn PenaltyTerm>, RcError> {
        Err(RcError::Unsupported {
            what: "test penalty",
        })
    }
}

/// A chain pinned at both ends (z = 0.5, 12 vertices), bulging up to z = 0.62 in the middle,
/// edge lengths fixed: gravity wants it to hang down to about z = 0.36.
fn hanging_chain() -> (CurveNetwork, Vec<P>, FlowOptions) {
    let n = 14;
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let s = i as f64 / (n - 1) as f64;
            [s, 0.02 * (PI * s).sin(), 0.5 + 0.12 * (PI * s).sin()]
        })
        .collect();
    let net = CurveNetwork::new(n, chain_edges(n)).unwrap();
    let opts = FlowOptions {
        constraints: vec![
            Constraint::Point {
                vertex: 0,
                target: None,
            },
            Constraint::Point {
                vertex: n as u32 - 1,
                target: None,
            },
            Constraint::EdgeLengths {
                edges: None,
                targets: None,
            },
        ],
        threads: Some(4),
        dense_threads: Some(2),
        record_path: true,
        max_iterations: 80,
        ..FlowOptions::default()
    };
    (net, pos, opts)
}

fn min_height_on_path(report: &ddg_repulsive_curves::FlowReport) -> f64 {
    let mut m = f64::INFINITY;
    for s in &report.steps {
        for state in s.path.iter().flatten() {
            for p in state {
                m = m.min(p[2]);
            }
        }
    }
    m
}

#[test]
fn guarded_flow_never_crosses_the_obstacle() {
    let plane_z = 0.45;
    // control: without an obstacle the chain hangs through the plane z = 0.45
    {
        let (net, pos, mut opts) = hanging_chain();
        opts.penalties = vec![Arc::new(Gravity(4.0))];
        let mut s = Session::new(net, pos, opts).unwrap();
        let rep = s.run(None).unwrap();
        let m = min_height_on_path(&rep);
        assert!(m < plane_z - 0.03, "control chain reaches only z = {m}");
    }
    // The CCD alone (weight 0: no potential, the obstacle comes from the penalty) and with the
    // barrier: neither run ever puts a vertex below the plane, on any state of the accepted
    // paths (vertices move linearly between states, so the sign of z at the states decides).
    for weight in [0.0, 1e-3] {
        let obs = Arc::new(quad_at(plane_z));
        let (net, pos, mut opts) = hanging_chain();
        opts.penalties = vec![
            Arc::new(Gravity(4.0)),
            Arc::new(SurfacePotential::new(obs.clone(), weight, 3.0, 0.25).unwrap()),
        ];
        let mut s = Session::new(net, pos, opts).unwrap();
        let rep = s.run(None).unwrap();
        assert!(!matches!(rep.status, Status::Error(_)), "{:?}", rep.status);
        let m = min_height_on_path(&rep);
        assert!(
            m >= plane_z,
            "weight {weight}: the chain reached z = {m} <= {plane_z}"
        );
        // it did come down towards the plane
        assert!(
            m < plane_z + 0.1,
            "weight {weight}: the chain stayed at z = {m}"
        );
        // independent audit of the final state with the CCD against the obstacle (zero motion
        // from the final state: the distance is positive) and brute-force distance
        let fin = s.positions();
        let d = fin
            .iter()
            .map(|&x| obs.distance_to_brute(x))
            .fold(f64::INFINITY, f64::min);
        assert!(d >= 0.0);
        if d < 1e-9 {
            continue; // the Zeno approach reached the gap floor: nothing left to audit
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        let prepared = [PreparedObstacle::new(obs.clone() as Arc<dyn CollisionObstacle>).unwrap()];
        let up = vec![[0.0, 0.0, 0.0]; fin.len()];
        let r = first_contact(
            &pool,
            s.network(),
            fin,
            &up,
            1.0,
            &prepared,
            &CcdOptions::default(),
        )
        .unwrap();
        assert!(r.contact.is_none());
        // and a motion straight down through the plane is a contact
        let down = vec![[0.0, 0.0, -1.0]; fin.len()];
        let r = first_contact(
            &pool,
            s.network(),
            fin,
            &down,
            1.0,
            &prepared,
            &CcdOptions::default(),
        )
        .unwrap();
        assert!(r.contact.is_some());
    }
}

// ---------------------------------------------------------------------------------------------
// Penalties in flows
// ---------------------------------------------------------------------------------------------

fn small_options() -> FlowOptions {
    FlowOptions {
        threads: Some(4),
        dense_threads: Some(2),
        ..FlowOptions::default()
    }
}

#[test]
fn field_penalty_aligns_a_chain_with_a_constant_field() {
    // A half circle pinned at (1, 0, 0) with fixed edge lengths: the field penalty turns it into
    // a straight segment along the x axis (the sign of the direction does not matter).
    let n = 12;
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let t = PI * i as f64 / (n - 1) as f64;
            [t.cos(), t.sin(), 0.0]
        })
        .collect();
    let net = CurveNetwork::new(n, chain_edges(n)).unwrap();
    let field = Arc::new(ConstantField::new([1.0, 0.0, 0.0]).unwrap());
    let misalignment = |net: &CurveNetwork, pos: &[P]| -> f64 {
        let g = net.edge_geometry(pos);
        g.tangent
            .iter()
            .map(|t| (t[1] * t[1] + t[2] * t[2]).sqrt())
            .fold(0.0, f64::max)
    };
    let before = misalignment(&net, &pos);
    assert!(before > 0.9);
    let opts = FlowOptions {
        constraints: vec![
            Constraint::Point {
                vertex: 0,
                target: None,
            },
            Constraint::EdgeLengths {
                edges: None,
                targets: None,
            },
        ],
        penalties: vec![Arc::new(FieldAlignment::new(field, 200.0).unwrap())],
        max_iterations: 150,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let e0 = s.energy();
    let rep = s.run(None).unwrap();
    assert!(!matches!(rep.status, Status::Error(_)));
    let after = misalignment(s.network(), s.positions());
    assert!(after < 0.15, "misalignment {before} -> {after}");
    assert!(rep.final_energy < e0);
    // lengths were kept
    let g = s.network().edge_geometry(s.positions());
    let l0 = PI / (n - 1) as f64;
    for &l in &g.length {
        assert!((l / (2.0 * (0.5 * l0).sin()) - 1.0).abs() < 1e-2);
    }
}

#[test]
fn length_difference_evens_out_edge_lengths() {
    let (net, pos) = wobbly_loop(40);
    let cv = |net: &CurveNetwork, pos: &[P]| -> f64 {
        let g = net.edge_geometry(pos);
        let mean = g.length.iter().sum::<f64>() / g.length.len() as f64;
        let var = g.length.iter().map(|l| (l - mean).powi(2)).sum::<f64>() / g.length.len() as f64;
        var.sqrt() / mean
    };
    let before = cv(&net, &pos);
    let total0: f64 = net.edge_geometry(&pos).length.iter().sum();
    assert!(before > 0.2, "{before}");
    let opts = FlowOptions {
        constraints: vec![
            Constraint::Barycenter { target: None },
            Constraint::TotalLength { target: None },
        ],
        penalties: vec![Arc::new(LengthDifference::new(30.0).unwrap())],
        max_iterations: 120,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let rep = s.run(None).unwrap();
    assert!(!matches!(rep.status, Status::Error(_)));
    let after = cv(s.network(), s.positions());
    assert!(
        after < 0.4 * before,
        "edge length variation {before} -> {after}"
    );
    let total: f64 = s.network().edge_geometry(s.positions()).length.iter().sum();
    assert!((total / total0 - 1.0).abs() < 1e-3);
}

#[test]
fn total_length_penalty_shortens_the_curve() {
    let (net, pos) = wobbly_loop(30);
    let l0: f64 = net.edge_geometry(&pos).length.iter().sum();
    let opts = FlowOptions {
        constraints: vec![Constraint::Barycenter { target: None }],
        penalties: vec![Arc::new(TotalLengthPenalty::new(5.0).unwrap())],
        max_iterations: 20,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    s.run(None).unwrap();
    let l1: f64 = s.network().edge_geometry(s.positions()).length.iter().sum();
    assert!(l1 < 0.95 * l0, "{l0} -> {l1}");
}

// ---------------------------------------------------------------------------------------------
// Implicit surfaces
// ---------------------------------------------------------------------------------------------

fn all_surfaces() -> Vec<(&'static str, Arc<dyn ImplicitSurface>)> {
    let s1: Arc<dyn ImplicitSurface> = Arc::new(Sphere::new([0.0, 0.0, 0.0], 1.0).unwrap());
    let s2: Arc<dyn ImplicitSurface> = Arc::new(Sphere::new([1.2, 0.3, 0.0], 0.7).unwrap());
    vec![
        (
            "sphere",
            Arc::new(Sphere::new([0.2, -0.3, 0.1], 1.3).unwrap()),
        ),
        (
            "torus",
            Arc::new(Torus::new([0.1, 0.2, -0.1], 1.0, 0.25).unwrap()),
        ),
        (
            "plane",
            Arc::new(Plane::new([0.0, 0.1, 0.0], [0.3, 1.0, -0.2]).unwrap()),
        ),
        ("yplane", Arc::new(Plane::y_zero())),
        ("doubletorus", Arc::new(DoubleTorus::new(0.04).unwrap())),
        (
            "union",
            Arc::new(Union::new(vec![s1.clone(), s2.clone()]).unwrap()),
        ),
        (
            "smooth-union",
            Arc::new(SmoothUnion::new(s1, s2, 0.4).unwrap()),
        ),
        (
            "rounded-box",
            Arc::new(RoundedBox::new([0.0, 0.1, 0.0], [1.0, 0.7, 0.5], 0.2).unwrap()),
        ),
        (
            "box",
            Arc::new(RoundedBox::cuboid([0.0, 0.0, 0.0], [1.0, 0.7, 0.5]).unwrap()),
        ),
        (
            "metaballs",
            Arc::new(
                Metaballs::new(
                    vec![([0.0, 0.0, 0.0], 1.0, 0.8), ([0.9, 0.2, 0.0], 0.8, 0.6)],
                    0.5,
                )
                .unwrap(),
            ),
        ),
    ]
}

#[test]
fn implicit_surface_gradients_match_finite_differences() {
    let mut rng = Rng(17);
    for (name, f) in all_surfaces() {
        let mut checked = 0;
        for _ in 0..200 {
            let x = [
                rng.range(-2.0, 2.5),
                rng.range(-1.8, 1.8),
                rng.range(-1.5, 1.5),
            ];
            let g = f.gradient(x);
            let h = 1e-6;
            let mut fd = [0.0; 3];
            let mut kink = false;
            for k in 0..3 {
                let (mut a, mut b) = (x, x);
                a[k] += h;
                b[k] -= h;
                fd[k] = (f.value(a) - f.value(b)) / (2.0 * h);
                // one-sided slopes differing: a crease (union, box); skip those points
                let fwd = (f.value(a) - f.value(x)) / h;
                let bwd = (f.value(x) - f.value(b)) / h;
                if (fwd - bwd).abs() > 1e-3 * (1.0 + fwd.abs()) {
                    kink = true;
                }
            }
            if kink {
                continue;
            }
            checked += 1;
            let scale = 1.0 + norm3(g);
            assert!(
                norm3(sub3(fd, g)) < 1e-5 * scale,
                "{name} at {x:?}: {fd:?} vs {g:?}"
            );
        }
        assert!(checked > 100, "{name}: only {checked} smooth points");
    }
}

#[test]
fn implicit_surface_values_and_signs() {
    let s = Sphere::new([1.0, 0.0, 0.0], 2.0).unwrap();
    assert!((s.value([1.0, 0.0, 0.0]) + 2.0).abs() < 1e-14);
    assert!((s.value([3.0, 0.0, 0.0])).abs() < 1e-14);
    let t = Torus::new([0.0; 3], 1.0, 0.25).unwrap();
    assert!(t.value([1.0, 0.0, 0.0]) < 0.0 && t.value([0.0, 0.0, 0.0]) > 0.0);
    assert!(t.value([1.25, 0.0, 0.0]).abs() < 1e-14 && t.value([0.0, 0.25, 1.0]).abs() < 1e-14);
    let p = Plane::new([0.0, 2.0, 0.0], [0.0, 3.0, 0.0]).unwrap();
    assert!((p.value([5.0, 3.0, 1.0]) - 1.0).abs() < 1e-14);
    let b = RoundedBox::cuboid([0.0; 3], [1.0, 2.0, 3.0]).unwrap();
    assert!((b.value([1.5, 0.0, 0.0]) - 0.5).abs() < 1e-14);
    assert!((b.value([0.0, 0.0, 0.0]) + 1.0).abs() < 1e-14);
    let m = Metaballs::new(vec![([0.0; 3], 1.0, 1.0)], 0.5).unwrap();
    assert!(m.value([0.0; 3]) < 0.0 && m.value([3.0, 0.0, 0.0]) > 0.0);
    // singular points give a zero gradient, never NaN
    assert_eq!(
        Sphere::new([0.0; 3], 1.0).unwrap().gradient([0.0; 3]),
        [0.0; 3]
    );
    assert!(t.gradient([0.0, 0.0, 0.0]).iter().all(|x| x.is_finite()));
    // invalid parameters
    assert!(Sphere::new([0.0; 3], 0.0).is_err());
    assert!(Sphere::new([f64::NAN, 0.0, 0.0], 1.0).is_err());
    assert!(Torus::new([0.0; 3], -1.0, 0.2).is_err());
    assert!(Plane::new([0.0; 3], [0.0; 3]).is_err());
    assert!(RoundedBox::new([0.0; 3], [1.0, 1.0, 1.0], 1.0).is_err());
    assert!(Union::new(vec![]).is_err());
    assert!(Metaballs::new(vec![], 0.5).is_err());
    assert!(DoubleTorus::new(0.0).is_err());
    assert!(ConstantField::unit([0.0; 3]).is_err());
}

/// Moves `x` onto `f = 0` by Newton steps along the gradient.
fn project_to_surface(f: &dyn ImplicitSurface, mut x: P) -> P {
    for _ in 0..60 {
        let v = f.value(x);
        let g = f.gradient(x);
        let gg = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        if v.abs() < 1e-15 || gg == 0.0 {
            break;
        }
        for k in 0..3 {
            x[k] -= v * g[k] / gg;
        }
    }
    x
}

#[test]
fn surface_constraint_keeps_vertices_on_each_surface() {
    // The curve starts on the surface; a short flow with Barycenter + Surface constraints
    // (the reference's "constrain_all") keeps every vertex within the projection tolerance.
    let n = 30;
    let cases: Vec<(&str, Arc<dyn ImplicitSurface>, Vec<P>)> = {
        let mut v: Vec<(&str, Arc<dyn ImplicitSurface>, Vec<P>)> = Vec::new();
        for (name, f) in all_surfaces() {
            if matches!(
                name,
                "union" | "smooth-union" | "box" | "doubletorus" | "yplane"
            ) {
                continue;
            }
            // a wobbly circle around the origin region, projected onto the surface
            let pts: Vec<P> = (0..n)
                .map(|i| {
                    let t = TAU * i as f64 / n as f64;
                    let x = [
                        1.1 * t.cos(),
                        0.15 * (3.0 * t).sin() + 0.1,
                        0.9 * t.sin() + 0.1 * (2.0 * t).cos(),
                    ];
                    project_to_surface(&*f, x)
                })
                .collect();
            v.push((name, f, pts));
        }
        // torus: a (2, 3) torus knot, exactly on the surface
        let tor: Arc<dyn ImplicitSurface> = Arc::new(Torus::new([0.0; 3], 1.0, 0.4).unwrap());
        let knot: Vec<P> = (0..60)
            .map(|i| {
                let t = TAU * i as f64 / 60.0;
                let r = 1.0 + 0.4 * (3.0 * t).cos();
                [
                    r * (2.0 * t).cos(),
                    0.4 * (3.0 * t).sin(),
                    r * (2.0 * t).sin(),
                ]
            })
            .collect();
        v.push(("torus-knot", tor, knot));
        // the y = 0 plane: a planar loop
        let yp: Arc<dyn ImplicitSurface> = Arc::new(Plane::y_zero());
        let ring: Vec<P> = (0..n)
            .map(|i| {
                let t = TAU * i as f64 / n as f64;
                [t.cos() + 0.2 * (3.0 * t).cos(), 0.0, t.sin()]
            })
            .collect();
        v.push(("yplane", yp, ring));
        // the union: a loop on the bigger sphere
        let u = all_surfaces()
            .into_iter()
            .find(|(n, _)| *n == "union")
            .unwrap()
            .1;
        let on_sphere: Vec<P> = (0..n)
            .map(|i| {
                let t = TAU * i as f64 / n as f64;
                project_to_surface(
                    &*u,
                    [0.6 * t.cos(), 0.9 * t.sin(), 0.5 + 0.3 * (2.0 * t).sin()],
                )
            })
            .collect();
        v.push(("union", u, on_sphere));
        v
    };
    for (name, f, pos) in cases {
        let nv = pos.len();
        let net = CurveNetwork::new(nv, loop_edges(nv)).unwrap();
        for p in &pos {
            assert!(f.value(*p).abs() < 1e-9, "{name} start {}", f.value(*p));
        }
        let lbar = net.mean_edge_length(&pos);
        let opts = FlowOptions {
            constraints: if name.contains("plane") {
                // a plane leaves the translation along its normal to the surface rows, so the
                // barycenter would be redundant: pin vertex 0 inside the plane instead
                vec![
                    Constraint::Point {
                        vertex: 0,
                        target: None,
                    },
                    Constraint::Surface {
                        vertices: Some((1..nv as u32).collect()),
                        surface: f.clone(),
                    },
                ]
            } else {
                vec![
                    Constraint::Barycenter { target: None },
                    Constraint::Surface {
                        vertices: None,
                        surface: f.clone(),
                    },
                ]
            },
            max_iterations: 12,
            ..small_options()
        };
        let mut s = Session::new(net, pos, opts).unwrap();
        let e0 = s.energy();
        for _ in 0..12 {
            let r = s.step().unwrap();
            for &x in s.positions() {
                let tol = 2e-4 * lbar * norm3(f.gradient(x)).max(1e-12);
                assert!(
                    f.value(x).abs() <= tol.max(1e-6 * lbar),
                    "{name}: |f| = {}",
                    f.value(x)
                );
            }
            if r.tau == 0.0 {
                break;
            }
        }
        assert!(s.energy() < e0, "{name}: energy {e0} -> {}", s.energy());
        let moved = s
            .positions()
            .iter()
            .map(|p| p[0].abs() + p[1].abs() + p[2].abs())
            .sum::<f64>();
        assert!(moved.is_finite());
    }
}

// ---------------------------------------------------------------------------------------------
// Subdivision of penalties
// ---------------------------------------------------------------------------------------------

#[test]
fn penalties_follow_a_subdivision() {
    let (net, pos) = wobbly_loop(16);
    let obs = Arc::new(uv_sphere([0.0, 0.0, 0.0], 0.5, 8, 12));
    let pens: Vec<Arc<dyn PenaltyTerm>> = vec![
        Arc::new(TotalLengthPenalty::new(1.5).unwrap()),
        Arc::new(LengthDifference::new(2.0).unwrap()),
        Arc::new(FieldAlignment::new(Arc::new(CircularField), 1.0).unwrap()),
        Arc::new(SurfacePotential::new(obs.clone(), 0.5, 3.0, 0.25).unwrap()),
        Arc::new(ddg_repulsive_curves::optimize::Hencky::at_rest(3.0, &net, &pos).unwrap()),
    ];
    let (net2, pos2, map) = subdivide(&net, &pos).unwrap();
    for p in &pens {
        let q = p.subdivided(&map).unwrap();
        q.check(&net2).unwrap();
        assert_eq!(q.name(), p.name());
        assert_eq!(q.obstacle().is_some(), p.obstacle().is_some());
        // finite value and a correct gradient on the new network
        let v = q.evaluate(&net2, &pos2, None).unwrap();
        assert!(v.is_finite(), "{}", q.name());
        check_gradient(&*q, &net2, &pos2, 1e-6);
        if q.name() == "TotalLength" {
            // splitting at midpoints keeps the length
            let v0 = p.evaluate(&net, &pos, None).unwrap();
            assert!((v - v0).abs() < 1e-12 * v0.abs());
        }
        if q.name() == "Hencky" {
            // rest lengths are halved with the edges: zero strain at the midpoint split
            assert!(v.abs() < 1e-24 + 1e-20);
        }
    }
    // In a session: the schedule grows the curve past twice its mean edge length; the penalties
    // are mapped at the subdivision and the flow continues.
    let (net, pos) = wobbly_loop(24);
    let opts = FlowOptions {
        schedule: Some(LengthSchedule {
            factor: 2.4,
            rate: 0.1,
        }),
        max_subdivisions: 1,
        penalties: vec![
            Arc::new(TotalLengthPenalty::new(0.1).unwrap()),
            Arc::new(LengthDifference::new(0.1).unwrap()),
            Arc::new(SurfacePotential::new(obs.clone(), 0.01, 3.0, 0.25).unwrap()),
        ],
        max_iterations: 60,
        ..small_options()
    };
    let mut s = Session::new(net, pos, opts).unwrap();
    let mut subdivided = false;
    for _ in 0..60 {
        let r = s.step().unwrap();
        if r.subdivided {
            subdivided = true;
            assert_eq!(s.network().n_vertices(), 48);
            break;
        }
    }
    assert!(subdivided);
    for p in s.penalties() {
        p.check(s.network()).unwrap();
        assert!(
            p.evaluate(s.network(), s.positions(), None)
                .unwrap()
                .is_finite()
        );
    }
    for _ in 0..3 {
        s.step().unwrap();
    }
    assert!(s.energy().is_finite());
}

// ---------------------------------------------------------------------------------------------
// Obstacle loading
// ---------------------------------------------------------------------------------------------

#[test]
fn obj_quad_mesh_is_fan_triangulated() {
    let text = "# a unit quad and a cube face above it\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\nf 1 2 3 4\nf 5 6 7 8\n";
    let o = Obstacle::parse_obj(text).unwrap();
    assert_eq!(o.n_triangles(), 4);
    assert_eq!(o.triangle_indices()[0], [0, 1, 2]);
    assert_eq!(o.triangle_indices()[1], [0, 2, 3]);
    assert_eq!(o.face_of_triangle(), &[0, 0, 1, 1]);
    assert!((o.areas().iter().sum::<f64>() - 2.0).abs() < 1e-14);
    assert_eq!(CollisionObstacle::triangles(&o).len(), 4);
    assert!((o.distance_to([0.5, 0.5, 0.25]).0 - 0.25).abs() < 1e-14);
}

#[test]
fn obj_ngons_and_non_manifold_meshes_load() {
    // a pentagon, a triangle, and three faces on one edge (non-manifold), negative indices
    let text = "v 0 0 0\nv 2 0 0\nv 3 1 0\nv 1 2 0\nv -1 1 0\nv 0 0 1\nv 0 0 -1\nv 0 1 0.5\n\
                f 1 2 3 4 5\nf 1 2 6\nf 1 2 7\nf 1 2 -1\n";
    let o = Obstacle::parse_obj(text).unwrap();
    assert_eq!(o.n_triangles(), 3 + 1 + 1 + 1);
    assert_eq!(o.face_of_triangle(), &[0, 0, 0, 1, 2, 3]);
    // the pentagon's fan has total area 6.0 - ...: check via shoelace of (0,0),(2,0),(3,1),(1,2),(-1,1)
    let pent = [[0.0, 0.0], [2.0, 0.0], [3.0, 1.0], [1.0, 2.0], [-1.0, 1.0]];
    let shoelace: f64 = 0.5
        * (0..5)
            .map(|i| pent[i][0] * pent[(i + 1) % 5][1] - pent[(i + 1) % 5][0] * pent[i][1])
            .sum::<f64>();
    let fan: f64 = o.areas()[..3].iter().sum();
    assert!((fan - shoelace.abs()).abs() < 1e-12, "{fan} vs {shoelace}");
    // a potential and a distance work on it
    assert!(o.distance_to([0.5, 0.5, 0.5]).0.is_finite());
}

#[test]
fn obj_errors_are_reported_not_panicked() {
    // no faces
    assert!(matches!(
        Obstacle::parse_obj("v 0 0 0\nv 1 0 0\n"),
        Err(RcError::InvalidOption { .. })
    ));
    // index out of range: the parser names the line
    match Obstacle::parse_obj("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 9\n") {
        Err(RcError::Parse { line, message }) => {
            assert_eq!(line, 4, "{message}");
        }
        other => panic!("{other:?}"),
    }
    // malformed number
    assert!(matches!(
        Obstacle::parse_obj("v 0 x 0\n"),
        Err(RcError::Parse { line: 1, .. })
    ));
    // a face with two corners
    assert!(matches!(
        Obstacle::parse_obj("v 0 0 0\nv 1 0 0\nf 1 2\n"),
        Err(RcError::Parse { .. })
    ));
    // a degenerate (collinear) triangle
    assert!(matches!(
        Obstacle::parse_obj("v 0 0 0\nv 1 0 0\nv 2 0 0\nf 1 2 3\n"),
        Err(RcError::DegenerateObstacleTriangle { triangle: 0 })
    ));
    // a quad with a repeated corner has a degenerate second fan triangle
    assert!(matches!(
        Obstacle::parse_obj("v 0 0 0\nv 1 0 0\nv 1 1 0\nf 1 2 3 3\n"),
        Err(RcError::DegenerateObstacleTriangle { triangle: 1 })
    ));
    // missing file
    assert!(matches!(
        Obstacle::read_obj("/nonexistent/dir/obstacle.obj"),
        Err(RcError::Io { .. })
    ));
    // direct construction
    let v = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(matches!(
        Obstacle::new(v.clone(), vec![[0, 1, 3]]),
        Err(RcError::ObstacleTriangleOutOfRange {
            triangle: 0,
            vertex: 3,
            n_vertices: 3
        })
    ));
    assert!(matches!(
        Obstacle::new(v.clone(), vec![]),
        Err(RcError::InvalidOption { .. })
    ));
    let mut bad = v.clone();
    bad[1][2] = f64::NAN;
    assert!(matches!(
        Obstacle::new(bad, vec![[0, 1, 2]]),
        Err(RcError::NonFinitePosition { vertex: 1 })
    ));
    assert!(Obstacle::from_soup(&[]).is_err());
}

#[test]
fn obj_file_round_trip_from_disk() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    let path = dir.join("rc_wp7_obstacle.obj");
    std::fs::write(&path, "v 0 0 0\nv 4 0 0\nv 4 4 0\nv 0 4 0\nf 1 2 3 4\n").unwrap();
    let o = Obstacle::read_obj(&path).unwrap();
    assert_eq!(o.n_triangles(), 2);
    assert!((o.areas().iter().sum::<f64>() - 16.0).abs() < 1e-12);
    let bad = dir.join("rc_wp7_bad.obj");
    std::fs::write(&bad, "v 0 0 0\nf 1 2 3\n").unwrap();
    match Obstacle::read_obj(&bad) {
        Err(RcError::Parse { message, .. }) => {
            assert!(message.contains("rc_wp7_bad.obj"), "{message}")
        }
        other => panic!("{other:?}"),
    }
    // the authors' scene meshes, if the reference checkout is next to the repository
    let bunny = std::path::Path::new(
        "/var/data/discrete-geometry/RepulsiveCurves/repulsive-curves/scenes/BunnyCurvePacking/bunny.obj",
    );
    if bunny.exists() {
        let o = Obstacle::read_obj(bunny).unwrap();
        assert!(o.n_triangles() > 100);
    }
}

// ---------------------------------------------------------------------------------------------
// Cost of the surface potential (ignored: run in release, alone)
// ---------------------------------------------------------------------------------------------

fn load_average() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .unwrap_or_default()
        .trim()
        .to_string()
}

#[test]
#[ignore = "timing; run with --release --ignored --nocapture"]
fn cost_of_the_surface_potential() {
    // 10k-triangle obstacle (UV sphere, radius 1) against a 1000-vertex curve in the shell
    // 1.2 < r < 2.4 (a (3, 7) torus-like knot).
    let obs = Arc::new(uv_sphere([0.0; 3], 1.0, 71, 71)); // 2*71*70 - ... ~ 9940 triangles
    let n = 1000;
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let t = TAU * i as f64 / n as f64;
            let r = 1.8 + 0.5 * (7.0 * t).cos();
            [
                r * (3.0 * t).cos(),
                r * (3.0 * t).sin(),
                0.5 * (7.0 * t).sin() + 0.4 * (2.0 * t).cos(),
            ]
        })
        .collect();
    let net = CurveNetwork::new(n, loop_edges(n)).unwrap();
    let pen = SurfacePotential::new(obs.clone(), 1.0, 3.0, 0.25).unwrap();
    eprintln!(
        "obstacle: {} triangles, curve: {n} vertices, load average at start: {}",
        obs.n_triangles(),
        load_average()
    );
    let mut rows = Vec::new();
    for threads in [1usize, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for (label, brute) in [("brute", true), ("bvh", false)] {
            let mut best = f64::INFINITY;
            let mut value = 0.0;
            for _ in 0..3 {
                let t0 = std::time::Instant::now();
                let mut g = vec![[0.0; 3]; n];
                let v = pool
                    .install(|| {
                        if brute {
                            pen.evaluate_brute(&net, &pos, Some(&mut g))
                        } else {
                            pen.evaluate(&net, &pos, Some(&mut g))
                        }
                    })
                    .unwrap();
                best = best.min(t0.elapsed().as_secs_f64());
                value = v;
            }
            rows.push((threads, label, best, value));
        }
    }
    let mut stats = (0u64, 0u64);
    for &x in &pos {
        let q = obs.potential_at(x, 3.0, 0.25, true);
        stats.0 += q.near_triangles as u64;
        stats.1 += q.far_nodes as u64;
    }
    eprintln!(
        "per vertex: {:.1} exact triangles, {:.1} lumped nodes",
        stats.0 as f64 / n as f64,
        stats.1 as f64 / n as f64
    );
    eprintln!("threads  method  best-of-3 seconds (value + gradient)  value");
    for (t, l, s, v) in &rows {
        eprintln!("{t:>7}  {l:<6}  {s:.5}  {v:.9e}");
    }
    eprintln!("load average at end: {}", load_average());
}
