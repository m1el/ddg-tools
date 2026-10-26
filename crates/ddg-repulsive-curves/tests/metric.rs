//! Tests of the metrics (`SPEC.md` §4, §13.4): the Gram increments against the quadratic
//! forms, Eq. 31 against Eq. 20–21, symmetry, positive semidefiniteness, null spaces (D5),
//! scaling (D3), determinism (D25), and the classic metrics (§4.6).

#![allow(clippy::needless_range_loop)]

use ddg_repulsive_curves::energy::exact::make_pool;
use ddg_repulsive_curves::metric::sobolev_dense::{GramGeometry, sobolev_gram_increments};
use ddg_repulsive_curves::metric::{
    DenseMatrix, GramQuadrature, Metric, apply_average, apply_average_transpose, apply_derivative,
    apply_derivative_transpose, bilaplacian, curve_laplacian, gram, kernel_matrices, mass_matrix,
    sobolev_gram,
};
use ddg_repulsive_curves::{CurveNetwork, Exponents};

type P = [f64; 3];

const QUADS: [GramQuadrature; 2] = [GramQuadrature::Midpoint, GramQuadrature::Trapezoid];

// ---------------------------------------------------------------- curves

fn closed_edges(n: usize, offset: u32) -> Vec<[u32; 2]> {
    (0..n as u32)
        .map(|k| [offset + k, offset + (k + 1) % n as u32])
        .collect()
}

fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            [
                t.sin() + 2.0 * (2.0 * t).sin(),
                t.cos() - 2.0 * (2.0 * t).cos(),
                -(3.0 * t).sin(),
            ]
        })
        .collect();
    (CurveNetwork::new(n, closed_edges(n, 0)).unwrap(), pos)
}

fn open_curve(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos = (0..n)
        .map(|k| {
            let t = 3.0 * k as f64 / n as f64;
            [t.cos(), t.sin(), 0.4 * t + 0.1 * (2.0 * t).sin()]
        })
        .collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn tripod(m: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = vec![[0.0, 0.0, 0.0]];
    let mut edges = Vec::new();
    let dirs = [[1.0, 0.0, 0.2], [-0.5, 0.9, -0.1], [-0.5, -0.8, 0.3]];
    for (a, d) in dirs.iter().enumerate() {
        let mut prev = 0u32;
        for k in 1..=m {
            let s = k as f64 / m as f64;
            let bend = 0.3 * s * s;
            pos.push([
                d[0] * s + bend * (a as f64 + 1.0).cos(),
                d[1] * s + bend * (a as f64 + 1.0).sin(),
                d[2] * s + 0.2 * (3.0 * s + a as f64).sin(),
            ]);
            let id = pos.len() as u32 - 1;
            // Alternate the orientation to exercise both signs.
            edges.push(if k % 2 == 0 { [prev, id] } else { [id, prev] });
            prev = id;
        }
    }
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

/// Two linked circles (Hopf link).
fn hopf(n: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = Vec::new();
    for k in 0..n {
        let t = std::f64::consts::TAU * k as f64 / n as f64;
        pos.push([t.cos(), t.sin(), 0.05 * (2.0 * t).sin()]);
    }
    for k in 0..n {
        let t = std::f64::consts::TAU * k as f64 / n as f64;
        pos.push([1.0 + t.cos(), 0.03 * (3.0 * t).cos(), t.sin()]);
    }
    let mut edges = closed_edges(n, 0);
    edges.extend(closed_edges(n, n as u32));
    (CurveNetwork::new(2 * n, edges).unwrap(), pos)
}

fn cases() -> Vec<(&'static str, CurveNetwork, Vec<P>)> {
    let (a, pa) = trefoil(24);
    let (b, pb) = open_curve(18);
    let (c, pc) = tripod(5);
    let (d, pd) = hopf(12);
    vec![
        ("trefoil", a, pa),
        ("open", b, pb),
        ("tripod", c, pc),
        ("hopf", d, pd),
    ]
}

/// Deterministic pseudo-random numbers in [−1, 1].
fn rand_vec(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        })
        .collect()
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The quadratic forms `uᵀBv`, `uᵀB⁰v` computed directly from Eq. 20–21 (independent of the
/// crate's weight code).
fn direct_forms(
    net: &CurveNetwork,
    pos: &[P],
    e: Exponents,
    q: GramQuadrature,
    u: &[f64],
    v: &[f64],
) -> (f64, f64) {
    let sigma = (e.beta() - 1.0) / e.alpha() - 1.0;
    let edges = net.edges();
    let m = edges.len();
    let len = |i: usize| {
        let d = sub(pos[edges[i][1] as usize], pos[edges[i][0] as usize]);
        dot(d, d).sqrt()
    };
    let tan = |i: usize| {
        let d = sub(pos[edges[i][1] as usize], pos[edges[i][0] as usize]);
        let l = dot(d, d).sqrt();
        [d[0] / l, d[1] / l, d[2] / l]
    };
    let mid = |i: usize| {
        let (a, b) = (pos[edges[i][0] as usize], pos[edges[i][1] as usize]);
        [
            (a[0] + b[0]) / 2.0,
            (a[1] + b[1]) / 2.0,
            (a[2] + b[2]) / 2.0,
        ]
    };
    let du = |f: &[f64], i: usize| {
        let s = (f[edges[i][1] as usize] - f[edges[i][0] as usize]) / len(i);
        let t = tan(i);
        [s * t[0], s * t[1], s * t[2]]
    };
    let eu = |f: &[f64], i: usize| 0.5 * (f[edges[i][0] as usize] + f[edges[i][1] as usize]);
    // k²_{2σ+5}(p, q, T) = |T × d|² / |d|^{2σ+5}
    let k2 = |p: P, qq: P, t: P| {
        let d = sub(p, qq);
        let c = [
            t[1] * d[2] - t[2] * d[1],
            t[2] * d[0] - t[0] * d[2],
            t[0] * d[1] - t[1] * d[0],
        ];
        dot(c, c) / dot(d, d).sqrt().powf(2.0 * sigma + 5.0)
    };
    let (mut fb, mut fb0) = (0.0, 0.0);
    for i in 0..m {
        for j in 0..m {
            if net.edges_adjacent(i, j) {
                continue;
            }
            let ll = len(i) * len(j);
            let (w, w0) = match q {
                GramQuadrature::Midpoint => {
                    let d = sub(mid(i), mid(j));
                    (
                        ll / dot(d, d).sqrt().powf(2.0 * sigma + 1.0),
                        ll * k2(mid(i), mid(j), tan(i)),
                    )
                }
                GramQuadrature::Trapezoid => {
                    let (mut w, mut w0) = (0.0, 0.0);
                    for a in edges[i] {
                        for b in edges[j] {
                            let (p, qq) = (pos[a as usize], pos[b as usize]);
                            let d = sub(p, qq);
                            w += 1.0 / dot(d, d).sqrt().powf(2.0 * sigma + 1.0);
                            w0 += k2(p, qq, tan(i));
                        }
                    }
                    (0.25 * ll * w, 0.25 * ll * w0)
                }
            };
            let gu = sub(du(u, i), du(u, j));
            let gv = sub(du(v, i), du(v, j));
            fb += w * dot(gu, gv);
            fb0 += w0 * (eu(u, i) - eu(u, j)) * (eu(v, i) - eu(v, j));
        }
    }
    (fb, fb0)
}

fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-300)
}

fn max_abs_diff(a: &DenseMatrix, b: &DenseMatrix) -> f64 {
    a.as_slice()
        .iter()
        .zip(b.as_slice())
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
}

fn eigenvalues(m: &DenseMatrix) -> Vec<f64> {
    let n = m.n();
    let a = faer::Mat::<f64>::from_fn(n, n, |i, j| m.get(i, j));
    let mut ev = a.self_adjoint_eigenvalues(faer::Side::Lower).unwrap();
    ev.sort_by(f64::total_cmp);
    ev
}

fn kernel_dim(m: &DenseMatrix, rel: f64) -> usize {
    let ev = eigenvalues(m);
    let max = ev.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    ev.iter().filter(|&&x| x.abs() <= rel * max).count()
}

// ---------------------------------------------------------------- tests

#[test]
fn increments_equal_quadratic_forms() {
    for e in [Exponents::default(), Exponents::new(2.0, 4.5).unwrap()] {
        for (name, net, pos) in cases() {
            for q in QUADS {
                let g = sobolev_gram_increments(&net, &pos, e, q);
                let n = net.n_vertices();
                for seed in 0..3 {
                    let u = rand_vec(n, 10 + seed);
                    let v = rand_vec(n, 100 + seed);
                    let (fb, fb0) = direct_forms(&net, &pos, e, q, &u, &v);
                    let (b, b0) = (g.b.bilinear(&u, &v), g.b0.bilinear(&u, &v));
                    assert!(rel_close(b, fb, 1e-12), "{name} {q:?} B: {b} vs {fb}");
                    assert!(rel_close(b0, fb0, 1e-12), "{name} {q:?} B0: {b0} vs {fb0}");
                    let a = g.a.bilinear(&u, &v);
                    assert!(rel_close(a, fb + fb0, 1e-12), "{name} {q:?} A");
                }
            }
        }
    }
}

#[test]
fn kernel_form_equals_increments() {
    let pool = make_pool(Some(4)).unwrap();
    for e in [Exponents::default(), Exponents::new(2.0, 4.5).unwrap()] {
        for (name, net, pos) in cases() {
            for q in QUADS {
                let fast = sobolev_gram(&pool, &net, &pos, e, q);
                let slow = sobolev_gram_increments(&net, &pos, e, q);
                for (what, f, s) in [("B", &fast.b, &slow.b), ("B0", &fast.b0, &slow.b0)] {
                    let scale = s.max_abs();
                    let diff = max_abs_diff(f, s);
                    assert!(
                        diff <= 1e-12 * scale,
                        "{name} {q:?} {what}: diff {diff:e}, scale {scale:e}"
                    );
                }
            }
        }
    }
}

#[test]
fn kernel_matrices_are_symmetric_and_vanish_on_adjacent_pairs() {
    let pool = make_pool(Some(3)).unwrap();
    for (name, net, pos) in cases() {
        for q in QUADS {
            let geom = GramGeometry::new(&net, &pos, Exponents::default());
            let k = kernel_matrices(&pool, &net, &geom, q);
            assert_eq!(k.kb.symmetry_defect(), 0.0, "{name}");
            assert_eq!(k.k0.symmetry_defect(), 0.0, "{name}");
            for i in 0..net.n_edges() {
                for j in 0..net.n_edges() {
                    if net.edges_adjacent(i, j) {
                        assert_eq!(k.kb.get(i, j), 0.0);
                        assert_eq!(k.k0.get(i, j), 0.0);
                    } else {
                        assert!(k.kb.get(i, j) > 0.0);
                        assert!(k.k0.get(i, j) >= 0.0);
                        // K = w + wᵀ (§4.4).
                        let (wij, w0ij) = geom.weights(i, j, q);
                        let (wji, w0ji) = geom.weights(j, i, q);
                        assert!(rel_close(k.kb.get(i, j), wij + wji, 1e-14));
                        assert!(rel_close(k.k0.get(i, j), w0ij + w0ji, 1e-13));
                    }
                }
            }
        }
    }
}

#[test]
fn derivative_operators_are_adjoint() {
    let (net, pos) = tripod(4);
    let geom = net.edge_geometry(&pos);
    let n = net.n_vertices();
    let m = net.n_edges();
    let u = rand_vec(n, 3);
    let y: Vec<P> = rand_vec(3 * m, 4)
        .chunks(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();
    let du = apply_derivative(&net, &geom, &u);
    let lhs: f64 = du.iter().zip(&y).map(|(a, b)| dot(*a, *b)).sum();
    let dty = apply_derivative_transpose(&net, &geom, &y);
    let rhs: f64 = u.iter().zip(&dty).map(|(a, b)| a * b).sum();
    assert!(rel_close(lhs, rhs, 1e-13));
    let z = rand_vec(m, 5);
    let eu = apply_average(&net, &u);
    let lhs: f64 = eu.iter().zip(&z).map(|(a, b)| a * b).sum();
    let etz = apply_average_transpose(&net, &z);
    let rhs: f64 = u.iter().zip(&etz).map(|(a, b)| a * b).sum();
    assert!(rel_close(lhs, rhs, 1e-13));
}

#[test]
fn symmetric_and_positive_semidefinite() {
    let pool = make_pool(Some(2)).unwrap();
    for (name, net, pos) in cases() {
        for q in QUADS {
            let g = sobolev_gram(&pool, &net, &pos, Exponents::default(), q);
            for (what, m) in [("B", &g.b), ("B0", &g.b0), ("A", &g.a)] {
                assert_eq!(m.symmetry_defect(), 0.0, "{name} {what}");
                let ev = eigenvalues(m);
                let max = ev[ev.len() - 1];
                assert!(ev[0] >= -1e-12 * max, "{name} {q:?} {what}: λmin {}", ev[0]);
            }
        }
    }
}

/// `+1/−1` on the two colour classes of each bipartite component (0 on the others): `E` kills
/// it, so it lies in `ker B⁰`. `None` if no component is bipartite.
fn bipartite_signs(net: &CurveNetwork) -> (usize, Vec<f64>) {
    let n = net.n_vertices();
    let mut colour = vec![0i8; n];
    let mut count = 0;
    let mut out = vec![0.0; n];
    for comp in net.components() {
        let s = comp.vertices[0] as usize;
        colour[s] = 1;
        let mut stack = vec![s];
        let mut ok = true;
        while let Some(v) = stack.pop() {
            for &e in net.incident_edges(v as u32) {
                let [a, b] = net.edges()[e as usize];
                let w = if a as usize == v { b } else { a } as usize;
                if colour[w] == 0 {
                    colour[w] = -colour[v];
                    stack.push(w);
                } else if colour[w] == colour[v] {
                    ok = false;
                }
            }
        }
        if ok {
            count += 1;
            for &v in &comp.vertices {
                out[v as usize] = colour[v as usize] as f64;
            }
        }
    }
    (count, out)
}

#[test]
fn null_spaces() {
    let pool = make_pool(Some(2)).unwrap();
    for q in QUADS {
        // One component: ker B = ker A = span{𝟙}; ker B⁰ = span{𝟙} plus the sign function of
        // a bipartite network (even loops, open curves, trees), which E kills.
        for (name, net, pos) in [trefoil(25), trefoil(24), open_curve(18), tripod(5)]
            .into_iter()
            .zip(["trefoil25", "trefoil24", "open", "tripod"])
            .map(|((n, p), s)| (s, n, p))
        {
            let g = sobolev_gram(&pool, &net, &pos, Exponents::default(), q);
            let (bip, sign) = bipartite_signs(&net);
            assert_eq!(kernel_dim(&g.b, 1e-10), 1, "{name} {q:?} B");
            assert_eq!(kernel_dim(&g.b0, 1e-10), 1 + bip, "{name} {q:?} B0");
            assert_eq!(kernel_dim(&g.a, 1e-10), 1, "{name} {q:?} A");
            let ones = vec![1.0; net.n_vertices()];
            let scale = g.a.max_abs();
            for x in g.a.mul_vec(&ones) {
                assert!(x.abs() <= 1e-12 * scale);
            }
            for x in g.b0.mul_vec(&sign) {
                assert!(x.abs() <= 1e-12 * g.b0.max_abs());
            }
        }
        // Two components: ker B = per-component constants (dim 2), ker A = span{𝟙}
        // (D5: B⁰ couples the components); odd loops, so ker B⁰ = span{𝟙}.
        for (n, bip) in [(13, 0), (12, 2)] {
            let (net, pos) = hopf(n);
            assert_eq!(bipartite_signs(&net).0, bip);
            let g = sobolev_gram(&pool, &net, &pos, Exponents::default(), q);
            assert_eq!(kernel_dim(&g.b, 1e-10), 2, "hopf B");
            assert_eq!(kernel_dim(&g.b0, 1e-10), 1 + bip, "hopf B0");
            assert_eq!(kernel_dim(&g.a, 1e-10), 1, "hopf A");
            let ind: Vec<f64> = (0..2 * n).map(|v| if v < n { 1.0 } else { 0.0 }).collect();
            let scale = g.b.max_abs();
            for x in g.b.mul_vec(&ind) {
                assert!(x.abs() <= 1e-12 * scale, "B kills a component constant");
            }
            let b0i = g.b0.mul_vec(&ind);
            assert!(b0i.iter().any(|x| x.abs() > 1e-6 * g.b0.max_abs()));
        }
    }
}

#[test]
fn scaling_exponent() {
    let pool = make_pool(Some(2)).unwrap();
    for e in [Exponents::default(), Exponents::new(2.0, 4.5).unwrap()] {
        let sigma = e.sigma();
        for q in QUADS {
            let (net, pos) = trefoil(20);
            let g = sobolev_gram(&pool, &net, &pos, e, q);
            for c in [1e-3, 0.7, 1e3] {
                let scaled: Vec<P> = pos.iter().map(|p| [c * p[0], c * p[1], c * p[2]]).collect();
                let gs = sobolev_gram(&pool, &net, &scaled, e, q);
                let f = c.powf(-(2.0 * sigma + 1.0));
                for (what, a, b) in [
                    ("B", &g.b, &gs.b),
                    ("B0", &g.b0, &gs.b0),
                    ("A", &g.a, &gs.a),
                ] {
                    let scale = a.max_abs();
                    let diff = a
                        .as_slice()
                        .iter()
                        .zip(b.as_slice())
                        .fold(0.0f64, |m, (x, y)| m.max((x * f - y).abs()));
                    assert!(
                        diff <= 1e-13 * scale * f,
                        "{what} {q:?} c={c}: {diff:e} vs {:e}",
                        scale * f
                    );
                }
            }
        }
    }
}

#[test]
fn bitwise_independent_of_thread_count() {
    let (net, pos) = trefoil(97);
    for q in QUADS {
        let reference = sobolev_gram(
            &make_pool(Some(1)).unwrap(),
            &net,
            &pos,
            Exponents::default(),
            q,
        );
        for t in [7, 64] {
            let g = sobolev_gram(
                &make_pool(Some(t)).unwrap(),
                &net,
                &pos,
                Exponents::default(),
                q,
            );
            assert!(g == reference, "{t} threads, {q:?}");
        }
    }
}

#[test]
fn quadratures_differ_on_coarse_curves_and_agree_under_refinement() {
    let pool = make_pool(Some(2)).unwrap();
    let rel_diff = |n: usize| {
        let (net, pos) = trefoil(n);
        let a = sobolev_gram(
            &pool,
            &net,
            &pos,
            Exponents::default(),
            GramQuadrature::Midpoint,
        );
        let b = sobolev_gram(
            &pool,
            &net,
            &pos,
            Exponents::default(),
            GramQuadrature::Trapezoid,
        );
        // Compare the quadratic forms on a smooth function.
        let u: Vec<f64> = pos.iter().map(|p| p[0] + 0.5 * p[2]).collect();
        let (fa, fb) = (a.a.bilinear(&u, &u), b.a.bilinear(&u, &u));
        (fa - fb).abs() / fa.abs()
    };
    let (coarse, fine) = (rel_diff(12), rel_diff(96));
    assert!(coarse > 1e-3, "coarse {coarse}");
    assert!(fine < coarse, "fine {fine} coarse {coarse}");
}

#[test]
fn gram_entry_point_validates_and_returns_classic_metrics() {
    let (net, pos) = tripod(4);
    let n = net.n_vertices();
    let e = Exponents::default();
    // Invalid positions are rejected.
    let mut bad = pos.clone();
    bad[3] = bad[2];
    assert!(gram(&net, &bad, e, &Metric::default(), Some(1)).is_err());
    assert!(gram(&net, &pos, e, &Metric::default(), Some(0)).is_err());
    // L²: diagonal dual lengths.
    let l2 = gram(&net, &pos, e, &Metric::L2, None).unwrap();
    let m = net.dual_lengths(&pos);
    for v in 0..n {
        assert_eq!(l2.a.get(v, v), m[v]);
        assert_eq!(l2.b0.get(v, v), m[v]);
    }
    // H¹: symmetric, kills constants, rank n − 1 (one component).
    let h1 = gram(&net, &pos, e, &Metric::H1, None).unwrap();
    assert_eq!(h1.a.symmetry_defect(), 0.0);
    for x in h1.a.mul_vec(&vec![1.0; n]) {
        assert!(x.abs() < 1e-12);
    }
    assert_eq!(kernel_dim(&h1.a, 1e-12), 1);
    // H² = L M⁻¹ L.
    let h2 = gram(&net, &pos, e, &Metric::H2, None).unwrap();
    for i in 0..n {
        for j in 0..n {
            let s: f64 = (0..n).map(|k| h1.a.get(i, k) * h1.a.get(k, j) / m[k]).sum();
            assert!((h2.a.get(i, j) - s).abs() <= 1e-12 * h2.a.max_abs());
        }
    }
    // Sparse entry points agree with the dense ones.
    let l = curve_laplacian(&net, &pos).unwrap();
    assert_eq!(l.nnz(), n + 2 * net.n_edges());
    assert_eq!(mass_matrix(&net, &pos).unwrap().diagonal(), m);
    assert!(bilaplacian(&net, &pos).unwrap().symmetry_defect().unwrap() <= 1e-12);
    // Two components: the per-component constants are in ker L.
    let (net, pos) = hopf(8);
    let h1 = gram(&net, &pos, e, &Metric::H1, None).unwrap();
    assert_eq!(kernel_dim(&h1.a, 1e-12), 2);
}

#[test]
fn collinear_network_has_vanishing_low_order_term() {
    let n = 8;
    let pos: Vec<P> = (0..n).map(|k| [k as f64 * 0.7, 0.0, 0.0]).collect();
    let edges = (0..n as u32 - 1).map(|k| [k, k + 1]).collect();
    let net = CurveNetwork::new(n, edges).unwrap();
    let g = gram(
        &net,
        &pos,
        Exponents::default(),
        &Metric::default(),
        Some(1),
    )
    .unwrap();
    // k² ≡ 0 up to rounding (|T×d|² = |d|² − (T·d)² is computed by cancellation).
    assert!(g.b0.max_abs() <= 1e-14 * g.b.max_abs());
    // ker B contains the linear functions along the line as well as the constants.
    assert_eq!(kernel_dim(&g.b, 1e-10), 2);
}
