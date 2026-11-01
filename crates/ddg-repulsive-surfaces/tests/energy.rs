//! Exponents, kernel, exact energy and its differential (`SPEC.md` §2.2, §3.1–3.3, §13.2).

#![allow(clippy::needless_range_loop)]

use ddg_core::generators::{SplitMix64, bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::energy::exact::{energy_value, make_pool};
use ddg_repulsive_surfaces::geometry::{FaceGeometry, mean_edge_length};
use ddg_repulsive_surfaces::kernel::{ipow, kernel_with_gradient};
use ddg_repulsive_surfaces::{EnergyEval, Exponents, RsError, Surface, energy, kernel};

/// A Möbius band with `n` quads (2 triangles each), one half twist.
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
        let (tj, bj) = if i + 1 == n {
            (b(0), t(0))
        } else {
            (t(i + 1), b(i + 1))
        };
        faces.push([ti, bi, bj]);
        faces.push([ti, bj, tj]);
    }
    (pos, faces)
}

/// The test meshes of `SPEC.md` §13.2, randomly perturbed.
fn test_meshes() -> Vec<(&'static str, Surface, Vec<[f64; 3]>)> {
    let mut out = Vec::new();
    let mut add = |name: &'static str, mut pos: Vec<[f64; 3]>, faces: Vec<[u32; 3]>, seed| {
        let s = Surface::new(pos.len(), faces).unwrap();
        let l = mean_edge_length(&s, &pos);
        perturb(&mut pos, [0.1 * l; 3], seed);
        out.push((name, s, pos));
    };
    let m = icosphere(1).unwrap();
    add("icosphere", m.positions, m.faces, 1);
    let m = torus(8, 5, 1.5, 0.6).unwrap();
    add("torus", m.positions, m.faces, 2);
    let m = bumpy_grid_disk(4, 4, 0.3).unwrap();
    add("disk", m.positions, m.faces, 3);
    let (p, f) = moebius(10);
    add("moebius", p, f, 4);
    // Two components: two icospheres, one shifted and scaled.
    let m = icosphere(0).unwrap();
    let mut p = m.positions.clone();
    let mut f = m.faces.clone();
    let off = p.len() as u32;
    p.extend(
        m.positions
            .iter()
            .map(|q| [0.7 * q[0] + 2.2, 0.7 * q[1] + 0.3, 0.7 * q[2]]),
    );
    f.extend(m.faces.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    add("two components", p, f, 5);
    out
}

fn exponent_sets() -> Vec<Exponents> {
    vec![
        Exponents::default(),
        Exponents::new(4.5, 9.5).unwrap(),
        Exponents::new(3.0, 6.0).unwrap(),
    ]
}

#[test]
fn exponents_are_validated() {
    let bad = |a, b| {
        assert!(
            matches!(Exponents::new(a, b), Err(RsError::InvalidExponents { .. })),
            "({a}, {b})"
        )
    };
    bad(1.0, 3.5);
    bad(0.5, 3.0);
    bad(6.0, 8.0); // β = α + 2: s = 1
    bad(6.0, 14.0); // β = 2α + 2: s = 2
    bad(f64::NAN, 12.0);
    bad(6.0, f64::INFINITY);
    let e = Exponents::default();
    assert_eq!((e.alpha(), e.beta()), (6.0, 12.0));
    assert_eq!(e, Exponents::new(6.0, 12.0).unwrap());
    assert_eq!(e, Exponents::from_p(6.0).unwrap());
    assert!((e.s() - 10.0 / 6.0).abs() < 1e-15);
    assert!((e.s() - (2.0 - 2.0 / 6.0)).abs() < 1e-15);
    assert!(!e.subcritical());
    for e in exponent_sets() {
        assert!(e.s() > 1.0 && e.s() < 2.0);
        assert!((e.sigma() - (e.s() - 1.0)).abs() < 1e-15);
        assert!((e.high_order_kernel_exponent() - 2.0 * e.s()).abs() < 1e-14);
        assert!((e.preconditioner_kernel_exponent() - (6.0 - 2.0 * e.s())).abs() < 1e-14);
        assert_eq!(e.energy_scaling(), 4.0 + e.alpha() - e.beta());
    }
    assert!(Exponents::from_p(4.0).unwrap().subcritical());
    assert!(Exponents::new(3.0, 6.0).unwrap().subcritical());
    let msg = Exponents::new(6.0, 14.0).unwrap_err().to_string();
    assert!(msg.contains("alpha 6") && msg.contains("beta 14"), "{msg}");
}

#[test]
fn integer_powers_agree_with_powf() {
    let mut rng = SplitMix64::new(11);
    for _ in 0..200 {
        let x = rng.uniform(0.01, 5.0);
        for n in 0..20u32 {
            let a = ipow(x, n);
            let b = x.powf(n as f64);
            assert!((a - b).abs() <= 1e-14 * b, "{x}^{n}");
        }
    }
    // The integer and the general path agree: (6, 12) vs a nearby non-integer pair, by
    // continuity in the exponent.
    let x = [0.3, -0.2, 0.5];
    let n = [0.0, 0.6, 0.8];
    let y = [-0.4, 0.1, 0.2];
    let k_int = kernel(&Exponents::default(), x, n, y);
    let k_near = kernel(&Exponents::new(6.0 + 1e-9, 12.0 + 1e-9).unwrap(), x, n, y);
    assert!((k_int - k_near).abs() <= 1e-7 * k_int);
}

#[test]
fn kernel_gradient_matches_finite_differences_and_is_continuous() {
    let h = 1e-6;
    let mut rng = SplitMix64::new(5);
    for e in exponent_sets() {
        for _ in 0..20 {
            let x: [f64; 3] = std::array::from_fn(|_| rng.uniform(-1.0, 1.0));
            let y: [f64; 3] = std::array::from_fn(|_| rng.uniform(-1.0, 1.0));
            let n: [f64; 3] = std::array::from_fn(|_| rng.uniform(-1.0, 1.0));
            let k = kernel_with_gradient(&e, x, n, y);
            assert_eq!(k.value, kernel(&e, x, n, y));
            let scale = k
                .d_v
                .iter()
                .chain(&k.d_n)
                .map(|v| v.abs())
                .fold(0.0, f64::max);
            for i in 0..3 {
                let mut xp = x;
                let mut xm = x;
                xp[i] += h;
                xm[i] -= h;
                let fd = (kernel(&e, xp, n, y) - kernel(&e, xm, n, y)) / (2.0 * h);
                assert!((fd - k.d_v[i]).abs() <= 1e-6 * scale, "dX {i}");
                let mut yp = y;
                let mut ym = y;
                yp[i] += h;
                ym[i] -= h;
                let fd = (kernel(&e, x, n, yp) - kernel(&e, x, n, ym)) / (2.0 * h);
                assert!((fd + k.d_v[i]).abs() <= 1e-6 * scale, "dY {i}");
                let mut np = n;
                let mut nm = n;
                np[i] += h;
                nm[i] -= h;
                let fd = (kernel(&e, x, np, y) - kernel(&e, x, nm, y)) / (2.0 * h);
                assert!((fd - k.d_n[i]).abs() <= 1e-6 * scale, "dN {i}");
            }
        }
        // ν = 0: value and the ν-terms vanish; tiny ν of both signs gives tiny gradients.
        let x = [1.0, 0.5, 0.0];
        let n = [0.0, 0.0, 1.0];
        let k0 = kernel_with_gradient(&e, x, n, [0.0; 3]);
        assert_eq!(k0.value, 0.0);
        assert!(k0.d_n.iter().all(|&v| v == 0.0));
        assert!(k0.d_v.iter().all(|&v| v == 0.0));
        for sgn in [1.0, -1.0] {
            let k = kernel_with_gradient(&e, [1.0, 0.5, sgn * 1e-9], n, [0.0; 3]);
            assert!(k.d_n.iter().chain(&k.d_v).all(|v| v.abs() < 1e-6), "{k:?}");
        }
    }
}

fn exact(s: &Surface, p: &[[f64; 3]], e: Exponents, threads: usize) -> (f64, Vec<[f64; 3]>) {
    let r = energy(s, p, e, &EnergyEval::Exact, &[], Some(threads)).unwrap();
    assert!(r.penalties.is_empty());
    assert_eq!(r.nodes_lumped, 0);
    let nf = s.n_faces() as u64;
    assert_eq!(r.pairs_exact, nf * (nf - 1));
    assert_eq!(r.subcritical, e.subcritical());
    (r.energy, r.differential)
}

#[test]
fn differential_matches_finite_differences() {
    // Every (mesh, exponents) case in parallel on its own pool.
    let cases: Vec<_> = test_meshes()
        .into_iter()
        .flat_map(|m| exponent_sets().into_iter().map(move |e| (m.clone(), e)))
        .collect();
    std::thread::scope(|scope| {
        let handles: Vec<_> = cases
            .iter()
            .map(|((name, s, pos), e)| {
                scope.spawn(move || {
                    let (_, g) = exact(s, pos, *e, 2);
                    let pool = make_pool(Some(2)).unwrap();
                    let h = 1e-6 * mean_edge_length(s, pos);
                    let val =
                        |p: &[[f64; 3]]| energy_value(&pool, &FaceGeometry::new(s, p).unwrap(), e);
                    let mut err2 = 0.0;
                    let mut norm2 = 0.0;
                    let mut p = pos.clone();
                    for v in 0..p.len() {
                        for k in 0..3 {
                            let x0 = p[v][k];
                            p[v][k] = x0 + h;
                            let ep = val(&p);
                            p[v][k] = x0 - h;
                            let em = val(&p);
                            p[v][k] = x0;
                            let fd = (ep - em) / (2.0 * h);
                            err2 += (fd - g[v][k]).powi(2);
                            norm2 += g[v][k].powi(2);
                        }
                    }
                    let rel = (err2 / norm2).sqrt();
                    eprintln!("FD {name} ({}, {}): {rel:.2e}", e.alpha(), e.beta());
                    assert!(
                        rel < 1e-6,
                        "{name}, ({}, {}): relative FD error {rel:e}",
                        e.alpha(),
                        e.beta()
                    );
                })
            })
            .collect();
        for h in handles {
            h.join().unwrap();
        }
    });
}

fn rotate(p: [f64; 3], r: &[[f64; 3]; 3], t: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| r[i][0] * p[0] + r[i][1] * p[1] + r[i][2] * p[2] + t[i])
}

#[test]
fn invariances_scaling_and_gradient_identities() {
    // A rotation (Rodrigues) about a generic axis.
    let axis = {
        let a: [f64; 3] = [0.3, -0.5, 0.81];
        let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        [a[0] / l, a[1] / l, a[2] / l]
    };
    let th: f64 = 0.83;
    let (c, s) = (th.cos(), th.sin());
    let k = [
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ];
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let kk: f64 = (0..3).map(|m| k[i][m] * k[m][j]).sum();
            r[i][j] = if i == j { 1.0 } else { 0.0 } + s * k[i][j] + (1.0 - c) * kk;
        }
    }
    let t = [0.4, -1.3, 2.2];
    for (name, surf, pos) in test_meshes() {
        for e in exponent_sets() {
            let (e0, g0) = exact(&surf, &pos, e, 4);
            // Rigid motion.
            let moved: Vec<_> = pos.iter().map(|&p| rotate(p, &r, t)).collect();
            let (e1, g1) = exact(&surf, &moved, e, 4);
            assert!((e1 - e0).abs() <= 1e-12 * e0, "{name}: rigid {e0} {e1}");
            let gs = g0.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
            for v in 0..pos.len() {
                let rg = rotate(g0[v], &r, [0.0; 3]);
                for i in 0..3 {
                    assert!(
                        (g1[v][i] - rg[i]).abs() <= 1e-10 * gs,
                        "{name}: rotated gradient"
                    );
                }
            }
            // Scaling.
            let cs = 1.7;
            let scaled: Vec<_> = pos.iter().map(|p| p.map(|x| cs * x)).collect();
            let (e2, _) = exact(&surf, &scaled, e, 4);
            let want = cs.powf(e.energy_scaling()) * e0;
            assert!(
                (e2 - want).abs() <= 1e-12 * want,
                "{name}: scaling {e2} {want}"
            );
            // Gradient identities.
            let mut sum = [0.0; 3];
            let mut torque = [0.0; 3];
            let mut euler = 0.0;
            let mut mag = 0.0;
            for v in 0..pos.len() {
                let (p, g) = (pos[v], g0[v]);
                for i in 0..3 {
                    sum[i] += g[i];
                }
                torque[0] += p[1] * g[2] - p[2] * g[1];
                torque[1] += p[2] * g[0] - p[0] * g[2];
                torque[2] += p[0] * g[1] - p[1] * g[0];
                euler += p[0] * g[0] + p[1] * g[1] + p[2] * g[2];
                let pn = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                let gn = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                mag += pn * gn + gn;
            }
            for i in 0..3 {
                assert!(sum[i].abs() <= 1e-12 * mag, "{name}: Σg = {sum:?}");
                assert!(torque[i].abs() <= 1e-12 * mag, "{name}: Σf×g = {torque:?}");
            }
            let want = e.energy_scaling() * e0;
            assert!(
                (euler - want).abs() <= 1e-12 * mag,
                "{name}: Euler {euler} {want}"
            );
        }
    }
}

#[test]
fn face_orientation_does_not_matter() {
    for (name, surf, pos) in test_meshes() {
        let mut rng = SplitMix64::new(9);
        let faces: Vec<[u32; 3]> = surf
            .faces()
            .iter()
            .map(|t| match (rng.next_u64() % 3, rng.next_u64() % 2) {
                (r, 0) => [
                    t[r as usize],
                    t[(r as usize + 1) % 3],
                    t[(r as usize + 2) % 3],
                ],
                (r, _) => [
                    t[r as usize],
                    t[(r as usize + 2) % 3],
                    t[(r as usize + 1) % 3],
                ],
            })
            .collect();
        let flipped = Surface::new(pos.len(), faces).unwrap();
        let e = Exponents::default();
        let (e0, g0) = exact(&surf, &pos, e, 3);
        let (e1, g1) = exact(&flipped, &pos, e, 3);
        assert!((e1 - e0).abs() <= 1e-13 * e0, "{name}: {e0} {e1}");
        let gs = g0.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
        for (a, b) in g0.iter().flatten().zip(g1.iter().flatten()) {
            assert!((a - b).abs() <= 1e-11 * gs, "{name}: gradient");
        }
    }
}

#[test]
fn two_triangles_by_hand() {
    let s = Surface::new(6, vec![[0, 1, 2], [3, 4, 5]]).unwrap();
    let e = Exponents::default();
    // Parallel triangles, X_1 − X_0 = (1, 0, 1): each direction |ν|^6 / r^12 = 1/2^6.
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [2.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    let (en, _) = exact(&s, &p, e, 1);
    assert_eq!(en, 2.0 * 0.25 / 64.0);
    // Perpendicular: T1 in the xz plane, X_0 − X_1 = (−2, 1/3, −1/3); |ν| = 1/3 both ways,
    // r² = 38/9.
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [2.0, 0.0, 1.0],
    ];
    let (en, _) = exact(&s, &p, e, 1);
    let want = 2.0 * 0.25 * (1.0f64 / 3.0).powi(6) * (9.0f64 / 38.0).powi(6);
    assert!((en - want).abs() <= 1e-14 * want, "{en} {want}");
    // Non-integer exponents.
    let e = Exponents::new(4.5, 9.5).unwrap();
    let (en, _) = exact(&s, &p, e, 1);
    let want = 2.0 * 0.25 * (1.0f64 / 3.0).powf(4.5) * (38.0f64 / 9.0).powf(-4.75);
    assert!((en - want).abs() <= 1e-14 * want, "{en} {want}");
}

/// The energy by the plain formula: per pair, cross products and `powf`.
fn naive_energy(s: &Surface, p: &[[f64; 3]], e: Exponents) -> f64 {
    let face = |t: &[u32; 3]| {
        let (a, b, c) = (p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]);
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let x = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        (0.5 * l, [n[0] / l, n[1] / l, n[2] / l], x)
    };
    let mut total = 0.0;
    for (i, ti) in s.faces().iter().enumerate() {
        for (j, tj) in s.faces().iter().enumerate() {
            if i == j {
                continue;
            }
            let (ai, ni, xi) = face(ti);
            let (aj, _, xj) = face(tj);
            let v = [xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]];
            let nu = (ni[0] * v[0] + ni[1] * v[1] + ni[2] * v[2]).abs();
            let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            total += nu.powf(e.alpha()) / r.powf(e.beta()) * ai * aj;
        }
    }
    total
}

#[test]
fn agrees_with_a_plain_double_loop() {
    for (name, surf, pos) in test_meshes() {
        for e in exponent_sets() {
            let (en, _) = exact(&surf, &pos, e, 2);
            let naive = naive_energy(&surf, &pos, e);
            assert!((en - naive).abs() <= 1e-12 * naive, "{name}: {en} {naive}");
        }
    }
}

#[test]
fn bitwise_identical_across_thread_counts() {
    let m = icosphere(3).unwrap();
    let mut pos = m.positions.clone();
    perturb(&mut pos, [0.01; 3], 17);
    let s = Surface::new(pos.len(), m.faces.clone()).unwrap();
    for e in [Exponents::default(), Exponents::new(4.5, 9.5).unwrap()] {
        let (e1, g1) = exact(&s, &pos, e, 1);
        for threads in [3, 8] {
            let (e2, g2) = exact(&s, &pos, e, threads);
            assert_eq!(e1.to_bits(), e2.to_bits(), "{threads} threads");
            for (a, b) in g1.iter().flatten().zip(g2.iter().flatten()) {
                assert_eq!(a.to_bits(), b.to_bits(), "{threads} threads");
            }
            let pool = make_pool(Some(threads)).unwrap();
            let v = energy_value(&pool, &FaceGeometry::new(&s, &pos).unwrap(), &e);
            assert!((v - e1).abs() <= 1e-13 * e1);
        }
    }
}

#[test]
fn errors_of_the_energy() {
    // Coincident barycenters of two distinct faces.
    let p = vec![
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [2.0, 2.0, 0.0],
        [-1.0, 2.0, 0.0],
        [2.0, -1.0, 0.0],
    ];
    let s = Surface::new(6, vec![[0, 1, 2], [3, 4, 5]]).unwrap();
    let e = energy(
        &s,
        &p,
        Exponents::default(),
        &EnergyEval::Exact,
        &[],
        Some(2),
    )
    .unwrap_err();
    assert_eq!(e, RsError::NonFiniteEnergy { faces: [0, 1] });
    assert!(e.to_string().contains("faces (0, 1)"));
    let p2: Vec<_> = p.iter().map(|q| [q[0], q[1], q[2] + 0.0]).collect();
    assert!(matches!(
        energy(
            &s,
            &p2[..5],
            Exponents::default(),
            &EnergyEval::Exact,
            &[],
            None
        ),
        Err(RsError::PositionCountMismatch { .. })
    ));
    let mut p3 = p.clone();
    p3[4][2] = 1.0;
    assert!(matches!(
        energy(
            &s,
            &p3,
            Exponents::default(),
            &EnergyEval::Exact,
            &[],
            Some(0)
        ),
        Err(RsError::InvalidOption {
            name: "threads",
            ..
        })
    ));
    let opts = ddg_repulsive_surfaces::BhOptions {
        theta: -1.0,
        leaf_size: 8,
    };
    assert!(matches!(
        energy(
            &s,
            &p3,
            Exponents::default(),
            &EnergyEval::BarnesHut(opts),
            &[],
            Some(1)
        ),
        Err(RsError::InvalidOption { name: "theta", .. })
    ));
}
