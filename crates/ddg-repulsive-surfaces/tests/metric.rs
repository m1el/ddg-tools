//! The dense metric (`SPEC.md` §4, §13.3): `D_f`, `L_τ`, `B`, `B⁰`, `A`, the kernel-matrix form,
//! the comparison metrics, and agreement with the oracle's `df` and `metric` fixtures.

#![allow(clippy::needless_range_loop)]

#[path = "wp2_fixture/mod.rs"]
mod fixture;

use ddg_core::generators::{bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::{cotan_stiffness, lumped_mass, mean_edge_length};
use ddg_repulsive_surfaces::metric::dense::{
    comparison_metric, dense_metric, fractional_kernel, fractional_operator, high_order_kernel,
    high_order_operator, kernel_matrix_apply, kernel_matrix_form, low_order_kernel,
    low_order_operator, pair_sum_matrix, to_dense, willmore_matrix,
};
use ddg_repulsive_surfaces::metric::{
    Channel, DenseMatrix, FaceDiagonal, FaceOperators, L2Kind, Metric, SobolevWeights,
    df_coefficients,
};
use ddg_repulsive_surfaces::{Exponents, RsError, Surface};

/// A Möbius band with `n` quads.
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

/// Perturbed test meshes: icosphere, torus, disk, Möbius band, two components.
fn meshes() -> Vec<(&'static str, Surface, Vec<[f64; 3]>)> {
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

fn max_abs_diff(a: &DenseMatrix, b: &DenseMatrix) -> f64 {
    a.as_slice()
        .iter()
        .zip(b.as_slice())
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
}

/// Eigenvalues (increasing) of a symmetric matrix.
fn eigenvalues(a: &DenseMatrix) -> Vec<f64> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::evd::{ComputeEigenvectors, self_adjoint_evd, self_adjoint_evd_scratch};
    let n = a.n();
    let m = faer::Mat::<f64>::from_fn(n, n, |i, j| 0.5 * (a.get(i, j) + a.get(j, i)));
    let mut s = faer::diag::Diag::<f64>::zeros(n);
    let par = faer::Par::Seq;
    let mut mem = MemBuffer::new(self_adjoint_evd_scratch::<f64>(
        n,
        ComputeEigenvectors::No,
        par,
        Default::default(),
    ));
    self_adjoint_evd(
        m.as_ref(),
        s.as_mut(),
        None,
        par,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .unwrap();
    let mut ev: Vec<f64> = (0..n).map(|i| s[i]).collect();
    ev.sort_by(f64::total_cmp);
    ev
}

#[test]
fn df_is_the_p1_gradient() {
    // On a linear function u = c·f + d, D_f u(S) is the tangential part of c: (I − N Nᵀ) c.
    for (name, s, p) in meshes() {
        let ops = FaceOperators::new(&s, &p).unwrap();
        let c = [0.3, -1.1, 0.7];
        let u: Vec<f64> = p
            .iter()
            .map(|q| c[0] * q[0] + c[1] * q[1] + c[2] * q[2] + 2.5)
            .collect();
        let g = ops.gradient(&u);
        let ones = ops.gradient(&vec![1.0; p.len()]);
        let mut err: f64 = 0.0;
        for (f, gf) in g.iter().enumerate() {
            let n = ops.normals[f];
            let cn = c[0] * n[0] + c[1] * n[1] + c[2] * n[2];
            for k in 0..3 {
                err = err.max((gf[k] - (c[k] - cn * n[k])).abs());
                assert!(ones[f][k].abs() < 1e-12, "{name}: D_f 1 != 0");
            }
        }
        assert!(err < 1e-12, "{name}: D_f on a linear function off by {err:e}");
        // D_f is orientation free: flipping a face gives the same coefficients, permuted.
        for (f, t) in s.faces().iter().enumerate() {
            let q = [p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]];
            let g0 = df_coefficients(q);
            let g1 = df_coefficients([q[0], q[2], q[1]]);
            for k in 0..3 {
                assert!((g0[1][k] - g1[2][k]).abs() < 1e-12, "{name}: face {f} flip");
                assert!((g0[0][k] - g1[0][k]).abs() < 1e-12, "{name}: face {f} flip");
            }
        }
    }
}

#[test]
fn kernel_matrix_form_equals_the_pair_sums() {
    let pool = make_pool(Some(4)).unwrap();
    for e in [Exponents::default(), Exponents::new(3.0, 7.5).unwrap()] {
        for (name, s, p) in meshes() {
            let ops = FaceOperators::new(&s, &p).unwrap();
            let rho2 = |a: usize, b: usize| {
                let x = ops.barycenters[a];
                let y = ops.barycenters[b];
                (0..3).map(|k| (x[k] - y[k]).powi(2)).sum::<f64>()
            };
            let ss = e.s();
            // L_τ, B and the ordered (asymmetric) B⁰ kernel of Eq. 13.
            let tau = 0.6;
            let lt = fractional_operator(&pool, &ops, tau).unwrap();
            let lt_ref = pair_sum_matrix(&pool, &ops, Channel::Average, &|a, b| {
                rho2(a, b).powf(-(tau + 1.0))
            })
            .unwrap();
            let b = high_order_operator(&pool, &ops, &e).unwrap();
            let b_ref = pair_sum_matrix(&pool, &ops, Channel::Gradient, &|a, b| {
                rho2(a, b).powf(-ss)
            })
            .unwrap();
            let b0 = low_order_operator(&pool, &ops, &e).unwrap();
            let b0_ref = pair_sum_matrix(&pool, &ops, Channel::Average, &|a, b| {
                let x = ops.barycenters[a];
                let y = ops.barycenters[b];
                let n = ops.normals[a];
                let nu: f64 = (0..3).map(|k| n[k] * (x[k] - y[k])).sum();
                let r2 = rho2(a, b);
                nu * nu / (r2 * r2) * r2.powf(-ss)
            })
            .unwrap();
            for (what, m, r) in [("L_tau", &lt, &lt_ref), ("B", &b, &b_ref), ("B0", &b0, &b0_ref)]
            {
                let err = max_abs_diff(m, r) / r.max_abs();
                assert!(err < 1e-12, "{name} {what}: kernel form vs pair sum {err:e}");
                assert!(m.symmetry_defect() <= 1e-13 * m.max_abs(), "{name} {what} symmetric");
            }
            // The 3-column product equals the matrix times the columns.
            let h = high_order_kernel(&pool, &ops, &e).unwrap();
            let u: Vec<[f64; 3]> = (0..s.n_vertices())
                .map(|v| [(v as f64).sin(), (0.3 * v as f64).cos(), 0.1 * v as f64])
                .collect();
            let prod = kernel_matrix_apply(&pool, &ops, &h, Channel::Gradient, &u).unwrap();
            let mut err: f64 = 0.0;
            for v in 0..s.n_vertices() {
                for k in 0..3 {
                    let x: f64 = (0..s.n_vertices()).map(|w| b.get(v, w) * u[w][k]).sum();
                    err = err.max((x - prod[v][k]).abs());
                }
            }
            assert!(err < 1e-11 * b.max_abs(), "{name}: apply vs matrix {err:e}");
        }
    }
}

#[test]
fn kernel_matrix_form_rejects_bad_kernels() {
    let pool = make_pool(Some(1)).unwrap();
    let (_, s, p) = meshes().remove(0);
    let ops = FaceOperators::new(&s, &p).unwrap();
    let mut h = low_order_kernel(&pool, &ops, &Exponents::default()).unwrap();
    let nf = s.n_faces();
    h.as_mut_slice()[1] += 1.0; // (0, 1) only
    assert!(matches!(
        kernel_matrix_form(&pool, &ops, &h, Channel::Average),
        Err(RsError::InvalidOption { .. })
    ));
    let mut h = fractional_kernel(&pool, &ops, 0.5).unwrap();
    h.as_mut_slice()[nf + 1] = 1.0; // diagonal (1, 1)
    assert!(matches!(
        kernel_matrix_form(&pool, &ops, &h, Channel::Average),
        Err(RsError::InvalidOption { .. })
    ));
}

#[test]
fn metrics_are_symmetric_psd_with_the_expected_kernel() {
    let pool = make_pool(Some(4)).unwrap();
    for (name, s, p) in meshes() {
        let nc = s.n_components();
        let e = Exponents::default();
        let dm = dense_metric(&pool, &s, &p, &e, &SobolevWeights::default(), None).unwrap();
        let dm_hi = dense_metric(
            &pool,
            &s,
            &p,
            &e,
            &SobolevWeights {
                low_order: false,
                ..Default::default()
            },
            None,
        )
        .unwrap();
        let ops = FaceOperators::new(&s, &p).unwrap();
        let lt = fractional_operator(&pool, &ops, 2.0 - e.s()).unwrap();
        // (matrix, expected kernel dimension or None (only PSD and constants), label)
        let cases: Vec<(&str, &DenseMatrix, Option<usize>)> = vec![
            // B: D_f u = 0 iff u is constant on every component.
            ("B", &dm.b, Some(nc)),
            // B⁰ and L_τ couple all face pairs, also across components: global constants.
            ("B0", &dm.b0, None),
            ("L_{2-s}", &lt, None),
            // A = B + B⁰: the per-component constants of B that B⁰ also kills: global constants.
            ("A", &dm.a, Some(1)),
            ("A (low_order off)", &dm_hi.a, Some(nc)),
        ];
        for (what, m, kdim) in cases {
            assert!(m.symmetry_defect() <= 1e-13 * m.max_abs(), "{name} {what}");
            let ev = eigenvalues(m);
            let top = ev[ev.len() - 1];
            assert!(ev[0] > -1e-11 * top, "{name} {what}: not PSD, {:e}", ev[0]);
            let ones: Vec<f64> = vec![1.0; s.n_vertices()];
            let r = m.mul_vec(&ones);
            let rn = r.iter().fold(0.0f64, |a, x| a.max(x.abs()));
            assert!(rn < 1e-11 * top, "{name} {what}: constants not in the kernel {rn:e}");
            if let Some(k) = kdim {
                assert!(
                    ev[k - 1] < 1e-11 * top,
                    "{name} {what}: kernel smaller than {k}: {:e}",
                    ev[k - 1] / top
                );
                assert!(
                    ev[k] > 1e-9 * top,
                    "{name} {what}: kernel larger than {k}: {:e}",
                    ev[k] / top
                );
                eprintln!(
                    "{name:>15} {what:<18} kernel {k}, smallest nonzero eigenvalue / largest = {:.3e}",
                    ev[k] / top
                );
            }
        }
    }
}

#[test]
fn metric_weights_and_obstacle_hook() {
    let pool = make_pool(Some(2)).unwrap();
    let (_, s, p) = meshes().remove(1);
    let e = Exponents::new(4.0, 8.5).unwrap();
    let w = SobolevWeights {
        tp_weight: 2.0,
        low_order: true,
        willmore_weight: 0.3,
    };
    let dm = dense_metric(&pool, &s, &p, &e, &w, None).unwrap();
    let wmw = to_dense(dm.wmw.as_ref().unwrap());
    let n = s.n_vertices();
    let mut err: f64 = 0.0;
    for i in 0..n * n {
        let x = 2.0 * (dm.b.as_slice()[i] + dm.b0.as_slice()[i]) + 0.3 * wmw.as_slice()[i];
        err = err.max((x - dm.a.as_slice()[i]).abs());
    }
    assert!(err <= 1e-12 * dm.a.max_abs());
    // W M⁻¹ W against the dense product.
    let wd = to_dense(&dm.w);
    let mut err: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x: f64 = (0..n).map(|k| wd.get(i, k) * wd.get(k, j) / dm.mass[k]).sum();
            err = err.max((x - wmw.get(i, j)).abs());
        }
    }
    assert!(err < 1e-12 * wmw.max_abs());
    // The obstacle hook adds 2 Xᵀ diag(d/a) X: a kernel-matrix form whose H is zero but whose
    // H a is d.
    let ops = FaceOperators::new(&s, &p).unwrap();
    let hi: Vec<f64> = (0..s.n_faces()).map(|f| 0.1 + 0.01 * f as f64).collect();
    let lo: Vec<f64> = (0..s.n_faces()).map(|f| 0.2 + 0.03 * f as f64).collect();
    let o = FaceDiagonal {
        hi: hi.clone(),
        lo: lo.clone(),
    };
    let dmo = dense_metric(&pool, &s, &p, &e, &w, Some(&o)).unwrap();
    let mut q = vec![0.0; n];
    q[3] = 1.0;
    q[7] = -0.5;
    let quad = |m: &DenseMatrix| -> f64 {
        let r = m.mul_vec(&q);
        r.iter().zip(&q).map(|(a, b)| a * b).sum()
    };
    let g = ops.gradient(&q);
    let ub = ops.average(&q);
    let expect: f64 = (0..s.n_faces())
        .map(|f| {
            let a = ops.areas[f];
            2.0 * hi[f] * a * (g[f][0].powi(2) + g[f][1].powi(2) + g[f][2].powi(2))
                + 2.0 * lo[f] * a * ub[f] * ub[f]
        })
        .sum();
    let got = quad(&dmo.a) - quad(&dm.a);
    assert!((got - expect).abs() < 1e-12 * expect.abs(), "{got} vs {expect}");
    // Wrong length.
    let bad = FaceDiagonal {
        hi: vec![1.0; 3],
        lo: Vec::new(),
    };
    assert!(dense_metric(&pool, &s, &p, &e, &w, Some(&bad)).is_err());
    // Invalid weights.
    let badw = SobolevWeights {
        tp_weight: 0.0,
        ..Default::default()
    };
    assert!(matches!(
        dense_metric(&pool, &s, &p, &e, &badw, None),
        Err(RsError::InvalidOption { .. })
    ));
}

#[test]
fn comparison_metrics() {
    for (name, s, p) in meshes() {
        let n = s.n_vertices();
        let m = lumped_mass(&s, &p).unwrap();
        let w = cotan_stiffness(&s, &p).unwrap();
        let l2 = to_dense(&comparison_metric(&s, &p, Metric::L2(L2Kind::Mass)).unwrap());
        let id = to_dense(&comparison_metric(&s, &p, Metric::L2(L2Kind::Identity)).unwrap());
        let h1 = to_dense(&comparison_metric(&s, &p, Metric::H1).unwrap());
        let h2 = to_dense(&comparison_metric(&s, &p, Metric::H2).unwrap());
        for i in 0..n {
            assert_eq!(l2.get(i, i), m[i]);
            assert_eq!(id.get(i, i), 1.0);
        }
        assert_eq!(max_abs_diff(&h1, &to_dense(&w)), 0.0);
        assert_eq!(
            max_abs_diff(&h2, &to_dense(&willmore_matrix(&w, &m).unwrap())),
            0.0
        );
        let nc = s.n_components();
        for (what, mm) in [("H1", &h1), ("H2", &h2)] {
            let ev = eigenvalues(mm);
            let top = ev[n - 1];
            assert!(ev[0] > -1e-11 * top, "{name} {what} PSD");
            assert!(ev[nc - 1] < 1e-11 * top, "{name} {what} kernel");
            assert!(ev[nc] > 1e-9 * top, "{name} {what} kernel = components");
        }
        assert!(matches!(
            comparison_metric(&s, &p, Metric::Sobolev { low_order: true }),
            Err(RsError::InvalidOption { .. })
        ));
    }
}

#[test]
fn dense_metric_is_bitwise_independent_of_threads() {
    let (_, s, p) = meshes().remove(1);
    let e = Exponents::default();
    let w = SobolevWeights::default();
    let a1 = dense_metric(&make_pool(Some(1)).unwrap(), &s, &p, &e, &w, None).unwrap();
    for t in [3, 8] {
        let at = dense_metric(&make_pool(Some(t)).unwrap(), &s, &p, &e, &w, None).unwrap();
        assert_eq!(a1, at, "threads {t}");
    }
}

#[test]
fn coincident_barycenters_are_an_error() {
    let pool = make_pool(Some(1)).unwrap();
    // Two triangles with the same barycenter (one rotated within its plane).
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0 / 3.0, 2.0 / 3.0, 0.0],
        [-1.0 / 3.0, 2.0 / 3.0, 0.0],
        [2.0 / 3.0, -1.0 / 3.0, 0.0],
    ];
    let s = Surface::new(6, vec![[0, 1, 2], [3, 4, 5]]).unwrap();
    let ops = FaceOperators::new(&s, &p).unwrap();
    assert!(matches!(
        high_order_kernel(&pool, &ops, &Exponents::default()),
        Err(RsError::NonFiniteEnergy { faces: [0, 1] })
    ));
}

// ---- the oracle's golden fixtures (relative 1e-9) ----

#[test]
fn oracle_df_fixtures() {
    let names = fixture::cases("df");
    assert!(names.len() >= 4);
    for name in names {
        let c = fixture::read_case(&name);
        let o = fixture::read_output(&name);
        let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
        let ops = FaceOperators::new(&s, &c.positions).unwrap();
        let mut got = Vec::new();
        let mut want = Vec::new();
        for l in o.get("df") {
            let k: usize = l[0].parse().unwrap();
            let f: usize = l[1].parse().unwrap();
            let u = &c.u[k];
            let ub = ops.average(u)[f];
            let g = ops.gradient(u)[f];
            got.extend([ub, g[0], g[1], g[2]]);
            want.extend(l[2..6].iter().map(|x| fixture::hex(x)));
        }
        let err = fixture::rel_err(&got, &want, 1e-300);
        eprintln!("{name}: df relative error {err:.2e}");
        assert!(err < 1e-9, "{name}: {err:e}");
    }
}

#[test]
fn oracle_metric_fixtures() {
    let names = fixture::cases("metric");
    assert!(names.len() >= 6);
    let pool = make_pool(Some(4)).unwrap();
    for name in names {
        let c = fixture::read_case(&name);
        let o = fixture::read_output(&name);
        let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
        let (al, be) = c.exponents.unwrap_or((6.0, 12.0));
        let e = Exponents::new(al, be).unwrap();
        let (wtp, clo, ww) = c.weights.unwrap_or((1.0, 1.0, 0.0));
        let w = SobolevWeights {
            tp_weight: wtp,
            low_order: clo != 0.0,
            willmore_weight: ww,
        };
        let dm = dense_metric(&pool, &s, &c.positions, &e, &w, None).unwrap();
        let flat = |m: &DenseMatrix| m.as_slice().to_vec();
        let flat_o = |k: &str| o.matrix(k).into_iter().flatten().collect::<Vec<f64>>();
        let mut checks: Vec<(&str, Vec<f64>, Vec<f64>)> = vec![("A", flat(&dm.a), flat_o("A"))];
        if !c.matrices_a_only {
            let ops = FaceOperators::new(&s, &c.positions).unwrap();
            let tau = c.tau.unwrap_or(0.75);
            let lt = fractional_operator(&pool, &ops, tau).unwrap();
            let wmw = willmore_matrix(&dm.w, &dm.mass).unwrap();
            checks.push(("L", flat(&lt), flat_o("L")));
            checks.push(("B", flat(&dm.b), flat_o("B")));
            checks.push(("B0", flat(&dm.b0), flat_o("B0")));
            checks.push(("W", flat(&to_dense(&dm.w)), flat_o("W")));
            checks.push(("WMW", flat(&to_dense(&wmw)), flat_o("WMW")));
            let m: Vec<f64> = o.get("M")[0].iter().map(|x| fixture::hex(x)).collect();
            checks.push(("M", dm.mass.clone(), m));
            assert!((fixture::hex(&o.get("s")[0][0]) - e.s()).abs() < 1e-15);
        }
        for (what, got, want) in checks {
            let err = fixture::rel_err(&got, &want, 1e-300);
            eprintln!("{name}: {what} relative error {err:.2e}");
            assert!(err < 1e-9, "{name} {what}: {err:e}");
        }
    }
}
