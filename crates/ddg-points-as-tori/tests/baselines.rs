//! Baselines and mesh ground truth (SPEC.md §9, WP5): point areas, direct and fast winding
//! numbers, the regularized winding number, SHC / SSPD / LSE with their gradients, the exact mesh
//! SDF against analytic SDFs, determinism across thread counts, and error cases.
//!
//! The ignored test `sanity_grid_errors` is the WP5 sanity run (64³ grid errors of PAT and the
//! baselines against the mesh truth on one golden cloud; see `NOTES-WP5.md`).

mod common;

use common::{ellipsoid_cloud, lcg};
use ddg_points_as_tori::baselines::fwn::triangle_winding;
use ddg_points_as_tori::baselines::{
    FastWinding, FwnOptions, KernelField, KernelMethod, KernelOptions, PointAreaOptions, point_areas, regularization_s,
    regularized_winding_number, voronoi_cell_area, winding_number, winding_numbers,
};
use ddg_points_as_tori::truth::{MeshTruth, SignMode, TruthOptions, check_mesh, closest_point_on_triangle};
use ddg_points_as_tori::{LambdaRule, Neighborhood, NormalPolicy, PatError, PointCloud};
use std::f64::consts::PI;

fn fib_sphere(n: usize, r: f64) -> PointCloud {
    let mut p = Vec::new();
    let mut nr = Vec::new();
    for i in 0..n {
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
        let s = (1.0 - z * z).sqrt();
        let th = PI * (1.0 + 5f64.sqrt()) * (i as f64 + 0.5);
        let u = [th.cos() * s, th.sin() * s, z];
        p.push([r * u[0], r * u[1], r * u[2]]);
        nr.push(u);
    }
    PointCloud::new(p, nr, NormalPolicy::RequireUnit { tol: 1e-9 }).unwrap().0
}

fn rand_queries(n: usize, seed: u64, half: f64) -> Vec<[f64; 3]> {
    let mut s = seed;
    (0..n).map(|_| [half * (2.0 * lcg(&mut s) - 1.0), half * (2.0 * lcg(&mut s) - 1.0), half * (2.0 * lcg(&mut s) - 1.0)]).collect()
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

// ---------------------------------------------------------------- meshes

fn icosphere(r: f64, levels: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0], [1.0, t, 0.0], [-1.0, -t, 0.0], [1.0, -t, 0.0],
        [0.0, -1.0, t], [0.0, 1.0, t], [0.0, -1.0, -t], [0.0, 1.0, -t],
        [t, 0.0, -1.0], [t, 0.0, 1.0], [-t, 0.0, -1.0], [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut mid = std::collections::BTreeMap::new();
        let mut nf = Vec::new();
        let mut m = |a: u32, b: u32, v: &mut Vec<[f64; 3]>| -> u32 {
            *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let (p, q) = (v[a as usize], v[b as usize]);
                v.push([(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0, (p[2] + q[2]) / 2.0]);
                (v.len() - 1) as u32
            })
        };
        for t in &f {
            let a = m(t[0], t[1], &mut v);
            let b = m(t[1], t[2], &mut v);
            let c = m(t[2], t[0], &mut v);
            nf.extend([[t[0], a, c], [t[1], b, a], [t[2], c, b], [a, b, c]]);
        }
        f = nf;
    }
    for p in &mut v {
        let l = norm(*p);
        *p = [r * p[0] / l, r * p[1] / l, r * p[2] / l];
    }
    (v, f)
}

fn torus_mesh(big: f64, small: f64, nu: usize, nv: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut v = Vec::new();
    for i in 0..nu {
        let u = 2.0 * PI * i as f64 / nu as f64;
        for j in 0..nv {
            let w = 2.0 * PI * j as f64 / nv as f64;
            v.push([(big + small * w.cos()) * u.cos(), (big + small * w.cos()) * u.sin(), small * w.sin()]);
        }
    }
    let id = |i: usize, j: usize| ((i % nu) * nv + (j % nv)) as u32;
    let mut f = Vec::new();
    for i in 0..nu {
        for j in 0..nv {
            // outward: (∂u × ∂v) points outward
            f.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            f.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    (v, f)
}

fn cube(h: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut v = Vec::new();
    for i in 0..8 {
        v.push([if i & 1 != 0 { h } else { -h }, if i & 2 != 0 { h } else { -h }, if i & 4 != 0 { h } else { -h }]);
    }
    // quads with outward orientation, then split
    let quads = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
    let mut f = Vec::new();
    for q in quads {
        f.push([q[0], q[1], q[2]]);
        f.push([q[0], q[2], q[3]]);
    }
    (v, f)
}

fn box_sdf(x: [f64; 3], h: f64) -> f64 {
    let q = [x[0].abs() - h, x[1].abs() - h, x[2].abs() - h];
    let o = norm([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]);
    o + q[0].max(q[1]).max(q[2]).min(0.0)
}

// ---------------------------------------------------------------- point areas

#[test]
fn voronoi_cell_hexagon_and_clipping() {
    // six neighbors on a hexagon of radius 1: the cell is the hexagon with apothem 1/2
    let q: Vec<[f64; 2]> = (0..6).map(|k| [(k as f64 * PI / 3.0).cos(), (k as f64 * PI / 3.0).sin()]).collect();
    let (a, clipped) = voronoi_cell_area(&q, 1.0);
    let hex = 6.0 * 0.5 * 0.5 * (PI / 6.0).tan();
    assert!((a - hex).abs() < 1e-14, "{a} vs {hex}");
    assert!(!clipped);
    // neighbors on one side only: unbounded, clipped by the square
    let (a, clipped) = voronoi_cell_area(&[[1.0, 0.0]], 1.0);
    assert!(clipped);
    assert!((a - 1.5 * 2.0).abs() < 1e-14, "{a}");
}

#[test]
fn point_areas_plane_grid_exact() {
    let m = 21;
    let h = 0.1;
    let mut p = Vec::new();
    for i in 0..m {
        for j in 0..m {
            p.push([i as f64 * h, j as f64 * h, 0.3]);
        }
    }
    let n = vec![[0.0, 0.0, 1.0]; p.len()];
    let cloud = PointCloud::new(p, n, NormalPolicy::default()).unwrap().0;
    let (a, rep) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
    for i in 3..m - 3 {
        for j in 3..m - 3 {
            let ai = a[i * m + j];
            assert!((ai - h * h).abs() < 1e-14, "interior point ({i}, {j}): {ai}");
        }
    }
    assert!(rep.clipped > 0 && rep.clipped < 4 * m, "only boundary cells are clipped: {rep:?}");
    assert_eq!(rep.isolated, 0);
}

#[test]
fn point_areas_sphere_total() {
    for (n, tol) in [(2000, 0.02), (8000, 0.01)] {
        let cloud = fib_sphere(n, 0.7);
        let (a, rep) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
        let exact = 4.0 * PI * 0.49;
        let rel = (rep.total_area - exact).abs() / exact;
        assert!(rel < tol, "n = {n}: Σ A = {} vs {exact} (rel {rel})", rep.total_area);
        assert_eq!(rep.clipped, 0, "{rep:?}");
        // uniform sampling: every area near the mean
        let mean = exact / n as f64;
        assert!(a.iter().all(|&x| (x / mean - 1.0).abs() < 0.25), "areas spread too much");
    }
}

#[test]
fn point_areas_duplicates_share_and_opposite_filtered() {
    // a plane grid, duplicated once: every cell is shared by two points
    let m = 11;
    let mut p = Vec::new();
    for i in 0..m {
        for j in 0..m {
            p.push([i as f64, j as f64, 0.0]);
        }
    }
    let mut pp = p.clone();
    pp.extend(p.iter().copied());
    let n = vec![[0.0, 0.0, 1.0]; pp.len()];
    let cloud = PointCloud::new(pp, n, NormalPolicy::default()).unwrap().0;
    let (a, rep) = point_areas(&cloud, &PointAreaOptions { k: 30, ..Default::default() }).unwrap();
    assert_eq!(rep.coincident, 2 * m * m);
    let c = 5 * m + 5;
    assert!((a[c] - 0.5).abs() < 1e-14 && (a[c + m * m] - 0.5).abs() < 1e-14, "{} {}", a[c], a[c + m * m]);
    // two parallel sheets with opposite normals close together: the other sheet is filtered
    let mut p = Vec::new();
    let mut n = Vec::new();
    for i in 0..m {
        for j in 0..m {
            p.push([i as f64, j as f64, 0.0]);
            n.push([0.0, 0.0, 1.0]);
            p.push([i as f64 + 0.5, j as f64 + 0.5, -0.01]);
            n.push([0.0, 0.0, -1.0]);
        }
    }
    let cloud = PointCloud::new(p, n, NormalPolicy::default()).unwrap().0;
    let (a, rep) = point_areas(&cloud, &PointAreaOptions { k: 40, ..Default::default() }).unwrap();
    assert!(rep.filtered > 0);
    assert!((a[2 * c] - 1.0).abs() < 1e-12, "{}", a[2 * c]);
}

// ---------------------------------------------------------------- winding numbers

#[test]
fn gwn_sphere_inside_outside() {
    let cloud = fib_sphere(4000, 0.5);
    let (a, _) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
    let (w_in, _) = winding_number(cloud.positions(), cloud.normals(), &a, [0.1, -0.05, 0.02]);
    let (w_out, _) = winding_number(cloud.positions(), cloud.normals(), &a, [0.9, 0.3, -0.2]);
    assert!((w_in - 1.0).abs() < 0.01, "{w_in}");
    assert!(w_out.abs() < 0.01, "{w_out}");
    // a query at a data point leaves that term out
    let (_, skipped) = winding_number(cloud.positions(), cloud.normals(), &a, cloud.positions()[7]);
    assert_eq!(skipped, 1);
}

fn direct(cloud: &PointCloud, a: &[f64], xs: &[[f64; 3]]) -> Vec<f64> {
    let mut w = vec![0.0; xs.len()];
    winding_numbers(cloud, a, xs, &mut w, None).unwrap();
    w
}

#[test]
fn fwn_points_matches_direct() {
    let cloud = ellipsoid_cloud(3000, 11, [0.8, 0.5, 0.3]);
    let (a, _) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
    let xs = rand_queries(2000, 5, 1.0);
    let exact = direct(&cloud, &a, &xs);
    // β → ∞: every leaf exactly (only the summation order differs)
    let all = FastWinding::from_cloud(&cloud, &a, &FwnOptions { beta: 1e12, ..Default::default() }).unwrap();
    let mut w = vec![0.0; xs.len()];
    all.evaluate(&xs, &mut w).unwrap();
    let d = w.iter().zip(&exact).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
    assert!(d < 1e-11, "β = ∞: {d}");
    // the error decreases with the expansion order; β = 2, order 2 is accurate to ~1e-3
    let mut errs = Vec::new();
    for order in 0..=2u8 {
        let f = FastWinding::from_cloud(&cloud, &a, &FwnOptions { order, ..Default::default() }).unwrap();
        f.evaluate(&xs, &mut w).unwrap();
        let rms = (w.iter().zip(&exact).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / xs.len() as f64).sqrt();
        errs.push(rms);
    }
    eprintln!("FWN point RMS error by order (β = 2): {errs:?}");
    assert!(errs[0] > errs[1] && errs[1] > errs[2], "{errs:?}");
    assert!(errs[2] < 2e-3, "{errs:?}");
    // a larger β is more accurate
    let f = FastWinding::from_cloud(&cloud, &a, &FwnOptions { beta: 4.0, ..Default::default() }).unwrap();
    f.evaluate(&xs, &mut w).unwrap();
    let rms4 = (w.iter().zip(&exact).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / xs.len() as f64).sqrt();
    assert!(rms4 < errs[2], "{rms4} vs {}", errs[2]);
}

#[test]
fn fwn_single_far_cluster_taylor_orders() {
    // a small cluster far away: the error of order p scales like (r/R)^{p+1}
    let mut s = 3;
    let mut p = Vec::new();
    let mut n = Vec::new();
    for _ in 0..40 {
        p.push([0.1 * lcg(&mut s), 0.1 * lcg(&mut s), 0.1 * lcg(&mut s)]);
        let v = [lcg(&mut s) - 0.5, lcg(&mut s) - 0.5, lcg(&mut s) - 0.5];
        let l = norm(v);
        n.push([v[0] / l, v[1] / l, v[2] / l]);
    }
    let cloud = PointCloud::new(p, n, NormalPolicy::default()).unwrap().0;
    let a: Vec<f64> = (0..40).map(|i| 0.5 + 0.01 * i as f64).collect();
    // truncation after order p leaves the order-(p+1) term, which decays like R^{−(p+3)}
    let dir = [1.3 / 1.8, -0.7 / 1.8, 0.9 / 1.8];
    for order in 0..=2u8 {
        // leaf_size 40: the root is a leaf, far by construction
        let f = FastWinding::from_cloud(&cloud, &a, &FwnOptions { order, leaf_size: Some(40), ..Default::default() }).unwrap();
        let err = |t: f64| {
            let x = [t * dir[0], t * dir[1], t * dir[2]];
            (f.winding_number(x).unwrap() - winding_number(cloud.positions(), cloud.normals(), &a, x).0).abs()
        };
        let ratio = err(8.0) / err(4.0);
        let expect = 0.5f64.powi(order as i32 + 3);
        eprintln!("order {order}: error ratio {ratio} (expected {expect}), error at R = 4: {}", err(4.0));
        assert!(ratio > 0.6 * expect && ratio < 1.4 * expect, "order {order}: ratio {ratio} vs {expect}");
    }
}

#[test]
fn fwn_triangles_matches_solid_angles() {
    let (v, f) = torus_mesh(0.5, 0.2, 48, 24);
    let tris: Vec<[[f64; 3]; 3]> = f.iter().map(|t| [v[t[0] as usize], v[t[1] as usize], v[t[2] as usize]]).collect();
    let fwn = FastWinding::from_triangles(&v, &f, &FwnOptions::default()).unwrap();
    let xs = rand_queries(1500, 9, 1.0);
    let mut w = vec![0.0; xs.len()];
    fwn.evaluate(&xs, &mut w).unwrap();
    let mut max: f64 = 0.0;
    for (x, wf) in xs.iter().zip(&w) {
        let we: f64 = tris.iter().map(|t| triangle_winding(t, *x)).sum();
        max = max.max((we - wf).abs());
        assert!((we.round() - we).abs() < 1e-9, "closed mesh: integer winding number, got {we}");
    }
    assert!(max < 1e-2, "FWN on triangles: max error {max}");
}

#[test]
fn regularized_winding_number_s_and_limit() {
    // S(t) = erf t − (2t/√π) e^{−t²}; erf(1), erf(2), erf(3) from tables
    let erf = [(1.0, 0.842_700_792_949_714_9), (2.0, 0.995_322_265_018_952_7), (3.0, 0.999_977_909_503_001_4)];
    for (t, e) in erf {
        let s_ref = e - 2.0 * t / PI.sqrt() * (-(t * t) as f64).exp();
        assert!((regularization_s(t) - s_ref).abs() < 1e-14, "S({t}) = {} vs {s_ref}", regularization_s(t));
    }
    // continuity across the series / continued-fraction switch at t² = 5/2
    let t0 = 2.5f64.sqrt();
    let (a, b) = (regularization_s(t0 * (1.0 - 1e-12)), regularization_s(t0 * (1.0 + 1e-12)));
    assert!((a - b).abs() < 1e-11, "{a} {b}");
    // small t: S(t) ≈ 4t³/(3√π)
    let t = 1e-4;
    assert!((regularization_s(t) / (4.0 * t * t * t / (3.0 * PI.sqrt())) - 1.0).abs() < 1e-7);
    // ε → 0 gives the winding number away from the points; finite at a data point
    let cloud = fib_sphere(2000, 0.5);
    let (a, _) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
    for x in [[0.05, 0.1, -0.1], [0.9, 0.1, 0.4]] {
        let (w, _) = winding_number(cloud.positions(), cloud.normals(), &a, x);
        let r = regularized_winding_number(cloud.positions(), cloud.normals(), &a, x, 1e-3);
        assert!((w - r).abs() < 1e-10, "{w} vs {r}");
    }
    let r = regularized_winding_number(cloud.positions(), cloud.normals(), &a, cloud.positions()[3], 0.02);
    assert!(r.is_finite() && (r - 0.5).abs() < 0.2, "on the surface: {r}");
}

// ---------------------------------------------------------------- SHC, SSPD, LSE

fn fd_check(field: &KernelField, method: KernelMethod, xs: &[[f64; 3]], tol: f64) {
    for &x in xs {
        let v = field.eval(method, x, true).unwrap();
        let h = 1e-6;
        for k in 0..3 {
            let mut a = x;
            let mut b = x;
            a[k] += h;
            b[k] -= h;
            let fd = (field.eval(method, a, false).unwrap().value - field.eval(method, b, false).unwrap().value) / (2.0 * h);
            let g = v.gradient[k];
            assert!((fd - g).abs() <= tol * (1.0 + g.abs()), "{method:?} at {x:?}, ∂{k}: fd {fd} vs {g}");
        }
    }
}

#[test]
fn kernel_gradients_by_finite_differences() {
    let cloud = ellipsoid_cloud(400, 21, [0.7, 0.5, 0.4]);
    let (a, _) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
    let xs = rand_queries(40, 77, 0.9);
    // constant λ: the frozen gradient is the gradient (the min-ρ shift cancels exactly)
    for nb in [Neighborhood::All, Neighborhood::Knn { k: 32 }] {
        let f = KernelField::new(&cloud, a.clone(), KernelOptions { neighborhood: nb, lambda: LambdaRule::Constant { lambda: 40.0 }, threads: Some(2) })
            .unwrap();
        for m in [KernelMethod::Shc, KernelMethod::Sspd, KernelMethod::Lse] {
            fd_check(&f, m, &xs, 2e-6);
        }
    }
    let f = KernelField::new(
        &cloud,
        a,
        KernelOptions { neighborhood: Neighborhood::All, lambda: LambdaRule::PaperHeuristic { c: 64.0, factor: 3.0, k_mean: 64 }, threads: None },
    )
    .unwrap();
    for m in [KernelMethod::Shc, KernelMethod::Sspd, KernelMethod::Lse] {
        fd_check(&f, m, &xs, 2e-6);
    }
}

#[test]
fn kernel_sphere_values_and_shift_invariance() {
    let cloud = fib_sphere(3000, 0.5);
    let (a, _) = point_areas(&cloud, &PointAreaOptions::default()).unwrap();
    let f = KernelField::new(&cloud, a.clone(), KernelOptions::default()).unwrap();
    // SSPD near the surface ≈ the signed distance
    for x in [[0.0, 0.0, 0.45], [0.0, 0.55, 0.0], [0.3, 0.3, 0.3]] {
        let exact = norm(x) - 0.5;
        let s = f.eval(KernelMethod::Sspd, x, false).unwrap().value;
        assert!((s - exact).abs() < 0.01, "SSPD {s} vs {exact} at {x:?}");
        let h = f.eval(KernelMethod::Shc, x, false).unwrap().value;
        assert!(h.signum() == exact.signum(), "SHC sign at {x:?}: {h}");
        let l = f.eval(KernelMethod::Lse, x, false).unwrap().value;
        assert!(l > 0.0 && l < exact.abs() + 0.05, "LSE {l} at {x:?}");
    }
    // the constant rule's shift (min ρ) cancels: a constant λ gives the same φ as the formula
    // without shift, computed by hand here for SSPD
    let lam = 30.0;
    let g = KernelField::new(&cloud, a.clone(), KernelOptions { neighborhood: Neighborhood::All, lambda: LambdaRule::Constant { lambda: lam }, threads: None })
        .unwrap();
    let x = [0.1, 0.2, 0.6];
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..cloud.len() {
        let p = cloud.positions()[i];
        let d = [x[0] - p[0], x[1] - p[1], x[2] - p[2]];
        let e = a[i] * (-lam * norm(d)).exp();
        let n = cloud.normals()[i];
        num += e * (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]);
        den += e;
    }
    let v = g.eval(KernelMethod::Sspd, x, false).unwrap().value;
    assert!((v - num / den).abs() < 1e-13, "{v} vs {}", num / den);
}

#[test]
fn kernel_errors() {
    // planar cloud, query in its plane: every ⟨d, n⟩ = 0, so w = 0
    let mut p = Vec::new();
    for i in 0..8 {
        for j in 0..8 {
            p.push([i as f64, j as f64, 0.0]);
        }
    }
    let n = vec![[0.0, 0.0, 1.0]; p.len()];
    let cloud = PointCloud::new(p, n, NormalPolicy::default()).unwrap().0;
    let a = vec![1.0; cloud.len()];
    let f = KernelField::new(&cloud, a.clone(), KernelOptions::default()).unwrap();
    assert!(matches!(f.eval(KernelMethod::Shc, [3.5, 3.5, 0.0], false), Err(PatError::DegenerateQuery { .. })));
    let mut out = vec![0.0; 2];
    match f.evaluate(KernelMethod::Shc, &[[3.5, 3.5, 1.0], [3.5, 2.5, 0.0]], &mut out, None) {
        Err(PatError::DegenerateQuery { index, .. }) => assert_eq!(index, 1),
        e => panic!("{e:?}"),
    }
    assert!(matches!(f.eval(KernelMethod::Sspd, [f64::NAN, 0.0, 0.0], false), Err(PatError::NonFinite { .. })));
    assert!(matches!(f.evaluate(KernelMethod::Sspd, &[[0.0; 3]], &mut [0.0; 2], None), Err(PatError::SizeMismatch { .. })));
    assert!(matches!(KernelField::new(&cloud, vec![1.0; 3], KernelOptions::default()), Err(PatError::SizeMismatch { .. })));
    assert!(matches!(
        KernelField::new(&cloud, a.clone(), KernelOptions { lambda: LambdaRule::Constant { lambda: -1.0 }, ..Default::default() }),
        Err(PatError::InvalidOption { .. })
    ));
    let mut neg = a.clone();
    neg[5] = -1.0;
    assert!(matches!(KernelField::new(&cloud, neg, KernelOptions::default()), Err(PatError::InvalidOption { .. })));
    assert!(matches!(FastWinding::from_cloud(&cloud, &a, &FwnOptions { order: 3, ..Default::default() }), Err(PatError::InvalidOption { .. })));
    assert!(matches!(FastWinding::from_triangles(&[[0.0; 3]; 3], &[[0, 1, 5]], &FwnOptions::default()), Err(PatError::Format { .. })));
    assert!(matches!(point_areas(&cloud, &PointAreaOptions { k: 0, ..Default::default() }), Err(PatError::InvalidOption { .. })));
}

// ---------------------------------------------------------------- mesh truth

#[test]
fn truth_cube_is_exact() {
    let (v, f) = cube(0.4);
    let (t, rep) = MeshTruth::new(&v, &f, &TruthOptions::default()).unwrap();
    assert!(rep.is_closed() && (rep.signed_volume - 0.512).abs() < 1e-14, "{rep:?}");
    let xs = rand_queries(5000, 1, 1.0);
    let mut sdf = vec![0.0; xs.len()];
    t.evaluate(&xs, &mut sdf, None).unwrap();
    for (x, s) in xs.iter().zip(&sdf) {
        let e = box_sdf(*x, 0.4);
        assert!((s - e).abs() < 1e-14, "{x:?}: {s} vs {e}");
    }
}

#[test]
fn truth_sphere_and_torus_vs_analytic() {
    // sphere: the inscribed polyhedron lies between the balls of radius r_in and r, so the SDFs
    // differ by at most r − r_in
    let r = 0.6;
    let (v, f) = icosphere(r, 4);
    let (t, rep) = MeshTruth::new(&v, &f, &TruthOptions::default()).unwrap();
    assert!(rep.is_closed());
    let r_in = f
        .iter()
        .map(|tr| norm(closest_point_on_triangle([0.0; 3], v[tr[0] as usize], v[tr[1] as usize], v[tr[2] as usize])))
        .fold(f64::INFINITY, f64::min);
    let xs = rand_queries(4000, 2, 1.0);
    let mut sdf = vec![0.0; xs.len()];
    let tr = t.evaluate(&xs, &mut sdf, None).unwrap();
    assert!(tr.inside > 0 && tr.ambiguous == 0, "{tr:?}");
    for (x, s) in xs.iter().zip(&sdf) {
        let e = norm(*x) - r;
        assert!((s - e).abs() <= (r - r_in) + 1e-12, "{x:?}: {s} vs {e} (bound {})", r - r_in);
    }
    // torus: the error is bounded by the chord sagittas of both circles
    let (big, small) = (0.5, 0.2);
    let (nu, nv) = (128, 64);
    let (v, f) = torus_mesh(big, small, nu, nv);
    let (t, rep) = MeshTruth::new(&v, &f, &TruthOptions::default()).unwrap();
    assert!(rep.is_closed() && rep.signed_volume > 0.0, "{rep:?}");
    let bound = (big + small) * (1.0 - (PI / nu as f64).cos()) + small * (1.0 - (PI / nv as f64).cos());
    let mut sdf = vec![0.0; xs.len()];
    t.evaluate(&xs, &mut sdf, None).unwrap();
    let mut max: f64 = 0.0;
    for (x, s) in xs.iter().zip(&sdf) {
        let q = ((x[0] * x[0] + x[1] * x[1]).sqrt() - big).hypot(x[2]) - small;
        max = max.max((s - q).abs());
    }
    eprintln!("torus mesh SDF vs analytic: max {max}, bound {bound}");
    assert!(max <= 2.0 * bound, "{max} vs {bound}");
}

#[test]
fn truth_closest_point_matches_brute_force_and_fast_sign() {
    let (v, f) = torus_mesh(0.5, 0.2, 40, 20);
    let (t, _) = MeshTruth::new(&v, &f, &TruthOptions { sign: SignMode::Exact, ..Default::default() }).unwrap();
    let (tf, _) = MeshTruth::new(&v, &f, &TruthOptions { sign: SignMode::Fast { beta: 2.0 }, ..Default::default() }).unwrap();
    let xs = rand_queries(3000, 4, 1.0);
    for x in &xs {
        let (d, _, face) = t.closest(*x);
        let mut best = (f64::INFINITY, u32::MAX);
        for (k, tr) in f.iter().enumerate() {
            let q = closest_point_on_triangle(*x, v[tr[0] as usize], v[tr[1] as usize], v[tr[2] as usize]);
            let dd = norm([x[0] - q[0], x[1] - q[1], x[2] - q[2]]);
            if dd < best.0 {
                best = (dd, k as u32);
            }
        }
        assert_eq!(d, best.0, "{x:?}");
        // ties within an ulp of d² may pick a neighbor of the brute-force face; its distance is the same
        let tr = f[face as usize];
        let q = closest_point_on_triangle(*x, v[tr[0] as usize], v[tr[1] as usize], v[tr[2] as usize]);
        assert_eq!(norm([x[0] - q[0], x[1] - q[1], x[2] - q[2]]), d, "face {face} vs {}", best.1);
        assert_eq!(t.signed_distance(*x).unwrap(), tf.signed_distance(*x).unwrap(), "fast sign differs at {x:?}");
    }
}

#[test]
fn truth_errors_and_soups() {
    let (v, mut f) = cube(0.5);
    // inward orientation
    let inv: Vec<[u32; 3]> = f.iter().map(|t| [t[0], t[2], t[1]]).collect();
    match MeshTruth::new(&v, &inv, &TruthOptions::default()) {
        Err(PatError::Format { reason }) => assert!(reason.contains("volume"), "{reason}"),
        e => panic!("{e:?}"),
    }
    // open mesh: rejected with the edge named, accepted as a soup
    f.pop();
    match MeshTruth::new(&v, &f, &TruthOptions::default()) {
        Err(PatError::Format { reason }) => assert!(reason.contains("edge"), "{reason}"),
        e => panic!("{e:?}"),
    }
    let (t, rep) = MeshTruth::new(&v, &f, &TruthOptions { require_closed: false, ..Default::default() }).unwrap();
    assert_eq!(rep.boundary_edges, 3);
    assert!(t.signed_distance([0.0, 0.0, 0.0]).unwrap() < 0.0);
    assert!(matches!(check_mesh(&v, &[[0, 1, 9]]), Err(PatError::Format { .. })));
    assert!(matches!(check_mesh(&v, &[[0, 1, 1]]), Err(PatError::Format { .. })));
    assert!(matches!(check_mesh(&[[f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], &[[0, 1, 2]]), Err(PatError::NonFinite { .. })));
    // a degenerate triangle does not break the closest point
    let q = closest_point_on_triangle([0.3, 1.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
    assert!((q[0] - 0.3).abs() < 1e-15 && q[1] == 0.0, "{q:?}");
    let (t, _) = MeshTruth::new(&v, &f, &TruthOptions { require_closed: false, ..Default::default() }).unwrap();
    assert!(matches!(t.signed_distance([f64::INFINITY, 0.0, 0.0]), Err(PatError::NonFinite { .. })));
}

// ---------------------------------------------------------------- determinism

#[test]
fn determinism_across_thread_counts() {
    let cloud = ellipsoid_cloud(1500, 31, [0.8, 0.6, 0.5]);
    let xs = rand_queries(3000, 8, 1.0);
    let (v, f) = torus_mesh(0.5, 0.2, 48, 24);
    let mut first: Option<Vec<u64>> = None;
    for threads in [1, 3, 8, 32] {
        let t = Some(threads);
        let (a, _) = point_areas(&cloud, &PointAreaOptions { threads: t, ..Default::default() }).unwrap();
        let mut bits: Vec<u64> = a.iter().map(|x| x.to_bits()).collect();
        let k = KernelField::new(&cloud, a.clone(), KernelOptions { threads: t, ..Default::default() }).unwrap();
        for m in [KernelMethod::Shc, KernelMethod::Sspd, KernelMethod::Lse] {
            let mut val = vec![0.0; xs.len()];
            let mut g = vec![[0.0; 3]; xs.len()];
            k.evaluate(m, &xs, &mut val, Some(&mut g)).unwrap();
            bits.extend(val.iter().map(|x| x.to_bits()));
            bits.extend(g.iter().flatten().map(|x| x.to_bits()));
        }
        let mut w = vec![0.0; xs.len()];
        winding_numbers(&cloud, &a, &xs, &mut w, t).unwrap();
        bits.extend(w.iter().map(|x| x.to_bits()));
        let fw = FastWinding::from_cloud(&cloud, &a, &FwnOptions { threads: t, ..Default::default() }).unwrap();
        fw.evaluate(&xs, &mut w).unwrap();
        bits.extend(w.iter().map(|x| x.to_bits()));
        let (tr, _) = MeshTruth::new(&v, &f, &TruthOptions { threads: t, sign: SignMode::Fast { beta: 2.0 }, ..Default::default() }).unwrap();
        let mut s = vec![0.0; xs.len()];
        tr.evaluate(&xs, &mut s, Some(&mut w)).unwrap();
        bits.extend(s.iter().chain(&w).map(|x| x.to_bits()));
        match &first {
            None => first = Some(bits),
            Some(b) => assert!(*b == bits, "results differ at {threads} threads"),
        }
    }
}

