//! Near-pair quadrature (WP16, `SPEC.md` §3.4, §17o): the inner integral, the graded outer
//! quadrature against an adaptive reference, the divergence rate, central differences of the
//! blended pair term and of whole-curve energies, the bitwise-unchanged guarantee when no pair
//! is near, thread independence, and the Barnes–Hut glue.

#![allow(clippy::needless_range_loop)]

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};
use std::sync::Arc;

use ddg_repulsive_curves::energy::barnes_hut::BhTree;
use ddg_repulsive_curves::energy::exact::make_pool;
use ddg_repulsive_curves::energy::near_pair::{
    InnerTables, NearPairOptions, blend_weight, closest_points, near_pairs, ordered_near_energy,
    pair_correction,
};
use ddg_repulsive_curves::step::{BarnesHutEnergy, EnergyModel, ExactEnergy, FlowOptions, Session};
use ddg_repulsive_curves::{BhOptions, CurveNetwork, EnergyEval, Exponents, energy};

type P = [f64; 3];

const GL8: [(f64, f64); 8] = [
    (-0.9602898564975362, 0.10122853629037669),
    (-0.7966664774136267, 0.22238103445337434),
    (-0.525532409916329, 0.31370664587788705),
    (-0.18343464249564978, 0.36268378337836177),
    (0.18343464249564978, 0.36268378337836177),
    (0.525532409916329, 0.31370664587788705),
    (0.7966664774136267, 0.22238103445337434),
    (0.9602898564975362, 0.10122853629037669),
];

fn gl(f: &dyn Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let (c, h) = (0.5 * (a + b), 0.5 * (b - a));
    GL8.iter().map(|&(x, w)| w * f(c + h * x)).sum::<f64>() * h
}

/// Adaptive Gauss–Legendre: split until two levels agree to `tol` (absolute).
fn adaptive(f: &dyn Fn(f64) -> f64, a: f64, b: f64, tol: f64, depth: u32) -> f64 {
    let whole = gl(f, a, b);
    let m = 0.5 * (a + b);
    let halves = gl(f, a, m) + gl(f, m, b);
    if depth == 0 || (whole - halves).abs() <= tol.max(1e-15 * halves.abs()) {
        halves
    } else {
        adaptive(f, a, m, 0.5 * tol, depth - 1) + adaptive(f, m, b, 0.5 * tol, depth - 1)
    }
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: P) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn lerp(a: P, b: P, t: f64) -> P {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}
fn kernel(x: P, y: P, t: P, al: f64, be: f64) -> f64 {
    let d = sub(x, y);
    let td = t[0] * d[0] + t[1] * d[1] + t[2] * d[2];
    let n = [d[0] - td * t[0], d[1] - td * t[1], d[2] - td * t[2]];
    norm(n).powf(al) / norm(d).powf(be)
}

fn exps(al: f64, be: f64) -> Exponents {
    Exponents::new(al, be).unwrap()
}

/// The reference pair integral `∫_I ∫_J k(x, y, T_I)`: the inner integral by
/// [`InnerTables::inner_value`] (checked separately), the outer one adaptively with breakpoints
/// at the closest point and the projections of I's ends, refined geometrically toward them.
fn reference(p: [P; 2], q: [P; 2], e: Exponents) -> f64 {
    let tab = InnerTables::new(e);
    let cl = closest_points(p[0], p[1], q[0], q[1]);
    let lj = norm(sub(q[1], q[0]));
    let f = |t: f64| tab.inner_value(lerp(q[0], q[1], t), p[0], p[1]);
    let proj = |x: P| {
        let b = sub(q[1], q[0]);
        let w = sub(x, q[0]);
        ((w[0] * b[0] + w[1] * b[1] + w[2] * b[2]) / (lj * lj)).clamp(0.0, 1.0)
    };
    let mut cuts = vec![0.0, 1.0];
    let scale = cl.d / lj;
    for c in [cl.t, proj(p[0]), proj(p[1])] {
        let mut s = 1e-3 * scale;
        while s < 1.0 {
            cuts.push(c - s);
            cuts.push(c + s);
            s *= 1.5;
        }
        cuts.push(c);
    }
    cuts.retain(|&c| (0.0..=1.0).contains(&c));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let mut sum = 0.0;
    for w in cuts.windows(2) {
        if w[1] > w[0] {
            // pieces of ratio 1.5 around every feature: 4 × 8-point Gauss–Legendre each
            let h = 0.25 * (w[1] - w[0]);
            for k in 0..4 {
                sum += gl(&f, w[0] + k as f64 * h, w[0] + (k + 1) as f64 * h);
            }
        }
    }
    sum * lj
}

fn near(p: [P; 2], q: [P; 2], e: Exponents, m: usize) -> f64 {
    let tab = InnerTables::new(e);
    let cl = closest_points(p[0], p[1], q[0], q[1]);
    ordered_near_energy(&tab, p, q, (cl.t, &cl.grad_t), (cl.d, &cl.grad_d), m, false).0
}

/// Crafted pairs at segment distance about `d` (I along x, J above).
fn configs(d: f64) -> Vec<(&'static str, [P; 2], [P; 2])> {
    let rot = |ang: f64| [ang.cos(), ang.sin(), 0.0];
    let i = [[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]];
    let sk = rot(0.3);
    let pa = rot(0.02);
    let tilt = [0.6f64.cos(), 0.0, 0.6f64.sin()];
    vec![
        ("perpendicular", i, [[0.0, -0.5, d], [0.0, 0.5, d]]),
        ("perpendicular-off", i, [[0.3, -0.7, d], [0.3, 0.3, d]]),
        (
            "skew",
            i,
            [
                [-0.45 * sk[0], -0.45 * sk[1], d],
                [0.55 * sk[0], 0.55 * sk[1], d],
            ],
        ),
        (
            "nearly-parallel",
            i,
            [[0.1, -0.005, d], [0.1 + pa[0], pa[1] - 0.005, d]],
        ),
        ("parallel", i, [[0.2, 0.0, d], [1.2, 0.0, d]]),
        ("t-mid", i, [[0.1, 0.0, d], [0.1, 0.0, d + 1.0]]),
        (
            "endpoint-closest",
            i,
            [[0.5 + d, 0.05, 0.0], [0.5 + d + tilt[0], 0.05, tilt[2]]],
        ),
    ]
}

#[test]
fn series_and_inner_integral_match_brute_force() {
    for (al, be) in [(2.0, 4.5), (3.0, 6.0), (2.5, 5.2)] {
        let e = exps(al, be);
        let tab = InnerTables::new(e);
        let m = be - 2.0;
        // S(θ) = ∫_0^θ cos^m.
        for th in [-1.55f64, -0.9, -0.3, 0.0, 0.2, FRAC_PI_4 - 1e-4, FRAC_PI_4 + 1e-4, 1.1, FRAC_PI_2 - 1e-4] {
            let f = |x: f64| x.cos().powf(m);
            let ex = if th >= 0.0 {
                adaptive(&f, 0.0, th, 1e-17, 24)
            } else {
                -adaptive(&f, th, 0.0, 1e-17, 24)
            };
            let s = tab.cos_power_integral(th.sin(), th.cos());
            assert!((s - ex).abs() < 2e-15, "S({th}) = {s}, quadrature {ex}");
        }
        // F(y) = ∫_I k(x, y, T_I) dx.
        let a = [0.1, -0.2, 0.05];
        let a1 = [1.0, 0.3, -0.1];
        let t = {
            let v = sub(a1, a);
            let l = norm(v);
            [v[0] / l, v[1] / l, v[2] / l]
        };
        for y in [
            [0.5, 0.4, 0.3],
            [0.55, 0.05, 0.0],
            [1.4, 0.5, -0.12],
            [-0.6, -0.5, 0.1],
            [0.4, 0.0, 0.03],
        ] {
            let f = |u: f64| kernel(lerp(a, a1, u), y, t, al, be);
            let ex = adaptive(&f, 0.0, 1.0, 1e-17, 24) * norm(sub(a1, a));
            let v = tab.inner_value(y, a, a1);
            assert!(
                ((v - ex) / ex).abs() < 1e-12,
                "({al}, {be}) y = {y:?}: {v} vs {ex}"
            );
            let g = tab.inner(y, a, a1);
            assert_eq!(g.value.to_bits(), v.to_bits());
        }
    }
}

#[test]
fn graded_quadrature_matches_the_adaptive_reference() {
    let mut worst = 0.0f64;
    for (al, be) in [(2.0, 4.5), (3.0, 6.0)] {
        let e = exps(al, be);
        for d in [0.5, 0.1, 1e-2, 1e-4, 1e-6] {
            for (name, p, q) in configs(d) {
                let r = reference(p, q, e);
                let v = near(p, q, e, 4);
                let rel = ((v - r) / r).abs();
                worst = worst.max(rel);
                assert!(
                    rel < 5e-6,
                    "({al},{be}) d = {d} {name}: {v} vs {r} ({rel:.1e})"
                );
                eprintln!("({al},{be}) d = {d:e} {name}: relative error {rel:.1e}");
            }
        }
    }
    eprintln!("graded quadrature, worst relative error {worst:.2e}");
    // The anchor of the prototype (SPEC §17o): perpendicular, d = 1e-6, (2, 4.5): 7531.4.
    let (_, p, q) = configs(1e-6)[0];
    let v = near(p, q, exps(2.0, 4.5), 4);
    assert!((v / 7531.4 - 1.0).abs() < 1e-4, "{v}");
}

#[test]
fn pair_integral_diverges_at_the_rate_of_the_polygon() {
    // Perpendicular pair crossing in the middle: e ~ C d^{2+α−β}.
    for (al, be) in [(3.0, 6.0), (2.0, 4.5)] {
        let e = exps(al, be);
        let rate = 2.0 + al - be;
        let mut prev: Option<f64> = None;
        for k in 5..10 {
            let d = 10f64.powi(-k);
            let (_, p, q) = configs(d)[0];
            let v = near(p, q, e, 4);
            if let Some(pv) = prev {
                let slope = (v / pv).log10() / -1.0;
                assert!(
                    (slope - rate).abs() < 0.01,
                    "({al}, {be}) d = {d}: slope {slope}, expected {rate}"
                );
            }
            prev = Some(v);
        }
    }
}

fn pair_fd(p: [P; 2], q: [P; 2], e: Exponents, o: &NearPairOptions) -> (f64, f64) {
    let tab = InnerTables::new(e);
    let f = |p: [P; 2], q: [P; 2]| pair_correction(&tab, e, o, p, q, true, true, true);
    let base = f(p, q);
    let scale = norm(sub(p[1], p[0]));
    let h = 1e-6 * scale;
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    let mut worst = 0.0f64;
    for k in 0..4 {
        for c in 0..3 {
            let mut pts = [p[0], p[1], q[0], q[1]];
            pts[k][c] += h;
            let up = f([pts[0], pts[1]], [pts[2], pts[3]]).value;
            pts[k][c] -= 2.0 * h;
            let dn = f([pts[0], pts[1]], [pts[2], pts[3]]).value;
            let fd = (up - dn) / (2.0 * h);
            let g = base.grad[k][c];
            num += (fd - g) * (fd - g);
            den += g * g;
            worst = worst.max((fd - g).abs());
        }
    }
    let _ = worst;
    ((num / den).sqrt(), base.weight)
}

#[test]
fn blended_pair_term_matches_central_differences() {
    let o = NearPairOptions::default();
    for (al, be) in [(2.0, 4.5), (3.0, 6.0)] {
        let e = exps(al, be);
        // Full weight (d ≪ c₁ ℓ), inside the blend (c₁ < r < c₂).
        for d in [3e-3, 0.08, 0.35, 0.5] {
            for (name, p, q) in configs(d) {
                let (err, w) = pair_fd(p, q, e, &o);
                // Exactly parallel overlapping edges: d_IJ is not differentiable there (the
                // closest points are not unique), so the blend has a kink unless w = 1.
                if w == 0.0 || (name == "parallel" && w < 1.0) {
                    continue;
                }
                assert!(
                    err < 2e-6,
                    "({al},{be}) d = {d} {name} (w = {w:.3}): relative FD error {err:.2e}"
                );
            }
        }
    }
}

#[test]
fn blend_weight_is_c2() {
    let (c1, c2) = (0.25, 0.75);
    assert_eq!(blend_weight(0.1, c1, c2), (1.0, 0.0));
    assert_eq!(blend_weight(0.25, c1, c2), (1.0, 0.0));
    assert_eq!(blend_weight(0.75, c1, c2), (0.0, 0.0));
    let h = 1e-6;
    for r in [0.3, 0.5, 0.7] {
        let fd = (blend_weight(r + h, c1, c2).0 - blend_weight(r - h, c1, c2).0) / (2.0 * h);
        assert!((fd - blend_weight(r, c1, c2).1).abs() < 1e-8);
    }
    assert!((blend_weight(0.5, c1, c2).0 - 0.5).abs() < 1e-15);
}

// ---------------------------------------------------------------------------------------------
// Whole curves
// ---------------------------------------------------------------------------------------------

/// A (2,3) torus knot of `n` vertices, with the strand near parameter `u₀` pushed toward the
/// strand near `u₁` so that two edges cross at distance about `gap` in their interiors.
fn trefoil(n: usize) -> (CurveNetwork, Vec<P>) {
    let pos: Vec<P> = (0..n)
        .map(|i| {
            let u = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            [
                (2.0 + (3.0 * u).cos()) * (2.0 * u).cos(),
                (2.0 + (3.0 * u).cos()) * (2.0 * u).sin(),
                (3.0 * u).sin(),
            ]
        })
        .collect();
    let edges = (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

/// Two straight strands crossing (one along x at z = 0, one along y at z = `gap`), a few
/// edges each, joined into one loop far away: a mid-edge near crossing.
fn crossing_loop(gap: f64, k: usize) -> (CurveNetwork, Vec<P>) {
    let mut pos = Vec::new();
    // Strand A: x from -k/2 to k/2 (offset by 0.5 so the crossing is mid-edge), z = 0.
    for i in 0..=k {
        pos.push([i as f64 - k as f64 / 2.0 + 0.5 - 0.37, 0.0, 0.0]);
    }
    // Arc up and over to strand B.
    let r = k as f64;
    pos.push([r, 0.0, 2.0]);
    pos.push([r * 0.7, r * 0.7, 4.0]);
    pos.push([0.31, r, gap + 3.0]);
    // Strand B: y from k/2 to -k/2 at z = gap, x = 0.31.
    for i in 0..=k {
        pos.push([0.31, k as f64 / 2.0 - i as f64 + 0.5 - 0.21, gap]);
    }
    pos.push([0.31, -r, -2.0]);
    pos.push([-r, -0.5 * r, -3.0]);
    let n = pos.len();
    let edges = (0..n as u32).map(|i| [i, (i + 1) % n as u32]).collect();
    (CurveNetwork::new(n, edges).unwrap(), pos)
}

fn pool(t: usize) -> Arc<rayon::ThreadPool> {
    Arc::new(make_pool(Some(t)).unwrap())
}

fn model_eval(model: &mut dyn EnergyModel, net: &CurveNetwork, pos: &[P]) -> (f64, Vec<P>) {
    model.prepare(net, pos).unwrap();
    let ev = model.energy_and_differential(net, pos).unwrap();
    (ev.energy, ev.differential)
}

fn bits(v: &[P]) -> Vec<u64> {
    v.iter().flatten().map(|x| x.to_bits()).collect()
}

#[test]
fn no_near_pair_means_bitwise_the_paper() {
    let (net, pos) = trefoil(120);
    assert!(near_pairs(&net, &pos, &NearPairOptions::default()).is_empty());
    for (al, be) in [(3.0, 6.0), (2.0, 4.5)] {
        let e = exps(al, be);
        let pl = pool(3);
        let mut a = ExactEnergy::new(pl.clone(), e);
        let mut b = ExactEnergy::new(pl.clone(), e).with_near_pair(Some(Default::default()));
        let (ea, ga) = model_eval(&mut a, &net, &pos);
        let (eb, gb) = model_eval(&mut b, &net, &pos);
        assert_eq!(ea.to_bits(), eb.to_bits());
        assert_eq!(bits(&ga), bits(&gb));
        assert_eq!(
            a.energy(&net, &pos).unwrap().to_bits(),
            b.energy(&net, &pos).unwrap().to_bits()
        );
        let bh = BhOptions::default();
        let mut a = BarnesHutEnergy::new(pl.clone(), e, bh);
        let mut b =
            BarnesHutEnergy::new(pl.clone(), e, bh).with_near_pair(Some(Default::default()));
        let (ea, ga) = model_eval(&mut a, &net, &pos);
        let (eb, gb) = model_eval(&mut b, &net, &pos);
        assert_eq!(ea.to_bits(), eb.to_bits());
        assert_eq!(bits(&ga), bits(&gb));
    }
    // A few flow steps with and without the option: the same curves, bit for bit.
    let run = |near: Option<NearPairOptions>| {
        let o = FlowOptions {
            exponents: exps(2.0, 4.5),
            near_pair: near,
            threads: Some(2),
            ..FlowOptions::default()
        };
        let mut s = Session::new(net.clone(), pos.clone(), o).unwrap();
        for _ in 0..3 {
            s.step().unwrap();
        }
        bits(s.positions())
    };
    assert_eq!(run(None), run(Some(NearPairOptions::default())));
}

#[test]
fn near_pairs_are_found_and_the_energy_sees_them() {
    let e = exps(2.0, 4.5);
    let (net, pos) = crossing_loop(1e-3, 6);
    let o = NearPairOptions::default();
    let pairs = near_pairs(&net, &pos, &o);
    assert!(!pairs.is_empty());
    let paper = energy(&net, &pos, e, &EnergyEval::Exact, Some(1))
        .unwrap()
        .energy;
    let mut m = ExactEnergy::new(pool(1), e).with_near_pair(Some(o));
    let (en, _) = model_eval(&mut m, &net, &pos);
    // The crossing pair alone carries about 7.5 d^{-1/2} · 2 (both orders) ≈ 470.
    assert!(en > paper + 300.0, "paper {paper}, with near pairs {en}");
}

fn curve_fd(model: &mut dyn EnergyModel, net: &CurveNetwork, pos: &[P], h: f64) -> f64 {
    let (_, g) = model_eval(model, net, pos);
    // Directional derivatives along a few fixed pseudo-random directions.
    let mut worst = 0.0f64;
    for seed in 0..3u64 {
        let v: Vec<P> = (0..pos.len())
            .map(|i| {
                let x = ((i as u64 * 2654435761 + seed * 40503) % 1000) as f64 / 500.0 - 1.0;
                let y = ((i as u64 * 97 + seed * 13 + 7) % 1000) as f64 / 500.0 - 1.0;
                let z = ((i as u64 * 31 + seed * 101 + 3) % 1000) as f64 / 500.0 - 1.0;
                [x, y, z]
            })
            .collect();
        let mut at = |s: f64| -> f64 {
            let p: Vec<P> = pos
                .iter()
                .zip(&v)
                .map(|(p, v)| [p[0] + s * v[0], p[1] + s * v[1], p[2] + s * v[2]])
                .collect();
            model.energy(net, &p).unwrap()
        };
        let fd = (at(h) - at(-h)) / (2.0 * h);
        let gv: f64 = g
            .iter()
            .zip(&v)
            .map(|(g, v)| g[0] * v[0] + g[1] * v[1] + g[2] * v[2])
            .sum();
        let gn: f64 = g.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
        let vn: f64 = v.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
        worst = worst.max((fd - gv).abs() / (gn * vn));
    }
    worst
}

#[test]
fn curve_differential_matches_central_differences() {
    for (al, be) in [(2.0, 4.5), (3.0, 6.0)] {
        let e = exps(al, be);
        for gap in [1e-4, 0.05, 0.4] {
            let (net, pos) = crossing_loop(gap, 6);
            let mut m = ExactEnergy::new(pool(2), e).with_near_pair(Some(Default::default()));
            let err = curve_fd(&mut m, &net, &pos, 1e-3 * gap.min(0.01));
            assert!(err < 1e-7, "({al},{be}) gap {gap}: {err:.2e}");
            let bh = BhOptions {
                theta_x: 0.5,
                theta_t: 0.5,
                leaf_size: 2,
            };
            let mut m =
                BarnesHutEnergy::new(pool(2), e, bh).with_near_pair(Some(Default::default()));
            let err = curve_fd(&mut m, &net, &pos, 1e-3 * gap.min(0.01));
            assert!(err < 1e-7, "BH ({al},{be}) gap {gap}: {err:.2e}");
        }
    }
}

#[test]
fn results_do_not_depend_on_the_thread_count() {
    let e = exps(2.0, 4.5);
    let (net, pos) = crossing_loop(2e-3, 8);
    let mut out = Vec::new();
    for t in [1, 3, 8] {
        let mut a = ExactEnergy::new(pool(t), e).with_near_pair(Some(Default::default()));
        let (ea, ga) = model_eval(&mut a, &net, &pos);
        let bh = BhOptions {
            theta_x: 0.5,
            theta_t: 0.5,
            leaf_size: 2,
        };
        let mut b = BarnesHutEnergy::new(pool(t), e, bh).with_near_pair(Some(Default::default()));
        let (eb, gb) = model_eval(&mut b, &net, &pos);
        out.push((ea.to_bits(), bits(&ga), eb.to_bits(), bits(&gb)));
    }
    assert!(out.windows(2).all(|w| w[0] == w[1]));
}

#[test]
fn barnes_hut_without_lumping_equals_exact() {
    let e = exps(3.0, 6.0);
    let (net, pos) = crossing_loop(1e-2, 6);
    let tiny = BhOptions {
        theta_x: 1e-12,
        theta_t: 1e-12,
        leaf_size: 4,
    };
    let tree = BhTree::build(&net, &pos, tiny).unwrap();
    let lists = tree.interaction_lists(&make_pool(Some(1)).unwrap());
    assert_eq!(lists.n_lumped(), 0);
    let mut a = ExactEnergy::new(pool(1), e).with_near_pair(Some(Default::default()));
    let mut b = BarnesHutEnergy::new(pool(1), e, tiny).with_near_pair(Some(Default::default()));
    let (ea, ga) = model_eval(&mut a, &net, &pos);
    let (eb, gb) = model_eval(&mut b, &net, &pos);
    assert!(((ea - eb) / ea).abs() < 1e-13, "{ea} vs {eb}");
    let gn: f64 = ga.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
    let dn: f64 = ga
        .iter()
        .flatten()
        .zip(gb.iter().flatten())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt();
    assert!(dn / gn < 1e-12, "{}", dn / gn);
}

#[test]
fn options_are_validated() {
    let bad = [
        NearPairOptions {
            blend_start: 0.0,
            ..Default::default()
        },
        NearPairOptions {
            blend_start: 0.8,
            blend_end: 0.75,
            ..Default::default()
        },
        NearPairOptions {
            blend_end: f64::NAN,
            ..Default::default()
        },
        NearPairOptions {
            subpanels: 0,
            ..Default::default()
        },
    ];
    for o in bad {
        assert!(o.validate().is_err(), "{o:?}");
        let f = FlowOptions {
            near_pair: Some(o),
            ..FlowOptions::default()
        };
        assert!(f.validate().is_err());
    }
    assert!(NearPairOptions::default().validate().is_ok());
}
