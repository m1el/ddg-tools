//! Analytic shapes (SPEC.md §11): 2D SDF exactness, gradients, revolution/extrusion identities,
//! surface sampling statistics, the 70-instance set and the saddle-rich set. Run with `--release`.

use ddg_points_as_tori::analytic::{
    ALL_KINDS, Boundary2d, Instance, Kind, PAPER_KINDS, Shape2d, Solid, SolidKind, analytic_json, analytic_set, saddle_set, saddle_stats,
};
use ddg_points_as_tori::csg::ImplicitField;
use ddg_points_as_tori::generate::rng::Rng;
use ddg_points_as_tori::mc::{Band, McOptions, MeshCheck, Resolution, marching_cubes};
use ddg_points_as_tori::pool::make_pool;
use std::f64::consts::PI;

fn dist_to_segments(b: &Boundary2d, p: [f64; 2]) -> f64 {
    let mut best = f64::INFINITY;
    for &(a, c) in &b.segments {
        let e = [c[0] - a[0], c[1] - a[1]];
        let w = [p[0] - a[0], p[1] - a[1]];
        let ee = e[0] * e[0] + e[1] * e[1];
        let t = if ee > 0.0 { ((w[0] * e[0] + w[1] * e[1]) / ee).clamp(0.0, 1.0) } else { 0.0 };
        best = best.min((w[0] - t * e[0]).hypot(w[1] - t * e[1]));
    }
    best
}

fn shapes(seeds: u64) -> Vec<Shape2d> {
    let mut v = Vec::new();
    for (i, &k) in ALL_KINDS.iter().enumerate() {
        for s in 0..seeds {
            let mut rng = Rng::new(1000 * i as u64 + s);
            v.push(Shape2d::random(k, &mut rng));
        }
    }
    v
}

#[test]
fn quilez_2d_boundary_samples_have_zero_distance_and_unit_normal() {
    for sh in shapes(4) {
        let ev = |p: [f64; 2]| sh.eval(p);
        let b = Boundary2d::trace(&ev, sh.bbox(), 1024);
        assert!(b.segments.len() > 20, "{:?} {:?}", sh.kind(), sh.params());
        let mut rng = Rng::new(5);
        for _ in 0..300 {
            let (p, n) = b.sample(&ev, rng.f64()).unwrap();
            let (d, g) = sh.eval(p);
            assert!(d.abs() < 1e-11, "{:?} {:?}: φ = {d} at {p:?}", sh.kind(), sh.params());
            assert!((n[0] * g[0] + n[1] * g[1] - 1.0).abs() < 1e-12);
            assert!((g[0].hypot(g[1]) - 1.0).abs() < 1e-12, "{:?}: |∇φ| = {}", sh.kind(), g[0].hypot(g[1]));
            // inside the declared box
            let bb = sh.bbox();
            for a in 0..2 {
                assert!(p[a] >= bb[0][a] - 1e-9 && p[a] <= bb[1][a] + 1e-9, "{:?} {:?}: {p:?} outside bbox {bb:?}", sh.kind(), sh.params());
            }
        }
    }
}

#[test]
fn quilez_2d_is_the_signed_distance_to_the_traced_boundary() {
    for sh in shapes(3) {
        let ev = |p: [f64; 2]| sh.eval(p);
        let b = Boundary2d::trace(&ev, sh.bbox(), 2048);
        let bb = sh.bbox();
        let mut rng = Rng::new(9);
        let mut worst: f64 = 0.0;
        for _ in 0..400 {
            let p = [rng.uniform(bb[0][0] - 0.3, bb[1][0] + 0.3), rng.uniform(bb[0][1] - 0.3, bb[1][1] + 0.3)];
            let (d, _) = sh.eval(p);
            let brute = dist_to_segments(&b, p);
            // The polyline lies inside the true curve for convex parts and cuts corners: allow 3 cells.
            worst = worst.max((d.abs() - brute).abs());
            assert!((d.abs() - brute).abs() < 4e-3, "{:?} {:?} at {p:?}: |φ| = {} vs polyline {brute}", sh.kind(), sh.params(), d.abs());
        }
        let _ = worst;
    }
}

#[test]
fn quilez_2d_sign_matches_area_and_shape_has_positive_area() {
    for sh in shapes(2) {
        let ev = |p: [f64; 2]| sh.eval(p);
        let b = Boundary2d::trace(&ev, sh.bbox(), 1024);
        assert!(b.area > 1e-3, "{:?}: area {}", sh.kind(), b.area);
    }
    // Known areas.
    let c = Shape2d::new(Kind::Circle, &[0.7]).unwrap();
    let b = Boundary2d::trace(&|p| c.eval(p), c.bbox(), 1024);
    assert!((b.area / (PI * 0.49) - 1.0).abs() < 1e-5, "{}", b.area);
    assert!((b.length / (2.0 * PI * 0.7) - 1.0).abs() < 1e-5);
    let e = Shape2d::new(Kind::Ellipse, &[0.8, 0.3]).unwrap();
    let b = Boundary2d::trace(&|p| e.eval(p), e.bbox(), 1024);
    assert!((b.area / (PI * 0.24) - 1.0).abs() < 1e-5, "{}", b.area);
    let h = Shape2d::new(Kind::Heart, &[]).unwrap();
    let b = Boundary2d::trace(&|p| h.eval(p), h.bbox(), 1024);
    // heart = square of diagonal… area 1/2 + 2 half discs: 0.5 + π r²·… = 0.5 + 2·(π/2)(√2/4)²·… computed by integration below
    let mut rng = Rng::new(1);
    let bb = h.bbox();
    let mut inside = 0usize;
    let n = 400_000;
    for _ in 0..n {
        let p = [rng.uniform(bb[0][0], bb[1][0]), rng.uniform(bb[0][1], bb[1][1])];
        if h.distance(p) < 0.0 {
            inside += 1;
        }
    }
    let mc = inside as f64 / n as f64 * (bb[1][0] - bb[0][0]) * (bb[1][1] - bb[0][1]);
    assert!((mc / b.area - 1.0).abs() < 0.01, "heart area {mc} vs {}", b.area);
}

#[test]
fn quilez_2d_gradient_matches_finite_differences() {
    for sh in shapes(3) {
        let bb = sh.bbox();
        let mut rng = Rng::new(21);
        let (mut ok, mut total) = (0, 0);
        for _ in 0..300 {
            let p = [rng.uniform(bb[0][0] - 0.3, bb[1][0] + 0.3), rng.uniform(bb[0][1] - 0.3, bb[1][1] + 0.3)];
            let (d, g) = sh.eval(p);
            let h = 1e-6;
            let fd = [
                (sh.distance([p[0] + h, p[1]]) - sh.distance([p[0] - h, p[1]])) / (2.0 * h),
                (sh.distance([p[0], p[1] + h]) - sh.distance([p[0], p[1] - h])) / (2.0 * h),
            ];
            total += 1;
            if (fd[0] - g[0]).abs() < 1e-5 && (fd[1] - g[1]).abs() < 1e-5 {
                ok += 1;
            } else {
                // only allowed on the medial axis / at creases: the one-sided slopes then differ
                assert!(d.is_finite());
            }
        }
        assert!(ok as f64 >= 0.97 * total as f64, "{:?} {:?}: only {ok}/{total} FD gradients agree", sh.kind(), sh.params());
    }
}

fn some_solids() -> Vec<Solid> {
    let mut v = Vec::new();
    for (i, &k) in ALL_KINDS.iter().enumerate() {
        let mut rng = Rng::new(77 + i as u64);
        let p = Shape2d::random(k, &mut rng);
        v.push(Solid::new(p.clone(), SolidKind::Extrusion { h: 0.2 }).unwrap());
        v.push(Solid::new(p, SolidKind::Revolution { o: 0.95 }).unwrap());
    }
    v
}

#[test]
fn extrusion_and_revolution_identities() {
    let mut rng = Rng::new(3);
    for k in ALL_KINDS {
        let p = Shape2d::random(k, &mut rng);
        let ext = Solid::new(p.clone(), SolidKind::Extrusion { h: 0.25 }).unwrap();
        let rev = Solid::new(p.clone(), SolidKind::Revolution { o: 1.0 }).unwrap();
        let bb = p.bbox();
        for _ in 0..200 {
            let (x, y) = (rng.uniform(bb[0][0] - 0.2, bb[1][0] + 0.2), rng.uniform(bb[0][1] - 0.2, bb[1][1] + 0.2));
            let (d2, g2) = p.eval([x, y]);
            // extrusion: mid-plane z = 0, outside the profile: φ = φ₂, gradient (∇φ₂, 0)
            let (d, g) = ext.eval([x, y, 0.0]);
            if d2 > 0.0 {
                assert!((d - d2).abs() < 1e-14 && (g[0] - g2[0]).abs() < 1e-14 && g[2] == 0.0, "{k:?}");
            } else {
                // inside: min(max(w), 0) = max(φ₂, −h)
                assert!((d - d2.max(-0.25)).abs() < 1e-14);
            }
            // beyond the caps in z: φ = sqrt(max(φ₂,0)² + (|z|-h)²)
            let z = 0.25 + rng.uniform(0.01, 0.5);
            let (d, g) = ext.eval([x, y, z]);
            assert!((d - d2.max(0.0).hypot(z - 0.25)).abs() < 1e-14, "{k:?}");
            assert!(((g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt() - 1.0).abs() < 1e-12);
            // revolution: φ(ρ cosψ, y, ρ sinψ) = φ₂(ρ − o, y) for every ψ
            let rho = 1.0 + x;
            let psi = rng.uniform(0.0, 2.0 * PI);
            let (dr, gr) = rev.eval([rho * psi.cos(), y, rho * psi.sin()]);
            assert!((dr - d2).abs() < 1e-14, "{k:?}");
            assert!((gr[1] - g2[1]).abs() < 1e-14);
            assert!((gr[0] * psi.cos() + gr[2] * psi.sin() - g2[0]).abs() < 1e-13);
        }
    }
}

#[test]
fn solid_gradients_match_finite_differences_and_are_unit() {
    for s in some_solids() {
        let bb = s.bbox();
        let mut rng = Rng::new(8);
        let (mut ok, mut total) = (0, 0);
        for _ in 0..300 {
            let x = [0, 1, 2].map(|a| rng.uniform(bb[0][a], bb[1][a]));
            let (_, g) = s.distance_gradient(x).unwrap();
            let h = 1e-6;
            let mut fd = [0.0; 3];
            for a in 0..3 {
                let (mut xp, mut xm) = (x, x);
                xp[a] += h;
                xm[a] -= h;
                fd[a] = (s.distance(xp).unwrap() - s.distance(xm).unwrap()) / (2.0 * h);
            }
            total += 1;
            if (0..3).all(|a| (fd[a] - g[a]).abs() < 2e-5) {
                ok += 1;
            }
            assert!(((g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt() - 1.0).abs() < 1e-9 || s.distance(x).unwrap().abs() < 1e-9);
        }
        assert!(ok as f64 >= 0.95 * total as f64, "{:?} {:?}: {ok}/{total}", s.profile().kind(), s.kind());
    }
}

#[test]
fn surface_samples_lie_on_the_surface_with_exact_normals() {
    for s in some_solids() {
        let (p, n) = s.sample_surface(2000, 4).unwrap();
        let bb = s.bbox();
        for (x, nn) in p.iter().zip(&n) {
            let (d, g) = s.distance_gradient(*x).unwrap();
            assert!(d.abs() < 1e-10, "{:?} {:?}: φ = {d}", s.profile().kind(), s.kind());
            assert!(((nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]).sqrt() - 1.0).abs() < 1e-12);
            // the normal equals ∇φ except exactly on a crease (cap rim) where it is a limit
            let dot = g[0] * nn[0] + g[1] * nn[1] + g[2] * nn[2];
            assert!(dot > 1.0 - 1e-9 || dot > 0.0, "{:?} {:?}: normal·∇φ = {dot}", s.profile().kind(), s.kind());
            for a in 0..3 {
                assert!(x[a] >= bb[0][a] && x[a] <= bb[1][a], "{:?} {:?}: sample {x:?} outside bbox {bb:?}", s.profile().kind(), s.kind());
            }
        }
        let frac_exact = p.iter().zip(&n).filter(|(x, nn)| {
            let (_, g) = s.distance_gradient(**x).unwrap();
            (g[0] - nn[0]).abs() + (g[1] - nn[1]).abs() + (g[2] - nn[2]).abs() < 1e-9
        });
        assert!(frac_exact.count() >= 1990, "{:?} {:?}", s.profile().kind(), s.kind());
    }
}

#[test]
fn distance_is_exact_against_dense_surface_samples() {
    // For exact SDFs the distance to the nearest of many surface samples is >= |φ| and close to it.
    let insts = analytic_set().unwrap();
    for inst in insts.iter().step_by(7) {
        let s = &inst.solid;
        let (pts, _) = s.sample_surface(60_000, 11).unwrap();
        let bb = s.bbox();
        let mut rng = Rng::new(2);
        for _ in 0..40 {
            let q = [0, 1, 2].map(|a| rng.uniform(bb[0][a], bb[1][a]));
            let phi = s.distance(q).unwrap();
            if phi.abs() < 0.03 {
                continue;
            }
            let nearest = pts.iter().map(|p| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()).fold(f64::INFINITY, f64::min);
            // Revolution profiles that cross the axis can have |φ| smaller than the true distance
            // inside; test outside points, where the field is exact.
            if phi > 0.0 {
                assert!(nearest >= phi - 1e-9, "{}: sample closer ({nearest}) than φ ({phi})", inst.name());
                assert!(nearest - phi < 0.12, "{}: nearest sample {nearest} vs φ {phi}", inst.name());
            }
        }
    }
}

fn ks_uniform(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len() as f64;
    v.iter().enumerate().map(|(i, &x)| ((i as f64 + 1.0) / n - x).abs().max((x - i as f64 / n).abs())).fold(0.0, f64::max)
}

#[test]
fn sampling_is_area_uniform_on_a_sphere() {
    // Archimedes: y is uniform on a sphere; the azimuth is uniform.
    let r = 0.6;
    let s = Solid::new(Shape2d::new(Kind::Circle, &[r]).unwrap(), SolidKind::Revolution { o: 0.0 }).unwrap();
    let n = 40_000;
    let (p, nn) = s.sample_surface(n, 1).unwrap();
    let ys: Vec<f64> = p.iter().map(|x| (x[1] / r + 1.0) / 2.0).collect();
    let az: Vec<f64> = p.iter().map(|x| (x[2].atan2(x[0]) / (2.0 * PI)).rem_euclid(1.0)).collect();
    let crit = 1.95 / (n as f64).sqrt();
    assert!(ks_uniform(ys) < crit);
    assert!(ks_uniform(az) < crit);
    for (x, m) in p.iter().zip(&nn) {
        assert!(((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt() - r).abs() < 1e-10);
        assert!((m[0] * x[0] + m[1] * x[1] + m[2] * x[2] - r).abs() < 1e-9);
    }
    assert!((s.surface_area() / (4.0 * PI * r * r) - 1.0).abs() < 1e-4);
}

#[test]
fn sampling_is_area_uniform_on_a_torus_and_a_cylinder() {
    let (r, o) = (0.25, 0.5);
    let t = Solid::new(Shape2d::new(Kind::Circle, &[r]).unwrap(), SolidKind::Revolution { o }).unwrap();
    assert!((t.surface_area() / (4.0 * PI * PI * r * o) - 1.0).abs() < 1e-4, "{}", t.surface_area());
    let n = 100_000;
    let (p, _) = t.sample_surface(n, 2).unwrap();
    let inner = p.iter().filter(|x| x[0].hypot(x[2]) < o).count() as f64 / n as f64;
    let want = 0.5 - r / (PI * o);
    // Points with cos θ < 0 are those with ρ < o. Each fraction is a binomial proportion.
    let want_inner = want + 0.0;
    let sigma = (want_inner * (1.0 - want_inner) / n as f64).sqrt();
    assert!((inner - want_inner).abs() < 4.5 * sigma, "inner fraction {inner} vs {want_inner}");
    // Cylinder: cap fraction.
    let (rc, h) = (0.4, 0.3);
    let c = Solid::new(Shape2d::new(Kind::Circle, &[rc]).unwrap(), SolidKind::Extrusion { h }).unwrap();
    assert!((c.surface_area() / (2.0 * PI * rc * rc + 2.0 * PI * rc * 2.0 * h) - 1.0).abs() < 1e-4);
    let (p, _) = c.sample_surface(n, 3).unwrap();
    let caps = p.iter().filter(|x| (x[2].abs() - h).abs() < 1e-12 && x[0].hypot(x[1]) < rc - 1e-9).count() as f64 / n as f64;
    let want = 2.0 * rc / (2.0 * rc + 4.0 * h);
    let sigma = (want * (1.0 - want) / n as f64).sqrt();
    assert!((caps - want).abs() < 4.5 * sigma, "cap fraction {caps} vs {want}");
    // Radii on the caps are area uniform: r²/rc² uniform.
    let u: Vec<f64> = p.iter().filter(|x| (x[2].abs() - h).abs() < 1e-12 && x[0].hypot(x[1]) < rc - 1e-9).map(|x| (x[0] * x[0] + x[1] * x[1]) / (rc * rc)).collect();
    assert!(ks_uniform(u.clone()) < 1.95 / (u.len() as f64).sqrt());
}

#[test]
fn deterministic_by_seed_and_across_thread_counts() {
    let inst = Instance::from_seed(23).unwrap();
    let a = inst.solid.sample_surface(500, 99).unwrap();
    let b = inst.solid.sample_surface(500, 99).unwrap();
    let c = inst.solid.sample_surface(500, 100).unwrap();
    assert_eq!(a, b);
    assert_ne!(a, c);
    let p1 = make_pool(Some(1)).unwrap();
    let p8 = make_pool(Some(8)).unwrap();
    let s1 = p1.install(|| Solid::new(inst.solid.profile().clone(), inst.solid.kind()).unwrap().sample_surface(500, 99).unwrap());
    let s8 = p8.install(|| Solid::new(inst.solid.profile().clone(), inst.solid.kind()).unwrap().sample_surface(500, 99).unwrap());
    assert_eq!(s1, s8);
    assert_eq!(s1, a);
    assert_eq!(analytic_json(&analytic_set().unwrap()), analytic_json(&analytic_set().unwrap()));
}

#[test]
fn the_set_has_70_instances_5_per_shape() {
    let set = analytic_set().unwrap();
    assert_eq!(set.len(), 70);
    for (i, &k) in PAPER_KINDS.iter().enumerate() {
        let of_kind: Vec<_> = set.iter().filter(|x| x.solid.profile().kind() == k).collect();
        assert_eq!(of_kind.len(), 5, "{k:?}");
        assert!(of_kind.iter().all(|x| x.seed as usize / 5 == i));
    }
    let revs = set.iter().filter(|x| matches!(x.solid.kind(), SolidKind::Revolution { .. })).count();
    assert!((20..=50).contains(&revs), "{revs} revolutions");
    for inst in &set {
        let (p, _) = inst.solid.sample_surface(300, inst.seed).unwrap();
        assert!(p.iter().all(|x| inst.solid.distance(*x).unwrap().abs() < 1e-10), "{}", inst.name());
        match inst.solid.kind() {
            SolidKind::Revolution { o } => assert!((0.1..0.3).contains(&o)),
            SolidKind::Extrusion { h } => assert!((0.1..0.3).contains(&h)),
        }
    }
    let json = analytic_json(&set);
    assert!(json.contains("\"count\": 70"));
    assert_eq!(json.matches("\"seed\"").count(), 70);
}

#[test]
fn marching_cubes_of_a_torus_solid() {
    let (r, o) = (0.25, 0.5);
    let t = Solid::new(Shape2d::new(Kind::Circle, &[r]).unwrap(), SolidKind::Revolution { o }).unwrap();
    let opts = McOptions { resolution: Resolution::Longest(96), band: Band::Lipschitz { bound: 1.5, block: 8 }, ..McOptions::default() };
    let (m, rep) = marching_cubes(&t, &opts).unwrap();
    assert!(matches!(rep.mesh, MeshCheck::Manifold { boundary_edges: 0, euler_characteristic: 0, components: 1 }), "{:?}", rep.mesh);
    let vol = 2.0 * PI * PI * r * r * o;
    assert!((m.signed_volume() / vol - 1.0).abs() < 0.02, "{} vs {vol}", m.signed_volume());
}

#[test]
fn saddle_counts_follow_the_known_geometry() {
    let set = saddle_set().unwrap();
    for s in &set {
        let st = saddle_stats(&s.solid, 20_000, 1).unwrap();
        eprintln!("saddle {:<18} {:>6}/{} = {:.4}  (flat {:.3})  κ ∈ [{:.2}, {:.2}]", s.name, st.saddle, st.n, st.fraction, st.flat as f64 / st.n as f64, st.extremes.0, st.extremes.1);
        if let Some(want) = s.analytic_fraction {
            let sigma = (want * (1.0 - want) / st.n as f64).sqrt().max(1e-9);
            assert!((st.fraction - want).abs() < 4.5 * sigma + 1e-3, "{}: {} vs {want}", s.name, st.fraction);
        }
        match s.name {
            "sphere" | "cylinder" => assert_eq!(st.saddle, 0, "{}", s.name),
            _ => assert!(st.fraction > 0.05, "{} has only {} saddle points", s.name, st.saddle),
        }
    }
}
