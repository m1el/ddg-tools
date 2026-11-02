//! Generator tests (WP12): valid manifold surfaces of the expected genus, no self-intersections,
//! linking numbers of the handcuff cores, determinism, face counts.

use ddg_repulsive_surfaces::Orientation;
use ddg_repulsive_surfaces::generate::check::{linking_number, self_intersections};
use ddg_repulsive_surfaces::generate::perturbed::PerturbedTorus;
use ddg_repulsive_surfaces::generate::primitives::*;
use ddg_repulsive_surfaces::generate::shapes::*;
use ddg_repulsive_surfaces::generate::Mesh;

/// Closed, one component, consistently oriented outward, genus `g`, no self-intersection.
fn check_closed(m: &Mesh, genus: i64, what: &str) {
    let s = m.surface().unwrap_or_else(|e| panic!("{what}: {e}"));
    assert_eq!(s.n_components(), 1, "{what}: components");
    assert!(s.boundary_loops().is_empty(), "{what}: boundary");
    assert_eq!(s.euler_characteristic(), 2 - 2 * genus, "{what}: chi");
    assert!(matches!(s.orientation(0), Orientation::Consistent), "{what}: orientation");
    assert!(m.volume() > 0.0, "{what}: outward");
    let x = self_intersections(m, 5);
    assert!(x.is_empty(), "{what}: self-intersections {x:?}");
    // no sliver triangles: the smallest area is a sane fraction of the mean
    let mean = m.area() / m.faces.len() as f64;
    let min = m
        .faces
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|i| m.positions[i as usize]);
            let u = [0, 1, 2].map(|k| b[k] - a[k]);
            let v = [0, 1, 2].map(|k| c[k] - a[k]);
            0.5 * ((u[1] * v[2] - u[2] * v[1]).powi(2)
                + (u[2] * v[0] - u[0] * v[2]).powi(2)
                + (u[0] * v[1] - u[1] * v[0]).powi(2))
            .sqrt()
        })
        .fold(f64::INFINITY, f64::min);
    assert!(min > 1e-4 * mean, "{what}: sliver, min area {min} vs mean {mean}");
}

#[test]
fn primitives() {
    check_closed(&sphere(2, 1.5, [1.0, 0.0, 0.0]).unwrap(), 0, "sphere");
    check_closed(&torus(24, 12, 2.0, 0.7).unwrap(), 1, "torus");
    check_closed(&box_mesh([1.0, 2.0, 0.5], 3).unwrap(), 0, "box");
    check_closed(&capsule(0.5, 2.0, 16, 6, 4).unwrap(), 0, "capsule");
    assert_eq!(box_mesh([1.0; 3], 3).unwrap().faces.len(), 12 * 9);
    assert!(box_mesh([1.0; 3], 0).is_err());
}

#[test]
fn punctured_torus_has_one_hole() {
    let m = punctured_torus(40, 20, 2.0, 0.7, 0.5).unwrap();
    let s = m.surface().unwrap();
    assert_eq!(s.boundary_loops().len(), 1);
    assert_eq!(s.euler_characteristic(), -1);
    assert!(self_intersections(&m, 1).is_empty());
}

#[test]
fn genus_series() {
    for g in 0..=11u32 {
        let m = genus_surface(g, 1500 + 700 * g as usize).unwrap();
        check_closed(&m, g as i64, &format!("genus {g}"));
    }
    let a = genus_surface(3, 3000).unwrap();
    let b = genus_surface(3, 3000).unwrap();
    assert_eq!(a, b, "deterministic");
}

#[test]
fn handcuffs_linking_and_genus() {
    for k in [1u32, 2, 6] {
        let t = handcuffs(k, 25_000).unwrap();
        check_closed(&t.mesh, 2, &format!("handcuffs {k}"));
        let lk = linking_number(&t.cores[0], &t.cores[1]);
        assert!((lk.abs() - k as f64).abs() < 0.02, "k={k}: linking number {lk}");
    }
}

#[test]
fn knotted_tunnels_are_tori() {
    for (knot, name) in [(Knot::Trefoil, "trefoil"), (Knot::FigureEight, "figure-8")] {
        let t = knotted_tunnel(knot, 30_000).unwrap();
        check_closed(&t.mesh, 1, name);
    }
}

#[test]
fn fig12_analogues_have_the_papers_face_counts() {
    let results: Vec<_> = (0..5)
        .map(|i| std::thread::spawn(move || fig12(i).unwrap()))
        .collect::<Vec<_>>()
        .into_iter()
        .map(|h| h.join().unwrap())
        .collect();
    for (i, t) in results.iter().enumerate() {
        let want = FIG12_FACES[i] as f64;
        let got = t.mesh.faces.len() as f64;
        eprintln!("fig12 #{i}: {got} faces (paper {want})");
        assert!((got - want).abs() / want < 0.05, "#{i}: {got} vs {want}");
        check_closed(&t.mesh, if i < 3 { 2 } else { 1 }, &format!("fig12 #{i}"));
    }
}

#[test]
fn perturbed_torus() {
    let f = PerturbedTorus::new(2.0, 0.8, 2, 0.15, 7).unwrap();
    let m = f.sample(60, 24, 0.3, 0.0).unwrap();
    check_closed(&m, 1, "perturbed torus");
    let m2 = PerturbedTorus::new(2.0, 0.8, 2, 0.15, 7).unwrap().sample(60, 24, 0.3, 0.0).unwrap();
    assert_eq!(m, m2);
    let m3 = PerturbedTorus::new(2.0, 0.8, 2, 0.15, 8).unwrap().sample(60, 24, 0.3, 0.0).unwrap();
    assert_ne!(m, m3);
    // periodic
    let p = f.point(0.3, 1.1);
    let q = f.point(0.3 + std::f64::consts::TAU, 1.1 - std::f64::consts::TAU);
    assert!((0..3).all(|k| (p[k] - q[k]).abs() < 1e-12));
}

#[test]
fn text_soup_cobblestones() {
    let h = h_shape_coarse();
    assert_eq!((h.positions.len(), h.faces.len()), (32, 60));
    let s = h.surface().unwrap();
    assert_eq!(s.euler_characteristic(), 2);
    assert!(h.volume() > 0.0);
    assert!((h.volume() - 9.0).abs() < 1e-12, "H area 9 x depth 1: {}", h.volume());

    let base = sphere(3, 1.0, [0.0; 3]).unwrap();
    let soup = soup_with_holes(&base, 0.2, 5).unwrap();
    let frac = 1.0 - soup.len() as f64 / base.faces.len() as f64;
    assert!((frac - 0.2).abs() < 0.05, "{frac}");
    assert_eq!(soup, soup_with_holes(&base, 0.2, 5).unwrap());

    let c = cobblestones(12, 0.4, 0.1, [2.0, 2.0, 1.0], 1, 3).unwrap();
    let s = c.spheres.surface().unwrap();
    assert_eq!(s.n_components(), 12);
    assert!(self_intersections(&c.spheres, 1).is_empty());
    check_closed(&c.slab, 0, "slab");
    assert!(cobblestones(500, 0.4, 0.1, [2.0, 2.0, 1.0], 1, 3).is_err());
}
