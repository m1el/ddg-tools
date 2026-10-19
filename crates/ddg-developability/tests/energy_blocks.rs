//! Finite-difference checks of the building blocks (SPEC §10.2): Eq. 7–9, `backprop_face` with
//! random sensitivities, the log-map adjoint (φ ∈ {1e-6, 1e-3, 1, 3}), the circumcenter
//! derivative and the altitude derivative.
#![allow(clippy::needless_range_loop)]

mod energy_common;
use energy_common::*;

use ddg_core::generators::SplitMix64;
use ddg_developability::energy::{
    circumcenter_value_gradient, log_map, log_map_adjoint, log_map_g, log_map_h,
    spherical_circumcenter, tangent_basis, width_direction,
};
use ddg_developability::geometry::{
    angle_gradients, area_gradient, backprop_face, corner_angles, face_frame, normal_jacobian,
    vector_area,
};

fn rv(r: &mut SplitMix64) -> [f64; 3] {
    [
        r.uniform(-1.0, 1.0),
        r.uniform(-1.0, 1.0),
        r.uniform(-1.0, 1.0),
    ]
}

fn random_triangle(r: &mut SplitMix64) -> [[f64; 3]; 3] {
    loop {
        let p = [rv(r), rv(r), rv(r)];
        let f = face_frame(p).unwrap();
        if f.area() > 0.1 {
            return p;
        }
    }
}

/// Central difference of `f` at `p` in direction (corner `m`, coordinate `d`).
fn fd<T>(
    p: [[f64; 3]; 3],
    m: usize,
    d: usize,
    h: f64,
    f: impl Fn([[f64; 3]; 3]) -> T,
    sub: impl Fn(T, T) -> T,
) -> T {
    let mut a = p;
    let mut b = p;
    a[m][d] += h;
    b[m][d] -= h;
    sub(f(a), f(b))
}

#[test]
fn eq7_eq8_eq9() {
    let mut r = SplitMix64::new(1);
    let h = 1e-6;
    let mut worst = [0.0f64; 3];
    for _ in 0..50 {
        let p = random_triangle(&mut r);
        let fr = face_frame(p).unwrap();
        for m in 0..3 {
            let ga = area_gradient(p, fr.normal, m);
            let jn = normal_jacobian(p, &fr, m);
            let gt: Vec<[[f64; 3]; 3]> = (0..3).map(|c| angle_gradients(p, fr.normal, c)).collect();
            for d in 0..3 {
                let da =
                    fd(p, m, d, h, |q| face_frame(q).unwrap().area(), |a, b| a - b) / (2.0 * h);
                worst[0] = worst[0].max((da - ga[d]).abs());
                let dn = fd(p, m, d, h, |q| face_frame(q).unwrap().normal, sub);
                for row in 0..3 {
                    worst[1] = worst[1].max((dn[row] / (2.0 * h) - jn[row][d]).abs());
                }
                let dt = fd(p, m, d, h, corner_angles, |a, b| {
                    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
                });
                for c in 0..3 {
                    // gradient of θ_c w.r.t. corner m: angle_gradients(c) is ordered [c, c+1, c+2]
                    let slot = (m + 3 - c) % 3;
                    worst[2] = worst[2].max((dt[c] / (2.0 * h) - gt[c][slot][d]).abs());
                }
            }
        }
    }
    println!(
        "Eq7 {:.1e}  Eq8 {:.1e}  Eq9 {:.1e}",
        worst[0], worst[1], worst[2]
    );
    assert!(worst.iter().all(|&w| w < 1e-8), "{worst:?}");
}

#[test]
fn backprop_face_matches_chain_rule() {
    let mut r = SplitMix64::new(2);
    let h = 1e-6;
    let mut worst = 0.0f64;
    for _ in 0..50 {
        let p = random_triangle(&mut r);
        let (y, z) = (rv(&mut r), rv(&mut r));
        let s = [
            r.uniform(-1.0, 1.0),
            r.uniform(-1.0, 1.0),
            r.uniform(-1.0, 1.0),
        ];
        // L = y·N + z·S (S = n/2 is this face's part of a vector area) + Σ s_m θ_m
        let l = |q: [[f64; 3]; 3]| {
            let fr = face_frame(q).unwrap();
            let t = corner_angles(q);
            dot(y, fr.normal)
                + dot(z, vector_area([&fr.n]))
                + s[0] * t[0]
                + s[1] * t[1]
                + s[2] * t[2]
        };
        let g = backprop_face(p, &face_frame(p).unwrap(), y, z, s);
        for m in 0..3 {
            for d in 0..3 {
                let dl = fd(p, m, d, h, l, |a, b| a - b) / (2.0 * h);
                worst = worst.max((dl - g[m][d]).abs());
            }
        }
    }
    println!("backprop_face {worst:.1e}");
    assert!(worst < 1e-8);
}

#[test]
fn log_map_series_and_adjoint() {
    // G, H continuous across the series thresholds; H = G'/sin
    for phi in [0.9e-3f64, 1e-3, 1.1e-3, 0.9e-2, 1e-2, 1.1e-2] {
        let gd = (phi / phi.sin() - log_map_g(phi)).abs();
        assert!(gd < 1e-15, "G at {phi}: {gd:e}");
        let (sn, cs) = phi.sin_cos();
        let hd = ((sn - phi * cs) / (sn * sn * sn) - log_map_h(phi)).abs();
        assert!(hd < 1e-8, "H at {phi}: {hd:e}"); // closed form loses ~eps/phi^2
    }
    for phi in [0.5f64, 1.0, 2.0] {
        let h = 1e-6f64;
        let dg = (phi + h) / (phi + h).sin() - (phi - h) / (phi - h).sin();
        assert!((dg / (2.0 * h) / phi.sin() - log_map_h(phi)).abs() < 1e-8);
    }
    let mut r = SplitMix64::new(3);
    let mut worst = 0.0f64;
    for &phi in &[1e-6f64, 1e-3, 1.0, 3.0] {
        for _ in 0..20 {
            let n0 = normalize(rv(&mut r));
            let (t1, _) = tangent_basis(n0);
            let rot_axis = normalize(cross(n0, add(t1, scale(rv(&mut r), 0.3))));
            // N at angle phi from n0
            let perp = normalize(cross(rot_axis, n0));
            let nn = add(scale(n0, phi.cos()), scale(perp, phi.sin()));
            let lm = log_map(n0, nn);
            assert!(
                (lm.phi - phi).abs() < 1e-12
                    && (norm(lm.tilde) - phi).abs() < 1e-12 * phi.max(1e-3)
            );
            assert!(dot(lm.tilde, n0).abs() < 1e-14);
            let z = rv(&mut r);
            let (gn, gn0) = log_map_adjoint(n0, nn, &lm, z);
            // tangential directional derivatives (unit vectors stay unit)
            let h = 1e-7;
            for which in 0..2 {
                let dir = rv(&mut r);
                let at = |t: f64| {
                    let (a, b) = if which == 0 {
                        (n0, normalize(add(nn, scale(dir, t))))
                    } else {
                        (normalize(add(n0, scale(dir, t))), nn)
                    };
                    dot(z, log_map(a, b).tilde)
                };
                let fdv = (at(h) - at(-h)) / (2.0 * h);
                let base = if which == 0 { nn } else { n0 };
                let tang = sub(dir, scale(base, dot(dir, base))); // d(normalize(base + t dir))/dt
                let an = dot(if which == 0 { gn } else { gn0 }, tang);
                worst = worst.max((fdv - an).abs() / (1.0 + norm(z)));
            }
        }
    }
    println!("log map adjoint {worst:.1e}");
    assert!(worst < 1e-7);
}

#[test]
fn circumcenter_and_its_derivative() {
    let mut r = SplitMix64::new(4);
    let h = 1e-6;
    let mut worst = 0.0f64;
    for _ in 0..50 {
        let s = [
            normalize(rv(&mut r)),
            normalize(rv(&mut r)),
            normalize(rv(&mut r)),
        ];
        let (v, len) = spherical_circumcenter(s[0], s[1], s[2]).unwrap();
        let g = dot(v, s[0]);
        assert!((dot(v, s[1]) - g).abs() < 1e-14 && (dot(v, s[2]) - g).abs() < 1e-14);
        let det = dot(s[0], cross(s[1], s[2]));
        assert!((g - det / len).abs() < 1e-14);
        let grad = circumcenter_value_gradient(s[0], s[1], s[2], v, len);
        let e = |q: [[f64; 3]; 3]| {
            let (v, _) = spherical_circumcenter(q[0], q[1], q[2]).unwrap();
            dot(v, q[0]).powi(2)
        };
        for m in 0..3 {
            for d in 0..3 {
                let de = fd(s, m, d, h, e, |a, b| a - b) / (2.0 * h);
                worst = worst.max((de - grad[m][d]).abs());
            }
        }
    }
    println!("circumcenter {worst:.1e}");
    assert!(worst < 1e-8);
}

#[test]
fn altitude_derivative() {
    let mut r = SplitMix64::new(5);
    let h = 1e-7;
    let mut worst = 0.0f64;
    for _ in 0..50 {
        let n0 = normalize(rv(&mut r));
        let (t1, t2) = tangent_basis(n0);
        let w = add(
            scale(t1, r.uniform(-1.0, 1.0)),
            scale(t2, r.uniform(-1.0, 1.0)),
        );
        let dw = add(
            scale(t1, r.uniform(-1.0, 1.0)),
            scale(t2, r.uniform(-1.0, 1.0)),
        );
        let (u, wh, wl) = width_direction(n0, w).unwrap();
        let up = width_direction(n0, add(w, scale(dw, h))).unwrap().0;
        let um = width_direction(n0, sub(w, scale(dw, h))).unwrap().0;
        let fdv = scale(sub(up, um), 0.5 / h);
        let an = scale(wh, -dot(u, dw) / wl);
        worst = worst.max(norm(sub(fdv, an)) * wl);
    }
    println!("altitude {worst:.1e}");
    assert!(worst < 1e-8);
}
