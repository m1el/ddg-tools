//! Exact values (SPEC §10.3): regular cones, hinges, flat stars, small boundary stars, valence 3.

mod energy_common;
use energy_common::*;

use std::f64::consts::PI;

use ddg_core::TriMesh;
use ddg_core::generators::SplitMix64;
use ddg_developability::{Energy, EnergyOptions, evaluate};

fn center_only(energy: Energy) -> EnergyOptions {
    EnergyOptions {
        energy,
        include_boundary: false,
        omit_valence_three: false,
        threads: Some(1),
    }
}

fn all_vertices(energy: Energy) -> EnergyOptions {
    EnergyOptions {
        energy,
        include_boundary: true,
        omit_valence_three: false,
        threads: Some(1),
    }
}

/// Closed forms of SPEC §10.3 for the regular cone of valence `k` and apex height `h`;
/// `None` where only oracle values exist (λ_max for odd k).
fn cone_closed_form(energy: Energy, k: usize, h: f64) -> Option<f64> {
    let kf = k as f64;
    let c = (PI / kf).cos();
    let sin2a = h * h / (h * h + c * c);
    let cos2a = 1.0 - sin2a;
    let alpha = h.atan2(c);
    let theta = 2.0 * ((PI / kf).sin() / (1.0 + h * h).sqrt()).asin();
    Some(match energy {
        Energy::Covariance => theta * kf * (sin2a / 2.0).min(cos2a),
        Energy::Combinatorial => {
            let rho = |m: usize| (m as f64 * PI / kf).sin() / (m as f64 * (PI / kf).sin());
            let m = (1..k)
                .map(|m| 2.0 - rho(m).powi(2) - rho(k - m).powi(2))
                .fold(f64::INFINITY, f64::min);
            sin2a * m
        }
        Energy::CovarianceMax => {
            if !k.is_multiple_of(2) {
                return None;
            }
            cos2a.min(sin2a * c * c)
        }
        Energy::Intrinsic => theta * kf * alpha * alpha / 2.0,
        Energy::IntrinsicMax => {
            if k.is_multiple_of(2) {
                (2.0 * alpha * c).powi(2)
            } else {
                (alpha * (1.0 + c)).powi(2)
            }
        }
        Energy::CombinatorialMax => return None,
    })
}

#[test]
fn regular_cones_match_closed_forms() {
    let mut checked = 0;
    for k in 3..=8 {
        for h in [0.1, 0.6, 2.0] {
            let (m, p) = regular_cone(k, h);
            for energy in Energy::ALL {
                let Some(want) = cone_closed_form(energy, k, h) else {
                    continue;
                };
                let e = evaluate(&m, &p, &center_only(energy)).unwrap();
                let got = e.vertex_energy[0];
                assert!(
                    (got - want).abs() <= 1e-12 * want.abs().max(1e-300),
                    "{energy:?} k={k} h={h}: got {got:.15e}, closed form {want:.15e}"
                );
                assert_eq!(e.energy, got);
                checked += 1;
            }
        }
    }
    assert!(checked > 70);
}

/// Values of the SPEC author's independent NumPy check (12 significant digits) for the cases
/// without a closed form: `λ_max` at odd valence.
#[test]
fn regular_cones_odd_covariance_max_match_numpy() {
    let table = [
        (5, 0.1, 0.012310175987),
        (5, 0.6, 0.288867938521),
        (5, 2.0, 0.140618176459),
        (7, 0.1, 0.010993598769),
        (7, 0.6, 0.277260164000),
        (7, 2.0, 0.168700734898),
    ];
    for (k, h, want) in table {
        let (m, p) = regular_cone(k, h);
        let got = evaluate(&m, &p, &center_only(Energy::CovarianceMax))
            .unwrap()
            .vertex_energy[0];
        assert!(
            (got - want).abs() <= 1e-12,
            "k={k} h={h}: {got:.12} vs {want:.12}"
        );
    }
}

/// `π_max` of a regular cone: the best split into two arcs of `m` and `k − m` faces has the
/// largest chord of the larger arc: `4 sin²α sin²((m_max − 1)π/k)` with `m_max = ⌈k/2⌉`
/// (chord between normals `j` apart: `|N_a − N_b|² = 4 sin²α sin²(jπ/k)`).
#[test]
fn regular_cones_combinatorial_max() {
    for k in 3..=8 {
        for h in [0.1, 0.6, 2.0] {
            let (m, p) = regular_cone(k, h);
            let c = (PI / k as f64).cos();
            let sin2a = h * h / (h * h + c * c);
            let j = k.div_ceil(2) - 1;
            let want = 4.0 * sin2a * (j as f64 * PI / k as f64).sin().powi(2);
            let got = evaluate(&m, &p, &center_only(Energy::CombinatorialMax))
                .unwrap()
                .vertex_energy[0];
            assert!(
                (got - want).abs() <= 1e-12 * want,
                "k={k} h={h}: {got} vs {want}"
            );
        }
    }
}

/// A star whose ring lies on two half-planes through the x axis: `z = 0` for `y ≥ 0` and
/// `z = y/2` for `y < 0`, with ring vertices on the x axis at both ends of each half. All
/// normals are exactly orthogonal to x (the factor 1/2 is exact), so every energy is ≤ 1e-28.
fn exact_hinge(n_up: usize, n_down: usize, seed: u64) -> Mesh {
    let mut r = SplitMix64::new(seed);
    let mut pos = vec![[0.0, 0.0, 0.0]];
    let mut ang = vec![0.0];
    for i in 1..n_up {
        ang.push(PI * (i as f64 + r.uniform(-0.3, 0.3)) / n_up as f64);
    }
    ang.push(PI);
    for i in 1..n_down {
        ang.push(PI + PI * (i as f64 + r.uniform(-0.3, 0.3)) / n_down as f64);
    }
    let k = ang.len();
    for (i, &t) in ang.iter().enumerate() {
        let rad = r.uniform(0.5, 1.5);
        let (x, y) = if i == 0 {
            (rad, 0.0)
        } else if t == PI {
            (-rad, 0.0)
        } else {
            (rad * t.cos(), rad * t.sin())
        };
        let z = if y < 0.0 { 0.5 * y } else { 0.0 };
        pos.push([x, y, z]);
    }
    let faces = (0..k)
        .map(|m| [0, 1 + m as u32, 1 + ((m + 1) % k) as u32])
        .collect();
    (TriMesh::new(k + 1, faces).unwrap(), pos)
}

#[test]
fn hinge_and_flat_stars_vanish() {
    for (n_up, n_down) in [(2, 2), (2, 3), (3, 4), (2, 5)] {
        for seed in 0..3 {
            let (m, p) = exact_hinge(n_up, n_down, 31 * seed + n_up as u64);
            for energy in Energy::ALL {
                let e = evaluate(&m, &p, &center_only(energy)).unwrap();
                assert!(
                    e.vertex_energy[0].abs() <= 1e-28,
                    "hinge {n_up}+{n_down} {energy:?}: {:e}",
                    e.vertex_energy[0]
                );
            }
        }
    }
    // rotated hinge: rounding level only
    let rot = random_rotation(5);
    let (m, p) = exact_hinge(3, 3, 9);
    let p: Vec<_> = p.iter().map(|&x| rotate(&rot, x)).collect();
    for energy in Energy::ALL {
        let e = evaluate(&m, &p, &center_only(energy)).unwrap();
        assert!(
            e.vertex_energy[0].abs() <= 1e-14,
            "rotated hinge {energy:?}: {:e}",
            e.vertex_energy[0]
        );
    }
    // flat stars and a flat grid (all vertices, boundary included)
    // (random rings for k >= 4 keep every sector below pi; k = 3 from the regular cone)
    for k in 3..=8 {
        let (m, mut p) = if k == 3 {
            regular_cone(3, 0.0)
        } else {
            interior_star(k, 77 + k as u64, 0.0)
        };
        p[0][2] = 0.0;
        for energy in Energy::ALL {
            let e = evaluate(&m, &p, &all_vertices(energy)).unwrap();
            for (v, x) in e.vertex_energy.iter().enumerate() {
                assert!(
                    x.abs() <= 1e-28,
                    "flat star k={k} {energy:?} vertex {v}: {x:e}"
                );
            }
        }
    }
}

#[test]
fn small_boundary_stars_vanish() {
    for k in 1..=2 {
        for seed in 0..5 {
            let (m, p) = boundary_fan(k, 500 + seed, 0.8);
            for energy in Energy::ALL {
                let e = evaluate(&m, &p, &all_vertices(energy)).unwrap();
                for (v, x) in e.vertex_energy.iter().enumerate() {
                    assert!(
                        x.abs() <= 1e-15,
                        "{k}-face fan {energy:?} vertex {v}: {x:e}"
                    );
                }
                assert!(e.vertex_margin.iter().all(|&g| g == f64::INFINITY));
                assert_eq!(e.near_ties, 0);
            }
        }
    }
}

#[test]
fn valence_three_is_omitted_by_default() {
    // A tetrahedron: all four vertices are interior with valence 3.
    let m = TriMesh::new(4, vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]]).unwrap();
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.1, 0.0],
        [0.2, 1.0, 0.1],
        [0.3, 0.2, 1.0],
    ];
    for energy in Energy::ALL {
        let e = evaluate(&m, &p, &EnergyOptions::new(energy)).unwrap();
        assert!(e.counted.iter().all(|&c| !c));
        assert_eq!(e.energy, 0.0);
        assert!(e.gradient.iter().all(|g| *g == [0.0; 3]));
        assert!(
            e.vertex_energy.iter().all(|x| x.is_finite() && *x > 0.0),
            "{energy:?}: {:?}",
            e.vertex_energy
        );
        let all = evaluate(&m, &p, &all_vertices(energy)).unwrap();
        assert_eq!(all.vertex_energy, e.vertex_energy);
        assert!(all.energy > 0.0);
    }
    // valence-3 interior vertex of a star mesh
    let (m, p) = interior_star(3, 4, 0.5);
    let e = evaluate(&m, &p, &EnergyOptions::new(Energy::Covariance)).unwrap();
    assert!(!e.counted[0] && e.counted[1]);
    assert!(e.vertex_energy[0] > 0.0);
}

/// `λ_max` with pruning equals the plain enumeration of SPEC §3.5 (bitwise, value and margin)
/// and is never beaten by a dense search over the sphere.
#[test]
fn covariance_max_matches_plain_enumeration_and_dense_search() {
    use ddg_developability::energy::{Scratch, Star, star_energy};
    let mut r = SplitMix64::new(77);
    for trial in 0..300 {
        let k = 3 + trial % 7;
        let h = [0.05, 0.4, 1.5][trial % 3];
        let normals: Vec<[f64; 3]> = (0..k)
            .map(|_| {
                normalize([
                    r.uniform(-h, h),
                    r.uniform(-h, h),
                    1.0 - 2.0 * (trial % 2) as f64 * r.uniform(0.0, 1.0),
                ])
            })
            .collect();
        let angles = vec![1.0; k];
        let star = Star {
            normals: &normals,
            angles: &angles,
            closed: true,
            vector_area: [0.0, 0.0, 1.0],
        };
        let (mut y, mut s) = (vec![[0.0; 3]; k], vec![0.0; k]);
        let out = star_energy(
            Energy::CovarianceMax,
            &star,
            &mut y,
            &mut s,
            &mut Scratch::default(),
        )
        .unwrap();
        // plain enumeration
        let phi = |v: [f64; 3]| {
            normals
                .iter()
                .map(|n| dot(v, *n).powi(2))
                .fold(0.0, f64::max)
        };
        let mut keys: Vec<(f64, f64)> = Vec::new();
        let mut cand = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            let n = cross(sub(b, a), sub(c, a));
            let len = norm(n);
            if len.is_finite() && len > 0.0 {
                let v = scale(n, 1.0 / len);
                let p = phi(v);
                keys.push((p, p - dot(v, a).powi(2)));
            }
        };
        for i in 0..k {
            for j in i + 1..k {
                for l in j + 1..k {
                    for (s2, s3) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                        cand(normals[i], scale(normals[j], s2), scale(normals[l], s3));
                    }
                }
            }
        }
        for i in 0..k {
            for j in i + 1..k {
                cand(normals[i], scale(normals[i], -1.0), normals[j]);
            }
        }
        let mut ps: Vec<f64> = keys.iter().map(|k| k.0).collect();
        ps.sort_by(f64::total_cmp);
        assert_eq!(out.value.to_bits(), ps[0].to_bits(), "trial {trial}");
        let margin = if ps.len() > 1 {
            ps[1] - ps[0]
        } else {
            f64::INFINITY
        };
        assert_eq!(out.margin.to_bits(), margin.to_bits(), "trial {trial}");
        // dense search: no direction is better
        let mut dense = f64::INFINITY;
        for a in 0..200 {
            for b in 0..100 {
                let (t, p) = (2.0 * PI * a as f64 / 200.0, PI * (b as f64 + 0.5) / 100.0);
                dense = dense.min(phi([p.sin() * t.cos(), p.sin() * t.sin(), p.cos()]));
            }
        }
        assert!(
            out.value <= dense + 1e-12,
            "trial {trial}: {} > dense {dense}",
            out.value
        );
    }
}
