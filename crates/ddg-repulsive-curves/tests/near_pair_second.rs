//! WP16b (`SPEC.md` §17p): the second-neighbour rules of the near-pair quadrature
//! ([`SecondNeighbourRule`]): the smooth positive part and the gate, central differences of the
//! pair terms, the one-sided bound, a whole curve, and that the default is untouched.

use std::sync::Arc;

use ddg_repulsive_curves::energy::near_pair::{
    InnerTables, NearPairOptions, SecondNeighbourRule, gate_factor, joining_ends,
    pair_correction_of_kind, smooth_positive_part,
};
use ddg_repulsive_curves::step::{EnergyModel, ExactEnergy};
use ddg_repulsive_curves::{CurveNetwork, Exponents};

type P = [f64; 3];

fn exps(al: f64, be: f64) -> Exponents {
    Exponents::new(al, be).unwrap()
}

#[test]
fn smooth_positive_part_is_c2_and_bounded() {
    assert_eq!(smooth_positive_part(-1.0), (0.0, 0.0));
    assert_eq!(smooth_positive_part(0.0), (0.0, 0.0));
    assert_eq!(smooth_positive_part(1.0), (1.0, 1.0));
    assert_eq!(smooth_positive_part(3.0), (3.0, 1.0));
    let h = 1e-6;
    let mut prev = 0.0;
    for k in 1..1000 {
        let z = k as f64 / 1000.0;
        let (v, d) = smooth_positive_part(z);
        let fd = (smooth_positive_part(z + h).0 - smooth_positive_part(z - h).0) / (2.0 * h);
        assert!((fd - d).abs() < 1e-8, "z {z}: fd {fd} d {d}");
        assert!(v >= 0.0 && v <= z && v >= prev, "z {z}: {v}");
        prev = v;
    }
    // Second derivative continuous at both junctions: ψ'' → 0 from inside.
    let d2 = |z: f64| (smooth_positive_part(z + h).1 - smooth_positive_part(z - h).1) / (2.0 * h);
    assert!(d2(1e-4).abs() < 1e-2 && d2(1.0 - 1e-4).abs() < 1e-2);
    assert_eq!(gate_factor(0.0), (0.0, 0.0));
    assert_eq!(gate_factor(2.0), (1.0, 0.0));
    for x in [0.1, 0.5, 0.9] {
        let fd = (gate_factor(x + h).0 - gate_factor(x - h).0) / (2.0 * h);
        assert!((fd - gate_factor(x).1).abs() < 1e-8);
    }
}

/// Second-neighbour configurations `I = p₀p₁`, `J = q₀q₁`, joined by `K = p₁q₀`:
/// a Z-fold (I and J cross in their interiors at height `d`), a corner (closest points at K's
/// ends), a hairpin (a U-turn that opens).
fn configs(d: f64) -> Vec<(&'static str, [P; 2], [P; 2])> {
    vec![
        (
            "zfold",
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.3, 0.05, d], [1.3, -0.07, d]],
        ),
        (
            "zfold-skew",
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.35, 0.2, d], [1.1, -0.4, d + 0.1]],
        ),
        (
            "corner",
            [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.1 + d, 0.0], [0.0, 0.1 + d, 1.0]],
        ),
        (
            "end-on-interior",
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[1.0, 0.4, d], [1.0, -0.6, d]],
        ),
        (
            "hairpin",
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[1.0, 0.05 + d, 0.0], [0.2, 0.35 + d, 0.0]],
        ),
    ]
}

fn rules() -> Vec<(&'static str, NearPairOptions)> {
    let base = NearPairOptions {
        second_neighbours: true,
        ..NearPairOptions::default()
    };
    vec![
        (
            "one-sided",
            NearPairOptions {
                second_neighbour_rule: SecondNeighbourRule::OneSided,
                ..base
            },
        ),
        (
            "gated",
            NearPairOptions {
                second_neighbour_rule: SecondNeighbourRule::Gated { width: 0.1 },
                ..base
            },
        ),
        (
            "one-sided-gated",
            NearPairOptions {
                second_neighbour_rule: SecondNeighbourRule::OneSidedGated { width: 0.1 },
                ..base
            },
        ),
        (
            "one-sided-either",
            NearPairOptions {
                second_neighbour_rule: SecondNeighbourRule::OneSidedGatedEither { width: 0.1 },
                ..base
            },
        ),
        (
            "narrow",
            NearPairOptions {
                second_neighbour_band: Some((0.05, 0.25)),
                ..base
            },
        ),
    ]
}

fn pair_fd(p: [P; 2], q: [P; 2], e: Exponents, o: &NearPairOptions) -> (f64, f64, f64) {
    let tab = InnerTables::new(e);
    let f = |p: [P; 2], q: [P; 2]| {
        pair_correction_of_kind(&tab, e, o, p, q, Some((1, 0)), true, true, true)
    };
    let base = f(p, q);
    let h = 1e-7;
    let (mut num, mut den) = (0.0f64, 0.0f64);
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
        }
    }
    let err = if den > 0.0 {
        (num / den).sqrt()
    } else {
        num.sqrt()
    };
    (err, base.weight, base.value)
}

#[test]
fn second_neighbour_terms_match_central_differences() {
    for (al, be) in [(2.0, 4.5), (3.0, 6.0)] {
        let e = exps(al, be);
        for d in [3e-3, 0.03, 0.12, 0.3] {
            for (cname, p, q) in configs(d) {
                for (rname, o) in rules() {
                    let (err, w, v) = pair_fd(p, q, e, &o);
                    if w == 0.0 {
                        continue;
                    }
                    eprintln!(
                        "FD ({al},{be}) {cname} d {d} {rname}: w {w:.3} value {v:.4e} err {err:.2e}"
                    );
                    assert!(
                        err < 1e-5,
                        "({al},{be}) {cname} d {d} {rname}: FD error {err:.2e}"
                    );
                }
            }
        }
    }
}

#[test]
fn one_sided_rules_never_lower_the_pair_and_see_the_fold() {
    let e = exps(2.0, 4.5);
    let tab = InnerTables::new(e);
    for d in [1e-6, 1e-3, 0.05, 0.3] {
        for (cname, p, q) in configs(d) {
            for (rname, o) in rules() {
                let c = pair_correction_of_kind(&tab, e, &o, p, q, Some((1, 0)), true, true, false);
                if rname.starts_with("one-sided") {
                    assert!(c.value >= 0.0, "{cname} {d} {rname}: {}", c.value);
                }
                if cname == "zfold" && d <= 1e-3 && rname != "narrow" {
                    assert!(c.value > 10.0, "{cname} {d} {rname}: {}", c.value);
                }
                if cname == "end-on-interior" && d <= 1e-3 && rname == "one-sided-either" {
                    // One closest point at an end of K, the other inside J: the either gate is open.
                    assert!(c.value > 10.0, "{cname} {d} {rname}: {}", c.value);
                }
                if cname == "corner" && rname != "narrow" && rname != "one-sided" {
                    // Closest points at the ends of K: the gate is shut.
                    assert_eq!(c.value, 0.0, "{cname} {d} {rname}");
                }
            }
        }
    }
}

/// An open polyline with an accordion fold in the middle: edges `I = (3, 4)`, `K = (4, 5)`,
/// `J = (5, 6)` with `J` passing over `I` at height `gap`.
fn accordion(gap: f64) -> (CurveNetwork, Vec<P>) {
    let mut pos: Vec<P> = vec![
        [-3.0, 0.0, 0.0],
        [-2.0, 0.1, 0.0],
        [-1.0, 0.0, 0.05],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.3, 0.05, gap],
        [1.3, -0.07, gap],
    ];
    for k in 0..4 {
        pos.push([2.3 + k as f64, -0.1 + 0.05 * k as f64, gap + 0.1 * k as f64]);
    }
    let edges: Vec<[u32; 2]> = (0..pos.len() as u32 - 1).map(|i| [i, i + 1]).collect();
    (CurveNetwork::new(pos.len(), edges).unwrap(), pos)
}

fn pool() -> Arc<rayon::ThreadPool> {
    Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap(),
    )
}

#[test]
fn curve_with_a_fold_matches_central_differences() {
    let e = exps(2.0, 4.5);
    for gap in [1e-3, 0.05, 0.2] {
        let (net, pos) = accordion(gap);
        assert_eq!(joining_ends(&net, 3, 5), Some((1, 0)));
        for (rname, o) in rules() {
            let mut m = ExactEnergy::new(pool(), e).with_near_pair(Some(o));
            m.prepare(&net, &pos).unwrap();
            let ev = m.energy_and_differential(&net, &pos).unwrap();
            let h = 1e-7;
            let mut worst = 0.0f64;
            let gn = ev
                .differential
                .iter()
                .flatten()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt();
            for v in 2..8 {
                for c in 0..3 {
                    let mut p = pos.clone();
                    p[v][c] += h;
                    let up = m.energy(&net, &p).unwrap();
                    p[v][c] -= 2.0 * h;
                    let dn = m.energy(&net, &p).unwrap();
                    let fd = (up - dn) / (2.0 * h);
                    worst = worst.max((fd - ev.differential[v][c]).abs() / gn);
                }
            }
            eprintln!(
                "curve gap {gap} {rname}: near pairs {} worst {worst:.2e}",
                ev.near_pairs
            );
            assert!(worst < 1e-6, "gap {gap} {rname}: {worst:.2e}");
        }
    }
}

#[test]
fn rules_are_inert_while_second_neighbours_are_excluded() {
    let e = exps(2.0, 4.5);
    let (net, pos) = accordion(0.01);
    let reference = {
        let mut m = ExactEnergy::new(pool(), e).with_near_pair(Some(NearPairOptions::default()));
        m.prepare(&net, &pos).unwrap();
        m.energy_and_differential(&net, &pos).unwrap()
    };
    for (_, o) in rules() {
        let o = NearPairOptions {
            second_neighbours: false,
            ..o
        };
        let mut m = ExactEnergy::new(pool(), e).with_near_pair(Some(o));
        m.prepare(&net, &pos).unwrap();
        let ev = m.energy_and_differential(&net, &pos).unwrap();
        assert_eq!(ev.energy.to_bits(), reference.energy.to_bits());
        for (a, b) in ev.differential.iter().zip(&reference.differential) {
            for c in 0..3 {
                assert_eq!(a[c].to_bits(), b[c].to_bits());
            }
        }
    }
}

#[test]
fn options_are_validated() {
    let bad = [
        NearPairOptions {
            second_neighbour_band: Some((0.3, 0.2)),
            ..Default::default()
        },
        NearPairOptions {
            second_neighbour_rule: SecondNeighbourRule::Gated { width: 0.0 },
            ..Default::default()
        },
        NearPairOptions {
            second_neighbour_rule: SecondNeighbourRule::OneSidedGated { width: 1.5 },
            ..Default::default()
        },
    ];
    for o in bad {
        assert!(o.validate().is_err(), "{o:?}");
    }
}
