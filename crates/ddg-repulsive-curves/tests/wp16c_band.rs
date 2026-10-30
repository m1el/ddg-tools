//! WP16c (`SPEC.md` §17r): the second-neighbour band near contact
//! (`NearPairOptions::second_neighbour_band = Some((0.001, 0.005))` with the one-sided rule), the
//! candidate against the hard-flow regression of the Q10 defaults. A hairpin at a short joining
//! edge (the geometry that slows the hard flows, §17r) is left to Eq. 18; a Z-fold near contact
//! (the geometry that ended the soft runs of §17p on the contact floor) is corrected as by the
//! Q10 default; inside the narrow band the term matches central differences.

use ddg_repulsive_curves::Exponents;
use ddg_repulsive_curves::energy::near_pair::{
    InnerTables, NearPairOptions, pair_correction_of_kind,
};

type P = [f64; 3];

/// The Q10 rule, spelled out (so the tests do not depend on the default).
fn q10_rule() -> NearPairOptions {
    NearPairOptions {
        second_neighbours: true,
        second_neighbour_rule:
            ddg_repulsive_curves::energy::near_pair::SecondNeighbourRule::OneSided,
        second_neighbour_band: None,
        ..NearPairOptions::default()
    }
}

fn near_contact() -> NearPairOptions {
    NearPairOptions {
        second_neighbour_band: Some((0.001, 0.005)),
        ..q10_rule()
    }
}

/// `I = p₀p₁` (length 3, along −x), `K = p₁q₀` (length ≈ 0.28, turning by ≈ 170° at `p₁`),
/// `J = q₀q₁` continuing `K` (length 1.9): `J` runs back along `I` at distance `gap` (§17r,
/// trefoil0073).
fn hairpin(gap: f64) -> ([P; 2], [P; 2]) {
    let q0 = [(0.28f64 * 0.28 - gap * gap).sqrt(), gap, 0.0];
    let n = (q0[0] * q0[0] + q0[1] * q0[1]).sqrt();
    let q1 = [q0[0] + 1.9 * q0[0] / n, q0[1] + 1.9 * q0[1] / n, 0.0];
    ([[3.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [q0, q1])
}

/// A Z-fold (`I` and `J` parallel and overlapping, `K` folded back between them) at height `h`.
fn zfold(h: f64) -> ([P; 2], [P; 2]) {
    (
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.3, 0.05, h], [1.3, -0.07, h]],
    )
}

#[test]
fn the_band_leaves_hairpins_to_eq18_and_catches_folds_at_contact() {
    let e = Exponents::new(2.0, 4.5).unwrap();
    let t = InnerTables::new(e);
    let q10 = q10_rule();
    let band = near_contact();
    // K = p₁q₀: I meets K at its end 1, J at its end 0.
    let joined = Some((1, 0));
    for gap in [0.02, 0.05, 0.1] {
        let (p, q) = hairpin(gap);
        let a = pair_correction_of_kind(&t, e, &q10, p, q, joined, true, true, true);
        let b = pair_correction_of_kind(&t, e, &band, p, q, joined, true, true, true);
        assert!(a.value > 0.0, "gap {gap}: the Q10 rule sees the hairpin");
        assert!(a.ratio > 0.005, "gap {gap}: r {}", a.ratio);
        assert_eq!(b.value, 0.0, "gap {gap}");
        assert!(b.grad.iter().flatten().all(|&g| g == 0.0));
    }
    for h in [1e-6, 1e-4, 5e-4] {
        let (p, q) = zfold(h);
        let a = pair_correction_of_kind(&t, e, &q10, p, q, joined, true, true, true);
        let b = pair_correction_of_kind(&t, e, &band, p, q, joined, true, true, true);
        assert!(a.ratio <= 0.001, "h {h}: r {}", a.ratio);
        assert!(
            b.value > 0.0 && b.value == a.value,
            "h {h}: {} {}",
            a.value,
            b.value
        );
    }
}

#[test]
fn the_narrow_band_matches_central_differences() {
    let e = Exponents::new(2.0, 4.5).unwrap();
    let t = InnerTables::new(e);
    let band = near_contact();
    let joined = Some((1, 0));
    // Heights inside the blend of the band (0.001 < r < 0.005).
    for h in [0.0015, 0.003, 0.0045] {
        let (p, q) = zfold(h);
        let c = pair_correction_of_kind(&t, e, &band, p, q, joined, true, true, true);
        assert!(c.weight > 0.0 && c.weight < 1.0, "h {h}: w {}", c.weight);
        let step = 1e-7;
        let (mut max_err, mut scale) = (0.0f64, 0.0f64);
        for v in 0..4 {
            for k in 0..3 {
                let shifted = |s: f64| {
                    let (mut p2, mut q2) = (p, q);
                    if v < 2 {
                        p2[v][k] += s;
                    } else {
                        q2[v - 2][k] += s;
                    }
                    pair_correction_of_kind(&t, e, &band, p2, q2, joined, true, true, false).value
                };
                let fd = (shifted(step) - shifted(-step)) / (2.0 * step);
                max_err = max_err.max((fd - c.grad[v][k]).abs());
                scale = scale.max(c.grad[v][k].abs());
            }
        }
        assert!(max_err <= 1e-5 * scale, "h {h}: error {max_err} of {scale}");
    }
}

fn fold_gated() -> NearPairOptions {
    NearPairOptions {
        second_neighbour_rule:
            ddg_repulsive_curves::energy::near_pair::SecondNeighbourRule::fold_gated(),
        ..NearPairOptions::default()
    }
}

/// A fold between a Z-fold and a hairpin: `J` leaves `K` at the angle `φ` to `I`
/// (oriented toward `K`), at height `h` above `I`.
fn bent(phi_deg: f64, h: f64) -> ([P; 2], [P; 2]) {
    let phi = phi_deg.to_radians();
    let q0 = [0.7, 0.02, h];
    let q1 = [
        q0[0] + phi.cos(),
        q0[1] + 0.3 * phi.sin(),
        q0[2] + 0.95 * phi.sin(),
    ];
    ([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [q0, q1])
}

#[test]
fn the_fold_gate_keeps_z_folds_drops_hairpins_and_matches_central_differences() {
    let e = Exponents::new(2.0, 4.5).unwrap();
    let t = InnerTables::new(e);
    let q10 = q10_rule();
    let fg = fold_gated();
    fg.validate().unwrap();
    let joined = Some((1, 0));
    // Hairpins outside the contact band: no correction; Z-folds: the Q10 term.
    for gap in [0.02, 0.05] {
        let (p, q) = hairpin(gap);
        let c = pair_correction_of_kind(&t, e, &fg, p, q, joined, true, true, true);
        assert_eq!(c.value, 0.0, "hairpin gap {gap}");
    }
    for h in [1e-4, 0.01, 0.05] {
        let (p, q) = zfold(h);
        let a = pair_correction_of_kind(&t, e, &q10, p, q, joined, true, true, false);
        let b = pair_correction_of_kind(&t, e, &fg, p, q, joined, true, true, false);
        assert!(a.value > 0.0);
        assert!(
            (b.value - a.value).abs() <= 1e-12 * a.value,
            "zfold {h}: {} {}",
            a.value,
            b.value
        );
    }
    // Inside the gate (60° < φ < 90°) and the contact band: never above the one-sided term,
    // never negative, and the differential matches central differences.
    let step = 1e-7;
    for (phi, h) in [(70.0, 0.02), (80.0, 0.003), (75.0, 0.004), (30.0, 0.02)] {
        let (p, q) = bent(phi, h);
        let a = pair_correction_of_kind(&t, e, &q10, p, q, joined, true, true, false);
        let c = pair_correction_of_kind(&t, e, &fg, p, q, joined, true, true, true);
        assert!(
            c.value >= 0.0 && c.value <= a.value * (1.0 + 1e-12),
            "φ {phi}"
        );
        let (mut max_err, mut scale) = (0.0f64, 0.0f64);
        for v in 0..4 {
            for k in 0..3 {
                let shifted = |s: f64| {
                    let (mut p2, mut q2) = (p, q);
                    if v < 2 {
                        p2[v][k] += s;
                    } else {
                        q2[v - 2][k] += s;
                    }
                    pair_correction_of_kind(&t, e, &fg, p2, q2, joined, true, true, false).value
                };
                let fd = (shifted(step) - shifted(-step)) / (2.0 * step);
                max_err = max_err.max((fd - c.grad[v][k]).abs());
                scale = scale.max(c.grad[v][k].abs());
            }
        }
        assert!(c.value > 0.0, "φ {phi} h {h}: the probe should be active");
        assert!(
            max_err <= 1e-5 * scale,
            "φ {phi} h {h}: error {max_err} of {scale}"
        );
    }
}
