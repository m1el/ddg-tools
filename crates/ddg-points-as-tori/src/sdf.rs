//! The closed-form SDF of one torus, its gradient and Laplacian (SPEC.md §7.2; paper Eq. 23–24;
//! the reference's `torus_signed_distance{,_gradient}` in `torus_distance.cpp`, points-as-tori,
//! MIT, Copyright (c) 2026 Nicole Feng). f64.
//!
//! With `x − c`, `v = ‖(x − c) × u‖`, `d = (v − R, ⟨x − c, u⟩)`:
//! `g(x) = sign⁺(r)(‖d‖ − |r|)`.
//!
//! **Singular sets** (decision D12): on the axis (`v = 0`) the gradient is the one-sided limit
//! `sign⁺(r) w₂ u` with `w = (−R, d₂)/‖·‖`; on the core circle (`‖d‖ = 0`) it is 0. Both are
//! reported by [`Singular`]; no NaN is produced.

use crate::Torus;

/// Which singular set a query hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Singular {
    /// A regular point.
    #[default]
    None,
    /// On the axis (`v = 0`).
    Axis,
    /// On the core circle (`‖d‖ = 0`).
    Core,
}

/// C++ `sign`: +1 for `x ≥ 0`.
#[inline]
fn sign_pos(x: f64) -> f64 {
    if x >= 0.0 { 1.0 } else { -1.0 }
}

#[inline]
fn parts(t: &Torus, x: [f64; 3]) -> ([f64; 3], f64, f64, f64) {
    let u = t.axis;
    let rv = [x[0] - t.center[0], x[1] - t.center[1], x[2] - t.center[2]];
    let cr = [rv[1] * u[2] - rv[2] * u[1], rv[2] * u[0] - rv[0] * u[2], rv[0] * u[1] - rv[1] * u[0]];
    let v = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
    let d2 = rv[0] * u[0] + rv[1] * u[1] + rv[2] * u[2];
    (rv, v, v - t.major, d2)
}

/// `g(x) = sign⁺(r)(‖d‖ − |r|)`: the signed distance to the torus for `R ≥ |r|`; for a spindle
/// torus (`R < |r|`) the distance to the tube of the generating circle.
///
/// [formula: pat-torus-sdf]
#[inline]
pub fn torus_sdf(t: &Torus, x: [f64; 3]) -> f64 {
    let (_, _, d1, d2) = parts(t, x);
    sign_pos(t.minor) * ((d1 * d1 + d2 * d2).sqrt() - t.minor.abs())
}

/// `g` and `∇g = sign⁺(r)(w₁ (I − uuᵀ)(x − c)/v + w₂ u)`, `w = d/‖d‖`, with the singular-set
/// rules of the module docs.
///
/// [formula: pat-torus-sdf-gradient]
#[inline]
pub fn torus_sdf_gradient(t: &Torus, x: [f64; 3]) -> (f64, [f64; 3], Singular) {
    let (rv, v, d1, d2) = parts(t, x);
    let u = t.axis;
    let nd = (d1 * d1 + d2 * d2).sqrt();
    let s = sign_pos(t.minor);
    let g = s * (nd - t.minor.abs());
    if nd == 0.0 {
        return (g, [0.0; 3], Singular::Core);
    }
    let (w1, w2) = (d1 / nd, d2 / nd);
    if v == 0.0 {
        return (g, [s * w2 * u[0], s * w2 * u[1], s * w2 * u[2]], Singular::Axis);
    }
    let gr = [
        s * (w1 * (rv[0] - d2 * u[0]) / v + w2 * u[0]),
        s * (w1 * (rv[1] - d2 * u[1]) / v + w2 * u[1]),
        s * (w1 * (rv[2] - d2 * u[2]) / v + w2 * u[2]),
    ];
    (g, gr, Singular::None)
}

/// `g`, `∇g` and `Δg = sign⁺(r)(1/‖d‖ + d₁/(‖d‖ v))` (ours: the 2D Laplacian of `‖d‖` in the
/// meridian plane plus the axisymmetric term `(1/v) ∂‖d‖/∂v`; replaces the reference's
/// placeholder).
///
/// Singular sets (SPEC resolution from WP3): on the core and on the axis `Δg` is returned as 0
/// (it is unbounded there); both are reported.
///
/// [formula: pat-torus-sdf-laplacian]
#[inline]
pub fn torus_sdf_laplacian(t: &Torus, x: [f64; 3]) -> (f64, [f64; 3], f64, Singular) {
    let (g, gr, sing) = torus_sdf_gradient(t, x);
    let (_, v, d1, d2) = parts(t, x);
    let nd = (d1 * d1 + d2 * d2).sqrt();
    let s = sign_pos(t.minor);
    let lap = match sing {
        Singular::Core | Singular::Axis => 0.0,
        Singular::None => s * (1.0 / nd + d1 / (nd * v)),
    };
    (g, gr, lap, sing)
}
