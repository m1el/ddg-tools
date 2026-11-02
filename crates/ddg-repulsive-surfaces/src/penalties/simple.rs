//! The geometric penalties of `SPEC.md` §5.2: total area and volume (and their deviations),
//! squared error, boundary length and boundary curvature. Each function returns the
//! unweighted value and its exact gradient (`rs-penalties`).

use crate::error::RsError;
use crate::geometry::v3::{add, axpy, cross, dot, norm, scale, sub};
use crate::geometry::{corners, face_area_gradient, face_quantities};
use crate::kernel::ipow;
use crate::surface::{Orientation, Surface};

/// Total area `A = Σ a_S` and its gradient `∂A/∂f_v = Σ_{S∋v} ½ N_S × e_opp`.
///
/// [formula: rs-penalties]
pub fn area_and_gradient(surface: &Surface, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
    let mut a = 0.0;
    let mut g = vec![[0.0; 3]; surface.n_vertices()];
    for t in surface.faces() {
        let p = corners(positions, t);
        a += face_quantities(p).0;
        let gr = face_area_gradient(p);
        for c in 0..3 {
            let v = t[c] as usize;
            g[v] = add(g[v], gr[c]);
        }
    }
    (a, g)
}

/// Checks that every component is consistently oriented (signed volume is defined).
///
/// # Errors
/// [`RsError::VolumeNeedsOrientation`] naming the component and an edge whose faces disagree.
pub fn require_orientation(surface: &Surface) -> Result<(), RsError> {
    for (c, comp) in surface.components().iter().enumerate() {
        if *comp.orientation() != Orientation::Consistent {
            return Err(RsError::VolumeNeedsOrientation {
                component: c as u32,
                edge: comp.inconsistent_edge().unwrap_or([0, 0]),
            });
        }
    }
    Ok(())
}

/// Signed volume `V = ⅙ Σ_S f_i·(f_j × f_k)` (of the cone over the origin for open
/// components) and its gradient (`∂V/∂f_i = ⅙ f_j × f_k`, cyclically). Needs consistent
/// orientation ([`require_orientation`]).
///
/// [formula: rs-volume-gradient]
pub fn volume_and_gradient(surface: &Surface, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
    let mut vol = 0.0;
    let mut g = vec![[0.0; 3]; surface.n_vertices()];
    for t in surface.faces() {
        let p = corners(positions, t);
        vol += dot(p[0], cross(p[1], p[2])) / 6.0;
        for c in 0..3 {
            let v = t[c] as usize;
            g[v] = axpy(g[v], 1.0 / 6.0, cross(p[(c + 1) % 3], p[(c + 2) % 3]));
        }
    }
    (vol, g)
}

/// `x^power` and `power x^{power − 1}` (integer powers by binary powering).
pub fn power_and_derivative(x: f64, power: f64) -> (f64, f64) {
    if power.fract() == 0.0 && (1.0..=64.0).contains(&power) {
        let n = power as u32;
        (ipow(x, n), power * ipow(x, n - 1))
    } else {
        (x.powf(power), power * x.powf(power - 1.0))
    }
}

/// `Σ_v |f_v − f_v⁰|²` and its gradient `2 (f − f⁰)`.
///
/// # Errors
/// [`RsError::PositionCountMismatch`] if `targets` has the wrong length.
pub fn squared_error_and_gradient(
    positions: &[[f64; 3]],
    targets: &[[f64; 3]],
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    if targets.len() != positions.len() {
        return Err(RsError::PositionCountMismatch {
            expected: positions.len(),
            found: targets.len(),
        });
    }
    let mut val = 0.0;
    let mut g = Vec::with_capacity(positions.len());
    for (p, q) in positions.iter().zip(targets) {
        let d = sub(*p, *q);
        val += dot(d, d);
        g.push(scale(d, 2.0));
    }
    Ok((val, g))
}

/// Total length `L_∂` of the boundary edges and its gradient (unit edge tangents; a
/// zero-length edge contributes a zero gradient).
pub fn boundary_length_and_gradient(
    surface: &Surface,
    positions: &[[f64; 3]],
) -> (f64, Vec<[f64; 3]>) {
    let mut l = 0.0;
    let mut g = vec![[0.0; 3]; surface.n_vertices()];
    for (e, &[a, b]) in surface.edges().iter().enumerate() {
        if !surface.is_boundary_edge(e as u32) {
            continue;
        }
        let d = sub(positions[a as usize], positions[b as usize]);
        let len = norm(d);
        l += len;
        if len > 0.0 {
            let u = scale(d, 1.0 / len);
            g[a as usize] = add(g[a as usize], u);
            g[b as usize] = sub(g[b as usize], u);
        }
    }
    (l, g)
}

/// The penalty `(L − L_∂)²` and its gradient.
pub fn boundary_length_penalty(
    surface: &Surface,
    positions: &[[f64; 3]],
    target: f64,
) -> (f64, Vec<[f64; 3]>) {
    let (l, mut g) = boundary_length_and_gradient(surface, positions);
    let d = target - l;
    for x in &mut g {
        *x = scale(*x, -2.0 * d);
    }
    (d * d, g)
}

/// `Σ_{v∈∂M} θ_v² / ℓ_v` (`dual_length`) or `Σ θ_v²`, with `θ_v` the turning angle between the
/// two boundary edges at `v` and `ℓ_v` half their total length, and its gradient. A straight
/// vertex (`θ = 0`) and a cusp (`θ = π`) have zero angle gradient.
///
/// [formula: rs-penalties]
pub fn boundary_curvature_and_gradient(
    surface: &Surface,
    positions: &[[f64; 3]],
    dual_length: bool,
) -> (f64, Vec<[f64; 3]>) {
    let mut val = 0.0;
    let mut g = vec![[0.0; 3]; surface.n_vertices()];
    for lp in surface.boundary_loops() {
        let m = lp.len();
        for i in 0..m {
            let (ip, iv, inx) = (
                lp[(i + m - 1) % m] as usize,
                lp[i] as usize,
                lp[(i + 1) % m] as usize,
            );
            let u = sub(positions[iv], positions[ip]);
            let w = sub(positions[inx], positions[iv]);
            let (lu, lw) = (norm(u), norm(w));
            let c = cross(u, w);
            let sin = norm(c);
            let theta = sin.atan2(dot(u, w));
            // Angle gradients with respect to u and w.
            let (mut gu, mut gw) = ([0.0; 3], [0.0; 3]);
            if sin > 0.0 {
                let n = scale(c, 1.0 / sin);
                gu = scale(cross(n, u), -1.0 / (lu * lu));
                gw = scale(cross(n, w), 1.0 / (lw * lw));
            }
            let (ell, mut dell_u, mut dell_w) = if dual_length {
                (
                    0.5 * (lu + lw),
                    scale(u, 0.5 / lu),
                    scale(w, 0.5 / lw),
                )
            } else {
                (1.0, [0.0; 3], [0.0; 3])
            };
            if !dual_length {
                dell_u = [0.0; 3];
                dell_w = [0.0; 3];
            }
            val += theta * theta / ell;
            let a = 2.0 * theta / ell;
            let b = theta * theta / (ell * ell);
            let du = sub(scale(gu, a), scale(dell_u, b));
            let dw = sub(scale(gw, a), scale(dell_w, b));
            // u = f_v − f_p ; w = f_n − f_v.
            g[iv] = add(g[iv], sub(du, dw));
            g[ip] = sub(g[ip], du);
            g[inx] = add(g[inx], dw);
        }
    }
    (val, g)
}
