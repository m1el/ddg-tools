//! Implicit-surface obstacles and attractors, and the Coulomb vertex energy (`SPEC.md` §3.5,
//! §5.2, D21, D22, D24).

use rayon::prelude::*;

use crate::error::RsError;
use crate::geometry::v3::{axpy, dot, scale, sub};
use crate::geometry::{FaceGeometry, corners, face_area_gradient, vertex_areas};
use crate::implicit::ImplicitSurface;
use crate::surface::Surface;

/// The weights `q_v` of the implicit potentials (D21).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Quadrature {
    /// `q_v = a_v`, the barycentric dual area: a discretization of `∫ d^{−p} dA` (default).
    #[default]
    DualArea,
    /// `q_v = 1` (the reference; scene import uses it).
    PerVertex,
}

/// `Σ_v ∂/∂a_v` pulled back to the vertices through `a_v = ⅓ Σ a_S`: `g_v` is the partial
/// derivative with respect to the dual area `a_v`.
fn dual_area_pull_back(
    surface: &Surface,
    positions: &[[f64; 3]],
    g_av: &[f64],
    out: &mut [[f64; 3]],
) {
    for t in surface.faces() {
        let da = (g_av[t[0] as usize] + g_av[t[1] as usize] + g_av[t[2] as usize]) / 3.0;
        let ga = face_area_gradient(corners(positions, t));
        for c in 0..3 {
            out[t[c] as usize] = axpy(out[t[c] as usize], da, ga[c]);
        }
    }
}

fn vertex_weights(
    surface: &Surface,
    positions: &[[f64; 3]],
    quadrature: Quadrature,
) -> Result<Vec<f64>, RsError> {
    Ok(match quadrature {
        Quadrature::PerVertex => vec![1.0; surface.n_vertices()],
        Quadrature::DualArea => {
            let g = FaceGeometry::new(surface, positions)?;
            vertex_areas(surface, &g.areas)
        }
    })
}

/// The obstacle potential `Σ_v q_v d(f_v)^{−power}` and its gradient (including the variation
/// of `q_v` for [`Quadrature::DualArea`]). If any `d(f_v) ≤ 0` the value is `+∞` and the
/// returned gradient is zero (a rejected trial).
///
/// [formula: rs-penalties]
///
/// # Errors
/// The errors of [`Surface::validate_positions`].
pub fn implicit_obstacle(
    surface: &Surface,
    positions: &[[f64; 3]],
    shape: &dyn ImplicitSurface,
    power: f64,
    quadrature: Quadrature,
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    let q = vertex_weights(surface, positions, quadrature)?;
    let n = surface.n_vertices();
    let mut value = 0.0;
    let mut g = vec![[0.0; 3]; n];
    let mut g_av = vec![0.0; n];
    for v in 0..n {
        let d = shape.value(positions[v]);
        if d.is_nan() || d <= 0.0 {
            return Ok((f64::INFINITY, vec![[0.0; 3]; n]));
        }
        let pot = d.powf(-power);
        value += q[v] * pot;
        g_av[v] = pot;
        let c = -power * q[v] * pot / d;
        g[v] = scale(shape.gradient(positions[v]), c);
    }
    if quadrature == Quadrature::DualArea {
        dual_area_pull_back(surface, positions, &g_av, &mut g);
    }
    Ok((value, g))
}

/// The attractor potential `Σ_{v∈set} q_v |d(f_v)|^{power}` (`vertices = None`: all vertices)
/// and its gradient; at `d = 0` the gradient of the term is zero.
///
/// [formula: rs-penalties]
///
/// # Errors
/// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] for a vertex out
/// of range or listed twice.
pub fn implicit_attractor(
    surface: &Surface,
    positions: &[[f64; 3]],
    shape: &dyn ImplicitSurface,
    power: f64,
    quadrature: Quadrature,
    vertices: Option<&[u32]>,
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    let q = vertex_weights(surface, positions, quadrature)?;
    let n = surface.n_vertices();
    let mut active = vec![vertices.is_none(); n];
    if let Some(vs) = vertices {
        for &v in vs {
            if v as usize >= n || active[v as usize] {
                return Err(RsError::InvalidOption {
                    name: "attractor vertices",
                    value: v.to_string(),
                    requirement: "distinct vertex indices below the vertex count",
                });
            }
            active[v as usize] = true;
        }
    }
    let mut value = 0.0;
    let mut g = vec![[0.0; 3]; n];
    let mut g_av = vec![0.0; n];
    for v in (0..n).filter(|&v| active[v]) {
        let d = shape.value(positions[v]);
        let pot = d.abs().powf(power);
        value += q[v] * pot;
        g_av[v] = pot;
        if d != 0.0 {
            let c = power * q[v] * pot / d;
            g[v] = scale(shape.gradient(positions[v]), c);
        }
    }
    if quadrature == Quadrature::DualArea {
        dual_area_pull_back(surface, positions, &g_av, &mut g);
    }
    Ok((value, g))
}

/// `Energy::Coulomb { alpha }` (`SPEC.md` §3.5): `Σ_{v≠w} a_v a_w |f_v − f_w|^{−α}` over ordered
/// vertex pairs with the dual areas `a_v`, dense `O(|V|²)`, and its exact differential
/// (including the variation of `a_v`). Coincident vertices give a non-finite energy; the caller
/// decides. Rows are summed in order, so the result is independent of the thread count.
///
/// [formula: rs-coulomb-energy]
///
/// # Errors
/// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] for a non-finite
/// or non-positive `alpha`.
pub fn coulomb_energy(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    alpha: f64,
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    if !(alpha.is_finite() && alpha > 0.0) {
        return Err(RsError::InvalidOption {
            name: "alpha",
            value: alpha.to_string(),
            requirement: "a finite positive number",
        });
    }
    let geom = FaceGeometry::new(surface, positions)?;
    let a = vertex_areas(surface, &geom.areas);
    let n = surface.n_vertices();
    // Per row v: (Σ_w a_w r^{-α}, ∂/∂f_v of Σ_w a_v a_w r^{-α} (one ordered copy)).
    let rows: Vec<(f64, [f64; 3])> = pool.install(|| {
        (0..n)
            .into_par_iter()
            .map(|v| {
                let mut s = 0.0;
                let mut d = [0.0; 3];
                for w in 0..n {
                    if w == v {
                        continue;
                    }
                    let x = sub(positions[v], positions[w]);
                    let r2 = dot(x, x);
                    let k = a[w] * r2.powf(-0.5 * alpha);
                    s += k;
                    d = axpy(d, -alpha * k / r2, x);
                }
                (s, scale(d, a[v]))
            })
            .collect()
    });
    let energy: f64 = (0..n).map(|v| a[v] * rows[v].0).sum();
    // Ordered pairs: both (v,w) and (w,v) depend on f_v → factor 2.
    let mut g: Vec<[f64; 3]> = rows.iter().map(|r| scale(r.1, 2.0)).collect();
    let g_av: Vec<f64> = rows.iter().map(|r| 2.0 * r.0).collect();
    dual_area_pull_back(surface, positions, &g_av, &mut g);
    Ok((energy, g))
}
