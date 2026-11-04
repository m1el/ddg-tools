//! A small saddle-rich analytic set for the torus sign-rule question (SPEC.md §18: WP1
//! resolution 1, WP5 resolution 6), and the curvature statistics that count saddle points under
//! the exact geometry.
//!
//! A saddle is a surface point with `κmin · κmax < 0` (principal curvatures with respect to the
//! outward normal, convex positive). Revolving a profile whose boundary faces the axis gives
//! saddles on that side (the parallel circles curve away from the normal, the meridian curves
//! towards it, or vice versa); a torus has them on its whole inner half. The set:
//!
//! | name | profile (`Shape2d`) | `o` | saddles |
//! |---|---|---|---|
//! | `torus` | circle `r = 0.25` | 0.5 | inner half; area fraction `1/2 − r/(π o)` exactly |
//! | `thin-torus` | circle `r = 0.12` | 0.6 | same, fraction → 1/2 |
//! | `fat-torus` | circle `r = 0.35` | 0.4 | same |
//! | `tall-ellipse-ring` | ellipse `(0.1, 0.35)` | 0.45 | hyperboloid-like waist on the inner side |
//! | `lens-ring` | vesica `r = 0.45, d = 0.2` | 0.4 | inner arc |
//! | `arc-ring` | arc `ra = 0.3, rb = 0.07, t = 2` | 0.45 | the bent band |
//! | `moon-ring` | moon `d = 0.3, ra = 0.4, rb = 0.3` | 0.5 | concave side of the crescent |
//! | `heart-ring` | heart | 0.65 | the cleft and the inner side |
//! | `washer` | pill `(0.2, 0.1)` | 0.4 | rounded inner rim; the flat faces are developable |
//! | `sphere` | circle `r = 0.6` | 0 | none (control) |
//! | `cylinder` | circle `r = 0.4` extruded, `h = 0.3` | — | none (control, K = 0 everywhere) |

use super::sdf2d::{Kind, Shape2d};
use super::solid::{Solid, SolidKind};
use crate::PatError;
use crate::csg::ImplicitField;
use rayon::prelude::*;

/// A named saddle-rich (or control) shape.
#[derive(Debug, Clone)]
pub struct SaddleShape {
    /// Name (see the table in the module docs).
    pub name: &'static str,
    /// The solid.
    pub solid: Solid,
    /// The exact saddle area fraction when it has a closed form.
    pub analytic_fraction: Option<f64>,
}

fn rev(kind: Kind, p: &[f64], o: f64, name: &'static str) -> Result<SaddleShape, PatError> {
    Ok(SaddleShape { name, solid: Solid::new(Shape2d::new(kind, p)?, SolidKind::Revolution { o })?, analytic_fraction: None })
}

fn torus(r: f64, o: f64, name: &'static str) -> Result<SaddleShape, PatError> {
    let mut s = rev(Kind::Circle, &[r], o, name)?;
    s.analytic_fraction = Some(0.5 - r / (std::f64::consts::PI * o));
    Ok(s)
}

/// The set (11 shapes; the last two are controls with no saddle points).
pub fn saddle_set() -> Result<Vec<SaddleShape>, PatError> {
    Ok(vec![
        torus(0.25, 0.5, "torus")?,
        torus(0.12, 0.6, "thin-torus")?,
        torus(0.35, 0.4, "fat-torus")?,
        rev(Kind::Ellipse, &[0.1, 0.35], 0.45, "tall-ellipse-ring")?,
        rev(Kind::Vesica, &[0.45, 0.2], 0.4, "lens-ring")?,
        rev(Kind::Arc, &[0.3, 0.07, 2.0], 0.45, "arc-ring")?,
        rev(Kind::Moon, &[0.3, 0.4, 0.3], 0.5, "moon-ring")?,
        rev(Kind::Heart, &[], 0.65, "heart-ring")?,
        rev(Kind::Pill, &[0.2, 0.1], 0.4, "washer")?,
        SaddleShape { analytic_fraction: Some(0.0), ..rev(Kind::Circle, &[0.6], 0.0, "sphere")? },
        SaddleShape {
            name: "cylinder",
            solid: Solid::new(Shape2d::new(Kind::Circle, &[0.4])?, SolidKind::Extrusion { h: 0.3 })?,
            analytic_fraction: None,
        },
    ])
}

/// Principal curvatures `(κmin, κmax)` of the level set through `x`, from the Hessian of the
/// field by central differences of its exact gradient (step `h`), projected on the tangent
/// plane; convex positive for `φ < 0` inside. Exact for an exact SDF (`‖∇φ‖ = 1`).
pub fn principal_curvatures(field: &dyn ImplicitField, x: [f64; 3], h: f64) -> Result<(f64, f64), PatError> {
    let (_, g) = field.distance_gradient(x)?;
    let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    if !(gl > 0.0) {
        return Err(PatError::DegenerateQuery { index: 0, reason: "zero gradient".into() });
    }
    let n = [g[0] / gl, g[1] / gl, g[2] / gl];
    let mut hess = [[0.0; 3]; 3];
    for a in 0..3 {
        let (mut xp, mut xm) = (x, x);
        xp[a] += h;
        xm[a] -= h;
        let (_, gp) = field.distance_gradient(xp)?;
        let (_, gm) = field.distance_gradient(xm)?;
        for b in 0..3 {
            hess[a][b] = (gp[b] - gm[b]) / (2.0 * h);
        }
    }
    // Tangent basis.
    let helper = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let t1 = cross(n, helper);
    let l1 = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
    let t1 = [t1[0] / l1, t1[1] / l1, t1[2] / l1];
    let t2 = cross(n, t1);
    let m = |u: [f64; 3], v: [f64; 3]| {
        let mut s = 0.0;
        for a in 0..3 {
            for b in 0..3 {
                s += u[a] * hess[a][b] * v[b];
            }
        }
        s
    };
    let (a, d) = (m(t1, t1), m(t2, t2));
    let b = 0.5 * (m(t1, t2) + m(t2, t1));
    let mean = 0.5 * (a + d);
    let disc = (0.25 * (a - d) * (a - d) + b * b).sqrt();
    Ok((mean - disc, mean + disc))
}

/// Curvature statistics of a surface sample.
#[derive(Debug, Clone, PartialEq)]
pub struct SaddleStats {
    /// Points examined.
    pub n: usize,
    /// Points with `κmin < −eps` and `κmax > eps` (saddle).
    pub saddle: usize,
    /// Points with both curvatures within `±eps` of zero (flat or developable).
    pub flat: usize,
    /// `saddle / n`.
    pub fraction: f64,
    /// Smallest `κmin` and largest `κmax` seen.
    pub extremes: (f64, f64),
}

/// Threshold below which a principal curvature counts as zero (units of 1/length).
pub const SADDLE_EPS: f64 = 1e-4;

/// Count the points of `solid`'s area-uniform surface sample (`n` points, stream `seed`) with
/// `κmin · κmax < 0` (both beyond [`SADDLE_EPS`] in magnitude).
pub fn saddle_stats(solid: &Solid, n: usize, seed: u64) -> Result<SaddleStats, PatError> {
    let (pts, _) = solid.sample_surface(n, seed)?;
    let ks: Result<Vec<(f64, f64)>, PatError> = pts.par_iter().map(|&p| principal_curvatures(solid, p, 1e-5)).collect();
    let ks = ks?;
    let saddle = ks.iter().filter(|k| k.0 < -SADDLE_EPS && k.1 > SADDLE_EPS).count();
    let flat = ks.iter().filter(|k| k.0.abs() <= SADDLE_EPS && k.1.abs() <= SADDLE_EPS).count();
    let lo = ks.iter().map(|k| k.0).fold(f64::INFINITY, f64::min);
    let hi = ks.iter().map(|k| k.1).fold(f64::NEG_INFINITY, f64::max);
    Ok(SaddleStats { n, saddle, flat, fraction: saddle as f64 / n.max(1) as f64, extremes: (lo, hi) })
}
