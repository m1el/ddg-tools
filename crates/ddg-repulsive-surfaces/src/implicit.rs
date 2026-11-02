//! Implicit surfaces for obstacles and attractors, including the infinite [`Cylinder`] and the
//! framed torus [`FramedTorus`] (`SPEC.md` §5.2, D21, D22).
//!
//! The generic [`ImplicitSurface`] trait and `Sphere`, `Torus` (fixed `y` axis), `Plane` and
//! unions come from `ddg-repulsive-curves` and are re-exported here. `f` is negative inside and
//! `f = 0` on the surface; for an obstacle `|f|` should be a distance (all surfaces here are
//! exact signed distances, away from the axes). Gradients at singular points (a cylinder's axis,
//! a torus' core circle and axis) are the zero vector, never NaN.

pub use ddg_repulsive_curves::constraints::ImplicitSurface;
pub use ddg_repulsive_curves::surfaces::{Plane, Sphere, Torus, Union};

use crate::error::RsError;
use crate::geometry::v3::{add, dot, norm, scale, sub};

fn check_finite(name: &'static str, v: &[f64]) -> Result<(), RsError> {
    match v.iter().find(|x| !x.is_finite()) {
        Some(x) => Err(RsError::InvalidOption {
            name,
            value: x.to_string(),
            requirement: "finite numbers",
        }),
        None => Ok(()),
    }
}

fn check_positive(name: &'static str, x: f64) -> Result<(), RsError> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(RsError::InvalidOption {
            name,
            value: x.to_string(),
            requirement: "a finite positive number",
        })
    }
}

fn unit_axis(name: &'static str, axis: [f64; 3]) -> Result<[f64; 3], RsError> {
    check_finite(name, &axis)?;
    let n = norm(axis);
    if n == 0.0 || !n.is_finite() {
        return Err(RsError::InvalidOption {
            name,
            value: format!("{axis:?}"),
            requirement: "a non-zero axis",
        });
    }
    Ok(scale(axis, 1.0 / n))
}

/// The infinite cylinder `f(x) = |(x − c)_⊥| − r` around the line through `point` with the
/// given axis (`SPEC.md` §5.2; the reference's `cylinder`). The gradient is the unit radial
/// direction, and zero on the axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    point: [f64; 3],
    axis: [f64; 3],
    radius: f64,
}

impl Cylinder {
    /// A cylinder through `point` along `axis` (normalized here) with the given radius.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for non-finite input, a zero axis or a non-positive radius.
    pub fn new(point: [f64; 3], axis: [f64; 3], radius: f64) -> Result<Cylinder, RsError> {
        check_finite("Cylinder.point", &point)?;
        let axis = unit_axis("Cylinder.axis", axis)?;
        check_positive("Cylinder.radius", radius)?;
        Ok(Cylinder {
            point,
            axis,
            radius,
        })
    }

    /// The radial part `(x − c) − ((x − c)·n) n`.
    fn radial(&self, x: [f64; 3]) -> [f64; 3] {
        let p = sub(x, self.point);
        sub(p, scale(self.axis, dot(p, self.axis)))
    }
}

impl ImplicitSurface for Cylinder {
    fn value(&self, x: [f64; 3]) -> f64 {
        norm(self.radial(x)) - self.radius
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let q = self.radial(x);
        let n = norm(q);
        if n > 0.0 { scale(q, 1.0 / n) } else { [0.0; 3] }
    }
}

/// A torus with a frame: centre `c`, unit axis `n`, major radius `R`, minor radius `r`;
/// `f(x) = √((ρ − R)² + h²) − r` with `h = (x − c)·n` and `ρ` the distance from the axis (the
/// signed distance to the torus; the reference's torus is the case `n = e_y`). The gradient
/// is zero on the core circle (`f = −r`) and the radial part is zero on the axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FramedTorus {
    center: [f64; 3],
    axis: [f64; 3],
    major: f64,
    minor: f64,
}

impl FramedTorus {
    /// A torus with the given frame and radii.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for non-finite input, a zero axis or non-positive radii.
    pub fn new(
        center: [f64; 3],
        axis: [f64; 3],
        major: f64,
        minor: f64,
    ) -> Result<FramedTorus, RsError> {
        check_finite("FramedTorus.center", &center)?;
        let axis = unit_axis("FramedTorus.axis", axis)?;
        check_positive("FramedTorus.major", major)?;
        check_positive("FramedTorus.minor", minor)?;
        Ok(FramedTorus {
            center,
            axis,
            major,
            minor,
        })
    }

    /// `(h, radial vector, ρ)` of `x`.
    fn coords(&self, x: [f64; 3]) -> (f64, [f64; 3], f64) {
        let p = sub(x, self.center);
        let h = dot(p, self.axis);
        let q = sub(p, scale(self.axis, h));
        (h, q, norm(q))
    }
}

impl ImplicitSurface for FramedTorus {
    fn value(&self, x: [f64; 3]) -> f64 {
        let (h, _, rho) = self.coords(x);
        (rho - self.major).hypot(h) - self.minor
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let (h, q, rho) = self.coords(x);
        let d = (rho - self.major).hypot(h);
        if d == 0.0 {
            return [0.0; 3];
        }
        let radial = if rho > 0.0 {
            scale(q, (rho - self.major) / (rho * d))
        } else {
            [0.0; 3]
        };
        add(radial, scale(self.axis, h / d))
    }
}
