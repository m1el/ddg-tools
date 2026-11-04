//! Revolutions and extrusions of 2D signed distance fields (Quilez; SPEC.md §11) with exact
//! gradients, and area-uniform surface sampling with exact normals.
//!
//! * Revolution about the y axis, offset `o`: `φ(x, y, z) = φ₂(‖(x, z)‖ − o, y)`
//!   `[formula: pat-revolution]`. Exact distance whenever the profile `{φ₂ = 0}` stays in
//!   `s ≥ −o`; the part of the profile with `s < −o` is not on the surface (the solid is the
//!   revolved half-plane `ρ ≥ 0` slice). On the axis the gradient is taken in the x direction.
//! * Extrusion by half-height `h`: with `w = (φ₂(x, y), |z| − h)`,
//!   `φ = min(max(w₀, w₁), 0) + ‖max(w, 0)‖` `[formula: pat-extrusion]`; gradient
//!   `(a/l)(∇φ₂, 0) + (b/l)(0, 0, sgn z)` for `a, b = max(w, 0)`, `l = ‖(a, b)‖ > 0`, and the
//!   gradient of the larger `wᵢ` inside `[formula: pat-extrusion-gradient]`.
//!
//! Surface samples are area-uniform by construction: extrusion sides are (boundary arclength) ×
//! (uniform z), the caps are uniform in the 2D region (rejection sampling in the profile's box);
//! revolution picks a boundary segment by `length × mean radius` and a point in it with density
//! proportional to the radius (exact inverse CDF of the linear radius), then a uniform angle.
//! Boundary vertices are projected exactly onto `φ₂ = 0`, so `|φ| < 1e-12` on all samples.

use super::boundary::{Boundary2d, project};
use super::sdf2d::Shape2d;
use crate::csg::ImplicitField;
use crate::PatError;
use crate::generate::rng::Rng;
use std::f64::consts::TAU;
use std::sync::OnceLock;

/// How the profile becomes a solid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolidKind {
    /// Revolution about the y axis; the profile is shifted by `o` in the radial direction.
    Revolution {
        /// Radial offset `o`.
        o: f64,
    },
    /// Extrusion along z to `[−h, h]`.
    Extrusion {
        /// Half height `h`.
        h: f64,
    },
}

/// A 3D shape from a 2D profile.
#[derive(Debug)]
pub struct Solid {
    profile: Shape2d,
    kind: SolidKind,
    boundary: OnceLock<Boundary2d>,
}

impl Clone for Solid {
    fn clone(&self) -> Self {
        Solid { profile: self.profile.clone(), kind: self.kind, boundary: OnceLock::new() }
    }
}

/// Grid cells (along the longest side) of the profile boundary trace.
pub const TRACE_CELLS: usize = 2048;

impl Solid {
    /// A revolution or extrusion (`o` may have either sign for a revolution; `h > 0`).
    pub fn new(profile: Shape2d, kind: SolidKind) -> Result<Solid, PatError> {
        match kind {
            SolidKind::Revolution { o } if !o.is_finite() => Err(PatError::InvalidOption { what: format!("revolution offset {o}") }),
            SolidKind::Extrusion { h } if !(h.is_finite() && h > 0.0) => Err(PatError::InvalidOption { what: format!("extrusion half height {h}") }),
            _ => Ok(Solid { profile, kind, boundary: OnceLock::new() }),
        }
    }

    /// The profile.
    pub fn profile(&self) -> &Shape2d {
        &self.profile
    }

    /// Revolution or extrusion.
    pub fn kind(&self) -> SolidKind {
        self.kind
    }

    fn eval2(&self) -> impl Fn([f64; 2]) -> (f64, [f64; 2]) + Sync + '_ {
        move |p| self.profile.eval(p)
    }

    /// The traced profile boundary (computed once).
    pub fn boundary(&self) -> &Boundary2d {
        self.boundary.get_or_init(|| Boundary2d::trace(&self.eval2(), self.profile.bbox(), TRACE_CELLS))
    }

    /// `φ` and `∇φ` (no error path: all inputs are finite or the result is non-finite).
    pub fn eval(&self, x: [f64; 3]) -> (f64, [f64; 3]) {
        match self.kind {
            SolidKind::Revolution { o } => {
                let rho = x[0].hypot(x[2]);
                let (d, g) = self.profile.eval([rho - o, x[1]]);
                if rho > 1e-300 {
                    (d, [g[0] * x[0] / rho, g[1], g[0] * x[2] / rho])
                } else {
                    (d, [g[0], g[1], 0.0])
                }
            }
            SolidKind::Extrusion { h } => {
                let (d, g) = self.profile.eval([x[0], x[1]]);
                let w1 = x[2].abs() - h;
                let sz = if x[2] >= 0.0 { 1.0 } else { -1.0 };
                let (a, b) = (d.max(0.0), w1.max(0.0));
                let l = a.hypot(b);
                let phi = d.max(w1).min(0.0) + l;
                let grad = if l > 0.0 {
                    [a / l * g[0], a / l * g[1], b / l * sz]
                } else if d > w1 {
                    [g[0], g[1], 0.0]
                } else {
                    [0.0, 0.0, sz]
                };
                (phi, grad)
            }
        }
    }

    /// Total surface area (from the traced boundary, relative accuracy about 1e-5).
    pub fn surface_area(&self) -> f64 {
        let b = self.boundary();
        match self.kind {
            SolidKind::Extrusion { h } => 2.0 * h * b.length + 2.0 * b.area,
            SolidKind::Revolution { o } => TAU * self.rev_segments(o).iter().map(|s| s.weight).sum::<f64>(),
        }
    }

    /// `n` area-uniform surface points with their exact unit normals, from the stream `seed`.
    pub fn sample_surface(&self, n: usize, seed: u64) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), PatError> {
        let mut rng = Rng::new(seed);
        self.sample_surface_with(n, &mut rng)
    }

    /// As [`Solid::sample_surface`] with the caller's generator.
    pub fn sample_surface_with(&self, n: usize, rng: &mut Rng) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), PatError> {
        let b = self.boundary();
        let ev = self.eval2();
        let mut pts = Vec::with_capacity(n);
        let mut nrm = Vec::with_capacity(n);
        let degenerate = || PatError::Format { reason: format!("the {} profile has an empty boundary", self.profile.kind().name()) };
        match self.kind {
            SolidKind::Extrusion { h } => {
                let side = 2.0 * h * b.length;
                let total = side + 2.0 * b.area;
                let bb = self.profile.bbox();
                for _ in 0..n {
                    if rng.f64() * total < side {
                        let (p, g) = b.sample(&ev, rng.f64()).ok_or_else(degenerate)?;
                        pts.push([p[0], p[1], rng.uniform(-h, h)]);
                        nrm.push([g[0], g[1], 0.0]);
                    } else {
                        let sz = if rng.f64() < 0.5 { 1.0 } else { -1.0 };
                        let mut found = None;
                        for _ in 0..100_000 {
                            let q = [rng.uniform(bb[0][0], bb[1][0]), rng.uniform(bb[0][1], bb[1][1])];
                            if self.profile.distance(q) < 0.0 {
                                found = Some(q);
                                break;
                            }
                        }
                        let q = found.ok_or_else(degenerate)?;
                        pts.push([q[0], q[1], sz * h]);
                        nrm.push([0.0, 0.0, sz]);
                    }
                }
            }
            SolidKind::Revolution { o } => {
                let segs = self.rev_segments(o);
                let total: f64 = segs.iter().map(|s| s.weight).sum();
                if !(total > 0.0) {
                    return Err(degenerate());
                }
                let mut cum = Vec::with_capacity(segs.len());
                let mut c = 0.0;
                for s in &segs {
                    c += s.weight;
                    cum.push(c);
                }
                for _ in 0..n {
                    let target = rng.f64() * total;
                    let i = cum.partition_point(|&c| c <= target).min(segs.len() - 1);
                    let s = &segs[i];
                    let (ra, rb) = (s.a[0] + o, s.b[0] + o);
                    let u = rng.f64();
                    let t = if (rb - ra).abs() < 1e-12 * (ra + rb).abs().max(1e-300) {
                        u
                    } else {
                        ((ra * ra + u * (rb * rb - ra * ra)).max(0.0).sqrt() - ra) / (rb - ra)
                    };
                    let p = project(&ev, [s.a[0] + t * (s.b[0] - s.a[0]), s.a[1] + t * (s.b[1] - s.a[1])], 3);
                    let (_, g) = ev(p);
                    let l = g[0].hypot(g[1]);
                    let (g, rho) = ([g[0] / l, g[1] / l], (p[0] + o).max(0.0));
                    let psi = rng.uniform(0.0, TAU);
                    let (sn, cs) = psi.sin_cos();
                    pts.push([rho * cs, p[1], rho * sn]);
                    nrm.push([g[0] * cs, g[1], g[0] * sn]);
                }
            }
        }
        Ok((pts, nrm))
    }

    /// Boundary segments clipped to the radial half-plane `s ≥ −o`, with `length × mean radius`.
    fn rev_segments(&self, o: f64) -> Vec<RevSeg> {
        let mut out = Vec::new();
        for &(a, b) in &self.boundary().segments {
            let (ra, rb) = (a[0] + o, b[0] + o);
            let (a, b) = if ra >= 0.0 && rb >= 0.0 {
                (a, b)
            } else if ra < 0.0 && rb < 0.0 {
                continue;
            } else {
                let t = ra / (ra - rb); // radius crosses 0
                let m = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])];
                if ra < 0.0 { (m, b) } else { (a, m) }
            };
            let len = (b[0] - a[0]).hypot(b[1] - a[1]);
            let weight = len * 0.5 * (a[0] + o + b[0] + o);
            if weight > 0.0 {
                out.push(RevSeg { a, b, weight });
            }
        }
        out
    }
}

struct RevSeg {
    a: [f64; 2],
    b: [f64; 2],
    weight: f64,
}

impl ImplicitField for Solid {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        self.distance_gradient(x).map(|(d, _)| d)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        if let Some(a) = x.iter().position(|c| !c.is_finite()) {
            return Err(PatError::NonFinite { index: a, what: "query" });
        }
        Ok(self.eval(x))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        let b = self.profile.bbox();
        let ext = (b[1][0] - b[0][0]).max(b[1][1] - b[0][1]);
        let pad = 0.05 * ext;
        match self.kind {
            SolidKind::Revolution { o } => {
                let r = (b[1][0] + o).max(b[0][0] + o).max(1e-9) + pad;
                [[-r, b[0][1] - pad, -r], [r, b[1][1] + pad, r]]
            }
            SolidKind::Extrusion { h } => [[b[0][0] - pad, b[0][1] - pad, -h - pad], [b[1][0] + pad, b[1][1] + pad, h + pad]],
        }
    }
}
