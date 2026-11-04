//! The 16 two-dimensional signed distance functions of Inigo Quilez (<https://iquilezles.org/articles/distfunctions2d/>),
//! with exact gradients. Port of the shapes of the reference's `shape_2d.{h,cpp}` (MIT, Copyright
//! (c) 2026 Nicole Feng), with these differences: the polygons (box, cross, triangles, hexagon,
//! pentagon, trapezoid, quad) share one exact point-to-polygon routine whose gradient is the edge
//! normal on an edge and `s (p − v)/|p − v|` at a vertex, so it is also exact *on* the boundary
//! (the reference's `q / |q|` is 0/0 there); the ellipse is an exact closest-point solve (Eberly's
//! robust bisection) instead of four Newton steps; the heart's straight edges are handled in
//! closed form. `[formula: pat-quilez-2d]` (exactness cited; Quilez 2019).
//!
//! Sign: negative inside. Gradient: unit outward normal on the boundary, `s (p − q)/|p − q|`
//! elsewhere (`q` the closest boundary point); at a point where it is not defined (medial axis
//! points, vertices of the shape itself) some unit vector of the one-sided limits is returned.

use crate::PatError;
use crate::generate::rng::Rng;
use std::f64::consts::PI;

/// The 16 classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
    /// Disc of radius `r`.
    Circle,
    /// Pie / pac-man wedge: radius `r`, half-angle `t`.
    Pie,
    /// Thick circular arc: radii `ra`, `rb` (thickness), half-angle.
    Arc,
    /// Capsule (the paper's "segment").
    Pill,
    /// Lens of two circles (vesica piscis); for `d < 0` the formula yields the union of the two discs.
    Vesica,
    /// Axis-aligned rectangle.
    Box,
    /// Plus sign.
    Cross,
    /// Isosceles triangle.
    Isosceles,
    /// Triangle.
    Triangle,
    /// Regular hexagon.
    Hexagon,
    /// Ellipse.
    Ellipse,
    /// Heart.
    Heart,
    /// Regular pentagon.
    Pentagon,
    /// Crescent.
    Moon,
    /// Trapezoid.
    Trapezoid,
    /// Quadrilateral.
    Quad,
}

/// All 16 kinds in the order of the table of NOTES-REFERENCE §12 (circle … quad).
pub const ALL_KINDS: [Kind; 16] = [
    Kind::Circle,
    Kind::Pie,
    Kind::Arc,
    Kind::Pill,
    Kind::Vesica,
    Kind::Box,
    Kind::Cross,
    Kind::Isosceles,
    Kind::Triangle,
    Kind::Hexagon,
    Kind::Ellipse,
    Kind::Heart,
    Kind::Pentagon,
    Kind::Moon,
    Kind::Trapezoid,
    Kind::Quad,
];

/// The 14 kinds of App. E (no heart, no isosceles triangle).
pub const PAPER_KINDS: [Kind; 14] = [
    Kind::Circle,
    Kind::Pie,
    Kind::Arc,
    Kind::Pill,
    Kind::Vesica,
    Kind::Box,
    Kind::Cross,
    Kind::Triangle,
    Kind::Hexagon,
    Kind::Ellipse,
    Kind::Pentagon,
    Kind::Moon,
    Kind::Trapezoid,
    Kind::Quad,
];

impl Kind {
    /// Lower-case name (as in the reference's `sample_SDF`).
    pub fn name(self) -> &'static str {
        match self {
            Kind::Circle => "circle",
            Kind::Pie => "pie",
            Kind::Arc => "arc",
            Kind::Pill => "pill",
            Kind::Vesica => "vesica",
            Kind::Box => "box",
            Kind::Cross => "cross",
            Kind::Isosceles => "isosceles",
            Kind::Triangle => "triangle",
            Kind::Hexagon => "hexagon",
            Kind::Ellipse => "ellipse",
            Kind::Heart => "heart",
            Kind::Pentagon => "pentagon",
            Kind::Moon => "moon",
            Kind::Trapezoid => "trapezoid",
            Kind::Quad => "quad",
        }
    }

    /// The kind with this name.
    pub fn from_name(name: &str) -> Option<Kind> {
        ALL_KINDS.iter().copied().find(|k| k.name() == name)
    }

    /// Number of parameters of [`Shape2d::new`].
    pub fn param_count(self) -> usize {
        match self {
            Kind::Circle | Kind::Hexagon | Kind::Pentagon => 1,
            Kind::Pie | Kind::Pill | Kind::Vesica | Kind::Box | Kind::Cross | Kind::Isosceles | Kind::Ellipse => 2,
            Kind::Arc | Kind::Moon | Kind::Trapezoid => 3,
            Kind::Heart => 0,
            Kind::Triangle => 6,
            Kind::Quad => 8,
        }
    }
}

#[derive(Debug, Clone)]
enum Geo {
    Circle { r: f64 },
    Pie { r: f64, sc: [f64; 2] },
    Arc { sc: [f64; 2], ra: f64, rb: f64 },
    Pill { a: [f64; 2], b: [f64; 2], r: f64 },
    Vesica { r: f64, d: f64, b: f64 },
    Ellipse { a: f64, b: f64 },
    Heart,
    Moon { d: f64, ra: f64, rb: f64 },
    /// CCW polygon.
    Poly { v: Vec<[f64; 2]> },
}

/// A 2D shape: its class, its parameters and the derived geometry.
#[derive(Debug, Clone)]
pub struct Shape2d {
    kind: Kind,
    params: Vec<f64>,
    geo: Geo,
}

fn sgn(x: f64) -> f64 {
    if x >= 0.0 { 1.0 } else { -1.0 }
}

fn norm2(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

fn bad(what: String) -> PatError {
    PatError::InvalidOption { what }
}

fn ccw(mut v: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    let mut a = 0.0;
    for i in 0..v.len() {
        let (p, q) = (v[i], v[(i + 1) % v.len()]);
        a += p[0] * q[1] - q[0] * p[1];
    }
    if a < 0.0 {
        v.reverse();
    }
    v
}

impl Shape2d {
    /// Build a shape. Parameters per kind: circle `[r]`; pie `[r, t]` (sector of half-angle `t`
    /// about +y); arc `[ra, rb, t]`; pill `[half_length, r]` (axis along x); vesica `[r, d]`
    /// (`|d| < r`); box `[bx, by]` (half extents); cross `[bx, by]` (arm half length / half
    /// width, either order); isosceles `[qx, qy]` (apex at the origin, base at `y = qy`);
    /// triangle `[v0, v1, v2]` flattened; hexagon `[r]` (apothem); ellipse `[a, b]`; heart `[]`;
    /// pentagon `[r]` (apothem); moon `[d, ra, rb]` (`|ra − rb| < d < ra + rb`); trapezoid
    /// `[ra, rb, he]` (half widths at `y = −he` and `y = +he`); quad `[v0, v1, v2, v3]` flattened.
    pub fn new(kind: Kind, p: &[f64]) -> Result<Shape2d, PatError> {
        if p.len() != kind.param_count() {
            return Err(bad(format!("{} takes {} parameters, got {}", kind.name(), kind.param_count(), p.len())));
        }
        if let Some(i) = p.iter().position(|x| !x.is_finite()) {
            return Err(bad(format!("{} parameter {i} is not finite", kind.name())));
        }
        let pos = |i: usize| -> Result<f64, PatError> {
            if p[i] > 0.0 { Ok(p[i]) } else { Err(bad(format!("{} parameter {i} = {} must be positive", kind.name(), p[i]))) }
        };
        let geo = match kind {
            Kind::Circle => Geo::Circle { r: pos(0)? },
            Kind::Pie => Geo::Pie { r: pos(0)?, sc: [p[1].sin(), p[1].cos()] },
            Kind::Arc => Geo::Arc { sc: [p[2].sin(), p[2].cos()], ra: pos(0)?, rb: pos(1)? },
            Kind::Pill => Geo::Pill { a: [-pos(0)?, 0.0], b: [p[0], 0.0], r: pos(1)? },
            Kind::Vesica => {
                let (r, d) = (pos(0)?, p[1]);
                if d.abs() >= r {
                    return Err(bad(format!("vesica needs |d| < r, got r = {r}, d = {d}")));
                }
                Geo::Vesica { r, d, b: (r * r - d * d).sqrt() }
            }
            Kind::Box => {
                let (bx, by) = (pos(0)?, pos(1)?);
                Geo::Poly { v: vec![[-bx, -by], [bx, -by], [bx, by], [-bx, by]] }
            }
            Kind::Cross => {
                let (a, b) = (pos(0)?, pos(1)?);
                let (l, t) = (a.max(b), a.min(b));
                Geo::Poly {
                    v: vec![
                        [l, t],
                        [t, t],
                        [t, l],
                        [-t, l],
                        [-t, t],
                        [-l, t],
                        [-l, -t],
                        [-t, -t],
                        [-t, -l],
                        [t, -l],
                        [t, -t],
                        [l, -t],
                    ],
                }
            }
            Kind::Isosceles => {
                let (qx, qy) = (pos(0)?, p[1]);
                if qy == 0.0 {
                    return Err(bad("isosceles needs qy ≠ 0".into()));
                }
                Geo::Poly { v: ccw(vec![[0.0, 0.0], [qx, qy], [-qx, qy]]) }
            }
            Kind::Triangle => {
                Geo::Poly { v: ccw(vec![[p[0], p[1]], [p[2], p[3]], [p[4], p[5]]]) }
            }
            Kind::Hexagon => {
                let r = pos(0)?;
                let k = r / 3f64.sqrt();
                Geo::Poly { v: vec![[2.0 * k, 0.0], [k, r], [-k, r], [-2.0 * k, 0.0], [-k, -r], [k, -r]] }
            }
            Kind::Ellipse => {
                let (a, b) = (pos(0)?, pos(1)?);
                Geo::Ellipse { a, b }
            }
            Kind::Heart => Geo::Heart,
            Kind::Pentagon => {
                let r = pos(0)?;
                let big = r / (PI / 5.0).cos();
                Geo::Poly {
                    v: (0..5)
                        .map(|k| {
                            let a = PI / 2.0 + 2.0 * PI * k as f64 / 5.0;
                            [big * a.cos(), big * a.sin()]
                        })
                        .collect(),
                }
            }
            Kind::Moon => {
                let (d, ra, rb) = (pos(0)?, pos(1)?, pos(2)?);
                if !((ra - rb).abs() < d && d < ra + rb) {
                    return Err(bad(format!("moon needs |ra − rb| < d < ra + rb, got d = {d}, ra = {ra}, rb = {rb}")));
                }
                Geo::Moon { d, ra, rb }
            }
            Kind::Trapezoid => {
                let (ra, rb, he) = (pos(0)?, pos(1)?, pos(2)?);
                Geo::Poly { v: vec![[-ra, -he], [ra, -he], [rb, he], [-rb, he]] }
            }
            Kind::Quad => Geo::Poly { v: ccw(vec![[p[0], p[1]], [p[2], p[3]], [p[4], p[5]], [p[6], p[7]]]) },
        };
        Ok(Shape2d { kind, params: p.to_vec(), geo })
    }

    /// Draw the parameters as App. E / the reference's `sample_SDF` (`r ∈ [0.1, 0.9]`,
    /// `t ∈ [π/6, 5π/6]`, radius ratios in `[0.1, 0.9]`), with independent draws for every
    /// parameter and perturbation directions uniform on the circle (SPEC.md §11; the reference
    /// draws them in the first quadrant and reuses keys). Invalid draws are redrawn.
    pub fn random(kind: Kind, rng: &mut Rng) -> Shape2d {
        let (rmin, rmax, ratio_min, ratio_max) = (0.1, 0.9, 0.1, 0.9);
        let (tmin, tmax) = (PI / 6.0, 5.0 * PI / 6.0);
        for _ in 0..1000 {
            let params: Vec<f64> = match kind {
                Kind::Circle | Kind::Hexagon | Kind::Pentagon => vec![rng.uniform(rmin, rmax)],
                Kind::Pie => vec![rng.uniform(rmin, rmax), rng.uniform(tmin, tmax)],
                Kind::Arc => {
                    let ra = rng.uniform(rmin, rmax);
                    let rb = rng.uniform(ratio_min, ratio_max) * ra;
                    vec![ra, rb, rng.uniform(tmin, tmax)]
                }
                Kind::Pill => {
                    let hl = rng.uniform(rmin, rmax);
                    vec![hl, rng.uniform(ratio_min, ratio_max) * hl]
                }
                Kind::Vesica => {
                    let r = rng.uniform(rmin, rmax);
                    vec![r, rng.uniform(ratio_min, ratio_max) * 2.0 * r - r]
                }
                Kind::Box | Kind::Cross | Kind::Ellipse => vec![rng.uniform(rmin, rmax), rng.uniform(rmin, rmax)],
                Kind::Isosceles => vec![rng.uniform(rmin, rmax), -rng.uniform(rmin, rmax)],
                Kind::Triangle => {
                    let h = rng.uniform(rmin, rmax);
                    let b = 2.0 * h / 3f64.sqrt();
                    let base = [[-b, -h], [b, -h], [0.0, h]];
                    let mut v = Vec::new();
                    for c in base {
                        let d = rng.unit_circle();
                        v.push(c[0] + ratio_max * h * d[0]);
                        v.push(c[1] + ratio_max * h * d[1]);
                    }
                    v
                }
                Kind::Heart => vec![],
                Kind::Moon => {
                    let ra = rng.uniform(rmin, rmax);
                    let rb = rng.uniform(ratio_min, ratio_max) * ra;
                    let d = rng.uniform((ra - rb) / ratio_max, (ra + rb) * ratio_max);
                    vec![d, ra, rb]
                }
                Kind::Trapezoid => vec![rng.uniform(rmin, rmax), rng.uniform(rmin, rmax), 2.0 * rng.uniform(rmin, rmax)],
                Kind::Quad => {
                    let (w, h) = (rng.uniform(rmin, rmax), rng.uniform(rmin, rmax));
                    let amt = (ratio_max * h).min(ratio_max * w);
                    let base = [[-w, -h], [w, -h], [w, h], [-w, h]];
                    let mut v = Vec::new();
                    for c in base {
                        let d = rng.unit_circle();
                        v.push(c[0] + amt * d[0]);
                        v.push(c[1] + amt * d[1]);
                    }
                    v
                }
            };
            if let Ok(s) = Shape2d::new(kind, &params) {
                if kind == Kind::Quad && !s.is_simple_polygon() {
                    continue;
                }
                return s;
            }
        }
        // Unreachable in practice (every kind has a valid set of draws); a fixed valid shape.
        Shape2d::new(Kind::Circle, &[0.5]).unwrap_or_else(|_| unreachable!())
    }

    /// The class.
    pub fn kind(&self) -> Kind {
        self.kind
    }

    /// The parameters as given to [`Shape2d::new`].
    pub fn params(&self) -> &[f64] {
        &self.params
    }

    /// Vertices (counter-clockwise) of the polygonal kinds.
    pub fn polygon(&self) -> Option<&[[f64; 2]]> {
        match &self.geo {
            Geo::Poly { v } => Some(v),
            _ => None,
        }
    }

    /// A polygon with no two non-adjacent edges crossing (true for non-polygons).
    pub fn is_simple_polygon(&self) -> bool {
        let Geo::Poly { v } = &self.geo else { return true };
        let n = v.len();
        let cross = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
        for i in 0..n {
            for j in i + 2..n {
                if i == 0 && j == n - 1 {
                    continue;
                }
                let (a, b, c, d) = (v[i], v[(i + 1) % n], v[j], v[(j + 1) % n]);
                if cross(a, b, c) * cross(a, b, d) < 0.0 && cross(c, d, a) * cross(c, d, b) < 0.0 {
                    return false;
                }
            }
        }
        true
    }

    /// An axis-aligned box `[min, max]` containing the shape (tight for most kinds, loose for
    /// the arc and the moon).
    pub fn bbox(&self) -> [[f64; 2]; 2] {
        match &self.geo {
            Geo::Circle { r } => [[-r, -r], [*r, *r]],
            Geo::Pie { r, sc } => {
                let t = sc[0].atan2(sc[1]);
                let xm = if t >= PI / 2.0 { *r } else { r * sc[0] };
                let y0 = if t > PI / 2.0 { r * sc[1] } else { 0.0 };
                [[-xm, y0.min(0.0)], [xm, *r]]
            }
            Geo::Arc { sc, ra, rb } => {
                let t = sc[0].atan2(sc[1]);
                let xm = if t >= PI / 2.0 { ra + rb } else { ra * sc[0] + rb };
                let y0 = if t < PI { ra * sc[1] - rb } else { -ra - rb };
                [[-xm, y0], [xm, ra + rb]]
            }
            Geo::Pill { a, b, r } => [[a[0].min(b[0]) - r, a[1].min(b[1]) - r], [a[0].max(b[0]) + r, a[1].max(b[1]) + r]],
            // d > 0: lens of the two discs (height 2b); d < 0 the formula gives their union (height 2r).
            Geo::Vesica { r, d, b } => {
                let h = if *d < 0.0 { *r } else { *b };
                [[-(r - d), -h], [r - d, h]]
            }
            Geo::Ellipse { a, b } => [[-a, -b], [*a, *b]],
            Geo::Heart => {
                let r = 2f64.sqrt() / 4.0;
                [[-0.25 - r, 0.0], [0.25 + r, 0.75 + r]]
            }
            Geo::Moon { ra, .. } => [[-ra, -ra], [*ra, *ra]],
            Geo::Poly { v } => {
                let mut lo = [f64::INFINITY; 2];
                let mut hi = [f64::NEG_INFINITY; 2];
                for p in v {
                    for a in 0..2 {
                        lo[a] = lo[a].min(p[a]);
                        hi[a] = hi[a].max(p[a]);
                    }
                }
                [lo, hi]
            }
        }
    }

    /// `φ(p)` and `∇φ(p)`.
    pub fn eval(&self, p: [f64; 2]) -> (f64, [f64; 2]) {
        match &self.geo {
            Geo::Circle { r } => {
                let l = norm2(p);
                (l - r, if l > 0.0 { [p[0] / l, p[1] / l] } else { [1.0, 0.0] })
            }
            Geo::Pie { r, sc } => eval_pie(p, *r, *sc),
            Geo::Arc { sc, ra, rb } => eval_arc(p, *sc, *ra, *rb),
            Geo::Pill { a, b, r } => {
                let ba = [b[0] - a[0], b[1] - a[1]];
                let pa = [p[0] - a[0], p[1] - a[1]];
                let h = ((pa[0] * ba[0] + pa[1] * ba[1]) / (ba[0] * ba[0] + ba[1] * ba[1])).clamp(0.0, 1.0);
                let q = [pa[0] - h * ba[0], pa[1] - h * ba[1]];
                let d = norm2(q);
                let g = if d > 0.0 {
                    [q[0] / d, q[1] / d]
                } else {
                    let l = norm2(ba);
                    [-ba[1] / l, ba[0] / l]
                };
                (d - r, g)
            }
            Geo::Vesica { r, d, b } => {
                let s = [sgn(p[0]), sgn(p[1])];
                let q = [p[0].abs(), p[1].abs()];
                if (q[1] - b) * d > q[0] * b {
                    let w = [q[0], q[1] - b];
                    let l = norm2(w);
                    let sd = sgn(*d);
                    if l > 0.0 {
                        (l * sd, [sd * s[0] * w[0] / l, sd * s[1] * w[1] / l])
                    } else {
                        (0.0, [0.0, s[1]])
                    }
                } else {
                    let w = [q[0] + d, q[1]];
                    let l = norm2(w);
                    (l - r, [s[0] * w[0] / l, s[1] * w[1] / l])
                }
            }
            Geo::Ellipse { a, b } => eval_ellipse(p, *a, *b),
            Geo::Heart => eval_heart(p),
            Geo::Moon { d, ra, rb } => eval_moon(p, *d, *ra, *rb),
            Geo::Poly { v } => eval_polygon(p, v),
        }
    }

    /// `φ(p)`.
    pub fn distance(&self, p: [f64; 2]) -> f64 {
        self.eval(p).0
    }
}

fn eval_pie(q: [f64; 2], r: f64, sc: [f64; 2]) -> (f64, [f64; 2]) {
    let s = sgn(q[0]);
    let p = [q[0].abs(), q[1]];
    let l = norm2(p);
    let n = l - r;
    let t = (p[0] * sc[0] + p[1] * sc[1]).clamp(0.0, r);
    let qv = [p[0] - sc[0] * t, p[1] - sc[1] * t];
    let qn = norm2(qv);
    let side = sgn(sc[1] * p[0] - sc[0] * p[1]);
    let m = qn * side;
    let edge_n = [sc[1], -sc[0]];
    let (phi, mut g) = if n > m {
        (n, if l > 0.0 { [p[0] / l, p[1] / l] } else { [0.0, 1.0] })
    } else if t > 1e-12 && t < r - 1e-12 {
        (m, edge_n)
    } else if qn > 0.0 {
        (m, [side * qv[0] / qn, side * qv[1] / qn])
    } else {
        (m, edge_n)
    };
    g[0] *= s;
    (phi, g)
}

fn eval_arc(q: [f64; 2], sc: [f64; 2], ra: f64, rb: f64) -> (f64, [f64; 2]) {
    let s = sgn(q[0]);
    let p = [q[0].abs(), q[1]];
    if sc[1] * p[0] > sc[0] * p[1] {
        let w = [p[0] - ra * sc[0], p[1] - ra * sc[1]];
        let d = norm2(w);
        (d - rb, if d > 0.0 { [s * w[0] / d, w[1] / d] } else { [0.0, 1.0] })
    } else {
        let l = norm2(q);
        let w = l - ra;
        let sw = sgn(w);
        (w.abs() - rb, if l > 0.0 { [sw * q[0] / l, sw * q[1] / l] } else { [0.0, 1.0] })
    }
}

fn eval_moon(q: [f64; 2], d: f64, ra: f64, rb: f64) -> (f64, [f64; 2]) {
    let s = sgn(q[1]);
    let p = [q[0], q[1].abs()];
    let a = (ra * ra - rb * rb + d * d) / (2.0 * d);
    let b = (ra * ra - a * a).max(0.0).sqrt();
    if d * (p[0] * b - p[1] * a) > d * d * (b - p[1]).max(0.0) {
        let w = [p[0] - a, p[1] - b];
        let l = norm2(w);
        (l, if l > 0.0 { [w[0] / l, s * w[1] / l] } else { [1.0, 0.0] })
    } else {
        let w1 = p;
        let l1 = norm2(w1);
        let d1 = l1 - ra;
        let w2 = [p[0] - d, p[1]];
        let l2 = norm2(w2);
        let d2 = rb - l2;
        if d1 > d2 {
            (d1, if l1 > 0.0 { [w1[0] / l1, s * w1[1] / l1] } else { [1.0, 0.0] })
        } else {
            (d2, if l2 > 0.0 { [-w2[0] / l2, -s * w2[1] / l2] } else { [-1.0, 0.0] })
        }
    }
}

fn eval_heart(q: [f64; 2]) -> (f64, [f64; 2]) {
    let sx = sgn(q[0]);
    let p = [q[0].abs(), q[1]];
    let fold = |g: [f64; 2]| [sx * g[0], g[1]];
    let r = 2f64.sqrt() / 4.0;
    if p[1] + p[0] > 1.0 {
        let w = [p[0] - 0.25, p[1] - 0.75];
        let l = norm2(w);
        return (l - r, if l > 0.0 { fold([w[0] / l, w[1] / l]) } else { [0.0, 1.0] });
    }
    let w1 = [p[0], p[1] - 1.0];
    let d1 = norm2(w1);
    // Distance to the line y = x clipped below at the origin (Quilez' `0.5 max(x + y, 0) (1, 1)`).
    let h = (p[0] + p[1]).max(0.0);
    let w2 = [p[0] - 0.5 * h, p[1] - 0.5 * h];
    let d2 = norm2(w2);
    if d1 < d2 {
        let sg = sgn(p[0] - p[1]);
        (sg * d1, if d1 > 0.0 { fold([sg * w1[0] / d1, sg * w1[1] / d1]) } else { [0.0, 1.0] })
    } else if p[0] + p[1] > 0.0 {
        let sg = sgn(p[0] - p[1]);
        let k = std::f64::consts::FRAC_1_SQRT_2;
        (sg * d2, fold([k, -k]))
    } else {
        let sg = sgn(p[0] - p[1]);
        let l = norm2(p);
        (sg * l, if l > 0.0 { fold([sg * p[0] / l, sg * p[1] / l]) } else { [0.0, -1.0] })
    }
}

/// Exact distance to an ellipse with semi-axes `(a, b)` (Eberly, "Distance from a point to an
/// ellipse", robust bisection) and the unit outward normal at the closest point.
fn eval_ellipse(p: [f64; 2], a: f64, b: f64) -> (f64, [f64; 2]) {
    // Work with e0 >= e1 in the first quadrant.
    let swap = a < b;
    let (e0, e1) = if swap { (b, a) } else { (a, b) };
    let (px, py) = if swap { (p[1], p[0]) } else { (p[0], p[1]) };
    let (sx, sy) = (sgn(px), sgn(py));
    let (x, y) = (px.abs(), py.abs());
    let (cx, cy);
    if y > 0.0 {
        if x > 0.0 {
            let (z0, z1) = (x / e0, y / e1);
            let g = z0 * z0 + z1 * z1 - 1.0;
            let r0 = (e0 / e1) * (e0 / e1);
            let f = |s: f64| (r0 * z0 / (s + r0)).powi(2) + (z1 / (s + 1.0)).powi(2) - 1.0;
            let mut s0 = z1 - 1.0;
            let mut s1 = if g < 0.0 { 0.0 } else { (r0 * z0).hypot(z1) - 1.0 };
            let mut s = 0.0;
            for _ in 0..200 {
                s = 0.5 * (s0 + s1);
                if s == s0 || s == s1 {
                    break;
                }
                let gs = f(s);
                if gs > 0.0 {
                    s0 = s;
                } else if gs < 0.0 {
                    s1 = s;
                } else {
                    break;
                }
            }
            cx = r0 * x / (s + r0);
            cy = y / (s + 1.0);
        } else {
            cx = 0.0;
            cy = e1;
        }
    } else {
        let numer0 = e0 * x;
        let denom0 = e0 * e0 - e1 * e1;
        if numer0 < denom0 {
            let xde0 = numer0 / denom0;
            cx = e0 * xde0;
            cy = e1 * (1.0 - xde0 * xde0).max(0.0).sqrt();
        } else {
            cx = e0;
            cy = 0.0;
        }
    }
    let inside = (px / e0).powi(2) + (py / e1).powi(2) < 1.0;
    let dist = (x - cx).hypot(y - cy);
    let n = [cx / (e0 * e0), cy / (e1 * e1)];
    let nl = norm2(n);
    let n = [sx * n[0] / nl, sy * n[1] / nl];
    let n = if swap { [n[1], n[0]] } else { n };
    (if inside { -dist } else { dist }, n)
}

/// Exact signed distance to a counter-clockwise polygon (even-odd inside test).
fn eval_polygon(p: [f64; 2], v: &[[f64; 2]]) -> (f64, [f64; 2]) {
    let n = v.len();
    let mut best = (f64::INFINITY, [0.0; 2], usize::MAX, 0.0);
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (v[i], v[(i + 1) % n]);
        let e = [b[0] - a[0], b[1] - a[1]];
        let ee = e[0] * e[0] + e[1] * e[1];
        if ee < 1e-28 {
            continue;
        }
        let w = [p[0] - a[0], p[1] - a[1]];
        let t = ((w[0] * e[0] + w[1] * e[1]) / ee).clamp(0.0, 1.0);
        let q = [a[0] + t * e[0], a[1] + t * e[1]];
        let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2);
        if d2 < best.0 {
            best = (d2, q, i, t);
        }
        if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]) {
            inside = !inside;
        }
    }
    let (d2, q, i, t) = best;
    if i == usize::MAX {
        return (f64::INFINITY, [1.0, 0.0]);
    }
    let d = d2.sqrt();
    let s = if inside { -1.0 } else { 1.0 };
    let edge_normal = |i: usize| {
        let (a, b) = (v[i], v[(i + 1) % n]);
        let e = [b[0] - a[0], b[1] - a[1]];
        let l = norm2(e);
        [e[1] / l, -e[0] / l]
    };
    let eps = 1e-12;
    let g = if t > eps && t < 1.0 - eps {
        edge_normal(i)
    } else if d > 1e-13 {
        [s * (p[0] - q[0]) / d, s * (p[1] - q[1]) / d]
    } else {
        // On a vertex: the mean of the outward normals of its two edges.
        let j = if t <= eps { (i + n - 1) % n } else { (i + 1) % n };
        let (n1, n2) = (edge_normal(i), edge_normal(j));
        let m = [n1[0] + n2[0], n1[1] + n2[1]];
        let l = norm2(m);
        if l > 1e-12 { [m[0] / l, m[1] / l] } else { n1 }
    };
    (s * d, g)
}
