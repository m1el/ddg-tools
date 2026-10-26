//! Implicit surfaces and vector fields (`SPEC.md` §5.1, §5.3).
//!
//! * [`ImplicitSurface`] (defined with the constraints) for [`Constraint::Surface`](crate::constraints::Constraint):
//!   [`Sphere`], [`Torus`], [`Plane`], [`DoubleTorus`], [`Union`], [`SmoothUnion`] (the set of
//!   the authors' `src/implicit_surface.cpp`, MIT, Copyright (c) 2019 Christopher Yu, ported)
//!   and the additions [`RoundedBox`] and [`Metaballs`]. `f` is negative inside and need not be
//!   a signed distance; the constraint scales its tolerance by `|∇f|`. Gradients at the
//!   singular points (a sphere's centre, the torus axis, the creases of a box) are returned as
//!   zero or as one of the one-sided gradients, never NaN.
//! * [`VectorField`] (value and Jacobian) for the field penalty: [`ConstantField`],
//!   [`CircularField`] and [`SphereTangentField`] (the reference's `interesting` field). The
//!   reference gives zero Jacobians for two of its fields and a non-Jacobian for the circular
//!   one (D20 family); ours are exact derivatives, checked against finite differences.

use std::fmt;
use std::sync::Arc;

pub use crate::constraints::ImplicitSurface;
use crate::error::RcError;
use crate::kernel::v3::{add, cross, dot, norm, scale, sub};

fn check_finite(name: &'static str, values: &[f64]) -> Result<(), RcError> {
    match values.iter().find(|x| !x.is_finite()) {
        Some(x) => Err(RcError::InvalidOption {
            name,
            value: x.to_string(),
            requirement: "finite numbers",
        }),
        None => Ok(()),
    }
}

fn check_positive(name: &'static str, x: f64) -> Result<(), RcError> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(RcError::InvalidOption {
            name,
            value: x.to_string(),
            requirement: "a finite positive number",
        })
    }
}

/// The sphere `|x − c| − r`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    center: [f64; 3],
    radius: f64,
}

impl Sphere {
    /// A sphere (the reference's default: centre 0, radius 1).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite centre or a non-positive radius.
    pub fn new(center: [f64; 3], radius: f64) -> Result<Sphere, RcError> {
        check_finite("Sphere.center", &center)?;
        check_positive("Sphere.radius", radius)?;
        Ok(Sphere { center, radius })
    }
}

impl ImplicitSurface for Sphere {
    fn value(&self, x: [f64; 3]) -> f64 {
        norm(sub(x, self.center)) - self.radius
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let d = sub(x, self.center);
        let n = norm(d);
        if n > 0.0 { scale(d, 1.0 / n) } else { [0.0; 3] }
    }
}

/// The torus around the `y` axis, `√((ρ − R)² + y²) − r` with `ρ = √(x² + z²)` (centred at `c`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    center: [f64; 3],
    major: f64,
    minor: f64,
}

impl Torus {
    /// A torus with major radius `R` and minor radius `r` (reference defaults: 1 and 0.25).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite centre or non-positive radii.
    pub fn new(center: [f64; 3], major: f64, minor: f64) -> Result<Torus, RcError> {
        check_finite("Torus.center", &center)?;
        check_positive("Torus.major", major)?;
        check_positive("Torus.minor", minor)?;
        Ok(Torus {
            center,
            major,
            minor,
        })
    }
}

impl ImplicitSurface for Torus {
    fn value(&self, x: [f64; 3]) -> f64 {
        let p = sub(x, self.center);
        let rho = p[0].hypot(p[2]);
        (rho - self.major).hypot(p[1]) - self.minor
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let p = sub(x, self.center);
        let rho = p[0].hypot(p[2]);
        let q = (rho - self.major).hypot(p[1]);
        if q == 0.0 {
            return [0.0; 3];
        }
        // On the axis (ρ = 0) the radial direction is undefined: take the radial part as 0.
        let k = if rho > 0.0 {
            (rho - self.major) / rho
        } else {
            0.0
        };
        [k * p[0] / q, p[1] / q, k * p[2] / q]
    }
}

/// The plane `n·(x − c)` with unit normal `n`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    point: [f64; 3],
    normal: [f64; 3],
}

impl Plane {
    /// The plane through `point` with normal `normal` (normalized here).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-finite input or a zero normal.
    pub fn new(point: [f64; 3], normal: [f64; 3]) -> Result<Plane, RcError> {
        check_finite("Plane.point", &point)?;
        check_finite("Plane.normal", &normal)?;
        let n = norm(normal);
        if n == 0.0 {
            return Err(RcError::InvalidOption {
                name: "Plane.normal",
                value: format!("{normal:?}"),
                requirement: "a non-zero normal",
            });
        }
        Ok(Plane {
            point,
            normal: scale(normal, 1.0 / n),
        })
    }

    /// The reference's `yplane`: `y = 0`.
    pub fn y_zero() -> Plane {
        Plane {
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
        }
    }
}

impl ImplicitSurface for Plane {
    fn value(&self, x: [f64; 3]) -> f64 {
        dot(self.normal, sub(x, self.point))
    }
    fn gradient(&self, _x: [f64; 3]) -> [f64; 3] {
        self.normal
    }
}

/// The reference's `doubletorus`: with `X = x + 2`, `f = (X(X−1)²(X−2) + y²)² + z² − r²`
/// (a topological double torus; the squared radius `r²` is the parameter).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleTorus {
    radius2: f64,
}

impl DoubleTorus {
    /// The double torus with `f = … − radius2`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-positive or non-finite `radius2`.
    pub fn new(radius2: f64) -> Result<DoubleTorus, RcError> {
        check_positive("DoubleTorus.radius2", radius2)?;
        Ok(DoubleTorus { radius2 })
    }
}

impl ImplicitSurface for DoubleTorus {
    fn value(&self, p: [f64; 3]) -> f64 {
        let x = p[0] + 2.0;
        let inner = x * (x - 1.0) * (x - 1.0) * (x - 2.0) + p[1] * p[1];
        inner * inner + p[2] * p[2] - self.radius2
    }
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let x = p[0] + 2.0;
        let y2 = p[1] * p[1];
        let polyn = x.powi(4) - 4.0 * x.powi(3) + 5.0 * x * x - 2.0 * x + y2;
        [
            2.0 * polyn * (4.0 * x.powi(3) - 12.0 * x * x + 10.0 * x - 2.0),
            4.0 * p[1] * polyn,
            2.0 * p[2],
        ]
    }
}

/// The union of surfaces, `min_i f_i`; the gradient is that of the minimizing surface (the
/// lowest index on ties). Away from the crease between two surfaces it is smooth.
#[derive(Clone)]
pub struct Union {
    parts: Vec<Arc<dyn ImplicitSurface>>,
}

impl fmt::Debug for Union {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Union({} surfaces)", self.parts.len())
    }
}

impl Union {
    /// The union of `parts`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for an empty list.
    pub fn new(parts: Vec<Arc<dyn ImplicitSurface>>) -> Result<Union, RcError> {
        if parts.is_empty() {
            return Err(RcError::InvalidOption {
                name: "Union.parts",
                value: "[]".into(),
                requirement: "at least one surface",
            });
        }
        Ok(Union { parts })
    }

    fn argmin(&self, x: [f64; 3]) -> (usize, f64) {
        let mut best = (0, self.parts[0].value(x));
        for (i, s) in self.parts.iter().enumerate().skip(1) {
            let v = s.value(x);
            if v < best.1 {
                best = (i, v);
            }
        }
        best
    }
}

impl ImplicitSurface for Union {
    fn value(&self, x: [f64; 3]) -> f64 {
        self.argmin(x).1
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        self.parts[self.argmin(x).0].gradient(x)
    }
}

/// The smooth union of two surfaces (polynomial smooth minimum of the reference):
/// `h = clamp(½ + ½ (f₂ − f₁)/k, 0, 1)`, `f = (1 − h) f₂ + h f₁ − k h (1 − h)`.
#[derive(Clone)]
pub struct SmoothUnion {
    a: Arc<dyn ImplicitSurface>,
    b: Arc<dyn ImplicitSurface>,
    k: f64,
}

impl fmt::Debug for SmoothUnion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SmoothUnion(k = {})", self.k)
    }
}

impl SmoothUnion {
    /// The smooth union with blend width `k`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for a non-positive `k`.
    pub fn new(
        a: Arc<dyn ImplicitSurface>,
        b: Arc<dyn ImplicitSurface>,
        k: f64,
    ) -> Result<SmoothUnion, RcError> {
        check_positive("SmoothUnion.k", k)?;
        Ok(SmoothUnion { a, b, k })
    }

    fn h(&self, d1: f64, d2: f64) -> f64 {
        (0.5 + 0.5 * (d2 - d1) / self.k).clamp(0.0, 1.0)
    }
}

impl ImplicitSurface for SmoothUnion {
    fn value(&self, x: [f64; 3]) -> f64 {
        let (d1, d2) = (self.a.value(x), self.b.value(x));
        let h = self.h(d1, d2);
        (1.0 - h) * d2 + h * d1 - self.k * h * (1.0 - h)
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let (d1, d2) = (self.a.value(x), self.b.value(x));
        let (g1, g2) = (self.a.gradient(x), self.b.gradient(x));
        let h = self.h(d1, d2);
        let raw = 0.5 + 0.5 * (d2 - d1) / self.k;
        // d clamp = 0 outside (0, 1)
        let dh = if raw > 0.0 && raw < 1.0 {
            scale(sub(g2, g1), 0.5 / self.k)
        } else {
            [0.0; 3]
        };
        // (1-h) g2 + h g1 + (d1 - d2 - k(1 - 2h)) dh
        let mut g = add(scale(g2, 1.0 - h), scale(g1, h));
        g = add(g, scale(dh, d1 - d2 - self.k * (1.0 - 2.0 * h)));
        g
    }
}

/// A rounded box: the signed distance of the box with half extents `b` rounded by `r`
/// (`r = 0`: a plain box), centred at `c`, axis aligned. Not smooth at the creases of the
/// underlying box when `r = 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoundedBox {
    center: [f64; 3],
    half: [f64; 3],
    radius: f64,
}

impl RoundedBox {
    /// A box with half extents `half` (each `> radius`), rounding radius `radius >= 0`.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for non-finite input, a negative radius, or half extents
    /// that are not larger than the radius.
    pub fn new(center: [f64; 3], half: [f64; 3], radius: f64) -> Result<RoundedBox, RcError> {
        check_finite("RoundedBox.center", &center)?;
        check_finite("RoundedBox.half", &half)?;
        check_finite("RoundedBox.radius", &[radius])?;
        if radius < 0.0 || half.iter().any(|&h| h <= radius) {
            return Err(RcError::InvalidOption {
                name: "RoundedBox.half",
                value: format!("{half:?} with radius {radius}"),
                requirement: "radius >= 0 and every half extent larger than the radius",
            });
        }
        Ok(RoundedBox {
            center,
            half,
            radius,
        })
    }

    /// A plain box.
    ///
    /// # Errors
    /// As [`RoundedBox::new`].
    pub fn cuboid(center: [f64; 3], half: [f64; 3]) -> Result<RoundedBox, RcError> {
        RoundedBox::new(center, half, 0.0)
    }

    fn q(&self, x: [f64; 3]) -> [f64; 3] {
        let p = sub(x, self.center);
        let mut q = [0.0; 3];
        for k in 0..3 {
            q[k] = p[k].abs() - (self.half[k] - self.radius);
        }
        q
    }
}

impl ImplicitSurface for RoundedBox {
    fn value(&self, x: [f64; 3]) -> f64 {
        let q = self.q(x);
        let out = norm([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]);
        out + q[0].max(q[1]).max(q[2]).min(0.0) - self.radius
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let p = sub(x, self.center);
        let q = self.q(x);
        let sign = |t: f64| if t < 0.0 { -1.0 } else { 1.0 };
        let o = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
        let n = norm(o);
        if n > 0.0 {
            return [
                sign(p[0]) * o[0] / n,
                sign(p[1]) * o[1] / n,
                sign(p[2]) * o[2] / n,
            ];
        }
        // inside the core box: the axis of the largest q (lowest index on ties)
        let k = if q[0] >= q[1] && q[0] >= q[2] {
            0
        } else if q[1] >= q[2] {
            1
        } else {
            2
        };
        let mut g = [0.0; 3];
        g[k] = sign(p[k]);
        g
    }
}

/// Metaballs: `f = T − Σ_i a_i exp(−|x − c_i|²/b_i²)`, negative where the sum of Gaussians
/// exceeds the threshold `T`.
#[derive(Debug, Clone, PartialEq)]
pub struct Metaballs {
    balls: Vec<([f64; 3], f64, f64)>,
    threshold: f64,
}

impl Metaballs {
    /// Metaballs from `(centre, amplitude a, width b)` triples.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for an empty list, non-finite numbers, or non-positive
    /// amplitudes, widths or threshold.
    pub fn new(balls: Vec<([f64; 3], f64, f64)>, threshold: f64) -> Result<Metaballs, RcError> {
        if balls.is_empty() {
            return Err(RcError::InvalidOption {
                name: "Metaballs.balls",
                value: "[]".into(),
                requirement: "at least one ball",
            });
        }
        check_positive("Metaballs.threshold", threshold)?;
        for &(c, a, b) in &balls {
            check_finite("Metaballs.center", &c)?;
            check_positive("Metaballs.amplitude", a)?;
            check_positive("Metaballs.width", b)?;
        }
        Ok(Metaballs { balls, threshold })
    }
}

impl ImplicitSurface for Metaballs {
    fn value(&self, x: [f64; 3]) -> f64 {
        let mut s = 0.0;
        for &(c, a, b) in &self.balls {
            let d = sub(x, c);
            s += a * (-dot(d, d) / (b * b)).exp();
        }
        self.threshold - s
    }
    fn gradient(&self, x: [f64; 3]) -> [f64; 3] {
        let mut g = [0.0; 3];
        for &(c, a, b) in &self.balls {
            let d = sub(x, c);
            let e = a * (-dot(d, d) / (b * b)).exp();
            g = add(g, scale(d, 2.0 * e / (b * b)));
        }
        g
    }
}

// =============================================================================================
// Vector fields
// =============================================================================================

/// A vector field `X` on `R³` for the field penalty (`SPEC.md` §5.3). `X` is not normalized by
/// the penalty.
pub trait VectorField: Send + Sync {
    /// `X(x)`.
    fn value(&self, x: [f64; 3]) -> [f64; 3];
    /// The Jacobian `J[i][j] = ∂X_i/∂x_j`.
    fn jacobian(&self, x: [f64; 3]) -> [[f64; 3]; 3];
}

/// The constant field `X(x) = v`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantField {
    v: [f64; 3],
}

impl ConstantField {
    /// The field `v`, used as given.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for non-finite components.
    pub fn new(v: [f64; 3]) -> Result<ConstantField, RcError> {
        check_finite("ConstantField.v", &v)?;
        Ok(ConstantField { v })
    }

    /// The unit field `v/|v|` (the reference normalizes its constant field).
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] for non-finite or zero `v`.
    pub fn unit(v: [f64; 3]) -> Result<ConstantField, RcError> {
        check_finite("ConstantField.v", &v)?;
        let n = norm(v);
        if n == 0.0 {
            return Err(RcError::InvalidOption {
                name: "ConstantField.v",
                value: format!("{v:?}"),
                requirement: "a non-zero vector",
            });
        }
        ConstantField::new(scale(v, 1.0 / n))
    }
}

impl VectorField for ConstantField {
    fn value(&self, _x: [f64; 3]) -> [f64; 3] {
        self.v
    }
    fn jacobian(&self, _x: [f64; 3]) -> [[f64; 3]; 3] {
        [[0.0; 3]; 3]
    }
}

/// The unit field circulating around the `y` axis, `(z, 0, −x)/√(x² + z²)` (the reference's
/// `circular`: zero on the axis, where the Jacobian is also set to zero).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CircularField;

impl VectorField for CircularField {
    fn value(&self, x: [f64; 3]) -> [f64; 3] {
        let rho = x[0].hypot(x[2]);
        if rho == 0.0 {
            return [0.0; 3];
        }
        [x[2] / rho, 0.0, -x[0] / rho]
    }
    fn jacobian(&self, x: [f64; 3]) -> [[f64; 3]; 3] {
        let rho = x[0].hypot(x[2]);
        if rho == 0.0 {
            return [[0.0; 3]; 3];
        }
        let r3 = rho * rho * rho;
        let (a, c) = (x[0], x[2]);
        [
            [-c * a / r3, 0.0, a * a / r3],
            [0.0; 3],
            [-c * c / r3, 0.0, a * c / r3],
        ]
    }
}

/// The reference's `interesting` field: on the unit sphere, the gradient of a sum of two
/// Gaussians projected to the tangent plane, rotated by 90 degrees about the normal `N = x/|x|`
/// and normalized. Defined away from the origin and from the zeros of the Gaussian gradient.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SphereTangentField;

const GAUSS: [(f64, f64, [f64; 3]); 2] = [
    (
        0.33570404546960386,
        0.45645123825263556,
        [0.6127842849933675, -0.7513894319788624, 0.24476384858808273],
    ),
    (
        0.6339815359591325,
        0.39270985340263576,
        [-0.5468629255159502, -0.8331733676192284, 0.082237948577103],
    ),
];

type Mat = [[f64; 3]; 3];

fn mat_mul(a: &Mat, b: &Mat) -> Mat {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    r
}

fn mat_sub(a: &Mat, b: &Mat) -> Mat {
    let mut r = *a;
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] -= b[i][j];
        }
    }
    r
}

fn mat_scale(a: &Mat, s: f64) -> Mat {
    let mut r = *a;
    for row in &mut r {
        for x in row {
            *x *= s;
        }
    }
    r
}

fn outer(a: [f64; 3], b: [f64; 3]) -> Mat {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i] * b[j];
        }
    }
    r
}

fn identity() -> Mat {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// The matrix `[v]×` with `[v]× u = v × u`.
fn skew(v: [f64; 3]) -> Mat {
    [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
}

fn transpose_vec(m: &Mat, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
        m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
        m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
    ]
}

impl SphereTangentField {
    /// `(X, J_X)` of the Gaussian sum.
    fn gaussians(x: [f64; 3]) -> ([f64; 3], Mat) {
        let mut v = [0.0; 3];
        let mut j = [[0.0; 3]; 3];
        for &(a, b, c) in &GAUSS {
            let d = sub(x, c);
            let b2 = b * b;
            let e = (-dot(d, d) / b2).exp();
            let k = -2.0 * a / b2 * e;
            v = add(v, scale(d, k));
            // d/dx [k d] = k (I - 2 d dᵀ / b²)
            let m = mat_sub(&identity(), &mat_scale(&outer(d, d), 2.0 / b2));
            for r in 0..3 {
                for s in 0..3 {
                    j[r][s] += k * m[r][s];
                }
            }
        }
        (v, j)
    }
}

impl VectorField for SphereTangentField {
    fn value(&self, x: [f64; 3]) -> [f64; 3] {
        self.eval(x).0
    }
    fn jacobian(&self, x: [f64; 3]) -> [[f64; 3]; 3] {
        self.eval(x).1
    }
}

impl SphereTangentField {
    fn eval(&self, x: [f64; 3]) -> ([f64; 3], Mat) {
        let r = norm(x);
        if r == 0.0 {
            return ([0.0; 3], [[0.0; 3]; 3]);
        }
        let n = scale(x, 1.0 / r);
        let jn = mat_scale(&mat_sub(&identity(), &outer(n, n)), 1.0 / r);
        let (xv, jx) = SphereTangentField::gaussians(x);
        // Y = X - (X.N) N
        let xn = dot(xv, n);
        let y = sub(xv, scale(n, xn));
        // dY = J_X - N (Nᵀ J_X + Xᵀ J_N) - (X.N) J_N
        let row = {
            let a = transpose_vec(&jx, n);
            let b = transpose_vec(&jn, xv);
            add(a, b)
        };
        let jy = mat_sub(&mat_sub(&jx, &outer(n, row)), &mat_scale(&jn, xn));
        // Z = Y x N; J_Z = -[N]x J_Y + [Y]x J_N
        let jz = {
            let a = mat_mul(&skew(n), &jy);
            let b = mat_mul(&skew(y), &jn);
            // b - a = -[N] J_Y + [Y] J_N
            mat_sub(&b, &a)
        };
        let z = cross(y, n);
        let zn = norm(z);
        if zn == 0.0 {
            return ([0.0; 3], [[0.0; 3]; 3]);
        }
        let w = scale(z, 1.0 / zn);
        let proj = mat_scale(&mat_sub(&identity(), &outer(w, w)), 1.0 / zn);
        (w, mat_mul(&proj, &jz))
    }
}

/// `J ᵀ v` for a Jacobian `J[i][j] = ∂X_i/∂x_j`.
pub fn jacobian_transpose_apply(j: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    transpose_vec(j, v)
}
