//! `ImplicitField`, Booleans and offsets (SPEC.md §8.1, §8.2). Owned by WP6.
//!
//! [`ImplicitField`] is the small interface every application takes (marching cubes, level-set
//! sampling, morphology): a signed distance (negative inside), its gradient, and a box that
//! contains the surface. [`TorusField`] implements it; so does [`Csg`], a tree of Booleans and
//! offsets over any fields (the baselines of WP5 can be leaves).
//!
//! * Union `min(φ_A, φ_B)`, intersection `max(φ_A, φ_B)`, difference `max(φ_A, −φ_B)`
//!   `[formula: pat-boolean]`. Gradients are taken from the active argument; at ties from `A`.
//!   Min/max of exact SDFs is exact outside the shape and a bound inside (cited).
//! * Offset by `c`: `φ − c` `[formula: pat-offset]` (no refit; exact for an exact SDF).
//! * Optional smooth variants (ours, off unless you build them): the polynomial smooth minimum of
//!   Quilez with radius `k`, `smin(a, b) = min(a, b) − h² k / 4`, `h = max(k − |a − b|, 0) / k`;
//!   `smax(a, b) = −smin(−a, −b)`. Smooth unions only ever make `φ` smaller than `min`.

use crate::{PatError, TorusField};
use rayon::prelude::*;
use std::sync::Arc;

/// Queries per parallel task of the default batch implementations.
const BATCH_CHUNK: usize = 256;

/// A scalar field `φ` on `R³` with `φ < 0` inside, `φ = 0` on the surface, and `‖∇φ‖ ≈ 1`.
///
/// Single queries return `Result` (a field may reject a query, e.g. `DegenerateQuery`). The batch
/// methods default to the single-point methods run in parallel on the *current* rayon pool (call
/// them inside `pool.install`; marching cubes does); fields with a faster batch path override them.
/// Results must not depend on the batch composition or the number of threads.
pub trait ImplicitField: Sync {
    /// `φ(x)`.
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError>;
    /// `φ(x)` and `∇φ(x)`.
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError>;
    /// An axis-aligned box `[min, max]` that contains the surface (padded; not tight).
    fn bbox(&self) -> [[f64; 3]; 2];
    /// Number of data points behind the field, if any (the default size of resampling).
    fn point_count_hint(&self) -> Option<usize> {
        None
    }
    /// `out[i] = φ(xs[i])`.
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        check_len(xs.len(), out.len(), "distance")?;
        let res: Result<Vec<()>, PatError> = xs
            .par_chunks(BATCH_CHUNK)
            .zip(out.par_chunks_mut(BATCH_CHUNK))
            .map(|(x, o)| {
                for (xi, oi) in x.iter().zip(o.iter_mut()) {
                    *oi = self.distance(*xi)?;
                }
                Ok(())
            })
            .collect();
        res.map(|_| ())
    }
    /// `d[i], g[i] = φ(xs[i]), ∇φ(xs[i])`.
    fn distance_gradients(&self, xs: &[[f64; 3]], d: &mut [f64], g: &mut [[f64; 3]]) -> Result<(), PatError> {
        check_len(xs.len(), d.len(), "distance")?;
        check_len(xs.len(), g.len(), "gradient")?;
        let res: Result<Vec<()>, PatError> = xs
            .par_chunks(BATCH_CHUNK)
            .zip(d.par_chunks_mut(BATCH_CHUNK).zip(g.par_chunks_mut(BATCH_CHUNK)))
            .map(|(x, (dd, gg))| {
                for ((xi, di), gi) in x.iter().zip(dd.iter_mut()).zip(gg.iter_mut()) {
                    let (a, b) = self.distance_gradient(*xi)?;
                    *di = a;
                    *gi = b;
                }
                Ok(())
            })
            .collect();
        res.map(|_| ())
    }
}

fn check_len(queries: usize, buffer: usize, name: &str) -> Result<(), PatError> {
    if queries != buffer {
        return Err(PatError::SizeMismatch { what: format!("{name} buffer has {buffer} entries for {queries} queries") });
    }
    Ok(())
}

impl<T: ImplicitField + ?Sized> ImplicitField for &T {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        (**self).distance(x)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        (**self).distance_gradient(x)
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        (**self).bbox()
    }
    fn point_count_hint(&self) -> Option<usize> {
        (**self).point_count_hint()
    }
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        (**self).distances(xs, out)
    }
    fn distance_gradients(&self, xs: &[[f64; 3]], d: &mut [f64], g: &mut [[f64; 3]]) -> Result<(), PatError> {
        (**self).distance_gradients(xs, d, g)
    }
}

impl ImplicitField for TorusField {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        TorusField::distance(self, x)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        TorusField::distance_gradient(self, x)
    }
    /// The cloud's box padded by 10 % of its largest extent and by `|isovalue|`: the surface of
    /// the blend stays within a fraction of the local spacing of the cloud, and the offset
    /// surface `φ = isovalue` within `|isovalue|` of it.
    fn bbox(&self) -> [[f64; 3]; 2] {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in self.positions() {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        let ext = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max);
        let pad = 0.1 * ext + self.options().isovalue.abs();
        for a in 0..3 {
            lo[a] -= pad;
            hi[a] += pad;
        }
        [lo, hi]
    }
    fn point_count_hint(&self) -> Option<usize> {
        Some(self.positions().len())
    }
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        self.evaluate(xs, crate::EvalOut { distance: Some(out), gradient: None, laplacian: None }).map(|_| ())
    }
    fn distance_gradients(&self, xs: &[[f64; 3]], d: &mut [f64], g: &mut [[f64; 3]]) -> Result<(), PatError> {
        self.evaluate(xs, crate::EvalOut { distance: Some(d), gradient: Some(g), laplacian: None }).map(|_| ())
    }
}

/// A tree of Booleans and offsets. Leaves are any [`ImplicitField`] (shared with `Arc`).
#[derive(Clone)]
pub enum Csg {
    /// A leaf field.
    Leaf(Arc<dyn ImplicitField + Send + Sync>),
    /// `min(φ_A, φ_B)`.
    Union(Box<Csg>, Box<Csg>),
    /// `max(φ_A, φ_B)`.
    Intersect(Box<Csg>, Box<Csg>),
    /// `max(φ_A, −φ_B)`.
    Subtract(Box<Csg>, Box<Csg>),
    /// `φ − c` (offset outward for `c > 0`).
    Offset(Box<Csg>, f64),
    /// Smooth minimum with radius `k > 0` (ours; off unless built).
    SmoothUnion(Box<Csg>, Box<Csg>, f64),
    /// Smooth maximum with radius `k > 0`.
    SmoothIntersect(Box<Csg>, Box<Csg>, f64),
    /// Smooth `max(φ_A, −φ_B)` with radius `k > 0`.
    SmoothSubtract(Box<Csg>, Box<Csg>, f64),
}

impl std::fmt::Debug for Csg {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Csg::Leaf(_) => write!(f, "Leaf"),
            Csg::Union(a, b) => write!(f, "Union({a:?}, {b:?})"),
            Csg::Intersect(a, b) => write!(f, "Intersect({a:?}, {b:?})"),
            Csg::Subtract(a, b) => write!(f, "Subtract({a:?}, {b:?})"),
            Csg::Offset(a, c) => write!(f, "Offset({a:?}, {c})"),
            Csg::SmoothUnion(a, b, k) => write!(f, "SmoothUnion({a:?}, {b:?}, {k})"),
            Csg::SmoothIntersect(a, b, k) => write!(f, "SmoothIntersect({a:?}, {b:?}, {k})"),
            Csg::SmoothSubtract(a, b, k) => write!(f, "SmoothSubtract({a:?}, {b:?}, {k})"),
        }
    }
}

/// Which argument a Boolean takes its value and gradient from, or the smooth blend weights.
#[derive(Clone, Copy)]
struct Pick {
    /// Value.
    v: f64,
    /// Weight of `∇φ_A`.
    wa: f64,
    /// Weight of `∇φ_B` (`σ_b` already folded in for subtraction).
    wb: f64,
}

/// `min(a, b)` with ties to `a`.
fn pick_min(a: f64, b: f64) -> Pick {
    if a <= b { Pick { v: a, wa: 1.0, wb: 0.0 } } else { Pick { v: b, wa: 0.0, wb: 1.0 } }
}

/// `max(a, b)` with ties to `a`.
fn pick_max(a: f64, b: f64) -> Pick {
    if a >= b { Pick { v: a, wa: 1.0, wb: 0.0 } } else { Pick { v: b, wa: 0.0, wb: 1.0 } }
}

/// Polynomial smooth minimum and the weights of the two gradients.
fn smooth_min(a: f64, b: f64, k: f64) -> Pick {
    let h = (k - (a - b).abs()).max(0.0) / k;
    let v = a.min(b) - h * h * k * 0.25;
    // For a <= b: d/da = 1 - h/2, d/db = h/2 (and the mirror image otherwise); at a = b both are 1/2.
    if a <= b { Pick { v, wa: 1.0 - 0.5 * h, wb: 0.5 * h } } else { Pick { v, wa: 0.5 * h, wb: 1.0 - 0.5 * h } }
}

fn smooth_max(a: f64, b: f64, k: f64) -> Pick {
    let p = smooth_min(-a, -b, k);
    Pick { v: -p.v, wa: p.wa, wb: p.wb }
}

impl Csg {
    /// A leaf from any field.
    pub fn leaf<F: ImplicitField + Send + Sync + 'static>(f: F) -> Csg {
        Csg::Leaf(Arc::new(f))
    }
    /// A leaf from a [`TorusField`].
    pub fn field(f: TorusField) -> Csg {
        Csg::leaf(f)
    }
    /// `min(self, other)`.
    pub fn union(self, other: Csg) -> Csg {
        Csg::Union(Box::new(self), Box::new(other))
    }
    /// `max(self, other)`.
    pub fn intersect(self, other: Csg) -> Csg {
        Csg::Intersect(Box::new(self), Box::new(other))
    }
    /// `max(self, −other)`.
    pub fn subtract(self, other: Csg) -> Csg {
        Csg::Subtract(Box::new(self), Box::new(other))
    }
    /// `self − c`.
    pub fn offset(self, c: f64) -> Csg {
        Csg::Offset(Box::new(self), c)
    }

    /// Check the parameters of the whole tree (finite offsets, smooth radii `> 0`).
    pub fn validate(&self) -> Result<(), PatError> {
        match self {
            Csg::Leaf(_) => Ok(()),
            Csg::Union(a, b) | Csg::Intersect(a, b) | Csg::Subtract(a, b) => {
                a.validate()?;
                b.validate()
            }
            Csg::Offset(a, c) => {
                if !c.is_finite() {
                    return Err(PatError::InvalidOption { what: format!("Csg::Offset c = {c}") });
                }
                a.validate()
            }
            Csg::SmoothUnion(a, b, k) | Csg::SmoothIntersect(a, b, k) | Csg::SmoothSubtract(a, b, k) => {
                if !(k.is_finite() && *k > 0.0) {
                    return Err(PatError::InvalidOption { what: format!("smooth Boolean radius k = {k} (must be finite and > 0)") });
                }
                a.validate()?;
                b.validate()
            }
        }
    }

    /// The combination rule of a binary node: `(pick, sign of B's gradient)`.
    fn combine(&self, a: f64, b: f64) -> Option<(Pick, f64)> {
        match self {
            Csg::Union(..) => Some((pick_min(a, b), 1.0)),
            Csg::Intersect(..) => Some((pick_max(a, b), 1.0)),
            Csg::Subtract(..) => Some((pick_max(a, -b), -1.0)),
            Csg::SmoothUnion(_, _, k) => Some((smooth_min(a, b, *k), 1.0)),
            Csg::SmoothIntersect(_, _, k) => Some((smooth_max(a, b, *k), 1.0)),
            Csg::SmoothSubtract(_, _, k) => {
                let p = smooth_max(a, -b, *k);
                Some((p, -1.0))
            }
            _ => None,
        }
    }

    fn children(&self) -> Option<(&Csg, &Csg)> {
        match self {
            Csg::Union(a, b)
            | Csg::Intersect(a, b)
            | Csg::Subtract(a, b)
            | Csg::SmoothUnion(a, b, _)
            | Csg::SmoothIntersect(a, b, _)
            | Csg::SmoothSubtract(a, b, _) => Some((a, b)),
            _ => None,
        }
    }

    /// `[formula: pat-boolean]` `[formula: pat-offset]` value and gradient of the tree at `x`.
    fn eval(&self, x: [f64; 3], want_grad: bool) -> Result<(f64, [f64; 3]), PatError> {
        match self {
            Csg::Leaf(f) => {
                if want_grad {
                    f.distance_gradient(x)
                } else {
                    Ok((f.distance(x)?, [0.0; 3]))
                }
            }
            Csg::Offset(a, c) => {
                let (v, g) = a.eval(x, want_grad)?;
                Ok((v - c, g))
            }
            _ => {
                let (a, b) = self.children().expect("binary node");
                let (va, ga) = a.eval(x, want_grad)?;
                let (vb, gb) = b.eval(x, want_grad)?;
                let (p, sb) = self.combine(va, vb).expect("binary node");
                let g = [0, 1, 2].map(|i| p.wa * ga[i] + sb * p.wb * gb[i]);
                Ok((p.v, g))
            }
        }
    }

    /// Batch evaluation that delegates to the children's batch methods (so a tree of
    /// `TorusField` leaves runs on the leaves' pools).
    fn eval_batch(&self, xs: &[[f64; 3]], d: &mut [f64], g: Option<&mut [[f64; 3]]>) -> Result<(), PatError> {
        match self {
            Csg::Leaf(f) => match g {
                Some(g) => f.distance_gradients(xs, d, g),
                None => f.distances(xs, d),
            },
            Csg::Offset(a, c) => {
                a.eval_batch(xs, d, g)?;
                d.iter_mut().for_each(|v| *v -= c);
                Ok(())
            }
            _ => {
                let (a, b) = self.children().expect("binary node");
                let n = xs.len();
                let mut db = vec![0.0; n];
                match g {
                    Some(g) => {
                        let mut gb = vec![[0.0; 3]; n];
                        a.eval_batch(xs, d, Some(&mut *g))?;
                        b.eval_batch(xs, &mut db, Some(&mut gb))?;
                        for i in 0..n {
                            let (p, sb) = self.combine(d[i], db[i]).expect("binary node");
                            d[i] = p.v;
                            g[i] = [0, 1, 2].map(|c| p.wa * g[i][c] + sb * p.wb * gb[i][c]);
                        }
                    }
                    None => {
                        a.eval_batch(xs, d, None)?;
                        b.eval_batch(xs, &mut db, None)?;
                        for i in 0..n {
                            d[i] = self.combine(d[i], db[i]).expect("binary node").0.v;
                        }
                    }
                }
                Ok(())
            }
        }
    }
}

impl ImplicitField for Csg {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        self.validate()?;
        Ok(self.eval(x, false)?.0)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        self.validate()?;
        self.eval(x, true)
    }
    /// Union: the union of the boxes; intersection: their overlap (a degenerate box if empty);
    /// difference: `A`'s box; offset by `c > 0`: the box grown by `c`; smooth blends: grown by `k/4`
    /// (the largest amount `smin` lowers `min`).
    fn bbox(&self) -> [[f64; 3]; 2] {
        match self {
            Csg::Leaf(f) => f.bbox(),
            Csg::Offset(a, c) => grow(a.bbox(), c.max(0.0)),
            Csg::Union(a, b) => merge(a.bbox(), b.bbox()),
            Csg::SmoothUnion(a, b, k) => grow(merge(a.bbox(), b.bbox()), 0.25 * k),
            Csg::Intersect(a, b) => overlap(a.bbox(), b.bbox()),
            Csg::SmoothIntersect(a, b, k) => grow(overlap(a.bbox(), b.bbox()), 0.25 * k),
            Csg::Subtract(a, _) => a.bbox(),
            Csg::SmoothSubtract(a, _, k) => grow(a.bbox(), 0.25 * k),
        }
    }
    fn point_count_hint(&self) -> Option<usize> {
        match self {
            Csg::Leaf(f) => f.point_count_hint(),
            _ => None,
        }
    }
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        self.validate()?;
        check_len(xs.len(), out.len(), "distance")?;
        self.eval_batch(xs, out, None)
    }
    fn distance_gradients(&self, xs: &[[f64; 3]], d: &mut [f64], g: &mut [[f64; 3]]) -> Result<(), PatError> {
        self.validate()?;
        check_len(xs.len(), d.len(), "distance")?;
        check_len(xs.len(), g.len(), "gradient")?;
        self.eval_batch(xs, d, Some(g))
    }
}

fn grow(b: [[f64; 3]; 2], r: f64) -> [[f64; 3]; 2] {
    [b[0].map(|v| v - r), b[1].map(|v| v + r)]
}

fn merge(a: [[f64; 3]; 2], b: [[f64; 3]; 2]) -> [[f64; 3]; 2] {
    [[0, 1, 2].map(|i| a[0][i].min(b[0][i])), [0, 1, 2].map(|i| a[1][i].max(b[1][i]))]
}

fn overlap(a: [[f64; 3]; 2], b: [[f64; 3]; 2]) -> [[f64; 3]; 2] {
    let lo = [0, 1, 2].map(|i| a[0][i].max(b[0][i]));
    let hi = [0, 1, 2].map(|i| a[1][i].min(b[1][i]));
    // an empty overlap collapses to a point at the midpoint of the gap
    let lo2 = [0, 1, 2].map(|i| if lo[i] <= hi[i] { lo[i] } else { 0.5 * (lo[i] + hi[i]) });
    let hi2 = [0, 1, 2].map(|i| if lo[i] <= hi[i] { hi[i] } else { 0.5 * (lo[i] + hi[i]) });
    [lo2, hi2]
}
