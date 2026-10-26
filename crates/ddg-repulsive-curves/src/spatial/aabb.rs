//! Axis-aligned bounding boxes in `R^D` (`SPEC.md` §8.1, §16).
//!
//! Generic in the dimension; no curve types. The BVH of [`super::bvh`] uses `D = 6` for the
//! tangent-points `(x, T)`, the endpoint guard uses `D = 3`.

/// A closed axis-aligned box `[lo, hi]` in `R^D`. The empty box has `lo = +∞`, `hi = −∞`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb<const D: usize> {
    /// Lower corner.
    pub lo: [f64; D],
    /// Upper corner.
    pub hi: [f64; D],
}

impl<const D: usize> Aabb<D> {
    /// The empty box (the identity of [`Aabb::union`]).
    pub const fn empty() -> Self {
        Aabb {
            lo: [f64::INFINITY; D],
            hi: [f64::NEG_INFINITY; D],
        }
    }

    /// The degenerate box `{p}`.
    pub const fn from_point(p: [f64; D]) -> Self {
        Aabb { lo: p, hi: p }
    }

    /// `true` for the empty box.
    pub fn is_empty(&self) -> bool {
        (0..D).any(|k| self.lo[k] > self.hi[k])
    }

    /// The smallest box containing `self` and the point `p`.
    #[inline]
    pub fn grow(&mut self, p: [f64; D]) {
        for k in 0..D {
            self.lo[k] = self.lo[k].min(p[k]);
            self.hi[k] = self.hi[k].max(p[k]);
        }
    }

    /// The smallest box containing both boxes.
    #[inline]
    pub fn union(&self, other: &Self) -> Self {
        let mut r = *self;
        for k in 0..D {
            r.lo[k] = r.lo[k].min(other.lo[k]);
            r.hi[k] = r.hi[k].max(other.hi[k]);
        }
        r
    }

    /// `true` if `p` lies in the closed box.
    #[inline]
    pub fn contains_point(&self, p: [f64; D]) -> bool {
        (0..D).all(|k| self.lo[k] <= p[k] && p[k] <= self.hi[k])
    }

    /// `true` if `other` lies in the closed box.
    pub fn contains_box(&self, other: &Self) -> bool {
        (0..D).all(|k| self.lo[k] <= other.lo[k] && other.hi[k] <= self.hi[k])
    }

    /// Squared length of the diagonal (`0` for the empty box).
    #[inline]
    pub fn diameter_sq(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        (0..D)
            .map(|k| {
                let w = self.hi[k] - self.lo[k];
                w * w
            })
            .sum()
    }

    /// Half the length of the diagonal, the radius `r` of the paper's admissibility criteria.
    #[inline]
    pub fn half_diagonal(&self) -> f64 {
        0.5 * self.diameter_sq().sqrt()
    }

    /// Half the diagonal of the projection onto the coordinates `range`.
    pub fn half_diagonal_of(&self, range: std::ops::Range<usize>) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        0.5 * range
            .map(|k| {
                let w = self.hi[k] - self.lo[k];
                w * w
            })
            .sum::<f64>()
            .sqrt()
    }
}

impl<const D: usize> Default for Aabb<D> {
    fn default() -> Self {
        Aabb::empty()
    }
}
