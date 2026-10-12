//! The scalar (per-vertex) formulas of the canonical problem.
//!
//! At every unknown vertex the cone angle `µ` is charged `g(µ) = t·|µ| + ι_[lo, hi](µ)` with
//! threshold `t = λ·wR ≥ 0` and bounds `lo ≤ 0 ≤ hi` (possibly infinite). A forbidden vertex
//! (`wR = +∞`) and a prescribed cone are both *fixed* sites: `µ = a` is forced at no cost, which
//! is `t = 0`, `lo = hi = a`.
//!
//! Everything the solver and the certificate need from `g` is here: the soft threshold, the
//! proximal map and its pieces, the regularized measure, the subdifferential, and the convex
//! conjugate.

/// The soft threshold `P_t(x) = max(0, x − t) + min(0, x + t)` (paper Eq. 18), for `t ≥ 0`.
///
/// It is the proximal map of `t·|·|`; `µ = P_t(φ + µ)` holds exactly when `φ ∈ t·∂|µ|`
/// (paper Eq. 19).
///
/// [formula: soft-threshold-fixed-point]
#[inline]
pub fn soft_threshold(t: f64, x: f64) -> f64 {
    (x - t).max(0.0) + (x + t).min(0.0)
}

/// What the problem allows and charges at one unknown vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Site {
    /// Threshold `t = λ·wR` (0 at a free site and at a fixed site).
    pub t: f64,
    /// Lower bound of the cone angle (`−∞` if unbounded; the fixed angle at a fixed site).
    pub lo: f64,
    /// Upper bound of the cone angle (`+∞` if unbounded; the fixed angle at a fixed site).
    pub hi: f64,
    /// Whether the cone angle is fixed (forbidden vertex or prescribed cone).
    pub fixed: bool,
}

/// The affine pieces of `x ↦ clamp(k·P_t(x), lo, hi)`, in the order of increasing `x`.
///
/// [formula: soft-threshold-active-sets]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Piece {
    /// Clamped at the lower bound: value `lo`, slope 0.
    Lower,
    /// Negative cone: value `k (x + t)`, slope `k`.
    Neg,
    /// No cone: value 0, slope 0.
    Zero,
    /// Positive cone: value `k (x − t)`, slope `k`.
    Pos,
    /// Clamped at the upper bound: value `hi`, slope 0.
    Upper,
    /// Fixed site: the constant fixed angle.
    Fixed,
}

impl Piece {
    /// Whether the piece has positive slope (the paper's active sets `A⁺ ∪ A⁻`, Eq. 24).
    #[inline]
    pub fn is_active(self) -> bool {
        matches!(self, Piece::Neg | Piece::Pos)
    }
}

impl Site {
    /// A site with cone cost `t·|µ|` and bounds `lo ≤ 0 ≤ hi`.
    pub fn penalized(t: f64, lo: f64, hi: f64) -> Site {
        Site {
            t,
            lo,
            hi,
            fixed: false,
        }
    }

    /// A site whose cone angle is fixed to `angle` at no cost (`angle = 0`: forbidden vertex).
    pub fn fixed(angle: f64) -> Site {
        Site {
            t: 0.0,
            lo: angle,
            hi: angle,
            fixed: true,
        }
    }

    /// The piece of `x ↦ clamp(k·P_t(x), lo, hi)` that `x` lies on (`k > 0`).
    ///
    /// **Tie convention.** The two sloped pieces are *closed* intervals,
    /// `Pos = [t, t + hi/k]` and `Neg = [−t + lo/k, −t]`: a point exactly on a kink counts as
    /// active. A sloped piece of zero width (`hi = 0` resp. `lo = 0`) does not exist. In
    /// particular a free site (`t = 0`) is active at `x = 0`, hence always.
    ///
    /// The paper (Eq. 24) and `Site.pieceReg` of `lean/proofs` use strict inequalities at the
    /// thresholds (`Pos = (t, t + hi/k]`), the Lean oracle's `coneState` gives both kinks to
    /// the flat pieces (`Pos = (t, t + hi/k)`). The three classifiers differ only on the
    /// kinks, where the two adjacent affine formulas have the same value; each satisfies
    /// "the formula of the piece of `x` gives the value at `x`" ([`Site::value_on`] of
    /// `piece(k, x)` is [`Site::measure`]), which is all the stopping rule "pieces unchanged ⟹
    /// exact solution" needs. For this classifier that property is checked by the unit test
    /// `measure_is_the_clamped_soft_threshold` on a grid that contains the kinks; it is not
    /// machine-checked (the Lean instances are for the other two classifiers).
    ///
    /// [formula: soft-threshold-active-sets]
    #[inline]
    pub fn piece(&self, k: f64, x: f64) -> Piece {
        if self.fixed {
            Piece::Fixed
        } else if self.is_two_sided_free() {
            // Positive and negative cones lie on one line: a single sloped piece.
            let y = k * x;
            if y > self.hi {
                Piece::Upper
            } else if y < self.lo {
                Piece::Lower
            } else {
                Piece::Pos
            }
        } else if x >= self.t && self.hi > 0.0 {
            if k * (x - self.t) <= self.hi {
                Piece::Pos
            } else {
                Piece::Upper
            }
        } else if x <= -self.t && self.lo < 0.0 {
            if k * (x + self.t) >= self.lo {
                Piece::Neg
            } else {
                Piece::Lower
            }
        } else {
            Piece::Zero
        }
    }

    /// A free site (`t = 0`) that admits cones of both signs. Its two sloped pieces are one
    /// line, which is labelled [`Piece::Pos`] throughout.
    #[inline]
    fn is_two_sided_free(&self) -> bool {
        !self.fixed && self.t == 0.0 && self.lo < 0.0 && self.hi > 0.0
    }

    /// The piece of a vertex after a Newton step that was computed with the vertex on the
    /// piece `current`.
    ///
    /// `x` is the new argument of `clamp(k·P_t(·))` and `mu` the new cone angle. For a vertex
    /// that was active, `x = ±t + mu/k` holds by the step's own equation, so the new piece is
    /// read off `mu`; evaluating `x − t` instead would lose all digits for large `k`. This is
    /// the same decision as `piece(k, x)` in exact arithmetic.
    ///
    /// [formula: soft-threshold-active-sets]
    #[inline]
    pub fn next_piece(&self, k: f64, current: Piece, x: f64, mu: f64) -> Piece {
        match current {
            Piece::Pos if self.is_two_sided_free() => {
                if mu > self.hi {
                    Piece::Upper
                } else if mu < self.lo {
                    Piece::Lower
                } else {
                    Piece::Pos
                }
            }
            Piece::Pos => {
                if mu > self.hi {
                    Piece::Upper
                } else if mu >= 0.0 {
                    Piece::Pos
                } else {
                    // Strictly below the kink at t, even if `t + mu/k` rounds to `t`.
                    match self.piece(k, self.t + mu / k) {
                        Piece::Pos => Piece::Zero,
                        piece => piece,
                    }
                }
            }
            Piece::Neg => {
                if mu < self.lo {
                    Piece::Lower
                } else if mu <= 0.0 {
                    Piece::Neg
                } else {
                    match self.piece(k, -self.t + mu / k) {
                        Piece::Neg => Piece::Zero,
                        piece => piece,
                    }
                }
            }
            _ => self.piece(k, x),
        }
    }

    /// The pieces that exist for this site, in the order of increasing argument.
    fn chain(&self) -> ([Piece; 5], usize) {
        let mut pieces = [Piece::Fixed; 5];
        let mut len = 0;
        let mut push = |piece| {
            pieces[len] = piece;
            len += 1;
        };
        if self.fixed {
            push(Piece::Fixed);
        } else if self.is_two_sided_free() {
            if self.lo.is_finite() {
                push(Piece::Lower);
            }
            push(Piece::Pos);
            if self.hi.is_finite() {
                push(Piece::Upper);
            }
        } else {
            if self.lo < 0.0 {
                if self.lo.is_finite() {
                    push(Piece::Lower);
                }
                push(Piece::Neg);
            }
            push(Piece::Zero);
            if self.hi > 0.0 {
                push(Piece::Pos);
                if self.hi.is_finite() {
                    push(Piece::Upper);
                }
            }
        }
        (pieces, len)
    }

    /// One step from the piece `current` towards the piece `target` along the chain
    /// lower bound – negative cone – no cone – positive cone – upper bound (restricted to the
    /// pieces that exist for the site). Equal to `current` exactly when `target` is.
    ///
    /// The exact Newton iteration uses this as a damping: a vertex changes to a neighbouring
    /// piece per step and never jumps, for example, from a positive to a negative cone.
    ///
    /// [formula: soft-threshold-active-sets]
    pub fn toward(&self, current: Piece, target: Piece) -> Piece {
        if current == target {
            return current;
        }
        let (chain, len) = self.chain();
        let position = |piece: Piece| chain[..len].iter().position(|&p| p == piece);
        match (position(current), position(target)) {
            (Some(a), Some(b)) if a < b => chain[a + 1],
            (Some(a), Some(b)) if a > b => chain[a - 1],
            _ => target,
        }
    }

    /// The value of `clamp(k·P_t(x), lo, hi)` on a given piece.
    ///
    /// [formula: soft-threshold-active-sets]
    #[inline]
    pub fn value_on(&self, piece: Piece, k: f64, x: f64) -> f64 {
        match piece {
            Piece::Lower => self.lo,
            Piece::Neg => k * (x + self.t),
            Piece::Zero => 0.0,
            Piece::Pos => k * (x - self.t),
            Piece::Upper => self.hi,
            Piece::Fixed => self.lo,
        }
    }

    /// The regularized measure `clamp(k·P_t(x), lo, hi)`. With `k = m/γ` and `x = φ` it is the
    /// cone angle of the regularized problem; with `k = 1` and `x = φ + µ` it is the proximal
    /// map whose fixed points are the solutions of the inclusion `φ ∈ ∂g(µ)`.
    ///
    /// [formula: regularized-cone-update]
    #[inline]
    pub fn measure(&self, k: f64, x: f64) -> f64 {
        self.value_on(self.piece(k, x), k, x)
    }

    /// The proximal map `clamp(P_t(z), lo, hi)` of `g`; `µ = prox(φ + µ)` is the cone equation
    /// of the optimality system.
    ///
    /// [formula: cone-fixed-point-bounded]
    #[inline]
    pub fn prox(&self, z: f64) -> f64 {
        if self.fixed {
            self.lo
        } else {
            soft_threshold(self.t, z).clamp(self.lo, self.hi)
        }
    }

    /// The kinks of `x ↦ clamp(k·P_t(x), lo, hi)`: the points where its slope changes. Up to
    /// four; unused slots are NaN.
    ///
    /// [formula: soft-threshold-active-sets]
    #[inline]
    pub fn kinks(&self, k: f64) -> [f64; 4] {
        let mut out = [f64::NAN; 4];
        if self.fixed {
            return out;
        }
        if self.hi > 0.0 {
            out[0] = self.t;
            if self.hi.is_finite() {
                out[1] = self.t + self.hi / k;
            }
        }
        if self.lo < 0.0 {
            out[2] = -self.t;
            if self.lo.is_finite() {
                out[3] = -self.t + self.lo / k;
            }
        }
        out
    }

    /// The cone cost `t·|µ|` (zero at a fixed site).
    #[inline]
    pub fn cost(&self, mu: f64) -> f64 {
        if self.fixed || self.t == 0.0 {
            0.0
        } else {
            self.t * mu.abs()
        }
    }

    /// Distance of `µ` from the admissible set (`[lo, hi]`, or the fixed angle).
    #[inline]
    pub fn bound_violation(&self, mu: f64) -> f64 {
        (self.lo - mu).max(mu - self.hi).max(0.0)
    }

    /// The subdifferential `∂g(µ)` of `g = t|·| + ι_[lo,hi]` at an admissible `µ`, as a closed
    /// interval `[a, b]` (ends may be infinite):
    /// `a = −∞` if `µ = lo`, else `t` if `µ > 0`, else `−t`;
    /// `b = +∞` if `µ = hi`, else `t` if `µ ≥ 0`, else `−t`.
    /// At a fixed site it is the whole line.
    ///
    /// [formula: abs-subdifferential-bounded]
    #[inline]
    pub fn subdifferential(&self, mu: f64) -> (f64, f64) {
        if self.fixed {
            return (f64::NEG_INFINITY, f64::INFINITY);
        }
        let a = if mu <= self.lo {
            f64::NEG_INFINITY
        } else if mu > 0.0 {
            self.t
        } else {
            -self.t
        };
        let b = if mu >= self.hi {
            f64::INFINITY
        } else if mu >= 0.0 {
            self.t
        } else {
            -self.t
        };
        (a, b)
    }

    /// The domain of the conjugate `g*` as a closed interval: `φ ≤ t` is required where the
    /// upper bound is infinite and `−t ≤ φ` where the lower bound is infinite.
    ///
    /// [formula: duality-gap-certificate]
    #[inline]
    pub fn dual_domain(&self) -> (f64, f64) {
        if self.fixed {
            return (f64::NEG_INFINITY, f64::INFINITY);
        }
        let a = if self.lo == f64::NEG_INFINITY {
            -self.t
        } else {
            f64::NEG_INFINITY
        };
        let b = if self.hi == f64::INFINITY {
            self.t
        } else {
            f64::INFINITY
        };
        (a, b)
    }

    /// The convex conjugate `g*(p) = sup_{lo ≤ µ ≤ hi} (p µ − t|µ|)` on its domain:
    /// `hi·max(0, p − t) + lo·min(0, p + t)`, where a term with an infinite bound is absent;
    /// `a·p` at a site fixed to `a`.
    ///
    /// [formula: duality-gap-certificate]
    #[inline]
    pub fn conjugate(&self, p: f64) -> f64 {
        if self.fixed {
            return self.lo * p;
        }
        let mut value = 0.0;
        if self.hi.is_finite() {
            value += self.hi * (p - self.t).max(0.0);
        }
        if self.lo.is_finite() {
            value += self.lo * (p + self.t).min(0.0);
        }
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const INF: f64 = f64::INFINITY;

    fn grid() -> Vec<f64> {
        let mut v: Vec<f64> = (-40..=40).map(|i| i as f64 * 0.125).collect();
        v.extend([1e-9, -1e-9, 3.0000001, -2.9999999]);
        v
    }

    fn sites() -> Vec<Site> {
        vec![
            Site::penalized(0.5, -INF, INF),
            Site::penalized(0.0, -INF, INF),
            Site::penalized(0.5, 0.0, INF),
            Site::penalized(0.5, -INF, 0.0),
            Site::penalized(0.5, -1.0, 2.0),
            Site::penalized(0.0, -1.0, 2.0),
            Site::penalized(1.5, 0.0, 0.0),
            Site::penalized(0.25, -0.5, 0.0),
            Site::fixed(0.0),
            Site::fixed(1.25),
            Site::fixed(-0.75),
        ]
    }

    #[test]
    fn soft_threshold_values() {
        assert_eq!(soft_threshold(1.0, 3.0), 2.0);
        assert_eq!(soft_threshold(1.0, -3.0), -2.0);
        assert_eq!(soft_threshold(1.0, 0.5), 0.0);
        assert_eq!(soft_threshold(1.0, 1.0), 0.0);
        assert_eq!(soft_threshold(0.0, -0.25), -0.25);
    }

    /// `measure` agrees with the literal `clamp(k P_t(x))` and with its piecewise description.
    #[test]
    fn measure_is_the_clamped_soft_threshold() {
        for s in sites() {
            for k in [1.0, 0.25, 8.0] {
                for x in grid() {
                    let literal = if s.fixed {
                        s.lo
                    } else {
                        (k * soft_threshold(s.t, x)).clamp(s.lo, s.hi)
                    };
                    let got = s.measure(k, x);
                    assert!((got - literal).abs() <= 1e-15, "{s:?} k={k} x={x}");
                    if k == 1.0 {
                        assert_eq!(s.prox(x), literal);
                    }
                }
            }
        }
    }

    /// The slope of `measure` changes only at the listed kinks.
    #[test]
    fn kinks_are_where_the_piece_changes() {
        for s in sites() {
            for k in [1.0, 4.0] {
                let kinks: Vec<f64> = s.kinks(k).into_iter().filter(|x| !x.is_nan()).collect();
                let mut xs = grid();
                xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
                for w in xs.windows(2) {
                    let (p, q) = (s.piece(k, w[0]), s.piece(k, w[1]));
                    if p != q {
                        assert!(
                            kinks.iter().any(|&c| w[0] <= c && c <= w[1]),
                            "{s:?} k={k}: piece changes in [{}, {}] without a kink",
                            w[0],
                            w[1]
                        );
                    }
                }
            }
        }
    }

    /// `µ = prox(φ + µ)` holds exactly when `φ` lies in the subdifferential at `µ`.
    #[test]
    fn prox_fixed_point_is_the_inclusion() {
        for s in sites() {
            for mu in grid() {
                if s.bound_violation(mu) > 0.0 {
                    continue;
                }
                let (a, b) = s.subdifferential(mu);
                for phi in grid() {
                    let inside = a <= phi && phi <= b;
                    // The grid contains values that are not dyadic, hence the tolerance.
                    let fixed_point = (s.prox(phi + mu) - mu).abs() <= 1e-12;
                    assert_eq!(inside, fixed_point, "{s:?} mu={mu} phi={phi}");
                }
            }
        }
    }

    /// Fenchel–Young: `p µ ≤ g(µ) + g*(p)` on the domain of `g*`, with equality exactly on the
    /// subdifferential.
    #[test]
    fn conjugate_satisfies_fenchel_young() {
        for s in sites() {
            let (d0, d1) = s.dual_domain();
            for p in grid() {
                if p < d0 || p > d1 {
                    continue;
                }
                let mut best = f64::NEG_INFINITY;
                for mu in grid() {
                    if s.bound_violation(mu) > 0.0 {
                        continue;
                    }
                    let slack = s.cost(mu) + s.conjugate(p) - p * mu;
                    assert!(slack >= -1e-12, "{s:?} p={p} mu={mu}: slack {slack}");
                    let (a, b) = s.subdifferential(mu);
                    if a <= p && p <= b {
                        assert!(slack.abs() <= 1e-12, "{s:?} p={p} mu={mu}: slack {slack}");
                    }
                    best = best.max(p * mu - s.cost(mu));
                }
                // The supremum is attained on the grid whenever the bounds are grid points.
                if s.fixed || (s.lo.is_finite() && s.hi.is_finite()) {
                    assert!((best - s.conjugate(p)).abs() <= 1e-12, "{s:?} p={p}");
                }
            }
        }
    }

    /// `next_piece` decides like `piece` at the point the step's equation implies.
    #[test]
    fn next_piece_agrees_with_piece() {
        for s in sites() {
            for k in [1.0, 4.0] {
                for mu in grid() {
                    for (current, sign) in [(Piece::Pos, 1.0), (Piece::Neg, -1.0)] {
                        // The piece must exist for the site.
                        if s.fixed || (sign > 0.0 && s.hi <= 0.0) || (sign < 0.0 && s.lo >= 0.0) {
                            continue;
                        }
                        if current == Piece::Neg && s.is_two_sided_free() {
                            continue;
                        }
                        let x = sign * s.t + mu / k;
                        let direct = s.piece(k, x);
                        let next = s.next_piece(k, current, x, mu);
                        // Equal values always; equal labels away from the kinks.
                        assert!(
                            (s.value_on(next, k, x) - s.measure(k, x)).abs() <= 1e-12,
                            "{s:?} k={k} mu={mu} {current:?}: {next:?} vs {direct:?}"
                        );
                        let at_kink = s
                            .kinks(k)
                            .iter()
                            .any(|c| !c.is_nan() && (c - x).abs() <= 1e-12);
                        if !at_kink {
                            assert_eq!(next, direct, "{s:?} k={k} mu={mu} {current:?}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn toward_moves_one_piece_at_a_time() {
        let s = Site::penalized(0.5, -1.0, 2.0);
        assert_eq!(s.toward(Piece::Pos, Piece::Neg), Piece::Zero);
        assert_eq!(s.toward(Piece::Zero, Piece::Upper), Piece::Pos);
        assert_eq!(s.toward(Piece::Upper, Piece::Lower), Piece::Pos);
        assert_eq!(s.toward(Piece::Lower, Piece::Neg), Piece::Neg);
        assert_eq!(s.toward(Piece::Zero, Piece::Zero), Piece::Zero);
        // Pieces that do not exist are skipped.
        let positive_only = Site::penalized(0.5, 0.0, INF);
        assert_eq!(positive_only.toward(Piece::Pos, Piece::Zero), Piece::Zero);
        let free = Site::penalized(0.0, -1.0, 2.0);
        assert_eq!(free.toward(Piece::Upper, Piece::Lower), Piece::Pos);
        assert_eq!(free.toward(Piece::Pos, Piece::Lower), Piece::Lower);
        assert_eq!(
            Site::fixed(1.0).toward(Piece::Fixed, Piece::Fixed),
            Piece::Fixed
        );
        // Every piece `piece()` can return is on the chain of its site.
        for s in sites() {
            let (chain, len) = s.chain();
            for k in [1.0, 4.0] {
                for x in grid() {
                    assert!(chain[..len].contains(&s.piece(k, x)), "{s:?} {x}");
                }
            }
        }
    }

    #[test]
    fn free_site_is_active_at_zero() {
        let s = Site::penalized(0.0, -INF, INF);
        assert_eq!(s.piece(1.0, 0.0), Piece::Pos);
        assert_eq!(s.piece(1.0, -3.0), Piece::Pos);
        assert_eq!(s.measure(2.0, -3.0), -6.0);
        let neg_only = Site::penalized(0.0, -INF, 0.0);
        assert_eq!(neg_only.piece(1.0, 0.0), Piece::Neg);
        assert_eq!(neg_only.piece(1.0, 1.0), Piece::Zero);
    }
}
