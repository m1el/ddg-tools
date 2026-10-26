//! L-BFGS (`SPEC.md` §6.4, D22; Liu & Nocedal 1989; Nocedal & Wright, *Numerical
//! Optimization*, 2nd ed., Alg. 7.4).
//!
//! Two-loop recursion with the dynamic initial inverse Hessian `H₀ = γ_k M_k⁻¹`: `M_k⁻¹ q` is the
//! mean-free Sobolev solve of the current curve (Eq. 22 with the rows `M_c`, Supp. 3.3), and
//! `γ_k = sᵀy / yᵀ M_k⁻¹ y` of the newest pair (Nocedal & Wright Eq. 7.20 with the
//! preconditioner) so that the unit step is natural. The line search is Armijo only, so curvature
//! is not guaranteed: a pair with `sᵀy ≤ 1e-12 |s||y|` is not stored (D22).
//!
//! [formula: rc-lbfgs-two-loop]

use std::collections::VecDeque;

use super::dotv;
use crate::error::RcError;

/// `q ↦ M⁻¹ q`, the preconditioner of the two-loop recursion.
pub type Preconditioner<'a> = dyn FnMut(&[f64]) -> Result<Vec<f64>, RcError> + 'a;

/// The L-BFGS memory: the last `memory` pairs `(s_k, y_k)` of position and gradient differences
/// (flat vectors in the session's layout).
#[derive(Debug, Clone, PartialEq)]
pub struct Lbfgs {
    memory: usize,
    pairs: VecDeque<(Vec<f64>, Vec<f64>, f64)>,
    skipped: usize,
}

impl Lbfgs {
    /// An empty memory of `memory` pairs.
    pub fn new(memory: usize) -> Lbfgs {
        Lbfgs {
            memory: memory.max(1),
            pairs: VecDeque::new(),
            skipped: 0,
        }
    }

    /// Number of stored pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// `true` without pairs.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Number of pairs rejected by the curvature test.
    pub fn skipped(&self) -> usize {
        self.skipped
    }

    /// Forgets all pairs (after a subdivision, or when the direction is not a descent
    /// direction).
    pub fn clear(&mut self) {
        self.pairs.clear();
    }

    /// Stores `(s, y)` if `sᵀy > 1e-12 |s||y|`; returns whether it was stored.
    pub fn push(&mut self, s: Vec<f64>, y: Vec<f64>) -> bool {
        let sy = dotv(&s, &y);
        let ns = dotv(&s, &s).sqrt();
        let ny = dotv(&y, &y).sqrt();
        if !sy.is_finite() || sy <= 1e-12 * ns * ny {
            self.skipped += 1;
            return false;
        }
        if self.pairs.len() == self.memory {
            self.pairs.pop_front();
        }
        self.pairs.push_back((s, y, 1.0 / sy));
        true
    }

    /// The quasi-Newton direction `p = H r` by the two-loop recursion, with `precondition(q)` =
    /// `M⁻¹ q` (the mean-free metric solve of the current curve). Trials are `γ − τ p`.
    ///
    /// [formula: rc-lbfgs-two-loop]
    ///
    /// # Errors
    /// The errors of `precondition`.
    pub fn direction(
        &self,
        r: &[f64],
        precondition: &mut Preconditioner<'_>,
    ) -> Result<Vec<f64>, RcError> {
        let mut q = r.to_vec();
        let mut alpha = vec![0.0; self.pairs.len()];
        for (i, (s, y, rho)) in self.pairs.iter().enumerate().rev() {
            let a = rho * dotv(s, &q);
            alpha[i] = a;
            for (qk, yk) in q.iter_mut().zip(y) {
                *qk -= a * yk;
            }
        }
        let mut z = precondition(&q)?;
        if let Some((s, y, _)) = self.pairs.back() {
            let hy = precondition(y)?;
            let yhy = dotv(y, &hy);
            let gamma = dotv(s, y) / yhy;
            if gamma.is_finite() && gamma > 0.0 {
                for zk in z.iter_mut() {
                    *zk *= gamma;
                }
            }
        }
        for (i, (s, y, rho)) in self.pairs.iter().enumerate() {
            let b = rho * dotv(y, &z);
            for (zk, sk) in z.iter_mut().zip(s) {
                *zk += (alpha[i] - b) * sk;
            }
        }
        Ok(z)
    }
}
