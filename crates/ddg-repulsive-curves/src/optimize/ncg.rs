//! Preconditioned nonlinear conjugate gradients, Polak–Ribière with restarts (`SPEC.md` §6.4,
//! D22; Polak & Ribière 1969, Nocedal & Wright §5.2).
//!
//! With the residual `r_k = dFᵀ` and the preconditioned gradient `z_k = M_k⁻¹ r_k` (the
//! mean-free Sobolev solve of the current curve):
//!
//! ```text
//! β_k = max(0, ⟨r_k, z_k − z_{k−1}⟩ / ⟨r_{k−1}, z_{k−1}⟩)        (PR+)
//! d_k = −z_k + β_k d_{k−1},     restart d_k = −z_k when dF·d_k ≥ 0
//! ```
//!
//! Trials are `γ + τ d_k` (in the session's convention `γ − τ p` with `p = −d_k`).
//!
//! [formula: rc-ncg-polak-ribiere]

use super::dotv;

/// The state of the NCG iteration.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Ncg {
    prev: Option<(Vec<f64>, Vec<f64>, Vec<f64>)>,
    restarts: usize,
}

/// One NCG direction.
#[derive(Debug, Clone, PartialEq)]
pub struct NcgDirection {
    /// `p = −d_k`; trials are `γ − τ p`.
    pub p: Vec<f64>,
    /// `β_k` (0 at a restart).
    pub beta: f64,
    /// Whether the direction was restarted (`d_k = −z_k`).
    pub restarted: bool,
}

impl Ncg {
    /// A fresh iteration.
    pub fn new() -> Ncg {
        Ncg::default()
    }

    /// Forgets the previous direction (after a subdivision).
    pub fn clear(&mut self) {
        self.prev = None;
    }

    /// Number of restarts so far (descent test).
    pub fn restarts(&self) -> usize {
        self.restarts
    }

    /// The direction for the residual `r` and preconditioned gradient `z`; remembers `r`, `z`
    /// and `d` for the next call.
    ///
    /// [formula: rc-ncg-polak-ribiere]
    pub fn direction(&mut self, r: &[f64], z: &[f64]) -> NcgDirection {
        let mut beta = 0.0;
        let mut d: Vec<f64> = z.iter().map(|x| -x).collect();
        if let Some((r0, z0, d0)) = &self.prev
            && r0.len() == r.len()
        {
            let num: f64 = r
                .iter()
                .zip(z.iter().zip(z0))
                .map(|(rk, (zk, zp))| rk * (zk - zp))
                .sum();
            let den = dotv(r0, z0);
            let b = num / den;
            if b.is_finite() && b > 0.0 {
                beta = b;
                for (dk, dp) in d.iter_mut().zip(d0) {
                    *dk += beta * dp;
                }
            }
        }
        let mut restarted = beta == 0.0;
        if dotv(r, &d) >= 0.0 && beta > 0.0 {
            d = z.iter().map(|x| -x).collect();
            beta = 0.0;
            restarted = true;
            self.restarts += 1;
        }
        self.prev = Some((r.to_vec(), z.to_vec(), d.clone()));
        NcgDirection {
            p: d.iter().map(|x| -x).collect(),
            beta,
            restarted,
        }
    }
}
