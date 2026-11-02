//! The state and the directions of AQP, H¹ L-BFGS and BCQN (see the module documentation).

// `!(x > tol)` is deliberate: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::collections::VecDeque;

use ddg_core::SparseMatrix;
use ddg_repulsive_curves::RcError;
use ddg_repulsive_curves::optimize::Lbfgs;

use crate::error::RsError;
use crate::solve::{SaddleSolver, SolveRecord};
use crate::step::{BcqnVariant, LbfgsVariant, Method, SparseSaddle};

/// What the comparison optimizer did so far (for the reports of the evaluation and the tests).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MethodDiagnostics {
    /// AQP: momentum resets (energy increase or remeshing).
    pub momentum_resets: usize,
    /// L-BFGS/BCQN: memory resets (remeshing, rejected pair, non-descent direction).
    pub memory_resets: usize,
    /// Steps whose remeshing changed the topology (all history was void).
    pub remeshes: usize,
    /// L-BFGS/BCQN: pairs stored now.
    pub pairs: usize,
    /// L-BFGS/BCQN: pairs rejected for non-positive curvature (standard variants).
    pub rejected_pairs: usize,
    /// AQP: the momentum `θ`.
    pub theta: Option<f64>,
    /// BCQN: the last blending factor `β`.
    pub last_beta: Option<f64>,
    /// BCQN: the last norm estimate of the Laplacian.
    pub last_norm: Option<f64>,
}

/// The optimizer state kept by the session between steps.
pub(crate) struct MethodState {
    method: Method,
    diag: MethodDiagnostics,
    start: Option<Vec<[f64; 3]>>,
    prev: Option<Vec<[f64; 3]>>,
    lbfgs: Lbfgs,
    ref_pairs: VecDeque<(Vec<f64>, Vec<f64>)>,
    last: Option<(Vec<f64>, Vec<f64>)>,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// `y = (K ⊗ I₃) x`.
pub fn apply_k3(k: &SparseMatrix, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; x.len()];
    let (cp, ri, val) = (k.col_ptr(), k.row_idx(), k.values());
    for j in 0..k.ncols() {
        for p in cp[j]..cp[j + 1] {
            let i = ri[p] as usize;
            for a in 0..3 {
                y[3 * i + a] += val[p] * x[3 * j + a];
            }
        }
    }
    y
}

/// The spectral norm of a symmetric positive semidefinite `K` by power iteration (fixed
/// deterministic start, at most 200 products, relative change 1e-6).
pub fn spectral_norm(k: &SparseMatrix) -> f64 {
    let n = k.ncols();
    let mut v: Vec<f64> = (0..n).map(|i| 1.0 + 0.5 * (i as f64).sin()).collect();
    let mut lambda = 0.0;
    let mut w = vec![0.0; n];
    for _ in 0..200 {
        let nv = dot(&v, &v).sqrt();
        if !(nv > 0.0) {
            return 0.0;
        }
        v.iter_mut().for_each(|x| *x /= nv);
        if k.mul_vec(&v, &mut w).is_err() {
            return 0.0;
        }
        let l = dot(&v, &w);
        std::mem::swap(&mut v, &mut w);
        let done = (l - lambda).abs() <= 1e-6 * l.abs();
        lambda = l;
        if done {
            break;
        }
    }
    lambda
}

/// Frobenius norm of the constrained Laplacian `[[K ⊗ I₃, C_Qᵀ], [C_Q, 0]]` (the reference's
/// `L.norm()`).
pub(crate) fn saddle_frobenius(k: &SparseMatrix, c: &SparseMatrix) -> f64 {
    let sq = |m: &SparseMatrix| m.values().iter().map(|v| v * v).sum::<f64>();
    (3.0 * sq(k) + 2.0 * sq(c)).sqrt()
}

/// `β = clamp(norm · yᵀ(Ws) / B, 0, 1)` of the blend (Zhu et al. Eq. 13 as coded by the
/// reference).
///
/// [formula: rs-bcqn]
pub(crate) fn blend_beta(norm: f64, y_ws: f64, b: f64) -> f64 {
    let beta = norm * y_ws / b;
    if beta.is_nan() { 0.0 } else { beta.clamp(0.0, 1.0) }
}

/// `q ↦ H₀ q`, the inverse of the constrained proxy.
pub type InverseMetric<'a> = dyn FnMut(&[f64]) -> Result<Vec<f64>, RsError> + 'a;

/// The reference's recursion: every inner product is taken in `W` (`ρ = 1/yᵀWs`,
/// `α = ρ sᵀWq`, `γ = sᵀWy / yᵀWy`, `z = γ H₀ q`, `β = ρ yᵀWz`). `None` for a degenerate pair.
pub fn two_loop_reference(
    pairs: &VecDeque<(Vec<f64>, Vec<f64>)>,
    g: &[f64],
    w: &dyn Fn(&[f64]) -> Vec<f64>,
    h0: &mut InverseMetric<'_>,
) -> Result<Option<Vec<f64>>, RsError> {
    let mut rho = Vec::with_capacity(pairs.len());
    for (s, y) in pairs {
        let r = 1.0 / dot(y, &w(s));
        if !r.is_finite() {
            return Ok(None);
        }
        rho.push(r);
    }
    let mut q = g.to_vec();
    let mut alpha = vec![0.0; pairs.len()];
    for i in (0..pairs.len()).rev() {
        let (s, y) = &pairs[i];
        alpha[i] = rho[i] * dot(s, &w(&q));
        for (qk, yk) in q.iter_mut().zip(y) {
            *qk -= alpha[i] * yk;
        }
    }
    let mut z = h0(&q)?;
    if let Some((s, y)) = pairs.back() {
        let wy = w(y);
        let gamma = dot(s, &wy) / dot(y, &wy);
        if gamma.is_finite() {
            z.iter_mut().for_each(|v| *v *= gamma);
        }
    }
    for (i, (s, y)) in pairs.iter().enumerate() {
        let b = rho[i] * dot(y, &w(&z));
        for (zk, sk) in z.iter_mut().zip(s) {
            *zk += (alpha[i] - b) * sk;
        }
    }
    Ok(Some(z))
}

fn flatten(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

/// `θ = (1 − 1/√κ) / (1 + 1/√κ)` (Kovalsky et al. 2016).
///
/// [formula: rs-aqp]
pub fn aqp_theta(kappa: f64) -> f64 {
    let r = 1.0 / kappa.sqrt();
    (1.0 - r) / (1.0 + r)
}

impl MethodState {
    pub fn new(method: Method) -> MethodState {
        let memory = match method {
            Method::H1Lbfgs { memory, .. } | Method::Bcqn { memory, .. } => memory,
            _ => 1,
        };
        let theta = match method {
            Method::Aqp { kappa } => Some(aqp_theta(kappa)),
            _ => None,
        };
        MethodState {
            method,
            diag: MethodDiagnostics { theta, ..Default::default() },
            start: None,
            prev: None,
            lbfgs: Lbfgs::new(memory),
            ref_pairs: VecDeque::new(),
            last: None,
        }
    }

    /// `true` when the direction is a quasi-Newton one with stored pairs: the natural first
    /// step size is 1.
    pub fn unit_step(&self) -> bool {
        !self.lbfgs.is_empty() || !self.ref_pairs.is_empty()
    }

    pub fn diagnostics(&self) -> MethodDiagnostics {
        let mut d = self.diag;
        d.pairs = self.lbfgs.len().max(self.ref_pairs.len());
        d
    }

    /// AQP: replaces the positions by the Nesterov point `y` (when there is a previous
    /// iterate); remembers `x_k`.
    pub fn begin_step(&mut self, positions: &mut [[f64; 3]]) {
        let Method::Aqp { kappa } = self.method else { return };
        self.start = Some(positions.to_vec());
        if let Some(prev) = &self.prev {
            let th = aqp_theta(kappa);
            for (x, p) in positions.iter_mut().zip(prev) {
                for a in 0..3 {
                    x[a] = (1.0 + th) * x[a] - th * p[a];
                }
            }
        }
    }

    /// After a failed step: gives the positions back.
    pub fn abort_step(&mut self, positions: &mut [[f64; 3]]) {
        if let Some(s) = self.start.take() {
            positions.copy_from_slice(&s);
        }
    }

    /// After a step: AQP keeps `x_k` as the previous iterate unless the energy increased or the
    /// mesh changed; the quasi-Newton methods forget everything when the mesh changed.
    pub fn end_step(&mut self, remeshed: bool, e_before: f64, e_after: f64) {
        if remeshed {
            self.on_remesh();
            self.start = None;
            return;
        }
        if let Method::Aqp { .. } = self.method {
            let start = self.start.take();
            if !(e_after <= e_before) {
                if self.prev.is_some() {
                    self.diag.momentum_resets += 1;
                }
                self.prev = None;
            } else {
                self.prev = start;
            }
        }
    }

    /// Remeshing changed the mesh: all history is void.
    pub fn on_remesh(&mut self) {
        self.diag.remeshes += 1;
        if self.prev.take().is_some() {
            self.diag.momentum_resets += 1;
        }
        if !self.lbfgs.is_empty() || !self.ref_pairs.is_empty() {
            self.diag.memory_resets += 1;
        }
        self.lbfgs.clear();
        self.ref_pairs.clear();
        self.last = None;
    }

    fn clear_memory(&mut self) {
        self.lbfgs.clear();
        self.ref_pairs.clear();
        self.diag.memory_resets += 1;
    }

    /// The direction `x` (the session steps `f − τ x`) at the current positions with the
    /// Euclidean differential `g`.
    pub fn direction(
        &mut self,
        solver: &mut SparseSaddle,
        positions: &[[f64; 3]],
        g: &[f64],
        c_q: &SparseMatrix,
        total_area: f64,
    ) -> Result<(Vec<f64>, SolveRecord), RsError> {
        let zeros = vec![0.0; c_q.nrows()];
        let k = solver.metric_matrix().clone();
        let mut record: Option<SolveRecord> = None;
        let mut h0 = |q: &[f64]| -> Result<Vec<f64>, RsError> {
            let sol = solver.solve(q, &zeros)?;
            record = Some(sol.record);
            Ok(sol.x)
        };
        let (variant_ref, blend) = match self.method {
            Method::Aqp { .. } => {
                let x = h0(g)?;
                return Ok((x, record.expect("one solve")));
            }
            Method::H1Lbfgs { variant, .. } => (variant == LbfgsVariant::Reference, None),
            Method::Bcqn { variant, .. } => (variant == BcqnVariant::Reference, Some(variant)),
            Method::ProjectedGradient => unreachable!("handled by the session"),
        };
        // 1. The pair of the last step.
        let pos = flatten(positions);
        if let Some((p0, g0)) = self.last.take() {
            let s: Vec<f64> = pos.iter().zip(&p0).map(|(a, b)| a - b).collect();
            let y: Vec<f64> = g.iter().zip(&g0).map(|(a, b)| a - b).collect();
            if dot(&s, &s) > 0.0 {
                let y = match blend {
                    Some(v) => {
                        let ws = apply_k3(&k, &s);
                        let norm = match v {
                            BcqnVariant::Standard => spectral_norm(&k),
                            BcqnVariant::Reference => saddle_frobenius(&k, c_q),
                        };
                        let b = total_area.powf(4.0 / 3.0);
                        let beta = blend_beta(norm, dot(&y, &ws), b);
                        self.diag.last_beta = Some(beta);
                        self.diag.last_norm = Some(norm);
                        y.iter().zip(&ws).map(|(yi, wi)| (1.0 - beta) * yi + beta * wi).collect()
                    }
                    None => y,
                };
                if variant_ref {
                    self.ref_pairs.push_back((s, y));
                    let cap = match self.method {
                        Method::H1Lbfgs { memory, .. } | Method::Bcqn { memory, .. } => memory,
                        _ => 1,
                    };
                    while self.ref_pairs.len() > cap {
                        self.ref_pairs.pop_front();
                    }
                } else if !self.lbfgs.push(s, y) {
                    self.diag.rejected_pairs += 1;
                    self.clear_memory();
                }
            }
        }
        // 2. The two-loop recursion.
        let wmul = |v: &[f64]| apply_k3(&k, v);
        let mut x = if variant_ref {
            two_loop_reference(&self.ref_pairs, g, &wmul, &mut h0)?
        } else {
            let mut err: Option<RsError> = None;
            let r = self.lbfgs.direction(g, &mut |q: &[f64]| {
                h0(q).map_err(|e| {
                    err = Some(e);
                    RcError::Cancelled
                })
            });
            match r {
                Ok(z) => Some(z),
                Err(_) => return Err(err.expect("the closure failed")),
            }
        };
        // 3. A descent direction, else forget and take the H1 gradient.
        if x.as_ref().is_none_or(|z| !(dot(g, z) > 0.0) || z.iter().any(|v| !v.is_finite())) {
            self.lbfgs.clear();
            self.ref_pairs.clear();
            self.diag.memory_resets += 1;
            x = Some(h0(g)?);
        }
        self.last = Some((pos, g.to_vec()));
        Ok((x.expect("set above"), record.expect("at least one solve")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aqp_momentum_resets_on_energy_increase_and_remeshing() {
        let mut m = MethodState::new(Method::Aqp { kappa: 100.0 });
        let mut x = vec![[1.0, 0.0, 0.0]; 2];
        m.begin_step(&mut x);
        assert_eq!(x, vec![[1.0, 0.0, 0.0]; 2], "no momentum without history");
        m.end_step(false, 2.0, 1.0);
        let mut x2 = vec![[2.0, 0.0, 0.0]; 2];
        m.begin_step(&mut x2);
        let th = aqp_theta(100.0);
        assert!((x2[0][0] - ((1.0 + th) * 2.0 - th)).abs() < 1e-15);
        // Energy increased: reset.
        m.end_step(false, 1.0, 1.5);
        assert_eq!(m.diagnostics().momentum_resets, 1);
        let mut x3 = vec![[3.0, 0.0, 0.0]; 2];
        m.begin_step(&mut x3);
        assert_eq!(x3[0][0], 3.0);
        m.end_step(false, 1.0, 0.5);
        // Remeshing: reset.
        m.begin_step(&mut x3);
        m.end_step(true, 1.0, 0.5);
        assert_eq!(m.diagnostics().momentum_resets, 2);
        let mut x4 = vec![[4.0, 0.0, 0.0]; 2];
        m.begin_step(&mut x4);
        assert_eq!(x4[0][0], 4.0);
    }
}
