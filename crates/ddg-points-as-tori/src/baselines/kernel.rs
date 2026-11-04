//! Kernel baselines on the PAT neighbor sets (SPEC.md §9; paper Eq. 28–29):
//!
//! * **SHC**, signed Hopf–Cole: `w = Σ A_i (λρ_i + 1)⟨x − p_i, n_i⟩/(2πρ_i³) e^{−λ(ρ_i − s)}`,
//!   `φ = sgn(w)(−log|w|/λ + s)`, `∇φ = −∇w/(λ|w|)` with `λ` and the shift frozen. Terms with
//!   `ρ = 0` are left out and counted; `w = 0` is [`PatError::DegenerateQuery`].
//! * **SSPD**, smoothed signed planar distance: `φ = Σ A⟨x − p, n⟩e / Σ A e`,
//!   `∇φ = Σ w (n − λ(g − φ) r̂)/W` (`r̂ = 0` at `ρ = 0`).
//! * **LSE**, the LogSumExp (Madan–Levin 2022) unsigned distance:
//!   `φ = −log(Σ e^{−λ(ρ − s)})/λ + s`, `∇φ = Σ e r̂ / Σ e`. No areas.
//!
//! The neighbor set, `σ_x`, `λ` and the shift `s` follow the PAT blend exactly
//! ([`crate::blend::shift`], [`LambdaRule`]): `s = σ_x` for `PerQuery`, `s = min ρ` for the
//! constant rules (SPEC resolution 1 after WP3). Sums run over the set in its order: kNN and
//! radius sets sorted by `(ρ², index)`, `All` (and kNN with `k ≥ n`) in index order (WP1
//! resolution 6). The operation order is the Lean oracle's (`lean/oracle/README.md`,
//! "Baselines").

use super::{CHUNK, check_len, check_queries, dot, norm, sub};
use crate::blend::{LambdaRule, shift};
use crate::field::{Neighborhood, mean_neighbor_distance};
use crate::spatial::{KdTree, KdTreeOptions, Neighbor};
use crate::{PatError, PointCloud};
use rayon::prelude::*;
use std::f64::consts::PI;

/// Which kernel baseline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelMethod {
    /// Signed Hopf–Cole (Eq. 28).
    Shc,
    /// Smoothed signed planar distance (Eq. 29).
    Sspd,
    /// LogSumExp unsigned distance (Madan–Levin 2022).
    Lse,
}

/// Options of [`KernelField`] (defaults: the PAT reference, `Knn { k: 32 }`, `PerQuery { c: 64 }`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct KernelOptions {
    /// Neighbor set of a query.
    pub neighborhood: Neighborhood,
    /// The λ rule.
    pub lambda: LambdaRule,
    /// Threads of batch evaluation (`None`: all cores).
    pub threads: Option<usize>,
    /// SHC where `w = 0` exactly (e.g. a query in the plane of a planar patch, where every
    /// `⟨x − p_i, n_i⟩ = 0`): `None` (default) returns [`PatError::DegenerateQuery`] as the SPEC
    /// requires; `Some(v)` returns `φ = v`, `∇φ = 0` and counts the query in
    /// [`KernelReport::degenerate`] (for error statistics over grids). `|φ| → ∞` as `w → 0`, so no
    /// value is right.
    pub shc_zero_value: Option<f64>,
}

/// One kernel evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct KernelValue {
    /// `φ(x)`.
    pub value: f64,
    /// `∇φ(x)` (λ and shift frozen); zero when not requested.
    pub gradient: [f64; 3],
    /// `σ_x = ½ max ρ`.
    pub sigma: f64,
    /// `λ` used.
    pub lambda: f64,
    /// The exponent shift `s`.
    pub shift: f64,
    /// Size of the neighbor set.
    pub terms: u32,
    /// SHC terms with `ρ = 0` left out; for SSPD and LSE, terms with `r̂ := 0`.
    pub skipped: u32,
    /// SHC with `w = 0` replaced by [`KernelOptions::shc_zero_value`].
    pub degenerate: bool,
}

/// What a batch evaluation did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct KernelReport {
    /// Queries.
    pub queries: usize,
    /// Neighbor terms summed.
    pub terms: u64,
    /// Terms with `ρ = 0` (left out for SHC, `r̂ := 0` otherwise).
    pub skipped: u64,
    /// Radius queries that fell back to kNN.
    pub radius_fallbacks: u64,
    /// SHC queries with `w = 0` given [`KernelOptions::shc_zero_value`].
    pub degenerate: u64,
    /// Wall time.
    pub seconds: f64,
}

/// A cloud with point areas, a kd-tree and the kernel options.
#[derive(Debug)]
pub struct KernelField {
    positions: Vec<[f64; 3]>,
    normals: Vec<[f64; 3]>,
    areas: Vec<f64>,
    tree: KdTree,
    options: KernelOptions,
    lambda_global: f64,
    pool: rayon::ThreadPool,
}

enum Fail {
    Degenerate(&'static str),
}

impl KernelField {
    /// A field over `cloud` with point areas `areas` (one per point, finite, `≥ 0`; e.g. from
    /// [`super::point_areas`]; LSE ignores them).
    pub fn new(cloud: &PointCloud, areas: Vec<f64>, options: KernelOptions) -> Result<Self, PatError> {
        let n = cloud.len();
        if areas.len() != n {
            return Err(PatError::SizeMismatch { what: format!("{} areas for {n} points", areas.len()) });
        }
        if let Some(i) = areas.iter().position(|a| !a.is_finite()) {
            return Err(PatError::NonFinite { index: i, what: "point area" });
        }
        if let Some(i) = areas.iter().position(|&a| a < 0.0) {
            return Err(PatError::InvalidOption { what: format!("point area {i} is negative ({})", areas[i]) });
        }
        if let Some(z) = options.shc_zero_value
            && !z.is_finite()
        {
            return Err(PatError::InvalidOption { what: format!("KernelOptions::shc_zero_value = {z} (must be finite)") });
        }
        match options.neighborhood {
            Neighborhood::Knn { k: 0 } => return Err(PatError::InvalidOption { what: "Neighborhood::Knn { k: 0 }".into() }),
            Neighborhood::Radius { r, fallback_k } if !(r.is_finite() && r >= 0.0) || fallback_k == 0 => {
                return Err(PatError::InvalidOption { what: format!("Neighborhood::Radius {{ r: {r}, fallback_k: {fallback_k} }}") });
            }
            _ => {}
        }
        match options.lambda {
            LambdaRule::PerQuery { c } if !(c.is_finite() && c > 0.0) => {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::PerQuery {{ c: {c} }}") });
            }
            LambdaRule::Constant { lambda } if !(lambda.is_finite() && lambda > 0.0) => {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::Constant {{ lambda: {lambda} }}") });
            }
            LambdaRule::PaperHeuristic { factor, k_mean, .. } if !(factor.is_finite() && factor > 0.0) || k_mean == 0 => {
                return Err(PatError::InvalidOption {
                    what: format!("LambdaRule::PaperHeuristic {{ factor: {factor}, k_mean: {k_mean} }}"),
                });
            }
            _ => {}
        }
        let pool = crate::pool::make_pool(options.threads)?;
        let positions = cloud.positions().to_vec();
        let (tree, lambda_global) = pool.install(|| -> Result<_, PatError> {
            let tree = KdTree::new(&positions, KdTreeOptions::default())?;
            let lg = match options.lambda {
                LambdaRule::PaperHeuristic { factor, k_mean, .. } => {
                    let d = mean_neighbor_distance(&tree, &positions, k_mean);
                    if d <= 0.0 {
                        return Err(PatError::InvalidOption {
                            what: "PaperHeuristic: the mean neighbor distance D is 0 (all points coincide)".into(),
                        });
                    }
                    factor / d
                }
                _ => 0.0,
            };
            Ok((tree, lg))
        })?;
        Ok(KernelField { positions, normals: cloud.normals().to_vec(), areas, tree, options, lambda_global, pool })
    }

    /// The options.
    pub fn options(&self) -> &KernelOptions {
        &self.options
    }

    /// The point areas.
    pub fn areas(&self) -> &[f64] {
        &self.areas
    }

    fn neighbors(&self, x: [f64; 3], buf: &mut Vec<Neighbor>, idx: &mut Vec<u32>) -> bool {
        let n = self.positions.len();
        idx.clear();
        let mut fell_back = false;
        match self.options.neighborhood {
            Neighborhood::All => idx.extend(0..n as u32),
            Neighborhood::Knn { k } if k >= n => idx.extend(0..n as u32),
            Neighborhood::Knn { k } => {
                self.tree.knn(x, k, buf);
                idx.extend(buf.iter().map(|nb| nb.index));
            }
            Neighborhood::Radius { r, fallback_k } => {
                self.tree.radius(x, r, buf);
                if buf.is_empty() {
                    fell_back = true;
                    self.tree.knn(x, fallback_k, buf);
                }
                idx.extend(buf.iter().map(|nb| nb.index));
            }
        }
        fell_back
    }

    fn lambda_shift(&self, x: [f64; 3], nbrs: &[u32]) -> Result<(f64, f64, f64), Fail> {
        let (sigma, rho_min, _) = shift(x, &self.positions, nbrs);
        match self.options.lambda {
            LambdaRule::PerQuery { c } => {
                if sigma == 0.0 {
                    return Err(Fail::Degenerate(
                        "every neighbor is at the query, so σ_x = 0 and the PerQuery λ = C/σ_x is undefined",
                    ));
                }
                Ok((sigma, c / sigma, sigma))
            }
            LambdaRule::Constant { lambda } => Ok((sigma, lambda, rho_min)),
            LambdaRule::PaperHeuristic { .. } => Ok((sigma, self.lambda_global, rho_min)),
        }
    }

    fn eval_set(&self, method: KernelMethod, x: [f64; 3], nbrs: &[u32], grad: bool) -> Result<KernelValue, Fail> {
        let (sigma, lambda, s) = self.lambda_shift(x, nbrs)?;
        let mut v = KernelValue { sigma, lambda, shift: s, terms: nbrs.len() as u32, ..Default::default() };
        match method {
            KernelMethod::Shc => {
                let (value, gradient, skipped) = shc_sum(&self.positions, &self.normals, &self.areas, x, nbrs, lambda, s, grad);
                v.skipped = skipped;
                if value == 0.0
                    && let Some(z) = self.options.shc_zero_value
                {
                    v.value = z;
                    v.degenerate = true;
                    return Ok(v);
                }
                if value == 0.0 || !value.is_finite() {
                    return Err(Fail::Degenerate("the SHC sum w is 0 (or not finite), so log|w| is undefined"));
                }
                v.value = value.signum() * (-value.abs().ln() / lambda + s);
                if grad {
                    let k = -1.0 / (lambda * value.abs());
                    v.gradient = [k * gradient[0], k * gradient[1], k * gradient[2]];
                }
            }
            KernelMethod::Sspd => {
                let (phi, g, skipped) = sspd_sum(&self.positions, &self.normals, &self.areas, x, nbrs, lambda, s, grad);
                if !phi.is_finite() {
                    return Err(Fail::Degenerate("the SSPD weights sum to 0 (all neighbor areas are 0)"));
                }
                v.value = phi;
                v.gradient = g;
                v.skipped = skipped;
            }
            KernelMethod::Lse => {
                let (phi, g, skipped) = lse_sum(&self.positions, x, nbrs, lambda, s, grad);
                v.value = phi;
                v.gradient = g;
                v.skipped = skipped;
            }
        }
        Ok(v)
    }

    /// One evaluation at `x`, with the gradient if `gradient`.
    pub fn eval(&self, method: KernelMethod, x: [f64; 3], gradient: bool) -> Result<KernelValue, PatError> {
        check_queries(&[x])?;
        let mut idx = Vec::new();
        self.neighbors(x, &mut Vec::new(), &mut idx);
        self.eval_set(method, x, &idx, gradient)
            .map_err(|Fail::Degenerate(r)| PatError::DegenerateQuery { index: 0, reason: r.into() })
    }

    /// Batch evaluation on the field's pool. `gradient`, if given, must have `xs.len()` entries.
    pub fn evaluate(
        &self,
        method: KernelMethod,
        xs: &[[f64; 3]],
        values: &mut [f64],
        mut gradient: Option<&mut [[f64; 3]]>,
    ) -> Result<KernelReport, PatError> {
        let t0 = std::time::Instant::now();
        check_queries(xs)?;
        check_len("value", values.len(), xs.len())?;
        if let Some(g) = gradient.as_deref() {
            check_len("gradient", g.len(), xs.len())?;
        }
        let want = gradient.is_some();
        let res: Vec<Result<(KernelValue, bool), Fail>> = self.pool.install(|| {
            xs.par_chunks(CHUNK)
                .flat_map_iter(|c| {
                    let mut buf = Vec::new();
                    let mut idx = Vec::new();
                    c.iter()
                        .map(|&x| {
                            let fb = self.neighbors(x, &mut buf, &mut idx);
                            self.eval_set(method, x, &idx, want).map(|v| (v, fb))
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        });
        let mut rep = KernelReport { queries: xs.len(), ..Default::default() };
        for (i, r) in res.into_iter().enumerate() {
            let (v, fb) = r.map_err(|Fail::Degenerate(reason)| PatError::DegenerateQuery { index: i, reason: reason.into() })?;
            values[i] = v.value;
            if let Some(g) = gradient.as_deref_mut() {
                g[i] = v.gradient;
            }
            rep.terms += v.terms as u64;
            rep.skipped += v.skipped as u64;
            rep.radius_fallbacks += fb as u64;
            rep.degenerate += v.degenerate as u64;
        }
        rep.seconds = t0.elapsed().as_secs_f64();
        Ok(rep)
    }
}

/// The SHC sum `w` and `∇w` (λ and shift frozen) over `nbrs`; terms with `ρ = 0` are left out
/// (their count is returned).
///
/// `∇w = Σ A e/(2π)[(λρ+1)ρ⁻³ n − ⟨d,n⟩(λ²ρ²+3λρ+3)ρ⁻⁵ d]`, `d = x − p`, `e = e^{−λ(ρ−s)}`.
///
/// [formula: pat-shc]
/// [formula: pat-shc-gradient]
#[allow(clippy::too_many_arguments)]
pub fn shc_sum(
    positions: &[[f64; 3]],
    normals: &[[f64; 3]],
    areas: &[f64],
    x: [f64; 3],
    nbrs: &[u32],
    lambda: f64,
    s: f64,
    grad: bool,
) -> (f64, [f64; 3], u32) {
    let mut w = 0.0;
    let mut gw = [0.0; 3];
    let mut skipped = 0;
    for &i in nbrs {
        let i = i as usize;
        let d = sub(x, positions[i]);
        let rho = norm(d);
        if rho == 0.0 {
            skipped += 1;
            continue;
        }
        let n = normals[i];
        let dn = dot(d, n);
        let lr = lambda * rho;
        let e = areas[i] * (-lambda * (rho - s)).exp();
        let r3 = rho * rho * rho;
        w += e * (lr + 1.0) * dn / (2.0 * PI * r3);
        if grad {
            let a = (lr + 1.0) / r3;
            let b = dn * (lr * lr + 3.0 * lr + 3.0) / (r3 * rho * rho);
            let k = e / (2.0 * PI);
            for t in 0..3 {
                gw[t] += k * (a * n[t] - b * d[t]);
            }
        }
    }
    (w, gw, skipped)
}

/// SSPD `φ = Σ A⟨d,n⟩e / Σ A e` and `∇φ = Σ A e (n − λ(g − φ) r̂)/W` (λ and shift frozen),
/// `r̂ = 0` at `ρ = 0` (counted).
///
/// [formula: pat-sspd]
/// [formula: pat-sspd-gradient]
#[allow(clippy::too_many_arguments)]
pub fn sspd_sum(
    positions: &[[f64; 3]],
    normals: &[[f64; 3]],
    areas: &[f64],
    x: [f64; 3],
    nbrs: &[u32],
    lambda: f64,
    s: f64,
    grad: bool,
) -> (f64, [f64; 3], u32) {
    let mut num = 0.0;
    let mut den = 0.0;
    let mut gn = [0.0; 3]; // Σ w (n − λ g r̂)
    let mut gr = [0.0; 3]; // Σ w r̂
    let mut skipped = 0;
    for &i in nbrs {
        let i = i as usize;
        let d = sub(x, positions[i]);
        let rho = norm(d);
        let n = normals[i];
        let g = dot(d, n);
        let w = areas[i] * (-lambda * (rho - s)).exp();
        num += w * g;
        den += w;
        if grad {
            let rhat = if rho > 0.0 {
                [d[0] / rho, d[1] / rho, d[2] / rho]
            } else {
                skipped += 1;
                [0.0; 3]
            };
            for t in 0..3 {
                gn[t] += w * (n[t] - lambda * g * rhat[t]);
                gr[t] += w * rhat[t];
            }
        } else if rho == 0.0 {
            skipped += 1;
        }
    }
    let phi = num / den;
    let mut gradient = [0.0; 3];
    if grad {
        for t in 0..3 {
            gradient[t] = (gn[t] + lambda * phi * gr[t]) / den;
        }
    }
    (phi, gradient, skipped)
}

/// LSE `φ = −log(Σ e^{−λ(ρ−s)})/λ + s` and `∇φ = Σ e r̂/Σ e` (λ and shift frozen), `r̂ = 0` at
/// `ρ = 0` (counted).
///
/// [formula: pat-lse]
pub fn lse_sum(positions: &[[f64; 3]], x: [f64; 3], nbrs: &[u32], lambda: f64, s: f64, grad: bool) -> (f64, [f64; 3], u32) {
    let mut sum = 0.0;
    let mut g = [0.0; 3];
    let mut skipped = 0;
    for &i in nbrs {
        let d = sub(x, positions[i as usize]);
        let rho = norm(d);
        let e = (-lambda * (rho - s)).exp();
        sum += e;
        if rho > 0.0 {
            if grad {
                for t in 0..3 {
                    g[t] += e * d[t] / rho;
                }
            }
        } else {
            skipped += 1;
        }
    }
    if grad {
        for t in &mut g {
            *t /= sum;
        }
    }
    (-sum.ln() / lambda + s, g, skipped)
}
