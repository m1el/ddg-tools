//! The 2D blend (Eq. 25-26 in 2D) and the facade [`Pat2d::build`].

use super::fit::{Circles, CirclePrimitive, Fit2dOptions, Fit2dReport, QuadraticPredictor, circle_sdf, circle_sdf_gradient, fit_circles};
use super::Cloud2;
use crate::blend::{GradientMode, LambdaRule};
use crate::pool::make_pool;
use crate::spatial::{KdTree, KdTreeOptions, Neighbor, fitting_neighborhoods};
use crate::{Neighborhood, PatError};
use rayon::prelude::*;
use std::time::Instant;

/// Options of [`Pat2dField`] (defaults as in 3D, decision D13).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Eval2dOptions {
    /// The neighbor set (default 32 nearest).
    pub neighborhood: Neighborhood,
    /// The λ rule (default `PerQuery { c: 64 }`).
    pub lambda: LambdaRule,
    /// Frozen (reference) or exact gradient.
    pub gradient: GradientMode,
    /// Subtracted from `φ`.
    pub isovalue: f64,
    /// Threads of batch evaluation (`None`: all cores).
    pub threads: Option<usize>,
}

/// What a batch evaluation did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Eval2dReport {
    /// Queries.
    pub queries: usize,
    /// Circle terms evaluated at their center (gradient 0).
    pub singular_center: u64,
    /// Queries coinciding with a data point (`r̂ := 0` for that term).
    pub at_data_point: u64,
    /// Radius queries that fell back to kNN.
    pub radius_fallbacks: u64,
    /// Seconds.
    pub seconds: f64,
}

/// A blended signed distance field of circles.
#[derive(Debug)]
pub struct Pat2dField {
    positions: Vec<[f64; 2]>,
    circles: Circles,
    tree: KdTree,
    options: Eval2dOptions,
    lambda_global: f64,
    pool: rayon::ThreadPool,
}

const CHUNK: usize = 512;

#[derive(Default, Clone, Copy)]
struct Flags {
    center: u64,
    at_point: u64,
    fallback: u64,
}

impl Pat2dField {
    /// Build the field from fitted circles.
    pub fn new(cloud: &Cloud2, circles: Circles, options: Eval2dOptions) -> Result<Pat2dField, PatError> {
        let n = cloud.len();
        if circles.primitives.len() != n {
            return Err(PatError::SizeMismatch { what: format!("{} circles for {n} points", circles.primitives.len()) });
        }
        match options.neighborhood {
            Neighborhood::Knn { k: 0 } => return Err(PatError::InvalidOption { what: "Neighborhood::Knn { k: 0 }".into() }),
            Neighborhood::Radius { r, fallback_k } if !(r.is_finite() && r > 0.0) || fallback_k == 0 => {
                return Err(PatError::InvalidOption { what: format!("Neighborhood::Radius {{ r: {r}, fallback_k: {fallback_k} }}") });
            }
            _ => {}
        }
        match options.lambda {
            LambdaRule::PerQuery { c } if !(c.is_finite() && c > 0.0) => return Err(PatError::InvalidOption { what: format!("LambdaRule::PerQuery c = {c}") }),
            LambdaRule::Constant { lambda } if !(lambda.is_finite() && lambda > 0.0) => {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::Constant lambda = {lambda}") });
            }
            LambdaRule::PaperHeuristic { c, factor, k_mean } if !(c > 0.0 && factor > 0.0 && k_mean > 0) => {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::PaperHeuristic c = {c}, factor = {factor}, k_mean = {k_mean}") });
            }
            _ => {}
        }
        let pool = make_pool(options.threads)?;
        let emb = cloud.embedded();
        let tree = pool.install(|| KdTree::new(&emb, KdTreeOptions::default()))?;
        let lambda_global = match options.lambda {
            LambdaRule::PaperHeuristic { factor, k_mean, .. } => factor / mean_neighbor_distance(&pool, &tree, &emb, k_mean),
            _ => f64::NAN,
        };
        Ok(Pat2dField { positions: cloud.positions().to_vec(), circles, tree, options, lambda_global, pool })
    }

    /// The paper's Sec. 4.2 options for `cloud` in 2D: `λ = 10³/D`, radius `2C/λ`, 32-NN fallback.
    pub fn paper_options(cloud: &Cloud2, threads: Option<usize>) -> Result<(Eval2dOptions, f64), PatError> {
        let (c, factor, k_mean) = (64.0, 1e3, 64);
        let pool = make_pool(threads)?;
        let emb = cloud.embedded();
        let tree = pool.install(|| KdTree::new(&emb, KdTreeOptions::default()))?;
        let d = mean_neighbor_distance(&pool, &tree, &emb, k_mean);
        let r = 2.0 * c / (factor / d);
        Ok((Eval2dOptions { neighborhood: Neighborhood::Radius { r, fallback_k: 32 }, lambda: LambdaRule::PaperHeuristic { c, factor, k_mean }, threads, ..Default::default() }, r))
    }

    /// The circles.
    pub fn circles(&self) -> &Circles {
        &self.circles
    }

    /// The data positions.
    pub fn positions(&self) -> &[[f64; 2]] {
        &self.positions
    }

    /// The options.
    pub fn options(&self) -> &Eval2dOptions {
        &self.options
    }

    fn neighbors(&self, x: [f64; 2], buf: &mut Vec<Neighbor>, idx: &mut Vec<u32>, fl: &mut Flags) {
        let n = self.positions.len();
        idx.clear();
        let q = [x[0], x[1], 0.0];
        match self.options.neighborhood {
            Neighborhood::All => idx.extend(0..n as u32),
            Neighborhood::Knn { k } if k >= n => idx.extend(0..n as u32),
            Neighborhood::Knn { k } => {
                self.tree.knn(q, k, buf);
                idx.extend(buf.iter().map(|nb| nb.index));
            }
            Neighborhood::Radius { r, fallback_k } => {
                self.tree.radius(q, r, buf);
                if buf.is_empty() {
                    fl.fallback += 1;
                    self.tree.knn(q, fallback_k, buf);
                }
                idx.extend(buf.iter().map(|nb| nb.index));
            }
        }
    }

    /// One blend over the neighbor set `idx`: `(φ − isovalue, ∇φ)`.
    ///
    /// [formula: pat2d-blend]
    fn blend(&self, x: [f64; 2], idx: &[u32], want_grad: bool, fl: &mut Flags) -> Option<(f64, [f64; 2])> {
        let m = idx.len();
        let mut rho = Vec::with_capacity(m);
        let mut rmax = (f64::NEG_INFINITY, 0usize);
        let mut rmin = f64::INFINITY;
        for (t, &j) in idx.iter().enumerate() {
            let p = self.positions[j as usize];
            let d = [x[0] - p[0], x[1] - p[1]];
            let r = (d[0] * d[0] + d[1] * d[1]).sqrt();
            if r >= rmax.0 {
                rmax = (r, t);
            }
            rmin = rmin.min(r);
            rho.push(r);
        }
        let sigma = 0.5 * rmax.0;
        let (lambda, shift, c) = match self.options.lambda {
            LambdaRule::PerQuery { c } => {
                if sigma <= 0.0 {
                    return None;
                }
                (c / sigma, sigma, c)
            }
            LambdaRule::Constant { lambda } => (lambda, rmin, f64::NAN),
            LambdaRule::PaperHeuristic { .. } => (self.lambda_global, rmin, f64::NAN),
        };
        let (mut w_sum, mut gw) = (0.0, 0.0);
        let mut ws = Vec::with_capacity(m);
        let mut gs = Vec::with_capacity(m);
        for (t, &j) in idx.iter().enumerate() {
            let w = (-lambda * (rho[t] - shift)).exp();
            let g = circle_sdf(&self.circles.primitives[j as usize], x);
            w_sum += w;
            gw += w * g;
            ws.push(w);
            gs.push(g);
        }
        let phi = gw / w_sum;
        if !want_grad {
            return Some((phi - self.options.isovalue, [0.0; 2]));
        }
        let (mut a, mut b) = ([0.0f64; 2], [0.0f64; 2]);
        let mut ds = 0.0;
        for (t, &j) in idx.iter().enumerate() {
            let p = self.positions[j as usize];
            let rhat = if rho[t] > 0.0 { [(x[0] - p[0]) / rho[t], (x[1] - p[1]) / rho[t]] } else {
                fl.at_point += 1;
                [0.0, 0.0]
            };
            let (dg, sing) = circle_sdf_gradient(&self.circles.primitives[j as usize], x);
            if sing {
                fl.center += 1;
            }
            for k in 0..2 {
                a[k] += ws[t] * (dg[k] - lambda * rhat[k] * gs[t]);
                b[k] += ws[t] * lambda * rhat[k];
            }
            ds += ws[t] * (sigma - rho[t]) * (gs[t] - phi);
        }
        let mut g = [a[0] / w_sum + phi * b[0] / w_sum, a[1] / w_sum + phi * b[1] / w_sum];
        if self.options.gradient == GradientMode::Exact && matches!(self.options.lambda, LambdaRule::PerQuery { .. }) {
            // ∂λ/∂x = −(C/σ²) · ½ r̂_argmax
            let p = self.positions[idx[rmax.1] as usize];
            let r = rho[rmax.1];
            if r > 0.0 {
                let rh = [(x[0] - p[0]) / r, (x[1] - p[1]) / r];
                let f = -(c / (sigma * sigma)) * 0.5 * ds / w_sum;
                g[0] += f * rh[0];
                g[1] += f * rh[1];
            }
        }
        Some((phi - self.options.isovalue, g))
    }

    fn check(x: [f64; 2], index: usize) -> Result<(), PatError> {
        if x.iter().all(|c| c.is_finite()) { Ok(()) } else { Err(PatError::NonFinite { index, what: "query" }) }
    }

    fn degenerate(index: usize) -> PatError {
        PatError::DegenerateQuery { index, reason: "every neighbor is at the query, so σ_x = 0 and the PerQuery λ = C/σ_x is undefined".into() }
    }

    /// `φ(x) − isovalue`.
    pub fn distance(&self, x: [f64; 2]) -> Result<f64, PatError> {
        Self::check(x, 0)?;
        let (mut idx, mut buf, mut fl) = (Vec::new(), Vec::new(), Flags::default());
        self.neighbors(x, &mut buf, &mut idx, &mut fl);
        self.blend(x, &idx, false, &mut fl).map(|v| v.0).ok_or_else(|| Self::degenerate(0))
    }

    /// `φ(x) − isovalue` and `∇φ(x)`.
    pub fn distance_gradient(&self, x: [f64; 2]) -> Result<(f64, [f64; 2]), PatError> {
        Self::check(x, 0)?;
        let (mut idx, mut buf, mut fl) = (Vec::new(), Vec::new(), Flags::default());
        self.neighbors(x, &mut buf, &mut idx, &mut fl);
        self.blend(x, &idx, true, &mut fl).ok_or_else(|| Self::degenerate(0))
    }

    /// Batch evaluation on the field's pool. Output buffers, if given, need one entry per query.
    /// The result does not depend on the thread count.
    pub fn evaluate(&self, xs: &[[f64; 2]], distance: Option<&mut [f64]>, gradient: Option<&mut [[f64; 2]]>) -> Result<Eval2dReport, PatError> {
        let t0 = Instant::now();
        let m = xs.len();
        if distance.as_ref().is_some_and(|d| d.len() != m) || gradient.as_ref().is_some_and(|g| g.len() != m) {
            return Err(PatError::SizeMismatch { what: format!("output buffers for {m} queries") });
        }
        if let Some(i) = xs.iter().position(|x| x.iter().any(|c| !c.is_finite())) {
            return Err(PatError::NonFinite { index: i, what: "query" });
        }
        let want_grad = gradient.is_some();
        let mut dbuf = vec![0.0; m];
        let mut gbuf = vec![[0.0; 2]; if want_grad { m } else { 0 }];
        let results: Vec<(Option<usize>, Flags)> = self.pool.install(|| {
            let chunks = dbuf.par_chunks_mut(CHUNK);
            if want_grad {
                chunks
                    .zip(gbuf.par_chunks_mut(CHUNK))
                    .enumerate()
                    .map(|(c, (d, g))| self.run_chunk(xs, c, d, Some(g)))
                    .collect()
            } else {
                chunks.enumerate().map(|(c, d)| self.run_chunk(xs, c, d, None)).collect()
            }
        });
        if let Some(i) = results.iter().filter_map(|r| r.0).min() {
            return Err(Self::degenerate(i));
        }
        let mut rep = Eval2dReport { queries: m, ..Default::default() };
        for (_, f) in &results {
            rep.singular_center += f.center;
            rep.at_data_point += f.at_point;
            rep.radius_fallbacks += f.fallback;
        }
        if let Some(d) = distance {
            d.copy_from_slice(&dbuf);
        }
        if let Some(g) = gradient {
            g.copy_from_slice(&gbuf);
        }
        rep.seconds = t0.elapsed().as_secs_f64();
        Ok(rep)
    }

    fn run_chunk(&self, xs: &[[f64; 2]], c: usize, d: &mut [f64], mut g: Option<&mut [[f64; 2]]>) -> (Option<usize>, Flags) {
        let (mut idx, mut buf, mut fl) = (Vec::new(), Vec::new(), Flags::default());
        for t in 0..d.len() {
            let i = c * CHUNK + t;
            self.neighbors(xs[i], &mut buf, &mut idx, &mut fl);
            match self.blend(xs[i], &idx, g.is_some(), &mut fl) {
                Some((phi, grad)) => {
                    d[t] = phi;
                    if let Some(g) = g.as_deref_mut() {
                        g[t] = grad;
                    }
                }
                None => return (Some(i), fl),
            }
        }
        (None, fl)
    }
}

fn mean_neighbor_distance(pool: &rayon::ThreadPool, tree: &KdTree, pts: &[[f64; 3]], k: usize) -> f64 {
    let k = k.min(pts.len());
    let sums: Vec<f64> = pool.install(|| {
        pts.par_chunks(256)
            .map(|c| {
                let mut buf = Vec::new();
                let mut s = 0.0;
                for &p in c {
                    tree.knn(p, k, &mut buf);
                    s += buf.iter().map(|nb| nb.d2.sqrt()).sum::<f64>() / buf.len() as f64;
                }
                s
            })
            .collect()
    });
    sums.iter().sum::<f64>() / pts.len() as f64
}

/// Options of [`Pat2d::build`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BuildOptions2 {
    /// Circles from coefficients.
    pub fit: Fit2dOptions,
    /// The blend.
    pub eval: Eval2dOptions,
}

impl BuildOptions2 {
    /// Set the thread count of every stage.
    pub fn with_threads(mut self, threads: Option<usize>) -> Self {
        self.fit.threads = threads;
        self.eval.threads = threads;
        self
    }
}

/// What [`Pat2d::build`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BuildReport2 {
    /// Points.
    pub n: usize,
    /// Seconds in the neighborhoods.
    pub knn_seconds: f64,
    /// Seconds in the predictor.
    pub predict_seconds: f64,
    /// Fitting counters.
    pub fit: Fit2dReport,
    /// Wall time.
    pub seconds: f64,
}

/// The 2D facade: neighborhoods, predictor, circles, field.
#[derive(Debug, Clone, Copy)]
pub struct Pat2d;

impl Pat2d {
    /// Build the circles of `cloud` with `predictor` and return the evaluable field.
    pub fn build(cloud: &Cloud2, predictor: &dyn QuadraticPredictor, opts: &BuildOptions2) -> Result<(Pat2dField, BuildReport2), PatError> {
        let t0 = Instant::now();
        let k = predictor.neighborhood_size();
        if k == 0 {
            return Err(PatError::InvalidOption { what: "predictor neighborhood size 0".into() });
        }
        let threads = opts.fit.threads;
        let pool = make_pool(threads)?;
        let emb = cloud.embedded();
        let tree = pool.install(|| KdTree::new(&emb, KdTreeOptions::default()))?;
        let nb = pool.install(|| fitting_neighborhoods(&tree, &emb, k));
        drop(pool);
        let knn_seconds = t0.elapsed().as_secs_f64();
        let tp = Instant::now();
        let pred = predictor.predict(cloud, &nb, threads)?;
        let predict_seconds = tp.elapsed().as_secs_f64();
        let (circles, mut fit) = fit_circles(cloud, &pred.coefficients, &pred.sigma, &opts.fit)?;
        fit.degenerate_fits = pred.fallbacks;
        fit.tokens = nb.k;
        let field = Pat2dField::new(cloud, circles, opts.eval.clone())?;
        Ok((field, BuildReport2 { n: cloud.len(), knn_seconds, predict_seconds, fit, seconds: t0.elapsed().as_secs_f64() }))
    }
}

impl CirclePrimitive {
    /// Signed distance of this primitive at `x`.
    pub fn distance(&self, x: [f64; 2]) -> f64 {
        circle_sdf(self, x)
    }
}
