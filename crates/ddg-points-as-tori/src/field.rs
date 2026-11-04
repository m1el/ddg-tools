//! The evaluable field and the `PointsAsTori` facade (SPEC.md §7.4–§7.5).

use crate::blend::{BlendFlags, BlendInputs, BlendValue, BlendWant, Degenerate, GradientMode, LambdaRule, blend};
use crate::fit::{FitOptions, FitReport, Tori, fit_tori_with};
use crate::network::{Device, Model, PredictOptions, PredictReport, predict_coefficients};
use crate::spatial::{KdTree, KdTreeOptions, Neighbor, fitting_neighborhoods};
use crate::{PatError, PointCloud};
use rayon::prelude::*;

/// The neighbor set of a query.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Neighborhood {
    /// The `k` nearest data points (the reference: 32). `k ≥ n` means all points (in index order,
    /// as the reference).
    Knn {
        /// Neighbors.
        k: usize,
    },
    /// All data points, in index order (Eq. 26 exactly).
    All,
    /// All points within `r` (sorted by `(d², index)`); the `fallback_k` nearest when empty.
    Radius {
        /// Radius.
        r: f64,
        /// Fallback neighbors.
        fallback_k: usize,
    },
}

impl Default for Neighborhood {
    fn default() -> Self {
        Neighborhood::Knn { k: 32 }
    }
}

/// Evaluation precision. The CPU is f64 Eq. 23 always; `F32Stable` is a GPU option (WP4).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EvalPrecision {
    /// f64.
    #[default]
    F64,
    /// f32 with the stable form of SPEC.md §7.3 (GPU only).
    F32Stable,
}

/// Options of [`TorusField`] (defaults: the released reference, decision D13).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct EvalOptions {
    /// The neighbor set.
    pub neighborhood: Neighborhood,
    /// The λ rule.
    pub lambda: LambdaRule,
    /// Frozen (reference) or exact gradient.
    pub gradient: GradientMode,
    /// Optional per-point weights `A_i` (D28; default off).
    pub area_weights: Option<Vec<f64>>,
    /// Subtracted from `φ` (offsets).
    pub isovalue: f64,
    /// Precision.
    pub precision: EvalPrecision,
    /// Device (the CPU here; WP4 adds CUDA).
    pub device: Device,
    /// Threads for batch evaluation (`None`: all cores).
    pub threads: Option<usize>,
}

/// Output buffers of [`TorusField::evaluate`]; each requested buffer has one entry per query.
#[derive(Debug, Default)]
pub struct EvalOut<'a> {
    /// `φ − isovalue`.
    pub distance: Option<&'a mut [f64]>,
    /// `∇φ` (mode as in the options).
    pub gradient: Option<&'a mut [[f64; 3]]>,
    /// `Δφ` (frozen λ).
    pub laplacian: Option<&'a mut [f64]>,
}

/// What an evaluation did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct EvalReport {
    /// Queries evaluated.
    pub queries: usize,
    /// Torus terms evaluated on an axis (one-sided gradient rule).
    pub singular_axis: u64,
    /// Torus terms evaluated on a core circle (gradient 0).
    pub singular_core: u64,
    /// Queries coinciding with a data point (`r̂ := 0` for that term).
    pub at_data_point: u64,
    /// Radius queries that fell back to kNN.
    pub radius_fallbacks: u64,
    /// Neighbor terms summed.
    pub terms: u64,
    /// Wall time.
    pub seconds: f64,
}

/// An analytic SDF approximation: cloud, tori, kd-tree and evaluation options.
#[derive(Debug)]
pub struct TorusField {
    positions: Vec<[f64; 3]>,
    tori: Tori,
    tree: KdTree,
    options: EvalOptions,
    lambda_global: f64,
    mean_distance: f64,
    pool: rayon::ThreadPool,
}

fn degenerate(index: usize) -> PatError {
    PatError::DegenerateQuery {
        index,
        reason: "every neighbor is at the query, so σ_x = 0 and the PerQuery λ = C/σ_x is undefined".into(),
    }
}

/// Queries per parallel task (a constant; every query is independent anyway).
const EVAL_CHUNK: usize = 512;

/// Mean over all points of the mean distance to their `k` nearest points (self included), the
/// `D` of Sec. 4.2. Runs on the current pool.
pub fn mean_neighbor_distance(tree: &KdTree, positions: &[[f64; 3]], k: usize) -> f64 {
    let sums: Vec<f64> = positions
        .par_chunks(EVAL_CHUNK)
        .map_init(Vec::new, |buf, ps| {
            ps.iter()
                .map(|p| {
                    tree.knn(*p, k, buf);
                    buf.iter().map(|nb| nb.d2.sqrt()).sum::<f64>() / buf.len() as f64
                })
                .sum::<f64>()
        })
        .collect();
    sums.iter().sum::<f64>() / positions.len() as f64
}

impl TorusField {
    /// Build a field from a cloud, its tori and options (builds the kd-tree).
    pub fn new(cloud: &PointCloud, tori: Tori, options: EvalOptions) -> Result<Self, PatError> {
        let n = cloud.len();
        if tori.tori.len() != n {
            return Err(PatError::SizeMismatch { what: format!("{} tori for {n} points", tori.tori.len()) });
        }
        for (i, t) in tori.tori.iter().enumerate() {
            let ok = t.center.iter().chain(t.axis.iter()).all(|x| x.is_finite()) && t.major.is_finite() && t.minor.is_finite();
            if !ok {
                return Err(PatError::NonFiniteTorus { indices: vec![i as u32] });
            }
        }
        if let Some(a) = &options.area_weights {
            if a.len() != n {
                return Err(PatError::SizeMismatch { what: format!("{} area weights for {n} points", a.len()) });
            }
            if let Some(i) = a.iter().position(|w| !(w.is_finite() && *w > 0.0)) {
                return Err(PatError::InvalidOption { what: format!("area weight {i} = {} (must be finite and > 0)", a[i]) });
            }
        }
        match options.neighborhood {
            Neighborhood::Knn { k: 0 } => return Err(PatError::InvalidOption { what: "Neighborhood::Knn k = 0".into() }),
            Neighborhood::Radius { r, fallback_k } if !(r.is_finite() && r > 0.0) || fallback_k == 0 => {
                return Err(PatError::InvalidOption { what: format!("Neighborhood::Radius r = {r}, fallback_k = {fallback_k}") });
            }
            _ => {}
        }
        match options.lambda {
            LambdaRule::PerQuery { c } if !(c.is_finite() && c > 0.0) => {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::PerQuery c = {c}") });
            }
            LambdaRule::Constant { lambda } if !(lambda.is_finite() && lambda > 0.0) => {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::Constant lambda = {lambda}") });
            }
            LambdaRule::PaperHeuristic { c, factor, k_mean }
                if !(c.is_finite() && c > 0.0 && factor.is_finite() && factor > 0.0 && k_mean > 0) =>
            {
                return Err(PatError::InvalidOption { what: format!("LambdaRule::PaperHeuristic {:?}", options.lambda) });
            }
            _ => {}
        }
        if !options.isovalue.is_finite() {
            return Err(PatError::InvalidOption { what: format!("isovalue = {}", options.isovalue) });
        }
        if options.precision != EvalPrecision::F64 || !matches!(options.device, Device::Cpu | Device::Auto) {
            return Err(PatError::Gpu { reason: "GPU evaluation (and F32Stable) is not built yet (WP4)".into() });
        }
        let pool = crate::pool::make_pool(options.threads)?;
        let positions = cloud.positions().to_vec();
        let tree = pool.install(|| KdTree::new(&positions, KdTreeOptions::default()))?;
        let (mut lambda_global, mut mean_distance) = (f64::NAN, f64::NAN);
        if let LambdaRule::PaperHeuristic { factor, k_mean, .. } = options.lambda {
            mean_distance = pool.install(|| mean_neighbor_distance(&tree, &positions, k_mean));
            lambda_global = factor / mean_distance;
            if !(lambda_global.is_finite() && lambda_global > 0.0) {
                return Err(PatError::InvalidOption { what: format!("PaperHeuristic gives λ = {lambda_global} (D = {mean_distance})") });
            }
        }
        Ok(TorusField { positions, tori, tree, options, lambda_global, mean_distance, pool })
    }

    /// The options of the paper's Sec. 4.2 for `cloud`: `λ = 10³/D` constant, radius
    /// `R_eval = 2C/λ` with `C = 64`, 32-NN fallback. Returns the options and `R_eval`.
    pub fn paper_options(cloud: &PointCloud, threads: Option<usize>) -> Result<(EvalOptions, f64), PatError> {
        let (c, factor, k_mean) = (64.0, 1e3, 64);
        let pool = crate::pool::make_pool(threads)?;
        let tree = pool.install(|| KdTree::new(cloud.positions(), KdTreeOptions::default()))?;
        let d = pool.install(|| mean_neighbor_distance(&tree, cloud.positions(), k_mean));
        let r_eval = 2.0 * c / (factor / d);
        let opts = EvalOptions {
            neighborhood: Neighborhood::Radius { r: r_eval, fallback_k: 32 },
            lambda: LambdaRule::PaperHeuristic { c, factor, k_mean },
            threads,
            ..Default::default()
        };
        Ok((opts, r_eval))
    }

    /// The tori.
    pub fn tori(&self) -> &Tori {
        &self.tori
    }

    /// The data positions.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// The options.
    pub fn options(&self) -> &EvalOptions {
        &self.options
    }

    /// The kd-tree over the data points.
    pub fn tree(&self) -> &KdTree {
        &self.tree
    }

    /// For [`LambdaRule::PaperHeuristic`]: `(D, λ)`; NaN otherwise.
    pub fn paper_constants(&self) -> (f64, f64) {
        (self.mean_distance, self.lambda_global)
    }

    fn inputs(&self) -> BlendInputs<'_> {
        BlendInputs {
            positions: &self.positions,
            tori: &self.tori.tori,
            area_weights: self.options.area_weights.as_deref(),
            lambda: self.options.lambda,
            lambda_global: self.lambda_global,
            isovalue: self.options.isovalue,
        }
    }

    /// The neighbor set of `x` (into `idx`, using `buf`); returns whether a radius query fell back.
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

    fn eval_one(
        &self,
        x: [f64; 3],
        want: BlendWant,
        buf: &mut Vec<Neighbor>,
        idx: &mut Vec<u32>,
    ) -> Result<(BlendValue, bool, usize), Degenerate> {
        let fb = self.neighbors(x, buf, idx);
        Ok((blend(&self.inputs(), x, idx, want)?, fb, idx.len()))
    }

    fn single(&self, x: [f64; 3], want: BlendWant) -> Result<BlendValue, PatError> {
        if !x.iter().all(|c| c.is_finite()) {
            return Err(PatError::NonFinite { index: 0, what: "query" });
        }
        self.eval_one(x, want, &mut Vec::new(), &mut Vec::new()).map(|v| v.0).map_err(|_| degenerate(0))
    }

    /// `φ(x) − isovalue`.
    pub fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok(self.single(x, BlendWant::Distance)?.distance)
    }

    /// `φ(x) − isovalue` and `∇φ(x)` (gradient mode of the options).
    pub fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let v = self.single(x, BlendWant::Gradient(self.options.gradient))?;
        Ok((v.distance, v.gradient))
    }

    /// `Δφ(x)` with frozen λ (App. D divergence with our `Δg`).
    pub fn laplacian(&self, x: [f64; 3]) -> Result<f64, PatError> {
        Ok(self.single(x, BlendWant::Laplacian)?.laplacian)
    }

    /// Everything one blend computes at `x` (value, gradient or Laplacian per `want`, `σ_x`, `λ`,
    /// the shift, counters). For tests and the oracle replay.
    pub fn blend_at(&self, x: [f64; 3], want: BlendWant) -> Result<(BlendValue, Vec<u32>), PatError> {
        if !x.iter().all(|c| c.is_finite()) {
            return Err(PatError::NonFinite { index: 0, what: "query" });
        }
        let mut idx = Vec::new();
        let (v, _, _) = self.eval_one(x, want, &mut Vec::new(), &mut idx).map_err(|_| degenerate(0))?;
        Ok((v, idx))
    }

    /// Batch evaluation on the field's pool. The requested buffers must have `xs.len()` entries.
    pub fn evaluate(&self, xs: &[[f64; 3]], out: EvalOut<'_>) -> Result<EvalReport, PatError> {
        let t0 = std::time::Instant::now();
        let q = xs.len();
        for (name, len) in [
            ("distance", out.distance.as_ref().map(|b| b.len())),
            ("gradient", out.gradient.as_ref().map(|b| b.len())),
            ("laplacian", out.laplacian.as_ref().map(|b| b.len())),
        ] {
            if let Some(l) = len
                && l != q
            {
                return Err(PatError::SizeMismatch { what: format!("{name} buffer has {l} entries for {q} queries") });
            }
        }
        if let Some(i) = xs.iter().position(|x| !x.iter().all(|c| c.is_finite())) {
            return Err(PatError::NonFinite { index: i, what: "query" });
        }
        let want = if out.laplacian.is_some() {
            BlendWant::Laplacian
        } else if out.gradient.is_some() {
            BlendWant::Gradient(self.options.gradient)
        } else {
            BlendWant::Distance
        };
        if out.laplacian.is_some() && out.gradient.is_some() && self.options.gradient == GradientMode::Exact {
            return Err(PatError::InvalidOption {
                what: "the Laplacian pass returns the frozen gradient; request the exact gradient separately".into(),
            });
        }
        let vals: Vec<Result<(BlendValue, bool, usize), Degenerate>> = self.pool.install(|| {
            xs.par_chunks(EVAL_CHUNK)
                .flat_map_iter(|c| {
                    let mut buf = Vec::new();
                    let mut idx = Vec::new();
                    c.iter().map(|x| self.eval_one(*x, want, &mut buf, &mut idx)).collect::<Vec<_>>()
                })
                .collect()
        });
        if let Some(i) = vals.iter().position(|v| v.is_err()) {
            return Err(degenerate(i));
        }
        let mut rep = EvalReport { queries: q, ..Default::default() };
        let EvalOut { mut distance, mut gradient, mut laplacian } = out;
        for (i, (v, fb, terms)) in vals.iter().flatten().enumerate() {
            let BlendFlags { singular_axis, singular_core, at_data_point } = v.flags;
            rep.singular_axis += singular_axis as u64;
            rep.singular_core += singular_core as u64;
            rep.at_data_point += at_data_point as u64;
            rep.radius_fallbacks += *fb as u64;
            rep.terms += *terms as u64;
            if let Some(d) = distance.as_deref_mut() {
                d[i] = v.distance;
            }
            if let Some(g) = gradient.as_deref_mut() {
                g[i] = v.gradient;
            }
            if let Some(l) = laplacian.as_deref_mut() {
                l[i] = v.laplacian;
            }
        }
        rep.seconds = t0.elapsed().as_secs_f64();
        Ok(rep)
    }
}

/// Options of [`PointsAsTori::build`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BuildOptions {
    /// Neighborhood size for the network (`None`: the model's `k`).
    pub k: Option<usize>,
    /// Prediction.
    pub predict: PredictOptions,
    /// Fitting.
    pub fit: FitOptions,
    /// Evaluation.
    pub eval: EvalOptions,
}

impl BuildOptions {
    /// Set the thread count of every stage.
    pub fn with_threads(mut self, threads: Option<usize>) -> Self {
        self.predict.threads = threads;
        self.fit.threads = threads;
        self.eval.threads = threads;
        self
    }
}

/// What [`PointsAsTori::build`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BuildReport {
    /// Points.
    pub n: usize,
    /// Tokens per neighborhood.
    pub tokens: usize,
    /// Seconds for the kd-tree and the fitting neighborhoods.
    pub knn_seconds: f64,
    /// Prediction.
    pub predict: PredictReport,
    /// Fitting.
    pub fit: FitReport,
    /// Wall time.
    pub seconds: f64,
}

/// The facade: kNN, network, fitting, field.
#[derive(Debug, Clone, Copy)]
pub struct PointsAsTori;

impl PointsAsTori {
    /// Build the tori of `cloud` with `model` and return the evaluable field.
    pub fn build(cloud: &PointCloud, model: &Model, opts: &BuildOptions) -> Result<(TorusField, BuildReport), PatError> {
        let t0 = std::time::Instant::now();
        let k = opts.k.unwrap_or(model.config().k);
        if k == 0 {
            return Err(PatError::InvalidOption { what: "BuildOptions::k = 0".into() });
        }
        let pool = crate::pool::make_pool(opts.predict.threads)?;
        let tk = std::time::Instant::now();
        let nb = pool.install(|| -> Result<_, PatError> {
            let tree = KdTree::new(cloud.positions(), KdTreeOptions::default())?;
            Ok(fitting_neighborhoods(&tree, cloud.positions(), k))
        })?;
        drop(pool);
        let knn_seconds = tk.elapsed().as_secs_f64();
        let (coeffs, sigma, predict) = predict_coefficients(model, cloud, &nb, &opts.predict)?;
        let tokens_lt_k = if nb.k < model.config().k { cloud.len() } else { 0 };
        let (tori, fit) = fit_tori_with(cloud, &coeffs, &sigma, tokens_lt_k, &opts.fit)?;
        let field = TorusField::new(cloud, tori, opts.eval.clone())?;
        let report = BuildReport { n: cloud.len(), tokens: nb.k, knn_seconds, predict, fit, seconds: t0.elapsed().as_secs_f64() };
        Ok((field, report))
    }
}
