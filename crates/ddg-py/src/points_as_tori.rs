//! Points as Tori (Feng, Gkioulekas, Crane, ACM TOG 45(4), 2026) for the Python module
//! `ddg_tools`: signed distance to an oriented point cloud by blending fitted tori, with
//! marching cubes, CSG, offsets, morphology, level-set sampling, normal orientation, the baselines
//! of the paper's evaluation, the exact mesh signed distance and the classical 2D variant.
//!
//! This file is plain Rust. Options arrive as a key/value list ([`Kw`], filled by the PyO3 glue
//! from the keyword arguments), arrays as flat `f64` vectors; reports go back as [`V`] trees.
//! Every numerical routine is the crate's (`ddg-points-as-tori`); nothing is computed here.
//! The Python-facing documentation is the crate README and `ddg_tools.pyi`.

use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex, OnceLock};

use ddg_points_as_tori as pat;
use pat::apps::orient::{Epsilon, OrientOptions, orient_normals};
use pat::apps::{Morph, SampleOptions, sample_level_set};
use pat::baselines::{FastWinding, FwnOptions, KernelField, KernelMethod, KernelOptions, PointAreaOptions, point_areas};
use pat::csg::{Csg, ImplicitField};
use pat::mc::{Band, McOptions, Resolution, marching_cubes};
use pat::pat2d;
use pat::truth::{MeshTruth, SignMode, TruthOptions};
use pat::fit::Torus;
use pat::{
    BuildOptions, BuildReport, Device, EigenRule, EvalOptions, EvalPrecision, EvalReport, FeatureOrder, GpuMatmul, GradientMode, LambdaRule, Model,
    Neighborhood, NormalPolicy, PatError, PointCloud, PointsAsTori, Progress, SignRule, Tori, TorusField,
};

// -------------------------------------------------------------------------------------------
// Values and errors
// -------------------------------------------------------------------------------------------

/// A report value (becomes a Python `None`, `bool`, `int`, `float`, `str`, `list` or `dict`).
#[derive(Debug, Clone, PartialEq)]
pub enum V {
    /// `None`.
    None,
    /// `bool`.
    B(bool),
    /// `int`.
    I(i64),
    /// `float`.
    F(f64),
    /// `str`.
    S(String),
    /// `list`.
    L(Vec<V>),
    /// `dict` (insertion order).
    D(Vec<(&'static str, V)>),
}

impl V {
    fn u(x: usize) -> V {
        V::I(x as i64)
    }
    fn s(x: &str) -> V {
        V::S(x.to_owned())
    }
    fn opt_f(x: Option<f64>) -> V {
        x.map_or(V::None, V::F)
    }
    fn opt_u(x: Option<usize>) -> V {
        x.map_or(V::None, V::u)
    }
    fn arr3(x: [f64; 3]) -> V {
        V::L(x.iter().map(|&c| V::F(c)).collect())
    }
}

/// An option value (a Python keyword argument).
#[derive(Debug, Clone, PartialEq)]
pub enum O {
    /// `None`.
    None,
    /// `bool`.
    Bool(bool),
    /// `int`.
    Int(i64),
    /// `float`.
    Float(f64),
    /// `str`.
    Str(String),
    /// A flat numeric buffer or sequence.
    Array(Vec<f64>),
    /// A progress/cancel handle (`PatProgress`).
    Progress(Progress),
}

/// Why a call failed. Both kinds are `PointsAsToriError` (a `ValueError`) in Python.
#[derive(Debug, Clone, PartialEq)]
pub struct PErr(pub String);

impl From<PatError> for PErr {
    fn from(e: PatError) -> Self {
        PErr(e.to_string())
    }
}

/// Result alias.
pub type R<T> = Result<T, PErr>;

fn bad<T>(msg: impl Into<String>) -> R<T> {
    Err(PErr(msg.into()))
}

/// Keyword options with use tracking: [`Kw::finish`] rejects the keys nobody read.
#[derive(Debug, Default)]
pub struct Kw {
    items: Vec<(String, O, bool)>,
}

impl Kw {
    /// From key/value pairs.
    pub fn new(items: Vec<(String, O)>) -> Kw {
        Kw { items: items.into_iter().map(|(k, v)| (k, v, false)).collect() }
    }
    fn take(&mut self, key: &str) -> Option<&O> {
        for (k, v, used) in self.items.iter_mut() {
            if k == key {
                *used = true;
                return Some(v);
            }
        }
        None
    }
    /// A string option (`None`/absent: `None`).
    pub fn opt_str(&mut self, key: &str) -> R<Option<String>> {
        match self.take(key) {
            None | Some(O::None) => Ok(None),
            Some(O::Str(s)) => Ok(Some(s.clone())),
            Some(o) => bad(format!("{key}: expected a string, got {o:?}")),
        }
    }
    /// A string option restricted to `allowed`.
    pub fn choice(&mut self, key: &str, allowed: &[&str]) -> R<Option<String>> {
        match self.opt_str(key)? {
            None => Ok(None),
            Some(s) if allowed.contains(&s.as_str()) => Ok(Some(s)),
            Some(s) => bad(format!("{key}: {s:?} is not one of {}", allowed.iter().map(|a| format!("{a:?}")).collect::<Vec<_>>().join(", "))),
        }
    }
    /// A float option.
    pub fn opt_f64(&mut self, key: &str) -> R<Option<f64>> {
        match self.take(key) {
            None | Some(O::None) => Ok(None),
            Some(O::Float(x)) => Ok(Some(*x)),
            Some(O::Int(i)) => Ok(Some(*i as f64)),
            Some(o) => bad(format!("{key}: expected a number, got {o:?}")),
        }
    }
    /// A non-negative integer option.
    pub fn opt_usize(&mut self, key: &str) -> R<Option<usize>> {
        match self.take(key) {
            None | Some(O::None) => Ok(None),
            Some(O::Int(i)) if *i >= 0 => Ok(Some(*i as usize)),
            Some(o) => bad(format!("{key}: expected a non-negative integer, got {o:?}")),
        }
    }
    /// A seed.
    pub fn opt_u64(&mut self, key: &str) -> R<Option<u64>> {
        self.opt_usize(key).map(|x| x.map(|v| v as u64))
    }
    /// A progress handle.
    pub fn opt_progress(&mut self, key: &str) -> R<Option<Progress>> {
        match self.take(key) {
            None | Some(O::None) => Ok(None),
            Some(O::Progress(p)) => Ok(Some(p.clone())),
            Some(o) => bad(format!("{key}: expected a PatProgress, got {o:?}")),
        }
    }
    /// A boolean option.
    pub fn opt_bool(&mut self, key: &str) -> R<Option<bool>> {
        match self.take(key) {
            None | Some(O::None) => Ok(None),
            Some(O::Bool(b)) => Ok(Some(*b)),
            Some(o) => bad(format!("{key}: expected a bool, got {o:?}")),
        }
    }
    /// An array option.
    pub fn opt_array(&mut self, key: &str) -> R<Option<Vec<f64>>> {
        match self.take(key) {
            None | Some(O::None) => Ok(None),
            Some(O::Array(a)) => Ok(Some(a.clone())),
            Some(O::Float(x)) => Ok(Some(vec![*x])),
            Some(O::Int(i)) => Ok(Some(vec![*i as f64])),
            Some(o) => bad(format!("{key}: expected an array of numbers, got {o:?}")),
        }
    }
    /// Reject the options that were not read.
    pub fn finish(&self, what: &str) -> R<()> {
        let unknown: Vec<&str> = self.items.iter().filter(|(_, _, u)| !u).map(|(k, _, _)| k.as_str()).collect();
        if unknown.is_empty() { Ok(()) } else { bad(format!("{what}: unknown option(s) {}", unknown.join(", "))) }
    }
}

fn to3(name: &str, v: &[f64]) -> R<Vec<[f64; 3]>> {
    if !v.len().is_multiple_of(3) {
        return bad(format!("{name}: {} numbers are not a whole number of 3-vectors", v.len()));
    }
    Ok(v.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

fn to2(name: &str, v: &[f64]) -> R<Vec<[f64; 2]>> {
    if !v.len().is_multiple_of(2) {
        return bad(format!("{name}: {} numbers are not a whole number of 2-vectors", v.len()));
    }
    Ok(v.chunks_exact(2).map(|c| [c[0], c[1]]).collect())
}

fn flat3(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

fn flat2(v: &[[f64; 2]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

// -------------------------------------------------------------------------------------------
// Options
// -------------------------------------------------------------------------------------------

/// How the point-area weights of the blend (and of the winding baselines) are given.
#[derive(Debug, Clone, PartialEq)]
pub enum AreaSpec {
    /// None (all weights 1).
    None,
    /// Given, one per point.
    Given(Vec<f64>),
    /// Barill's point areas of the cloud.
    Barill,
}

/// Everything the constructor needs besides the cloud.
#[derive(Debug, Clone)]
pub struct BuildSpec {
    /// Normal validation.
    pub policy: NormalPolicy,
    /// Weights file (`None`: the bundled ones).
    pub weights: Option<String>,
    /// Build options.
    pub opts: BuildOptions,
    /// Use the paper's evaluation options as the base.
    pub paper: bool,
    /// Area weights of the blend.
    pub areas: AreaSpec,
    /// Explicit evaluation overrides, applied after `paper`.
    pub eval: EvalOverrides,
}

/// Evaluation options that were given explicitly.
#[derive(Debug, Clone, Default)]
pub struct EvalOverrides {
    neighborhood: Option<Neighborhood>,
    lambda: Option<LambdaRule>,
    gradient: Option<GradientMode>,
    isovalue: Option<f64>,
    precision: Option<EvalPrecision>,
}

fn parse_device(kw: &mut Kw) -> R<Device> {
    let ord = kw.opt_usize("cuda_ordinal")?.unwrap_or(0);
    Ok(match kw.choice("device", &["cpu", "cuda", "auto"])?.as_deref() {
        None | Some("cpu") => Device::Cpu,
        Some("cuda") => Device::Cuda { ordinal: ord },
        _ => Device::Auto,
    })
}

fn parse_threads(kw: &mut Kw) -> R<Option<usize>> {
    match kw.opt_usize("threads")? {
        Some(0) => bad("threads: must be at least 1 (or None for all cores)"),
        t => Ok(t),
    }
}

fn parse_eval_overrides(kw: &mut Kw) -> R<EvalOverrides> {
    let mut o = EvalOverrides::default();
    let knn_k = kw.opt_usize("knn_k")?;
    let radius = kw.opt_f64("radius")?;
    let fb = kw.opt_usize("radius_fallback_k")?;
    o.neighborhood = match kw.choice("neighborhood", &["knn", "all", "radius"])?.as_deref() {
        None if radius.is_some() => Some(Neighborhood::Radius { r: radius.unwrap_or(0.0), fallback_k: fb.unwrap_or(32) }),
        None if knn_k.is_some() => Some(Neighborhood::Knn { k: knn_k.unwrap_or(32) }),
        None => None,
        Some("knn") => Some(Neighborhood::Knn { k: knn_k.unwrap_or(32) }),
        Some("all") => Some(Neighborhood::All),
        Some(_) => match radius {
            Some(r) => Some(Neighborhood::Radius { r, fallback_k: fb.unwrap_or(32) }),
            None => return bad("radius: neighborhood=\"radius\" needs radius=<float> (or paper=True)"),
        },
    };
    let c = kw.opt_f64("lambda_c")?;
    let lambda_value = kw.opt_f64("lambda_value")?;
    let factor = kw.opt_f64("paper_factor")?;
    let k_mean = kw.opt_usize("paper_k_mean")?;
    o.lambda = match kw.choice("lambda_rule", &["per_query", "constant", "paper"])?.as_deref() {
        None if c.is_none() && lambda_value.is_none() && factor.is_none() && k_mean.is_none() => None,
        None | Some("per_query") => Some(LambdaRule::PerQuery { c: c.unwrap_or(64.0) }),
        Some("constant") => match lambda_value {
            Some(lambda) => Some(LambdaRule::Constant { lambda }),
            None => return bad("lambda_value: lambda_rule=\"constant\" needs lambda_value=<float>"),
        },
        Some(_) => Some(LambdaRule::PaperHeuristic { c: c.unwrap_or(64.0), factor: factor.unwrap_or(1e3), k_mean: k_mean.unwrap_or(64) }),
    };
    o.gradient = kw.choice("gradient", &["frozen", "exact"])?.map(|g| if g == "exact" { GradientMode::Exact } else { GradientMode::Frozen });
    o.isovalue = kw.opt_f64("isovalue")?;
    o.precision = kw.choice("precision", &["f64", "f32_stable"])?.map(|p| if p == "f64" { EvalPrecision::F64 } else { EvalPrecision::F32Stable });
    Ok(o)
}

/// Parse the constructor options (consumes the keys it knows).
pub fn parse_build(kw: &mut Kw, paper_default: bool) -> R<BuildSpec> {
    let policy = match kw.choice("normal_policy", &["require_unit", "normalize"])?.as_deref() {
        Some("normalize") => NormalPolicy::Normalize,
        _ => NormalPolicy::RequireUnit { tol: kw.opt_f64("normal_tol")?.unwrap_or(1e-6) },
    };
    let _ = kw.opt_f64("normal_tol")?;
    let weights = kw.opt_str("weights")?;
    let device = parse_device(kw)?;
    let threads = parse_threads(kw)?;
    let matmul = match kw.choice("gpu_matmul", &["tf32", "split3"])?.as_deref() {
        Some("split3") => GpuMatmul::Split3,
        _ => GpuMatmul::Tf32,
    };
    let feature_order = match kw.choice("feature_order", &["difference_first", "reference"])?.as_deref() {
        Some("reference") => FeatureOrder::Reference,
        _ => FeatureOrder::DifferenceFirst,
    };
    let sign_rule = match kw.choice("sign_rule", &["reference", "geometric"])?.as_deref() {
        Some("geometric") => SignRule::Geometric,
        _ => SignRule::Reference,
    };
    let eigen_rule = match kw.choice("eigen_rule", &["reference", "robust"])?.as_deref() {
        Some("robust") => EigenRule::Robust,
        _ => EigenRule::Reference,
    };
    let k = kw.opt_usize("k")?;
    let paper = kw.opt_bool("paper")?.unwrap_or(paper_default);
    let areas = match kw.take("area_weights") {
        None | Some(O::None) => AreaSpec::None,
        Some(O::Str(s)) if s == "barill" => AreaSpec::Barill,
        Some(O::Array(a)) => AreaSpec::Given(a.clone()),
        Some(o) => return bad(format!("area_weights: expected None, \"barill\" or one number per point, got {o:?}")),
    };
    let eval = parse_eval_overrides(kw)?;
    let progress = kw.opt_progress("progress")?;
    let mut opts = BuildOptions { k, progress: progress.clone(), ..Default::default() };
    opts.eval.progress = progress;
    opts.predict.device = device;
    opts.predict.threads = threads;
    opts.predict.feature_order = feature_order;
    opts.predict.gpu_matmul = matmul;
    opts.fit.sign_rule = sign_rule;
    opts.fit.eigen_rule = eigen_rule;
    opts.fit.threads = threads;
    opts.eval.device = device;
    opts.eval.threads = threads;
    Ok(BuildSpec { policy, weights, opts, paper, areas, eval })
}

static REFERENCE_MODEL: OnceLock<Result<Arc<Model>, String>> = OnceLock::new();
static CUSTOM_MODELS: OnceLock<Mutex<HashMap<String, Arc<Model>>>> = OnceLock::new();

/// The bundled network, or a safetensors file (cached by path).
pub fn model(weights: Option<&str>) -> R<Arc<Model>> {
    match weights {
        None => REFERENCE_MODEL.get_or_init(|| Model::reference().map(Arc::new).map_err(|e| e.to_string())).clone().map_err(PErr),
        Some(path) => {
            let cache = CUSTOM_MODELS.get_or_init(Default::default);
            if let Some(m) = cache.lock().map_err(|_| PErr("model cache poisoned".into()))?.get(path) {
                return Ok(m.clone());
            }
            let bytes = std::fs::read(path).map_err(|e| PErr(format!("weights: cannot read {path}: {e}")))?;
            let m = Arc::new(Model::from_safetensors(&bytes)?);
            cache.lock().map_err(|_| PErr("model cache poisoned".into()))?.insert(path.to_owned(), m.clone());
            Ok(m)
        }
    }
}

fn apply_eval(spec: &BuildSpec, cloud: &PointCloud, base: &mut EvalOptions) -> R<()> {
    if spec.paper {
        let (paper, _) = TorusField::paper_options(cloud, spec.opts.eval.threads)?;
        base.neighborhood = paper.neighborhood;
        base.lambda = paper.lambda;
    }
    let o = &spec.eval;
    if let Some(v) = o.neighborhood {
        base.neighborhood = v;
    }
    if let Some(v) = o.lambda {
        base.lambda = v;
    }
    if let Some(v) = o.gradient {
        base.gradient = v;
    }
    if let Some(v) = o.isovalue {
        base.isovalue = v;
    }
    if let Some(v) = o.precision {
        base.precision = v;
    }
    base.area_weights = match &spec.areas {
        AreaSpec::None => None,
        AreaSpec::Given(a) => Some(a.clone()),
        AreaSpec::Barill => Some(point_areas(cloud, &PointAreaOptions { threads: spec.opts.eval.threads, ..Default::default() })?.0),
    };
    Ok(())
}

// -------------------------------------------------------------------------------------------
// The field
// -------------------------------------------------------------------------------------------

/// A built field with its cloud and the reports of how it came to be.
pub struct Pat {
    /// The validated cloud.
    pub cloud: PointCloud,
    /// The field (shared with CSG trees).
    pub field: Arc<TorusField>,
    /// What building did (absent for fields loaded from tori).
    pub build_report: Option<BuildReport>,
    /// What the cloud validation found.
    pub cloud_report: pat::CloudReport,
    /// The report of the last evaluation.
    pub last_eval: Mutex<Option<EvalReport>>,
}

impl Pat {
    /// Build from flat arrays.
    pub fn build(points: &[f64], normals: &[f64], spec: &BuildSpec) -> R<Pat> {
        let (cloud, creport) = PointCloud::new(to3("points", points)?, to3("normals", normals)?, spec.policy)?;
        let m = model(spec.weights.as_deref())?;
        let mut opts = spec.opts.clone();
        apply_eval(spec, &cloud, &mut opts.eval)?;
        let (field, report) = PointsAsTori::build(&cloud, &m, &opts)?;
        Ok(Pat { cloud, field: Arc::new(field), build_report: Some(report), cloud_report: creport, last_eval: Mutex::new(None) })
    }

    /// A field from saved tori.
    pub fn from_tori(points: &[f64], normals: &[f64], tori: Tori, sign_rule: Option<SignRule>, spec: &BuildSpec) -> R<Pat> {
        let (cloud, creport) = PointCloud::new(to3("points", points)?, to3("normals", normals)?, spec.policy)?;
        let mut eval = EvalOptions { device: spec.opts.eval.device, threads: spec.opts.eval.threads, progress: spec.opts.eval.progress.clone(), ..Default::default() };
        apply_eval(spec, &cloud, &mut eval)?;
        let field = TorusField::new(&cloud, tori, eval)?.with_sign_rule(sign_rule);
        Ok(Pat { cloud, field: Arc::new(field), build_report: None, cloud_report: creport, last_eval: Mutex::new(None) })
    }

    /// The same tori and cloud with other evaluation options.
    pub fn with_options(&self, kw: &mut Kw) -> R<Pat> {
        let mut spec = parse_build(kw, false)?;
        let cur = self.field.options();
        // keep the current threads and device unless given
        if !kw_given(kw) {
            spec.opts.eval.device = cur.device;
            spec.opts.eval.threads = cur.threads;
        }
        let mut eval = cur.clone();
        eval.device = spec.opts.eval.device;
        eval.threads = spec.opts.eval.threads;
        // the progress handle is never inherited: `progress=None` (the default) detaches it
        eval.progress = spec.opts.eval.progress.clone();
        if matches!(spec.areas, AreaSpec::None) {
            spec.areas = match &cur.area_weights {
                Some(a) => AreaSpec::Given(a.clone()),
                None => AreaSpec::None,
            };
        }
        apply_eval(&spec, &self.cloud, &mut eval)?;
        let field = TorusField::new(&self.cloud, self.field.tori().clone(), eval)?.with_sign_rule(self.field.sign_rule());
        Ok(Pat { cloud: self.cloud.clone(), field: Arc::new(field), build_report: self.build_report.clone(), cloud_report: self.cloud_report.clone(), last_eval: Mutex::new(None) })
    }

    /// Batch evaluation. Returns `(distance, gradient, laplacian)` as requested.
    pub fn evaluate(&self, queries: &[f64], d: bool, g: bool, l: bool) -> R<(Vec<f64>, Vec<f64>, Vec<f64>)> {
        let xs = to3("queries", queries)?;
        let n = xs.len();
        let mut dist = vec![0.0; if d { n } else { 0 }];
        let mut grad = vec![[0.0; 3]; if g { n } else { 0 }];
        let mut lap = vec![0.0; if l { n } else { 0 }];
        let out = pat::field::EvalOut {
            distance: d.then_some(&mut dist[..]),
            gradient: g.then_some(&mut grad[..]),
            laplacian: l.then_some(&mut lap[..]),
        };
        let rep = self.field.evaluate(&xs, out)?;
        if let Ok(mut s) = self.last_eval.lock() {
            *s = Some(rep);
        }
        Ok((dist, flat3(&grad), lap))
    }

    /// The tori in our layout and the reference's.
    pub fn tori(&self) -> &Tori {
        self.field.tori()
    }

    /// Save our binary tori file (tied to the cloud by its hash, and marked with the sign rule the
    /// tori were fitted with when that is known).
    pub fn save_tori(&self, path: &str) -> R<()> {
        std::fs::write(path, self.field.to_bytes()).map_err(|e| PErr(format!("cannot write {path}: {e}")))
    }

    /// The evaluation options as a report.
    pub fn options(&self) -> V {
        eval_options_v(self.field.options())
    }
}

fn kw_given(_kw: &Kw) -> bool {
    false
}

/// Read a tori file (checked against the cloud's hash).
pub fn load_tori_file(path: &str) -> R<(Tori, Option<SignRule>)> {
    let bytes = std::fs::read(path).map_err(|e| PErr(format!("cannot read {path}: {e}")))?;
    let (t, _, rule) = Tori::from_bytes_with(&bytes, None)?;
    Ok((t, rule))
}

/// Check a tori file against the cloud that is meant to go with it, then read it.
pub fn load_tori_for(path: &str, points: &[f64], normals: &[f64], policy: NormalPolicy) -> R<(Tori, Option<SignRule>)> {
    let bytes = std::fs::read(path).map_err(|e| PErr(format!("cannot read {path}: {e}")))?;
    let (cloud, _) = PointCloud::new(to3("points", points)?, to3("normals", normals)?, policy)?;
    let (t, _, rule) = Tori::from_bytes_with(&bytes, Some(cloud.hash()))?;
    Ok((t, rule))
}

/// The tori as a dict of flat arrays: the reference's layout (`centers`, `axes` as `(theta, phi)`,
/// `major`, `minor`) plus ours (`axis_vectors`, `sigma`, `coefficients`).
pub struct ToriArrays {
    /// `3n`.
    pub centers: Vec<f64>,
    /// `2n` angles.
    pub axes: Vec<f64>,
    /// `n`.
    pub major: Vec<f64>,
    /// `n` (signed).
    pub minor: Vec<f64>,
    /// `3n` unit axes.
    pub axis_vectors: Vec<f64>,
    /// `n`.
    pub sigma: Vec<f64>,
    /// `6n` (f32 values widened).
    pub coefficients: Vec<f64>,
}

/// Export.
pub fn export_tori(t: &Tori) -> ToriArrays {
    let mut a = ToriArrays { centers: vec![], axes: vec![], major: vec![], minor: vec![], axis_vectors: vec![], sigma: t.sigma.clone(), coefficients: vec![] };
    for x in &t.tori {
        a.centers.extend(x.center);
        a.axes.extend(x.axis_angles());
        a.major.push(x.major);
        a.minor.push(x.minor);
        a.axis_vectors.extend(x.axis);
    }
    for c in &t.coefficients {
        a.coefficients.extend(c.iter().map(|&v| f64::from(v)));
    }
    a
}

/// Import the reference's layout `(centers (n,3), axes (n,2) angles, R (n), r (n))`.
pub fn import_tori(centers: &[f64], axes: &[f64], major: &[f64], minor: &[f64]) -> R<Tori> {
    let c = to3("centers", centers)?;
    let a = to2("axes", axes)?;
    let n = c.len();
    if a.len() != n || major.len() != n || minor.len() != n {
        return bad(format!("tori: centers has {n} rows, axes {}, major {}, minor {}", a.len(), major.len(), minor.len()));
    }
    let tori = (0..n).map(|i| Torus { center: c[i], axis: Torus::axis_from_angles(a[i][0], a[i][1]), major: major[i], minor: minor[i] }).collect();
    Ok(Tori { tori, coefficients: vec![[f32::NAN; 6]; n], sigma: vec![f64::NAN; n] })
}

// -------------------------------------------------------------------------------------------
// Reports
// -------------------------------------------------------------------------------------------

fn dev_v(d: Device) -> V {
    match d {
        Device::Cpu => V::s("cpu"),
        Device::Cuda { ordinal } => V::S(format!("cuda:{ordinal}")),
        Device::Auto => V::s("auto"),
    }
}

/// `EvalOptions` as a dict.
pub fn eval_options_v(o: &EvalOptions) -> V {
    let (nb, knn_k, radius, fb) = match o.neighborhood {
        Neighborhood::Knn { k } => ("knn", V::u(k), V::None, V::None),
        Neighborhood::All => ("all", V::None, V::None, V::None),
        Neighborhood::Radius { r, fallback_k } => ("radius", V::None, V::F(r), V::u(fallback_k)),
    };
    let (rule, c, lv, factor, kmean) = match o.lambda {
        LambdaRule::PerQuery { c } => ("per_query", V::F(c), V::None, V::None, V::None),
        LambdaRule::Constant { lambda } => ("constant", V::None, V::F(lambda), V::None, V::None),
        LambdaRule::PaperHeuristic { c, factor, k_mean } => ("paper", V::F(c), V::None, V::F(factor), V::u(k_mean)),
    };
    V::D(vec![
        ("neighborhood", V::s(nb)),
        ("knn_k", knn_k),
        ("radius", radius),
        ("radius_fallback_k", fb),
        ("lambda_rule", V::s(rule)),
        ("lambda_c", c),
        ("lambda_value", lv),
        ("paper_factor", factor),
        ("paper_k_mean", kmean),
        ("gradient", V::s(if o.gradient == GradientMode::Exact { "exact" } else { "frozen" })),
        ("isovalue", V::F(o.isovalue)),
        ("precision", V::s(if o.precision == EvalPrecision::F64 { "f64" } else { "f32_stable" })),
        ("area_weights", V::B(o.area_weights.is_some())),
        ("device", dev_v(o.device)),
        ("threads", V::opt_u(o.threads)),
    ])
}

/// `BuildReport` as a dict.
pub fn build_report_v(r: &BuildReport) -> V {
    V::D(vec![
        ("n", V::u(r.n)),
        ("tokens", V::u(r.tokens)),
        ("knn_seconds", V::F(r.knn_seconds)),
        ("knn_device", dev_v(r.knn_device)),
        ("fit_device", dev_v(r.fit_device)),
        ("predict_device", dev_v(r.predict.device)),
        ("predict_chunks", V::u(r.predict.chunks)),
        ("features_seconds", V::F(r.predict.features_seconds)),
        ("forward_seconds", V::F(r.predict.forward_seconds)),
        ("predict_seconds", V::F(r.predict.seconds)),
        ("clamped_disc", V::u(r.fit.clamped_disc)),
        ("clamped_eigvec", V::u(r.fit.clamped_eigvec)),
        ("clamped_denominator", V::u(r.fit.clamped_denominator)),
        ("clamped_curvature", V::u(r.fit.clamped_curvature)),
        ("umbilic", V::u(r.fit.umbilic)),
        ("saddle", V::u(r.fit.saddle)),
        ("zero_curvature", V::u(r.fit.zero_curvature)),
        ("eigen_switched", V::u(r.fit.eigen_switched)),
        ("tokens_lt_k", V::u(r.fit.tokens_lt_k)),
        ("fit_seconds", V::F(r.fit.seconds)),
        ("seconds", V::F(r.seconds)),
    ])
}

/// `EvalReport` as a dict.
pub fn eval_report_v(r: &EvalReport) -> V {
    V::D(vec![
        ("queries", V::u(r.queries)),
        ("singular_axis", V::I(r.singular_axis as i64)),
        ("singular_core", V::I(r.singular_core as i64)),
        ("at_data_point", V::I(r.at_data_point as i64)),
        ("radius_fallbacks", V::I(r.radius_fallbacks as i64)),
        ("terms", V::I(r.terms as i64)),
        ("device", dev_v(r.device)),
        ("seconds", V::F(r.seconds)),
    ])
}

/// `CloudReport` as a dict.
pub fn cloud_report_v(r: &pat::CloudReport) -> V {
    V::D(vec![
        ("n", V::u(r.n)),
        ("normalized", V::u(r.normalized)),
        ("exact_duplicates", V::u(r.exact_duplicates)),
        ("bbox_min", V::arr3(r.bbox[0])),
        ("bbox_max", V::arr3(r.bbox[1])),
    ])
}

fn mc_report_v(r: &pat::mc::McReport) -> V {
    use pat::mc::MeshCheck;
    let mesh = match &r.mesh {
        MeshCheck::Skipped => V::D(vec![("status", V::s("skipped"))]),
        MeshCheck::Manifold { boundary_edges, euler_characteristic, components } => V::D(vec![
            ("status", V::s("manifold")),
            ("boundary_edges", V::u(*boundary_edges)),
            ("euler_characteristic", V::I(*euler_characteristic)),
            ("components", V::u(*components)),
        ]),
        MeshCheck::Invalid(m) => V::D(vec![("status", V::s("invalid")), ("message", V::S(m.clone()))]),
    };
    V::D(vec![
        ("cells", V::L(r.cells.iter().map(|&c| V::u(c)).collect())),
        ("spacing", V::arr3(r.spacing)),
        ("origin", V::arr3(r.origin)),
        ("blocks", V::u(r.blocks)),
        ("blocks_grown", V::u(r.blocks_grown)),
        ("band_rounds", V::u(r.band_rounds)),
        ("blocks_active", V::u(r.blocks_active)),
        ("cells_visited", V::u(r.cells_visited)),
        ("cells_crossing", V::u(r.cells_crossing)),
        ("ambiguous_cell_faces", V::u(r.ambiguous_cell_faces)),
        ("grid_values", V::u(r.grid_values)),
        ("grid_points", V::u(r.grid_points)),
        ("extra_evaluations", V::u(r.extra_evaluations)),
        ("vertices", V::u(r.vertices)),
        ("triangles", V::u(r.triangles)),
        ("lipschitz_bound", V::opt_f(r.lipschitz_bound)),
        ("max_edge_slope", V::F(r.max_edge_slope)),
        ("bound_violated", V::B(r.bound_violated)),
        ("zero_normals", V::u(r.zero_normals)),
        ("mesh", mesh),
        ("seconds", V::F(r.seconds[4])),
        ("eval", r.eval.as_ref().map_or(V::None, eval_report_v)),
    ])
}

// -------------------------------------------------------------------------------------------
// Fields as arguments: a Pat or a CSG tree
// -------------------------------------------------------------------------------------------

/// A shared implicit field (a torus field or a CSG tree).
pub type Shared = Arc<dyn ImplicitField + Send + Sync>;

/// Booleans.
pub fn boolean(a: Shared, b: Shared, op: &str, smooth: Option<f64>) -> R<Csg> {
    let (x, y) = (Csg::Leaf(a), Csg::Leaf(b));
    let c = match (op, smooth) {
        ("union", None) => x.union(y),
        ("intersect" | "intersection", None) => x.intersect(y),
        ("subtract" | "difference", None) => x.subtract(y),
        ("union", Some(k)) => Csg::SmoothUnion(Box::new(x), Box::new(y), k),
        ("intersect" | "intersection", Some(k)) => Csg::SmoothIntersect(Box::new(x), Box::new(y), k),
        ("subtract" | "difference", Some(k)) => Csg::SmoothSubtract(Box::new(x), Box::new(y), k),
        _ => return bad(format!("op: {op:?} is not one of \"union\", \"intersect\", \"subtract\"")),
    };
    c.validate()?;
    Ok(c)
}

/// An offset (`φ − c`).
pub fn offset(a: Shared, c: f64) -> R<Csg> {
    let t = Csg::Leaf(a).offset(c);
    t.validate()?;
    Ok(t)
}

/// Batch evaluation of any shared field (CSG trees).
pub fn eval_shared(f: &dyn ImplicitField, queries: &[f64], grad: bool) -> R<(Vec<f64>, Vec<f64>)> {
    let xs = to3("queries", queries)?;
    if let Some(i) = xs.iter().position(|x| x.iter().any(|c| !c.is_finite())) {
        return Err(PatError::NonFinite { index: i, what: "query" }.into());
    }
    let mut d = vec![0.0; xs.len()];
    let pool = pat::pool::make_pool(None)?;
    if grad {
        let mut g = vec![[0.0; 3]; xs.len()];
        pool.install(|| f.distance_gradients(&xs, &mut d, &mut g))?;
        Ok((d, flat3(&g)))
    } else {
        pool.install(|| f.distances(&xs, &mut d))?;
        Ok((d, vec![]))
    }
}

/// The box of a field.
pub fn bbox(f: &dyn ImplicitField) -> [f64; 6] {
    let b = f.bbox();
    [b[0][0], b[0][1], b[0][2], b[1][0], b[1][1], b[1][2]]
}

fn parse_mc(kw: &mut Kw) -> R<McOptions> {
    let mut m = McOptions::default();
    if let Some(b) = kw.opt_array("bbox")? {
        if b.len() != 6 {
            return bad(format!("bbox: expected 6 numbers (min x y z, max x y z), got {}", b.len()));
        }
        m.bbox = Some([[b[0], b[1], b[2]], [b[3], b[4], b[5]]]);
    }
    if let Some(r) = kw.opt_array("resolution")? {
        let ok = r.iter().all(|&x| x >= 1.0 && x.fract() == 0.0 && x < 4096.0);
        m.resolution = match r.len() {
            1 if ok => Resolution::Longest(r[0] as usize),
            3 if ok => Resolution::Cells([r[0] as usize, r[1] as usize, r[2] as usize]),
            _ => return bad(format!("resolution: expected a positive integer or 3 of them (below 4096), got {r:?}")),
        };
    }
    if let Some(l) = kw.opt_f64("level")? {
        m.isovalue = l;
    }
    let bound = kw.opt_f64("lipschitz")?;
    let block = kw.opt_usize("block")?;
    if kw.choice("band", &["lipschitz", "dense"])?.as_deref() == Some("dense") {
        m.band = Band::Dense;
    } else if bound.is_some() || block.is_some() {
        m.band = Band::Lipschitz { bound: bound.unwrap_or(2.0), block: block.unwrap_or(8) };
    }
    if let Some(r) = kw.opt_usize("refine_steps")? {
        m.refine_steps = r;
    }
    if let Some(n) = kw.opt_bool("normals")? {
        m.normals = n;
    }
    if let Some(c) = kw.opt_bool("check_mesh")? {
        m.check_mesh = c;
    }
    m.threads = parse_threads(kw)?;
    Ok(m)
}

/// A reconstructed mesh.
pub struct MeshOut {
    /// `3V`.
    pub positions: Vec<f64>,
    /// `3F`.
    pub faces: Vec<u32>,
    /// `3V` unit normals when requested.
    pub normals: Option<Vec<f64>>,
    /// The report.
    pub report: V,
    /// What the fitted fields evaluated during the run (`None` for analytic fields).
    pub eval: Option<EvalReport>,
}

/// Marching cubes on `field`.
pub fn reconstruct(field: &dyn ImplicitField, kw: &mut Kw) -> R<MeshOut> {
    let m = parse_mc(kw)?;
    kw.finish("reconstruct")?;
    let (mesh, rep) = marching_cubes(&field, &m)?;
    Ok(MeshOut {
        positions: flat3(&mesh.positions),
        faces: mesh.faces.iter().flatten().copied().collect(),
        normals: mesh.normals.as_ref().map(|n| flat3(n)),
        report: mc_report_v(&rep),
        eval: rep.eval.clone(),
    })
}

/// A sphere-traced image (`pat::render`): RGBA bytes, row 0 at the top, and the report.
pub struct ImageOut {
    /// `4 · width · height` bytes (alpha 0 on the plain background).
    pub rgba: Vec<u8>,
    /// Report.
    pub report: V,
}

/// Camera and image size of [`render`].
pub struct CameraSpec {
    /// Eye.
    pub eye: [f64; 3],
    /// Target.
    pub target: [f64; 3],
    /// Up.
    pub up: [f64; 3],
    /// Vertical field of view in degrees.
    pub fov: f64,
    /// Width.
    pub width: usize,
    /// Height.
    pub height: usize,
}

fn render_report_v(r: &pat::render::RenderReport) -> V {
    V::D(vec![
        ("rays", V::I(r.rays as i64)),
        ("clipped", V::I(r.clipped as i64)),
        ("hits", V::I(r.hits as i64)),
        ("misses", V::I(r.misses as i64)),
        ("step_limit", V::I(r.step_limit as i64)),
        ("non_finite", V::I(r.non_finite as i64)),
        ("steps_total", V::I(r.steps_total as i64)),
        ("steps_mean", V::F(r.steps_mean)),
        ("steps_max", V::I(r.steps_max as i64)),
        ("shadow_rays", V::I(r.shadow_rays as i64)),
        ("shadow_hits", V::I(r.shadow_hits as i64)),
        ("shadow_step_limit", V::I(r.shadow_step_limit as i64)),
        ("cut_plane_pixels", V::I(r.cut_plane_pixels as i64)),
        ("seconds", V::F(r.seconds)),
        ("threads", V::u(r.threads)),
        ("device", dev_v(r.device)),
    ])
}

/// Sphere-trace `field` (`pat::render::Tracer`, CPU). Options: `level`, `safety`, `g_min`,
/// `max_steps`, `eps_hit` (relative to the bounding-box diagonal), `sky` (`"forest"`, `"white"`),
/// `ground` (bool), `ground_up` (3 numbers), `demo` (the demo's constants), `base_color` (3
/// numbers in [0, 1]), `cut_axis` (0, 1, 2) with `cut_offset`, `cut_rotation_x`, `cut_rotation_y`
/// (radians), `threads`.
pub fn render(field: &dyn ImplicitField, cam: &CameraSpec, kw: &mut Kw) -> R<ImageOut> {
    use pat::render::{Camera, CutPlane, Environment, Tracer, TraceOptions};
    let mut o = if kw.opt_bool("demo")?.unwrap_or(false) { TraceOptions::demo() } else { TraceOptions::default() };
    if let Some(v) = kw.opt_f64("level")? {
        o.level = v;
    }
    if let Some(v) = kw.opt_f64("safety")? {
        o.safety = v;
    }
    if let Some(v) = kw.opt_f64("g_min")? {
        o.g_min = v;
    }
    if let Some(v) = kw.opt_usize("max_steps")? {
        o.max_steps = u32::try_from(v).map_err(|_| PErr(format!("max_steps = {v} is too large")))?;
    }
    if let Some(v) = kw.opt_f64("eps_hit")? {
        o.eps_hit = pat::render::Epsilon::Relative(v);
    }
    if let Some(s) = kw.choice("sky", &["forest", "white"])? {
        o.environment = if s == "white" { Environment::Constant([1.0, 1.0, 1.0]) } else { Environment::Forest };
    }
    if let Some(false) = kw.opt_bool("ground")? {
        o.ground = None;
    }
    if let Some(u) = kw.opt_array("ground_up")? {
        let u: [f64; 3] = u.try_into().map_err(|_| PErr("ground_up: expected 3 numbers".into()))?;
        if let Some(g) = o.ground.as_mut() {
            g.up = u;
        }
    }
    if let Some(c) = kw.opt_array("base_color")? {
        o.base_color = c.try_into().map_err(|_| PErr("base_color: expected 3 numbers".into()))?;
    }
    if let Some(axis) = kw.opt_usize("cut_axis")? {
        o.cut_plane = Some(CutPlane {
            axis,
            offset: kw.opt_f64("cut_offset")?.unwrap_or(0.0),
            rotation_x: kw.opt_f64("cut_rotation_x")?.unwrap_or(0.0),
            rotation_y: kw.opt_f64("cut_rotation_y")?.unwrap_or(0.0),
            ..Default::default()
        });
    }
    o.threads = parse_threads(kw)?;
    kw.finish("render")?;
    let camera = Camera::look_at(cam.eye, cam.target, cam.up, cam.fov, cam.width, cam.height)?;
    let (img, rep) = Tracer::new(o)?.render(field, &camera)?;
    let mut rgba = Vec::with_capacity(4 * img.alpha.len());
    for (c, a) in img.rgb.chunks_exact(3).zip(&img.alpha) {
        rgba.extend_from_slice(&[c[0], c[1], c[2], *a]);
    }
    Ok(ImageOut { rgba, report: render_report_v(&rep) })
}

fn parse_sample(kw: &mut Kw) -> R<SampleOptions> {
    let mut s = SampleOptions { count: kw.opt_usize("count")?, ..Default::default() };
    s.mc = parse_mc(kw)?;
    if let Some(v) = kw.opt_usize("newton_iters")? {
        s.newton_iters = v;
    }
    if let Some(v) = kw.opt_f64("residual_tol_rel")? {
        s.residual_tol_rel = v;
    }
    if let Some(v) = kw.opt_f64("max_move_cells")? {
        s.max_move_cells = v;
    }
    if let Some(v) = kw.opt_u64("seed")? {
        s.seed = v;
    }
    if let Some(v) = kw.opt_usize("rounds")? {
        s.rounds = v;
    }
    s.threads = s.mc.threads;
    Ok(s)
}

fn sample_report_v(r: &pat::apps::SampleReport) -> V {
    V::D(vec![
        ("requested", V::u(r.requested)),
        ("produced", V::u(r.produced)),
        ("rejected", V::u(r.rejected)),
        ("residual_max", V::F(r.residual_max)),
        ("residual_mean", V::F(r.residual_mean)),
        ("residual_before_max", V::F(r.residual_before_max)),
        ("rounds", V::u(r.rounds)),
        ("mc", mc_report_v(&r.mc)),
        ("seconds", V::F(r.seconds)),
    ])
}

/// Points and unit normals on `{φ = level}`.
pub fn level_set(field: &dyn ImplicitField, level: f64, kw: &mut Kw) -> R<(Vec<f64>, Vec<f64>, V)> {
    let s = parse_sample(kw)?;
    kw.finish("sample_level_set")?;
    let (p, n, rep) = sample_level_set(field, level, &s)?;
    Ok((flat3(&p), flat3(&n), sample_report_v(&rep)))
}

/// Erode, dilate, open or close; the result is a new field with its cloud.
pub fn morphology(field: &dyn ImplicitField, op: &str, radius: f64, kw: &mut Kw) -> R<(Pat, V)> {
    let morph = match op {
        "erode" => Morph::Erode(radius),
        "dilate" => Morph::Dilate(radius),
        "open" => Morph::Open(radius),
        "close" => Morph::Close(radius),
        _ => return bad(format!("op: {op:?} is not one of \"erode\", \"dilate\", \"open\", \"close\"")),
    };
    if !(radius.is_finite() && radius >= 0.0) {
        return bad(format!("radius: {radius} (finite and >= 0)"));
    }
    let levels = morph.levels()?;
    let mut sample = parse_sample(kw)?;
    let mut spec = parse_build(kw, true)?;
    kw.finish("morphology")?;
    spec.policy = NormalPolicy::Normalize;
    let th = spec.opts.eval.threads.or(sample.threads);
    spec.opts = spec.opts.clone().with_threads(th);
    sample.threads = sample.threads.or(spec.opts.eval.threads);
    if sample.count.is_none() {
        sample.count = field.point_count_hint();
    }
    let m = model(spec.weights.as_deref())?;
    let t0 = std::time::Instant::now();
    let mut steps = Vec::new();
    let mut current: Option<Pat> = None;
    for level in levels {
        let src: &dyn ImplicitField = match &current {
            Some(p) => &*p.field,
            None => field,
        };
        let (pts, nrm, srep) = sample_level_set(src, level, &sample)?;
        let (cloud, creport) = PointCloud::new(pts, nrm, NormalPolicy::Normalize)?;
        let mut opts = spec.opts.clone();
        apply_eval(&spec, &cloud, &mut opts.eval)?;
        let (f, brep) = PointsAsTori::build(&cloud, &m, &opts)?;
        steps.push(V::D(vec![("level", V::F(level)), ("sample", sample_report_v(&srep)), ("build", build_report_v(&brep))]));
        current = Some(Pat { cloud, field: Arc::new(f), build_report: Some(brep), cloud_report: creport, last_eval: Mutex::new(None) });
    }
    match current {
        Some(p) => Ok((p, V::D(vec![("steps", V::L(steps)), ("seconds", V::F(t0.elapsed().as_secs_f64()))]))),
        None => bad("morphology without steps"),
    }
}

/// Fig. 27: flip normals until the field agrees with them.
pub fn orient(points: &[f64], normals: &[f64], kw: &mut Kw) -> R<(Vec<f64>, V)> {
    let mut o = OrientOptions::default();
    if let Some(k) = kw.opt_usize("k")? {
        o.k = k;
    }
    // `k` of the build group is the network's; the orientation vote has its own name.
    if let Some(k) = kw.opt_usize("vote_k")? {
        o.k = k;
    }
    let scale = kw.opt_f64("epsilon_scale")?;
    let constant = kw.opt_f64("epsilon")?;
    o.epsilon = match (scale, constant) {
        (Some(_), Some(_)) => return bad("epsilon: give either epsilon_scale or epsilon, not both"),
        (_, Some(c)) => Epsilon::Constant(c),
        (s, None) => Epsilon::InverseSigma { scale: s.unwrap_or(1.0) },
    };
    if let Some(i) = kw.opt_usize("max_iterations")? {
        o.max_iterations = i;
    }
    let spec = parse_build(kw, false)?;
    kw.finish("orient_normals")?;
    o.normals = spec.policy;
    o.build = spec.opts.clone();
    o.build.k = None;
    o.threads = spec.opts.eval.threads;
    let m = model(spec.weights.as_deref())?;
    let (n, rep) = orient_normals(&to3("points", points)?, &to3("normals", normals)?, &m, &o)?;
    let report = V::D(vec![
        ("n", V::u(rep.n)),
        ("flips_per_iteration", V::L(rep.flips_per_iteration.iter().map(|&x| V::u(x)).collect())),
        ("converged", V::B(rep.converged)),
        ("flipped_total", V::u(rep.flipped_total)),
        ("negative_final", V::u(rep.negative_final)),
        ("negative_per_iteration", V::L(rep.negative_per_iteration.iter().map(|&x| V::u(x)).collect())),
        ("undecided", V::u(rep.undecided)),
        ("seconds", V::F(rep.seconds)),
    ]);
    Ok((flat3(&n), report))
}

// -------------------------------------------------------------------------------------------
// Baselines
// -------------------------------------------------------------------------------------------

fn cloud_of(points: &[f64], normals: &[f64], kw: &mut Kw) -> R<PointCloud> {
    let policy = match kw.choice("normal_policy", &["require_unit", "normalize"])?.as_deref() {
        Some("normalize") => NormalPolicy::Normalize,
        _ => NormalPolicy::RequireUnit { tol: kw.opt_f64("normal_tol")?.unwrap_or(1e-6) },
    };
    let _ = kw.opt_f64("normal_tol")?;
    Ok(PointCloud::new(to3("points", points)?, to3("normals", normals)?, policy)?.0)
}

/// Barill's point areas.
pub fn areas(points: &[f64], normals: &[f64], kw: &mut Kw) -> R<(Vec<f64>, V)> {
    let cloud = cloud_of(points, normals, kw)?;
    let mut o = PointAreaOptions::default();
    if let Some(k) = kw.opt_usize("k")? {
        o.k = k;
    }
    if let Some(f) = kw.opt_bool("filter_opposite")? {
        o.filter_opposite = f;
    }
    o.threads = parse_threads(kw)?;
    kw.finish("point_areas")?;
    let (a, r) = point_areas(&cloud, &o)?;
    Ok((a, V::D(vec![
        ("n", V::u(r.n)),
        ("clipped", V::u(r.clipped)),
        ("filtered", V::I(r.filtered as i64)),
        ("coincident", V::u(r.coincident)),
        ("isolated", V::u(r.isolated)),
        ("total_area", V::F(r.total_area)),
        ("seconds", V::F(r.seconds)),
    ])))
}

fn area_arg(cloud: &PointCloud, kw: &mut Kw, threads: Option<usize>) -> R<Vec<f64>> {
    match kw.take("areas") {
        None | Some(O::None) => Ok(point_areas(cloud, &PointAreaOptions { threads, ..Default::default() })?.0),
        Some(O::Array(a)) => Ok(a.clone()),
        Some(o) => bad(format!("areas: expected None or one number per point, got {o:?}")),
    }
}

/// Generalized winding number (direct sum) or its regularization (`regularization = eps`).
pub fn gwn(points: &[f64], normals: &[f64], queries: &[f64], kw: &mut Kw) -> R<(Vec<f64>, V)> {
    let cloud = cloud_of(points, normals, kw)?;
    let threads = parse_threads(kw)?;
    let eps = kw.opt_f64("regularization")?;
    let areas = area_arg(&cloud, kw, threads)?;
    kw.finish("gwn")?;
    let xs = to3("queries", queries)?;
    let mut out = vec![0.0; xs.len()];
    let rep = match eps {
        None => pat::baselines::winding_numbers(&cloud, &areas, &xs, &mut out, threads)?,
        Some(e) => pat::baselines::regularized_winding_numbers(&cloud, &areas, e, &xs, &mut out, threads)?,
    };
    Ok((out, V::D(vec![("queries", V::u(rep.queries)), ("skipped", V::I(rep.skipped as i64)), ("seconds", V::F(rep.seconds))])))
}

fn parse_fwn(kw: &mut Kw) -> R<FwnOptions> {
    let mut o = FwnOptions::default();
    if let Some(b) = kw.opt_f64("beta")? {
        o.beta = b;
    }
    if let Some(v) = kw.opt_usize("order")? {
        o.order = u8::try_from(v).map_err(|_| PErr(format!("order: {v} (0, 1 or 2)")))?;
    }
    o.leaf_size = kw.opt_usize("leaf_size")?;
    o.threads = parse_threads(kw)?;
    Ok(o)
}

fn fwn_v(r: &pat::baselines::FwnReport) -> V {
    V::D(vec![
        ("queries", V::u(r.queries)),
        ("far_nodes", V::I(r.far_nodes as i64)),
        ("near_elements", V::I(r.near_elements as i64)),
        ("skipped", V::I(r.skipped as i64)),
        ("seconds", V::F(r.seconds)),
    ])
}

/// Barill's fast winding number of a point cloud.
pub fn fwn(points: &[f64], normals: &[f64], queries: &[f64], kw: &mut Kw) -> R<(Vec<f64>, V)> {
    let cloud = cloud_of(points, normals, kw)?;
    let o = parse_fwn(kw)?;
    let areas = area_arg(&cloud, kw, o.threads)?;
    kw.finish("fwn")?;
    let tree = FastWinding::from_cloud(&cloud, &areas, &o)?;
    let xs = to3("queries", queries)?;
    let mut out = vec![0.0; xs.len()];
    let rep = tree.evaluate(&xs, &mut out)?;
    Ok((out, fwn_v(&rep)))
}

/// Fast winding number of a triangle soup.
pub fn fwn_mesh(vertices: &[f64], faces: &[u32], queries: &[f64], kw: &mut Kw) -> R<(Vec<f64>, V)> {
    let o = parse_fwn(kw)?;
    kw.finish("fwn_mesh")?;
    if !faces.len().is_multiple_of(3) {
        return bad(format!("faces: {} indices are not whole triangles", faces.len()));
    }
    let f: Vec<[u32; 3]> = faces.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let tree = FastWinding::from_triangles(&to3("vertices", vertices)?, &f, &o)?;
    let xs = to3("queries", queries)?;
    let mut out = vec![0.0; xs.len()];
    let rep = tree.evaluate(&xs, &mut out)?;
    Ok((out, fwn_v(&rep)))
}

/// SHC, SSPD or LSE at the queries. Returns `(values, gradient)`.
pub fn kernel(method: &str, points: &[f64], normals: &[f64], queries: &[f64], gradient: bool, kw: &mut Kw) -> R<(Vec<f64>, Vec<f64>, V)> {
    let (m, name) = match method {
        "shc" => (KernelMethod::Shc, "shc"),
        "sspd" => (KernelMethod::Sspd, "sspd"),
        "lse" => (KernelMethod::Lse, "lse"),
        _ => return bad(format!("method: {method:?} is not one of \"shc\", \"sspd\", \"lse\"")),
    };
    let cloud = cloud_of(points, normals, kw)?;
    let threads = parse_threads(kw)?;
    let ev = parse_eval_overrides(kw)?;
    let zero = kw.opt_f64("shc_zero_value")?;
    let mut o = KernelOptions { threads, shc_zero_value: zero, ..Default::default() };
    if let Some(v) = ev.neighborhood {
        o.neighborhood = v;
    }
    if let Some(v) = ev.lambda {
        o.lambda = v;
    }
    if kw.opt_bool("paper")?.unwrap_or(false) {
        let (p, _) = TorusField::paper_options(&cloud, threads)?;
        o.neighborhood = ev.neighborhood.unwrap_or(p.neighborhood);
        o.lambda = ev.lambda.unwrap_or(p.lambda);
    }
    let areas = if m == KernelMethod::Lse && matches!(kw.take("areas"), None | Some(O::None)) {
        vec![1.0; cloud.len()]
    } else {
        area_arg(&cloud, kw, threads)?
    };
    kw.finish(name)?;
    let field = KernelField::new(&cloud, areas, o)?;
    let xs = to3("queries", queries)?;
    let mut vals = vec![0.0; xs.len()];
    let mut g = vec![[0.0; 3]; if gradient { xs.len() } else { 0 }];
    let rep = field.evaluate(m, &xs, &mut vals, gradient.then_some(&mut g[..]))?;
    Ok((vals, flat3(&g), V::D(vec![
        ("queries", V::u(rep.queries)),
        ("terms", V::I(rep.terms as i64)),
        ("skipped", V::I(rep.skipped as i64)),
        ("radius_fallbacks", V::I(rep.radius_fallbacks as i64)),
        ("degenerate", V::I(rep.degenerate as i64)),
        ("seconds", V::F(rep.seconds)),
    ])))
}

// -------------------------------------------------------------------------------------------
// Mesh truth
// -------------------------------------------------------------------------------------------

/// The exact signed distance to a closed triangle mesh.
pub struct Truth {
    /// The BVH.
    pub truth: MeshTruth,
    /// Topology of the mesh.
    pub report: V,
}

/// Build the mesh truth.
pub fn truth(vertices: &[f64], faces: &[u32], kw: &mut Kw) -> R<Truth> {
    let beta = kw.opt_f64("beta")?;
    let sign = match kw.choice("sign", &["auto", "exact", "fast"])?.as_deref() {
        Some("exact") => SignMode::Exact,
        Some("fast") => SignMode::Fast { beta: beta.unwrap_or(2.0) },
        _ => SignMode::Auto,
    };
    let mut o = TruthOptions { sign, ..Default::default() };
    if let Some(t) = kw.opt_f64("threshold")? {
        o.threshold = t;
    }
    if let Some(c) = kw.opt_bool("require_closed")? {
        o.require_closed = c;
    }
    o.threads = parse_threads(kw)?;
    kw.finish("MeshTruth")?;
    if !faces.len().is_multiple_of(3) {
        return bad(format!("faces: {} indices are not whole triangles", faces.len()));
    }
    let f: Vec<[u32; 3]> = faces.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    let (t, r) = MeshTruth::new(&to3("vertices", vertices)?, &f, &o)?;
    let report = V::D(vec![
        ("vertices", V::u(r.vertices)),
        ("faces", V::u(r.faces)),
        ("degenerate_faces", V::u(r.degenerate_faces)),
        ("boundary_edges", V::u(r.boundary_edges)),
        ("nonmanifold_edges", V::u(r.nonmanifold_edges)),
        ("inconsistent_edges", V::u(r.inconsistent_edges)),
        ("signed_volume", V::F(r.signed_volume)),
        ("closed", V::B(r.is_closed())),
    ]);
    Ok(Truth { truth: t, report })
}

/// Signed distance (and winding number) of the mesh at the queries.
pub fn truth_eval(t: &MeshTruth, queries: &[f64], winding: bool) -> R<(Vec<f64>, Vec<f64>, V)> {
    let xs = to3("queries", queries)?;
    let mut d = vec![0.0; xs.len()];
    let mut w = vec![0.0; if winding { xs.len() } else { 0 }];
    let r = t.evaluate(&xs, &mut d, winding.then_some(&mut w[..]))?;
    Ok((d, w, V::D(vec![
        ("queries", V::u(r.queries)),
        ("inside", V::u(r.inside)),
        ("on_surface", V::u(r.on_surface)),
        ("ambiguous", V::u(r.ambiguous)),
        ("seconds", V::F(r.seconds)),
    ])))
}

// -------------------------------------------------------------------------------------------
// The classical 2D variant
// -------------------------------------------------------------------------------------------

/// A 2D field with its report.
pub struct Pat2 {
    /// The field.
    pub field: pat2d::Pat2dField,
    /// Build report.
    pub report: V,
}

/// Build the classical 2D variant (a circle per point from a weighted quadratic fit).
pub fn build_2d(points: &[f64], normals: &[f64], kw: &mut Kw) -> R<Pat2> {
    let policy = match kw.choice("normal_policy", &["require_unit", "normalize"])?.as_deref() {
        Some("normalize") => NormalPolicy::Normalize,
        _ => NormalPolicy::RequireUnit { tol: kw.opt_f64("normal_tol")?.unwrap_or(1e-6) },
    };
    let _ = kw.opt_f64("normal_tol")?;
    let threads = parse_threads(kw)?;
    let mut fit = pat2d::ClassicalFit::default();
    if let Some(k) = kw.opt_usize("fit_k")? {
        fit.k = k;
    }
    if let Some(b) = kw.opt_f64("bandwidth")? {
        fit.bandwidth = b;
    }
    let ev = parse_eval_overrides(kw)?;
    let paper = kw.opt_bool("paper")?.unwrap_or(false);
    kw.finish("Pat2d")?;
    let (cloud, _) = pat2d::Cloud2::new(to2("points", points)?, to2("normals", normals)?, policy)?;
    let mut eval = if paper { pat2d::Pat2dField::paper_options(&cloud, threads)?.0 } else { pat2d::Eval2dOptions::default() };
    eval.threads = threads;
    if let Some(v) = ev.neighborhood {
        eval.neighborhood = v;
    }
    if let Some(v) = ev.lambda {
        eval.lambda = v;
    }
    if let Some(v) = ev.gradient {
        eval.gradient = v;
    }
    if let Some(v) = ev.isovalue {
        eval.isovalue = v;
    }
    let mut opts = pat2d::BuildOptions2 { eval, ..Default::default() };
    opts.fit.threads = threads;
    let (field, r) = pat2d::Pat2d::build(&cloud, &fit, &opts)?;
    let report = V::D(vec![
        ("n", V::u(r.n)),
        ("knn_seconds", V::F(r.knn_seconds)),
        ("predict_seconds", V::F(r.predict_seconds)),
        ("lines", V::u(r.fit.lines)),
        ("circles", V::u(r.fit.circles)),
        ("degenerate_fits", V::u(r.fit.degenerate_fits)),
        ("tokens", V::u(r.fit.tokens)),
        ("seconds", V::F(r.seconds)),
    ]);
    Ok(Pat2 { field, report })
}

/// 2D evaluation: `(distance, gradient)`.
pub fn eval_2d(f: &pat2d::Pat2dField, queries: &[f64], grad: bool) -> R<(Vec<f64>, Vec<f64>)> {
    let xs = to2("queries", queries)?;
    let mut d = vec![0.0; xs.len()];
    let mut g = vec![[0.0; 2]; if grad { xs.len() } else { 0 }];
    f.evaluate(&xs, Some(&mut d), grad.then_some(&mut g[..]))?;
    Ok((d, flat2(&g)))
}

/// Contour polylines `{φ = level}` of the 2D field on an `nx × ny` grid over `bbox`
/// `[xmin, ymin, xmax, ymax]`: `(flat points, closed)` per polyline.
pub fn contours_2d(f: &pat2d::Pat2dField, bbox: [f64; 4], nx: usize, ny: usize, level: f64) -> R<Vec<(Vec<f64>, bool)>> {
    let g = pat2d::evaluate_grid(f, [[bbox[0], bbox[1]], [bbox[2], bbox[3]]], nx, ny)?;
    let segs = pat2d::marching_squares(&g, level);
    Ok(pat2d::chain_segments(&segs).into_iter().map(|p| (flat2(&p.points), p.closed)).collect())
}

/// Whether a file exists (tests).
pub fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

/// The defaults of the constructor options, read from the library.
pub fn defaults() -> V {
    let b = BuildOptions::default();
    let mut items = match eval_options_v(&b.eval) {
        V::D(d) => d,
        _ => vec![],
    };
    items.retain(|(k, _)| !matches!(*k, "area_weights" | "device" | "threads"));
    let g = PredictOptionsDefaults::get();
    items.extend([
        ("normal_policy", V::s("require_unit")),
        ("normal_tol", V::F(1e-6)),
        ("weights", V::None),
        ("device", V::s("cpu")),
        ("cuda_ordinal", V::I(0)),
        ("threads", V::None),
        ("gpu_matmul", V::s(g.0)),
        ("feature_order", V::s(g.1)),
        ("sign_rule", V::s(if b.fit.sign_rule == SignRule::Reference { "reference" } else { "geometric" })),
        ("eigen_rule", V::s(if b.fit.eigen_rule == EigenRule::Reference { "reference" } else { "robust" })),
        ("k", V::None),
        ("paper", V::B(false)),
        ("area_weights", V::None),
        ("progress", if b.progress.is_none() { V::None } else { V::s("set") }),
    ]);
    V::D(items)
}

struct PredictOptionsDefaults;

impl PredictOptionsDefaults {
    /// The binding's defaults for the GPU matmul (SPEC.md section 18: TF32) and the feature order.
    fn get() -> (&'static str, &'static str) {
        ("tf32", if FeatureOrder::default() == FeatureOrder::DifferenceFirst { "difference_first" } else { "reference" })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unknown_and_mistyped_options_are_rejected() {
        let mut kw = Kw::new(vec![("device".into(), O::Str("cpu".into())), ("bogus".into(), O::Int(1))]);
        parse_build(&mut kw, false).unwrap();
        assert!(kw.finish("x").unwrap_err().0.contains("bogus"));
        let mut kw = Kw::new(vec![("sign_rule".into(), O::Str("nope".into()))]);
        assert!(parse_build(&mut kw, false).unwrap_err().0.contains("not one of"));
        let mut kw = Kw::new(vec![("threads".into(), O::Int(0))]);
        assert!(parse_build(&mut kw, false).is_err());
    }

    #[test]
    fn a_small_cloud_builds_and_round_trips_the_tori() {
        let n = 200;
        let (mut p, mut nr) = (vec![], vec![]);
        for i in 0..n {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let r = (1.0 - z * z).sqrt();
            let a = 2.399963 * i as f64;
            let d = [r * a.cos(), r * a.sin(), z];
            nr.extend(d);
            p.extend(d.iter().map(|c| 0.5 * c));
        }
        let spec = parse_build(&mut Kw::default(), false).unwrap();
        let f = Pat::build(&p, &nr, &spec).unwrap();
        let (d, _, _) = f.evaluate(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], true, false, false).unwrap();
        assert!(d[0] < -0.3 && (d[1] - 0.5).abs() < 0.05, "{d:?}");
        let t = export_tori(f.tori());
        let back = import_tori(&t.centers, &t.axes, &t.major, &t.minor).unwrap();
        let g = Pat::from_tori(&p, &nr, back, None, &spec).unwrap();
        let (e, _, _) = g.evaluate(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], true, false, false).unwrap();
        assert!((e[0] - d[0]).abs() < 1e-8);
        assert!(matches!(defaults(), V::D(_)));
    }
}
