//! The CUDA back end (feature `cuda`; SPEC.md §5.4, §7.5, §18.2 D4; `NOTES-GPU.md`).
//!
//! Two parts:
//! - **the network** through burn-cuda ([`crate::network::cuda`]): selected with
//!   [`crate::Device`] in [`crate::PredictOptions`];
//! - **our NVRTC kernels** (`kernels.cu`, compiled at run time by NVRTC through `cudarc`, nothing
//!   CUDA at build time): exact kNN on the flattened kd-tree, fitting neighborhoods, the torus
//!   fit, and the blend with its gradient and Laplacian in f64 (default) or in f32 with the stable
//!   form ([`crate::EvalPrecision::F32Stable`], SPEC.md §7.3).
//!
//! **Determinism** (the workspace GPU rule). Every kernel computes one output per thread in the
//! CPU's operation order, with gathers only (no atomics, no cross-thread reductions), compiled
//! with `--fmad=false`, IEEE division and square root. Hence:
//! - level 1, bitwise run to run, for everything (the network included, measured);
//! - level 2, bitwise equal to the CPU, for kNN, fitting neighborhoods, the fit and the per-torus
//!   SDF; the blend uses the device's `exp` instead of the host libm's, so it matches the CPU
//!   within a measured tolerance (NOTES-GPU.md), bitwise wherever the two `exp` agree.
//!
//! **Devices.** [`crate::Device::Cuda`] is an error when the device cannot be used (no driver,
//! no device of that ordinal, no NVRTC); [`crate::Device::Auto`] uses device 0 when usable and
//! the CPU otherwise. Every report says which device ran each stage. A [`GpuContext`] is opened
//! explicitly per field or build; the compiled PTX is memoized per compute capability (a pure
//! function of the source and the options, so it carries no state that could change a result).

use crate::field::{EvalOut, EvalReport, Neighborhood};
use crate::fit::{EigenRule, FitOptions, FitReport, SignRule, Tori, Torus};
use crate::spatial::{KdTree, Neighborhoods};
use crate::{Device, GradientMode, LambdaRule, PatError, PointCloud};
use cudarc::driver::{CudaContext, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchConfig, PushKernelArg, ValidAsZeroBits};
use cudarc::nvrtc::{CompileOptions, Ptx, compile_ptx_with_opts};
use std::sync::{Arc, Mutex};

const SOURCE: &str = include_str!("kernels.cu");
const BLOCK: u32 = 128;
/// Queries per launch group (bounds device memory; results do not depend on it).
const QUERY_GROUP: usize = 1 << 20;
/// The largest `k` of the kNN kernels (larger `k < n` stays on the CPU).
pub const KMAX: usize = 128;

fn gpu_err(what: &'static str) -> impl FnOnce(cudarc::driver::DriverError) -> PatError {
    move |e| PatError::Gpu { reason: format!("{what}: {e:?}") }
}

fn launch_cfg(n: u64) -> Result<LaunchConfig, PatError> {
    let blocks = n.div_ceil(BLOCK as u64);
    if blocks > i32::MAX as u64 {
        return Err(PatError::Gpu { reason: format!("{n} threads exceed the grid limit") });
    }
    Ok(LaunchConfig { grid_dim: (blocks.max(1) as u32, 1, 1), block_dim: (BLOCK, 1, 1), shared_mem_bytes: 0 })
}

/// One CUDA device as the driver reports it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GpuDeviceInfo {
    /// Ordinal.
    pub ordinal: usize,
    /// Name.
    pub name: String,
    /// Compute capability.
    pub compute_capability: (i32, i32),
    /// Total memory in bytes.
    pub total_memory: u64,
}

fn open(ordinal: usize) -> Result<Arc<CudaContext>, PatError> {
    let n = std::panic::catch_unwind(CudaContext::device_count)
        .map_err(|_| PatError::Gpu { reason: "the CUDA driver library (libcuda.so) could not be loaded".into() })?
        .map_err(gpu_err("device count"))?;
    if ordinal >= usize::try_from(n).unwrap_or(0) {
        return Err(PatError::Gpu { reason: format!("CUDA device {ordinal} does not exist ({n} devices)") });
    }
    CudaContext::new(ordinal).map_err(gpu_err("context"))
}

fn info_of(ctx: &CudaContext, ordinal: usize) -> Result<GpuDeviceInfo, PatError> {
    Ok(GpuDeviceInfo {
        ordinal,
        name: ctx.name().map_err(gpu_err("name"))?,
        compute_capability: ctx.compute_capability().map_err(gpu_err("compute capability"))?,
        total_memory: ctx.total_mem().map_err(gpu_err("memory"))? as u64,
    })
}

/// Open device `ordinal` with the driver (no kernels compiled). An error, never a panic, when the
/// driver library is missing or the ordinal does not exist.
pub fn probe(ordinal: usize) -> Result<GpuDeviceInfo, PatError> {
    let ctx = open(ordinal)?;
    info_of(&ctx, ordinal)
}

/// The CUDA devices the driver reports.
pub fn devices() -> Result<Vec<GpuDeviceInfo>, PatError> {
    let n = std::panic::catch_unwind(CudaContext::device_count)
        .map_err(|_| PatError::Gpu { reason: "the CUDA driver library (libcuda.so) could not be loaded".into() })?
        .map_err(gpu_err("device count"))?;
    (0..usize::try_from(n).unwrap_or(0)).map(probe).collect()
}

/// The ordinal `device` asks for: `None` for the CPU (or `Auto` without a usable GPU).
pub fn resolve_ordinal(device: Device) -> Result<Option<usize>, PatError> {
    match device {
        Device::Cpu => Ok(None),
        Device::Cuda { ordinal } => probe(ordinal).map(|_| Some(ordinal)),
        Device::Auto => Ok(probe(0).ok().map(|_| 0)),
    }
}

/// The GPU context `device` asks for: `None` for the CPU (or `Auto` without a usable GPU).
pub fn resolve(device: Device) -> Result<Option<Arc<GpuContext>>, PatError> {
    match device {
        Device::Cpu => Ok(None),
        Device::Cuda { ordinal } => GpuContext::new(ordinal).map(Some),
        Device::Auto => Ok(GpuContext::new(0).ok()),
    }
}

struct Kernels {
    knn_batch: [CudaFunction; 3],
    fitting: [CudaFunction; 3],
    eval_knn: [CudaFunction; 3],
    eval_knn_f32: [CudaFunction; 3],
    eval_all: CudaFunction,
    eval_all_f32: CudaFunction,
    radius_counts: CudaFunction,
    eval_radius: CudaFunction,
    eval_radius_f32: CudaFunction,
    fit_tori: CudaFunction,
    exp_values: CudaFunction,
}

/// The template instance (32, 64 or 128) for `k ≤ KMAX`.
fn kslot(k: usize) -> usize {
    if k <= 32 {
        0
    } else if k <= 64 {
        1
    } else {
        2
    }
}

static PTX_CACHE: Mutex<Vec<((i32, i32), Ptx)>> = Mutex::new(Vec::new());

fn compiled(cc: (i32, i32)) -> Result<Ptx, PatError> {
    if let Ok(cache) = PTX_CACHE.lock()
        && let Some((_, p)) = cache.iter().find(|(c, _)| *c == cc)
    {
        return Ok(p.clone());
    }
    let opts = CompileOptions {
        fmad: Some(false),
        prec_div: Some(true),
        prec_sqrt: Some(true),
        ftz: Some(false),
        options: vec![format!("--gpu-architecture=sm_{}{}", cc.0, cc.1)],
        name: Some("ddg_points_as_tori_kernels.cu".into()),
        ..Default::default()
    };
    let ptx = std::panic::catch_unwind(|| compile_ptx_with_opts(SOURCE, opts))
        .map_err(|_| PatError::Gpu { reason: "NVRTC (libnvrtc.so.12) could not be loaded".into() })?
        .map_err(|e| PatError::Gpu { reason: format!("NVRTC: {e:?}") })?;
    if let Ok(mut cache) = PTX_CACHE.lock() {
        cache.push((cc, ptx.clone()));
    }
    Ok(ptx)
}

/// A CUDA context, its stream and our compiled kernels.
pub struct GpuContext {
    _ctx: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    k: Kernels,
    info: GpuDeviceInfo,
}

impl std::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuContext").field("device", &self.info).finish()
    }
}

impl GpuContext {
    /// Open device `ordinal` and load the kernels (NVRTC compiles them on the first call per
    /// compute capability).
    pub fn new(ordinal: usize) -> Result<Arc<GpuContext>, PatError> {
        let ctx = open(ordinal)?;
        let info = info_of(&ctx, ordinal)?;
        // One stream, used in order.
        // SAFETY: every buffer of this context is used on `stream` only.
        unsafe { ctx.disable_event_tracking() };
        let module = ctx.load_module(compiled(info.compute_capability)?).map_err(gpu_err("module"))?;
        let f = |n: &str| module.load_function(n).map_err(gpu_err("module"));
        let three = |p: &str| -> Result<[CudaFunction; 3], PatError> { Ok([f(&format!("{p}_32"))?, f(&format!("{p}_64"))?, f(&format!("{p}_128"))?]) };
        let k = Kernels {
            knn_batch: three("knn_batch")?,
            fitting: three("fitting")?,
            eval_knn: three("eval_knn")?,
            eval_knn_f32: three("eval_knn_f32")?,
            eval_all: f("eval_all")?,
            eval_all_f32: f("eval_all_f32")?,
            radius_counts: f("radius_counts")?,
            eval_radius: f("eval_radius")?,
            eval_radius_f32: f("eval_radius_f32")?,
            fit_tori: f("fit_tori")?,
            exp_values: f("exp_values")?,
        };
        let stream = ctx.default_stream();
        Ok(Arc::new(GpuContext { _ctx: ctx, stream, k, info }))
    }

    /// The device.
    pub fn info(&self) -> &GpuDeviceInfo {
        &self.info
    }

    /// The [`Device`] value of this context.
    pub fn device(&self) -> Device {
        Device::Cuda { ordinal: self.info.ordinal }
    }

    fn up<T: DeviceRepr + Copy>(&self, v: &[T], dummy: T) -> Result<CudaSlice<T>, PatError> {
        if v.is_empty() {
            return self.stream.clone_htod(&[dummy]).map_err(gpu_err("upload"));
        }
        self.stream.clone_htod(v).map_err(gpu_err("upload"))
    }

    fn zeros<T: DeviceRepr + ValidAsZeroBits>(&self, n: usize) -> Result<CudaSlice<T>, PatError> {
        self.stream.alloc_zeros::<T>(n.max(1)).map_err(gpu_err("alloc"))
    }

    fn down<T: DeviceRepr + Default + Clone>(&self, d: &CudaSlice<T>, n: usize) -> Result<Vec<T>, PatError> {
        let mut out = vec![T::default(); n];
        if n > 0 {
            self.stream.memcpy_dtoh(&d.slice(0..n), &mut out).map_err(gpu_err("download"))?;
        }
        self.stream.synchronize().map_err(gpu_err("synchronize"))?;
        Ok(out)
    }

    /// CUDA's `exp` of `x` (for the libm comparison of NOTES-GPU.md).
    pub fn exp_values(&self, x: &[f64]) -> Result<Vec<f64>, PatError> {
        let n = x.len() as u64;
        let dx = self.up(x, 0.0)?;
        let mut dy = self.zeros::<f64>(x.len())?;
        unsafe { self.stream.launch_builder(&self.k.exp_values).arg(&n).arg(&dx).arg(&mut dy).launch(launch_cfg(n)?) }
            .map_err(gpu_err("exp"))?;
        self.down(&dy, x.len())
    }
}

/// A kd-tree on the device.
pub(crate) struct GpuTree {
    n: u32,
    lo: CudaSlice<f64>,
    hi: CudaSlice<f64>,
    start: CudaSlice<u32>,
    end: CudaSlice<u32>,
    right: CudaSlice<u32>,
    pts: CudaSlice<f64>,
    idx: CudaSlice<u32>,
}

impl GpuTree {
    pub(crate) fn upload(g: &GpuContext, tree: &KdTree) -> Result<GpuTree, PatError> {
        let (lo, hi, start, end, right, pts, idx) = tree.flat();
        Ok(GpuTree {
            n: tree.len() as u32,
            lo: g.up(&lo, 0.0)?,
            hi: g.up(&hi, 0.0)?,
            start: g.up(&start, 0)?,
            end: g.up(&end, 0)?,
            right: g.up(&right, 0)?,
            pts: g.up(&pts, 0.0)?,
            idx: g.up(idx, 0)?,
        })
    }
}

macro_rules! tree_args {
    ($b:expr, $t:expr) => {
        $b.arg(&$t.lo).arg(&$t.hi).arg(&$t.start).arg(&$t.end).arg(&$t.right).arg(&$t.pts).arg(&$t.idx)
    };
}

fn flat3(v: &[[f64; 3]]) -> Vec<f64> {
    v.iter().flatten().copied().collect()
}

/// The `k` nearest data points of each query by `(d², index)` on the GPU: indices and squared
/// distances, `queries.len() × min(k, n)` each (as [`crate::spatial::knn_batch`]; bitwise equal).
/// `k ≤` [`KMAX`].
pub fn knn_batch(g: &GpuContext, tree: &KdTree, queries: &[[f64; 3]], k: usize) -> Result<(Vec<u32>, Vec<f64>), PatError> {
    if k == 0 || k > KMAX {
        return Err(PatError::InvalidOption { what: format!("GPU kNN with k = {k} (1..={KMAX})") });
    }
    if let Some(i) = queries.iter().position(|q| !q.iter().all(|c| c.is_finite())) {
        return Err(PatError::NonFinite { index: i, what: "query" });
    }
    let t = GpuTree::upload(g, tree)?;
    let kk = k.min(tree.len());
    let mut idx = Vec::with_capacity(queries.len() * kk);
    let mut d2 = Vec::with_capacity(queries.len() * kk);
    let f = &g.k.knn_batch[kslot(k)];
    for qs in queries.chunks(QUERY_GROUP) {
        let nq = qs.len() as u64;
        let dq = g.up(&flat3(qs), 0.0)?;
        let mut oi = g.zeros::<u32>(qs.len() * kk)?;
        let mut od = g.zeros::<f64>(qs.len() * kk)?;
        let ki = k as i32;
        let mut b = g.stream.launch_builder(f);
        b.arg(&nq).arg(&dq).arg(&t.n).arg(&ki);
        tree_args!(b, t);
        b.arg(&mut oi).arg(&mut od);
        unsafe { b.launch(launch_cfg(nq)?) }.map_err(gpu_err("knn"))?;
        idx.extend(g.down(&oi, qs.len() * kk)?);
        d2.extend(g.down(&od, qs.len() * kk)?);
    }
    Ok((idx, d2))
}

/// The fitting neighborhoods on the GPU (as [`crate::spatial::fitting_neighborhoods`]; bitwise
/// equal). `k ≤` [`KMAX`].
pub fn fitting_neighborhoods(g: &GpuContext, tree: &KdTree, points: &[[f64; 3]], k: usize) -> Result<Neighborhoods, PatError> {
    let n = points.len();
    if n != tree.len() {
        return Err(PatError::SizeMismatch { what: format!("{n} points for a tree of {}", tree.len()) });
    }
    if k > KMAX {
        return Err(PatError::InvalidOption { what: format!("GPU fitting neighborhoods with k = {k} (≤ {KMAX})") });
    }
    let k = k.min(n).max(1);
    let t = GpuTree::upload(g, tree)?;
    let dp = g.up(&flat3(points), 0.0)?;
    let mut rows = g.zeros::<u32>(n * k)?;
    let (nn, ki) = (n as u32, k as i32);
    let mut b = g.stream.launch_builder(&g.k.fitting[kslot(k)]);
    b.arg(&nn).arg(&dp).arg(&ki);
    tree_args!(b, t);
    b.arg(&mut rows);
    unsafe { b.launch(launch_cfg(n as u64)?) }.map_err(gpu_err("fitting neighborhoods"))?;
    Ok(Neighborhoods { k, indices: g.down(&rows, n * k)? })
}

/// Fit one torus per point on the GPU (as [`crate::fit::fit_tori_with`]; bitwise equal).
pub fn fit_tori(
    g: &GpuContext,
    cloud: &PointCloud,
    coefficients: &[[f32; 6]],
    sigma: &[f64],
    tokens_lt_k: usize,
    opts: &FitOptions,
) -> Result<(Tori, FitReport), PatError> {
    let t0 = std::time::Instant::now();
    let n = cloud.len();
    if coefficients.len() != n {
        return Err(PatError::SizeMismatch { what: format!("{} coefficient rows for {n} points", coefficients.len()) });
    }
    if sigma.len() != n {
        return Err(PatError::SizeMismatch { what: format!("{} sigma values for {n} points", sigma.len()) });
    }
    if !(opts.eig_eps.is_finite() && opts.eig_eps > 0.0) {
        return Err(PatError::InvalidOption { what: format!("FitOptions::eig_eps = {}", opts.eig_eps) });
    }
    if !(opts.curvature_floor.is_finite() && opts.curvature_floor > 0.0) {
        return Err(PatError::InvalidOption { what: format!("FitOptions::curvature_floor = {}", opts.curvature_floor) });
    }
    if let Some(i) = coefficients.iter().position(|a| !a.iter().all(|x| x.is_finite())) {
        return Err(PatError::NonFinite { index: i, what: "coefficient row" });
    }
    let a: Vec<f32> = coefficients.iter().flatten().copied().collect();
    let dp = g.up(&flat3(cloud.positions()), 0.0)?;
    let dn = g.up(&flat3(cloud.normals()), 0.0)?;
    let da = g.up(&a, 0.0)?;
    let mut c = g.zeros::<f64>(3 * n)?;
    let mut u = g.zeros::<f64>(3 * n)?;
    let mut big = g.zeros::<f64>(n)?;
    let mut small = g.zeros::<f64>(n)?;
    let mut cl = g.zeros::<u32>(n)?;
    let nn = n as u32;
    let geometric = i32::from(opts.sign_rule == SignRule::Geometric);
    let robust = i32::from(opts.eigen_rule == EigenRule::Robust);
    unsafe {
        g.stream
            .launch_builder(&g.k.fit_tori)
            .arg(&nn)
            .arg(&dp)
            .arg(&dn)
            .arg(&da)
            .arg(&opts.eig_eps)
            .arg(&opts.curvature_floor)
            .arg(&geometric)
            .arg(&robust)
            .arg(&mut c)
            .arg(&mut u)
            .arg(&mut big)
            .arg(&mut small)
            .arg(&mut cl)
            .launch(launch_cfg(n as u64)?)
    }
    .map_err(gpu_err("fit"))?;
    let (c, u, big, small, cl) = (g.down(&c, 3 * n)?, g.down(&u, 3 * n)?, g.down(&big, n)?, g.down(&small, n)?, g.down(&cl, n)?);
    let mut report = FitReport { tokens_lt_k, ..Default::default() };
    let mut tori = Vec::with_capacity(n);
    for i in 0..n {
        let f = cl[i];
        report.clamped_disc += (f & 1 != 0) as usize;
        report.umbilic += (f & 1 != 0) as usize;
        report.clamped_eigvec += (f & 2 != 0) as usize;
        report.clamped_denominator += (f & 4 != 0) as usize;
        report.clamped_curvature += (f & 8 != 0) as usize;
        report.saddle += (f & 16 != 0) as usize;
        report.zero_curvature += (f & 32 != 0) as usize;
        report.eigen_switched += (f & 64 != 0) as usize;
        let t = Torus { center: [c[3 * i], c[3 * i + 1], c[3 * i + 2]], axis: [u[3 * i], u[3 * i + 1], u[3 * i + 2]], major: big[i], minor: small[i] };
        let finite = t.center.iter().chain(t.axis.iter()).all(|x| x.is_finite()) && t.major.is_finite() && t.minor.is_finite();
        if !finite {
            report.nonfinite.push(i as u32);
        }
        tori.push(t);
    }
    report.seconds = t0.elapsed().as_secs_f64();
    if !report.nonfinite.is_empty() {
        return Err(PatError::NonFiniteTorus { indices: report.nonfinite });
    }
    Ok((Tori { tori, coefficients: coefficients.to_vec(), sigma: sigma.to_vec() }, report))
}

/// What the blend kernels compute (`WANT_*` of `kernels.cu`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Want {
    Distance = 0,
    Frozen = 1,
    Exact = 2,
    Laplacian = 3,
}

impl Want {
    pub(crate) fn of(gradient: bool, laplacian: bool, mode: GradientMode) -> Want {
        if laplacian {
            Want::Laplacian
        } else if gradient {
            match mode {
                GradientMode::Frozen => Want::Frozen,
                GradientMode::Exact => Want::Exact,
            }
        } else {
            Want::Distance
        }
    }
}

/// The f32 data of the stable form (SPEC.md §7.3, NOTES-GPU.md), relative to the expansion point
/// `p_i` of every torus. [formula: pat-torus-sdf-stable]
struct StableData {
    origin: [f64; 3],
    pos: CudaSlice<f32>,
    w: CudaSlice<f32>,
    axis: CudaSlice<f32>,
    vp: CudaSlice<f32>,
    d1p: CudaSlice<f32>,
    d2p: CudaSlice<f32>,
    k2: CudaSlice<f32>,
    minor: CudaSlice<f32>,
    area: Option<CudaSlice<f32>>,
}

/// Offsets, counts, distances and indices of the radius sets of one query group, `r²`, `fallback_k`.
type RadiusScratch = (CudaSlice<u64>, CudaSlice<u32>, CudaSlice<f64>, CudaSlice<u32>, f64, i32);

/// A field's data on the device.
pub(crate) struct GpuField {
    g: Arc<GpuContext>,
    tree: GpuTree,
    n: u32,
    pos: CudaSlice<f64>,
    center: CudaSlice<f64>,
    axis: CudaSlice<f64>,
    major: CudaSlice<f64>,
    minor: CudaSlice<f64>,
    area: Option<CudaSlice<f64>>,
    stable: Option<StableData>,
}

impl std::fmt::Debug for GpuField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuField").field("device", &self.g.info).field("n", &self.n).field("f32_stable", &self.stable.is_some()).finish()
    }
}

/// The stable form's per-torus constants, in f64 on the host, relative to `p`: `w = p − c`,
/// `v_p = ‖w × u‖`, `d_p = (v_p − R, ⟨w, u⟩)`, `K₂ = ‖d_p‖² − r²` (= `(‖d_p‖ − |r|)(‖d_p‖ + |r|)`,
/// no cancellation).
///
/// [formula: pat-torus-sdf-stable]
pub fn stable_constants(p: [f64; 3], t: &Torus) -> ([f64; 3], f64, f64, f64, f64) {
    let w = [p[0] - t.center[0], p[1] - t.center[1], p[2] - t.center[2]];
    let u = t.axis;
    let cr = [w[1] * u[2] - w[2] * u[1], w[2] * u[0] - w[0] * u[2], w[0] * u[1] - w[1] * u[0]];
    let vp = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
    let d1 = vp - t.major;
    let d2 = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
    let nd = (d1 * d1 + d2 * d2).sqrt();
    let r = t.minor.abs();
    (w, vp, d1, d2, (nd - r) * (nd + r))
}

impl GpuField {
    pub(crate) fn new(
        g: Arc<GpuContext>,
        tree: &KdTree,
        positions: &[[f64; 3]],
        tori: &[Torus],
        area: Option<&[f64]>,
        f32_stable: bool,
    ) -> Result<GpuField, PatError> {
        let n = positions.len();
        let tr = GpuTree::upload(&g, tree)?;
        let center: Vec<f64> = tori.iter().flat_map(|t| t.center).collect();
        let axis: Vec<f64> = tori.iter().flat_map(|t| t.axis).collect();
        let major: Vec<f64> = tori.iter().map(|t| t.major).collect();
        let minor: Vec<f64> = tori.iter().map(|t| t.minor).collect();
        let stable = if f32_stable {
            // origin: the center of the bounding box (f32 coordinates relative to it)
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for p in positions {
                for a in 0..3 {
                    lo[a] = lo[a].min(p[a]);
                    hi[a] = hi[a].max(p[a]);
                }
            }
            let origin = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
            let mut pos = Vec::with_capacity(3 * n);
            let (mut w, mut ax) = (Vec::with_capacity(3 * n), Vec::with_capacity(3 * n));
            let (mut vp, mut d1p, mut d2p, mut k2, mut mi) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
            for (p, t) in positions.iter().zip(tori) {
                let (wv, v, d1, d2, k) = stable_constants(*p, t);
                for a in 0..3 {
                    pos.push((p[a] - origin[a]) as f32);
                    w.push(wv[a] as f32);
                    ax.push(t.axis[a] as f32);
                }
                vp.push(v as f32);
                d1p.push(d1 as f32);
                d2p.push(d2 as f32);
                k2.push(k as f32);
                mi.push(t.minor as f32);
            }
            let area = match area {
                Some(a) => Some(g.up(&a.iter().map(|x| *x as f32).collect::<Vec<_>>(), 0.0)?),
                None => None,
            };
            Some(StableData {
                origin,
                pos: g.up(&pos, 0.0)?,
                w: g.up(&w, 0.0)?,
                axis: g.up(&ax, 0.0)?,
                vp: g.up(&vp, 0.0)?,
                d1p: g.up(&d1p, 0.0)?,
                d2p: g.up(&d2p, 0.0)?,
                k2: g.up(&k2, 0.0)?,
                minor: g.up(&mi, 0.0)?,
                area,
            })
        } else {
            None
        };
        let area = match area {
            Some(a) => Some(g.up(a, 0.0)?),
            None => None,
        };
        Ok(GpuField {
            tree: tr,
            n: n as u32,
            pos: g.up(&flat3(positions), 0.0)?,
            center: g.up(&center, 0.0)?,
            axis: g.up(&axis, 0.0)?,
            major: g.up(&major, 0.0)?,
            minor: g.up(&minor, 0.0)?,
            area,
            stable,
            g,
        })
    }

    pub(crate) fn device(&self) -> Device {
        self.g.device()
    }

    /// Whether this neighborhood runs on the GPU (`Knn` with `k ≤ KMAX` or `k ≥ n`, `All`,
    /// `Radius` with `fallback_k ≤ KMAX`).
    pub(crate) fn supports(nb: Neighborhood, n: usize) -> bool {
        match nb {
            Neighborhood::Knn { k } => k <= KMAX || k >= n,
            Neighborhood::All => true,
            Neighborhood::Radius { fallback_k, .. } => fallback_k <= KMAX,
        }
    }

    /// Batch evaluation (the checks of `TorusField::evaluate` are done by the caller). Returns
    /// the report, or `Err(index)` of the first degenerate query.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate(
        &self,
        xs: &[[f64; 3]],
        nb: Neighborhood,
        rule: LambdaRule,
        lambda_global: f64,
        isovalue: f64,
        want: Want,
        out: &mut EvalOut<'_>,
    ) -> Result<Result<EvalReport, usize>, PatError> {
        let g = &*self.g;
        let n = self.n as usize;
        let (ri, c, lc) = match rule {
            LambdaRule::PerQuery { c } => (0i32, c, 0.0),
            LambdaRule::Constant { lambda } => (1, 0.0, lambda),
            LambdaRule::PaperHeuristic { c, .. } => (2, c, lambda_global),
        };
        let wi = want as i32;
        let want_grad = !matches!(want, Want::Distance);
        let want_lap = matches!(want, Want::Laplacian);
        let mut rep = EvalReport { queries: xs.len(), ..Default::default() };
        let null = 0u64;
        for (gi, qs) in xs.chunks(QUERY_GROUP).enumerate() {
            let base = gi * QUERY_GROUP;
            let m = qs.len();
            let nq = m as u64;
            let mut dist = g.zeros::<f64>(m)?;
            let mut grad = if want_grad { Some(g.zeros::<f64>(3 * m)?) } else { None };
            let mut lap = if want_lap { Some(g.zeros::<f64>(m)?) } else { None };
            let mut flags = g.zeros::<u32>(5 * m)?;
            // queries: f64 (f64 path; also the f32 path's kNN), and f32 relative to the origin
            let dq = g.up(&flat3(qs), 0.0)?;
            let dq32 = match &self.stable {
                Some(s) => {
                    let o = s.origin;
                    let v: Vec<f32> = qs.iter().flat_map(|q| [(q[0] - o[0]) as f32, (q[1] - o[1]) as f32, (q[2] - o[2]) as f32]).collect();
                    Some(g.up(&v, 0.0)?)
                }
                None => None,
            };
            // radius scratch
            let mut radius: Option<RadiusScratch> = None;
            if let Neighborhood::Radius { r, fallback_k } = nb {
                let r2 = r * r;
                let mut counts = g.zeros::<u32>(m)?;
                let mut b = g.stream.launch_builder(&g.k.radius_counts);
                b.arg(&nq).arg(&dq).arg(&r2);
                tree_args!(b, self.tree);
                b.arg(&mut counts);
                unsafe { b.launch(launch_cfg(nq)?) }.map_err(gpu_err("radius count"))?;
                let hc = g.down(&counts, m)?;
                let mut offs = Vec::with_capacity(m);
                let mut tot = 0u64;
                for &cnt in &hc {
                    offs.push(tot);
                    tot += cnt as u64;
                }
                let need = tot as usize * 12;
                if need > g.info.total_memory as usize / 2 {
                    return Err(PatError::Gpu { reason: format!("radius sets of {tot} entries ({need} bytes) do not fit; use fewer queries per call") });
                }
                radius = Some((g.up(&offs, 0)?, counts, g.zeros::<f64>(tot as usize)?, g.zeros::<u32>(tot as usize)?, r2, fallback_k as i32));
            }
            macro_rules! common_args {
                ($b:expr, $q:expr) => {
                    $b.arg(&nq).arg($q).arg(&self.n).arg(&ri).arg(&c).arg(&lc).arg(&isovalue).arg(&wi)
                };
            }
            macro_rules! outs {
                ($b:expr) => {{
                    $b.arg(&mut dist);
                    match grad.as_mut() {
                        Some(x) => $b.arg(x),
                        None => $b.arg(&null),
                    };
                    match lap.as_mut() {
                        Some(x) => $b.arg(x),
                        None => $b.arg(&null),
                    };
                    $b.arg(&mut flags);
                }};
            }
            let all = matches!(nb, Neighborhood::All) || matches!(nb, Neighborhood::Knn { k } if k >= n);
            match (&self.stable, dq32.as_ref()) {
                (None, _) => {
                    let f = if all {
                        &g.k.eval_all
                    } else {
                        match nb {
                            Neighborhood::Knn { k } => &g.k.eval_knn[kslot(k)],
                            _ => &g.k.eval_radius,
                        }
                    };
                    let mut b = g.stream.launch_builder(f);
                    let ki = match nb {
                        Neighborhood::Knn { k } => k as i32,
                        _ => 0,
                    };
                    if !all {
                        if let Some((offs, counts, sd, si, r2, fk)) = radius.as_mut() {
                            b.arg(&*r2).arg(&*fk).arg(&*offs).arg(&*counts).arg(sd).arg(si);
                        } else {
                            b.arg(&ki);
                        }
                    }
                    common_args!(b, &dq);
                    b.arg(&self.pos).arg(&self.center).arg(&self.axis).arg(&self.major).arg(&self.minor);
                    match &self.area {
                        Some(a) => b.arg(a),
                        None => b.arg(&null),
                    };
                    tree_args!(b, self.tree);
                    outs!(b);
                    unsafe { b.launch(launch_cfg(nq)?) }.map_err(gpu_err("blend"))?;
                }
                (Some(s), Some(dq32)) => {
                    let f = if all {
                        &g.k.eval_all_f32
                    } else {
                        match nb {
                            Neighborhood::Knn { k } => &g.k.eval_knn_f32[kslot(k)],
                            _ => &g.k.eval_radius_f32,
                        }
                    };
                    let ki = match nb {
                        Neighborhood::Knn { k } => k as i32,
                        _ => 0,
                    };
                    let mut b = g.stream.launch_builder(f);
                    if !all {
                        if let Some((offs, counts, sd, si, r2, fk)) = radius.as_mut() {
                            b.arg(&*r2).arg(&*fk).arg(&*offs).arg(&*counts).arg(sd).arg(si);
                        } else {
                            b.arg(&ki);
                        }
                    }
                    common_args!(b, &dq);
                    b.arg(dq32).arg(&s.pos).arg(&s.w).arg(&s.axis).arg(&s.vp).arg(&s.d1p).arg(&s.d2p).arg(&s.k2).arg(&s.minor);
                    match &s.area {
                        Some(a) => b.arg(a),
                        None => b.arg(&null),
                    };
                    tree_args!(b, self.tree);
                    outs!(b);
                    unsafe { b.launch(launch_cfg(nq)?) }.map_err(gpu_err("blend (f32)"))?;
                }
                (Some(_), None) => unreachable!("f32 queries are uploaded with the stable data"),
            }
            let hf = g.down(&flags, 5 * m)?;
            for (j, f) in hf.chunks_exact(5).enumerate() {
                if f[3] >> 31 != 0 {
                    return Ok(Err(base + j));
                }
                rep.singular_axis += f[0] as u64;
                rep.singular_core += f[1] as u64;
                rep.at_data_point += f[2] as u64;
                rep.radius_fallbacks += (f[3] & 1) as u64;
                rep.terms += f[4] as u64;
            }
            if let Some(d) = out.distance.as_deref_mut() {
                d[base..base + m].copy_from_slice(&g.down(&dist, m)?);
            } else {
                // the distance buffer is always computed; nothing to copy
            }
            if let (Some(o), Some(dg)) = (out.gradient.as_deref_mut(), grad.as_ref()) {
                let h = g.down(dg, 3 * m)?;
                for (j, v) in h.chunks_exact(3).enumerate() {
                    o[base + j] = [v[0], v[1], v[2]];
                }
            }
            if let (Some(o), Some(dl)) = (out.laplacian.as_deref_mut(), lap.as_ref()) {
                o[base..base + m].copy_from_slice(&g.down(dl, m)?);
            }
        }
        Ok(Ok(rep))
    }
}

/// Whether [`crate::PointsAsTori::build`] computes the fitting neighborhoods on the GPU when the
/// precompute device is a GPU (NOTES-GPU.md: measured).
pub const USE_GPU_KNN: bool = true;
/// Whether [`crate::PointsAsTori::build`] fits the tori on the GPU when the precompute device is a
/// GPU (NOTES-GPU.md: measured).
pub const USE_GPU_FIT: bool = true;
