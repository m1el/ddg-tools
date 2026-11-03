//! The GPU back end of the hierarchical products (`NOTES-GPU.md`; optional, cargo feature `cuda`).
//!
//! What runs on the GPU: the block-cluster-tree products of the metric (`B̃`, `B̃⁰` with the
//! obstacle diagonal, [`crate::metric::hierarchical::FaceMetricTree`]) and of the vertex-pair
//! fractional operator of the preconditioner
//! ([`crate::metric::hierarchical::VertexFractionalTree`]). The trees are built on the CPU (as
//! before) and uploaded once per step; everything else (GMRES, the sparse `L′` solves, energy,
//! line search, remeshing) stays on the CPU. The CPU path is the reference and is unchanged.
//!
//! **Choosing the device** is explicit: [`Device`] in
//! [`crate::solve::iterative::HierarchicalOptions::device`] (default [`Device::Cpu`]).
//! [`Device::Cuda`] fails with [`RsError::Gpu`] if the crate was built without the `cuda` feature
//! or no usable device/driver/NVRTC is found; [`Device::Auto`] uses the GPU when one is usable and
//! the CPU otherwise. A [`GpuContext`] (driver context, compiled kernels) is created by
//! [`Device::resolve`] and passed on explicitly ([`crate::Session`] creates one per session);
//! there is no global state.
//!
//! **Determinism** ([`Determinism`]). Every GPU kernel computes one output value per thread in the
//! operation order of the CPU loop it replaces, with gathers only (no atomics, no
//! launch-shape-dependent reductions), so results are bitwise identical from run to run on any
//! device ([`Determinism::RunToRun`], level 1, the default: the compiler may contract `a·b + c`
//! into FMAs, so results differ from the CPU's by round-off). [`Determinism::MatchCpu`]
//! (level 2) compiles the kernels with `--fmad=false`, which makes every product bitwise equal to
//! the CPU's (the CPU loops are compiled without FMA, `block_tree::simd`), and therefore a whole
//! flow bitwise equal to the CPU flow.

use std::sync::Arc;

use crate::error::RsError;

#[cfg(feature = "cuda")]
pub(crate) mod cuda;

/// How bitwise reproducible GPU results are (module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Determinism {
    /// Level 1: bitwise identical from run to run (any device, any launch configuration); FMA
    /// contraction allowed, so within round-off of the CPU.
    #[default]
    RunToRun,
    /// Level 2: additionally bitwise equal to the CPU path (`--fmad=false`).
    MatchCpu,
}

/// Options of the CUDA device.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CudaOptions {
    /// The device ordinal (default 0).
    pub ordinal: usize,
    /// The determinism level (default [`Determinism::RunToRun`]).
    pub determinism: Determinism,
}

/// Where the hierarchical products run (module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Device {
    /// The CPU (the reference; default).
    #[default]
    Cpu,
    /// A CUDA device; an error if it is not usable.
    Cuda(CudaOptions),
    /// A CUDA device if one is usable, the CPU otherwise.
    Auto(CudaOptions),
}

impl Device {
    /// The GPU context this choice asks for: `None` for the CPU (or [`Device::Auto`] without a
    /// usable GPU), a new context otherwise (compiles the kernels: about a second).
    ///
    /// # Errors
    /// [`RsError::Gpu`] for [`Device::Cuda`] when the crate was built without the `cuda` feature
    /// or the device cannot be used.
    pub fn resolve(&self) -> Result<Option<Arc<GpuContext>>, RsError> {
        match self {
            Device::Cpu => Ok(None),
            Device::Cuda(o) => GpuContext::new(o).map(Some),
            Device::Auto(o) => Ok(GpuContext::new(o).ok()),
        }
    }

    /// Whether this is the CPU.
    pub fn is_cpu(&self) -> bool {
        matches!(self, Device::Cpu)
    }
}

/// Whether the crate was built with the `cuda` feature.
pub const fn cuda_compiled() -> bool {
    cfg!(feature = "cuda")
}

/// A CUDA context with the compiled kernels (module docs). Shared through `Arc`.
pub struct GpuContext {
    #[cfg(feature = "cuda")]
    pub(crate) inner: cuda::CudaBackend,
    #[cfg(not(feature = "cuda"))]
    _never: std::convert::Infallible,
}

impl std::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuContext")
            .field("device", &self.device_name())
            .field("determinism", &self.determinism())
            .finish()
    }
}

impl GpuContext {
    /// Opens the device `options.ordinal` and compiles the kernels for it.
    ///
    /// # Errors
    /// [`RsError::Gpu`]: built without the `cuda` feature, no driver, no such device, NVRTC
    /// missing or failing.
    pub fn new(options: &CudaOptions) -> Result<Arc<GpuContext>, RsError> {
        #[cfg(feature = "cuda")]
        {
            Ok(Arc::new(GpuContext {
                inner: cuda::CudaBackend::new(options)?,
            }))
        }
        #[cfg(not(feature = "cuda"))]
        {
            let _ = options;
            Err(RsError::Gpu {
                operation: "context",
                detail: "ddg-repulsive-surfaces was built without the `cuda` feature \
                         (cargo build --features cuda)"
                    .into(),
            })
        }
    }

    /// The device's name and compute capability.
    pub fn device_name(&self) -> String {
        #[cfg(feature = "cuda")]
        {
            self.inner.name.clone()
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }

    /// The determinism level the kernels were compiled for.
    pub fn determinism(&self) -> Determinism {
        #[cfg(feature = "cuda")]
        {
            self.inner.determinism
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }

    /// Sum of `a_i b_i` by the fixed-shape pairwise tree ([`pairwise_sum`]) on the device (a
    /// determinism check and a building block; the flow's GMRES keeps its vectors on the CPU).
    ///
    /// # Errors
    /// [`RsError::Gpu`] on a device failure; [`RsError::InvalidOption`] for different lengths.
    pub fn dot(&self, a: &[f64], b: &[f64]) -> Result<f64, RsError> {
        if a.len() != b.len() {
            return Err(RsError::InvalidOption {
                name: "gpu dot",
                value: format!("{} and {} entries", a.len(), b.len()),
                requirement: "vectors of equal length",
            });
        }
        #[cfg(feature = "cuda")]
        {
            self.inner.dot(a, b)
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }
}

/// The fixed-shape pairwise sum of `v` (the reduction of determinism levels 1–2): pad `v` with
/// `+0.0` to the next power of two `2^k` and add adjacent pairs level by level,
/// `s_{l+1}[i] = s_l[2i] + s_l[2i+1]`. The shape depends only on `v.len()`; the GPU computes the
/// same tree level by level (`pairwise_level` kernel), so the two are bitwise equal.
pub fn pairwise_sum(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut cur: Vec<f64> = v.to_vec();
    while cur.len() > 1 {
        let n_out = cur.len().div_ceil(2);
        let next: Vec<f64> = (0..n_out)
            .map(|i| {
                let a = cur[2 * i];
                let b = cur.get(2 * i + 1).copied().unwrap_or(0.0);
                a + b
            })
            .collect();
        cur = next;
    }
    cur[0]
}

/// `Σ a_i b_i` by [`pairwise_sum`] (products rounded first, then the tree).
pub fn pairwise_dot(a: &[f64], b: &[f64]) -> f64 {
    let p: Vec<f64> = a.iter().zip(b).map(|(x, y)| x * y).collect();
    pairwise_sum(&p)
}

/// The GPU copy of a [`crate::metric::hierarchical::FaceMetricTree`] with its obstacle diagonal
/// (built by [`crate::metric::hierarchical::FaceMetricTree::upload`]).
pub struct GpuFaceMetric {
    #[cfg(feature = "cuda")]
    pub(crate) inner: cuda::DeviceFaceMetric,
    #[cfg(not(feature = "cuda"))]
    pub(crate) _never: std::convert::Infallible,
}

impl std::fmt::Debug for GpuFaceMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("GpuFaceMetric")
    }
}

/// The GPU copy of a [`crate::metric::hierarchical::VertexFractionalTree`] (built by
/// [`crate::metric::hierarchical::VertexFractionalTree::upload`]).
pub struct GpuVertexFractional {
    #[cfg(feature = "cuda")]
    pub(crate) inner: cuda::DeviceVertexFractional,
    #[cfg(not(feature = "cuda"))]
    pub(crate) _never: std::convert::Infallible,
}

impl std::fmt::Debug for GpuVertexFractional {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("GpuVertexFractional")
    }
}

fn check_batch(xs: &[&[f64]], ys: &[&mut [f64]], n: usize) -> Result<(), RsError> {
    if xs.len() != ys.len() || xs.iter().any(|x| x.len() != n) || ys.iter().any(|y| y.len() != n) {
        return Err(RsError::InvalidOption {
            name: "gpu product batch",
            value: format!("{} inputs, {} outputs", xs.len(), ys.len()),
            requirement: "one output per input, 3|V| interleaved entries each",
        });
    }
    Ok(())
}

impl GpuFaceMetric {
    /// `|V|` of the uploaded tree.
    pub fn n_vertices(&self) -> usize {
        #[cfg(feature = "cuda")]
        {
            self.inner.n_vertices()
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }

    /// Bytes of the block cluster tree on the device.
    pub fn device_bytes(&self) -> usize {
        #[cfg(feature = "cuda")]
        {
            self.inner.tree_bytes()
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }

    /// [`crate::metric::hierarchical::FaceMetricTree::apply_batch`] on the device, with the
    /// diagonal given at the upload: `y = (w_hi B̃ + w_lo B̃⁰ + O_diag) x` for each vector.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths or non-finite weights; [`RsError::Gpu`].
    pub fn apply_batch(&self, w_hi: f64, w_lo: f64, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        check_batch(xs, ys, 3 * self.n_vertices())?;
        if !(w_hi.is_finite() && w_lo.is_finite()) {
            return Err(RsError::InvalidOption {
                name: "metric weights",
                value: format!("({w_hi}, {w_lo})"),
                requirement: "finite weights",
            });
        }
        #[cfg(feature = "cuda")]
        {
            self.inner.apply_batch(w_hi, w_lo, xs, ys)
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }
}

impl GpuVertexFractional {
    /// `|V|`.
    pub fn n_vertices(&self) -> usize {
        #[cfg(feature = "cuda")]
        {
            self.inner.n_vertices()
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }

    /// [`crate::metric::hierarchical::VertexFractionalTree::apply_batch`] on the device.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] for wrong lengths; [`RsError::Gpu`].
    pub fn apply_batch(&self, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        check_batch(xs, ys, 3 * self.n_vertices())?;
        #[cfg(feature = "cuda")]
        {
            self.inner.apply_batch(xs, ys)
        }
        #[cfg(not(feature = "cuda"))]
        {
            match self._never {}
        }
    }
}

/// A Barnes–Hut tree's frozen lists on the device ([`crate::energy::barnes_hut::BhTree::attach_gpu`]).
pub(crate) struct GpuBh {
    #[cfg(feature = "cuda")]
    pub(crate) inner: cuda::DeviceBh,
    #[cfg(not(feature = "cuda"))]
    _never: std::convert::Infallible,
}

impl std::fmt::Debug for GpuBh {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("GpuBh")
    }
}
