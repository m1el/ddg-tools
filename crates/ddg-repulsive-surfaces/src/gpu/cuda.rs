//! The CUDA implementation (cargo feature `cuda`): `cudarc` with the driver API and NVRTC, both
//! opened at run time. Kernels: `kernels.cu` (module docs of [`crate::gpu`] for the
//! determinism argument).

use std::sync::Arc;

use cudarc::driver::{CudaContext, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};

use super::{CudaOptions, Determinism, GpuContext};
use crate::error::RsError;
use crate::spatial::block_tree::BlockTreeParts;

const SOURCE: &str = include_str!("kernels.cu");
const BLOCK: u32 = 256;

fn gpu_err(operation: &'static str) -> impl FnOnce(cudarc::driver::DriverError) -> RsError {
    move |e| RsError::Gpu {
        operation,
        detail: format!("{e:?}"),
    }
}

fn cfg(n: u64) -> Result<LaunchConfig, RsError> {
    let blocks = n.div_ceil(BLOCK as u64);
    if blocks > i32::MAX as u64 {
        return Err(RsError::Gpu {
            operation: "launch",
            detail: format!("{n} threads exceed the grid limit"),
        });
    }
    Ok(LaunchConfig {
        grid_dim: (blocks.max(1) as u32, 1, 1),
        block_dim: (BLOCK, 1, 1),
        shared_mem_bytes: 0,
    })
}

/// The compiled kernels.
pub(crate) struct Kernels {
    up_leaves: CudaFunction,
    up_internal: CudaFunction,
    far: CudaFunction,
    rows: CudaFunction,
    face_m: CudaFunction,
    face_z: CudaFunction,
    face_gather: CudaFunction,
    vertex_psi: CudaFunction,
    vertex_out: CudaFunction,
    pairwise_level: CudaFunction,
    mul_elems: CudaFunction,
}

/// A device, its stream and the kernels.
pub(crate) struct CudaBackend {
    #[allow(dead_code)]
    ctx: Arc<CudaContext>,
    pub(crate) stream: Arc<CudaStream>,
    k: Kernels,
    pub(crate) name: String,
    pub(crate) determinism: Determinism,
}

impl CudaBackend {
    pub(crate) fn new(o: &CudaOptions) -> Result<CudaBackend, RsError> {
        let ctx = std::panic::catch_unwind(|| CudaContext::new(o.ordinal))
            .map_err(|_| RsError::Gpu {
                operation: "context",
                detail: "the CUDA driver library (libcuda.so) could not be loaded".into(),
            })?
            .map_err(gpu_err("context"))?;
        // One stream, used in order: the per-buffer event tracking of cudarc is not needed.
        // SAFETY: every buffer of this back end is used on `stream` only.
        unsafe { ctx.disable_event_tracking() };
        let (major, minor) = ctx.compute_capability().map_err(gpu_err("context"))?;
        let name = format!(
            "{} (compute capability {major}.{minor})",
            ctx.name().map_err(gpu_err("context"))?
        );
        let opts = CompileOptions {
            fmad: Some(o.determinism == Determinism::RunToRun),
            prec_div: Some(true),
            prec_sqrt: Some(true),
            ftz: Some(false),
            options: vec![format!("--gpu-architecture=sm_{major}{minor}")],
            name: Some("ddg_repulsive_surfaces_kernels.cu".into()),
            ..Default::default()
        };
        let ptx = std::panic::catch_unwind(|| compile_ptx_with_opts(SOURCE, opts))
            .map_err(|_| RsError::Gpu {
                operation: "compile",
                detail: "NVRTC (libnvrtc.so.12) could not be loaded".into(),
            })?
            .map_err(|e| RsError::Gpu {
                operation: "compile",
                detail: format!("{e:?}"),
            })?;
        let module = ctx.load_module(ptx).map_err(gpu_err("module"))?;
        let f = |n: &str| module.load_function(n).map_err(gpu_err("module"));
        let k = Kernels {
            up_leaves: f("bct_up_leaves")?,
            up_internal: f("bct_up_internal")?,
            far: f("bct_far")?,
            rows: f("bct_rows")?,
            face_m: f("face_m")?,
            face_z: f("face_z")?,
            face_gather: f("face_gather")?,
            vertex_psi: f("vertex_psi")?,
            vertex_out: f("vertex_out")?,
            pairwise_level: f("pairwise_level")?,
            mul_elems: f("mul_elems")?,
        };
        let stream = ctx.default_stream();
        Ok(CudaBackend {
            ctx,
            stream,
            k,
            name,
            determinism: o.determinism,
        })
    }

    /// Uploads `v` (a one-element `dummy` if empty: zero-size allocations are not allowed).
    fn upload<T: DeviceRepr + Copy>(&self, v: &[T], dummy: T) -> Result<CudaSlice<T>, RsError> {
        if v.is_empty() {
            return self.stream.clone_htod(&[dummy]).map_err(gpu_err("upload"));
        }
        self.stream.clone_htod(v).map_err(gpu_err("upload"))
    }

    fn alloc(&self, n: usize) -> Result<CudaSlice<f64>, RsError> {
        // SAFETY: every buffer allocated here is fully written by a kernel before it is read.
        unsafe { self.stream.alloc::<f64>(n.max(1)) }.map_err(gpu_err("alloc"))
    }

    fn download(&self, d: &CudaSlice<f64>, out: &mut [f64]) -> Result<(), RsError> {
        self.stream.memcpy_dtoh(d, out).map_err(gpu_err("download"))?;
        self.stream.synchronize().map_err(gpu_err("download"))
    }

    /// The pairwise-tree dot product (`super::pairwise_dot`).
    pub(crate) fn dot(&self, a: &[f64], b: &[f64]) -> Result<f64, RsError> {
        if a.is_empty() {
            return Ok(0.0);
        }
        let n = a.len() as u64;
        let da = self.upload(a, 0.0)?;
        let db = self.upload(b, 0.0)?;
        let mut cur = self.alloc(a.len())?;
        unsafe {
            self.stream
                .launch_builder(&self.k.mul_elems)
                .arg(&n)
                .arg(&da)
                .arg(&db)
                .arg(&mut cur)
                .launch(cfg(n)?)
        }
        .map_err(gpu_err("dot"))?;
        let s = self.tree_sum(cur, a.len())?;
        Ok(s)
    }

    /// The fixed-shape pairwise sum of the first `n` entries of `cur` (consumed).
    fn tree_sum(&self, mut cur: CudaSlice<f64>, n: usize) -> Result<f64, RsError> {
        let mut len = n as u64;
        while len > 1 {
            let n_out = len.div_ceil(2);
            let mut next = self.alloc(n_out as usize)?;
            unsafe {
                self.stream
                    .launch_builder(&self.k.pairwise_level)
                    .arg(&len)
                    .arg(&n_out)
                    .arg(&cur)
                    .arg(&mut next)
                    .launch(cfg(n_out)?)
            }
            .map_err(gpu_err("reduce"))?;
            cur = next;
            len = n_out;
        }
        let mut out = [0.0f64];
        self.stream.memcpy_dtoh(&cur.slice(0..1), &mut out).map_err(gpu_err("download"))?;
        self.stream.synchronize().map_err(gpu_err("download"))?;
        Ok(out[0])
    }
}

/// A [`crate::spatial::block_tree::BlockTree`] on the device.
pub(crate) struct DeviceBlockTree {
    n: u32,
    n_nodes: u32,
    nk: u32,
    n_leaves: u32,
    start: CudaSlice<u32>,
    end: CudaSlice<u32>,
    right: CudaSlice<u32>,
    leaves: CudaSlice<u32>,
    /// Internal nodes grouped by depth, deepest group first, and the group boundaries.
    levels: CudaSlice<u32>,
    level_ptr: Vec<usize>,
    leaf_of: CudaSlice<u32>,
    chain_ptr: CudaSlice<u32>,
    chain: CudaSlice<u32>,
    weight: CudaSlice<f64>,
    far_ptr: CudaSlice<u32>,
    far_col: CudaSlice<u32>,
    far_coef: CudaSlice<f64>,
    near_ptr: CudaSlice<u32>,
    near_col: CudaSlice<u32>,
    near_off: CudaSlice<u64>,
    near_val: CudaSlice<f64>,
    /// Bytes on the device.
    pub(crate) bytes: usize,
}

fn to_u32(v: &[usize], what: &'static str) -> Result<Vec<u32>, RsError> {
    v.iter()
        .map(|&x| {
            u32::try_from(x).map_err(|_| RsError::Gpu {
                operation: "upload",
                detail: format!("{what}: {x} does not fit 32 bits"),
            })
        })
        .collect()
}

impl DeviceBlockTree {
    pub(crate) fn upload<const N: usize>(g: &CudaBackend, p: &BlockTreeParts<'_, N>, n: usize) -> Result<Self, RsError> {
        let n_nodes = p.start.len();
        // Depths (pre-order: children after parents) and the internal nodes per depth.
        let mut depth = vec![0u32; n_nodes];
        let mut max_depth = 0;
        for a in 0..n_nodes {
            let r = p.right[a] as usize;
            if r != 0 {
                depth[a + 1] = depth[a] + 1;
                depth[r] = depth[a] + 1;
                max_depth = max_depth.max(depth[a] + 1);
            }
        }
        let mut levels = Vec::new();
        let mut level_ptr = vec![0usize];
        for d in (0..=max_depth).rev() {
            levels.extend((0..n_nodes as u32).filter(|&a| p.right[a as usize] != 0 && depth[a as usize] == d));
            level_ptr.push(levels.len());
        }
        let mut leaf_of = vec![0u32; n];
        for (li, &l) in p.leaves.iter().enumerate() {
            for q in p.start[l as usize]..p.end[l as usize] {
                leaf_of[q as usize] = li as u32;
            }
        }
        let far_coef: Vec<f64> = p.far_coef.iter().flatten().copied().collect();
        let near_val: Vec<f64> = p.near_val.iter().flatten().copied().collect();
        let near_off: Vec<u64> = p.near_off.iter().map(|&x| x as u64).collect();
        let bytes = 8 * (far_coef.len() + near_val.len() + near_off.len() + n)
            + 4 * (5 * n_nodes + p.chain.len() + p.far_col.len() + p.near_col.len() + 2 * n);
        Ok(DeviceBlockTree {
            n: n as u32,
            n_nodes: n_nodes as u32,
            nk: N as u32,
            n_leaves: p.leaves.len() as u32,
            start: g.upload(p.start, 0)?,
            end: g.upload(p.end, 0)?,
            right: g.upload(p.right, 0)?,
            leaves: g.upload(p.leaves, 0)?,
            levels: g.upload(&levels, 0)?,
            level_ptr,
            leaf_of: g.upload(&leaf_of, 0)?,
            chain_ptr: g.upload(&to_u32(p.chain_ptr, "chain pointer")?, 0)?,
            chain: g.upload(p.chain, 0)?,
            weight: g.upload(p.weight, 0.0)?,
            far_ptr: g.upload(&to_u32(p.far_ptr, "far pointer")?, 0)?,
            far_col: g.upload(p.far_col, 0)?,
            far_coef: g.upload(&far_coef, 0.0)?,
            near_ptr: g.upload(&to_u32(p.near_ptr, "near pointer")?, 0)?,
            near_col: g.upload(p.near_col, 0)?,
            near_off: g.upload(&near_off, 0)?,
            near_val: g.upload(&near_val, 0.0)?,
            bytes,
        })
    }

    /// `H̃ ψ` for `C` channels (channels `< h` on kernel 0, the others on kernel `N − 1`), all in
    /// BVH order (`BlockTree::multiply_split`).
    pub(crate) fn multiply(&self, g: &CudaBackend, c: u32, h: u32, psi: &CudaSlice<f64>) -> Result<CudaSlice<f64>, RsError> {
        let s = &g.stream;
        let e = gpu_err("product");
        let cc = c as u64;
        let mut up = g.alloc(self.n_nodes as usize * c as usize)?;
        unsafe {
            s.launch_builder(&g.k.up_leaves)
                .arg(&self.n_leaves)
                .arg(&c)
                .arg(&self.leaves)
                .arg(&self.start)
                .arg(&self.end)
                .arg(&self.weight)
                .arg(psi)
                .arg(&mut up)
                .launch(cfg(self.n_leaves as u64 * cc)?)
        }
        .map_err(e)?;
        for w in self.level_ptr.windows(2) {
            let count = (w[1] - w[0]) as u32;
            if count == 0 {
                continue;
            }
            let nodes = self.levels.slice(w[0]..w[1]);
            unsafe {
                s.launch_builder(&g.k.up_internal)
                    .arg(&count)
                    .arg(&c)
                    .arg(&nodes)
                    .arg(&self.right)
                    .arg(&mut up)
                    .launch(cfg(count as u64 * cc)?)
            }
            .map_err(gpu_err("product"))?;
        }
        let mut far = g.alloc(self.n_nodes as usize * c as usize)?;
        unsafe {
            s.launch_builder(&g.k.far)
                .arg(&self.n_nodes)
                .arg(&c)
                .arg(&h)
                .arg(&self.nk)
                .arg(&self.far_ptr)
                .arg(&self.far_col)
                .arg(&self.far_coef)
                .arg(&up)
                .arg(&mut far)
                .launch(cfg(self.n_nodes as u64 * cc)?)
        }
        .map_err(gpu_err("product"))?;
        let mut out = g.alloc(self.n as usize * c as usize)?;
        unsafe {
            s.launch_builder(&g.k.rows)
                .arg(&self.n)
                .arg(&c)
                .arg(&h)
                .arg(&self.nk)
                .arg(&self.leaf_of)
                .arg(&self.chain_ptr)
                .arg(&self.chain)
                .arg(&self.start)
                .arg(&self.end)
                .arg(&self.near_ptr)
                .arg(&self.near_col)
                .arg(&self.near_off)
                .arg(&self.near_val)
                .arg(&self.weight)
                .arg(&far)
                .arg(psi)
                .arg(&mut out)
                .launch(cfg(self.n as u64 * cc)?)
        }
        .map_err(gpu_err("product"))?;
        Ok(out)
    }
}

/// The host data of a face metric tree for its upload.
pub(crate) struct FaceMetricHost<'a> {
    pub order: &'a [u32],
    pub position: &'a [u32],
    pub faces: &'a [[u32; 3]],
    pub df: &'a [[[f64; 3]; 3]],
    pub areas: &'a [f64],
    pub ha_perm: &'a [[f64; 2]],
    pub inc_ptr: &'a [usize],
    pub inc: &'a [(u32, u8)],
    pub n_vertices: usize,
    pub d_hi: &'a [f64],
    pub d_lo: &'a [f64],
}

/// [`crate::metric::hierarchical::FaceMetricTree`] on the device.
pub(crate) struct DeviceFaceMetric {
    gpu: Arc<GpuContext>,
    tree: DeviceBlockTree,
    nf: u32,
    nv: u32,
    order: CudaSlice<u32>,
    position: CudaSlice<u32>,
    faces: CudaSlice<u32>,
    df: CudaSlice<f64>,
    areas: CudaSlice<f64>,
    ha_perm: CudaSlice<f64>,
    d_hi: CudaSlice<f64>,
    d_lo: CudaSlice<f64>,
    has_d_hi: bool,
    has_d_lo: bool,
    inc_ptr: CudaSlice<u32>,
    inc_face: CudaSlice<u32>,
    inc_corner: CudaSlice<u8>,
}

impl DeviceFaceMetric {
    pub(crate) fn new(gpu: &Arc<GpuContext>, tree: &BlockTreeParts<'_, 2>, h: &FaceMetricHost<'_>) -> Result<Self, RsError> {
        let g = &gpu.inner;
        let nf = h.faces.len();
        let faces: Vec<u32> = h.faces.iter().flatten().copied().collect();
        let df: Vec<f64> = h.df.iter().flatten().flatten().copied().collect();
        let ha: Vec<f64> = h.ha_perm.iter().flatten().copied().collect();
        let inc_face: Vec<u32> = h.inc.iter().map(|e| e.0).collect();
        let inc_corner: Vec<u8> = h.inc.iter().map(|e| e.1).collect();
        Ok(DeviceFaceMetric {
            gpu: gpu.clone(),
            tree: DeviceBlockTree::upload(g, tree, nf)?,
            nf: nf as u32,
            nv: h.n_vertices as u32,
            order: g.upload(h.order, 0)?,
            position: g.upload(h.position, 0)?,
            faces: g.upload(&faces, 0)?,
            df: g.upload(&df, 0.0)?,
            areas: g.upload(h.areas, 0.0)?,
            ha_perm: g.upload(&ha, 0.0)?,
            d_hi: g.upload(h.d_hi, 0.0)?,
            d_lo: g.upload(h.d_lo, 0.0)?,
            has_d_hi: !h.d_hi.is_empty(),
            has_d_lo: !h.d_lo.is_empty(),
            inc_ptr: g.upload(&to_u32(h.inc_ptr, "incidence pointer")?, 0)?,
            inc_face: g.upload(&inc_face, 0)?,
            inc_corner: g.upload(&inc_corner, 0)?,
        })
    }

    pub(crate) fn n_vertices(&self) -> usize {
        self.nv as usize
    }

    /// Bytes of the block tree on the device.
    pub(crate) fn tree_bytes(&self) -> usize {
        self.tree.bytes
    }

    /// `FaceMetricTree::apply_batch` with the uploaded diagonal.
    pub(crate) fn apply_batch(&self, w_hi: f64, w_lo: f64, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        let g = &self.gpu.inner;
        let s = &g.stream;
        let nb = xs.len();
        if nb == 0 {
            return Ok(());
        }
        let nd = 3 * self.nv as usize;
        let hi = w_hi != 0.0 || self.has_d_hi;
        let lo = w_lo != 0.0 || self.has_d_lo;
        if !hi && !lo {
            ys.iter_mut().for_each(|y| y.iter_mut().for_each(|v| *v = 0.0));
            return Ok(());
        }
        let c = (nb * (9 * hi as usize + 3 * lo as usize)) as u32;
        let off = if hi { 9 * nb as u32 } else { 0 };
        let mut xh = Vec::with_capacity(nb * nd);
        for x in xs {
            xh.extend_from_slice(x);
        }
        let dx = g.upload(&xh, 0.0)?;
        let total = self.nf as u64 * c as u64;
        let mut m = g.alloc(total as usize)?;
        let mut psi = g.alloc(total as usize)?;
        let ndu = nd as u64;
        unsafe {
            s.launch_builder(&g.k.face_m)
                .arg(&self.nf)
                .arg(&c)
                .arg(&off)
                .arg(&ndu)
                .arg(&self.order)
                .arg(&self.faces)
                .arg(&self.df)
                .arg(&self.areas)
                .arg(&dx)
                .arg(&mut m)
                .arg(&mut psi)
                .launch(cfg(total)?)
        }
        .map_err(gpu_err("product"))?;
        let mut hz = self.tree.multiply(g, c, off, &psi)?;
        let (has_hi, has_lo) = (self.has_d_hi as u32, self.has_d_lo as u32);
        unsafe {
            s.launch_builder(&g.k.face_z)
                .arg(&self.nf)
                .arg(&c)
                .arg(&off)
                .arg(&w_hi)
                .arg(&w_lo)
                .arg(&self.order)
                .arg(&self.ha_perm)
                .arg(&self.d_hi)
                .arg(&self.d_lo)
                .arg(&has_hi)
                .arg(&has_lo)
                .arg(&m)
                .arg(&mut hz)
                .launch(cfg(total)?)
        }
        .map_err(gpu_err("product"))?;
        let mut dy = g.alloc(nb * nd)?;
        let (nbu, hiu, lou) = (nb as u32, hi as u32, lo as u32);
        unsafe {
            s.launch_builder(&g.k.face_gather)
                .arg(&self.nv)
                .arg(&nbu)
                .arg(&c)
                .arg(&off)
                .arg(&hiu)
                .arg(&lou)
                .arg(&ndu)
                .arg(&self.inc_ptr)
                .arg(&self.inc_face)
                .arg(&self.inc_corner)
                .arg(&self.areas)
                .arg(&self.df)
                .arg(&self.position)
                .arg(&hz)
                .arg(&mut dy)
                .launch(cfg(self.nv as u64 * nb as u64)?)
        }
        .map_err(gpu_err("product"))?;
        let mut yh = vec![0.0; nb * nd];
        g.download(&dy, &mut yh)?;
        for (q, y) in ys.iter_mut().enumerate() {
            y.copy_from_slice(&yh[q * nd..(q + 1) * nd]);
        }
        Ok(())
    }
}

/// [`crate::metric::hierarchical::VertexFractionalTree`] on the device.
pub(crate) struct DeviceVertexFractional {
    gpu: Arc<GpuContext>,
    tree: DeviceBlockTree,
    n: u32,
    order: CudaSlice<u32>,
    position: CudaSlice<u32>,
    mass: CudaSlice<f64>,
    ha: CudaSlice<f64>,
}

impl DeviceVertexFractional {
    pub(crate) fn new(
        gpu: &Arc<GpuContext>,
        tree: &BlockTreeParts<'_, 1>,
        order: &[u32],
        position: &[u32],
        mass: &[f64],
        ha: &[f64],
    ) -> Result<Self, RsError> {
        let g = &gpu.inner;
        Ok(DeviceVertexFractional {
            gpu: gpu.clone(),
            tree: DeviceBlockTree::upload(g, tree, mass.len())?,
            n: mass.len() as u32,
            order: g.upload(order, 0)?,
            position: g.upload(position, 0)?,
            mass: g.upload(mass, 0.0)?,
            ha: g.upload(ha, 0.0)?,
        })
    }

    pub(crate) fn n_vertices(&self) -> usize {
        self.n as usize
    }

    /// `VertexFractionalTree::apply_batch`.
    pub(crate) fn apply_batch(&self, xs: &[&[f64]], ys: &mut [&mut [f64]]) -> Result<(), RsError> {
        let g = &self.gpu.inner;
        let s = &g.stream;
        let nb = xs.len();
        if nb == 0 {
            return Ok(());
        }
        let nd = 3 * self.n as usize;
        let c = 3 * nb as u32;
        let mut xh = Vec::with_capacity(nb * nd);
        for x in xs {
            xh.extend_from_slice(x);
        }
        let dx = g.upload(&xh, 0.0)?;
        let total = self.n as u64 * c as u64;
        let mut psi = g.alloc(total as usize)?;
        let ndu = nd as u64;
        unsafe {
            s.launch_builder(&g.k.vertex_psi)
                .arg(&self.n)
                .arg(&c)
                .arg(&ndu)
                .arg(&self.order)
                .arg(&self.mass)
                .arg(&dx)
                .arg(&mut psi)
                .launch(cfg(total)?)
        }
        .map_err(gpu_err("product"))?;
        let hx = self.tree.multiply(g, c, c, &psi)?;
        let mut dy = g.alloc(nb * nd)?;
        let nbu = nb as u32;
        unsafe {
            s.launch_builder(&g.k.vertex_out)
                .arg(&self.n)
                .arg(&nbu)
                .arg(&c)
                .arg(&ndu)
                .arg(&self.position)
                .arg(&self.mass)
                .arg(&self.ha)
                .arg(&dx)
                .arg(&hx)
                .arg(&mut dy)
                .launch(cfg(nb as u64 * nd as u64)?)
        }
        .map_err(gpu_err("product"))?;
        let mut yh = vec![0.0; nb * nd];
        g.download(&dy, &mut yh)?;
        for (q, y) in ys.iter_mut().enumerate() {
            y.copy_from_slice(&yh[q * nd..(q + 1) * nd]);
        }
        Ok(())
    }
}
