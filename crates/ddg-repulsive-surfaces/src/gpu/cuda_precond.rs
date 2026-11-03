//! The preconditioner of `src/solve/precond.rs` on the device (GPU-2, `NOTES-GPU2.md`): the `L′`
//! solves by the chain schedule of [`crate::solve::ldlt::LdltSchedule`] (kernels `ldl_*`), the
//! vertex-pair fractional product (GPU-1's tree, on device vectors) and the augmented middle factor
//! (`prec_aug`), with one upload and one download per application.

use std::sync::Arc;

use cudarc::driver::{CudaSlice, LaunchConfig, PushKernelArg};

use super::GpuContext;
use super::cuda::{CudaBackend, DeviceVertexFractional, cfg, gpu_err};
use crate::error::RsError;
use crate::solve::ldlt::{LdltSchedule, SimplicialLdlt, Sweep};

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

struct DeviceSweep {
    level_ptr: Vec<usize>,
    level_block: Vec<u32>,
    slot_ptr: Vec<usize>,
    chains: CudaSlice<u32>,
    slot_row: CudaSlice<u32>,
    slot_chunk: CudaSlice<u32>,
    slot_len: CudaSlice<u32>,
    nfull: CudaSlice<u32>,
    nmix: CudaSlice<u32>,
    slot_base: CudaSlice<u32>,
}

impl DeviceSweep {
    fn upload(g: &CudaBackend, s: &Sweep) -> Result<Self, RsError> {
        Ok(DeviceSweep {
            level_ptr: s.level_ptr.clone(),
            level_block: s.level_block.clone(),
            slot_ptr: s.slot_ptr.clone(),
            chains: g.upload(&s.chains, 0)?,
            slot_row: g.upload(&s.slot_row, 0)?,
            slot_chunk: g.upload(&s.slot_chunk, 0)?,
            slot_len: g.upload(&s.slot_len, 0)?,
            nfull: g.upload(&s.nfull, 0)?,
            nmix: g.upload(&s.nmix, 0)?,
            slot_base: g.upload(&s.slot_base, 0)?,
        })
    }
}

/// A [`SimplicialLdlt`] with its schedule on the device.
pub(crate) struct DeviceLdlt {
    n: u32,
    per_axis: bool,
    nv: u32,
    perm: CudaSlice<u32>,
    rp: CudaSlice<u32>,
    rj: CudaSlice<u32>,
    rx: CudaSlice<f64>,
    lp: CudaSlice<u32>,
    ci: CudaSlice<u32>,
    cx: CudaSlice<f64>,
    d: CudaSlice<f64>,
    cstart: CudaSlice<u32>,
    cend: CudaSlice<u32>,
    nsub: u32,
    sub_start: CudaSlice<u32>,
    sub_end: CudaSlice<u32>,
    fwd: DeviceSweep,
    bwd: DeviceSweep,
    max_slots: usize,
}

impl DeviceLdlt {
    pub(crate) fn upload(
        g: &CudaBackend,
        f: &SimplicialLdlt,
        s: &LdltSchedule,
        per_axis: bool,
        nv: usize,
    ) -> Result<Self, RsError> {
        let (rp, rj, rx) = f.rows();
        let (lp, _, _) = f.columns();
        Ok(DeviceLdlt {
            n: f.dim() as u32,
            per_axis,
            nv: nv as u32,
            perm: g.upload(f.perm(), 0)?,
            rp: g.upload(&to_u32(rp, "L row pointers")?, 0)?,
            rj: g.upload(rj, 0)?,
            rx: g.upload(rx, 0.0)?,
            lp: g.upload(&to_u32(lp, "L column pointers")?, 0)?,
            ci: g.upload(&s.ci, 0)?,
            cx: g.upload(&s.cx, 0.0)?,
            d: g.upload(f.diagonal(), 0.0)?,
            cstart: g.upload(&s.chain_start, 0)?,
            cend: g.upload(&s.chain_end, 0)?,
            nsub: s.sub_start.len() as u32,
            sub_start: g.upload(&s.sub_start, 0)?,
            sub_end: g.upload(&s.sub_end, 0)?,
            fwd: DeviceSweep::upload(g, &s.fwd)?,
            bwd: DeviceSweep::upload(g, &s.bwd)?,
            max_slots: s.n_slots.0.max(s.n_slots.1),
        })
    }

    /// `dst = L′⁻¹ src` for `nb` interleaved vectors of dimension `dim` (the order of
    /// [`SimplicialLdlt::solve_rows`] on the gathered columns).
    pub(crate) fn solve(
        &self,
        g: &CudaBackend,
        src: &CudaSlice<f64>,
        nb: usize,
        dim: u64,
        dst: &mut CudaSlice<f64>,
        mut timing: Option<&mut [f64; 4]>,
    ) -> Result<(), RsError> {
        let st = &g.stream;
        let e = || gpu_err("preconditioner solve");
        // Diagnostics (`timing`): synchronize after every launch and add its wall time to
        // [gather + scatter, partials, forward chains, backward chains].
        let mut t0 = std::time::Instant::now();
        let mut tick = |slot: usize| -> Result<(), RsError> {
            if let Some(t) = timing.as_deref_mut() {
                st.synchronize().map_err(gpu_err("preconditioner solve"))?;
                let dt = t0.elapsed().as_secs_f64();
                t[slot] += dt;
                t0 = std::time::Instant::now();
            }
            Ok(())
        };
        tick(0)?;
        let m = if self.per_axis { 3 * nb as u32 } else { nb as u32 };
        let total = self.n as u64 * m as u64;
        let pa = u32::from(self.per_axis);
        let mut x = g.alloc(total as usize)?;
        let mut z = g.alloc(total as usize)?;
        let mut part = g.alloc(self.max_slots * m as usize)?;
        unsafe {
            st.launch_builder(&g.k.ldl_gather)
                .arg(&self.n)
                .arg(&m)
                .arg(&pa)
                .arg(&self.nv)
                .arg(&dim)
                .arg(&self.perm)
                .arg(src)
                .arg(&mut x)
                .launch(cfg(total)?)
        }
        .map_err(e())?;
        tick(0)?;
        for fwd in [true, false] {
            let (sw, ptr, idx, val) = if fwd {
                (&self.fwd, &self.rp, &self.rj, &self.rx)
            } else {
                (&self.bwd, &self.lp, &self.ci, &self.cx)
            };
            let fwd_u = u32::from(fwd);
            let leaves = |data: &mut CudaSlice<f64>, zz: &mut CudaSlice<f64>| -> Result<(), RsError> {
                if self.nsub == 0 {
                    return Ok(());
                }
                unsafe {
                    st.launch_builder(&g.k.ldl_leaves)
                        .arg(&fwd_u)
                        .arg(&self.nsub)
                        .arg(&m)
                        .arg(&self.sub_start)
                        .arg(&self.sub_end)
                        .arg(ptr)
                        .arg(idx)
                        .arg(val)
                        .arg(&self.d)
                        .arg(data)
                        .arg(zz)
                        .launch(cfg(self.nsub as u64 * m as u64)?)
                }
                .map_err(gpu_err("preconditioner solve"))?;
                Ok(())
            };
            if fwd {
                leaves(&mut x, &mut z)?;
                tick(0)?;
            }
            for l in 0..sw.level_ptr.len() - 1 {
                let count = (sw.slot_ptr[l + 1] - sw.slot_ptr[l]) as u32;
                let data = if fwd { &mut x } else { &mut z };
                if count > 0 {
                    let off = sw.slot_ptr[l] as u32;
                    unsafe {
                        st.launch_builder(&g.k.ldl_partials)
                            .arg(&count)
                            .arg(&m)
                            .arg(&off)
                            .arg(&sw.slot_row)
                            .arg(&sw.slot_chunk)
                            .arg(&sw.slot_len)
                            .arg(ptr)
                            .arg(idx)
                            .arg(val)
                            .arg(&*data)
                            .arg(&mut part)
                            .launch(cfg(count as u64 * m as u64)?)
                    }
                    .map_err(e())?;
                    tick(1)?;
                }
                let nch = (sw.level_ptr[l + 1] - sw.level_ptr[l]) as u32;
                if nch == 0 {
                    continue;
                }
                let chains = sw.chains.slice(sw.level_ptr[l]..sw.level_ptr[l + 1]);
                let lc = LaunchConfig {
                    grid_dim: (nch, m, 1),
                    block_dim: (sw.level_block[l], 1, 1),
                    shared_mem_bytes: 0,
                };
                if fwd {
                    unsafe {
                        st.launch_builder(&g.k.ldl_chains)
                            .arg(&fwd_u)
                            .arg(&m)
                            .arg(&chains)
                            .arg(&self.cstart)
                            .arg(&self.cend)
                            .arg(ptr)
                            .arg(idx)
                            .arg(val)
                            .arg(&sw.nfull)
                            .arg(&sw.nmix)
                            .arg(&sw.slot_base)
                            .arg(&part)
                            .arg(&self.d)
                            .arg(&mut x)
                            .arg(&mut z)
                            .launch(lc)
                    }
                    .map_err(e())?;
                    tick(2)?;
                } else {
                    unsafe {
                        st.launch_builder(&g.k.ldl_chains)
                            .arg(&fwd_u)
                            .arg(&m)
                            .arg(&chains)
                            .arg(&self.cstart)
                            .arg(&self.cend)
                            .arg(ptr)
                            .arg(idx)
                            .arg(val)
                            .arg(&sw.nfull)
                            .arg(&sw.nmix)
                            .arg(&sw.slot_base)
                            .arg(&part)
                            .arg(&self.d)
                            .arg(&mut z)
                            .arg(&mut x)
                            .launch(lc)
                    }
                    .map_err(e())?;
                    tick(3)?;
                }
            }
            if !fwd {
                leaves(&mut z, &mut x)?;
                tick(3)?;
            }
        }
        unsafe {
            st.launch_builder(&g.k.ldl_scatter)
                .arg(&self.n)
                .arg(&m)
                .arg(&pa)
                .arg(&self.nv)
                .arg(&dim)
                .arg(&self.perm)
                .arg(&z)
                .arg(dst)
                .launch(cfg(total)?)
        }
        .map_err(e())?;
        tick(0)?;
        Ok(())
    }
}

/// The host data of [`DevicePreconditioner::new`].
pub(crate) struct PrecHost<'a> {
    pub factor: &'a SimplicialLdlt,
    pub schedule: &'a LdltSchedule,
    pub per_axis: bool,
    pub n_vertices: usize,
    pub n_constraints: usize,
    /// The Laplacian block `K` by columns (symmetric).
    pub k_ptr: &'a [usize],
    pub k_idx: &'a [u32],
    pub k_val: &'a [f64],
    /// Augmented form (else congruence).
    pub augmented: bool,
    pub scale: f64,
}

/// The fine preconditioner `P⁻¹` (without the coarse correction) on the device.
pub(crate) struct DevicePreconditioner {
    gpu: Arc<GpuContext>,
    ldlt: DeviceLdlt,
    nv: u32,
    dim: usize,
    kp: CudaSlice<u32>,
    ki: CudaSlice<u32>,
    kx: CudaSlice<f64>,
    augmented: bool,
    pub(crate) scale: f64,
    pub(crate) bytes: usize,
}

impl DevicePreconditioner {
    pub(crate) fn new(gpu: &Arc<GpuContext>, h: &PrecHost<'_>) -> Result<Self, RsError> {
        let g = &gpu.inner;
        let ldlt = DeviceLdlt::upload(g, h.factor, h.schedule, h.per_axis, h.n_vertices)?;
        let (rp, _, _) = h.factor.rows();
        let nnz = rp[rp.len() - 1];
        let bytes = 2 * 12 * nnz + 16 * h.factor.dim() + 12 * h.k_val.len();
        Ok(DevicePreconditioner {
            gpu: gpu.clone(),
            ldlt,
            nv: h.n_vertices as u32,
            dim: 3 * h.n_vertices + h.n_constraints,
            kp: g.upload(&to_u32(h.k_ptr, "K column pointers")?, 0)?,
            ki: g.upload(h.k_idx, 0)?,
            kx: g.upload(h.k_val, 0.0)?,
            augmented: h.augmented,
            scale: h.scale,
            bytes,
        })
    }

    /// Diagnostics: seconds of one `L′` solve of `rs` split as
    /// [gather + scatter, partials, forward chains, backward chains] (synchronized after every
    /// launch, so the sum exceeds an unsynchronized solve), and of the unsynchronized solve.
    pub(crate) fn profile_solve(&self, rs: &[&[f64]]) -> Result<([f64; 4], f64), RsError> {
        let g = &self.gpu.inner;
        let nb = rs.len();
        let dim = self.dim;
        let mut h = Vec::with_capacity(nb * dim);
        for r in rs {
            h.extend_from_slice(r);
        }
        let dr = g.upload(&h, 0.0)?;
        let mut y = g.alloc(nb * dim)?;
        let mut t = [0.0; 4];
        self.ldlt.solve(g, &dr, nb, dim as u64, &mut y, Some(&mut t))?;
        g.stream.synchronize().map_err(gpu_err("profile"))?;
        let t0 = std::time::Instant::now();
        self.ldlt.solve(g, &dr, nb, dim as u64, &mut y, None)?;
        g.stream.synchronize().map_err(gpu_err("profile"))?;
        Ok((t, t0.elapsed().as_secs_f64()))
    }

    /// `P⁻¹ r` for every `r` of `rs` (the fine part of `SparsePreconditioner::apply_many`).
    pub(crate) fn apply_many(&self, frac: &DeviceVertexFractional, rs: &[&[f64]]) -> Result<Vec<Vec<f64>>, RsError> {
        let g = &self.gpu.inner;
        let st = &g.stream;
        let nb = rs.len();
        if nb == 0 {
            return Ok(Vec::new());
        }
        let dim = self.dim;
        let mut h = Vec::with_capacity(nb * dim);
        for r in rs {
            h.extend_from_slice(r);
        }
        let dr = g.upload(&h, 0.0)?;
        let total = (nb * dim) as u64;
        let mut y = g.alloc(nb * dim)?;
        self.ldlt.solve(g, &dr, nb, dim as u64, &mut y, None)?;
        let mut v = g.alloc(nb * dim)?;
        unsafe {
            st.launch_builder(&g.k.vec_zero)
                .arg(&total)
                .arg(&mut v)
                .launch(cfg(total)?)
        }
        .map_err(gpu_err("preconditioner"))?;
        frac.apply_device(&y, dim as u64, nb, &mut v)?;
        if self.augmented {
            let nbu = nb as u32;
            let dimu = dim as u64;
            let scaled = u32::from(self.scale != 1.0);
            unsafe {
                st.launch_builder(&g.k.prec_aug)
                    .arg(&self.nv)
                    .arg(&nbu)
                    .arg(&dimu)
                    .arg(&scaled)
                    .arg(&self.scale)
                    .arg(&self.kp)
                    .arg(&self.ki)
                    .arg(&self.kx)
                    .arg(&y)
                    .arg(&mut v)
                    .launch(cfg(nb as u64 * 3 * self.nv as u64)?)
            }
            .map_err(gpu_err("preconditioner"))?;
        }
        let mut z = g.alloc(nb * dim)?;
        self.ldlt.solve(g, &v, nb, dim as u64, &mut z, None)?;
        if self.augmented {
            unsafe {
                st.launch_builder(&g.k.vec_add)
                    .arg(&total)
                    .arg(&y)
                    .arg(&mut z)
                    .launch(cfg(total)?)
            }
            .map_err(gpu_err("preconditioner"))?;
        }
        let mut out = vec![0.0; nb * dim];
        g.download(&z, &mut out)?;
        Ok(out.chunks(dim).map(<[f64]>::to_vec).collect())
    }
}
