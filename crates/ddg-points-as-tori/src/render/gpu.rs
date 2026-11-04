//! The CUDA tracer (feature `cuda`): rays are traced by the kernels `trace_knn_*` / `trace_all` of
//! `src/gpu/kernels.cu`, which repeat [`super::trace_ray`] operation for operation on the f64
//! blend (no FMA, crate-owned `exp`, one ray per thread). Shading stays on the host
//! ([`super::Tracer::render_with`]), so only the hit positions and gradients depend on the device.
//!
//! Supported: a single [`crate::TorusField`] on a CUDA device with a kNN (`k ≤ 128`) or
//! all-points neighborhood; the level `ℓ` and offsets are options of the tracer / the field's
//! isovalue. Booleans (a [`crate::csg::Csg`] of several fields) render on the CPU only.

use super::{Backend, RayJob, RayResult, TraceParams};
use crate::{Device, PatError, TorusField};

/// The CUDA backend over one field.
pub struct GpuBackend<'a> {
    field: &'a TorusField,
}

impl<'a> GpuBackend<'a> {
    /// A backend for `field` (which must have been built with a CUDA evaluation device).
    pub fn new(field: &'a TorusField) -> Result<GpuBackend<'a>, PatError> {
        if matches!(field.device(), Device::Cpu) {
            return Err(PatError::Gpu {
                reason: "the field was built for the CPU; set EvalOptions::device to a CUDA device"
                    .into(),
            });
        }
        Ok(GpuBackend { field })
    }
}

impl Backend for GpuBackend<'_> {
    fn trace(&self, jobs: &[RayJob], p: &TraceParams) -> Result<Vec<RayResult>, PatError> {
        let flat: Vec<f64> = jobs
            .iter()
            .flat_map(|j| {
                [
                    j.ro[0], j.ro[1], j.ro[2], j.rd[0], j.rd[1], j.rd[2], j.t0, j.t1,
                ]
            })
            .collect();
        let (out, st) = self.field.gpu_trace(&flat, p)?;
        Ok((0..jobs.len())
            .map(|i| RayResult {
                status: st[2 * i],
                steps: st[2 * i + 1],
                t: out[4 * i],
                grad: [out[4 * i + 1], out[4 * i + 2], out[4 * i + 3]],
            })
            .collect())
    }

    fn phi(&self, pts: &[[f64; 3]]) -> Result<Vec<f64>, PatError> {
        use crate::csg::ImplicitField;
        let mut out = vec![0.0; pts.len()];
        self.field.distances(pts, &mut out)?;
        Ok(out)
    }

    fn device(&self) -> Device {
        self.field.device()
    }

    fn threads(&self) -> usize {
        1
    }
}
