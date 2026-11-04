//! Sphere tracer and shading (Fig. 12, SPEC.md §8.5). WP11.
//!
//! [`Tracer::render`] renders any [`ImplicitField`] (a [`crate::TorusField`], a [`crate::csg::Csg`] of
//! offsets and Booleans, an analytic field) from an explicit [`Camera`]: the level set `φ = ℓ` is
//! sphere-traced, shaded with the demo's Lambert + IBL-GGX + Fresnel terms, over a white
//! background with the demo's shadow-only ground plane and an optional cut plane coloured by `φ`.
//! The shading is a port of the authors' demo shader (`demo/shaders_3d.py` of `points-as-tori`,
//! MIT, Copyright (c) 2026 Nicole Feng; see `shade.rs`, `cubemap.rs`, `assets/`).
//!
//! # Traversal
//!
//! `[formula: pat-render-step]` Rays are clipped to the inflated bounding sphere (centre of the
//! field's bounding box, radius `0.8 · diag + max(ℓ, 0)` as the demo). Then
//! `t ← t + safety · |φ(x) − ℓ| / max(‖∇φ(x)‖, g_min)`, a hit when `|φ − ℓ| < ε_hit`, at most
//! `max_steps` evaluations. For an `L`-Lipschitz `φ` with `L ≤ max(‖∇φ‖, g_min)` along the ray
//! (the safety bound, SPEC D17) the step never crosses the level set, since `|φ(x+sd) − φ(x)| ≤ L s`.
//!
//! # Determinism
//!
//! Every pixel is computed independently in f64 with a fixed operation order, so the image and
//! the report counters are bitwise identical at any thread count and from run to run
//! (`tests/render.rs`). The CUDA path (`gpu.rs`, kernel `trace_*` in `src/gpu/kernels.cu`) traces
//! with the CPU's step arithmetic (no FMA, crate-owned `exp`); shading stays on the host.

pub mod cubemap;
#[cfg(feature = "cuda")]
pub mod gpu;
mod shade;

pub use cubemap::CubeMap;

use crate::PatError;
use crate::csg::ImplicitField;
use rayon::prelude::*;
use shade::{
    Lights, V3, add, cross, cutplane_color, dot, length, normalize, scale, shade_surface, sub,
};

/// The largest image side and pixel count accepted.
pub const MAX_SIDE: usize = 16384;
/// The largest number of pixels accepted (`2^28`).
pub const MAX_PIXELS: usize = 1 << 28;

fn bad(what: impl Into<String>) -> PatError {
    PatError::InvalidOption { what: what.into() }
}

fn finite3(v: V3) -> bool {
    v.iter().all(|x| x.is_finite())
}

/// A pinhole camera.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Eye position.
    pub eye: [f64; 3],
    /// Point looked at.
    pub target: [f64; 3],
    /// Up direction (not parallel to the view direction).
    pub up: [f64; 3],
    /// Vertical field of view in radians, in `(0, π)`.
    pub fov_y: f64,
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
}

impl Camera {
    /// `look_at(eye, target, up)` with a vertical field of view in degrees; validated.
    pub fn look_at(
        eye: [f64; 3],
        target: [f64; 3],
        up: [f64; 3],
        fov_y_degrees: f64,
        width: usize,
        height: usize,
    ) -> Result<Camera, PatError> {
        let c = Camera {
            eye,
            target,
            up,
            fov_y: fov_y_degrees.to_radians(),
            width,
            height,
        };
        c.validate()?;
        Ok(c)
    }

    /// The demo's camera for a scene of the given centre and scale (`scale = 1` for the demo's
    /// `[-1, 1]³` normalization): eye at `(10 cos 5.4, 2.4, 10 sin 5.4)·scale`, looking at the
    /// centre, y up, focal length 4.5 (`fov_y = 2 atan(1/4.5)` = 25.1°).
    pub fn demo(
        center: [f64; 3],
        scale: f64,
        width: usize,
        height: usize,
    ) -> Result<Camera, PatError> {
        let eye = [
            center[0] + scale * 10.0 * 5.4f64.cos(),
            center[1] + scale * 2.4,
            center[2] + scale * 10.0 * 5.4f64.sin(),
        ];
        let c = Camera {
            eye,
            target: center,
            up: [0.0, 1.0, 0.0],
            fov_y: 2.0 * (1.0f64 / 4.5).atan(),
            width,
            height,
        };
        c.validate()?;
        Ok(c)
    }

    /// A camera orbiting `center` at `distance`: `azimuth` about the `up` axis (radians, from the
    /// `+x`/`+z` plane's x direction when `up` is y), `elevation` above the horizon.
    #[allow(clippy::too_many_arguments)]
    pub fn orbit(
        center: [f64; 3],
        distance: f64,
        azimuth: f64,
        elevation: f64,
        up_axis: usize,
        fov_y_degrees: f64,
        width: usize,
        height: usize,
    ) -> Result<Camera, PatError> {
        if up_axis > 2 {
            return Err(bad(format!("up_axis = {up_axis} (0, 1 or 2)")));
        }
        let (a, b) = ((up_axis + 1) % 3, (up_axis + 2) % 3);
        let mut dir = [0.0; 3];
        dir[a] = elevation.cos() * azimuth.cos();
        dir[b] = elevation.cos() * azimuth.sin();
        dir[up_axis] = elevation.sin();
        let mut up = [0.0; 3];
        up[up_axis] = 1.0;
        let eye = add(center, scale(dir, distance));
        Camera::look_at(eye, center, up, fov_y_degrees, width, height)
    }

    /// Check the camera.
    pub fn validate(&self) -> Result<(), PatError> {
        if !(finite3(self.eye) && finite3(self.target) && finite3(self.up)) {
            return Err(bad("camera eye, target and up must be finite"));
        }
        if !(self.fov_y.is_finite() && self.fov_y > 0.0 && self.fov_y < std::f64::consts::PI) {
            return Err(bad(format!(
                "camera fov_y = {} rad (must be in (0, pi))",
                self.fov_y
            )));
        }
        if self.width == 0
            || self.height == 0
            || self.width > MAX_SIDE
            || self.height > MAX_SIDE
            || self.width * self.height > MAX_PIXELS
        {
            return Err(bad(format!(
                "resolution {}x{} (each side 1..={MAX_SIDE}, at most {MAX_PIXELS} pixels)",
                self.width, self.height
            )));
        }
        let w = sub(self.target, self.eye);
        if length(w) <= 0.0 {
            return Err(bad("camera eye equals target"));
        }
        if length(cross(normalize(w), self.up)) < 1e-9 * length(self.up).max(1e-300)
            || length(self.up) == 0.0
        {
            return Err(bad("camera up is zero or parallel to the view direction"));
        }
        Ok(())
    }

    /// Orthonormal frame `(right, up, forward)` of the camera (`setCamera` of the demo).
    fn basis(&self) -> (V3, V3, V3) {
        let cw = normalize(sub(self.target, self.eye));
        let cu = normalize(cross(cw, self.up));
        let cv = cross(cu, cw);
        (cu, cv, cw)
    }
}

/// How a distance is given: relative to the bounding-box diagonal or absolute.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Epsilon {
    /// A fraction of the bounding-box diagonal.
    Relative(f64),
    /// In world units.
    Absolute(f64),
}

impl Epsilon {
    fn resolve(self, diag: f64) -> f64 {
        match self {
            Epsilon::Relative(r) => r * diag,
            Epsilon::Absolute(a) => a,
        }
    }
    fn validate(self, what: &str) -> Result<(), PatError> {
        let v = match self {
            Epsilon::Relative(r) | Epsilon::Absolute(r) => r,
        };
        if v.is_finite() && v > 0.0 {
            Ok(())
        } else {
            Err(bad(format!("{what} = {self:?} (finite and > 0)")))
        }
    }
}

/// The environment of the IBL terms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Environment {
    /// The demo's blurred forest cube map (`assets/ForestBlurred*.png`).
    Forest,
    /// A uniform sky of this colour.
    Constant([f64; 3]),
}

/// The ground plane of the demo: invisible except for the shadow it receives (a white background
/// pixel whose ground point is shadowed is drawn at `shadow_value`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ground {
    /// Normal (the up direction; the demo's is `+y`).
    pub up: [f64; 3],
    /// Distance below the lowest bounding-box corner, in units of the scene scale `s = diag/(2√3)`
    /// (the demo's 1.0 in its `[-1, 1]³` normalization).
    pub gap: f64,
    /// Colour factor of a shadowed ground pixel (0.8).
    pub shadow_value: f64,
}

impl Default for Ground {
    fn default() -> Self {
        Ground {
            up: [0.0, 1.0, 0.0],
            gap: 1.0,
            shadow_value: 0.8,
        }
    }
}

/// The cut plane of the demo: an axis-aligned rectangle (rotated about two of its axes) coloured
/// by `φ − ℓ` with the demo's colour map and isobands.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutPlane {
    /// Axis of the unrotated normal (0, 1, 2).
    pub axis: usize,
    /// Offset of the plane along that axis from the bounding-box centre, in world units.
    pub offset: f64,
    /// Rotation about the first in-plane axis (radians).
    pub rotation_x: f64,
    /// Rotation about the second in-plane axis (radians).
    pub rotation_y: f64,
    /// The colour map covers `φ − ℓ ∈ [−range, range]`; `None`: `2 s` (the demo's 2).
    pub range: Option<f64>,
    /// Isoband frequency (`sin(−f φ)`); `None`: `60/s` (the demo's 60).
    pub isoband_frequency: Option<f64>,
    /// Saturate `φ − ℓ` into the colour-map range (default). The demo's map is black outside
    /// `[−range, range)`; `false` reproduces that.
    pub clamp_colormap: bool,
}

impl Default for CutPlane {
    fn default() -> Self {
        CutPlane {
            axis: 2,
            offset: 0.0,
            rotation_x: 0.0,
            rotation_y: 0.0,
            range: None,
            isoband_frequency: None,
            clamp_colormap: true,
        }
    }
}

/// Options of [`Tracer`]. [`TraceOptions::default`] is ours (SPEC D17); [`TraceOptions::demo`] has
/// the demo's constants.
#[derive(Debug, Clone, PartialEq)]
pub struct TraceOptions {
    /// The level `ℓ` traced (`φ = ℓ`).
    pub level: f64,
    /// Step factor (0.9; the demo's 1.0).
    pub safety: f64,
    /// Lower bound of `‖∇φ‖` in the step (0.1; the demo has none, i.e. 0).
    pub g_min: f64,
    /// Maximum field evaluations per primary ray (128; the demo's 40).
    pub max_steps: u32,
    /// Hit threshold on `|φ − ℓ|` (1e-3 of the bounding-box diagonal; the demo's 1e-2 units).
    pub eps_hit: Epsilon,
    /// Maximum evaluations per shadow ray (128; the demo's 80).
    pub shadow_max_steps: u32,
    /// Hit threshold of shadow rays (1e-3 of the diagonal; the demo's 1e-4 units).
    pub shadow_eps: Epsilon,
    /// Bounding-sphere radius factor on the bounding-box diagonal (0.8, as the demo).
    pub bounding_factor: f64,
    /// Start rays at the entry of the bounding sphere (the demo starts at the camera).
    pub start_at_sphere: bool,
    /// Surface colour (the demo's contour purple `(68, 74, 145)/255`).
    pub base_color: [f64; 3],
    /// Background colour (white).
    pub background: [f64; 3],
    /// The environment of the IBL terms.
    pub environment: Environment,
    /// The diffuse point light, as an offset from the scene centre in units of `s` (demo `(0, 10, 2)`).
    pub light: [f64; 3],
    /// Direction of the specular light (demo `(-0.5, 1, 0)`).
    pub specular_dir: [f64; 3],
    /// The shadow-receiving ground (`None`: none).
    pub ground: Option<Ground>,
    /// The cut plane (`None`: none).
    pub cut_plane: Option<CutPlane>,
    /// Render threads (`None`: all cores).
    pub threads: Option<usize>,
    /// Where the traversal runs (`Cpu`; `Cuda`/`Auto` through [`Tracer::render_torus`]).
    pub device: crate::Device,
}

impl Default for TraceOptions {
    fn default() -> Self {
        TraceOptions {
            level: 0.0,
            safety: 0.9,
            g_min: 0.1,
            max_steps: 128,
            eps_hit: Epsilon::Relative(1e-3),
            shadow_max_steps: 128,
            shadow_eps: Epsilon::Relative(1e-3),
            bounding_factor: 0.8,
            start_at_sphere: true,
            base_color: [68.0 / 255.0, 74.0 / 255.0, 145.0 / 255.0],
            background: [1.0, 1.0, 1.0],
            environment: Environment::Forest,
            light: [0.0, 10.0, 2.0],
            specular_dir: [-0.5, 1.0, 0.0],
            ground: Some(Ground::default()),
            cut_plane: None,
            threads: None,
            device: crate::Device::Cpu,
        }
    }
}

impl TraceOptions {
    /// The demo's constants: safety 1, no gradient floor, 40 primary steps with `ε = 1e-2 s`,
    /// 80 shadow steps with `ε = 1e-4 s` (`s = diag/(2√3)`), rays started at the camera.
    pub fn demo() -> TraceOptions {
        let two_sqrt3 = 2.0 * 3.0f64.sqrt();
        TraceOptions {
            safety: 1.0,
            g_min: 0.0,
            max_steps: 40,
            eps_hit: Epsilon::Relative(1e-2 / two_sqrt3),
            shadow_max_steps: 80,
            shadow_eps: Epsilon::Relative(1e-4 / two_sqrt3),
            start_at_sphere: false,
            ..Default::default()
        }
    }

    fn validate(&self) -> Result<(), PatError> {
        if !self.level.is_finite() {
            return Err(bad(format!("level = {}", self.level)));
        }
        if !(self.safety.is_finite() && self.safety > 0.0 && self.safety <= 1.0) {
            return Err(bad(format!("safety = {} (in (0, 1])", self.safety)));
        }
        if !(self.g_min.is_finite() && self.g_min >= 0.0) {
            return Err(bad(format!("g_min = {} (finite, >= 0)", self.g_min)));
        }
        if self.max_steps == 0 || self.shadow_max_steps == 0 {
            return Err(bad("max_steps and shadow_max_steps must be > 0"));
        }
        self.eps_hit.validate("eps_hit")?;
        self.shadow_eps.validate("shadow_eps")?;
        if !(self.bounding_factor.is_finite() && self.bounding_factor > 0.0) {
            return Err(bad(format!("bounding_factor = {}", self.bounding_factor)));
        }
        for (name, v) in [
            ("base_color", self.base_color),
            ("background", self.background),
        ] {
            if !(finite3(v) && v.iter().all(|x| *x >= 0.0)) {
                return Err(bad(format!("{name} = {v:?} (finite, >= 0)")));
            }
        }
        if !finite3(self.light) || !finite3(self.specular_dir) || length(self.specular_dir) == 0.0 {
            return Err(bad(
                "light must be finite and specular_dir finite and nonzero",
            ));
        }
        if let Some(g) = &self.ground
            && (!finite3(g.up)
                || length(g.up) == 0.0
                || !g.gap.is_finite()
                || !g.shadow_value.is_finite())
        {
            return Err(bad(format!("ground {g:?}")));
        }
        if let Some(c) = &self.cut_plane {
            let ok = c.axis <= 2
                && c.offset.is_finite()
                && c.rotation_x.is_finite()
                && c.rotation_y.is_finite()
                && c.range.is_none_or(|r| r.is_finite() && r > 0.0)
                && c.isoband_frequency.is_none_or(|f| f.is_finite());
            if !ok {
                return Err(bad(format!("cut plane {c:?}")));
            }
        }
        if let Environment::Constant(c) = self.environment {
            CubeMap::constant(c)?;
        }
        Ok(())
    }
}

/// An RGB image (row 0 at the top) with an alpha channel (255 where the scene or a shadow is
/// drawn, 0 on the plain background).
#[derive(Debug, Clone, PartialEq)]
pub struct Image {
    /// Width.
    pub width: usize,
    /// Height.
    pub height: usize,
    /// `3 · width · height` bytes.
    pub rgb: Vec<u8>,
    /// `width · height` bytes.
    pub alpha: Vec<u8>,
}

impl Image {
    /// PNG bytes, RGB or (`with_alpha`) RGBA.
    pub fn to_png_bytes(&self, with_alpha: bool) -> Result<Vec<u8>, PatError> {
        let mut out = Vec::new();
        {
            let mut enc = png::Encoder::new(&mut out, self.width as u32, self.height as u32);
            enc.set_color(if with_alpha {
                png::ColorType::Rgba
            } else {
                png::ColorType::Rgb
            });
            enc.set_depth(png::BitDepth::Eight);
            let mut w = enc.write_header().map_err(|e| PatError::Format {
                reason: format!("png: {e}"),
            })?;
            let data: Vec<u8> = if with_alpha {
                self.rgb
                    .chunks_exact(3)
                    .zip(&self.alpha)
                    .flat_map(|(c, a)| [c[0], c[1], c[2], *a])
                    .collect()
            } else {
                self.rgb.clone()
            };
            w.write_image_data(&data).map_err(|e| PatError::Format {
                reason: format!("png: {e}"),
            })?;
        }
        Ok(out)
    }

    /// Write a PNG file.
    pub fn write_png(&self, path: &std::path::Path, with_alpha: bool) -> Result<(), PatError> {
        let bytes = self.to_png_bytes(with_alpha)?;
        std::fs::write(path, bytes).map_err(|e| PatError::Format {
            reason: format!("{}: {e}", path.display()),
        })
    }
}

/// What a render did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RenderReport {
    /// Primary rays (one per pixel).
    pub rays: u64,
    /// Rays that miss the bounding sphere or the range before the ground (no traversal).
    pub clipped: u64,
    /// Primary rays that hit the level set.
    pub hits: u64,
    /// Traversed primary rays that left the range without a hit.
    pub misses: u64,
    /// Rays that ran out of steps (`max_steps`).
    pub step_limit: u64,
    /// Rays stopped by a non-finite `φ` or gradient (treated as misses).
    pub non_finite: u64,
    /// Field evaluations of all primary rays.
    pub steps_total: u64,
    /// Mean evaluations per traversed primary ray.
    pub steps_mean: f64,
    /// Maximum evaluations of a primary ray.
    pub steps_max: u32,
    /// Shadow rays cast.
    pub shadow_rays: u64,
    /// Shadow rays that hit (shadowed ground pixels).
    pub shadow_hits: u64,
    /// Shadow rays that ran out of steps.
    pub shadow_step_limit: u64,
    /// Field evaluations of all shadow rays.
    pub shadow_steps_total: u64,
    /// Pixels coloured by the cut plane.
    pub cut_plane_pixels: u64,
    /// Wall time in seconds.
    pub seconds: f64,
    /// Threads used (CPU path).
    pub threads: usize,
    /// The device that traced the rays.
    pub device: crate::Device,
}

/// A ray segment to trace.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayJob {
    /// Origin.
    pub ro: [f64; 3],
    /// Unit direction.
    pub rd: [f64; 3],
    /// Start parameter.
    pub t0: f64,
    /// End parameter (exclusive).
    pub t1: f64,
}

/// Status codes of [`RayResult`] (the same values in the CUDA kernels).
pub mod status {
    /// Left the range without a hit.
    pub const MISS: u32 = 0;
    /// `|φ − ℓ| < ε`.
    pub const HIT: u32 = 1;
    /// Out of steps.
    pub const STEP_LIMIT: u32 = 2;
    /// `φ` or `∇φ` not finite.
    pub const NON_FINITE: u32 = 3;
    /// The field could not be evaluated (degenerate query).
    pub const EVAL_ERROR: u32 = 4;
}

/// Result of tracing one ray.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayResult {
    /// One of [`status`].
    pub status: u32,
    /// Field evaluations.
    pub steps: u32,
    /// Parameter of the hit.
    pub t: f64,
    /// `∇φ` at the hit (not normalized).
    pub grad: [f64; 3],
}

/// Traversal parameters, resolved against the scene.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceParams {
    /// Level.
    pub level: f64,
    /// Step factor.
    pub safety: f64,
    /// Gradient floor.
    pub g_min: f64,
    /// Hit threshold (absolute).
    pub eps: f64,
    /// Step budget.
    pub max_steps: u32,
}

/// The traversal of one ray, `[formula: pat-render-step]`. The GPU kernel repeats this
/// arithmetic operation for operation.
pub fn trace_ray(field: &dyn ImplicitField, job: &RayJob, p: &TraceParams) -> RayResult {
    let mut t = job.t0;
    let mut steps = 0u32;
    let mut res = RayResult {
        status: status::MISS,
        steps: 0,
        t: 0.0,
        grad: [0.0; 3],
    };
    while steps < p.max_steps && t < job.t1 {
        let x = [
            job.ro[0] + t * job.rd[0],
            job.ro[1] + t * job.rd[1],
            job.ro[2] + t * job.rd[2],
        ];
        steps += 1;
        let (phi, g) = match field.distance_gradient(x) {
            Ok(v) => v,
            Err(_) => {
                res.status = status::EVAL_ERROR;
                res.steps = steps;
                return res;
            }
        };
        let f = phi - p.level;
        let gl = ((g[0] * g[0] + g[1] * g[1]) + g[2] * g[2]).sqrt();
        if !(f.is_finite() && gl.is_finite()) {
            res.status = status::NON_FINITE;
            res.steps = steps;
            return res;
        }
        let af = f.abs();
        if af < p.eps {
            res.status = status::HIT;
            res.steps = steps;
            res.t = t;
            res.grad = g;
            return res;
        }
        let denom = if gl > p.g_min { gl } else { p.g_min };
        t += af / denom * p.safety;
    }
    res.status = if t >= job.t1 {
        status::MISS
    } else {
        status::STEP_LIMIT
    };
    res.steps = steps;
    res
}

/// Where rays are traced (CPU here, CUDA in `gpu.rs`).
pub trait Backend {
    /// Trace every job (results in job order).
    fn trace(&self, jobs: &[RayJob], p: &TraceParams) -> Result<Vec<RayResult>, PatError>;
    /// `φ` at each point.
    fn phi(&self, pts: &[[f64; 3]]) -> Result<Vec<f64>, PatError>;
    /// Which device this is.
    fn device(&self) -> crate::Device;
    /// Threads of the host pool.
    fn threads(&self) -> usize;
}

struct CpuBackend<'a> {
    field: &'a dyn ImplicitField,
    pool: rayon::ThreadPool,
}

impl Backend for CpuBackend<'_> {
    fn trace(&self, jobs: &[RayJob], p: &TraceParams) -> Result<Vec<RayResult>, PatError> {
        Ok(self.pool.install(|| {
            jobs.par_iter()
                .map(|j| trace_ray(self.field, j, p))
                .collect()
        }))
    }
    fn phi(&self, pts: &[[f64; 3]]) -> Result<Vec<f64>, PatError> {
        let mut out = vec![0.0; pts.len()];
        self.pool.install(|| self.field.distances(pts, &mut out))?;
        Ok(out)
    }
    fn device(&self) -> crate::Device {
        crate::Device::Cpu
    }
    fn threads(&self) -> usize {
        self.pool.current_num_threads()
    }
}

/// The sphere tracer.
#[derive(Debug, Clone)]
pub struct Tracer {
    options: TraceOptions,
    env: CubeMap,
}

/// Per-pixel setup.
struct Prep {
    ro: V3,
    rd: V3,
    job: RayJob,
    in_sphere: bool,
    ground: Option<f64>,
    plane: Option<f64>,
}

impl Tracer {
    /// A tracer with these options (validated).
    pub fn new(options: TraceOptions) -> Result<Tracer, PatError> {
        options.validate()?;
        let env = match options.environment {
            Environment::Forest => CubeMap::forest()?,
            Environment::Constant(c) => CubeMap::constant(c)?,
        };
        Ok(Tracer { options, env })
    }

    /// The options.
    pub fn options(&self) -> &TraceOptions {
        &self.options
    }

    /// Render `field` on the CPU (the pool of `options.threads`). `options.device` must be `Cpu`
    /// or `Auto` here; the CUDA path is [`Tracer::render_torus`].
    pub fn render(
        &self,
        field: &dyn ImplicitField,
        camera: &Camera,
    ) -> Result<(Image, RenderReport), PatError> {
        if let crate::Device::Cuda { .. } = self.options.device {
            return Err(PatError::Gpu {
                reason: "Tracer::render runs on the CPU; use render_torus for a CUDA field".into(),
            });
        }
        let pool = crate::pool::make_pool(self.options.threads)?;
        let backend = CpuBackend { field, pool };
        self.render_with(field, camera, &backend)
    }

    /// Render a [`crate::TorusField`] on `options.device`: `Cpu` is [`Tracer::render`]; `Cuda`
    /// needs the `cuda` feature and a field built with a CUDA evaluation device (kNN or all-points
    /// neighborhood); `Auto` uses the GPU when the field is on one.
    pub fn render_torus(
        &self,
        field: &crate::TorusField,
        camera: &Camera,
    ) -> Result<(Image, RenderReport), PatError> {
        match self.options.device {
            crate::Device::Cpu => self.render(field, camera),
            #[cfg(feature = "cuda")]
            crate::Device::Cuda { .. } | crate::Device::Auto => {
                if matches!(field.device(), crate::Device::Cpu) {
                    if matches!(self.options.device, crate::Device::Auto) {
                        return self.render(field, camera);
                    }
                    return Err(PatError::Gpu { reason: "the field was built for the CPU (EvalOptions::device); build it with a CUDA device to trace on the GPU".into() });
                }
                let backend = gpu::GpuBackend::new(field)?;
                self.render_with(field, camera, &backend)
            }
            #[cfg(not(feature = "cuda"))]
            crate::Device::Auto => self.render(field, camera),
            #[cfg(not(feature = "cuda"))]
            crate::Device::Cuda { .. } => Err(PatError::Gpu {
                reason: "tracing on the GPU needs the `cuda` feature".into(),
            }),
        }
    }

    /// The three-phase render over any backend.
    pub fn render_with(
        &self,
        field: &dyn ImplicitField,
        camera: &Camera,
        backend: &dyn Backend,
    ) -> Result<(Image, RenderReport), PatError> {
        let t_start = std::time::Instant::now();
        camera.validate()?;
        let o = &self.options;
        let bb = field.bbox();
        if !(finite3(bb[0]) && finite3(bb[1])) || (0..3).any(|a| bb[0][a] > bb[1][a]) {
            return Err(bad(format!(
                "field bounding box {bb:?} is not a finite box"
            )));
        }
        let center = [
            0.5 * (bb[0][0] + bb[1][0]),
            0.5 * (bb[0][1] + bb[1][1]),
            0.5 * (bb[0][2] + bb[1][2]),
        ];
        let diag = length(sub(bb[1], bb[0]));
        if diag <= 0.0 {
            return Err(bad("field bounding box has zero extent"));
        }
        let s = diag / (2.0 * 3.0f64.sqrt());
        let sphere_r = o.bounding_factor * diag + o.level.max(0.0);
        let params = TraceParams {
            level: o.level,
            safety: o.safety,
            g_min: o.g_min,
            eps: o.eps_hit.resolve(diag),
            max_steps: o.max_steps,
        };
        let shadow_params = TraceParams {
            eps: o.shadow_eps.resolve(diag),
            max_steps: o.shadow_max_steps,
            ..params
        };
        let lights = Lights {
            position: add(center, scale(o.light, s)),
            specular_dir: o.specular_dir,
        };
        let (w, h) = (camera.width, camera.height);
        let npix = w * h;
        let (cu, cv, cw) = camera.basis();
        let focal = 1.0 / (0.5 * camera.fov_y).tan();
        let ground = o.ground.map(|g| {
            let n = normalize(g.up);
            let mut lowest = 0.0;
            for a in 0..3 {
                lowest += (n[a] * bb[0][a]).min(n[a] * bb[1][a]);
            }
            (n, lowest - g.gap * s, g)
        });
        let plane = o.cut_plane.map(|c| CutGeom::new(&c, center, bb, s));

        // The sphere clip of a ray: parameters of entry and exit, or None.
        let sphere = |ro: V3, rd: V3| -> Option<(f64, f64)> {
            let oc = sub(ro, center);
            let b = dot(oc, rd);
            let c = dot(oc, oc) - sphere_r * sphere_r;
            let disc = b * b - c;
            if disc < 0.0 {
                return None;
            }
            let q = disc.sqrt();
            Some((-b - q, -b + q))
        };

        // Phase 0: rays and clip ranges.
        let preps: Vec<Prep> = (0..npix)
            .into_par_iter()
            .map(|i| {
                let (px, py) = ((i % w) as f64, (i / w) as f64);
                let ndc_x = (2.0 * (px + 0.5) - w as f64) / h as f64;
                let ndc_y = (2.0 * ((h as f64 - py) - 0.5) - h as f64) / h as f64;
                let rd = normalize(add(
                    add(scale(cu, ndc_x), scale(cv, ndc_y)),
                    scale(cw, focal),
                ));
                let ro = camera.eye;
                let mut tmax = f64::INFINITY;
                let mut tg = None;
                if let Some((n, d, _)) = &ground {
                    let denom = dot(*n, rd);
                    let t = (d - dot(*n, ro)) / denom;
                    if denom < 0.0 && t >= 0.0 && t.is_finite() {
                        tg = Some(t);
                        tmax = t;
                    }
                }
                let mut tp = None;
                if let Some(pg) = &plane
                    && let Some(t) = pg.intersect(ro, rd)
                    && t < tmax
                {
                    tp = Some(t);
                    tmax = t;
                }
                let (mut t0, mut t1, mut in_sphere) = (0.0, 0.0, false);
                if let Some((ts0, ts1)) = sphere(ro, rd) {
                    t0 = if o.start_at_sphere { ts0.max(0.0) } else { 0.0 };
                    t1 = ts1.min(tmax);
                    in_sphere = ts1 > 0.0 && t0 < t1;
                }
                Prep {
                    ro,
                    rd,
                    job: RayJob { ro, rd, t0, t1 },
                    in_sphere,
                    ground: tg,
                    plane: tp,
                }
            })
            .collect();

        // Phase A: primary rays.
        let jobs: Vec<RayJob> = preps.iter().map(|p| p.job).collect();
        let prim = backend.trace(&jobs, &params)?;
        drop(jobs);
        let mut rep = RenderReport {
            rays: npix as u64,
            threads: backend.threads(),
            device: backend.device(),
            ..Default::default()
        };
        for (p, r) in preps.iter().zip(&prim) {
            if !p.in_sphere {
                rep.clipped += 1;
                continue;
            }
            rep.steps_total += r.steps as u64;
            rep.steps_max = rep.steps_max.max(r.steps);
            match r.status {
                status::HIT => rep.hits += 1,
                status::MISS => rep.misses += 1,
                status::STEP_LIMIT => {
                    rep.misses += 1;
                    rep.step_limit += 1;
                }
                status::NON_FINITE => {
                    rep.misses += 1;
                    rep.non_finite += 1;
                }
                _ => {
                    return Err(PatError::DegenerateQuery {
                        index: rep.rays as usize,
                        reason: "the field could not be evaluated along a primary ray".into(),
                    });
                }
            }
        }
        let traversed = rep.hits + rep.misses;
        rep.steps_mean = if traversed > 0 {
            rep.steps_total as f64 / traversed as f64
        } else {
            0.0
        };

        // Phase B: shadow rays of the ground pixels that nothing else covers.
        let hit_of = |i: usize| preps[i].in_sphere && prim[i].status == status::HIT;
        let mut shadow_pix = Vec::new();
        let mut shadow_jobs = Vec::new();
        if let Some((_, _, _)) = &ground {
            for (i, pr) in preps.iter().enumerate() {
                if hit_of(i) || pr.plane.is_some() {
                    continue;
                }
                if let Some(tg) = pr.ground {
                    let gp = add(pr.ro, scale(pr.rd, tg));
                    let to_light = sub(lights.position, gp);
                    let dist = length(to_light);
                    let ld = scale(to_light, 1.0 / dist);
                    let ro = add(gp, scale(ld, 0.01 * s));
                    let (t0, t1, ok) = match sphere(ro, ld) {
                        Some((a, b)) => (a.max(0.0), b.min(dist), true),
                        None => (0.0, 0.0, false),
                    };
                    shadow_pix.push(i);
                    shadow_jobs.push(RayJob {
                        ro,
                        rd: ld,
                        t0,
                        t1: if ok { t1 } else { 0.0 },
                    });
                }
            }
        }
        let shadow_res = backend.trace(&shadow_jobs, &shadow_params)?;
        let mut shadowed = vec![false; npix];
        for (&i, r) in shadow_pix.iter().zip(&shadow_res) {
            rep.shadow_rays += 1;
            rep.shadow_steps_total += r.steps as u64;
            match r.status {
                status::HIT => {
                    rep.shadow_hits += 1;
                    shadowed[i] = true;
                }
                status::STEP_LIMIT => rep.shadow_step_limit += 1,
                status::EVAL_ERROR => {
                    return Err(PatError::DegenerateQuery {
                        index: i,
                        reason: "the field could not be evaluated along a shadow ray".into(),
                    });
                }
                _ => {}
            }
        }
        drop(shadow_jobs);

        // Phase D: φ at the visible cut-plane points.
        let mut plane_pix = Vec::new();
        let mut plane_pts = Vec::new();
        for (i, p) in preps.iter().enumerate() {
            if !hit_of(i)
                && let Some(t) = p.plane
            {
                plane_pix.push(i);
                plane_pts.push(add(p.ro, scale(p.rd, t)));
            }
        }
        let plane_phi = if plane_pts.is_empty() {
            Vec::new()
        } else {
            backend.phi(&plane_pts)?
        };
        let mut plane_val = vec![f64::NAN; npix];
        for (&i, &v) in plane_pix.iter().zip(&plane_phi) {
            plane_val[i] = v - o.level;
        }
        rep.cut_plane_pixels = plane_pix.len() as u64;

        // Phase C: colours.
        let (range, freq) = match &o.cut_plane {
            Some(c) => (
                c.range.unwrap_or(2.0 * s),
                c.isoband_frequency.unwrap_or(60.0 / s),
            ),
            None => (1.0, 1.0),
        };
        let to_u8 = |c: V3| [0, 1, 2].map(|a| (c[a].clamp(0.0, 1.0) * 255.0 + 0.5).floor() as u8);
        let shadow_value = o.ground.map(|g| g.shadow_value).unwrap_or(1.0);
        let px: Vec<([u8; 3], u8)> = (0..npix)
            .into_par_iter()
            .map(|i| {
                let p = &preps[i];
                if hit_of(i) {
                    let r = &prim[i];
                    let pos = add(p.ro, scale(p.rd, r.t));
                    let gl = length(r.grad);
                    let mut n = if gl > 0.0 {
                        scale(r.grad, 1.0 / gl)
                    } else {
                        scale(p.rd, -1.0)
                    };
                    if dot(n, p.rd) > 0.0 {
                        n = scale(n, -1.0);
                    }
                    (
                        to_u8(shade_surface(
                            &self.env,
                            &lights,
                            pos,
                            p.rd,
                            n,
                            o.base_color,
                        )),
                        255,
                    )
                } else if p.plane.is_some() {
                    let phi = plane_val[i];
                    // the demo's colour map over [-range, range] (non-finite φ is drawn as the lowest colour)
                    let mut v = if phi.is_finite() { phi } else { -range };
                    if o.cut_plane.is_some_and(|c| c.clamp_colormap) {
                        v = v.clamp(-range, range * (1.0 - 1e-9));
                    }
                    (to_u8(cutplane_color(v, -range, range, freq)), 255)
                } else if shadowed[i] {
                    (to_u8(scale(o.background, shadow_value)), 255)
                } else {
                    (to_u8(o.background), 0)
                }
            })
            .collect();
        let mut rgb = Vec::with_capacity(3 * npix);
        let mut alpha = Vec::with_capacity(npix);
        for (c, a) in px {
            rgb.extend_from_slice(&c);
            alpha.push(a);
        }
        rep.seconds = t_start.elapsed().as_secs_f64();
        Ok((
            Image {
                width: w,
                height: h,
                rgb,
                alpha,
            },
            rep,
        ))
    }
}

/// The geometry of the cut plane (`getCutPlaneBasis`, `intersectCutPlane` of the demo).
struct CutGeom {
    normal: V3,
    u: V3,
    v: V3,
    center: V3,
    half: [f64; 2],
    eps: f64,
}

fn rotate_about(v: V3, axis: V3, theta: f64) -> V3 {
    let (c, s) = (theta.cos(), theta.sin());
    let t = 1.0 - c;
    add(
        add(scale(v, c), scale(cross(axis, v), s)),
        scale(axis, dot(axis, v) * t),
    )
}

impl CutGeom {
    fn new(c: &CutPlane, bbox_center: V3, bb: [V3; 2], s: f64) -> CutGeom {
        let (n0, u0, v0) = match c.axis {
            0 => ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
            1 => ([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            _ => ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        };
        let mut normal = rotate_about(n0, u0, c.rotation_x);
        let mut u = rotate_about(u0, u0, c.rotation_x);
        let v = rotate_about(v0, u0, c.rotation_x);
        normal = rotate_about(normal, v, c.rotation_y);
        u = rotate_about(u, v, c.rotation_y);
        let (normal, u, v) = (normalize(normal), normalize(u), normalize(v));
        let diag = sub(bb[1], bb[0]);
        let half = match c.axis {
            0 => [diag[1], diag[2]],
            1 => [diag[0], diag[2]],
            _ => [diag[0], diag[1]],
        };
        CutGeom {
            normal,
            u,
            v,
            center: add(bbox_center, scale(n0, c.offset)),
            half,
            eps: 1e-2 * s,
        }
    }

    fn intersect(&self, ro: V3, rd: V3) -> Option<f64> {
        let eps = self.eps;
        let denom = dot(self.normal, rd);
        if denom.abs() < 1e-2 {
            return None;
        }
        let t = (dot(self.normal, self.center) - dot(self.normal, ro)) / denom;
        if !t.is_finite() || t <= eps {
            return None;
        }
        let local = sub(add(ro, scale(rd, t)), self.center);
        if dot(local, self.u).abs() > self.half[0] || dot(local, self.v).abs() > self.half[1] {
            return None;
        }
        Some(t)
    }
}
