//! The blob generator (NOTES-REFERENCE §11; `mesh_data.py` `generate_blob_mesh`): the implicit
//! function `F(x) = ‖x‖² − (R + α·2·n(scale·x))²` with `n` fractal Perlin noise
//! ([`Perlin`](super::perlin::Perlin), our own implementation, `base = seed % 256`), sampled on a
//! `64³` grid over `±1.2 R (1 + α)` and contoured by our marching cubes ([`crate::mc`], dense,
//! `F = 0`), then simplified to at most 2048 faces ([`super::simplify`], standing in for
//! `pyfqmr`), checked for watertightness, centered and scaled to `[-1, 1]³`.
//!
//! `F` is an implicit function, not a distance; [`BlobField`]'s gradient is exact in `x` and the
//! noise (central differences of the noise, step `1e-5`).

use super::perlin::Perlin;
use super::rng::Rng;
use super::sampling::center_and_scale;
use super::simplify::simplify_closed;
use crate::PatError;
use crate::csg::ImplicitField;
use crate::mc::{Band, McOptions, MeshCheck, Resolution, marching_cubes};
use crate::truth::check_mesh;

/// Parameters of one blob.
#[derive(Debug, Clone, PartialEq)]
pub struct BlobParams {
    /// Deformation amplitude `α`.
    pub alpha: f64,
    /// Noise frequency scale.
    pub scale: f64,
    /// Noise octaves (1..=4 in the reference).
    pub octaves: usize,
    /// Octave amplitude ratio.
    pub persistence: f64,
    /// Octave frequency ratio.
    pub lacunarity: f64,
    /// Base radius `R`.
    pub base_radius: f64,
    /// Noise field number (the reference uses `seed % 256`).
    pub noise_base: u64,
    /// Marching-cubes grid points per axis (64).
    pub grid: usize,
}

impl BlobParams {
    /// Draw parameters as the reference's dataset generator: `α ∈ (0.1, 0.7)`, `scale ∈ (0.5, 2)`,
    /// `octaves ∈ 1..=4`, persistence 0.5, lacunarity 2, `R = 1`, 64 grid points.
    pub fn draw(rng: &mut Rng, noise_base: u64) -> BlobParams {
        BlobParams {
            alpha: rng.uniform(0.1, 0.7),
            scale: rng.uniform(0.5, 2.0),
            octaves: 1 + rng.below(4),
            persistence: 0.5,
            lacunarity: 2.0,
            base_radius: 1.0,
            noise_base,
            grid: 64,
        }
    }

    fn extent(&self) -> f64 {
        self.base_radius * (1.0 + self.alpha) * 1.2
    }
}

/// The implicit function `F` of a blob.
#[derive(Debug, Clone)]
pub struct BlobField {
    params: BlobParams,
    noise: Perlin,
}

impl BlobField {
    /// The field of `params`.
    pub fn new(params: BlobParams) -> BlobField {
        let noise = Perlin::new(params.noise_base);
        BlobField { params, noise }
    }

    fn n(&self, x: [f64; 3]) -> f64 {
        let p = &self.params;
        self.noise.fbm3([x[0] * p.scale, x[1] * p.scale, x[2] * p.scale], p.octaves, p.persistence, p.lacunarity)
    }

    /// `F(x)`.
    pub fn value(&self, x: [f64; 3]) -> f64 {
        let r = self.params.base_radius + self.params.alpha * 2.0 * self.n(x);
        x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - r * r
    }
}

impl ImplicitField for BlobField {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        if let Some(a) = x.iter().position(|c| !c.is_finite()) {
            return Err(PatError::NonFinite { index: a, what: "query" });
        }
        Ok(self.value(x))
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let f = self.distance(x)?;
        let h = 1e-5;
        let r = self.params.base_radius + self.params.alpha * 2.0 * self.n(x);
        let mut g = [0.0; 3];
        for a in 0..3 {
            let (mut xp, mut xm) = (x, x);
            xp[a] += h;
            xm[a] -= h;
            let dn = (self.n(xp) - self.n(xm)) / (2.0 * h);
            g[a] = 2.0 * x[a] - 2.0 * r * 2.0 * self.params.alpha * dn;
        }
        Ok((f, g))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        let e = self.params.extent();
        [[-e; 3], [e; 3]]
    }
}

/// A watertight triangle mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshData {
    /// Vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Outward-oriented triangles.
    pub faces: Vec<[u32; 3]>,
}

/// Options of [`generate_blob_mesh`].
#[derive(Debug, Clone, PartialEq)]
pub struct BlobOptions {
    /// Faces after simplification (2048).
    pub max_faces: usize,
    /// Reject meshes with more than one connected component (the reference accepts them).
    pub require_single_component: bool,
    /// Threads of the marching cubes (`Some(1)` when blobs are generated in parallel).
    pub threads: Option<usize>,
}

impl Default for BlobOptions {
    fn default() -> Self {
        BlobOptions { max_faces: 2048, require_single_component: false, threads: Some(1) }
    }
}

/// Contour, simplify and validate one blob (no centering). Fails with `PatError::Format` naming
/// the reason when the mesh is not a closed oriented manifold with positive volume.
pub fn generate_blob_mesh(params: &BlobParams, opts: &BlobOptions) -> Result<MeshData, PatError> {
    if !(params.alpha.is_finite() && params.alpha >= 0.0 && params.scale.is_finite() && params.base_radius > 0.0 && params.grid >= 8) {
        return Err(PatError::InvalidOption { what: format!("blob parameters {params:?}") });
    }
    let field = BlobField::new(params.clone());
    let e = params.extent();
    let mc = McOptions {
        bbox: Some([[-e; 3], [e; 3]]),
        resolution: Resolution::Cells([params.grid - 1; 3]),
        band: Band::Dense,
        check_mesh: true,
        threads: opts.threads,
        ..McOptions::default()
    };
    let (mesh, rep) = marching_cubes(&field, &mc)?;
    let fail = |why: String| PatError::Format { reason: format!("blob {params:?}: {why}") };
    match &rep.mesh {
        MeshCheck::Manifold { boundary_edges: 0, components, .. } => {
            if opts.require_single_component && *components != 1 {
                return Err(fail(format!("{components} connected components")));
            }
        }
        MeshCheck::Manifold { boundary_edges, .. } => return Err(fail(format!("{boundary_edges} boundary edges"))),
        MeshCheck::Invalid(m) => return Err(fail(format!("not a manifold: {m}"))),
        MeshCheck::Skipped => return Err(fail("mesh check skipped".into())),
    }
    let (vertices, faces) = if mesh.faces.len() > opts.max_faces {
        simplify_closed(&mesh.positions, &mesh.faces, opts.max_faces)
    } else {
        (mesh.positions.clone(), mesh.faces.clone())
    };
    if faces.len() > opts.max_faces {
        return Err(fail(format!("simplification stopped at {} faces (> {})", faces.len(), opts.max_faces)));
    }
    let r = check_mesh(&vertices, &faces)?;
    if !r.is_closed() || r.degenerate_faces > 0 || !(r.signed_volume > 0.0) {
        return Err(fail(format!("simplified mesh is not a closed outward manifold: {r:?}")));
    }
    Ok(MeshData { vertices, faces })
}

/// Blob number `index` of the run `seed`: parameters from `Rng::for_item(seed, index)` (noise base
/// `index % 256`), retried with fresh draws (up to 20 times) when the mesh fails validation.
/// Returns the mesh centered and scaled to `[-1, 1]³` and the parameters used.
pub fn generate_blob(seed: u64, index: u64, opts: &BlobOptions) -> Result<(MeshData, BlobParams), PatError> {
    let mut last = None;
    for attempt in 0..20u64 {
        let mut rng = Rng::for_item(seed, index.wrapping_mul(64).wrapping_add(attempt));
        let params = BlobParams::draw(&mut rng, (index.wrapping_add(attempt)) % 256);
        match generate_blob_mesh(&params, opts) {
            Ok(mut m) => {
                m.vertices = center_and_scale(&m.vertices);
                return Ok((m, params));
            }
            Err(e) => last = Some(e),
        }
    }
    Err(last.unwrap_or(PatError::Format { reason: "no blob attempt was made".into() }))
}
