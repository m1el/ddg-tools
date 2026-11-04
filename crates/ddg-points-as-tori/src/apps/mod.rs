//! Applications (SPEC.md §8): level-set sampling and morphology (WP6), normal orientation (WP9).
//!
//! [`sample_level_set`] is shared by morphology and orientation: marching cubes at a resolution,
//! area-uniform points on the triangles, then Newton projections onto `{φ = ℓ}` with residual
//! checks (SPEC.md §8.3, D14).

pub mod morph;
pub mod orient;

use crate::csg::ImplicitField;
use crate::mc::{McOptions, McReport, marching_cubes};
use crate::pool::make_pool;
use crate::PatError;
use rayon::prelude::*;
use std::time::Instant;

pub use morph::{Morph, MorphOptions, MorphReport, MorphStep, morphology};

/// Options of [`sample_level_set`].
#[derive(Debug, Clone, PartialEq)]
pub struct SampleOptions {
    /// Number of points (`None`: the field's [`ImplicitField::point_count_hint`], else 4096).
    pub count: Option<usize>,
    /// Marching cubes for the triangles (its `isovalue` is overridden by the level).
    pub mc: McOptions,
    /// Newton projections `x ← x − (φ(x) − ℓ) ∇φ / ‖∇φ‖²` (1 to 8).
    pub newton_iters: usize,
    /// A projected point is accepted if `|φ(x) − ℓ| ≤ residual_tol_rel · h` (`h` the grid spacing).
    pub residual_tol_rel: f64,
    /// A projected point must stay within this many grid spacings of its start.
    pub max_move_cells: f64,
    /// Seed of the (counter-based, thread-independent) random numbers.
    pub seed: u64,
    /// Rounds that replace rejected points (≥ 1).
    pub rounds: usize,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
}

impl Default for SampleOptions {
    fn default() -> Self {
        SampleOptions {
            count: None,
            mc: McOptions::default(),
            newton_iters: 3,
            residual_tol_rel: 1e-3,
            max_move_cells: 2.0,
            seed: 0,
            rounds: 4,
            threads: None,
        }
    }
}

/// What [`sample_level_set`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleReport {
    /// Marching cubes.
    pub mc: McReport,
    /// Points asked for.
    pub requested: usize,
    /// Points returned.
    pub produced: usize,
    /// Candidate points rejected (residual, gradient, or moved too far).
    pub rejected: usize,
    /// Largest `|φ − ℓ|` among the returned points, after projection.
    pub residual_max: f64,
    /// Mean `|φ − ℓ|` among the returned points.
    pub residual_mean: f64,
    /// Largest `|φ − ℓ|` of the candidates before projection (the MC interpolation error).
    pub residual_before_max: f64,
    /// Rounds used.
    pub rounds: usize,
    /// Seconds.
    pub seconds: f64,
}

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Uniform in `[0, 1)` from a counter.
fn unit(seed: u64, stream: u64, index: u64, k: u64) -> f64 {
    let h = splitmix(splitmix(splitmix(seed ^ stream.wrapping_mul(0xA24B_AED4_963E_E407)) ^ index) ^ k);
    (h >> 11) as f64 / (1u64 << 53) as f64
}

/// Points, unit normals and the report of [`sample_level_set`].
pub type LevelSetSamples = (Vec<[f64; 3]>, Vec<[f64; 3]>, SampleReport);

/// Points on `{φ = level}` and their unit normals `∇φ/‖∇φ‖`.
///
/// Marching cubes gives the triangles; points are drawn area-uniformly (stratified along the
/// cumulative area, so the result is deterministic and independent of the thread count), projected
/// by Newton steps, and kept only if the residual is below tolerance. Rejected points are replaced
/// in further rounds; `produced < requested` is possible and reported.
pub fn sample_level_set(
    field: &dyn ImplicitField,
    level: f64,
    opts: &SampleOptions,
) -> Result<LevelSetSamples, PatError> {
    let t0 = Instant::now();
    if !level.is_finite() {
        return Err(PatError::InvalidOption { what: format!("level = {level}") });
    }
    if !(1..=8).contains(&opts.newton_iters) || opts.rounds == 0 || opts.residual_tol_rel.is_nan() || opts.residual_tol_rel <= 0.0 || opts.max_move_cells.is_nan() || opts.max_move_cells <= 0.0 {
        return Err(PatError::InvalidOption {
            what: format!(
                "SampleOptions newton_iters = {}, rounds = {}, residual_tol_rel = {}, max_move_cells = {}",
                opts.newton_iters, opts.rounds, opts.residual_tol_rel, opts.max_move_cells
            ),
        });
    }
    let count = opts.count.or_else(|| field.point_count_hint()).unwrap_or(4096);
    if count == 0 {
        return Err(PatError::InvalidOption { what: "SampleOptions::count = 0".into() });
    }
    let mut mc_opts = opts.mc.clone();
    mc_opts.isovalue = level;
    mc_opts.normals = false;
    if mc_opts.threads.is_none() {
        mc_opts.threads = opts.threads;
    }
    let (mesh, mc) = marching_cubes(field, &mc_opts)?;
    if mesh.faces.is_empty() {
        return Err(PatError::InvalidOption { what: format!("the level set φ = {level} is empty on the grid {:?} of box origin {:?}", mc.cells, mc.origin) });
    }
    let h = mc.spacing.iter().cloned().fold(0.0, f64::max);
    let tol = opts.residual_tol_rel * h;
    let pool = make_pool(opts.threads)?;
    let areas = mesh.areas();
    let mut cum = Vec::with_capacity(areas.len());
    let mut acc = 0.0;
    for a in &areas {
        acc += a;
        cum.push(acc);
    }
    let total = acc;
    if total.is_nan() || total <= 0.0 {
        return Err(PatError::InvalidOption { what: "the level set has zero area".into() });
    }
    let mut pts: Vec<[f64; 3]> = Vec::with_capacity(count);
    let mut nrm: Vec<[f64; 3]> = Vec::with_capacity(count);
    let (mut rejected, mut rounds, mut res_before) = (0usize, 0usize, 0.0f64);
    let (mut res_max, mut res_sum) = (0.0f64, 0.0f64);
    for round in 0..opts.rounds {
        let need = count - pts.len();
        if need == 0 {
            break;
        }
        rounds = round + 1;
        // draw
        let cand: Vec<[f64; 3]> = pool.install(|| {
            (0..need)
                .into_par_iter()
                .map(|s| {
                    let u = (s as f64 + unit(opts.seed, round as u64, s as u64, 0)) / need as f64 * total;
                    let t = cum.partition_point(|c| *c <= u).min(cum.len() - 1);
                    let (r1, r2) = (unit(opts.seed, round as u64, s as u64, 1), unit(opts.seed, round as u64, s as u64, 2));
                    let sq = r1.sqrt();
                    let (w0, w1, w2) = (1.0 - sq, sq * (1.0 - r2), sq * r2);
                    let f = mesh.faces[t];
                    let (a, b, c) = (mesh.positions[f[0] as usize], mesh.positions[f[1] as usize], mesh.positions[f[2] as usize]);
                    [0, 1, 2].map(|i| w0 * a[i] + w1 * b[i] + w2 * c[i])
                })
                .collect()
        });
        let start = cand.clone();
        let mut x = cand;
        let mut d = vec![0.0; need];
        let mut g = vec![[0.0; 3]; need];
        let mut ok = vec![true; need];
        for it in 0..=opts.newton_iters {
            pool.install(|| field.distance_gradients(&x, &mut d, &mut g))?;
            if it == 0 {
                res_before = res_before.max(d.iter().map(|v| (v - level).abs()).fold(0.0, f64::max));
            }
            if it == opts.newton_iters {
                break;
            }
            for i in 0..need {
                let g2 = g[i][0] * g[i][0] + g[i][1] * g[i][1] + g[i][2] * g[i][2];
                if g2.is_nan() || g2 <= 1e-24 || !d[i].is_finite() {
                    ok[i] = false;
                    continue;
                }
                let s = (d[i] - level) / g2;
                for c in 0..3 {
                    x[i][c] -= s * g[i][c];
                }
            }
        }
        for i in 0..need {
            let g2 = g[i][0] * g[i][0] + g[i][1] * g[i][1] + g[i][2] * g[i][2];
            let moved = ((x[i][0] - start[i][0]).powi(2) + (x[i][1] - start[i][1]).powi(2) + (x[i][2] - start[i][2]).powi(2)).sqrt();
            let r = (d[i] - level).abs();
            if ok[i] && g2 > 1e-24 && r <= tol && moved <= opts.max_move_cells * h {
                let l = g2.sqrt();
                pts.push(x[i]);
                nrm.push([g[i][0] / l, g[i][1] / l, g[i][2] / l]);
                res_max = res_max.max(r);
                res_sum += r;
            } else {
                rejected += 1;
            }
        }
    }
    let produced = pts.len();
    let report = SampleReport {
        mc,
        requested: count,
        produced,
        rejected,
        residual_max: res_max,
        residual_mean: if produced > 0 { res_sum / produced as f64 } else { 0.0 },
        residual_before_max: res_before,
        rounds,
        seconds: t0.elapsed().as_secs_f64(),
    };
    Ok((pts, nrm, report))
}
