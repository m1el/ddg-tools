//! Orienting inconsistently oriented point clouds (Sec. 5.3, Fig. 27; SPEC.md §8.4, decision D16).
//!
//! Not in the authors' code; implemented from the paper's description and the SPEC:
//!
//! 1. fit the tori to the current normals (the network, then [`fit_tori_with`]),
//! 2. for every point `i` sum the field gradients at its neighbors with anisotropic Mahalanobis
//!    weights built from the fitted principal curvatures and directions,
//!    `ñ_i = Σ_{j ∈ N_k(i), j ≠ i} exp(−M_i(p_i − p_j)) ∇φ(p_j)`,
//!    `M_i(r) = (κmin² ⟨r, v̂min⟩² + κmax² ⟨r, v̂max⟩² + ε² ⟨r, n_i⟩²)^{1/2}`,
//! 3. flip `n_i` where `⟨n_i, ñ_i⟩ < 0` (all points at once: Jacobi), refit, repeat until no flip
//!    or `max_iterations`.
//!
//! The field `φ` is the blended torus SDF of the whole cloud with the *current* normals, so the
//! gradient at `p_j` is dominated by the torus of `p_j` and its neighbors; where the normals are
//! consistent it is `≈ n_j` and the vote of the neighborhood confirms `n_i`, where a point is
//! flipped against its neighbors the neighbors outvote it.
//!
//! The paper does not give `k`, `ε` or the flip schedule; ours are `k = 64`, `ε = 1/σ_i`
//! ([`Epsilon::InverseSigma`], `σ_i` the median neighbor distance of the network), Jacobi flips,
//! and the `ε` choice is measured in `NOTES-WP9a.md` (SPEC.md §16.3).

use crate::fit::{FitDiagnostics, fit_torus};
use crate::pool::make_pool;
use crate::spatial::{KdTree, KdTreeOptions, fitting_neighborhoods};
use crate::{BuildOptions, Model, NormalPolicy, PatError, PointCloud, PointsAsTori, TorusField};
use crate::field::EvalOut;
use rayon::prelude::*;
use std::time::Instant;

/// The weight `ε` of the normal direction in the metric `M_i`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Epsilon {
    /// `ε_i = scale / σ_i` with `σ_i` the median neighbor distance of the point's network
    /// neighborhood. `scale = 1` is the default of SPEC.md §8.4.
    InverseSigma {
        /// The factor in front of `1/σ_i`.
        scale: f64,
    },
    /// The same `ε` for all points (in the cloud's units).
    Constant(f64),
}

impl Default for Epsilon {
    fn default() -> Self {
        Epsilon::InverseSigma { scale: 1.0 }
    }
}

/// Options of [`orient_normals`].
#[derive(Debug, Clone, PartialEq)]
pub struct OrientOptions {
    /// Neighbors `k` of the vote (self excluded); default 64.
    pub k: usize,
    /// The normal-direction weight `ε`.
    pub epsilon: Epsilon,
    /// Maximum iterations (each: build, vote, flip); default 10 as in Fig. 27.
    pub max_iterations: usize,
    /// How the input normals are validated; default `RequireUnit`.
    pub normals: NormalPolicy,
    /// Network, fit and evaluation options of every refit. Its thread count is overridden by
    /// [`OrientOptions::threads`] when that is set.
    pub build: BuildOptions,
    /// Threads of every stage (`None`: all cores, or `build`'s own settings if given).
    pub threads: Option<usize>,
}

impl Default for OrientOptions {
    fn default() -> Self {
        OrientOptions {
            k: 64,
            epsilon: Epsilon::default(),
            max_iterations: 10,
            normals: NormalPolicy::default(),
            build: BuildOptions::default(),
            threads: None,
        }
    }
}

/// What [`orient_normals`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct OrientReport {
    /// Points.
    pub n: usize,
    /// Flips made in each iteration (the last entry is 0 if the iteration converged).
    pub flips_per_iteration: Vec<usize>,
    /// True if an iteration made no flips.
    pub converged: bool,
    /// Points whose final normal differs from the input normal.
    pub flipped_total: usize,
    /// Final count of `⟨n_i, ∇φ(p_i)⟩ < 0` for the returned normals and their own field.
    pub negative_final: usize,
    /// Per iteration, the count of `⟨n_i, ∇φ(p_i)⟩ < 0` before the flips of that iteration.
    pub negative_per_iteration: Vec<usize>,
    /// Points with `ñ_i = 0` or non-finite, which were not flipped (summed over iterations).
    pub undecided: usize,
    /// Wall time.
    pub seconds: f64,
}

/// `M_i(r)` of [`orient_normals`]: `(κmin² ⟨r, v̂min⟩² + κmax² ⟨r, v̂max⟩² + ε² ⟨r, n⟩²)^{1/2}`. A
/// norm when `κmin, κmax, ε ≠ 0` and `v̂min, v̂max, n` are independent (Lean `orientMetric_norm`).
///
/// [formula: pat-orientation-metric]
pub fn orientation_metric(r: [f64; 3], d: &FitDiagnostics, normal: [f64; 3], epsilon: f64) -> f64 {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let (a, b, c) = (d.kappa_min * dot(r, d.v_min), d.kappa_max * dot(r, d.v_max), epsilon * dot(r, normal));
    (a * a + b * b + c * c).sqrt()
}

/// The vote `ñ_i` of every point for the field `field` that was built from `cloud` (one
/// iteration's step 1, exposed for tests and diagnostics). Returns `(ñ, ∇φ(p_i))`.
///
/// `neighbors` are the rows `[i, j₁, j₂, …]` of [`fitting_neighborhoods`] (self first, skipped).
pub fn orientation_votes(
    cloud: &PointCloud,
    field: &TorusField,
    neighbors: &crate::Neighborhoods,
    epsilon: Epsilon,
    opts: &BuildOptions,
    threads: Option<usize>,
) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), PatError> {
    let n = cloud.len();
    if neighbors.len() != n {
        return Err(PatError::SizeMismatch { what: format!("{} neighborhood rows for {n} points", neighbors.len()) });
    }
    match epsilon {
        Epsilon::InverseSigma { scale } if !(scale.is_finite() && scale >= 0.0) => {
            return Err(PatError::InvalidOption { what: format!("Epsilon::InverseSigma scale = {scale}") });
        }
        Epsilon::Constant(e) if !(e.is_finite() && e >= 0.0) => {
            return Err(PatError::InvalidOption { what: format!("Epsilon::Constant({e})") });
        }
        _ => {}
    }
    let pos = cloud.positions();
    let mut grad = vec![[0.0f64; 3]; n];
    let mut eo = EvalOut::default();
    eo.gradient = Some(&mut grad);
    // the field's own pool is used (its thread count was set by the caller)
    field.evaluate(pos, eo)?;
    let tori = field.tori();
    let pool = make_pool(threads.or(opts.fit.threads))?;
    let votes: Vec<[f64; 3]> = pool.install(|| {
        (0..n)
            .into_par_iter()
            .map(|i| {
                let a = tori.coefficients[i].map(|x| x as f64);
                let (_, diag) = fit_torus(pos[i], cloud.normals()[i], a, &opts.fit);
                let eps = match epsilon {
                    Epsilon::InverseSigma { scale } => {
                        let s = tori.sigma[i];
                        if s > 0.0 && s.is_finite() { scale / s } else { 0.0 }
                    }
                    Epsilon::Constant(e) => e,
                };
                let row = neighbors.row(i);
                // weights exp(-M) with the smallest M factored out (a positive factor: only the sign of the vote is used)
                let m: Vec<f64> = row.iter().skip(1).map(|&j| orientation_metric(sub(pos[i], pos[j as usize]), &diag, cloud.normals()[i], eps)).collect();
                let m0 = m.iter().cloned().fold(f64::INFINITY, f64::min);
                let mut acc = [0.0; 3];
                for (&j, &mj) in row.iter().skip(1).zip(&m) {
                    let w = (-(mj - m0)).exp();
                    let g = grad[j as usize];
                    for c in 0..3 {
                        acc[c] += w * g[c];
                    }
                }
                acc
            })
            .collect()
    });
    let own: Vec<[f64; 3]> = grad;
    Ok((votes, own))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Orient the normals of a point cloud by the iteration of Fig. 27 (SPEC.md §8.4).
///
/// `positions` and `normals` are validated like [`PointCloud::new`]; the normals may be oriented
/// arbitrarily (their *lines* should be roughly right: this routine only flips, it never
/// re-estimates a direction). Returns the oriented normals and the report.
pub fn orient_normals(
    positions: &[[f64; 3]],
    normals: &[[f64; 3]],
    model: &Model,
    opts: &OrientOptions,
) -> Result<(Vec<[f64; 3]>, OrientReport), PatError> {
    let t0 = Instant::now();
    if opts.k == 0 {
        return Err(PatError::InvalidOption { what: "OrientOptions::k = 0".into() });
    }
    let (cloud0, _) = PointCloud::new(positions.to_vec(), normals.to_vec(), opts.normals)?;
    let mut build = opts.build.clone();
    if opts.threads.is_some() {
        build = build.with_threads(opts.threads);
    }
    let threads = build.eval.threads;
    let n = cloud0.len();
    let mut cur: Vec<[f64; 3]> = cloud0.normals().to_vec();
    let start = cur.clone();
    let tree = make_pool(threads)?.install(|| KdTree::new(cloud0.positions(), KdTreeOptions::default()))?;
    let nb = make_pool(threads)?.install(|| fitting_neighborhoods(&tree, cloud0.positions(), opts.k + 1));
    let mut report = OrientReport { n, ..Default::default() };
    let mut negative_final = None;
    for _ in 0..opts.max_iterations {
        let (cloud, _) = PointCloud::new(positions.to_vec(), cur.clone(), NormalPolicy::Normalize)?;
        let (field, _) = PointsAsTori::build(&cloud, model, &build)?;
        let (votes, own) = orientation_votes(&cloud, &field, &nb, opts.epsilon, &build, threads)?;
        let neg = (0..n).filter(|&i| dot(cur[i], own[i]) < 0.0).count();
        report.negative_per_iteration.push(neg);
        let mut flips = 0;
        for i in 0..n {
            let s = dot(cur[i], votes[i]);
            if !s.is_finite() || votes[i] == [0.0; 3] {
                report.undecided += 1;
            } else if s < 0.0 {
                cur[i] = [-cur[i][0], -cur[i][1], -cur[i][2]];
                flips += 1;
            }
        }
        report.flips_per_iteration.push(flips);
        if flips == 0 {
            report.converged = true;
            negative_final = Some(neg);
            break;
        }
    }
    let neg = match negative_final {
        Some(v) => v,
        None => {
            // the normals changed in the last iteration: one more fit to count
            let (cloud, _) = PointCloud::new(positions.to_vec(), cur.clone(), NormalPolicy::Normalize)?;
            let (field, _) = PointsAsTori::build(&cloud, model, &build)?;
            let mut grad = vec![[0.0f64; 3]; n];
            let mut eo = EvalOut::default();
            eo.gradient = Some(&mut grad);
            field.evaluate(cloud.positions(), eo)?;
            (0..n).filter(|&i| dot(cur[i], grad[i]) < 0.0).count()
        }
    };
    report.negative_final = neg;
    report.flipped_total = (0..n).filter(|&i| cur[i] != start[i]).count();
    report.seconds = t0.elapsed().as_secs_f64();
    Ok((cur, report))
}

/// Unoriented normal estimates: the smallest principal axis of the covariance of the `k` nearest
/// points (self included), with a **random sign** from a seeded generator. This is the input of
/// the "unoriented estimate" experiment of `NOTES-WP9a.md`: the method is then given lines but no
/// orientation.
pub fn pca_normals_random_sign(positions: &[[f64; 3]], k: usize, seed: u64, threads: Option<usize>) -> Result<Vec<[f64; 3]>, PatError> {
    if k < 3 {
        return Err(PatError::InvalidOption { what: format!("pca_normals_random_sign: k = {k} < 3") });
    }
    if positions.len() < 3 {
        return Err(PatError::TooFewPoints { n: positions.len(), need: 3 });
    }
    if let Some(i) = positions.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
        return Err(PatError::NonFinite { index: i, what: "position" });
    }
    let pool = make_pool(threads)?;
    let tree = pool.install(|| KdTree::new(positions, KdTreeOptions::default()))?;
    let nb = pool.install(|| fitting_neighborhoods(&tree, positions, k));
    Ok(pool.install(|| {
        (0..positions.len())
            .into_par_iter()
            .map(|i| {
                let row = nb.row(i);
                let mut c = [0.0; 3];
                for &j in row {
                    for a in 0..3 {
                        c[a] += positions[j as usize][a] / row.len() as f64;
                    }
                }
                let mut cov = [[0.0f64; 3]; 3];
                for &j in row {
                    let d = sub(positions[j as usize], c);
                    for a in 0..3 {
                        for b in 0..3 {
                            cov[a][b] += d[a] * d[b];
                        }
                    }
                }
                let mut v = smallest_eigvec(cov);
                // sign from a counter-based hash of (seed, i)
                let mut z = seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                if z & 1 == 1 {
                    v = [-v[0], -v[1], -v[2]];
                }
                v
            })
            .collect()
    }))
}

/// Eigenvector of the smallest eigenvalue of a symmetric 3x3 matrix by shifted inverse-free power
/// iteration on `σI − C` (deterministic; plenty for a normal estimate).
fn smallest_eigvec(c: [[f64; 3]; 3]) -> [f64; 3] {
    let tr = c[0][0] + c[1][1] + c[2][2];
    let shift = tr.max(1e-300);
    let mut v = [0.577_350_269_189_626, 0.577_350_269_189_626 + 1e-3, 0.577_350_269_189_626 - 1e-3];
    for _ in 0..200 {
        let mut w = [0.0; 3];
        for a in 0..3 {
            w[a] = shift * v[a] - (c[a][0] * v[0] + c[a][1] * v[1] + c[a][2] * v[2]);
        }
        let l = dot(w, w).sqrt();
        if l == 0.0 {
            break;
        }
        v = [w[0] / l, w[1] / l, w[2] / l];
    }
    v
}
