//! The exact `O(|F|²)` tangent-point energy and its exact differential (`SPEC.md` §3.1–3.2).
//!
//! Parallel over faces on an explicit rayon pool. Every face's task evaluates its row
//! (`S` as the source, `K(X_S, N_S; X_T)`) and its column (`S` as the target,
//! `K(X_T, N_T; X_S)`) in increasing `T`, sharing `r^{−β}` between the two; the per-face
//! results are combined in face order. The result is therefore bitwise independent of the
//! number of threads.

use rayon::prelude::*;

use crate::error::RsError;
use crate::geometry::v3::{axpy, sub};
use crate::geometry::{FaceGeometry, pull_back};
use crate::kernel::{Exponents, kernel_with_gradient_of_v, pair_terms};
use crate::surface::Surface;

/// An explicit rayon pool with `threads` threads (`None`: all cores). Equivalent to
/// `ddg_repulsive_curves::energy::exact::make_pool`, with this crate's error.
///
/// # Errors
/// [`RsError::InvalidOption`] for `Some(0)` or if the pool cannot be started.
pub fn make_pool(threads: Option<usize>) -> Result<rayon::ThreadPool, RsError> {
    let n = match threads {
        Some(0) => {
            return Err(RsError::InvalidOption {
                name: "threads",
                value: "0".into(),
                requirement: "at least 1 thread",
            });
        }
        Some(n) => n,
        None => std::thread::available_parallelism().map_or(1, |n| n.get()),
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .map_err(|_| RsError::InvalidOption {
            name: "threads",
            value: n.to_string(),
            requirement: "the thread pool could not be started",
        })
}

/// Per-face partials of `ℰ̂` (`∂/∂a_S`, `∂/∂X_S`, `∂/∂N_S`) and the energy.
#[derive(Debug, Clone, PartialEq)]
pub struct FacePartials {
    /// `ℰ̂`.
    pub energy: f64,
    /// `∂ℰ̂/∂a_S`.
    pub d_area: Vec<f64>,
    /// `∂ℰ̂/∂X_S`.
    pub d_bary: Vec<[f64; 3]>,
    /// `∂ℰ̂/∂N_S` (unconstrained partial in `R³`).
    pub d_normal: Vec<[f64; 3]>,
}

/// `ℰ̂(f) = Σ_S Σ_{T≠S} K(X_S, N_S; X_T) a_S a_T` (Eq. 10; ordered pairs, neighbours included,
/// no `2^p`, no ½). Value only. May be non-finite (coincident barycenters); callers check.
///
/// [formula: rs-energy]
pub fn energy_value(pool: &rayon::ThreadPool, geom: &FaceGeometry, e: &Exponents) -> f64 {
    let nf = geom.areas.len();
    let rows: Vec<f64> = pool.install(|| {
        (0..nf)
            .into_par_iter()
            .map(|s| {
                let xs = geom.barycenters[s];
                let ns = geom.normals[s];
                let mut acc = 0.0;
                for t in 0..nf {
                    if t != s {
                        let v = sub(xs, geom.barycenters[t]);
                        acc += geom.areas[t] * crate::kernel::kernel_of_v(e, v, ns);
                    }
                }
                geom.areas[s] * acc
            })
            .collect()
    });
    rows.iter().sum()
}

/// `ℰ̂` and its per-face partials (see [`FacePartials`]). May be non-finite; callers check.
///
/// [formula: rs-energy]
/// [formula: rs-energy-gradient]
pub fn face_partials(pool: &rayon::ThreadPool, geom: &FaceGeometry, e: &Exponents) -> FacePartials {
    let nf = geom.areas.len();
    let per_face: Vec<(f64, f64, [f64; 3], [f64; 3])> = pool.install(|| {
        (0..nf)
            .into_par_iter()
            .map(|s| {
                let xs = geom.barycenters[s];
                let ns = geom.normals[s];
                let a_s = geom.areas[s];
                let mut e_row = 0.0;
                let mut d_a = 0.0;
                let mut d_x = [0.0; 3];
                let mut d_n = [0.0; 3];
                for t in 0..nf {
                    if t == s {
                        continue;
                    }
                    let a_t = geom.areas[t];
                    let v = sub(xs, geom.barycenters[t]);
                    let (k_st, k_ts, dv_st, dv_ts, dn_st) = pair_terms(e, v, ns, geom.normals[t]);
                    e_row += a_t * k_st;
                    d_a += a_t * (k_st + k_ts);
                    d_x = axpy(d_x, a_t, sub(dv_st, dv_ts));
                    d_n = axpy(d_n, a_t, dn_st);
                }
                (
                    a_s * e_row,
                    d_a,
                    crate::geometry::v3::scale(d_x, a_s),
                    crate::geometry::v3::scale(d_n, a_s),
                )
            })
            .collect()
    });
    let mut out = FacePartials {
        energy: 0.0,
        d_area: Vec::with_capacity(nf),
        d_bary: Vec::with_capacity(nf),
        d_normal: Vec::with_capacity(nf),
    };
    for (en, da, dx, dn) in per_face {
        out.energy += en;
        out.d_area.push(da);
        out.d_bary.push(dx);
        out.d_normal.push(dn);
    }
    out
}

/// The exact energy and its differential `g_v = ∂ℰ̂/∂f_v`, by the chain rule of `SPEC.md`
/// §3.2 through `(a_S, X_S, N_S)` (`∂X_S/∂f_i = I/3`, [`crate::geometry::face_pull_back`]).
/// Positions must be valid ([`Surface::validate_positions`]).
///
/// [formula: rs-energy]
/// [formula: rs-energy-gradient]
///
/// # Errors
/// [`RsError::NonFiniteEnergy`] naming a pair of faces when the energy or a differential
/// entry is not finite.
pub fn energy_and_differential(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    geom: &FaceGeometry,
    e: &Exponents,
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    let p = face_partials(pool, geom, e);
    let g = pull_back(surface, positions, &p.d_area, &p.d_bary, &p.d_normal);
    if !p.energy.is_finite() || !g.iter().flatten().all(|x| x.is_finite()) {
        return Err(RsError::NonFiniteEnergy {
            faces: culprit_pair(geom, e),
        });
    }
    Ok((p.energy, g))
}

/// The first ordered pair `(S, T)` (in row-major order) whose term `K a_S a_T` or its
/// gradient is not finite; if there is none, the pair with the largest term.
pub(crate) fn culprit_pair(geom: &FaceGeometry, e: &Exponents) -> [u32; 2] {
    let nf = geom.areas.len();
    let mut best = ([0u32, 0u32], f64::NEG_INFINITY);
    for s in 0..nf {
        for t in 0..nf {
            if s == t {
                continue;
            }
            let v = sub(geom.barycenters[s], geom.barycenters[t]);
            let k = kernel_with_gradient_of_v(e, v, geom.normals[s]);
            let term = k.value * geom.areas[s] * geom.areas[t];
            let finite = term.is_finite() && k.d_v.iter().chain(&k.d_n).all(|x| x.is_finite());
            if !finite {
                return [s as u32, t as u32];
            }
            if term > best.1 {
                best = ([s as u32, t as u32], term);
            }
        }
    }
    best.0
}

/// Number of ordered pairs `(S, T)`, `S ≠ T`, evaluated by the exact energy.
pub fn exact_pair_count(n_faces: usize) -> u64 {
    let n = n_faces as u64;
    n * n.saturating_sub(1)
}
