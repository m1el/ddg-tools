//! Exact `O(|E|²)` energy and differential (`SPEC.md` §3.1–3.2), parallel over edges with
//! results independent of the thread count (D25).
//!
//! Evaluation order (SPEC §3.2): per edge `I` the row `R_I = Σ_J (e_IJ + e_JI)` over `J` with
//! `(I, J) ∈ 𝒫` in increasing `J`, and its derivative with respect to the two endpoints of `I`
//! only; then `∂Ê/∂γ_v = Σ_{I ∋ v} G_{I,a(v)}` in increasing edge index, and `Ê = Σ_I Σ_J e_IJ`.
//! Rows are computed independently (any thread, any order); every reduction runs in a fixed
//! sequential order, hence the bitwise independence of the thread count.

use rayon::prelude::*;

use crate::error::RcError;
use crate::kernel::v3::{axpy, scale};
use crate::kernel::{Exponents, kernel, kernel_with_gradient};
use crate::network::CurveNetwork;

/// Builds an explicit rayon pool with `threads` workers (`None`: the number of available cores).
///
/// # Errors
/// [`RcError::InvalidOption`] for `Some(0)` or if the pool cannot be started.
pub fn make_pool(threads: Option<usize>) -> Result<rayon::ThreadPool, RcError> {
    let n = match threads {
        Some(0) => {
            return Err(RcError::InvalidOption {
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
        .map_err(|_| RcError::InvalidOption {
            name: "threads",
            value: n.to_string(),
            requirement: "the thread pool could not be started",
        })
}

/// Everything about the geometry that rows need.
struct Geom<'a> {
    network: &'a CurveNetwork,
    /// Endpoint positions per edge.
    ends: Vec<[[f64; 3]; 2]>,
    length: Vec<f64>,
    tangent: Vec<[f64; 3]>,
}

impl<'a> Geom<'a> {
    fn new(network: &'a CurveNetwork, positions: &[[f64; 3]]) -> Self {
        let g = network.edge_geometry(positions);
        let ends = network
            .edges()
            .iter()
            .map(|&[a, b]| [positions[a as usize], positions[b as usize]])
            .collect();
        Geom {
            network,
            ends,
            length: g.length,
            tangent: g.tangent,
        }
    }
}

/// `(k̂)_IJ = ¼ Σ_{i ∈ I} Σ_{j ∈ J} k(γ_i, γ_j, T_I)`.
///
/// [formula: rc-discrete-kernel]
#[inline]
fn discrete_kernel(pi: &[[f64; 3]; 2], pj: &[[f64; 3]; 2], t: [f64; 3], e: Exponents) -> f64 {
    let mut s = 0.0;
    for a in pi {
        for b in pj {
            s += kernel(*a, *b, t, e);
        }
    }
    0.25 * s
}

/// Row energy `Σ_J e_IJ` of edge `i` (J over `(I, J) ∈ 𝒫`, increasing).
///
/// [formula: rc-energy]
fn row_energy(g: &Geom, i: usize, e: Exponents) -> f64 {
    let mut sum = 0.0;
    for j in 0..g.length.len() {
        if g.network.edges_adjacent(i, j) {
            continue;
        }
        let khat = discrete_kernel(&g.ends[i], &g.ends[j], g.tangent[i], e);
        sum += khat * g.length[i] * g.length[j];
    }
    sum
}

/// Row energy `Σ_J e_IJ` and the derivative of `Σ_J (e_IJ + e_JI)` with respect to the two
/// endpoints of `I`.
///
/// [formula: rc-energy-differential]
fn row_with_gradient(g: &Geom, i: usize, e: Exponents) -> (f64, [[f64; 3]; 2]) {
    let (pi, ti, li) = (&g.ends[i], g.tangent[i], g.length[i]);
    let signs = [-1.0, 1.0];
    let mut energy = 0.0;
    let mut grad = [[0.0; 3]; 2];
    for j in 0..g.length.len() {
        if g.network.edges_adjacent(i, j) {
            continue;
        }
        let (pj, tj, lj) = (&g.ends[j], g.tangent[j], g.length[j]);
        // e_IJ = k̂_IJ ℓ_I ℓ_J with k evaluated at (γ_i, γ_j, T_I).
        let mut khat_ij = 0.0;
        let mut gp = [[0.0; 3]; 2]; // Σ_{j∈J} ∇_p k(γ_{i_a}, γ_j, T_I)
        let mut gt = [0.0; 3]; // Σ_i Σ_j ∇_T k
        for a in 0..2 {
            for b in 0..2 {
                let ke = kernel_with_gradient(pi[a], pj[b], ti, e);
                khat_ij += ke.value;
                gp[a] = axpy(gp[a], 1.0, ke.grad_p);
                gt = axpy(gt, 1.0, ke.grad_t);
            }
        }
        khat_ij *= 0.25;
        // e_JI = k̂_JI ℓ_J ℓ_I with k evaluated at (γ_j, γ_i, T_J); I is the second edge.
        let mut khat_ji = 0.0;
        let mut gq = [[0.0; 3]; 2]; // Σ_{j∈J} ∇_q k(γ_j, γ_{i_a}, T_J)
        for a in 0..2 {
            for b in 0..2 {
                let ke = kernel_with_gradient(pj[b], pi[a], tj, e);
                khat_ji += ke.value;
                gq[a] = axpy(gq[a], -1.0, ke.grad_p);
            }
        }
        khat_ji *= 0.25;

        energy += khat_ij * li * lj;
        let ll = 0.25 * li * lj;
        for a in 0..2 {
            // ∂e_IJ/∂γ_{i_a}
            let first = axpy(scale(ti, khat_ij), 0.25, gt);
            let mut v = scale(first, signs[a] * lj);
            v = axpy(v, ll, gp[a]);
            // ∂e_JI/∂γ_{i_a}
            v = axpy(v, signs[a] * lj * khat_ji, ti);
            v = axpy(v, ll, gq[a]);
            grad[a] = axpy(grad[a], 1.0, v);
        }
    }
    (energy, grad)
}

/// The energy `Ê` only, without the differential.
///
/// [formula: rc-energy]
pub fn energy_only(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
) -> f64 {
    let g = Geom::new(network, positions);
    let rows: Vec<f64> = pool.install(|| {
        (0..network.n_edges())
            .into_par_iter()
            .map(|i| row_energy(&g, i, exponents))
            .collect()
    });
    rows.iter().sum()
}

/// The energy `Ê` and its differential `∂Ê/∂γ_v` (one `[f64; 3]` per vertex), by the hand-derived
/// formulas of `SPEC.md` §3.2, with the exact `ν = 0` limit (D24).
///
/// The caller must have validated the positions
/// ([`CurveNetwork::validate_positions`]); coincident vertices give infinities.
///
/// [formula: rc-energy-differential]
pub fn energy_and_differential(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
) -> (f64, Vec<[f64; 3]>) {
    let g = Geom::new(network, positions);
    let rows: Vec<(f64, [[f64; 3]; 2])> = pool.install(|| {
        (0..network.n_edges())
            .into_par_iter()
            .map(|i| row_with_gradient(&g, i, exponents))
            .collect()
    });
    let mut energy = 0.0;
    let mut diff = vec![[0.0; 3]; network.n_vertices()];
    for (i, (e, gr)) in rows.iter().enumerate() {
        energy += e;
        let [a, b] = network.edges()[i];
        diff[a as usize] = axpy(diff[a as usize], 1.0, gr[0]);
        diff[b as usize] = axpy(diff[b as usize], 1.0, gr[1]);
    }
    (energy, diff)
}
