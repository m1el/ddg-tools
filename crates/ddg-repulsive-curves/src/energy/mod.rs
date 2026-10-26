//! The discrete tangent-point energy and its differential (`SPEC.md` §3).
//!
//! ```text
//! Ê(γ) = Σ_{(I,J) ∈ 𝒫} (k̂)_IJ ℓ_I ℓ_J,     (k̂)_IJ = ¼ Σ_{i ∈ I} Σ_{j ∈ J} k(γ_i, γ_j, T_I)
//! ```
//!
//! # Invariances and scaling
//!
//! `Ê` is invariant under translations, rotations, reflections and edge re-orientation and
//! homogeneous of degree `2 + α − β`: `Ê(cγ) = c^{2+α−β} Ê(γ)` for `c > 0`
//! (see [`Exponents::energy_scaling`]).
//!
//! [formula: rc-energy-scaling]
//!
//! Consequently the differential `g_v = ∂Ê/∂γ_v` satisfies
//! `Σ_v g_v = 0`, `Σ_v γ_v × g_v = 0` and the Euler identity `Σ_v γ_v · g_v = (2 + α − β) Ê`.
//!
//! [formula: rc-gradient-identities]

use std::time::Instant;

use crate::error::RcError;
use crate::kernel::Exponents;
use crate::network::CurveNetwork;

pub mod barnes_hut;
pub mod exact;

/// Parameters of the Barnes–Hut evaluation (`SPEC.md` §8.2; implemented by WP4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BhOptions {
    /// Spatial admissibility parameter (bound on `r_x / d`).
    pub theta_x: f64,
    /// Tangent admissibility parameter.
    pub theta_t: f64,
    /// Maximal number of edges in a leaf.
    pub leaf_size: usize,
}

impl Default for BhOptions {
    /// `theta_x 0.125, theta_t 0.25, leaf_size 8` (`SPEC.md` §10).
    fn default() -> Self {
        BhOptions {
            theta_x: 0.125,
            theta_t: 0.25,
            leaf_size: 8,
        }
    }
}

/// How the energy and its differential are evaluated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnergyEval {
    /// All pairs, `O(|E|²)`.
    Exact,
    /// Barnes–Hut with the given parameters ([`barnes_hut`]).
    BarnesHut(BhOptions),
}

/// Result of [`energy`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyReport {
    /// `Ê(γ)`.
    pub energy: f64,
    /// Values of the penalties (empty: penalties are added by work package WP7).
    pub penalties: Vec<f64>,
    /// `∂Ê/∂γ_v` for every vertex.
    pub differential: Vec<[f64; 3]>,
    /// Number of ordered edge pairs evaluated exactly.
    pub pairs_exact: u64,
    /// Number of lumped nodes used (0 for the exact evaluation).
    pub nodes_lumped: u64,
    /// Wall time of the evaluation.
    pub seconds: f64,
}

/// Energy and differential of the curve `positions` on `network`.
///
/// `threads` is the size of the explicit rayon pool used (`None`: all cores); the result does not
/// depend on it, bit for bit (D25). Positions are validated with
/// [`CurveNetwork::validate_positions`].
///
/// Penalties (`SPEC.md` §5.3) belong to work package WP7; the signature of §9.1 gets its
/// `penalties` argument there.
///
/// # Errors
/// The errors of [`CurveNetwork::validate_positions`]; [`RcError::InvalidOption`] for
/// `threads = Some(0)`; [`RcError::InvalidOption`] for bad [`BhOptions`].
pub fn energy(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    eval: &EnergyEval,
    threads: Option<usize>,
) -> Result<EnergyReport, RcError> {
    let t0 = Instant::now();
    network.validate_positions(positions)?;
    let pool = exact::make_pool(threads)?;
    if let EnergyEval::BarnesHut(opts) = eval {
        let (energy, differential, stats) =
            barnes_hut::energy_and_differential(&pool, network, positions, exponents, *opts)?;
        return Ok(EnergyReport {
            energy,
            penalties: Vec::new(),
            differential,
            pairs_exact: stats.pairs_exact,
            nodes_lumped: stats.nodes_lumped,
            seconds: t0.elapsed().as_secs_f64(),
        });
    }
    let (energy, differential) =
        exact::energy_and_differential(&pool, network, positions, exponents);
    let m = network.n_edges() as u64;
    let adjacent: u64 = (0..network.n_vertices() as u32)
        .map(|v| {
            let d = network.degree(v) as u64;
            d * (d - 1)
        })
        .sum();
    Ok(EnergyReport {
        energy,
        penalties: Vec::new(),
        differential,
        pairs_exact: m * (m - 1) - adjacent,
        nodes_lumped: 0,
        seconds: t0.elapsed().as_secs_f64(),
    })
}
