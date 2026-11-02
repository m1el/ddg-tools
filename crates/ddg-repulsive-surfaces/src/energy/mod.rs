//! The discrete tangent-point energy and its differential (`SPEC.md` §3).
//!
//! ```text
//! ℰ̂(f) = Σ_{S ∈ F} Σ_{T ∈ F, T ≠ S} K(X_S, N_S; X_T) a_S a_T,     K(X, N; Y) = |N·(X − Y)|^α |X − Y|^{−β}
//! ```
//!
//! over ordered pairs (both orientations), `S = T` excluded, neighbours not excluded; no `2^p`,
//! no ½ (paper Eq. 10, 16).
//!
//! # Invariances and scaling (`SPEC.md` §3.3)
//!
//! `ℰ̂` is invariant under rigid motions and under flipping any face's orientation (the kernel
//! uses `|N·v|`), and homogeneous: `ℰ̂(cf) = c^{4+α−β} ℰ̂(f)` for `c > 0`
//! ([`Exponents::energy_scaling`]).
//!
//! [formula: rs-energy-scaling]
//!
//! Consequently the differential `g_v = ∂ℰ̂/∂f_v` satisfies `Σ_v g_v = 0`,
//! `Σ_v f_v × g_v = 0` and the Euler identity `Σ_v f_v · g_v = (4 + α − β) ℰ̂`.
//!
//! [formula: rs-gradient-identities]
//!
//! # Work packages
//!
//! WP1 wrote the exact evaluation ([`exact`]) and the entry point [`energy`]. The Barnes–Hut
//! evaluation lives in [`barnes_hut`] (WP3) behind the hook
//! `barnes_hut::energy_and_differential`; penalties live in [`crate::penalties`] (WP7)
//! behind the hook `penalties::evaluate`. Both hooks are called from [`energy`], so the later
//! packages fill their own files only.

use std::time::Instant;

use crate::error::RsError;
use crate::geometry::FaceGeometry;
use crate::kernel::Exponents;
use crate::penalties::Penalty;
use crate::surface::Surface;

pub mod barnes_hut;
pub mod exact;

/// Which energy the flow minimizes (`SPEC.md` §3.5). [`energy`] evaluates the tangent-point
/// energy; the other two are evaluation-only variants filled by WP7.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Energy {
    /// The tangent-point energy `ℰ̂` (default).
    #[default]
    TangentPoint,
    /// `Σ_{v≠w} a_v a_w |f_v − f_w|^{−α}` over ordered vertex pairs (paper Fig. 2; exact
    /// dense, small meshes).
    Coulomb {
        /// The exponent.
        alpha: f64,
    },
    /// The Willmore term of `SPEC.md` §5.2 in place of `ℰ̂` (paper Fig. 4, with the H²
    /// metric).
    Willmore,
}

/// Parameters of the Barnes–Hut evaluation (`SPEC.md` §3.4, §8.2; implemented by WP3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BhOptions {
    /// Admissibility parameter: node `C` is far from leaf `ℓ` when
    /// `max(r_ℓ, r_C) < θ · d_box(ℓ, C)`. `0` gives the exact energy.
    pub theta: f64,
    /// Maximal number of faces in a BVH leaf.
    pub leaf_size: usize,
}

impl Default for BhOptions {
    /// `θ = 0.5`, leaf size 8 (`SPEC.md` §3.4).
    fn default() -> Self {
        BhOptions {
            theta: 0.5,
            leaf_size: 8,
        }
    }
}

impl BhOptions {
    /// Checks `θ ≥ 0` finite and `leaf_size ≥ 1`.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] naming the option.
    pub fn validate(&self) -> Result<(), RsError> {
        if !(self.theta.is_finite() && self.theta >= 0.0) {
            return Err(RsError::InvalidOption {
                name: "theta",
                value: self.theta.to_string(),
                requirement: "finite and non-negative",
            });
        }
        if self.leaf_size == 0 {
            return Err(RsError::InvalidOption {
                name: "leaf_size",
                value: "0".into(),
                requirement: "at least 1",
            });
        }
        Ok(())
    }
}

/// How the energy and its differential are evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum EnergyEval {
    /// All ordered pairs, `O(|F|²)` ([`exact`]).
    #[default]
    Exact,
    /// Barnes–Hut with the given parameters ([`barnes_hut`], WP3).
    BarnesHut(BhOptions),
}

/// Energy, differential and counts of one evaluation without penalties (what the exact and
/// the Barnes–Hut evaluations return).
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyEvaluation {
    /// `ℰ̂(f)` or `Ẽ(f)`.
    pub energy: f64,
    /// Its differential `∂/∂f_v`.
    pub differential: Vec<[f64; 3]>,
    /// Number of ordered face pairs evaluated exactly.
    pub pairs_exact: u64,
    /// Number of (leaf, node) interactions lumped (0 for the exact evaluation).
    pub nodes_lumped: u64,
}

/// Result of [`energy`] (`SPEC.md` §9.1).
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyReport {
    /// `ℰ̂(f)` (or its Barnes–Hut approximation `Ẽ`), without penalties.
    pub energy: f64,
    /// The weighted value `w_i P_i(f)` of every penalty, in the order given.
    pub penalties: Vec<f64>,
    /// `∂F/∂f_v` for every vertex, `F = energy + Σ penalties` (the differential the flow
    /// descends).
    pub differential: Vec<[f64; 3]>,
    /// [`Exponents::subcritical`] of the exponents used.
    pub subcritical: bool,
    /// Number of ordered face pairs evaluated exactly.
    pub pairs_exact: u64,
    /// Number of (leaf, node) interactions lumped by Barnes–Hut (0 for the exact evaluation).
    pub nodes_lumped: u64,
    /// Wall time of the evaluation.
    pub seconds: f64,
}

impl EnergyReport {
    /// `energy + Σ penalties`.
    pub fn total(&self) -> f64 {
        self.energy + self.penalties.iter().sum::<f64>()
    }
}

/// Energy, penalties and differential of `positions` on `surface` (`SPEC.md` §9.1).
///
/// `threads` is the size of the explicit rayon pool (`None`: all cores); the result does not
/// depend on it, bit for bit. Positions are validated with [`Surface::validate_positions`].
///
/// [formula: rs-energy]
/// [formula: rs-energy-gradient]
///
/// # Errors
/// The errors of [`Surface::validate_positions`]; [`RsError::InvalidOption`] for
/// `threads = Some(0)` or bad [`BhOptions`]; [`RsError::NonFiniteEnergy`] naming a pair of
/// faces (e.g. coincident barycenters); the errors of the Barnes–Hut evaluation and of the
/// penalties.
pub fn energy(
    surface: &Surface,
    positions: &[[f64; 3]],
    exponents: Exponents,
    eval: &EnergyEval,
    penalties: &[Penalty],
    threads: Option<usize>,
) -> Result<EnergyReport, RsError> {
    let t0 = Instant::now();
    let geom = FaceGeometry::new(surface, positions)?;
    let pool = exact::make_pool(threads)?;
    let EnergyEvaluation {
        energy: value,
        mut differential,
        pairs_exact,
        nodes_lumped,
    } = match eval {
        EnergyEval::Exact => {
            let (energy, differential) =
                exact::energy_and_differential(&pool, surface, positions, &geom, &exponents)?;
            EnergyEvaluation {
                energy,
                differential,
                pairs_exact: exact::exact_pair_count(surface.n_faces()),
                nodes_lumped: 0,
            }
        }
        EnergyEval::BarnesHut(opts) => {
            opts.validate()?;
            barnes_hut::energy_and_differential(&pool, surface, positions, &geom, &exponents, opts)?
        }
    };
    let penalty_values =
        crate::penalties::evaluate(
            &pool,
            surface,
            positions,
            penalties,
            &exponents,
            eval,
            &mut differential,
        )?;
    Ok(EnergyReport {
        energy: value,
        penalties: penalty_values,
        differential,
        subcritical: exponents.subcritical(),
        pairs_exact,
        nodes_lumped,
        seconds: t0.elapsed().as_secs_f64(),
    })
}
