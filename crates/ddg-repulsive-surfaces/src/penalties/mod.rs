//! Penalties of `SPEC.md` §5.2 (area/volume deviation, total area/volume, squared error,
//! boundary length and curvature, Willmore, implicit obstacles and attractors, mesh and
//! point-cloud obstacles), each `w · P(f)` with an exact gradient (D24).
//!
//! `evaluate` is the hook of [`crate::energy::energy`]; [`evaluate_one`] evaluates a single
//! penalty (weighted value and weighted gradient). **A penalty whose value is not finite
//! (an implicit obstacle with `d ≤ 0` gives `+∞`; coincident points of an obstacle) adds no
//! gradient**: the caller (line search) rejects the trial.
//!
//! Penalties that need the exponents or the `θ` of the energy (tangent-point obstacles) take
//! them from the `energy()` call (`θ = 0` for the exact evaluation).

use std::fmt;
use std::sync::Arc;

use crate::energy::{BhOptions, EnergyEval};
use crate::error::RsError;
use crate::geometry::v3::axpy;
use crate::implicit::ImplicitSurface;
use crate::kernel::Exponents;
use crate::obstacle::{CrossEvaluator, ObstacleKind, StaticSet};
use crate::surface::Surface;

pub mod potentials;
pub mod simple;
pub mod willmore;

pub use potentials::{Quadrature, coulomb_energy};
pub use willmore::{willmore_and_gradient, willmore_matrix};

/// A soft penalty `w · P(f)` (`SPEC.md` §5.2). Targets that are "the initial value" in the
/// constraints are explicit numbers here (the session resolves them).
#[derive(Clone)]
#[non_exhaustive]
pub enum Penalty {
    /// `w (A/A₀ − 1)²`.
    AreaDeviation {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// `A₀ ≠ 0`.
        target: f64,
    },
    /// `w (V/V₀ − 1)²` (signed volume; needs consistent orientation).
    VolumeDeviation {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// `V₀ ≠ 0`.
        target: f64,
    },
    /// `w A^power`.
    TotalArea {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// The power (1: minimal surfaces; 2: the reference's `minimize area`).
        power: f64,
    },
    /// `w V^power` (signed volume; needs consistent orientation).
    TotalVolume {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// The power.
        power: f64,
    },
    /// `w Σ_v |f_v − f_v⁰|²`.
    SquaredError {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// `f_v⁰` per vertex.
        targets: Vec<[f64; 3]>,
    },
    /// `w (L − Σ_{e∈∂M} ℓ_e)²`.
    BoundaryLength {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// `L`.
        target: f64,
    },
    /// `w Σ_{v∈∂M} θ_v² / ℓ_v` (`dual_length`) or `w Σ θ_v²`.
    BoundaryCurvature {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// Divide by the dual length `ℓ_v`.
        dual_length: bool,
    },
    /// `w fᵀ W M⁻¹ W f` (interior rows).
    Willmore {
        /// Weight `w ≥ 0`.
        weight: f64,
    },
    /// `w Σ_v q_v d(f_v)^{−power}`; `+∞` for `d ≤ 0`.
    ImplicitObstacle {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// The surface (negative inside).
        surface: Arc<dyn ImplicitSurface>,
        /// The power `p > 0`.
        power: f64,
        /// The weights `q_v`.
        quadrature: Quadrature,
    },
    /// `w Σ_{v∈set} q_v |d(f_v)|^{power}`.
    ImplicitAttractor {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// The surface.
        surface: Arc<dyn ImplicitSurface>,
        /// The power `p > 0`.
        power: f64,
        /// The weights `q_v`.
        quadrature: Quadrature,
        /// The attracted vertices (`None`: all).
        vertices: Option<Vec<u32>>,
    },
    /// A triangle-soup obstacle (`SPEC.md` §5.3).
    MeshObstacle {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// The obstacle ([`StaticSet::from_obstacle`]).
        obstacle: Arc<StaticSet>,
        /// Tangent-point cross energy or the printed Coulomb potential.
        kind: ObstacleKind,
    },
    /// A point cloud, one-sided tangent-point energy `Σ_S Σ_p K(X_S, N_S; y_p) a_S`.
    PointCloudObstacle {
        /// Weight `w ≥ 0`.
        weight: f64,
        /// The cloud ([`StaticSet::from_points`]).
        points: Arc<StaticSet>,
    },
}

impl Penalty {
    /// The weight `w`.
    pub fn weight(&self) -> f64 {
        match self {
            Penalty::AreaDeviation { weight, .. }
            | Penalty::VolumeDeviation { weight, .. }
            | Penalty::TotalArea { weight, .. }
            | Penalty::TotalVolume { weight, .. }
            | Penalty::SquaredError { weight, .. }
            | Penalty::BoundaryLength { weight, .. }
            | Penalty::BoundaryCurvature { weight, .. }
            | Penalty::Willmore { weight }
            | Penalty::ImplicitObstacle { weight, .. }
            | Penalty::ImplicitAttractor { weight, .. }
            | Penalty::MeshObstacle { weight, .. }
            | Penalty::PointCloudObstacle { weight, .. } => *weight,
        }
    }

    /// The name of the variant.
    pub fn name(&self) -> &'static str {
        match self {
            Penalty::AreaDeviation { .. } => "AreaDeviation",
            Penalty::VolumeDeviation { .. } => "VolumeDeviation",
            Penalty::TotalArea { .. } => "TotalArea",
            Penalty::TotalVolume { .. } => "TotalVolume",
            Penalty::SquaredError { .. } => "SquaredError",
            Penalty::BoundaryLength { .. } => "BoundaryLength",
            Penalty::BoundaryCurvature { .. } => "BoundaryCurvature",
            Penalty::Willmore { .. } => "Willmore",
            Penalty::ImplicitObstacle { .. } => "ImplicitObstacle",
            Penalty::ImplicitAttractor { .. } => "ImplicitAttractor",
            Penalty::MeshObstacle { .. } => "MeshObstacle",
            Penalty::PointCloudObstacle { .. } => "PointCloudObstacle",
        }
    }

    /// Checks the options of the penalty against the surface.
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] naming the option;
    /// [`RsError::VolumeNeedsOrientation`]; [`RsError::PositionCountMismatch`] for
    /// `SquaredError` targets of the wrong length.
    pub fn validate(&self, surface: &Surface) -> Result<(), RsError> {
        let bad = |name: &'static str, x: f64, requirement: &'static str| {
            Err(RsError::InvalidOption {
                name,
                value: x.to_string(),
                requirement,
            })
        };
        let w = self.weight();
        if !(w.is_finite() && w >= 0.0) {
            return bad("weight", w, "finite and non-negative");
        }
        let positive = |name: &'static str, x: f64| {
            if x.is_finite() && x > 0.0 {
                Ok(())
            } else {
                bad(name, x, "a finite positive number")
            }
        };
        let nonzero = |name: &'static str, x: f64| {
            if x.is_finite() && x != 0.0 {
                Ok(())
            } else {
                bad(name, x, "finite and non-zero")
            }
        };
        match self {
            Penalty::AreaDeviation { target, .. } => nonzero("target area", *target),
            Penalty::VolumeDeviation { target, .. } => {
                nonzero("target volume", *target)?;
                simple::require_orientation(surface)
            }
            Penalty::TotalArea { power, .. } => positive("power", *power),
            Penalty::TotalVolume { power, .. } => {
                positive("power", *power)?;
                simple::require_orientation(surface)
            }
            Penalty::SquaredError { targets, .. } => {
                if targets.len() != surface.n_vertices() {
                    return Err(RsError::PositionCountMismatch {
                        expected: surface.n_vertices(),
                        found: targets.len(),
                    });
                }
                match targets.iter().flatten().find(|x| !x.is_finite()) {
                    Some(x) => bad("targets", *x, "finite coordinates"),
                    None => Ok(()),
                }
            }
            Penalty::BoundaryLength { target, .. } => {
                if target.is_finite() && *target >= 0.0 {
                    Ok(())
                } else {
                    bad("target length", *target, "finite and non-negative")
                }
            }
            Penalty::BoundaryCurvature { .. } | Penalty::Willmore { .. } => Ok(()),
            Penalty::ImplicitObstacle { power, .. } | Penalty::ImplicitAttractor { power, .. } => {
                positive("power", *power)
            }
            Penalty::MeshObstacle { obstacle, kind, .. } => {
                if !obstacle.has_normals() {
                    return Err(RsError::ObstacleError {
                        obstacle: 0,
                        message: "a mesh obstacle needs triangles (use PointCloudObstacle)".into(),
                    });
                }
                match kind {
                    ObstacleKind::Coulomb { power } => positive("power", *power),
                    ObstacleKind::TangentPoint => Ok(()),
                }
            }
            Penalty::PointCloudObstacle { points, .. } => {
                if points.has_normals() {
                    Err(RsError::ObstacleError {
                        obstacle: 0,
                        message: "a point-cloud obstacle must be built with StaticSet::from_points"
                            .into(),
                    })
                } else {
                    Ok(())
                }
            }
        }
    }
}

impl fmt::Debug for Penalty {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct(self.name())
            .field("weight", &self.weight())
            .finish_non_exhaustive()
    }
}

fn chain(value: f64, d: f64, mut grad: Vec<[f64; 3]>) -> (f64, Vec<[f64; 3]>) {
    for g in &mut grad {
        *g = crate::geometry::v3::scale(*g, d);
    }
    (value, grad)
}

/// The weighted value `w P(f)` and weighted gradient `w ∇P(f)` of one penalty (the gradient is
/// zero when the value is not finite). `exponents` and `eval` are those of the energy, used by
/// tangent-point obstacles (`θ = 0` for [`EnergyEval::Exact`]).
///
/// [formula: rs-penalties]
/// [formula: rs-willmore]
///
/// # Errors
/// The errors of [`Penalty::validate`] and [`Surface::validate_positions`].
pub fn evaluate_one(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    penalty: &Penalty,
    exponents: &Exponents,
    eval: &EnergyEval,
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    penalty.validate(surface)?;
    surface.validate_positions(positions)?;
    let w = penalty.weight();
    let (value, grad) = match penalty {
        Penalty::AreaDeviation { target, .. } => {
            let (a, g) = simple::area_and_gradient(surface, positions);
            let dev = a / target - 1.0;
            chain(dev * dev, 2.0 * dev / target, g)
        }
        Penalty::VolumeDeviation { target, .. } => {
            let (v, g) = simple::volume_and_gradient(surface, positions);
            let dev = v / target - 1.0;
            chain(dev * dev, 2.0 * dev / target, g)
        }
        Penalty::TotalArea { power, .. } => {
            let (a, g) = simple::area_and_gradient(surface, positions);
            let (p, dp) = simple::power_and_derivative(a, *power);
            chain(p, dp, g)
        }
        Penalty::TotalVolume { power, .. } => {
            let (v, g) = simple::volume_and_gradient(surface, positions);
            let (p, dp) = simple::power_and_derivative(v, *power);
            chain(p, dp, g)
        }
        Penalty::SquaredError { targets, .. } => {
            simple::squared_error_and_gradient(positions, targets)?
        }
        Penalty::BoundaryLength { target, .. } => {
            simple::boundary_length_penalty(surface, positions, *target)
        }
        Penalty::BoundaryCurvature { dual_length, .. } => {
            simple::boundary_curvature_and_gradient(surface, positions, *dual_length)
        }
        Penalty::Willmore { .. } => willmore::willmore_and_gradient(surface, positions)?,
        Penalty::ImplicitObstacle {
            surface: s,
            power,
            quadrature,
            ..
        } => potentials::implicit_obstacle(surface, positions, s.as_ref(), *power, *quadrature)?,
        Penalty::ImplicitAttractor {
            surface: s,
            power,
            quadrature,
            vertices,
            ..
        } => potentials::implicit_attractor(
            surface,
            positions,
            s.as_ref(),
            *power,
            *quadrature,
            vertices.as_deref(),
        )?,
        Penalty::MeshObstacle { obstacle, kind, .. } => {
            cross(pool, surface, positions, obstacle, *kind, exponents, eval)?
        }
        Penalty::PointCloudObstacle { points, .. } => cross(
            pool,
            surface,
            positions,
            points,
            ObstacleKind::TangentPoint,
            exponents,
            eval,
        )?,
    };
    let wv = w * value;
    if !wv.is_finite() || !grad.iter().flatten().all(|x| x.is_finite()) {
        return Ok((wv, vec![[0.0; 3]; surface.n_vertices()]));
    }
    Ok((wv, chain(0.0, w, grad).1))
}

fn cross(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    set: &Arc<StaticSet>,
    kind: ObstacleKind,
    exponents: &Exponents,
    eval: &EnergyEval,
) -> Result<(f64, Vec<[f64; 3]>), RsError> {
    let opts = match eval {
        EnergyEval::Exact => BhOptions {
            theta: 0.0,
            ..BhOptions::default()
        },
        EnergyEval::BarnesHut(o) => *o,
    };
    let ev = CrossEvaluator::build(pool, surface, positions, set.clone(), kind, exponents, opts)?;
    Ok(ev.energy_and_differential(pool))
}

/// Hook for [`crate::energy::energy`]: returns the weighted value `w_i P_i(f)` of every
/// penalty (in order) and adds the weighted gradients to `differential`. Positions are
/// validated.
///
/// # Errors
/// The errors of [`evaluate_one`].
pub(crate) fn evaluate(
    pool: &rayon::ThreadPool,
    surface: &Surface,
    positions: &[[f64; 3]],
    penalties: &[Penalty],
    exponents: &Exponents,
    eval: &EnergyEval,
    differential: &mut [[f64; 3]],
) -> Result<Vec<f64>, RsError> {
    let mut values = Vec::with_capacity(penalties.len());
    for p in penalties {
        let (v, g) = evaluate_one(pool, surface, positions, p, exponents, eval)?;
        for (d, x) in differential.iter_mut().zip(&g) {
            *d = axpy(*d, 1.0, *x);
        }
        values.push(v);
    }
    Ok(values)
}
