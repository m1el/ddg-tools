//! The problem description: which energy is minimized, with which weights, and the options of one
//! evaluation. Everything here is plain data that the caller may change between evaluations.

use std::fmt;
use std::sync::Arc;

use ddg_core::sparse::Parallelism;

use crate::error::CutError;

/// The number of regions this version supports.
pub const SUPPORTED_REGIONS: usize = 2;

/// The smoothing constant `c = 0.02` the paper uses for the curvature-based weight (Sec. 5.6).
pub const PAPER_FEATURE_SMOOTHING: f64 = 0.02;

/// The scaling of the input surface to unit total area (SPEC §2).
///
/// The method runs on the surface scaled by `length_scale = 1/√(total area)`. Every parameter of
/// a [`CutProblem`] and every number in an [`Evaluation`](crate::Evaluation) refers to that
/// scaled surface unless its documentation says otherwise; this struct converts between the two
/// units. Dimensionless quantities (the log scale factor `u`, the Dirichlet energy, integrated
/// normal derivatives) are the same in both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitAreaScale {
    /// Total area of the input surface, in the caller's unit.
    pub total_area: f64,
    /// The factor `s = 1/√(total area)`: a length of the caller times `s` is a length on the
    /// unit-area surface.
    pub length_scale: f64,
}

impl UnitAreaScale {
    /// The scaling for a surface of the given total area.
    ///
    /// # Errors
    /// [`CutError::InvalidTotalArea`] unless the area and the resulting scale are positive and
    /// finite.
    pub fn from_total_area(total_area: f64) -> Result<UnitAreaScale, CutError> {
        let length_scale = 1.0 / total_area.sqrt();
        if total_area > 0.0 && total_area.is_finite() && length_scale.is_finite() {
            Ok(UnitAreaScale {
                total_area,
                length_scale,
            })
        } else {
            Err(CutError::InvalidTotalArea { total_area })
        }
    }

    /// Converts a length (or a level-set value) from the caller's unit to the unit-area surface.
    #[inline]
    pub fn length_to_unit_area(&self, length: f64) -> f64 {
        length * self.length_scale
    }

    /// Converts a length (or a level-set value) from the unit-area surface to the caller's unit.
    #[inline]
    pub fn length_from_unit_area(&self, length: f64) -> f64 {
        length / self.length_scale
    }

    /// Converts an area from the caller's unit to the unit-area surface.
    #[inline]
    pub fn area_to_unit_area(&self, area: f64) -> f64 {
        area / self.total_area
    }

    /// Converts an area from the unit-area surface to the caller's unit.
    #[inline]
    pub fn area_from_unit_area(&self, area: f64) -> f64 {
        area * self.total_area
    }
}

/// A pointwise energy density `b(u)` and its derivative `b′(u)` (paper Sec. 3.5.2):
/// `J = ∫ b(u) dA`.
#[derive(Clone)]
pub struct PointwiseEnergy {
    value: Arc<dyn Fn(f64) -> f64 + Send + Sync>,
    derivative: Arc<dyn Fn(f64) -> f64 + Send + Sync>,
}

impl PointwiseEnergy {
    /// Wraps the density `b` and its derivative `b′`. The two must be consistent; nothing checks
    /// that. Both must return finite values for every log scale factor that occurs.
    pub fn new(
        value: impl Fn(f64) -> f64 + Send + Sync + 'static,
        derivative: impl Fn(f64) -> f64 + Send + Sync + 'static,
    ) -> PointwiseEnergy {
        PointwiseEnergy {
            value: Arc::new(value),
            derivative: Arc::new(derivative),
        }
    }

    /// `b(u)`.
    #[inline]
    pub fn value(&self, u: f64) -> f64 {
        (self.value)(u)
    }

    /// `b′(u)`.
    #[inline]
    pub fn derivative(&self, u: f64) -> f64 {
        (self.derivative)(u)
    }
}

impl fmt::Debug for PointwiseEnergy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PointwiseEnergy")
            .field("b(0)", &self.value(0.0))
            .field("b'(0)", &self.derivative(0.0))
            .finish_non_exhaustive()
    }
}

/// The weight `a > 0` of the scaled Dirichlet energy `∫ a |∇u|²` (paper Eq. 17), one value per
/// vertex of the input mesh, interpolated linearly.
#[derive(Debug, Clone, PartialEq)]
pub enum ScaleWeight {
    /// The curvature-based weight of the paper (Sec. 3.5.1, 5.6; SPEC §3.6): `a = r³` with the
    /// smoothed feature size `r`. The evaluator computes it once and caches it.
    FeatureSize {
        /// The cap `r_ε > 0` on the feature size in flat regions, a length on the unit-area
        /// surface. The paper gives no value (SPEC D9f).
        r_eps: f64,
        /// The smoothing constant `c ≥ 0` of `(W_V + c L) r = W_V r⁰`, on the unit-area surface.
        /// The paper uses [`PAPER_FEATURE_SMOOTHING`].
        smoothing: f64,
    },
    /// One positive value per vertex, used as it is on the unit-area surface.
    PerVertex(Vec<f64>),
}

/// The distortion energy `J` of one region (SPEC §1).
#[derive(Debug, Clone)]
pub enum Energy {
    /// `J = ∫ |∇u|²` (paper Eq. 2). Needs no adjoint solve.
    Dirichlet,
    /// `J = ∫ a |∇u|²` (paper Eq. 17).
    ScaledDirichlet(ScaleWeight),
    /// `J = ∫ u²` (paper Eq. 19), the pointwise energy with `b(u) = u²`.
    Hencky,
    /// `J = ∫ b(u)` for a caller-supplied `b` (paper Sec. 3.5.2).
    Pointwise(PointwiseEnergy),
}

/// Where the right-hand side `K` of the Yamabe equation `Δu = −K` comes from (SPEC §3.1).
#[derive(Debug, Clone, PartialEq)]
pub enum CurvatureSource {
    /// The Gaussian curvature of the mesh, distributed to faces (paper Sec. 5.2).
    Gaussian,
    /// The same density everywhere (paper Sec. 7.4, "uniform density curves"), per unit area of
    /// the unit-area surface.
    Uniform(f64),
    /// One density per face of the input mesh (paper Fig. 18), per unit area of the unit-area
    /// surface.
    PerFace(Vec<f64>),
}

/// The weights of one region in the energy (SPEC §1, paper Fig. 14).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegionWeights {
    /// `w_r ≥ 0`: the weight of the distortion `J_r` of the region. Default 1.
    pub distortion: f64,
    /// `β_r`: the weight of the area of the region (on the unit-area surface). Default 0.
    pub area: f64,
}

impl Default for RegionWeights {
    fn default() -> Self {
        RegionWeights {
            distortion: 1.0,
            area: 0.0,
        }
    }
}

/// The energy whose gradient the flow follows (SPEC §1, §6):
///
/// ```text
/// E = Σ_r [ w_r · J_r(u) + β_r · Area(S_r) ]  +  ∫_γ (α_L + η) ds
/// ```
///
/// on the surface scaled to unit area ([`UnitAreaScale`]). The caller may change any field
/// between two evaluations.
#[derive(Debug, Clone)]
pub struct CutProblem {
    /// The distortion energy `J`.
    pub energy: Energy,
    /// The source term of the Yamabe equation.
    pub curvature: CurvatureSource,
    /// The weights `w_r`, `β_r` of every region; its length is the number of regions and must be
    /// [`SUPPORTED_REGIONS`]. Region 0 is where the level-set function is `≤ 0`.
    pub regions: Vec<RegionWeights>,
    /// The length weight `α_L ≥ 0` (paper Eq. 5), per unit length on the unit-area surface.
    pub length_weight: f64,
    /// The penalty `η_i ≥ 0` per vertex of the input mesh (paper Sec. 7.1), per unit length on
    /// the unit-area surface; `None` means 0 everywhere.
    pub penalty: Option<Vec<f64>>,
}

impl CutProblem {
    /// The given energy with Gaussian curvature, two regions with default weights, the given
    /// length weight and no penalty.
    pub fn new(energy: Energy, length_weight: f64) -> CutProblem {
        CutProblem {
            energy,
            curvature: CurvatureSource::Gaussian,
            regions: vec![RegionWeights::default(); SUPPORTED_REGIONS],
            length_weight,
            penalty: None,
        }
    }
}

/// How the normal derivatives `∂u/∂n` and `∂p/∂n` at the cut vertices are obtained from the
/// integrated normal derivatives `v`, `q` and the dual lengths `ℓ` (SPEC §3.5, D12).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FluxRecovery {
    /// `∂u/∂n = v_p / ℓ_p` at every cut vertex on its own, as in the paper (Sec. 6.2). Where two
    /// or more consecutive crossings lie close to one mesh vertex, `ℓ_p` is of the order of
    /// their distance while the error of `v_p` is of the order of the mesh size, and the
    /// quotient is wrong by a factor of up to `1/min_crossing_offset` (measured: errors of
    /// `10⁴ %`). The paper repairs this with sliver smoothing
    /// ([`EvalOptions::smoothing_passes`]), which does not help when the neighbours are slivers
    /// too.
    Nodal,
    /// The nodal fluxes are smoothed along the cut before they are divided: in every region the
    /// values `d_p` at its cut vertices solve
    ///
    /// ```text
    /// ℓ_p d_p + Σ_{segments σ = pq at p} κ_σ (d_p − d_q) = v_p,     κ_σ = (c · h_f)² / |pq|
    /// ```
    ///
    /// with `h_f²` the mean squared edge length of the face that contains `σ`. This is one
    /// backward Euler step of the 1D heat flow on the cut with the lumped mass `ℓ`, applied to
    /// `v/ℓ`, with the local time `(c h_f)²`: it conserves `Σ ℓ_p d_p = Σ v_p` on every piece of
    /// the cut, keeps `d` within the range of the nodal values, and pools cut vertices that are
    /// much closer to each other than `c h_f`. `length = 0` is [`FluxRecovery::Nodal`].
    Smoothed {
        /// The smoothing length `c ≥ 0` in units of the local edge length. Default `0.5`.
        length: f64,
    },
}

impl Default for FluxRecovery {
    fn default() -> Self {
        FluxRecovery::Smoothed { length: 0.5 }
    }
}

/// Options of one evaluation.
#[derive(Debug, Clone, Copy)]
pub struct EvalOptions<'a> {
    /// `δ` of SPEC §3.2: a crossing closer than this (as a fraction of its edge) to a vertex is
    /// moved to that distance, and counted in
    /// [`EvalReport::n_clamped_crossings`](crate::EvalReport::n_clamped_crossings). Must lie in
    /// `[0, ½)`; with 0 nothing is moved and a crossing exactly at a vertex is an error.
    /// Default `10⁻³`.
    pub min_crossing_offset: f64,
    /// How `∂u/∂n` and `∂p/∂n` on the cut are recovered from the integrated normal derivatives
    /// (SPEC §3.5). Default: [`FluxRecovery::Smoothed`] with length `0.5`.
    pub flux_recovery: FluxRecovery,
    /// `ρ₀` of SPEC §3.7: a cut vertex with an incident cut-mesh triangle of aspect ratio above
    /// this is a sliver vertex. Must be positive (`+∞` means no vertex is a sliver; every
    /// triangle has an aspect ratio of at least 1). Default 5 (paper Sec. 6.2).
    pub sliver_aspect_ratio: f64,
    /// Number of rounds of the paper's sliver smoothing of the speed (SPEC §3.7); 0 switches it
    /// off. Default 0: with [`FluxRecovery::Smoothed`] it brings nothing measurable.
    /// [`EvalOptions::paper`] sets 1 (the paper does not give a number, SPEC D7).
    pub smoothing_passes: usize,
    /// Where factorizations and solves run. Default: on the calling thread.
    pub parallelism: Parallelism<'a>,
}

impl Default for EvalOptions<'_> {
    fn default() -> Self {
        EvalOptions {
            min_crossing_offset: 1e-3,
            flux_recovery: FluxRecovery::default(),
            sliver_aspect_ratio: 5.0,
            smoothing_passes: 0,
            parallelism: Parallelism::Sequential,
        }
    }
}

impl EvalOptions<'_> {
    /// The procedure of the paper (Sec. 6.2): [`FluxRecovery::Nodal`] and one round of sliver
    /// smoothing. Not robust when the cut passes close to a mesh vertex; see [`FluxRecovery`].
    pub fn paper() -> Self {
        EvalOptions {
            flux_recovery: FluxRecovery::Nodal,
            smoothing_passes: 1,
            ..EvalOptions::default()
        }
    }
}

impl EvalOptions<'_> {
    /// Checks the ranges documented on the fields.
    pub(crate) fn validate(&self) -> Result<(), CutError> {
        let d = self.min_crossing_offset;
        if !(d >= 0.0 && d < 0.5) {
            return Err(CutError::InvalidOptions {
                reason: format!("min_crossing_offset must lie in [0, 0.5), got {d}"),
            });
        }
        let r = self.sliver_aspect_ratio;
        if !(r > 0.0) {
            return Err(CutError::InvalidOptions {
                reason: format!("sliver_aspect_ratio must be positive, got {r}"),
            });
        }
        if let FluxRecovery::Smoothed { length } = self.flux_recovery
            && !(length.is_finite() && length >= 0.0)
        {
            return Err(CutError::InvalidOptions {
                reason: format!(
                    "the smoothing length of the flux recovery must be finite and non-negative, \
                     got {length}"
                ),
            });
        }
        Ok(())
    }
}

/// Checks that `array` has one entry per `per`.
pub(crate) fn check_length(
    what: &'static str,
    per: &'static str,
    got: usize,
    expected: usize,
) -> Result<(), CutError> {
    if got == expected {
        Ok(())
    } else {
        Err(CutError::WrongLength {
            what,
            per,
            expected,
            got,
        })
    }
}

impl CutProblem {
    /// Validates everything that does not depend on the level-set function.
    pub(crate) fn validate(&self, n_vertices: usize, n_faces: usize) -> Result<(), CutError> {
        if self.regions.len() != SUPPORTED_REGIONS {
            return Err(CutError::UnsupportedRegionCount {
                got: self.regions.len(),
            });
        }
        for (r, w) in self.regions.iter().enumerate() {
            if !(w.distortion.is_finite() && w.distortion >= 0.0) {
                return Err(CutError::InvalidRegionWeight {
                    region: r as u32,
                    what: "distortion",
                    value: w.distortion,
                });
            }
            if !w.area.is_finite() {
                return Err(CutError::InvalidRegionWeight {
                    region: r as u32,
                    what: "area",
                    value: w.area,
                });
            }
        }
        if !(self.length_weight.is_finite() && self.length_weight >= 0.0) {
            return Err(CutError::InvalidLengthWeight {
                value: self.length_weight,
            });
        }
        if let Some(eta) = &self.penalty {
            check_length("penalty", "vertices", eta.len(), n_vertices)?;
            if let Some(v) = eta.iter().position(|x| !(x.is_finite() && *x >= 0.0)) {
                return Err(CutError::InvalidPenalty {
                    vertex: v as u32,
                    value: eta[v],
                });
            }
        }
        match &self.curvature {
            CurvatureSource::Gaussian => {}
            CurvatureSource::Uniform(c) => {
                if !c.is_finite() {
                    return Err(CutError::InvalidCurvature {
                        face: None,
                        value: *c,
                    });
                }
            }
            CurvatureSource::PerFace(k) => {
                check_length("curvature override", "faces", k.len(), n_faces)?;
                if let Some(f) = k.iter().position(|x| !x.is_finite()) {
                    return Err(CutError::InvalidCurvature {
                        face: Some(f as u32),
                        value: k[f],
                    });
                }
            }
        }
        match &self.energy {
            Energy::ScaledDirichlet(ScaleWeight::PerVertex(a)) => {
                check_length("scaled Dirichlet weight", "vertices", a.len(), n_vertices)?;
                if let Some(v) = a.iter().position(|x| !(x.is_finite() && *x > 0.0)) {
                    return Err(CutError::InvalidScaleWeight {
                        vertex: v as u32,
                        value: a[v],
                    });
                }
            }
            Energy::ScaledDirichlet(ScaleWeight::FeatureSize { r_eps, smoothing }) => {
                validate_feature_size(*r_eps, *smoothing)?;
            }
            Energy::Dirichlet | Energy::Hencky | Energy::Pointwise(_) => {}
        }
        Ok(())
    }
}

/// Checks the parameters of the curvature-based weight.
pub(crate) fn validate_feature_size(r_eps: f64, smoothing: f64) -> Result<(), CutError> {
    if !(r_eps.is_finite() && r_eps > 0.0) {
        return Err(CutError::InvalidFeatureSize {
            what: "r_eps must be positive and finite",
            value: r_eps,
        });
    }
    if !(smoothing.is_finite() && smoothing >= 0.0) {
        return Err(CutError::InvalidFeatureSize {
            what: "the smoothing constant must be finite and non-negative",
            value: smoothing,
        });
    }
    Ok(())
}
