//! Morphology (SPEC.md §8.3, Fig. 24): erosion, dilation, opening and closing of the shape of a
//! field by refitting tori to its offset level sets. Owned by WP6.
//!
//! Erosion by `c` samples the `−c` level set (`φ = −c`, an exact SDF's inward offset), takes the
//! normals `∇φ/‖∇φ‖` and refits tori (network + fit) to the new cloud; dilation uses `+c`.
//! Opening is erosion then dilation, closing is dilation then erosion; every step refits.

use super::{SampleOptions, SampleReport, sample_level_set};
use crate::csg::ImplicitField;
use crate::{BuildOptions, BuildReport, Model, NormalPolicy, PatError, PointCloud, PointsAsTori, TorusField};
use std::time::Instant;

/// A morphological operation with radius `c ≥ 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Morph {
    /// Shrink by `c`.
    Erode(f64),
    /// Grow by `c`.
    Dilate(f64),
    /// Erode, then dilate.
    Open(f64),
    /// Dilate, then erode.
    Close(f64),
}

impl Morph {
    /// The levels sampled in order (`−c` for erosion, `+c` for dilation, two for opening and
    /// closing), after checking the radius (finite, `≥ 0`).
    pub fn levels(self) -> Result<Vec<f64>, PatError> {
        let (c, l) = match self {
            Morph::Erode(c) => (c, vec![-c]),
            Morph::Dilate(c) => (c, vec![c]),
            Morph::Open(c) => (c, vec![-c, c]),
            Morph::Close(c) => (c, vec![c, -c]),
        };
        if !(c.is_finite() && c >= 0.0) {
            return Err(PatError::InvalidOption { what: format!("morphology radius c = {c} (finite and >= 0)") });
        }
        Ok(l)
    }
}

/// Options of [`morphology`].
#[derive(Debug, Clone, PartialEq)]
pub struct MorphOptions {
    /// Level-set sampling (`count: None` keeps the input's point count when the field knows it).
    pub sample: SampleOptions,
    /// Network, fitting and (unless `paper_eval`) evaluation options of every refit.
    pub build: BuildOptions,
    /// Use the paper's evaluation options (`λ = 10³/D`, radius `2C/λ`) for each refit cloud.
    pub paper_eval: bool,
}

impl Default for MorphOptions {
    fn default() -> Self {
        MorphOptions { sample: SampleOptions::default(), build: BuildOptions::default(), paper_eval: true }
    }
}

/// One refit.
#[derive(Debug, Clone, PartialEq)]
pub struct MorphStep {
    /// The level set sampled (`−c` or `+c`).
    pub level: f64,
    /// Sampling.
    pub sample: SampleReport,
    /// The refit.
    pub build: BuildReport,
}

/// What [`morphology`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct MorphReport {
    /// One entry per refit.
    pub steps: Vec<MorphStep>,
    /// Seconds.
    pub seconds: f64,
}

/// Apply `op` to the shape of `field` and refit tori to the result. Returns the new field, the
/// cloud that was sampled from the morphed surface (the one the field was built on), and a report.
pub fn morphology(
    field: &dyn ImplicitField,
    model: &Model,
    op: Morph,
    opts: &MorphOptions,
) -> Result<(TorusField, PointCloud, MorphReport), PatError> {
    let t0 = Instant::now();
    let levels = op.levels()?;
    let mut sample = opts.sample.clone();
    if sample.count.is_none() {
        sample.count = field.point_count_hint();
    }
    let mut steps = Vec::new();
    let mut current: Option<(TorusField, PointCloud)> = None;
    for level in levels {
        let src: &dyn ImplicitField = match &current {
            Some((f, _)) => f,
            None => field,
        };
        let (pts, nrm, srep) = sample_level_set(src, level, &sample)?;
        let (cloud, _) = PointCloud::new(pts, nrm, NormalPolicy::Normalize)?;
        let mut build = opts.build.clone();
        if opts.paper_eval {
            build.eval = TorusField::paper_options(&cloud, build.eval.threads)?.0;
        }
        let (f, brep) = PointsAsTori::build(&cloud, model, &build)?;
        steps.push(MorphStep { level, sample: srep, build: brep });
        current = Some((f, cloud));
    }
    match current {
        Some((f, cloud)) => Ok((f, cloud, MorphReport { steps, seconds: t0.elapsed().as_secs_f64() })),
        None => Err(PatError::InvalidOption { what: "morphology without steps".into() }),
    }
}
