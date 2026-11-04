//! The 2D variant of Points as Tori (Fig. 10; SPEC.md §8.7, decision D15): a circle per point
//! instead of a torus, the same distance-weighted blend (Eq. 25-26 in 2D) of their signed
//! distance functions.
//!
//! Not in the authors' code. Two ways to get the per-point circle share one interface,
//! [`QuadraticPredictor`]: a *predictor* maps a point and its nearest neighbors to the three
//! coefficients `(a0, a1, a2)` of a local height function `Q(s) = a0 + a1 s + a2 s²` (the 2D
//! analogue of the network's six coefficients, SPEC.md §4.2), and [`fit_circles`] turns them into
//! circles. [`ClassicalFit`] is the classical predictor (a weighted least-squares quadratic over
//! the `k` nearest points, WP9a); the trained 2D network of WP9b implements the same trait.
//!
//! ```no_run
//! use ddg_points_as_tori::pat2d::*;
//! # fn main() -> Result<(), ddg_points_as_tori::PatError> {
//! # let (positions, normals): (Vec<[f64; 2]>, Vec<[f64; 2]>) = (vec![], vec![]);
//! let (cloud, _) = Cloud2::new(positions, normals, ddg_points_as_tori::NormalPolicy::default())?;
//! let (field, _report) = Pat2d::build(&cloud, &ClassicalFit::default(), &BuildOptions2::default())?;
//! let (phi, grad) = field.distance_gradient([0.1, 0.2])?;
//! # let _ = (phi, grad); Ok(()) }
//! ```

mod field;
mod fit;
mod grid;
mod output;

pub use field::{BuildOptions2, BuildReport2, Eval2dOptions, Eval2dReport, Pat2d, Pat2dField};
pub use fit::{
    Circles, CirclePrimitive, ClassicalFit, Fit2dOptions, Fit2dReport, Predicted2, Quad2, QuadraticPredictor, circle_sdf, circle_sdf_gradient, fit_circles,
    local_frame,
};
pub use grid::{Grid2d, Polyline, Segment2, chain_segments, evaluate_grid, marching_squares};
pub use output::{Raster, write_png, write_svg_contours};

use crate::{NormalPolicy, PatError};

/// A validated oriented 2D point cloud (finite, unit normals, `n ≥ 2`).
#[derive(Debug, Clone, PartialEq)]
pub struct Cloud2 {
    positions: Vec<[f64; 2]>,
    normals: Vec<[f64; 2]>,
}

impl Cloud2 {
    /// Validate like [`crate::PointCloud::new`]: lengths, finiteness, zero normals, unit length
    /// under [`NormalPolicy::RequireUnit`] (or normalized, counted in the returned number, under
    /// [`NormalPolicy::Normalize`]).
    pub fn new(positions: Vec<[f64; 2]>, mut normals: Vec<[f64; 2]>, policy: NormalPolicy) -> Result<(Cloud2, usize), PatError> {
        if positions.len() != normals.len() {
            return Err(PatError::LengthMismatch { positions: positions.len(), normals: normals.len() });
        }
        if positions.len() < 2 {
            return Err(PatError::TooFewPoints { n: positions.len(), need: 2 });
        }
        if let Some(i) = positions.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
            return Err(PatError::NonFinite { index: i, what: "position" });
        }
        if let Some(i) = normals.iter().position(|p| p.iter().any(|c| !c.is_finite())) {
            return Err(PatError::NonFinite { index: i, what: "normal" });
        }
        let mut normalized = 0;
        for (i, n) in normals.iter_mut().enumerate() {
            let l = (n[0] * n[0] + n[1] * n[1]).sqrt();
            if l == 0.0 {
                return Err(PatError::ZeroNormal { index: i });
            }
            match policy {
                NormalPolicy::RequireUnit { tol } => {
                    if (l - 1.0).abs() > tol {
                        return Err(PatError::NonUnitNormal { index: i, norm: l });
                    }
                }
                NormalPolicy::Normalize => {
                    if l != 1.0 {
                        *n = [n[0] / l, n[1] / l];
                        normalized += 1;
                    }
                }
            }
        }
        Ok((Cloud2 { positions, normals }, normalized))
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// True when empty (never, after validation).
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Positions.
    pub fn positions(&self) -> &[[f64; 2]] {
        &self.positions
    }

    /// Unit normals.
    pub fn normals(&self) -> &[[f64; 2]] {
        &self.normals
    }

    /// The positions embedded in the plane `z = 0` (the kd-tree of the crate is 3D).
    pub(crate) fn embedded(&self) -> Vec<[f64; 3]> {
        self.positions.iter().map(|p| [p[0], p[1], 0.0]).collect()
    }
}
