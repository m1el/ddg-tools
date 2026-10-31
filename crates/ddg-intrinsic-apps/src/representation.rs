//! Which representation of an intrinsic triangulation an application runs on (Integer
//! Coordinates contract, WP11; Gillespie, Sharp, Crane 2021).
//!
//! The applications that need only the lengths and the connectivity are generic over
//! [`ddg_intrinsic::Triangulation`]; the ones that also locate points of the input mesh take
//! an [`InputPoints`]. [`Representation::DEFAULT`] is **signposts** for now. After the
//! evaluation of the integer representation (WP10) one commit changes this constant to
//! [`Representation::Integer`]; nothing else here depends on it.

use ddg_integer_coords::IntegerTriangulation;
use ddg_intrinsic::{IntrinsicError, IntrinsicTriangulation, SurfacePoint, Triangulation};

/// The representation of an intrinsic triangulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Representation {
    /// Signposts (Sharp, Soliman, Crane 2019): `ddg_intrinsic::IntrinsicTriangulation`.
    Signposts,
    /// Integer coordinates (Gillespie, Sharp, Crane 2021): `ddg_integer_coords::IntegerTriangulation`.
    Integer,
}

impl Representation {
    /// The default of the consumers. Signposts until the evaluation (WP10) passes; then
    /// `Integer` (Integer Coordinates SPEC §13 Q3).
    pub const DEFAULT: Representation = Representation::Signposts;
}

impl Default for Representation {
    fn default() -> Self {
        Representation::DEFAULT
    }
}

/// A triangulation that can say where a point of the input mesh is on it (`T⁰ → T¹`).
pub trait InputPoints: Triangulation {
    /// The place on this triangulation of a point of the input mesh (face or edge or
    /// vertex), as `IntrinsicTriangulation::point_on_triangulation`.
    ///
    /// # Errors
    /// [`IntrinsicError::InvalidPoint`] for a point that is not one of the input mesh, and the
    /// representation's refusals.
    fn locate_on_triangulation(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError>;

    /// The place on the input mesh of a point of this triangulation (`T¹ → T⁰`).
    ///
    /// # Errors
    /// As [`locate_on_triangulation`](Self::locate_on_triangulation).
    fn locate_on_input(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError>;
}

impl InputPoints for IntrinsicTriangulation {
    fn locate_on_triangulation(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError> {
        self.point_on_triangulation(point)
    }
    fn locate_on_input(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError> {
        self.point_on_input(point)
    }
}

impl InputPoints for IntegerTriangulation {
    fn locate_on_triangulation(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError> {
        self.point_on_triangulation(point)
    }
    fn locate_on_input(&self, point: SurfacePoint) -> Result<SurfacePoint, IntrinsicError> {
        self.point_on_input(point)
            .map_err(|e| ddg_integer_coords::to_intrinsic_error(e, "point_on_input"))
    }
}
