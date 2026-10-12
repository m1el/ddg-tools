//! The options of the pipeline: one plain struct, every policy explicit.

use std::f64::consts::TAU;

use crate::error::UvError;

/// Which signs cone angles may have.
///
/// A cone of negative angle has more than a full turn around it. Where the cut ends at such a
/// cone the layout overlaps itself (measured in the `ddg-bff` README: 1 % to 15 % of the
/// footprint on nine of sixteen dataset meshes with negative cones, none with positive cones
/// only). On a closed piece the cone angles must sum to `2πχ`, so positive cones alone are
/// impossible for genus 2 and above.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConeSigns {
    /// `µ ≥ 0` on every piece where that is feasible: pieces with boundary, and closed pieces of
    /// genus 0 and 1. On a closed piece of genus 1 it leaves no cone at all (`Σµ = 0`). On a
    /// closed piece of genus 2 or more any sign is allowed, and the report carries a warning.
    /// The default.
    #[default]
    PositiveWhereFeasible,
    /// `µ ≥ 0` everywhere; a closed piece of genus 2 or more is an error
    /// ([`UvError::NegativeConesRequired`]).
    PositiveOnly,
    /// No restriction of the sign.
    Any,
}

/// Where an automatic slit is put on a closed piece.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SlitPlacement {
    /// In the middle of a longest shortest path: from the vertex of largest absolute angle
    /// defect the farthest vertex `a` is found, from `a` the farthest vertex `b` (distances
    /// along polygon edges), and the slit is the middle part of the shortest path from `a` to
    /// `b`. The default.
    #[default]
    DiameterMiddle,
    /// Starting at the vertex of largest absolute angle defect (the site the paper suggests
    /// for a free cone, Sec. 8.1), along the shortest path to the vertex farthest from it.
    AtLargestCurvature,
}

/// What to do with a connected piece of the surface that has no boundary after the seams were
/// cut.
///
/// The variants other than [`ClosedSurfaces::Refuse`] concern closed pieces of **genus 0**.
/// There, with positive cones, the cone cost is the constant `4πλ`, `λ` steers nothing, and on
/// a round or convex shape the optimum is a cone at every vertex. Closed pieces of higher
/// genus are solved as they are (see [`ConeSigns`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClosedSurfaces {
    /// Return [`UvError::ClosedSurface`] for every closed piece, of any genus: the caller has
    /// to mark a seam.
    Refuse,
    /// Solve the closed problem as it stands and use its optimum (with cluster merging as
    /// [`MergeClusters`] says).
    RawOptimum,
    /// One free cone site (paper Sec. 8.1): a vertex where a cone costs nothing. The site is
    /// the vertex of largest absolute angle defect, unless the caller pinned vertices on the
    /// piece (pinned vertices are free sites already; no further site is added then).
    FreeSite {
        /// Upper bound of the cone angle at the site. In the solver's linear model a free site
        /// can take `2π` or more, which no cone can have; such a result is an error
        /// ([`UvError::UnrealizableCone`]). The bound of all other vertices is
        /// [`UnwrapOptions::max_cone_angle`].
        max_angle: f64,
    },
    /// Cut a slit along polygon edges, so that the piece gets a small boundary. The scale factor
    /// is zero along it, the two sides of the slit keep their lengths, and `λ` trades cones
    /// against distortion as on any surface with boundary. The slit is part of the reported
    /// cut.
    Slit {
        /// Where.
        placement: SlitPlacement,
        /// Its length as a fraction of the longest shortest path it is taken from (see
        /// [`SlitPlacement`]), between 0 and 1. The slit always has at least two edges.
        relative_length: f64,
    },
    /// [`ClosedSurfaces::RawOptimum`], and if its result is not usable — more than `max_cones`
    /// cones are left after merging — a slit instead (a second solve).
    Auto {
        /// The largest number of cones accepted from the raw optimum.
        max_cones: usize,
        /// Placement of the fallback slit.
        placement: SlitPlacement,
        /// Length of the fallback slit, see [`ClosedSurfaces::Slit`].
        relative_length: f64,
    },
}

impl Default for ClosedSurfaces {
    fn default() -> Self {
        ClosedSurfaces::Slit {
            placement: SlitPlacement::DiameterMiddle,
            relative_length: 0.5,
        }
    }
}

/// Whether edge-connected clusters of cones are merged into single cones (paper Sec. 6.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MergeClusters {
    /// Use the optimal measure as it is.
    Never,
    /// Merge every cluster of at most `max_cluster_size` cones, whatever it does to the
    /// distortion (which is reported).
    Always {
        /// Clusters with more cones are left alone.
        max_cluster_size: usize,
    },
    /// Merge, then compare the distortion `½ Σ m u²` of the merged cones with that of the
    /// optimum. The merged cones are kept if
    /// `after ≤ max(max_factor · before, before + max_increase)`; otherwise the unmerged
    /// optimum is used and the report carries a warning.
    Guarded {
        /// Clusters with more cones are left alone.
        max_cluster_size: usize,
        /// Largest accepted ratio of the distortion after merging to the distortion before.
        max_factor: f64,
        /// An absolute increase of the distortion that is always accepted (the distortion of
        /// the optimum can be zero, which makes any ratio infinite).
        max_increase: f64,
    },
}

impl Default for MergeClusters {
    fn default() -> Self {
        MergeClusters::Guarded {
            max_cluster_size: usize::MAX,
            max_factor: 2.0,
            max_increase: 0.01,
        }
    }
}

/// What to do when polygons of one connected component are wound inconsistently.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orientation {
    /// Flip the minority of every component (`ddg_core::topology::orient_faces`), report how
    /// many polygons were flipped, and warn. The UVs of a flipped polygon are mirrored with
    /// respect to its own winding. The default.
    #[default]
    Repair,
    /// Return [`UvError::InconsistentOrientation`].
    Reject,
}

/// What to do with triangles of the triangulation that have no area (collinear or coincident
/// corners; Blender's tessellation of polygons with collinear vertices produces them).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Degenerate {
    /// Return [`UvError::DegenerateTriangles`].
    Reject,
    /// If, and only if, some triangle is degenerate: lengthen **all** edge lengths by the
    /// smallest common amount `ε` that makes every triangle satisfy the triangle inequality
    /// with the margin `relative_delta · (mean edge length)`
    /// (`ddg_bff::mollify_edge_lengths`, after Sharp and Crane 2020). A mesh without
    /// degenerate triangles is not touched. `ε`, the margin and the number of triangles that
    /// missed it are reported, with a warning. The default, with the reference
    /// implementation's margin `1e-6`.
    MollifyIfNeeded {
        /// The margin relative to the mean edge length.
        relative_delta: f64,
    },
    /// Mollify always, as the reference implementation of Boundary First Flattening does. This
    /// also changes meshes whose thinnest triangles are valid but miss the margin.
    MollifyAlways {
        /// The margin relative to the mean edge length.
        relative_delta: f64,
    },
}

impl Default for Degenerate {
    fn default() -> Self {
        Degenerate::MollifyIfNeeded {
            relative_delta: ddg_bff::REFERENCE_MOLLIFY_DELTA,
        }
    }
}

/// Options of the chart packer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackOptions {
    /// Free space around every chart, in units of the unit square: charts keep `margin` from
    /// the border and `2 · margin` from each other (their bounding boxes do). With many charts
    /// it is capped at `1 / (4 ⌈√n⌉)`; the report says which margin was used.
    pub margin: f64,
    /// Whether a chart may be rotated to the orientation in which its bounding box is
    /// smallest. Rotation is a rigid motion; it does not change the map's distortion.
    pub rotate: bool,
}

impl Default for PackOptions {
    fn default() -> Self {
        PackOptions {
            margin: 1.0 / 256.0,
            rotate: true,
        }
    }
}

/// The options of [`crate::unwrap`]. `UnwrapOptions::default()` is the configuration the crate
/// is tested with; the crate README gives the measurements behind every default.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnwrapOptions {
    /// The trade-off `λ ≥ 0` between area distortion and total cone angle, in the cone
    /// solver's convention: the area of **every piece** (connected part of the surface after
    /// the seams were cut) is normalized to 1. Larger means fewer cones.
    pub lambda: f64,
    /// Which signs cone angles may have.
    pub cone_signs: ConeSigns,
    /// Upper bound of every cone angle, at most `2π`. A cone of `2π` or more does not exist;
    /// the solver's linear model can ask for one (at free sites, or when clusters are merged).
    pub max_cone_angle: f64,
    /// What to do with closed pieces.
    pub closed_surfaces: ClosedSurfaces,
    /// Whether clusters of cones are merged.
    pub merge: MergeClusters,
    /// What to do with inconsistently wound polygons.
    pub orientation: Orientation,
    /// What to do with triangles without area.
    pub degenerate: Degenerate,
    /// Packing of the charts into the unit square, or `None` to leave every chart where the
    /// flattener put it (first boundary vertex at the origin, at the scale `e^u`; charts then
    /// lie on top of each other).
    pub packing: Option<PackOptions>,
    /// Resolution of the overlap estimate of every chart: about this many sample points per
    /// triangle of average size. `0.0` switches the estimate off.
    pub overlap_samples_per_triangle: f64,
    /// Whether cones that the solver could not certify as optimal are an error
    /// ([`UvError::ConesNotCertified`]). With `false` they are used and the report warns.
    pub require_certified_cones: bool,
}

impl Default for UnwrapOptions {
    fn default() -> Self {
        UnwrapOptions {
            lambda: 0.5,
            cone_signs: ConeSigns::default(),
            max_cone_angle: 0.95 * TAU,
            closed_surfaces: ClosedSurfaces::default(),
            merge: MergeClusters::default(),
            orientation: Orientation::default(),
            degenerate: Degenerate::default(),
            packing: Some(PackOptions::default()),
            overlap_samples_per_triangle: ddg_bff::DEFAULT_OVERLAP_SAMPLES,
            require_certified_cones: false,
        }
    }
}

fn invalid(name: &'static str, value: impl ToString, requirement: &'static str) -> UvError {
    UvError::InvalidOption {
        name,
        value: value.to_string(),
        requirement,
    }
}

impl UnwrapOptions {
    /// Checks every option.
    ///
    /// # Errors
    /// [`UvError::InvalidOption`] naming the option, its value and the requirement.
    pub fn validate(&self) -> Result<(), UvError> {
        if !(self.lambda.is_finite() && self.lambda >= 0.0) {
            return Err(invalid(
                "lambda",
                self.lambda,
                "must be a finite number that is not negative",
            ));
        }
        if !(self.max_cone_angle > 0.0 && self.max_cone_angle <= TAU) {
            return Err(invalid(
                "max_cone_angle",
                self.max_cone_angle,
                "must be positive and at most 2π",
            ));
        }
        let check_slit = |length: f64| {
            if length.is_finite() && (0.0..=1.0).contains(&length) {
                Ok(())
            } else {
                Err(invalid(
                    "closed_surfaces.relative_length",
                    length,
                    "must be between 0 and 1",
                ))
            }
        };
        match self.closed_surfaces {
            ClosedSurfaces::Refuse | ClosedSurfaces::RawOptimum => {}
            ClosedSurfaces::FreeSite { max_angle } => {
                // Infinity is allowed: it is the paper's free site as it stands.
                if max_angle.is_nan() || max_angle <= 0.0 {
                    return Err(invalid(
                        "closed_surfaces.max_angle",
                        max_angle,
                        "must be positive",
                    ));
                }
            }
            ClosedSurfaces::Slit {
                relative_length, ..
            } => check_slit(relative_length)?,
            ClosedSurfaces::Auto {
                relative_length, ..
            } => check_slit(relative_length)?,
        }
        match self.merge {
            MergeClusters::Never | MergeClusters::Always { .. } => {}
            MergeClusters::Guarded {
                max_factor,
                max_increase,
                ..
            } => {
                if max_factor.is_nan() || max_factor < 1.0 {
                    return Err(invalid(
                        "merge.max_factor",
                        max_factor,
                        "must be at least 1",
                    ));
                }
                if max_increase.is_nan() || max_increase < 0.0 {
                    return Err(invalid(
                        "merge.max_increase",
                        max_increase,
                        "must not be negative",
                    ));
                }
            }
        }
        match self.degenerate {
            Degenerate::Reject => {}
            Degenerate::MollifyIfNeeded { relative_delta }
            | Degenerate::MollifyAlways { relative_delta } => {
                if !(relative_delta.is_finite() && relative_delta > 0.0) {
                    return Err(invalid(
                        "degenerate.relative_delta",
                        relative_delta,
                        "must be positive and finite",
                    ));
                }
            }
        }
        if let Some(pack) = self.packing
            && !(pack.margin.is_finite() && (0.0..0.5).contains(&pack.margin))
        {
            return Err(invalid(
                "packing.margin",
                pack.margin,
                "must be at least 0 and below 0.5",
            ));
        }
        if !(self.overlap_samples_per_triangle.is_finite()
            && self.overlap_samples_per_triangle >= 0.0)
        {
            return Err(invalid(
                "overlap_samples_per_triangle",
                self.overlap_samples_per_triangle,
                "must be finite and not negative",
            ));
        }
        Ok(())
    }
}
