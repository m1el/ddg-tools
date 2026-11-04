//! Baselines of Sec. 5.2 (SPEC.md §9): point areas, winding numbers (direct, fast, regularized),
//! signed Hopf–Cole (SHC, Eq. 28), smoothed signed planar distance (SSPD, Eq. 29) and the
//! LogSumExp distance.
//!
//! * [`areas`]: Barill point areas (Barill, Dickson, Schmidt, Levin, Jacobson, *Fast Winding
//!   Numbers for Soups and Clouds*, ACM TOG 37(4), 2018, Sec. 3.1): the Voronoi cell of each
//!   point among its `k = 20` nearest neighbors projected to its tangent plane. Written from the
//!   paper's text; libigl's `copyleft/cgal/point_areas` (GPL) was not consulted.
//! * [`winding`]: the generalized winding number of a point cloud by direct summation of area
//!   weighted dipoles (Barill Eq. 5), the regularized winding number (Chen et al. 2024), and the
//!   winding-number-signed distance.
//! * [`fwn`]: Barill's hierarchical approximation (Algorithm 1, the order-0/1/2 Taylor expansion
//!   of Eq. 13–15 and 18–20) for point clouds (octree) and triangle soups (AABB tree).
//! * [`kernel`]: SHC, SSPD and LSE on the neighbor sets and λ rules of the PAT blend
//!   (the reference's `signed_log_sum_exp`, `self_normalized_signed_log_sum_exp` and
//!   `log_sum_exp`, points-as-tori, MIT, Copyright (c) 2026 Nicole Feng).
//!
//! Sign convention everywhere: the winding number is 1 inside a closed, outward oriented surface,
//! signed distances are negative inside. Every batch evaluation runs on an explicit pool and every
//! query is computed independently in a fixed order, so results are bitwise independent of the
//! thread count.

pub mod areas;
pub mod fwn;
pub mod kernel;
pub mod winding;

pub use areas::{PointAreaOptions, PointAreaReport, point_areas, voronoi_cell_area};
pub use fwn::{FastWinding, FwnOptions, FwnReport};
pub use kernel::{KernelField, KernelMethod, KernelOptions, KernelReport, KernelValue};
pub use winding::{
    WindingReport, regularization_s, regularized_kernel, regularized_winding_number, regularized_winding_numbers, winding_number,
    winding_numbers, winding_signed_distance,
};

#[inline]
pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
pub(crate) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// The first non-finite query, as an error.
pub(crate) fn check_queries(xs: &[[f64; 3]]) -> Result<(), crate::PatError> {
    match xs.iter().position(|x| !x.iter().all(|c| c.is_finite())) {
        Some(i) => Err(crate::PatError::NonFinite { index: i, what: "query" }),
        None => Ok(()),
    }
}

/// Check an output buffer length.
pub(crate) fn check_len(name: &str, len: usize, q: usize) -> Result<(), crate::PatError> {
    if len != q {
        return Err(crate::PatError::SizeMismatch { what: format!("{name} buffer has {len} entries for {q} queries") });
    }
    Ok(())
}

/// Queries per parallel task.
pub(crate) const CHUNK: usize = 256;
