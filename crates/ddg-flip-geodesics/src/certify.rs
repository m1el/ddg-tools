//! The geodesic certificate on the triangulation (SPEC §7 layer 1a): computed on any state of
//! a network, never assumed.

use std::f64::consts::PI;

use ddg_intrinsic::Triangulation;

use crate::error::FlipError;
use crate::network::{Joint, Network};
use crate::wedge::{JointGeometry, Side};

/// Options of [`certify`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifyOptions {
    /// A joint is straight when `α^min ≥ π − tolerance` (radians). Default `1e-4`, the
    /// paper's criterion (Sec. 5.3; F3).
    pub tolerance: f64,
}

impl Default for CertifyOptions {
    fn default() -> Self {
        CertifyOptions { tolerance: 1e-4 }
    }
}

/// How the certificate classified a joint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum JointStatus {
    /// At a pinned vertex: exempt.
    Pinned,
    /// Both wedge angles are at least `π − tolerance`.
    Straight,
    /// Every side below `π − tolerance` is blocked by a path end: locally shortest in its
    /// isotopy class (Sec. 4.2, case (b)).
    BlockedByEnd,
    /// Some side below `π − tolerance` is not blocked by a path end.
    NotStraight,
}

/// A joint that fails the certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JointFailure {
    /// The joint.
    pub joint: Joint,
    /// Its vertex.
    pub vertex: usize,
    /// The side below `π − tolerance` that is not blocked by a path end (the smaller one if
    /// both are).
    pub side: Side,
    /// Its angle.
    pub angle: f64,
    /// Whether that wedge is flexible (FlipOut could shorten it right away); if not, it is
    /// blocked by another part of the network, which is then not geodesic either.
    pub flexible: bool,
}

/// The result of [`certify`].
#[derive(Debug, Clone, PartialEq)]
pub struct GeodesicCertificate {
    /// Whether every joint is pinned, straight or blocked by an end.
    pub geodesic: bool,
    /// The tolerance used.
    pub tolerance: f64,
    /// Number of joints.
    pub n_joints: usize,
    /// Joints at pinned vertices (exempt).
    pub n_pinned: usize,
    /// Straight joints.
    pub n_straight: usize,
    /// Joints that are not straight but blocked by a path end (exempt).
    pub n_blocked: usize,
    /// The smallest `α^min` over the joints that are neither pinned nor blocked; `+∞` if
    /// there are none.
    pub min_angle: f64,
    /// The joints that fail, in path order.
    pub failures: Vec<JointFailure>,
}

/// Certifies a network on `tri` (SPEC §7 layer 1a): the network is valid
/// ([`Network::check`]), and every joint not at a pinned vertex is straight to the tolerance
/// (both wedge angles `≥ π − tolerance`) or every side below that is blocked by a path end.
///
/// Stricter than "the small wedge is blocked" of SPEC §1.1: at a vertex with `Θ_b < 2π` both
/// sides can be below `π`, and an unblocked large side is a shortening the algorithm makes.
///
/// # Errors
/// [`FlipError::InvalidOptions`] for a tolerance that is negative or not finite; the errors
/// of [`Network::check`] for an invalid network.
pub fn certify<T: Triangulation>(
    tri: &T,
    network: &Network,
    options: &CertifyOptions,
) -> Result<GeodesicCertificate, FlipError> {
    let tol = options.tolerance;
    if !(tol.is_finite() && tol >= 0.0) {
        return Err(FlipError::InvalidOptions {
            what: format!("the certificate tolerance must be finite and >= 0, got {tol}"),
        });
    }
    network.check(tri.complex())?;
    let mut certificate = GeodesicCertificate {
        geodesic: true,
        tolerance: tol,
        n_joints: 0,
        n_pinned: 0,
        n_straight: 0,
        n_blocked: 0,
        min_angle: f64::INFINITY,
        failures: Vec::new(),
    };
    for joint in network.joints() {
        certificate.n_joints += 1;
        let geometry = JointGeometry::new(tri, network, joint)?;
        match classify(&geometry, network, tol) {
            (JointStatus::Pinned, _) => certificate.n_pinned += 1,
            (JointStatus::Straight, _) => {
                certificate.n_straight += 1;
                certificate.min_angle = certificate.min_angle.min(geometry.min_angle());
            }
            (JointStatus::BlockedByEnd, _) => certificate.n_blocked += 1,
            (JointStatus::NotStraight, failing) => {
                certificate.min_angle = certificate.min_angle.min(geometry.min_angle());
                let side = failing.expect("a failing side");
                certificate.geodesic = false;
                certificate.failures.push(JointFailure {
                    joint,
                    vertex: geometry.vertex(),
                    side,
                    angle: geometry.angle(side),
                    flexible: geometry.is_flexible(side),
                });
            }
        }
    }
    Ok(certificate)
}

/// The status of one joint for the certificate, and the first failing side (smaller first).
pub fn classify(
    geometry: &JointGeometry,
    network: &Network,
    tolerance: f64,
) -> (JointStatus, Option<Side>) {
    if network.is_pinned(geometry.vertex()) {
        return (JointStatus::Pinned, None);
    }
    let min = geometry.min_side();
    let bent: Vec<Side> = [min, min.opposite()]
        .into_iter()
        .filter(|&s| geometry.angle(s) < PI - tolerance)
        .collect();
    if bent.is_empty() {
        return (JointStatus::Straight, None);
    }
    match bent
        .iter()
        .copied()
        .find(|&s| !geometry.is_blocked_by_end(network, s))
    {
        Some(side) => (JointStatus::NotStraight, Some(side)),
        None => (JointStatus::BlockedByEnd, None),
    }
}
