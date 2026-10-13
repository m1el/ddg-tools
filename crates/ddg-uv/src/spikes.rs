//! The spike check: do a handful of vertices carry almost all of the distortion that the
//! surface would have without cones?
//!
//! # Why
//!
//! A vertex pulled far out of the surface (a thin needle) has a very large log scale factor in
//! the linearized model of the cone paper: for a spike of height `h` over a base of radius `ρ`
//! the cone-free `u` at the tip is of order `−h/ρ` on a mass of order `ρ h`, so the vertex alone
//! adds about `h³/ρ` to `½ Σ m u²`. The optimum then spends its cones on the spikes and the map
//! is poor (`ddg-optimal-cones`, README, "Needle-like spikes"; on `camelshotnoise` the map is
//! unusable). The paper ignores such outliers; the model does not.
//!
//! # The criterion
//!
//! For every piece (connected part of the surface after the seams were cut) the pipeline
//! evaluates the state **without cones** with [`ddg_optimal_cones::evaluate_cones`] (only the
//! cones the caller prescribed are kept, since the solver cannot change them). That gives the log
//! scale `u_i` and the distortion `½ Σ m_i u_i²` (masses normalized to sum 1). With
//! `c_i = m_i u_i²`, sorted in decreasing order, let `k` be the smallest number such that the `k`
//! largest `c_i` add up to at least a share `s` of `Σ c_i`. The piece has spikes if
//!
//! 1. `Σ c_i` is at least `min_distortion` (below it the `u_i` are rounding noise),
//! 2. `k ≤ max_vertices`, and
//! 3. `k ≤ max_fraction · (number of vertices of the piece)`.
//!
//! The spikes are those `k` vertices. The defaults (`s = 0.99`, 64 vertices, 1 %) are in
//! [`SpikeCheck`]; the crate README has the measurements behind them.
//!
//! The cost is one factorization of the cotan Laplacian and two solves per piece, much less
//! than the cone solve (see the README).

use ddg_core::io::Cone;
use ddg_core::operators::lumped_mass_normalized;
use ddg_optimal_cones::{ConeProblem, Parallelism, evaluate_cones};

use crate::options::SpikeCheck;
use crate::prepare::{Piece, Prepared};

/// What the spike check found on one piece.
#[derive(Debug, Clone, PartialEq)]
pub struct SpikeReport {
    /// The distortion `½ Σ m u²` of the piece without cones (only prescribed cones), with the
    /// area of the piece normalized to 1.
    pub cone_free_distortion: f64,
    /// The smallest number of vertices that carry the configured share of that distortion.
    pub n_carrying: usize,
    /// The spikes, as vertices of the caller's mesh, ascending: empty unless the criterion
    /// holds.
    pub vertices: Vec<u32>,
    /// The share of the cone-free distortion that the spikes carry (0 if there are none).
    pub share: f64,
    /// The spike vertices that were given zero distortion weight in the cone problem, as
    /// vertices of the caller's mesh, ascending. Empty if the action is warn only, or the
    /// piece has no spikes.
    pub excluded: Vec<u32>,
    /// Spike vertices that were not excluded because they lie on the boundary of the piece.
    pub not_excluded_boundary: Vec<u32>,
}

/// Runs the check on one piece (cone-free state with `u = 0` on the boundary).
///
/// # Errors
/// The message of the error of the Laplacian solve; the caller reports it as a warning and
/// carries on without the check.
pub(crate) fn check_piece(
    prepared: &Prepared,
    piece: &Piece,
    check: &SpikeCheck,
) -> Result<Option<SpikeReport>, String> {
    let mesh = &piece.mesh;
    let n = mesh.n_vertices();
    // Prescribed cones at interior vertices are fixed whatever the solver does.
    let prescribed: Vec<Cone> = (0..n)
        .filter(|&v| !mesh.is_boundary_vertex(v))
        .filter(|&v| prepared.steering_role[piece.vertex_input[v] as usize] == 2)
        .map(|v| Cone {
            vertex: v as u32,
            angle: prepared.prescribed_angle[piece.vertex_input[v] as usize],
        })
        .collect();
    if (0..n).all(|v| mesh.is_boundary_vertex(v)) {
        return Ok(None);
    }
    let problem = ConeProblem::new(mesh, &piece.geometry, 0.0);
    // Always with `u = 0` on the boundary, also when the pipeline optimizes the boundary
    // values: the cone-free state with free boundary spreads the curvature of a needle over
    // the whole piece (measured: on a disk with two needles 317 of 625 vertices carry 99 %
    // of it), whereas the local, fixed-boundary state is what identifies the needle.
    let log_scale = evaluate_cones(&problem, &prescribed, Parallelism::Sequential)
        .map_err(|e| e.to_string())?
        .log_scale;
    let mass = lumped_mass_normalized(mesh, &piece.geometry).map_err(|e| e.to_string())?;
    let mut parts: Vec<(f64, usize)> = (0..n)
        .map(|v| (mass[v] * log_scale[v] * log_scale[v], v))
        .collect();
    let total: f64 = parts.iter().map(|p| p.0).sum();
    if !total.is_finite() {
        return Err("the cone-free scale factors are not finite".to_owned());
    }
    let mut report = SpikeReport {
        cone_free_distortion: 0.5 * total,
        n_carrying: 0,
        vertices: Vec::new(),
        share: 0.0,
        excluded: Vec::new(),
        not_excluded_boundary: Vec::new(),
    };
    if total <= 0.0 {
        return Ok(Some(report));
    }
    // Largest first; ties by vertex, so that the result does not depend on the sort.
    parts.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut carried = 0.0;
    let mut k = 0;
    while k < n && carried < check.share * total {
        carried += parts[k].0;
        k += 1;
    }
    report.n_carrying = k;
    if 0.5 * total >= check.min_distortion
        && k <= check.max_vertices
        && (k as f64) <= check.max_fraction * n as f64
    {
        report.share = carried / total;
        report.vertices = parts[..k].iter().map(|p| piece.vertex_input[p.1]).collect();
        report.vertices.sort_unstable();
    }
    Ok(Some(report))
}

/// The warning for a piece with spikes.
pub(crate) fn warning(piece: &Piece, spikes: &SpikeReport) -> String {
    let handled = if spikes.excluded.is_empty() {
        String::new()
    } else {
        format!(
            ". {} of them were excluded from the distortion measure of the cone problem \
             (zero distortion weight){}",
            spikes.excluded.len(),
            if spikes.not_excluded_boundary.is_empty() {
                String::new()
            } else {
                format!(
                    "; {} on the boundary were not",
                    spikes.not_excluded_boundary.len()
                )
            }
        )
    };
    const LISTED: usize = 16;
    let mut list: Vec<String> = spikes
        .vertices
        .iter()
        .take(LISTED)
        .map(u32::to_string)
        .collect();
    if spikes.vertices.len() > LISTED {
        list.push(format!("and {} more", spikes.vertices.len() - LISTED));
    }
    format!(
        "{} of the {} vertices of the piece of the surface that contains face {} carry {:.1} % \
         of the distortion it would have without cones ({:.3e}): vertices {}. This is what a thin \
         spike (a vertex pulled far out of the surface) looks like to this method; the spikes \
         dominate the distortion measure and the map is likely to be poor, with many cones at \
         the largest allowed angle and a UV map that overlaps itself. Repair the mesh first: \
         move these vertices back onto the surface or remove them{}",
        spikes.vertices.len(),
        piece.mesh.n_vertices(),
        piece.first_polygon,
        100.0 * spikes.share,
        spikes.cone_free_distortion,
        list.join(", "),
        handled
    )
}
