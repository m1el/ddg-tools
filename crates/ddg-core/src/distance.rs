//! Upwind (Hopf–Lax) distance on a triangle mesh: a first-order solver for the eikonal equation
//! `|∇T| = 1` with prescribed values at seed vertices, computed from the edge lengths alone.
//!
//! The scheme defines the value `T_i` of a non-seed vertex as the minimum, over the triangles
//! `(i, j, k)` at `i`, of the cheapest way to reach `i` from the opposite edge `jk`, where the
//! value on the edge is the linear interpolant of `T_j` and `T_k`:
//!
//! ```text
//! T_i = min over triangles (i, j, k) of  min_{λ ∈ [0,1]}  λ T_j + (1−λ) T_k + d_i(λ)
//! ```
//!
//! with `d_i(λ)` the distance inside the flat triangle from `i` to the point
//! `λ x_j + (1−λ) x_k` (see [`triangle_update`] for the closed form of the inner minimum). If
//! only one of `j`, `k` is *usable*, the edge term `T_j + ℓ_ij` is taken instead. A neighbour is
//! usable when its value is finite (or, in the single pass, when it is accepted) and, if a
//! per-vertex class array is given, its class is the class of `i`: the fronts of different
//! classes never see each other. This is the discrete Hopf–Lax formula of the local variational
//! principle of Bornemann and Rasch; Kimmel and Sethian use the same update inside fast
//! marching (and unfold obtuse triangles, which is not done here).
//!
//! # Two modes ([`DistanceMode`])
//!
//! * **Single pass (the default).** Fast marching: seeds are accepted at the start; the vertex
//!   with the smallest tentative value (ties: lower vertex index) is accepted next and is final;
//!   a tentative value is the minimum of the update over triangles whose other two vertices are
//!   both accepted, and of the edge terms over accepted neighbours. `O(n log n)` whatever the
//!   angles, deterministic, and independent of the vertex labeling when there are no ties.
//! * **Fixed point.** Iterates `T_i = min over all incident triangles` to its fixed point with
//!   a Dijkstra-ordered pass that re-opens a vertex whenever its value can still strictly
//!   decrease. Independent of the processing order, always at most the single-pass value, but
//!   expensive on obtuse triangles (below).
//!
//! The two agree (to rounding) on meshes without obtuse angles, where a vertex never needs a
//! neighbour with a larger value. On obtuse triangles the single pass is larger: on sheared
//! unit grids with three seeds, by a mean of 0.2% (largest angle 117°), 2% (135°), 4% (146°)
//! and 6% (153°) of the value, with maxima 2%, 11%, 19% and 25% (`tests/distance.rs`). It is
//! never below the fixed point (by induction over the order of acceptance: each accepted value
//! is the update of some triangles at values at least the fixed-point values, and the update
//! is monotone).
//!
//! # Cost of the fixed-point mode on obtuse triangles
//!
//! On a mesh of non-obtuse triangles the Dijkstra order settles every vertex once. With obtuse
//! triangles the scheme is a stochastic-shortest-path iteration (each update is a convex
//! combination of two neighbour values plus a positive length), and the relaxation after a
//! re-opening converges geometrically, with a rate that tends to 1 as the obtuse angle tends
//! to π. Measured on sheared grids: about 0.1 to 0.3 re-openings per vertex at a largest angle
//! of 117°, 1 to 3 at 135°, 100 to 400 at 153°, and thousands at 162°.
//! [`DistanceMode::FixedPoint`] bounds the work and [`FixedPointReport::converged`] reports
//! whether the bound was hit.
//!
//! With the first-order scheme the error against the exact distance on a smooth surface is
//! `O(h)` (a front that crosses an edge sees the linear interpolant of the values at its ends).
//!
//! # References
//!
//! * R. Kimmel, J. A. Sethian, *Computing geodesic paths on manifolds*, Proceedings of the
//!   National Academy of Sciences 95(15), 8431–8435, 1998, <https://doi.org/10.1073/pnas.95.15.8431>.
//! * F. Bornemann, C. Rasch, *Finite-element discretization of static Hamilton–Jacobi equations
//!   based on a local variational principle*, Computing and Visualization in Science 9(2),
//!   57–69, 2006 (they solve the system by an adaptive Gauss–Seidel iteration), <https://doi.org/10.1007/s00791-006-0016-y>
//!   (arXiv:math/0403517; DOI, volume, issue and pages confirmed against the Crossref record).
//!
//! # Example
//!
//! ```
//! use ddg_core::distance::upwind_distance;
//! use ddg_core::generators::grid_disk;
//! use ddg_core::Geometry;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let grid = grid_disk(4, 4)?;
//! let mesh = grid.topology()?;
//! let geometry = Geometry::from_positions(&mesh, &grid.positions)?;
//! let result = upwind_distance(&mesh, &geometry, &[(0, 0.0)], None)?;
//! assert_eq!(result.n_unreached, 0);
//! assert!((result.values[5] - 0.25).abs() < 1e-14); // the neighbour along an edge
//! # Ok(())
//! # }
//! ```

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fmt;

use crate::geometry::{Geometry, GeometryError, corner_cotan, triangle_area};
use crate::topology::TriMesh;

/// Why a distance computation was refused. Nothing is computed in these cases.
#[derive(Debug, Clone, PartialEq)]
pub enum DistanceError {
    /// The geometry does not belong to the mesh.
    Geometry(GeometryError),
    /// A seed names a vertex that the mesh does not have.
    SeedOutOfRange {
        /// Position of the seed in the seed list.
        seed: usize,
        /// The vertex index it names.
        vertex: u32,
        /// Number of vertices of the mesh.
        n_vertices: usize,
    },
    /// Two seeds name the same vertex.
    DuplicateSeed {
        /// The vertex.
        vertex: u32,
        /// Position of the first seed naming it.
        first: usize,
        /// Position of the second seed naming it.
        second: usize,
    },
    /// A seed value is negative, NaN or infinite.
    InvalidSeedValue {
        /// Position of the seed in the seed list.
        seed: usize,
        /// The seed vertex.
        vertex: u32,
        /// The offending value.
        value: f64,
    },
    /// The class array does not have one entry per vertex.
    ClassCountMismatch {
        /// Number of vertices of the mesh.
        expected: usize,
        /// Length of the class array.
        got: usize,
    },
}

impl fmt::Display for DistanceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DistanceError::Geometry(e) => write!(f, "{e}"),
            DistanceError::SeedOutOfRange {
                seed,
                vertex,
                n_vertices,
            } => write!(
                f,
                "seed {seed} names vertex {vertex}, but the mesh has only {n_vertices} vertices"
            ),
            DistanceError::DuplicateSeed {
                vertex,
                first,
                second,
            } => write!(f, "seeds {first} and {second} both name vertex {vertex}"),
            DistanceError::InvalidSeedValue {
                seed,
                vertex,
                value,
            } => write!(
                f,
                "seed {seed} at vertex {vertex} has value {value}, which is negative or not finite"
            ),
            DistanceError::ClassCountMismatch { expected, got } => write!(
                f,
                "the mesh has {expected} vertices but {got} vertex classes were given"
            ),
        }
    }
}

impl std::error::Error for DistanceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            DistanceError::Geometry(e) => Some(e),
            _ => None,
        }
    }
}

impl From<GeometryError> for DistanceError {
    fn from(e: GeometryError) -> Self {
        DistanceError::Geometry(e)
    }
}

/// The solution of the upwind scheme and how it was obtained.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceResult {
    /// The value `T_i` of every vertex; seeds keep their seed value. `f64::INFINITY` where no
    /// front arrives (no seed of the vertex's class in its connected component).
    pub values: Vec<f64>,
    /// Number of vertices whose value is infinite.
    pub n_unreached: usize,
    /// Statistics of the iteration; `Some` only in [`DistanceMode::FixedPoint`], `None` for the
    /// single pass (which has nothing to converge and never re-opens).
    pub fixed_point: Option<FixedPointReport>,
}

/// What the fixed-point iteration did; see [`DistanceMode::FixedPoint`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FixedPointReport {
    /// How many times a vertex that had already been settled received a strictly smaller value
    /// and was put back on the heap. Zero on meshes where the Dijkstra order is already
    /// consistent (no obtuse angles); positive when obtuse triangles let a later, larger value
    /// improve an earlier one.
    pub n_reopened: usize,
    /// Whether the heap ran empty. `false` means `max_reopenings` was reached; the values are
    /// then still upper bounds of the fixed point (every value is the update of one triangle at
    /// some earlier state), but not the fixed point.
    pub converged: bool,
}

/// The local update of the upwind scheme at vertex `i` of one triangle `(i, j, k)`: the minimum
/// over `λ ∈ [0, 1]` of `λ T_j + (1−λ) T_k + d_i(λ)`, where `d_i(λ)` is the distance from `i`
/// to `λ x_j + (1−λ) x_k` in the flat triangle with side lengths `l_ij`, `l_ik`, `l_jk`.
///
/// # Closed form
///
/// Write `c = ℓ_jk`, let `h` be the height of `i` over the line `jk`, and let `u` be the
/// coordinate along that line measured from the foot of the height. Then `j` is at
/// `u = −x₀` and `k` at `u = c − x₀`, with `x₀ = h cot θ_j` and `c − x₀ = h cot θ_k`, where
/// `θ_j`, `θ_k` are the triangle's angles at `j` and `k` (a cotangent is negative for an
/// obtuse angle, which puts the foot outside the edge). The point at fraction `t = 1 − λ` from
/// `j` towards `k` has `u = t c − x₀`. With the slope `g = (T_k − T_j)/c` the objective is
///
/// ```text
/// f(u) = T_j + g x₀ + g u + √(u² + h²).
/// ```
///
/// It is convex in `u`. For `|g| ≥ 1` it is monotone and the minimum is at an end of the
/// interval. For `|g| < 1`, `f′(u) = g + u/√(u²+h²) = 0` gives `u* = −g h/√(1−g²)` and
/// `√(u*²+h²) = h/√(1−g²)`, hence
///
/// ```text
/// f(u*) = T_j + g x₀ + h √(1−g²) = T_j + h ( g cot θ_j + √(1−g²) ).
/// ```
///
/// By convexity the constrained minimum is `f(u*)` if `u*` lies in the interval and the value at
/// the nearer end otherwise. The test `u* > −h cot θ_j`, in terms of
/// `q = √(1−g²)`, reads `g < cot θ_j · q` (otherwise the minimum is `T_j + ℓ_ij`), and
/// `u* < h cot θ_k` reads `−g < cot θ_k · q` (otherwise the minimum is `T_k + ℓ_ik`). No
/// division and no angle appears, `1 − g²` is evaluated as `(1−g)(1+g)`, and `T_j = T_k`
/// (`g = 0`) gives `T_j + h` when the foot lies on the edge, `min(T_j + ℓ_ij, T_k + ℓ_ik)`
/// otherwise. The cotangents and `h = 2A/c` come from the stable formulas of
/// [`crate::geometry`], so needle and obtuse triangles lose no accuracy. The result is capped by
/// `min(T_j + ℓ_ij, T_k + ℓ_ik)` so that it never exceeds the two edge terms, whatever the
/// rounding.
///
/// Returns `None` if the three lengths do not form a triangle of positive area or if `T_j` or
/// `T_k` is not finite. The update is nondecreasing in `T_j` and `T_k` (it is a minimum of
/// increasing functions), which the solver relies on.
///
/// [formula: cut-redistance]
pub fn triangle_update(t_j: f64, t_k: f64, l_ij: f64, l_ik: f64, l_jk: f64) -> Option<f64> {
    if !(t_j.is_finite() && t_k.is_finite()) {
        return None;
    }
    // Corner j is opposite the side ik, corner k opposite the side ij.
    let area = triangle_area(l_ij, l_ik, l_jk)?;
    let cot_j = corner_cotan(l_ik, l_ij, l_jk)?;
    let cot_k = corner_cotan(l_ij, l_ik, l_jk)?;
    Some(update_core(
        t_j,
        t_k,
        [l_ij, l_ik, l_jk],
        [cot_j, cot_k],
        2.0 * area / l_jk,
    ))
}

/// The closed form of [`triangle_update`] with `cot θ_j`, `cot θ_k` and the height `h` already
/// known. The sides are `[ℓ_ij, ℓ_ik, ℓ_jk]`, the cotangents `[cot θ_j, cot θ_k]`.
#[inline]
fn update_core(tj: f64, tk: f64, [a, b, c]: [f64; 3], [cot_j, cot_k]: [f64; 2], h: f64) -> f64 {
    let ends = (tj + a).min(tk + b);
    let g = (tk - tj) / c;
    if g >= 1.0 || g <= -1.0 {
        return ends;
    }
    let q = ((1.0 - g) * (1.0 + g)).sqrt();
    if g >= cot_j * q || -g >= cot_k * q {
        return ends;
    }
    (tj + h * (g * cot_j + q)).min(ends)
}

/// Which scheme [`upwind_distance_with`] solves.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum DistanceMode {
    /// Fast marching with the Hopf–Lax update, one pass (the default). Exactly:
    ///
    /// * all seeds are accepted from the start and never change;
    /// * the tentative value of a non-seed vertex `i` is the minimum of [`triangle_update`]
    ///   over the incident triangles `(i, j, k)` whose other two vertices are both accepted and
    ///   of `i`'s class, and of `T_j + ℓ_ij` over the accepted neighbours `j` of `i`'s class;
    /// * repeatedly the not yet accepted vertex with the smallest tentative value is accepted,
    ///   ties going to the lower vertex index; an accepted value is final.
    ///
    /// Cost: `O(n log n)`, no re-openings, whatever the angles. On a mesh without obtuse angles
    /// it equals [`DistanceMode::FixedPoint`]; on obtuse triangles it can be larger (it never
    /// uses a neighbour with a larger value than the vertex being accepted).
    #[default]
    SinglePass,
    /// Iterate the scheme `T_i = min over all incident triangles` to its fixed point: a
    /// Dijkstra-ordered pass that re-opens a vertex whenever its value can still strictly
    /// decrease. Order independent, but on obtuse triangles the number of re-openings grows
    /// without bound as the angle tends to π (see the module documentation).
    FixedPoint {
        /// Stop after this many re-openings (see [`FixedPointReport::converged`]).
        max_reopenings: usize,
        /// A settled vertex is re-opened only if its value would drop by more than this
        /// fraction of itself. `0.0` re-opens for every strict decrease; a small positive value
        /// (say `1e-14`) stops the solver from chasing differences in the last bits between
        /// triangles that agree mathematically. Must be finite and in `[0, 1)`; otherwise
        /// `0.0` is used.
        reopen_tolerance: f64,
    },
}

impl DistanceMode {
    /// The fixed-point mode with the budget `1000 n + 10 000` re-openings and strict decrease.
    pub fn fixed_point(n_vertices: usize) -> Self {
        DistanceMode::FixedPoint {
            max_reopenings: n_vertices.saturating_mul(1000).saturating_add(10_000),
            reopen_tolerance: 0.0,
        }
    }
}

/// Heap key: values are non-negative finite, so the bit pattern orders like the number.
type Entry = Reverse<(u64, u32)>;

/// Solves the upwind scheme of the module documentation with [`DistanceMode::SinglePass`].
///
/// `seeds` lists `(vertex, value)` pairs: seed values are kept as they are (they may differ, as
/// when a front starts from a curve with known distances). `class`, if given, has one entry per
/// vertex; a vertex only ever takes information from neighbours of its own class.
///
/// # Errors
/// [`DistanceError`] for a geometry of another mesh, a seed out of range, duplicated, negative
/// or non-finite, and a class array of the wrong length.
pub fn upwind_distance(
    mesh: &TriMesh,
    geometry: &Geometry,
    seeds: &[(u32, f64)],
    class: Option<&[u32]>,
) -> Result<DistanceResult, DistanceError> {
    upwind_distance_with(mesh, geometry, seeds, class, &DistanceMode::SinglePass)
}

/// [`upwind_distance`] with an explicit [`DistanceMode`].
///
/// **Fixed-point mode.** A binary heap releases the vertex with the smallest value, and the two
/// other vertices of each incident triangle are updated from it with [`triangle_update`] (or
/// `T_j + ℓ_ij` if the third vertex is not usable, i.e. infinite or of another class). A vertex
/// whose value can still strictly decrease is *re-opened*: it goes back on the heap. On
/// termination every non-seed value is the minimum of the updates of its triangles evaluated
/// with the final values, so the result is the fixed point of the scheme reached from above
/// (it is unique: every update exceeds the smallest value it uses by a positive amount, so a
/// fixed point has no cycles) and relabelling the vertices changes nothing but the order of
/// processing. Proof sketch: a value is only ever set to the update of one triangle evaluated
/// with values that are at least the current ones, and the update is nondecreasing in them, so
/// a value is never below the update with the current values; conversely whenever a value
/// `T_j` decreases, every triangle at `j` is re-evaluated, so on termination no update is
/// below the value it feeds.
///
/// # Errors
/// As [`upwind_distance`].
pub fn upwind_distance_with(
    mesh: &TriMesh,
    geometry: &Geometry,
    seeds: &[(u32, f64)],
    class: Option<&[u32]>,
    mode: &DistanceMode,
) -> Result<DistanceResult, DistanceError> {
    geometry.check_mesh(mesh)?;
    let n = mesh.n_vertices();
    if let Some(class) = class
        && class.len() != n
    {
        return Err(DistanceError::ClassCountMismatch {
            expected: n,
            got: class.len(),
        });
    }
    let mut first_seed = vec![usize::MAX; n];
    for (s, &(vertex, value)) in seeds.iter().enumerate() {
        if vertex as usize >= n {
            return Err(DistanceError::SeedOutOfRange {
                seed: s,
                vertex,
                n_vertices: n,
            });
        }
        if !(value.is_finite() && value >= 0.0) {
            return Err(DistanceError::InvalidSeedValue {
                seed: s,
                vertex,
                value,
            });
        }
        if first_seed[vertex as usize] != usize::MAX {
            return Err(DistanceError::DuplicateSeed {
                vertex,
                first: first_seed[vertex as usize],
                second: s,
            });
        }
        first_seed[vertex as usize] = s;
    }

    let faces = mesh.faces();
    let cotans = geometry.corner_cotans();
    let areas = geometry.face_areas();
    // Side lengths per face: side[m] joins corner m to corner m + 1 (halfedge 3f + m).
    let sides: Vec<[f64; 3]> = (0..mesh.n_faces())
        .map(|f| [0, 1, 2].map(|m| geometry.edge_lengths()[mesh.edge(3 * f + m)]))
        .collect();

    let same_class = |a: u32, b: u32| class.is_none_or(|c| c[a as usize] == c[b as usize]);
    // The value offered to corner `ci` of face `f` by the accepted vertex `v` (corner `m`):
    // the triangle update if the third vertex is usable, the edge term from `v` otherwise.
    let offer = |values: &[f64], f: usize, ci: usize, m: usize, v: u32, other_usable: bool| {
        // Corner i has j at ci + 1 and k at ci + 2; one of them is v.
        let (cj, ck) = ((ci + 1) % 3, (ci + 2) % 3);
        let l = &sides[f];
        let (a, b, c) = (l[ci], l[ck], l[cj]); // ij, ik, jk
        if other_usable {
            update_core(
                values[faces[f][cj] as usize],
                values[faces[f][ck] as usize],
                [a, b, c],
                [cotans[f][cj], cotans[f][ck]],
                2.0 * areas[f] / c,
            )
        } else if cj == m {
            values[v as usize] + a
        } else {
            values[v as usize] + b
        }
    };

    let mut values = vec![f64::INFINITY; n];
    let mut heap: BinaryHeap<Entry> = BinaryHeap::with_capacity(seeds.len().max(16));
    for &(vertex, value) in seeds {
        values[vertex as usize] = value;
    }

    let report = match *mode {
        DistanceMode::SinglePass => {
            let mut accepted = vec![false; n];
            for &(vertex, _) in seeds {
                accepted[vertex as usize] = true;
            }
            // Offers `v` (just accepted) makes to the other corners of its faces.
            let relax =
                |v: u32, values: &mut Vec<f64>, accepted: &[bool], heap: &mut BinaryHeap<Entry>| {
                    for h in mesh.outgoing_halfedges(v as usize) {
                        let (f, m) = (h / 3, h % 3);
                        for ci in [(m + 1) % 3, (m + 2) % 3] {
                            let i = faces[f][ci];
                            if accepted[i as usize] || !same_class(i, v) {
                                continue;
                            }
                            let other = faces[f][if (ci + 1) % 3 == m {
                                (ci + 2) % 3
                            } else {
                                (ci + 1) % 3
                            }];
                            let usable = accepted[other as usize] && same_class(i, other);
                            let candidate = offer(values, f, ci, m, v, usable);
                            if candidate < values[i as usize] {
                                values[i as usize] = candidate;
                                heap.push(Reverse((candidate.to_bits(), i)));
                            }
                        }
                    }
                };
            for &(vertex, _) in seeds {
                relax(vertex, &mut values, &accepted, &mut heap);
            }
            while let Some(Reverse((bits, v))) = heap.pop() {
                if accepted[v as usize] || values[v as usize].to_bits() != bits {
                    continue; // stale entry
                }
                accepted[v as usize] = true;
                relax(v, &mut values, &accepted, &mut heap);
            }
            None
        }
        DistanceMode::FixedPoint {
            max_reopenings,
            reopen_tolerance,
        } => {
            let tolerance =
                if reopen_tolerance.is_finite() && (0.0..1.0).contains(&reopen_tolerance) {
                    reopen_tolerance
                } else {
                    0.0
                };
            let mut settled = vec![false; n];
            for &(vertex, value) in seeds {
                heap.push(Reverse((value.to_bits(), vertex)));
            }
            let mut n_reopened = 0usize;
            let mut converged = true;

            // Termination: every push is caused by a strict decrease of one value, and a value
            // only takes values from a finite set (the update of a triangle is a deterministic
            // function of the two other values and the lengths; the values are doubles in
            // [0, max seed + diameter]), so the number of pushes is finite. In exact arithmetic
            // the decreases stop long before that: a vertex can only be improved through a
            // strictly shorter polygonal path, and the paths that matter (through mesh edges
            // and across triangles) form a finite family. The budget bounds the work anyway.
            while let Some(Reverse((bits, v))) = heap.pop() {
                if values[v as usize].to_bits() != bits {
                    continue; // stale entry: the vertex received a smaller value afterwards
                }
                if n_reopened >= max_reopenings {
                    converged = false;
                    break;
                }
                settled[v as usize] = true;
                for h in mesh.outgoing_halfedges(v as usize) {
                    let (f, m) = (h / 3, h % 3);
                    // Update the two other corners of the face from v (corner m).
                    for ci in [(m + 1) % 3, (m + 2) % 3] {
                        let i = faces[f][ci];
                        if first_seed[i as usize] != usize::MAX || !same_class(i, v) {
                            continue;
                        }
                        let (cj, ck) = ((ci + 1) % 3, (ci + 2) % 3);
                        let other = faces[f][if cj == m { ck } else { cj }];
                        let usable = values[other as usize].is_finite() && same_class(i, other);
                        let candidate = offer(&values, f, ci, m, v, usable);
                        let current = values[i as usize];
                        // A settled vertex is re-opened only for a decrease beyond rounding
                        // noise (see `reopen_tolerance`).
                        let improves = if settled[i as usize] {
                            candidate < current * (1.0 - tolerance)
                        } else {
                            candidate < current
                        };
                        if improves {
                            values[i as usize] = candidate;
                            if settled[i as usize] {
                                settled[i as usize] = false;
                                n_reopened += 1;
                            }
                            heap.push(Reverse((candidate.to_bits(), i)));
                        }
                    }
                }
            }
            Some(FixedPointReport {
                n_reopened,
                converged,
            })
        }
    };
    let n_unreached = values.iter().filter(|x| x.is_infinite()).count();
    Ok(DistanceResult {
        values,
        n_unreached,
        fixed_point: report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Brute-force minimum over a fine grid of `t`, for checking the closed form.
    fn brute(tj: f64, tk: f64, a: f64, b: f64, c: f64) -> f64 {
        // Place j = (0,0), k = (c,0), i above: i = (x, y) with a = |i-j|, b = |i-k|.
        let x = (a * a + c * c - b * b) / (2.0 * c);
        let y = (a * a - x * x).max(0.0).sqrt();
        let mut best = f64::INFINITY;
        let n = 200_000;
        for s in 0..=n {
            let t = s as f64 / n as f64;
            let d = ((t * c - x).powi(2) + y * y).sqrt();
            best = best.min(tj + t * (tk - tj) + d);
        }
        best
    }

    #[test]
    fn closed_form_matches_brute_force() {
        let tris = [
            (1.0, 1.0, 1.0),
            (1.0, 2f64.sqrt(), 1.0),
            (3.0, 1.0, 2.5),
            (0.5, 3.0, 3.2),
            (1.0, 1.0, 1.9),
        ];
        for &(a, b, c) in &tris {
            for &(tj, tk) in &[
                (0.0, 0.0),
                (0.0, 0.3),
                (0.4, 0.0),
                (0.0, 5.0),
                (2.0, 2.0),
                (1.0, 0.2),
            ] {
                let v = triangle_update(tj, tk, a, b, c).unwrap();
                let w = brute(tj, tk, a, b, c);
                assert!((v - w).abs() < 1e-8, "{a} {b} {c} {tj} {tk}: {v} vs {w}");
            }
        }
    }

    #[test]
    fn degenerate_triangle_and_infinite_values_give_none() {
        assert!(triangle_update(0.0, 0.0, 1.0, 1.0, 2.0).is_none());
        assert!(triangle_update(f64::INFINITY, 0.0, 1.0, 1.0, 1.0).is_none());
    }
}
