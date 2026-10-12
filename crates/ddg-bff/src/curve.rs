//! Step II of Boundary First Flattening: from boundary lengths and exterior angles to a closed
//! plane curve (Sec. 4.5 of the paper, `BFF::constructBestFitCurve` and `BFF::closeLengths` of
//! the reference implementation).
//!
//! A boundary loop is described by its vertices `0, …, n − 1` in order, **with the surface on
//! the left** (counter-clockwise, the order of `TriMesh::boundary_loops`). Edge `i` runs from
//! vertex `i` to vertex `i + 1 (mod n)`. The curve starts at the origin and its first edge
//! points along `+x`.
//!
//! The reference implementation walks the boundary clockwise and mirrors the result; walking
//! counter-clockwise gives the same curve up to a rigid motion without the mirror step.

use std::f64::consts::PI;

/// Target lengths of the boundary edges: `ℓ*_i = e^{(u_i + u_{i+1}) / 2} ℓ_i` for edge `i` from
/// boundary vertex `i` to `i + 1 (mod n)`, with `u` the log scale factors at the boundary
/// vertices and `ℓ` the edge lengths in the input metric (paper Alg. 1, line 9).
///
/// [formula: bff-target-boundary-length]
pub fn target_boundary_lengths(scale_factors: &[f64], lengths: &[f64]) -> Vec<f64> {
    let n = lengths.len();
    (0..n)
        .map(|i| (0.5 * (scale_factors[i] + scale_factors[(i + 1) % n])).exp() * lengths[i])
        .collect()
}

/// Dual lengths of the boundary vertices: `ℓ_j = ½ (ℓ_{j−1} + ℓ_j)`, half of the two boundary
/// edges at vertex `j` (paper Sec. 3.3).
///
/// [formula: bff-dual-boundary-length]
pub fn dual_boundary_lengths(lengths: &[f64]) -> Vec<f64> {
    let n = lengths.len();
    (0..n)
        .map(|j| 0.5 * (lengths[(j + n - 1) % n] + lengths[j]))
        .collect()
}

/// Unit tangents of the boundary edges from the exterior angles: with `φ_0 = 0` and
/// `φ_i = φ_{i−1} + k̃_i`, edge `i` has direction `T_i = (cos φ_i, sin φ_i)` (paper Sec. 4.5.1).
/// The exterior angle at vertex `i` is the turn from edge `i − 1` to edge `i`; positive angles
/// turn left.
///
/// [formula: bff-boundary-tangents]
pub fn boundary_tangents(exterior_angles: &[f64]) -> Vec<[f64; 2]> {
    let mut phi = 0.0;
    exterior_angles
        .iter()
        .enumerate()
        .map(|(i, &k)| {
            if i > 0 {
                phi += k;
            }
            [phi.cos(), phi.sin()]
        })
        .collect()
}

/// Adds to each exterior angle its share of the defect `2π − Σ k̃`, in proportion to the dual
/// boundary lengths: `k̃_j ← k̃_j + (ℓ_j / Σ ℓ)(2π − Σ k̃)`. Afterwards the angles sum to `2π`
/// (`BFF::closeCurvatures` of the reference implementation). Returns the defect that was
/// distributed.
///
/// [formula: bff-close-exterior-angles]
pub fn close_exterior_angles(exterior_angles: &mut [f64], dual_lengths: &[f64]) -> f64 {
    let total_length: f64 = dual_lengths.iter().sum();
    let defect = 2.0 * PI - exterior_angles.iter().sum::<f64>();
    for (k, &l) in exterior_angles.iter_mut().zip(dual_lengths) {
        *k += l / total_length * defect;
    }
    defect
}

/// Vertex positions of a boundary polygon from edge lengths and tangents:
/// `γ_0 = 0`, `γ_{p} = Σ_{i < p} ℓ̃_i T_i` (paper Alg. 5, lines 9–11). Returns the positions of
/// the `n` vertices and the end point `γ_n`, which equals `γ_0 = 0` for a closed curve.
///
/// [formula: bff-boundary-curve]
pub fn integrate_curve(lengths: &[f64], tangents: &[[f64; 2]]) -> (Vec<[f64; 2]>, [f64; 2]) {
    let mut position = [0.0, 0.0];
    let mut out = Vec::with_capacity(lengths.len());
    for (&l, t) in lengths.iter().zip(tangents) {
        out.push(position);
        position[0] += l * t[0];
        position[1] += l * t[1];
    }
    (out, position)
}

/// How a boundary curve was closed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosureReport {
    /// `Σ k̃ − 2π`: how far the exterior angles are from one full turn. Zero up to rounding
    /// when the angles come from the discrete Cherrier formula (paper Proposition 1).
    pub turning_error: f64,
    /// `|Σ ℓ*_i T_i| / Σ ℓ*_i`: the gap of the curve built from the target lengths, relative to
    /// its length. This is the discretization error that closing has to remove.
    pub gap_before: f64,
    /// `|Σ ℓ̃_i T_i| / Σ ℓ̃_i`: the relative gap after the lengths were adjusted. Rounding
    /// level unless closing directions were dropped (see `dropped_directions`). This remaining
    /// gap is spread linearly along the loop, so the returned polygon is exactly closed.
    pub gap_after: f64,
    /// `max_i |ℓ̃_i / ℓ*_i − 1|`: the largest relative change of a boundary length.
    pub max_length_change: f64,
    /// `min_i ℓ̃_i / ℓ*_i`. A value `<= 0` means an edge of the boundary polygon was reversed.
    pub min_length_ratio: f64,
    /// Number of loops whose closing problems are coupled with this one through shared lengths
    /// (cut edges between different charts), including this loop.
    pub coupled_loops: usize,
    /// Number of closing constraints of that coupled problem (two per loop) that were not
    /// enforced: either they are redundant (the two sides of all cut edges are parallel, as on
    /// a flat torus) or enforcing them would change a length by more than the allowed
    /// fraction.
    pub dropped_directions: usize,
}

/// One boundary loop to lay out.
#[derive(Debug, Clone)]
pub(crate) struct LoopData {
    /// Per boundary edge: the id of its length degree of freedom. The two sides of a cut edge
    /// share one id; at most two edges (of all loops together) may have the same id.
    pub dof: Vec<u32>,
    /// Per boundary edge: its length in the input metric.
    pub length: Vec<f64>,
    /// Per boundary edge: its target length `ℓ*`.
    pub target: Vec<f64>,
    /// Per boundary vertex: its exterior angle `k̃`.
    pub exterior: Vec<f64>,
}

/// A closed boundary polygon.
#[derive(Debug, Clone)]
pub(crate) struct Curve {
    /// Per boundary vertex: its position.
    pub positions: Vec<[f64; 2]>,
    /// How the curve was closed.
    pub report: ClosureReport,
}

/// Options of [`best_fit_curves`].
#[derive(Debug, Clone, Copy)]
pub(crate) struct CurveOptions {
    /// Largest allowed relative change of a boundary length.
    pub max_length_adjustment: f64,
    /// Largest number of loops whose closing problems are solved together.
    pub max_coupled_loops: usize,
}

/// Pivots below this fraction of the total boundary length are treated as zero.
const RANK_TOLERANCE: f64 = 1e-12;

/// Builds closed polygons with the given exterior angles whose edge lengths are as close as
/// possible to the target lengths.
///
/// The lengths solve (paper Eq. 19, with one unknown per pair of cut sides as in Sec. 6.4)
///
/// ```text
/// minimize  ½ Σ_d N_d (ℓ̃_d − ℓ*_d)²   subject to   Σ_{i ∈ loop c} ℓ̃_{d(i)} T_i = 0  for every loop c,
/// ```
///
/// where `d(i)` is the degree of freedom of boundary edge `i` and `N_d = Σ 1/ℓ_i` over the one
/// or two edges with that degree of freedom. With `T` the `2C × D` matrix of the constraints
/// the minimizer is `ℓ̃ = ℓ* − N⁻¹ Tᵀ (T N⁻¹ Tᵀ)⁺ T ℓ*` (paper Eq. 20). For a single loop
/// without shared lengths this is exactly the formula of the paper.
///
/// [formula: bff-closing-lengths]
pub(crate) fn best_fit_curves(loops: &[LoopData], options: &CurveOptions) -> Vec<Curve> {
    let tangents: Vec<Vec<[f64; 2]>> = loops
        .iter()
        .map(|l| boundary_tangents(&l.exterior))
        .collect();

    // The sides of every degree of freedom, as (loop, edge).
    let n_dofs = loops
        .iter()
        .flat_map(|l| l.dof.iter())
        .map(|&d| d as usize + 1)
        .max()
        .unwrap_or(0);
    let mut sides: Vec<[(u32, u32); 2]> = vec![[(u32::MAX, 0); 2]; n_dofs];
    for (c, l) in loops.iter().enumerate() {
        for (i, &d) in l.dof.iter().enumerate() {
            let slot = &mut sides[d as usize];
            let k = usize::from(slot[0].0 != u32::MAX);
            debug_assert!(
                k == 0 || slot[1].0 == u32::MAX,
                "a length is shared by 3 edges"
            );
            slot[k] = (c as u32, i as u32);
        }
    }

    // Loops coupled through shared lengths.
    let mut group = Groups::new(loops.len());
    for slot in &sides {
        if slot[0].0 != u32::MAX && slot[1].0 != u32::MAX {
            group.union(slot[0].0, slot[1].0);
        }
    }
    let mut members: Vec<Vec<u32>> = vec![Vec::new(); loops.len()];
    for c in 0..loops.len() as u32 {
        members[group.find(c) as usize].push(c);
    }

    // Per boundary edge: the adjusted length.
    let mut adjusted: Vec<Vec<f64>> = loops.iter().map(|l| l.target.clone()).collect();
    let mut coupled = vec![1usize; loops.len()];
    let mut dropped = vec![0usize; loops.len()];
    let mut local = vec![u32::MAX; loops.len()];
    for group_members in members.iter().filter(|m| !m.is_empty()) {
        // Too large a group is solved loop by loop, without sharing lengths between loops.
        let together = group_members.len() <= options.max_coupled_loops.max(1);
        let problems: Vec<&[u32]> = if together {
            vec![group_members.as_slice()]
        } else {
            group_members.chunks(1).collect()
        };
        for problem in problems {
            for (k, &c) in problem.iter().enumerate() {
                local[c as usize] = k as u32;
            }
            let n_dropped = close_group(
                loops,
                &tangents,
                &sides,
                problem,
                &local,
                options.max_length_adjustment,
                &mut adjusted,
            );
            for &c in problem {
                coupled[c as usize] = problem.len();
                dropped[c as usize] = n_dropped;
                local[c as usize] = u32::MAX;
            }
        }
    }

    loops
        .iter()
        .enumerate()
        .map(|(c, l)| {
            let t = &tangents[c];
            let total_target: f64 = l.target.iter().sum();
            let (_, end_before) = integrate_curve(&l.target, t);
            let (mut positions, end) = integrate_curve(&adjusted[c], t);
            let total: f64 = adjusted[c].iter().sum();
            // Spread what is left of the gap linearly along the loop.
            let mut arc = 0.0;
            for (p, &len) in positions.iter_mut().zip(&adjusted[c]) {
                let s = if total > 0.0 { arc / total } else { 0.0 };
                p[0] -= s * end[0];
                p[1] -= s * end[1];
                arc += len;
            }
            let mut max_length_change = 0.0f64;
            let mut min_length_ratio = f64::INFINITY;
            for (&new, &target) in adjusted[c].iter().zip(&l.target) {
                let ratio = new / target;
                max_length_change = max_length_change.max((ratio - 1.0).abs());
                min_length_ratio = min_length_ratio.min(ratio);
            }
            Curve {
                positions,
                report: ClosureReport {
                    turning_error: l.exterior.iter().sum::<f64>() - 2.0 * PI,
                    gap_before: end_before[0].hypot(end_before[1]) / total_target,
                    gap_after: end[0].hypot(end[1]) / total,
                    max_length_change,
                    min_length_ratio,
                    coupled_loops: coupled[c],
                    dropped_directions: dropped[c],
                },
            }
        })
        .collect()
}

/// Solves the closing problem of the loops in `problem` together and writes the adjusted
/// lengths of their edges. A length shared with a loop outside `problem` is treated as a
/// length of its own. Returns the number of dropped constraint directions.
fn close_group(
    loops: &[LoopData],
    tangents: &[Vec<[f64; 2]>],
    sides: &[[(u32, u32); 2]],
    problem: &[u32],
    local: &[u32],
    max_length_adjustment: f64,
    adjusted: &mut [Vec<f64>],
) -> usize {
    let m = 2 * problem.len();
    // One entry per degree of freedom of the problem: its sides inside the problem.
    struct Dof {
        /// `1 / N_d`.
        weight: f64,
        target: f64,
        /// Up to two (row block, tangent, loop, edge).
        sides: [(u32, [f64; 2], u32, u32); 2],
        n_sides: usize,
    }
    let mut dofs: Vec<Dof> = Vec::new();
    for &c in problem {
        let l = &loops[c as usize];
        for (i, &d) in l.dof.iter().enumerate() {
            let slot = &sides[d as usize];
            let mine = (c, i as u32);
            let other = if slot[0] == mine { slot[1] } else { slot[0] };
            let other_inside = other.0 != u32::MAX && local[other.0 as usize] != u32::MAX;
            if other_inside && other < mine {
                continue; // already added from the other side
            }
            let mut dof = Dof {
                weight: 0.0,
                target: l.target[i],
                sides: [(0, [0.0; 2], 0, 0); 2],
                n_sides: 0,
            };
            let mut inverse_weight = 0.0;
            let mut add = |(c, i): (u32, u32), dof: &mut Dof| {
                inverse_weight += 1.0 / loops[c as usize].length[i as usize];
                dof.sides[dof.n_sides] =
                    (local[c as usize], tangents[c as usize][i as usize], c, i);
                dof.n_sides += 1;
            };
            add(mine, &mut dof);
            if other_inside {
                add(other, &mut dof);
            }
            dof.weight = 1.0 / inverse_weight;
            dofs.push(dof);
        }
    }

    // S = T N⁻¹ Tᵀ and r = T ℓ*.
    let mut s = vec![0.0; m * m];
    let mut r = vec![0.0; m];
    let mut scale = 0.0;
    let column = |dof: &Dof| -> ([usize; 4], [f64; 4], usize) {
        // The non-zero entries of the column of T that belongs to this degree of freedom.
        let mut rows = [0usize; 4];
        let mut values = [0.0; 4];
        let mut n = 0;
        for &(block, t, _, _) in &dof.sides[..dof.n_sides] {
            for (axis, &component) in t.iter().enumerate() {
                let row = 2 * block as usize + axis;
                match rows[..n].iter().position(|&x| x == row) {
                    Some(k) => values[k] += component,
                    None => {
                        rows[n] = row;
                        values[n] = component;
                        n += 1;
                    }
                }
            }
        }
        (rows, values, n)
    };
    for dof in &dofs {
        let (rows, values, n) = column(dof);
        for a in 0..n {
            r[rows[a]] += dof.target * values[a];
            for b in 0..n {
                s[rows[a] * m + rows[b]] += dof.weight * values[a] * values[b];
            }
        }
        scale += dof.weight * dof.n_sides as f64;
    }

    // Solve, dropping the weakest directions until the lengths change by an acceptable amount.
    let mut tolerance = RANK_TOLERANCE;
    loop {
        let (lambda, rank) = solve_semidefinite(&s, &r, m, tolerance * scale);
        let mut worst = 0.0f64;
        for dof in &dofs {
            let (rows, values, n) = column(dof);
            let dot: f64 = (0..n).map(|a| values[a] * lambda[rows[a]]).sum();
            let new = dof.target - dof.weight * dot;
            worst = worst.max((new / dof.target - 1.0).abs());
            for &(_, _, c, i) in &dof.sides[..dof.n_sides] {
                adjusted[c as usize][i as usize] = new;
            }
        }
        if rank == 0 || worst <= max_length_adjustment {
            return m - rank;
        }
        tolerance *= 100.0;
        if tolerance > 1.0 {
            // Give up closing by lengths altogether.
            for dof in &dofs {
                for &(_, _, c, i) in &dof.sides[..dof.n_sides] {
                    adjusted[c as usize][i as usize] = dof.target;
                }
            }
            return m;
        }
    }
}

/// Solves `S x = r` for a symmetric positive semidefinite `S` (`n × n`, row major) by Cholesky
/// factorization with diagonal pivoting, stopping when the largest remaining pivot is at most
/// `threshold`. The unknowns that were not eliminated are set to zero, which gives a solution
/// whenever `r` is in the range of the eliminated block. Returns `x` and the rank found.
fn solve_semidefinite(s: &[f64], r: &[f64], n: usize, threshold: f64) -> (Vec<f64>, usize) {
    let mut a = s.to_vec();
    let mut order: Vec<usize> = (0..n).collect();
    let mut rank = 0;
    // In-place factorization of the permuted matrix: row/column k of `a` belongs to unknown
    // order[k]; the lower triangle holds L.
    for k in 0..n {
        let (pivot, value) =
            (k..n)
                .map(|i| (i, a[i * n + i]))
                .fold(
                    (k, f64::NEG_INFINITY),
                    |best, x| if x.1 > best.1 { x } else { best },
                );
        if value.is_nan() || value <= threshold {
            break;
        }
        if pivot != k {
            order.swap(k, pivot);
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            for i in 0..n {
                a.swap(i * n + k, i * n + pivot);
            }
        }
        let diagonal = value.sqrt();
        a[k * n + k] = diagonal;
        for i in k + 1..n {
            a[i * n + k] /= diagonal;
        }
        for i in k + 1..n {
            let lik = a[i * n + k];
            if lik != 0.0 {
                for j in k + 1..=i {
                    a[i * n + j] -= lik * a[j * n + k];
                }
            }
        }
        // Keep the trailing block symmetric for the swaps of later steps.
        for i in k + 1..n {
            for j in k + 1..i {
                a[j * n + i] = a[i * n + j];
            }
        }
        rank += 1;
    }
    // Forward and backward substitution on the leading rank × rank block.
    let mut y: Vec<f64> = (0..rank).map(|k| r[order[k]]).collect();
    for i in 0..rank {
        for j in 0..i {
            y[i] -= a[i * n + j] * y[j];
        }
        y[i] /= a[i * n + i];
    }
    for i in (0..rank).rev() {
        for j in i + 1..rank {
            y[i] -= a[j * n + i] * y[j];
        }
        y[i] /= a[i * n + i];
    }
    let mut x = vec![0.0; n];
    for k in 0..rank {
        x[order[k]] = y[k];
    }
    (x, rank)
}

/// Minimal union-find for grouping loops.
struct Groups {
    parent: Vec<u32>,
}

impl Groups {
    fn new(n: usize) -> Self {
        Groups {
            parent: (0..n as u32).collect(),
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let grandparent = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = grandparent;
            x = grandparent;
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.parent[a.max(b) as usize] = a.min(b);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const OPTIONS: CurveOptions = CurveOptions {
        max_length_adjustment: 0.5,
        max_coupled_loops: 64,
    };

    fn square(target: [f64; 4]) -> LoopData {
        LoopData {
            dof: vec![0, 1, 2, 3],
            length: vec![1.0; 4],
            target: target.to_vec(),
            exterior: vec![PI / 2.0; 4],
        }
    }

    #[test]
    fn a_closed_curve_is_reproduced() {
        let curve = &best_fit_curves(&[square([1.0; 4])], &OPTIONS)[0];
        let expect = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        for (p, q) in curve.positions.iter().zip(&expect) {
            assert!((p[0] - q[0]).abs() < 1e-15 && (p[1] - q[1]).abs() < 1e-15);
        }
        assert!(curve.report.gap_before < 1e-15);
        assert!(curve.report.max_length_change < 1e-15);
        assert_eq!(curve.report.dropped_directions, 0);
        assert!(curve.report.turning_error.abs() < 1e-15);
    }

    #[test]
    fn lengths_are_adjusted_minimally_to_close() {
        // A rectangle whose opposite sides disagree: the closest closed one averages them.
        let curve = &best_fit_curves(&[square([1.0, 2.0, 1.2, 2.0])], &OPTIONS)[0];
        let p = &curve.positions;
        assert!((p[1][0] - 1.1).abs() < 1e-14, "{p:?}");
        assert!((p[2][1] - 2.0).abs() < 1e-14, "{p:?}");
        assert!(
            (p[3][0]).abs() < 1e-14 && (p[3][1] - 2.0).abs() < 1e-14,
            "{p:?}"
        );
        assert!(curve.report.gap_before > 0.03);
        assert!(curve.report.gap_after < 1e-15);
        assert!((curve.report.max_length_change - 0.1).abs() < 1e-14);
    }

    #[test]
    fn matches_the_formula_of_the_paper_for_one_loop() {
        // ℓ̃ = ℓ* − N⁻¹ Tᵀ (T N⁻¹ Tᵀ)⁻¹ T ℓ*, written out for a random pentagon.
        let exterior = vec![1.1, 1.4, 1.2, 1.5, 2.0 * PI - 5.2];
        let length = vec![1.0, 0.7, 1.3, 0.9, 1.1];
        let target = vec![1.05, 0.66, 1.4, 0.95, 1.0];
        let data = LoopData {
            dof: (0..5).collect(),
            length: length.clone(),
            target: target.clone(),
            exterior: exterior.clone(),
        };
        let t = boundary_tangents(&exterior);
        let mut m = [[0.0; 2]; 2];
        let mut rhs = [0.0; 2];
        for i in 0..5 {
            for a in 0..2 {
                rhs[a] += t[i][a] * target[i];
                for b in 0..2 {
                    m[a][b] += length[i] * t[i][a] * t[i][b];
                }
            }
        }
        let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        let lambda = [
            (m[1][1] * rhs[0] - m[0][1] * rhs[1]) / det,
            (m[0][0] * rhs[1] - m[1][0] * rhs[0]) / det,
        ];
        let expect: Vec<f64> = (0..5)
            .map(|i| target[i] - length[i] * (t[i][0] * lambda[0] + t[i][1] * lambda[1]))
            .collect();
        let (positions, end) = integrate_curve(&expect, &t);
        assert!(end[0].hypot(end[1]) < 1e-14);
        let curve = &best_fit_curves(&[data], &OPTIONS)[0];
        for (p, q) in curve.positions.iter().zip(&positions) {
            assert!((p[0] - q[0]).abs() < 1e-13 && (p[1] - q[1]).abs() < 1e-13);
        }
    }

    #[test]
    fn shared_lengths_stay_equal() {
        // A square whose edges 0 and 2 are the two sides of one cut edge (a cylinder cut open).
        let data = LoopData {
            dof: vec![0, 1, 0, 2],
            length: vec![1.0; 4],
            target: vec![1.0, 1.0, 1.0, 1.3],
            exterior: vec![PI / 2.0; 4],
        };
        let curve = &best_fit_curves(&[data], &OPTIONS)[0];
        let p = &curve.positions;
        let len = |a: usize, b: usize| (p[b][0] - p[a][0]).hypot(p[b][1] - p[a][1]);
        assert!((len(0, 1) - len(2, 3)).abs() < 1e-14);
        assert!((len(1, 2) - 1.15).abs() < 1e-14);
        assert!((len(3, 0) - 1.15).abs() < 1e-14);
        // The pair of parallel sides closes for any length: one direction is redundant.
        assert_eq!(curve.report.dropped_directions, 1);
        assert!(curve.report.gap_after < 1e-15);
    }

    #[test]
    fn two_loops_sharing_an_edge_agree_on_its_length() {
        // Two unit squares glued along one edge whose target lengths disagree with the rest.
        let a = LoopData {
            dof: vec![0, 1, 2, 3],
            length: vec![1.0; 4],
            target: vec![1.0, 1.2, 1.0, 1.0],
            exterior: vec![PI / 2.0; 4],
        };
        let b = LoopData {
            dof: vec![4, 5, 6, 1],
            length: vec![1.0; 4],
            target: vec![1.0, 0.9, 1.0, 1.2],
            exterior: vec![PI / 2.0; 4],
        };
        let curves = best_fit_curves(&[a.clone(), b.clone()], &OPTIONS);
        let len = |c: usize, i: usize| {
            let p = &curves[c].positions;
            let j = (i + 1) % 4;
            (p[j][0] - p[i][0]).hypot(p[j][1] - p[i][1])
        };
        assert!((len(0, 1) - len(1, 3)).abs() < 1e-14);
        assert_eq!(curves[0].report.coupled_loops, 2);
        assert!(curves.iter().all(|c| c.report.gap_after < 1e-15));
        // Solved one by one, the shared edge gets two lengths.
        let alone = best_fit_curves(
            &[a, b],
            &CurveOptions {
                max_coupled_loops: 1,
                ..OPTIONS
            },
        );
        let len = |c: usize, i: usize| {
            let p = &alone[c].positions;
            let j = (i + 1) % 4;
            (p[j][0] - p[i][0]).hypot(p[j][1] - p[i][1])
        };
        assert!((len(0, 1) - len(1, 3)).abs() > 0.01);
        assert_eq!(alone[0].report.coupled_loops, 1);
    }

    #[test]
    fn hopeless_closing_is_left_open_and_spread() {
        // All edges are paired and parallel, but the angles do not describe a closed curve:
        // lengths cannot close it, so the gap is spread and the polygon still returned closed.
        let data = LoopData {
            dof: vec![0, 1, 0, 1],
            length: vec![1.0; 4],
            target: vec![1.0; 4],
            exterior: vec![PI / 2.0, PI / 2.0 + 0.01, PI / 2.0 - 0.01, PI / 2.0],
        };
        let curve = &best_fit_curves(
            &[data],
            &CurveOptions {
                max_length_adjustment: 0.1,
                ..OPTIONS
            },
        )[0];
        assert!(curve.report.dropped_directions >= 1);
        assert!(curve.report.gap_after > 1e-3);
        assert!(curve.report.max_length_change <= 0.1);
        assert!(
            curve
                .positions
                .iter()
                .all(|p| p[0].is_finite() && p[1].is_finite())
        );
    }

    #[test]
    fn helpers() {
        assert_eq!(dual_boundary_lengths(&[1.0, 2.0, 4.0]), vec![2.5, 1.5, 3.0]);
        let l = target_boundary_lengths(&[0.0, 2.0f64.ln(), 0.0], &[1.0, 1.0, 3.0]);
        assert!((l[0] - 2.0f64.sqrt()).abs() < 1e-15);
        assert!((l[2] - 3.0).abs() < 1e-15);
        let mut k = vec![1.0, 2.0, 3.0];
        let defect = close_exterior_angles(&mut k, &[1.0, 1.0, 2.0]);
        assert!((defect - (2.0 * PI - 6.0)).abs() < 1e-15);
        assert!((k.iter().sum::<f64>() - 2.0 * PI).abs() < 1e-15);
        assert!((k[2] - 3.0 - 0.5 * defect).abs() < 1e-15);
    }
}
