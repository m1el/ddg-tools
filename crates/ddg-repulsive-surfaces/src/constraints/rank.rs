//! The translation and rank checks of a constraint set (`SPEC.md` §5.1, D10).
//!
//! **Rank (D10).** Rows are taken greedily, each kept only if it is linearly independent of the
//! rows kept before it: the relative residual `|r⊥|² / |r|²` of the row after projecting out the
//! span of the kept rows (Euclidean, in the `3|V|` coordinates) must exceed [`RANK_TOLERANCE`]
//! (the row-normalized test of the curves crate, `ddg-repulsive-curves` `constraints::rank`).
//! Order of preference: pins (a vertex pinned twice keeps its first pin; different targets are
//! [`RsError::InfeasibleConstraint`]), then barycenters, then area and volume rows, each in list
//! order. A dependent positional row is dropped with its whole group of three (positional rows
//! act on each axis alike, so their dependencies are the same on all three axes); so a barycenter
//! of a component whose vertices are all pinned is dropped, and a second `Area` of the same
//! component is dropped. A barycenter with some pins on its component is **not** dependent (it
//! constrains the weighted mean of the unpinned vertices) and is kept.
//!
//! **Translations.** `A₃` has the translations of every component in its kernel (§4.4); each
//! component needs an active barycenter or at least one active pin
//! ([`RsError::TranslationNotFixed`]). In [`super::BarycenterMode::Global`] a global barycenter
//! counts for every component (the reference's convention).

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use super::{BarycenterMode, ConstraintSet, RowKind};
use crate::error::RsError;

/// A row whose relative residual `|r⊥|²/|r|²` against the rows kept before it is at most this is
/// linearly dependent (the curves crate's `TRANSLATION_RANK_TOLERANCE`).
pub const RANK_TOLERANCE: f64 = 1e-10;

/// Result of the rank check.
#[derive(Debug, Clone, PartialEq)]
pub struct RankReport {
    /// Rows dropped as linearly dependent on rows kept before them (indices into
    /// [`ConstraintSet::rows`]).
    pub dropped_rows: Vec<usize>,
    /// Number of active rows.
    pub n_active_rows: usize,
    /// The smallest relative residual `|r⊥|²/|r|²` among the kept rows (1 without rows): how
    /// close the kept rows are to a dependency.
    pub min_relative_residual: f64,
}

impl Default for RankReport {
    fn default() -> Self {
        RankReport {
            dropped_rows: Vec::new(),
            n_active_rows: 0,
            min_relative_residual: 1.0,
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Projects `r` onto the orthogonal complement of the orthonormal `basis` (modified
/// Gram–Schmidt, twice for stability).
fn orthogonalize(r: &mut [f64], basis: &[Vec<f64>]) {
    for _ in 0..2 {
        for q in basis {
            let c = dot(r, q);
            for (x, y) in r.iter_mut().zip(q) {
                *x -= c * y;
            }
        }
    }
}

/// The greedy rank check (module docs) at `positions`; all rows of `set` are candidates.
pub(super) fn analyze(set: &ConstraintSet, positions: &[[f64; 3]]) -> Result<RankReport, RsError> {
    let nv = set.n_vertices;
    let rows = &set.rows;
    let mut dropped: Vec<usize> = Vec::new();
    let mut min_res: f64 = 1.0;
    // 1. Pins (groups of three rows).
    let mut pin_of: Vec<Option<usize>> = vec![None; nv];
    let mut r = 0;
    while r < rows.len() {
        if let RowKind::Pin { vertex, axis: 0 } = rows[r].kind {
            let v = vertex as usize;
            match pin_of[v] {
                Some(first) => {
                    let mut diff: f64 = 0.0;
                    for k in 0..3 {
                        diff = diff.max((rows[first + k].target - rows[r + k].target).abs());
                    }
                    if diff > 0.0 {
                        return Err(RsError::InfeasibleConstraint {
                            constraint: rows[r].constraint,
                            violation: diff,
                        });
                    }
                    dropped.extend(r..r + 3);
                }
                None => pin_of[v] = Some(r),
            }
            r += 3;
        } else {
            r += 1;
        }
    }
    let pinned: Vec<bool> = pin_of.iter().map(Option::is_some).collect();
    // 2. Barycenters: scalar rows over the vertices (the same on each axis).
    let va = set.vertex_areas(positions);
    let mut bary_basis: Vec<Vec<f64>> = Vec::new();
    let mut r = 0;
    while r < rows.len() {
        if let RowKind::Barycenter { component, axis: 0 } = rows[r].kind {
            let w: Vec<f64> = (0..nv)
                .map(|v| if set.in_component(v, component) { va[v] } else { 0.0 })
                .collect();
            let n2 = dot(&w, &w);
            let mut res: Vec<f64> = w
                .iter()
                .zip(&pinned)
                .map(|(x, &p)| if p { 0.0 } else { *x })
                .collect();
            orthogonalize(&mut res, &bary_basis);
            let rel = dot(&res, &res) / n2;
            if !(rel > RANK_TOLERANCE) {
                dropped.extend(r..r + 3);
            } else {
                min_res = min_res.min(rel);
                let inv = 1.0 / dot(&res, &res).sqrt();
                res.iter_mut().for_each(|x| *x *= inv);
                bary_basis.push(res);
            }
            r += 3;
        } else {
            r += 1;
        }
    }
    // 3. Area and volume rows, in 3|V| coordinates, against the positional span and each other.
    let nonlinear: Vec<usize> = (0..rows.len())
        .filter(|&r| !rows[r].kind.is_positional())
        .collect();
    if !nonlinear.is_empty() {
        let mut basis: Vec<Vec<f64>> = Vec::new();
        for &r in &nonlinear {
            let g = match rows[r].kind {
                RowKind::Area { component } => set.area_gradient(positions, component),
                RowKind::Volume { component } => set.volume_gradient(positions, component),
                _ => unreachable!("positional rows filtered out"),
            };
            let full: Vec<f64> = g.iter().flatten().copied().collect();
            let n2 = dot(&full, &full);
            // Positional span, per axis: pinned coordinates and the barycenter basis.
            let mut res = vec![0.0; 3 * nv];
            for k in 0..3 {
                let mut ax: Vec<f64> = (0..nv)
                    .map(|v| if pinned[v] { 0.0 } else { g[v][k] })
                    .collect();
                orthogonalize(&mut ax, &bary_basis);
                for v in 0..nv {
                    res[3 * v + k] = ax[v];
                }
            }
            orthogonalize(&mut res, &basis);
            let rel = if n2 > 0.0 { dot(&res, &res) / n2 } else { 0.0 };
            if !(rel > RANK_TOLERANCE) {
                dropped.push(r);
            } else {
                min_res = min_res.min(rel);
                let inv = 1.0 / dot(&res, &res).sqrt();
                res.iter_mut().for_each(|x| *x *= inv);
                basis.push(res);
            }
        }
    }
    dropped.sort_unstable();
    Ok(RankReport {
        n_active_rows: rows.len() - dropped.len(),
        dropped_rows: dropped,
        min_relative_residual: min_res,
    })
}

/// Every component needs an active barycenter or an active pin (module docs).
pub(super) fn check_translations(set: &ConstraintSet) -> Result<(), RsError> {
    let mut fixed = vec![false; set.n_components];
    let vc = set.vertex_component();
    for (r, row) in set.rows.iter().enumerate() {
        if !set.active[r] {
            continue;
        }
        match row.kind {
            RowKind::Barycenter {
                component: Some(c), ..
            } => fixed[c as usize] = true,
            RowKind::Barycenter {
                component: None, ..
            } => {
                debug_assert_eq!(set.mode, BarycenterMode::Global);
                fixed.iter_mut().for_each(|f| *f = true);
            }
            RowKind::Pin { vertex, .. } => fixed[vc[vertex as usize] as usize] = true,
            _ => {}
        }
    }
    match fixed.iter().position(|f| !f) {
        Some(c) => Err(RsError::TranslationNotFixed {
            component: c as u32,
        }),
        None => Ok(()),
    }
}
