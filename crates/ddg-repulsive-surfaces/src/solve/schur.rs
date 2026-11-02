//! Schur complement direction and Newton correction (`SPEC.md` §6.1–6.2, Sec. 5.3, 7.2).
//!
//! The positional rows `C_Q` are inside `A′` (any [`SaddleSolver`]); the `k` nonlinear rows
//! `C_N` (area, volume) go through the dense `k × k` Schur complement. With
//! `Z = A′⁻¹ C_Nᵀ` (position block, `k` solves; `C_N` padded with zeros on the multiplier
//! columns) and `S = −C_N Z`:
//!
//! ```text
//! direction   x = y + Z S⁻¹ C_N y,     y = A′⁻¹ g                (C_N x = 0, C_Q x = 0)
//! correction  h = −Z S⁻¹ Φ_N                                       (C_N h = Φ_N, C_Q h = 0)
//! ```
//!
//! (`rs-schur-direction`, `rs-schur-correction`; together they solve the full saddle with all
//! rows, `rs-nested-saddle`). The step is `f ← f − τ x`; the correction is applied as
//! `f ← f − h`, reusing `Z` and `S` of the step (old Jacobian).
//!
//! `N = −S = C_N (A′⁻¹)₁₁ C_Nᵀ` is positive semidefinite (the position block of `A′⁻¹` is the
//! inverse of `A₃` on `ker C_Q`). It is symmetrized, equilibrated to unit diagonal and factored
//! by Cholesky with diagonal pivoting (largest remaining pivot first, ties by row order); a row
//! whose pivot falls to [`SCHUR_DROP_TOLERANCE`] or below is dropped for this step and reported
//! (§6.2, D10): area and volume together on a near-sphere keep one row instead of producing a
//! huge correction.

// `!(x > tol)` is deliberate throughout: a NaN must fail the test.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use super::{LinearReport, SaddleSolver, SolveRecord};
use crate::constraints::ConstraintSet;
use crate::error::RsError;
use ddg_core::SparseMatrix;

/// Equilibrated Schur pivots at or below this drop their row for the step (§6.2).
pub const SCHUR_DROP_TOLERANCE: f64 = 1e-10;

/// The Schur complement of the nonlinear rows of one step.
#[derive(Debug, Clone)]
pub struct SchurComplement {
    /// `3|V|`.
    n_dofs: usize,
    /// The constraint-set rows of `C_N`, in order (all candidates).
    rows: Vec<usize>,
    /// `C_N` (`k × 3|V|`).
    c: SparseMatrix,
    /// `Z = A′⁻¹ C_Nᵀ`, position block, one vector per row.
    z: Vec<Vec<f64>>,
    /// `S = −C_N Z` as computed (`k × k`, row-major, before symmetrization).
    s: Vec<f64>,
    /// Positions (into `rows`) of the kept rows, in pivot order.
    kept: Vec<usize>,
    /// Lower Cholesky factor of the equilibrated `N` restricted to `kept` (pivot order),
    /// row-major `kept.len()²`.
    l: Vec<f64>,
    /// Equilibration `1/√N_ii` per row.
    d: Vec<f64>,
    /// Smallest kept pivot (1 without rows).
    min_pivot: f64,
    /// `max |S − Sᵀ| / max |S|` (inexact solves).
    asymmetry: f64,
    /// The solves that built `Z`.
    records: Vec<SolveRecord>,
}

impl SchurComplement {
    /// Builds `Z` (`k` solves with `solver`) and `S` for the nonlinear rows `rows` of `set`
    /// with the Jacobian `c` (`rows.len() × 3|V|`, e.g. `set.jacobian(positions, &rows)`).
    ///
    /// [formula: rs-schur-direction]
    ///
    /// # Errors
    /// [`RsError::InvalidOption`] if `c` does not match; the solver's errors.
    pub fn new(
        solver: &mut dyn SaddleSolver,
        rows: Vec<usize>,
        c: SparseMatrix,
    ) -> Result<SchurComplement, RsError> {
        let nd = solver.n_dofs();
        let k = rows.len();
        if c.nrows() != k || c.ncols() != nd {
            return Err(RsError::InvalidOption {
                name: "Schur Jacobian",
                value: format!("{} x {}", c.nrows(), c.ncols()),
                requirement: "one row per nonlinear row, 3|V| columns",
            });
        }
        let ct = c.transpose();
        let zero_b = vec![0.0; solver.n_constraints()];
        let cols = Self::columns(&ct, k, nd);
        let rhs: Vec<(&[f64], &[f64])> = cols.iter().map(|a| (a.as_slice(), zero_b.as_slice())).collect();
        let sols = solver.solve_many(&rhs)?;
        Ok(Self::from_solutions(nd, rows, c, ct, sols))
    }

    /// [`SchurComplement::new`] and [`SchurComplement::direction`] for `g` together: the `k`
    /// Schur columns and `y = A′⁻¹ g` are solved as one batch ([`SaddleSolver::solve_many`]),
    /// with the same results.
    ///
    /// [formula: rs-schur-direction]
    ///
    /// # Errors
    /// As [`SchurComplement::new`] and [`SchurComplement::direction`].
    pub fn with_direction(
        solver: &mut dyn SaddleSolver,
        rows: Vec<usize>,
        c: SparseMatrix,
        g: &[f64],
    ) -> Result<(SchurComplement, Direction), RsError> {
        let nd = solver.n_dofs();
        let k = rows.len();
        if c.nrows() != k || c.ncols() != nd {
            return Err(RsError::InvalidOption {
                name: "Schur Jacobian",
                value: format!("{} x {}", c.nrows(), c.ncols()),
                requirement: "one row per nonlinear row, 3|V| columns",
            });
        }
        let ct = c.transpose();
        let zero_b = vec![0.0; solver.n_constraints()];
        let cols = Self::columns(&ct, k, nd);
        let mut rhs: Vec<(&[f64], &[f64])> = cols.iter().map(|a| (a.as_slice(), zero_b.as_slice())).collect();
        rhs.push((g, zero_b.as_slice()));
        let mut sols = solver.solve_many(&rhs)?;
        let ysol = sols.pop().expect("the direction's solve");
        let schur = Self::from_solutions(nd, rows, c, ct, sols);
        let dir = schur.direction_from(ysol);
        Ok((schur, dir))
    }

    /// The right-hand sides `C_Nᵀ e_r` (dense, `3|V|`).
    fn columns(ct: &SparseMatrix, k: usize, nd: usize) -> Vec<Vec<f64>> {
        let (cp, ri, val) = (ct.col_ptr(), ct.row_idx(), ct.values());
        (0..k)
            .map(|r| {
                let mut a = vec![0.0; nd];
                for p in cp[r]..cp[r + 1] {
                    a[ri[p] as usize] = val[p];
                }
                a
            })
            .collect()
    }

    fn from_solutions(
        nd: usize,
        rows: Vec<usize>,
        c: SparseMatrix,
        ct: SparseMatrix,
        sols: Vec<super::SaddleSolution>,
    ) -> SchurComplement {
        let k = rows.len();
        let (cp, ri, val) = (ct.col_ptr(), ct.row_idx(), ct.values());
        let mut z = Vec::with_capacity(k);
        let mut records = Vec::with_capacity(k);
        for sol in sols {
            z.push(sol.x);
            records.push(sol.record);
        }
        // S = −C_N Z.
        let mut s = vec![0.0; k * k];
        for r in 0..k {
            for q in 0..k {
                let mut acc = 0.0;
                for p in cp[r]..cp[r + 1] {
                    acc += val[p] * z[q][ri[p] as usize];
                }
                s[r * k + q] = -acc;
            }
        }
        let mut out = SchurComplement {
            n_dofs: nd,
            rows,
            c,
            z,
            s,
            kept: Vec::new(),
            l: Vec::new(),
            d: vec![1.0; k],
            min_pivot: 1.0,
            asymmetry: 0.0,
            records,
        };
        out.factor();
        out
    }

    /// Equilibrated, pivoted Cholesky of `N = −½(S + Sᵀ)` with row dropping (module docs).
    fn factor(&mut self) {
        let k = self.rows.len();
        let s = &self.s;
        let smax = s.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let mut asym: f64 = 0.0;
        for i in 0..k {
            for j in 0..i {
                asym = asym.max((s[i * k + j] - s[j * k + i]).abs());
            }
        }
        self.asymmetry = if smax > 0.0 { asym / smax } else { 0.0 };
        let n = |i: usize, j: usize| -0.5 * (s[i * k + j] + s[j * k + i]);
        for i in 0..k {
            let nii = n(i, i);
            self.d[i] = if nii > 0.0 && nii.is_finite() {
                1.0 / nii.sqrt()
            } else {
                0.0
            };
        }
        // Working matrix (equilibrated), the Schur complement of the kept rows so far.
        let mut w: Vec<f64> = (0..k * k)
            .map(|idx| {
                let (i, j) = (idx / k, idx % k);
                n(i, j) * self.d[i] * self.d[j]
            })
            .collect();
        let mut remaining: Vec<usize> = (0..k).filter(|&i| self.d[i] > 0.0).collect();
        let mut kept: Vec<usize> = Vec::new();
        // Columns of L (in the original row indexing) for the kept pivots.
        let mut lcols: Vec<Vec<f64>> = Vec::new();
        let mut min_pivot: f64 = 1.0;
        while !remaining.is_empty() {
            // Largest remaining pivot, ties by row order.
            let (pos, &p) = remaining
                .iter()
                .enumerate()
                .fold(None::<(usize, &usize)>, |best, (pos, i)| match best {
                    Some((_, b)) if w[b * k + b] >= w[i * k + i] => best,
                    _ => Some((pos, i)),
                })
                .expect("non-empty");
            let piv = w[p * k + p];
            if !(piv > SCHUR_DROP_TOLERANCE) {
                break;
            }
            min_pivot = min_pivot.min(piv);
            remaining.remove(pos);
            let sq = piv.sqrt();
            let mut col = vec![0.0; k];
            col[p] = sq;
            for &i in &remaining {
                col[i] = w[i * k + p] / sq;
            }
            for &i in &remaining {
                for &j in &remaining {
                    w[i * k + j] -= col[i] * col[j];
                }
            }
            kept.push(p);
            lcols.push(col);
        }
        let m = kept.len();
        let mut l = vec![0.0; m * m];
        for (b, col) in lcols.iter().enumerate() {
            for (a, &i) in kept.iter().enumerate() {
                l[a * m + b] = col[i];
            }
        }
        self.kept = kept;
        self.l = l;
        self.min_pivot = if m == 0 && k > 0 { 0.0 } else { min_pivot };
    }

    /// `N_kept⁻¹ v` for `v` indexed by `rows` (entries of dropped rows ignored); result indexed
    /// by `rows` (0 on dropped rows).
    fn solve_n(&self, v: &[f64]) -> Vec<f64> {
        let k = self.rows.len();
        let m = self.kept.len();
        // Ñ = D N D, N⁻¹ = D Ñ⁻¹ D; Ñ_kept = L Lᵀ in pivot order.
        let mut y: Vec<f64> = self.kept.iter().map(|&i| self.d[i] * v[i]).collect();
        for a in 0..m {
            let mut s = y[a];
            for b in 0..a {
                s -= self.l[a * m + b] * y[b];
            }
            y[a] = s / self.l[a * m + a];
        }
        for a in (0..m).rev() {
            let mut s = y[a];
            for b in a + 1..m {
                s -= self.l[b * m + a] * y[b];
            }
            y[a] = s / self.l[a * m + a];
        }
        let mut out = vec![0.0; k];
        for (a, &i) in self.kept.iter().enumerate() {
            out[i] = self.d[i] * y[a];
        }
        out
    }

    /// `Σ_r c_r Z_r`.
    fn combine(&self, c: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.n_dofs];
        for (cr, zr) in c.iter().zip(&self.z) {
            if *cr != 0.0 {
                for (o, x) in out.iter_mut().zip(zr) {
                    *o += cr * x;
                }
            }
        }
        out
    }

    /// `C_N v` for a position vector.
    pub fn apply_c(&self, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.rows.len()];
        if self.c.mul_vec(v, &mut out).is_err() {
            out.iter_mut().for_each(|x| *x = f64::NAN);
        }
        out
    }

    /// The direction `x = y + Z S⁻¹ C_N y`, `y = A′⁻¹ g` (one more solve; `g` interleaved,
    /// `3|V|`). Dropped rows are not enforced.
    ///
    /// [formula: rs-schur-direction]
    /// [formula: rs-nested-saddle]
    ///
    /// # Errors
    /// The solver's errors.
    pub fn direction(&self, solver: &mut dyn SaddleSolver, g: &[f64]) -> Result<Direction, RsError> {
        let zero_b = vec![0.0; solver.n_constraints()];
        let sol = solver.solve(g, &zero_b)?;
        Ok(self.direction_from(sol))
    }

    /// The direction from the solution `y = A′⁻¹ g`.
    fn direction_from(&self, sol: super::SaddleSolution) -> Direction {
        let y = sol.x;
        // S⁻¹ = −N⁻¹: x = y − Z N⁻¹ C_N y.
        let cy = self.apply_c(&y);
        let c = self.solve_n(&cy);
        let zc = self.combine(&c);
        let x: Vec<f64> = y.iter().zip(&zc).map(|(a, b)| a - b).collect();
        let cnx = self.apply_c(&x);
        let max_cnx = self
            .kept
            .iter()
            .fold(0.0f64, |m, &i| m.max(cnx[i].abs()));
        Direction {
            x,
            lambda: sol.lambda,
            record: sol.record,
            max_nonlinear_residual: max_cnx,
        }
    }

    /// The Newton correction `h = −Z S⁻¹ Φ_N` (`C_N h = Φ_N` on the kept rows, `C_Q h = 0`), to
    /// be applied as `f ← f − h`. `phi` is indexed like the Schur rows (dropped rows ignored).
    ///
    /// [formula: rs-schur-correction]
    pub fn correction(&self, phi: &[f64]) -> Vec<f64> {
        // −S⁻¹ = N⁻¹.
        self.combine(&self.solve_n(phi))
    }

    /// The constraint-set rows of `C_N` (all candidates), in order.
    pub fn rows(&self) -> &[usize] {
        &self.rows
    }

    /// The constraint-set rows dropped for this step (pivot at or below
    /// [`SCHUR_DROP_TOLERANCE`]), increasing.
    pub fn dropped_rows(&self) -> Vec<usize> {
        let mut out: Vec<usize> = (0..self.rows.len())
            .filter(|i| !self.kept.contains(i))
            .map(|i| self.rows[i])
            .collect();
        out.sort_unstable();
        out
    }

    /// Positions (into [`SchurComplement::rows`]) of the kept rows, increasing.
    pub fn kept(&self) -> Vec<usize> {
        let mut k = self.kept.clone();
        k.sort_unstable();
        k
    }

    /// `S = −C_N A′⁻¹ C_Nᵀ` as computed (`k × k`, row-major).
    pub fn matrix(&self) -> &[f64] {
        &self.s
    }

    /// `Z = A′⁻¹ C_Nᵀ` (position block), one vector per row.
    pub fn z(&self) -> &[Vec<f64>] {
        &self.z
    }

    /// The smallest kept equilibrated pivot (1 without rows, 0 if every row was dropped).
    pub fn min_pivot(&self) -> f64 {
        self.min_pivot
    }

    /// `max |S − Sᵀ| / max |S|` (nonzero with inexact solves).
    pub fn asymmetry(&self) -> f64 {
        self.asymmetry
    }

    /// The solves that built `Z`.
    pub fn records(&self) -> &[SolveRecord] {
        &self.records
    }

    /// Adds the Schur data to a step's report.
    pub fn fill_report(&self, report: &mut LinearReport) {
        report.solves.extend(self.records.iter().cloned());
        report.dropped_rows = self.dropped_rows();
        report.schur_pivot = self.min_pivot;
    }
}

/// The constrained direction of §6.1.
#[derive(Debug, Clone, PartialEq)]
pub struct Direction {
    /// `x` (interleaved, `3|V|`); the step is `f ← f − τ x`.
    pub x: Vec<f64>,
    /// The positional multipliers of the solve `y = A′⁻¹ g`.
    pub lambda: Vec<f64>,
    /// The solve `y = A′⁻¹ g`.
    pub record: SolveRecord,
    /// `max |C_N x|` over the kept rows (round-off or solver tolerance).
    pub max_nonlinear_residual: f64,
}

/// How many Newton rounds the projection may take (`SPEC.md` §6.2, D12).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonSteps {
    /// Exactly this many rounds (the reference: 1), each still subject to the safeguard.
    Fixed(u32),
    /// Repeat while a kept row is above its tolerance, at most `max` rounds (default 5).
    UntilTolerance {
        /// The cap.
        max: u32,
    },
}

impl Default for NewtonSteps {
    fn default() -> Self {
        NewtonSteps::UntilTolerance { max: 5 }
    }
}

/// Options of [`project`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ProjectionOptions {
    /// The Newton rounds.
    pub newton: NewtonSteps,
}

/// How a projection ended (`SPEC.md` §6.2). Not an error in any case.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionStatus {
    /// Every kept nonlinear row within its tolerance.
    Converged,
    /// The rounds were used up (`Fixed` or the cap) with a row above its tolerance.
    MaxRounds,
    /// A round did not lower `max_r |Φ_r|/tol_r`; it was undone and the rounds stopped.
    NotImproving,
    /// The projection would have crossed an obstacle or the surface itself (or made the
    /// objective infinite); only part of its motion was applied (the step's
    /// `projection_fraction`; set by the session, `SPEC.md` §6.2).
    Blocked,
}

/// Result of [`project`].
#[derive(Debug, Clone, PartialEq)]
pub struct Projection {
    /// The projected positions.
    pub positions: Vec<[f64; 3]>,
    /// How it ended.
    pub status: ProjectionStatus,
    /// Newton rounds kept.
    pub rounds: u32,
    /// `max_r |Φ_r| / tol_r` over the kept nonlinear rows: after the positional projection of
    /// the input, then after every kept round, then (if one was undone) the rejected value.
    pub history: Vec<f64>,
    /// `max_r |Φ_r|` over all active rows at the end (positional and nonlinear, kept or not).
    pub max_violation: f64,
    /// `Φ` of every row of the set at the end.
    pub phi: Vec<f64>,
    /// Rows dropped from the Schur complement for this step.
    pub dropped_rows: Vec<usize>,
}

/// The projection of §6.2: the positional rows exactly, then Newton rounds
/// `f ← f − h`, `h = −Z S⁻¹ Φ_N` with the step's `schur` (old Jacobian), each followed by the
/// positional projection. A round is kept only if it lowers `max_r |Φ_r|/tol_r` over the kept
/// rows (the safeguard); otherwise it is undone and the rounds stop (`NotImproving`).
/// `UntilTolerance` stops as soon as every kept row is within its tolerance; `Fixed(n)` takes
/// `n` rounds. A violation above tolerance at the end is reported, not an error.
///
/// [formula: rs-schur-correction]
///
/// # Errors
/// [`RsError::PositionCountMismatch`]; [`RsError::InvalidOption`] if `schur` was not built
/// for rows of `set`.
pub fn project(
    set: &ConstraintSet,
    schur: Option<&SchurComplement>,
    positions: &[[f64; 3]],
    options: &ProjectionOptions,
) -> Result<Projection, RsError> {
    let mut pos = positions.to_vec();
    set.project_positional(&mut pos)?;
    let rows: Vec<usize> = schur.map_or_else(Vec::new, |s| s.rows.clone());
    if rows.iter().any(|&r| r >= set.n_rows()) {
        return Err(RsError::InvalidOption {
            name: "schur",
            value: format!("{rows:?}"),
            requirement: "rows of the constraint set",
        });
    }
    let kept: Vec<usize> = schur.map_or_else(Vec::new, |s| s.kept());
    let tol: Vec<f64> = rows.iter().map(|&r| set.tolerance(r)).collect();
    let measure = |phi: &[f64]| -> f64 {
        kept.iter().fold(0.0f64, |m, &i| {
            let v = phi[rows[i]].abs() / tol[i];
            if v.is_nan() || m.is_nan() { f64::NAN } else { m.max(v) }
        })
    };
    let mut phi = set.values(&pos)?;
    let mut m = measure(&phi);
    let mut history = vec![m];
    let mut rounds = 0u32;
    let (cap, until) = match options.newton {
        NewtonSteps::Fixed(n) => (n, false),
        NewtonSteps::UntilTolerance { max } => (max, true),
    };
    let mut stopped = false;
    if let Some(s) = schur {
        while rounds < cap && !kept.is_empty() {
            if until && m <= 1.0 {
                break;
            }
            let phin: Vec<f64> = rows.iter().map(|&r| phi[r]).collect();
            let h = s.correction(&phin);
            let mut trial = pos.clone();
            for (v, p) in trial.iter_mut().enumerate() {
                for k in 0..3 {
                    p[k] -= h[3 * v + k];
                }
            }
            set.project_positional(&mut trial)?;
            let phi_t = set.values(&trial)?;
            let m_t = measure(&phi_t);
            if !(m_t < m) {
                history.push(m_t);
                stopped = true;
                break;
            }
            pos = trial;
            phi = phi_t;
            m = m_t;
            history.push(m);
            rounds += 1;
        }
    }
    let status = if m <= 1.0 || kept.is_empty() {
        ProjectionStatus::Converged
    } else if stopped {
        ProjectionStatus::NotImproving
    } else {
        ProjectionStatus::MaxRounds
    };
    let max_violation = (0..set.n_rows())
        .filter(|&r| set.is_active(r))
        .fold(0.0f64, |acc, r| acc.max(phi[r].abs()));
    Ok(Projection {
        positions: pos,
        status,
        rounds,
        history,
        max_violation,
        phi,
        dropped_rows: schur.map_or_else(Vec::new, SchurComplement::dropped_rows),
    })
}
