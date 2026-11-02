//! A coarse correction of the saddle preconditioner on the per-component translations
//! (`SPEC.md` §6.3, `NOTES-SOLVERFIX.md`).
//!
//! With one global barycenter on several components (`BarycenterMode::Global`, the default
//! with `B⁰`), the translations of each component are not fixed by rows of their own, so the
//! preconditioner's `L′ = [[W + ε_multi M, C_Qᵀ], [C_Q, 0]]` sees them only through `ε_multi`:
//! `L′⁻¹` scales them by `1/ε_multi` and the preconditioned operator has `3 (n_c − 1)` outlying
//! eigenvalues anywhere between `10²` and `10⁹` (∝ `ε_multi⁻²`; measured on 14 spheres:
//! 39 of them, the rest of the spectrum in `[0.36, 1]`). GMRES then needs 50–60 iterations per
//! solve instead of the 4–9 of one sphere or of per-component barycenters.
//!
//! The correction is the balancing (deflation) preconditioner on the space `Z` spanned by the
//! translations of the free components (positions) and the multiplier coordinates of the rows
//! that touch them:
//!
//! ```text
//! P_c⁻¹ = Z E⁺ Zᵀ + (I − Z E⁺ Zᵀ A′) P⁻¹ (I − A′ Z E⁺ Zᵀ),      E = Zᵀ A′ Z
//! ```
//!
//! `P_c⁻¹ A′ Z = Z`: the coarse directions get the eigenvalue 1 exactly, whatever `ε_multi`.
//! `A′ Z` is computed once (`⌈n_c/3⌉` products: `A₃ = A ⊗ I₃`, three components in the three
//! axes of one product) and kept dense; `Zᵀ A′ z = (A′ Z)ᵀ z` uses the symmetry of `A′`, so an
//! application costs two dense `n × m` products besides `P⁻¹`. `E` is inverted by a symmetric
//! eigendecomposition with a relative cutoff (pseudo-inverse; `E` is singular only when the
//! coarse rows are dependent).

use ddg_core::SparseMatrix;
use faer::{Mat, Side};

use super::gmres::LinearOperator;
use crate::error::RsError;
use crate::surface::Surface;

/// The largest `dim · m` (entries of the dense `A′ Z`) for which a coarse space is built
/// (2²⁵ entries = 256 MB).
pub const COARSE_MAX_ENTRIES: usize = 1 << 25;

/// The coarse space of [`super::precond::SparsePreconditioner`] (module docs).
#[derive(Debug, Clone)]
pub struct CoarseSpace {
    /// `3|V|`.
    nd: usize,
    /// For each position column: (component, axis, 1/√|c|); then the multiplier rows.
    columns: Vec<CoarseColumn>,
    /// The vertices of every free component (indexed as in `columns`).
    members: Vec<Vec<u32>>,
    /// `A′ Z`, column-major, `dim × m`.
    az: Vec<f64>,
    /// `E⁺` (`m × m`, row-major).
    e_pinv: Vec<f64>,
    /// The smallest kept `|λ(E)|` over the largest.
    conditioning: f64,
    /// Dropped eigenvalues of `E` (rank deficiency).
    dropped: usize,
}

#[derive(Debug, Clone, Copy)]
enum CoarseColumn {
    /// Translation of free component `slot` (index into `members`) along `axis`, scaled.
    Translation { slot: usize, axis: usize, scale: f64 },
    /// The multiplier of positional row `row`.
    Multiplier { row: usize },
}

impl CoarseSpace {
    /// The components whose translations are not fixed by positional rows supported on that
    /// component alone, and the positional rows that touch one of them.
    pub fn free_components(surface: &Surface, c: &SparseMatrix) -> (Vec<u32>, Vec<usize>) {
        let vc = surface.vertex_components();
        let nc = surface.n_components();
        let mut fixed = vec![false; nc];
        let ct = c.transpose();
        let (cp, ri) = (ct.col_ptr(), ct.row_idx());
        let mut touched: Vec<Vec<u32>> = Vec::with_capacity(c.nrows());
        for r in 0..c.nrows() {
            let mut comps: Vec<u32> = (cp[r]..cp[r + 1]).map(|p| vc[ri[p] as usize / 3]).collect();
            comps.sort_unstable();
            comps.dedup();
            if comps.len() == 1 {
                fixed[comps[0] as usize] = true;
            }
            touched.push(comps);
        }
        let free: Vec<u32> = (0..nc as u32).filter(|&k| !fixed[k as usize]).collect();
        let rows = touched
            .iter()
            .enumerate()
            .filter(|(_, cs)| cs.iter().any(|k| !fixed[*k as usize]))
            .map(|(r, _)| r)
            .collect();
        (free, rows)
    }

    /// Builds the coarse space of the free components `free` and the positional rows `rows`
    /// (from [`CoarseSpace::free_components`]) for the saddle operator `op` (dimension
    /// `3|V| + k_Q`, `A₃ = A ⊗ I₃`) with the positional Jacobian `c`. `None` when there is no
    /// free component or `A′ Z` would exceed [`COARSE_MAX_ENTRIES`].
    ///
    /// # Errors
    /// Operator errors; [`RsError::InvalidOption`] for mismatched sizes.
    pub fn build(
        op: &dyn LinearOperator,
        surface: &Surface,
        c: &SparseMatrix,
        free: &[u32],
        rows: &[usize],
    ) -> Result<Option<CoarseSpace>, RsError> {
        let nv = surface.n_vertices();
        let nd = 3 * nv;
        let dim = op.dim();
        if dim != nd + c.nrows() || c.ncols() != nd {
            return Err(RsError::InvalidOption {
                name: "coarse space",
                value: format!("operator {dim}, rows {} x {}", c.nrows(), c.ncols()),
                requirement: "an operator of dimension 3|V| + k_Q",
            });
        }
        if free.is_empty() {
            return Ok(None);
        }
        let m = 3 * free.len() + rows.len();
        if dim.saturating_mul(m) > COARSE_MAX_ENTRIES {
            return Ok(None);
        }
        let vc = surface.vertex_components();
        let mut slot_of = vec![usize::MAX; surface.n_components()];
        for (k, &comp) in free.iter().enumerate() {
            slot_of[comp as usize] = k;
        }
        let mut members: Vec<Vec<u32>> = vec![Vec::new(); free.len()];
        for (v, &comp) in vc.iter().enumerate() {
            let k = slot_of[comp as usize];
            if k != usize::MAX {
                members[k].push(v as u32);
            }
        }
        // A 1_c for every free component, three components per product.
        let mut a1: Vec<Vec<f64>> = vec![Vec::new(); free.len()];
        let mut x = vec![0.0; dim];
        let mut y = vec![0.0; dim];
        for chunk in (0..free.len()).collect::<Vec<_>>().chunks(3) {
            x.iter_mut().for_each(|v| *v = 0.0);
            for (axis, &k) in chunk.iter().enumerate() {
                for &v in &members[k] {
                    x[3 * v as usize + axis] = 1.0;
                }
            }
            op.apply(&x, &mut y)?;
            for (axis, &k) in chunk.iter().enumerate() {
                a1[k] = (0..nv).map(|v| y[3 * v + axis]).collect();
            }
        }
        let mut columns = Vec::with_capacity(m);
        for (k, mem) in members.iter().enumerate() {
            let scale = 1.0 / (mem.len() as f64).sqrt();
            for axis in 0..3 {
                columns.push(CoarseColumn::Translation { slot: k, axis, scale });
            }
        }
        columns.extend(rows.iter().map(|&row| CoarseColumn::Multiplier { row }));
        let ct = c.transpose();
        let err = |e: ddg_core::SparseError| RsError::InvalidOption {
            name: "coarse space",
            value: e.to_string(),
            requirement: "matching sparse dimensions",
        };
        let mut az = vec![0.0; dim * m];
        let mut zcol = vec![0.0; nd];
        for (j, col) in columns.iter().enumerate() {
            let out = &mut az[j * dim..(j + 1) * dim];
            match *col {
                CoarseColumn::Translation { slot, axis, scale } => {
                    for v in 0..nv {
                        out[3 * v + axis] = scale * a1[slot][v];
                    }
                    zcol.iter_mut().for_each(|v| *v = 0.0);
                    for &v in &members[slot] {
                        zcol[3 * v as usize + axis] = scale;
                    }
                    c.mul_vec(&zcol, &mut out[nd..]).map_err(err)?;
                }
                CoarseColumn::Multiplier { row } => {
                    let (cp, ri, val) = (ct.col_ptr(), ct.row_idx(), ct.values());
                    for p in cp[row]..cp[row + 1] {
                        out[ri[p] as usize] = val[p];
                    }
                }
            }
        }
        let mut space = CoarseSpace {
            nd,
            columns,
            members,
            az,
            e_pinv: Vec::new(),
            conditioning: 0.0,
            dropped: 0,
        };
        // E = Zᵀ A′ Z, symmetrized.
        let mut e = Mat::<f64>::zeros(m, m);
        for j in 0..m {
            let col = &space.az[j * dim..(j + 1) * dim];
            let t = space.zt(col);
            for i in 0..m {
                e[(i, j)] = t[i];
            }
        }
        for j in 0..m {
            for i in 0..j {
                let s = 0.5 * (e[(i, j)] + e[(j, i)]);
                e[(i, j)] = s;
                e[(j, i)] = s;
            }
        }
        let evd = e.self_adjoint_eigen(Side::Lower).map_err(|e| RsError::InvalidOption {
            name: "coarse space",
            value: format!("{e:?}"),
            requirement: "a symmetric eigendecomposition of the coarse matrix",
        })?;
        let (u, s) = (evd.U(), evd.S());
        let lam: Vec<f64> = (0..m).map(|i| s[i]).collect();
        let lmax = lam.iter().fold(0.0f64, |a, l| a.max(l.abs()));
        let cut = 1e-12 * lmax;
        let mut pinv = vec![0.0; m * m];
        let mut lmin = f64::INFINITY;
        for (k, &l) in lam.iter().enumerate() {
            if l.is_nan() || l.abs() <= cut {
                space.dropped += 1;
                continue;
            }
            lmin = lmin.min(l.abs());
            for i in 0..m {
                let ui = u[(i, k)] / l;
                for jj in 0..m {
                    pinv[i * m + jj] += ui * u[(jj, k)];
                }
            }
        }
        space.e_pinv = pinv;
        space.conditioning = if lmax > 0.0 { lmin / lmax } else { 0.0 };
        Ok(Some(space))
    }

    /// `m`, the coarse dimension.
    pub fn dim(&self) -> usize {
        self.columns.len()
    }

    /// `(smallest kept |λ(E)| / largest, dropped eigenvalues)`.
    pub fn conditioning(&self) -> (f64, usize) {
        (self.conditioning, self.dropped)
    }

    /// `Zᵀ r`.
    fn zt(&self, r: &[f64]) -> Vec<f64> {
        self.columns
            .iter()
            .map(|col| match *col {
                CoarseColumn::Translation { slot, axis, scale } => {
                    scale * self.members[slot].iter().map(|&v| r[3 * v as usize + axis]).sum::<f64>()
                }
                CoarseColumn::Multiplier { row } => r[self.nd + row],
            })
            .collect()
    }

    /// `out += Z t`.
    fn add_z(&self, t: &[f64], out: &mut [f64]) {
        for (col, &tj) in self.columns.iter().zip(t) {
            match *col {
                CoarseColumn::Translation { slot, axis, scale } => {
                    for &v in &self.members[slot] {
                        out[3 * v as usize + axis] += scale * tj;
                    }
                }
                CoarseColumn::Multiplier { row } => out[self.nd + row] += tj,
            }
        }
    }

    fn e_solve(&self, t: &[f64]) -> Vec<f64> {
        let m = self.dim();
        (0..m).map(|i| (0..m).map(|j| self.e_pinv[i * m + j] * t[j]).sum()).collect()
    }

    /// `(A′Z)ᵀ z`.
    fn azt(&self, z: &[f64]) -> Vec<f64> {
        let dim = z.len();
        (0..self.dim())
            .map(|j| self.az[j * dim..(j + 1) * dim].iter().zip(z).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// The balancing preconditioner (module docs) around the fine `P⁻¹` (`fine`, applied to
    /// all right sides in one batch).
    ///
    /// # Errors
    /// The errors of `fine`.
    pub fn apply_many<F>(&self, rs: &[&[f64]], fine: F) -> Result<Vec<Vec<f64>>, RsError>
    where
        F: FnOnce(&[&[f64]]) -> Result<Vec<Vec<f64>>, RsError>,
    {
        let mut ycs = Vec::with_capacity(rs.len());
        let mut rps = Vec::with_capacity(rs.len());
        for r in rs {
            let dim = r.len();
            let yc = self.e_solve(&self.zt(r));
            let mut rp = r.to_vec();
            for (j, &yj) in yc.iter().enumerate() {
                if yj != 0.0 {
                    for (ri, a) in rp.iter_mut().zip(&self.az[j * dim..(j + 1) * dim]) {
                        *ri -= yj * a;
                    }
                }
            }
            ycs.push(yc);
            rps.push(rp);
        }
        let rpr: Vec<&[f64]> = rps.iter().map(Vec::as_slice).collect();
        let mut zs = fine(&rpr)?;
        for (z, yc) in zs.iter_mut().zip(&ycs) {
            let t2 = self.e_solve(&self.azt(z));
            let corr: Vec<f64> = yc.iter().zip(&t2).map(|(a, b)| a - b).collect();
            self.add_z(&corr, z);
        }
        Ok(zs)
    }
}
