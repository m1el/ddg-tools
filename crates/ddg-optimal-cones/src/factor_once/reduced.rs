//! The reduced problem: everything the factor-once method needs from one Cholesky factor of
//! the Laplacian block.
//!
//! With `A` the Laplacian block of the unknowns (`L_II`, or all of `L` when every vertex is an
//! unknown), `Ω` the angle defects, `W = diag(w)`, `w_i = wE_i m_i`, and the cone angles of the
//! fixed sites (`µ_F`, prescribed or forbidden) moved to the right-hand side, the state equation
//! is `A u = E µ_S − (Ω − µ_F)` with `E` the columns of the free sites `S`.
//!
//! * `A` positive definite: `u = Z µ_S − w0` with `Z = A⁻¹ E` (columns `z_j = A⁻¹ e_j`) and
//!   `w0 = A⁻¹ (Ω − µ_F)`.
//! * `A` singular (annihilates constants): one unknown, the *pin* `p ∈ S`, loses its equation
//!   and its `u` is set to zero before shifting; then `A_p` (`A` without row and column `p`) is
//!   positive definite, the columns are `z_j = P A_p⁻¹ e_j` (zero at `p`) and `w0 = P A_p⁻¹
//!   (Ω − µ_F)_p̂`, with the projection `P x = x − (Σ w_i x_i / Σ w_i)·1`, which chooses the
//!   additive constant of `u` that minimizes the distortion. The pin's equation holds by the sum
//!   constraint `Σ_S µ = Σ (Ω − µ_F)`; `z_p = 0`.
//!
//! The distortion `½ uᵀ W u` is then a quadratic in `µ_S`.

use std::sync::atomic::{AtomicUsize, Ordering};

use ddg_core::SparseMatrix;
use ddg_core::sparse::{Factorization, Method, Parallelism, RefineOptions, Symbolic, ThreadPool};
use faer::linalg::matmul::matmul;
use faer::{Accum, MatMut, MatRef, Par};

use super::{FallbackReason, Unsupported};
use crate::problem::Discrete;

pub(super) struct Reduced<'a> {
    pub d: &'a Discrete,
    pub n: usize,
    /// The pin (unknown index) when `A` is singular.
    pub pin: Option<usize>,
    /// Unknown index of each row of the factored matrix (all unknowns but the pin).
    pub rows: Vec<u32>,
    /// Row of each unknown, `u32::MAX` for the pin.
    pub row_of: Vec<u32>,
    pub matrix: SparseMatrix,
    pub fac: Factorization,
    /// The free sites (unknown indices); with a pin, the pin first.
    pub slots: Vec<u32>,
    /// `Σ_S µ = c` on a singular `A`.
    pub sum: Option<f64>,
    pub w_sum: f64,
    /// `w0` (projected), `W w0`.
    pub w0: Vec<f64>,
    pub ww0: Vec<f64>,
    /// The cone angles of the fixed sites, zero elsewhere.
    pub mu_fixed: Vec<f64>,
    /// Factor entries per row: one solve touches about this many numbers per row.
    pub factor_per_row: f64,
    pub solves: AtomicUsize,
}

impl<'a> Reduced<'a> {
    /// Whether the problem is within what the method supports; on a singular `A` the pin and
    /// the required sum of the free cone angles.
    pub fn check(d: &Discrete, lambda: f64) -> Result<Option<(usize, f64)>, Unsupported> {
        let n = d.n;
        if n == 0 {
            return Err(Unsupported::NoUnknowns);
        }
        if lambda.is_nan() || lambda <= 0.0 {
            return Err(Unsupported::ZeroLambda);
        }
        if let Some(i) = d.w.iter().position(|&w| w.is_nan() || w <= 0.0) {
            return Err(Unsupported::ZeroDistortionWeight {
                vertex: d.unknowns[i],
            });
        }
        if !d.closed {
            return Ok(None);
        }
        let c: f64 = d
            .omega
            .iter()
            .zip(&d.sites)
            .map(|(o, s)| o - if s.fixed { s.lo } else { 0.0 })
            .sum();
        let pin = (0..n).find(|&i| {
            let s = &d.sites[i];
            !s.fixed && s.lo <= c && c <= s.hi
        });
        match pin {
            Some(p) => Ok(Some((p, c))),
            None => Err(Unsupported::NoFeasiblePin),
        }
    }

    /// Checks that the problem is within what the method supports, factors, and computes `w0`.
    ///
    /// [formula: reduced-lasso-problem]
    /// [formula: reduced-sum-constraint-pin]
    pub fn new(d: &'a Discrete, lambda: f64) -> Result<Reduced<'a>, FallbackReason> {
        let n = d.n;
        let pin_sum = Self::check(d, lambda)?;
        let (pin, sum) = match pin_sum {
            Some((p, c)) => (Some(p), Some(c)),
            None => (None, None),
        };
        let mu_fixed: Vec<f64> = d
            .sites
            .iter()
            .map(|s| if s.fixed { s.lo } else { 0.0 })
            .collect();
        let rhs0: Vec<f64> = d.omega.iter().zip(&mu_fixed).map(|(o, m)| o - m).collect();
        let rows: Vec<u32> = (0..n as u32).filter(|&i| Some(i as usize) != pin).collect();
        let mut row_of = vec![u32::MAX; n];
        for (r, &i) in rows.iter().enumerate() {
            row_of[i as usize] = r as u32;
        }
        let matrix = if pin.is_some() {
            d.a.submatrix(&rows, &rows).map_err(numeric)?
        } else {
            d.a.clone()
        };
        let symbolic = Symbolic::analyze(&matrix, Method::Llt).map_err(numeric)?;
        let factor_per_row =
            symbolic.info().factor_values.unwrap_or(0) as f64 / rows.len().max(1) as f64;
        let fac = symbolic
            .factor(&matrix, Parallelism::Sequential)
            .map_err(numeric)?;
        let mut slots: Vec<u32> = Vec::new();
        if let Some(p) = pin {
            slots.push(p as u32);
        }
        slots.extend(
            (0..n as u32).filter(|&i| !d.sites[i as usize].fixed && Some(i as usize) != pin),
        );
        let w_sum = d.w.iter().sum();
        let mut red = Reduced {
            d,
            n,
            pin,
            rows,
            row_of,
            matrix,
            fac,
            slots,
            sum,
            w_sum,
            w0: Vec::new(),
            ww0: Vec::new(),
            mu_fixed,
            factor_per_row,
            solves: AtomicUsize::new(0),
        };
        let mut x = red.gather(&rhs0);
        red.solve_rows(&mut x, 1).map_err(numeric)?;
        let mut w0 = red.scatter(&x);
        red.project(&mut w0);
        red.ww0 = w0.iter().zip(&d.w).map(|(a, b)| a * b).collect();
        red.w0 = w0;
        Ok(red)
    }

    /// Number of free sites (including the pin).
    pub fn n_slots(&self) -> usize {
        self.slots.len()
    }

    /// The entries of an unknown-indexed vector that are rows of the factored matrix.
    pub fn gather(&self, x: &[f64]) -> Vec<f64> {
        self.rows.iter().map(|&i| x[i as usize]).collect()
    }

    /// A row vector as an unknown-indexed one (zero at the pin).
    pub fn scatter(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.n];
        for (r, &i) in self.rows.iter().enumerate() {
            out[i as usize] = x[r];
        }
        out
    }

    /// Removes the `w`-weighted mean (singular `A` only): the additive constant of `u` that
    /// minimizes the distortion.
    ///
    /// [formula: reduced-sum-constraint-pin]
    pub fn project(&self, x: &mut [f64]) {
        if self.pin.is_some() {
            let mean: f64 = x.iter().zip(&self.d.w).map(|(a, b)| a * b).sum::<f64>() / self.w_sum;
            x.iter_mut().for_each(|v| *v -= mean);
        }
    }

    pub fn solve_rows(
        &self,
        buf: &mut [f64],
        nrhs: usize,
    ) -> Result<(), ddg_core::sparse::SparseError> {
        self.solves.fetch_add(nrhs, Ordering::Relaxed);
        self.fac
            .solve_many_in_place(buf, nrhs, Parallelism::Sequential)
    }

    /// The accurate solve `A x = b` with iterative refinement (final answers only).
    pub fn solve_refined(&self, b: &[f64]) -> Result<Vec<f64>, ddg_core::sparse::SparseError> {
        self.solves.fetch_add(1, Ordering::Relaxed);
        let mut x = vec![0.0; b.len()];
        self.fac.solve_refined(
            &self.matrix,
            b,
            &mut x,
            RefineOptions {
                max_steps: 2,
                tolerance: 1e-15,
            },
            Parallelism::Sequential,
        )?;
        Ok(x)
    }

    /// `u = P A⁻¹ (E µ_S − Ω + µ_F)` for the full cone vector `mu` over the unknowns (fixed
    /// sites at their angles), with refinement.
    pub fn u_of_mu(&self, mu: &[f64]) -> Result<Vec<f64>, ddg_core::sparse::SparseError> {
        let rhs: Vec<f64> = mu.iter().zip(&self.d.omega).map(|(m, o)| m - o).collect();
        let x = self.solve_refined(&self.gather(&rhs))?;
        let mut u = self.scatter(&x);
        self.project(&mut u);
        Ok(u)
    }

    /// `y = A⁻¹ (W u)` as a row vector (zero at the pin): the gradient of the distortion with
    /// respect to the cone angles of the free sites; one solve.
    ///
    /// [formula: reduced-optimality-conditions]
    pub fn gradient(
        &self,
        u: &[f64],
        refined: bool,
    ) -> Result<Vec<f64>, ddg_core::sparse::SparseError> {
        let wu: Vec<f64> = u.iter().zip(&self.d.w).map(|(a, b)| a * b).collect();
        let b = self.gather(&wu);
        if refined {
            self.solve_refined(&b)
        } else {
            let mut b = b;
            self.solve_rows(&mut b, 1)?;
            Ok(b)
        }
    }

    /// Gradient entry of a slot from the row vector `y`.
    #[inline]
    pub fn grad_of_slot(&self, y: &[f64], slot: usize) -> f64 {
        let r = self.row_of[self.slots[slot] as usize];
        if r == u32::MAX { 0.0 } else { y[r as usize] }
    }

    /// The projected columns `z_j` of the given slots (none of them the pin), column major,
    /// `n × k`. Chunks of `chunk` right-hand sides are solved together, on the pool if any.
    pub fn columns(
        &self,
        slots: &[usize],
        pool: Option<&ThreadPool>,
        chunk: usize,
    ) -> Result<Vec<f64>, ddg_core::sparse::SparseError> {
        let n = self.n;
        let nr = self.rows.len();
        let chunk = chunk.max(1);
        let nchunks = slots.len().div_ceil(chunk);
        let work = |ci: usize| -> Result<Vec<f64>, ddg_core::sparse::SparseError> {
            let part = &slots[ci * chunk..((ci + 1) * chunk).min(slots.len())];
            let cnt = part.len();
            let mut buf = vec![0.0; nr * cnt];
            for (t, &s) in part.iter().enumerate() {
                buf[t * nr + self.row_of[self.slots[s] as usize] as usize] = 1.0;
            }
            self.solve_rows(&mut buf, cnt)?;
            let mut out = vec![0.0; n * cnt];
            for t in 0..cnt {
                let col = &mut out[t * n..(t + 1) * n];
                for (r, &i) in self.rows.iter().enumerate() {
                    col[i as usize] = buf[t * nr + r];
                }
                self.project(col);
            }
            Ok(out)
        };
        let parts: Vec<Result<Vec<f64>, _>> = match pool {
            Some(p) if nchunks > 1 => p.map(nchunks, work),
            _ => (0..nchunks).map(work).collect(),
        };
        let mut data = Vec::with_capacity(n * slots.len());
        for part in parts {
            data.extend_from_slice(&part?);
        }
        Ok(data)
    }

    /// Gram cross terms of a new block of `k` columns `zn` against the stored blocks (each
    /// `(data, columns)`, the first dense index being `offset`): returns `(cross, q)` with
    /// `cross` of size `(offset + stored + k) × k`, column major, rows indexed like the dense
    /// problem (the pin's row is zero), and `q_s = −z_sᵀ W w0`.
    ///
    /// [formula: reduced-gram-matrix]
    pub fn gram_cross(
        &self,
        blocks: &[(Vec<f64>, usize)],
        offset: usize,
        zn: &[f64],
        k: usize,
        pool: Option<&ThreadPool>,
        nthreads: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let stored: usize = blocks.iter().map(|b| b.1).sum();
        let dim = offset + stored + k;
        let partial = |r0: usize, r1: usize| -> Vec<f64> {
            let m = r1 - r0;
            let mut mz = vec![0.0; m * k];
            for s in 0..k {
                for i in 0..m {
                    mz[s * m + i] = zn[s * n + r0 + i] * self.d.w[r0 + i];
                }
            }
            let rhs = MatRef::from_column_major_slice(&mz[..], m, k);
            let mut part = vec![0.0; dim * k + k];
            for s in 0..k {
                let mut acc = 0.0;
                for i in 0..m {
                    acc += zn[s * n + r0 + i] * self.ww0[r0 + i];
                }
                part[dim * k + s] = -acc;
            }
            let mut put = |zb: &[f64], kb: usize, off: usize| {
                let lhs = MatRef::from_column_major_slice(zb, n, kb).subrows(r0, m);
                let mut out = vec![0.0; kb * k];
                matmul(
                    MatMut::from_column_major_slice_mut(&mut out, kb, k),
                    Accum::Replace,
                    lhs.transpose(),
                    rhs,
                    1.0,
                    Par::Seq,
                );
                for s in 0..k {
                    for i in 0..kb {
                        part[s * dim + off + i] = out[s * kb + i];
                    }
                }
            };
            let mut off = offset;
            for (zb, kb) in blocks {
                put(zb, *kb, off);
                off += kb;
            }
            put(zn, k, off);
            part
        };
        let nchunks = match pool {
            Some(_) if nthreads > 1 => nthreads.min(n / 2048).max(1),
            _ => 1,
        };
        let parts: Vec<Vec<f64>> = match pool {
            Some(p) if nchunks > 1 => {
                p.map(nchunks, |c| partial(c * n / nchunks, (c + 1) * n / nchunks))
            }
            _ => vec![partial(0, n)],
        };
        let mut cross = vec![0.0; dim * k];
        let mut q = vec![0.0; k];
        for part in parts {
            for (c, p) in cross.iter_mut().zip(&part[..dim * k]) {
                *c += p;
            }
            for (c, p) in q.iter_mut().zip(&part[dim * k..]) {
                *c += p;
            }
        }
        (cross, q)
    }

    /// `u = Σ µ_j z_j − w0` from stored columns. `terms` are `(block, column in block, µ)`.
    pub fn u_from_columns(
        &self,
        blocks: &[(Vec<f64>, usize)],
        terms: &[(usize, usize, f64)],
        pool: Option<&ThreadPool>,
    ) -> Vec<f64> {
        const CH: usize = 4096;
        let n = self.n;
        let nch = n.div_ceil(CH);
        let body = |ci: usize| -> Vec<f64> {
            let r0 = ci * CH;
            let r1 = (r0 + CH).min(n);
            let mut out: Vec<f64> = self.w0[r0..r1].iter().map(|w| -w).collect();
            for &(b, c, mu) in terms {
                let col = &blocks[b].0[c * n + r0..c * n + r1];
                for (o, z) in out.iter_mut().zip(col) {
                    *o += mu * z;
                }
            }
            out
        };
        let parts: Vec<Vec<f64>> = match pool {
            Some(p) if nch > 1 => p.map(nch, body),
            _ => (0..nch).map(body).collect(),
        };
        parts.concat()
    }
}

fn numeric(e: ddg_core::sparse::SparseError) -> FallbackReason {
    FallbackReason::FactorizationFailed {
        message: e.to_string(),
    }
}
