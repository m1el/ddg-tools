//! The reduced problem. With `A` the Laplacian restricted to the unknown rows (interior
//! vertices; on a closed mesh all vertices but a pinned one, whose `u` is fixed to 0 and
//! shifted afterwards), `u = P A⁻¹ (µ − Ω)` where `P` removes the `M`-weighted mean on a closed
//! mesh and is the identity on a mesh with boundary. One Cholesky factor of `A` is all the
//! linear algebra the solver needs.

use std::time::Instant;

use ddg_core::operators::{angle_defects, cotan_laplacian, lumped_mass_normalized};
use ddg_core::sparse::{Factorization, Method, Parallelism, Symbolic};
use faer::linalg::matmul::matmul;
use faer::{Accum, MatMut, MatRef, Par};
use rayon::prelude::*;

use crate::meshes::Case;

/// Per-vertex steering, as in `ConeProblem`: `cost` = `wR` (0 free site, `INFINITY` forbidden),
/// `nonneg` = `lo = 0`, `nonpos` = `hi = 0`.
#[derive(Clone, Default)]
pub struct Constraints {
    pub cost: Option<Vec<f64>>,
    pub nonneg: Option<Vec<bool>>,
    pub nonpos: Option<Vec<bool>>,
}

pub struct Prepared<'a> {
    pub case: &'a Case,
    pub closed: bool,
    pub n_unk: usize,
    /// vertex of each row: unknown rows first, then (closed) the pin
    pub rows: Vec<u32>,
    /// vertex of each slot (variables of µ): bounded: interior vertices; closed: pin first, then the rest
    pub slot_vertex: Vec<u32>,
    pub mrow: Vec<f64>,
    pub omega_rows: Vec<f64>,
    pub fac: Factorization,
    /// projected `A⁻¹ Ω` over rows
    pub w0: Vec<f64>,
    pub mw0: Vec<f64>,
    /// closed: `Σ Ω = 2πχ`
    pub sum: Option<f64>,
    /// cost weight and allowed signs (bit 0 positive, bit 1 negative) of each slot
    pub slot_w: Vec<f64>,
    pub slot_allow: Vec<u8>,
    slot_unk_tab: Vec<u32>,
    /// neighbours of each slot, as slots
    pub adj_ptr: Vec<usize>,
    pub adj_idx: Vec<u32>,
    /// factor values per row: a solve touches about this many numbers per row
    pub factor_per_row: f64,
    pub setup_seconds: f64,
    pub factor_seconds: f64,
    pub n_solves: std::sync::atomic::AtomicUsize,
}

impl<'a> Prepared<'a> {
    pub fn new(case: &'a Case) -> Result<Prepared<'a>, String> {
        Prepared::with_constraints(case, &Constraints::default())
    }

    pub fn with_constraints(case: &'a Case, cons: &Constraints) -> Result<Prepared<'a>, String> {
        let t0 = Instant::now();
        let mesh = &case.mesh;
        let nv = mesh.n_vertices();
        let closed = mesh.n_boundary_edges() == 0;
        let l = cotan_laplacian(mesh, &case.geometry).map_err(|e| e.to_string())?;
        let mass = lumped_mass_normalized(mesh, &case.geometry).map_err(|e| e.to_string())?;
        let omega = angle_defects(mesh, &case.geometry).map_err(|e| e.to_string())?;
        let cost_of = |v: usize| cons.cost.as_ref().map_or(1.0, |c| c[v]);
        let allow_of = |v: usize| -> u8 {
            let mut a = 3u8;
            if cons.nonneg.as_ref().is_some_and(|c| c[v]) { a &= !2; }
            if cons.nonpos.as_ref().is_some_and(|c| c[v]) { a &= !1; }
            a
        };
        let total_omega: f64 = omega.iter().sum();
        let (unk, rows, slot_vertex): (Vec<u32>, Vec<u32>, Vec<u32>) = if closed {
            let need = if total_omega >= 0.0 { 1u8 } else { 2u8 };
            let pin = (0..nv)
                .find(|&v| cost_of(v).is_finite() && allow_of(v) & need != 0)
                .ok_or("no vertex can carry the total curvature (all forbidden or sign-restricted)")?;
            let unk: Vec<u32> = (0..nv as u32).filter(|&v| v as usize != pin).collect();
            let mut rows = unk.clone();
            rows.push(pin as u32);
            let mut slots = vec![pin as u32];
            slots.extend(unk.iter().copied().filter(|&v| cost_of(v as usize).is_finite()));
            (unk, rows, slots)
        } else {
            let unk = mesh.interior_vertices();
            let slots: Vec<u32> = unk.iter().copied().filter(|&v| cost_of(v as usize).is_finite()).collect();
            (unk.clone(), unk, slots)
        };
        let mut unk_pos = vec![u32::MAX; nv];
        for (i, &v) in unk.iter().enumerate() {
            unk_pos[v as usize] = i as u32;
        }
        let slot_unk_tab: Vec<u32> = slot_vertex.iter().map(|&v| unk_pos[v as usize]).collect();
        let slot_w: Vec<f64> = slot_vertex.iter().map(|&v| cost_of(v as usize)).collect();
        let slot_allow: Vec<u8> = slot_vertex.iter().map(|&v| allow_of(v as usize)).collect();
        let n_unk = unk.len();
        let a = l.submatrix(&unk, &unk).map_err(|e| e.to_string())?;
        let t1 = Instant::now();
        let sym = Symbolic::analyze(&a, Method::Llt).map_err(|e| e.to_string())?;
        let factor_per_row = sym.info().factor_values.unwrap_or(0) as f64 / n_unk.max(1) as f64;
        let fac = sym.factor(&a, Parallelism::Sequential).map_err(|e| e.to_string())?;
        let factor_seconds = t1.elapsed().as_secs_f64();
        let mrow: Vec<f64> = rows.iter().map(|&v| mass[v as usize]).collect();
        let omega_rows: Vec<f64> = rows.iter().map(|&v| omega[v as usize]).collect();
        let sum = closed.then_some(total_omega);
        // slot adjacency
        let mut slot_of = vec![u32::MAX; nv];
        for (s, &v) in slot_vertex.iter().enumerate() {
            slot_of[v as usize] = s as u32;
        }
        let mut adj_ptr = vec![0usize];
        let mut adj_idx = Vec::new();
        for &v in &slot_vertex {
            for w in mesh.vertex_neighbors(v as usize) {
                if slot_of[w] != u32::MAX {
                    adj_idx.push(slot_of[w]);
                }
            }
            adj_ptr.push(adj_idx.len());
        }
        let mut p = Prepared {
            case,
            closed,
            n_unk,
            rows,
            slot_vertex,
            mrow,
            omega_rows,
            fac,
            w0: vec![],
            mw0: vec![],
            sum,
            slot_w,
            slot_allow,
            slot_unk_tab,
            adj_ptr,
            adj_idx,
            factor_per_row,
            setup_seconds: 0.0,
            factor_seconds,
            n_solves: 0.into(),
        };
        let mut rhs = p.omega_rows[..n_unk].to_vec();
        p.solve_unk(&mut rhs, 1);
        p.w0 = p.embed_project(rhs);
        p.mw0 = p.w0.iter().zip(&p.mrow).map(|(a, b)| a * b).collect();
        p.setup_seconds = t0.elapsed().as_secs_f64();
        Ok(p)
    }

    pub fn n_slots(&self) -> usize {
        self.slot_vertex.len()
    }

    /// Position among the unknown rows of a slot, `None` for the pin.
    #[inline]
    pub fn slot_unk(&self, s: usize) -> Option<usize> {
        let t = self.slot_unk_tab[s];
        (t != u32::MAX).then_some(t as usize)
    }

    #[inline]
    pub fn slot_row(&self, s: usize) -> usize {
        match self.slot_unk(s) {
            Some(i) => i,
            None => self.n_unk,
        }
    }

    pub fn solve_unk(&self, buf: &mut [f64], nrhs: usize) {
        self.n_solves.fetch_add(nrhs, std::sync::atomic::Ordering::Relaxed);
        self.fac
            .solve_many_in_place(buf, nrhs, Parallelism::Sequential)
            .expect("solve");
    }

    /// Unknown-row vector to row vector (pin gets 0), then the `M`-weighted mean removed (closed).
    pub fn embed_project(&self, mut x: Vec<f64>) -> Vec<f64> {
        if self.closed {
            x.push(0.0);
            let mean: f64 = x.iter().zip(&self.mrow).map(|(a, b)| a * b).sum();
            x.iter_mut().for_each(|v| *v -= mean);
        }
        x
    }

    /// `u` over rows for a sparse `µ` given as `(slot, value)`; one solve.
    pub fn u_of_mu(&self, mu: &[(usize, f64)]) -> Vec<f64> {
        let mut rhs: Vec<f64> = self.omega_rows[..self.n_unk].iter().map(|o| -o).collect();
        for &(s, v) in mu {
            if let Some(i) = self.slot_unk(s) {
                rhs[i] += v;
            }
        }
        self.solve_unk(&mut rhs, 1);
        self.embed_project(rhs)
    }

    /// `y = A⁻¹ (M u)` over unknown rows: the gradient of the distortion with respect to `µ`
    /// (up to the multiplier of the sum constraint); one solve.
    pub fn gradient(&self, u: &[f64]) -> Vec<f64> {
        let mut rhs: Vec<f64> = (0..self.n_unk).map(|i| self.mrow[i] * u[i]).collect();
        self.solve_unk(&mut rhs, 1);
        rhs
    }

    /// Gradient entry of a slot.
    #[inline]
    pub fn grad_of_slot(&self, y: &[f64], s: usize) -> f64 {
        self.slot_unk(s).map_or(0.0, |i| y[i])
    }

    /// Projected columns `P A⁻¹ e_s` for the given slots, column major `n_rows × k`.
    /// `pool` runs chunks of `chunk` right-hand sides on separate threads.
    pub fn columns(&self, slots: &[usize], pool: Option<&rayon::ThreadPool>, chunk: usize) -> Vec<f64> {
        let nr = self.rows.len();
        let k = slots.len();
        let mut data = vec![0.0; nr * k];
        let chunk = chunk.max(1);
        let work = |(ci, out): (usize, &mut [f64])| {
            let cnt = out.len() / nr;
            let mut buf = vec![0.0; self.n_unk * cnt];
            for t in 0..cnt {
                if let Some(i) = self.slot_unk(slots[ci * chunk + t]) {
                    buf[t * self.n_unk + i] = 1.0;
                }
            }
            self.solve_unk(&mut buf, cnt);
            for t in 0..cnt {
                let col = &mut out[t * nr..(t + 1) * nr];
                col[..self.n_unk].copy_from_slice(&buf[t * self.n_unk..(t + 1) * self.n_unk]);
                if self.closed {
                    col[self.n_unk] = 0.0;
                    let mean: f64 = col.iter().zip(&self.mrow).map(|(a, b)| a * b).sum();
                    col.iter_mut().for_each(|v| *v -= mean);
                }
            }
        };
        match pool {
            Some(pool) => pool.install(|| {
                data.par_chunks_mut(nr * chunk).enumerate().for_each(work);
            }),
            None => data.chunks_mut(nr * chunk).enumerate().for_each(work),
        }
        data
    }

    /// Gram cross terms for a new block of columns `zn` (`k` columns) against the existing blocks
    /// (each `(data, columns)`): returns `(cross, q)` with `cross` of size `(dim + k) × k`,
    /// column major, and `q_t = −z_tᵀ M w0`.
    pub fn gram_cross(
        &self,
        blocks: &[(Vec<f64>, usize)],
        zn: &[f64],
        k: usize,
        pool: Option<&rayon::ThreadPool>,
        nthreads: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        let nr = self.rows.len();
        let dim_old: usize = blocks.iter().map(|b| b.1).sum();
        let dim = dim_old + k;
        let mut cross = vec![0.0; dim * k];
        // partial products over a range of rows, all blocks sequentially
        let partial = |r0: usize, r1: usize| -> Vec<f64> {
            let m = r1 - r0;
            // mass-weighted copy of the new columns for this row range
            let mut mz = vec![0.0; m * k];
            for t in 0..k {
                for i in 0..m {
                    mz[t * m + i] = zn[t * nr + r0 + i] * self.mrow[r0 + i];
                }
            }
            let rhs = MatRef::from_column_major_slice(&mz[..], m, k);
            let mut part = vec![0.0; dim * k + k];
            for t in 0..k {
                let mut acc = 0.0;
                for i in 0..m {
                    acc += zn[t * nr + r0 + i] * self.mw0[r0 + i];
                }
                part[dim * k + t] = -acc;
            }
            let mut off = 0;
            let mut put = |zb: &[f64], kb: usize, off: usize| {
                let lhs = MatRef::from_column_major_slice(zb, nr, kb).subrows(r0, m);
                let mut out = vec![0.0; kb * k];
                matmul(
                    MatMut::from_column_major_slice_mut(&mut out, kb, k),
                    Accum::Replace,
                    lhs.transpose(),
                    rhs,
                    1.0,
                    Par::Seq,
                );
                for t in 0..k {
                    for i in 0..kb {
                        part[t * dim + off + i] = out[t * kb + i];
                    }
                }
            };
            for (zb, kb) in blocks {
                put(zb, *kb, off);
                off += kb;
            }
            put(zn, k, off);
            part
        };
        let nchunks = if pool.is_some() && nthreads > 1 { nthreads.min(nr / 2048).max(1) } else { 1 };
        let parts: Vec<Vec<f64>> = match pool {
            Some(p) if nchunks > 1 => p.install(|| {
                (0..nchunks)
                    .into_par_iter()
                    .map(|c| partial(c * nr / nchunks, (c + 1) * nr / nchunks))
                    .collect()
            }),
            _ => vec![partial(0, nr)],
        };
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
}

impl Prepared<'_> {
    /// `u = Σ µ_j z_j − w0` from stored columns. `terms` are `(block, column in block, µ)`.
    pub fn u_from_columns(
        &self,
        blocks: &[(Vec<f64>, usize)],
        terms: &[(usize, usize, f64)],
        pool: Option<&rayon::ThreadPool>,
    ) -> Vec<f64> {
        let nr = self.rows.len();
        let mut u: Vec<f64> = self.w0.iter().map(|w| -w).collect();
        let body = |(ci, out): (usize, &mut [f64])| {
            let r0 = ci * CH;
            for &(b, c, mu) in terms {
                let col = &blocks[b].0[c * nr + r0..c * nr + r0 + out.len()];
                for (o, z) in out.iter_mut().zip(col) {
                    *o += mu * z;
                }
            }
        };
        match pool {
            Some(p) => p.install(|| u.par_chunks_mut(CH).enumerate().for_each(body)),
            None => u.chunks_mut(CH).enumerate().for_each(body),
        }
        u
    }
}

const CH: usize = 4096;
