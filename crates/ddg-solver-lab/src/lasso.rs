//! The factor-once solver: working set outer loop around the dense active-set method.
//!
//! Round: dense solve on the working set W (columns `z_j = P A⁻¹ e_j` known, Gram matrix
//! `H = ZᵀMZ` dense), then the exact gradient of the distortion over all vertices by two sparse
//! solves (`u = P A⁻¹(µ − Ω)`, `y = A⁻¹ M u`), then new variables whose optimality condition
//! `|y_j + ν| ≤ λ` is violated enter W, with their columns computed by one solve each (in
//! parallel). Stops when nothing is violated by more than `tol·λ`.

use std::time::Instant;

use crate::dense::{Dense, DenseStats};
use crate::prepare::Prepared;

#[derive(Debug, Clone)]
pub struct Options {
    pub threads: usize,
    pub chunk: usize,
    /// continuation: start at `lambda · factor^steps` and divide by `factor` per stage (steps 0 = none)
    pub cont_factor: f64,
    pub cont_steps: usize,
    pub batch_min: usize,
    pub batch_frac: f64,
    pub nms_radius: usize,
    pub tol: f64,
    pub max_rounds: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            threads: 1,
            chunk: 4,
            cont_factor: 2.0,
            cont_steps: 0,
            batch_min: 16,
            batch_frac: 0.5,
            nms_radius: 3,
            tol: 1e-10,
            max_rounds: 500,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Stats {
    pub rounds: usize,
    pub stages: usize,
    pub column_solves: usize,
    pub gradient_solves: usize,
    pub total_solves: usize,
    pub peak_active: usize,
    pub final_active: usize,
    pub working_set: usize,
    pub t_setup: f64,
    pub t_columns: f64,
    pub t_gram: f64,
    pub t_dense: f64,
    pub t_gradient: f64,
    pub t_total: f64,
    pub dense: DenseStats,
    pub converged: bool,
    pub max_excess: f64,
}

pub struct Solution {
    /// µ per mesh vertex
    pub mu: Vec<f64>,
    /// u per mesh vertex
    pub u: Vec<f64>,
    /// adjoint per mesh vertex
    pub phi: Vec<f64>,
    pub objective: f64,
    pub distortion: f64,
    pub stats: Stats,
}

pub fn solve(prep: &Prepared, lambda: f64, opts: &Options) -> Solution {
    let t_start = Instant::now();
    let pool = (opts.threads > 1).then(|| {
        rayon::ThreadPoolBuilder::new().num_threads(opts.threads).build().expect("pool")
    });
    let pool = pool.as_ref();
    let ns = prep.n_slots();
    let nr = prep.rows.len();
    let solves_before = prep.n_solves.load(std::sync::atomic::Ordering::Relaxed);
    let mut st = Stats::default();
    st.t_setup = prep.setup_seconds;

    let mut dense = Dense::new(lambda, prep.sum);
    let mut w_slots: Vec<usize> = Vec::new();
    let mut in_w = vec![usize::MAX; ns];
    let mut blocks: Vec<(Vec<f64>, usize)> = Vec::new();
    if prep.closed {
        w_slots.push(0);
        in_w[0] = 0;
        blocks.push((vec![0.0; nr], 1));
    }

    // schedule of lambdas
    let mut schedule: Vec<f64> = Vec::new();
    for k in (0..=opts.cont_steps).rev() {
        schedule.push(lambda * opts.cont_factor.powi(k as i32));
    }

    let mut nu_final = 0.0;
    let mut y_final: Vec<f64> = Vec::new();
    let mut converged = false;
    let mut round_total = 0;
    'stages: for (si, &lam) in schedule.iter().enumerate() {
        dense.lambda = lam;
        st.stages += 1;
        loop {
            round_total += 1;
            st.rounds = round_total;
            if round_total > opts.max_rounds {
                break 'stages;
            }
            let t = Instant::now();
            let ok = dense.solve(opts.tol * 0.1, 20 * (dense.dim() + 100));
            st.t_dense += t.elapsed().as_secs_f64();
            let _ = ok;
            // exact gradient
            let t = Instant::now();
            let mu_sparse: Vec<(usize, f64)> =
                dense.active().into_iter().map(|j| (w_slots[j], dense.mu[j])).collect();
            let u = prep.u_of_mu(&mu_sparse);
            let y = prep.gradient(&u);
            st.t_gradient += t.elapsed().as_secs_f64();
            let nu = dense.nu;
            // violators
            let mut cand: Vec<(usize, f64)> = Vec::new();
            let mut max_in_w = 0.0f64;
            let mut max_all = 0.0f64;
            for s in 0..ns {
                let in_active = in_w[s] != usize::MAX && dense.is_active(in_w[s]);
                if in_active {
                    continue;
                }
                let ex = (prep.grad_of_slot(&y, s) + nu).abs() - lam;
                max_all = max_all.max(ex);
                if ex > opts.tol * lam {
                    if in_w[s] == usize::MAX {
                        cand.push((s, ex));
                    } else {
                        max_in_w = max_in_w.max(ex);
                    }
                }
            }
            st.peak_active = st.peak_active.max(dense.n_active());
            let last = si + 1 == schedule.len();
            if cand.is_empty() {
                if last {
                    st.max_excess = max_all.max(0.0) / lam;
                    nu_final = nu;
                    y_final = y;
                    converged = max_in_w <= 1e-6 * lam;
                    break 'stages;
                }
                break;
            }
            // choose a batch
            cand.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let k = ((opts.batch_frac * dense.n_active() as f64) as usize).max(opts.batch_min);
            let mut chosen: Vec<usize> = Vec::new();
            if opts.nms_radius == 0 {
                chosen.extend(cand.iter().take(k).map(|c| c.0));
            } else {
                let mut blocked = vec![false; ns];
                for &(s, _) in &cand {
                    if chosen.len() >= k {
                        break;
                    }
                    if blocked[s] {
                        continue;
                    }
                    chosen.push(s);
                    // block the radius-ball around s
                    let mut frontier = vec![s];
                    blocked[s] = true;
                    for _ in 0..opts.nms_radius {
                        let mut next = Vec::new();
                        for &f in &frontier {
                            for &w in &prep.adj_idx[prep.adj_ptr[f]..prep.adj_ptr[f + 1]] {
                                if !blocked[w as usize] {
                                    blocked[w as usize] = true;
                                    next.push(w as usize);
                                }
                            }
                        }
                        frontier = next;
                    }
                }
            }
            // columns
            let t = Instant::now();
            let zn = prep.columns(&chosen, pool, opts.chunk);
            st.t_columns += t.elapsed().as_secs_f64();
            st.column_solves += chosen.len();
            let t = Instant::now();
            let (cross, q) = prep.gram_cross(&blocks, &zn, chosen.len(), pool, opts.threads);
            st.t_gram += t.elapsed().as_secs_f64();
            dense.extend(&cross, chosen.len(), &q);
            for &s in &chosen {
                in_w[s] = w_slots.len();
                w_slots.push(s);
            }
            blocks.push((zn, chosen.len()));
        }
    }
    st.converged = converged;
    st.working_set = w_slots.len();
    st.final_active = dense.n_active();
    st.dense = dense.stats.clone();

    // final answer
    let mu_sparse: Vec<(usize, f64)> =
        dense.active().into_iter().map(|j| (w_slots[j], dense.mu[j])).collect();
    let u_rows = prep.u_of_mu(&mu_sparse);
    let nv = prep.case.mesh.n_vertices();
    let mut mu = vec![0.0; nv];
    for &(s, v) in &mu_sparse {
        mu[prep.slot_vertex[s] as usize] = v;
    }
    let mut u = vec![0.0; nv];
    for (r, &v) in prep.rows.iter().enumerate() {
        u[v as usize] = u_rows[r];
    }
    let mut phi = vec![0.0; nv];
    if !y_final.is_empty() {
        for (i, y) in y_final.iter().enumerate() {
            phi[prep.rows[i] as usize] = -y - nu_final;
        }
        if prep.closed {
            phi[prep.rows[prep.n_unk] as usize] = -nu_final;
        }
    }
    let distortion: f64 = 0.5 * u_rows.iter().zip(&prep.mrow).map(|(u, m)| m * u * u).sum::<f64>();
    let cost: f64 = lambda * mu.iter().map(|m| m.abs()).sum::<f64>();
    st.t_total = t_start.elapsed().as_secs_f64() + prep.setup_seconds;
    st.total_solves = prep.n_solves.load(std::sync::atomic::Ordering::Relaxed) - solves_before;
    st.gradient_solves = st.total_solves - st.column_solves;
    Solution { mu, u, phi, objective: distortion + cost, distortion, stats: st }
}
