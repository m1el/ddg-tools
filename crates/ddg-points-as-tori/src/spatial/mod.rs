//! Neighbor search: our own deterministic kd-tree (SPEC.md §3, decision D10) and the fitting
//! neighborhoods `N_k(i)` of §3.2.

mod kdtree;

pub use kdtree::{KdTree, KdTreeOptions, Neighbor};

use rayon::prelude::*;

/// Fitting neighborhoods: row `i` holds `N_k(i) = [i] ++ (the k − 1 nearest j ≠ i by
/// (‖p_j − p_i‖², j))`. With `n < k` every row has `n` entries (`k` is then `n`).
#[derive(Debug, Clone, PartialEq)]
pub struct Neighborhoods {
    /// Entries per row.
    pub k: usize,
    /// `n × k` indices, row major.
    pub indices: Vec<u32>,
}

impl Neighborhoods {
    /// Number of rows.
    pub fn len(&self) -> usize {
        self.indices.len().checked_div(self.k).unwrap_or(0)
    }

    /// True when there are no rows.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Row `i`.
    pub fn row(&self, i: usize) -> &[u32] {
        &self.indices[i * self.k..(i + 1) * self.k]
    }
}

/// Rows per parallel task in batch queries (a constant, so the work split never depends on the
/// thread count; each row is computed independently anyway).
const QUERY_CHUNK: usize = 256;

/// The fitting neighborhoods of every point of `points` (SPEC.md §3.2, decision D11): self
/// first, then the `k − 1` nearest other points by `(d², index)`. Runs on the current rayon pool.
pub fn fitting_neighborhoods(tree: &KdTree, points: &[[f64; 3]], k: usize) -> Neighborhoods {
    let n = points.len();
    let k = k.min(n).max(1);
    let mut indices = vec![0u32; n * k];
    indices.par_chunks_mut(k * QUERY_CHUNK).enumerate().for_each_init(Vec::new, |buf, (c, rows)| {
        for (r, row) in rows.chunks_mut(k).enumerate() {
            let i = (c * QUERY_CHUNK + r) as u32;
            tree.knn(points[i as usize], k, buf);
            row[0] = i;
            let mut w = 1;
            for nb in buf.iter() {
                if nb.index != i && w < k {
                    row[w] = nb.index;
                    w += 1;
                }
            }
        }
    });
    Neighborhoods { k, indices }
}

/// `k` nearest points of each query, ordered by `(d², index)`: indices and squared distances,
/// `queries.len() × min(k, n)` each, row major. Runs on the current rayon pool.
pub fn knn_batch(tree: &KdTree, queries: &[[f64; 3]], k: usize) -> (Vec<u32>, Vec<f64>) {
    let k = k.min(tree.len());
    let mut idx = vec![0u32; queries.len() * k];
    let mut d2 = vec![0f64; queries.len() * k];
    if k == 0 {
        return (idx, d2);
    }
    idx.par_chunks_mut(k * QUERY_CHUNK)
        .zip(d2.par_chunks_mut(k * QUERY_CHUNK))
        .enumerate()
        .for_each_init(Vec::new, |buf, (c, (ri, rd))| {
            for (r, (oi, od)) in ri.chunks_mut(k).zip(rd.chunks_mut(k)).enumerate() {
                tree.knn(queries[c * QUERY_CHUNK + r], k, buf);
                for (t, nb) in buf.iter().enumerate() {
                    oi[t] = nb.index;
                    od[t] = nb.d2;
                }
            }
        });
    (idx, d2)
}
