//! The smallest eigenpair of `A x = λ M x` (`A` symmetric positive semidefinite, `M` symmetric
//! positive definite) by block inverse iteration with Rayleigh–Ritz.
//!
//! GODF Alg. 2 is inverse iteration with one vector. It converges with the ratio of the two
//! smallest eigenvalues, which is close to one when they nearly coincide (a round sphere:
//! the smoothest 4-direction fields form an almost degenerate cluster). A block of `p`
//! vectors converges with the ratio of the smallest and the `(p+1)`-th, and its Ritz vectors
//! are the best approximations within the block; with `p = 1` this is Alg. 2.

use ddg_core::sparse::Factorization;
use ddg_core::{SparseError, SparseMatrix};
use ddg_core::generators::SplitMix64;

use super::linsolve::{dot, mul, norm, solve};

pub(crate) struct Eigenpair {
    pub vector: Vec<f64>,
    pub value: f64,
    pub residual: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// `J x`: multiplication by `i` in the stored complex form.
fn times_i(x: &[f64]) -> Vec<f64> {
    x.chunks_exact(2).flat_map(|p| [-p[1], p[0]]).collect()
}

/// Orthonormalizes the columns of `y` in the `M` inner product (modified Gram–Schmidt, twice);
/// drops columns that become dependent.
fn m_orthonormalize(m: &SparseMatrix, y: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(y.len());
    let mut m_out: Vec<Vec<f64>> = Vec::with_capacity(y.len());
    for mut v in y {
        let original = dot(&v, &mul(m, &v)).sqrt();
        for _ in 0..2 {
            for (q, mq) in out.iter().zip(&m_out) {
                let c = dot(&v, mq);
                for (x, qx) in v.iter_mut().zip(q) {
                    *x -= c * qx;
                }
            }
        }
        let mv = mul(m, &v);
        let len = dot(&v, &mv).sqrt();
        if !(len > 1e-10 * original) {
            continue;
        }
        out.push(v.iter().map(|x| x / len).collect());
        m_out.push(mv.iter().map(|x| x / len).collect());
    }
    out
}

/// Eigen-decomposition of a small symmetric matrix by cyclic Jacobi rotations: eigenvalues
/// ascending, eigenvectors as columns (`vectors[i][k]` is entry `i` of vector `k`).
pub(crate) fn jacobi(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let scale: f64 = a.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum::<f64>()
            .sqrt();
        if off <= 1e-15 * scale.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vp, vq) = (row[p], row[q]);
                    row[p] = c * vp - s * vq;
                    row[q] = s * vp + c * vq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let values = order.iter().map(|&i| a[i][i]).collect();
    let vectors = (0..n)
        .map(|r| order.iter().map(|&k| v[r][k]).collect())
        .collect();
    (values, vectors)
}

/// Block inverse iteration with `block` complex start vectors (each with its multiple by
/// `i`, so the block is a complex subspace), Rayleigh–Ritz after every step.
#[allow(clippy::too_many_arguments)]
pub(crate) fn smallest(
    a: &SparseMatrix,
    m: &SparseMatrix,
    factor: &Factorization,
    block: usize,
    max_iterations: usize,
    tolerance: f64,
    eps: f64,
    seed: u64,
) -> Result<Eigenpair, SparseError> {
    let n = a.nrows();
    let mut rng = SplitMix64::new(seed);
    let mut y: Vec<Vec<f64>> = Vec::new();
    for _ in 0..block.max(1) {
        let x: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
        let ix = times_i(&x);
        y.push(x);
        y.push(ix);
    }
    let mut iterations = 0;
    loop {
        let q = m_orthonormalize(m, y);
        let aq: Vec<Vec<f64>> = q.iter().map(|x| mul(a, x)).collect();
        let p = q.len();
        let b: Vec<Vec<f64>> = (0..p)
            .map(|i| (0..p).map(|j| 0.5 * (dot(&q[i], &aq[j]) + dot(&q[j], &aq[i]))).collect())
            .collect();
        let (values, vectors) = jacobi(b);
        // Ritz vectors, ascending.
        let ritz: Vec<Vec<f64>> = (0..p)
            .map(|k| {
                let mut x = vec![0.0; n];
                for (i, qi) in q.iter().enumerate() {
                    let c = vectors[i][k];
                    for (xi, qv) in x.iter_mut().zip(qi) {
                        *xi += c * qv;
                    }
                }
                x
            })
            .collect();
        let x = &ritz[0];
        let ax = mul(a, x);
        let mx = mul(m, x);
        let value = dot(x, &ax) / dot(x, &mx);
        let r: Vec<f64> = ax.iter().zip(&mx).map(|(u, v)| u - value * v).collect();
        let residual =
            norm(&r) / (norm(&ax) + (value.abs() + eps) * norm(&mx)).max(f64::MIN_POSITIVE);
        let _ = values;
        if residual <= tolerance || iterations >= max_iterations {
            return Ok(Eigenpair {
                vector: x.clone(),
                value,
                residual,
                iterations,
                converged: residual <= tolerance,
            });
        }
        y = Vec::with_capacity(p);
        for x in &ritz {
            let mut z = mul(m, x);
            solve(factor, &mut z)?;
            y.push(z);
        }
        iterations += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_diagonalizes() {
        let a = vec![
            vec![4.0, 1.0, 0.5],
            vec![1.0, 3.0, 0.2],
            vec![0.5, 0.2, 1.0],
        ];
        let (values, vectors) = jacobi(a.clone());
        for k in 0..3 {
            for i in 0..3 {
                let av: f64 = (0..3).map(|j| a[i][j] * vectors[j][k]).sum();
                assert!((av - values[k] * vectors[i][k]).abs() < 1e-13);
            }
        }
        assert!(values[0] <= values[1] && values[1] <= values[2]);
    }
}
