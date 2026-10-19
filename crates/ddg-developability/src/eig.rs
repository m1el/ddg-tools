//! Accurate eigensolvers for symmetric 3x3 and 2x2 real matrices.
//!
//! The developability energy `E_λ` is the smallest eigenvalue `λ` of
//! `A = Σ θ_k N_k N_kᵀ` (symmetric positive semidefinite, often nearly rank 2), and its gradient
//! uses the matching unit eigenvector `x` (`∇_A λ = x xᵀ`). Stein, Grinspun and Crane (2018,
//! Sec. 4.3) report that an accurate 3x3 eigensolver (Kopp 2008) is needed for convergence. What
//! matters is therefore the *absolute* accuracy of the smallest eigenvalue relative to `‖A‖`,
//! also when `λ ≪ ‖A‖`, and the eigenvector accuracy when the gap to the next eigenvalue is small.
//!
//! # Algorithm
//!
//! [`sym_eig3`] is cyclic Jacobi (Golub & Van Loan, *Matrix Computations*, 4th ed., Sec. 8.5;
//! Rutishauser's stable rotation formulas; Demmel & Veselić 1992 for the accuracy of Jacobi).
//! Kopp (2008, "Efficient numerical diagonalization of hermitian 3x3 matrices", arXiv
//! physics/0610206) compares analytic and iterative 3x3 methods and hybrids of them; nothing here
//! is taken from that paper or its code. A 3x3 Jacobi converges in 3-5 sweeps and costs only a few
//! hundred flops, so we use Jacobi alone: there is no accuracy check and no fallback to get wrong.
//!
//! Details: the matrix is first scaled by a power of two (exact) so its largest entry is in
//! `[1, 2)`; this makes overflow and underflow impossible, also for entries near `1e±300`.
//! Rotations update the diagonal as `a_pp -= t a_pq`, which keeps small eigenvalues accurate;
//! an off-diagonal entry is set to zero once it cannot change either of its diagonal neighbours
//! (`|a_pp| + 100 |a_pq| == |a_pp|`, the Numerical Recipes test) and iteration stops when all
//! three vanish (at most 50 sweeps, never reached in practice).
//!
//! # Accuracy
//!
//! Eigenvalues: `|λ̂_i − λ_i| ≤ c ε ‖A‖` with `c` a small constant (Jacobi is in practice even
//! relatively accurate for well-scaled definite matrices). Eigenvectors: the computed `x̂_i`
//! satisfies `‖A x̂_i − λ̂_i x̂_i‖ ≤ c ε ‖A‖` and the three vectors are orthonormal to a few
//! `ε`; but the *direction* of an eigenvector is only determined up to the gap: by the
//! Davis–Kahan theorem the angle to the exact eigenvector is about `ε ‖A‖ / gap_i`, where
//! `gap_i` is the distance of `λ_i` to the other eigenvalues. For `E_λ` this means the gradient
//! `x xᵀ` is only meaningful when the two smallest eigenvalues are separated by well more than
//! `ε ‖A‖`. The integration tests check these bounds against a double-double reference.
//!
//! # Input conventions
//!
//! Only the upper triangle (`a[i][j]` with `i ≤ j`) is read; the matrix is treated as symmetric.
//! Non-finite entries in the upper triangle give [`EigError::NonFinite`]; nothing panics.
//! Eigenvalues are ascending, the eigenvectors form a right-handed orthonormal basis
//! (determinant `+1`). For repeated eigenvalues any orthonormal basis of the eigenspace may be
//! returned. The zero matrix returns zeros and the identity basis.
//!
//! # Speed
//!
//! `cargo bench -p ddg-developability --bench eig` (release, one thread): `sym_eig3` about 260 ns
//! per call on matrices of the form `Σ θ N Nᵀ` (6 nearly planar normals), `sym_eig2` about 50 ns.
//! A dedicated smallest-eigenpair routine was not added: the full Jacobi is already this cheap
//! and a cheaper analytic route would give up the accuracy this module exists for.

use std::fmt;

/// Error returned by the eigensolvers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EigError {
    /// An entry of the (upper triangle of the) input is NaN or infinite.
    NonFinite,
}

impl fmt::Display for EigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EigError::NonFinite => write!(f, "matrix has a non-finite (NaN or infinite) entry"),
        }
    }
}

impl std::error::Error for EigError {}

/// Eigendecomposition of a symmetric 3x3 matrix: `A = Σ values[k] · vectors[k] vectors[k]ᵀ`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymEig3 {
    /// Eigenvalues, ascending.
    pub values: [f64; 3],
    /// `vectors[k]` is the unit eigenvector of `values[k]`; the three form a right-handed
    /// orthonormal basis (`vectors[0] × vectors[1] = vectors[2]`).
    pub vectors: [[f64; 3]; 3],
}

/// Eigendecomposition of a symmetric 2x2 matrix: `A = Σ values[k] · vectors[k] vectors[k]ᵀ`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymEig2 {
    /// Eigenvalues, ascending.
    pub values: [f64; 2],
    /// `vectors[k]` is the unit eigenvector of `values[k]`; `vectors[1]` is `vectors[0]`
    /// rotated by +90 degrees (right-handed).
    pub vectors: [[f64; 2]; 2],
}

/// `x · 2^k`, exact (barring subnormal results), for any `k` in the range of `f64` exponents.
fn ldexp(x: f64, k: i32) -> f64 {
    let h = k / 2;
    x * 2f64.powi(h) * 2f64.powi(k - h)
}

/// Exponent `k` with `2^k <= m < 2^(k+1)` for finite `m > 0` (also subnormal).
fn exponent(m: f64) -> i32 {
    let bits = m.to_bits();
    let e = ((bits >> 52) & 0x7ff) as i32;
    if e == 0 {
        // subnormal: scale up first
        exponent(m * 2f64.powi(100)) - 100
    } else {
        e - 1023
    }
}

/// Jacobi rotation parameters `(t, c)` zeroing `a_pq`, given `a_pp`, `a_qq`, `a_pq != 0`.
/// Rutishauser's formulas: `t = tan φ` is the smaller root of `t² + 2θt − 1 = 0`.
fn rotation(app: f64, aqq: f64, apq: f64) -> (f64, f64) {
    let theta = (aqq - app) / (2.0 * apq);
    let t = if theta.abs() > 1e150 {
        0.5 / theta
    } else {
        let s = if theta >= 0.0 { 1.0 } else { -1.0 };
        s / (theta.abs() + (theta * theta + 1.0).sqrt())
    };
    let c = 1.0 / (t * t + 1.0).sqrt();
    (t, c)
}

/// Eigendecomposition of the symmetric matrix `a` (upper triangle read), eigenvalues ascending.
///
/// See the [module documentation](self) for the algorithm, accuracy and conventions.
///
/// # Errors
/// [`EigError::NonFinite`] if an entry of the upper triangle is NaN or infinite.
pub fn sym_eig3(a: [[f64; 3]; 3]) -> Result<SymEig3, EigError> {
    let mut m = [
        [a[0][0], a[0][1], a[0][2]],
        [a[0][1], a[1][1], a[1][2]],
        [a[0][2], a[1][2], a[2][2]],
    ];
    let mut mx = 0.0f64;
    for row in &m {
        for &x in row {
            if !x.is_finite() {
                return Err(EigError::NonFinite);
            }
            mx = mx.max(x.abs());
        }
    }
    if mx == 0.0 {
        return Ok(SymEig3 {
            values: [0.0; 3],
            vectors: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        });
    }
    let k = exponent(mx);
    for row in &mut m {
        for x in row {
            *x = ldexp(*x, -k);
        }
    }
    // v[i][j]: component i of the j-th column (eigenvector candidate).
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    for _sweep in 0..50 {
        let mut off = 0.0;
        for (p, q, r) in [(0usize, 1usize, 2usize), (0, 2, 1), (1, 2, 0)] {
            let apq = m[p][q];
            if apq == 0.0 {
                continue;
            }
            let (app, aqq) = (m[p][p], m[q][q]);
            if app.abs() + 100.0 * apq.abs() == app.abs()
                && aqq.abs() + 100.0 * apq.abs() == aqq.abs()
            {
                m[p][q] = 0.0;
                m[q][p] = 0.0;
                continue;
            }
            off += apq.abs();
            let (t, c) = rotation(app, aqq, apq);
            let s = t * c;
            m[p][p] = app - t * apq;
            m[q][q] = aqq + t * apq;
            m[p][q] = 0.0;
            m[q][p] = 0.0;
            let (arp, arq) = (m[r][p], m[r][q]);
            let (np, nq) = (c * arp - s * arq, s * arp + c * arq);
            m[r][p] = np;
            m[p][r] = np;
            m[r][q] = nq;
            m[q][r] = nq;
            for row in &mut v {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
        if off == 0.0 {
            break;
        }
    }

    let mut idx = [0usize, 1, 2];
    let d = [m[0][0], m[1][1], m[2][2]];
    if d[idx[0]] > d[idx[1]] {
        idx.swap(0, 1);
    }
    if d[idx[1]] > d[idx[2]] {
        idx.swap(1, 2);
    }
    if d[idx[0]] > d[idx[1]] {
        idx.swap(0, 1);
    }
    let values = [
        ldexp(d[idx[0]], k),
        ldexp(d[idx[1]], k),
        ldexp(d[idx[2]], k),
    ];
    let col = |j: usize| [v[0][j], v[1][j], v[2][j]];
    let (e0, e1) = (col(idx[0]), col(idx[1]));
    let mut e2 = col(idx[2]);
    let det = e0[0] * (e1[1] * e2[2] - e1[2] * e2[1]) - e0[1] * (e1[0] * e2[2] - e1[2] * e2[0])
        + e0[2] * (e1[0] * e2[1] - e1[1] * e2[0]);
    if det < 0.0 {
        e2 = [-e2[0], -e2[1], -e2[2]];
    }
    Ok(SymEig3 {
        values,
        vectors: [e0, e1, e2],
    })
}

/// Eigendecomposition of the symmetric 2x2 matrix `a` (upper triangle read), ascending.
///
/// Uses the Jacobi rotation that annihilates the off-diagonal entry, so the eigenvalues are
/// `a_00 − t a_01` and `a_11 + t a_01` with `|t| ≤ 1` (no cancellation between nearly equal
/// eigenvalues, and the rotation angle is well defined even if they coincide).
///
/// # Errors
/// [`EigError::NonFinite`] if an entry of the upper triangle is NaN or infinite.
pub fn sym_eig2(a: [[f64; 2]; 2]) -> Result<SymEig2, EigError> {
    let (a00, a01, a11) = (a[0][0], a[0][1], a[1][1]);
    if !(a00.is_finite() && a01.is_finite() && a11.is_finite()) {
        return Err(EigError::NonFinite);
    }
    let mx = a00.abs().max(a01.abs()).max(a11.abs());
    if mx == 0.0 {
        return Ok(SymEig2 {
            values: [0.0; 2],
            vectors: [[1.0, 0.0], [0.0, 1.0]],
        });
    }
    let k = exponent(mx);
    let (p, b, q) = (ldexp(a00, -k), ldexp(a01, -k), ldexp(a11, -k));
    // (value, unit vector) of the two eigenpairs
    let (mut l0, mut v0, mut l1);
    if b == 0.0 || (p.abs() + 100.0 * b.abs() == p.abs() && q.abs() + 100.0 * b.abs() == q.abs()) {
        l0 = p;
        v0 = [1.0, 0.0];
        l1 = q;
    } else {
        let (t, c) = rotation(p, q, b);
        let s = t * c;
        l0 = p - t * b;
        v0 = [c, -s];
        l1 = q + t * b;
    }
    if l0 > l1 {
        std::mem::swap(&mut l0, &mut l1);
        // the other eigenvector is v0 rotated by +90 degrees
        v0 = [-v0[1], v0[0]];
    }
    Ok(SymEig2 {
        values: [ldexp(l0, k), ldexp(l1, k)],
        vectors: [v0, [-v0[1], v0[0]]],
    })
}
