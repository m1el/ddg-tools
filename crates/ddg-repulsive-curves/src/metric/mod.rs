//! Gram matrices of the Sobolev–Slobodeckij metric and the classic L², H¹, H² metrics
//! (`SPEC.md` §4).
//!
//! * [`sobolev_dense`]: the dense `B` (Eq. 21), `B⁰` (Eq. 20) and `A = B + B⁰` (Eq. 19) with
//!   either quadrature (D2), assembled in the kernel-matrix form of Eq. 31–32 in parallel, and a
//!   serial reference assembly by the increments of the paper's Sec. 5.2.2.
//! * [`sparse_classic`]: the sparse comparison metrics of §4.6 (D23).
//! * [`DofLayout`]: the interleaved coordinate layout of `Ā = A ⊗ I₃` (D29), with the planar mode
//!   of §2.2 (D18) dropping the held coordinate.
//!
//! # Properties
//!
//! `B`, `B⁰`, `A` are symmetric positive semidefinite; `ker B` contains the functions that are
//! constant on each component, `ker B⁰` the global constants, `ker A = ker B ∩ ker B⁰`, which is
//! `span{𝟙}` for a generic network (so `Ā` has exactly the translations in its kernel, D5).
//!
//! [formula: rc-gram-psd-kernel]
//!
//! Under `γ → cγ` (`c > 0`) both terms scale alike:
//! `B(cγ) = c^{−(2σ+1)} B(γ)`, `B⁰(cγ) = c^{−(2σ+1)} B⁰(γ)` (the printed exponents of the
//! paper's A.4.2 are wrong, D3).
//!
//! [formula: rc-gram-scaling]

use crate::energy::exact::make_pool;
use crate::error::RcError;
use crate::kernel::Exponents;
use crate::network::{Axis, CurveNetwork, EdgeGeometry};

pub mod hierarchical;
pub mod sobolev_dense;
pub mod sparse_classic;

pub use sobolev_dense::{GramGeometry, KernelMatrices, kernel_matrices, sobolev_gram};
pub use sparse_classic::{bilaplacian, classic_metric, curve_laplacian, mass_matrix};

/// The quadrature of the Sobolev–Slobodeckij Gram matrices (D2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GramQuadrature {
    /// One point at the edge midpoints (Eq. 32; the reference, dense and BCT). The default.
    #[default]
    Midpoint,
    /// The trapezoid rule over the four endpoint pairs (Eq. 20–21 as printed).
    Trapezoid,
}

/// The metric of the descent direction (`SPEC.md` §9.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Metric {
    /// The fractional Sobolev–Slobodeckij metric `A = B + B⁰` (Sec. 4–5).
    Sobolev {
        /// Quadrature of `B` and `B⁰`.
        quadrature: GramQuadrature,
    },
    /// The lumped mass matrix `diag(m_v)` (evaluation only, §4.6).
    L2,
    /// The curve Laplacian with weights `1/ℓ` (evaluation only, §4.6).
    H1,
    /// `L M⁻¹ L` (evaluation only, §4.6).
    H2,
}

impl Default for Metric {
    /// `Sobolev { quadrature: Midpoint }` (`SPEC.md` §10).
    fn default() -> Self {
        Metric::Sobolev {
            quadrature: GramQuadrature::Midpoint,
        }
    }
}

/// The layout of a vector field on the vertices as one flat vector (D29): the coordinates of a
/// vertex are consecutive, `(γ₁ₓ, γ₁ᵧ, γ₁_z, γ₂ₓ, …)`. In planar mode (§2.2, D18) the held
/// coordinate is dropped, so each vertex has two entries (in increasing coordinate order).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DofLayout {
    planar: Option<Axis>,
}

impl DofLayout {
    /// The layout for `planar` (`None`: all three coordinates).
    pub fn new(planar: Option<Axis>) -> Self {
        DofLayout { planar }
    }

    /// The held axis in planar mode.
    pub fn planar(&self) -> Option<Axis> {
        self.planar
    }

    /// Number of coordinates per vertex `d` (3, or 2 in planar mode).
    pub fn dim(&self) -> usize {
        self.coords().len()
    }

    /// The active coordinate indices, increasing.
    pub fn coords(&self) -> &'static [usize] {
        match self.planar {
            None => &[0, 1, 2],
            Some(Axis::X) => &[1, 2],
            Some(Axis::Y) => &[0, 2],
            Some(Axis::Z) => &[0, 1],
        }
    }

    /// The slot (position within a vertex block) of coordinate `coord`, `None` if it is held.
    pub fn slot(&self, coord: usize) -> Option<usize> {
        self.coords().iter().position(|&c| c == coord)
    }

    /// Length `d |V|` of a flat vector.
    pub fn n_dofs(&self, n_vertices: usize) -> usize {
        self.dim() * n_vertices
    }

    /// Flat index of slot `slot` of vertex `v`.
    #[inline]
    pub fn index(&self, v: usize, slot: usize) -> usize {
        v * self.dim() + slot
    }

    /// Flattens per-vertex 3-vectors, dropping the held coordinate.
    pub fn gather(&self, x: &[[f64; 3]]) -> Vec<f64> {
        let cs = self.coords();
        let mut out = Vec::with_capacity(cs.len() * x.len());
        for p in x {
            for &c in cs {
                out.push(p[c]);
            }
        }
        out
    }

    /// Unflattens to per-vertex 3-vectors; the held coordinate is 0.
    pub fn scatter(&self, dofs: &[f64]) -> Vec<[f64; 3]> {
        let cs = self.coords();
        let d = cs.len();
        dofs.chunks_exact(d)
            .map(|ch| {
                let mut p = [0.0; 3];
                for (s, &c) in cs.iter().enumerate() {
                    p[c] = ch[s];
                }
                p
            })
            .collect()
    }

    /// Adds `scale · dofs` to `positions` (the held coordinate is unchanged).
    pub fn add_scaled(&self, positions: &mut [[f64; 3]], scale: f64, dofs: &[f64]) {
        let cs = self.coords();
        let d = cs.len();
        for (p, ch) in positions.iter_mut().zip(dofs.chunks_exact(d)) {
            for (s, &c) in cs.iter().enumerate() {
                p[c] += scale * ch[s];
            }
        }
    }
}

/// A dense square matrix in row-major order (the Gram matrices are symmetric, so the order only
/// matters for non-symmetric uses).
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    n: usize,
    data: Vec<f64>,
}

impl DenseMatrix {
    /// The `n × n` zero matrix.
    pub fn zeros(n: usize) -> Self {
        DenseMatrix {
            n,
            data: vec![0.0; n * n],
        }
    }

    /// Wraps `n · n` row-major values.
    ///
    /// # Errors
    /// [`RcError::InvalidOption`] if `data.len() != n²`.
    pub fn from_row_major(n: usize, data: Vec<f64>) -> Result<Self, RcError> {
        if data.len() != n * n {
            return Err(RcError::InvalidOption {
                name: "DenseMatrix data",
                value: data.len().to_string(),
                requirement: "exactly n * n values",
            });
        }
        Ok(DenseMatrix { n, data })
    }

    /// Dimension `n`.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Entry `(i, j)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    /// Row `i`.
    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.n..(i + 1) * self.n]
    }

    /// The row-major values.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// The row-major values, mutable.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// `y = M x` (sequential, fixed summation order).
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|i| self.row(i).iter().zip(x).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// `uᵀ M v`.
    pub fn bilinear(&self, u: &[f64], v: &[f64]) -> f64 {
        u.iter().zip(self.mul_vec(v)).map(|(a, b)| a * b).sum()
    }

    /// `max_ij |M_ij − M_ji|`.
    pub fn symmetry_defect(&self) -> f64 {
        let mut d: f64 = 0.0;
        for i in 0..self.n {
            for j in 0..i {
                d = d.max((self.get(i, j) - self.get(j, i)).abs());
            }
        }
        d
    }

    /// `max_ij |M_ij|`.
    pub fn max_abs(&self) -> f64 {
        self.data.iter().fold(0.0f64, |m, x| m.max(x.abs()))
    }
}

/// The Gram matrices of a metric on vertex functions, `|V| × |V|` (the scalar matrices; the
/// metric on positions is `Ā = A ⊗ I_d`, see [`DofLayout`]).
///
/// For [`Metric::Sobolev`]: `b = B` (Eq. 21), `b0 = B⁰` (Eq. 20), `a = B + B⁰` (Eq. 19). For the
/// classic metrics (§4.6) the matrix is put where its order belongs: L² in `b0` (`b = 0`), H¹
/// and H² in `b` (`b0 = 0`); `a` is the metric in all cases.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseGram {
    /// High-order term.
    pub b: DenseMatrix,
    /// Low-order term.
    pub b0: DenseMatrix,
    /// The metric `A = B + B⁰`.
    pub a: DenseMatrix,
}

/// The dense Gram matrices of `metric` at `positions` (for tests, the oracle and small
/// problems; the solvers use [`sobolev_gram`] and [`classic_metric`] directly).
///
/// `threads` is the size of the explicit pool (`None`: all cores); the result is bitwise
/// independent of it (D25).
///
/// # Errors
/// The errors of [`CurveNetwork::validate_positions`]; [`RcError::InvalidOption`] for
/// `threads = Some(0)`.
pub fn gram(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    metric: &Metric,
    threads: Option<usize>,
) -> Result<DenseGram, RcError> {
    network.validate_positions(positions)?;
    let n = network.n_vertices();
    match *metric {
        Metric::Sobolev { quadrature } => {
            let pool = make_pool(threads)?;
            Ok(sobolev_gram(
                &pool, network, positions, exponents, quadrature,
            ))
        }
        Metric::L2 | Metric::H1 | Metric::H2 => {
            let s = classic_metric(network, positions, metric)?;
            let mut m = DenseMatrix::zeros(n);
            for (r, c, v) in s.triplets() {
                m.data[r as usize * n + c as usize] += v;
            }
            let zero = DenseMatrix::zeros(n);
            Ok(if matches!(metric, Metric::L2) {
                DenseGram {
                    b: zero,
                    b0: m.clone(),
                    a: m,
                }
            } else {
                DenseGram {
                    b: m.clone(),
                    b0: zero,
                    a: m,
                }
            })
        }
    }
}

/// `(D u)_I = (u_{i₂} − u_{i₁}) T_I / ℓ_I` for a vertex function `u` (Sec. 5.2.1).
///
/// [formula: rc-derivative-operator]
pub fn apply_derivative(network: &CurveNetwork, geom: &EdgeGeometry, u: &[f64]) -> Vec<[f64; 3]> {
    network
        .edges()
        .iter()
        .enumerate()
        .map(|(i, &[a, b])| {
            let s = (u[b as usize] - u[a as usize]) / geom.length[i];
            let t = geom.tangent[i];
            [s * t[0], s * t[1], s * t[2]]
        })
        .collect()
}

/// `Dᵀ y` for an edge field `y ∈ (R³)^{|E|}` (summed in increasing edge order).
///
/// [formula: rc-derivative-operator]
pub fn apply_derivative_transpose(
    network: &CurveNetwork,
    geom: &EdgeGeometry,
    y: &[[f64; 3]],
) -> Vec<f64> {
    let mut out = vec![0.0; network.n_vertices()];
    for (i, &[a, b]) in network.edges().iter().enumerate() {
        let t = geom.tangent[i];
        let s = (t[0] * y[i][0] + t[1] * y[i][1] + t[2] * y[i][2]) / geom.length[i];
        out[b as usize] += s;
        out[a as usize] -= s;
    }
    out
}

/// `(E u)_I = (u_{i₁} + u_{i₂}) / 2` (B.3.4).
///
/// [formula: rc-derivative-operator]
pub fn apply_average(network: &CurveNetwork, u: &[f64]) -> Vec<f64> {
    network
        .edges()
        .iter()
        .map(|&[a, b]| 0.5 * (u[a as usize] + u[b as usize]))
        .collect()
}

/// `Eᵀ y` for an edge function `y` (summed in increasing edge order).
///
/// [formula: rc-derivative-operator]
pub fn apply_average_transpose(network: &CurveNetwork, y: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; network.n_vertices()];
    for (i, &[a, b]) in network.edges().iter().enumerate() {
        out[a as usize] += 0.5 * y[i];
        out[b as usize] += 0.5 * y[i];
    }
    out
}
