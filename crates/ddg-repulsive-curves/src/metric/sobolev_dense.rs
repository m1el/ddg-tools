//! Dense Sobolev–Slobodeckij Gram matrices `B`, `B⁰`, `A` for both quadratures
//! (`SPEC.md` §4.2–4.4).
//!
//! The fast path ([`sobolev_gram`]) uses the kernel-matrix form (Eq. 31–32):
//!
//! ```text
//! B  = Dᵀ ((diag(K^B 𝟙) − K^B) ⊗ I₃) D,      B⁰ = Eᵀ (diag(K⁰ 𝟙) − K⁰) E,
//! K^B_IJ = w_IJ + w_JI,   K⁰_IJ = w⁰_IJ + w⁰_JI   for (I, J) ∈ 𝒫, 0 otherwise.
//! ```
//!
//! It costs `O(|E|²)` kernel evaluations (one `powf` per pair with the midpoint rule, four with
//! the trapezoid rule) and `O(|V| |E|)` for the products, all parallel over rows with fixed
//! summation orders, so the result is bitwise independent of the thread count (D25).
//!
//! [`sobolev_gram_increments`] assembles the same matrices by the paper's increments
//! (Sec. 5.2.2; §4.2–4.3), serially; it is the reference the tests compare with.
//!
//! The paper's Eq. 32 and its Eq. 13/20 do not disagree in their exponents
//! (`k²₄(p, q, T) |p − q|^{−(2σ+1)} = k²_{2σ+5}(p, q, T)`); they differ in the quadrature, which is
//! an option ([`GramQuadrature`], D2). The midpoint rule is the one the reference uses
//! (`ref: src/sobo_slobo.cpp`, `IntegrateLocalMatrix` with the single point 0.5; MIT,
//! Copyright (c) 2019 Christopher Yu).

use rayon::prelude::*;

use super::{DenseGram, DenseMatrix, GramQuadrature};
use crate::kernel::Exponents;
use crate::kernel::v3::{dot, sub};
use crate::network::CurveNetwork;

/// Per-edge data needed by the Gram weights: endpoints, lengths, tangents, midpoints.
///
/// Built once per geometry; shared by the dense assembly and (for the exact near-field blocks)
/// the block cluster tree, so that both use the same `K` (D2).
#[derive(Debug, Clone, PartialEq)]
pub struct GramGeometry {
    /// Endpoint positions `(γ_{i₁}, γ_{i₂})` per edge.
    pub ends: Vec<[[f64; 3]; 2]>,
    /// `ℓ_I`.
    pub length: Vec<f64>,
    /// `T_I`.
    pub tangent: Vec<[f64; 3]>,
    /// `x_I`.
    pub midpoint: Vec<[f64; 3]>,
    /// Exponent `−(2σ + 1)/2` applied to squared distances.
    half_exponent: f64,
}

impl GramGeometry {
    /// The edge data of `positions` (assumed valid, see
    /// [`CurveNetwork::validate_positions`]).
    pub fn new(network: &CurveNetwork, positions: &[[f64; 3]], exponents: Exponents) -> Self {
        let g = network.edge_geometry(positions);
        let ends = network
            .edges()
            .iter()
            .map(|&[a, b]| [positions[a as usize], positions[b as usize]])
            .collect();
        let sigma = exponents.sigma();
        GramGeometry {
            ends,
            length: g.length,
            tangent: g.tangent,
            midpoint: g.midpoint,
            half_exponent: -(2.0 * sigma + 1.0) / 2.0,
        }
    }

    /// Number of edges.
    pub fn n_edges(&self) -> usize {
        self.length.len()
    }

    /// `(|d|^{−(2σ+1)}, |d|^{−(2σ+5)})` for `d² = r2`.
    #[inline(always)]
    fn powers(&self, r2: f64) -> (f64, f64) {
        let p = r2.powf(self.half_exponent);
        (p, p / (r2 * r2))
    }

    /// The ordered weights `(w_IJ, w⁰_IJ)` of §4.2–4.3 for an ordered pair `(I, J) ∈ 𝒫`:
    ///
    /// ```text
    /// midpoint:  w_IJ = ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)},          w⁰_IJ = ℓ_I ℓ_J k²_{2σ+5}(x_I, x_J, T_I)
    /// trapezoid: w_IJ = ¼ ℓ_I ℓ_J Σ_{i∈I, j∈J} |γ_i − γ_j|^{−(2σ+1)},
    ///            w⁰_IJ = ¼ ℓ_I ℓ_J Σ_{i∈I, j∈J} k²_{2σ+5}(γ_i, γ_j, T_I)
    /// ```
    ///
    /// with `k²_{2σ+5}(p, q, T) = |T × (p − q)|² / |p − q|^{2σ+5}`. The caller must exclude pairs
    /// with a common vertex.
    ///
    /// [formula: rc-high-order-gram]
    /// [formula: rc-low-order-gram]
    pub fn weights(&self, i: usize, j: usize, quadrature: GramQuadrature) -> (f64, f64) {
        let ll = self.length[i] * self.length[j];
        let t = self.tangent[i];
        match quadrature {
            GramQuadrature::Midpoint => {
                let d = sub(self.midpoint[i], self.midpoint[j]);
                let r2 = dot(d, d);
                let td = dot(t, d);
                let (p1, p5) = self.powers(r2);
                (ll * p1, ll * (r2 - td * td) * p5)
            }
            GramQuadrature::Trapezoid => {
                let (mut w, mut w0) = (0.0, 0.0);
                for a in &self.ends[i] {
                    for b in &self.ends[j] {
                        let d = sub(*a, *b);
                        let r2 = dot(d, d);
                        let td = dot(t, d);
                        let (p1, p5) = self.powers(r2);
                        w += p1;
                        w0 += (r2 - td * td) * p5;
                    }
                }
                (0.25 * ll * w, 0.25 * ll * w0)
            }
        }
    }

    /// The symmetric kernel entries `[K^B_IJ, K⁰_IJ] = [w_IJ + w_JI, w⁰_IJ + w⁰_JI]` of Eq. 31 for
    /// `(I, J) ∈ 𝒫` (the caller excludes adjacent pairs and `I = J`). The computation is done
    /// for `(min, max)` of the pair, so `K_IJ` and `K_JI` are bitwise equal. With the midpoint
    /// rule this is Eq. 32 (with `k⁰_{2σ+1}` for `K^B`, as B.3.4 says).
    ///
    /// [formula: rc-kernel-matrix-form]
    #[inline]
    pub fn pair_kernels(&self, i: usize, j: usize, quadrature: GramQuadrature) -> [f64; 2] {
        let (i, j) = if i <= j { (i, j) } else { (j, i) };
        let ll = self.length[i] * self.length[j];
        let (ti, tj) = (self.tangent[i], self.tangent[j]);
        match quadrature {
            GramQuadrature::Midpoint => {
                let d = sub(self.midpoint[i], self.midpoint[j]);
                let r2 = dot(d, d);
                let (a, b) = (dot(ti, d), dot(tj, d));
                let (p1, p5) = self.powers(r2);
                [2.0 * ll * p1, ll * (((r2 - a * a) + (r2 - b * b)) * p5)]
            }
            GramQuadrature::Trapezoid => {
                let (mut kb, mut k0) = (0.0, 0.0);
                for p in &self.ends[i] {
                    for q in &self.ends[j] {
                        let d = sub(*p, *q);
                        let r2 = dot(d, d);
                        let (a, b) = (dot(ti, d), dot(tj, d));
                        let (p1, p5) = self.powers(r2);
                        kb += p1;
                        k0 += ((r2 - a * a) + (r2 - b * b)) * p5;
                    }
                }
                [0.5 * ll * kb, 0.25 * ll * k0]
            }
        }
    }
}

/// The symmetric `|E| × |E|` kernel matrices `K^B`, `K⁰` of Eq. 31 and their row sums.
#[derive(Debug, Clone, PartialEq)]
pub struct KernelMatrices {
    /// `K^B` (zero on the diagonal and on adjacent pairs).
    pub kb: DenseMatrix,
    /// `K⁰` (zero on the diagonal and on adjacent pairs).
    pub k0: DenseMatrix,
    /// `K^B 𝟙`, each row summed in increasing column order.
    pub kb_rowsum: Vec<f64>,
    /// `K⁰ 𝟙`.
    pub k0_rowsum: Vec<f64>,
}

/// `K^B` and `K⁰` (Eq. 31–32), parallel over rows on `pool`, bitwise independent of the thread
/// count. Positions must be valid ([`CurveNetwork::validate_positions`]).
///
/// [formula: rc-kernel-matrix-form]
pub fn kernel_matrices(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    geom: &GramGeometry,
    quadrature: GramQuadrature,
) -> KernelMatrices {
    let m = geom.n_edges();
    let mut kb = vec![0.0; m * m];
    let mut k0 = vec![0.0; m * m];
    let (kb_rowsum, k0_rowsum): (Vec<f64>, Vec<f64>) = pool.install(|| {
        kb.par_chunks_mut(m)
            .zip(k0.par_chunks_mut(m))
            .enumerate()
            .map(|(i, (rb, r0))| {
                for j in 0..m {
                    if network.edges_adjacent(i, j) {
                        continue;
                    }
                    let [b, z] = geom.pair_kernels(i, j, quadrature);
                    rb[j] = b;
                    r0[j] = z;
                }
                (rb.iter().sum::<f64>(), r0.iter().sum::<f64>())
            })
            .unzip()
    });
    KernelMatrices {
        kb: DenseMatrix { n: m, data: kb },
        k0: DenseMatrix { n: m, data: k0 },
        kb_rowsum,
        k0_rowsum,
    }
}

/// `B`, `B⁰` and `A = B + B⁰` by the kernel-matrix form (Eq. 31), parallel on `pool` and
/// bitwise independent of its size (D25). The matrices are exactly symmetric (each is
/// averaged with its transpose after the row-wise products).
///
/// Positions must be valid ([`CurveNetwork::validate_positions`]); [`super::gram`] is the
/// checked entry point.
///
/// [formula: rc-kernel-matrix-form]
pub fn sobolev_gram(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    quadrature: GramQuadrature,
) -> DenseGram {
    let geom = GramGeometry::new(network, positions, exponents);
    let k = kernel_matrices(pool, network, &geom, quadrature);
    gram_from_kernels(pool, network, &geom, &k)
}

/// The products of Eq. 31 for given kernel matrices.
///
/// Row `u` of `B` is `Σ_{I ∋ u} (s_{I,u}/ℓ_I) Σ_J L^B_IJ ⟨T_I, T_J⟩/ℓ_J (e_{j₂} − e_{j₁})` with
/// `L = diag(K𝟙) − K`, and row `u` of `B⁰` is `Σ_{I ∋ u} ¼ Σ_J L⁰_IJ (e_{j₁} + e_{j₂})`.
///
/// [formula: rc-kernel-matrix-form]
pub fn gram_from_kernels(
    pool: &rayon::ThreadPool,
    network: &CurveNetwork,
    geom: &GramGeometry,
    k: &KernelMatrices,
) -> DenseGram {
    let n = network.n_vertices();
    let m = geom.n_edges();
    let edges = network.edges();
    let t_over_l: Vec<[f64; 3]> = (0..m)
        .map(|j| {
            let (t, l) = (geom.tangent[j], geom.length[j]);
            [t[0] / l, t[1] / l, t[2] / l]
        })
        .collect();
    let mut b = vec![0.0; n * n];
    let mut b0 = vec![0.0; n * n];
    pool.install(|| {
        b.par_chunks_mut(n)
            .zip(b0.par_chunks_mut(n))
            .enumerate()
            .for_each(|(u, (rb, r0))| {
                for &ie in network.incident_edges(u as u32) {
                    let i = ie as usize;
                    let s = if edges[i][1] as usize == u { 1.0 } else { -1.0 };
                    let ti = t_over_l[i];
                    let si = [s * ti[0], s * ti[1], s * ti[2]];
                    let (kb_row, k0_row) = (k.kb.row(i), k.k0.row(i));
                    for j in 0..m {
                        let (lb, l0) = if j == i {
                            (k.kb_rowsum[i], k.k0_rowsum[i])
                        } else {
                            (-kb_row[j], -k0_row[j])
                        };
                        if lb == 0.0 && l0 == 0.0 {
                            continue;
                        }
                        let [j1, j2] = edges[j];
                        let c = lb * dot(si, t_over_l[j]);
                        rb[j2 as usize] += c;
                        rb[j1 as usize] -= c;
                        let c0 = 0.25 * l0;
                        r0[j1 as usize] += c0;
                        r0[j2 as usize] += c0;
                    }
                }
            });
    });
    let (b, b0, a) = symmetrize(pool, n, &b, &b0);
    DenseGram { b, b0, a }
}

/// `(M + Mᵀ)/2` for both matrices, and their sum; parallel over rows, exact symmetry.
fn symmetrize(
    pool: &rayon::ThreadPool,
    n: usize,
    b: &[f64],
    b0: &[f64],
) -> (DenseMatrix, DenseMatrix, DenseMatrix) {
    let mut sb = vec![0.0; n * n];
    let mut sb0 = vec![0.0; n * n];
    let mut sa = vec![0.0; n * n];
    pool.install(|| {
        sb.par_chunks_mut(n)
            .zip(sb0.par_chunks_mut(n))
            .zip(sa.par_chunks_mut(n))
            .enumerate()
            .for_each(|(u, ((rb, r0), ra))| {
                for v in 0..n {
                    let x = 0.5 * (b[u * n + v] + b[v * n + u]);
                    let y = 0.5 * (b0[u * n + v] + b0[v * n + u]);
                    rb[v] = x;
                    r0[v] = y;
                    ra[v] = x + y;
                }
            });
    });
    (
        DenseMatrix { n, data: sb },
        DenseMatrix { n, data: sb0 },
        DenseMatrix { n, data: sa },
    )
}

/// `B`, `B⁰`, `A` by the increments of the paper's Sec. 5.2.2 (§4.2–4.3), summed over the
/// ordered pairs `(I, J) ∈ 𝒫` in lexicographic order, serially. For every pair and
/// `a, b ∈ {1, 2}`:
///
/// ```text
/// B_{i_a i_b} += (−1)^{a+b} w_IJ/ℓ_I²,     B_{i_a j_b} −= (−1)^{a+b} w_IJ ⟨T_I, T_J⟩/(ℓ_I ℓ_J),
/// B_{j_a j_b} += (−1)^{a+b} w_IJ/ℓ_J²,     B_{j_a i_b} −= (−1)^{a+b} w_IJ ⟨T_J, T_I⟩/(ℓ_J ℓ_I),
/// B⁰_{i_a i_b} += ¼ w⁰_IJ,  B⁰_{i_a j_b} −= ¼ w⁰_IJ,  B⁰_{j_a i_b} −= ¼ w⁰_IJ,  B⁰_{j_a j_b} += ¼ w⁰_IJ.
/// ```
///
/// `O(|E|²)` and single-threaded: a reference for tests, not for production.
///
/// [formula: rc-high-order-gram]
/// [formula: rc-low-order-gram]
pub fn sobolev_gram_increments(
    network: &CurveNetwork,
    positions: &[[f64; 3]],
    exponents: Exponents,
    quadrature: GramQuadrature,
) -> DenseGram {
    let n = network.n_vertices();
    let geom = GramGeometry::new(network, positions, exponents);
    let m = geom.n_edges();
    let edges = network.edges();
    let mut b = DenseMatrix::zeros(n);
    let mut b0 = DenseMatrix::zeros(n);
    let sign = |a: usize, c: usize| if (a + c).is_multiple_of(2) { 1.0 } else { -1.0 };
    for i in 0..m {
        for j in 0..m {
            if network.edges_adjacent(i, j) {
                continue;
            }
            let (w, w0) = geom.weights(i, j, quadrature);
            let (li, lj) = (geom.length[i], geom.length[j]);
            let tt = dot(geom.tangent[i], geom.tangent[j]);
            let (ei, ej) = (edges[i], edges[j]);
            for a in 0..2 {
                for c in 0..2 {
                    let s = sign(a, c);
                    let (ia, ic) = (ei[a] as usize, ei[c] as usize);
                    let (ja, jc) = (ej[a] as usize, ej[c] as usize);
                    b.data[ia * n + ic] += s * w / (li * li);
                    b.data[ia * n + jc] -= s * w * tt / (li * lj);
                    b.data[ja * n + jc] += s * w / (lj * lj);
                    b.data[ja * n + ic] -= s * w * tt / (lj * li);
                    b0.data[ia * n + ic] += 0.25 * w0;
                    b0.data[ia * n + jc] -= 0.25 * w0;
                    b0.data[ja * n + ic] -= 0.25 * w0;
                    b0.data[ja * n + jc] += 0.25 * w0;
                }
            }
        }
    }
    let mut a = DenseMatrix::zeros(n);
    for k in 0..n * n {
        a.data[k] = b.data[k] + b0.data[k];
    }
    DenseGram { b, b0, a }
}
