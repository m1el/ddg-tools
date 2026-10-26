//! Near-contacts in the multigrid hierarchy (`SPEC.md` §8.4, §17g): the Galerkin correction of
//! the rediscretized coarse metrics and the block preconditioner of the smoother.
//!
//! The coarse metric `Ā_{ℓ+1}` is rediscretized on the coarse curve (§8.4). That works when the
//! coarse curve resolves the couplings of the fine one, and fails at **near-contacts**: edge
//! pairs whose midpoints are closer than their length (two branches that nearly cross, as in the
//! authors' CapillaryNetwork scene) carry weights `ℓ_I ℓ_J |x_I − x_J|^{−(2σ+1)}` and
//! `k²_{2σ+5}` far above their neighbours', which the coarse edges, whose midpoints are elsewhere,
//! do not reproduce: on that scene `(Πx)ᵀ H_ℓ (Πx) / xᵀ H_{ℓ+1} x` is 53 for a smooth `x`, the
//! coarse correction overshoots by that factor, and the cycles diverge.
//!
//! **Correction.** For the coarse edge pairs `(I', J')` that contain a near pair of the finer
//! level (or a pair corrected there), the level's metric replaces the rediscretized pair term
//! `R(I', J')` by the Galerkin sum `S(I', J') = Σ_{(I, J) → (I', J')} Jᵀ E(I, J) J` over the
//! finer pairs inside it, where `E` is the finer level's effective pair term (its own `S` if the
//! pair was corrected there, else its exact term): the level's metric is
//! `Ā_ℓ + Σ (S − R) ⊗ I_d`. Pair terms are scalar `4 × 4` blocks on the vertices
//! `(i₁, i₂, j₁, j₂)` (Eq. 31: `K^B (aaᵀ + bbᵀ − ⟨T_I, T_J⟩ (abᵀ + baᵀ)) + K⁰ eeᵀ`,
//! `a = (−1, 1, 0, 0)/ℓ_I`, `b = (0, 0, −1, 1)/ℓ_J`, `e = (1, 1, −1, −1)/2`), so the correction
//! is local and costs `O(#pairs)` per product. It is a partial Galerkin coarse operator; the
//! constraint projectors between the levels are not included (the operators remain
//! `P_ℓ Ā_ℓ P_ℓ`).
//!
//! **Smoother.** The near pairs also make the metric stiff locally (condition numbers grow as
//! the gap shrinks); [`BlockSmoother`] preconditions the CG smoother with the diagonal plus the
//! exact near-pair blocks, which solves those couplings exactly.

use std::collections::BTreeMap;

use crate::kernel::v3::{dot, sub};
use crate::metric::GramQuadrature;
use crate::metric::sobolev_dense::GramGeometry;
use crate::network::CurveNetwork;

use super::coarsen::{Coarsening, ProlongationRow};

/// A scalar `4 × 4` block on the vertices `(i₁, i₂, j₁, j₂)` of an edge pair.
pub type Block = [[f64; 4]; 4];

/// The corrected pair `(I, J)` (`I ≤ J`, possibly equal or adjacent) of one level.
#[derive(Debug, Clone, PartialEq)]
pub struct PairCorrection {
    /// The edge pair, `I ≤ J`.
    pub pair: [u32; 2],
    /// The vertices `(i₁, i₂, j₁, j₂)` of the block.
    pub vertices: [u32; 4],
    /// The Galerkin pair term `S(I, J)`.
    pub galerkin: Block,
    /// `S(I, J) − R(I, J)`, what is added to the rediscretized metric.
    pub delta: Block,
}

/// The corrections of one level, sorted by pair.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NearCorrection {
    pairs: Vec<PairCorrection>,
}

/// The exact pair term `R(I, J)` of Eq. 31 on `(i₁, i₂, j₁, j₂)`; zero for `I = J` and for
/// adjacent edges (they are not in `𝒫`).
///
/// [formula: rc-kernel-matrix-form]
pub fn pair_block(
    network: &CurveNetwork,
    geom: &GramGeometry,
    i: usize,
    j: usize,
    quadrature: GramQuadrature,
) -> Block {
    let mut m = [[0.0; 4]; 4];
    let edges = network.edges();
    let (ei, ej) = (edges[i], edges[j]);
    if i == j || ei.iter().any(|v| ej.contains(v)) {
        return m;
    }
    let [kb, k0] = geom.pair_kernels(i, j, quadrature);
    let (li, lj) = (geom.length[i], geom.length[j]);
    let c = dot(geom.tangent[i], geom.tangent[j]);
    let a = [-1.0 / li, 1.0 / li, 0.0, 0.0];
    let b = [0.0, 0.0, -1.0 / lj, 1.0 / lj];
    let e = [0.5, 0.5, -0.5, -0.5];
    for r in 0..4 {
        for s in 0..4 {
            m[r][s] = kb * (a[r] * a[s] + b[r] * b[s] - c * (a[r] * b[s] + b[r] * a[s]))
                + k0 * e[r] * e[s];
        }
    }
    m
}

/// The non-adjacent edge pairs `I < J` whose midpoints are closer than
/// `kappa · max(ℓ_I, ℓ_J)`, sorted (uniform grid of the midpoints, cell `kappa · max ℓ`).
pub fn near_pairs(network: &CurveNetwork, geom: &GramGeometry, kappa: f64) -> Vec<[u32; 2]> {
    let m = geom.n_edges();
    if kappa.is_nan() || kappa <= 0.0 || m < 2 {
        return Vec::new();
    }
    let lmax = geom.length.iter().copied().fold(0.0, f64::max);
    let h = kappa * lmax;
    if !h.is_finite() || h <= 0.0 {
        return Vec::new();
    }
    let cell = |p: [f64; 3]| -> [i64; 3] {
        [
            (p[0] / h).floor() as i64,
            (p[1] / h).floor() as i64,
            (p[2] / h).floor() as i64,
        ]
    };
    let mut grid: BTreeMap<[i64; 3], Vec<u32>> = BTreeMap::new();
    for e in 0..m {
        grid.entry(cell(geom.midpoint[e]))
            .or_default()
            .push(e as u32);
    }
    let edges = network.edges();
    let mut out = Vec::new();
    for i in 0..m {
        let c = cell(geom.midpoint[i]);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(list) = grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) else {
                        continue;
                    };
                    for &j in list {
                        let j = j as usize;
                        if j <= i || edges[i].iter().any(|v| edges[j].contains(v)) {
                            continue;
                        }
                        let d = sub(geom.midpoint[i], geom.midpoint[j]);
                        let r = kappa * geom.length[i].max(geom.length[j]);
                        if dot(d, d) < r * r {
                            out.push([i as u32, j as u32]);
                        }
                    }
                }
            }
        }
    }
    out.sort_unstable();
    out
}

impl NearCorrection {
    /// No correction (the finest level).
    pub fn none() -> Self {
        NearCorrection::default()
    }

    /// The corrected pairs.
    pub fn pairs(&self) -> &[PairCorrection] {
        &self.pairs
    }

    /// Whether there is nothing to correct.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// The corrections of the coarse level of `coarsening`, from the fine level's network,
    /// geometry and corrections: the coarse pairs that contain a fine near pair
    /// ([`near_pairs`] with `kappa`) or a corrected fine pair.
    pub fn coarsen(
        fine_network: &CurveNetwork,
        fine_geom: &GramGeometry,
        fine: &NearCorrection,
        coarsening: &Coarsening,
        coarse_geom: &GramGeometry,
        quadrature: GramQuadrature,
        kappa: f64,
    ) -> NearCorrection {
        let e2c = coarsening.edge_to_coarse();
        let coarse = coarsening.coarse();
        let key = |a: u32, b: u32| if a <= b { [a, b] } else { [b, a] };
        let mut keys: Vec<[u32; 2]> = Vec::new();
        for p in near_pairs(fine_network, fine_geom, kappa) {
            keys.push(key(e2c[p[0] as usize], e2c[p[1] as usize]));
        }
        for c in &fine.pairs {
            keys.push(key(e2c[c.pair[0] as usize], e2c[c.pair[1] as usize]));
        }
        keys.sort_unstable();
        keys.dedup();
        if keys.is_empty() {
            return NearCorrection::none();
        }
        let fine_index: BTreeMap<[u32; 2], usize> = fine
            .pairs
            .iter()
            .enumerate()
            .map(|(k, c)| (c.pair, k))
            .collect();
        let fine_edges = fine_network.edges();
        let coarse_edges = coarse.edges();
        let parts = coarsening.coarse_edge_parts();
        let rows = coarsening.rows();
        let mut pairs = Vec::with_capacity(keys.len());
        for k in keys {
            let (ci, cj) = (k[0] as usize, k[1] as usize);
            let verts = [
                coarse_edges[ci][0],
                coarse_edges[ci][1],
                coarse_edges[cj][0],
                coarse_edges[cj][1],
            ];
            // (slot, weight) of a fine vertex of a fine edge inside coarse edge `ce` at `base`.
            let slots = |v: u32, ce: usize, base: usize| -> [(usize, f64); 2] {
                let ends = coarse_edges[ce];
                let slot = |c: u32| if c == ends[0] { base } else { base + 1 };
                match rows[v as usize] {
                    ProlongationRow::Black(c) => [(slot(c), 1.0), (slot(c), 0.0)],
                    ProlongationRow::White(a, b) => [(slot(a), 0.5), (slot(b), 0.5)],
                }
            };
            let mut s: Block = [[0.0; 4]; 4];
            let mut add = |fi: u32, fj: u32| {
                let fk = key(fi, fj);
                let block = match fine_index.get(&fk) {
                    Some(&x) => fine.pairs[x].galerkin,
                    None => pair_block(
                        fine_network,
                        fine_geom,
                        fk[0] as usize,
                        fk[1] as usize,
                        quadrature,
                    ),
                };
                // Fine block vertices in the key's orientation: the edge inside coarse edge ci
                // first.
                let (first, second) = if fi == fk[0] { (0, 2) } else { (2, 0) };
                let fv = [
                    fine_edges[fk[0] as usize][0],
                    fine_edges[fk[0] as usize][1],
                    fine_edges[fk[1] as usize][0],
                    fine_edges[fk[1] as usize][1],
                ];
                // map[r] = coarse slots of fine block row r
                let mut map = [[(0usize, 0.0f64); 2]; 4];
                for r in 0..2 {
                    map[first + r] = slots(fv[first + r], ci, 0);
                    map[second + r] = slots(fv[second + r], cj, 2);
                }
                for r in 0..4 {
                    for q in 0..4 {
                        let x = block[r][q];
                        if x == 0.0 {
                            continue;
                        }
                        for &(sr, wr) in &map[r] {
                            for &(sq, wq) in &map[q] {
                                s[sr][sq] += wr * wq * x;
                            }
                        }
                    }
                }
            };
            if ci == cj {
                let p = &parts[ci];
                for (x, &a) in p.iter().enumerate() {
                    for &b in &p[x..] {
                        if a == b && !fine_index.contains_key(&[a, a]) {
                            continue;
                        }
                        add(a, b);
                    }
                }
            } else {
                for &a in &parts[ci] {
                    for &b in &parts[cj] {
                        add(a, b);
                    }
                }
            }
            let r = pair_block(coarse, coarse_geom, ci, cj, quadrature);
            let mut delta = [[0.0; 4]; 4];
            for x in 0..4 {
                for y in 0..4 {
                    delta[x][y] = s[x][y] - r[x][y];
                }
            }
            pairs.push(PairCorrection {
                pair: k,
                vertices: verts,
                galerkin: s,
                delta,
            });
        }
        NearCorrection { pairs }
    }

    /// `y += (Σ δ) ⊗ I_d x` (sequential, fixed order).
    pub fn apply_add(&self, x: &[f64], y: &mut [f64], d: usize) {
        for c in &self.pairs {
            for r in 0..4 {
                let vr = c.vertices[r] as usize;
                for q in 0..4 {
                    let w = c.delta[r][q];
                    if w == 0.0 {
                        continue;
                    }
                    let vq = c.vertices[q] as usize;
                    for s in 0..d {
                        y[vr * d + s] += w * x[vq * d + s];
                    }
                }
            }
        }
    }

    /// `A += Σ δ` for a dense scalar `n × n` matrix in row-major order.
    pub fn add_to_dense(&self, a: &mut [f64], n: usize) {
        for c in &self.pairs {
            for r in 0..4 {
                for q in 0..4 {
                    a[c.vertices[r] as usize * n + c.vertices[q] as usize] += c.delta[r][q];
                }
            }
        }
    }
}

/// The preconditioner of the smoother of one level (`SPEC.md` §17g): the scalar matrix
/// `M = D_rest + Σ_{(I, J) near} E(I, J)` (applied as `M ⊗ I_d`), with `E` the effective pair
/// terms of the level's near pairs (exact, or the Galerkin term where corrected) and `D_rest`
/// the diagonal of the level's metric without those pairs. It solves the stiff local couplings
/// of near-contacts exactly, which a plain CG smoother cannot.
pub struct BlockSmoother {
    factor: ddg_core::sparse::Factorization,
    n: usize,
}

impl std::fmt::Debug for BlockSmoother {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BlockSmoother").field("n", &self.n).finish()
    }
}

impl BlockSmoother {
    /// The preconditioner of a level from the diagonal `diag` of its metric (corrections
    /// included), its near pairs (`kappa`) and its corrections. Returns `None` (the smoother
    /// stays plain CG) when the level has no near pair and no correction, or if the matrix
    /// cannot be factored.
    pub fn new(
        network: &CurveNetwork,
        geom: &GramGeometry,
        correction: &NearCorrection,
        diag: &[f64],
        kappa: f64,
        quadrature: GramQuadrature,
    ) -> Option<BlockSmoother> {
        use ddg_core::sparse::{Method, Parallelism, SparseMatrix, Symbolic};
        let n = network.n_vertices();
        if diag.len() != n || n == 0 {
            return None;
        }
        let edges = network.edges();
        let corrected: BTreeMap<[u32; 2], usize> = correction
            .pairs
            .iter()
            .enumerate()
            .map(|(k, c)| (c.pair, k))
            .collect();
        let mut keys = near_pairs(network, geom, kappa);
        keys.extend(correction.pairs.iter().map(|c| c.pair));
        if keys.is_empty() {
            // No near-contact on this level: plain CG smoothing, as on smooth closed curves.
            return None;
        }
        keys.sort_unstable();
        keys.dedup();
        let mut rest = diag.to_vec();
        let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
        for k in keys {
            let (block, verts) = match corrected.get(&k) {
                Some(&x) => (correction.pairs[x].galerkin, correction.pairs[x].vertices),
                None => {
                    let (i, j) = (k[0] as usize, k[1] as usize);
                    (
                        pair_block(network, geom, i, j, quadrature),
                        [edges[i][0], edges[i][1], edges[j][0], edges[j][1]],
                    )
                }
            };
            for r in 0..4 {
                rest[verts[r] as usize] -= block[r][r];
                for q in 0..4 {
                    if block[r][q] != 0.0 {
                        triplets.push((verts[r], verts[q], block[r][q]));
                    }
                }
            }
        }
        for (v, (&r, &dv)) in rest.iter().zip(diag).enumerate() {
            // The rest of the diagonal is a sum of non-negative pair diagonals; rounding (and
            // the far-field approximation) can make it slightly smaller: keep M definite.
            triplets.push((v as u32, v as u32, r.max(1e-2 * dv)));
        }
        let m = SparseMatrix::from_triplets(n, n, &triplets).ok()?;
        let sym = Symbolic::analyze(&m, Method::Llt).ok()?;
        let factor = sym.factor(&m, Parallelism::Sequential).ok()?;
        Some(BlockSmoother { factor, n })
    }

    /// `(M⁻¹ ⊗ I_d) r`.
    pub fn solve(&self, r: &[f64], d: usize) -> Option<Vec<f64>> {
        use ddg_core::sparse::Parallelism;
        if d == 0 || r.len() != self.n * d {
            return None;
        }
        // channels as columns
        let mut rhs = vec![0.0; r.len()];
        for v in 0..self.n {
            for s in 0..d {
                rhs[s * self.n + v] = r[v * d + s];
            }
        }
        self.factor
            .solve_many_in_place(&mut rhs, d, Parallelism::Sequential)
            .ok()?;
        let mut out = vec![0.0; r.len()];
        for v in 0..self.n {
            for s in 0..d {
                out[v * d + s] = rhs[s * self.n + v];
            }
        }
        Some(out)
    }
}

/// The diagonal of the scalar metric `A` from the kernel row sums `[K^B 𝟙, K⁰ 𝟙]` per edge:
/// `A_uu = Σ_{I ∋ u} ((K^B 𝟙)_I / ℓ_I² + (K⁰ 𝟙)_I / 4)` (Eq. 31: adjacent pairs and `I = J`
/// have `K = 0`), plus the diagonal of the corrections.
pub fn metric_diagonal(
    network: &CurveNetwork,
    lengths: &[f64],
    row_sums: &[[f64; 2]],
    correction: &NearCorrection,
) -> Vec<f64> {
    let mut diag = vec![0.0; network.n_vertices()];
    for (i, &[a, b]) in network.edges().iter().enumerate() {
        let v = row_sums[i][0] / (lengths[i] * lengths[i]) + 0.25 * row_sums[i][1];
        diag[a as usize] += v;
        diag[b as usize] += v;
    }
    for c in &correction.pairs {
        for r in 0..4 {
            diag[c.vertices[r] as usize] += c.delta[r][r];
        }
    }
    diag
}
