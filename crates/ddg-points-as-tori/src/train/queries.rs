//! Per-neighborhood query sampling and ground truth (SPEC.md §10.1, D20, D22; NOTES-REFERENCE.md
//! §9.1), ported from `generate_queries_for_neighborhood` of the reference's
//! `training/mesh_sampling_jax.py` (points-as-tori, MIT, Copyright (c) 2026 Nicole Feng).
//!
//! For every row (neighborhood) `3·m` queries (`m = 40`):
//! 1. **zero set**: faces whose bounding box meets the row's regularized box (each extent at least
//!    `min_extent`), chosen with probability ∝ area, a uniform point on the face;
//! 2. **narrow band**: the same face rule with fresh draws, offset along the unit face normal by
//!    `U(−band, band)`;
//! 3. **cube**: uniform in a cube around the center of the row's (unregularized) box with half-side
//!    `h = c·max(max extent, min_side)`, `c = ½√3` (code, default) or `½` (paper).
//!
//! Randomness is ChaCha8, one stream per (seed, epoch, shape, row) (D22). The reference's JAX keys
//! (threefry, `fold_in`) are not reproduced; distributions are the same.

use super::truth::{GwnPrecision, MeshTruth};
use crate::PatError;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::{RngCore, SeedableRng};
use rayon::prelude::*;

/// The box sizes (D20).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QueryBoxes {
    /// The code's (default; the weights were trained with them): cube half-side `½√3·extent`,
    /// band `±0.2`.
    #[default]
    Code,
    /// The paper's (Sec. 5.1): cube side = extent, band width 0.2 (`±0.1`).
    Paper,
}

/// Options of the query sampling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QueryOptions {
    /// Queries of each kind per row (40; 120 in all).
    pub per_kind: usize,
    /// Box sizes.
    pub boxes: QueryBoxes,
    /// Minimum extent of the regularized face-selection box (0.1).
    pub min_extent: f64,
    /// Minimum cube side before scaling (0.01).
    pub min_side: f64,
    /// Precision of the ground truth's sign test.
    pub gwn: GwnPrecision,
}

impl Default for QueryOptions {
    fn default() -> Self {
        QueryOptions { per_kind: 40, boxes: QueryBoxes::Code, min_extent: 0.1, min_side: 0.01, gwn: GwnPrecision::F64 }
    }
}

impl QueryOptions {
    /// Narrow-band half-width.
    pub fn band(&self) -> f64 {
        match self.boxes {
            QueryBoxes::Code => 0.2,
            QueryBoxes::Paper => 0.1,
        }
    }

    /// Cube half-side factor `c` (`h = c·side`).
    pub fn cube_factor(&self) -> f64 {
        match self.boxes {
            QueryBoxes::Code => 0.5 * 3f64.sqrt(),
            QueryBoxes::Paper => 0.5,
        }
    }

    /// Queries per row.
    pub fn per_row(&self) -> usize {
        3 * self.per_kind
    }
}

/// Queries and ground truth of a set of rows: `rows.len() × per_row` points, row major.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RowQueries {
    /// The rows (indices into the shape's points).
    pub rows: Vec<u32>,
    /// Query positions.
    pub points: Vec<[f64; 3]>,
    /// Signed distances to the mesh.
    pub truth: Vec<f64>,
    /// Rows where no face met the regularized box: their zero-set and band queries were drawn from
    /// all faces instead (the reference draws face 0 then).
    pub rows_without_faces: usize,
}

/// The ChaCha8 stream of one (seed, epoch, shape, row) (D22): the key is SplitMix64 of
/// `(seed, epoch, shape)`, the stream number the row.
pub fn row_rng(seed: u64, epoch: u64, shape: u64, row: u64) -> ChaCha8Rng {
    let mut x = seed ^ 0x5041_545f_5452_4149; // "PAT_TRAI"
    let mut next = |v: u64| {
        x = x.wrapping_add(v).wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = x;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    let mut key = [0u8; 32];
    let words = [next(epoch), next(shape), next(0x51), next(0x52)];
    for (i, w) in words.iter().enumerate() {
        key[8 * i..8 * i + 8].copy_from_slice(&w.to_le_bytes());
    }
    let mut r = ChaCha8Rng::from_seed(key);
    r.set_stream(row);
    r
}

/// Uniform in `[0, 1)` with 53 random bits.
#[inline]
pub fn uniform(r: &mut ChaCha8Rng) -> f64 {
    (r.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

/// Per-face data for sampling.
struct Faces {
    lo: Vec<[f64; 3]>,
    hi: Vec<[f64; 3]>,
    area: Vec<f64>,
    normal: Vec<[f64; 3]>,
    tri: Vec<[[f64; 3]; 3]>,
}

fn faces_of(vertices: &[[f64; 3]], faces: &[[u32; 3]]) -> Faces {
    let mut f = Faces { lo: vec![], hi: vec![], area: vec![], normal: vec![], tri: vec![] };
    for t in faces {
        let [a, b, c] = t.map(|j| vertices[j as usize]);
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        f.area.push(0.5 * l);
        let ls = l.max(1e-10);
        f.normal.push([n[0] / ls, n[1] / ls, n[2] / ls]);
        f.lo.push(std::array::from_fn(|k| a[k].min(b[k]).min(c[k])));
        f.hi.push(std::array::from_fn(|k| a[k].max(b[k]).max(c[k])));
        f.tri.push([a, b, c]);
    }
    f
}

/// Draw a face index from the cumulative areas `cum` (∝ area), then a uniform point on it with the
/// reference's barycentrics `(√u(1−v), √u·v, 1 − √u)`.
///
/// [formula: pat-triangle-sampling]
fn sample_on(f: &Faces, sel: &[u32], cum: &[f64], r: &mut ChaCha8Rng) -> (u32, [f64; 3]) {
    let total = *cum.last().unwrap_or(&0.0);
    let x = uniform(r) * total;
    let k = cum.partition_point(|&c| c <= x).min(sel.len() - 1);
    let fi = sel[k];
    let (u, v) = (uniform(r), uniform(r));
    let su = u.sqrt();
    let (al, be) = (su * (1.0 - v), v * su);
    let ga = 1.0 - al - be;
    let [a, b, c] = f.tri[fi as usize];
    (fi, std::array::from_fn(|i| al * a[i] + be * b[i] + ga * c[i]))
}

/// The bounding box of a row's neighbor positions (`cols` valid columns).
pub fn row_box(positions: &[[f64; 3]], cols: &[usize]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &j in cols {
        for a in 0..3 {
            lo[a] = lo[a].min(positions[j][a]);
            hi[a] = hi[a].max(positions[j][a]);
        }
    }
    (lo, hi)
}

/// Sample the queries of one row. `cols` are the row's neighbor positions (indices into
/// `positions`).
///
/// [formula: pat-query-sampling]
fn sample_row(f: &Faces, all: &[u32], all_cum: &[f64], positions: &[[f64; 3]], cols: &[usize], opts: &QueryOptions, r: &mut ChaCha8Rng, out: &mut Vec<[f64; 3]>) -> bool {
    let (lo, hi) = row_box(positions, cols);
    let center: [f64; 3] = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
    let ext: [f64; 3] = std::array::from_fn(|a| (hi[a] - lo[a]).abs());
    let rlo: [f64; 3] = std::array::from_fn(|a| center[a] - 0.5 * ext[a].max(opts.min_extent));
    let rhi: [f64; 3] = std::array::from_fn(|a| center[a] + 0.5 * ext[a].max(opts.min_extent));
    let mut sel = Vec::new();
    let mut cum = Vec::new();
    let mut acc = 0.0;
    for (i, (flo, fhi)) in f.lo.iter().zip(&f.hi).enumerate() {
        if (0..3).all(|a| flo[a] <= rhi[a] && fhi[a] >= rlo[a]) && f.area[i] > 0.0 {
            sel.push(i as u32);
            acc += f.area[i];
            cum.push(acc);
        }
    }
    let found = !sel.is_empty();
    let (sel, cum) = if found { (&sel[..], &cum[..]) } else { (all, all_cum) };
    for _ in 0..opts.per_kind {
        out.push(sample_on(f, sel, cum, r).1);
    }
    let band = opts.band();
    for _ in 0..opts.per_kind {
        let (fi, p) = sample_on(f, sel, cum, r);
        let o = -band + 2.0 * band * uniform(r);
        let n = f.normal[fi as usize];
        out.push([p[0] + o * n[0], p[1] + o * n[1], p[2] + o * n[2]]);
    }
    let side = ext[0].max(ext[1]).max(ext[2]).max(opts.min_side);
    let h = opts.cube_factor() * side;
    for _ in 0..opts.per_kind {
        let p: [f64; 3] = std::array::from_fn(|a| center[a] - h + 2.0 * h * uniform(r));
        out.push(p);
    }
    found
}

/// Sample queries and ground truth for `rows` of a shape. `columns(i)` lists the positions the box
/// of row `i` is taken over. Runs on the current rayon pool; the result does not depend on it.
#[allow(clippy::too_many_arguments)]
pub fn sample_rows(
    positions: &[[f64; 3]],
    rows: &[u32],
    columns: &(dyn Fn(usize) -> Vec<usize> + Sync),
    mesh: &MeshTruth,
    vertices: &[[f64; 3]],
    faces: &[[u32; 3]],
    opts: &QueryOptions,
    key: (u64, u64, u64),
) -> Result<RowQueries, PatError> {
    if opts.per_kind == 0 {
        return Err(PatError::InvalidOption { what: "QueryOptions::per_kind = 0".into() });
    }
    let f = faces_of(vertices, faces);
    let all: Vec<u32> = (0..f.area.len() as u32).filter(|&i| f.area[i as usize] > 0.0).collect();
    if all.is_empty() {
        return Err(PatError::Format { reason: "the shape's mesh has no face of positive area".into() });
    }
    let mut acc = 0.0;
    let all_cum: Vec<f64> = all
        .iter()
        .map(|&i| {
            acc += f.area[i as usize];
            acc
        })
        .collect();
    let m = opts.per_row();
    let per: Vec<(Vec<[f64; 3]>, Vec<f64>, bool)> = rows
        .par_iter()
        .map(|&row| {
            let mut r = row_rng(key.0, key.1, key.2, row as u64);
            let mut pts = Vec::with_capacity(m);
            let found = sample_row(&f, &all, &all_cum, positions, &columns(row as usize), opts, &mut r, &mut pts);
            let truth = pts.iter().map(|&p| mesh.signed_distance(p, opts.gwn)).collect();
            (pts, truth, found)
        })
        .collect();
    let mut out = RowQueries { rows: rows.to_vec(), ..Default::default() };
    for (p, t, found) in per {
        out.points.extend(p);
        out.truth.extend(t);
        out.rows_without_faces += (!found) as usize;
    }
    Ok(out)
}
