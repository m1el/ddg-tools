//! Barill point areas (Barill et al. 2018, Sec. 3.1; SPEC.md §9, D26).
//!
//! The paper: "we project each point and its k-nearest neighbors (filtering out those with
//! opposite normals) onto the best-fit plane, and compute that point's associated 2D Voronoi
//! area", with `k = 20`. Ours, written from that text (libigl's GPL `point_areas` was not
//! consulted):
//!
//! 1. The `k` nearest points of `p_i` other than itself, by the deterministic kd-tree
//!    (`(d², index)` order). Those with `⟨n_j, n_i⟩ < 0` are dropped when
//!    [`PointAreaOptions::filter_opposite`] (counted).
//! 2. They are projected to the tangent plane of `p_i` (normal `n_i`), giving `q_j ∈ R²` with
//!    `p_i` at the origin.
//! 3. The area is that of the Voronoi cell of the origin among `{0} ∪ {q_j}`, i.e. the
//!    intersection of the half-planes `⟨y, q_j⟩ ≤ ‖q_j‖²/2`. This is exactly the Voronoi cell of
//!    the center vertex in the planar Delaunay triangulation (its dual), computed by clipping a
//!    convex polygon, with no predicates. When the origin is on the convex hull of the `q_j`
//!    the cell is unbounded; every cell is therefore intersected with the square of half-width
//!    `h = max_j ‖q_j‖` (the projected k-NN radius) and the cells that touch the square are
//!    counted as `clipped`.
//! 4. Neighbors whose projection is exactly the origin (duplicates, or points straight along the
//!    normal) have no bisector. They are left out of the clipping, and the cell is shared among
//!    the coincident points: its area is divided by `1 + their number` (counted).

use super::{CHUNK, dot, sub};
use crate::spatial::{KdTree, KdTreeOptions};
use crate::{PatError, PointCloud};
use rayon::prelude::*;

/// Options of [`point_areas`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PointAreaOptions {
    /// Neighbors besides the point itself (Barill: 20).
    pub k: usize,
    /// Drop neighbors with opposite normals (`⟨n_j, n_i⟩ < 0`), as Barill.
    pub filter_opposite: bool,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
}

impl Default for PointAreaOptions {
    fn default() -> Self {
        PointAreaOptions { k: 20, filter_opposite: true, threads: None }
    }
}

/// What [`point_areas`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointAreaReport {
    /// Points.
    pub n: usize,
    /// Cells that touch the clipping square (unbounded or very elongated Voronoi cells).
    pub clipped: usize,
    /// Neighbors dropped for an opposite normal (summed over points).
    pub filtered: u64,
    /// Points that share their cell with coincident projected neighbors.
    pub coincident: usize,
    /// Points with no usable neighbor (all filtered or coincident): area 0.
    pub isolated: usize,
    /// `Σ A_i`.
    pub total_area: f64,
    /// Wall time.
    pub seconds: f64,
}

/// The area of the Voronoi cell of the origin among `{0} ∪ q` intersected with the square
/// `[−h, h]²`, and whether the square bounds the result (`q_j ≠ 0`, `h > 0`).
///
/// [formula: pat-point-area]
pub fn voronoi_cell_area(q: &[[f64; 2]], h: f64) -> (f64, bool) {
    // polygon vertices with the label of the edge that starts at each (usize::MAX: square side)
    const BOX: usize = usize::MAX;
    let mut poly: Vec<([f64; 2], usize)> = vec![([-h, -h], BOX), ([h, -h], BOX), ([h, h], BOX), ([-h, h], BOX)];
    let mut next = Vec::with_capacity(poly.len() + 4);
    for (j, &qj) in q.iter().enumerate() {
        let c = 0.5 * (qj[0] * qj[0] + qj[1] * qj[1]);
        let f = |v: [f64; 2]| v[0] * qj[0] + v[1] * qj[1] - c;
        next.clear();
        let m = poly.len();
        for t in 0..m {
            let (a, la) = poly[t];
            let (b, _) = poly[(t + 1) % m];
            let (fa, fb) = (f(a), f(b));
            let ain = fa <= 0.0;
            let bin = fb <= 0.0;
            if ain {
                next.push((a, la));
            }
            if ain != bin {
                let s = fa / (fa - fb);
                let p = [a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1])];
                // leaving: the new edge runs along the bisector j; entering: the rest of edge la
                next.push((p, if ain { j } else { la }));
            }
        }
        std::mem::swap(&mut poly, &mut next);
        if poly.len() < 3 {
            return (0.0, false);
        }
    }
    let m = poly.len();
    let mut a2 = 0.0;
    for t in 0..m {
        let (u, v) = (poly[t].0, poly[(t + 1) % m].0);
        a2 += u[0] * v[1] - u[1] * v[0];
    }
    (0.5 * a2.abs(), poly.iter().any(|&(_, l)| l == BOX))
}

/// A unit tangent basis `(s, t)` of the plane with unit normal `n` (any fixed rule; the area does
/// not depend on it).
fn tangent_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let a = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let d = dot(a, n);
    let s = [a[0] - d * n[0], a[1] - d * n[1], a[2] - d * n[2]];
    let l = dot(s, s).sqrt();
    let s = [s[0] / l, s[1] / l, s[2] / l];
    let t = super::cross(n, s);
    (s, t)
}

#[derive(Default, Clone, Copy)]
struct One {
    area: f64,
    clipped: bool,
    filtered: u32,
    coincident: bool,
    isolated: bool,
}

/// Barill point areas of `cloud` (one per point), on an explicit pool.
///
/// [formula: pat-point-area]
pub fn point_areas(cloud: &PointCloud, opts: &PointAreaOptions) -> Result<(Vec<f64>, PointAreaReport), PatError> {
    let t0 = std::time::Instant::now();
    if opts.k == 0 {
        return Err(PatError::InvalidOption { what: "PointAreaOptions::k = 0".into() });
    }
    let pool = crate::pool::make_pool(opts.threads)?;
    let (p, n) = (cloud.positions(), cloud.normals());
    let ones: Vec<One> = pool.install(|| -> Result<Vec<One>, PatError> {
        let tree = KdTree::new(p, KdTreeOptions::default())?;
        Ok((0..p.len())
            .into_par_iter()
            .with_min_len(CHUNK)
            .map_init(
                || (Vec::new(), Vec::new()),
                |(buf, q), i| {
                    tree.knn(p[i], opts.k + 1, buf);
                    let (s, t) = tangent_basis(n[i]);
                    q.clear();
                    let mut one = One::default();
                    let mut h: f64 = 0.0;
                    let mut same = 0usize;
                    let mut seen_self = false;
                    for nb in buf.iter() {
                        let j = nb.index as usize;
                        if j == i && !seen_self {
                            seen_self = true;
                            continue;
                        }
                        if q.len() + same + one.filtered as usize >= opts.k {
                            break;
                        }
                        if opts.filter_opposite && dot(n[j], n[i]) < 0.0 {
                            one.filtered += 1;
                            continue;
                        }
                        let d = sub(p[j], p[i]);
                        let qq = [dot(d, s), dot(d, t)];
                        if qq[0] == 0.0 && qq[1] == 0.0 {
                            same += 1;
                            continue;
                        }
                        h = h.max((qq[0] * qq[0] + qq[1] * qq[1]).sqrt());
                        q.push(qq);
                    }
                    if q.is_empty() || h == 0.0 {
                        one.isolated = true;
                        return one;
                    }
                    let (a, clipped) = voronoi_cell_area(q, h);
                    one.area = a / (1 + same) as f64;
                    one.clipped = clipped;
                    one.coincident = same > 0;
                    one
                },
            )
            .collect())
    })?;
    let mut rep = PointAreaReport { n: p.len(), ..Default::default() };
    let areas: Vec<f64> = ones.iter().map(|o| o.area).collect();
    for o in &ones {
        rep.clipped += o.clipped as usize;
        rep.filtered += o.filtered as u64;
        rep.coincident += o.coincident as usize;
        rep.isolated += o.isolated as usize;
    }
    rep.total_area = areas.iter().sum();
    rep.seconds = t0.elapsed().as_secs_f64();
    Ok((areas, rep))
}
