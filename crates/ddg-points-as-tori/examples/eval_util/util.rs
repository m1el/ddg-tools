//! Shared helpers of the WP14 examples `eval` and `density` (SPEC.md §16): `.npy` I/O, grids,
//! error metrics, the jump statistic, a small JSON writer, the tapered-kNN field (ours, for the
//! continuity study only) and the self-intersection test used to select Thingi10K meshes.

#![allow(dead_code)]

use ddg_points_as_tori::csg::ImplicitField;
use ddg_points_as_tori::npy::read_npy;
use ddg_points_as_tori::sdf::torus_sdf;
use ddg_points_as_tori::spatial::Neighbor;
use ddg_points_as_tori::truth::MeshTruth;
use ddg_points_as_tori::{PatError, TorusField};
use rayon::prelude::*;
use std::fmt::Write as _;
use std::path::Path;

pub type R<T> = Result<T, String>;

pub fn err<E: std::fmt::Display>(e: E) -> String {
    e.to_string()
}

// ------------------------------------------------------------------ npy

fn npy_header(descr: &str, shape: &[usize]) -> Vec<u8> {
    let sh = match shape.len() {
        1 => format!("({},)", shape[0]),
        _ => format!("({})", shape.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", ")),
    };
    let mut h = format!("{{'descr': '{descr}', 'fortran_order': False, 'shape': {sh}, }}");
    let total = 10 + h.len() + 1;
    let pad = (64 - total % 64) % 64;
    h.push_str(&" ".repeat(pad));
    h.push('\n');
    let mut out = b"\x93NUMPY\x01\x00".to_vec();
    out.extend_from_slice(&(h.len() as u16).to_le_bytes());
    out.extend_from_slice(h.as_bytes());
    out
}

pub fn write_npy_f64(path: &Path, shape: &[usize], data: &[f64]) -> R<()> {
    let mut b = npy_header("<f8", shape);
    for v in data {
        b.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(path, b).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn write_npy_f32(path: &Path, shape: &[usize], data: &[f64]) -> R<()> {
    let mut b = npy_header("<f4", shape);
    for v in data {
        b.extend_from_slice(&(*v as f32).to_le_bytes());
    }
    std::fs::write(path, b).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn write_npy_i32(path: &Path, shape: &[usize], data: &[i64]) -> R<()> {
    let mut b = npy_header("<i4", shape);
    for v in data {
        b.extend_from_slice(&(*v as i32).to_le_bytes());
    }
    std::fs::write(path, b).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn write3(path: &Path, v: &[[f64; 3]]) -> R<()> {
    let flat: Vec<f64> = v.iter().flatten().copied().collect();
    write_npy_f64(path, &[v.len(), 3], &flat)
}

pub fn write_faces(path: &Path, f: &[[u32; 3]]) -> R<()> {
    let flat: Vec<i64> = f.iter().flatten().map(|&x| x as i64).collect();
    write_npy_i32(path, &[f.len(), 3], &flat)
}

pub fn read_flat(path: &Path) -> R<(Vec<usize>, Vec<f64>)> {
    let a = read_npy(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok((a.shape, a.data))
}

pub fn read3(path: &Path) -> R<Vec<[f64; 3]>> {
    let (s, d) = read_flat(path)?;
    if s.len() != 2 || s[1] != 3 {
        return Err(format!("{}: expected shape (n, 3), got {s:?}", path.display()));
    }
    Ok(d.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect())
}

pub fn read_faces(path: &Path) -> R<Vec<[u32; 3]>> {
    let (s, d) = read_flat(path)?;
    if s.len() != 2 || s[1] != 3 {
        return Err(format!("{}: expected shape (m, 3), got {s:?}", path.display()));
    }
    d.chunks_exact(3)
        .enumerate()
        .map(|(i, c)| {
            if c.iter().all(|&x| x >= 0.0 && x < u32::MAX as f64 && x.fract() == 0.0) {
                Ok([c[0] as u32, c[1] as u32, c[2] as u32])
            } else {
                Err(format!("{}: face {i} has invalid indices {c:?}", path.display()))
            }
        })
        .collect()
}

// ------------------------------------------------------------------ grids

/// `linspace(-1, 1, n)` (the golden grid axis for n = 64).
pub fn axis(n: usize) -> Vec<f64> {
    (0..n).map(|i| -1.0 + 2.0 * i as f64 / (n - 1) as f64).collect()
}

/// The `n³` grid over `[-1, 1]³` in C order of `(ix, iy, iz)` (`meshgrid(..., indexing='ij')`).
pub fn grid(n: usize) -> Vec<[f64; 3]> {
    let a = axis(n);
    let mut g = Vec::with_capacity(n * n * n);
    for &x in &a {
        for &y in &a {
            for &z in &a {
                g.push([x, y, z]);
            }
        }
    }
    g
}

#[inline]
pub fn gi(n: usize, i: usize, j: usize, k: usize) -> usize {
    (i * n + j) * n + k
}

// ------------------------------------------------------------------ metrics

/// Error metrics of one method on one grid.
#[derive(Debug, Clone, Default)]
pub struct Metrics {
    pub mae: f64,
    pub max: f64,
    pub band_mae: f64,
    pub band_n: usize,
    pub band2_mae: f64,
    pub eik: f64,
    pub eik_band: f64,
    pub sign_err: f64,
    pub nonfinite: usize,
}

pub const BAND: f64 = 0.05;
pub const BAND2: f64 = 0.1;
pub const SIGN_MARGIN: f64 = 0.05;

/// Grid MAE, max, narrow-band MAE (`|φ_true| < BAND`, `< BAND2`), eikonality by central
/// differences (mean `| ‖∇_h φ‖ − 1 |` over interior grid points; also over the band), and the
/// fraction of grid points with `|φ_true| > SIGN_MARGIN` whose sign is wrong.
pub fn metrics(n: usize, phi: &[f64], truth: &[f64]) -> Metrics {
    let h = 2.0 / (n - 1) as f64;
    let mut m = Metrics::default();
    let (mut s, mut sb, mut nb, mut sb2, mut nb2, mut se, mut ns) = (0.0, 0.0, 0usize, 0.0, 0usize, 0usize, 0usize);
    for (a, b) in phi.iter().zip(truth) {
        if !a.is_finite() {
            m.nonfinite += 1;
            continue;
        }
        let e = (a - b).abs();
        s += e;
        m.max = m.max.max(e);
        if b.abs() < BAND {
            sb += e;
            nb += 1;
        }
        if b.abs() < BAND2 {
            sb2 += e;
            nb2 += 1;
        }
        if b.abs() > SIGN_MARGIN {
            ns += 1;
            if (a < &0.0) != (b < &0.0) {
                se += 1;
            }
        }
    }
    let q = phi.len() - m.nonfinite;
    m.mae = s / q.max(1) as f64;
    m.band_mae = sb / nb.max(1) as f64;
    m.band_n = nb;
    m.band2_mae = sb2 / nb2.max(1) as f64;
    m.sign_err = se as f64 / ns.max(1) as f64;
    let (eik, eikb) = eikonality(n, h, phi, truth);
    m.eik = eik;
    m.eik_band = eikb;
    m
}

pub fn eikonality(n: usize, h: f64, phi: &[f64], truth: &[f64]) -> (f64, f64) {
    let rows: Vec<(f64, usize, f64, usize)> = (1..n - 1)
        .into_par_iter()
        .map(|i| {
            let (mut s, mut c, mut sb, mut cb) = (0.0, 0, 0.0, 0);
            for j in 1..n - 1 {
                for k in 1..n - 1 {
                    let gx = phi[gi(n, i + 1, j, k)] - phi[gi(n, i - 1, j, k)];
                    let gy = phi[gi(n, i, j + 1, k)] - phi[gi(n, i, j - 1, k)];
                    let gz = phi[gi(n, i, j, k + 1)] - phi[gi(n, i, j, k - 1)];
                    let g = (gx * gx + gy * gy + gz * gz).sqrt() / (2.0 * h);
                    if !g.is_finite() {
                        continue;
                    }
                    let e = (g - 1.0).abs();
                    s += e;
                    c += 1;
                    if truth[gi(n, i, j, k)].abs() < BAND {
                        sb += e;
                        cb += 1;
                    }
                }
            }
            (s, c, sb, cb)
        })
        .collect();
    let (s, c, sb, cb) = rows.iter().fold((0.0, 0, 0.0, 0), |a, r| (a.0 + r.0, a.1 + r.1, a.2 + r.2, a.3 + r.3));
    (s / c.max(1) as f64, sb / cb.max(1) as f64)
}

/// Jump statistic: `|Δφ|/h` along every grid edge (all three axes).
#[derive(Debug, Clone, Default)]
pub struct Jumps {
    pub edges: usize,
    pub p50: f64,
    pub p90: f64,
    pub p99: f64,
    pub p999: f64,
    pub max: f64,
    pub over2: usize,
    pub over5: usize,
    pub over20: usize,
    /// The largest edges: (slope, grid index a, axis).
    pub top: Vec<(f64, usize, usize)>,
}

pub fn jumps(n: usize, phi: &[f64], keep_top: usize) -> Jumps {
    let h = 2.0 / (n - 1) as f64;
    let mut s: Vec<(f64, usize, usize)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            let mut v = Vec::with_capacity(3 * n * n);
            for j in 0..n {
                for k in 0..n {
                    let a = gi(n, i, j, k);
                    if i + 1 < n {
                        v.push(((phi[gi(n, i + 1, j, k)] - phi[a]).abs() / h, a, 0));
                    }
                    if j + 1 < n {
                        v.push(((phi[gi(n, i, j + 1, k)] - phi[a]).abs() / h, a, 1));
                    }
                    if k + 1 < n {
                        v.push(((phi[gi(n, i, j, k + 1)] - phi[a]).abs() / h, a, 2));
                    }
                }
            }
            v
        })
        .collect();
    let m = s.len();
    let q = |s: &[(f64, usize, usize)], p: f64| s[((p * (m - 1) as f64).round() as usize).min(m - 1)].0;
    s.par_sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut j = Jumps {
        edges: m,
        p50: q(&s, 0.5),
        p90: q(&s, 0.9),
        p99: q(&s, 0.99),
        p999: q(&s, 0.999),
        max: s[m - 1].0,
        over2: s.iter().filter(|x| x.0 > 2.0).count(),
        over5: s.iter().filter(|x| x.0 > 5.0).count(),
        over20: s.iter().filter(|x| x.0 > 20.0).count(),
        top: Vec::new(),
    };
    j.top = s.iter().rev().take(keep_top).copied().collect();
    j
}

pub fn quantile(v: &mut [f64], p: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.total_cmp(b));
    v[((p * (v.len() - 1) as f64).round() as usize).min(v.len() - 1)]
}

// ------------------------------------------------------------------ JSON

/// A minimal JSON object writer (no dependency).
#[derive(Default)]
pub struct Json {
    s: String,
    first: bool,
}

pub fn jnum(x: f64) -> String {
    if x.is_finite() { format!("{x:.6e}") } else { "null".into() }
}

pub fn jstr(s: &str) -> String {
    let mut o = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c if (c as u32) < 0x20 => {
                let _ = write!(o, "\\u{:04x}", c as u32);
            }
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

impl Json {
    pub fn new() -> Json {
        Json { s: "{".into(), first: true }
    }
    fn key(&mut self, k: &str) {
        if !self.first {
            self.s.push_str(", ");
        }
        self.first = false;
        self.s.push_str(&jstr(k));
        self.s.push_str(": ");
    }
    pub fn num(&mut self, k: &str, v: f64) -> &mut Self {
        self.key(k);
        self.s.push_str(&jnum(v));
        self
    }
    pub fn int(&mut self, k: &str, v: i64) -> &mut Self {
        self.key(k);
        let _ = write!(self.s, "{v}");
        self
    }
    pub fn str(&mut self, k: &str, v: &str) -> &mut Self {
        self.key(k);
        self.s.push_str(&jstr(v));
        self
    }
    pub fn raw(&mut self, k: &str, v: &str) -> &mut Self {
        self.key(k);
        self.s.push_str(v);
        self
    }
    pub fn bool(&mut self, k: &str, v: bool) -> &mut Self {
        self.key(k);
        self.s.push_str(if v { "true" } else { "false" });
        self
    }
    pub fn finish(&self) -> String {
        format!("{}}}", self.s)
    }
}

pub fn jarr(items: &[String]) -> String {
    format!("[{}]", items.join(", "))
}

pub fn metrics_json(m: &Metrics) -> String {
    Json::new()
        .num("mae", m.mae)
        .num("max", m.max)
        .num("band_mae", m.band_mae)
        .int("band_n", m.band_n as i64)
        .num("band2_mae", m.band2_mae)
        .num("eik", m.eik)
        .num("eik_band", m.eik_band)
        .num("sign_err", m.sign_err)
        .int("nonfinite", m.nonfinite as i64)
        .finish()
}

pub fn jumps_json(j: &Jumps) -> String {
    Json::new()
        .int("edges", j.edges as i64)
        .num("p50", j.p50)
        .num("p90", j.p90)
        .num("p99", j.p99)
        .num("p999", j.p999)
        .num("max", j.max)
        .int("over2", j.over2 as i64)
        .int("over5", j.over5 as i64)
        .int("over20", j.over20 as i64)
        .finish()
}

// ------------------------------------------------------------------ fields

/// The mesh truth as an [`ImplicitField`] (for marching cubes of the truth).
pub struct TruthField<'a>(pub &'a MeshTruth);

impl ImplicitField for TruthField<'_> {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        self.0.signed_distance(x)
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let d = self.0.signed_distance(x)?;
        Ok((d, [0.0; 3]))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.1; 3], [1.1; 3]]
    }
}

/// The tapered-kNN blend (ours; WP6 resolution 1, an option for the continuity comparison only).
///
/// For a query `x`: the `k + 1` nearest points (`ρ_1 ≤ … ≤ ρ_{k+1}`), `σ = ½ ρ_k` (as the default
/// rule: half the largest distance in the k-set), `λ = C/σ`, weights
/// `w_i = exp(−λ(ρ_i − σ)) · τ(ρ_i / ρ_{k+1})` over the `k` nearest with `τ(t) = (1 − t²)²`.
/// A point enters or leaves the k-set only when its distance equals `ρ_{k+1}`, where `τ = 0`, so
/// `φ = Σ w g / Σ w` is continuous wherever `ρ_{k+1} > ρ_1` (`ρ_k` and `ρ_{k+1}` are continuous).
/// Near the surface (`ρ_i ≪ ρ_{k+1}`) `τ ≈ 1` and the blend equals the default one.
pub struct TaperField<'a> {
    pub field: &'a TorusField,
    pub k: usize,
    pub c: f64,
}

impl TaperField<'_> {
    pub fn eval(&self, x: [f64; 3], buf: &mut Vec<Neighbor>) -> Result<f64, PatError> {
        let pos = self.field.positions();
        let tori = &self.field.tori().tori;
        let kk = (self.k + 1).min(pos.len());
        self.field.tree().knn(x, kk, buf);
        let rk1 = buf[kk - 1].d2.sqrt();
        let used = if kk > self.k { self.k } else { kk };
        let sigma = 0.5 * buf[used - 1].d2.sqrt();
        if !(sigma > 0.0) || !(rk1 > 0.0) {
            return Err(PatError::DegenerateQuery { index: 0, reason: "tapered kNN: all neighbors at the query".into() });
        }
        let lambda = self.c / sigma;
        let (mut ws, mut gs) = (0.0, 0.0);
        for nb in &buf[..used] {
            let rho = nb.d2.sqrt();
            let t = rho / rk1;
            let tau = if t < 1.0 { (1.0 - t * t) * (1.0 - t * t) } else { 0.0 };
            let w = (-lambda * (rho - sigma)).exp() * tau;
            ws += w;
            gs += w * torus_sdf(&tori[nb.index as usize], x);
        }
        if !(ws > 0.0) {
            // every weight tapered to 0 (all k at distance ρ_{k+1}): fall back to the plain mean
            let mut s = 0.0;
            for nb in &buf[..used] {
                s += torus_sdf(&tori[nb.index as usize], x);
            }
            return Ok(s / used as f64);
        }
        Ok(gs / ws)
    }

    pub fn evaluate(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        let r: Result<Vec<()>, PatError> = xs
            .par_chunks(512)
            .zip(out.par_chunks_mut(512))
            .map(|(x, o)| {
                let mut buf = Vec::new();
                for (xi, oi) in x.iter().zip(o.iter_mut()) {
                    *oi = self.eval(*xi, &mut buf)?;
                }
                Ok(())
            })
            .collect();
        r.map(|_| ())
    }
}

impl ImplicitField for TaperField<'_> {
    fn distance(&self, x: [f64; 3]) -> Result<f64, PatError> {
        self.eval(x, &mut Vec::new())
    }
    fn distance_gradient(&self, x: [f64; 3]) -> Result<(f64, [f64; 3]), PatError> {
        let h = 1e-6;
        let d = self.distance(x)?;
        let mut g = [0.0; 3];
        for a in 0..3 {
            let mut p = x;
            let mut m = x;
            p[a] += h;
            m[a] -= h;
            g[a] = (self.distance(p)? - self.distance(m)?) / (2.0 * h);
        }
        Ok((d, g))
    }
    fn bbox(&self) -> [[f64; 3]; 2] {
        [[-1.1; 3], [1.1; 3]]
    }
    fn distances(&self, xs: &[[f64; 3]], out: &mut [f64]) -> Result<(), PatError> {
        self.evaluate(xs, out)
    }
}

// ------------------------------------------------------------------ self-intersection

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn orient(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    dot(sub(b, a), cross(sub(c, a), sub(d, a)))
}

/// Does the open segment `pq` cross the triangle `abc` (strictly: both endpoints strictly on
/// opposite sides of the plane, and the crossing strictly inside the triangle)? Determinants with
/// a magnitude below `eps` count as zero (no crossing): touching contacts are not reported.
fn segment_crosses(p: [f64; 3], q: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3], eps: f64) -> bool {
    let sp = orient(a, b, c, p);
    let sq = orient(a, b, c, q);
    if !(sp > eps && sq < -eps || sp < -eps && sq > eps) {
        return false;
    }
    let s1 = orient(p, q, a, b);
    let s2 = orient(p, q, b, c);
    let s3 = orient(p, q, c, a);
    (s1 > eps && s2 > eps && s3 > eps) || (s1 < -eps && s2 < -eps && s3 < -eps)
}

/// Count pairs of faces that intersect other than in their shared vertices/edge (proper
/// crossings of an edge of one through the interior of the other), and folded edge pairs
/// (two faces on an edge with opposite normals, dihedral ≈ 0). Candidate pairs come from a
/// uniform grid on the triangle bounding boxes. Coplanar overlaps of non-adjacent faces are not
/// detected by this test (counted separately as `coplanar_close` when their boxes overlap and
/// they are coplanar within `eps`).
pub struct SelfIntersection {
    pub crossings: usize,
    pub folds: usize,
    pub coplanar_close: usize,
    pub pairs_tested: usize,
}

pub fn self_intersections(v: &[[f64; 3]], f: &[[u32; 3]]) -> SelfIntersection {
    let m = f.len();
    let tri = |i: usize| [v[f[i][0] as usize], v[f[i][1] as usize], v[f[i][2] as usize]];
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in v {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    let ext = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max).max(1e-300);
    let scale = ext;
    let eps = 1e-14 * scale * scale * scale;
    let cells = ((m as f64).cbrt() * 2.0).clamp(1.0, 256.0) as usize;
    let cs = ext / cells as f64 * 1.000001;
    let cell_of = |x: f64, a: usize| (((x - lo[a]) / cs).floor().max(0.0) as usize).min(cells - 1);
    let mut buckets: std::collections::HashMap<(usize, usize, usize), Vec<u32>> = std::collections::HashMap::new();
    let boxes: Vec<([f64; 3], [f64; 3])> = (0..m)
        .map(|i| {
            let t = tri(i);
            let mut l = [f64::INFINITY; 3];
            let mut h = [f64::NEG_INFINITY; 3];
            for p in t {
                for a in 0..3 {
                    l[a] = l[a].min(p[a]);
                    h[a] = h[a].max(p[a]);
                }
            }
            (l, h)
        })
        .collect();
    for (i, (l, h)) in boxes.iter().enumerate() {
        for x in cell_of(l[0], 0)..=cell_of(h[0], 0) {
            for y in cell_of(l[1], 1)..=cell_of(h[1], 1) {
                for z in cell_of(l[2], 2)..=cell_of(h[2], 2) {
                    buckets.entry((x, y, z)).or_default().push(i as u32);
                }
            }
        }
    }
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    for b in buckets.values() {
        for (s, &i) in b.iter().enumerate() {
            for &j in &b[s + 1..] {
                let (a, c) = if i < j { (i, j) } else { (j, i) };
                let (la, ha) = boxes[a as usize];
                let (lb, hb) = boxes[c as usize];
                if (0..3).all(|k| la[k] <= hb[k] && lb[k] <= ha[k]) {
                    pairs.push((a, c));
                }
            }
        }
    }
    pairs.par_sort_unstable();
    pairs.dedup();
    let res: Vec<(usize, usize, usize)> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let (fi, fj) = (f[i as usize], f[j as usize]);
            let shared: Vec<u32> = fi.iter().filter(|x| fj.contains(x)).copied().collect();
            let (ti, tj) = (tri(i as usize), tri(j as usize));
            let ni = cross(sub(ti[1], ti[0]), sub(ti[2], ti[0]));
            let nj = cross(sub(tj[1], tj[0]), sub(tj[2], tj[0]));
            let li = dot(ni, ni).sqrt();
            let lj = dot(nj, nj).sqrt();
            match shared.len() {
                0 => {
                    let mut hit = false;
                    for e in 0..3 {
                        hit |= segment_crosses(ti[e], ti[(e + 1) % 3], tj[0], tj[1], tj[2], eps);
                        hit |= segment_crosses(tj[e], tj[(e + 1) % 3], ti[0], ti[1], ti[2], eps);
                    }
                    let coplanar = li > 0.0 && lj > 0.0 && (0..3).all(|k| orient(ti[0], ti[1], ti[2], tj[k]).abs() <= eps);
                    (hit as usize, 0, (coplanar && !hit) as usize)
                }
                1 => {
                    // the edges opposite the shared vertex
                    let s = shared[0];
                    let opp = |t: [u32; 3], p: [[f64; 3]; 3]| {
                        let k = t.iter().position(|&x| x == s).unwrap_or(0);
                        (p[(k + 1) % 3], p[(k + 2) % 3])
                    };
                    let (a, b) = opp(fi, ti);
                    let (c, d) = opp(fj, tj);
                    let hit = segment_crosses(a, b, tj[0], tj[1], tj[2], eps) || segment_crosses(c, d, ti[0], ti[1], ti[2], eps);
                    (hit as usize, 0, 0)
                }
                2 => {
                    let fold = li > 0.0 && lj > 0.0 && dot(ni, nj) / (li * lj) < -1.0 + 1e-9;
                    (0, fold as usize, 0)
                }
                _ => (0, 0, 0),
            }
        })
        .collect();
    SelfIntersection {
        crossings: res.iter().map(|r| r.0).sum(),
        folds: res.iter().map(|r| r.1).sum(),
        coplanar_close: res.iter().map(|r| r.2).sum(),
        pairs_tested: pairs.len(),
    }
}
