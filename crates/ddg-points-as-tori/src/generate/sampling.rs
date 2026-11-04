//! Uniform area sampling of triangle meshes and farthest-point sampling (SPEC.md §10.3;
//! NOTES-REFERENCE §11: `sample_uniform_point_cloud`, `subsample_point_cloud`).

use super::rng::Rng;
use crate::PatError;
use rayon::prelude::*;

/// Translate by the vertex mean and scale by the largest absolute coordinate so that the mesh
/// fits `[-1, 1]³` (the reference's `center_and_scale`). A degenerate mesh is only centered.
pub fn center_and_scale(vertices: &[[f64; 3]]) -> Vec<[f64; 3]> {
    if vertices.is_empty() {
        return Vec::new();
    }
    let n = vertices.len() as f64;
    let mut c = [0.0; 3];
    for v in vertices {
        for a in 0..3 {
            c[a] += v[a] / n;
        }
    }
    let mut m: f64 = 0.0;
    for v in vertices {
        for a in 0..3 {
            m = m.max((v[a] - c[a]).abs());
        }
    }
    let s = if m > 1e-10 { 1.0 / m } else { 1.0 };
    vertices.iter().map(|v| [(v[0] - c[0]) * s, (v[1] - c[1]) * s, (v[2] - c[2]) * s]).collect()
}

/// `n` points uniform in area on a triangle mesh, with the unit normals of their faces (right-hand
/// rule). Draws sequentially from `rng`, so the result depends only on the generator state.
pub fn sample_mesh_uniform(vertices: &[[f64; 3]], faces: &[[u32; 3]], n: usize, rng: &mut Rng) -> Result<(Vec<[f64; 3]>, Vec<[f64; 3]>), PatError> {
    let mut cum = Vec::with_capacity(faces.len());
    let mut normals = Vec::with_capacity(faces.len());
    let mut total = 0.0;
    for (f, t) in faces.iter().enumerate() {
        if let Some(&v) = t.iter().find(|&&v| v as usize >= vertices.len()) {
            return Err(PatError::Format { reason: format!("face {f} has vertex index {v}, but the mesh has {} vertices", vertices.len()) });
        }
        let (a, b, c) = (vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]);
        let (u, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let x = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        total += 0.5 * l;
        cum.push(total);
        normals.push(if l > 0.0 { [x[0] / l, x[1] / l, x[2] / l] } else { [0.0; 3] });
    }
    if !(total > 0.0) || !total.is_finite() {
        return Err(PatError::Format { reason: "mesh has no area to sample".into() });
    }
    let mut pts = Vec::with_capacity(n);
    let mut nrm = Vec::with_capacity(n);
    for _ in 0..n {
        let target = rng.f64() * total;
        let f = cum.partition_point(|&c| c <= target).min(faces.len() - 1);
        let t = faces[f];
        let (a, b, c) = (vertices[t[0] as usize], vertices[t[1] as usize], vertices[t[2] as usize]);
        let (r1, r2) = (rng.f64(), rng.f64());
        let s = r1.sqrt();
        let (wa, wb, wc) = (1.0 - s, s * (1.0 - r2), s * r2);
        pts.push([0, 1, 2].map(|k| wa * a[k] + wb * b[k] + wc * c[k]));
        nrm.push(normals[f]);
    }
    Ok((pts, nrm))
}

/// Farthest-point sampling: `m` indices, starting from `start`, each next one maximizing the
/// distance to the chosen set (ties: smallest index). Returns all indices when `m ≥ n`. Runs on
/// the current rayon pool; the result is independent of the thread count.
pub fn farthest_point_sampling(points: &[[f64; 3]], m: usize, start: usize) -> Vec<usize> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    if m >= n {
        return (0..n).collect();
    }
    const CHUNK: usize = 4096;
    let mut mind = vec![f64::INFINITY; n];
    let mut chosen = Vec::with_capacity(m);
    let mut cur = start.min(n - 1);
    chosen.push(cur);
    while chosen.len() < m {
        let q = points[cur];
        let best: (f64, usize) = mind
            .par_chunks_mut(CHUNK)
            .enumerate()
            .map(|(c, chunk)| {
                let mut best = (-1.0f64, usize::MAX);
                for (i, d) in chunk.iter_mut().enumerate() {
                    let p = points[c * CHUNK + i];
                    let dd = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
                    if dd < *d {
                        *d = dd;
                    }
                    if *d > best.0 {
                        best = (*d, c * CHUNK + i);
                    }
                }
                best
            })
            .collect::<Vec<_>>()
            .into_iter()
            .fold((-1.0, usize::MAX), |a, b| if b.0 > a.0 { b } else { a });
        cur = best.1;
        chosen.push(cur);
        mind[cur] = 0.0;
    }
    chosen
}
