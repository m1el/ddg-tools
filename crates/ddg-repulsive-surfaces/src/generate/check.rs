//! Checks used by the generator tests and available to callers: self-intersections of a
//! triangle mesh and the linking number of closed polylines.

use std::f64::consts::PI;

use super::Mesh;
use crate::geometry::v3::{cross, dot, norm, sub};

type P = [f64; 3];

/// Does the segment `a..b` pierce the triangle (proper crossing, strict interior)?
fn segment_triangle(a: P, b: P, t: [P; 3]) -> bool {
    let e1 = sub(t[1], t[0]);
    let e2 = sub(t[2], t[0]);
    let d = sub(b, a);
    let pv = cross(d, e2);
    let det = dot(e1, pv);
    if det.abs() < 1e-14 * norm(e1) * norm(e2) * norm(d) {
        return false; // parallel or coplanar
    }
    let tv = sub(a, t[0]);
    let u = dot(tv, pv) / det;
    let qv = cross(tv, e1);
    let v = dot(d, qv) / det;
    let s = dot(e2, qv) / det;
    let eps = 1e-12;
    u > eps && v > eps && u + v < 1.0 - eps && s > eps && s < 1.0 - eps
}

/// Do the two triangles cross (an edge of one pierces the other; coplanar overlaps are not
/// detected)?
pub fn triangles_cross(a: [P; 3], b: [P; 3]) -> bool {
    (0..3).any(|i| segment_triangle(a[i], a[(i + 1) % 3], b))
        || (0..3).any(|i| segment_triangle(b[i], b[(i + 1) % 3], a))
}

/// Pairs of faces that intersect each other (faces sharing a vertex are not compared), at most
/// `limit` of them. Sweep over the `x` extents.
pub fn self_intersections(mesh: &Mesh, limit: usize) -> Vec<[u32; 2]> {
    let tri = |f: usize| mesh.faces[f].map(|i| mesh.positions[i as usize]);
    let boxes: Vec<([f64; 3], [f64; 3])> = (0..mesh.faces.len())
        .map(|f| {
            let t = tri(f);
            let mut lo = t[0];
            let mut hi = t[0];
            for p in &t {
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            (lo, hi)
        })
        .collect();
    let mut order: Vec<usize> = (0..mesh.faces.len()).collect();
    order.sort_by(|&a, &b| boxes[a].0[0].total_cmp(&boxes[b].0[0]));
    let mut out = Vec::new();
    for (pi, &i) in order.iter().enumerate() {
        for &j in &order[pi + 1..] {
            if boxes[j].0[0] > boxes[i].1[0] {
                break;
            }
            if (1..3).any(|k| boxes[j].0[k] > boxes[i].1[k] || boxes[i].0[k] > boxes[j].1[k]) {
                continue;
            }
            let (fi, fj) = (mesh.faces[i], mesh.faces[j]);
            if fi.iter().any(|v| fj.contains(v)) {
                continue;
            }
            if triangles_cross(tri(i), tri(j)) {
                out.push([i.min(j) as u32, i.max(j) as u32]);
                if out.len() >= limit {
                    return out;
                }
            }
        }
    }
    out
}

/// The Gauss linking number of two closed polylines (midpoint rule; accurate when the curves are
/// well separated compared with their segment lengths). Round to the nearest integer.
pub fn linking_number(a: &[P], b: &[P]) -> f64 {
    let seg = |c: &[P], i: usize| {
        let p = c[i];
        let q = c[(i + 1) % c.len()];
        ([0, 1, 2].map(|k| 0.5 * (p[k] + q[k])), sub(q, p))
    };
    let mut s = 0.0;
    for i in 0..a.len() {
        let (ma, da) = seg(a, i);
        for j in 0..b.len() {
            let (mb, db) = seg(b, j);
            let r = sub(ma, mb);
            let n = norm(r);
            s += dot(r, cross(da, db)) / (n * n * n);
        }
    }
    s / (4.0 * PI)
}

/// Smallest distance between two sample points of closed polylines that are at least `exclude`
/// indices apart along the curve (`a` and `b` the same slice: self distance); a lower bound of
/// nothing in particular, an estimate of the clearance of the curves.
pub fn clearance(a: &[P], b: &[P], same: bool, exclude: usize) -> f64 {
    let mut m = f64::INFINITY;
    for i in 0..a.len() {
        for j in 0..b.len() {
            if same {
                let d = i.abs_diff(j);
                if d.min(a.len() - d) <= exclude {
                    continue;
                }
            }
            m = m.min(norm(sub(a[i], b[j])));
        }
    }
    m
}
