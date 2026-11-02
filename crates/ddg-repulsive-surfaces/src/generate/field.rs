//! Distance-like fields (negative inside) for the tube mesher: spheres, boxes, infinite
//! cylinders, tubes around polylines and curve networks, smooth union and difference.

use std::sync::Arc;

use crate::geometry::v3::{add, dot, norm, scale, sub};

/// A scalar field `f(x)`, negative inside; roughly 1-Lipschitz.
pub type Field = Arc<dyn Fn([f64; 3]) -> f64 + Send + Sync>;

/// Distance from `p` to the segment `ab`.
pub fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    let t = if l2 > 0.0 {
        (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    norm(sub(p, add(a, scale(ab, t))))
}

/// Distance to a polyline (`closed`: with the segment from the last to the first point).
pub fn polyline_distance(p: [f64; 3], pts: &[[f64; 3]], closed: bool) -> f64 {
    let n = pts.len();
    let mut d = f64::INFINITY;
    for i in 0..n.saturating_sub(1) {
        d = d.min(segment_distance(p, pts[i], pts[i + 1]));
    }
    if closed && n > 2 {
        d = d.min(segment_distance(p, pts[n - 1], pts[0]));
    }
    d
}

/// Smooth minimum (polynomial, blend width `k`).
pub fn smooth_min(a: f64, b: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return a.min(b);
    }
    let h = (k - (a - b).abs()).max(0.0) / k;
    a.min(b) - h * h * k * 0.25
}

/// Smooth maximum.
pub fn smooth_max(a: f64, b: f64, k: f64) -> f64 {
    -smooth_min(-a, -b, k)
}

/// The ball.
pub fn sphere(center: [f64; 3], radius: f64) -> Field {
    Arc::new(move |p| norm(sub(p, center)) - radius)
}

/// The box with half extents `half` rounded by `round` (`round = 0`: sharp).
pub fn rounded_box(center: [f64; 3], half: [f64; 3], round: f64) -> Field {
    Arc::new(move |p| {
        let q = sub(p, center);
        let d = [
            q[0].abs() - (half[0] - round),
            q[1].abs() - (half[1] - round),
            q[2].abs() - (half[2] - round),
        ];
        let out = norm([d[0].max(0.0), d[1].max(0.0), d[2].max(0.0)]);
        out + d[0].max(d[1]).max(d[2]).min(0.0) - round
    })
}

/// The infinite cylinder around the line through `point` along the unit vector `axis`.
pub fn cylinder(point: [f64; 3], axis: [f64; 3], radius: f64) -> Field {
    let n = norm(axis);
    let axis = scale(axis, 1.0 / n);
    Arc::new(move |p| {
        let q = sub(p, point);
        norm(sub(q, scale(axis, dot(q, axis)))) - radius
    })
}

/// The half space `{ (p − point)·normal ≤ 0 }`.
pub fn half_space(point: [f64; 3], normal: [f64; 3]) -> Field {
    let n = norm(normal);
    let normal = scale(normal, 1.0 / n);
    Arc::new(move |p| dot(sub(p, point), normal))
}

/// A tube of the given radius around a polyline.
pub fn tube(points: Vec<[f64; 3]>, closed: bool, radius: f64) -> Field {
    Arc::new(move |p| polyline_distance(p, &points, closed) - radius)
}

/// Union of tubes around several polylines `(points, closed, radius)`, blended smoothly with
/// width `blend` (0: sharp union).
pub fn tube_network(curves: Vec<(Vec<[f64; 3]>, bool, f64)>, blend: f64) -> Field {
    Arc::new(move |p| {
        let mut f = f64::INFINITY;
        for (pts, closed, r) in &curves {
            let d = polyline_distance(p, pts, *closed) - r;
            f = if f.is_finite() { smooth_min(f, d, blend) } else { d };
        }
        f
    })
}

/// Smooth union.
pub fn union(a: Field, b: Field, blend: f64) -> Field {
    Arc::new(move |p| smooth_min(a(p), b(p), blend))
}

/// Smooth difference `a \ b`.
pub fn difference(a: Field, b: Field, blend: f64) -> Field {
    Arc::new(move |p| smooth_max(a(p), -b(p), blend))
}

/// Smooth intersection.
pub fn intersection(a: Field, b: Field, blend: f64) -> Field {
    Arc::new(move |p| smooth_max(a(p), b(p), blend))
}
