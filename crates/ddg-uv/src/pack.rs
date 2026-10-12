//! Stage 4: packing the charts into the unit square.
//!
//! All charts are scaled by **one** factor, so their relative sizes stay what the flattening
//! made them (lengths in the plane are `e^u` times lengths on the surface in every chart
//! alike, and `u = 0` along seams and boundaries). A chart is moved, optionally rotated, never
//! mirrored and never scaled on its own.
//!
//! 1. Optionally every chart is rotated to the orientation in which its bounding box has the
//!    smallest area: the directions of the edges of its convex hull are tried (at most 720 of
//!    them, evenly spaced along the hull), ties go to the direction tried first.
//! 2. All boxes are laid on their long side, or all are stood on their short side; both are
//!    tried and the better one is kept (ties: on the long side).
//! 3. Every box is padded by `margin` on all four sides, in units of the unit square. With `n`
//!    charts the margin is capped at `1 / (4 ⌈√n⌉)`.
//! 4. The padded boxes go on shelves: tallest first (ties: lower chart index first), from left
//!    to right, a new shelf above whenever the next box does not fit the remaining width.
//! 5. The scale is the largest for which the shelves fit, found by bisection (60 steps) below
//!    the scale at which the largest box alone would fill the square.
//!
//! The layout is deterministic. It is a shelf layout of bounding boxes: it wastes the space
//! between a chart and its box and the ends of the shelves.

use ddg_bff::quality::signed_double_area;

use crate::error::UvError;
use crate::options::PackOptions;

/// What the packer did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackReport {
    /// Number of charts that have at least one polygon.
    pub n_charts: usize,
    /// The scale applied to every chart.
    pub scale: f64,
    /// The margin that was used: the requested one unless it had to be capped.
    pub margin: f64,
    /// Number of charts that were rotated (by more than a quarter turn counts too).
    pub n_rotated: usize,
    /// The fraction of the unit square covered by the bounding boxes of the charts.
    pub box_fill: f64,
    /// The fraction of the unit square covered by the charts themselves: the sum of the
    /// absolute areas of the polygons in the plane after scaling. Overlap within a chart is
    /// counted as if it were not there.
    pub fill: f64,
}

/// Largest number of hull edge directions tried per chart.
const MAX_DIRECTIONS: usize = 720;

/// Convex hull by Andrew's monotone chain, counter-clockwise, without collinear points.
fn convex_hull(points: &mut Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    points.sort_unstable_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    points.dedup();
    let n = points.len();
    if n < 3 {
        return points.clone();
    }
    let turn = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| signed_double_area(o, a, b);
    let mut hull: Vec<[f64; 2]> = Vec::with_capacity(2 * n);
    for &p in points.iter() {
        while hull.len() >= 2 && turn(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
            hull.pop();
        }
        hull.push(p);
    }
    let lower = hull.len() + 1;
    for &p in points.iter().rev().skip(1) {
        while hull.len() >= lower && turn(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
            hull.pop();
        }
        hull.push(p);
    }
    hull.pop();
    hull
}

/// The bounding box of `points` in the frame rotated by `(cos, sin)`: `(min, size)` of the
/// rotated coordinates `(c x + s y, −s x + c y)`.
fn rotated_box(points: &[[f64; 2]], cos: f64, sin: f64) -> ([f64; 2], [f64; 2]) {
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in points {
        let q = [cos * p[0] + sin * p[1], -sin * p[0] + cos * p[1]];
        for k in 0..2 {
            lo[k] = lo[k].min(q[k]);
            hi[k] = hi[k].max(q[k]);
        }
    }
    (lo, [hi[0] - lo[0], hi[1] - lo[1]])
}

/// The rotation `(cos, sin)` that minimizes the area of the bounding box of the hull, among
/// the directions of the hull's edges.
fn best_rotation(hull: &[[f64; 2]]) -> (f64, f64) {
    let (_, size) = rotated_box(hull, 1.0, 0.0);
    let mut best = (size[0] * size[1], 1.0, 0.0);
    if hull.len() < 3 {
        return (1.0, 0.0);
    }
    let step = hull.len().div_ceil(MAX_DIRECTIONS);
    for i in (0..hull.len()).step_by(step) {
        let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let length = dx.hypot(dy);
        if !(length > 0.0 && length.is_finite()) {
            continue;
        }
        let (cos, sin) = (dx / length, dy / length);
        let (_, size) = rotated_box(hull, cos, sin);
        let area = size[0] * size[1];
        // A strictly smaller box only: the unrotated chart wins ties, then earlier edges.
        if area < best.0 * (1.0 - 1e-12) {
            best = (area, cos, sin);
        }
    }
    (best.1, best.2)
}

#[derive(Debug, Clone, Copy)]
struct ChartBox {
    chart: usize,
    /// Rotation `(cos, sin)` applied to the chart.
    rotation: (f64, f64),
    /// Lower left corner and size of the rotated chart.
    min: [f64; 2],
    size: [f64; 2],
}

/// Shelf layout at scale `s`: the lower left corner of every padded box, in the order of
/// `boxes`, or `None` if they do not fit a square of side `1 + tolerance`.
fn shelve(boxes: &[ChartBox], s: f64, margin: f64, tolerance: f64) -> Option<Vec<[f64; 2]>> {
    let limit = 1.0 + tolerance;
    let mut corners = Vec::with_capacity(boxes.len());
    let (mut x, mut y, mut shelf_height) = (0.0f64, 0.0f64, 0.0f64);
    let mut shelf_is_empty = true;
    for b in boxes {
        let w = s * b.size[0] + 2.0 * margin;
        let h = s * b.size[1] + 2.0 * margin;
        if w > limit {
            return None;
        }
        if !shelf_is_empty && x + w > limit {
            y += shelf_height;
            x = 0.0;
            shelf_height = 0.0;
        }
        shelf_height = shelf_height.max(h);
        if y + h > limit {
            return None;
        }
        corners.push([x, y]);
        x += w;
        shelf_is_empty = false;
    }
    Some(corners)
}

/// The largest scale at which the boxes fit, and their corners. Boxes must be sorted.
fn fit(boxes: &[ChartBox], margin: f64) -> (f64, Vec<[f64; 2]>) {
    const TOLERANCE: f64 = 1e-12;
    let largest = boxes
        .iter()
        .map(|b| b.size[0].max(b.size[1]))
        .fold(0.0f64, f64::max);
    if !(largest > 0.0 && largest.is_finite()) {
        // Every chart is a point: any scale will do.
        let corners = shelve(boxes, 1.0, margin, TOLERANCE)
            .unwrap_or_else(|| vec![[0.0, 0.0]; boxes.len()]);
        return (1.0, corners);
    }
    let upper = (1.0 - 2.0 * margin) / largest;
    if let Some(corners) = shelve(boxes, upper, margin, TOLERANCE) {
        return (upper, corners);
    }
    let (mut lo, mut hi) = (0.0f64, upper);
    let mut best = shelve(boxes, 0.0, margin, TOLERANCE);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        match shelve(boxes, mid, margin, TOLERANCE) {
            Some(corners) => {
                lo = mid;
                best = Some(corners);
            }
            None => hi = mid,
        }
    }
    // With the capped margin the padding alone always fits, so `best` is set.
    (lo, best.unwrap_or_else(|| vec![[0.0, 0.0]; boxes.len()]))
}

/// Packs a UV map into the unit square, in place. `corner_chart[c]` is the chart of polygon
/// corner `c` (any `u32`; charts without corners do not exist for the packer), and
/// `polygon_offsets` delimits the polygons (for the fill figure of the report only).
pub(crate) fn pack_in_place(
    uvs: &mut [[f64; 2]],
    corner_chart: &[u32],
    polygon_offsets: &[u32],
    options: &PackOptions,
) -> PackReport {
    // Unsigned polygon areas, before anything moves (rigid motions keep them).
    let mut area = 0.0;
    for w in polygon_offsets.windows(2) {
        let (start, end) = (w[0] as usize, w[1] as usize);
        let mut twice = 0.0;
        for i in start + 1..end.saturating_sub(1) {
            twice += signed_double_area(uvs[start], uvs[i], uvs[i + 1]);
        }
        area += 0.5 * twice.abs();
    }

    let n_ids = corner_chart.iter().map(|&c| c as usize + 1).max().unwrap_or(0);
    let mut points: Vec<Vec<[f64; 2]>> = vec![Vec::new(); n_ids];
    for (c, &chart) in corner_chart.iter().enumerate() {
        points[chart as usize].push(uvs[c]);
    }
    let mut boxes: Vec<ChartBox> = Vec::new();
    let mut n_rotated = 0usize;
    for (chart, chart_points) in points.iter_mut().enumerate() {
        if chart_points.is_empty() {
            continue;
        }
        let hull = convex_hull(chart_points);
        let rotation = if options.rotate {
            best_rotation(&hull)
        } else {
            (1.0, 0.0)
        };
        let (min, size) = rotated_box(&hull, rotation.0, rotation.1);
        boxes.push(ChartBox {
            chart,
            rotation,
            min,
            size,
        });
    }
    let n_charts = boxes.len();
    if n_charts == 0 {
        return PackReport {
            n_charts: 0,
            scale: 1.0,
            margin: options.margin,
            n_rotated: 0,
            box_fill: 0.0,
            fill: 0.0,
        };
    }
    let cap = 1.0 / (4.0 * (n_charts as f64).sqrt().ceil());
    let margin = options.margin.min(cap);

    // Two candidates: every box on its long side, every box on its short side.
    let candidate = |stand: bool| {
        let mut oriented: Vec<ChartBox> = boxes
            .iter()
            .map(|b| {
                let wide = b.size[0] >= b.size[1];
                if !options.rotate || wide != stand {
                    *b
                } else {
                    // A further quarter turn of the rotated chart, (x, y) -> (−y, x): in the
                    // convention of `rotated_box` that is the rotation (sin, −cos).
                    let (c, s) = b.rotation;
                    ChartBox {
                        chart: b.chart,
                        rotation: (s, -c),
                        min: [-(b.min[1] + b.size[1]), b.min[0]],
                        size: [b.size[1], b.size[0]],
                    }
                }
            })
            .collect();
        oriented.sort_by(|a, b| b.size[1].total_cmp(&a.size[1]).then(a.chart.cmp(&b.chart)));
        let (scale, corners) = fit(&oriented, margin);
        (scale, oriented, corners)
    };
    let (mut scale, mut placed, mut corners) = candidate(false);
    if options.rotate {
        let (other_scale, other_placed, other_corners) = candidate(true);
        if other_scale > scale {
            (scale, placed, corners) = (other_scale, other_placed, other_corners);
        }
    }

    // Per chart id: rotation and the offset to add after rotating and scaling.
    let mut transform = vec![((1.0f64, 0.0f64), [0.0f64; 2]); n_ids];
    let mut box_area = 0.0;
    for (b, corner) in placed.iter().zip(&corners) {
        transform[b.chart] = (
            b.rotation,
            [
                corner[0] + margin - scale * b.min[0],
                corner[1] + margin - scale * b.min[1],
            ],
        );
        box_area += scale * b.size[0] * scale * b.size[1];
        if b.rotation != (1.0, 0.0) {
            n_rotated += 1;
        }
    }
    for (uv, &chart) in uvs.iter_mut().zip(corner_chart) {
        let ((cos, sin), offset) = transform[chart as usize];
        let rotated = [cos * uv[0] + sin * uv[1], -sin * uv[0] + cos * uv[1]];
        *uv = [
            scale * rotated[0] + offset[0],
            scale * rotated[1] + offset[1],
        ];
    }
    PackReport {
        n_charts,
        scale,
        margin,
        n_rotated,
        box_fill: box_area,
        fill: area * scale * scale,
    }
}

/// Packs the charts of a UV map into the unit square, in place: see the
/// [module documentation](self).
///
/// `uvs` holds one UV per polygon corner, `polygon_offsets` the `n_polygons + 1` corner
/// offsets of the polygons and `polygon_charts` the chart of every polygon.
///
/// # Errors
/// [`UvError::PolygonOffsets`] or [`UvError::NonFinite`] if the arrays do not fit each other
/// or a UV is not finite, [`UvError::InvalidOption`] for a margin outside `[0, ½)`.
pub fn pack_charts(
    uvs: &mut [[f64; 2]],
    polygon_offsets: &[u32],
    polygon_charts: &[u32],
    options: &PackOptions,
) -> Result<PackReport, UvError> {
    if !(options.margin.is_finite() && (0.0..0.5).contains(&options.margin)) {
        return Err(UvError::InvalidOption {
            name: "packing.margin",
            value: options.margin.to_string(),
            requirement: "must be at least 0 and below 0.5",
        });
    }
    let n_polygons = polygon_charts.len();
    if polygon_offsets.len() != n_polygons + 1 {
        return Err(UvError::PolygonOffsets {
            position: polygon_offsets.len().saturating_sub(1),
            value: polygon_offsets.last().copied().unwrap_or(0),
            requirement: format!(
                "there must be one more offset than faces ({n_polygons} faces, {} offsets)",
                polygon_offsets.len()
            ),
        });
    }
    let mut corner_chart = vec![0u32; uvs.len()];
    let mut previous = 0u32;
    for (p, w) in polygon_offsets.windows(2).enumerate() {
        if w[0] != previous || w[1] < w[0] || w[1] as usize > uvs.len() {
            return Err(UvError::PolygonOffsets {
                position: p + 1,
                value: w[1],
                requirement: format!(
                    "offsets must start at 0, not decrease and end at the number of UVs, {}",
                    uvs.len()
                ),
            });
        }
        corner_chart[w[0] as usize..w[1] as usize].fill(polygon_charts[p]);
        previous = w[1];
    }
    if previous as usize != uvs.len() {
        return Err(UvError::PolygonOffsets {
            position: n_polygons,
            value: previous,
            requirement: format!("the last offset must be the number of UVs, {}", uvs.len()),
        });
    }
    if let Some(c) = uvs.iter().position(|uv| !(uv[0].is_finite() && uv[1].is_finite())) {
        return Err(UvError::NonFinite {
            what: "uvs",
            position: c,
            value: if uvs[c][0].is_finite() {
                uvs[c][1]
            } else {
                uvs[c][0]
            },
        });
    }
    Ok(pack_in_place(uvs, &corner_chart, polygon_offsets, options))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(x: f64, y: f64, w: f64, h: f64) -> [[f64; 2]; 4] {
        [[x, y], [x + w, y], [x + w, y + h], [x, y + h]]
    }

    #[test]
    fn one_chart_spans_the_square_minus_the_margin() {
        let mut uvs = square(5.0, -3.0, 4.0, 2.0).to_vec();
        let options = PackOptions {
            margin: 0.1,
            rotate: false,
        };
        let report = pack_charts(&mut uvs, &[0, 4], &[0], &options).unwrap();
        assert!((report.scale - 0.2).abs() < 1e-12);
        assert!((uvs[0][0] - 0.1).abs() < 1e-12 && (uvs[0][1] - 0.1).abs() < 1e-12);
        assert!((uvs[2][0] - 0.9).abs() < 1e-12 && (uvs[2][1] - 0.5).abs() < 1e-12);
        assert!((report.fill - 0.32).abs() < 1e-12);
        assert!((report.box_fill - 0.32).abs() < 1e-12);
    }

    #[test]
    fn a_tilted_chart_is_rotated_into_its_smallest_box() {
        // A 4 x 1 rectangle turned by 30 degrees.
        let (c, s) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
        let mut uvs: Vec<[f64; 2]> = square(0.0, 0.0, 4.0, 1.0)
            .iter()
            .map(|p| [c * p[0] - s * p[1], s * p[0] + c * p[1]])
            .collect();
        let report = pack_charts(&mut uvs, &[0, 4], &[0], &PackOptions {
            margin: 0.0,
            rotate: true,
        })
        .unwrap();
        assert_eq!(report.n_rotated, 1);
        assert!((report.scale - 0.25).abs() < 1e-12, "{}", report.scale);
        assert!((report.box_fill - 0.25).abs() < 1e-12);
        // Lengths are scaled uniformly: it is still a rectangle of sides 1 and 1/4.
        let d = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).hypot(a[1] - b[1]);
        assert!((d(uvs[0], uvs[1]) - 1.0).abs() < 1e-12);
        assert!((d(uvs[1], uvs[2]) - 0.25).abs() < 1e-12);
        // Not mirrored.
        assert!(signed_double_area(uvs[0], uvs[1], uvs[2]) > 0.0);
    }

    #[test]
    fn charts_share_one_scale_and_do_not_touch() {
        let mut uvs = Vec::new();
        uvs.extend(square(0.0, 0.0, 2.0, 1.0));
        uvs.extend(square(0.0, 0.0, 1.0, 1.0));
        uvs.extend(square(7.0, 7.0, 0.5, 3.0));
        let original = uvs.clone();
        let options = PackOptions {
            margin: 0.01,
            rotate: true,
        };
        let report = pack_charts(&mut uvs, &[0, 4, 8, 12], &[0, 1, 2], &options).unwrap();
        assert_eq!(report.n_charts, 3);
        let d = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).hypot(a[1] - b[1]);
        for chart in 0..3 {
            for k in 0..4 {
                let (i, j) = (4 * chart + k, 4 * chart + (k + 1) % 4);
                let ratio = d(uvs[i], uvs[j]) / d(original[i], original[j]);
                assert!((ratio - report.scale).abs() < 1e-12);
            }
        }
        let bbox = |chart: usize| {
            let p = &uvs[4 * chart..4 * chart + 4];
            let lo = [
                p.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min),
                p.iter().map(|q| q[1]).fold(f64::INFINITY, f64::min),
            ];
            let hi = [
                p.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max),
                p.iter().map(|q| q[1]).fold(f64::NEG_INFINITY, f64::max),
            ];
            (lo, hi)
        };
        for a in 0..3 {
            let (lo, hi) = bbox(a);
            assert!(lo[0] >= 0.01 - 1e-12 && lo[1] >= 0.01 - 1e-12);
            assert!(hi[0] <= 0.99 + 1e-12 && hi[1] <= 0.99 + 1e-12);
            for b in a + 1..3 {
                let (lo2, hi2) = bbox(b);
                let apart = lo2[0] >= hi[0] + 0.02 - 1e-12
                    || lo[0] >= hi2[0] + 0.02 - 1e-12
                    || lo2[1] >= hi[1] + 0.02 - 1e-12
                    || lo[1] >= hi2[1] + 0.02 - 1e-12;
                assert!(apart, "charts {a} and {b} are closer than twice the margin");
            }
        }
    }

    #[test]
    fn bad_arrays_are_refused() {
        let mut uvs = square(0.0, 0.0, 1.0, 1.0).to_vec();
        let options = PackOptions::default();
        assert!(pack_charts(&mut uvs, &[0, 3], &[0], &options).is_err());
        assert!(pack_charts(&mut uvs, &[0, 4], &[0, 1], &options).is_err());
        assert!(pack_charts(&mut uvs, &[1, 4], &[0], &options).is_err());
        uvs[2][1] = f64::NAN;
        let err = pack_charts(&mut uvs, &[0, 4], &[0], &options).unwrap_err();
        assert_eq!(err.to_string(), "uvs[2] = NaN is not a finite number");
        let mut uvs = square(0.0, 0.0, 1.0, 1.0).to_vec();
        let bad = PackOptions {
            margin: 0.5,
            rotate: true,
        };
        assert!(pack_charts(&mut uvs, &[0, 4], &[0], &bad).is_err());
    }

    #[test]
    fn many_tiny_charts_get_a_capped_margin() {
        let mut uvs = Vec::new();
        let mut offsets = vec![0u32];
        let mut charts = Vec::new();
        for i in 0..100u32 {
            uvs.extend(square(0.0, 0.0, 1.0, 1.0));
            offsets.push(uvs.len() as u32);
            charts.push(i);
        }
        let report = pack_charts(&mut uvs, &offsets, &charts, &PackOptions {
            margin: 0.2,
            rotate: true,
        })
        .unwrap();
        assert!((report.margin - 1.0 / 40.0).abs() < 1e-15);
        assert!(report.scale > 0.0);
        assert!(uvs.iter().all(|p| p.iter().all(|x| (0.0..=1.0).contains(x))));
    }
}
