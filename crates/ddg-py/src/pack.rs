//! Packing charts into the unit square: simple, predictable, no rotation.
//!
//! A pipeline returns every chart in its own coordinate system, so charts overlap. Blender wants
//! one layout in `[0, 1]²`. [`pack_charts`] produces it:
//!
//! 1. Every chart gets the axis-aligned bounding box of its UVs. Charts are never rotated or
//!    mirrored, and all charts are scaled by the **same** factor `s`, so relative sizes (texel
//!    density) are exactly what the pipeline produced.
//! 2. Each box is padded by `margin` on all four sides, in units of the final unit square.
//! 3. The padded boxes are put on shelves: sorted by height, tallest first (ties: lower chart
//!    index first), placed left to right from the lower left corner, a new shelf above the
//!    previous one whenever the next box does not fit in the remaining width of the unit square.
//! 4. `s` is the largest scale for which the shelves fit in the unit square, found by 60 steps
//!    of bisection below the scale at which the largest box alone would fill the square.
//!
//! A single chart therefore ends up in the lower left corner with its longer side spanning
//! exactly `[margin, 1 − margin]`. With several charts, neighbouring charts are at least
//! `2 · margin` apart and every chart is at least `margin` away from the border.
//!
//! With many charts a fixed margin would leave no room for the charts themselves, so the margin
//! is capped at `1 / (4 · ⌈√n⌉)` for `n` charts: a quarter of the cell width of a `⌈√n⌉ × ⌈√n⌉`
//! grid. The margin that was used is in [`PackReport::margin`].
//!
//! The layout wastes space compared with a real packer. That is accepted for now; the pipeline
//! crate may take packing over later.

use crate::error::InputError;

/// What [`pack_charts`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackReport {
    /// Number of charts that have at least one polygon.
    pub n_charts: usize,
    /// The uniform scale applied to all UVs.
    pub scale: f64,
    /// The margin that was used: the requested one unless it had to be reduced.
    pub margin: f64,
    /// Fraction of the unit square covered by the (unpadded) chart bounding boxes.
    pub coverage: f64,
}

#[derive(Debug, Clone, Copy)]
struct ChartBox {
    chart: u32,
    min: [f64; 2],
    size: [f64; 2],
}

/// Shelf layout of the padded boxes at scale `s`. Returns the lower left corner of every padded
/// box (in the order of `boxes`), or `None` if they do not fit in a square of side
/// `1 + tolerance`.
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
            // Next shelf. Boxes are sorted by height, so the first box of a shelf is its tallest.
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

/// Packs the charts of a UV map into the unit square, in place. See the
/// [module documentation](self) for the layout rule.
///
/// * `uvs`: one UV per polygon corner.
/// * `polygon_offsets`: `n_polygons + 1` corner offsets, as in
///   [`crate::contract::UnwrapInput::polygon_offsets`].
/// * `polygon_charts`: the chart of every polygon. Chart ids need not be contiguous.
/// * `margin`: free space around every chart's bounding box, as a fraction of the unit square;
///   `0 ≤ margin < 0.5`.
///
/// The result is deterministic. UVs that are already NaN or infinite are rejected.
///
/// # Errors
/// [`InputError::BadOption`] for a margin outside `[0, 0.5)`, [`InputError::LengthMismatch`] when
/// the arrays do not describe the same mesh, [`InputError::NonFinite`] for a non-finite UV.
pub fn pack_charts(
    uvs: &mut [[f64; 2]],
    polygon_offsets: &[u32],
    polygon_charts: &[u32],
    margin: f64,
) -> Result<PackReport, InputError> {
    if !(margin.is_finite() && (0.0..0.5).contains(&margin)) {
        return Err(InputError::BadOption {
            name: "margin",
            value: format!("{margin}"),
            requirement: "the margin must be at least 0 and below 0.5",
        });
    }
    let n_polygons = polygon_charts.len();
    if polygon_offsets.len() != n_polygons + 1 {
        return Err(InputError::LengthMismatch {
            name: "polygon_charts",
            got: n_polygons,
            expected: polygon_offsets.len().saturating_sub(1),
            reason: "one chart id per polygon",
        });
    }
    let n_corners = polygon_offsets.last().copied().unwrap_or(0) as usize;
    if uvs.len() != n_corners {
        return Err(InputError::LengthMismatch {
            name: "uvs",
            got: uvs.len(),
            expected: n_corners,
            reason: "one UV per polygon corner",
        });
    }
    if let Some(c) = uvs
        .iter()
        .position(|uv| !(uv[0].is_finite() && uv[1].is_finite()))
    {
        let k = usize::from(uvs[c][0].is_finite());
        return Err(InputError::NonFinite {
            name: "uvs",
            position: 2 * c + k,
            value: uvs[c][k],
        });
    }

    // Bounding box per chart. Chart ids are compacted through a sorted list so that sparse ids
    // cost nothing.
    let mut ids: Vec<u32> = polygon_charts.to_vec();
    ids.sort_unstable();
    ids.dedup();
    let slot = |chart: u32| ids.binary_search(&chart).expect("id was collected above");
    let mut lo = vec![[f64::INFINITY; 2]; ids.len()];
    let mut hi = vec![[f64::NEG_INFINITY; 2]; ids.len()];
    for (p, &chart) in polygon_charts.iter().enumerate() {
        let k = slot(chart);
        for uv in &uvs[polygon_offsets[p] as usize..polygon_offsets[p + 1] as usize] {
            for axis in 0..2 {
                lo[k][axis] = lo[k][axis].min(uv[axis]);
                hi[k][axis] = hi[k][axis].max(uv[axis]);
            }
        }
    }
    // A chart whose polygons have no corners cannot occur with validated offsets; skip it anyway.
    let mut boxes: Vec<ChartBox> = (0..ids.len())
        .filter(|&k| lo[k][0] <= hi[k][0])
        .map(|k| ChartBox {
            chart: ids[k],
            min: lo[k],
            size: [hi[k][0] - lo[k][0], hi[k][1] - lo[k][1]],
        })
        .collect();
    if boxes.is_empty() {
        return Ok(PackReport {
            n_charts: 0,
            scale: 1.0,
            margin,
            coverage: 0.0,
        });
    }
    // Tallest first; total order, so the layout does not depend on the sort algorithm.
    boxes.sort_by(|a, b| b.size[1].total_cmp(&a.size[1]).then(a.chart.cmp(&b.chart)));

    // Cap the margin so that the padding leaves room for the charts (see the module docs).
    let grid = (boxes.len() as f64).sqrt().ceil();
    let mut margin = margin.min(0.25 / grid);
    // The padding alone fits by construction of the cap; if rounding says otherwise, pack
    // without a margin, which always fits at scale 0.
    let padding_only = match shelve(&boxes, 0.0, margin, 0.0) {
        Some(corners) => corners,
        None => {
            margin = 0.0;
            vec![[0.0, 0.0]; boxes.len()]
        }
    };

    let extent = boxes
        .iter()
        .map(|b| b.size[0].max(b.size[1]))
        .fold(0.0f64, f64::max);
    let (scale, corners) = if extent > 0.0 {
        // At `upper`, the largest box alone spans the square. The tolerance absorbs the rounding
        // of `upper * extent + 2 * margin`, which is 1 in exact arithmetic.
        let upper = (1.0 - 2.0 * margin) / extent;
        match shelve(&boxes, upper, margin, 1e-12) {
            Some(corners) => (upper, corners),
            None => {
                let (mut feasible, mut infeasible) = (0.0f64, upper);
                let mut best = padding_only;
                for _ in 0..60 {
                    let middle = 0.5 * (feasible + infeasible);
                    match shelve(&boxes, middle, margin, 0.0) {
                        Some(corners) => {
                            feasible = middle;
                            best = corners;
                        }
                        None => infeasible = middle,
                    }
                }
                (feasible, best)
            }
        }
    } else {
        // Every chart is a single point.
        (1.0, padding_only)
    };

    // Per chart: uv -> (uv - min) * scale + corner + margin.
    let mut offset = vec![[0.0f64; 2]; ids.len()];
    let mut origin = vec![[0.0f64; 2]; ids.len()];
    for (b, corner) in boxes.iter().zip(&corners) {
        let k = slot(b.chart);
        origin[k] = b.min;
        offset[k] = [corner[0] + margin, corner[1] + margin];
    }
    for (p, &chart) in polygon_charts.iter().enumerate() {
        let k = slot(chart);
        for uv in &mut uvs[polygon_offsets[p] as usize..polygon_offsets[p + 1] as usize] {
            for axis in 0..2 {
                uv[axis] = (uv[axis] - origin[k][axis]) * scale + offset[k][axis];
            }
        }
    }
    let coverage = boxes.iter().map(|b| b.size[0] * b.size[1]).sum::<f64>() * scale * scale;
    Ok(PackReport {
        n_charts: boxes.len(),
        scale,
        margin,
        coverage,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One quad per chart, with the given lower left corner and size.
    fn quads(rects: &[([f64; 2], [f64; 2])]) -> (Vec<[f64; 2]>, Vec<u32>, Vec<u32>) {
        let mut uvs = Vec::new();
        let mut offsets = vec![0u32];
        let mut charts = Vec::new();
        for (i, &([x, y], [w, h])) in rects.iter().enumerate() {
            uvs.extend([[x, y], [x + w, y], [x + w, y + h], [x, y + h]]);
            offsets.push(uvs.len() as u32);
            charts.push(i as u32);
        }
        (uvs, offsets, charts)
    }

    fn bbox(uvs: &[[f64; 2]]) -> ([f64; 2], [f64; 2]) {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for uv in uvs {
            for k in 0..2 {
                lo[k] = lo[k].min(uv[k]);
                hi[k] = hi[k].max(uv[k]);
            }
        }
        (lo, hi)
    }

    #[test]
    fn single_chart_fills_the_square_along_its_longer_side() {
        let (mut uvs, offsets, charts) = quads(&[([-3.0, 10.0], [4.0, 2.0])]);
        let report = pack_charts(&mut uvs, &offsets, &charts, 0.0).unwrap();
        assert_eq!(report.n_charts, 1);
        assert_eq!(report.scale, 0.25);
        assert_eq!(uvs, vec![[0.0, 0.0], [1.0, 0.0], [1.0, 0.5], [0.0, 0.5]]);
        assert_eq!(report.coverage, 0.5);

        let (mut uvs, offsets, charts) = quads(&[([-3.0, 10.0], [2.0, 4.0])]);
        let report = pack_charts(&mut uvs, &offsets, &charts, 0.125).unwrap();
        assert_eq!(report.scale, 0.1875);
        assert_eq!(
            uvs,
            vec![[0.125, 0.125], [0.5, 0.125], [0.5, 0.875], [0.125, 0.875]]
        );
        assert_eq!(report.margin, 0.125);
    }

    #[test]
    fn charts_do_not_overlap_and_keep_their_relative_size() {
        let rects = [
            ([0.0, 0.0], [1.0, 1.0]),
            ([0.0, 0.0], [2.0, 0.5]),
            ([5.0, 5.0], [0.5, 3.0]),
            ([-1.0, -1.0], [1.5, 1.5]),
            ([0.0, 0.0], [0.25, 0.25]),
        ];
        let margin = 0.01;
        let (mut uvs, offsets, charts) = quads(&rects);
        let report = pack_charts(&mut uvs, &offsets, &charts, margin).unwrap();
        assert_eq!(report.n_charts, 5);
        assert!(report.scale > 0.0);
        let boxes: Vec<_> = uvs.chunks(4).map(bbox).collect();
        for (i, (lo, hi)) in boxes.iter().enumerate() {
            // Same scale for every chart.
            for k in 0..2 {
                assert!(
                    ((hi[k] - lo[k]) - rects[i].1[k] * report.scale).abs() < 1e-12,
                    "chart {i} axis {k}"
                );
                // Inside the square, a margin away from the border.
                assert!(lo[k] >= margin - 1e-12 && hi[k] <= 1.0 - margin + 1e-9);
            }
            // At least two margins away from every other chart along some axis.
            for (j, (lo2, hi2)) in boxes.iter().enumerate().skip(i + 1) {
                let separated = (0..2).any(|k| {
                    lo[k] - hi2[k] >= 2.0 * margin - 1e-12 || lo2[k] - hi[k] >= 2.0 * margin - 1e-12
                });
                assert!(separated, "charts {i} and {j} are too close");
            }
        }
        // The layout uses the full width or the full height (the bisection converged).
        let (lo, hi) = bbox(&uvs);
        assert!((hi[0] - lo[0]).max(hi[1] - lo[1]) > 1.0 - 2.0 * margin - 1e-6);
        // Tallest chart first: chart 2 sits in the lower left corner.
        assert!((boxes[2].0[0] - margin).abs() < 1e-12 && (boxes[2].0[1] - margin).abs() < 1e-12);
    }

    #[test]
    fn packing_is_deterministic_and_handles_sparse_chart_ids() {
        let rects = [
            ([0.0, 0.0], [1.0, 1.0]),
            ([0.0, 0.0], [1.0, 1.0]),
            ([0.0, 0.0], [1.0, 1.0]),
        ];
        let (uvs, offsets, _) = quads(&rects);
        let charts = vec![7u32, 1_000_000, 7];
        let mut first = uvs.clone();
        let report = pack_charts(&mut first, &offsets, &charts, 0.0).unwrap();
        // Polygons 0 and 2 form one chart (identical boxes), polygon 1 another.
        assert_eq!(report.n_charts, 2);
        assert_eq!(first[0..4], first[8..12]);
        assert_eq!(report.scale, 0.5);
        // Equal heights: the lower chart id (7) comes first.
        assert_eq!(first[0], [0.0, 0.0]);
        assert_eq!(first[4], [0.5, 0.0]);
        let mut second = uvs.clone();
        pack_charts(&mut second, &offsets, &charts, 0.0).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn degenerate_charts_are_placed_without_dividing_by_zero() {
        // A chart that is a single point and one that is a horizontal segment.
        let mut uvs = vec![[3.0, 3.0]; 3];
        uvs.extend([[0.0, 1.0], [2.0, 1.0], [1.0, 1.0]]);
        let report = pack_charts(&mut uvs, &[0, 3, 6], &[0, 1], 0.05).unwrap();
        assert!(uvs.iter().flatten().all(|x| x.is_finite()));
        assert!(uvs.iter().flatten().all(|&x| (0.0..=1.0).contains(&x)));
        assert_eq!(report.n_charts, 2);
        // All points coincide: scale 1, nothing to fit.
        let mut uvs = vec![[3.0, 3.0]; 3];
        let report = pack_charts(&mut uvs, &[0, 3], &[0], 0.25).unwrap();
        assert_eq!(report.scale, 1.0);
        assert_eq!(uvs, vec![[0.25, 0.25]; 3]);
    }

    #[test]
    fn the_margin_is_capped_for_many_charts() {
        // 100 charts: the cap is 1 / (4 * 10).
        let rects: Vec<_> = (0..100).map(|_| ([0.0, 0.0], [1.0, 1.0])).collect();
        let (mut uvs, offsets, charts) = quads(&rects);
        let report = pack_charts(&mut uvs, &offsets, &charts, 0.2).unwrap();
        assert_eq!(report.margin, 0.025);
        // Padded boxes of side s + 0.05: 10 per shelf on 10 shelves fit for s up to 0.05.
        assert!((report.scale - 0.05).abs() < 1e-9, "{}", report.scale);
        assert!(uvs.iter().flatten().all(|&x| (0.0..=1.0).contains(&x)));
        // A margin below the cap is used as given.
        let (mut uvs, offsets, charts) = quads(&rects);
        let report = pack_charts(&mut uvs, &offsets, &charts, 0.01).unwrap();
        assert_eq!(report.margin, 0.01);
    }

    #[test]
    fn bad_arguments_are_rejected() {
        let (mut uvs, offsets, charts) = quads(&[([0.0, 0.0], [1.0, 1.0])]);
        assert_eq!(
            pack_charts(&mut uvs, &offsets, &charts, 0.5)
                .unwrap_err()
                .to_string(),
            "margin = 0.5: the margin must be at least 0 and below 0.5"
        );
        assert!(pack_charts(&mut uvs, &offsets, &charts, f64::NAN).is_err());
        assert_eq!(
            pack_charts(&mut uvs, &offsets, &[0, 0], 0.0)
                .unwrap_err()
                .to_string(),
            "polygon_charts: got 2 rows, expected 1 (one chart id per polygon)"
        );
        assert_eq!(
            pack_charts(&mut uvs[..3], &offsets, &charts, 0.0)
                .unwrap_err()
                .to_string(),
            "uvs: got 3 rows, expected 4 (one UV per polygon corner)"
        );
        uvs[2][1] = f64::NAN;
        assert_eq!(
            pack_charts(&mut uvs, &offsets, &charts, 0.0)
                .unwrap_err()
                .to_string(),
            "uvs[5] = NaN is not a finite number"
        );
    }
}
