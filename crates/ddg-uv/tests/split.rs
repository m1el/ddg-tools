//! Automatic splitting of charts that overlap themselves (`UnwrapOptions::split_overlaps`).

mod common;

use std::f64::consts::PI;

use common::{check, without_times};
use ddg_core::generators::{icosphere, open_cylinder};
use ddg_uv::{PolygonMesh, SplitOptions, SplitOutcome, UnwrapOptions, unwrap};

/// A flat annular sector `1 <= r <= 2`, `0 <= θ <= angle`, as an abstract surface: vertices
/// `(i, j)`, never identified, so that for `angle > 2π` the sector lies on top of itself in
/// space (the pipeline only looks at lengths). Developable with no cone to place, and its
/// flattening is the sector itself, which overlaps itself by `angle - 2π`.
fn spiral(nr: usize, nt: usize, angle: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let id = |i: usize, j: usize| (j * (nr + 1) + i) as u32;
    let mut positions = Vec::new();
    for j in 0..=nt {
        for i in 0..=nr {
            let r = 1.0 + i as f64 / nr as f64;
            let t = angle * j as f64 / nt as f64;
            positions.push([r * t.cos(), r * t.sin(), 0.0]);
        }
    }
    let mut faces = Vec::new();
    for j in 0..nt {
        for i in 0..nr {
            faces.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            faces.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    (positions, faces)
}

fn on() -> UnwrapOptions {
    UnwrapOptions {
        split_overlaps: Some(SplitOptions::default()),
        ..UnwrapOptions::default()
    }
}

#[test]
fn the_option_is_off_by_default() {
    assert!(UnwrapOptions::default().split_overlaps.is_none());
    let (p, f) = spiral(4, 24, 3.0 * PI);
    let mesh = PolygonMesh::from_triangles(p, &f);
    let map = unwrap(&mesh.input(), &UnwrapOptions::default()).unwrap();
    assert!(map.report.splitting.is_none());
    assert_eq!(map.report.charts.len(), 1);
    // Measured: the half of the footprint that is covered twice.
    let o = map.report.max_overlapped_fraction.unwrap();
    assert!((0.3..0.7).contains(&o), "{o}");
}

#[test]
fn a_spiral_that_overlaps_itself_is_split() {
    let (p, f) = spiral(4, 36, 3.0 * PI);
    let mesh = PolygonMesh::from_triangles(p, &f);
    let input = mesh.input();
    let options = on();
    let before = unwrap(&input, &UnwrapOptions::default()).unwrap();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let split = map.report.splitting.as_ref().unwrap();
    println!("{split:#?}");
    assert!(split.n_kept >= 1 && split.n_rejected == 0, "{split:?}");
    assert!(map.report.charts.len() >= 2);
    assert!(map.cones.is_empty());
    assert_eq!(map.report.flipped_triangles, 0);
    assert!(map.report.max_overlapped_fraction.unwrap() <= 0.001);
    assert!(split.overlap_before.unwrap() > 0.3);
    assert_eq!(
        split.overlap_before,
        before.report.overlapped_fraction,
        "the figure before the first split is the unsplit result"
    );
    let kept = split
        .splits
        .iter()
        .find(|s| s.outcome == SplitOutcome::Kept)
        .unwrap();
    assert!(kept.overlap_before > 0.3 && kept.overlap_after < 0.5 * kept.overlap_before);
    // The cut is short: it runs across the strip (width 1), not along it.
    assert!(split.cut_length < 4.0, "{}", split.cut_length);
    assert!(map.cut_edges.len() < 12);
    // The seams of the unsplit result are all still cut.
    assert!(before.cut_edges.iter().all(|e| map.cut_edges.contains(e)));
    // Deterministic.
    let again = unwrap(&input, &options).unwrap();
    assert_eq!(without_times(map), without_times(again));
}

#[test]
fn a_spiral_of_quads_is_cut_along_polygon_edges_only() {
    let (nr, nt) = (3, 24);
    let id = |i: usize, j: usize| (j * (nr + 1) + i) as u32;
    let (p, _) = spiral(nr, nt, 2.6 * PI);
    let mut polygons = Vec::new();
    for j in 0..nt {
        for i in 0..nr {
            polygons.push(vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    let mesh = PolygonMesh::from_polygons(p, &polygons);
    let input = mesh.input();
    let options = on();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let split = map.report.splitting.as_ref().unwrap();
    assert!(split.n_kept >= 1, "{split:?}");
    assert_eq!(map.report.flipped_triangles, 0);
    assert!(map.report.max_overlapped_fraction.unwrap() <= 0.001);
}

#[test]
fn a_map_without_overlap_is_unchanged() {
    let cylinder = open_cylinder(12, 4, 1.0, 2.0).unwrap();
    let sphere = icosphere(2).unwrap();
    for (name, g) in [("cylinder", cylinder), ("sphere", sphere)] {
        let mesh = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
        let input = mesh.input();
        let plain = unwrap(&input, &UnwrapOptions::default()).unwrap();
        let split = unwrap(&input, &on()).unwrap();
        let report = split.report.splitting.clone().unwrap();
        let o = plain.report.max_overlapped_fraction.unwrap();
        if o <= SplitOptions::default().min_overlap_fraction {
            assert_eq!(report.n_kept + report.n_rejected + report.n_no_cut, 0, "{name}");
            assert_eq!(report.runs, 1);
            assert_eq!(plain.uvs, split.uvs, "{name}");
            assert_eq!(plain.cut_edges, split.cut_edges);
            assert_eq!(plain.polygon_charts, split.polygon_charts);
        } else {
            // Never worse, never a mirrored triangle that was not there before.
            assert!(split.report.max_overlapped_fraction.unwrap() <= o, "{name}");
            assert!(split.report.flipped_triangles <= plain.report.flipped_triangles);
        }
    }
}

#[test]
fn the_budget_is_respected_and_reported() {
    let (p, f) = spiral(4, 36, 5.5 * PI);
    let mesh = PolygonMesh::from_triangles(p, &f);
    let options = UnwrapOptions {
        split_overlaps: Some(SplitOptions {
            max_splits: 1,
            ..SplitOptions::default()
        }),
        ..UnwrapOptions::default()
    };
    let map = unwrap(&mesh.input(), &options).unwrap();
    let split = map.report.splitting.as_ref().unwrap();
    assert_eq!(split.n_kept, 1);
    assert!(split.exhausted);
    assert!(map.report.warnings.iter().any(|w| w.contains("max_splits")));
}

#[test]
fn bad_options_are_refused() {
    let (p, f) = spiral(2, 8, 1.0);
    let mesh = PolygonMesh::from_triangles(p, &f);
    let input = mesh.input();
    let no_samples = UnwrapOptions {
        overlap_samples_per_triangle: 0.0,
        ..on()
    };
    assert!(unwrap(&input, &no_samples).is_err());
    for split in [
        SplitOptions { band: 1.0, ..SplitOptions::default() },
        SplitOptions { min_gain: -0.1, ..SplitOptions::default() },
        SplitOptions { crease_bias: f64::NAN, ..SplitOptions::default() },
        SplitOptions { min_overlap_fraction: 1.5, ..SplitOptions::default() },
    ] {
        let options = UnwrapOptions {
            split_overlaps: Some(split),
            ..UnwrapOptions::default()
        };
        assert!(unwrap(&input, &options).is_err(), "{split:?}");
    }
}
