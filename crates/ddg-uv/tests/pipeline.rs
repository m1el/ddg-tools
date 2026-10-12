//! End-to-end tests of `ddg_uv::unwrap` on small generated meshes.

use ddg_core::generators::{bumpy_grid_disk, icosphere, open_cylinder, torus};
use ddg_uv::{ErrorKind, MeshInput, PolygonMesh, UnwrapOptions, UvMap, unwrap};

/// The checks every result must pass, whatever the input.
fn check(mesh: &PolygonMesh, input: &MeshInput<'_>, map: &UvMap, options: &UnwrapOptions) {
    assert_eq!(map.uvs.len(), mesh.corner_vertices.len());
    assert!(map.uvs.iter().all(|q| q[0].is_finite() && q[1].is_finite()));
    assert_eq!(map.polygon_charts.len(), mesh.n_polygons());
    let n_charts = map.report.charts.len() as u32;
    assert!(map.polygon_charts.iter().all(|&c| c < n_charts));
    if options.packing.is_some() {
        let eps = 1e-9;
        assert!(
            map.uvs
                .iter()
                .all(|q| q.iter().all(|&x| (-eps..=1.0 + eps).contains(&x))),
            "packed UVs must lie in the unit square"
        );
    }
    // Cuts run along polygon edges only, sorted and unique; every seam that is an interior
    // polygon edge is cut.
    let edges = mesh.edges();
    assert!(map.cut_edges.windows(2).all(|w| w[0] < w[1]));
    for e in &map.cut_edges {
        assert!(e[0] < e[1] && edges.binary_search(e).is_ok(), "cut {e:?}");
    }
    // Across an edge that is not cut, the two polygons agree on the UVs of its end points.
    let mut first_uv: std::collections::BTreeMap<[u32; 2], ([f64; 2], [f64; 2])> =
        std::collections::BTreeMap::new();
    for p in 0..mesh.n_polygons() {
        let corners: Vec<usize> = mesh.polygon_corners(p).collect();
        for (i, &c) in corners.iter().enumerate() {
            let d = corners[(i + 1) % corners.len()];
            let (a, b) = (mesh.corner_vertices[c], mesh.corner_vertices[d]);
            let key = [a.min(b), a.max(b)];
            if map.cut_edges.binary_search(&key).is_ok() {
                continue;
            }
            let here = if a < b {
                (map.uvs[c], map.uvs[d])
            } else {
                (map.uvs[d], map.uvs[c])
            };
            if let Some(there) = first_uv.insert(key, here) {
                assert_eq!(there, here, "edge {key:?} is not cut but its UVs differ");
            }
        }
    }
    assert!(map.cones.windows(2).all(|w| w[0].vertex < w[1].vertex));
    assert!(
        map.cones
            .iter()
            .all(|c| c.angle.is_finite() && c.angle < std::f64::consts::TAU)
    );
    for &v in input.forbidden_vertices {
        assert!(map.cones.iter().all(|c| c.vertex != v), "forbidden {v}");
    }
}

#[test]
fn a_bumpy_disk_is_unwrapped() {
    let g = bumpy_grid_disk(16, 16, 0.4).unwrap();
    let mesh = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
    let options = UnwrapOptions {
        lambda: 0.02,
        ..UnwrapOptions::default()
    };
    let input = mesh.input();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    assert_eq!(map.report.charts.len(), 1);
    assert!(map.report.warnings.is_empty(), "{:?}", map.report.warnings);
    assert_eq!(map.report.flipped_triangles, 0);
    assert!(!map.cones.is_empty());
    assert!(map.cones.iter().all(|c| c.angle > 0.0));
    let chart = &map.report.charts[0];
    println!(
        "bumpy disk: {} cones, Q mean {:.4} max {:.4}, distortion {:.5}, overlap {:?}",
        map.cones.len(),
        chart.flatten.quasi_conformal_mean,
        chart.flatten.quasi_conformal_max,
        map.report.area_distortion.scale_free,
        chart.flatten.overlapped_fraction
    );
    // The quasi-conformal error is discretization error of this coarse, steep grid. Measured
    // mean (max): 1.269 (1.979) with 16 x 16 cells, 1.091 (1.257) with 48 x 48, 1.045 (1.123)
    // with 96 x 96 (`examples/scratch.rs`).
    assert!(chart.flatten.quasi_conformal_mean < 1.35);
    // The README measure of the map is about twice the solver's `½ Σ m u²` (measured 0.1336
    // against 2 x 0.0841 here; 0.1720 against 2 x 0.0871 with 96 x 96 cells).
    let solver = chart.solve.as_ref().unwrap().distortion;
    assert!(map.report.area_distortion.scale_free < 2.0 * solver);
    assert!(map.report.area_distortion.scale_free > solver);
    // Identical on a second run, apart from the clock.
    let mut again = unwrap(&input, &options).unwrap();
    let mut first = map.clone();
    first.report.timings = Default::default();
    again.report.timings = Default::default();
    for report in [&mut first.report, &mut again.report] {
        for chart in &mut report.charts {
            if let Some(solve) = &mut chart.solve {
                solve.seconds = 0.0;
            }
        }
    }
    assert_eq!(first, again);
}

#[test]
fn closed_and_open_generated_meshes_are_unwrapped() {
    let meshes = [
        ("icosphere", icosphere(2).unwrap()),
        ("torus", torus(16, 8, 1.0, 0.35).unwrap()),
        ("cylinder", open_cylinder(12, 4, 1.0, 2.0).unwrap()),
    ];
    for (name, g) in meshes {
        let mesh = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
        let options = UnwrapOptions::default();
        let input = mesh.input();
        let map = unwrap(&input, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
        check(&mesh, &input, &map, &options);
        assert_eq!(map.report.charts.len(), 1, "{name}");
        println!(
            "{name}: {} cones, {} cut edges, distortion {:.4}, warnings {:?}",
            map.cones.len(),
            map.cut_edges.len(),
            map.report.area_distortion.scale_free,
            map.report.warnings
        );
    }
}

#[test]
fn malformed_arrays_are_input_errors() {
    let g = bumpy_grid_disk(2, 2, 0.0).unwrap();
    let mesh = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
    let mut bad = mesh.clone();
    bad.corner_vertices[4] = 99;
    let err = unwrap(&bad.input(), &UnwrapOptions::default()).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::Input);
    assert_eq!(
        err.to_string(),
        "corner_vertices[4] = 99 is out of range: there are 9 vertices (valid indices 0..=8)"
    );
}
