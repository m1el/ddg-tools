//! Checks shared by the integration tests.
#![allow(dead_code)]

use std::collections::BTreeMap;

use ddg_uv::{MeshInput, PolygonMesh, UnwrapOptions, UvMap};

/// The checks every result must pass, whatever the input.
pub fn check(mesh: &PolygonMesh, input: &MeshInput<'_>, map: &UvMap, options: &UnwrapOptions) {
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
    let mut first_uv: BTreeMap<[u32; 2], ([f64; 2], [f64; 2])> = BTreeMap::new();
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


/// The map with every wall-clock time set to zero, for comparing two runs.
pub fn without_times(mut map: UvMap) -> UvMap {
    map.report.timings = Default::default();
    for chart in &mut map.report.charts {
        if let Some(solve) = &mut chart.solve {
            solve.seconds = 0.0;
        }
    }
    map
}

/// For every polygon edge `[min, max]`, the number of polygons that use it.
pub fn edge_uses(mesh: &PolygonMesh) -> BTreeMap<[u32; 2], usize> {
    let mut uses = BTreeMap::new();
    for p in 0..mesh.n_polygons() {
        let polygon = mesh.polygon(p);
        for (i, &a) in polygon.iter().enumerate() {
            let b = polygon[(i + 1) % polygon.len()];
            *uses.entry([a.min(b), a.max(b)]).or_insert(0) += 1;
        }
    }
    uses
}
