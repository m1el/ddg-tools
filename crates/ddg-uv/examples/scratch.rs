//! Scratch measurements (not part of the documented examples).
use ddg_core::generators::bumpy_grid_disk;
use ddg_uv::{PolygonMesh, UnwrapOptions, unwrap};

fn main() {
    for (n, amplitude, lambda) in [
        (16, 0.4, 0.02),
        (48, 0.4, 0.02),
        (96, 0.4, 0.02),
        (16, 0.1, 0.02),
        (48, 0.1, 0.02),
        (48, 0.1, 0.002),
    ] {
        let g = bumpy_grid_disk(n, n, amplitude).unwrap();
        let mesh = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
        let options = UnwrapOptions {
            lambda,
            ..UnwrapOptions::default()
        };
        let map = unwrap(&mesh.input(), &options).unwrap();
        let chart = &map.report.charts[0];
        println!(
            "n {n} amplitude {amplitude} lambda {lambda}: {} cones, Q mean {:.4} max {:.4}, \
             distortion {:.5}, flipped {}, overlap {:?}, solver distortion {:?}",
            map.cones.len(),
            chart.flatten.quasi_conformal_mean,
            chart.flatten.quasi_conformal_max,
            map.report.area_distortion.scale_free,
            map.report.flipped_triangles,
            chart.flatten.overlapped_fraction,
            chart.solve.as_ref().map(|s| s.distortion),
        );
    }
}
