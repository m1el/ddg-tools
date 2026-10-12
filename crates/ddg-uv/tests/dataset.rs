//! The pipeline on the optimal-cones dataset (`data/optimal-cones/optimal-cones-data/data`, or
//! `DDG_DATASET`). Without the data the test prints a notice and passes. A release build runs
//! all 23 meshes (about a minute, most of it the 306k-triangle brain). A debug build runs only
//! `bump` (6.5 s): unoptimized, the cone solver needs 32 s for `airplane`, 85 s for `hand` and
//! 137 s for `axes`.

mod common;

use std::path::PathBuf;

use common::check;
use ddg_core::io::read_obj;
use ddg_uv::{ClosedHandling, PolygonMesh, UnwrapOptions, area_distortion, unwrap};

fn dataset_dir() -> PathBuf {
    std::env::var_os("DDG_DATASET").map_or_else(
        || {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../data/optimal-cones/optimal-cones-data/data")
        },
        PathBuf::from,
    )
}

fn names() -> Option<Vec<String>> {
    let entries = std::fs::read_dir(dataset_dir()).ok()?;
    let mut names: Vec<String> = entries
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let name = e.file_name().into_string().ok()?;
            name.strip_suffix(".obj").map(str::to_owned)
        })
        .collect();
    names.sort();
    (!names.is_empty()).then_some(names)
}

#[test]
fn every_dataset_mesh_is_unwrapped() {
    let Some(mut names) = names() else {
        println!("dataset not found at {}: skipped", dataset_dir().display());
        return;
    };
    if cfg!(debug_assertions) {
        names.retain(|n| n == "bump");
        println!("debug build: only bump is run; use --release for all 23 meshes");
    }
    let options = UnwrapOptions::default();
    for name in &names {
        let obj = read_obj(dataset_dir().join(format!("{name}.obj"))).unwrap();
        let mesh = PolygonMesh::from_triangles(obj.positions.clone(), &obj.triangles().unwrap());
        let input = mesh.input();
        let map = unwrap(&input, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
        check(&mesh, &input, &map, &options);
        let report = &map.report;
        assert_eq!(report.charts.len(), 1, "{name}");
        assert!(report.cones.n_negative == 0, "{name}: negative cones with the defaults");
        let chart = &report.charts[0];
        if chart.closed != ClosedHandling::NotClosed {
            // A closed surface of genus 0: the cone angles sum to 4π.
            assert!(
                (report.cones.angle_sum - 2.0 * std::f64::consts::TAU).abs() < 1e-6,
                "{name}: {}",
                report.cones.angle_sum
            );
        }
        // The spiked camel is the one mesh on which the result is not usable; it must say so.
        let spiked = name == "camelshotnoise";
        assert_eq!(
            report.warnings.iter().any(|w| w.contains("upper bound")),
            spiked,
            "{name}: {:?}",
            report.warnings
        );
        if !spiked {
            // Measured with the defaults: at most 0.028 % of a footprint overlaps
            // (camel_noise2) and at most 0.48 % of the triangles are mirrored (camel_noise2:
            // 108 of 22 704; brain 38, camel_noise1 19, mannequin 2, all others none). The
            // bounds are the thresholds of the pipeline's own warnings.
            assert!(report.overlapped_fraction.unwrap() < 0.01, "{name}");
            assert!(100 * report.flipped_triangles < mesh.n_polygons(), "{name}");
            assert!(report.warnings.is_empty(), "{name}: {:?}", report.warnings);
        }
        // The published map, by the same measure, for the record.
        let published = obj.corner_uvs().map(|uvs| {
            area_distortion(&mesh.positions, &mesh.corner_vertices, &mesh.triangles, &uvs)
                .scale_free
        });
        println!(
            "{name}: {} cones, distortion {:.5} (published {:?}), overlap {:?}, flipped {}, {:.2} s",
            report.cones.n_cones,
            report.area_distortion.scale_free,
            published,
            report.overlapped_fraction,
            report.flipped_triangles,
            report.timings.total
        );
    }
}
