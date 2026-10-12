//! The spike check: a handful of vertices that carry almost all of the distortion of the
//! surface without cones. The dataset mesh `camelshotnoise` is in `tests/dataset.rs`.

mod common;

use common::check;
use ddg_core::generators::{bumpy_grid_disk, icosphere, open_cylinder, torus};
use ddg_uv::primitives::{
    capsule, cone, cube, cylinder, grid_cube, holed_slab, torus as quad_torus, uv_sphere,
};
use ddg_uv::{ErrorKind, PolygonMesh, SpikeCheck, UnwrapOptions, Value, unwrap};

/// A 24 x 24 grid with two vertices pulled 2.0 out of a sheet of side 1, behind three loose
/// vertices (so that the numbers of the caller differ from those of the pipeline's meshes).
fn spiked_disk() -> (PolygonMesh, [u32; 2]) {
    let g = bumpy_grid_disk(24, 24, 0.02).unwrap();
    let at = |i: usize, j: usize| (i * 25 + j) as u32;
    let spikes = [at(12, 12), at(5, 18)];
    let mut positions = vec![[7.0, 7.0, 7.0], [8.0, 8.0, 8.0], [9.0, 9.0, 9.0]];
    positions.extend_from_slice(&g.positions);
    for s in spikes {
        positions[s as usize + 3][2] += 2.0;
    }
    let faces: Vec<[u32; 3]> = g.faces.iter().map(|f| f.map(|v| v + 3)).collect();
    let mesh = PolygonMesh::from_triangles(positions, &faces);
    (mesh, spikes.map(|s| s + 3))
}

#[test]
fn spikes_are_reported_in_the_callers_numbering() {
    let (mesh, spikes) = spiked_disk();
    let options = UnwrapOptions::default();
    let input = mesh.input();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let report = &map.report;
    assert_eq!(report.validation.n_loose_vertices, 3);
    let found = report.charts[0].spikes.as_ref().expect("the check ran");
    let mut expected = spikes;
    expected.sort_unstable();
    assert_eq!(found.vertices, expected, "{found:?}");
    assert!(found.n_carrying == 2 && found.share >= 0.99);
    assert!(found.cone_free_distortion > 1.0, "{found:?}");

    let warnings: Vec<&String> = report.warnings.iter().filter(|w| w.contains("thin spike")).collect();
    assert_eq!(warnings.len(), 1, "{:?}", report.warnings);
    let text = warnings[0];
    assert!(text.contains(&format!("vertices {}, {}", expected[0], expected[1])), "{text}");
    assert!(text.contains("Repair the mesh"), "{text}");

    let summary = report.summary();
    let get = |key: &str| summary.iter().find(|(k, _)| *k == key).map(|(_, v)| v.clone());
    assert_eq!(get("n_spike_vertices"), Some(Value::Int(2)));
    assert_eq!(
        get("spike_vertices"),
        Some(Value::Text(format!("{}, {}", expected[0], expected[1])))
    );
    assert!(matches!(get("cone_free_distortion"), Some(Value::Float(d)) if d > 1.0));
    assert!(matches!(get("seconds_spike_check"), Some(Value::Float(_))));
}

#[test]
fn the_check_can_be_switched_off_and_does_not_change_the_map() {
    let (mesh, _) = spiked_disk();
    let on = unwrap(&mesh.input(), &UnwrapOptions::default()).unwrap();
    let options = UnwrapOptions {
        spike_check: None,
        ..UnwrapOptions::default()
    };
    let off = unwrap(&mesh.input(), &options).unwrap();
    assert!(off.report.charts[0].spikes.is_none());
    assert!(off.report.warnings.iter().all(|w| !w.contains("thin spike")));
    assert!(off.report.summary().iter().all(|(k, _)| !k.contains("spike_vertices")));
    assert_eq!(on.uvs, off.uvs);
    assert_eq!(on.cones, off.cones);
}

#[test]
fn a_looser_limit_finds_more_and_a_stricter_one_nothing() {
    let (mesh, _) = spiked_disk();
    let strict = UnwrapOptions {
        spike_check: Some(SpikeCheck {
            max_vertices: 1,
            ..SpikeCheck::default()
        }),
        ..UnwrapOptions::default()
    };
    let map = unwrap(&mesh.input(), &strict).unwrap();
    let found = map.report.charts[0].spikes.as_ref().unwrap();
    assert!(found.vertices.is_empty() && found.n_carrying == 2, "{found:?}");
    assert!(map.report.warnings.iter().all(|w| !w.contains("thin spike")));
}

#[test]
fn bad_limits_are_refused() {
    let (mesh, _) = spiked_disk();
    for (name, check) in [
        ("spike_check.share", SpikeCheck { share: 0.0, ..SpikeCheck::default() }),
        ("spike_check.share", SpikeCheck { share: 1.5, ..SpikeCheck::default() }),
        ("spike_check.max_fraction", SpikeCheck { max_fraction: -0.1, ..SpikeCheck::default() }),
        ("spike_check.min_distortion", SpikeCheck { min_distortion: f64::NAN, ..SpikeCheck::default() }),
    ] {
        let options = UnwrapOptions {
            spike_check: Some(check),
            ..UnwrapOptions::default()
        };
        let error = unwrap(&mesh.input(), &options).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::Input);
        assert!(error.to_string().contains(name), "{error}");
    }
}

/// None of the generated meshes has spikes. Most have more than the 100 vertices below which the
/// check cannot fire at all.
#[test]
fn generated_meshes_have_no_spikes() {
    let mut meshes: Vec<(String, PolygonMesh)> = vec![
        ("uv_sphere 32x16".into(), uv_sphere(32, 16, 1.0).unwrap()),
        ("cube".into(), cube(2.0)),
        ("grid_cube 4".into(), grid_cube(4, 2.0)),
        ("grid_cube 8".into(), grid_cube(8, 2.0)),
        ("cylinder 32".into(), cylinder(32, 1, 1.0, 2.0).unwrap()),
        ("cone 32".into(), cone(32, 1.0, 2.0).unwrap()),
        ("capsule 32x8".into(), capsule(32, 8, 0.5, 1.0).unwrap()),
        ("torus 48x12".into(), quad_torus(48, 12, 1.0, 0.25).unwrap()),
        ("genus2".into(), holed_slab(2, 2, 1.0)),
    ];
    for (name, g) in [
        ("icosphere 2", icosphere(2).unwrap()),
        ("tri torus", torus(16, 8, 1.0, 0.35).unwrap()),
        ("open cylinder", open_cylinder(12, 4, 1.0, 2.0).unwrap()),
        ("bumpy disk", bumpy_grid_disk(24, 24, 0.4).unwrap()),
        ("flat disk", bumpy_grid_disk(24, 24, 0.0).unwrap()),
    ] {
        meshes.push((name.into(), PolygonMesh::from_triangles(g.positions.clone(), &g.faces)));
    }
    for (name, mesh) in &meshes {
        let options = UnwrapOptions::default();
        let input = mesh.input();
        let map = unwrap(&input, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(
            map.report.warnings.iter().all(|w| !w.contains("thin spike")),
            "{name}: {:?}",
            map.report.warnings
        );
        for chart in &map.report.charts {
            if let Some(found) = &chart.spikes {
                assert!(found.vertices.is_empty(), "{name}: {found:?}");
                println!(
                    "{name}: {} of {} vertices carry 99 % of the cone-free distortion {:.3e}",
                    found.n_carrying,
                    mesh.positions.len(),
                    found.cone_free_distortion
                );
            }
        }
    }
}
