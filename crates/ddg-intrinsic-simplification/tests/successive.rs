//! The successive self-parameterization of the QEM baseline (Liu et al. 2021 Sec. 4; SPEC
//! §10.4, S23): valid barycentric locations, closeness to the surface, a closed-mesh requirement.

use ddg_core::generators::{grid_disk, icosphere, torus};
use ddg_intrinsic_simplification::baseline::{QemError, QemOptions, qem_simplify};

fn check(positions: &[[f64; 3]], faces: &[[u32; 3]], target: usize) {
    let o = QemOptions {
        successive_map: true,
        ..QemOptions::new(target)
    };
    let r = qem_simplify(positions, faces, &o).unwrap();
    let map = r.successive.as_ref().unwrap();
    let st = r.successive_stats.unwrap();
    assert_eq!(map.len(), positions.len());
    assert_eq!(st.n_lscm + st.n_tutte, r.report.n_collapses);
    for (v, p) in map.iter().enumerate() {
        assert!((p.face as usize) < r.faces.len(), "vertex {v}");
        assert!(p.bary.iter().all(|&b| b >= 0.0 && b.is_finite()));
        assert!((p.bary.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }
    // The image of a fine vertex stays near it (closest-point projection is the lower bound).
    let mean: f64 = map.iter().map(|p| p.distance).sum::<f64>() / map.len() as f64;
    let proj = r.report.mean_projection_distance;
    let size = r.report.fine_area.sqrt();
    assert!(
        mean < 0.2 * size,
        "mean distance {mean} (projection {proj}), size {size}, stats {st:?}"
    );
    eprintln!(
        "target {target}: {} collapses, stats {st:?}, mean dist {mean:.3e} vs projection {proj:.3e}, rejected_map {}",
        r.report.n_collapses, r.report.n_rejected_map
    );
}

#[test]
fn sphere_and_torus_maps_are_valid() {
    let s = icosphere(4).unwrap();
    check(&s.positions, &s.faces, 400);
    check(&s.positions, &s.faces, 60);
    let t = torus(24, 12, 2.0, 0.7).unwrap();
    check(&t.positions, &t.faces, 100);
}

#[test]
fn an_open_mesh_is_refused_with_a_message() {
    let g = grid_disk(6, 6).unwrap();
    let o = QemOptions {
        successive_map: true,
        ..QemOptions::new(10)
    };
    match qem_simplify(&g.positions, &g.faces, &o) {
        Err(QemError::InvalidOptions(m)) => assert!(m.contains("closed mesh"), "{m}"),
        other => panic!("{other:?}"),
    }
}
