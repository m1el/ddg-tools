//! End-to-end tests of `ddg_uv::unwrap` on small generated meshes.

mod common;

use common::{check, without_times};
use ddg_core::generators::{bumpy_grid_disk, icosphere, open_cylinder, torus};
use ddg_uv::{ErrorKind, MeshInput, PolygonMesh, UnwrapOptions, unwrap};

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
    let again = without_times(unwrap(&input, &options).unwrap());
    let first = without_times(map.clone());
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

use ddg_uv::primitives::{
    capsule, cone, cube, cylinder, grid_cube, holed_slab, torus as quad_torus, uv_sphere,
};
use ddg_uv::{ClosedHandling, Cone};

fn primitives() -> Vec<(&'static str, PolygonMesh)> {
    vec![
        ("uv_sphere", uv_sphere(16, 8, 1.0).unwrap()),
        ("cube", cube(2.0)),
        ("grid_cube", grid_cube(4, 2.0)),
        ("cylinder", cylinder(16, 2, 1.0, 2.0).unwrap()),
        ("cone", cone(16, 1.0, 2.0).unwrap()),
        ("capsule", capsule(16, 4, 0.5, 1.0).unwrap()),
        ("torus", quad_torus(16, 8, 1.0, 0.25).unwrap()),
        ("genus2", holed_slab(2, 2, 1.0)),
    ]
}

#[test]
fn polygonal_primitives_are_unwrapped_with_the_defaults() {
    for (name, mesh) in primitives() {
        let options = UnwrapOptions::default();
        let input = mesh.input();
        let map = unwrap(&input, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
        check(&mesh, &input, &map, &options);
        assert_eq!(map.report.charts.len(), 1, "{name}");
        println!(
            "{name}: {} cones (total {:.3}), {} cut edges, distortion {:.4}, flipped {}, \
             overlap {:?}, closed {:?}, warnings {}",
            map.cones.len(),
            map.report.cones.total_angle,
            map.cut_edges.len(),
            map.report.area_distortion.scale_free,
            map.report.flipped_triangles,
            map.report.overlapped_fraction,
            map.report.charts[0].closed,
            map.report.warnings.len()
        );
    }
}

/// The edges of a mesh that join two vertices at the given heights.
fn edges_at(mesh: &PolygonMesh, keep: impl Fn([f64; 3], [f64; 3]) -> bool) -> Vec<[u32; 2]> {
    mesh.edges()
        .into_iter()
        .filter(|e| keep(mesh.positions[e[0] as usize], mesh.positions[e[1] as usize]))
        .collect()
}

#[test]
fn a_cylinder_cut_along_its_rims_and_one_side_edge_is_developed_exactly() {
    let mesh = cylinder(16, 2, 1.0, 2.0).unwrap();
    // Both rims, and the side edges above the first vertex (x = 1, y = 0).
    let seams = edges_at(&mesh, |a, b| {
        let rim = a[2] == b[2] && a[2].abs() == 1.0;
        let side = a[2] != b[2] && a[0] == 1.0 && b[0] == 1.0 && a[1] == 0.0 && b[1] == 0.0;
        rim || side
    });
    assert_eq!(seams.len(), 16 + 16 + 2);
    let input = MeshInput {
        seam_edges: &seams,
        ..mesh.input()
    };
    let options = UnwrapOptions::default();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    assert_eq!(map.report.charts.len(), 3);
    assert!(map.cones.is_empty());
    // Exactly the seams are cut: the three pieces are disks without cones.
    assert_eq!(map.cut_edges, seams);
    assert_eq!(map.report.validation.seams.n_cut, 34);
    assert_eq!(map.report.validation.seams.n_pieces, 3);
    assert!(map.report.charts.iter().all(|c| c.closed == ClosedHandling::NotClosed));
    // The pieces are flat, so the map is an isometry up to the common scale of the packing.
    assert!(map.report.area_distortion.scale_free < 1e-20, "{:?}", map.report.area_distortion);
    assert!(map.report.charts.iter().all(|c| c.flatten.quasi_conformal_max < 1.0 + 1e-9));
    assert_eq!(map.report.flipped_triangles, 0);
    assert_eq!(map.report.overlapped_fraction, Some(0.0));
    let scale = map.report.packing.unwrap().scale;
    for p in 0..mesh.n_polygons() {
        let corners: Vec<usize> = mesh.polygon_corners(p).collect();
        for (i, &c) in corners.iter().enumerate() {
            let d = corners[(i + 1) % corners.len()];
            let (a, b) = (
                mesh.positions[mesh.corner_vertices[c] as usize],
                mesh.positions[mesh.corner_vertices[d] as usize],
            );
            let length = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2))
                .sqrt();
            let (q, r) = (map.uvs[c], map.uvs[d]);
            let flat = (q[0] - r[0]).hypot(q[1] - r[1]);
            assert!((flat - scale * length).abs() < 1e-12, "face {p}");
        }
    }
}

#[test]
fn seams_that_are_not_edges_are_ignored_and_one_edge_seams_extended() {
    let mesh = grid_cube(4, 2.0);
    let edges = mesh.edges();
    // An edge in the middle of a face: both ends are interior vertices of the closed cube.
    let one = edges[edges.len() / 2];
    // A diagonal of the first quad, and a pair of vertices far apart.
    let quad = mesh.polygon(0);
    let diagonal = [quad[0].min(quad[2]), quad[0].max(quad[2])];
    let far = [0, mesh.positions.len() as u32 - 1];
    let seams = [one, diagonal, far, [one[1], one[0]]];
    let input = MeshInput {
        seam_edges: &seams,
        ..mesh.input()
    };
    let options = UnwrapOptions::default();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let report = &map.report.validation.seams;
    assert_eq!((report.n_given, report.n_cut, report.n_not_edges), (3, 1, 2));
    assert_eq!(report.n_extension_edges, 1);
    assert!(map.cut_edges.binary_search(&one).is_ok());
    assert!(map.cut_edges.binary_search(&diagonal).is_err());
    assert_eq!(map.report.warnings.len(), 2, "{:?}", map.report.warnings);
    assert!(map.report.warnings[0].contains("2 seam edges are not edges of a face"));
    assert!(map.report.warnings[1].contains("extended by the shortest neighbouring edge"));
    // The seam opened the cube: it is not closed any more, so no slit is cut.
    assert_eq!(map.report.charts[0].closed, ClosedHandling::NotClosed);
}

#[test]
fn steering_is_respected_where_it_can_be() {
    let mesh = grid_cube(4, 2.0);
    // Cut the cube open along the four vertical edges and the top rim (a box without lid,
    // unfolded into a cross plus a lid): corners become boundary vertices.
    // (`grid_cube(4, 2.0)` spans [0, 2]³.)
    let side = |x: f64| x == 0.0 || x == 2.0;
    let seams = edges_at(&mesh, |a, b| {
        let vertical = side(a[0]) && side(a[1]) && a[0] == b[0] && a[1] == b[1];
        let rim = a[2] == 2.0 && b[2] == 2.0 && (side(a[0]) || side(a[1]))
            && (side(b[0]) || side(b[1]))
            && (a[0] == b[0] && side(a[0]) || a[1] == b[1] && side(a[1]));
        vertical || rim
    });
    assert_eq!(seams.len(), 4 * 4 + 4 * 4);
    // Interior vertices of the bottom face (z = 0).
    let bottom: Vec<u32> = (0..mesh.positions.len() as u32)
        .filter(|&v| {
            let p = mesh.positions[v as usize];
            p[2] == 0.0 && !side(p[0]) && !side(p[1])
        })
        .collect();
    assert_eq!(bottom.len(), 9);
    let on_seam = seams[0][0];
    let prescribed = [
        Cone {
            vertex: bottom[4],
            angle: 0.25,
        },
        Cone {
            vertex: on_seam,
            angle: 0.5,
        },
    ];
    let input = MeshInput {
        seam_edges: &seams,
        pinned_vertices: &bottom[0..2],
        forbidden_vertices: &bottom[5..9],
        prescribed_cones: &prescribed,
        ..mesh.input()
    };
    let options = UnwrapOptions::default();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    let cone = map.cones.iter().find(|c| c.vertex == bottom[4]).unwrap();
    assert_eq!(cone.angle, 0.25);
    assert!(map.cones.iter().all(|c| c.vertex != on_seam));
    let steering = &map.report.steering;
    assert_eq!(steering.n_prescribed_applied, 1);
    assert_eq!(steering.prescribed_not_applied, vec![on_seam]);
    assert_eq!(steering.n_pinned_applied, 2);
    assert_eq!(steering.n_forbidden_applied, 4);
    assert!(map.report.warnings.iter().any(|w| w.contains("1 prescribed cones could not")));
}

#[test]
fn loose_vertices_are_ignored() {
    let g = bumpy_grid_disk(3, 3, 0.1).unwrap();
    let mut positions = g.positions.clone();
    // Two vertices that no face uses, one of them steered.
    positions.push([5.0, 5.0, 5.0]);
    positions.push([6.0, 6.0, 6.0]);
    let loose = positions.len() as u32 - 1;
    let mesh = PolygonMesh::from_triangles(positions, &g.faces);
    let input = MeshInput {
        pinned_vertices: &[loose],
        seam_edges: &[[loose, 0]],
        ..mesh.input()
    };
    let options = UnwrapOptions::default();
    let map = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &map, &options);
    assert_eq!(map.report.validation.n_loose_vertices, 2);
    assert_eq!(map.report.steering.n_pinned_without_effect, 1);
    assert_eq!(map.report.validation.seams.n_not_edges, 1);
    // The same map as without them.
    let plain = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
    assert_eq!(map.uvs, unwrap(&plain.input(), &options).unwrap().uvs);
}

#[test]
fn the_merge_guard_falls_back_to_the_unmerged_cones() {
    use ddg_uv::{ClosedSurfaces, MergeClusters, MergeOutcome, SlitPlacement};
    // A cube with a slit across it at a small λ: the optimum has neighbouring cones.
    let mesh = cube(2.0);
    let base = UnwrapOptions {
        lambda: 0.02,
        closed_surfaces: ClosedSurfaces::Slit {
            placement: SlitPlacement::DiameterMiddle,
            relative_length: 1.0,
        },
        ..UnwrapOptions::default()
    };
    let run = |merge| {
        let options = UnwrapOptions { merge, ..base };
        let map = unwrap(&mesh.input(), &options).unwrap();
        check(&mesh, &mesh.input(), &map, &options);
        map
    };
    let always = run(MergeClusters::Always {
        max_cluster_size: usize::MAX,
    });
    let solve = always.report.charts[0].solve.as_ref().unwrap();
    assert_eq!(solve.merge, MergeOutcome::Merged);
    assert!(always.cones.len() < solve.n_cones_optimal);
    assert!(solve.distortion_merged > solve.distortion_optimal);
    assert_eq!(solve.distortion, solve.distortion_merged);

    // A guard that accepts no increase at all rejects that merge.
    let strict = run(MergeClusters::Guarded {
        max_cluster_size: usize::MAX,
        max_factor: 1.0,
        max_increase: 0.0,
    });
    let solve = strict.report.charts[0].solve.as_ref().unwrap();
    assert_eq!(solve.merge, MergeOutcome::Rejected);
    assert_eq!(strict.cones.len(), solve.n_cones_optimal);
    assert_eq!(solve.distortion, solve.distortion_optimal);
    assert!(
        strict.report.warnings.iter().any(|w| w.contains("merging the cone clusters")),
        "{:?}",
        strict.report.warnings
    );
    // Without merging the same cones come out, without the warning.
    let never = run(MergeClusters::Never);
    assert_eq!(never.cones, strict.cones);
    assert_eq!(
        never.report.charts[0].solve.as_ref().unwrap().merge,
        MergeOutcome::NotRequested
    );
    assert!(never.report.warnings.is_empty());
}

#[test]
fn closed_round_shapes_fall_back_to_free_sites_and_reoptimization_is_reported() {
    let mesh = uv_sphere(16, 8, 1.0).unwrap();
    let options = UnwrapOptions::default();
    let map = unwrap(&mesh.input(), &options).unwrap();
    // The raw optimum is a cone at every one of the 114 vertices.
    assert_eq!(
        map.report.charts[0].closed,
        ClosedHandling::FreeSites {
            n_automatic: 8,
            n_pinned: 0,
            fallback_from_cones: Some(114)
        }
    );
    assert_eq!(map.cones.len(), 8);
    let sum: f64 = map.cones.iter().map(|c| c.angle).sum();
    assert!((sum - 2.0 * std::f64::consts::TAU).abs() < 1e-9, "{sum}");
    assert_eq!(map.report.flipped_triangles, 0);
    assert_eq!(map.report.overlapped_fraction, Some(0.0));

    // Re-optimized angles on a disk: reported, and not worse in the solver's measure.
    let g = bumpy_grid_disk(16, 16, 0.4).unwrap();
    let disk = PolygonMesh::from_triangles(g.positions.clone(), &g.faces);
    let reoptimized = UnwrapOptions {
        lambda: 0.02,
        reoptimize_angles: true,
        ..UnwrapOptions::default()
    };
    let map = unwrap(&disk.input(), &reoptimized).unwrap();
    check(&disk, &disk.input(), &map, &reoptimized);
    let solve = map.report.charts[0].solve.as_ref().unwrap();
    let (before, n_before) = solve.before_reoptimization.unwrap();
    assert!(solve.distortion <= before, "{} > {before}", solve.distortion);
    assert!(map.cones.len() <= n_before);
}
