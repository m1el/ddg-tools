//! Analytic and structural tests of the cone flattening.

// Index loops mirror the subscripts of the formulas being checked.
#![allow(clippy::needless_range_loop)]

mod common;

use std::f64::consts::PI;

use common::{compare_layouts, flat_polyhedron, layout_angle_sums, max_length_deviation};
use ddg_bff::{
    BffError, Cone, ConeFlattenInput, ConeFlattenOptions, ConeFlattening, CutSource, DiskBff,
    DiskOptions, GaussBonnetPolicy, ScaleFactorOptions, component_topology, flatten_with_cones,
    mollify_edge_lengths, solve_scale_factors,
};
use ddg_core::generators;
use ddg_core::geometry::edge_length;
use ddg_core::operators::angle_defects;
use ddg_core::{Geometry, TriMesh};

const PLAN: CutSource<'static> = CutSource::Plan {
    may_cut: None,
    must_cut: &[],
};

fn flatten(
    mesh: &TriMesh,
    geometry: &Geometry,
    cones: &[Cone],
    cut: CutSource<'_>,
) -> Result<ConeFlattening, BffError> {
    flatten_with_cones(
        &ConeFlattenInput {
            mesh,
            geometry,
            cones,
            scale_factors: None,
            cut,
            anchors: &[],
        },
        &ConeFlattenOptions::default(),
    )
}

fn setup(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> (TriMesh, Geometry) {
    let mesh = TriMesh::new(positions.len(), faces.to_vec()).unwrap();
    let geometry = Geometry::from_positions(&mesh, positions).unwrap();
    (mesh, geometry)
}

/// Every chart of a flattening is a disk: connected, Euler characteristic 1, one boundary loop.
fn assert_disk_charts(flat: &ConeFlattening) {
    let topology = component_topology(&flat.cut.mesh);
    assert_eq!(topology.len(), flat.n_charts());
    for (chart, t) in topology.iter().enumerate() {
        assert_eq!(t.euler_characteristic, 1, "chart {chart}");
        assert_eq!(t.boundary_loops, 1, "chart {chart}");
        assert_eq!(t.genus, 0, "chart {chart}");
    }
    assert_eq!(flat.cut.mesh.boundary_loops().len(), flat.n_charts());
    for (f, &chart) in flat.face_chart.iter().enumerate() {
        assert_eq!(chart, flat.cut.mesh.face_components()[f]);
    }
}

/// The lengths of the two sides of every cut edge agree, measured on the corner UVs.
fn max_seam_mismatch(mesh: &TriMesh, flat: &ConeFlattening) -> f64 {
    let side = |h: usize| {
        let (a, b) = (flat.corner_uvs[h], flat.corner_uvs[mesh.next(h)]);
        (b[0] - a[0]).hypot(b[1] - a[1])
    };
    flat.cut_edges()
        .iter()
        .map(|&e| {
            let h = mesh.edge_halfedge(e as usize);
            let t = mesh.twin(h).unwrap();
            (side(h) - side(t)).abs() / side(h).max(side(t))
        })
        .fold(0.0, f64::max)
}

#[test]
fn flat_polyhedra_unfold_isometrically() {
    // A polyhedron with flat faces is a cone metric already: with its own angle defects as
    // cone angles the scale factor vanishes and the flattening is an isometric unfolding.
    for (n_faces, level) in [(4, 3), (8, 3), (20, 2)] {
        let (positions, faces) = flat_polyhedron(n_faces, level);
        let (mesh, geometry) = setup(&positions, &faces);
        let omega = angle_defects(&mesh, &geometry).unwrap();
        let cones: Vec<Cone> = (0..mesh.n_vertices())
            .filter(|&v| omega[v].abs() > 1e-9)
            .map(|v| Cone {
                vertex: v as u32,
                angle: omega[v],
            })
            .collect();
        assert_eq!(
            cones.len(),
            [4, 6, 12][[4, 8, 20].iter().position(|&n| n == n_faces).unwrap()]
        );
        let flat = flatten(&mesh, &geometry, &cones, PLAN).unwrap();
        assert_disk_charts(&flat);
        assert_eq!(flat.n_charts(), 1);
        let u_max = flat
            .scale_factors
            .iter()
            .fold(0.0f64, |m, u| m.max(u.abs()));
        assert!(u_max < 1e-10, "{n_faces} faces: |u| up to {u_max}");
        let deviation = max_length_deviation(&mesh, &geometry, &flat.corner_uvs);
        assert!(
            deviation < 1e-9,
            "{n_faces} faces: lengths off by {deviation}"
        );
        let report = &flat.charts[0];
        assert_eq!(report.distortion.flipped_triangles, 0);
        assert!(report.distortion.quasi_conformal.max < 1.0 + 1e-9);
        assert!(report.distortion.log_scale.max.abs() < 1e-9);
        assert!(report.closure.gap_before < 1e-10);
        assert!(report.closure.turning_error.abs() < 1e-10);
        assert!(report.is_injective(), "{report:?}");
        assert!(flat.seams.max_length_mismatch < 1e-12);
        for cone in &flat.cones {
            assert!((cone.realized_angle - cone.angle).abs() < 1e-9, "{cone:?}");
            assert!((cone.implied_angle - cone.angle).abs() < 1e-9, "{cone:?}");
        }
    }
}

#[test]
fn angle_sums_are_two_pi_minus_cone_angle() {
    // The twelve vertices of the icosahedron under a sphere, each with a cone of π/3.
    let g = generators::icosphere(3).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let cones: Vec<Cone> = (0..12)
        .map(|vertex| Cone {
            vertex,
            angle: PI / 3.0,
        })
        .collect();
    let flat = flatten(&mesh, &geometry, &cones, PLAN).unwrap();
    assert_disk_charts(&flat);
    assert_eq!(flat.charts[0].distortion.flipped_triangles, 0);
    assert!(flat.charts[0].is_injective());
    let sums = layout_angle_sums(&mesh, &flat.corner_uvs);
    for v in 0..mesh.n_vertices() {
        let mu = if v < 12 { PI / 3.0 } else { 0.0 };
        assert!(
            (sums[v] - (2.0 * PI - mu)).abs() < 1e-9,
            "vertex {v}: angle sum {} instead of {}",
            sums[v],
            2.0 * PI - mu
        );
    }
    // Seamless: equal lengths across every cut edge.
    assert!(max_seam_mismatch(&mesh, &flat) < 1e-12);
    assert!(flat.seams.max_length_mismatch < 1e-12);
    assert_eq!(flat.seams.n_cut_edges, flat.cut_edges().len());
    // Nearly conformal, with the scale factor it was built from.
    let q = flat.charts[0].distortion.quasi_conformal;
    assert!(q.mean < 1.05 && q.max < 1.5, "{q:?}");
    let solve = flat.scale_factor_solve.as_ref().unwrap();
    assert!(solve.max_residual < 1e-11 && flat.max_interior_curvature < 1e-11);
    assert!(solve.closed_components[0].discrepancy.abs() < 1e-12);
    assert!(flat.extension_residual < 1e-11);
    // The cut passes through all cones and nothing else was cut than the plan says.
    let plan = flat.plan.as_ref().unwrap();
    assert_eq!(plan.cut_edges, flat.cut_edges());
    for cone in &flat.cones {
        assert!(cone.n_wedges >= 1);
        assert!(cone.max_wedge_angle < 2.0 * PI);
    }
}

#[test]
fn a_negative_cone_at_the_end_of_a_cut_overlaps_and_seams_can_fix_it() {
    // A flat square with a cone of +π/2 and one of −π/2.
    let g = generators::grid_disk(12, 12).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let at = |i: usize, j: usize| (i * 13 + j) as u32;
    let cones = [
        Cone {
            vertex: at(4, 4),
            angle: PI / 2.0,
        },
        Cone {
            vertex: at(8, 8),
            angle: -PI / 2.0,
        },
    ];
    let flat = flatten(&mesh, &geometry, &cones, PLAN).unwrap();
    assert_disk_charts(&flat);
    let negative = &flat.cones[1];
    assert!(
        (negative.realized_angle - negative.angle).abs() < 1e-9 || {
            // With overlap the signed sum can be off by a full turn only if triangles flip.
            flat.charts[0].distortion.flipped_triangles > 0
        }
    );
    assert_eq!(
        negative.n_wedges, 1,
        "the planner joins each cone to the boundary"
    );
    // One wedge has to hold 2π + π/2: the layout overlaps itself around the cone.
    assert!(negative.max_wedge_angle > 2.0 * PI);
    assert!(flat.charts[0].overlap.wedges_over_2pi >= 1);
    assert!(!flat.charts[0].is_injective());
    let positive = &flat.cones[0];
    assert!((positive.realized_angle - PI / 2.0).abs() < 1e-9);
    assert!(positive.max_wedge_angle < 2.0 * PI);

    // The caller's remedy: a seam straight through the negative cone, from boundary to
    // boundary. The cone then has two wedges of 5π/4 each.
    let seam: Vec<u32> = (0..12)
        .map(|i| {
            mesh.find_edge(at(i, 8) as usize, at(i + 1, 8) as usize)
                .unwrap() as u32
        })
        .collect();
    let flat = flatten(
        &mesh,
        &geometry,
        &cones,
        CutSource::Plan {
            may_cut: None,
            must_cut: &seam,
        },
    )
    .unwrap();
    assert_disk_charts(&flat);
    assert_eq!(flat.n_charts(), 2);
    let negative = &flat.cones[1];
    assert_eq!(negative.n_wedges, 2);
    assert!((negative.realized_angle - negative.angle).abs() < 1e-9);
    assert!(negative.max_wedge_angle < 2.0 * PI);
    assert_eq!(
        flat.charts
            .iter()
            .map(|c| c.overlap.wedges_over_2pi)
            .sum::<usize>(),
        0
    );
    assert_eq!(flat.seams.n_between_charts, 12);
    assert!(flat.seams.max_length_mismatch < 1e-12);
}

#[test]
fn flat_annulus_and_cylinder_are_unfolded_isometrically() {
    // A flat annulus: u = 0, so cutting it open must reproduce it.
    let g = generators::annulus(24, 4, 1.0, 2.0).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let flat = flatten(&mesh, &geometry, &[], PLAN).unwrap();
    assert_disk_charts(&flat);
    assert_eq!(flat.cut_edges().len(), 4);
    let original: Vec<[f64; 2]> = mesh
        .faces()
        .iter()
        .flat_map(|f| f.map(|v| [g.positions[v as usize][0], g.positions[v as usize][1]]))
        .collect();
    let agreement = compare_layouts(&flat.face_chart, 1, &flat.corner_uvs, &original);
    assert!(agreement.max < 1e-10, "{agreement:?}");
    assert!((agreement.scale - 1.0).abs() < 1e-10);
    assert!(flat.seams.max_length_mismatch < 1e-12);
    assert_eq!(flat.charts[0].distortion.flipped_triangles, 0);

    // A cylinder is flat intrinsically: it unrolls to a rectangle without stretch.
    let g = generators::open_cylinder(16, 5, 1.0, 2.0).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let flat = flatten(&mesh, &geometry, &[], PLAN).unwrap();
    assert_disk_charts(&flat);
    assert!(max_length_deviation(&mesh, &geometry, &flat.corner_uvs) < 1e-10);
    assert!(flat.charts[0].distortion.quasi_conformal.max < 1.0 + 1e-10);
    assert!(flat.charts[0].is_injective());
}

/// A torus with the metric of a flat rectangle with opposite sides identified.
fn flat_torus(n_major: usize, n_minor: usize, a: f64, b: f64) -> (TriMesh, Geometry) {
    let g = generators::torus(n_major, n_minor, 2.0, 0.5).unwrap();
    let mesh = g.topology().unwrap();
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[p, q]| {
            let (ip, jp) = (p as usize / n_minor, p as usize % n_minor);
            let (iq, jq) = (q as usize / n_minor, q as usize % n_minor);
            let di = if ip == iq { 0.0 } else { a };
            let dj = if jp == jq { 0.0 } else { b };
            di.hypot(dj)
        })
        .collect();
    let geometry = Geometry::from_edge_lengths(&mesh, lengths).unwrap();
    (mesh, geometry)
}

#[test]
fn tori_without_cones() {
    // Flat torus: the unfolding is an isometry onto a parallelogram. The two sides of every
    // cut edge are parallel, so the closing constraints are redundant.
    let (mesh, geometry) = flat_torus(12, 8, 0.3, 0.2);
    let flat = flatten(&mesh, &geometry, &[], PLAN).unwrap();
    assert_disk_charts(&flat);
    assert_eq!(flat.plan.as_ref().unwrap().n_loop_edges, 2);
    assert!(max_length_deviation(&mesh, &geometry, &flat.corner_uvs) < 1e-10);
    let closure = &flat.charts[0].closure;
    assert_eq!(closure.dropped_directions, 2, "{closure:?}");
    assert!(closure.gap_after < 1e-12 && closure.max_length_change < 1e-12);
    assert!(flat.seams.max_length_mismatch < 1e-11);
    assert!(flat.charts[0].is_injective());

    // A torus of revolution: its flat metric is not the given one, and the discrete boundary
    // curve closes only up to discretization error, which lengths shared across parallel cut
    // sides cannot remove. The gap is spread along the curve and reported.
    let g = generators::torus(32, 16, 2.0, 0.7).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let flat = flatten(&mesh, &geometry, &[], PLAN).unwrap();
    assert_disk_charts(&flat);
    assert_eq!(flat.n_charts(), 1);
    let report = &flat.charts[0];
    assert_eq!(report.distortion.flipped_triangles, 0);
    assert!(report.distortion.quasi_conformal.mean < 1.1, "{report:?}");
    assert!(report.closure.max_length_change <= 0.5);
    assert!(report.closure.turning_error.abs() < 1e-10);
    // Seams agree up to the spread gap.
    assert!(
        flat.seams.max_length_mismatch <= 4.0 * report.closure.gap_after + 1e-10,
        "mismatch {} gap {}",
        flat.seams.max_length_mismatch,
        report.closure.gap_after
    );
    assert!(
        flat.corner_uvs
            .iter()
            .all(|p| p[0].is_finite() && p[1].is_finite())
    );
}

#[test]
fn seams_split_the_surface_into_charts_with_matching_lengths() {
    let mut g = generators::grid_disk(8, 6).unwrap();
    generators::perturb(&mut g.positions, [0.02, 0.02, 0.0], 3);
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let at = |i: usize, j: usize| i * 7 + j;
    let seam: Vec<u32> = (0..6)
        .map(|j| mesh.find_edge(at(4, j), at(4, j + 1)).unwrap() as u32)
        .collect();
    let flat = flatten(
        &mesh,
        &geometry,
        &[],
        CutSource::Plan {
            may_cut: None,
            must_cut: &seam,
        },
    )
    .unwrap();
    assert_disk_charts(&flat);
    assert_eq!(flat.n_charts(), 2);
    assert_eq!(flat.cut_edges(), seam.as_slice());
    assert_eq!(flat.seams.n_between_charts, 6);
    // Flat input: both charts are reproduced, at the same scale.
    let original: Vec<[f64; 2]> = mesh
        .faces()
        .iter()
        .flat_map(|f| f.map(|v| [g.positions[v as usize][0], g.positions[v as usize][1]]))
        .collect();
    assert!(max_length_deviation(&mesh, &geometry, &flat.corner_uvs) < 1e-10);
    let agreement = compare_layouts(&flat.face_chart, 2, &flat.corner_uvs, &original);
    assert!(agreement.max < 1e-10, "{agreement:?}");
    // Each chart starts at the origin: nothing is packed.
    for chart in 0..2 {
        let start =
            flat.cut.mesh.boundary_loops().iter().find(|l| {
                flat.face_chart[flat.cut.mesh.outgoing_halfedge(l[0] as usize) / 3] == chart
            });
        let v = start.unwrap()[0] as usize;
        assert_eq!(flat.vertex_uvs[v], [0.0, 0.0]);
    }
    // A curved surface with a seam between two charts: lengths across it still agree.
    let g = generators::bumpy_grid_disk(8, 6, 0.2).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let flat = flatten(
        &mesh,
        &geometry,
        &[],
        CutSource::Plan {
            may_cut: None,
            must_cut: &seam,
        },
    )
    .unwrap();
    assert_eq!(flat.n_charts(), 2);
    assert_eq!(flat.charts[0].closure.coupled_loops, 2);
    assert!(max_seam_mismatch(&mesh, &flat) < 1e-12);
    // Without coupling the two sides differ by the closing adjustment.
    let uncoupled = flatten_with_cones(
        &ConeFlattenInput {
            mesh: &mesh,
            geometry: &geometry,
            cones: &[],
            scale_factors: None,
            cut: CutSource::Plan {
                may_cut: None,
                must_cut: &seam,
            },
            anchors: &[],
        },
        &ConeFlattenOptions {
            max_coupled_charts: 1,
            ..ConeFlattenOptions::default()
        },
    )
    .unwrap();
    assert!(uncoupled.seams.max_length_mismatch > 1e-6);
    assert!(uncoupled.seams.max_length_mismatch < 0.1);
}

#[test]
fn precomputed_scale_factors_and_explicit_cuts() {
    let g = generators::icosphere(2).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let cones: Vec<Cone> = [0u32, 3, 5, 9]
        .iter()
        .map(|&vertex| Cone { vertex, angle: PI })
        .collect();
    let planned = flatten(&mesh, &geometry, &cones, PLAN).unwrap();
    // The same scale factors and the same cut, passed in: the same layout.
    let again = flatten_with_cones(
        &ConeFlattenInput {
            mesh: &mesh,
            geometry: &geometry,
            cones: &cones,
            scale_factors: Some(&planned.scale_factors),
            cut: CutSource::Edges(planned.cut_edges()),
            anchors: &[],
        },
        &ConeFlattenOptions::default(),
    )
    .unwrap();
    assert!(again.scale_factor_solve.is_none() && again.plan.is_none());
    for (a, b) in again.corner_uvs.iter().zip(&planned.corner_uvs) {
        assert!((a[0] - b[0]).abs() < 1e-12 && (a[1] - b[1]).abs() < 1e-12);
    }
    for cone in &again.cones {
        assert!((cone.implied_angle - PI).abs() < 1e-10);
    }
    // solve_scale_factors gives the same u.
    let solved =
        solve_scale_factors(&mesh, &geometry, &cones, &ScaleFactorOptions::default()).unwrap();
    assert_eq!(solved.u, planned.scale_factors);
    // Scale factors with a cone that is not listed: the planner cannot know about it, and
    // the flattening refuses to hide the curvature inside a chart.
    let err = flatten_with_cones(
        &ConeFlattenInput {
            mesh: &mesh,
            geometry: &geometry,
            cones: &cones[..3],
            scale_factors: Some(&planned.scale_factors),
            cut: PLAN,
            anchors: &[],
        },
        &ConeFlattenOptions::default(),
    )
    .unwrap_err();
    match err {
        BffError::ScaleFactorsNotFlat {
            vertex, curvature, ..
        } => {
            assert_eq!(vertex, 9);
            assert!((curvature - PI).abs() < 1e-9);
        }
        other => panic!("unexpected error: {other}"),
    }
    // An anchor moves the layout rigidly.
    let corner = (0..mesh.n_halfedges())
        .find(|&h| mesh.from_vertex(h) == 5)
        .unwrap() as u32;
    let anchored = flatten_with_cones(
        &ConeFlattenInput {
            mesh: &mesh,
            geometry: &geometry,
            cones: &cones,
            scale_factors: None,
            cut: CutSource::Edges(planned.cut_edges()),
            anchors: &[corner],
        },
        &ConeFlattenOptions::default(),
    )
    .unwrap();
    assert_eq!(anchored.corner_uvs[corner as usize], [0.0, 0.0]);
    let agreement = compare_layouts(
        &anchored.face_chart,
        1,
        &anchored.corner_uvs,
        &planned.corner_uvs,
    );
    assert!(
        agreement.max < 1e-10 && (agreement.scale - 1.0).abs() < 1e-10,
        "{agreement:?}"
    );
}

#[test]
fn gauss_bonnet_is_checked_on_closed_surfaces() {
    let g = generators::icosphere(2).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let excess = 1e-3;
    let cones: Vec<Cone> = (0..12)
        .map(|vertex| Cone {
            vertex,
            angle: PI / 3.0 + excess,
        })
        .collect();
    let run = |policy| {
        flatten_with_cones(
            &ConeFlattenInput {
                mesh: &mesh,
                geometry: &geometry,
                cones: &cones,
                scale_factors: None,
                cut: PLAN,
                anchors: &[],
            },
            &ConeFlattenOptions {
                scale_factors: ScaleFactorOptions {
                    gauss_bonnet: policy,
                    ..ScaleFactorOptions::default()
                },
                ..ConeFlattenOptions::default()
            },
        )
    };
    // Default: the discrepancy is reported, not hidden.
    let err = run(GaussBonnetPolicy::Reject).unwrap_err();
    match &err {
        BffError::GaussBonnetViolated {
            n_cones,
            cone_sum,
            expected,
            euler_characteristic,
            ..
        } => {
            assert_eq!((*n_cones, *euler_characteristic), (12, 2));
            assert!((cone_sum - expected - 12.0 * excess).abs() < 1e-12);
        }
        other => panic!("unexpected error: {other}"),
    }
    assert!(err.to_string().contains("Gauss-Bonnet"), "{err}");
    // On request the cones absorb it ...
    let flat = run(GaussBonnetPolicy::DistributeOverCones).unwrap();
    let total: f64 = flat.cones.iter().map(|c| c.angle).sum();
    assert!((total - 4.0 * PI).abs() < 1e-12);
    for cone in &flat.cones {
        assert!((cone.angle - PI / 3.0).abs() < 1e-12);
        assert!((cone.realized_angle - cone.angle).abs() < 1e-9);
    }
    let report = &flat.scale_factor_solve.as_ref().unwrap().closed_components[0];
    assert!((report.discrepancy - 12.0 * excess).abs() < 1e-12);
    assert!(flat.charts[0].closure.turning_error.abs() < 1e-10);
    // ... or all vertices do: the cone angles stay, the charts are slightly curved inside
    // and the boundary curve misses a full turn by the part of the discrepancy that lies
    // inside the chart.
    let flat = run(GaussBonnetPolicy::SpreadOverVertices).unwrap();
    let n = mesh.n_vertices() as f64;
    let per_vertex = 12.0 * excess / n;
    assert!((flat.max_interior_curvature - per_vertex).abs() < 1e-10);
    for cone in &flat.cones {
        assert!((cone.angle - (PI / 3.0 + excess)).abs() < 1e-15);
    }
    let interior = flat.cut.mesh.interior_vertices().len() as f64;
    assert!(
        (flat.charts[0].closure.turning_error.abs() - per_vertex * interior).abs() < 1e-9,
        "{:?}",
        flat.charts[0].closure
    );
    // A surface with boundary has no such condition.
    let g = generators::bumpy_grid_disk(6, 6, 0.1).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let flat = flatten(
        &mesh,
        &geometry,
        &[Cone {
            vertex: 24,
            angle: 0.7,
        }],
        PLAN,
    )
    .unwrap();
    assert!(
        flat.scale_factor_solve
            .unwrap()
            .closed_components
            .is_empty()
    );
}

#[test]
fn bad_input_gives_descriptive_errors() {
    let g = generators::icosphere(1).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let cone = |vertex, angle| Cone { vertex, angle };
    let four = [cone(0, PI), cone(3, PI), cone(5, PI), cone(9, PI)];

    let err = flatten(&mesh, &geometry, &[cone(1000, 1.0)], PLAN).unwrap_err();
    assert_eq!(
        err,
        BffError::ConeVertexOutOfRange {
            position: 0,
            vertex: 1000,
            n_vertices: 42
        }
    );
    assert!(err.to_string().contains("vertex 1000"), "{err}");
    let err = flatten(&mesh, &geometry, &[cone(1, 1.0), cone(1, 2.0)], PLAN).unwrap_err();
    assert!(
        matches!(
            err,
            BffError::DuplicateCone {
                vertex: 1,
                first: 0,
                second: 1
            }
        ),
        "{err}"
    );
    let err = flatten(&mesh, &geometry, &[cone(1, 2.0 * PI)], PLAN).unwrap_err();
    assert!(matches!(err, BffError::InvalidConeAngle { .. }), "{err}");
    let err = flatten(&mesh, &geometry, &[cone(1, f64::NAN)], PLAN).unwrap_err();
    assert!(matches!(err, BffError::InvalidConeAngle { .. }), "{err}");
    // Two cones cannot carry 4π.
    let err = flatten(&mesh, &geometry, &[cone(1, 3.0), cone(2, 3.0)], PLAN).unwrap_err();
    assert!(matches!(err, BffError::GaussBonnetViolated { .. }), "{err}");
    // An explicit cut that does not make a disk, or misses a cone.
    let err = flatten(&mesh, &geometry, &four, CutSource::Edges(&[])).unwrap_err();
    assert!(
        matches!(
            err,
            BffError::ChartNotDisk {
                genus: 0,
                boundary_loops: 0,
                ..
            }
        ),
        "{err}"
    );
    let planned = flatten(&mesh, &geometry, &four[..3], PLAN);
    assert!(matches!(planned, Err(BffError::GaussBonnetViolated { .. })));
    let cut_through_three = ddg_bff::plan_cut(&ddg_bff::CutPlanInput {
        mesh: &mesh,
        edge_lengths: geometry.edge_lengths(),
        cone_vertices: &[0, 3, 5],
        may_cut: None,
        must_cut: &[],
    })
    .unwrap();
    let err = flatten(
        &mesh,
        &geometry,
        &four,
        CutSource::Edges(&cut_through_three.cut_edges),
    )
    .unwrap_err();
    assert!(
        matches!(
            err,
            BffError::ConeNotOnCut {
                vertex: 9,
                chart: 0
            }
        ) || cut_through_three
            .cut_edges
            .iter()
            .any(|&e| mesh.edge_vertices(e as usize).contains(&9)),
        "{err}"
    );
    let err = flatten(&mesh, &geometry, &four, CutSource::Edges(&[100_000])).unwrap_err();
    assert!(matches!(err, BffError::EdgeOutOfRange { .. }), "{err}");
    // Scale factors of the wrong length or with NaN.
    for (u, expect_length) in [(vec![0.0; 5], true), (vec![f64::NAN; 42], false)] {
        let err = flatten_with_cones(
            &ConeFlattenInput {
                mesh: &mesh,
                geometry: &geometry,
                cones: &four,
                scale_factors: Some(&u),
                cut: PLAN,
                anchors: &[],
            },
            &ConeFlattenOptions::default(),
        )
        .unwrap_err();
        if expect_length {
            assert!(matches!(err, BffError::LengthMismatch { .. }), "{err}");
        } else {
            assert!(matches!(err, BffError::NonFiniteInput { .. }), "{err}");
        }
    }
    // A geometry of another mesh.
    let other = generators::icosphere(2).unwrap();
    let (_, other_geometry) = setup(&other.positions, &other.faces);
    let err = flatten(&mesh, &other_geometry, &four, PLAN).unwrap_err();
    assert!(matches!(err, BffError::Geometry(_)), "{err}");
    // Bad anchors.
    let anchored = |anchors: &[u32]| {
        flatten_with_cones(
            &ConeFlattenInput {
                mesh: &mesh,
                geometry: &geometry,
                cones: &four,
                scale_factors: None,
                cut: PLAN,
                anchors,
            },
            &ConeFlattenOptions::default(),
        )
    };
    let err = anchored(&[10_000]).unwrap_err();
    assert!(matches!(err, BffError::AnchorOutOfRange { .. }), "{err}");
    let ok = anchored(&[]).unwrap();
    let interior_corner = (0..mesh.n_halfedges())
        .find(|&h| !ok.cut.mesh.is_boundary_vertex(ok.cut.mesh.from_vertex(h)))
        .unwrap() as u32;
    let err = anchored(&[interior_corner]).unwrap_err();
    assert!(matches!(err, BffError::AnchorNotOnBoundary { .. }), "{err}");
    let on_cut: Vec<u32> = (0..mesh.n_halfedges() as u32)
        .filter(|&h| {
            ok.cut
                .mesh
                .is_boundary_vertex(ok.cut.mesh.from_vertex(h as usize))
        })
        .take(2)
        .collect();
    let err = anchored(&on_cut).unwrap_err();
    assert!(matches!(err, BffError::DuplicateAnchor { .. }), "{err}");
    // A cone on the boundary of a disk.
    let g = generators::grid_disk(4, 4).unwrap();
    let (mesh, geometry) = setup(&g.positions, &g.faces);
    let err = flatten(&mesh, &geometry, &[cone(0, 0.5)], PLAN).unwrap_err();
    assert!(
        matches!(err, BffError::ConeOnBoundary { vertex: 0, .. }),
        "{err}"
    );
}

#[test]
fn a_degenerate_mesh_is_flattened_after_mollification() {
    // A grid with one vertex moved onto the line through two neighbours: a zero-area face.
    let mut g = generators::grid_disk(6, 6).unwrap();
    let at = |i: usize, j: usize| i * 7 + j;
    // Face (3,3) (4,3) (4,4) or its mirror: put (4,4) on the line of the other two.
    let (a, b) = (g.positions[at(3, 3)], g.positions[at(4, 3)]);
    g.positions[at(4, 4)] = [b[0] + (b[0] - a[0]) * 0.5, b[1] + (b[1] - a[1]) * 0.5, 0.0];
    let mesh = g.topology().unwrap();
    assert!(Geometry::from_positions(&mesh, &g.positions).is_err());
    let lengths: Vec<f64> = mesh
        .edges()
        .iter()
        .map(|&[p, q]| edge_length(g.positions[p as usize], g.positions[q as usize]))
        .collect();
    let mollified = mollify_edge_lengths(&mesh, &lengths, 1e-4).unwrap();
    assert!(mollified.epsilon > 0.0);
    let geometry = Geometry::from_edge_lengths(&mesh, mollified.edge_lengths).unwrap();
    let mut bff = DiskBff::new(&mesh, &geometry, DiskOptions::default()).unwrap();
    let flat = bff.flatten_automatic().unwrap();
    assert!(flat.uv.iter().all(|p| p[0].is_finite() && p[1].is_finite()));
    assert_eq!(flat.report.distortion.n_faces, mesh.n_faces());
}
