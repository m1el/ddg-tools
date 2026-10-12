//! Property tests of cutting, cut planning and flattening on generated meshes.

// Index loops mirror the subscripts of the formulas being checked.
#![allow(clippy::needless_range_loop)]

mod common;

use std::f64::consts::PI;

use common::{layout_angle_sums, max_length_deviation};
use ddg_bff::{
    BffError, Cone, ConeFlattenInput, ConeFlattenOptions, ConeFlattening, CutPlanInput, CutSource,
    component_topology, cut_euler_characteristic, cut_mesh, flatten_with_cones, plan_cut,
};
use ddg_core::generators::{self, GeneratedMesh, SplitMix64};
use ddg_core::{Geometry, TriMesh};
use proptest::prelude::*;

/// A small mesh of one of the topology classes: sphere, torus, disk, annulus.
fn arbitrary_mesh() -> impl Strategy<Value = (GeneratedMesh, &'static str)> {
    prop_oneof![
        (0u32..3).prop_map(|level| (generators::icosphere(level).unwrap(), "sphere")),
        (3usize..9, 3usize..7)
            .prop_map(|(a, b)| (generators::torus(a, b, 2.0, 0.7).unwrap(), "torus")),
        (1usize..7, 1usize..7, 0.0..0.3f64)
            .prop_map(|(a, b, amp)| (generators::bumpy_grid_disk(a, b, amp).unwrap(), "disk")),
        (3usize..10, 1usize..4)
            .prop_map(|(a, b)| (generators::annulus(a, b, 1.0, 2.0).unwrap(), "annulus")),
    ]
}

fn setup(g: &GeneratedMesh) -> (TriMesh, Geometry) {
    let mesh = g.topology().unwrap();
    let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
    (mesh, geometry)
}

/// `k` distinct interior vertices, chosen by a seeded shuffle.
fn pick_interior(mesh: &TriMesh, k: usize, rng: &mut SplitMix64) -> Vec<u32> {
    let mut interior = mesh.interior_vertices();
    for i in (1..interior.len()).rev() {
        interior.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
    }
    interior.truncate(k);
    interior
}

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

/// The structural facts every successful flattening has to satisfy.
fn check_structure(mesh: &TriMesh, flat: &ConeFlattening) -> Result<(), TestCaseError> {
    let topology = component_topology(&flat.cut.mesh);
    prop_assert_eq!(topology.len(), flat.n_charts());
    for t in &topology {
        prop_assert!(t.is_disk(), "{:?}", t);
        prop_assert_eq!(t.euler_characteristic, 1);
        prop_assert_eq!(t.boundary_loops, 1);
    }
    prop_assert_eq!(flat.corner_uvs.len(), 3 * mesh.n_faces());
    prop_assert!(
        flat.corner_uvs
            .iter()
            .all(|p| p[0].is_finite() && p[1].is_finite())
    );
    // Corner UVs are the vertex UVs of the cut mesh, whose faces are the input faces.
    for (f, face) in flat.cut.mesh.faces().iter().enumerate() {
        for k in 0..3 {
            prop_assert_eq!(
                flat.corner_uvs[3 * f + k],
                flat.vertex_uvs[face[k] as usize]
            );
            prop_assert_eq!(
                flat.cut.vertex_to_original[face[k] as usize],
                mesh.faces()[f][k]
            );
        }
    }
    // Every cone is on the cut.
    for cone in &flat.cones {
        prop_assert!(cone.n_wedges >= 1);
    }
    for (v, &original) in flat.cut.vertex_to_original.iter().enumerate() {
        if flat.cones.iter().any(|c| c.vertex == original) {
            prop_assert!(flat.cut.mesh.is_boundary_vertex(v));
        }
    }
    // Lengths agree across the cut whenever closing did not have to give up.
    if flat
        .charts
        .iter()
        .all(|c| c.closure.dropped_directions == 0)
    {
        prop_assert!(
            flat.seams.max_length_mismatch < 1e-9,
            "seam mismatch {}",
            flat.seams.max_length_mismatch
        );
        for chart in &flat.charts {
            prop_assert!(chart.closure.gap_after < 1e-9, "{:?}", chart.closure);
        }
    }
    Ok(())
}

/// Number of cases per property: 48, or `PROPTEST_CASES` from the environment.
fn cases() -> u32 {
    std::env::var("PROPTEST_CASES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(48)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(cases()))]

    /// Cutting along any edge set keeps faces and corners, and changes the Euler
    /// characteristic by (added vertices) − (cut edges).
    #[test]
    fn cutting_obeys_the_euler_bookkeeping(
        (g, _) in arbitrary_mesh(),
        seed in any::<u64>(),
        density in 0.0..0.6f64,
    ) {
        let (mesh, _) = setup(&g);
        let mut rng = SplitMix64::new(seed);
        let edges: Vec<u32> = (0..mesh.n_edges() as u32)
            .filter(|_| rng.next_f64() < density)
            .collect();
        match cut_mesh(&mesh, &edges) {
            Err(BffError::IsolatedCutEdge { edge, .. }) => {
                // Really isolated: no other requested edge touches its end points, and they
                // are interior.
                let [a, b] = mesh.edge_vertices(edge as usize);
                for v in [a, b] {
                    prop_assert!(!mesh.is_boundary_vertex(v));
                    let others = edges
                        .iter()
                        .filter(|&&e| e != edge && mesh.edge_vertices(e as usize).contains(&v))
                        .count();
                    prop_assert_eq!(others, 0);
                }
            }
            Err(other) => prop_assert!(false, "unexpected error: {}", other),
            Ok(cut) => {
                prop_assert_eq!(cut.mesh.n_faces(), mesh.n_faces());
                let interior: Vec<u32> = edges
                    .iter()
                    .copied()
                    .filter(|&e| !mesh.is_boundary_edge(e as usize))
                    .collect();
                prop_assert_eq!(&cut.cut_edges, &interior);
                prop_assert_eq!(cut.ignored_boundary_edges, edges.len() - interior.len());
                prop_assert_eq!(
                    cut.mesh.euler_characteristic(),
                    cut_euler_characteristic(
                        mesh.euler_characteristic(),
                        cut.n_added_vertices(),
                        cut.cut_edges.len()
                    )
                );
                prop_assert_eq!(cut.mesh.n_edges(), mesh.n_edges() + interior.len());
                prop_assert_eq!(
                    cut.mesh.n_boundary_edges(),
                    mesh.n_boundary_edges() + 2 * interior.len()
                );
                for (f, face) in cut.mesh.faces().iter().enumerate() {
                    for k in 0..3 {
                        prop_assert_eq!(cut.vertex_to_original[face[k] as usize], mesh.faces()[f][k]);
                    }
                }
                for (e, &original) in cut.edge_to_original.iter().enumerate() {
                    let [a, b] = cut.mesh.edge_vertices(e);
                    let mut mapped = [cut.vertex_to_original[a], cut.vertex_to_original[b]];
                    mapped.sort_unstable();
                    prop_assert_eq!(mapped, mesh.edges()[original as usize]);
                    // An edge is on the boundary of the cut mesh iff it was cut or already
                    // on the boundary.
                    let was = mesh.is_boundary_edge(original as usize)
                        || interior.binary_search(&original).is_ok();
                    prop_assert_eq!(cut.mesh.is_boundary_edge(e), was);
                }
                let copies = cut.copies_per_vertex();
                prop_assert_eq!(copies.iter().map(|&c| c as usize).sum::<usize>(), cut.mesh.n_vertices());
            }
        }
    }

    /// Whatever the constraints, a plan is either a valid cut that honours them or one of the
    /// documented refusals.
    #[test]
    fn planned_cuts_are_disks_and_honour_the_constraints(
        (g, kind) in arbitrary_mesh(),
        seed in any::<u64>(),
        n_cones in 0usize..6,
        forbid in prop_oneof![Just(0.0), 0.0..0.5f64],
        require in prop_oneof![Just(0.0), 0.0..0.2f64],
    ) {
        let (mesh, geometry) = setup(&g);
        let mut rng = SplitMix64::new(seed);
        let cones = pick_interior(&mesh, n_cones, &mut rng);
        let may: Vec<bool> = (0..mesh.n_edges()).map(|_| rng.next_f64() >= forbid).collect();
        let must: Vec<u32> = (0..mesh.n_edges() as u32)
            .filter(|&e| may[e as usize] && rng.next_f64() < require)
            .collect();
        let input = CutPlanInput {
            mesh: &mesh,
            edge_lengths: geometry.edge_lengths(),
            cone_vertices: &cones,
            may_cut: Some(&may),
            must_cut: &must,
        };
        match plan_cut(&input) {
            Ok(plan) => {
                // Deterministic.
                prop_assert_eq!(&plan_cut(&input).unwrap(), &plan);
                // Constraints.
                for &e in &plan.cut_edges {
                    prop_assert!(may[e as usize], "forbidden edge {} was cut", e);
                    prop_assert!(!mesh.is_boundary_edge(e as usize));
                }
                for &e in &must {
                    prop_assert!(
                        mesh.is_boundary_edge(e as usize) || plan.cut_edges.binary_search(&e).is_ok(),
                        "required edge {} was not cut", e
                    );
                }
                prop_assert_eq!(plan.cut_edges.len(), plan.n_required + plan.n_chosen);
                // Verified independently: cut and look.
                let cut = cut_mesh(&mesh, &plan.cut_edges).unwrap();
                let topology = component_topology(&cut.mesh);
                prop_assert_eq!(topology.len(), plan.charts.len());
                for (t, chart) in topology.iter().zip(&plan.charts) {
                    prop_assert!(t.is_disk(), "{:?}", t);
                    prop_assert_eq!(t, &chart.after);
                }
                for (v, &original) in cut.vertex_to_original.iter().enumerate() {
                    if cones.contains(&original) {
                        prop_assert!(cut.mesh.is_boundary_vertex(v), "cone {} not on the cut", original);
                    }
                }
                // The planner never splits a piece: the charts are the pieces of the required
                // cuts (which the planner may have extended by an edge to make them cuttable).
                if must.is_empty() {
                    prop_assert_eq!(plan.charts.len(), mesh.n_components());
                }
            }
            Err(
                BffError::SphereNeedsTwoCones { .. }
                | BffError::ConeEnclosedByUncuttableEdges { .. }
                | BffError::UncuttableRegionNotDisk { .. }
                | BffError::IsolatedCutEdge { .. },
            ) => {
                // Without constraints only a sphere with fewer than two cones can be refused.
                if forbid == 0.0 && require == 0.0 {
                    prop_assert!(kind == "sphere" && n_cones < 2, "{} with {} cones", kind, n_cones);
                }
            }
            Err(other) => prop_assert!(false, "unexpected error: {}", other),
        }
    }

    /// Spheres with positive cones that satisfy Gauss–Bonnet: one disk, seamless, and the
    /// angle sums in the plane are 2π − µ at the cones and 2π elsewhere.
    #[test]
    fn spheres_with_cones_flatten_seamlessly(
        level in 1u32..3,
        seed in any::<u64>(),
        n_cones in 3usize..9,
        noise in 0.0..0.05f64,
    ) {
        let mut g = generators::icosphere(level).unwrap();
        generators::perturb(&mut g.positions, [noise; 3], seed ^ 1);
        let (mesh, geometry) = setup(&g);
        let mut rng = SplitMix64::new(seed);
        let vertices = pick_interior(&mesh, n_cones, &mut rng);
        // Random positive angles that sum to 4π, each below 2π.
        let weights: Vec<f64> = (0..n_cones).map(|_| rng.uniform(0.6, 1.4)).collect();
        let total: f64 = weights.iter().sum();
        let cones: Vec<Cone> = vertices
            .iter()
            .zip(&weights)
            .map(|(&vertex, w)| Cone { vertex, angle: 4.0 * PI * w / total })
            .collect();
        prop_assume!(cones.iter().all(|c| c.angle < 1.9 * PI));
        let flat = flatten(&mesh, &geometry, &cones, CutSource::Plan { may_cut: None, must_cut: &[] });
        let flat = match flat {
            Ok(flat) => flat,
            Err(e) => return Err(TestCaseError::fail(format!("flatten failed: {e}"))),
        };
        check_structure(&mesh, &flat)?;
        prop_assert_eq!(flat.n_charts(), 1);
        prop_assert!(flat.charts[0].closure.turning_error.abs() < 1e-9);
        // The boundary curve has exactly the prescribed angles unless closing it by lengths
        // had to be given up (extreme cone angles on a coarse mesh), in which case the
        // remaining gap was spread along the curve.
        let closed_by_lengths = flat.charts[0].closure.dropped_directions == 0;
        for cone in &flat.cones {
            prop_assert!((cone.implied_angle - cone.angle).abs() < 1e-9);
            // Modulo 2π in any case; exactly when no triangle is flipped.
            let off = (cone.realized_angle - cone.angle) / (2.0 * PI);
            prop_assert!(!closed_by_lengths || (off - off.round()).abs() < 1e-8, "{:?}", cone);
        }
        if closed_by_lengths && flat.charts[0].distortion.flipped_triangles == 0 {
            let sums = layout_angle_sums(&mesh, &flat.corner_uvs);
            for v in 0..mesh.n_vertices() {
                let mu = cones.iter().find(|c| c.vertex as usize == v).map_or(0.0, |c| c.angle);
                prop_assert!((sums[v] - (2.0 * PI - mu)).abs() < 1e-7, "vertex {}: {}", v, sums[v]);
            }
        }
    }

    /// Bumpy disks with a few cones of either sign.
    #[test]
    fn bumpy_disks_with_cones(
        nx in 4usize..10,
        ny in 4usize..10,
        amplitude in 0.0..0.25f64,
        seed in any::<u64>(),
        n_cones in 0usize..4,
    ) {
        let g = generators::bumpy_grid_disk(nx, ny, amplitude).unwrap();
        let (mesh, geometry) = setup(&g);
        let mut rng = SplitMix64::new(seed);
        let cones: Vec<Cone> = pick_interior(&mesh, n_cones, &mut rng)
            .into_iter()
            .map(|vertex| Cone { vertex, angle: rng.uniform(-1.0, 1.5) })
            .collect();
        let flat = flatten(&mesh, &geometry, &cones, CutSource::Plan { may_cut: None, must_cut: &[] });
        let flat = match flat {
            Ok(flat) => flat,
            Err(e) => return Err(TestCaseError::fail(format!("flatten failed: {e}"))),
        };
        check_structure(&mesh, &flat)?;
        prop_assert_eq!(flat.n_charts(), 1);
        prop_assert!(flat.max_interior_curvature < 1e-9);
        let closed_by_lengths = flat.charts[0].closure.dropped_directions == 0;
        for cone in &flat.cones {
            let off = (cone.realized_angle - cone.angle) / (2.0 * PI);
            prop_assert!(!closed_by_lengths || (off - off.round()).abs() < 1e-8, "{:?}", cone);
            // A negative cone at the end of a cut is reported as overlapping.
            if closed_by_lengths
                && cone.angle < -1e-6
                && cone.n_wedges == 1
                && flat.charts[0].distortion.flipped_triangles == 0
            {
                prop_assert!(cone.max_wedge_angle > 2.0 * PI);
                prop_assert!(flat.charts[0].overlap.wedges_over_2pi >= 1);
            }
        }
        if cones.is_empty() {
            prop_assert!(flat.cut_edges().is_empty());
            if amplitude == 0.0 {
                prop_assert!(max_length_deviation(&mesh, &geometry, &flat.corner_uvs) < 1e-9);
            }
        }
    }

    /// Tori without cones become one disk with two handle loops; flat annuli are unfolded
    /// isometrically.
    #[test]
    fn tori_and_annuli_without_cones(
        n_major in 5usize..16,
        n_minor in 4usize..10,
        seed in any::<u64>(),
        noise in 0.0..0.03f64,
    ) {
        let mut g = generators::torus(n_major, n_minor, 2.0, 0.7).unwrap();
        generators::perturb(&mut g.positions, [noise; 3], seed);
        let (mesh, geometry) = setup(&g);
        let flat = flatten(&mesh, &geometry, &[], CutSource::Plan { may_cut: None, must_cut: &[] });
        let flat = match flat {
            Ok(flat) => flat,
            Err(e) => return Err(TestCaseError::fail(format!("torus: {e}"))),
        };
        check_structure(&mesh, &flat)?;
        prop_assert_eq!(flat.n_charts(), 1);
        prop_assert_eq!(flat.plan.as_ref().unwrap().n_loop_edges, 2);
        prop_assert!(flat.charts[0].closure.turning_error.abs() < 1e-9);
        prop_assert!(flat.charts[0].closure.max_length_change <= 0.5);

        let mut g = generators::annulus(n_major, n_minor - 3, 1.0, 2.5).unwrap();
        generators::perturb(&mut g.positions, [0.3 * noise, 0.3 * noise, 0.0], seed);
        let (mesh, geometry) = setup(&g);
        let flat = flatten(&mesh, &geometry, &[], CutSource::Plan { may_cut: None, must_cut: &[] });
        let flat = match flat {
            Ok(flat) => flat,
            Err(e) => return Err(TestCaseError::fail(format!("annulus: {e}"))),
        };
        check_structure(&mesh, &flat)?;
        prop_assert_eq!(flat.n_charts(), 1);
        prop_assert!(max_length_deviation(&mesh, &geometry, &flat.corner_uvs) < 1e-8);
        prop_assert_eq!(flat.charts[0].distortion.flipped_triangles, 0);
    }

    /// Whatever is thrown at it — cones of any angle and sign, on any vertex, repeated or out
    /// of range, random masks, seams and anchors, every Gauss–Bonnet policy, with and without
    /// the diagonal shift — the flattening returns a finite result or an error with a
    /// message. It never panics and never returns NaN. (Of 20 000 cases run once, about a
    /// quarter ended in a layout and three quarters in an error.)
    #[test]
    fn nothing_panics_and_nothing_is_nan(
        (g, _) in arbitrary_mesh(),
        doubled in any::<bool>(),
        seed in any::<u64>(),
        n_cones in 0usize..7,
        policy in 0usize..3,
        shifted in any::<bool>(),
        wild in any::<bool>(),
    ) {
        // Half of the time the connected sum with the mirror image: twice the genus. (Where
        // the removed face touches the boundary the sum can be pinched; then the mesh itself.)
        let (mut positions, mut faces) = (g.positions.clone(), g.faces.clone());
        if doubled {
            let (p, f) = common::mirror_double(&g.positions, &g.faces, 0);
            if TriMesh::new(p.len(), f.clone()).is_ok() {
                (positions, faces) = (p, f);
            }
        }
        let mesh = TriMesh::new(positions.len(), faces).unwrap();
        let geometry = Geometry::from_positions(&mesh, &positions).unwrap();
        let mut rng = SplitMix64::new(seed);
        let n = mesh.n_vertices() as u64;
        // Tame input: distinct interior vertices. Wild input: any vertex, boundary or not,
        // repeated or (rarely) one that does not exist.
        let tame = pick_interior(&mesh, n_cones, &mut rng);
        let cones: Vec<Cone> = (0..n_cones.min(if wild { n_cones } else { tame.len() }))
            .map(|k| Cone {
                vertex: if wild { (rng.next_u64() % (n + 1)) as u32 } else { tame[k] },
                angle: match rng.next_u64() % 8 {
                    0 if wild => f64::NAN,
                    1 if wild => 7.0,
                    2 => rng.uniform(-12.0, -3.0),
                    3 => rng.uniform(5.0, 2.0 * PI - 1e-3),
                    _ => rng.uniform(-2.0, 3.0),
                },
            })
            .collect();
        let forbid = if rng.next_u64().is_multiple_of(2) { 0.0 } else { rng.uniform(0.0, 0.4) };
        let may: Vec<bool> = (0..mesh.n_edges()).map(|_| rng.next_f64() >= forbid).collect();
        let require = if rng.next_u64().is_multiple_of(2) { 0.0 } else { rng.uniform(0.0, 0.2) };
        let must: Vec<u32> = (0..mesh.n_edges() as u32)
            .filter(|&e| rng.next_f64() < require && (wild || may[e as usize]))
            .collect();
        let anchors: Vec<u32> = (0..rng.next_u64() % 3)
            .map(|_| (rng.next_u64() % (mesh.n_halfedges() as u64 + u64::from(wild))) as u32)
            .collect();
        let options = ConeFlattenOptions {
            scale_factors: ddg_bff::ScaleFactorOptions {
                gauss_bonnet: [
                    ddg_bff::GaussBonnetPolicy::Reject,
                    ddg_bff::GaussBonnetPolicy::DistributeOverCones,
                    ddg_bff::GaussBonnetPolicy::SpreadOverVertices,
                ][policy],
                diagonal_shift: if shifted { 1e-8 } else { 0.0 },
                ..ddg_bff::ScaleFactorOptions::default()
            },
            ..ConeFlattenOptions::default()
        };
        let result = flatten_with_cones(
            &ConeFlattenInput {
                mesh: &mesh,
                geometry: &geometry,
                cones: &cones,
                scale_factors: None,
                cut: CutSource::Plan { may_cut: Some(&may), must_cut: &must },
                anchors: &anchors,
            },
            &options,
        );
        match result {
            Err(error) => prop_assert!(!error.to_string().is_empty()),
            Ok(flat) => {
                prop_assert!(flat.corner_uvs.iter().all(|p| p[0].is_finite() && p[1].is_finite()));
                prop_assert!(flat.scale_factors.iter().all(|u| u.is_finite()));
                prop_assert!(flat.max_interior_curvature.is_finite());
                prop_assert!(flat.extension_residual.is_finite());
                prop_assert!(flat.seams.max_length_mismatch.is_finite());
                for chart in &flat.charts {
                    let (c, d) = (&chart.closure, &chart.distortion);
                    for x in [
                        c.turning_error, c.gap_before, c.gap_after, c.max_length_change,
                        c.min_length_ratio, d.surface_area, d.plane_area, d.quasi_conformal.min,
                        d.quasi_conformal.max, d.quasi_conformal.mean, d.log_scale.min,
                        d.log_scale.max, d.log_scale.mean, d.log_scale_deviation, d.flipped_area,
                    ] {
                        prop_assert!(x.is_finite(), "{:?}", chart);
                    }
                    if let Some(coverage) = chart.overlap.coverage {
                        prop_assert!(coverage.overlapped_fraction().is_finite());
                        prop_assert!(coverage.overlapped_area <= coverage.covered_area);
                    }
                }
                for cone in &flat.cones {
                    prop_assert!(
                        cone.angle.is_finite()
                            && cone.implied_angle.is_finite()
                            && cone.realized_angle.is_finite()
                            && cone.max_wedge_angle.is_finite(),
                        "{:?}", cone
                    );
                }
                // A result is always made of disks with the cones on their boundaries.
                check_structure(&mesh, &flat)?;
            }
        }
    }

    /// Random seams that split the surface into several charts: every chart is flattened and
    /// lengths agree across all seams.
    #[test]
    fn random_seams_split_the_surface(
        closed in any::<bool>(),
        seed in any::<u64>(),
        density in 0.05..0.35f64,
    ) {
        let g = if closed {
            generators::icosphere(2).unwrap()
        } else {
            generators::bumpy_grid_disk(8, 7, 0.15).unwrap()
        };
        let (mesh, geometry) = setup(&g);
        let mut rng = SplitMix64::new(seed);
        let must: Vec<u32> = (0..mesh.n_edges() as u32)
            .filter(|_| rng.next_f64() < density)
            .collect();
        // Cones: on the sphere four of angle π (Gauss–Bonnet), on the disk two small ones.
        let cones: Vec<Cone> = if closed {
            pick_interior(&mesh, 4, &mut rng)
                .into_iter()
                .map(|vertex| Cone { vertex, angle: PI })
                .collect()
        } else {
            pick_interior(&mesh, 2, &mut rng)
                .into_iter()
                .map(|vertex| Cone { vertex, angle: rng.uniform(-0.5, 1.0) })
                .collect()
        };
        let flat = flatten(&mesh, &geometry, &cones, CutSource::Plan { may_cut: None, must_cut: &must });
        let flat = match flat {
            Ok(flat) => flat,
            Err(e) => return Err(TestCaseError::fail(format!("flatten failed: {e}"))),
        };
        check_structure(&mesh, &flat)?;
        for &e in &must {
            prop_assert!(
                mesh.is_boundary_edge(e as usize) || flat.cut_edges().binary_search(&e).is_ok()
            );
        }
        // The charts are exactly the pieces that the seams (as extended by the planner) make.
        let plan = flat.plan.as_ref().unwrap();
        prop_assert_eq!(plan.charts.len(), flat.n_charts());
        prop_assert!(flat.max_interior_curvature < 1e-9);
        prop_assert!(flat.seams.n_between_charts <= flat.seams.n_cut_edges);
        for (chart, report) in flat.charts.iter().enumerate() {
            prop_assert!(report.closure.turning_error.abs() < 1e-8, "chart {}: {:?}", chart, report.closure);
            prop_assert_eq!(report.distortion.n_faces, plan.charts[chart].after.n_faces);
        }
    }
}
