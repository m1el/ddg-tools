//! What the cut planner guarantees, on every class of topology it claims to handle: closed
//! surfaces and surfaces with one or several boundary loops, of genus 0, 1, 2 and 4, connected
//! or not, with forbidden edges, required edges and cones.
//!
//! Three kinds of test:
//!
//! 1. Without constraints the planner must succeed on every surface of the zoo, and the result
//!    is checked from scratch (cut, then look at the pieces).
//! 2. With random forbidden and required edges a plan must honour them, and a refusal must be
//!    one of the documented ones *and be true*: the stated cause is verified independently.
//! 3. On meshes small enough to try every subset of the cuttable edges, the planner is compared
//!    with the truth: it finds a cut exactly when one exists.

mod common;

use common::{mirror_double, punch_holes};
use ddg_bff::{
    BffError, CutPlan, CutPlanInput, component_topology, cut_euler_characteristic, cut_mesh,
    plan_cut,
};
use ddg_core::generators::{self, SplitMix64};
use ddg_core::{Geometry, TriMesh};

/// A surface of the zoo with what is known about its topology by construction.
struct Surface {
    name: String,
    mesh: TriMesh,
    lengths: Vec<f64>,
    components: usize,
    /// Sum of the genera of the components.
    genus: usize,
    boundary_loops: usize,
    /// Whether some component is a closed surface of genus 0 (it needs two cones).
    has_closed_sphere: bool,
}

fn surface(
    name: &str,
    positions: &[[f64; 3]],
    faces: &[[u32; 3]],
    components: usize,
    genus: usize,
    boundary_loops: usize,
) -> Surface {
    let mesh =
        TriMesh::new(positions.len(), faces.to_vec()).unwrap_or_else(|e| panic!("{name}: {e}"));
    let geometry =
        Geometry::from_positions(&mesh, positions).unwrap_or_else(|e| panic!("{name}: {e}"));
    // The zoo is what it says: checked against ddg-core's own counts.
    assert_eq!(mesh.n_components(), components, "{name}: components");
    assert_eq!(mesh.genus(), genus, "{name}: genus");
    assert_eq!(mesh.boundary_loops().len(), boundary_loops, "{name}: loops");
    let has_closed_sphere = component_topology(&mesh)
        .iter()
        .any(|t| t.is_closed() && t.genus == 0);
    Surface {
        name: name.to_string(),
        mesh,
        lengths: geometry.edge_lengths().to_vec(),
        components,
        genus,
        boundary_loops,
        has_closed_sphere,
    }
}

/// Closed surfaces of genus 0, 1, 2, 4; the same with one to four holes; disks with holes; an
/// annulus; a surface with two components.
fn zoo() -> Vec<Surface> {
    let mut out = Vec::new();
    let sphere = generators::icosphere(2).unwrap();
    out.push(surface("sphere", &sphere.positions, &sphere.faces, 1, 0, 0));
    for holes in [1, 2, 4] {
        let faces = punch_holes(sphere.positions.len(), &sphere.faces, holes, 7).unwrap();
        out.push(surface(
            &format!("sphere with {holes} holes"),
            &sphere.positions,
            &faces,
            1,
            0,
            holes,
        ));
    }
    let torus = generators::torus(10, 7, 2.0, 0.7).unwrap();
    out.push(surface("torus", &torus.positions, &torus.faces, 1, 1, 0));
    for holes in [1, 3] {
        let faces = punch_holes(torus.positions.len(), &torus.faces, holes, 3).unwrap();
        out.push(surface(
            &format!("torus with {holes} holes"),
            &torus.positions,
            &faces,
            1,
            1,
            holes,
        ));
    }
    let (p2, f2) = mirror_double(&torus.positions, &torus.faces, 0);
    out.push(surface("genus 2", &p2, &f2, 1, 2, 0));
    let holes2 = punch_holes(p2.len(), &f2, 2, 5).unwrap();
    out.push(surface("genus 2 with 2 holes", &p2, &holes2, 1, 2, 2));
    let (p4, f4) = mirror_double(&p2, &f2, f2.len() / 3);
    out.push(surface("genus 4", &p4, &f4, 1, 4, 0));
    // A torus with a hole, doubled: genus 2 with two holes, by another route.
    let one_hole = punch_holes(torus.positions.len(), &torus.faces, 1, 11).unwrap();
    let (pd, fd) = mirror_double(&torus.positions, &one_hole, 5);
    out.push(surface("doubled torus with a hole", &pd, &fd, 1, 2, 2));
    let disk = generators::bumpy_grid_disk(7, 6, 0.1).unwrap();
    out.push(surface("disk", &disk.positions, &disk.faces, 1, 0, 1));
    let grid = generators::grid_disk(10, 9).unwrap();
    for holes in [1, 3] {
        let faces = punch_holes(grid.positions.len(), &grid.faces, holes, 2).unwrap();
        out.push(surface(
            &format!("disk with {holes} holes"),
            &grid.positions,
            &faces,
            1,
            0,
            1 + holes,
        ));
    }
    let annulus = generators::annulus(12, 3, 1.0, 2.0).unwrap();
    out.push(surface(
        "annulus",
        &annulus.positions,
        &annulus.faces,
        1,
        0,
        2,
    ));
    // Two components in one mesh: a torus and a sphere with two holes.
    let mut positions = torus.positions.clone();
    let mut faces = torus.faces.clone();
    let offset = positions.len() as u32;
    positions.extend(sphere.positions.iter().map(|p| [p[0] + 10.0, p[1], p[2]]));
    let punched = punch_holes(sphere.positions.len(), &sphere.faces, 2, 9).unwrap();
    faces.extend(punched.iter().map(|f| f.map(|v| v + offset)));
    out.push(surface(
        "torus and sphere with 2 holes",
        &positions,
        &faces,
        2,
        1,
        2,
    ));
    out
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

/// Checks a plan from scratch: cuts the mesh along it and looks at the result. Returns the
/// number of charts.
fn check_plan(
    s: &Surface,
    plan: &CutPlan,
    cones: &[u32],
    may: Option<&[bool]>,
    must: &[u32],
) -> usize {
    let mesh = &s.mesh;
    let name = &s.name;
    let cut = cut_mesh(mesh, &plan.cut_edges).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(cut.cut_edges, plan.cut_edges, "{name}");
    // Every chart is a disk: connected (by definition of a component), Euler characteristic 1,
    // one boundary loop, genus 0.
    let topology = component_topology(&cut.mesh);
    for (chart, t) in topology.iter().enumerate() {
        assert_eq!(t.euler_characteristic, 1, "{name}: chart {chart}: {t:?}");
        assert_eq!(t.boundary_loops, 1, "{name}: chart {chart}: {t:?}");
        assert_eq!(t.genus, 0, "{name}: chart {chart}: {t:?}");
    }
    assert_eq!(cut.mesh.boundary_loops().len(), topology.len(), "{name}");
    assert_eq!(cut.mesh.euler_characteristic(), topology.len() as i64);
    assert_eq!(
        cut.mesh.euler_characteristic(),
        cut_euler_characteristic(
            mesh.euler_characteristic(),
            cut.n_added_vertices(),
            cut.cut_edges.len()
        ),
        "{name}"
    );
    // What the plan reports about the charts is what cutting gives.
    assert_eq!(plan.charts.len(), topology.len(), "{name}");
    for (chart, t) in plan.charts.iter().zip(&topology) {
        assert_eq!(&chart.after, t, "{name}");
    }
    // Every copy of every cone is on the boundary of its chart.
    for (v, &original) in cut.vertex_to_original.iter().enumerate() {
        if cones.contains(&original) {
            assert!(
                cut.mesh.is_boundary_vertex(v),
                "{name}: cone {original} is inside a chart"
            );
        }
    }
    // No forbidden edge is cut, no boundary edge is listed, every required edge is cut.
    for &e in &plan.cut_edges {
        assert!(!mesh.is_boundary_edge(e as usize), "{name}: edge {e}");
        if let Some(may) = may {
            assert!(may[e as usize], "{name}: forbidden edge {e} was cut");
        }
    }
    for &e in must {
        assert!(
            mesh.is_boundary_edge(e as usize) || plan.cut_edges.binary_search(&e).is_ok(),
            "{name}: required edge {e} was not cut"
        );
    }
    assert_eq!(
        plan.cut_edges.len(),
        plan.n_required + plan.n_chosen,
        "{name}"
    );
    topology.len()
}

#[test]
fn every_topology_class_is_cut_into_disks() {
    for s in zoo() {
        let mesh = &s.mesh;
        for seed in 0..12u64 {
            let mut rng = SplitMix64::new(seed);
            // 0, 1, 2, 3 or 6 cones; a closed sphere needs two.
            let mut n_cones = [0, 1, 2, 3, 6][(seed % 5) as usize];
            if s.has_closed_sphere {
                n_cones = n_cones.max(2);
            }
            let cones = pick_interior(mesh, n_cones, &mut rng);
            let input = CutPlanInput {
                mesh,
                edge_lengths: &s.lengths,
                cone_vertices: &cones,
                may_cut: None,
                must_cut: &[],
            };
            let plan =
                plan_cut(&input).unwrap_or_else(|e| panic!("{} with {cones:?}: {e}", s.name));
            // Deterministic.
            assert_eq!(plan_cut(&input).unwrap(), plan, "{}", s.name);
            let n_charts = check_plan(&s, &plan, &cones, None, &[]);
            // Without required edges nothing is split: one chart per component.
            assert_eq!(n_charts, s.components, "{}", s.name);
            assert_eq!(plan.n_required, 0);
            // Two loops per handle, and the plan knows the topology it started from.
            assert_eq!(plan.n_loop_edges, 2 * s.genus, "{}", s.name);
            assert_eq!(
                plan.charts.iter().map(|c| c.before.genus).sum::<usize>(),
                s.genus
            );
            assert_eq!(
                plan.charts
                    .iter()
                    .map(|c| c.before.boundary_loops)
                    .sum::<usize>(),
                s.boundary_loops
            );
            assert_eq!(
                plan.charts
                    .iter()
                    .map(|c| c.n_interior_cones)
                    .sum::<usize>(),
                cones.len()
            );
            // Nothing superfluous: a cut edge that ends at an interior vertex without another
            // cut edge ends at a cone (or is the extension of a one-edge slit).
            let mut degree = vec![0u32; mesh.n_vertices()];
            for &e in &plan.cut_edges {
                for v in mesh.edge_vertices(e as usize) {
                    degree[v] += 1;
                }
            }
            let stray = (0..mesh.n_vertices())
                .filter(|&v| {
                    degree[v] == 1 && !mesh.is_boundary_vertex(v) && !cones.contains(&(v as u32))
                })
                .count();
            assert!(
                stray <= plan.n_extension_edges,
                "{}: {stray} cut ends at vertices that are not cones",
                s.name
            );
            // A disk without cones needs no cut at all.
            if s.genus == 0 && s.boundary_loops == 1 && s.components == 1 && cones.is_empty() {
                assert!(plan.cut_edges.is_empty());
            }
        }
    }
}

/// Whether cutting along `edges` gives disks with all cones on their boundary. `None` if the
/// cut cannot be represented (an isolated one-edge slit).
fn is_valid_cut(mesh: &TriMesh, edges: &[u32], cones: &[u32]) -> Option<bool> {
    let cut = match cut_mesh(mesh, edges) {
        Ok(cut) => cut,
        Err(BffError::IsolatedCutEdge { .. }) => return None,
        Err(e) => panic!("{e}"),
    };
    let disks = component_topology(&cut.mesh).iter().all(|t| t.is_disk());
    let cones_on_cut = cut
        .vertex_to_original
        .iter()
        .enumerate()
        .all(|(v, original)| !cones.contains(original) || cut.mesh.is_boundary_vertex(v));
    Some(disks && cones_on_cut)
}

/// Verifies that a refusal of the planner states a true cause. Returns a label for counting.
fn check_refusal(
    s: &Surface,
    error: &BffError,
    cones: &[u32],
    may: &[bool],
    must: &[u32],
) -> &'static str {
    let mesh = &s.mesh;
    let name = &s.name;
    let is_must = |e: usize| must.contains(&(e as u32));
    match *error {
        BffError::ConeEnclosedByUncuttableEdges { vertex } => {
            // An interior cone none of whose edges may be cut, and no seam passes through it.
            assert!(cones.contains(&vertex), "{name}: {error}");
            assert!(!mesh.is_boundary_vertex(vertex as usize), "{name}: {error}");
            for h in mesh.outgoing_halfedges(vertex as usize) {
                let e = mesh.edge(h);
                assert!(
                    !may[e] && !is_must(e),
                    "{name}: {error}: edge {e} is usable"
                );
            }
            "cone enclosed"
        }
        BffError::SphereNeedsTwoCones {
            first_face,
            n_cones,
            ..
        } => {
            assert!(n_cones < 2, "{name}: {error}");
            // The piece of the surface that the required cuts leave around `first_face` is
            // closed, of genus 0, and has that many cones inside. (Checked when the required
            // cuts can be made as they are.)
            if let Ok(pieces) = cut_mesh(mesh, must) {
                let component = pieces.mesh.face_components()[first_face as usize];
                let t = component_topology(&pieces.mesh)[component as usize];
                assert!(t.is_closed() && t.genus == 0, "{name}: {error}: {t:?}");
                let vertex_component = pieces.mesh.vertex_components();
                let inside = (0..pieces.mesh.n_vertices())
                    .filter(|&v| {
                        vertex_component[v] == component
                            && cones.contains(&pieces.vertex_to_original[v])
                    })
                    .count();
                assert_eq!(inside, n_cones, "{name}: {error}");
            }
            "sphere needs two cones"
        }
        BffError::UncuttableRegionNotDisk {
            first_face,
            n_faces,
            euler_characteristic,
            boundary_loops,
        } => {
            // With every usable edge cut, the piece around `first_face` is what the error says.
            let everything: Vec<u32> = (0..mesh.n_edges() as u32)
                .filter(|&e| !mesh.is_boundary_edge(e as usize) && may[e as usize])
                .collect();
            if let Ok(pieces) = cut_mesh(mesh, &everything) {
                let component = pieces.mesh.face_components()[first_face as usize];
                let t = component_topology(&pieces.mesh)[component as usize];
                assert_eq!(t.n_faces, n_faces, "{name}: {error}");
                assert_eq!(
                    t.euler_characteristic, euler_characteristic,
                    "{name}: {error}"
                );
                assert_eq!(t.boundary_loops, boundary_loops, "{name}: {error}");
                assert!(!t.is_disk(), "{name}: {error}");
            }
            "uncuttable region"
        }
        BffError::IsolatedCutEdge { edge, vertices } => {
            // A one-edge slit that cannot be extended: every other edge at its two interior
            // end points is forbidden.
            assert_eq!(
                mesh.edge_vertices(edge as usize),
                vertices.map(|v| v as usize)
            );
            for v in vertices {
                assert!(!mesh.is_boundary_vertex(v as usize), "{name}: {error}");
                for h in mesh.outgoing_halfedges(v as usize) {
                    let e = mesh.edge(h);
                    assert!(
                        e == edge as usize || (!may[e] && !is_must(e)),
                        "{name}: {error}: edge {e} could extend the slit"
                    );
                }
            }
            "slit cannot be extended"
        }
        _ => panic!("{name}: undocumented refusal: {error}"),
    }
}

#[test]
fn constraints_are_honoured_and_refusals_are_true() {
    let mut outcomes: std::collections::BTreeMap<&'static str, usize> = Default::default();
    for s in zoo() {
        let mesh = &s.mesh;
        for seed in 0..40u64 {
            let mut rng = SplitMix64::new(1000 + seed);
            let forbid = [0.0, 0.1, 0.3, 0.6][(seed % 4) as usize];
            let require = [0.0, 0.02, 0.1, 0.3][((seed / 4) % 4) as usize];
            let n_cones = (seed % 6) as usize;
            let cones = pick_interior(mesh, n_cones, &mut rng);
            let may: Vec<bool> = (0..mesh.n_edges())
                .map(|_| rng.next_f64() >= forbid)
                .collect();
            let must: Vec<u32> = (0..mesh.n_edges() as u32)
                .filter(|&e| may[e as usize] && rng.next_f64() < require)
                .collect();
            let input = CutPlanInput {
                mesh,
                edge_lengths: &s.lengths,
                cone_vertices: &cones,
                may_cut: Some(&may),
                must_cut: &must,
            };
            let label = match plan_cut(&input) {
                Ok(plan) => {
                    assert_eq!(plan_cut(&input).unwrap(), plan, "{}", s.name);
                    let n_charts = check_plan(&s, &plan, &cones, Some(&may), &must);
                    if must.is_empty() {
                        assert_eq!(n_charts, s.components, "{}", s.name);
                    }
                    "planned"
                }
                Err(error) => check_refusal(&s, &error, &cones, &may, &must),
            };
            *outcomes.entry(label).or_default() += 1;
        }
    }
    println!("outcomes: {outcomes:?}");
    // The test exercised both successes and every kind of refusal.
    for label in [
        "planned",
        "cone enclosed",
        "sphere needs two cones",
        "uncuttable region",
    ] {
        assert!(
            outcomes.get(label).copied().unwrap_or(0) > 0,
            "no case of: {label}"
        );
    }
}

/// Tiny meshes on which every subset of the cuttable edges can be tried.
fn tiny_meshes() -> Vec<(&'static str, TriMesh)> {
    let tetrahedron = vec![[0u32, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
    let octahedron = vec![
        [0u32, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    let mut out = vec![
        ("tetrahedron", TriMesh::new(4, tetrahedron).unwrap()),
        ("octahedron", TriMesh::new(6, octahedron).unwrap()),
    ];
    for (name, g) in [
        ("grid 2x2", generators::grid_disk(2, 2).unwrap()),
        ("grid 3x2", generators::grid_disk(3, 2).unwrap()),
        ("grid 3x3", generators::grid_disk(3, 3).unwrap()),
        ("annulus 3x1", generators::annulus(3, 1, 1.0, 2.0).unwrap()),
        ("annulus 4x2", generators::annulus(4, 2, 1.0, 2.0).unwrap()),
        ("torus 3x3", generators::torus(3, 3, 2.0, 0.7).unwrap()),
        ("icosahedron", generators::icosphere(0).unwrap()),
    ] {
        out.push((name, g.topology().unwrap()));
    }
    out
}

#[test]
fn the_planner_finds_a_cut_exactly_when_one_exists() {
    // Most free edges for which all subsets are tried.
    const MAX_FREE: usize = 13;
    let mut compared = 0usize;
    let mut feasible_count = 0usize;
    let mut planned_count = 0usize;
    let mut by_design = 0usize;
    let mut refusals: std::collections::BTreeMap<&'static str, usize> = Default::default();
    for (name, mesh) in tiny_meshes() {
        let interior_edges: Vec<u32> = (0..mesh.n_edges() as u32)
            .filter(|&e| !mesh.is_boundary_edge(e as usize))
            .collect();
        let lengths: Vec<f64> = (0..mesh.n_edges())
            .map(|e| 1.0 + 0.37 * ((e * 7919) % 13) as f64)
            .collect();
        for seed in 0..160u64 {
            let mut rng = SplitMix64::new(77 + seed);
            // Forbid enough edges that at most MAX_FREE remain free.
            let base = [0.0, 0.25, 0.5, 0.75][(seed % 4) as usize];
            let mut may: Vec<bool> = (0..mesh.n_edges())
                .map(|_| rng.next_f64() >= base)
                .collect();
            let require = [0.0, 0.1, 0.3][((seed / 4) % 3) as usize];
            let must: Vec<u32> = interior_edges
                .iter()
                .copied()
                .filter(|&e| may[e as usize] && rng.next_f64() < require)
                .collect();
            let mut free: Vec<u32> = interior_edges
                .iter()
                .copied()
                .filter(|&e| may[e as usize] && !must.contains(&e))
                .collect();
            while free.len() > MAX_FREE {
                let drop = free.swap_remove((rng.next_u64() % free.len() as u64) as usize);
                may[drop as usize] = false;
            }
            let n_cones = (seed % 4) as usize;
            let cones = pick_interior(&mesh, n_cones, &mut rng);

            let planned = plan_cut(&CutPlanInput {
                mesh: &mesh,
                edge_lengths: &lengths,
                cone_vertices: &cones,
                may_cut: Some(&may),
                must_cut: &must,
            });
            // The truth: try every subset of the free edges.
            let mut feasible = false;
            let mut edges = must.clone();
            for subset in 0u32..(1 << free.len()) {
                edges.truncate(must.len());
                edges.extend(
                    free.iter()
                        .enumerate()
                        .filter(|&(k, _)| subset >> k & 1 == 1)
                        .map(|(_, &e)| e),
                );
                if is_valid_cut(&mesh, &edges, &cones) == Some(true) {
                    feasible = true;
                    break;
                }
            }
            match &planned {
                Ok(plan) => {
                    planned_count += 1;
                    assert!(feasible, "{name} seed {seed}: planned an impossible cut");
                    assert_eq!(
                        is_valid_cut(&mesh, &plan.cut_edges, &cones),
                        Some(true),
                        "{name} seed {seed}"
                    );
                    assert!(plan.cut_edges.iter().all(|&e| may[e as usize]));
                    assert!(must.iter().all(|e| plan.cut_edges.contains(e)));
                }
                // By design the planner does not invent the end points of a cut on a closed
                // surface of genus 0: it asks for two cones even though a cut exists.
                Err(BffError::SphereNeedsTwoCones { .. }) => by_design += 1,
                Err(error) => {
                    assert!(
                        !feasible,
                        "{name} seed {seed}: refused with \"{error}\" although a cut exists \
                         (cones {cones:?}, must {must:?}, forbidden {:?})",
                        (0..mesh.n_edges()).filter(|&e| !may[e]).collect::<Vec<_>>()
                    );
                    let kind = match error {
                        BffError::ConeEnclosedByUncuttableEdges { .. } => "cone enclosed",
                        BffError::UncuttableRegionNotDisk { .. } => "uncuttable region",
                        BffError::IsolatedCutEdge { .. } => "slit cannot be extended",
                        _ => panic!("{name} seed {seed}: undocumented refusal: {error}"),
                    };
                    *refusals.entry(kind).or_default() += 1;
                }
            }
            compared += 1;
            feasible_count += usize::from(feasible);
        }
    }
    println!(
        "{compared} instances compared with exhaustive search: {feasible_count} have a cut; \
         {planned_count} planned, {by_design} refused by design (closed genus 0 with fewer than \
         two cones), the rest refused although no cut exists: {refusals:?}"
    );
    // Every instance with a cut was planned, except the refusals by design.
    assert!(planned_count + by_design >= feasible_count);
    assert!(planned_count > compared / 5);
    assert!(refusals.values().sum::<usize>() > compared / 20);
}

#[test]
fn malformed_planner_input_is_rejected_with_a_cause() {
    let g = generators::icosphere(1).unwrap();
    let mesh = g.topology().unwrap();
    let lengths = vec![1.0; mesh.n_edges()];
    let plan = |lengths: &[f64], cones: &[u32], may: Option<&[bool]>, must: &[u32]| {
        plan_cut(&CutPlanInput {
            mesh: &mesh,
            edge_lengths: lengths,
            cone_vertices: cones,
            may_cut: may,
            must_cut: must,
        })
    };
    let err = plan(&lengths[1..], &[0, 1], None, &[]).unwrap_err();
    assert!(matches!(err, BffError::LengthMismatch { .. }), "{err}");
    let mut bad = lengths.clone();
    bad[3] = f64::NAN;
    let err = plan(&bad, &[0, 1], None, &[]).unwrap_err();
    assert!(
        matches!(err, BffError::NonFiniteInput { index: 3, .. }),
        "{err}"
    );
    bad[3] = -1.0;
    let err = plan(&bad, &[0, 1], None, &[]).unwrap_err();
    assert!(err.to_string().contains("non-negative"), "{err}");
    let err = plan(&lengths, &[0, 4000], None, &[]).unwrap_err();
    assert!(
        matches!(err, BffError::ConeVertexOutOfRange { vertex: 4000, .. }),
        "{err}"
    );
    let err = plan(&lengths, &[0, 1], Some(&[true; 5]), &[]).unwrap_err();
    assert!(matches!(err, BffError::LengthMismatch { .. }), "{err}");
    let err = plan(&lengths, &[0, 1], None, &[9999]).unwrap_err();
    assert!(
        matches!(err, BffError::EdgeOutOfRange { edge: 9999, .. }),
        "{err}"
    );
    let mut may = vec![true; mesh.n_edges()];
    may[5] = false;
    let err = plan(&lengths, &[0, 1], Some(&may), &[5]).unwrap_err();
    assert!(
        matches!(err, BffError::RequiredEdgeForbidden { edge: 5, .. }),
        "{err}"
    );
    // Zero lengths are allowed (they are costs), and repeated cones count once.
    let zero = vec![0.0; mesh.n_edges()];
    let plan_zero = plan(&zero, &[0, 1, 1, 0], None, &[]).unwrap();
    assert_eq!(
        is_valid_cut(&mesh, &plan_zero.cut_edges, &[0, 1]),
        Some(true)
    );
}
