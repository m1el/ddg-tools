//! Property tests: random meshes of every topology, polygonal input with a supplied
//! triangulation, random seams, steering and options.

mod common;

use common::{check, edge_uses, without_times};
use ddg_core::generators::SplitMix64;
use ddg_uv::primitives::{
    End, cone, cylinder, grid_cube, holed_slab, jittered, revolution, torus, uv_sphere,
};
use ddg_uv::{
    BoundaryValues, ClosedSurfaces, Cone, ConeSigns, MeshInput, MergeClusters, PackOptions, PolygonMesh,
    SlitPlacement, UnwrapOptions, UvError, unwrap,
};
use proptest::prelude::*;

/// A bumpy grid of `nx × ny` quads.
fn quad_grid(nx: u32, ny: u32) -> PolygonMesh {
    let mut positions = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let (x, y) = (f64::from(i) / f64::from(nx), f64::from(j) / f64::from(ny));
            positions.push([x, y, 0.2 * (4.0 * x).sin() * (3.0 * y).cos()]);
        }
    }
    let at = |i: u32, j: u32| j * (nx + 1) + i;
    let mut polygons = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            polygons.push(vec![at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    PolygonMesh::from_polygons(positions, &polygons)
}

/// Two meshes as one, the second moved away.
fn union(a: &PolygonMesh, b: &PolygonMesh) -> PolygonMesh {
    let offset = a.positions.len() as u32;
    let mut polygons: Vec<Vec<u32>> = (0..a.n_polygons()).map(|p| a.polygon(p).to_vec()).collect();
    polygons.extend((0..b.n_polygons()).map(|p| b.polygon(p).iter().map(|v| v + offset).collect()));
    let mut positions = a.positions.clone();
    positions.extend(b.positions.iter().map(|p| [p[0] + 10.0, p[1], p[2]]));
    PolygonMesh::from_polygons(positions, &polygons)
}

/// The mesh of a case: its kind decides the topology.
fn mesh_of(kind: usize, rng: &mut SplitMix64) -> PolygonMesh {
    let base = match kind {
        0 => quad_grid(5, 4),
        1 => uv_sphere(8, 5, 1.0).unwrap(),
        2 => grid_cube(3, 2.0),
        3 => cylinder(8, 2, 1.0, 2.0).unwrap(),
        4 => cone(8, 1.0, 2.0).unwrap(),
        5 => torus(8, 6, 1.0, 0.4).unwrap(),
        6 => holed_slab(2, 1, 1.0),
        7 => revolution(&[(1.0, 1.0), (1.0, 0.3), (0.8, -0.2)], 7, End::Open, End::Cap).unwrap(),
        8 => revolution(&[(0.5, 0.3), (1.0, 0.2), (1.5, 0.0)], 9, End::Open, End::Open).unwrap(),
        _ => union(&grid_cube(2, 1.0), &quad_grid(3, 3)),
    };
    let mut mesh = jittered(&base, 0.02, rng.next_u64());
    // A supplied triangulation that is not the fan: some quads are split along their other
    // diagonal.
    let mut t = 0;
    for p in 0..mesh.n_polygons() {
        let corners = mesh.polygon_corners(p);
        if corners.len() == 4 && rng.next_u64().is_multiple_of(2) {
            let c = corners.start as u32;
            mesh.triangles[t] = [c, c + 1, c + 3];
            mesh.triangles[t + 1] = [c + 1, c + 2, c + 3];
        }
        t += corners.len() - 2;
    }
    mesh
}

fn options_of(rng: &mut SplitMix64) -> UnwrapOptions {
    let pick = |rng: &mut SplitMix64, n: u64| (rng.next_u64() % n) as usize;
    UnwrapOptions {
        lambda: [0.5, 0.1, 0.02][pick(rng, 3)],
        cone_signs: [
            ConeSigns::PositiveWhereFeasible,
            ConeSigns::Any,
            ConeSigns::PositiveOnly,
        ][pick(rng, 3)],
        closed_surfaces: [
            UnwrapOptions::default().closed_surfaces,
            ClosedSurfaces::RawOptimum,
            ClosedSurfaces::Slit {
                placement: SlitPlacement::DiameterMiddle,
                relative_length: 0.5,
            },
            ClosedSurfaces::Slit {
                placement: SlitPlacement::AtLargestCurvature,
                relative_length: 0.0,
            },
            ClosedSurfaces::FreeSites {
                count: 6,
                max_angle: 3.5,
                forbid_elsewhere: true,
            },
        ][pick(rng, 5)],
        merge: [
            UnwrapOptions::default().merge,
            MergeClusters::Never,
            MergeClusters::Always {
                max_cluster_size: 8,
            },
        ][pick(rng, 3)],
        packing: [
            Some(PackOptions::default()),
            Some(PackOptions {
                margin: 0.02,
                rotate: false,
            }),
            None,
        ][pick(rng, 3)],
        reoptimize_angles: pick(rng, 3) == 0,
        boundary: [BoundaryValues::Zero, BoundaryValues::Optimized][pick(rng, 2)],
        ..UnwrapOptions::default()
    }
}

/// Everything a case is made of.
struct Case {
    mesh: PolygonMesh,
    seams: Vec<[u32; 2]>,
    pinned: Vec<u32>,
    forbidden: Vec<u32>,
    prescribed: Vec<Cone>,
    options: UnwrapOptions,
}

fn case_of(kind: usize, seed: u64) -> Case {
    let mut rng = SplitMix64::new(seed);
    let mesh = mesh_of(kind, &mut rng);
    let edges = mesh.edges();
    let n_vertices = mesh.positions.len() as u32;
    // Seams: a few random polygon edges, in either orientation, one of them twice, and
    // sometimes a pair of vertices that is no edge.
    let mut seams: Vec<[u32; 2]> = Vec::new();
    for _ in 0..rng.next_u64() % 7 {
        let e = edges[(rng.next_u64() % edges.len() as u64) as usize];
        seams.push(if rng.next_u64().is_multiple_of(2) {
            e
        } else {
            [e[1], e[0]]
        });
    }
    if let Some(&first) = seams.first() {
        seams.push(first);
    }
    if rng.next_u64().is_multiple_of(4) {
        seams.push([0, n_vertices - 1]);
    }
    // Steering: disjoint random vertices.
    let mut order: Vec<u32> = (0..n_vertices).collect();
    for i in (1..order.len()).rev() {
        order.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
    }
    let (n_pinned, n_forbidden, n_prescribed) = (
        (rng.next_u64() % 3) as usize,
        (rng.next_u64() % 4) as usize,
        (rng.next_u64() % 3) as usize,
    );
    let pinned = order[..n_pinned].to_vec();
    let forbidden = order[n_pinned..n_pinned + n_forbidden].to_vec();
    let prescribed: Vec<Cone> = order[n_pinned + n_forbidden..]
        .iter()
        .take(n_prescribed)
        .map(|&vertex| Cone {
            vertex,
            angle: rng.uniform(0.05, 0.6),
        })
        .collect();
    let options = options_of(&mut rng);
    Case {
        mesh,
        seams,
        pinned,
        forbidden,
        prescribed,
        options,
    }
}

/// How a case ended, for the census.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Outcome {
    Unwrapped,
    NegativeConesRequired,
    InfeasibleSteering,
}

/// Runs one case and checks everything; `Err` is a failed property.
fn run_case(kind: usize, seed: u64) -> Result<(Outcome, usize, usize, usize), String> {
    let Case {
        mesh,
        seams,
        pinned,
        forbidden,
        prescribed,
        options,
    } = case_of(kind, seed);
    let untouched = mesh.clone();
    let (seams_before, pinned_before) = (seams.clone(), pinned.clone());
    let input = MeshInput {
        seam_edges: &seams,
        pinned_vertices: &pinned,
        forbidden_vertices: &forbidden,
        prescribed_cones: &prescribed,
        ..mesh.input()
    };
    let result = unwrap(&input, &options);
    // The input is not modified.
    if mesh != untouched || seams != seams_before || pinned != pinned_before {
        return Err("the input was modified".into());
    }
    let map = match result {
        Ok(map) => map,
        // Refusals that these random inputs can legitimately run into.
        Err(UvError::NegativeConesRequired { .. })
            if options.cone_signs == ConeSigns::PositiveOnly && kind == 6 =>
        {
            return Ok((Outcome::NegativeConesRequired, 0, 0, 0));
        }
        Err(UvError::InfeasibleSteering { .. })
            if !prescribed.is_empty() || !forbidden.is_empty() =>
        {
            return Ok((Outcome::InfeasibleSteering, 0, 0, 0));
        }
        Err(e) => return Err(format!("{e}")),
    };
    check(&mesh, &input, &map, &options);

    // Every seam that is an interior polygon edge is cut.
    let uses = edge_uses(&mesh);
    let mut n_seams_cut = 0;
    for seam in &seams {
        let key = [seam[0].min(seam[1]), seam[0].max(seam[1])];
        if uses.get(&key) == Some(&2) {
            if map.cut_edges.binary_search(&key).is_err() {
                return Err(format!("the seam {key:?} was not cut"));
            }
            n_seams_cut += 1;
        }
    }
    // Prescribed cones are imposed exactly, or reported as not imposable.
    for cone in &prescribed {
        let found = map.cones.iter().find(|c| c.vertex == cone.vertex);
        let skipped = map.report.steering.prescribed_not_applied.contains(&cone.vertex);
        if skipped != found.is_none() || (!skipped && found.map(|c| c.angle) != Some(cone.angle)) {
            return Err(format!("prescribed {cone:?} became {found:?} (skipped: {skipped})"));
        }
    }
    let steering = &map.report.steering;
    if steering.n_prescribed_applied + steering.prescribed_not_applied.len() != prescribed.len()
        || steering.n_pinned_applied + steering.n_pinned_without_effect != pinned.len()
        || steering.n_forbidden_applied + steering.n_forbidden_without_effect != forbidden.len()
    {
        return Err(format!("the steering counts do not add up: {steering:?}"));
    }
    // Positive cones where they were asked for (prescribed ones are the caller's).
    if options.cone_signs != ConeSigns::Any && kind != 6 {
        for cone in &map.cones {
            let own = prescribed.iter().any(|c| c.vertex == cone.vertex);
            if !(own || cone.angle > 0.0) {
                return Err(format!("a negative cone {cone:?} where none was allowed"));
            }
        }
    }
    // The same again.
    let n = (map.cones.len(), map.report.charts.len());
    let again = unwrap(&input, &options).map_err(|e| format!("second run: {e}"))?;
    if without_times(map) != without_times(again) {
        return Err("the second run differs from the first".into());
    }
    Ok((Outcome::Unwrapped, n_seams_cut, n.0, n.1))
}

proptest! {
    // A failure prints its kind and seed, which reproduce it; no regression file is kept.
    #![proptest_config(ProptestConfig {
        cases: 64,
        failure_persistence: None,
        ..ProptestConfig::default()
    })]

    #[test]
    fn random_meshes_seams_and_steering(kind in 0usize..10, seed in any::<u64>()) {
        if let Err(message) = run_case(kind, seed) {
            return Err(TestCaseError::fail(format!("kind {kind} seed {seed}: {message}")));
        }
    }
}

/// What the random cases exercise: `cargo test -p ddg-uv --release --test properties census
/// -- --ignored --nocapture` prints how 2000 fixed cases end.
#[test]
#[ignore = "a census of the property test's cases, run on request"]
fn census() {
    let mut outcomes = std::collections::BTreeMap::new();
    let (mut with_seams, mut with_cones, mut several_charts, mut steered) = (0, 0, 0, 0);
    let mut seeds = SplitMix64::new(2026);
    let total = 2000;
    for i in 0..total {
        let (kind, seed) = (i % 10, seeds.next_u64());
        let (outcome, n_seams_cut, n_cones, n_charts) =
            run_case(kind, seed).unwrap_or_else(|e| panic!("kind {kind} seed {seed}: {e}"));
        *outcomes.entry(outcome).or_insert(0usize) += 1;
        with_seams += usize::from(n_seams_cut > 0);
        with_cones += usize::from(n_cones > 0);
        several_charts += usize::from(n_charts > 1);
        let case = case_of(kind, seed);
        steered += usize::from(
            !(case.pinned.is_empty() && case.forbidden.is_empty() && case.prescribed.is_empty()),
        );
    }
    println!(
        "{total} cases: {outcomes:?}; unwrapped with a seam cut {with_seams}, with cones \
         {with_cones}, with several charts {several_charts}; cases with steering {steered}"
    );
}
