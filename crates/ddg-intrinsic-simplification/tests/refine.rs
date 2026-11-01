//! Pre- and post-refinement (SPEC §6; WP5): the coarse mesh of a refined triangulation, the
//! refinement of a coarse mesh, and the tracked map through both.

mod common;

use std::f64::consts::PI;

use common::{Shape, crane_models, dist, planar_grid, reference_meshes, small_shapes};
use ddg_core::geometry::triangle_area;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::RefineRules;
use ddg_intrinsic_simplification::refine::{RefineError, RefinementOptions, Representation};
use ddg_intrinsic_simplification::{CoarseMesh, CoarseMeshOptions, CoarseOutput, CoarsePoint};

fn tri_mesh(shape: &Shape) -> (TriMesh, Geometry) {
    let mesh = TriMesh::new(shape.positions.len(), shape.faces.clone()).unwrap();
    let geometry = Geometry::from_positions(&mesh, &shape.positions).unwrap();
    (mesh, geometry)
}

fn options(representation: Representation) -> RefinementOptions {
    RefinementOptions {
        representation,
        ..RefinementOptions::default()
    }
}

/// Removes vertices in passes by ascending fine id (inserted vertices included) until `target`
/// vertices are left or nothing more can be removed; certifies every `every` removals.
fn remove_in_passes(mesh: &mut CoarseMesh, target: usize, every: usize, name: &str) -> usize {
    let mut n = 0usize;
    loop {
        let mut progress = false;
        let mut order = mesh.vertices().to_vec();
        order.sort_unstable();
        for v in order {
            if mesh.n_vertices() <= target {
                return n;
            }
            if mesh.remove_vertex(v).is_ok() {
                progress = true;
                n += 1;
                if n.is_multiple_of(every) {
                    mesh.certify(false)
                        .unwrap_or_else(|e| panic!("{name}: after {n} removals: {e}"));
                }
            }
        }
        if !progress {
            return n;
        }
    }
}

/// The planar positions of the vertices of a flat coarse mesh (a disk with flat interior),
/// laid out face by face from the lengths, then moved rigidly so that the vertices `a` and `b`
/// (compacted ids) sit at `pa` and `pb`.
fn layout(out: &CoarseOutput, a: u32, b: u32, pa: [f64; 2], pb: [f64; 2]) -> Vec<Option<[f64; 2]>> {
    let mut pos: Vec<Option<[f64; 2]>> = vec![None; out.n_vertices()];
    let mut seen = vec![false; out.n_faces()];
    let l = &out.face_lengths[0];
    pos[out.faces[0][0] as usize] = Some([0.0, 0.0]);
    pos[out.faces[0][1] as usize] = Some([l[0], 0.0]);
    let x = (l[0] * l[0] + l[2] * l[2] - l[1] * l[1]) / (2.0 * l[0]);
    pos[out.faces[0][2] as usize] =
        Some([x, 2.0 * triangle_area(l[0], l[1], l[2]).unwrap() / l[0]]);
    seen[0] = true;
    let mut queue = std::collections::VecDeque::from([0usize]);
    while let Some(f) = queue.pop_front() {
        for k in 0..3 {
            let Some(t) = out.twin[3 * f + k] else {
                continue;
            };
            let g = t as usize / 3;
            if seen[g] {
                continue;
            }
            seen[g] = true;
            // The twin goes from slot k' to k'+1 of g: from the head to the tail of ours.
            let kk = t as usize % 3;
            let (vb, va, vw) = (
                out.faces[g][kk] as usize,
                out.faces[g][(kk + 1) % 3] as usize,
                out.faces[g][(kk + 2) % 3] as usize,
            );
            let (pb_, pa_) = (pos[vb].unwrap(), pos[va].unwrap());
            let lg = &out.face_lengths[g];
            let (len, m, n) = (lg[kk], lg[(kk + 1) % 3], lg[(kk + 2) % 3]);
            let x = (len * len + n * n - m * m) / (2.0 * len);
            let y = 2.0 * triangle_area(lg[0], lg[1], lg[2]).unwrap() / len;
            let u = [(pa_[0] - pb_[0]) / len, (pa_[1] - pb_[1]) / len];
            let w = [pb_[0] + x * u[0] - y * u[1], pb_[1] + x * u[1] + y * u[0]];
            if pos[vw].is_none() {
                pos[vw] = Some(w);
            }
            queue.push_back(g);
        }
    }
    // Rigid motion: rotate the direction a -> b onto pa -> pb.
    let (qa, qb) = (pos[a as usize].unwrap(), pos[b as usize].unwrap());
    let ang = (pb[1] - pa[1]).atan2(pb[0] - pa[0]) - (qb[1] - qa[1]).atan2(qb[0] - qa[0]);
    let (s, c) = ang.sin_cos();
    pos.into_iter()
        .map(|p| {
            p.map(|p| {
                let d = [p[0] - qa[0], p[1] - qa[1]];
                [pa[0] + c * d[0] - s * d[1], pa[1] + s * d[0] + c * d[1]]
            })
        })
        .collect()
}

/// The planar position of every tracked point of a flat coarse mesh, by the layout.
fn point_positions(shape: &Shape, mesh: &CoarseMesh, corners: [u32; 2]) -> Vec<[f64; 2]> {
    let out = mesh.compact().unwrap();
    let ca = out.coarse_of_fine[corners[0] as usize].unwrap();
    let cb = out.coarse_of_fine[corners[1] as usize].unwrap();
    let p = |v: u32| {
        let q = shape.positions[v as usize];
        [q[0], q[1]]
    };
    let pos = layout(&out, ca, cb, p(corners[0]), p(corners[1]));
    out.points
        .iter()
        .map(|pt| match *pt {
            CoarsePoint::Vertex(v) => pos[v as usize].unwrap(),
            CoarsePoint::Face { face, bary } => {
                let mut r = [0.0; 2];
                for k in 0..3 {
                    let q = pos[out.faces[face as usize][k] as usize].unwrap();
                    r[0] += bary[k] * q[0];
                    r[1] += bary[k] * q[1];
                }
                r
            }
        })
        .collect()
}

fn check_prolongation(mesh: &CoarseMesh, name: &str) {
    let out = mesh.compact().unwrap();
    let p = out.scalar_prolongation().unwrap();
    let mut y = vec![0.0; p.nrows()];
    p.mul_vec(&vec![1.0; p.ncols()], &mut y).unwrap();
    let worst = y.iter().map(|s| (s - 1.0).abs()).fold(0.0, f64::max);
    assert!(worst < 1e-13, "{name}: rows of P sum to 1 within {worst:e}");
    let q = out.point_prolongation().unwrap();
    let mut y = vec![0.0; q.nrows()];
    q.mul_vec(&vec![1.0; q.ncols()], &mut y).unwrap();
    let worst = y.iter().map(|s| (s - 1.0).abs()).fold(0.0, f64::max);
    assert!(
        worst < 1e-13,
        "{name}: point rows sum to 1 within {worst:e}"
    );
}

/// Pre-refinement, removal of everything but the four corners, post-refinement on a planar
/// grid with extra points: the metric stays planar throughout, so every tracked point is
/// reproduced (SPEC §12.2).
#[test]
fn planar_grid_through_pre_refinement_removals_and_post_refinement() {
    for representation in [Representation::Integer, Representation::Signposts] {
        for (n, jitter, seed) in [(6usize, 0.0, 1u64), (10, 0.3, 2), (16, 0.5, 3)] {
            let shape = planar_grid(n, jitter, seed);
            let name = format!("{} {representation:?}", shape.name);
            let points = shape.random_points(seed);
            let (mesh, geometry) = tri_mesh(&shape);
            let (mut coarse, pre) = CoarseMesh::from_refined(
                &mesh,
                &geometry,
                &points,
                CoarseMeshOptions::default(),
                &options(representation),
            )
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            coarse.certify(true).unwrap();
            let n_fine = shape.positions.len();
            assert_eq!(coarse.n_fine_vertices(), n_fine);
            assert_eq!(coarse.n_points(), n_fine + points.len());
            assert_eq!(coarse.n_vertices(), n_fine + pre.n_added);
            // Inserted vertices: ids from |V|, free, not tracked, removable.
            for v in n_fine as u32..(n_fine + pre.n_added) as u32 {
                assert!(coarse.is_vertex(v) && !coarse.is_fixed(v));
            }
            check_prolongation(&coarse, &name);
            // Planar already after pre-refinement.
            let corners = [
                0u32,
                n as u32,
                (n * (n + 1)) as u32,
                ((n + 1) * (n + 1) - 1) as u32,
            ];
            // The corners and a scattering of other input vertices stay (a mesh with thin
            // triangles for the post-refinement to repair).
            let keep = |v: u32| {
                corners.contains(&v)
                    || (v as u64).wrapping_mul(2654435761).wrapping_add(seed) % 7 == 0
            };
            for v in 0..n_fine as u32 {
                coarse.set_fixed(v, keep(v)).unwrap();
            }
            let n_kept = (0..n_fine as u32).filter(|&v| keep(v)).count();
            let truth = |k: usize| -> [f64; 2] {
                if k < n_fine {
                    [shape.positions[k][0], shape.positions[k][1]]
                } else {
                    let q = shape.input_position(&points[k - n_fine]);
                    [q[0], q[1]]
                }
            };
            let worst_of = |mesh: &CoarseMesh| -> f64 {
                let got = point_positions(&shape, mesh, [corners[0], corners[1]]);
                (0..got.len())
                    .map(|k| dist([got[k][0], got[k][1], 0.0], [truth(k)[0], truth(k)[1], 0.0]))
                    .fold(0.0, f64::max)
            };
            let w0 = worst_of(&coarse);
            assert!(
                w0 < 1e-9,
                "{name}: after pre-refinement, worst point error {w0:e}"
            );
            // The edges of T¹ between two input vertices are the straight segments.
            for e in 0..coarse.n_edges() {
                let [a, b] = coarse.complex().edge_vertices(e);
                let (fa, fb) = (
                    coarse.fine_vertex(a) as usize,
                    coarse.fine_vertex(b) as usize,
                );
                if fa < n_fine && fb < n_fine {
                    let d = dist(shape.positions[fa], shape.positions[fb]);
                    assert!(
                        (d - coarse.edge_lengths()[e]).abs() < 1e-12,
                        "{name}: edge {e}"
                    );
                }
            }
            let removed = remove_in_passes(&mut coarse, n_kept, 25, &name);
            assert_eq!(
                coarse.n_vertices(),
                n_kept,
                "{name}: only the kept vertices remain"
            );
            coarse.certify(true).unwrap();
            let w1 = worst_of(&coarse);
            assert!(
                w1 < 1e-9,
                "{name}: after {removed} removals, worst point error {w1:e}"
            );
            // Post-refinement of the 4-vertex square.
            let before = coarse.clone();
            let (refined, post) = match coarse.post_refine(&options(representation)) {
                Ok(x) => x,
                Err(RefineError::Locate { point, what })
                    if representation == Representation::Signposts =>
                {
                    // The signposts' point query cannot reach a point on an input edge past a
                    // non-flat vertex (reported with the point, never placed wrongly).
                    eprintln!("{name}: signposts could not locate point {point}: {what}");
                    continue;
                }
                Err(e) => panic!("{name}: {e}"),
            };
            assert!(coarse.bitwise_eq(&before), "post_refine changed its input");
            refined.certify(true).unwrap();
            assert_eq!(refined.n_points(), coarse.n_points());
            assert_eq!(refined.n_fine_vertices(), n_fine);
            assert_eq!(refined.n_vertices(), n_kept + post.n_added);
            for &c in &corners {
                assert!(refined.is_vertex(c) && refined.is_fixed(c));
            }
            // New ids are above everything used so far.
            let used = n_fine + pre.n_added;
            for &v in refined.vertices() {
                assert!(
                    (v as usize) >= used || keep(v),
                    "{name}: vertex {v}"
                );
            }
            let w2 = worst_of(&refined);
            assert!(
                w2 < 1e-9,
                "{name}: after post-refinement, worst point error {w2:e}"
            );
            // Quality: a flat disk, every input angle sum is 90° (< 60° never), not closed.
            assert!(!post.closed && !post.guarantee_applies);
            eprintln!(
                "{name}: +{} (pre) -> {removed} removed -> +{} (post), min angle {:.2} deg, \
                 errors {w0:e} {w1:e} {w2:e}",
                pre.n_added,
                post.n_added,
                post.min_angle.to_degrees()
            );
            assert!(
                post.min_angle >= 24.9f64.to_radians()
                    || post.refine.min_angle_reachable.is_infinite(),
                "{name}: {}",
                post.min_angle
            );
            check_prolongation(&refined, &name);
        }
    }
}

/// Models: pre-refinement of the input, every extra point located, inserted vertices are
/// ordinary removable vertices; coarsening; post-refinement where the coarse mesh is simplicial.
#[test]
fn models_pre_and_post_refinement() {
    let mut shapes = small_shapes();
    shapes.extend(
        crane_models()
            .into_iter()
            .filter(|s| s.positions.len() < 8_000),
    );
    shapes.extend(reference_meshes(8_000));
    let mut n_post = 0usize;
    let mut n_refused = 0usize;
    let mut n_unlocated = 0usize;
    for shape in &shapes {
        for representation in [Representation::Integer, Representation::Signposts] {
            let name = format!("{} {representation:?}", shape.name);
            let points = shape.random_points(5);
            let (mesh, geometry) = tri_mesh(shape);
            let (mut coarse, pre) = CoarseMesh::from_refined(
                &mesh,
                &geometry,
                &points,
                CoarseMeshOptions::default(),
                &options(representation),
            )
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            let n_fine = shape.positions.len();
            coarse
                .certify(true)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(coarse.n_points(), n_fine + points.len());
            assert_eq!(coarse.n_vertices(), n_fine + pre.n_added);
            assert!(pre.certificate.n_not_delaunay == 0, "{name}");
            check_prolongation(&coarse, &name);
            // The pre-refined mesh reaches the bound (or says why not).
            let reachable = pre.refine.min_angle_reachable;
            assert!(
                reachable >= 24.9f64.to_radians() || pre.refine.min_angle_reachable.is_infinite(),
                "{name}: smallest reachable angle {:.3} deg",
                reachable.to_degrees()
            );
            // Coarsen to a quarter of the vertices in a fixed order (inserted ones included).
            let target = (coarse.n_vertices() / 4).max(8);
            let removed = remove_in_passes(&mut coarse, target, 40, &name);
            coarse
                .certify(false)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            check_prolongation(&coarse, &name);
            let kept_inserted = coarse
                .vertices()
                .iter()
                .filter(|&&v| v as usize >= n_fine)
                .count();
            if kept_inserted > 0 {
                // Inserted vertices that stay are columns of P.
                let out = coarse.compact().unwrap();
                let q = out.point_prolongation().unwrap();
                assert_eq!(q.ncols(), out.n_vertices());
                let mut colsum = vec![0.0; q.ncols()];
                for (_, c, v) in q.triplets() {
                    colsum[c as usize] += v;
                }
                assert!(
                    out.vertices
                        .iter()
                        .enumerate()
                        .any(|(i, &f)| f as usize >= n_fine && colsum[i] > 0.0),
                    "{name}: no kept inserted vertex has a column"
                );
            }
            match coarse.post_refine(&options(representation)) {
                Ok((refined, post)) => {
                    n_post += 1;
                    refined
                        .certify(true)
                        .unwrap_or_else(|e| panic!("{name}: post: {e}"));
                    assert_eq!(refined.n_points(), coarse.n_points());
                    check_prolongation(&refined, &name);
                    // Vertices of the coarse mesh stay.
                    for &v in coarse.vertices() {
                        assert!(refined.is_vertex(v), "{name}: vertex {v} lost");
                    }
                    assert!(post.certificate.n_not_delaunay == 0);
                    let guarantee = post.guarantee_applies;
                    if guarantee {
                        assert!(post.min_angle >= 24.9f64.to_radians(), "{name}: {post:?}");
                    }
                    eprintln!(
                        "{name}: pre +{} ({} kept after {removed} removals), post +{} min angle \
                         {:.2} deg (guarantee {guarantee}), located {}",
                        pre.n_added,
                        kept_inserted,
                        post.n_added,
                        post.min_angle.to_degrees(),
                        post.n_located
                    );
                }
                Err(RefineError::NotSimplicial { what }) => {
                    n_refused += 1;
                    eprintln!("{name}: post-refinement refused: {what}");
                }
                Err(RefineError::Locate { point, what })
                    if representation == Representation::Signposts =>
                {
                    n_unlocated += 1;
                    eprintln!("{name}: signposts could not locate point {point}: {what}");
                }
                Err(e) => panic!("{name}: {e}"),
            }
        }
    }
    eprintln!(
        "post-refinement: {n_post} done, {n_refused} refused (not simplicial), {n_unlocated} signpost \
         queries failed"
    );
    assert!(n_post > 0);
}

/// The Δ-complex refusal names what is wrong.
#[test]
fn delta_complex_is_refused_with_a_reason() {
    // A closed tetrahedron coarsened to the minimum has a bigon or a self-edge at some point.
    let mut refused = 0;
    for shape in small_shapes() {
        let mut coarse = shape.coarse(&[]);
        coarse.flip_to_delaunay().unwrap();
        remove_in_passes(&mut coarse, 3, 1000, &shape.name);
        match coarse.post_refine(&RefinementOptions::default()) {
            Err(RefineError::NotSimplicial { what }) => {
                assert!(!what.is_empty());
                refused += 1;
                eprintln!("{}: {what}", shape.name);
            }
            Ok(_) => {}
            Err(e) => panic!("{}: {e}", shape.name),
        }
    }
    assert!(refused > 0, "no fully coarsened shape was a Δ-complex");
}

/// Both rules and both representations on a closed model, and the guarantee flag.
#[test]
fn closed_surface_reaches_the_bound() {
    let shape = small_shapes()
        .into_iter()
        .find(|s| s.name == "torus 12x6")
        .unwrap();
    let (mesh, geometry) = tri_mesh(&shape);
    for representation in [Representation::Integer, Representation::Signposts] {
        for rules in [RefineRules::Chew, RefineRules::Ruppert] {
            let opts = RefinementOptions {
                representation,
                rules: Some(rules),
                ..RefinementOptions::default()
            };
            let (coarse, rep) = CoarseMesh::from_refined(
                &mesh,
                &geometry,
                &[],
                CoarseMeshOptions::default(),
                &opts,
            )
            .unwrap();
            coarse.certify(true).unwrap();
            assert!(rep.closed);
            assert_eq!(rep.guarantee_applies, rep.min_input_angle_sum >= PI / 3.0);
            if rep.guarantee_applies {
                assert!(rep.min_angle >= 24.9f64.to_radians(), "{rep:?}");
            }
        }
    }
}

/// A coarse Δ-complex that is refused becomes acceptable after `make_simplicial` flips (the
/// metric and the tracked map are kept), and the caller's mesh is not changed by the option.
#[test]
fn flips_make_a_refused_mesh_acceptable() {
    use ddg_core::generators::SplitMix64;
    let Some(shape) = common::load_data("intrinsic-simplification/boxpart.obj") else {
        return;
    };
    let (mesh, geometry) = tri_mesh(&shape);
    let points = shape.random_points(9);
    let (pre, _) = CoarseMesh::from_refined(
        &mesh,
        &geometry,
        &points,
        CoarseMeshOptions::default(),
        &RefinementOptions::default(),
    )
    .unwrap();
    let mut repaired = 0;
    for fraction in [0.5, 0.25] {
        let mut coarse = pre.clone();
        let mut order = coarse.vertices().to_vec();
        let mut rng = SplitMix64::new(42);
        for i in (1..order.len()).rev() {
            order.swap(i, (rng.next_u64() % (i as u64 + 1)) as usize);
        }
        let target = (pre.n_vertices() as f64 * fraction) as usize;
        for v in order {
            if coarse.n_vertices() <= target {
                break;
            }
            let _ = coarse.remove_vertex(v);
        }
        let plain = coarse.post_refine(&RefinementOptions::default());
        let Err(RefineError::NotSimplicial { .. }) = plain else {
            continue;
        };
        let before = coarse.clone();
        let opts = RefinementOptions {
            make_simplicial_flips: Some(100_000),
            ..RefinementOptions::default()
        };
        let (refined, report) = coarse.post_refine(&opts).unwrap();
        assert!(coarse.bitwise_eq(&before), "post_refine changed its input");
        let repair = report.repair.expect("the mesh needed the repair");
        assert!(repair.defects_before > 0 && repair.defects_after == 0 && repair.n_flips > 0);
        refined.certify(true).unwrap();
        assert!((refined.total_area() - coarse.total_area()).abs() < 1e-9 * coarse.total_area());
        check_prolongation(&refined, "boxpart repaired");
        repaired += 1;
        // The flips alone keep the metric: the repaired copy certifies with the same area.
        let mut copy = coarse.clone();
        let again = copy.make_simplicial(100_000).unwrap();
        assert_eq!(again.defects_after, 0);
        copy.certify(false).unwrap();
        assert!((copy.total_area() - coarse.total_area()).abs() < 1e-9 * coarse.total_area());
    }
    assert!(repaired > 0, "no coarse boxpart needed the repair");
}
