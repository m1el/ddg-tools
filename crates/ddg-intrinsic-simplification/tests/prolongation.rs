//! The compacted output, the scalar and vector prolongation and the colouring of the coarse
//! faces (SPEC §7; WP4). The coarse meshes are made by removing vertices in a fixed order, so
//! nothing here depends on the ICE metric.

#![allow(clippy::needless_range_loop)]

mod common;

use std::f64::consts::{PI, TAU};

use common::{Shape, crane_models, planar_grid, reference_meshes, small_shapes};
use ddg_core::generators::SplitMix64;
use ddg_core::{Geometry, TriMesh};
use ddg_intrinsic::IntrinsicTriangulation;
use ddg_intrinsic_apps::vector::{DirectionFieldOptions, smoothest_field};
use ddg_intrinsic_simplification::{
    C64, CoarseMesh, CoarsePoint, FineFrames, PointLocation, ProlongationError,
    VectorProlongationReport, normalize, unnormalize,
};

/// Removes vertices in the order of `order` (passes until `target` vertices remain or nothing
/// changes), skipping those for which `allow` is false and those that cannot be removed.
/// Returns the number removed.
fn coarsen(
    mesh: &mut CoarseMesh,
    order: &[u32],
    target: usize,
    allow: impl Fn(&CoarseMesh, u32) -> bool,
) -> usize {
    let mut removed = 0;
    loop {
        let mut progress = false;
        for &v in order {
            if mesh.n_vertices() <= target {
                return removed;
            }
            if !mesh.is_vertex(v) || !allow(mesh, v) {
                continue;
            }
            if mesh.remove_vertex(v).is_ok() {
                removed += 1;
                progress = true;
            }
        }
        if !progress {
            return removed;
        }
    }
}

fn shuffled(n: usize, seed: u64) -> Vec<u32> {
    let mut rng = SplitMix64::new(seed);
    let mut v: Vec<u32> = (0..n as u32).collect();
    for i in (1..n).rev() {
        let j = (rng.next_f64() * (i + 1) as f64) as usize % (i + 1);
        v.swap(i, j);
    }
    v
}

fn fine_frames(shape: &Shape) -> (TriMesh, Geometry, FineFrames) {
    let mesh = TriMesh::new(shape.positions.len(), shape.faces.clone()).unwrap();
    let geometry = Geometry::from_positions(&mesh, &shape.positions).unwrap();
    let frames = FineFrames::from_mesh(&mesh, &geometry).unwrap();
    (mesh, geometry, frames)
}

/// A flat grid coarsened by removing only vertices of zero curvature (everything but the
/// corners), so that no flattening happens and the map stays the identity of the plane.
fn flat_grid_coarse(n: usize, jitter: f64, keep_fraction: f64) -> (Shape, CoarseMesh) {
    let shape = planar_grid(n, jitter, 5);
    let points = shape.centroids();
    let mut mesh = shape.coarse(&points);
    mesh.flip_to_delaunay().unwrap();
    let target = (shape.positions.len() as f64 * keep_fraction) as usize;
    let order = shuffled(shape.positions.len(), 11);
    coarsen(&mut mesh, &order, target, |m, v| {
        m.curvature(m.complex_vertex(v).unwrap()).abs() < 1e-9
    });
    (shape, mesh)
}

fn count_rows(r: &VectorProlongationReport) -> usize {
    r.n_same_face
        + r.n_adjacent_face
        + r.n_kept_from_images
        + r.n_averaged
        + r.n_kept_averaged
        + r.n_default
}

#[test]
fn compact_output_agrees_with_the_mesh() {
    for shape in small_shapes() {
        let points = shape.random_points(2);
        let mut mesh = shape.coarse(&points);
        mesh.flip_to_delaunay().unwrap();
        let order = shuffled(shape.positions.len(), 3);
        coarsen(&mut mesh, &order, shape.positions.len() / 3, |_, _| true);
        let out = mesh.compact().unwrap();
        assert_eq!(out.n_vertices(), mesh.n_vertices(), "{}", shape.name);
        assert_eq!(out.n_faces(), mesh.n_faces(), "{}", shape.name);
        assert_eq!(out.faces_fine(), mesh.faces(), "{}", shape.name);
        assert_eq!(out.points.len(), mesh.n_points());
        for p in 0..mesh.n_points() {
            match (out.points[p], mesh.point_location(p).unwrap()) {
                (CoarsePoint::Vertex(v), PointLocation::Vertex(fine)) => {
                    assert_eq!(out.vertices[v as usize], fine);
                }
                (
                    CoarsePoint::Face { face, bary },
                    PointLocation::Face {
                        face: f2,
                        corners,
                        bary: b2,
                    },
                ) => {
                    assert_eq!(face, f2);
                    assert_eq!(bary, b2);
                    assert_eq!(
                        out.faces[face as usize].map(|v| out.vertices[v as usize]),
                        corners
                    );
                    assert!((bary[0] + bary[1] + bary[2] - 1.0).abs() < 1e-12);
                }
                (a, b) => panic!("{}: point {p}: {a:?} vs {b:?}", shape.name),
            }
        }
        // Lengths of the sides of a face are those of its edges.
        for f in 0..out.n_faces() {
            for k in 0..3 {
                assert_eq!(
                    out.face_lengths[f][k],
                    out.edge_lengths[out.face_edges[f][k] as usize]
                );
            }
        }
        // Fine -> coarse ids are inverse.
        for (c, &fine) in out.vertices.iter().enumerate() {
            assert_eq!(out.coarse_of_fine[fine as usize], Some(c as u32));
        }
    }
}

#[test]
fn compact_refuses_a_pending_removal() {
    let shape = planar_grid(4, 0.3, 1);
    let mut mesh = shape.coarse(&[]);
    mesh.flip_to_delaunay().unwrap();
    let _prepared = mesh.begin_removal(12).unwrap();
    assert!(mesh.compact().is_err());
    mesh.revert_removal().unwrap();
    assert!(mesh.compact().is_ok());
}

#[test]
fn scalar_prolongation_reproduces_linear_functions_on_a_flat_grid() {
    for (n, jitter, keep) in [(10, 0.0, 0.3), (14, 0.4, 0.15), (20, 0.4, 0.05)] {
        let (shape, mesh) = flat_grid_coarse(n, jitter, keep);
        let out = mesh.compact().unwrap();
        let p = out.scalar_prolongation().unwrap();
        assert_eq!(p.nrows(), shape.positions.len());
        assert_eq!(p.ncols(), out.n_vertices());
        // Rows sum to one; P 1 = 1.
        let ones = vec![1.0; p.ncols()];
        let mut y = vec![0.0; p.nrows()];
        p.mul_vec(&ones, &mut y).unwrap();
        let max_row = y.iter().map(|s| (s - 1.0).abs()).fold(0.0, f64::max);
        assert!(max_row < 1e-14, "row sums off by {max_row}");
        // Coordinate functions.
        let mut max_err: f64 = 0.0;
        for c in 0..2 {
            let coarse: Vec<f64> = out
                .vertices
                .iter()
                .map(|&v| shape.positions[v as usize][c])
                .collect();
            p.mul_vec(&coarse, &mut y).unwrap();
            for (i, pos) in shape.positions.iter().enumerate() {
                max_err = max_err.max((y[i] - pos[c]).abs());
            }
        }
        assert!(max_err < 1e-7, "{}: linear error {max_err}", shape.name);
        // The extra points (centroids) too.
        let pp = out.point_prolongation().unwrap();
        assert_eq!(pp.nrows(), out.points.len());
        let mut yy = vec![0.0; pp.nrows()];
        let mut max_pt: f64 = 0.0;
        for c in 0..2 {
            let coarse: Vec<f64> = out
                .vertices
                .iter()
                .map(|&v| shape.positions[v as usize][c])
                .collect();
            pp.mul_vec(&coarse, &mut yy).unwrap();
            for (f, pt) in shape.centroids().iter().enumerate() {
                let want = shape.input_position(pt)[c];
                max_pt = max_pt.max((yy[shape.positions.len() + f] - want).abs());
            }
        }
        assert!(max_pt < 1e-7, "{}: centroid error {max_pt}", shape.name);
        eprintln!(
            "{}: {} -> {} vertices, linear reproduction error {:e} (vertices), {:e} (centroids), \
             row sums {:e}",
            shape.name,
            shape.positions.len(),
            out.n_vertices(),
            max_err,
            max_pt,
            max_row
        );
    }
}

#[test]
fn scalar_prolongation_on_curved_models_has_unit_row_sums() {
    for shape in small_shapes()
        .into_iter()
        .chain(crane_models())
        .chain(reference_meshes(3_000))
    {
        let mut mesh = shape.coarse(&[]);
        mesh.flip_to_delaunay().unwrap();
        let order = shuffled(shape.positions.len(), 17);
        coarsen(&mut mesh, &order, shape.positions.len() / 8, |_, _| true);
        let out = mesh.compact().unwrap();
        let p = out.scalar_prolongation().unwrap();
        let mut y = vec![0.0; p.nrows()];
        p.mul_vec(&vec![1.0; p.ncols()], &mut y).unwrap();
        let worst = y.iter().map(|s| (s - 1.0).abs()).fold(0.0, f64::max);
        assert!(worst < 1e-13, "{}: row sums off by {worst}", shape.name);
        assert!(p.values().iter().all(|&b| (0.0..=1.0 + 1e-12).contains(&b)));
    }
}

/// The orientation of a vertex frame: global direction of a halfedge minus its unnormalized
/// angle in the frame. Returns the first and the largest disagreement over the halfedges.
fn frame_orientation(angles: &[(f64, f64)]) -> (f64, f64) {
    let r0 = angles[0].0 - angles[0].1;
    let spread = angles
        .iter()
        .map(|&(g, a)| {
            let d = (g - a - r0).rem_euclid(TAU);
            d.min(TAU - d)
        })
        .fold(0.0, f64::max);
    (r0, spread)
}

/// Prolongation of constant (n-direction) fields on a flat mesh without flattening.
fn constant_field_error(shape: &Shape, mesh: &CoarseMesh, gamma: f64, n_sym: u32) -> (f64, usize) {
    let out = mesh.compact().unwrap();
    let (tmesh, geometry, frames) = fine_frames(shape);
    let (p, report) = out.vector_prolongation(&frames, n_sym).unwrap();
    assert_eq!(p.n_empty_rows(), 0);
    assert_eq!(count_rows(&report), shape.positions.len());
    if n_sym == 1 && gamma == 0.0 {
        eprintln!("  rows of {}: {report:?}", shape.name);
    }
    let n = n_sym as f64;
    let pos = |fine: u32| shape.positions[fine as usize];
    let mut halfedges: Vec<Vec<(f64, f64)>> = vec![Vec::new(); out.n_vertices()];
    for f in 0..out.n_faces() {
        for k in 0..3 {
            let v = out.faces[f][k] as usize;
            let w = out.faces[f][(k + 1) % 3];
            let (a, b) = (pos(out.vertices[v]), pos(out.vertices[w as usize]));
            let g = (b[1] - a[1]).atan2(b[0] - a[0]);
            let s = if out.boundary[v] { PI } else { TAU } / out.angle_sum[v];
            halfedges[v].push((g, out.phi[3 * f + k] / s));
        }
    }
    // The coarse unnormalized field: n (gamma - orientation of the frame).
    let mut coarse = Vec::new();
    for h in &halfedges {
        let (r0, spread) = frame_orientation(h);
        assert!(spread < 1e-7, "coarse frame inconsistent by {spread}");
        coarse.push(C64::from_angle(n * (gamma - r0)));
    }
    // The fine expectation from the signposts (scale one away from the corners).
    let tri = IntrinsicTriangulation::new(&tmesh, &geometry).unwrap();
    let c = tri.complex();
    let mut max_err: f64 = 0.0;
    let mut checked = 0;
    let fine = p.mul_vec(&coarse).unwrap();
    for i in 0..shape.positions.len() {
        let theta = tri.vertex_angle_sum(i);
        let boundary = c.is_boundary_vertex(i);
        if (theta - if boundary { PI } else { TAU }).abs() > 1e-9 {
            continue; // a corner: its frame is normalized with a scale other than one
        }
        let h = c.outgoing_halfedge(i);
        let (a, b) = (pos(i as u32), pos(c.to_vertex(h) as u32));
        let g = (b[1] - a[1]).atan2(b[0] - a[0]);
        let r0 = g - tri.signpost(h);
        let want = C64::from_angle(n * (gamma - r0));
        max_err = max_err.max((fine[i] - want).abs());
        checked += 1;
    }
    (max_err, checked)
}

#[test]
fn constant_fields_are_reproduced_on_a_flat_mesh() {
    for (n, jitter, keep) in [(8, 0.0, 0.3), (12, 0.4, 0.15), (18, 0.4, 0.04)] {
        let (shape, mesh) = flat_grid_coarse(n, jitter, keep);
        for n_sym in [1u32, 2, 4] {
            for gamma in [0.0, 0.7, 2.3] {
                let (err, checked) = constant_field_error(&shape, &mesh, gamma, n_sym);
                eprintln!(
                    "{} ({} coarse vertices): n={n_sym} gamma={gamma}: max error {err:e} at {checked} vertices",
                    shape.name,
                    mesh.n_vertices()
                );
                assert!(err < 1e-6, "{}: error {err}", shape.name);
            }
        }
    }
}

#[test]
fn prolongation_commutes_with_rotating_the_field() {
    let (shape, mesh) = flat_grid_coarse(12, 0.4, 0.12);
    let out = mesh.compact().unwrap();
    let (_, _, frames) = fine_frames(&shape);
    let (p, _) = out.vector_prolongation(&frames, 1).unwrap();
    let mut rng = SplitMix64::new(9);
    let delta = 1.1;
    let u0: Vec<C64> = (0..out.n_vertices())
        .map(|_| C64::from_angle(TAU * rng.next_f64()).scale(0.5 + rng.next_f64()))
        .collect();
    let rot = C64::from_angle(delta);
    let u = out.unnormalize_field(&u0).unwrap();
    let ur: Vec<C64> = u.iter().map(|&x| x * rot).collect();
    let a = p.mul_vec(&u).unwrap();
    let b = p.mul_vec(&ur).unwrap();
    let err = a
        .iter()
        .zip(&b)
        .map(|(&x, &y)| (x * rot - y).abs())
        .fold(0.0, f64::max);
    assert!(err < 1e-12, "rotation consistency {err}");
    // At a vertex of cone angle 2 pi (or pi on the boundary) the frames are metric, and
    // rotating the unnormalized field rotates the normalized one.
    let unrot = out.normalize_field(&ur).unwrap();
    let mut worst = 0.0f64;
    for v in 0..out.n_vertices() {
        if (out.angle_sum[v] - if out.boundary[v] { PI } else { TAU }).abs() < 1e-9 {
            worst = worst.max((unrot[v] - u0[v] * rot).abs());
        }
    }
    assert!(
        worst < 1e-12,
        "normalize/unnormalize on flat vertices {worst}"
    );
}

#[test]
fn normalize_and_unnormalize_are_inverse() {
    let mut rng = SplitMix64::new(5);
    for _ in 0..200 {
        let theta = 0.5 + 5.5 * rng.next_f64();
        let u = C64::from_angle(TAU * rng.next_f64()).scale(0.1 + rng.next_f64());
        let back = normalize(unnormalize(u, theta, false), theta, false);
        assert!((back - u).abs() < 1e-12);
        // Boundary vectors point into the half plane: arg in [0, pi].
        let ub = C64::from_angle(PI * rng.next_f64()).scale(0.1 + rng.next_f64());
        let back = normalize(unnormalize(ub, theta, true), theta, true);
        assert!((back - ub).abs() < 1e-12);
        // The normalized direction pi (half way around) is the unnormalized direction theta / 2.
        let half = unnormalize(C64::from_angle(PI), theta, false);
        assert!((C64::from_angle(theta / 2.0) - half).abs() < 1e-12);
    }
}

#[test]
fn bad_arguments_are_errors() {
    let (shape, mesh) = flat_grid_coarse(6, 0.0, 0.5);
    let out = mesh.compact().unwrap();
    let (_, _, frames) = fine_frames(&shape);
    assert_eq!(
        out.vector_prolongation(&frames, 0).unwrap_err(),
        ProlongationError::InvalidSymmetry(0)
    );
    let other = planar_grid(3, 0.0, 1);
    let (_, _, wrong) = fine_frames(&other);
    assert!(matches!(
        out.vector_prolongation(&wrong, 1),
        Err(ProlongationError::Mismatch { .. })
    ));
    let (p, _) = out.vector_prolongation(&frames, 1).unwrap();
    assert!(p.mul_vec(&[C64::default()]).is_err());
    assert!(out.unnormalize_field(&[]).is_err());
}

#[test]
fn no_empty_rows_and_fallback_counts() {
    let mut total = VectorProlongationReport::default();
    let mut rows_total = 0usize;
    for shape in small_shapes()
        .into_iter()
        .chain(crane_models())
        .chain(reference_meshes(3_000))
    {
        let (_, _, frames) = fine_frames(&shape);
        for keep in [4usize, 10, 40] {
            let mut mesh = shape.coarse(&[]);
            mesh.flip_to_delaunay().unwrap();
            let order = shuffled(shape.positions.len(), 23);
            coarsen(&mut mesh, &order, shape.positions.len() / keep, |_, _| true);
            let out = mesh.compact().unwrap();
            for n_sym in [1u32, 4] {
                let (p, r) = out.vector_prolongation(&frames, n_sym).unwrap();
                assert_eq!(p.n_empty_rows(), 0, "{}", shape.name);
                assert_eq!(
                    count_rows(&r),
                    shape.positions.len(),
                    "{}: {r:?}",
                    shape.name
                );
                assert!(
                    p.values
                        .iter()
                        .all(|a| a.re.is_finite() && a.im.is_finite())
                );
                if n_sym == 1 {
                    eprintln!(
                        "{} 1/{keep}: {} -> {} vertices: {r:?}",
                        shape.name,
                        shape.positions.len(),
                        out.n_vertices()
                    );
                    rows_total += shape.positions.len();
                    total.n_same_face += r.n_same_face;
                    total.n_adjacent_face += r.n_adjacent_face;
                    total.n_kept_from_images += r.n_kept_from_images;
                    total.n_averaged += r.n_averaged;
                    total.n_kept_averaged += r.n_kept_averaged;
                    total.n_default += r.n_default;
                    total.n_cancelled += r.n_cancelled;
                    total.n_coincident_neighbours += r.n_coincident_neighbours;
                }
            }
        }
    }
    eprintln!("all runs, {rows_total} rows: {total:?}");
    assert_eq!(total.n_default, 0);
}

/// Percentile of a vector (sorted in place).
fn percentile(v: &mut [f64], q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[((v.len() - 1) as f64 * q) as usize]
}

/// The smoothest field on the coarse mesh, prolonged, against the smoothest field on the fine
/// mesh: the mean, median and 90th percentile of the angle between the two directions (after
/// the best global phase), in degrees, and the coarse vertex count.
fn smooth_field_agreement(
    shape: &Shape,
    keep: usize,
    n_sym: u32,
) -> Option<(f64, f64, f64, usize)> {
    let (tmesh, geometry, frames) = fine_frames(shape);
    let mut mesh = shape.coarse(&[]);
    mesh.flip_to_delaunay().unwrap();
    let order = shuffled(shape.positions.len(), 31);
    coarsen(&mut mesh, &order, shape.positions.len() / keep, |_, _| true);
    let out = mesh.compact().unwrap();
    // The coarse mesh as an input mesh of its own (needs a simplicial complex).
    let cmesh = match TriMesh::new(out.n_vertices(), out.faces.clone()) {
        Ok(m) => m,
        Err(e) => {
            eprintln!(
                "{}: coarse mesh is not simplicial ({e}); skipped",
                shape.name
            );
            return None;
        }
    };
    // Lengths per edge of the TriMesh from the halfedges (3f + k is the same halfedge).
    let mut lengths = vec![f64::NAN; cmesh.n_edges()];
    for f in 0..out.n_faces() {
        for k in 0..3 {
            let e = cmesh.edge(3 * f + k);
            let l = out.face_lengths[f][k];
            assert!(lengths[e].is_nan() || (lengths[e] - l).abs() < 1e-12 * (1.0 + l));
            lengths[e] = l;
        }
    }
    let cgeom = Geometry::from_edge_lengths(&cmesh, lengths).ok()?;
    let ctri = IntrinsicTriangulation::new(&cmesh, &cgeom).unwrap();
    let ftri = IntrinsicTriangulation::new(&tmesh, &geometry).unwrap();
    let opts = DirectionFieldOptions {
        n_sym,
        ..DirectionFieldOptions::default()
    };
    let (cf, _) = smoothest_field(&ctri, &opts).unwrap();
    let (ff, _) = smoothest_field(&ftri, &opts).unwrap();
    // The frame of the apps' coarse triangulation -> the frame of the coarse mesh: a rotation
    // per (interior) vertex, n times for an n-field.
    let c = ctri.complex();
    let mut u0 = Vec::new();
    for v in 0..out.n_vertices() {
        if out.boundary[v] {
            return None;
        }
        let h = c.outgoing_halfedge(v);
        let delta = out.phi[h] - ctri.signpost(h);
        let x = cf.values[v];
        u0.push(C64::new(x.re, x.im) * C64::from_angle(n_sym as f64 * delta));
    }
    let u = out.unnormalize_field(&u0).unwrap();
    let (p, _) = out.vector_prolongation(&frames, n_sym).unwrap();
    let pu = p.mul_vec(&u).unwrap();
    // Best global phase, then the angle between the directions (n-th root of the n-th powers).
    let mut inner = C64::default();
    for (a, b) in pu.iter().zip(&ff.values) {
        inner += a.conj() * C64::new(b.re, b.im);
    }
    let phase = C64::from_angle(inner.arg());
    let mut errs: Vec<f64> = pu
        .iter()
        .zip(&ff.values)
        .map(|(&a, b)| ((a * phase) / C64::new(b.re, b.im)).arg().abs() / n_sym as f64)
        .collect();
    let mean = errs.iter().sum::<f64>() / errs.len() as f64;
    let med = percentile(&mut errs, 0.5);
    let p90 = percentile(&mut errs, 0.9);
    Some((
        mean.to_degrees(),
        med.to_degrees(),
        p90.to_degrees(),
        out.n_vertices(),
    ))
}

#[test]
fn smooth_fields_prolong_close_to_the_fine_solution() {
    let mut any = false;
    for shape in crane_models()
        .into_iter()
        .take(1)
        .chain(small_shapes().into_iter().take(2))
    {
        for n_sym in [1u32, 2] {
            for keep in [2usize, 4, 10] {
                if let Some((mean, med, p90, nc)) = smooth_field_agreement(&shape, keep, n_sym) {
                    eprintln!(
                        "{}: n={n_sym}, {} -> {nc} vertices: angle to the fine smoothest field: \
                         mean {mean:.2} deg, median {med:.2}, 90th percentile {p90:.2}",
                        shape.name,
                        shape.positions.len()
                    );
                    any = true;
                    if shape.name.starts_with("spot") {
                        assert!(mean < 30.0, "{}: mean angle error {mean}", shape.name);
                    }
                }
            }
        }
    }
    assert!(any, "no mesh could be compared");
}

#[test]
fn colouring_is_valid_and_deterministic() {
    for shape in small_shapes()
        .into_iter()
        .chain(crane_models())
        .chain(reference_meshes(3_000))
    {
        let points = shape.centroids();
        let mut mesh = shape.coarse(&points);
        mesh.flip_to_delaunay().unwrap();
        let order = shuffled(shape.positions.len(), 41);
        coarsen(&mut mesh, &order, shape.positions.len() / 10, |_, _| true);
        let out = mesh.compact().unwrap();
        let c1 = out.colour_faces();
        let c2 = mesh.compact().unwrap().colour_faces();
        assert_eq!(c1, c2, "{}: colouring is not deterministic", shape.name);
        assert_eq!(c1.colours.len(), out.n_faces());
        let mut by_vertex: Vec<Vec<u32>> = vec![Vec::new(); out.n_vertices()];
        for (f, corners) in out.faces.iter().enumerate() {
            for &v in corners {
                by_vertex[v as usize].push(f as u32);
            }
        }
        for list in &by_vertex {
            for (a, &f) in list.iter().enumerate() {
                for &g in &list[a + 1..] {
                    if f != g {
                        assert_ne!(
                            c1.colours[f as usize], c1.colours[g as usize],
                            "{}: faces {f} and {g} share a vertex",
                            shape.name
                        );
                    }
                }
            }
        }
        let max_degree = by_vertex.iter().map(|l| l.len()).max().unwrap();
        assert!(c1.colours.iter().all(|&c| c < c1.n_colours));
        assert_eq!(c1.n_colours, c1.colours.iter().max().map_or(0, |m| m + 1));
        // Every point has a face; a point inside a face is coloured as the face.
        let faces = out.face_of_points();
        let colours = out.point_colours();
        assert!(faces.iter().all(|f| f.is_some()), "{}", shape.name);
        for (p, loc) in out.points.iter().enumerate() {
            if let CoarsePoint::Face { face, .. } = loc {
                assert_eq!(faces[p], Some(*face));
                assert_eq!(colours[p], Some(c1.colours[*face as usize]));
            }
            assert_eq!(out.face_of_point(p), faces[p]);
        }
        eprintln!(
            "{}: {} coarse faces, {} colours (max vertex degree {max_degree})",
            shape.name,
            out.n_faces(),
            c1.n_colours
        );
    }
}
