//! The option `boundary: BoundaryValues::Optimized`: closed pieces unchanged, the flattening
//! realizes the solver's boundary values on well-behaved meshes, seams are honoured, the seam
//! length mismatch is reported (and checked against an independent computation), determinism.

mod common;

use common::{check, without_times};
use ddg_uv::primitives::{End, revolution, uv_sphere};
use ddg_uv::{
    BoundaryValues, MeshInput, PolygonMesh, UnwrapOptions, UvMap, flatten, prepare, solve_cones,
    unwrap,
};

fn vase() -> PolygonMesh {
    let profile: Vec<(f64, f64)> = (0..16)
        .map(|i| {
            let t = f64::from(i) / 15.0;
            (1.0 + 0.6 * (3.0 * t).sin() + 0.3 * t, 2.0 * t)
        })
        .collect();
    revolution(&profile, 32, End::Open, End::Open).unwrap()
}

fn ring(i: u32, s: u32) -> Vec<[u32; 2]> {
    (0..s).map(|j| [i * s + j, i * s + (j + 1) % s]).collect()
}

fn optimized() -> UnwrapOptions {
    UnwrapOptions {
        boundary: BoundaryValues::Optimized,
        ..UnwrapOptions::default()
    }
}

/// The mismatch of the two images of every edge in `edges`, from the UVs per corner.
fn mismatch(mesh: &PolygonMesh, map: &UvMap, edges: &[[u32; 2]]) -> f64 {
    let mut images: Vec<([u32; 2], f64)> = Vec::new();
    for p in 0..mesh.n_polygons() {
        let corners: Vec<usize> = mesh.polygon_corners(p).collect();
        for (i, &c) in corners.iter().enumerate() {
            let d = corners[(i + 1) % corners.len()];
            let (a, b) = (mesh.corner_vertices[c], mesh.corner_vertices[d]);
            let key = [a.min(b), a.max(b)];
            if edges.iter().any(|&e| e == key || e == [key[1], key[0]]) {
                let (u, v) = (map.uvs[c], map.uvs[d]);
                images.push((key, (u[0] - v[0]).hypot(u[1] - v[1])));
            }
        }
    }
    let mut worst = 0.0f64;
    for (i, (ka, la)) in images.iter().enumerate() {
        for (kb, lb) in &images[i + 1..] {
            if ka == kb {
                worst = worst.max((la - lb).abs() / la.max(*lb));
            }
        }
    }
    worst
}

#[test]
fn closed_pieces_are_exactly_as_without_the_option() {
    let mesh = uv_sphere(16, 8, 1.0).unwrap();
    let plain = unwrap(&mesh.input(), &UnwrapOptions::default()).unwrap();
    let free = unwrap(&mesh.input(), &optimized()).unwrap();
    assert_eq!(without_times(plain), without_times(free));
}

#[test]
fn a_seam_loop_gives_finite_seam_honouring_deterministic_results() {
    let mesh = vase();
    let seams = ring(4, 32);
    let input = MeshInput {
        seam_edges: &seams,
        ..mesh.input()
    };
    let options = optimized();
    let a = unwrap(&input, &options).unwrap();
    check(&mesh, &input, &a, &options);
    for e in &seams {
        let key = [e[0].min(e[1]), e[0].max(e[1])];
        assert!(a.cut_edges.binary_search(&key).is_ok(), "{key:?} not cut");
    }
    assert_eq!(a.report.charts.len(), 2);
    assert!(a.report.charts.iter().all(|c| c.boundary_optimized));
    let b = unwrap(&input, &options).unwrap();
    assert_eq!(without_times(a.clone()), without_times(b));

    // The seam length mismatch is real (the option does something), reported, and equals an
    // independent computation from the UVs.
    let reported = a.report.seam_lengths;
    assert_eq!(reported.n_seam_edges, 32);
    assert!(reported.max_seam_mismatch > 0.01, "{reported:?}");
    let own = mismatch(&mesh, &a, &seams);
    assert!(
        (own - reported.max_seam_mismatch).abs() < 1e-9,
        "{own} against {}",
        reported.max_seam_mismatch
    );
    // Without the option the sides of the seam keep their lengths.
    let plain = unwrap(&input, &UnwrapOptions::default()).unwrap();
    assert!(plain.report.seam_lengths.max_seam_mismatch < 1e-3);
    assert!(plain.report.charts.iter().all(|c| !c.boundary_optimized));
    // Lower distortion by the solver's own measure on both pieces.
    assert!(a.report.area_distortion.scale_free < plain.report.area_distortion.scale_free);
}

#[test]
fn the_flattening_realizes_the_solvers_boundary_values() {
    let mesh = vase();
    for seams in [Vec::new(), ring(4, 32), ring(9, 32)] {
        let input = MeshInput {
            seam_edges: &seams,
            ..mesh.input()
        };
        let options = optimized();
        let prepared = prepare(&input, &options).unwrap();
        let cones = solve_cones(&prepared, &options).unwrap();
        let flat = flatten(&prepared, &cones, &options).unwrap();
        for (k, piece) in prepared.pieces.iter().enumerate() {
            let (pc, pf) = (&cones.pieces[k], &flat.pieces[k]);
            assert!(pc.boundary_optimized && pf.boundary_values_used);
            let u = pc.log_scale.as_ref().unwrap();
            let on = pc.slit_piece.as_ref().unwrap_or(piece);
            let mut nonzero = false;
            let mut worst = 0.0f64;
            for h in (0..on.mesh.n_halfedges()).filter(|&h| on.mesh.twin(h).is_none()) {
                let (a, b) = (on.mesh.from_vertex(h), on.mesh.to_vertex(h));
                let (p, q) = (
                    mesh.positions[on.vertex_input[a] as usize],
                    mesh.positions[on.vertex_input[b] as usize],
                );
                let length =
                    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
                let ratio = (0.5 * (u[a] + u[b])).exp();
                nonzero |= (ratio - 1.0).abs() > 1e-3;
                let (x, y) = (pf.corner_uvs[h], pf.corner_uvs[on.mesh.next(h)]);
                let got = (x[0] - y[0]).hypot(x[1] - y[1]);
                worst = worst.max((got - ratio * length).abs() / (ratio * length));
            }
            assert!(
                nonzero,
                "piece {k}: the solver's boundary values are all zero"
            );
            assert!(worst < 1e-3, "piece {k}: layout off by {worst}");
            // The pipeline's own figure agrees.
            assert!((pf.report.boundary_length_error - worst).abs() < 1e-9 + 0.01 * worst);
        }
    }
}
