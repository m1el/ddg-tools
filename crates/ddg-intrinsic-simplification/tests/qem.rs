//! The QEM baseline (SPEC §10.4, S23): planar meshes, a sphere, boundaries, determinism, the
//! counts, the maps, bad input.

use ddg_core::TriMesh;
use ddg_core::generators::{grid_disk, icosphere, open_cylinder, torus};
use ddg_intrinsic_simplification::baseline::{
    QemError, QemOptions, Quadric, optimal_point, optimal_point_on_segment, plane_quadric,
    qem_simplify,
};

fn norm(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

#[test]
fn plane_quadric_is_the_squared_distance_and_the_optimum_solves_the_system() {
    // The planes x = 1, y = 2, z = 3: the optimum is (1, 2, 3) with error 0.
    let q = [
        plane_quadric([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),
        plane_quadric([0.0, 2.0, 0.0], [0.0, 1.0, 0.0], 1.0),
        plane_quadric([0.0, 0.0, 3.0], [0.0, 0.0, 1.0], 1.0),
    ]
    .iter()
    .fold(Quadric::ZERO, |s, q| s.plus(q));
    let p = optimal_point(&q, 1e-8).unwrap();
    for (x, e) in p.iter().zip([1.0, 2.0, 3.0]) {
        assert!((x - e).abs() < 1e-14);
    }
    assert!(q.error(p).abs() < 1e-14);
    assert!((q.error([2.0, 2.0, 3.0]) - 1.0).abs() < 1e-14);
    // One plane: squared distance times the weight; the system is singular.
    let n = [0.0, 0.6, 0.8];
    let one = plane_quadric([1.0, 1.0, 1.0], n, 2.5);
    let x = [3.0, -1.0, 4.0];
    let dist = n[0] * (x[0] - 1.0) + n[1] * (x[1] - 1.0) + n[2] * (x[2] - 1.0);
    assert!((one.error(x) - 2.5 * dist * dist).abs() < 1e-12);
    assert_eq!(optimal_point(&one, 1e-8), None);
    // A random-ish non-orthogonal triple: the gradient vanishes at the optimum.
    let q = [
        plane_quadric([0.3, 0.1, 0.0], [0.6, 0.0, 0.8], 1.0),
        plane_quadric([0.0, 0.7, 0.2], [0.0, 1.0, 0.0], 2.0),
        plane_quadric([0.5, 0.5, 0.5], [0.8, 0.6, 0.0], 0.5),
    ]
    .iter()
    .fold(Quadric::ZERO, |s, q| s.plus(q));
    let p = optimal_point(&q, 1e-12).unwrap();
    let e0 = q.error(p);
    for k in 0..3 {
        for s in [-1e-4, 1e-4] {
            let mut r = p;
            r[k] += s;
            assert!(q.error(r) > e0);
        }
    }
}

#[test]
fn planar_mesh_collapses_with_zero_error_and_keeps_its_area() {
    let g = grid_disk(12, 12).unwrap();
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(30)).unwrap();
    let rep = r.report;
    assert!(rep.reached_target, "{rep:?}");
    assert_eq!(rep.n_output_vertices, 30);
    assert!(rep.total_cost < 1e-12, "cost {}", rep.total_cost);
    assert!(
        rep.area_error_percent < 1e-9,
        "area error {}",
        rep.area_error_percent
    );
    assert!(
        rep.max_projection_distance < 1e-9,
        "{}",
        rep.max_projection_distance
    );
    for p in &r.positions {
        assert!(p[2].abs() < 1e-12);
    }
}

#[test]
fn sphere_keeps_its_shape() {
    let g = icosphere(3).unwrap();
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(100)).unwrap();
    assert!(r.report.reached_target);
    println!("sphere 642 -> 100: {:?}", r.report);
    // The coarse vertices stay on (or just inside) the unit sphere.
    for p in &r.positions {
        assert!((norm(*p) - 1.0).abs() < 0.03, "radius {}", norm(*p));
    }
    // Hausdorff distance of the fine vertices from the coarse surface, and the other way round
    // through the coarse vertices (the sphere is the exact reference).
    assert!(
        r.report.max_projection_distance < 0.05,
        "{}",
        r.report.max_projection_distance
    );
    // The area of a piecewise linear sphere with 100 vertices is within a few percent.
    assert!(
        r.report.area_error_percent < 8.0,
        "{}",
        r.report.area_error_percent
    );
    let coarse = TriMesh::new(r.positions.len(), r.faces.clone()).unwrap();
    assert_eq!(coarse.euler_characteristic(), 2);
}

#[test]
fn boundary_stays_on_the_boundary_of_the_square() {
    let g = grid_disk(10, 10).unwrap();
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(25)).unwrap();
    let coarse = TriMesh::new(r.positions.len(), r.faces.clone()).unwrap();
    let mut n_boundary = 0;
    for v in 0..coarse.n_vertices() {
        if coarse.is_boundary_vertex(v) {
            n_boundary += 1;
            let p = r.positions[v];
            let on = p[0].abs() < 1e-9
                || (p[0] - 1.0).abs() < 1e-9
                || p[1].abs() < 1e-9
                || (p[1] - 1.0).abs() < 1e-9;
            assert!(on, "boundary vertex {v} at {p:?}");
        }
    }
    assert!(n_boundary >= 4);
    // The four corners survive: the area is that of the square.
    assert!(
        r.report.area_error_percent < 1e-9,
        "{}",
        r.report.area_error_percent
    );
    // The topology is that of the disk.
    assert_eq!(coarse.euler_characteristic(), 1);
}

#[test]
fn deterministic() {
    let g = torus(16, 9, 1.0, 0.4).unwrap();
    let o = QemOptions::new(40);
    let a = qem_simplify(&g.positions, &g.faces, &o).unwrap();
    let b = qem_simplify(&g.positions, &g.faces, &o).unwrap();
    assert_eq!(a, b);
}

#[test]
fn counts_maps_and_topology() {
    for (g, target) in [
        (icosphere(3).unwrap(), 200),
        (icosphere(2).unwrap(), 60),
        (torus(16, 9, 1.0, 0.4).unwrap(), 50),
    ] {
        let n = g.positions.len();
        let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(target)).unwrap();
        assert!(r.report.reached_target);
        assert_eq!(r.positions.len(), target);
        assert_eq!(r.report.n_collapses, n - target);
        let coarse = TriMesh::new(r.positions.len(), r.faces.clone()).unwrap();
        // A closed surface: V - E + F = chi, and the input's chi is kept.
        let fine = TriMesh::new(n, g.faces.clone()).unwrap();
        assert_eq!(coarse.euler_characteristic(), fine.euler_characteristic());
        assert_eq!(
            r.faces.len(),
            2 * (target as i64 - fine.euler_characteristic()) as usize
        );
        // The vertex map is onto and idempotent on survivors.
        assert_eq!(r.vertex_map.len(), n);
        let mut hit = vec![false; target];
        for &c in &r.vertex_map {
            hit[c as usize] = true;
        }
        assert!(hit.iter().all(|&h| h));
        // The projection: valid faces and barycentric coordinates, distance consistent.
        assert_eq!(r.projection.len(), n);
        for (v, pr) in r.projection.iter().enumerate() {
            assert!((pr.face as usize) < r.faces.len());
            let sum: f64 = pr.bary.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12 && pr.bary.iter().all(|&b| b >= -1e-12));
            let t = r.faces[pr.face as usize].map(|x| r.positions[x as usize]);
            let q: [f64; 3] = std::array::from_fn(|k| {
                pr.bary[0] * t[0][k] + pr.bary[1] * t[1][k] + pr.bary[2] * t[2][k]
            });
            let d = norm([
                q[0] - g.positions[v][0],
                q[1] - g.positions[v][1],
                q[2] - g.positions[v][2],
            ]);
            assert!((d - pr.distance).abs() < 1e-12);
        }
    }
}

#[test]
fn bad_and_degenerate_input_is_an_error_or_survives() {
    let o = QemOptions::new(3);
    let tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(matches!(
        qem_simplify(&tri, &[], &o),
        Err(QemError::InvalidMesh(_))
    ));
    assert!(matches!(
        qem_simplify(&tri, &[[0, 1, 5]], &o),
        Err(QemError::InvalidMesh(_))
    ));
    let mut nan = tri;
    nan[1][0] = f64::NAN;
    assert_eq!(
        qem_simplify(&nan, &[[0, 1, 2]], &o).unwrap_err(),
        QemError::NonFinitePosition { vertex: 1 }
    );
    let four = [tri[0], tri[1], tri[2], [5.0, 5.0, 5.0]];
    assert!(matches!(
        qem_simplify(&four, &[[0, 1, 2]], &o),
        Err(QemError::InvalidMesh(m)) if m.contains("not referenced")
    ));
    // Three faces on one edge.
    let fan = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    assert!(matches!(
        qem_simplify(&fan, &[[0, 1, 2], [1, 0, 3], [0, 1, 4]], &o),
        Err(QemError::InvalidMesh(_))
    ));
    let mut bad = QemOptions::new(3);
    bad.boundary_weight = f64::NAN;
    assert!(matches!(
        qem_simplify(&tri, &[[0, 1, 2]], &bad),
        Err(QemError::InvalidOptions(_))
    ));
    // A single triangle: nothing to collapse, the target is not reached, no panic.
    let r = qem_simplify(&tri, &[[0, 1, 2]], &QemOptions::new(2)).unwrap();
    assert!(!r.report.reached_target);
    assert_eq!(r.positions.len(), 3);
    // Zero-area faces (collinear points) and coincident points: no panic.
    let flat = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
    ];
    let _ = qem_simplify(
        &flat,
        &[[0, 1, 2], [2, 1, 3], [2, 3, 4]],
        &QemOptions::new(3),
    )
    .unwrap();
    let _ = qem_simplify(
        &flat,
        &[[0, 1, 2], [2, 1, 3], [2, 3, 4]],
        &QemOptions::new(0),
    )
    .unwrap();
    // The target above the vertex count and target 0 on a closed mesh (stops at a tetrahedron
    // or earlier, never panics).
    let g = icosphere(1).unwrap();
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(1000)).unwrap();
    assert!(r.report.reached_target && r.report.n_collapses == 0);
    assert_eq!(r.positions, g.positions);
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(0)).unwrap();
    assert!(!r.report.reached_target);
    assert!(r.positions.len() >= 4);
    TriMesh::new(r.positions.len(), r.faces.clone()).unwrap();
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The distance from `p` to a triangle: the plane distance when the projection is inside,
/// else the smallest distance to the three segments (independent of the Voronoi-region code).
fn triangle_distance(p: [f64; 3], t: [[f64; 3]; 3]) -> f64 {
    let seg = |a: [f64; 3], b: [f64; 3]| {
        let ab = sub(b, a);
        let s = (dot(sub(p, a), ab) / dot(ab, ab)).clamp(0.0, 1.0);
        let q = [a[0] + s * ab[0], a[1] + s * ab[1], a[2] + s * ab[2]];
        norm(sub(p, q))
    };
    let mut best = seg(t[0], t[1]).min(seg(t[1], t[2])).min(seg(t[2], t[0]));
    let (e1, e2) = (sub(t[1], t[0]), sub(t[2], t[0]));
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    let n2 = dot(n, n);
    let w = sub(p, t[0]);
    let h = dot(w, n) / n2;
    let q = sub(w, [h * n[0], h * n[1], h * n[2]]);
    // Barycentric coordinates of the projection by Cramer's rule.
    let (d00, d01, d11) = (dot(e1, e1), dot(e1, e2), dot(e2, e2));
    let (d20, d21) = (dot(q, e1), dot(q, e2));
    let den = d00 * d11 - d01 * d01;
    let v = (d11 * d20 - d01 * d21) / den;
    let u = (d00 * d21 - d01 * d20) / den;
    if v >= 0.0 && u >= 0.0 && v + u <= 1.0 {
        best = best.min(h.abs() * n2.sqrt());
    }
    best
}

#[test]
fn projection_equals_brute_force() {
    let g = torus(14, 8, 1.0, 0.4).unwrap();
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(40)).unwrap();
    for (v, pr) in r.projection.iter().enumerate() {
        let best = r
            .faces
            .iter()
            .map(|f| triangle_distance(g.positions[v], f.map(|x| r.positions[x as usize])))
            .fold(f64::INFINITY, f64::min);
        assert!(
            (pr.distance - best).abs() < 1e-10,
            "vertex {v}: {} != {best}",
            pr.distance
        );
    }
}

// ---- review against Garland and Heckbert 1997 ----

#[test]
fn defaults_follow_the_1997_paper() {
    let o = QemOptions::new(10);
    // Sec. 5: every plane has weight 1 (area weighting is a later variant).
    assert!(!o.area_weighted);
    // Sec. 4: the optimal point is not bounded.
    assert!(o.placement_radius.is_infinite());
    // Sec. 6: a large penalty on the boundary constraint planes.
    assert!(o.boundary_weight >= 100.0);
    // Infinity is allowed for the radius, NaN and negative values are not.
    let g = icosphere(2).unwrap();
    assert!(qem_simplify(&g.positions, &g.faces, &o).is_ok());
    for bad in [f64::NAN, -1.0] {
        let mut b = o;
        b.placement_radius = bad;
        assert!(matches!(
            qem_simplify(&g.positions, &g.faces, &b),
            Err(QemError::InvalidOptions(_))
        ));
    }
}

#[test]
fn weighting_changes_the_costs() {
    let g = icosphere(3).unwrap();
    let plain = QemOptions::new(60);
    let weighted = QemOptions {
        area_weighted: true,
        ..plain
    };
    let a = qem_simplify(&g.positions, &g.faces, &plain).unwrap();
    let b = qem_simplify(&g.positions, &g.faces, &weighted).unwrap();
    // Unit weights against face areas (about 1e-2 on this sphere): the costs differ by the
    // order of magnitude of the areas.
    assert!(
        a.report.total_cost > 10.0 * b.report.total_cost,
        "{} vs {}",
        a.report.total_cost,
        b.report.total_cost
    );
    assert!(a.report.reached_target && b.report.reached_target);
}

#[test]
fn optimum_on_the_segment() {
    // The plane z = 0 only: the 3x3 system is singular; on a segment crossing the plane the
    // minimizer is the crossing.
    let q = plane_quadric([0.0; 3], [0.0, 0.0, 1.0], 1.0);
    assert_eq!(optimal_point(&q, 1e-8), None);
    let p = optimal_point_on_segment(&q, [1.0, 2.0, 1.0], [1.0, 2.0, -3.0], 1e-8).unwrap();
    for (x, e) in p.iter().zip([1.0, 2.0, 0.0]) {
        assert!((x - e).abs() < 1e-14, "{p:?}");
    }
    // Clamped to the segment: the plane is beyond the end.
    let p = optimal_point_on_segment(&q, [0.0, 0.0, 2.0], [0.0, 0.0, 1.0], 1e-8).unwrap();
    assert_eq!(p, [0.0, 0.0, 1.0]);
    // A segment inside the plane: the quadric is constant along it.
    assert_eq!(
        optimal_point_on_segment(&q, [0.0; 3], [1.0, 1.0, 0.0], 1e-8),
        None
    );
    // The result is never worse than the endpoints or the midpoint.
    let q2 = plane_quadric([0.0; 3], [0.6, 0.0, 0.8], 1.0).plus(&plane_quadric(
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        1.0,
    ));
    let (a, b) = ([2.0, 3.0, -1.0], [-1.0, 0.5, 2.5]);
    let p = optimal_point_on_segment(&q2, a, b, 1e-8).unwrap();
    let mid = [0.5, 1.75, 0.75];
    assert!(
        q2.error(p) <= q2.error(a) && q2.error(p) <= q2.error(b) && q2.error(p) <= q2.error(mid)
    );
}

#[test]
fn singular_systems_use_the_segment_placement() {
    // A prism over a wavy closed curve: every face plane is parallel to the axis, so A has rank
    // 2 and the 3x3 system is singular (no boundary penalty here). Collapses along a ring edge
    // go to the optimum on the segment, which lies strictly inside it.
    let mut g = open_cylinder(30, 6, 1.0, 4.0).unwrap();
    for p in &mut g.positions {
        let theta = p[1].atan2(p[0]);
        let r = 1.0 + 0.2 * (5.0 * theta).sin();
        p[0] *= r;
        p[1] *= r;
    }
    let o = QemOptions {
        boundary_weight: 0.0,
        ..QemOptions::new(30)
    };
    let r = qem_simplify(&g.positions, &g.faces, &o).unwrap();
    assert!(r.report.n_collapses > 0);
    assert!(r.report.n_placed_segment > 0, "{:?}", r.report);
}

/// A vertex of degree three collapses into any neighbour (its three faces become two and one that
/// is not doubled); the link test must not reject it, or the run stalls far above a small target.
#[test]
fn a_sphere_collapses_to_a_handful_of_vertices() {
    let g = icosphere(3).unwrap();
    let r = qem_simplify(&g.positions, &g.faces, &QemOptions::new(8)).unwrap();
    assert!(r.report.reached_target, "{:?}", r.report);
    assert_eq!(r.positions.len(), 8);
    // Euler: V - E + F = 2 with E = 3F/2.
    assert_eq!(r.faces.len(), 2 * 8 - 4);
}
