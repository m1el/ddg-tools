//! Every refusal of the pipeline, with its message in the caller's numbering.

use std::f64::consts::TAU;

use ddg_uv::primitives::{cube, holed_slab};
use ddg_uv::{
    ClosedSurfaces, Cone, ConeSigns, Degenerate, ErrorKind, MeshInput, Orientation, PackOptions,
    PolygonMesh, UnwrapOptions, UvError, unwrap,
};

fn flat(points: &[[f64; 2]]) -> Vec<[f64; 3]> {
    points.iter().map(|p| [p[0], p[1], 0.0]).collect()
}

/// Two triangles sharing the edge (1, 2), in the plane.
fn two_triangles() -> PolygonMesh {
    PolygonMesh::from_triangles(
        flat(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        &[[0, 1, 2], [1, 3, 2]],
    )
}

fn message(mesh: &PolygonMesh, options: &UnwrapOptions) -> (ErrorKind, String) {
    let e = unwrap(&mesh.input(), options).unwrap_err();
    (e.kind(), e.to_string())
}

fn input_error(input: &MeshInput<'_>) -> String {
    let e = unwrap(input, &UnwrapOptions::default()).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Input, "{e}");
    e.to_string()
}

#[test]
fn malformed_polygon_arrays() {
    let good = two_triangles();
    let d = UnwrapOptions::default();

    let mut m = good.clone();
    m.polygon_offsets = vec![0];
    assert_eq!(message(&m, &d), (ErrorKind::Input, "the mesh has no faces".into()));

    let mut m = good.clone();
    m.polygon_offsets = vec![1, 3, 6];
    assert_eq!(
        message(&m, &d).1,
        "polygon_offsets[0] = 1: the first face must start at corner 0"
    );

    let mut m = good.clone();
    m.polygon_offsets = vec![0, 3, 5];
    assert_eq!(
        message(&m, &d).1,
        "face 1 has 2 corners; a face needs at least 3"
    );

    let mut m = good.clone();
    m.polygon_offsets = vec![0, 4, 3];
    assert_eq!(
        message(&m, &d).1,
        "polygon_offsets[2] = 3: offsets must not decrease, but the previous one is 4"
    );

    let mut m = good.clone();
    m.polygon_offsets = vec![0, 3, 7];
    assert_eq!(
        message(&m, &d).1,
        "polygon_offsets[2] = 7: the last offset must be the number of corners, 6"
    );

    let mut m = good.clone();
    m.positions[2][1] = f64::NAN;
    assert_eq!(message(&m, &d).1, "positions[2] = NaN is not a finite number");

    let mut m = good.clone();
    m.triangles[1] = [3, 4, 6];
    assert_eq!(
        message(&m, &d).1,
        "triangles[1] = 6 is out of range: there are 6 polygon corners (valid indices 0..=5)"
    );

    let mut m = good.clone();
    m.triangles.pop();
    assert_eq!(
        message(&m, &d).1,
        "triangles has 1 triangles, but the faces need 2 (corners minus twice the number of \
         faces)"
    );
}

#[test]
fn triangles_that_do_not_triangulate() {
    let d = UnwrapOptions::default();
    let good = two_triangles();

    let mut m = good.clone();
    m.triangles[1] = [3, 4, 2];
    assert_eq!(
        message(&m, &d).1,
        "the triangles of face 1 do not triangulate it: triangle 1 = corners (3, 4, 2) spans \
         the faces 1, 1 and 0"
    );

    let mut m = good.clone();
    m.triangles[0] = [0, 1, 1];
    assert_eq!(
        message(&m, &d).1,
        "the triangles of face 0 do not triangulate it: triangle 0 = corners (0, 1, 1) repeats \
         a corner"
    );

    let mut m = good.clone();
    m.triangles[0] = [0, 2, 1];
    assert_eq!(
        message(&m, &d).1,
        "the triangles of face 0 do not triangulate it: triangle 0 = corners (0, 2, 1) is \
         wound against the face"
    );

    // A quad whose two triangles use the same diagonal in the same direction.
    let quad = PolygonMesh {
        positions: flat(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        corner_vertices: vec![0, 1, 2, 3],
        polygon_offsets: vec![0, 4],
        triangles: vec![[0, 1, 2], [0, 1, 2]],
    };
    assert_eq!(
        message(&quad, &d).1,
        "the triangles of face 0 do not triangulate it: the side that starts at corner 0 is in \
         2 triangles instead of 1"
    );
    // A pentagon with a triangle too many on one side and a missing one elsewhere.
    let pentagon = PolygonMesh {
        positions: flat(&[[0.0, 0.0], [1.0, 0.0], [1.5, 1.0], [0.5, 1.5], [-0.5, 1.0]]),
        corner_vertices: vec![0, 1, 2, 3, 4],
        polygon_offsets: vec![0, 5],
        triangles: vec![[0, 1, 2], [0, 2, 3], [0, 2, 4]],
    };
    assert_eq!(
        message(&pentagon, &d).1,
        "the triangles of face 0 do not triangulate it: the side that starts at corner 3 is in \
         0 triangles instead of 1"
    );
}

#[test]
fn malformed_seams_and_steering() {
    let mesh = two_triangles();
    let seams = [[1, 2], [3, 3]];
    assert_eq!(
        input_error(&MeshInput {
            seam_edges: &seams,
            ..mesh.input()
        }),
        "seam_edges[1] joins vertex 3 to itself"
    );
    assert_eq!(
        input_error(&MeshInput {
            seam_edges: &[[1, 4]],
            ..mesh.input()
        }),
        "seam_edges[0] = 4 is out of range: there are 4 vertices (valid indices 0..=3)"
    );
    assert_eq!(
        input_error(&MeshInput {
            pinned_vertices: &[0, 9],
            ..mesh.input()
        }),
        "pinned_vertices[1] = 9 is out of range: there are 4 vertices (valid indices 0..=3)"
    );
    assert_eq!(
        input_error(&MeshInput {
            forbidden_vertices: &[4],
            ..mesh.input()
        }),
        "forbidden_vertices[0] = 4 is out of range: there are 4 vertices (valid indices 0..=3)"
    );
    let cone = |vertex, angle| Cone { vertex, angle };
    assert_eq!(
        input_error(&MeshInput {
            prescribed_cones: &[cone(7, 0.1)],
            ..mesh.input()
        }),
        "prescribed_cones[0] = 7 is out of range: there are 4 vertices (valid indices 0..=3)"
    );
    assert_eq!(
        input_error(&MeshInput {
            prescribed_cones: &[cone(1, f64::INFINITY)],
            ..mesh.input()
        }),
        "prescribed_cones[0] = inf is not a finite number"
    );
    assert_eq!(
        input_error(&MeshInput {
            prescribed_cones: &[cone(1, TAU)],
            ..mesh.input()
        }),
        format!(
            "the cone angle {TAU} prescribed at vertex 1 is not below 2π ({TAU}): no flattening \
             has such a cone"
        )
    );
    assert_eq!(
        input_error(&MeshInput {
            prescribed_cones: &[cone(1, 0.1), cone(1, 0.2)],
            ..mesh.input()
        }),
        "vertex 1 has more than one prescribed cone"
    );
    assert_eq!(
        input_error(&MeshInput {
            pinned_vertices: &[2],
            forbidden_vertices: &[2],
            ..mesh.input()
        }),
        "vertex 2 is both pinned and forbidden"
    );
    assert_eq!(
        input_error(&MeshInput {
            prescribed_cones: &[cone(3, 0.1)],
            forbidden_vertices: &[3],
            ..mesh.input()
        }),
        "vertex 3 is both prescribed and forbidden"
    );
}

#[test]
fn options_out_of_range() {
    let mesh = two_triangles();
    let d = UnwrapOptions::default();
    let cases: Vec<(UnwrapOptions, &str)> = vec![
        (
            UnwrapOptions { lambda: -1.0, ..d },
            "option lambda = -1: must be a finite number that is not negative",
        ),
        (
            UnwrapOptions {
                max_cone_angle: 7.0,
                ..d
            },
            "option max_cone_angle = 7: must be positive and at most 2π",
        ),
        (
            UnwrapOptions {
                packing: Some(PackOptions {
                    margin: 0.5,
                    rotate: true,
                }),
                ..d
            },
            "option packing.margin = 0.5: must be at least 0 and below 0.5",
        ),
        (
            UnwrapOptions {
                degenerate: Degenerate::MollifyIfNeeded {
                    relative_delta: 0.0,
                },
                ..d
            },
            "option degenerate.relative_delta = 0: must be positive and finite",
        ),
        (
            UnwrapOptions {
                closed_surfaces: ClosedSurfaces::FreeSites {
                    count: 2,
                    max_angle: 3.0,
                    forbid_elsewhere: true,
                },
                ..d
            },
            "option closed_surfaces.count = 2: with cones forbidden elsewhere, count times \
             max_angle must be at least 4π, the sum of the cone angles of a closed surface of \
             genus 0",
        ),
        (
            UnwrapOptions {
                overlap_samples_per_triangle: f64::NAN,
                ..d
            },
            "option overlap_samples_per_triangle = NaN: must be finite and not negative",
        ),
    ];
    for (options, expected) in cases {
        assert_eq!(message(&mesh, &options), (ErrorKind::Input, expected.to_string()));
    }
}

#[test]
fn meshes_that_are_not_manifold() {
    let d = UnwrapOptions::default();
    // Three triangles on the edge (0, 1).
    let fins = PolygonMesh::from_triangles(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 1, 2], [1, 0, 3], [0, 1, 4]],
    );
    assert_eq!(
        message(&fins, &d),
        (ErrorKind::Mesh, "edge (0, 1) is shared by 3 faces: 0, 1, 2".into())
    );

    // Two triangles that touch in one vertex only.
    let bowtie = PolygonMesh::from_triangles(
        flat(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]),
        &[[0, 1, 2], [0, 3, 4]],
    );
    assert_eq!(
        message(&bowtie, &d),
        (
            ErrorKind::Mesh,
            "vertex 0 is non-manifold: its incident faces form 2 separate fans".into()
        )
    );

    // The same triangle twice (vertex 4 is loose and does not matter).
    let twice = PolygonMesh::from_triangles(
        flat(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [5.0, 5.0], [9.0, 9.0]]),
        &[[0, 1, 2], [2, 0, 1]],
    );
    assert_eq!(
        message(&twice, &d).1,
        "faces 0 and 1 are duplicates or overlap: both cover the triangle of vertices (0, 1, 2)"
    );

    // A face that uses a vertex twice.
    let repeated = PolygonMesh::from_polygons(
        flat(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        &[vec![0, 1, 2, 1]],
    );
    assert_eq!(
        message(&repeated, &d),
        (ErrorKind::Mesh, "face 0 uses vertex 1 more than once".into())
    );

    // A quad split along (0, 2), and a triangle that uses (0, 2) as an edge on each side.
    let diagonal = PolygonMesh::from_polygons(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ],
        &[vec![0, 1, 2, 3], vec![2, 0, 4]],
    );
    assert_eq!(
        message(&diagonal, &d),
        (
            ErrorKind::Mesh,
            "edge (0, 2) is shared by 3 triangles of the triangulation, in faces 0, 1: a \
             diagonal of a face coincides with another edge"
                .into()
        )
    );
}

/// A strip of `n` quads closed with a half twist.
fn moebius(n: u32) -> PolygonMesh {
    let mut positions = Vec::new();
    for i in 0..n {
        let phi = TAU * f64::from(i) / f64::from(n);
        for side in [-0.3, 0.3] {
            let r = 1.0 + side * (0.5 * phi).cos();
            positions.push([r * phi.cos(), r * phi.sin(), side * (0.5 * phi).sin()]);
        }
    }
    let mut polygons = Vec::new();
    for i in 0..n {
        let (a, b) = (2 * i, 2 * i + 1);
        let (c, d) = if i + 1 < n { (2 * i + 2, 2 * i + 3) } else { (1, 0) };
        polygons.push(vec![a, c, d, b]);
    }
    PolygonMesh::from_polygons(positions, &polygons)
}

#[test]
fn orientation() {
    let d = UnwrapOptions::default();
    let (kind, text) = message(&moebius(8), &d);
    assert_eq!(kind, ErrorKind::Mesh);
    // Where the search for an orientation closes its loop depends on its order; here it is
    // inside face 1, across the diagonal of its fan.
    assert_eq!(
        text,
        "the surface is not orientable: the faces connected to face 1 have no consistent \
         orientation (the contradiction closes at the diagonal (2, 5) of that face)"
    );
    // On a strip of triangles the two faces are named.
    let mut strip = moebius(8);
    let faces: Vec<[u32; 3]> = strip
        .triangles
        .iter()
        .map(|t| t.map(|c| strip.corner_vertices[c as usize]))
        .collect();
    strip = PolygonMesh::from_triangles(strip.positions.clone(), &faces);
    let (kind, text) = message(&strip, &d);
    assert_eq!(kind, ErrorKind::Mesh);
    assert!(
        text.starts_with("the surface is not orientable: faces ")
            && text.contains(" cannot be oriented consistently across edge ("),
        "{text}"
    );

    // One of two triangles is wound the other way.
    let mut flipped = two_triangles();
    flipped.corner_vertices = vec![0, 1, 2, 1, 2, 3];
    let reject = UnwrapOptions {
        orientation: Orientation::Reject,
        ..d
    };
    assert_eq!(
        message(&flipped, &reject),
        (
            ErrorKind::Mesh,
            "faces 0 and 1 traverse their shared edge (1, 2) in the same direction: the mesh is \
             not consistently oriented (one of them is flipped)"
                .into()
        )
    );
    // By default it is repaired, counted and warned about.
    let map = unwrap(&flipped.input(), &d).unwrap();
    assert_eq!(map.report.validation.n_flipped_polygons, 1);
    assert_eq!(map.report.warnings.len(), 1);
    assert!(map.report.warnings[0].starts_with(
        "the faces were not consistently oriented: 1 of 2 faces were flipped (the first is \
         face 1)"
    ));
    // The map is the same as that of the consistently wound mesh, corner by corner.
    let good = unwrap(&two_triangles().input(), &d).unwrap();
    assert_eq!(map.uvs[0..3], good.uvs[0..3]);
    assert_eq!(
        [map.uvs[3], map.uvs[4], map.uvs[5]],
        [good.uvs[3], good.uvs[5], good.uvs[4]]
    );
}

#[test]
fn degenerate_triangles() {
    // A quad with three collinear corners: the fan from corner 0 has a triangle without area.
    let mesh = PolygonMesh::from_polygons(
        flat(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        &[vec![0, 1, 2, 3], vec![0, 3, 4]],
    );
    let d = UnwrapOptions::default();
    let reject = UnwrapOptions {
        degenerate: Degenerate::Reject,
        ..d
    };
    assert_eq!(
        message(&mesh, &reject),
        (
            ErrorKind::Mesh,
            "1 triangle of the triangulation has no area; the first is in face 0 (vertices 0, \
             1, 2; edge lengths 1, 1, 2)"
                .into()
        )
    );
    // The default mollifies, and says so.
    let map = unwrap(&mesh.input(), &d).unwrap();
    let report = &map.report.validation.degenerate;
    assert!(report.mollified);
    assert_eq!(report.n_degenerate_triangles, 1);
    assert_eq!(report.worst_polygon, Some(0));
    assert!(report.epsilon > 0.0 && report.epsilon < 1e-5, "{}", report.epsilon);
    assert_eq!(map.report.warnings.len(), 1);
    assert!(map.report.warnings[0].contains("intrinsic mollification"));
    assert!(map.uvs.iter().all(|q| q[0].is_finite() && q[1].is_finite()));
    // A mesh without degenerate triangles is not touched.
    let clean = unwrap(&two_triangles().input(), &d).unwrap();
    assert!(!clean.report.validation.degenerate.mollified);
    assert_eq!(clean.report.validation.degenerate.epsilon, 0.0);
}

#[test]
fn hostile_geometry_is_refused_without_a_panic() {
    let d = UnwrapOptions::default();
    // Every vertex in the same place: mollification cannot help (the mean edge length is 0).
    let mut collapsed = two_triangles();
    collapsed.positions = vec![[1.0, 2.0, 3.0]; 4];
    assert_eq!(
        message(&collapsed, &d),
        (
            ErrorKind::Mesh,
            "2 triangles of the triangulation have no area even after mollifying the edge \
             lengths; the first is in face 0 (vertices 0, 1, 2; edge lengths 0, 0, 0)"
                .into()
        )
    );
    // Coordinates whose differences overflow.
    let mut huge = two_triangles();
    huge.positions[3] = [1e308, -1e308, 1e308];
    huge.positions[0] = [-1e308, 1e308, -1e308];
    let e = unwrap(&huge.input(), &d).unwrap_err();
    assert_ne!(e.kind(), ErrorKind::Solve, "{e}");
    // The summary of a report never contains a NaN, no key twice, and none of the keys the
    // Python binding sets itself.
    let map = unwrap(&two_triangles().input(), &d).unwrap();
    let summary = map.report.summary();
    let reserved = [
        "pipeline", "placeholder", "n_vertices", "n_polygons", "n_corners", "n_triangles",
        "n_seam_edges", "n_pinned", "n_prescribed", "n_forbidden", "lambda",
        "positive_cones_only", "n_charts", "n_cones", "n_cut_edges", "n_degenerate_triangles",
        "n_collapsed_uv_triangles", "n_flipped_uv_triangles", "measured_triangles",
        "seconds_pipeline", "seconds_measure", "seconds_decode", "seconds_validate",
        "seconds_encode", "seconds_total", "warnings",
    ];
    for (i, (key, value)) in summary.iter().enumerate() {
        if let ddg_uv::Value::Float(x) = value {
            assert!(!x.is_nan(), "{key}");
        }
        assert!(!reserved.contains(key) && !key.starts_with("measured_"), "{key}");
        assert!(summary[..i].iter().all(|(other, _)| other != key), "{key} twice");
    }
}

#[test]
fn closed_surfaces_and_infeasible_cones() {
    let d = UnwrapOptions::default();
    let refuse = UnwrapOptions {
        closed_surfaces: ClosedSurfaces::Refuse,
        ..d
    };
    assert_eq!(
        message(&cube(2.0), &refuse),
        (
            ErrorKind::Mesh,
            "the piece of the surface that contains face 0 (6 faces, genus 0) is closed: mark a \
             seam on it"
                .into()
        )
    );
    let positive = UnwrapOptions {
        cone_signs: ConeSigns::PositiveOnly,
        ..d
    };
    let (kind, text) = message(&holed_slab(2, 1, 1.0), &positive);
    assert_eq!(kind, ErrorKind::Mesh);
    assert_eq!(
        text,
        format!(
            "the piece of the surface that contains face 0 is closed with genus 2: its cone \
             angles must sum to {} (2π times the Euler characteristic), which positive cones \
             cannot do; allow negative cones or mark seams",
            -2.0 * TAU
        )
    );
    // Cones forbidden at every vertex of a closed surface.
    let mesh = cube(2.0);
    let all: Vec<u32> = (0..8).collect();
    let raw = UnwrapOptions {
        closed_surfaces: ClosedSurfaces::RawOptimum,
        ..d
    };
    let e = unwrap(
        &MeshInput {
            forbidden_vertices: &all,
            ..mesh.input()
        },
        &raw,
    )
    .unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Mesh);
    assert!(
        e.to_string().starts_with(
            "the cone constraints on the piece of the surface that contains face 0 cannot be \
             met: infeasible: on this closed mesh (Euler characteristic 2)"
        ),
        "{e}"
    );
    // The paper's single free site takes all of 4π on a cube.
    let one_site = UnwrapOptions {
        closed_surfaces: ClosedSurfaces::FreeSites {
            count: 1,
            max_angle: f64::INFINITY,
            forbid_elsewhere: false,
        },
        max_cone_angle: TAU,
        ..d
    };
    let e = unwrap(&mesh.input(), &one_site).unwrap_err();
    assert_eq!(e.kind(), ErrorKind::Solve);
    assert!(matches!(e, UvError::UnrealizableCone { vertex: 0, angle } if angle > TAU), "{e}");
}

#[test]
fn messages_of_the_failures_that_no_test_input_triggers() {
    // These four are raised when the solver or the flattener fails or a final check does;
    // no input of the tests makes that happen, so only their wording is checked here.
    let cases = [
        (
            UvError::ConeSolver {
                polygon: 3,
                message: "x".into(),
            },
            "the cone solver failed on the piece of the surface that contains face 3: x",
        ),
        (
            UvError::ConesNotCertified {
                polygon: 3,
                message: "status NotConverged".into(),
            },
            "the cones of the piece of the surface that contains face 3 are not a certified \
             optimum: status NotConverged",
        ),
        (
            UvError::Flatten {
                polygon: 4,
                message: "y".into(),
            },
            "flattening failed on the piece of the surface that contains face 4: y",
        ),
        (
            UvError::Internal {
                message: "z".into(),
            },
            "internal error: z",
        ),
    ];
    for (error, expected) in cases {
        assert_eq!(error.kind(), ErrorKind::Solve);
        assert_eq!(error.to_string(), expected);
    }
    assert_eq!(
        UvError::TooLarge {
            what: "positions",
            count: 5_000_000_000
        }
        .to_string(),
        "positions has 5000000000 entries, more than 32-bit ids can address (limit 4294967294)"
    );
}
