//! Tests of the conformal map to the sphere (`map_to_sphere`).

// Index loops mirror the subscripts of the matrix formulas.
#![allow(clippy::needless_range_loop)]

mod common;

use std::path::PathBuf;

use common::flat_polyhedron;
use ddg_bff::sphere::{inverse_stereographic, sphere_inversion};
use ddg_bff::{
    BffError, CenteringOptions, CenteringStatus, PoleChoice, SphereMap, SphereOptions,
    map_to_sphere, mobius_center, sphere_distortion,
};
use ddg_core::generators::{self, GeneratedMesh};
use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};

fn setup(positions: &[[f64; 3]], faces: &[[u32; 3]]) -> (TriMesh, Geometry) {
    let mesh = TriMesh::new(positions.len(), faces.to_vec()).unwrap();
    let geometry = Geometry::from_positions(&mesh, positions).unwrap();
    (mesh, geometry)
}

fn setup_generated(g: &GeneratedMesh) -> (TriMesh, Geometry) {
    setup(&g.positions, &g.faces)
}

/// The models of `data/models` that are closed of genus 0 (spot, blub), if present.
/// A named mesh: positions and faces.
type Case = (String, Vec<[f64; 3]>, Vec<[u32; 3]>);

fn genus_zero_models() -> Vec<Case> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models");
    [
        "spot/spot/spot_triangulated.obj",
        "blub/blub/blub_triangulated.obj",
    ]
    .iter()
    .filter_map(|p| {
        let path = root.join(p);
        let obj = read_obj(&path).ok()?;
        let faces = obj.triangles()?;
        Some((p.to_string(), obj.positions, faces))
    })
    .collect()
}

fn norm(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
}

fn mat_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn inverse_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let cof =
        |r0: usize, r1: usize, c0: usize, c1: usize| m[r0][c0] * m[r1][c1] - m[r0][c1] * m[r1][c0];
    let c = [
        [cof(1, 2, 1, 2), -cof(1, 2, 0, 2), cof(1, 2, 0, 1)],
        [-cof(0, 2, 1, 2), cof(0, 2, 0, 2), -cof(0, 2, 0, 1)],
        [cof(0, 1, 1, 2), -cof(0, 1, 0, 2), cof(0, 1, 0, 1)],
    ];
    let det = m[0][0] * c[0][0] + m[0][1] * c[0][1] + m[0][2] * c[0][2];
    c.map(|row| row.map(|x| x / det))
}

/// The rotation `R` that best maps the points `p` to `q` (polar factor of `Σ q pᵀ`), and the
/// largest `|R p_i − q_i|`.
fn rotation_fit(p: &[[f64; 3]], q: &[[f64; 3]]) -> ([[f64; 3]; 3], f64) {
    let mut m = [[0.0; 3]; 3];
    for (a, b) in p.iter().zip(q) {
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += b[i] * a[j];
            }
        }
    }
    // Newton iteration for the polar factor: X ← ½ (X + X⁻ᵀ).
    let mut x = m;
    for _ in 0..100 {
        let it = inverse_transpose(&x);
        let next: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (x[i][j] + it[i][j])));
        let change: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| (i, j)))
            .map(|(i, j)| (next[i][j] - x[i][j]).abs())
            .fold(0.0, f64::max);
        x = next;
        if change < 1e-15 {
            break;
        }
    }
    // Orthogonal, and a rotation, not a reflection.
    let xt: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| x[j][i]));
    let xtx = mat_mul(&xt, &x);
    for i in 0..3 {
        for j in 0..3 {
            let id = if i == j { 1.0 } else { 0.0 };
            assert!(
                (xtx[i][j] - id).abs() < 1e-9,
                "polar factor is not orthogonal"
            );
        }
    }
    let det = x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1])
        - x[0][1] * (x[1][0] * x[2][2] - x[1][2] * x[2][0])
        + x[0][2] * (x[1][0] * x[2][1] - x[1][1] * x[2][0]);
    assert!(
        (det - 1.0).abs() < 1e-9,
        "best fit is not a rotation: det {det}"
    );
    let worst = p
        .iter()
        .zip(q)
        .map(|(a, b)| {
            let r: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| x[i][j] * a[j]).sum());
            norm([r[0] - b[0], r[1] - b[1], r[2] - b[2]])
        })
        .fold(0.0, f64::max);
    (x, worst)
}

/// The checks every successful map must pass.
fn check_map(name: &str, mesh: &TriMesh, map: &SphereMap) {
    let r = &map.report;
    assert_eq!(map.positions.len(), mesh.n_vertices(), "{name}");
    for (v, p) in map.positions.iter().enumerate() {
        assert!(
            (norm(*p) - 1.0).abs() < 1e-12,
            "{name}: vertex {v} at radius {}",
            norm(*p)
        );
    }
    assert!(r.distortion.max_radius_error < 1e-12, "{name}");
    let d = &r.distortion.faces;
    assert_eq!(d.flipped_triangles, 0, "{name}: {d:?}");
    assert_eq!(d.degenerate_triangles, 0, "{name}: {d:?}");
    assert!(
        (r.distortion.covering_degree - 1.0).abs() < 1e-9,
        "{name}: {:?}",
        r.distortion
    );
    assert!(r.is_injective(), "{name}");
    let c = r.centering.expect("centering is on");
    assert_eq!(c.status, CenteringStatus::Converged, "{name}: {c:?}");
    assert!(c.final_center <= 1e-10, "{name}: {c:?}");
    // The pole is a vertex of the input and the disk flattening of the rest is clean.
    assert!((r.pole as usize) < mesh.n_vertices());
    assert_eq!(r.disk.distortion.flipped_triangles, 0, "{name}");
    assert!(r.disk_solve_residual < 1e-9, "{name}");
}

#[test]
fn icosphere_maps_to_itself_up_to_a_rotation() {
    let g = generators::icosphere(4).unwrap();
    let (mesh, geometry) = setup_generated(&g);
    let map = map_to_sphere(&mesh, &geometry, &SphereOptions::default()).unwrap();
    check_map("icosphere 4", &mesh, &map);
    let d = &map.report.distortion.faces;
    // Angle distortion: about 1 everywhere except next to the removed star.
    assert!(d.quasi_conformal.mean < 1.005, "{d:?}");
    assert!(d.quasi_conformal.max < 1.3, "{d:?}");
    // Area distortion after centering: the scale factors are nearly constant.
    assert!(d.log_scale_deviation < 0.01, "{d:?}");
    assert!(d.log_scale.min > -0.1 && d.log_scale.max < 0.1, "{d:?}");
    // The input points are on the unit sphere and the map is close to a rotation of them.
    let (_, worst) = rotation_fit(&g.positions, &map.positions);
    assert!(worst < 0.02, "largest deviation from a rotation: {worst}");
}

#[test]
fn without_centering_the_pole_is_at_the_south_pole() {
    let g = generators::icosphere(2).unwrap();
    let (mesh, geometry) = setup_generated(&g);
    let options = SphereOptions {
        pole: PoleChoice::Vertex(17),
        centering: None,
        ..SphereOptions::default()
    };
    let map = map_to_sphere(&mesh, &geometry, &options).unwrap();
    assert_eq!(map.report.pole, 17);
    assert_eq!(map.positions[17], [0.0, 0.0, -1.0]);
    assert!(map.report.centering.is_none());
    assert_eq!(map.report.distortion.faces.flipped_triangles, 0);
    // The radius search moved the disk off the northern hemisphere.
    let radius = map.report.radius.unwrap();
    assert!(radius.radius > 1.0 && radius.iterations > 0, "{radius:?}");
}

#[test]
fn centering_puts_the_area_center_at_the_origin() {
    let g = generators::icosphere(3).unwrap();
    let (mesh, geometry) = setup_generated(&g);
    // Start from a badly off-center map: no radius search, no centering.
    let options = SphereOptions {
        radius_search: None,
        centering: None,
        ..SphereOptions::default()
    };
    let mut positions = map_to_sphere(&mesh, &geometry, &options).unwrap().positions;
    let before = sphere_distortion(&mesh, &geometry, &positions).unwrap();
    let report = mobius_center(
        &mesh,
        &geometry,
        &mut positions,
        &CenteringOptions::default(),
    )
    .unwrap();
    assert_eq!(report.status, CenteringStatus::Converged, "{report:?}");
    assert!(report.initial_center > 0.5, "{report:?}");
    assert!(report.final_center <= 1e-10);
    // The center of mass of the surface area, carried to the sphere, is at the origin.
    let total: f64 = geometry.face_areas().iter().sum();
    let mut mu = [0.0; 3];
    for (face, area) in mesh.faces().iter().zip(geometry.face_areas()) {
        let s: [f64; 3] =
            std::array::from_fn(|i| face.iter().map(|&v| positions[v as usize][i]).sum());
        let n = norm(s);
        for i in 0..3 {
            mu[i] += area / total * s[i] / n;
        }
    }
    assert!(norm(mu) < 1e-10, "{mu:?}");
    let after = sphere_distortion(&mesh, &geometry, &positions).unwrap();
    assert!(
        after.faces.log_scale_deviation < 0.1 * before.faces.log_scale_deviation,
        "area distortion {} -> {}",
        before.faces.log_scale_deviation,
        after.faces.log_scale_deviation
    );
    // Inversions are conformal; what changes is the angle distortion of the flat triangles
    // spanned by the image points, which is largest where triangles on the sphere are large.
    assert!(after.faces.quasi_conformal.mean < before.faces.quasi_conformal.mean);
    // A step limit that is hit is counted.
    let mut positions = map_to_sphere(&mesh, &geometry, &options).unwrap().positions;
    let short = CenteringOptions {
        max_step: 0.05,
        ..CenteringOptions::default()
    };
    let report = mobius_center(&mesh, &geometry, &mut positions, &short).unwrap();
    assert!(
        report.shortened_steps > 0 && report.converged(),
        "{report:?}"
    );
}

#[test]
fn genus_zero_meshes_map_injectively() {
    let mut cases: Vec<Case> = Vec::new();
    for (n, level) in [(4, 0), (4, 2), (8, 0), (8, 2), (20, 1)] {
        let (p, f) = flat_polyhedron(n, level);
        cases.push((format!("polyhedron {n} level {level}"), p, f));
    }
    {
        let mut g = generators::icosphere(3).unwrap();
        for p in &mut g.positions {
            p[0] *= 2.0;
            p[2] *= 0.5;
        }
        generators::perturb(&mut g.positions, [0.01; 3], 3);
        cases.push(("noisy ellipsoid".into(), g.positions, g.faces));
    }
    {
        // A bumpy sphere: radius 1 + 0.3 sin(3x) sin(3y) sin(3z).
        let mut g = generators::icosphere(3).unwrap();
        for p in &mut g.positions {
            let r = 1.0 + 0.3 * (3.0 * p[0]).sin() * (3.0 * p[1]).sin() * (3.0 * p[2]).sin();
            *p = p.map(|x| x * r);
        }
        cases.push(("bumpy sphere".into(), g.positions, g.faces));
    }
    cases.extend(genus_zero_models());
    for (name, positions, faces) in &cases {
        let (mesh, geometry) = setup(positions, faces);
        let map = map_to_sphere(&mesh, &geometry, &SphereOptions::default())
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        check_map(name, &mesh, &map);
        // The reference's pole works as well.
        let options = SphereOptions {
            pole: PoleChoice::Vertex(0),
            ..SphereOptions::default()
        };
        let map = map_to_sphere(&mesh, &geometry, &options).unwrap();
        check_map(name, &mesh, &map);
        assert_eq!(map.report.pole, 0);
    }
}

#[test]
fn rigid_motion_and_scale_do_not_change_the_map() {
    let mut g = generators::icosphere(3).unwrap();
    generators::perturb(&mut g.positions, [0.05; 3], 11);
    let (mesh, geometry) = setup_generated(&g);
    let map = map_to_sphere(&mesh, &geometry, &SphereOptions::default()).unwrap();
    let (s, c) = (0.6f64.sin(), 0.6f64.cos());
    for scale in [1.0, 1e-3, 250.0] {
        let moved: Vec<[f64; 3]> = g
            .positions
            .iter()
            .map(|p| {
                let r = [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]];
                [scale * r[0] + 3.0, scale * r[1] - 7.0, scale * r[2] + 0.5]
            })
            .collect();
        let (mesh2, geometry2) = setup(&moved, &g.faces);
        let map2 = map_to_sphere(&mesh2, &geometry2, &SphereOptions::default()).unwrap();
        assert_eq!(map.report.pole, map2.report.pole);
        // The map depends on edge lengths only; what remains is rounding, and the rotation
        // left free by the centering.
        let (_, worst) = rotation_fit(&map.positions, &map2.positions);
        assert!(worst < 1e-8, "scale {scale}: {worst}");
        let (a, b) = (&map.report.distortion.faces, &map2.report.distortion.faces);
        assert!((a.quasi_conformal.mean - b.quasi_conformal.mean).abs() < 1e-10);
        assert!((a.log_scale_deviation - b.log_scale_deviation).abs() < 1e-10);
    }
}

#[test]
fn different_poles_give_the_same_map_up_to_rotation() {
    let g = generators::icosphere(4).unwrap();
    let (mesh, geometry) = setup_generated(&g);
    let maps: Vec<SphereMap> = [0u32, 1000, 2000]
        .iter()
        .map(|&v| {
            let options = SphereOptions {
                pole: PoleChoice::Vertex(v),
                ..SphereOptions::default()
            };
            map_to_sphere(&mesh, &geometry, &options).unwrap()
        })
        .collect();
    for other in &maps[1..] {
        let (_, worst) = rotation_fit(&maps[0].positions, &other.positions);
        // Discretization error of the map near the two removed stars.
        assert!(worst < 0.05, "{worst}");
    }
}

#[test]
fn surfaces_that_are_not_spheres_are_refused() {
    let disk = generators::grid_disk(4, 4).unwrap();
    let (mesh, geometry) = setup_generated(&disk);
    let err = map_to_sphere(&mesh, &geometry, &SphereOptions::default()).unwrap_err();
    assert!(
        matches!(
            err,
            BffError::NotASphere {
                components: 1,
                euler_characteristic: 1,
                boundary_loops: 1
            }
        ),
        "{err}"
    );

    let torus = generators::torus(12, 8, 2.0, 0.5).unwrap();
    let (mesh, geometry) = setup_generated(&torus);
    let err = map_to_sphere(&mesh, &geometry, &SphereOptions::default()).unwrap_err();
    assert!(
        matches!(
            err,
            BffError::NotASphere {
                euler_characteristic: 0,
                boundary_loops: 0,
                ..
            }
        ),
        "{err}"
    );
    assert!(err.to_string().contains("genus 1"), "{err}");

    // Two spheres.
    let a = generators::icosphere(1).unwrap();
    let n = a.positions.len() as u32;
    let mut positions = a.positions.clone();
    positions.extend(a.positions.iter().map(|p| [p[0] + 3.0, p[1], p[2]]));
    let mut faces = a.faces.clone();
    faces.extend(a.faces.iter().map(|f| f.map(|v| v + n)));
    let (mesh, geometry) = setup(&positions, &faces);
    let err = map_to_sphere(&mesh, &geometry, &SphereOptions::default()).unwrap_err();
    assert!(
        matches!(err, BffError::NotASphere { components: 2, .. }),
        "{err}"
    );
    assert!(err.to_string().contains("one by one"), "{err}");

    // bob is a torus.
    let bob =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models/bob/bob/bob_tri.obj");
    if let Ok(obj) = read_obj(&bob) {
        let (mesh, geometry) = setup(&obj.positions, &obj.triangles().unwrap());
        let err = map_to_sphere(&mesh, &geometry, &SphereOptions::default()).unwrap_err();
        assert!(
            matches!(
                err,
                BffError::NotASphere {
                    euler_characteristic: 0,
                    ..
                }
            ),
            "{err}"
        );
    }
}

#[test]
fn bad_options_are_refused() {
    let g = generators::icosphere(1).unwrap();
    let (mesh, geometry) = setup_generated(&g);
    let err = map_to_sphere(
        &mesh,
        &geometry,
        &SphereOptions {
            pole: PoleChoice::Vertex(42),
            ..SphereOptions::default()
        },
    )
    .unwrap_err();
    assert!(
        matches!(
            err,
            BffError::PoleOutOfRange {
                vertex: 42,
                n_vertices: 42
            }
        ),
        "{err}"
    );
    for options in [
        SphereOptions {
            diagonal_shift: f64::NAN,
            ..SphereOptions::default()
        },
        SphereOptions {
            centering: Some(CenteringOptions {
                max_step: 1.0,
                ..CenteringOptions::default()
            }),
            ..SphereOptions::default()
        },
        SphereOptions {
            radius_search: Some(ddg_bff::RadiusSearch {
                min: 0.0,
                ..Default::default()
            }),
            ..SphereOptions::default()
        },
    ] {
        let err = map_to_sphere(&mesh, &geometry, &options).unwrap_err();
        assert!(matches!(err, BffError::NonFiniteInput { .. }), "{err}");
    }
    let mut positions = vec![[0.0, 0.0, 1.0]; 3];
    let err = mobius_center(
        &mesh,
        &geometry,
        &mut positions,
        &CenteringOptions::default(),
    )
    .unwrap_err();
    assert!(matches!(err, BffError::LengthMismatch { .. }), "{err}");
}

#[test]
fn the_formulas_hold_numerically() {
    // Inverse stereographic projection is conformal: the Jacobian columns are orthogonal and of
    // equal length 2 / (1 + |w|²), and it preserves orientation (outward normal).
    for w in [[0.0, 0.0], [0.4, -1.1], [3.0, 2.0]] {
        let h = 1e-6;
        let p = inverse_stereographic(w);
        let dx = inverse_stereographic([w[0] + h, w[1]]);
        let dy = inverse_stereographic([w[0], w[1] + h]);
        let jx: [f64; 3] = std::array::from_fn(|i| (dx[i] - p[i]) / h);
        let jy: [f64; 3] = std::array::from_fn(|i| (dy[i] - p[i]) / h);
        let expected = 2.0 / (1.0 + w[0] * w[0] + w[1] * w[1]);
        assert!((norm(jx) / expected - 1.0).abs() < 1e-5);
        assert!((norm(jy) / expected - 1.0).abs() < 1e-5);
        let dot: f64 = (0..3).map(|i| jx[i] * jy[i]).sum();
        assert!(dot.abs() < 1e-5 * expected * expected);
        let n = [
            jx[1] * jy[2] - jx[2] * jy[1],
            jx[2] * jy[0] - jx[0] * jy[2],
            jx[0] * jy[1] - jx[1] * jy[0],
        ];
        assert!(n[0] * p[0] + n[1] * p[1] + n[2] * p[2] > 0.0);
    }
    // The derivative of the inversion at c = 0 is 2 (I − x xᵀ).
    let x = inverse_stereographic([0.3, 0.8]);
    let h = 1e-7;
    for k in 0..3 {
        let mut c = [0.0; 3];
        c[k] = h;
        let y = sphere_inversion(c, x);
        for i in 0..3 {
            let id = if i == k { 1.0 } else { 0.0 };
            let expected = 2.0 * (id - x[i] * x[k]);
            assert!(((y[i] - x[i]) / h - expected).abs() < 1e-5);
        }
    }
}
