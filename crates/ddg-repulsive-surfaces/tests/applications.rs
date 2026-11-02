//! Smoke tests of the applications of Sec. 9 (`SPEC.md` §14.4) at reduced size: a few steps of
//! each, checking that the flow runs without error and keeps its constraints. The full runs are
//! `examples/apps` (see `showcase/repulsive-surfaces/`). Tests that need the reference scenes
//! (`scripts/rs-fetch-scenes.sh`) are skipped when `data/repulsive-surfaces/scenes` is missing.

#[path = "../examples/apps/apps.rs"]
mod apps;

use apps::{Cfg, Outcome};

fn cfg() -> Cfg {
    Cfg::smoke(2)
}

fn have_scenes() -> bool {
    let ok = apps::scenes_dir().join("FORMAT.txt").exists();
    if !ok {
        eprintln!("skipped: reference scenes not fetched");
    }
    ok
}

fn check(o: &Outcome) {
    println!("{}", o.line());
    assert!(!o.status.starts_with("Error"), "{}: {}", o.name, o.status);
    assert!(o.e1.is_finite() && o.e0.is_finite(), "{}: energy {} -> {}", o.name, o.e0, o.e1);
    assert!(o.steps > 0, "{}: no step", o.name);
    assert!(o.max_violation < 1e-3, "{}: constraint violation {}", o.name, o.max_violation);
}

#[test]
fn genus_surfaces() {
    for g in [0, 2] {
        let o = apps::genus(&cfg(), g, None).unwrap();
        check(&o);
        assert!((o.area1 / o.area0 - 1.0).abs() < 1e-3);
        assert_eq!(o.intersects, Some(false));
        assert!(o.extra.contains(&format!("genus={g}")), "{}", o.extra);
    }
}

#[test]
fn genus_with_attracting_plane() {
    check(&apps::genus(&cfg(), 2, Some(5.0)).unwrap());
}

#[test]
fn inertia_indicators_of_a_sphere_and_a_torus() {
    let s = ddg_repulsive_surfaces::generate::primitives::sphere(3, 1.0, [0.0; 3]).unwrap();
    let (_, mult, _) = apps::inertia_indicators(&s.faces, &s.positions);
    assert_eq!(mult, 3);
    let t = ddg_repulsive_surfaces::generate::primitives::torus(48, 24, 2.0, 0.7).unwrap();
    let (e, mult, _) = apps::inertia_indicators(&t.faces, &t.positions);
    assert_eq!(mult, 2);
    assert!(e[2] / e[0] > 1.5 || e[0] / e[2] < 0.7);
}

/// `eig3` on non-diagonal input (the torus test above is already diagonal): Q diag(1, 2, 5) Qᵀ
/// for a generic rotation Q, and a tilted torus keeps its double moment. A swapped `atan2`
/// in the Jacobi angle left the off-diagonal entries unreduced (NOTES-FIDELITY.md).
#[test]
fn eig3_and_inertia_of_rotated_input() {
    let (c1, s1, c2, s2) = (0.3f64.cos(), 0.3f64.sin(), 1.1f64.cos(), 1.1f64.sin());
    let rx = [[1.0, 0.0, 0.0], [0.0, c1, -s1], [0.0, s1, c1]];
    let rz = [[c2, -s2, 0.0], [s2, c2, 0.0], [0.0, 0.0, 1.0]];
    let mut q = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                q[i][j] += rz[i][k] * rx[k][j];
            }
        }
    }
    let d = [1.0, 2.0, 5.0];
    let mut a = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                a[i][j] += q[i][k] * d[k] * q[j][k];
            }
        }
    }
    let e = apps::eig3(a);
    for k in 0..3 {
        assert!((e[k] - d[k]).abs() < 1e-12, "{e:?}");
    }
    let t = ddg_repulsive_surfaces::generate::primitives::torus(48, 24, 2.0, 0.7).unwrap();
    let x: Vec<[f64; 3]> = t
        .positions
        .iter()
        .map(|p| [0, 1, 2].map(|i| q[i][0] * p[0] + q[i][1] * p[1] + q[i][2] * p[2]))
        .collect();
    let (_, mult, desc) = apps::inertia_indicators(&t.faces, &x);
    assert_eq!(mult, 2, "{desc}");
}

#[test]
fn handcuffs_linked_and_unlinked() {
    let v = apps::handcuffs(&cfg()).unwrap();
    assert_eq!(v.len(), 2);
    for o in &v {
        check(o);
        assert!(o.extra.contains("genus=2"), "{}", o.extra);
    }
}

#[test]
fn punctured_torus_eversion() {
    let o = apps::eversion(&cfg()).unwrap();
    check(&o);
    assert!(o.extra.contains("V_mid"));
}

#[test]
fn h_shape_variants() {
    for v in apps::HVariant::ALL {
        let o = apps::h_shape(&cfg(), v).unwrap();
        check(&o);
        assert!(o.extra.contains("max_pin_drift"));
    }
}

#[test]
fn subcritical_torus_and_klein() {
    for o in apps::subcritical_small(&cfg()).unwrap() {
        check(&o);
    }
}

#[test]
fn scene_applications() {
    if !have_scenes() {
        return;
    }
    check(&apps::walnut(&cfg()).unwrap());
    check(&apps::klein(&cfg()).unwrap());
    check(&apps::cobblestones(&cfg()).unwrap());
    let env = apps::nested_envelopes(&cfg(), 2).unwrap();
    check(&env[0]);
    // Level 2 starts on the coarsened level 1 in contact with it: the projection guard
    // (`LineSearch::guard_projection`) may keep the volume behind its target rather than push
    // the envelope into its obstacle; the violation is reported, not hidden.
    println!("{}", env[1].line());
    assert!(!env[1].status.starts_with("Error"), "{}", env[1].status);
    assert!(env[1].max_violation < 0.1, "{}", env[1].max_violation);
    assert_eq!(env[1].intersects, Some(false));
}

#[test]
#[ignore]
fn heavy_scene_applications() {
    // 10k to 20k faces: a few steps each, run on demand (`cargo test -- --ignored`)
    if !have_scenes() {
        return;
    }
    for o in apps::handcuffs_pole(&cfg()).unwrap() {
        check(&o);
    }
    for o in apps::slide_off_bar(&cfg()).unwrap() {
        check(&o);
    }
    check(&apps::shrink_wrap_points(&cfg(), "bunny").unwrap());
    check(&apps::shrink_wrap_soup(&cfg()).unwrap());
}
