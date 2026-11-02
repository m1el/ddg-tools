//! Constraints (`SPEC.md` §5.1, §13.2): values, Jacobians by finite differences, the exact
//! positional projection, the rank check and its drops (D10), the translation check, errors,
//! and agreement with the oracle's `constraints` fixtures.

#![allow(clippy::needless_range_loop)]

#[path = "wp2_fixture/mod.rs"]
mod fixture;

use ddg_core::generators::{bumpy_grid_disk, icosphere, perturb, torus};
use ddg_repulsive_surfaces::constraints::{BarycenterMode, Constraint, ConstraintSet, RowKind};
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::{RsError, Surface};

fn two_spheres() -> (Surface, Vec<[f64; 3]>) {
    let m = icosphere(1).unwrap();
    let mut p = m.positions.clone();
    let mut f = m.faces.clone();
    let off = p.len() as u32;
    p.extend(
        m.positions
            .iter()
            .map(|q| [0.7 * q[0] + 2.2, 0.7 * q[1] + 0.3, 0.7 * q[2]]),
    );
    f.extend(m.faces.iter().map(|t| [t[0] + off, t[1] + off, t[2] + off]));
    let s = Surface::new(p.len(), f).unwrap();
    let l = mean_edge_length(&s, &p);
    perturb(&mut p, [0.1 * l; 3], 11);
    (s, p)
}

fn perturbed(m: ddg_core::generators::GeneratedMesh, seed: u64) -> (Surface, Vec<[f64; 3]>) {
    let mut p = m.positions;
    let s = Surface::new(p.len(), m.faces).unwrap();
    let l = mean_edge_length(&s, &p);
    perturb(&mut p, [0.1 * l; 3], seed);
    (s, p)
}

/// Central differences of the area and volume rows against the Jacobian.
#[test]
fn jacobians_match_finite_differences() {
    let meshes = vec![
        ("icosphere", perturbed(icosphere(1).unwrap(), 1)),
        ("torus", perturbed(torus(8, 5, 1.5, 0.6).unwrap(), 2)),
        ("disk", perturbed(bumpy_grid_disk(4, 4, 0.3).unwrap(), 3)),
        ("two", two_spheres()),
    ];
    for (name, (s, p)) in meshes {
        let mut cons = vec![
            Constraint::Barycenter {
                component: None,
                target: None,
            },
            Constraint::Area {
                component: None,
                target: Some(1.0),
            },
            Constraint::Volume {
                component: None,
                target: Some(0.5),
            },
        ];
        if s.n_components() > 1 {
            cons.push(Constraint::Area {
                component: Some(1),
                target: None,
            });
            cons.push(Constraint::Volume {
                component: Some(0),
                target: None,
            });
        }
        let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
        assert!(set.dropped_rows().is_empty(), "{name}: {:?}", set.dropped_rows());
        let rows = set.nonlinear_rows();
        let j = set.jacobian(&p, &rows).unwrap();
        let jd = {
            let mut d = vec![vec![0.0; 3 * s.n_vertices()]; rows.len()];
            for (r, c, x) in j.triplets() {
                d[r as usize][c as usize] += x;
            }
            d
        };
        let h = 1e-6 * mean_edge_length(&s, &p);
        let mut worst: f64 = 0.0;
        for v in 0..s.n_vertices() {
            for k in 0..3 {
                let mut pp = p.clone();
                pp[v][k] += h;
                let mut pm = p.clone();
                pm[v][k] -= h;
                let fp = set.values(&pp).unwrap();
                let fm = set.values(&pm).unwrap();
                for (i, &r) in rows.iter().enumerate() {
                    let fd = (fp[r] - fm[r]) / (2.0 * h);
                    let scale = jd[i].iter().fold(0.0f64, |m, x| m.max(x.abs()));
                    worst = worst.max((fd - jd[i][3 * v + k]).abs() / scale);
                }
            }
        }
        eprintln!("{name}: area/volume Jacobian vs central differences: {worst:.2e}");
        assert!(worst < 1e-6, "{name}: {worst:e}");
        // The barycenter rows: frozen weights a_v / A_c; a translation moves Φ by exactly t.
        let pos_rows = set.positional_rows();
        let jq = set.jacobian(&p, &pos_rows).unwrap();
        let t = [0.3, -0.2, 0.7];
        let tv: Vec<f64> = (0..s.n_vertices()).flat_map(|_| t).collect();
        let mut jt = vec![0.0; pos_rows.len()];
        jq.mul_vec(&tv, &mut jt).unwrap();
        let shifted: Vec<[f64; 3]> = p
            .iter()
            .map(|q| [q[0] + t[0], q[1] + t[1], q[2] + t[2]])
            .collect();
        let f0 = set.values(&p).unwrap();
        let f1 = set.values(&shifted).unwrap();
        for (i, &r) in pos_rows.iter().enumerate() {
            let RowKind::Barycenter { axis, .. } = set.rows()[r].kind else {
                panic!()
            };
            assert!((jt[i] - t[axis as usize]).abs() < 1e-14);
            assert!((f1[r] - f0[r] - t[axis as usize]).abs() < 1e-12);
        }
    }
}

#[test]
fn volume_gradient_is_a_third_of_the_area_weighted_normal_on_closed_surfaces() {
    let (s, p) = perturbed(icosphere(2).unwrap(), 4);
    let set = ConstraintSet::new(
        &s,
        &p,
        &[Constraint::Barycenter {
            component: None,
            target: None,
        }],
        BarycenterMode::PerComponent,
    )
    .unwrap();
    let g = set.volume_gradient(&p, None);
    let mut an = vec![[0.0; 3]; s.n_vertices()];
    for t in s.faces() {
        let a = p[t[0] as usize];
        let b = p[t[1] as usize];
        let c = p[t[2] as usize];
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        // a_S N_S = ½ e1 × e2.
        let n = [
            0.5 * (e1[1] * e2[2] - e1[2] * e2[1]),
            0.5 * (e1[2] * e2[0] - e1[0] * e2[2]),
            0.5 * (e1[0] * e2[1] - e1[1] * e2[0]),
        ];
        for &v in t {
            for k in 0..3 {
                an[v as usize][k] += n[k] / 3.0;
            }
        }
    }
    for v in 0..s.n_vertices() {
        for k in 0..3 {
            assert!((g[v][k] - an[v][k]).abs() < 1e-14, "vertex {v}");
        }
    }
}

#[test]
fn positional_projection_is_exact() {
    let (s, p) = two_spheres();
    let cons = vec![
        Constraint::Barycenter {
            component: Some(0),
            target: Some([0.1, 0.2, -0.3]),
        },
        Constraint::Pins {
            vertices: vec![50, 60],
            targets: Some(vec![[2.0, 0.5, 0.1], [2.5, 0.0, 0.2]]),
        },
        Constraint::Barycenter {
            component: Some(1),
            target: None,
        },
        Constraint::Area {
            component: None,
            target: None,
        },
    ];
    let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
    let mut q = p.clone();
    perturb(&mut q, [0.05; 3], 99);
    set.project_positional(&mut q).unwrap();
    let phi = set.values(&q).unwrap();
    // Component 1 has pins and a barycenter: translation and reset are repeated, both hold.
    for (r, row) in set.rows().iter().enumerate() {
        match row.kind {
            RowKind::Pin { .. } => assert_eq!(phi[r], 0.0, "pin row {r}"),
            RowKind::Barycenter { .. } => {
                assert!(phi[r].abs() < 1e-13, "barycenter row {r}: {:e}", phi[r])
            }
            _ => {}
        }
    }
}

#[test]
fn rank_check_drops_dependent_rows() {
    // A vertex pinned twice: the second pin is dropped.
    let (s, p) = perturbed(icosphere(1).unwrap(), 5);
    let cons = vec![
        Constraint::Pins {
            vertices: vec![3, 4],
            targets: None,
        },
        Constraint::Pins {
            vertices: vec![4],
            targets: None,
        },
        Constraint::Area {
            component: None,
            target: None,
        },
        Constraint::Area {
            component: Some(0),
            target: None,
        },
        Constraint::Volume {
            component: None,
            target: None,
        },
    ];
    let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
    // Rows: pins 0..6, pin 6..9, area 9, area 10, volume 11.
    assert_eq!(set.dropped_rows(), &[6, 7, 8, 10]);
    assert_eq!(set.positional_rows(), (0..6).collect::<Vec<_>>());
    assert_eq!(set.nonlinear_rows(), vec![9, 11]);
    eprintln!("rank report: {:?}", set.rank_report());
    // Pinned twice to different targets: infeasible.
    let bad = vec![
        Constraint::Pins {
            vertices: vec![3],
            targets: None,
        },
        Constraint::Pins {
            vertices: vec![3],
            targets: Some(vec![[9.0, 9.0, 9.0]]),
        },
    ];
    assert!(matches!(
        ConstraintSet::new(&s, &p, &bad, BarycenterMode::PerComponent),
        Err(RsError::InfeasibleConstraint { constraint: 1, .. })
    ));
    // A barycenter of a component whose vertices are all pinned is dropped; an area row of a
    // fully pinned mesh is in the span of the pins and dropped too.
    let all: Vec<u32> = (0..s.n_vertices() as u32).collect();
    let cons = vec![
        Constraint::Barycenter {
            component: None,
            target: None,
        },
        Constraint::Pins {
            vertices: all,
            targets: None,
        },
        Constraint::Area {
            component: None,
            target: None,
        },
    ];
    let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
    assert_eq!(set.dropped_rows(), &[0, 1, 2, 3 + 3 * s.n_vertices()]);
    // A barycenter with some pins is independent and kept.
    let cons = vec![
        Constraint::Barycenter {
            component: None,
            target: None,
        },
        Constraint::Pins {
            vertices: vec![0, 1, 2],
            targets: None,
        },
    ];
    let set = ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent).unwrap();
    assert!(set.dropped_rows().is_empty());
}

#[test]
fn translation_and_input_errors() {
    let (s, p) = two_spheres();
    // Component 1 has neither a barycenter nor a pin.
    let cons = vec![
        Constraint::Barycenter {
            component: Some(0),
            target: None,
        },
        Constraint::Area {
            component: None,
            target: None,
        },
    ];
    assert!(matches!(
        ConstraintSet::new(&s, &p, &cons, BarycenterMode::PerComponent),
        Err(RsError::TranslationNotFixed { component: 1 })
    ));
    // A pin on component 1 fixes it.
    let mut ok = cons.clone();
    ok.push(Constraint::Pins {
        vertices: vec![45],
        targets: None,
    });
    ConstraintSet::new(&s, &p, &ok, BarycenterMode::PerComponent).unwrap();
    // Global mode: one barycenter group, accepted for every component.
    let g = vec![Constraint::Barycenter {
        component: None,
        target: None,
    }];
    let set = ConstraintSet::new(&s, &p, &g, BarycenterMode::Global).unwrap();
    assert_eq!(set.n_rows(), 3);
    let set = ConstraintSet::new(&s, &p, &g, BarycenterMode::PerComponent).unwrap();
    assert_eq!(set.n_rows(), 6);
    // Out of range.
    let bad = vec![Constraint::Area {
        component: Some(2),
        target: None,
    }];
    assert!(matches!(
        ConstraintSet::new(&s, &p, &bad, BarycenterMode::PerComponent),
        Err(RsError::ConstraintOutOfRange {
            constraint: 0,
            index: 2
        })
    ));
    let bad = vec![Constraint::Pins {
        vertices: vec![1000],
        targets: None,
    }];
    assert!(matches!(
        ConstraintSet::new(&s, &p, &bad, BarycenterMode::PerComponent),
        Err(RsError::ConstraintOutOfRange {
            constraint: 0,
            index: 1000
        })
    ));
    // Volume on an inconsistently ordered (orientable) component: refused, naming an edge.
    let m = icosphere(1).unwrap();
    let mut f = m.faces.clone();
    f[5] = [f[5][0], f[5][2], f[5][1]];
    let s2 = Surface::new(m.positions.len(), f).unwrap();
    let vol = vec![
        Constraint::Barycenter {
            component: None,
            target: None,
        },
        Constraint::Volume {
            component: None,
            target: None,
        },
    ];
    match ConstraintSet::new(&s2, &m.positions, &vol, BarycenterMode::PerComponent) {
        Err(RsError::VolumeNeedsOrientation { component: 0, edge }) => {
            let t = s2.faces()[5];
            assert!(t.contains(&edge[0]) && t.contains(&edge[1]));
        }
        other => panic!("{other:?}"),
    }
    // Non-finite target.
    let bad = vec![Constraint::Area {
        component: None,
        target: Some(f64::NAN),
    }];
    assert!(matches!(
        ConstraintSet::new(&s, &p, &bad, BarycenterMode::PerComponent),
        Err(RsError::InvalidOption { .. })
    ));
}

// ---- the oracle's golden fixtures (relative 1e-9) ----

#[test]
fn oracle_constraint_fixtures() {
    let names = fixture::cases("constraints");
    assert!(names.len() >= 6);
    for name in names {
        let c = fixture::read_case(&name);
        let o = fixture::read_output(&name);
        let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
        let ncomp: usize = o.get("mesh")[0][2].parse().unwrap();
        assert_eq!(s.n_components(), ncomp);
        let vcomp: Vec<u32> = o.get("vcomp")[0].iter().map(|x| x.parse().unwrap()).collect();
        assert_eq!(s.vertex_components(), vcomp.as_slice(), "{name}: component numbering");
        // The oracle's translation check is the saddle's; constraint fixtures may leave a
        // component free, so add a pin on every unfixed component only if needed.
        let set = match ConstraintSet::new(&s, &c.positions, &c.constraints, Default::default()) {
            Ok(set) => set,
            Err(RsError::TranslationNotFixed { .. }) => {
                let mut cons = c.constraints.clone();
                for comp in s.components() {
                    cons.push(Constraint::Pins {
                        vertices: vec![comp.vertices()[0]],
                        targets: None,
                    });
                }
                let set =
                    ConstraintSet::new(&s, &c.positions, &cons, Default::default()).unwrap();
                eprintln!("{name}: pins added for the translation check");
                set
            }
            Err(e) => panic!("{name}: {e}"),
        };
        let rows_out = o.get("row");
        let nr = rows_out.len();
        assert!(set.dropped_rows().iter().all(|&r| r >= nr), "{name}: drops");
        let phi = set.values(&c.positions).unwrap();
        let rows: Vec<usize> = (0..nr).collect();
        let j = set.jacobian(&c.positions, &rows).unwrap();
        let mut jd = vec![0.0; nr * 3 * s.n_vertices()];
        for (r, col, x) in j.triplets() {
            jd[r as usize * 3 * s.n_vertices() + col as usize] += x;
        }
        let mut want_phi = Vec::new();
        for (i, l) in rows_out.iter().enumerate() {
            assert_eq!(l[0].parse::<usize>().unwrap(), i);
            assert_eq!(l[1], set.rows()[i].kind.name(), "{name} row {i}");
            assert_eq!(
                l[2] == "positional",
                set.rows()[i].kind.is_positional(),
                "{name} row {i}"
            );
            want_phi.push(fixture::hex(&l[3]));
        }
        let want_j: Vec<f64> = o.matrix("J").into_iter().flatten().collect();
        let e_phi = fixture::rel_err(&phi[..nr], &want_phi, 1e-300);
        let e_j = fixture::rel_err(&jd, &want_j, 1e-300);
        eprintln!("{name}: phi {e_phi:.2e}, J {e_j:.2e}");
        // Φ is a difference (value − target): compare relative to the values' scale.
        let scale = set.rows()[..nr]
            .iter()
            .fold(1.0f64, |m, r| m.max(r.target.abs()));
        let e_phi_abs = phi[..nr]
            .iter()
            .zip(&want_phi)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
            / scale;
        assert!(e_phi_abs < 1e-9, "{name}: phi {e_phi_abs:e}");
        assert!(e_j < 1e-9, "{name}: J {e_j:e}");
    }
}

#[test]
fn oracle_constraint_error_fixtures() {
    for (name, expect) in [
        ("rs-error-volume-mobius", "VolumeNeedsOrientation"),
        ("rs-error-pin-range", "ConstraintOutOfRange"),
        ("rs-error-component", "ConstraintOutOfRange"),
        ("rs-error-translation", "TranslationNotFixed"),
    ] {
        let c = fixture::read_case(name);
        let o = fixture::read_output(name);
        assert!(o.status.starts_with("error"), "{name}");
        let s = Surface::new(c.positions.len(), c.faces.clone()).unwrap();
        let e = ConstraintSet::new(&s, &c.positions, &c.constraints, Default::default())
            .expect_err(name);
        let got = format!("{e:?}");
        assert!(got.starts_with(expect), "{name}: {got}");
    }
}
