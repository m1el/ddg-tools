//! Tests of the flow driver (SPEC §5.1, §10.7). The fast tests run by default; the tests on the
//! meshes of `data/models/` are `#[ignore]`d (they need `--release` to be quick):
//!
//! ```text
//! cargo test --release -p ddg-developability --test flow -- --ignored --nocapture
//! ```
//!
//! (`data/models/` is unpacked from `/var/data/discrete-geometry/models/*.zip`; tests skip with a
//! message if a mesh is missing.)

mod energy_common;
use energy_common::*;

use std::f64::consts::PI;

use ddg_core::TriMesh;
use ddg_core::generators::{grid_disk, icosphere, perturb};
use ddg_core::topology::NONE;
use ddg_developability::geometry::corner_angles;
use ddg_developability::{
    Control, DevError, Energy, EnergyOptions, FlowOptions, Round, evaluate_energy, flow, RemeshOptions,
};

fn energy_of(mesh: &TriMesh, pos: &[[f64; 3]], e: Energy) -> f64 {
    evaluate_energy(mesh, pos, &EnergyOptions::new(e)).unwrap().0
}

fn n_bad(mesh: &TriMesh, pos: &[[f64; 3]], small: f64) -> usize {
    mesh.faces()
        .iter()
        .filter(|f| corner_angles(f.map(|v| pos[v as usize])).iter().any(|&a| a < small))
        .count()
}

/// A developable sheet: the grid mapped to a cylinder arc (`cone = false`) or a cone, with noise.
fn developable_sheet(n: usize, cone: bool, noise: f64, seed: u64) -> Mesh {
    let g = grid_disk(n, n).unwrap();
    let mut pos: Vec<[f64; 3]> = g
        .positions
        .iter()
        .map(|p| {
            let (x, y) = (p[0], p[1]);
            if cone {
                let s = 1.0 + y; // distance from the apex along the cone
                let (beta, phi) = (0.6f64, 1.5 * x);
                [s * beta.sin() * phi.cos(), s * beta.sin() * phi.sin(), s * beta.cos()]
            } else {
                let r = 0.5;
                [r * (2.0 * x).cos(), r * (2.0 * x).sin(), y]
            }
        })
        .collect();
    perturb(&mut pos, [noise; 3], seed);
    (g.topology().unwrap(), pos)
}

fn no_remesh(opts: &mut FlowOptions) {
    opts.remesh = None;
}

#[test]
fn noisy_developable_sheets_become_developable() {
    for (cone, seed) in [(false, 1u64), (true, 2)] {
        let n = 10;
        let (mesh, pos) = developable_sheet(n, cone, 0.15 / n as f64, seed);
        let e0 = energy_of(&mesh, &pos, Energy::Covariance);
        let opts = FlowOptions::uniform(EnergyOptions::default(), 0, 600);
        let res = flow(&mesh, &pos, &opts, None).unwrap();
        let e1 = energy_of(&res.mesh, &res.positions, Energy::Covariance);
        eprintln!(
            "sheet cone={cone}: E {e0:.3e} -> {e1:.3e}, iterations {}, status {:?}",
            res.report.rounds[0].optimize.iterations, res.report.rounds[0].optimize.status
        );
        assert!(e1 < 1e-3 * e0, "cone={cone}: energy {e0} -> {e1}");
        assert!(res.report.rounds[0].optimize.energy_final.is_finite());
    }
}

/// D34: a flow run until the line search cannot improve the energy (no tolerance, a large
/// budget) ends `Stationary` at an energy at rounding level, not `LineSearchFailed`.
#[test]
fn converged_flow_is_stationary() {
    let (mesh, pos) = developable_sheet(8, false, 0.02, 3);
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 0, 20_000);
    opts.optimize.gradient_tolerance = 0.0;
    opts.optimize.energy_tolerance = 0.0;
    no_remesh(&mut opts);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    let o = &res.report.rounds[0].optimize;
    eprintln!("converged sheet: E {:.3e} -> {:.3e}, {} iterations, status {:?}, |g| {:.2e}", o.energy_initial, o.energy_final, o.iterations, o.status, o.gradient_norm_history.last().unwrap());
    assert_eq!(o.status, ddg_developability::OptimizeStatus::Stationary);
    assert!(o.energy_final < 1e-12, "E = {}", o.energy_final);
}

/// D36: the default flow is the paper's recipe of Figs. 23-24 (E_P coarse, E_lambda fine).
#[test]
fn default_flow_is_coarse_to_fine() {
    let o = FlowOptions::default();
    assert_eq!(o, FlowOptions::coarse_to_fine(2, 500));
    let e: Vec<(Energy, bool, usize)> = o.rounds.iter().map(|r| (r.energy.energy, r.subdivide_before, r.max_iterations)).collect();
    assert_eq!(e, vec![(Energy::Combinatorial, false, 500), (Energy::Covariance, true, 500), (Energy::Covariance, true, 500)]);
    assert!(o.remesh.is_some() && o.pinned.is_empty());
}

#[test]
fn sphere_energy_concentrates() {
    let g = icosphere(2).unwrap();
    let mesh = g.topology().unwrap();
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 0, 400);
    no_remesh(&mut opts);
    let res = flow(&mesh, &g.positions, &opts, None).unwrap();
    let r = &res.report.rounds[0].optimize;
    assert!(r.energy_final < r.energy_initial);
    let share = |mesh: &TriMesh, pos: &[[f64; 3]]| {
        let (e, ve) = evaluate_energy(mesh, pos, &EnergyOptions::new(Energy::Covariance)).unwrap();
        let mut s = ve;
        s.sort_by(|a, b| b.total_cmp(a));
        let top: f64 = s.iter().take(s.len() / 5).sum();
        (e, top / s.iter().sum::<f64>())
    };
    let (_, before) = share(&mesh, &g.positions);
    let (e, after) = share(&res.mesh, &res.positions);
    eprintln!("sphere: E {:.3e} -> {:.3e}; top 20% of vertices hold {:.1}% (before {:.1}%)", r.energy_initial, e, 100.0 * after, 100.0 * before);
    assert!(after > 0.6 && after > before + 0.2, "energy should concentrate on few vertices (seams, cones)");
}

#[test]
fn subdivision_preserves_energy_and_rounds_decrease() {
    let (mesh, pos) = developable_sheet(5, false, 0.03, 3);
    for energy in Energy::ALL {
        let mut opts = FlowOptions::uniform(EnergyOptions::new(energy), 2, 40);
        no_remesh(&mut opts);
        let res = flow(&mesh, &pos, &opts, None).unwrap();
        let rr = &res.report.rounds;
        assert_eq!(rr.len(), 3);
        for r in rr {
            // monotone within a round (strict decrease of accepted steps, L&O with c1 = 0)
            let h = &r.optimize.energy_history;
            assert!(h.windows(2).all(|w| w[1] <= w[0] * (1.0 + 1e-12) + 1e-300), "{energy:?}: {h:?}");
            assert!(r.optimize.energy_final <= r.optimize.energy_initial);
        }
        for w in rr.windows(2) {
            let (a, b) = (w[0].optimize.energy_final, w[1].energy_after_subdivision);
            assert!((a - b).abs() <= 1e-9 * a.abs().max(1e-6), "{energy:?}: subdivision changed the energy {a} -> {b}");
        }
        assert_eq!(rr[2].n_faces, 16 * mesh.n_faces());
    }
}

#[test]
fn vertex_source_and_pinned_vertices() {
    let (mesh, pos) = developable_sheet(6, false, 0.03, 4);
    let pinned: Vec<u32> = (0..mesh.n_vertices() as u32).filter(|v| v % 5 == 0).collect();
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 1, 30);
    opts.remesh = None;
    opts.pinned = pinned.clone();
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    let nv0 = mesh.n_vertices();
    assert_eq!(res.vertex_source.len(), res.mesh.n_vertices());
    for v in 0..res.mesh.n_vertices() {
        if v < nv0 {
            assert_eq!(res.vertex_source[v], v as u32);
        } else {
            assert_eq!(res.vertex_source[v], NONE);
        }
    }
    for &p in &pinned {
        assert_eq!(res.positions[p as usize], pos[p as usize], "pinned vertex {p} moved");
    }
    assert!(res.pinned.len() >= pinned.len());
    for &p in &res.pinned {
        let s = res.vertex_source[p as usize];
        // a pinned midpoint is created between two pinned vertices and stays at their midpoint
        if s == NONE {
            continue;
        }
        assert!(pinned.contains(&s));
    }
    assert_ne!(res.positions, pos[..].to_vec());
}

#[test]
fn remeshing_in_a_flow_removes_needles_and_tracks_sources() {
    // a needle: pull a vertex of a sphere almost onto a neighbour
    let g = icosphere(2).unwrap();
    let mesh = g.topology().unwrap();
    let mut pos = g.positions.clone();
    let nb = mesh.vertex_neighbors(7).next().unwrap();
    let target = pos[nb];
    for (x, t) in pos[7].iter_mut().zip(target) {
        *x = t + 1e-3 * (*x - t);
    }
    let small = 1f64.to_radians();
    assert!(n_bad(&mesh, &pos, small) > 0);
    let opts = FlowOptions::uniform(EnergyOptions::default(), 0, 30);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    let rr = &res.report.rounds[0];
    assert!(rr.remesh[0].collapses + rr.remesh[0].flips > 0);
    assert!(res.mesh.n_vertices() < mesh.n_vertices());
    assert_eq!(res.vertex_source.len(), res.mesh.n_vertices());
    // sources are distinct input ids
    let mut s = res.vertex_source.clone();
    s.sort_unstable();
    s.dedup();
    assert_eq!(s.len(), res.mesh.n_vertices());
    assert!(rr.remesh[0].energy_before.is_some() && rr.remesh[0].energy_after.is_some());
    assert_eq!(n_bad(&res.mesh, &res.positions, small), 0);
    assert_eq!(res.mesh.euler_characteristic(), 2);
}

#[test]
fn cancel_from_the_progress_callback() {
    let (mesh, pos) = developable_sheet(8, false, 0.05, 5);
    let opts = FlowOptions::uniform(EnergyOptions::default(), 1, 100);
    let mut calls = 0usize;
    let mut cb = |p: &ddg_developability::Progress| {
        calls += 1;
        assert_eq!(p.round, 0);
        assert_eq!(p.iteration, calls);
        if calls == 7 { Control::Break } else { Control::Continue }
    };
    let res = flow(&mesh, &pos, &opts, Some(&mut cb)).unwrap();
    assert!(res.report.cancelled);
    assert_eq!(res.report.rounds.len(), 1);
    assert_eq!(res.report.rounds[0].optimize.iterations, 7);
    assert_eq!(calls, 7);
}

#[test]
fn result_does_not_depend_on_threads() {
    let (mesh, pos) = developable_sheet(8, true, 0.05, 6);
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 1, 25);
    let base = flow(&mesh, &pos, &opts, None).unwrap();
    let mut iters = 0;
    for t in [1usize, 3] {
        opts.threads = Some(t);
        let mut cb = |_: &ddg_developability::Progress| {
            iters += 1;
            Control::Continue
        };
        let r = flow(&mesh, &pos, &opts, Some(&mut cb)).unwrap();
        assert_eq!(r.positions, base.positions, "threads = {t}");
        assert_eq!(r.mesh, base.mesh);
    }
    assert!(iters > 0);
}

#[test]
fn invalid_input_is_an_error() {
    let (mesh, pos) = developable_sheet(4, false, 0.0, 7);
    let ok = FlowOptions::uniform(EnergyOptions::default(), 0, 2);
    let mut o = ok.clone();
    o.rounds.clear();
    assert!(matches!(flow(&mesh, &pos, &o, None), Err(DevError::InvalidOption { name: "rounds", .. })));
    let mut o = ok.clone();
    o.optimize.pinned = vec![0];
    assert!(matches!(flow(&mesh, &pos, &o, None), Err(DevError::InvalidOption { .. })));
    let mut o = ok.clone();
    o.pinned = vec![10_000];
    assert!(matches!(flow(&mesh, &pos, &o, None), Err(DevError::VertexOutOfRange { .. })));
    let mut o = ok.clone();
    o.threads = Some(0);
    assert!(matches!(flow(&mesh, &pos, &o, None), Err(DevError::InvalidOption { .. })));
    let mut o = ok.clone();
    o.remesh = Some(RemeshOptions { small_angle: 2.0, ..Default::default() });
    assert!(matches!(flow(&mesh, &pos, &o, None), Err(DevError::InvalidOption { .. })));
    assert!(matches!(flow(&mesh, &pos[1..], &ok, None), Err(DevError::LengthMismatch { .. })));
    let mut bad = pos.clone();
    bad[3][1] = f64::NAN;
    assert!(matches!(flow(&mesh, &bad, &ok, None), Err(DevError::NonFinitePosition { vertex: 3 })));
    // zero iterations returns the input
    let mut z = ok.clone();
    z.rounds[0].max_iterations = 0;
    z.remesh = None;
    let r = flow(&mesh, &pos, &z, None).unwrap();
    assert_eq!(r.positions, pos);
}

#[test]
fn coarse_then_fine_rounds() {
    // E_P on the coarse round, E_lambda after a subdivision (as for Figs. 23-24)
    let (mesh, pos) = developable_sheet(5, true, 0.04, 8);
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 1, 60);
    opts.rounds = vec![
        Round { energy: EnergyOptions::new(Energy::Combinatorial), max_iterations: 60, subdivide_before: false },
        Round { energy: EnergyOptions::new(Energy::Covariance), max_iterations: 60, subdivide_before: true },
    ];
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    assert_eq!(res.report.rounds.len(), 2);
    assert!(res.report.rounds[1].n_faces >= 3 * mesh.n_faces());
    for r in &res.report.rounds {
        assert!(r.optimize.energy_final < r.optimize.energy_initial);
    }
}

// ---------------------------------------------------------------------------------------------
// Meshes of data/models (long; run with --release --ignored)

fn print_round_table(name: &str, res: &ddg_developability::FlowResult) {
    for (i, r) in res.report.rounds.iter().enumerate() {
        eprintln!(
            "{name} round {i}: {} f, E {:.4e} -> {:.4e}, {} it ({:?}), segments {}, remesh flips/collapses {}/{}, {:.2}s",
            r.n_faces,
            r.optimize.energy_initial,
            r.optimize.energy_final,
            r.optimize.iterations,
            r.optimize.status,
            r.segments.len(),
            r.remesh.iter().map(|m| m.flips).sum::<usize>(),
            r.remesh.iter().map(|m| m.collapses).sum::<usize>(),
            r.seconds_optimize + r.seconds_remesh + r.seconds_subdivide
        );
    }
}

fn check_rounds(res: &ddg_developability::FlowResult) {
    for r in &res.report.rounds {
        let o = &r.optimize;
        assert!(o.energy_final.is_finite());
        assert!(o.energy_final < o.energy_initial, "no decrease in a round: {} -> {}", o.energy_initial, o.energy_final);
        for s in &r.segments {
            assert!(s.energy_history.windows(2).all(|w| w[1] <= w[0] * (1.0 + 1e-12)));
        }
    }
    assert_eq!(res.vertex_source.len(), res.mesh.n_vertices());
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test flow -- --ignored"]
fn crumpled_developable_energy_falls_by_90_percent() {
    let Some((mesh, pos)) = load_model("crumpleddevelopable/CrumpledDevelopable.obj") else { return };
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 0, 500);
    opts.remesh = None; // topology fixed: dihedral angles comparable
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    print_round_table("crumpled", &res);
    let r = &res.report.rounds[0].optimize;
    eprintln!("crumpled: displacement {:.3} lbar", r.max_vertex_displacement);
    assert!(r.energy_final <= 0.1 * r.energy_initial);
    check_rounds(&res);
    let (m0, m1) = (mean_abs_defect(&mesh, &pos), mean_abs_defect(&res.mesh, &res.positions));
    eprintln!("crumpled: mean |angle defect| {m0:.4} -> {m1:.4}");
    let (h0, h1) = (creases(&mesh, &pos), creases(&res.mesh, &res.positions));
    eprintln!("crumpled: dihedral > 30deg: {h0} -> {h1} edges of {}", mesh.n_edges());
    assert!(h1 <= h0);
}

fn mean_abs_defect(mesh: &TriMesh, pos: &[[f64; 3]]) -> f64 {
    let mut sum = vec![0.0; mesh.n_vertices()];
    for f in mesh.faces() {
        let a = corner_angles(f.map(|v| pos[v as usize]));
        for k in 0..3 {
            sum[f[k] as usize] += a[k];
        }
    }
    let int = mesh.interior_vertices();
    int.iter().map(|&v| (2.0 * PI - sum[v as usize]).abs()).sum::<f64>() / int.len() as f64
}

fn creases(mesh: &TriMesh, pos: &[[f64; 3]]) -> usize {
    let n: Vec<[f64; 3]> = mesh
        .faces()
        .iter()
        .map(|f| {
            let p = f.map(|v| pos[v as usize]);
            normalize(cross(sub(p[1], p[0]), sub(p[2], p[0])))
        })
        .collect();
    (0..mesh.n_edges())
        .filter(|&e| {
            if mesh.is_boundary_edge(e) {
                return false;
            }
            let h = mesh.edge_halfedge(e);
            let (a, b) = (mesh.face(h), mesh.face(mesh.twin(h).unwrap()));
            dot(n[a], n[b]).clamp(-1.0, 1.0).acos() > 30f64.to_radians()
        })
        .count()
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test flow -- --ignored"]
fn spot_control_mesh_two_subdivisions_ep_then_elambda() {
    let Some((mesh, pos)) = load_model("spot/spot/spot_control_mesh.obj") else { return };
    let mut opts = FlowOptions::uniform(EnergyOptions::default(), 2, 300);
    opts.rounds[0].energy = EnergyOptions::new(Energy::Combinatorial);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    print_round_table("spot", &res);
    check_rounds(&res);
    assert_eq!(res.mesh.genus(), 0);
    assert_eq!(res.mesh.euler_characteristic(), 2);
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test flow -- --ignored"]
fn bob_and_blub_control_meshes() {
    for rel in ["bob/bob/bob_controlmesh.obj", "blub/blub/blub_control_mesh.obj"] {
        let Some((mesh, pos)) = load_model(rel) else { continue };
        let opts = FlowOptions::uniform(EnergyOptions::default(), 1, 200);
        let res = flow(&mesh, &pos, &opts, None).unwrap();
        print_round_table(rel, &res);
        check_rounds(&res);
        assert_eq!(res.mesh.genus(), mesh.genus(), "{rel}: topology changed");
    }
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test flow -- --ignored"]
fn dense_models_run() {
    for rel in ["spot/spot/spot_triangulated.obj", "bob/bob/bob_tri.obj", "blub/blub/blub_triangulated.obj"] {
        let Some((mesh, pos)) = load_model(rel) else { continue };
        let opts = FlowOptions::uniform(EnergyOptions::default(), 0, 50);
        let res = flow(&mesh, &pos, &opts, None).unwrap();
        print_round_table(rel, &res);
        check_rounds(&res);
    }
}
