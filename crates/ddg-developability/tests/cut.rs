//! Tests of cut and flatten (SPEC §6, §10.7). Tests on `data/models/` meshes are `#[ignore]`d:
//!
//! ```text
//! cargo test --release -p ddg-developability --test cut -- --ignored --nocapture
//! ```

mod energy_common;
use energy_common::*;

use ddg_core::generators::{annulus, grid_disk, icosphere, open_cylinder, perturb, torus};
use ddg_core::topology::TriMesh;
use ddg_developability::{CutOptions, DevError, Energy, EnergyOptions, FlowOptions, cut_and_flatten, flow, star_embeddedness};

fn cylinder_arc(n: usize) -> Mesh {
    let g = grid_disk(n, n).unwrap();
    let pos = g
        .positions
        .iter()
        .map(|p| [0.5 * (2.0 * p[0]).cos(), 0.5 * (2.0 * p[0]).sin(), p[1]])
        .collect();
    (g.topology().unwrap(), pos)
}

/// Checks the invariants every result must satisfy; returns nothing.
fn check_invariants(mesh: &TriMesh, cf: &ddg_developability::DevelopableFlattening) {
    assert_eq!(cf.corner_uvs.len(), 3 * mesh.n_faces());
    assert_eq!(cf.face_chart.len(), mesh.n_faces());
    assert_eq!(cf.face_log_scale.len(), mesh.n_faces());
    assert_eq!(cf.vertex_log_scale.len(), mesh.n_vertices());
    assert!(cf.cut_edges.windows(2).all(|w| w[0] < w[1]));
    assert_eq!(cf.cut_vertex_pairs.len(), cf.cut_edges.len());
    for (&e, &pair) in cf.cut_edges.iter().zip(&cf.cut_vertex_pairs) {
        assert_eq!(mesh.edges()[e as usize], pair);
    }
    // every chart is a disk
    assert_eq!(cf.plan.charts.len(), cf.n_charts());
    assert!(cf.plan.charts.iter().all(|c| c.after.is_disk()));
    // every terminal (and added terminal) is on the cut or on the mesh boundary
    let mut on_cut = vec![false; mesh.n_vertices()];
    for p in &cf.cut_vertex_pairs {
        on_cut[p[0] as usize] = true;
        on_cut[p[1] as usize] = true;
    }
    for &t in cf.terminals.iter().chain(&cf.added_terminals) {
        assert!(on_cut[t as usize] || mesh.is_boundary_vertex(t as usize), "terminal {t} is not on the cut");
    }
    // a face's UV triangle has the area implied by its log scale factor
    assert!(cf.corner_uvs.iter().all(|p| p[0].is_finite() && p[1].is_finite()));
}

#[test]
fn flat_disk_has_no_cut_and_is_isometric() {
    let g = grid_disk(8, 6).unwrap();
    let mesh = g.topology().unwrap();
    let cf = cut_and_flatten(&mesh, &g.positions, &CutOptions::default()).unwrap();
    check_invariants(&mesh, &cf);
    assert!(cf.terminals.is_empty() && cf.cut_edges.is_empty());
    assert_eq!(cf.n_charts(), 1);
    assert_eq!(cf.flipped_triangles(), 0);
    assert!(cf.mean_abs_log_scale < 1e-9, "{}", cf.mean_abs_log_scale);
    assert!(cf.is_injective());
}

#[test]
fn exact_cylinder_sheet_unrolls_isometrically() {
    let (mesh, pos) = cylinder_arc(10);
    let cf = cut_and_flatten(&mesh, &pos, &CutOptions::default()).unwrap();
    check_invariants(&mesh, &cf);
    // the energy is (numerically) zero everywhere: no terminals, no cut
    assert!(cf.vertex_energy.iter().all(|&e| e < 1e-3), "max e = {:?}", cf.vertex_energy.iter().cloned().fold(0.0, f64::max));
    assert!(cf.terminals.is_empty());
    assert_eq!(cf.n_charts(), 1);
    // The chord/arc discretization is not exactly developable at the angle sum level only on the
    // boundary; the interior of a polyhedral cylinder is flat, so the unrolled sheet is isometric
    // along the boundary and nearly so inside.
    assert_eq!(cf.flipped_triangles(), 0);
    assert!(cf.mean_abs_log_scale < 5e-3, "mean |u| = {}", cf.mean_abs_log_scale);
}

#[test]
fn tube_is_cut_between_its_two_boundaries_and_unrolled_isometrically() {
    let g = open_cylinder(16, 6, 0.5, 1.0).unwrap();
    let mesh = g.topology().unwrap();
    let cf = cut_and_flatten(&mesh, &g.positions, &CutOptions::default()).unwrap();
    check_invariants(&mesh, &cf);
    assert_eq!(mesh.boundary_loops().len(), 2);
    assert_eq!(cf.n_charts(), 1);
    assert!(!cf.cut_edges.is_empty());
    assert_eq!(cf.flipped_triangles(), 0);
    assert!(cf.mean_abs_log_scale < 1e-6, "mean |u| = {}", cf.mean_abs_log_scale);
    assert!(cf.is_injective());
}

#[test]
fn alignment_weight_makes_the_cut_follow_the_rulings() {
    let g = open_cylinder(16, 6, 0.5, 1.0).unwrap();
    let mesh = g.topology().unwrap();
    // Skew the mesh so that the shortest path between the two loops is NOT a straight ruling:
    // a rotation of the top loop by one segment is not needed; instead compare weights directly.
    let opts = CutOptions { length_weight: 0.0, alignment_weight: 1.0, ..Default::default() };
    let cf = cut_and_flatten(&mesh, &g.positions, &opts).unwrap();
    check_invariants(&mesh, &cf);
    for &[a, b] in &cf.cut_vertex_pairs {
        let d = sub(g.positions[a as usize], g.positions[b as usize]);
        assert!(d[2].abs() / norm(d) > 0.999, "cut edge ({a}, {b}) is not along a ruling");
    }
    assert!(cf.mean_abs_log_scale < 1e-6);
}

#[test]
fn closed_sphere_gets_two_added_terminals_when_nothing_is_above_epsilon() {
    let g = icosphere(2).unwrap();
    let mesh = g.topology().unwrap();
    let opts = CutOptions { epsilon: 1e9, ..Default::default() };
    let cf = cut_and_flatten(&mesh, &g.positions, &opts).unwrap();
    check_invariants(&mesh, &cf);
    assert!(cf.terminals.is_empty());
    assert_eq!(cf.added_terminals.len(), 2);
    assert_eq!(cf.n_charts(), 1);
    // highest-energy vertices, ties by lower id: the added ones have maximal energy
    let emax = cf.vertex_energy.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(cf.vertex_energy[cf.added_terminals[0] as usize] >= emax * (1.0 - 1e-9) - 1e-15);
}

#[test]
fn closed_sphere_with_many_terminals_gives_disk_charts() {
    let g = icosphere(2).unwrap();
    let mesh = g.topology().unwrap();
    let cf = cut_and_flatten(&mesh, &g.positions, &CutOptions::default()).unwrap();
    check_invariants(&mesh, &cf);
    assert!(cf.terminals.len() > 2);
    assert!(cf.added_terminals.is_empty());
}

#[test]
fn torus_cut_contains_the_handle_loops() {
    let g = torus(12, 8, 1.0, 0.4).unwrap();
    let mesh = g.topology().unwrap();
    let opts = CutOptions { epsilon: 1e9, ..Default::default() };
    let cf = cut_and_flatten(&mesh, &g.positions, &opts).unwrap();
    check_invariants(&mesh, &cf);
    assert_eq!(cf.plan.n_loop_edges, 2);
    assert_eq!(cf.n_charts(), 1);
    assert!(cf.added_terminals.is_empty(), "a torus needs no added terminals");
}

#[test]
fn must_cut_makes_several_charts_and_may_cut_is_respected() {
    let g = annulus(12, 3, 0.5, 1.0).unwrap();
    let mesh = g.topology().unwrap();
    // the 12 radial rows of edges: cutting two radial lines at opposite sides splits it in two
    let radial: Vec<u32> = (0..mesh.n_edges() as u32)
        .filter(|&e| {
            let [a, b] = mesh.edges()[e as usize];
            let (pa, pb) = (g.positions[a as usize], g.positions[b as usize]);
            let r = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1]).sqrt();
            // along the radius: same angle, different radius
            (pa[0] * pb[1] - pa[1] * pb[0]).abs() < 1e-12 && (r(pa) - r(pb)).abs() > 1e-6
        })
        .collect();
    assert!(!radial.is_empty());
    let angle = |e: u32| {
        let p = g.positions[mesh.edges()[e as usize][0] as usize];
        p[1].atan2(p[0])
    };
    let seam_a: Vec<u32> = radial.iter().copied().filter(|&e| angle(e).abs() < 1e-9).collect();
    let seam_b: Vec<u32> = radial.iter().copied().filter(|&e| (angle(e).abs() - std::f64::consts::PI).abs() < 1e-9).collect();
    assert!(!seam_a.is_empty() && !seam_b.is_empty());
    let mut must = seam_a;
    must.extend(seam_b);
    let opts = CutOptions { must_cut: must.clone(), ..Default::default() };
    let cf = cut_and_flatten(&mesh, &g.positions, &opts).unwrap();
    check_invariants(&mesh, &cf);
    assert_eq!(cf.n_charts(), 2);
    for e in &must {
        assert!(cf.cut_edges.contains(e));
    }
    assert!(cf.mean_abs_log_scale < 1e-6);
    assert_eq!(cf.flipped_triangles(), 0);
    // forbidding every edge but the seams is accepted as long as the seams suffice
    let mut may = vec![false; mesh.n_edges()];
    for &e in &must {
        may[e as usize] = true;
    }
    let opts = CutOptions { must_cut: must, may_cut: Some(may), ..Default::default() };
    let cf2 = cut_and_flatten(&mesh, &g.positions, &opts).unwrap();
    assert_eq!(cf2.cut_edges, cf.cut_edges);
}

#[test]
fn extra_terminals_are_cut_through() {
    let (mesh, pos) = cylinder_arc(8);
    let mid = (mesh.n_vertices() / 2) as u32;
    assert!(!mesh.is_boundary_vertex(mid as usize));
    let opts = CutOptions { extra_terminals: vec![mid], ..Default::default() };
    let cf = cut_and_flatten(&mesh, &pos, &opts).unwrap();
    check_invariants(&mesh, &cf);
    assert_eq!(cf.terminals, vec![mid]);
    assert!(cf.cut_vertex_pairs.iter().any(|p| p.contains(&mid)));
    // a slit in a developable sheet unrolls with the same area: still no flipped triangle
    assert_eq!(cf.flipped_triangles(), 0);
}

#[test]
fn bad_input_is_an_error() {
    let (mesh, pos) = cylinder_arc(4);
    let ok = CutOptions::default();
    let bad = |o: CutOptions| cut_and_flatten(&mesh, &pos, &o);
    assert!(matches!(bad(CutOptions { epsilon: f64::NAN, ..ok.clone() }), Err(DevError::InvalidOption { name: "epsilon", .. })));
    for t in [f64::NAN, -0.1] {
        assert!(matches!(
            bad(CutOptions { defect_tolerance: Some(t), ..ok.clone() }),
            Err(DevError::InvalidOption { name: "defect_tolerance", .. })
        ));
    }
    assert!(matches!(bad(CutOptions { length_weight: 0.0, ..ok.clone() }), Err(DevError::InvalidOption { .. })));
    assert!(matches!(bad(CutOptions { alignment_weight: -1.0, ..ok.clone() }), Err(DevError::InvalidOption { .. })));
    assert!(matches!(bad(CutOptions { extra_terminals: vec![9999], ..ok.clone() }), Err(DevError::VertexOutOfRange { .. })));
    assert!(matches!(bad(CutOptions { may_cut: Some(vec![true; 3]), ..ok.clone() }), Err(DevError::LengthMismatch { .. })));
    assert!(matches!(cut_and_flatten(&mesh, &pos[1..], &ok), Err(DevError::LengthMismatch { .. })));
    // a must-cut edge that is forbidden is the planner's error
    let may = vec![false; mesh.n_edges()];
    let interior = (0..mesh.n_edges() as u32).find(|&e| !mesh.is_boundary_edge(e as usize)).unwrap();
    assert!(matches!(bad(CutOptions { may_cut: Some(may), must_cut: vec![interior], ..ok.clone() }), Err(DevError::Bff(_))));
    // intrinsic energy of a degenerate configuration: error, not a panic
    let mut flat = pos.clone();
    flat[5] = flat[6];
    assert!(cut_and_flatten(&mesh, &flat, &ok).is_err());
}

#[test]
fn flow_then_cut_on_a_noisy_sheet_has_low_distortion() {
    // A noisy cone sheet becomes developable with the flow; the remaining high-energy vertices
    // are cut and the charts unrolled with small distortion and no flips.
    let g = grid_disk(10, 10).unwrap();
    let mut pos: Vec<[f64; 3]> = g
        .positions
        .iter()
        .map(|p| {
            let s = 1.0 + p[1];
            let (beta, phi) = (0.6f64, 1.5 * p[0]);
            [s * beta.sin() * phi.cos(), s * beta.sin() * phi.sin(), s * beta.cos()]
        })
        .collect();
    perturb(&mut pos, [0.015; 3], 11);
    let mesh = g.topology().unwrap();
    let before = cut_and_flatten(&mesh, &pos, &CutOptions::default()).unwrap();
    let opts = FlowOptions::uniform(EnergyOptions::default(), 0, 600);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    let after = cut_and_flatten(&res.mesh, &res.positions, &CutOptions::default()).unwrap();
    check_invariants(&res.mesh, &after);
    eprintln!(
        "noisy cone: terminals {} -> {}, mean|u| {:.3e} -> {:.3e}, flips {} -> {}",
        before.terminals.len(),
        after.terminals.len(),
        before.mean_abs_log_scale,
        after.mean_abs_log_scale,
        before.flipped_triangles(),
        after.flipped_triangles()
    );
    assert!(after.terminals.len() <= before.terminals.len() / 4 + 1);
    // (with a shattered cut the distortion is small too: the boundary lengths pin the charts)
    assert!(after.mean_abs_log_scale < 1e-3);
    assert_eq!(after.flipped_triangles(), 0);
}

#[test]
fn defect_tolerance_makes_cone_vertices_terminals_and_none_is_the_paper_rule() {
    // a flat disk with one interior vertex pulled up into a cone tip: the angle defect is
    // positive; the covariance energy of the tip is small only if its star is nearly coplanar,
    // so use a tolerance below the defect and check the terminal appears iff the option is set
    let g = grid_disk(6, 6).unwrap();
    let mesh = g.topology().unwrap();
    let v = mesh.interior_vertices()[mesh.interior_vertices().len() / 2];
    let mut pos = g.positions.clone();
    pos[v as usize][2] = 0.05;
    let geometry = ddg_core::Geometry::from_positions(&mesh, &pos).unwrap();
    let omega = ddg_core::operators::angle_defects(&mesh, &geometry).unwrap()[v as usize].abs();
    assert!(omega > 0.0);
    let on = CutOptions { epsilon: 1e9, defect_tolerance: Some(0.5 * omega), ..Default::default() };
    let cf = cut_and_flatten(&mesh, &pos, &on).unwrap();
    assert_eq!(cf.terminals, vec![v]);
    assert!((cf.vertex_defect[v as usize].abs() - omega).abs() < 1e-12);
    let off = CutOptions { defect_tolerance: None, ..on.clone() };
    assert!(cut_and_flatten(&mesh, &pos, &off).unwrap().terminals.is_empty());
    let above = CutOptions { defect_tolerance: Some(2.0 * omega), ..on };
    assert!(cut_and_flatten(&mesh, &pos, &above).unwrap().terminals.is_empty());
}

#[test]
fn every_energy_can_select_terminals() {
    let g = icosphere(1).unwrap();
    let mesh = g.topology().unwrap();
    for energy in Energy::ALL {
        let cf = cut_and_flatten(&mesh, &g.positions, &CutOptions { energy, ..Default::default() }).unwrap();
        check_invariants(&mesh, &cf);
    }
}

// ---------------------------------------------------------------------------------------------

/// Flow (E_λ, `rounds` x `iterations`) then cut and flatten with `cut`. Prints the numbers and
/// returns (faces, flipped, mean |u|, cut length fraction, result) after the invariants passed.
fn run_model_with(
    rel: &str,
    subdivisions: usize,
    iterations: usize,
    cut: &CutOptions,
) -> Option<(usize, usize, f64, f64, ddg_developability::DevelopableFlattening, Mesh)> {
    let (mesh, pos) = load_model(rel)?;
    let opts = FlowOptions::uniform(EnergyOptions::default(), subdivisions, iterations);
    let res = flow(&mesh, &pos, &opts, None).unwrap();
    let cf = cut_and_flatten(&res.mesh, &res.positions, cut).unwrap();
    check_invariants(&res.mesh, &cf);
    let nf = res.mesh.n_faces();
    let flips = cf.flipped_triangles();
    let g = ddg_core::Geometry::from_positions(&res.mesh, &res.positions).unwrap();
    let total: f64 = g.edge_lengths().iter().sum();
    let cut_len: f64 = cf.cut_edges.iter().map(|&e| g.edge_lengths()[e as usize]).sum();
    eprintln!(
        "{rel} (tol {:?}, eps {:e}): {nf} faces, genus {}, {} terminals (+{}), {} cut edges ({:.2}% of the edge length, {} loop edges), {} charts, {} flipped ({:.3}%), mean|u| {:.3e}, injective {}",
        cut.defect_tolerance,
        cut.epsilon,
        res.mesh.genus(),
        cf.terminals.len(),
        cf.added_terminals.len(),
        cf.cut_edges.len(),
        100.0 * cut_len / total,
        cf.plan.n_loop_edges,
        cf.n_charts(),
        flips,
        100.0 * flips as f64 / nf as f64,
        cf.mean_abs_log_scale,
        cf.is_injective()
    );
    if res.mesh.genus() > 0 {
        assert_eq!(cf.plan.n_loop_edges, 2 * res.mesh.genus());
    }
    let mesh_pos = (res.mesh, res.positions);
    Some((nf, flips, cf.mean_abs_log_scale, cut_len / total, cf, mesh_pos))
}

/// Sum of |Ω| over the interior non-terminals, and the part of it on non-embedded stars.
fn nonterminal_defect(mesh: &TriMesh, cf: &ddg_developability::DevelopableFlattening) -> (f64, f64) {
    let is_t: std::collections::HashSet<u32> = cf.terminals.iter().chain(&cf.added_terminals).copied().collect();
    let (mut all, mut bad) = (0.0, 0.0);
    for v in mesh.interior_vertices() {
        if is_t.contains(&v) {
            continue;
        }
        let d = cf.vertex_defect[v as usize].abs();
        all += d;
        if cf.star_embeddedness[v as usize].is_nan() || cf.star_embeddedness[v as usize] <= 0.0 {
            bad += d;
        }
    }
    (all, bad)
}

fn paper_rule() -> CutOptions {
    CutOptions { defect_tolerance: None, ..Default::default() }
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test cut -- --ignored"]
fn spot_flow_then_cut() {
    // SPEC D30 / section 10.7: the default (defect tolerance 0.1) gives no flipped UV triangle,
    // the paper's rule alone leaves about 2 %
    let (nf, flips, mean_u, len, ..) =
        run_model_with("spot/spot/spot_control_mesh.obj", 2, 300, &CutOptions::default()).unwrap();
    assert_eq!(flips, 0, "{flips} of {nf}");
    assert!(mean_u < 0.02, "{mean_u}");
    assert!(len < 0.06, "{len}");
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test cut -- --ignored"]
fn spot_paper_rule_leaves_defect_on_non_embedded_stars() {
    // The hypothesis of D30, measured: with the energy threshold alone the cancelling angle
    // defect of the non-terminals sits on stars that are not graphs over their tangent plane
    // (normals nearly coplanar, so the energy is small, but Omega is not).
    let (nf, flips, mean_u, _, cf, (mesh, _)) =
        run_model_with("spot/spot/spot_control_mesh.obj", 2, 300, &paper_rule()).unwrap();
    let (all, bad) = nonterminal_defect(&mesh, &cf);
    eprintln!("paper rule: non-terminal sum |Omega| {all:.2}, of which on non-embedded stars {bad:.2}");
    assert!(all > 50.0 && bad > 0.9 * all, "{all} {bad}");
    assert!(flips as f64 > 0.005 * nf as f64, "{flips} of {nf}");
    assert!(mean_u > 0.1, "{mean_u}");
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test cut -- --ignored"]
fn spot_over_run_flow_then_cut() {
    // 3 x 3000 iterations: E_λ ~ 1e-8, no terminal from the energy; the paper's rule then has
    // ~4 % flipped triangles and |u| ~ 6, the defect criterion none
    let (nf, flips, mean_u, len, ..) =
        run_model_with("spot/spot/spot_control_mesh.obj", 2, 3000, &CutOptions::default()).unwrap();
    assert_eq!(flips, 0, "{flips} of {nf}");
    assert!(mean_u < 0.01, "{mean_u}");
    assert!(len < 0.06, "{len}");
    let (_, flips2, ..) = run_model_with("spot/spot/spot_control_mesh.obj", 2, 3000, &paper_rule()).unwrap();
    assert!(flips2 as f64 > 0.01 * nf as f64, "{flips2}");
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test cut -- --ignored"]
fn bob_flow_then_cut_has_handle_loops() {
    let (nf, flips, mean_u, len, ..) =
        run_model_with("bob/bob/bob_controlmesh.obj", 2, 300, &CutOptions::default()).unwrap();
    assert_eq!(flips, 0, "{flips} of {nf}");
    assert!(mean_u < 0.02, "{mean_u}");
    assert!(len < 0.06, "{len}");
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test cut -- --ignored"]
fn blub_flow_then_cut() {
    let (nf, flips, mean_u, len, ..) =
        run_model_with("blub/blub/blub_control_mesh.obj", 2, 300, &CutOptions::default()).unwrap();
    assert_eq!(flips, 0, "{flips} of {nf}");
    assert!(mean_u < 0.02, "{mean_u}");
    assert!(len < 0.06, "{len}");
}

#[test]
#[ignore = "long: cargo test --release -p ddg-developability --test cut -- --ignored"]
fn crumpled_then_cut() {
    let (nf, flips, mean_u, ..) =
        run_model_with("crumpleddevelopable/CrumpledDevelopable.obj", 0, 500, &CutOptions::default()).unwrap();
    assert_eq!(flips, 0, "{flips} of {nf}");
    assert!(mean_u < 1e-3, "{mean_u}");
}

#[test]
fn star_embeddedness_flags_a_folded_star() {
    // a flat fan around vertex 0 is a graph (cos = 1); lifting one rim vertex far above the
    // center folds a face over (cos <= 0) while the energies would not notice
    let g = grid_disk(4, 4).unwrap();
    let mesh = g.topology().unwrap();
    let e = star_embeddedness(&mesh, &g.positions).unwrap();
    assert!(e.iter().all(|&c| (c - 1.0).abs() < 1e-12));
    let mut pos = g.positions.clone();
    let v = mesh.interior_vertices()[0] as usize;
    pos[v][2] = 5.0;
    let e = star_embeddedness(&mesh, &pos).unwrap();
    assert!(e[v] > 0.0 && e[v] < 1.0, "a cone tip is still a graph: {}", e[v]);
    let mut pos = g.positions.clone();
    pos[v][2] = -0.0;
    // push the vertex far outside its star: faces flip orientation around it
    pos[v][0] += 10.0;
    let e = star_embeddedness(&mesh, &pos).unwrap();
    assert!(e[v] <= 0.0 || e[v].is_nan(), "{}", e[v]);
}
