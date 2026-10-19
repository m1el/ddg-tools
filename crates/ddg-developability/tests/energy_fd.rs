//! Finite-difference checks of the energies' gradients (SPEC §10.2).
//!
//! Central differences with `h = 1e-6 ℓ̄` per coordinate of sampled vertices; the difference
//! quotient uses the per-vertex energies of the affected stars only (less cancellation). A vertex
//! is perturbed only if every counted star containing it has a margin above `1e-4 · max(1, e_v)`
//! (stricter than the SPEC's `1e-6`, so that the active piece cannot switch at `±h`), and the
//! margins are re-checked at `±h`.

mod energy_common;
use energy_common::*;

use ddg_core::TriMesh;
use ddg_developability::{Energy, EnergyOptions, evaluate, evaluate_energy};

const MARGIN: f64 = 1e-4;

struct FdStats {
    max_rel: f64,
    scale: f64,
    checked: usize,
}

fn fd_check(mesh: &TriMesh, pos: &[[f64; 3]], opts: &EnergyOptions, max_samples: usize) -> FdStats {
    let ev = evaluate(mesh, pos, opts).unwrap();
    let nv = mesh.n_vertices();
    let affected: Vec<Vec<usize>> = (0..nv)
        .map(|u| std::iter::once(u).chain(mesh.vertex_neighbors(u)).collect())
        .collect();
    let ok = |e: &ddg_developability::Evaluation, w: usize, thr: f64| {
        !e.counted[w] || e.vertex_margin[w] > thr * e.vertex_energy[w].max(1.0)
    };
    let eligible: Vec<usize> = (0..nv)
        .filter(|&u| affected[u].iter().all(|&w| ok(&ev, w, MARGIN)))
        .collect();
    let step = (eligible.len() / max_samples.max(1)).max(1);
    let lbar = mean_edge_length(mesh, pos);
    // Relative to the largest gradient entry, floored at 1e-3/ℓ̄ (stars whose energy vanishes
    // identically, e.g. boundary stars with <= 2 faces, have gradients at rounding level).
    let gscale = ev
        .gradient
        .iter()
        .map(|g| norm(*g))
        .fold(1e-3 / lbar, f64::max);
    let h = 1e-6 * lbar;
    let mut max_err = 0.0f64;
    let mut checked = 0;
    for &u in eligible.iter().step_by(step) {
        for d in 0..3 {
            let mut vals = [0.0; 2];
            for (i, sgn) in [1.0, -1.0].into_iter().enumerate() {
                let mut p = pos.to_vec();
                p[u][d] += sgn * h;
                let e = evaluate(mesh, &p, opts).unwrap();
                for &w in &affected[u] {
                    assert!(
                        ok(&e, w, 1e-6),
                        "margin of star {w} collapsed at +-h: {} -> {} (e {} -> {}), {:?}",
                        ev.vertex_margin[w],
                        e.vertex_margin[w],
                        ev.vertex_energy[w],
                        e.vertex_energy[w],
                        opts.energy
                    );
                }
                let (_, ve) = evaluate_energy(mesh, &p, opts).unwrap();
                assert_eq!(ve, e.vertex_energy);
                vals[i] = affected[u]
                    .iter()
                    .filter(|&&w| ev.counted[w])
                    .map(|&w| ve[w])
                    .sum();
            }
            let fd = (vals[0] - vals[1]) / (2.0 * h);
            max_err = max_err.max((fd - ev.gradient[u][d]).abs() / gscale);
        }
        checked += 1;
    }
    FdStats {
        max_rel: max_err,
        scale: gscale * lbar,
        checked,
    }
}

fn opts_all(energy: Energy) -> EnergyOptions {
    EnergyOptions {
        energy,
        include_boundary: true,
        omit_valence_three: false,
        threads: None,
    }
}

fn run(name: &str, meshes: &[Mesh], max_samples: usize, opts_of: impl Fn(Energy) -> EnergyOptions) {
    for energy in Energy::ALL {
        let mut worst = 0.0f64;
        let mut total = 0;
        for (m, p) in meshes {
            let st = fd_check(m, p, &opts_of(energy), max_samples);
            if st.max_rel > 1e-8 {
                println!(
                    "  {name} {energy:?}: mesh with {} faces: rel err {:.2e} (scale {:.2e})",
                    m.n_faces(),
                    st.max_rel,
                    st.scale
                );
            }
            worst = worst.max(st.max_rel);
            total += st.checked;
        }
        println!("FD {name:>18} {energy:?}: max rel err {worst:.2e} over {total} vertices");
        assert!(total > 0, "{name} {energy:?}: no vertex had a margin");
        assert!(worst <= 1e-6, "{name} {energy:?}: FD rel err {worst:e}");
    }
}

#[test]
fn fd_random_interior_stars() {
    let meshes: Vec<Mesh> = (3..=9)
        .flat_map(|k| (0..4).map(move |s| interior_star(k, 100 * k as u64 + s, 0.5)))
        .collect();
    run("interior stars", &meshes, 100, opts_all);
}

#[test]
fn fd_boundary_fans() {
    let meshes: Vec<Mesh> = (1..=7)
        .flat_map(|k| (0..3).map(move |s| boundary_fan(k, 1000 + 10 * k as u64 + s, 0.5)))
        .collect();
    run("boundary fans", &meshes, 100, opts_all);
}

#[test]
fn fd_small_meshes() {
    let meshes = vec![
        perturbed_icosphere(1, 0.08, 7),
        bumpy_grid(6, 0.15, 3),
        perturbed_torus(11),
    ];
    run("small meshes", &meshes, 30, opts_all);
    // default domain (valence-3 omission, boundary counted) and no boundary
    run("small meshes (def)", &meshes, 15, EnergyOptions::new);
    run("small meshes (nob)", &meshes, 15, |e| EnergyOptions {
        include_boundary: false,
        ..EnergyOptions::new(e)
    });
}

#[test]
fn fd_crumpled_and_spot() {
    let mut meshes = Vec::new();
    meshes.extend(load_model("crumpleddevelopable/CrumpledDevelopable.obj"));
    meshes.extend(load_model("spot/spot/spot_control_mesh.obj"));
    if meshes.is_empty() {
        return;
    }
    run("crumpled+spot ctrl", &meshes, 25, opts_all);
}
