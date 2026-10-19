//! Invariances and determinism (SPEC §10.4).

mod energy_common;
use energy_common::*;

use ddg_core::TriMesh;
use ddg_developability::{Energy, EnergyOptions, Evaluation, evaluate};

fn opts(energy: Energy) -> EnergyOptions {
    EnergyOptions {
        energy,
        include_boundary: true,
        omit_valence_three: false,
        threads: None,
    }
}

fn meshes() -> Vec<Mesh> {
    vec![
        perturbed_icosphere(2, 0.02, 3),
        bumpy_grid(10, 0.2, 5),
        perturbed_torus(8),
    ]
}

/// Vertices whose gradient is determined away from ties: every star containing it has a margin.
fn stable_vertices(mesh: &TriMesh, e: &Evaluation) -> Vec<usize> {
    (0..mesh.n_vertices())
        .filter(|&u| {
            std::iter::once(u)
                .chain(mesh.vertex_neighbors(u))
                .all(|w| !e.counted[w] || e.vertex_margin[w] > 1e-8 * e.vertex_energy[w].max(1.0))
        })
        .collect()
}

fn gmax(e: &Evaluation) -> f64 {
    e.gradient.iter().map(|g| norm(*g)).fold(0.0, f64::max)
}

#[test]
fn rigid_motion_and_scale() {
    let rot = random_rotation(42);
    let shift = [3.0, -1.5, 0.25];
    for (m, p) in meshes() {
        for energy in Energy::ALL {
            let e0 = evaluate(&m, &p, &opts(energy)).unwrap();
            let stable = stable_vertices(&m, &e0);
            assert!(
                stable.len() > m.n_vertices() / 4,
                "{energy:?}: too few stable vertices"
            );
            let g0 = gmax(&e0);
            // rigid motion
            let pr: Vec<_> = p.iter().map(|&x| add(rotate(&rot, x), shift)).collect();
            let er = evaluate(&m, &pr, &opts(energy)).unwrap();
            assert!(
                (er.energy - e0.energy).abs() <= 1e-12 * e0.energy.abs().max(1e-300),
                "{energy:?} rigid"
            );
            for v in 0..m.n_vertices() {
                let (a, b) = (e0.vertex_energy[v], er.vertex_energy[v]);
                assert!(
                    (a - b).abs() <= 1e-12 * a.abs().max(1.0),
                    "{energy:?} rigid vertex {v}: {a} vs {b}"
                );
            }
            for &v in &stable {
                let want = rotate(&rot, e0.gradient[v]);
                let err = norm(sub(er.gradient[v], want));
                assert!(
                    err <= 1e-9 * g0,
                    "{energy:?} rigid gradient {v}: {err:e} (|g| {g0:e})"
                );
            }
            // uniform scaling
            for s in [1e-3, 1e3] {
                let ps: Vec<_> = p.iter().map(|&x| scale(x, s)).collect();
                let es = evaluate(&m, &ps, &opts(energy)).unwrap();
                assert!(
                    (es.energy - e0.energy).abs() <= 1e-12 * e0.energy.abs(),
                    "{energy:?} scale {s}"
                );
                for &v in &stable {
                    let want = scale(e0.gradient[v], 1.0 / s);
                    let err = norm(sub(es.gradient[v], want));
                    assert!(
                        err <= 1e-9 * g0 / s,
                        "{energy:?} scale {s} gradient {v}: {err:e}"
                    );
                }
            }
        }
    }
}

#[test]
fn gradient_is_orthogonal_to_similarity_motions() {
    for (m, p) in meshes() {
        for energy in Energy::ALL {
            for o in [opts(energy), EnergyOptions::new(energy)] {
                let e = evaluate(&m, &p, &o).unwrap();
                let mut sg = [0.0; 3];
                let mut sfg = [0.0; 3];
                let mut sdot = 0.0;
                let mut ng = 0.0;
                let mut nfg = 0.0;
                for (f, g) in p.iter().zip(&e.gradient) {
                    sg = add(sg, *g);
                    sfg = add(sfg, cross(*f, *g));
                    sdot += dot(*f, *g);
                    ng += norm(*g);
                    nfg += norm(*f) * norm(*g);
                }
                assert!(
                    norm(sg) <= 1e-10 * ng,
                    "{energy:?}: |sum g| = {:e} vs {ng:e}",
                    norm(sg)
                );
                assert!(
                    norm(sfg) <= 1e-10 * nfg,
                    "{energy:?}: |sum f x g| = {:e}",
                    norm(sfg)
                );
                assert!(sdot.abs() <= 1e-10 * nfg, "{energy:?}: sum f.g = {sdot:e}");
            }
        }
    }
}

/// Split the interior edge `(a, b)` at its midpoint `m` (new vertex) and each of its two faces
/// into two coplanar faces.
fn split_edge(mesh: &TriMesh, pos: &[[f64; 3]], e: usize) -> Mesh {
    let [a, b] = mesh.edges()[e];
    let mid = scale(add(pos[a as usize], pos[b as usize]), 0.5);
    let nm = mesh.n_vertices() as u32;
    let mut faces = Vec::new();
    for &t in mesh.faces() {
        let has = |x: u32| t.contains(&x);
        if has(a) && has(b) {
            // rotate so that t = (p, q, r) with edge (p, q) = (a, b) or (b, a)
            let i = (0..3).find(|&i| {
                let (p, q) = (t[i], t[(i + 1) % 3]);
                (p == a && q == b) || (p == b && q == a)
            });
            let i = i.unwrap();
            let (p, q, r) = (t[i], t[(i + 1) % 3], t[(i + 2) % 3]);
            faces.push([p, nm, r]);
            faces.push([nm, q, r]);
        } else {
            faces.push(t);
        }
    }
    let mut p2 = pos.to_vec();
    p2.push(mid);
    (TriMesh::new(mesh.n_vertices() + 1, faces).unwrap(), p2)
}

#[test]
fn retessellation_keeps_covariance_energies() {
    let (m, p) = perturbed_icosphere(2, 0.02, 9);
    for e in [3, 50, 200] {
        let (m2, p2) = split_edge(&m, &p, e);
        for energy in [Energy::Covariance, Energy::CovarianceMax] {
            let o = opts(energy);
            let a = evaluate(&m, &p, &o).unwrap();
            let b = evaluate(&m2, &p2, &o).unwrap();
            let rel = (a.energy - b.energy).abs() / a.energy;
            assert!(
                rel <= 1e-12,
                "{energy:?} edge {e}: {} vs {} (rel {rel:e})",
                a.energy,
                b.energy
            );
            // the new vertex is a hinge
            assert!(b.vertex_energy[m.n_vertices()].abs() <= 1e-14);
        }
    }
}

#[test]
fn deterministic_for_any_thread_count() {
    let (m, p) = bumpy_grid(40, 0.2, 17); // 1681 vertices: several work chunks
    for energy in Energy::ALL {
        let run = |t: usize| {
            evaluate(
                &m,
                &p,
                &EnergyOptions {
                    threads: Some(t),
                    ..opts(energy)
                },
            )
            .unwrap()
        };
        let base = run(1);
        let bits = |e: &Evaluation| {
            let mut v = vec![e.energy.to_bits(), e.near_ties as u64];
            v.extend(e.vertex_energy.iter().map(|x| x.to_bits()));
            v.extend(e.vertex_margin.iter().map(|x| x.to_bits()));
            v.extend(e.gradient.iter().flat_map(|g| g.map(f64::to_bits)));
            v
        };
        let b0 = bits(&base);
        for t in [2, 7, 64] {
            assert!(bits(&run(t)) == b0, "{energy:?}: {t} threads differ from 1");
        }
        let none = evaluate(&m, &p, &opts(energy)).unwrap();
        assert!(bits(&none) == b0, "{energy:?}: default pool differs");
        let (en, ve) = ddg_developability::evaluate_energy(&m, &p, &opts(energy)).unwrap();
        assert_eq!(en.to_bits(), base.energy.to_bits());
        assert_eq!(ve, base.vertex_energy);
    }
}
