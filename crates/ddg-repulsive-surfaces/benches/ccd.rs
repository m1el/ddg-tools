//! CCD cost against the energy + gradient cost (`SPEC.md` §7.3, WP8). A plain `main` printing
//! CSV: for each mesh (5k–50k faces) and motion, the median time of `ccd::first_contact` and of
//! the Barnes–Hut (θ = 0.5) energy + differential, at 1 and 8 threads, and their ratio.
//!
//! Motions: `grad` = minus the L² gradient scaled so that the farthest vertex moves one mean
//! edge length at `t = 1` (`tau_end = 1`); `approach_x1` / `approach_x25` = the upper sheet of a
//! sheet pair moves onto the lower one by 1 / 2.5 mean edge lengths (the second runs into the
//! gap); `squeeze` = the torus is pushed towards its axis by one mean edge length.
//!
//! Run alone for timings: `cargo bench -p ddg-repulsive-surfaces --bench ccd`. The argument
//! `quick` runs the smallest sizes only.

use std::time::Instant;

use ddg_core::generators::{GeneratedMesh, grid_disk, icosphere, torus};
use ddg_repulsive_surfaces::ccd::{CcdOptions, first_contact};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::geometry::mean_edge_length;
use ddg_repulsive_surfaces::{BhOptions, EnergyEval, Exponents, Surface, energy};

type Mesh = (Surface, Vec<[f64; 3]>);

fn mesh(g: GeneratedMesh) -> Mesh {
    (
        Surface::new(g.positions.len(), g.faces).unwrap(),
        g.positions,
    )
}

fn sheets(n: usize, sep_edges: f64) -> Mesh {
    let g = grid_disk(n, n).unwrap();
    let a = Surface::new(g.positions.len(), g.faces.clone()).unwrap();
    let l = mean_edge_length(&a, &g.positions);
    let nv = g.positions.len() as u32;
    let mut pos = g.positions.clone();
    pos.extend(
        g.positions
            .iter()
            .map(|p| [p[0], p[1], p[2] + sep_edges * l]),
    );
    let mut faces = g.faces.clone();
    faces.extend(g.faces.iter().map(|f| f.map(|v| v + nv)));
    (Surface::new(pos.len(), faces).unwrap(), pos)
}

fn scaled_to(vel: &mut [[f64; 3]], length: f64) {
    let vmax = vel
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .fold(0.0, f64::max);
    for v in vel.iter_mut() {
        for c in v.iter_mut() {
            *c *= length / vmax;
        }
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

struct Case {
    name: &'static str,
    mesh: Mesh,
    kinds: Vec<&'static str>,
}

fn main() {
    let quick = std::env::args().any(|a| a == "quick");
    let ex = Exponents::default();
    let eval = EnergyEval::BarnesHut(BhOptions {
        theta: 0.5,
        leaf_size: 8,
    });
    let case = |name, mesh, kinds: &[&'static str]| Case {
        name,
        mesh,
        kinds: kinds.to_vec(),
    };
    let tor = |a, b, r| mesh(torus(a, b, 1.0, r).unwrap());
    let mut cases = vec![
        case("icosphere4_5k", mesh(icosphere(4).unwrap()), &["grad"]),
        case("torus_loose_10k", tor(100, 50, 0.35), &["grad", "squeeze"]),
        case("torus_tight_10k", tor(100, 50, 0.8), &["grad", "squeeze"]),
        case(
            "sheets_5k",
            sheets(35, 3.0),
            &["approach_x1", "approach_x25", "grad"],
        ),
    ];
    if !quick {
        cases.extend([
            case("icosphere5_20k", mesh(icosphere(5).unwrap()), &["grad"]),
            case("torus_loose_20k", tor(140, 70, 0.35), &["grad", "squeeze"]),
            case("torus_tight_20k", tor(140, 70, 0.8), &["grad", "squeeze"]),
            case("torus_tight_50k", tor(220, 110, 0.8), &["grad", "squeeze"]),
            case(
                "sheets_20k",
                sheets(70, 3.0),
                &["approach_x1", "approach_x25", "grad"],
            ),
            case(
                "sheets_50k",
                sheets(112, 3.0),
                &["approach_x1", "approach_x25", "grad"],
            ),
        ]);
    }
    println!(
        "case,motion,faces,threads,t_energy_s,t_ccd_s,ratio,tau_max,windows,pairs,iterations,close_pairs"
    );
    for c in &cases {
        let (s, x) = (&c.mesh.0, &c.mesh.1);
        let l = mean_edge_length(s, x);
        let opts = CcdOptions {
            reference_length: Some(l),
            ..Default::default()
        };
        let nv = x.len();
        for &kind in &c.kinds {
            for threads in [1usize, 8] {
                let rep = energy(s, x, ex, &eval, &[], Some(threads)).unwrap();
                let mut vel: Vec<[f64; 3]> = match kind {
                    "grad" => rep
                        .differential
                        .iter()
                        .map(|d| [-d[0], -d[1], -d[2]])
                        .collect(),
                    "squeeze" => x
                        .iter()
                        .map(|p| {
                            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
                            [-p[0] / r, -p[1] / r, 0.0]
                        })
                        .collect(),
                    _ => (0..nv)
                        .map(|i| if i >= nv / 2 { [0.0, 0.0, -1.0] } else { [0.0; 3] })
                        .collect(),
                };
                let len = if kind == "approach_x25" { 2.5 * l } else { l };
                scaled_to(&mut vel, len);
                let pool = make_pool(Some(threads)).unwrap();
                let reps = if s.n_faces() > 30000 { 3 } else { 5 };
                let te = median(
                    (0..reps)
                        .map(|_| energy(s, x, ex, &eval, &[], Some(threads)).unwrap().seconds)
                        .collect(),
                );
                let mut last = None;
                let tc = median(
                    (0..reps)
                        .map(|_| {
                            let t0 = Instant::now();
                            let r = first_contact(&pool, s, x, &vel, 1.0, &[], &opts).unwrap();
                            let dt = t0.elapsed().as_secs_f64();
                            last = Some(r);
                            dt
                        })
                        .collect(),
                );
                let r = last.unwrap();
                println!(
                    "{},{},{},{},{:.4},{:.4},{:.4},{:.4},{},{},{},{}",
                    c.name,
                    kind,
                    s.n_faces(),
                    threads,
                    te,
                    tc,
                    tc / te,
                    r.tau_max,
                    r.windows,
                    r.pairs_tested,
                    r.iterations,
                    r.close_pairs
                );
            }
        }
    }
}
