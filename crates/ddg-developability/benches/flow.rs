//! Benchmark of the flow (SPEC §10.8): `cargo bench -p ddg-developability --bench flow`.
//!
//! Part 1: `N` iterations (default 50, env `FLOW_BENCH_ITERATIONS`) on the dense models of
//! `data/models` (spot_triangulated, bob_tri, blub_triangulated), `E_λ` and the intrinsic energy,
//! no subdivision, remeshing on; time per iteration and per evaluation.
//! Part 2: the setups of SPEC §10.8 (spot control mesh + 2 subdivisions, crumpled, bob control mesh)
//! with a fixed iteration budget per round (env `FLOW_BENCH_BUDGET`, default 300).
//!
//! Every job runs with 1 thread and with all threads. Single-thread numbers are only meaningful
//! when nothing else is running: run this alone. The paper reports "seconds to a few minutes on a
//! 3.4GHz Intel Core i7 laptop ... we did not use multithreading".

use std::path::PathBuf;

use ddg_core::TriMesh;
use ddg_developability::{Energy, EnergyOptions, FlowOptions, flow};

fn load(rel: &str) -> Option<(TriMesh, Vec<[f64; 3]>)> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/models").join(rel);
    let obj = ddg_core::io::read_obj(&path).ok()?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let c: Vec<u32> = obj.face_corners(f).map(|i| obj.corner_verts[i]).collect();
        for i in 1..c.len() - 1 {
            faces.push([c[0], c[i], c[i + 1]]);
        }
    }
    Some((TriMesh::new(obj.positions.len(), faces).ok()?, obj.positions))
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn run(label: &str, mesh: &TriMesh, pos: &[[f64; 3]], energy: Energy, subdivisions: usize, iters: usize, threads: usize) {
    let mut opts = FlowOptions::uniform(EnergyOptions::new(energy), subdivisions, iters);
    opts.threads = Some(threads);
    // warm up the allocator and the pool once on a short run
    let mut warm = opts.clone();
    for r in &mut warm.rounds {
        r.max_iterations = 2;
    }
    let _ = flow(mesh, pos, &warm, None);
    let res = flow(mesh, pos, &opts, None);
    match res {
        Ok(r) => {
            let it: usize = r.report.rounds.iter().map(|x| x.optimize.iterations).sum();
            let ev: usize = r.report.rounds.iter().map(|x| x.optimize.evaluations).sum();
            let last = r.report.rounds.last().unwrap();
            let secs = r.report.seconds;
            println!(
                "{label:<34} {energy:<16?} thr {threads:>2}  {:>6} f  {it:>5} it {ev:>5} ev  {secs:>8.2} s  {:>8.2} ms/it  {:>7.2} ms/ev  E {:.3e} -> {:.3e}",
                last.n_faces,
                1e3 * secs / it.max(1) as f64,
                1e3 * secs / ev.max(1) as f64,
                r.report.rounds[0].optimize.energy_initial,
                last.optimize.energy_final
            );
        }
        Err(e) => println!("{label:<34} {energy:?} thr {threads}: error {e}"),
    }
}

fn main() {
    let all = std::thread::available_parallelism().map_or(1, |n| n.get());
    let n = env_usize("FLOW_BENCH_ITERATIONS", 50);
    let budget = env_usize("FLOW_BENCH_BUDGET", 300);
    println!("# part 1: {n} iterations, no subdivision, remeshing on");
    for rel in ["spot/spot/spot_triangulated.obj", "bob/bob/bob_tri.obj", "blub/blub/blub_triangulated.obj"] {
        let Some((mesh, pos)) = load(rel) else {
            println!("skipping {rel}: not found");
            continue;
        };
        for energy in [Energy::Covariance, Energy::Intrinsic] {
            for threads in [1, all] {
                run(rel.rsplit('/').next().unwrap(), &mesh, &pos, energy, 0, n, threads);
            }
        }
    }
    println!("# part 2: SPEC 10.8 setups, {budget} iterations per round, E_lambda");
    for (rel, subs) in [
        ("spot/spot/spot_control_mesh.obj", 2),
        ("crumpleddevelopable/CrumpledDevelopable.obj", 0),
        ("bob/bob/bob_controlmesh.obj", 2),
    ] {
        let Some((mesh, pos)) = load(rel) else {
            println!("skipping {rel}: not found");
            continue;
        };
        for threads in [1, all] {
            run(&format!("{} +{subs} sub", rel.rsplit('/').next().unwrap()), &mesh, &pos, Energy::Covariance, subs, budget, threads);
        }
    }
}
