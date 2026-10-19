//! Energy + gradient evaluation (SPEC §10.8 `energy`):
//! `cargo bench -p ddg-developability --bench energy [-- <mesh-name-filter>...]`.
//!
//! For each of the six energies, `evaluate` (value and gradient) on spot (5856 faces), bob
//! (10688), blub (14208) and spot subdivided twice (93696), with one thread and with all
//! threads; reports the median time per call and per vertex. The pool is built once per
//! measurement (`evaluate` reuses the pool it runs in). Meshes come from `data/models/` (unzip
//! `/var/data/discrete-geometry/models/{spot,bob,blub}.zip` there); missing ones are skipped.
//! Single-thread numbers are meaningful only when nothing else runs on the machine.
//! `DEV_BENCH_THREADS=n` replaces "all threads" by `n` threads (`1`: single-thread run only).

use std::path::PathBuf;
use std::time::{Duration, Instant};

use ddg_core::TriMesh;
use ddg_developability::{Energy, EnergyOptions, evaluate, subdivide};

fn load(rel: &str) -> Option<(TriMesh, Vec<[f64; 3]>)> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data/models")
        .join(rel);
    let Ok(obj) = ddg_core::io::read_obj(&path) else {
        eprintln!("skipping {}: not found or unreadable", path.display());
        return None;
    };
    let faces = obj.triangles()?;
    Some((
        TriMesh::new(obj.positions.len(), faces).ok()?,
        obj.positions,
    ))
}

/// Median seconds per call of `f`, after one warm-up call; at least 3 calls and ~1 s.
fn time(mut f: impl FnMut()) -> f64 {
    f();
    let mut samples = Vec::new();
    let start = Instant::now();
    while samples.len() < 3 || (start.elapsed() < Duration::from_secs(1) && samples.len() < 200) {
        let t = Instant::now();
        f();
        samples.push(t.elapsed().as_secs_f64());
    }
    samples.sort_by(f64::total_cmp);
    samples[samples.len() / 2]
}

fn main() {
    let filters: Vec<String> = std::env::args()
        .skip(1)
        .filter(|a| !a.starts_with('-'))
        .collect();
    let mut meshes = Vec::new();
    for (name, rel) in [
        ("spot", "spot/spot/spot_triangulated.obj"),
        ("bob", "bob/bob/bob_tri.obj"),
        ("blub", "blub/blub/blub_triangulated.obj"),
    ] {
        if let Some(m) = load(rel) {
            meshes.push((name.to_string(), m));
        }
    }
    if let Some((_, (m, p))) = meshes.iter().find(|(n, _)| n == "spot") {
        let s1 = subdivide(m, p).expect("subdivide");
        let s2 = subdivide(&s1.mesh, &s1.positions).expect("subdivide");
        meshes.push(("spot-sub2".to_string(), (s2.mesh, s2.positions)));
    }
    let all = std::env::var("DEV_BENCH_THREADS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    println!("threads available: {all}");
    println!(
        "{:<10} {:>6} {:>6}  {:<17} {:>11} {:>10} {:>11} {:>10} {:>7}",
        "mesh", "|V|", "|F|", "energy", "1T ms", "1T ns/v", "allT ms", "allT ns/v", "speedup"
    );
    for (name, (m, p)) in &meshes {
        if !filters.is_empty() && !filters.iter().any(|f| name.contains(f.as_str())) {
            continue;
        }
        let nv = m.n_vertices() as f64;
        for energy in Energy::ALL {
            let mut t = [0.0; 2];
            for (i, threads) in [1, all].into_iter().enumerate() {
                if i == 1 && all == 1 {
                    t[1] = t[0];
                    continue;
                }
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let opts = EnergyOptions {
                    threads: Some(threads),
                    ..EnergyOptions::new(energy)
                };
                t[i] = pool.install(|| {
                    time(|| {
                        let e = evaluate(m, p, &opts).expect("evaluate");
                        std::hint::black_box(e);
                    })
                });
            }
            println!(
                "{:<10} {:>6} {:>6}  {:<17} {:>11.3} {:>10.0} {:>11.3} {:>10.0} {:>7.1}",
                name,
                m.n_vertices(),
                m.n_faces(),
                format!("{energy:?}"),
                t[0] * 1e3,
                t[0] * 1e9 / nv,
                t[1] * 1e3,
                t[1] * 1e9 / nv,
                t[0] / t[1]
            );
        }
    }
}
