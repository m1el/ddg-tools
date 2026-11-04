//! CPU network throughput (SPEC.md §5.3 re-measured in the crate).
//!
//! `cargo run --release -p ddg-points-as-tori --example throughput -- <points.npy> <normals.npy>
//!  <threads> [reps] [limit]`
//!
//! Prints neighborhoods/s per repetition and an FNV hash of the coefficient bits. Run alone; the
//! caller records `uptime` (outside load).

use ddg_points_as_tori::network::{PredictOptions, predict_coefficients};
use ddg_points_as_tori::npy::read_npy;
use ddg_points_as_tori::spatial::{KdTree, KdTreeOptions, fitting_neighborhoods};
use ddg_points_as_tori::{Model, NormalPolicy, PointCloud};
use std::path::Path;

fn to3(v: &[f64]) -> Vec<[f64; 3]> {
    v.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: throughput <points.npy> <normals.npy> <threads> [reps] [limit]");
        std::process::exit(2);
    }
    let mut p = to3(&read_npy(Path::new(&a[1]))?.data);
    let mut n = to3(&read_npy(Path::new(&a[2]))?.data);
    let threads: usize = a[3].parse()?;
    let reps: usize = a.get(4).map(|s| s.parse()).transpose()?.unwrap_or(3);
    if let Some(l) = a.get(5) {
        let l: usize = l.parse()?;
        p.truncate(l);
        n.truncate(l);
    }
    let (cloud, _) = PointCloud::new(p, n, NormalPolicy::Normalize)?;
    let model = Model::reference()?;
    let tree = KdTree::new(cloud.positions(), KdTreeOptions::default())?;
    let nb = fitting_neighborhoods(&tree, cloud.positions(), 64);
    let opts = PredictOptions { threads: Some(threads), ..Default::default() };
    for r in 0..reps {
        let t0 = std::time::Instant::now();
        let (co, _, rep) = predict_coefficients(&model, &cloud, &nb, &opts)?;
        let dt = t0.elapsed().as_secs_f64();
        let mut h: u64 = 0xcbf29ce484222325;
        for c in &co {
            for x in c {
                h = (h ^ x.to_bits() as u64).wrapping_mul(0x100000001b3);
            }
        }
        println!(
            "rep {r} threads {threads} n {} : {:.3} s, {:.0} neighborhoods/s (features {:.2} s, forward {:.2} s summed), hash {h:016x}",
            cloud.len(),
            dt,
            cloud.len() as f64 / dt,
            rep.features_seconds,
            rep.forward_seconds
        );
    }
    Ok(())
}
