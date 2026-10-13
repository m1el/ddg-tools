//! Wall-clock time of the single-pass upwind distance on icospheres (one thread).
//!
//! `cargo run --release -p ddg-core --example distance_timings`

use std::time::Instant;

use ddg_core::distance::upwind_distance;
use ddg_core::generators::icosphere;
use ddg_core::{Geometry, TriMesh};

fn main() {
    println!("level  vertices  faces     setup ms  solve ms (best of 3)");
    for level in 6..=8u32 {
        let t0 = Instant::now();
        let m = icosphere(level).expect("icosphere");
        let mesh = TriMesh::new(m.positions.len(), m.faces.clone()).expect("mesh");
        let geometry = Geometry::from_positions(&mesh, &m.positions).expect("geometry");
        let setup = t0.elapsed().as_secs_f64() * 1e3;
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t = Instant::now();
            let r = upwind_distance(&mesh, &geometry, &[(0, 0.0)], None).expect("solve");
            best = best.min(t.elapsed().as_secs_f64() * 1e3);
            assert_eq!(r.n_unreached, 0);
        }
        println!(
            "{level:>5}  {:>8}  {:>7}  {setup:>8.0}  {best:>8.1}",
            m.positions.len(),
            m.faces.len()
        );
    }
}
