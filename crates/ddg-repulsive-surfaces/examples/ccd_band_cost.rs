//! Cost of the start-of-step intersection test (`ccd::IntersectionBand::new`, SPEC 7.3) on
//! embedded and intersecting tori of 10k to 50k faces, with 1 and 8 threads (best of 7), next to
//! the CCD itself (`first_contact` on a push of `0.3 ℓ̄`).
//!
//! ```text
//! cargo run --release -p ddg-repulsive-surfaces --example ccd_band_cost
//! ```

use std::time::Instant;

use ddg_repulsive_surfaces::ccd::{CcdOptions, IntersectionBand, first_contact_with_band};
use ddg_repulsive_surfaces::energy::exact::make_pool;
use ddg_repulsive_surfaces::generate::primitives::torus;
use ddg_repulsive_surfaces::geometry::mean_edge_length;

fn best<T>(mut f: impl FnMut() -> T) -> (f64, T) {
    let mut out = f();
    let mut t = f64::INFINITY;
    for _ in 0..7 {
        let s = Instant::now();
        out = f();
        t = t.min(s.elapsed().as_secs_f64());
    }
    (t, out)
}

fn main() {
    println!("case,faces,threads,band_ms,band_faces_pairs,ccd_ms");
    for (nu, nv) in [(100usize, 50usize), (160, 80), (250, 100)] {
        for crossing in [false, true] {
            let a = torus(nu, nv, 1.0, 0.3).unwrap();
            let mesh = if crossing {
                let mut b = torus(nu, nv, 1.0, 0.3).unwrap();
                for p in &mut b.positions {
                    *p = [p[0] + 0.5, -p[2], p[1]];
                }
                a.merged(&b)
            } else {
                a
            };
            let s = mesh.surface().unwrap();
            let p = mesh.positions.clone();
            let l = mean_edge_length(&s, &p);
            // A push along the vertex normal-ish direction (towards the origin axis).
            let vel: Vec<[f64; 3]> = p.iter().map(|q| [-0.3 * l * q[0], -0.3 * l * q[1], -0.3 * l * q[2]]).collect();
            for threads in [1usize, 8] {
                let pool = make_pool(Some(threads)).unwrap();
                let (tb, band) = best(|| IntersectionBand::new(&pool, &s, &p).unwrap());
                let o = CcdOptions { reference_length: Some(l), ..Default::default() };
                let (tc, _) = best(|| {
                    first_contact_with_band(&pool, &s, &p, &vel, 1.0, &[], &o, Some(&band)).unwrap()
                });
                println!(
                    "{},{},{},{:.2},{},{:.2}",
                    if crossing { "crossing" } else { "embedded" },
                    s.n_faces(),
                    threads,
                    1e3 * tb,
                    band.intersecting_pairs(),
                    1e3 * tc
                );
            }
        }
    }
}
