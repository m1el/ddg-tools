//! Micro-benchmark of the eigensolvers: `cargo bench -p ddg-developability --bench eig`.
#![allow(clippy::needless_range_loop)]
use ddg_developability::eig::{sym_eig2, sym_eig3};
use std::hint::black_box;
use std::time::Instant;

fn main() {
    // A = sum theta N N^T for nearly planar normals, plus a generic matrix.
    let mut s = 12345u64;
    let mut rnd = || {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((s >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    };
    let n = 4096;
    let mats: Vec<[[f64; 3]; 3]> = (0..n)
        .map(|_| {
            let mut a = [[0.0; 3]; 3];
            for _ in 0..6 {
                let (x, y, z) = (rnd(), rnd(), 1e-3 * rnd());
                let l = (x * x + y * y + z * z).sqrt();
                let nn = [x / l, y / l, z / l];
                for i in 0..3 {
                    for j in 0..3 {
                        a[i][j] += 1.0 * nn[i] * nn[j];
                    }
                }
            }
            a
        })
        .collect();
    let reps = 500;
    let t = Instant::now();
    let mut acc = 0.0;
    for _ in 0..reps {
        for m in &mats {
            acc += sym_eig3(black_box(*m)).unwrap().values[0];
        }
    }
    let ns = t.elapsed().as_secs_f64() * 1e9 / (reps * n) as f64;
    println!("sym_eig3: {ns:.1} ns/call (checksum {acc:e})");
    let t = Instant::now();
    for _ in 0..reps {
        for m in &mats {
            acc += sym_eig2(black_box([[m[0][0], m[0][1]], [m[0][1], m[1][1]]])).unwrap().values[0];
        }
    }
    let ns = t.elapsed().as_secs_f64() * 1e9 / (reps * n) as f64;
    println!("sym_eig2: {ns:.1} ns/call (checksum {acc:e})");
}
