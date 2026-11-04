//! Accuracy and cost of the crate's `exp` (WP4b): `cargo run --release -p ddg-points-as-tori
//! --example exp_check -- [args_in_millions=100] [threads=16] [dump_n]`.
//!
//! Compares `exp::exp` with the host libm (error < 0.52 ulp) over uniform arguments in the blend's
//! range [-800, 0] and in [-745, 709], and times both single threaded on an array of blend-range
//! arguments. With `dump_n`, writes `dump_n` arguments and our results as hex to stdout
//! (for a high-precision cross-check).

use ddg_points_as_tori::exp::exp;
use rayon::prelude::*;

fn lcg(s: &mut u64) -> f64 {
    *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    (*s >> 11) as f64 / (1u64 << 53) as f64
}

fn ulp_diff(a: f64, b: f64) -> u64 {
    if a == b {
        0
    } else {
        (a.to_bits() as i64 - b.to_bits() as i64).unsigned_abs()
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let millions: usize = a.get(1).map(|s| s.parse().unwrap()).unwrap_or(100);
    let threads: usize = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(16);
    if let Some(n) = a.get(3) {
        let n: usize = n.parse().unwrap();
        let mut s = 99;
        for _ in 0..n {
            let x = lcg(&mut s) * 1500.0 - 740.0;
            println!("{:016x} {:016x}", x.to_bits(), exp(x).to_bits());
        }
        return;
    }
    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
    for (name, lo, hi) in [("blend range [-800, 0]", -800.0, 0.0), ("full range [-745, 709.78]", -745.0, 709.78)] {
        let total = millions * 1_000_000;
        let chunks = 1000;
        let (maxu, hist, nd) = pool.install(|| {
            (0..chunks)
                .into_par_iter()
                .map(|c| {
                    let mut s = 1000 + c as u64;
                    let mut maxu = 0;
                    let mut hist = [0u64; 4];
                    for _ in 0..total / chunks {
                        let x = lo + (hi - lo) * lcg(&mut s);
                        let u = ulp_diff(exp(x), x.exp());
                        maxu = maxu.max(u);
                        hist[u.min(3) as usize] += 1;
                    }
                    (maxu, hist, 0u64)
                })
                .reduce(|| (0, [0; 4], 0), |a, b| (a.0.max(b.0), [a.1[0] + b.1[0], a.1[1] + b.1[1], a.1[2] + b.1[2], a.1[3] + b.1[3]], a.2 + b.2))
        });
        let _ = nd;
        println!("{name}: {total} arguments, max |difference to libm| {maxu} ulp; differing by 0/1/2/>=3 ulp: {:?}", hist);
    }
    // edge cases
    for x in [0.0, -0.0, 1e-300, 709.782712893384, 709.79, -745.1332191019412, -745.14, -708.3964185322641, -708.4, f64::NEG_INFINITY, f64::INFINITY, f64::NAN] {
        println!("exp({x:e}) = {:e} (libm {:e})", exp(x), x.exp());
    }
    // cost, single threaded: 1M blend-range arguments, summed so nothing is optimized away
    let mut s = 5;
    let xs: Vec<f64> = (0..1_000_000).map(|_| -128.0 * lcg(&mut s)).collect();
    for rep in 0..3 {
        let t0 = std::time::Instant::now();
        let mut acc = 0.0;
        for _ in 0..20 {
            for &x in &xs {
                acc += std::hint::black_box(x).exp();
            }
        }
        let t_libm = t0.elapsed().as_secs_f64() / 20e6 * 1e9;
        let t0 = std::time::Instant::now();
        let mut acc2 = 0.0;
        for _ in 0..20 {
            for &x in &xs {
                acc2 += exp(std::hint::black_box(x));
            }
        }
        let t_own = t0.elapsed().as_secs_f64() / 20e6 * 1e9;
        println!("rep {rep}: libm {t_libm:.2} ns/exp, crate {t_own:.2} ns/exp (sums {acc:.6e} {acc2:.6e})");
    }
}
