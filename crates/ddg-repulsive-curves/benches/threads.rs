//! Thread scaling (`SPEC.md` §14.6, §17k): every component at 1, 2, 4, 8, 16, 32, 64 threads on
//! curves of 2k, 8k, 32k vertices (Freedman–He–Wang 2400 resampled by arclength), with
//! efficiency `T₁ / (t · T_t)` and a bitwise check of the results across thread counts.
//! Prints CSV on stdout (the `eff` column) and `MISMATCH` lines on stderr.
//!
//! ```text
//! cargo bench --bench threads -- --sizes 2048,8192,32768 --threads 1,2,4,8,16,32,64
//!     [--components all|a,b] [--reps 3] [--dense-max 8192]
//! ```

mod common;

use common::components::*;
use common::*;
use ddg_repulsive_curves::Exponents;

fn main() {
    let a = args();
    if a.contains_key("header") {
        println!("n,component,threads,reps,best,mean,eff,hash,same_as_1,load");
        return;
    }
    let list = |k: &str, d: &str| -> Vec<usize> {
        a.get(k)
            .map(String::as_str)
            .unwrap_or(d)
            .split(',')
            .filter_map(|x| x.parse().ok())
            .collect()
    };
    let sizes = list("sizes", "2048,8192,32768");
    let threads = list("threads", "1,2,4,8,16,32,64");
    let reps: usize = get(&a, "reps", 3);
    let dense_max: usize = get(&a, "dense-max", 8192);
    let comps: Vec<String> = match a.get("components").map(String::as_str) {
        None | Some("all") => COMPONENTS.iter().map(|s| s.to_string()).collect(),
        Some(s) => s.split(',').map(String::from).collect(),
    };
    for &n in &sizes {
        let (net, pos) = curve("fhw-r", n).expect("curve");
        let inp = Inputs::new(&net, &pos, Exponents::new(2.0, 4.5).unwrap());
        for c in &comps {
            if matches!(c.as_str(), "gram" | "lu") && n > dense_max {
                continue;
            }
            let dense = (c == "lu").then(|| dense_data(&inp));
            let mut t1: Option<(f64, u64)> = None;
            for &t in &threads {
                let mut prep = prepare(c, &inp, t, dense.as_ref());
                let mut times = Vec::new();
                let mut hash = 0;
                let mut failed = false;
                // warm-up run; the number of timed runs depends on its time
                let mut todo = reps;
                let mut r = 0;
                while r <= todo {
                    match run(c, &inp, t, dense.as_ref(), &mut prep) {
                        Ok((s, h)) => {
                            hash = h;
                            if r == 0 {
                                todo = if s > 5.0 {
                                    0
                                } else if s < 0.05 {
                                    reps.max(20)
                                } else if s < 0.5 {
                                    reps.max(10)
                                } else {
                                    reps
                                };
                                if todo == 0 {
                                    times.push(s);
                                }
                            } else {
                                times.push(s);
                            }
                        }
                        Err(e) => {
                            eprintln!("{c} n={n} t={t}: {e}");
                            failed = true;
                            break;
                        }
                    }
                    r += 1;
                }
                if failed || times.is_empty() {
                    continue;
                }
                let best = times.iter().cloned().fold(f64::MAX, f64::min);
                let mean = times.iter().sum::<f64>() / times.len() as f64;
                let base = *t1.get_or_insert((best * t as f64, hash));
                // efficiency relative to the first (smallest) thread count in the list
                let eff = base.0 / (best * t as f64);
                let same = hash == base.1;
                if !same {
                    eprintln!("MISMATCH {c} n={n} t={t}: {hash:016x} vs {:016x}", base.1);
                }
                println!(
                    "{n},{c},{t},{},{best:.5},{mean:.5},{eff:.3},{hash:016x},{same},{:.1}",
                    times.len(),
                    load_avg()
                );
            }
        }
    }
}
