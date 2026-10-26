//! One component on one curve at given thread counts (`SPEC.md` §14.6, §17k); one CSV row each.
//!
//! ```text
//! cargo bench --bench components -- --curve fhw-r --size 2048 --component gram --threads 1,8,64 [--reps 3]
//! ```
//! With `--component all`, every component. Columns: `best` and `mean` seconds over `reps`
//! runs (after one warm-up run if the first run took under 2 s), `hash` of the result.

mod common;

use common::components::*;
use common::*;
use ddg_repulsive_curves::Exponents;

fn main() {
    let a = args();
    if a.contains_key("header") {
        println!("curve,n,component,threads,reps,best,mean,hash,load");
        return;
    }
    let curve_name = a.get("curve").cloned().unwrap_or("fhw-r".into());
    let size: usize = get(&a, "size", 2048);
    let comp = a.get("component").cloned().unwrap_or("all".into());
    let threads: Vec<usize> = a
        .get("threads")
        .map(|s| s.split(',').filter_map(|x| x.parse().ok()).collect())
        .unwrap_or(vec![1, 8]);
    let reps: usize = get(&a, "reps", 3);
    let (net, pos) = curve(&curve_name, size).expect("curve");
    let inp = Inputs::new(&net, &pos, Exponents::new(2.0, 4.5).unwrap());
    let comps: Vec<&str> = if comp == "all" {
        COMPONENTS.to_vec()
    } else {
        comp.split(',').collect()
    };
    for c in comps {
        let dense = matches!(c, "lu" | "schur").then(|| dense_data(&inp));
        for &t in &threads {
            let mut prep = prepare(c, &inp, t, dense.as_ref());
            let mut times = Vec::new();
            let mut hash = 0;
            let mut err = None;
            for r in 0..reps + 1 {
                match run(c, &inp, t, dense.as_ref(), &mut prep) {
                    Ok((s, h)) => {
                        if r == 0 && s > 2.0 && reps > 0 {
                            // slow: the first run counts
                            times.push(s);
                            hash = h;
                            break;
                        }
                        if r > 0 || reps == 0 {
                            times.push(s);
                        }
                        hash = h;
                    }
                    Err(e) => {
                        err = Some(e);
                        break;
                    }
                }
            }
            if let Some(e) = err {
                eprintln!("{c} t={t}: {e}");
                continue;
            }
            let best = times.iter().cloned().fold(f64::MAX, f64::min);
            let mean = times.iter().sum::<f64>() / times.len() as f64;
            println!(
                "{curve_name},{size},{c},{t},{},{best:.5},{mean:.5},{hash:016x},{:.1}",
                times.len(),
                load_avg()
            );
        }
    }
}
