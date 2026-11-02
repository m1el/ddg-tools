//! Applications of Sec. 9 of "Repulsive Surfaces" (Yu, Brakensiek, Schumacher, Crane 2021).
//!
//! ```text
//! cargo run --release -p ddg-repulsive-surfaces --example apps -- list
//! cargo run --release -p ddg-repulsive-surfaces --example apps -- NAME [--out DIR] [--threads T]
//!        [--frame-every K] [--steps N] [--smoke]
//! ```
//! Prints `RESULT ...` lines (steps, time, energy, area/volume, constraint violation,
//! self-intersection) and writes OBJ frames to `DIR/NAME/` when `--out` is given.
mod apps;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("list") || args.is_empty() {
        for n in apps::all_names() {
            println!("{n}");
        }
        return;
    }
    let name = args[0].clone();
    let mut cfg = apps::Cfg { smoke: false, threads: 4, out: None, frame_every: 0, steps: None, scenes: apps::scenes_dir() };
    let mut i = 1;
    while i < args.len() {
        let v = |i: usize| args.get(i + 1).cloned().expect("option needs a value");
        match args[i].as_str() {
            "--out" => {
                cfg.out = Some(v(i).into());
                i += 1;
            }
            "--threads" => {
                cfg.threads = v(i).parse().unwrap();
                i += 1;
            }
            "--frame-every" => {
                cfg.frame_every = v(i).parse().unwrap();
                i += 1;
            }
            "--steps" => {
                cfg.steps = Some(v(i).parse().unwrap());
                i += 1;
            }
            "--smoke" => cfg.smoke = true,
            o => panic!("unknown option {o}"),
        }
        i += 1;
    }
    match apps::by_name(&cfg, &name) {
        Ok(v) => {
            for o in v {
                println!("{}", o.line());
            }
        }
        Err(e) => {
            println!("RESULT name={name} status=ERROR error={e}");
            std::process::exit(1);
        }
    }
}
