//! Tracer `safety` / `g_min` study (SPEC.md §16.3, D17): for one cloud, render at 256² with each
//! setting and compare with a reference render (safety 0.5, 2048 steps): hits that differ from the
//! reference (an overshoot or a missed thin part), step statistics and step-limit rays.
//!
//! ```text
//! cargo run --release -p ddg-points-as-tori --example render_study -- <cloud.ply> [threads]
//! ```

use ddg_points_as_tori::render::{Camera, TraceOptions, Tracer};
use ddg_points_as_tori::{BuildOptions, Model, NormalPolicy, PointCloud, PointsAsTori};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let path = args
        .get(1)
        .ok_or("usage: render_study <cloud.ply> [threads]")?;
    let threads: usize = args.get(2).map(|s| s.parse()).transpose()?.unwrap_or(8);
    let bytes = std::fs::read(path)?;
    let end = bytes
        .windows(10)
        .position(|w| w == b"end_header")
        .ok_or("no end_header")?;
    let header = std::str::from_utf8(&bytes[..end])?;
    let n: usize = header
        .lines()
        .find_map(|l| l.strip_prefix("element vertex "))
        .ok_or("no vertex element")?
        .trim()
        .parse()?;
    let start = bytes[end..]
        .iter()
        .position(|b| *b == b'\n')
        .ok_or("bad header")?
        + end
        + 1;
    let f = |i: usize| {
        f32::from_le_bytes(bytes[start + 4 * i..start + 4 * i + 4].try_into().unwrap()) as f64
    };
    let mut p: Vec<[f64; 3]> = (0..n)
        .map(|i| [f(6 * i), f(6 * i + 1), f(6 * i + 2)])
        .collect();
    let nr: Vec<[f64; 3]> = (0..n)
        .map(|i| [f(6 * i + 3), f(6 * i + 4), f(6 * i + 5)])
        .collect();
    let mean = [0, 1, 2].map(|a| p.iter().map(|x| x[a]).sum::<f64>() / n as f64);
    let scale = p
        .iter()
        .flat_map(|x| (0..3).map(move |a| (x[a] - mean[a]).abs()))
        .fold(0.0, f64::max);
    for x in p.iter_mut() {
        for a in 0..3 {
            x[a] = (x[a] - mean[a]) / scale;
        }
    }
    let (cloud, _) = PointCloud::new(p, nr, NormalPolicy::Normalize)?;
    let (field, _) = PointsAsTori::build(
        &cloud,
        &Model::reference()?,
        &BuildOptions::default().with_threads(Some(threads)),
    )?;
    let cam = Camera::demo([0.0; 3], 1.0, 256, 256)?;
    let render = |safety: f64, g_min: f64, steps: u32| -> Result<_, Box<dyn std::error::Error>> {
        let o = TraceOptions {
            safety,
            g_min,
            max_steps: steps,
            shadow_max_steps: steps,
            threads: Some(threads),
            ..Default::default()
        };
        Ok(Tracer::new(o)?.render(&field, &cam)?)
    };
    let (rimg, rrep) = render(0.5, 0.1, 2048)?;
    println!(
        "{path}: reference (safety 0.5, g_min 0.1, 2048 steps): hits {}, mean steps {:.2}, max {}, step-limit {}",
        rrep.hits, rrep.steps_mean, rrep.steps_max, rrep.step_limit
    );
    println!(
        "safety g_min steps | hits | alpha-mismatch px | rgb-mismatch px | mean steps | max steps | step-limit | seconds"
    );
    for (safety, g_min, steps) in [
        (0.5, 0.1, 128),
        (0.7, 0.1, 128),
        (0.9, 0.1, 128),
        (1.0, 0.1, 128),
        (0.9, 0.0, 128),
        (0.9, 0.3, 128),
        (0.9, 0.1, 64),
        (0.9, 0.1, 40),
        (1.0, 0.0, 40),
    ] {
        let (img, rep) = render(safety, g_min, steps)?;
        let am = img
            .alpha
            .iter()
            .zip(&rimg.alpha)
            .filter(|(a, b)| a != b)
            .count();
        let cm = img
            .rgb
            .chunks_exact(3)
            .zip(rimg.rgb.chunks_exact(3))
            .filter(|(a, b)| a != b)
            .count();
        println!(
            "{safety:>6} {g_min:>5} {steps:>5} | {:>6} | {am:>6} | {cm:>6} | {:>6.2} | {:>4} | {:>5} | {:.2}",
            rep.hits, rep.steps_mean, rep.steps_max, rep.step_limit, rep.seconds
        );
    }
    Ok(())
}
