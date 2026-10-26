//! Runs a scene of the applications (`SPEC.md` §14.4) and writes OBJ frames and a JSON report.
//!
//! ```text
//! scenes <scene> [--out DIR] [--steps N] [--backend auto|exact|bh|mg|emg] [--threads T]
//!                [--every K] [--audit-every K] [--size N] [--curve FILE]
//! ```
//!
//! `<scene>` is a scene file of the authors' format (`scenes/FORMATS.md`), a directory with a
//! `scene.txt`, or `builtin:<name>` (see `scene::builtin::NAMES`). `--curve FILE` replaces the
//! scene's curve (e.g. the 600-vertex Freedman curve under the `FreedmanUnknot` scene).
//! Output: `DIR/frame00000.obj ...` (every K steps, the first and the last always) and
//! `DIR/report.json`.

use std::path::PathBuf;

use ddg_repulsive_curves::io::read_curve_obj;
use ddg_repulsive_curves::scene::{
    BackendChoice, RunOptions, Scene, builtin, read_scene, run_scene,
};

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let target = args.next().ok_or("usage: scenes <scene> [options]")?;
    let mut o = RunOptions::default();
    let mut size = 0usize;
    let mut curve: Option<PathBuf> = None;
    while let Some(a) = args.next() {
        let mut val = || args.next().ok_or(format!("{a} needs a value"));
        match a.as_str() {
            "--out" => o.out_dir = Some(PathBuf::from(val()?)),
            "--steps" => o.steps = Some(val()?.parse()?),
            "--backend" => {
                let v = val()?;
                o.backend = BackendChoice::parse(&v).ok_or("backend: auto, exact, bh or mg")?;
            }
            "--threads" => o.threads = Some(val()?.parse()?),
            "--every" => o.frame_every = val()?.parse()?,
            "--audit-every" => o.audit_every = val()?.parse()?,
            "--size" => size = val()?.parse()?,
            "--curve" => curve = Some(PathBuf::from(val()?)),
            other => return Err(format!("unknown option {other}").into()),
        }
    }
    let mut scene: Scene = if let Some(name) = target.strip_prefix("builtin:") {
        builtin::by_name(name, size)?
    } else {
        let mut p = PathBuf::from(&target);
        if p.is_dir() {
            p = p.join("scene.txt");
        }
        read_scene(&p)?
    };
    if let Some(c) = curve {
        let (net, pos) = read_curve_obj(&c)?;
        scene = Scene::new(scene.name.clone(), net, pos, scene.spec.clone())?;
    }
    eprintln!(
        "scene {}: {} vertices, {} edges",
        scene.name,
        scene.network.n_vertices(),
        scene.network.n_edges()
    );
    for n in &scene.notes {
        eprintln!("note: {n}");
    }
    let out = run_scene(&scene, &o)?;
    let r = &out.report;
    eprintln!(
        "{} steps ({}), {:.1} s, vertices {} -> {}, energy {:.6e} -> {:.6e}, length {:.4} -> {:.4}, \
         min edge distance {:.3e} -> {:.3e} (over run {:.3e}), residual {:?}",
        r.steps,
        r.status,
        r.seconds,
        r.vertices.0,
        r.vertices.1,
        r.energy.0,
        r.energy.1,
        r.length.0,
        r.length.1,
        r.audit_start.min_edge_distance,
        r.audit_end.min_edge_distance,
        r.min_distance_over_run,
        r.audit_end.max_constraint_residual
    );
    if o.out_dir.is_none() {
        println!("{}", r.to_json());
    }
    Ok(())
}
