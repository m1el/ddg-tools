//! Per-step time of the flow with its breakdown (`SPEC.md` §14.5, Fig. 13; WP17). A plain `main`
//! printing one CSV row per mesh (`--header` prints the column names first).
//!
//! ```text
//! cargo bench -p ddg-repulsive-surfaces --bench step -- [--header] [--suite fig13|wp18b]
//!     [--mesh SPEC] [--threads N] [--warmup N] [--steps N] [--compat] [--cache DIR]
//!     [--gmres-tol X] [--chi X]
//! ```
//!
//! `SPEC`: `genus:G:FACES`, `handcuffs:K:FACES`, `tunnel8:FACES`, `tunnel3:FACES`, `fig12:I`,
//! `punctured:FACES`, `torus:FACES`, `obj:PATH` (barycenter + area constraints, as the paper's
//! runs), or `scene:PATH` (the scene's own constraints and penalties). `--suite fig13` runs the 14
//! Fig. 13 analogues (`fig13_cases`), `--suite wp18b` the scenes of the black-box speed comparison
//! (`COMPARISON.md` §2.3). Generated meshes are cached as OBJ in `--cache DIR` if given.
//!
//! Each case runs `warmup` steps, then `steps` timed steps (remeshing and CCD on, the stop rule
//! off). The columns are means over the timed steps (`*_ms`), the median step, and the counts.
//! Run alone: the numbers are wall times.

use std::path::{Path, PathBuf};

use ddg_repulsive_surfaces::generate::Mesh;
use ddg_repulsive_surfaces::generate::primitives::{punctured_torus, torus};
use ddg_repulsive_surfaces::generate::shapes::{Knot, fig12, genus_surface, handcuffs, knotted_tunnel};
use ddg_repulsive_surfaces::scene::{SceneFile, SceneSpec};
use ddg_repulsive_surfaces::{
    Backend, Constraint, FlowOptions, Session, StepReport, StopRule, Surface, read_obj, scene_options, write_obj,
};

const SCENES: &str = "data/repulsive-surfaces/scenes";

struct Args {
    header: bool,
    suite: Option<String>,
    meshes: Vec<String>,
    threads: usize,
    warmup: usize,
    steps: usize,
    compat: bool,
    cache: Option<PathBuf>,
    gmres_tol: Option<f64>,
    chi: Option<f64>,
    restart: Option<usize>,
    no_lo: bool,
}

fn parse_args() -> Args {
    let mut a = Args {
        header: false,
        suite: None,
        meshes: Vec::new(),
        threads: 4,
        warmup: 2,
        steps: 5,
        compat: false,
        cache: None,
        gmres_tol: None,
        chi: None,
        restart: None,
        no_lo: false,
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let next = |i: usize| argv.get(i + 1).cloned().unwrap_or_else(|| panic!("{} needs a value", argv[i]));
        match argv[i].as_str() {
            "--header" => a.header = true,
            "--compat" => a.compat = true,
            "--bench" => {}
            "--suite" => {
                a.suite = Some(next(i));
                i += 1;
            }
            "--mesh" => {
                a.meshes.push(next(i));
                i += 1;
            }
            "--threads" => {
                a.threads = next(i).parse().expect("--threads N");
                i += 1;
            }
            "--warmup" => {
                a.warmup = next(i).parse().expect("--warmup N");
                i += 1;
            }
            "--steps" => {
                a.steps = next(i).parse().expect("--steps N");
                i += 1;
            }
            "--cache" => {
                a.cache = Some(PathBuf::from(next(i)));
                i += 1;
            }
            "--gmres-tol" => {
                a.gmres_tol = Some(next(i).parse().expect("--gmres-tol X"));
                i += 1;
            }
            "--no-lo" => a.no_lo = true,
            "--restart" => {
                a.restart = Some(next(i).parse().expect("--restart N"));
                i += 1;
            }
            "--chi" => {
                a.chi = Some(next(i).parse().expect("--chi X"));
                i += 1;
            }
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    a
}

/// The 14 models of Fig. 13 (paper face counts in brackets) as analogues: generated surfaces
/// (WP12) and reference scenes at about the same size.
fn fig13_cases() -> Vec<(&'static str, String)> {
    let s = |p: &str| format!("scene:{SCENES}/{p}");
    vec![
        ("01 triple torus [6k]", "genus:3:6000".into()),
        ("02 moving bars [9k] ~ HandcuffsPole scene", s("HandcuffsPole/scene.txt")),
        ("03 broken figure-8 [10k]", "obj:".to_string() + SCENES + "/BrokenFigureEight/BrokenFigureEightTubeMesh_11740F.obj"),
        ("04 punctured torus [10k]", "punctured:10000".into()),
        ("05 handcuffs [11k]", "handcuffs:1:11000".into()),
        ("06 2x-wound handcuffs [14k]", "handcuffs:2:16000".into()),
        ("07 point cloud [20k] ~ sphere20k + bunny cloud", "pointcloud".into()),
        ("08 nested cages [20k] ~ sphere20k + mesh obstacle", "nested".into()),
        ("09 walnut [23k] ~ genus-0 slab", "genus:0:23000".into()),
        ("10 trefoil tunnel [23k]", "fig12:4".into()),
        ("11 6x-wound handcuffs [22k]", "fig12:2".into()),
        ("12 double torus [36k]", s("EntangledDoubleTorus/scene-bench.txt")),
        ("13 keys [37k] ~ genus-4 slab", "genus:4:37000".into()),
        ("14 figure-8 tunnel [45k]", "fig12:3".into()),
    ]
}

/// The scenes of the 1-thread comparison with the reference (`COMPARISON.md` §2.3).
fn wp18b_cases() -> Vec<(&'static str, String)> {
    let s = |p: &str| format!("scene:{SCENES}/{p}");
    vec![
        ("EDT", s("EntangledDoubleTorus/scene.txt")),
        ("Hopf", s("HopfLink/scene.txt")),
        ("Handcuffs", s("HandcuffsPole/scene.txt")),
        ("Walnut", s("Walnut/scene.txt")),
        ("Lung", s("LungGrowing/scene.txt")),
        ("Bunny", s("Bunny/scene.txt")),
        ("Recon", s("Reconstruction/scene_bunny.txt")),
        ("kbb", s("KleinBottle/scene.txt")),
    ]
}

fn generated(spec: &str) -> Mesh {
    let p: Vec<&str> = spec.split(':').collect();
    let n = |i: usize| -> usize { p[i].parse().expect("a number") };
    match p[0] {
        "genus" => genus_surface(n(1) as u32, n(2)).expect("genus surface"),
        "handcuffs" => handcuffs(n(1) as u32, n(2)).expect("handcuffs").mesh,
        "tunnel8" => knotted_tunnel(Knot::FigureEight, n(1)).expect("tunnel").mesh,
        "tunnel3" => knotted_tunnel(Knot::Trefoil, n(1)).expect("tunnel").mesh,
        "fig12" => fig12(n(1)).expect("fig12").mesh,
        "torus" | "punctured" => {
            // n_major ≈ 2.5 n_minor, 2 n_major n_minor faces.
            let f = n(1) as f64;
            let nm = (f / 5.0).sqrt().round() as usize;
            let nmaj = (f / (2.0 * nm as f64)).round() as usize;
            if p[0] == "torus" {
                torus(nmaj, nm, 1.5, 0.6).expect("torus")
            } else {
                punctured_torus(nmaj, nm, 1.5, 0.6, 0.4).expect("punctured torus")
            }
        }
        other => panic!("unknown mesh kind {other}"),
    }
}

fn cached(spec: &str, cache: Option<&Path>) -> (Surface, Vec<[f64; 3]>) {
    let file = cache.map(|d| d.join(format!("{}.obj", spec.replace(':', "_"))));
    if let Some(f) = &file
        && f.exists()
    {
        return read_obj(f).expect("cached OBJ");
    }
    let m = generated(spec);
    let s = m.surface().expect("generated surface");
    if let Some(f) = &file {
        std::fs::create_dir_all(f.parent().unwrap()).ok();
        write_obj(f, &s, &m.positions).ok();
    }
    (s, m.positions)
}

fn paper_options() -> FlowOptions {
    FlowOptions {
        constraints: vec![
            Constraint::Barycenter { component: None, target: None },
            Constraint::Area { component: None, target: None },
        ],
        ..FlowOptions::default()
    }
}

fn scene_text(dir: &str, text: &str) -> (Surface, Vec<[f64; 3]>, FlowOptions) {
    let path = Path::new(SCENES).join(dir).join("wp17-bench.txt");
    let file = SceneFile::parse(text, &path).expect("scene text");
    let spec = SceneSpec::resolve(file).expect("scene");
    let o = scene_options(&spec).expect("scene options");
    (spec.surface, spec.positions, o)
}

fn setup(spec: &str, cache: Option<&Path>) -> (Surface, Vec<[f64; 3]>, FlowOptions) {
    if let Some(path) = spec.strip_prefix("scene:") {
        let sc = SceneSpec::load(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let o = scene_options(&sc).expect("scene options");
        return (sc.surface, sc.positions, o);
    }
    if let Some(path) = spec.strip_prefix("obj:") {
        let (s, p) = read_obj(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        return (s, p, paper_options());
    }
    match spec {
        "pointcloud" => scene_text(
            "Reconstruction",
            "repel_mesh sphere20k.obj 6 12\nconstrain volume 0.11 200\nconstrain barycenter\n\
             point_cloud_obstacle bunny-scaled.obj 1\nallow_barycenter_shift\n",
        ),
        "nested" => scene_text(
            "Reconstruction",
            "repel_mesh sphere20k.obj 6 12\nconstrain barycenter\nconstrain area\n\
             obstacle bunny-scaled.obj 0.1\n",
        ),
        _ => {
            let (s, p) = cached(spec, cache);
            (s, p, paper_options())
        }
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

const HEADER: &str = "case,spec,threads,faces0,faces_end,steps,step_ms,step_med_ms,bvh_ms,energy_ms,factor_ms,\
bct_ms,fractional_ms,lprime_ms,solve_ms,products,product_ms,precond_ms,gmres_its,solves,line_search_ms,trials,\
ccd_ms,projection_ms,remesh_ms,max_residual,energy_end,status";

fn run_case(name: &str, spec: &str, a: &Args) {
    let (surface, positions, mut o) = setup(spec, a.cache.as_deref());
    let faces0 = surface.n_faces();
    o.threads = Some(a.threads);
    o.stop = StopRule::Off;
    o.max_iterations = usize::MAX;
    o.reference_compat = a.compat;
    if let Backend::Auto { hierarchical: h, .. } | Backend::Hierarchical(h) = &mut o.backend {
        if let Some(t) = a.gmres_tol {
            h.gmres.tolerance = t;
        }
        if let Some(c) = a.chi {
            h.cluster.chi = c;
        }
        if let Some(r) = a.restart {
            h.gmres.restart = r;
        }
    }
    if a.no_lo {
        o.metric = ddg_repulsive_surfaces::Metric::Sobolev { low_order: false };
    }
    let mut session = Session::new(surface, positions, o).unwrap_or_else(|e| panic!("{name}: {e}"));
    let mut reps: Vec<StepReport> = Vec::new();
    let mut status = "ok".to_string();
    for k in 0..a.warmup + a.steps {
        match session.step() {
            Ok(r) => {
                if std::env::var_os("WP17_HISTORY").is_some() {
                    for (i, s) in r.linear.solves.iter().enumerate() {
                        if let Some(g) = &s.gmres {
                            let h: Vec<String> = g.history.iter().map(|x| format!("{x:.1e}")).collect();
                            eprintln!(
                                "step {} solve {i}: {} its, true {:.2e}, cycles {:?}: {}",
                                k + 1,
                                g.iterations,
                                g.relative_residual,
                                g.cycle_residuals,
                                h.join(" ")
                            );
                        }
                    }
                }
                if k >= a.warmup {
                    reps.push(r);
                }
            }
            Err(e) => {
                status = format!("error at step {}: {e}", k + 1).replace(',', ";");
                break;
            }
        }
    }
    let n = reps.len().max(1) as f64;
    let mean = |f: &dyn Fn(&StepReport) -> f64| 1e3 * reps.iter().map(f).sum::<f64>() / n;
    let cnt = |f: &dyn Fn(&StepReport) -> f64| reps.iter().map(f).sum::<f64>() / n;
    let s = |r: &StepReport| r.seconds;
    println!(
        "{name},{spec},{},{faces0},{},{},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.2},{:.1},{:.1},{:.1},{:.2e},{:.17e},{status}",
        a.threads,
        session.surface().n_faces(),
        reps.len(),
        mean(&|r| s(r).total),
        1e3 * median(reps.iter().map(|r| s(r).total).collect()),
        mean(&|r| s(r).bvh),
        mean(&|r| s(r).energy),
        mean(&|r| s(r).factor),
        mean(&|r| s(r).metric_tree),
        mean(&|r| s(r).fractional_tree),
        mean(&|r| s(r).lprime),
        mean(&|r| s(r).solve),
        cnt(&|r| r.products.map_or(0.0, |p| p.0 as f64)),
        mean(&|r| r.products.map_or(0.0, |p| p.1)),
        mean(&|r| s(r).preconditioner),
        cnt(&|r| r.linear.total_iterations() as f64),
        cnt(&|r| r.linear.solves.len() as f64),
        mean(&|r| s(r).line_search),
        cnt(&|r| r.trials.len() as f64),
        mean(&|r| s(r).ccd),
        mean(&|r| s(r).projection),
        mean(&|r| s(r).remesh),
        reps.iter().map(|r| r.linear.max_relative_residual()).fold(0.0, f64::max),
        session.energy(),
    );
}

fn main() {
    let a = parse_args();
    if a.header {
        println!("{HEADER}");
    }
    let mut cases: Vec<(String, String)> = match a.suite.as_deref() {
        Some("fig13") => fig13_cases().into_iter().map(|(n, s)| (n.to_string(), s)).collect(),
        Some("wp18b") => wp18b_cases().into_iter().map(|(n, s)| (n.to_string(), s)).collect(),
        Some(other) => panic!("unknown suite {other}"),
        None => Vec::new(),
    };
    cases.extend(a.meshes.iter().map(|m| (m.clone(), m.clone())));
    for (name, spec) in &cases {
        run_case(name, spec, &a);
    }
}
