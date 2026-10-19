//! Command line driver of the developability flow (SPEC §5.1) and of cut and flatten (SPEC §6).
//!
//! ```text
//! cargo run --release -p ddg-developability --example develop -- IN.obj [options]
//!   --energy NAME        covariance (default) | combinatorial | combinatorial-max |
//!                        covariance-max | intrinsic | intrinsic-max
//!   --subdivisions N     4-1 subdivision rounds after the first (default 1: two rounds)
//!   --iterations N       iterations per round (default 200)
//!   --no-remesh          switch small-angle remeshing off
//!   --small-angle DEG    remeshing threshold in degrees (default 1)
//!   --max-displacement X step cap in mean edge lengths (default 0.5)
//!   --threads N          threads of the flow (default: all)
//!   --out FILE.obj       write the developable mesh; the per-vertex energy goes to FILE.energy.txt
//!   --cut-out FILE.obj   cut and flatten the result, write an OBJ with per-corner UVs (charts unpacked)
//!   --epsilon X          terminal threshold of the cut (default 1e-3)
//!   --cut-only           skip the flow, cut and flatten the input
//!   --no-cut             skip the cut and flatten stage (by default it always runs)
//! ```
//!
//! Polygons of the input are fan-triangulated. The last line of stdout is a `RESULT key=value ...`
//! summary for scripts.

use std::process::ExitCode;
use std::time::Instant;

use ddg_core::TriMesh;
use ddg_core::io::{read_obj, write_obj};
use ddg_developability::{
    CutOptions, Energy, EnergyOptions, FlowOptions, OptimizeStatus, RemeshOptions, cut_and_flatten, evaluate_energy,
    flow,
};

fn parse_energy(s: &str) -> Option<Energy> {
    Some(match s {
        "covariance" | "e_lambda" => Energy::Covariance,
        "combinatorial" | "e_p" => Energy::Combinatorial,
        "combinatorial-max" => Energy::CombinatorialMax,
        "covariance-max" => Energy::CovarianceMax,
        "intrinsic" => Energy::Intrinsic,
        "intrinsic-max" => Energy::IntrinsicMax,
        _ => return None,
    })
}

struct Args {
    input: String,
    energy: Energy,
    subdivisions: usize,
    iterations: usize,
    remesh: bool,
    small_angle_deg: f64,
    max_displacement: f64,
    threads: Option<usize>,
    out: Option<String>,
    cut_out: Option<String>,
    epsilon: f64,
    cut_only: bool,
    no_cut: bool,
}

fn parse_args() -> Result<Args, String> {
    let mut it = std::env::args().skip(1);
    let mut a = Args {
        input: String::new(),
        energy: Energy::Covariance,
        subdivisions: 1,
        iterations: 200,
        remesh: true,
        small_angle_deg: 1.0,
        max_displacement: 0.5,
        threads: None,
        out: None,
        cut_out: None,
        epsilon: 1e-3,
        cut_only: false,
        no_cut: false,
    };
    while let Some(arg) = it.next() {
        let mut val = |name: &str| it.next().ok_or_else(|| format!("{name} needs a value"));
        match arg.as_str() {
            "--energy" => {
                let v = val("--energy")?;
                a.energy = parse_energy(&v).ok_or_else(|| format!("unknown energy {v}"))?;
            }
            "--subdivisions" => a.subdivisions = val("--subdivisions")?.parse().map_err(|e| format!("{e}"))?,
            "--iterations" => a.iterations = val("--iterations")?.parse().map_err(|e| format!("{e}"))?,
            "--no-remesh" => a.remesh = false,
            "--small-angle" => a.small_angle_deg = val("--small-angle")?.parse().map_err(|e| format!("{e}"))?,
            "--max-displacement" => a.max_displacement = val("--max-displacement")?.parse().map_err(|e| format!("{e}"))?,
            "--threads" => a.threads = Some(val("--threads")?.parse().map_err(|e| format!("{e}"))?),
            "--out" => a.out = Some(val("--out")?),
            "--cut-out" => a.cut_out = Some(val("--cut-out")?),
            "--epsilon" => a.epsilon = val("--epsilon")?.parse().map_err(|e| format!("{e}"))?,
            "--cut-only" => a.cut_only = true,
            "--no-cut" => a.no_cut = true,
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            s => a.input = s.to_string(),
        }
    }
    if a.input.is_empty() {
        return Err("usage: develop IN.obj [options] (see the source header)".into());
    }
    Ok(a)
}

fn status_name(s: OptimizeStatus) -> &'static str {
    match s {
        OptimizeStatus::GradientTolerance => "gradient",
        OptimizeStatus::EnergyStalled => "stalled",
        OptimizeStatus::MaxIterations => "max_iter",
        OptimizeStatus::LineSearchFailed => "linesearch_failed",
        OptimizeStatus::Cancelled => "cancelled",
        OptimizeStatus::NonDescent => "nondescent",
    }
}

fn run() -> Result<(), String> {
    let a = parse_args()?;
    let obj = read_obj(&a.input).map_err(|e| format!("{}: {e}", a.input))?;
    let positions = obj.positions.clone();
    let mut faces: Vec<[u32; 3]> = Vec::new();
    for f in 0..obj.n_faces() {
        let r = obj.face_corners(f);
        let c = &obj.corner_verts[r];
        for k in 1..c.len() - 1 {
            faces.push([c[0], c[k], c[k + 1]]);
        }
    }
    let mesh = TriMesh::new(positions.len(), faces).map_err(|e| format!("{}: {e}", a.input))?;
    let mut summary: Vec<(String, String)> = vec![
        ("mesh".into(), a.input.rsplit('/').next().unwrap_or(&a.input).to_string()),
        ("energy".into(), format!("{:?}", a.energy)),
        ("faces_in".into(), mesh.n_faces().to_string()),
        ("genus".into(), mesh.genus().to_string()),
        ("boundary_loops".into(), mesh.boundary_loops().len().to_string()),
    ];
    let (mut cur_mesh, mut cur_pos) = (mesh.clone(), positions.clone());
    let eopts = EnergyOptions::new(a.energy);

    if !a.cut_only {
        let mut fo = FlowOptions::uniform(eopts.clone(), a.subdivisions, a.iterations);
        fo.optimize.max_displacement = a.max_displacement;
        fo.remesh = a.remesh.then(|| RemeshOptions {
            small_angle: a.small_angle_deg.to_radians(),
            ..Default::default()
        });
        fo.threads = a.threads;
        let t = Instant::now();
        let res = flow(&mesh, &positions, &fo, None).map_err(|e| format!("flow failed: {e}"))?;
        let secs = t.elapsed().as_secs_f64();
        let mut flips = 0;
        let mut collapses = 0;
        let mut iters = 0;
        let mut evals = 0;
        let mut undefined = 0;
        let mut resets = 0;
        for (i, r) in res.report.rounds.iter().enumerate() {
            let o = &r.optimize;
            println!(
                "round {i}: {} v / {} f, E {:.6e} -> {:.6e}, {} iterations in {} segments ({}), {} evals, remesh {:?}, {:.2}s",
                r.n_vertices,
                r.n_faces,
                o.energy_initial,
                o.energy_final,
                o.iterations,
                r.segments.len(),
                status_name(o.status),
                o.evaluations,
                r.remesh.iter().map(|m| (m.flips, m.collapses)).collect::<Vec<_>>(),
                r.seconds_optimize + r.seconds_remesh + r.seconds_subdivide
            );
            flips += r.remesh.iter().map(|m| m.flips).sum::<usize>();
            collapses += r.remesh.iter().map(|m| m.collapses).sum::<usize>();
            iters += o.iterations;
            evals += o.evaluations;
            undefined += o.undefined_trials;
            resets += o.memory_resets;
        }
        let first = &res.report.rounds[0];
        let last = res.report.rounds.last().expect("at least one round");
        let (e_final, ve) = evaluate_energy(&res.mesh, &res.positions, &eopts).map_err(|e| format!("final energy: {e}"))?;
        summary.extend([
            ("e_initial".into(), format!("{:.6e}", first.optimize.energy_initial)),
            ("e_final".into(), format!("{:.6e}", e_final)),
            ("seconds".into(), format!("{secs:.3}")),
            ("status".into(), res.report.rounds.iter().map(|r| status_name(r.optimize.status)).collect::<Vec<_>>().join("/")),
            ("iterations".into(), iters.to_string()),
            ("evaluations".into(), evals.to_string()),
            ("undefined_trials".into(), undefined.to_string()),
            ("memory_resets".into(), resets.to_string()),
            ("flips".into(), flips.to_string()),
            ("collapses".into(), collapses.to_string()),
            ("faces_out".into(), last.n_faces.to_string()),
            ("finite".into(), (e_final.is_finite() && res.positions.iter().all(|p| p.iter().all(|c| c.is_finite()))).to_string()),
        ]);
        if let Some(out) = &a.out {
            let cv: Vec<u32> = res.mesh.faces().iter().flatten().copied().collect();
            let fs: Vec<u32> = (0..res.mesh.n_faces() as u32).map(|f| 3 * f).collect();
            write_obj(out, &res.positions, &cv, &fs, None).map_err(|e| format!("{out}: {e}"))?;
            let txt: String = ve.iter().map(|e| format!("{e:e}\n")).collect();
            let ep = format!("{}.energy.txt", out.trim_end_matches(".obj"));
            std::fs::write(&ep, txt).map_err(|e| format!("{ep}: {e}"))?;
        }
        cur_mesh = res.mesh;
        cur_pos = res.positions;
    }

    if !a.no_cut {
        let t = Instant::now();
        let co = CutOptions { energy: a.energy, epsilon: a.epsilon, threads: a.threads, ..Default::default() };
        match cut_and_flatten(&cur_mesh, &cur_pos, &co) {
            Ok(cf) => {
                let secs = t.elapsed().as_secs_f64();
                let overlapped: f64 = cf.charts.iter().filter_map(|c| c.overlap.overlapped_fraction()).fold(0.0, f64::max);
                println!(
                    "cut: {} terminals (+{} added), {} cut edges, {} charts, {} flipped uv triangles, mean |u| {:.4e}",
                    cf.terminals.len(),
                    cf.added_terminals.len(),
                    cf.cut_edges.len(),
                    cf.n_charts(),
                    cf.flipped_triangles(),
                    cf.mean_abs_log_scale
                );
                summary.extend([
                    ("cut".into(), "ok".into()),
                    ("terminals".into(), cf.terminals.len().to_string()),
                    ("added_terminals".into(), cf.added_terminals.len().to_string()),
                    ("cut_edges".into(), cf.cut_edges.len().to_string()),
                    ("charts".into(), cf.n_charts().to_string()),
                    ("flipped_uv".into(), cf.flipped_triangles().to_string()),
                    ("injective".into(), cf.is_injective().to_string()),
                    ("max_chart_overlap".into(), format!("{overlapped:.4}")),
                    ("mean_abs_u".into(), format!("{:.6e}", cf.mean_abs_log_scale)),
                    ("cut_seconds".into(), format!("{secs:.3}")),
                ]);
                if let Some(out) = &a.cut_out {
                    let cv: Vec<u32> = cur_mesh.faces().iter().flatten().copied().collect();
                    let fs: Vec<u32> = (0..cur_mesh.n_faces() as u32).map(|f| 3 * f).collect();
                    write_obj(out, &cur_pos, &cv, &fs, Some(&cf.corner_uvs)).map_err(|e| format!("{out}: {e}"))?;
                }
            }
            Err(e) => {
                println!("cut failed: {e}");
                summary.push(("cut".into(), format!("failed:{}", e.to_string().replace(' ', "_"))));
            }
        }
    }
    let line: Vec<String> = summary.iter().map(|(k, v)| format!("{k}={v}")).collect();
    println!("RESULT {}", line.join(" "));
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            println!("RESULT error={}", e.replace(' ', "_"));
            ExitCode::FAILURE
        }
    }
}
