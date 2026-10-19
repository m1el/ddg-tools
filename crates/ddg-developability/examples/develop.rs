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
//!   --min-face-angle X   admissibility floor on the smallest corner angle in radians (default 1e-4; 0 = off)
//!   --threads N          threads of the flow (default: all)
//!   --rounds SPEC        explicit rounds, comma separated `[+]energy:iterations` (`+` = 4-1 subdivide
//!                        first), e.g. `combinatorial:1000,+combinatorial:1000,+covariance:1000`;
//!                        overrides --energy/--subdivisions/--iterations (the energy reported is the last round's)
//!   --gradient-tol X     stop a round when mean edge length * |grad|_inf <= X (default 1e-8)
//!   --memory N           L-BFGS memory (default 10)
//!   --initial-step X     steepest-descent trial step in mean edge lengths (default 0.1)
//!   --stall-iters N      consecutive stalled iterations that end a round (default 10)
//!   --energy-tol X       relative decrease that counts as stalled (default 1e-10)
//!   --history FILE       write `round iteration energy gradient_norm` per accepted iterate
//!   --out FILE.obj       write the developable mesh; the per-vertex energy goes to FILE.energy.txt
//!   --cut-out FILE.obj   cut and flatten the result, write an OBJ with per-corner UVs (charts unpacked)
//!   --epsilon X          terminal threshold of the cut (default 1e-3)
//!   --cut-energy NAME    energy that selects the terminals (default: the CutOptions default)
//!   --defect-tol X       also make every vertex with |angle defect| > X a terminal (default: off)
//!   --cut-only           skip the flow, cut and flatten the input
//!   --conjugate          BFF conjugate extension instead of the harmonic one
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
    CutOptions, Energy, EnergyOptions, FlowOptions, Method, OptimizeStatus, RemeshOptions, Round, cut_and_flatten,
    evaluate_energy, flow,
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
    min_face_angle: Option<f64>,
    threads: Option<usize>,
    out: Option<String>,
    cut_out: Option<String>,
    epsilon: f64,
    cut_energy: Option<Energy>,
    defect_tol: Option<f64>,
    cut_only: bool,
    no_cut: bool,
    conjugate: bool,
    rounds: Option<Vec<(bool, Energy, usize)>>,
    gradient_tol: Option<f64>,
    memory: Option<usize>,
    initial_step: Option<f64>,
    stall_iters: Option<usize>,
    energy_tol: Option<f64>,
    history: Option<String>,
}

fn parse_rounds(s: &str) -> Result<Vec<(bool, Energy, usize)>, String> {
    s.split(',')
        .map(|r| {
            let (sub, r) = match r.strip_prefix('+') {
                Some(rest) => (true, rest),
                None => (false, r),
            };
            let (e, it) = r.split_once(':').ok_or_else(|| format!("round {r}: expected energy:iterations"))?;
            let e = parse_energy(e).ok_or_else(|| format!("unknown energy {e}"))?;
            Ok((sub, e, it.parse().map_err(|x| format!("round {r}: {x}"))?))
        })
        .collect()
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
        min_face_angle: None,
        threads: None,
        out: None,
        cut_out: None,
        epsilon: 1e-3,
        cut_energy: None,
        defect_tol: None,
        cut_only: false,
        no_cut: false,
        conjugate: false,
        rounds: None,
        gradient_tol: None,
        memory: None,
        initial_step: None,
        stall_iters: None,
        energy_tol: None,
        history: None,
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
            "--min-face-angle" => a.min_face_angle = Some(val("--min-face-angle")?.parse().map_err(|e| format!("{e}"))?),
            "--threads" => a.threads = Some(val("--threads")?.parse().map_err(|e| format!("{e}"))?),
            "--out" => a.out = Some(val("--out")?),
            "--cut-out" => a.cut_out = Some(val("--cut-out")?),
            "--epsilon" => a.epsilon = val("--epsilon")?.parse().map_err(|e| format!("{e}"))?,
            "--cut-energy" => {
                let v = val("--cut-energy")?;
                a.cut_energy = Some(parse_energy(&v).ok_or_else(|| format!("unknown energy {v}"))?);
            }
            "--defect-tol" => a.defect_tol = Some(val("--defect-tol")?.parse().map_err(|e| format!("{e}"))?),
            "--cut-only" => a.cut_only = true,
            "--no-cut" => a.no_cut = true,
            "--conjugate" => a.conjugate = true,
            "--rounds" => a.rounds = Some(parse_rounds(&val("--rounds")?)?),
            "--gradient-tol" => a.gradient_tol = Some(val("--gradient-tol")?.parse().map_err(|e| format!("{e}"))?),
            "--memory" => a.memory = Some(val("--memory")?.parse().map_err(|e| format!("{e}"))?),
            "--initial-step" => a.initial_step = Some(val("--initial-step")?.parse().map_err(|e| format!("{e}"))?),
            "--stall-iters" => a.stall_iters = Some(val("--stall-iters")?.parse().map_err(|e| format!("{e}"))?),
            "--history" => a.history = Some(val("--history")?),
            "--energy-tol" => a.energy_tol = Some(val("--energy-tol")?.parse().map_err(|e| format!("{e}"))?),
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
        OptimizeStatus::Stationary => "stationary",
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
    let eopts = EnergyOptions::new(a.rounds.as_ref().and_then(|r| r.last()).map_or(a.energy, |r| r.1));

    if !a.cut_only {
        let mut fo = FlowOptions::uniform(eopts.clone(), a.subdivisions, a.iterations);
        if let Some(r) = &a.rounds {
            fo.rounds = r
                .iter()
                .map(|&(subdivide_before, e, max_iterations)| Round {
                    energy: EnergyOptions::new(e),
                    max_iterations,
                    subdivide_before,
                })
                .collect();
        }
        if let Some(x) = a.gradient_tol {
            fo.optimize.gradient_tolerance = x;
        }
        if let Some(m) = a.memory {
            fo.optimize.method = Method::Lbfgs { memory: m };
        }
        if let Some(x) = a.initial_step {
            fo.optimize.initial_step = x;
        }
        if let Some(x) = a.stall_iters {
            fo.optimize.stall_iterations = x;
        }
        if let Some(x) = a.energy_tol {
            fo.optimize.energy_tolerance = x;
        }
        fo.optimize.max_displacement = a.max_displacement;
        if let Some(x) = a.min_face_angle {
            fo.optimize.min_face_angle = x;
        }
        fo.remesh = a.remesh.then(|| RemeshOptions {
            small_angle: a.small_angle_deg.to_radians(),
            ..Default::default()
        });
        if let (Some(x), Some(r)) = (a.min_face_angle, fo.remesh.as_mut()) {
            r.min_face_angle = x;
        }
        fo.threads = a.threads;
        let t = Instant::now();
        let res = flow(&mesh, &positions, &fo, None).map_err(|e| format!("flow failed: {e}"))?;
        let secs = t.elapsed().as_secs_f64();
        let mut flips = 0;
        let mut collapses = 0;
        let mut iters = 0;
        let mut evals = 0;
        let mut undefined = 0;
        let mut inadmissible = 0;
        let mut resets = 0;
        for (i, r) in res.report.rounds.iter().enumerate() {
            let o = &r.optimize;
            println!(
                "round {i}: {} v / {} f, E {:.6e} -> {:.6e}, |g| {:.3e} -> {:.3e}, {} iterations in {} segments ({}), {} evals, remesh {:?}, {:.2}s",
                r.n_vertices,
                r.n_faces,
                o.energy_initial,
                o.energy_final,
                o.gradient_norm_history.first().copied().unwrap_or(f64::NAN),
                o.gradient_norm_history.last().copied().unwrap_or(f64::NAN),
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
            inadmissible += o.inadmissible_trials;
            resets += o.memory_resets;
        }
        if let Some(h) = &a.history {
            let mut txt = String::new();
            for (i, r) in res.report.rounds.iter().enumerate() {
                for (k, (e, g)) in r.optimize.energy_history.iter().zip(&r.optimize.gradient_norm_history).enumerate() {
                    txt.push_str(&format!("{i} {k} {e:e} {g:e}\n"));
                }
            }
            std::fs::write(h, txt).map_err(|e| format!("{h}: {e}"))?;
        }
        let first = &res.report.rounds[0];
        let Some(last) = res.report.rounds.last() else { return Err(format!("flow ran no round: {:?}", res.report.stopped)) };
        let (e_final, ve) = evaluate_energy(&res.mesh, &res.positions, &eopts).map_err(|e| format!("final energy: {e}"))?;
        summary.extend([
            ("e_initial".into(), format!("{:.6e}", first.optimize.energy_initial)),
            ("e_final".into(), format!("{:.6e}", e_final)),
            ("seconds".into(), format!("{secs:.3}")),
            ("status".into(), res.report.rounds.iter().map(|r| status_name(r.optimize.status)).collect::<Vec<_>>().join("/")),
            ("gnorm_final".into(), res.report.rounds.iter().map(|r| format!("{:.2e}", r.optimize.gradient_norm_history.last().copied().unwrap_or(f64::NAN))).collect::<Vec<_>>().join("/")),
            ("round_iterations".into(), res.report.rounds.iter().map(|r| r.optimize.iterations.to_string()).collect::<Vec<_>>().join("/")),
            ("iterations".into(), iters.to_string()),
            ("evaluations".into(), evals.to_string()),
            ("undefined_trials".into(), undefined.to_string()),
            ("inadmissible_trials".into(), inadmissible.to_string()),
            ("memory_resets".into(), resets.to_string()),
            ("flips".into(), flips.to_string()),
            ("collapses".into(), collapses.to_string()),
            ("faces_out".into(), last.n_faces.to_string()),
            ("min_angle".into(), format!("{:.3e}", res.report.min_angle)),
            ("bad_faces".into(), res.report.bad_faces.to_string()),
            ("stopped".into(), res.report.stopped.as_ref().map_or("none".into(), |s| format!("round{}:v{}", s.round, s.vertex))),
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
        let co = CutOptions {
            energy: a.cut_energy.unwrap_or(CutOptions::default().energy),
            defect_tolerance: a.defect_tol,
            epsilon: a.epsilon,
            threads: a.threads,
            extension: if a.conjugate { ddg_bff::Extension::Conjugate } else { ddg_bff::Extension::Harmonic },
            ..Default::default()
        };
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
                // diagnostics (SPEC §12.3): cut length, defects, embeddedness, flips
                let lens: Vec<f64> = cur_mesh
                    .edges()
                    .iter()
                    .map(|&[p, q]| {
                        let (p, q) = (cur_pos[p as usize], cur_pos[q as usize]);
                        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
                    })
                    .collect();
                let total_len: f64 = lens.iter().sum();
                let cut_len: f64 = cf.cut_edges.iter().map(|&e| lens[e as usize]).sum();
                let mut is_term = vec![false; cur_mesh.n_vertices()];
                for &t in cf.terminals.iter().chain(&cf.added_terminals) {
                    is_term[t as usize] = true;
                }
                let nv = cur_mesh.n_vertices();
                let (mut om_term, mut om_non, mut non_emb, mut non_emb_non, mut big_non) = (0.0, 0.0, 0, 0, 0);
                let (mut om_non_bad, mut big_non_bad) = (0.0, 0);
                let mut max_om_non: f64 = 0.0;
                #[allow(clippy::needless_range_loop)]
                for v in 0..nv {
                    if cur_mesh.is_boundary_vertex(v) {
                        continue;
                    }
                    let (om, emb) = (cf.vertex_defect[v].abs(), cf.star_embeddedness[v]);
                    if is_term[v] {
                        om_term += om;
                    } else {
                        om_non += om;
                        max_om_non = max_om_non.max(om);
                        if om > 0.1 {
                            big_non += 1;
                        }
                        if emb.is_nan() || emb <= 0.0 {
                            om_non_bad += om;
                            if om > 0.1 {
                                big_non_bad += 1;
                            }
                        }
                    }
                    if emb.is_nan() || emb <= 0.0 {
                        non_emb += 1;
                        if !is_term[v] {
                            non_emb_non += 1;
                        }
                    }
                }
                // flipped uv triangles (negative signed area) and how many touch a non-embedded star
                let (mut flips, mut flips_near) = (0, 0);
                for (f, tri) in cur_mesh.faces().iter().enumerate() {
                    let (a, b, c) = (cf.corner_uvs[3 * f], cf.corner_uvs[3 * f + 1], cf.corner_uvs[3 * f + 2]);
                    if (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) < 0.0 {
                        flips += 1;
                        if tri.iter().any(|&v| { let e = cf.star_embeddedness[v as usize]; e.is_nan() || e <= 0.0 }) {
                            flips_near += 1;
                        }
                    }
                }
                println!(
                    "diag: cut length {:.4} of total edge length, sum|Omega| terminals {om_term:.3} / non-terminals {om_non:.3} (max {max_om_non:.3}, {big_non} above 0.1), non-embedded stars {non_emb} of {} interior ({non_emb_non} non-terminal, their sum|Omega| {om_non_bad:.3}, {big_non_bad} above 0.1), flipped {flips} ({flips_near} touch one)",
                    cut_len / total_len,
                    (0..nv).filter(|&v| !cur_mesh.is_boundary_vertex(v)).count()
                );
                summary.extend([
                    ("cut_len_frac".into(), format!("{:.5}", cut_len / total_len)),
                    ("om_terminals".into(), format!("{om_term:.4}")),
                    ("om_nonterminals".into(), format!("{om_non:.4}")),
                    ("max_om_nonterminal".into(), format!("{max_om_non:.4}")),
                    ("nonembedded".into(), non_emb.to_string()),
                    ("nonembedded_nonterminal".into(), non_emb_non.to_string()),
                    ("om_nonterminal_nonembedded".into(), format!("{om_non_bad:.4}")),
                    ("big_nonterminal".into(), big_non.to_string()),
                    ("big_nonterminal_nonembedded".into(), big_non_bad.to_string()),
                    ("flips_geom".into(), flips.to_string()),
                    ("flips_near_nonembedded".into(), flips_near.to_string()),
                    ("faces_cut".into(), cur_mesh.n_faces().to_string()),
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
