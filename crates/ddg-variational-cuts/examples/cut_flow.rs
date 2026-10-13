//! Runs the cut flow on an OBJ file and writes the history and the result.
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example cut_flow -- MESH.obj [options]
//!
//! --out DIR                   output directory (default runs_tmp/cut_flow)
//! --energy dirichlet | hencky | scaled[:R_EPS[:SMOOTHING]]     (R_EPS defaults to 1)
//! --alpha A | auto[:L]        length weight; auto sets it from the roughening rate of the
//!                             initial cut so that ripples longer than L grow (default auto:0.25)
//! --weights W0,W1             distortion weights of the two regions (default 1,1)
//! --area B0,B1                area weights of the two regions (default 0,0)
//! --curvature gaussian | uniform:K
//! --camera X,Y,Z --penalty W  visibility penalty W * max(n.e, 0), every vertex visible
//! --penalty-rule gradient | paper
//! --init halfspace[:NX,NY,NZ[,OFFSET]] | circle:VERTEX:RADIUS | phi:FILE | labels:FILE
//!        | face-labels:FILE   (default: halfspace across the longest axis of the bounding box,
//!                             through the centroid; files hold one number per line)
//! --raw-init                  do not redistance a halfspace or phi start
//! --symmetry FILE             one class label per vertex
//! --step adaptive:TIME_SCALE:MOTION | capped:TAU_MAX:MOTION | fixed:TAU
//!                             (default adaptive:0.05:0.25)
//! --steps N                   at most N steps (default 1000)
//! --tolerance T               stop when the cut moves less than T per step (default 1e-5)
//! --energy-tolerance R        stop when the energy fell by less than R of its size over 25
//!                             steps (default 1e-4)
//! --paper                     the stages of the paper as printed
//! --length-step cut | mesh    --extension local | harmonic    --speed unit | mean
//! --flux nodal | smoothed:C   --delta D (smallest crossing offset)
//! --threads N                 (default 1)
//! ```
//!
//! Lengths on the command line are in the unit of the file, except `TAU`, `A`, `B`, `W`, `K`,
//! `R_EPS`, `SMOOTHING` and `T`, which refer to the surface scaled to unit area like every
//! parameter of the library.
//!
//! Output: `steps.csv` (one line per step), `phi.txt` (the final level-set function in the unit
//! of the file), `cut.obj` (the cut as polylines), `cut_mesh.obj` (the cut mesh with one group
//! per region).

use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;

use ddg_core::io::read_obj;
use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::length_weight_for_ripples;
use ddg_variational_cuts::{
    CurvatureSource, CutProblem, Energy, EvalOptions, Extension, Flow, FlowOptions, FluxRecovery,
    LengthStep, PAPER_FEATURE_SMOOTHING, Parallelism, PenaltyRule, ScaleWeight, SpeedScaling,
    StepSize, StopRule, ThreadPool, init,
};

type AnyError = Box<dyn std::error::Error>;

fn numbers(text: &str) -> Result<Vec<f64>, AnyError> {
    text.split(',')
        .map(|x| {
            x.trim()
                .parse::<f64>()
                .map_err(|e| format!("{x:?}: {e}").into())
        })
        .collect()
}

fn column(path: &str) -> Result<Vec<f64>, AnyError> {
    std::fs::read_to_string(path)
        .map_err(|e| format!("{path}: {e}"))?
        .split_whitespace()
        .map(|x| {
            x.parse::<f64>()
                .map_err(|e| format!("{path}: {x:?}: {e}").into())
        })
        .collect()
}

fn labels(path: &str) -> Result<Vec<u32>, AnyError> {
    Ok(column(path)?.iter().map(|x| *x as u32).collect())
}

fn run() -> Result<(), AnyError> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut input = None;
    let mut out = PathBuf::from("runs_tmp/cut_flow");
    let mut problem = CutProblem::new(Energy::Dirichlet, 1.0);
    let mut step_size = StepSize::default();
    let mut stop = StopRule::default();
    let mut paper = false;
    let (mut length_step, mut extension, mut speed) = (None, None, None);
    let (mut flux, mut delta) = (None, None);
    let mut penalty_rule = PenaltyRule::Gradient;
    let mut init_spec = String::from("halfspace");
    let mut raw_init = false;
    let (mut camera, mut penalty_weight) = (None, 1.0);
    let mut symmetry = None;
    let mut threads = 1usize;
    let mut ripple = Some(0.25);

    let mut k = 0;
    while k < args.len() {
        let flag = args[k].as_str();
        let mut value = || -> Result<&str, AnyError> {
            k += 1;
            args.get(k)
                .map(String::as_str)
                .ok_or_else(|| format!("{flag} needs a value").into())
        };
        match flag {
            "--out" => out = PathBuf::from(value()?),
            "--energy" => {
                let v = value()?;
                let parts: Vec<&str> = v.split(':').collect();
                problem.energy = match parts[0] {
                    "dirichlet" => Energy::Dirichlet,
                    "hencky" => Energy::Hencky,
                    "scaled" => Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
                        r_eps: parts.get(1).map_or(Ok(1.0), |x| x.parse())?,
                        smoothing: parts
                            .get(2)
                            .map_or(Ok(PAPER_FEATURE_SMOOTHING), |x| x.parse())?,
                    }),
                    other => return Err(format!("unknown energy {other:?}").into()),
                };
            }
            "--alpha" => {
                let v = value()?;
                match v.split_once(':') {
                    Some(("auto", wavelength)) => ripple = Some(wavelength.parse()?),
                    None if v == "auto" => ripple = Some(0.25),
                    _ => {
                        ripple = None;
                        problem.length_weight = v.parse()?;
                    }
                }
            }
            "--weights" => {
                let w = numbers(value()?)?;
                for (region, w) in problem.regions.iter_mut().zip(w) {
                    region.distortion = w;
                }
            }
            "--area" => {
                let b = numbers(value()?)?;
                for (region, b) in problem.regions.iter_mut().zip(b) {
                    region.area = b;
                }
            }
            "--curvature" => {
                let v = value()?;
                problem.curvature = match v.split_once(':') {
                    None if v == "gaussian" => CurvatureSource::Gaussian,
                    Some(("uniform", density)) => CurvatureSource::Uniform(density.parse()?),
                    _ => return Err(format!("unknown curvature {v:?}").into()),
                };
            }
            "--camera" => camera = Some(numbers(value()?)?),
            "--penalty" => penalty_weight = value()?.parse()?,
            "--penalty-rule" => {
                penalty_rule = match value()? {
                    "gradient" => PenaltyRule::Gradient,
                    "paper" => PenaltyRule::PaperDiagonal,
                    other => return Err(format!("unknown penalty rule {other:?}").into()),
                }
            }
            "--init" => init_spec = value()?.to_string(),
            "--raw-init" => raw_init = true,
            "--symmetry" => symmetry = Some(labels(value()?)?),
            "--step" => {
                let v = value()?;
                let parts: Vec<&str> = v.split(':').collect();
                step_size = match parts.as_slice() {
                    ["fixed", tau] => StepSize::Fixed(tau.parse()?),
                    ["capped", tau_max, motion] => StepSize::Capped {
                        tau_max: tau_max.parse()?,
                        max_motion: motion.parse()?,
                    },
                    ["adaptive", time_scale, motion] => StepSize::Adaptive {
                        time_scale: time_scale.parse()?,
                        max_motion: motion.parse()?,
                        tau_max: 1.0,
                    },
                    _ => return Err(format!("unknown step rule {v:?}").into()),
                };
            }
            "--steps" => stop.max_steps = value()?.parse()?,
            "--tolerance" => stop.motion_tolerance = value()?.parse()?,
            "--energy-tolerance" => stop.energy_tolerance = value()?.parse()?,
            "--paper" => paper = true,
            "--length-step" => {
                length_step = Some(match value()? {
                    "cut" => LengthStep::OnCut,
                    "mesh" => LengthStep::OnMesh,
                    other => return Err(format!("unknown length step {other:?}").into()),
                })
            }
            "--extension" => {
                extension = Some(match value()? {
                    "local" => Extension::Local,
                    "harmonic" => Extension::Harmonic,
                    other => return Err(format!("unknown extension {other:?}").into()),
                })
            }
            "--speed" => {
                speed = Some(match value()? {
                    "unit" => SpeedScaling::Unit,
                    "mean" => SpeedScaling::MeanGradient,
                    other => return Err(format!("unknown speed scaling {other:?}").into()),
                })
            }
            "--flux" => {
                let v = value()?;
                flux = Some(match v.split_once(':') {
                    None if v == "nodal" => FluxRecovery::Nodal,
                    Some(("smoothed", c)) => FluxRecovery::Smoothed { length: c.parse()? },
                    _ => return Err(format!("unknown flux recovery {v:?}").into()),
                });
            }
            "--delta" => delta = Some(value()?.parse::<f64>()?),
            "--threads" => threads = value()?.parse()?,
            other if other.starts_with("--") => {
                return Err(format!("unknown option {other}").into());
            }
            path => input = Some(path.to_string()),
        }
        k += 1;
    }
    let input = input.ok_or("usage: cut_flow MESH.obj [options]; see the top of the example")?;

    // The mesh: polygons are split into fans.
    let obj = read_obj(&input)?;
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let positions = obj.positions;
    let mesh = TriMesh::new(positions.len(), faces)?;
    let geometry = Geometry::from_positions(&mesh, &positions)?;
    println!(
        "{input}: {} vertices, {} faces, {} boundary loops, area {:.6}",
        mesh.n_vertices(),
        mesh.n_faces(),
        mesh.boundary_loops().len(),
        geometry.total_area()
    );

    // The initial function.
    let parts: Vec<&str> = init_spec.split(':').collect();
    let phi0 = match parts.as_slice() {
        ["halfspace"] | ["halfspace", _] => {
            let n = positions.len() as f64;
            let centroid = [0, 1, 2].map(|k| positions.iter().map(|p| p[k]).sum::<f64>() / n);
            let given = match parts.get(1) {
                Some(text) => numbers(text)?,
                None => {
                    let extent = [0, 1, 2].map(|k| {
                        let (lo, hi) =
                            positions.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
                                (lo.min(p[k]), hi.max(p[k]))
                            });
                        hi - lo
                    });
                    let axis = (0..3)
                        .max_by(|a, b| extent[*a].total_cmp(&extent[*b]))
                        .unwrap();
                    let mut normal = vec![0.0; 3];
                    normal[axis] = 1.0;
                    normal
                }
            };
            if given.len() < 3 {
                return Err("halfspace needs NX,NY,NZ[,OFFSET]".into());
            }
            let offset = given
                .get(3)
                .copied()
                .unwrap_or_else(|| (0..3).map(|k| given[k] * centroid[k]).sum());
            let phi: Vec<f64> = positions
                .iter()
                .map(|p| (0..3).map(|k| given[k] * p[k]).sum::<f64>() - offset)
                .collect();
            if raw_init {
                phi
            } else {
                init::redistanced(&mesh, &geometry, &phi)?.0
            }
        }
        ["circle", vertex, radius] => {
            init::geodesic_circle(&mesh, &geometry, vertex.parse()?, radius.parse()?)?
        }
        ["phi", path] => {
            let phi = column(path)?;
            if raw_init {
                phi
            } else {
                init::redistanced(&mesh, &geometry, &phi)?.0
            }
        }
        ["labels", path] => init::from_vertex_labels(&mesh, &geometry, &labels(path)?)?,
        ["face-labels", path] => init::from_face_labels(&mesh, &geometry, &labels(path)?)?,
        _ => return Err(format!("unknown initial cut {init_spec:?}").into()),
    };

    if let Some(camera) = camera {
        if camera.len() != 3 {
            return Err("--camera needs X,Y,Z".into());
        }
        let visible = vec![penalty_weight; positions.len()];
        problem.penalty = Some(init::visibility_penalty(
            &mesh,
            &positions,
            [camera[0], camera[1], camera[2]],
            &visible,
        )?);
    }

    let pool = if threads > 1 {
        Some(ThreadPool::new(threads)?)
    } else {
        None
    };
    let mut options = if paper {
        FlowOptions::paper(step_size)
    } else {
        FlowOptions::new(step_size)
    };
    options.eval = EvalOptions {
        parallelism: pool
            .as_ref()
            .map_or(Parallelism::Sequential, Parallelism::Pool),
        ..options.eval
    };
    options.penalty_rule = penalty_rule;
    options.symmetry = symmetry;
    if let Some(x) = length_step {
        options.length_step = x;
    }
    if let Some(x) = extension {
        options.extension = x;
    }
    if let Some(x) = speed {
        options.speed_scaling = x;
    }
    if let Some(x) = flux {
        options.eval.flux_recovery = x;
    }
    if let Some(x) = delta {
        options.eval.min_crossing_offset = x;
    }

    let mut flow = Flow::new(&mesh, &geometry, &phi0)?;
    if let Some(wavelength) = ripple {
        let omega = flow.roughening_rate(&problem, &options.eval)?;
        problem.length_weight = length_weight_for_ripples(omega, wavelength);
        println!(
            "roughening rate of the initial cut {omega:.4e}: alpha_L = {:.4e} (ripples longer than {wavelength} grow)",
            problem.length_weight
        );
    }
    let start = std::time::Instant::now();
    let report = flow.run(&problem, &options, &stop)?;
    let seconds = start.elapsed().as_secs_f64();

    std::fs::create_dir_all(&out)?;
    let mut csv = String::from(
        "step,time,tau,energy,distortion_0,distortion_1,area_0,area_1,length,penalty,max_speed,\
         max_length_speed,roughening_rate,cut_motion,cut_points,clamped,unreached,\
         uncut_components,grad_min,grad_mean,grad_max,t_evaluate,t_extension,t_length,\
         t_redistance\n",
    );
    for (k, r) in report.steps.iter().enumerate() {
        let e = &r.energy;
        let g = r
            .gradient_norm
            .unwrap_or(ddg_variational_cuts::GradientNorm {
                min: f64::NAN,
                max: f64::NAN,
                mean: f64::NAN,
            });
        writeln!(
            csv,
            "{k},{:e},{:e},{:.12e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{:e},{},{},{},{},{:.4},{:.4},{:.4},{:.3e},{:.3e},{:.3e},{:.3e}",
            r.time,
            r.tau,
            e.total,
            e.distortion[0],
            e.distortion[1],
            e.area[0],
            e.area[1],
            e.length,
            e.penalty,
            r.max_speed,
            r.max_length_speed,
            r.roughening_rate,
            r.cut_motion,
            r.n_cut_points,
            r.n_clamped_crossings,
            r.n_unreached,
            r.n_uncut_closed_components,
            g.min,
            g.mean,
            g.max,
            r.timings.evaluate,
            r.timings.extension,
            r.timings.length,
            r.timings.redistance
        )?;
    }
    std::fs::write(out.join("steps.csv"), csv)?;

    let phi = flow.level_set();
    let mut text = String::new();
    for x in &phi {
        writeln!(text, "{x:e}")?;
    }
    std::fs::write(out.join("phi.txt"), text)?;

    let evaluation = flow.evaluate(&problem, &options.eval)?;
    let mut text = String::new();
    let mut first = 1usize;
    for (k, line) in evaluation
        .cut
        .polylines(&mesh, &positions)?
        .iter()
        .enumerate()
    {
        writeln!(
            text,
            "o cut_{k}_{}",
            if line.closed { "closed" } else { "open" }
        )?;
        for p in &line.positions {
            writeln!(text, "v {} {} {}", p[0], p[1], p[2])?;
        }
        let n = line.positions.len();
        let mut indices: Vec<String> = (first..first + n).map(|i| i.to_string()).collect();
        if line.closed {
            indices.push(first.to_string());
        }
        writeln!(text, "l {}", indices.join(" "))?;
        first += n;
    }
    std::fs::write(out.join("cut.obj"), text)?;

    let cut_mesh = evaluation.cut_mesh()?;
    let at = cut_mesh.positions(&evaluation.cut, &mesh, &positions)?;
    let mut text = String::new();
    for p in &at {
        writeln!(text, "v {} {} {}", p[0], p[1], p[2])?;
    }
    for region in 0..problem.regions.len() as u32 {
        writeln!(text, "g region_{region}")?;
        for (face, _) in cut_mesh
            .mesh
            .faces()
            .iter()
            .zip(&cut_mesh.face_region)
            .filter(|(_, r)| **r == region)
        {
            writeln!(text, "f {} {} {}", face[0] + 1, face[1] + 1, face[2] + 1)?;
        }
    }
    std::fs::write(out.join("cut_mesh.obj"), text)?;

    let first_energy = report.steps.first().map_or(f64::NAN, |r| r.energy.total);
    let e = &report.final_energy;
    println!(
        "{} steps in {seconds:.2} s, stopped: {:?}; flow time {:.4e}",
        report.steps.len(),
        report.stop,
        flow.time()
    );
    println!(
        "energy {first_energy:.6} -> {:.6}   distortion {:.6} + {:.6}   areas {:.4} / {:.4}   cut length {:.5} (unit area), penalty {:.5}",
        e.total, e.distortion[0], e.distortion[1], e.area[0], e.area[1], e.length, e.penalty
    );
    println!("wrote {}", out.display());
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("cut_flow: {e}");
            ExitCode::FAILURE
        }
    }
}
