//! Which stage of a step for several regions changes the energy by how much (SPEC §4, D26).
//!
//! ```text
//! cargo run --release -p ddg-variational-cuts --example rise_study -- MESH ENERGY ALPHA [KEY=VALUE ...]
//!
//! MESH      an OBJ file, or icosphere:LEVEL
//! ENERGY    dirichlet | hencky | scaled
//! ALPHA     the length weight, or auto:L for the ripple rule
//! clusters=8                  clusters of normals as the start
//! stages=paper | default      FlowOptions::for_regions | FlowOptions::new
//! length=mesh | cut | none    the length step (none: the step is taken with α_L = 0, the
//!                             energy is still measured with α_L)
//! ext=harmonic | local        speed=unit | mean        redist=keep | off
//! stats=robust | extremes
//! step=adaptive:TS:MOTION:TAUMAX | capped:TAUMAX:MOTION | fixed:TAU
//! steps=40                    number of steps
//! time=T                      stop when the flow time reaches T
//! quiet=1                     only the summary
//! ```
//!
//! For every step one line: the energy of the state and its two parts, the number of patches,
//! and the change of the energy by each stage, measured by evaluating the intermediate states
//! (`φ̃` after the explicit step, `φ′` after the length step on the mesh, the state after the
//! local limit at the vertices, after the projection of Losasso et al., and the new state after
//! the redistancing).

use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CutProblem, Energy, Evaluation, Evaluator, Extension, FlowOptions, LengthStep,
    PAPER_FEATURE_SMOOTHING, Redistance, RegionFlow, ScaleWeight, SpeedScaling, StepSize,
    StepStatistics, init, length_weight_for_ripples,
};

struct Model {
    mesh: TriMesh,
    geometry: Geometry,
    positions: Vec<[f64; 3]>,
}

fn load(path: &str) -> Model {
    if let Some(level) = path.strip_prefix("icosphere:") {
        let g = ddg_core::generators::icosphere(level.parse().unwrap()).unwrap();
        let mesh = TriMesh::new(g.positions.len(), g.faces).unwrap();
        let geometry = Geometry::from_positions(&mesh, &g.positions).unwrap();
        return Model {
            mesh,
            geometry,
            positions: g.positions,
        };
    }
    let obj = ddg_core::io::read_obj(path).unwrap();
    let mut faces = Vec::new();
    for f in 0..obj.n_faces() {
        let corners = &obj.corner_verts[obj.face_corners(f)];
        for k in 1..corners.len().saturating_sub(1) {
            faces.push([corners[0], corners[k], corners[k + 1]]);
        }
    }
    let mesh = TriMesh::new(obj.positions.len(), faces).unwrap();
    let geometry = Geometry::from_positions(&mesh, &obj.positions).unwrap();
    Model {
        mesh,
        geometry,
        positions: obj.positions,
    }
}

fn patches(e: &Evaluation) -> usize {
    e.cut_mesh().map_or(0, |m| m.mesh.n_components())
}

fn distortion(e: &Evaluation, problem: &CutProblem) -> f64 {
    e.energy
        .distortion
        .iter()
        .zip(&problem.regions)
        .map(|(j, w)| w.distortion * j)
        .sum()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (path, energy_name, alpha) = (&args[1], &args[2], &args[3]);
    let key = |name: &str, default: &str| -> String {
        args[4..]
            .iter()
            .find_map(|a| a.strip_prefix(&format!("{name}=")).map(str::to_string))
            .unwrap_or(default.to_string())
    };
    let model = load(path);
    let n: usize = key("clusters", "8").parse().unwrap();
    let energy = match energy_name.as_str() {
        "dirichlet" => Energy::Dirichlet,
        "hencky" => Energy::Hencky,
        _ => Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
            r_eps: 1.0,
            smoothing: PAPER_FEATURE_SMOOTHING,
        }),
    };
    let floats = |text: &str| -> Vec<f64> { text.split(':').map(|x| x.parse().unwrap()).collect() };
    let step = key("step", "adaptive:0.05:0.25:1");
    let step_size = match step.split_once(':').unwrap() {
        ("adaptive", rest) => {
            let v = floats(rest);
            StepSize::Adaptive {
                time_scale: v[0],
                max_motion: v[1],
                tau_max: v[2],
            }
        }
        ("capped", rest) => {
            let v = floats(rest);
            StepSize::Capped {
                tau_max: v[0],
                max_motion: v[1],
            }
        }
        ("fixed", rest) => StepSize::Fixed(rest.parse().unwrap()),
        other => panic!("unknown step rule {other:?}"),
    };
    let mut options = match key("stages", "paper").as_str() {
        "paper" => FlowOptions::for_regions(step_size),
        "default" => FlowOptions::new(step_size),
        other => panic!("unknown stage set {other}"),
    };
    let length = key("length", "");
    match length.as_str() {
        "mesh" => options.length_step = LengthStep::OnMesh,
        "cut" => options.length_step = LengthStep::OnCut,
        _ => {}
    }
    match key("ext", "").as_str() {
        "harmonic" => options.extension = Extension::Harmonic,
        "local" => options.extension = Extension::Local,
        _ => {}
    }
    match key("speed", "").as_str() {
        "unit" => options.speed_scaling = SpeedScaling::Unit,
        "mean" => options.speed_scaling = SpeedScaling::MeanGradient,
        _ => {}
    }
    match key("redist", "").as_str() {
        "off" => options.redistance = Redistance::Off,
        "keep" => options.redistance = Redistance::KeepAdjacent,
        _ => {}
    }
    match key("stats", "").as_str() {
        "extremes" => options.step_statistics = StepStatistics::Extremes,
        "robust" => options.step_statistics = StepStatistics::Robust,
        _ => {}
    }
    let max_steps: usize = key("steps", "40").parse().unwrap();
    let max_time: f64 = key("time", "inf").parse().unwrap();
    let quiet = key("quiet", "0") == "1";

    let labels = init::normal_clusters(&model.mesh, &model.positions, n, 100).unwrap();
    let phis = init::regions_from_vertex_labels(&model.mesh, &model.geometry, &labels, n).unwrap();
    let mut flow = RegionFlow::new(&model.mesh, &model.geometry, &phis).unwrap();
    let mut problem = CutProblem::for_regions(energy, 0.0, n);
    problem.length_weight = match alpha.strip_prefix("auto:") {
        Some(ripple) => length_weight_for_ripples(
            flow.roughening_rate(&problem, &options.eval).unwrap(),
            ripple.parse().unwrap(),
        ),
        None => alpha.parse().unwrap(),
    };
    // The problem that the step sees: without the length weight for `length=none`.
    let mut stepped = problem.clone();
    if length == "none" {
        stepped.length_weight = 0.0;
    }
    let mut evaluator = Evaluator::new(&model.mesh, &model.geometry).unwrap();
    let mut measure = |phis: &[Vec<f64>]| -> Evaluation {
        evaluator
            .evaluate_regions(phis, &problem, &options.eval)
            .unwrap()
    };

    println!(
        "{path}: {} vertices, {n} clusters, energy {energy_name}, alpha_L {:.4e}; {}",
        model.mesh.n_vertices(),
        problem.length_weight,
        args[4..].join(" ")
    );
    if !quiet {
        println!(
            "step time tau | E distortion length_term patches cut_points | dE: explicit length_step \
             limit projection redistance total | patches after: explicit length_step | max|s| \
             omega omega_max limited motion"
        );
    }
    let start = measure(flow.level_sets_unit_area());
    let mut state = measure(flow.level_sets_unit_area());
    let mut by_stage = [0.0f64; 5];
    let (mut rises, mut largest_rise) = (0usize, 0.0f64);
    let (mut lowest, mut lowest_step) = (start.energy.total, 0usize);
    let mut energies = vec![start.energy.total];
    let mut taken = 0;
    let clock = std::time::Instant::now();
    for k in 0..max_steps {
        if flow.time() >= max_time {
            break;
        }
        let (report, stages) = flow.step_with_stages(&stepped, &options).unwrap();
        taken += 1;
        if report.cut_vanished {
            println!("the cut vanished in step {k}");
            break;
        }
        let explicit = measure(&stages.explicit);
        let after_length = measure(&stages.length);
        let limited = measure(&stages.limited);
        let projected = if stages.projected.is_empty() {
            None
        } else {
            Some(measure(&stages.projected))
        };
        let new = measure(flow.level_sets_unit_area());
        let e0 = state.energy.total;
        let before_redistance = projected
            .as_ref()
            .map_or(limited.energy.total, |e| e.energy.total);
        let changes = [
            explicit.energy.total - e0,
            after_length.energy.total - explicit.energy.total,
            limited.energy.total - after_length.energy.total,
            before_redistance - limited.energy.total,
            new.energy.total - before_redistance,
        ];
        for (sum, change) in by_stage.iter_mut().zip(changes) {
            *sum += change;
        }
        if !quiet {
            println!(
                "{k:4} {:.4e} {:.3e} | {:.6e} {:.4e} {:.4e} {:3} {:5} | {:+.2e} {:+.2e} {:+.2e} {:+.2e} {:+.2e} {:+.2e} | {:3} {:3} | {:.3e} {:.3e} {:.3e} {:4} {:.2e}",
                report.time,
                report.tau,
                e0,
                distortion(&state, &problem),
                e0 - distortion(&state, &problem),
                patches(&state),
                state.cut.points.len(),
                changes[0],
                changes[1],
                changes[2],
                changes[3],
                changes[4],
                new.energy.total - e0,
                patches(&explicit),
                patches(&after_length),
                report.max_speed,
                report.roughening_rate,
                report.roughening_rate_max,
                report.n_limited,
                report.cut_motion,
            );
        }
        if new.energy.total > e0 {
            rises += 1;
            largest_rise = largest_rise.max((new.energy.total - e0) / e0.abs());
        }
        if new.energy.total < lowest {
            lowest = new.energy.total;
            lowest_step = k + 1;
        }
        energies.push(new.energy.total);
        state = new;
    }
    let seconds = clock.elapsed().as_secs_f64();
    println!(
        "summary: {taken} steps, flow time {:.4e}, {seconds:.1} s (with the measurements); energy {:.6e} -> {:.6e} \
         (lowest {:.6e} after {lowest_step} steps; rose in {rises} steps, at most {:.2e} of the energy); \
         distortion {:.4e} -> {:.4e}; cut length {:.4} -> {:.4}; patches {} -> {}; \
         sum of the changes by stage: explicit {:+.3e}, length step {:+.3e}, limit {:+.3e}, projection {:+.3e}, redistancing {:+.3e}",
        flow.time(),
        start.energy.total,
        state.energy.total,
        lowest,
        largest_rise,
        distortion(&start, &problem),
        distortion(&state, &problem),
        start.energy.length,
        state.energy.length,
        patches(&start),
        patches(&state),
        by_stage[0],
        by_stage[1],
        by_stage[2],
        by_stage[3],
        by_stage[4],
    );
    if !quiet {
        println!("region: distortion at the start -> at the end; area at the start -> at the end");
        for r in 0..n {
            println!(
                "  {r}: {:.4e} -> {:.4e}; {:.4} -> {:.4}",
                start.energy.distortion[r],
                state.energy.distortion[r],
                start.energy.area[r],
                state.energy.area[r]
            );
        }
    }
}
