//! Differential test against the Lean oracle (`lean/oracle`, an independent dense implementation
//! written from `SPEC.md` alone): the committed fixtures `lean/oracle/fixtures/<set>/*.in|out`
//! are replayed, read by the oracle's README alone.
//!
//! | set | what |
//! |---|---|
//! | `cuts` | two regions, the paper's stages (format 1) |
//! | `cuts-default` | two regions, every stage choice, the three step rules (format 2) |
//! | `cuts-robust` | two regions, the statistics of the step rule of SPEC D22 (format 3) |
//! | `regions` | `N` regions, both stage sets |
//! | `regions-robust` | `N` regions, the statistics of the step rule (format 2 of `regions`) |
//!
//! **Evaluation** (SPEC §3.1–§3.7, §6, §8.2–§8.3): the scaling, the face curvatures, the
//! crossings and triple points, every triangle of the cut mesh, the segments, `u`, `p`, `v`,
//! `ℓ`, `q`, the recovered `d`, `e` at every node, `G`, the speed `s` before and after sliver
//! smoothing, `J` and the areas per region, and the energy. Every case states its flux
//! recovery (`flux nodal` or `flux smoothed C`).
//!
//! **One step** (SPEC §3.8–§3.12, §4, §8.4–§8.5), with the stages and the statistics of the
//! step rule that the fixture names: `τ` and its three limits, the roughening rates, the length
//! force and the rates on the cut, `ḡ`, the extension, `φ̃`, `φ′`, the symmetry projection, the
//! projection of Losasso et al., the cut-adjacent vertices, the new level-set functions, the
//! number of unreached vertices, and the energy and the cut of the new state. With the robust
//! statistics also: the quantile `ω`, both speeds `m̂`, `h_c`, `ℓ_c`, and at every cut point or
//! boundary node the factor of the local limit and the rate before and after it.

use std::collections::HashMap;
use std::path::PathBuf;

use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CurvatureSource, CutProblem, Energy, EvalOptions, Evaluation, Evaluator, Extension, Flow,
    FlowOptions, FluxRecovery, LengthStep, Parallelism, PenaltyRule, ROUGHENING_QUANTILE,
    Redistance, RegionWeights, ScaleWeight, SpeedScaling, StepSize, StepStatistics,
    face_curvatures, weighted_quantile,
};

type Lines = HashMap<String, Vec<Vec<String>>>;

fn float(token: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(token, 16).unwrap_or_else(|_| panic!("bad float {token}")))
}

fn int(token: &str) -> usize {
    token
        .parse()
        .unwrap_or_else(|_| panic!("bad integer {token}"))
}

/// The lines of a block, grouped by their first word.
fn lines_by_key(text: &str) -> Lines {
    let mut out: Lines = HashMap::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        let mut words = line.split_whitespace();
        if let Some(key) = words.next() {
            out.entry(key.to_string())
                .or_default()
                .push(words.map(str::to_string).collect());
        }
    }
    out
}

struct Case {
    mesh: TriMesh,
    geometry: Geometry,
    phi: Vec<f64>,
    problem: CutProblem,
    /// `(r_eps, c)` if the energy is the curvature-scaled Dirichlet energy.
    feature: Option<(f64, f64)>,
    eval: EvalOptions<'static>,
    step_size: StepSize,
    penalty_rule: PenaltyRule,
    classes: Option<Vec<u32>>,
    /// Format 2: the stages; `None` in format 1 (the stages of the paper).
    stages: Option<(LengthStep, Extension, SpeedScaling)>,
    /// Format 3: the statistics of the step rule; the extremes where the line is absent.
    statistics: StepStatistics,
}

/// The line `statistics extremes | robust` of an input (absent: the extremes).
fn statistics_of(lines: &Lines) -> StepStatistics {
    match lines.get("statistics").map(|l| l[0][0].as_str()) {
        None | Some("extremes") => StepStatistics::Extremes,
        Some("robust") => StepStatistics::Robust,
        Some(other) => panic!("unknown statistics {other}"),
    }
}

/// Compares a number that the oracle may print as `+inf` (an absent limit).
fn close_or_infinite(report: &mut Report, what: &str, ours: f64, oracle: f64) {
    if oracle.is_infinite() || ours.is_infinite() {
        report.same(|| what.to_string(), ours, oracle);
    } else {
        report.close(|| what.to_string(), ours, oracle);
    }
}

/// The three limits of the step rule of SPEC §4 as the oracle prints them: `tau_max`, the
/// limit by the roughening rate and the limit by the motion (`+inf` where a rule has none or
/// the denominator is not positive).
fn step_limits(step: &StepSize, omega: f64, motion_speed: f64, edge: f64) -> [f64; 3] {
    let over = |numerator: f64, denominator: f64| {
        if denominator > 0.0 {
            numerator / denominator
        } else {
            f64::INFINITY
        }
    };
    match *step {
        StepSize::Fixed(tau) => [tau, f64::INFINITY, f64::INFINITY],
        StepSize::Capped {
            tau_max,
            max_motion,
        } => [
            tau_max,
            f64::INFINITY,
            over(max_motion * edge, motion_speed),
        ],
        StepSize::Adaptive {
            time_scale,
            max_motion,
            tau_max,
        } => [
            tau_max,
            over(time_scale, omega),
            over(max_motion * edge, motion_speed),
        ],
    }
}

fn parse_case(text: &str) -> Case {
    let lines = lines_by_key(text);
    let one = |key: &str| -> &Vec<String> {
        let all = lines
            .get(key)
            .unwrap_or_else(|| panic!("missing line {key}"));
        assert_eq!(all.len(), 1, "line {key} occurs {} times", all.len());
        &all[0]
    };
    let floats = |words: &[String]| -> Vec<f64> { words.iter().map(|w| float(w)).collect() };
    let nv = int(&one("mesh")[0]);
    let faces: Vec<[u32; 3]> = lines["f"]
        .iter()
        .map(|w| [0, 1, 2].map(|k| int(&w[k]) as u32))
        .collect();
    let mesh = TriMesh::new(nv, faces).unwrap();
    let mut lengths = vec![f64::NAN; mesh.n_edges()];
    for w in &lines["e"] {
        let e = mesh
            .find_edge(int(&w[0]), int(&w[1]))
            .expect("an edge of the mesh");
        lengths[e] = float(&w[2]);
    }
    let geometry = Geometry::from_edge_lengths(&mesh, lengths).unwrap();
    let energy_line = one("energy");
    let mut feature = None;
    let energy = match energy_line[0].as_str() {
        "dirichlet" => Energy::Dirichlet,
        "hencky" => Energy::Hencky,
        "scaled_a" => Energy::ScaledDirichlet(ScaleWeight::PerVertex(floats(&energy_line[1..]))),
        "scaled_curv" => {
            let (r_eps, smoothing) = (float(&energy_line[1]), float(&energy_line[2]));
            feature = Some((r_eps, smoothing));
            Energy::ScaledDirichlet(ScaleWeight::FeatureSize { r_eps, smoothing })
        }
        other => panic!("unknown energy {other}"),
    };
    let curvature_line = one("curvature");
    let curvature = match curvature_line[0].as_str() {
        "gaussian" => CurvatureSource::Gaussian,
        "uniform" => CurvatureSource::Uniform(float(&curvature_line[1])),
        "per_face" => CurvatureSource::PerFace(floats(&curvature_line[1..])),
        other => panic!("unknown curvature {other}"),
    };
    let (w, beta) = (floats(one("w")), floats(one("beta")));
    let problem = CutProblem {
        energy,
        curvature,
        regions: (0..2)
            .map(|r| RegionWeights {
                distortion: w[r],
                area: beta[r],
            })
            .collect(),
        length_weight: float(&one("alpha")[0]),
        penalty: lines.get("eta").map(|l| floats(&l[0])),
    };
    let flux = one("flux");
    let flux_recovery = match flux[0].as_str() {
        "nodal" => FluxRecovery::Nodal,
        "smoothed" => FluxRecovery::Smoothed {
            length: float(&flux[1]),
        },
        other => panic!("unknown flux recovery {other}"),
    };
    let step = one("step");
    let step_size = match step[0].as_str() {
        "fixed" => StepSize::Fixed(float(&step[1])),
        "capped" => StepSize::Capped {
            tau_max: float(&step[1]),
            max_motion: float(&step[2]),
        },
        "adaptive" => StepSize::Adaptive {
            time_scale: float(&step[1]),
            max_motion: float(&step[2]),
            tau_max: float(&step[3]),
        },
        other => panic!("unknown step rule {other}"),
    };
    let stages = lines.get("stages").map(|l| {
        let w = &l[0];
        (
            match w[0].as_str() {
                "onmesh" => LengthStep::OnMesh,
                "oncut" => LengthStep::OnCut,
                other => panic!("unknown length step {other}"),
            },
            match w[1].as_str() {
                "harmonic" => Extension::Harmonic,
                "local" => Extension::Local,
                other => panic!("unknown extension {other}"),
            },
            match w[2].as_str() {
                "unit" => SpeedScaling::Unit,
                "mean_gradient" => SpeedScaling::MeanGradient,
                other => panic!("unknown speed scaling {other}"),
            },
        )
    });
    let penalty_rule = match one("penalty")[0].as_str() {
        "gradient" => PenaltyRule::Gradient,
        "paper" => PenaltyRule::PaperDiagonal,
        other => panic!("unknown penalty rule {other}"),
    };
    Case {
        mesh,
        geometry,
        phi: floats(one("phi")),
        problem,
        feature,
        eval: EvalOptions {
            min_crossing_offset: float(&one("delta")[0]),
            flux_recovery,
            sliver_aspect_ratio: float(&one("rho0")[0]),
            smoothing_passes: int(&one("passes")[0]),
            parallelism: Parallelism::Sequential,
        },
        step_size,
        penalty_rule,
        classes: lines
            .get("classes")
            .map(|l| l[0].iter().map(|w| int(w) as u32).collect()),
        stages,
        statistics: statistics_of(&lines),
    }
}

/// The flow as SPEC §3.8–§3.12 first specified it and as the oracle implements it.
fn paper_flow(case: &Case) -> FlowOptions<'static> {
    let (length_step, extension, speed_scaling) =
        case.stages
            .unwrap_or((LengthStep::OnMesh, Extension::Harmonic, SpeedScaling::Unit));
    FlowOptions {
        eval: case.eval,
        step_size: case.step_size,
        length_step,
        penalty_rule: case.penalty_rule,
        symmetry: case.classes.clone(),
        speed_scaling,
        extension,
        redistance: Redistance::KeepAdjacent,
        step_statistics: case.statistics,
    }
}

/// Collects disagreements instead of stopping at the first.
struct Report {
    tolerance: f64,
    problems: Vec<String>,
    worst: f64,
    compared: usize,
}

impl Report {
    fn new(tolerance: f64) -> Report {
        Report {
            tolerance,
            problems: Vec::new(),
            worst: 0.0,
            compared: 0,
        }
    }

    fn note(&mut self, message: String) {
        if self.problems.len() < 8 {
            self.problems.push(message);
        } else if self.problems.len() == 8 {
            self.problems.push("...".to_string());
        }
    }

    fn close(&mut self, what: impl Fn() -> String, ours: f64, oracle: f64) {
        let error = (ours - oracle).abs() / oracle.abs().max(1.0);
        self.compared += 1;
        if error.is_nan() || error > self.tolerance {
            self.note(format!("{}: {ours:e} vs oracle {oracle:e}", what()));
        }
        self.worst = if error.is_finite() {
            self.worst.max(error)
        } else {
            f64::INFINITY
        };
    }

    fn same<T: PartialEq + std::fmt::Debug>(
        &mut self,
        what: impl Fn() -> String,
        ours: T,
        oracle: T,
    ) {
        if ours != oracle {
            self.note(format!("{}: {ours:?} vs oracle {oracle:?}", what()));
        }
    }
}

/// The vertex of our cut mesh for the oracle's node `node` of region `region`.
fn our_vertex(e: &Evaluation, region: u32, node: usize) -> usize {
    let nv = e.triangles.n_original;
    if node < nv {
        node
    } else {
        let copies = e.triangles.point_vertices(node - nv);
        copies
            .into_iter()
            .find(|&v| e.triangles.vertex_region[v] == region)
            .expect("a cut point has a copy in each of its regions")
    }
}

/// Everything of the evaluation that is well conditioned.
fn compare_evaluation(
    case: &Case,
    evaluator: &mut Evaluator,
    e: &Evaluation,
    out: &Lines,
    report: &mut Report,
) {
    let nv = case.mesh.n_vertices();
    let tri = &e.triangles;
    let line = |key: &str| -> &Vec<String> { &out[key][0] };

    let scale = line("scale");
    report.close(|| "scale S".into(), e.scale.length_scale, float(&scale[0]));
    report.close(|| "total area".into(), e.scale.total_area, float(&scale[1]));
    let omega = face_curvatures(&case.mesh, &case.geometry).unwrap();
    report.close(
        || "sum_omega".into(),
        omega.iter().sum(),
        float(&line("sum_omega")[0]),
    );
    let counts: Vec<usize> = line("counts").iter().map(|w| int(w)).collect();
    report.same(|| "cut points".into(), e.cut.points.len(), counts[0]);
    report.same(|| "segments".into(), e.cut.segments.len(), counts[1]);
    report.same(|| "triangles".into(), tri.n_faces(), counts[2]);
    report.same(
        || "moved crossings".into(),
        e.report.n_clamped_crossings,
        counts[3],
    );
    report.same(
        || "sliver vertices".into(),
        e.report.n_sliver_points,
        counts[4],
    );
    report.same(
        || "near ties".into(),
        e.report.n_diagonal_near_ties,
        counts[6],
    );

    for w in &out["f"] {
        let f = int(&w[0]);
        report.close(|| format!("Omega of face {f}"), omega[f], float(&w[1]));
        report.close(
            || format!("k of face {f}"),
            evaluator.gaussian_curvature()[f],
            float(&w[2]),
        );
        let area = e.scale.area_to_unit_area(case.geometry.face_areas()[f]);
        report.close(|| format!("area of face {f}"), area, float(&w[4]));
    }
    let weight = case.feature.map(|(r_eps, c)| {
        evaluator
            .feature_size_weight(r_eps, c, Parallelism::Sequential)
            .unwrap()
            .to_vec()
    });
    for w in &out["v"] {
        let i = int(&w[0]);
        report.same(
            || format!("region of vertex {i}"),
            e.cut.vertex_region[i] as usize,
            int(&w[1]),
        );
        report.close(
            || format!("area of vertex {i}"),
            evaluator.vertex_areas()[i],
            float(&w[3]),
        );
        if let Some(a) = &weight {
            report.close(|| format!("weight a of vertex {i}"), a[i], float(&w[5]));
        }
    }
    for w in out.get("x").map_or(&[][..], |l| &l[..]) {
        let k = int(&w[0]);
        let point = e.cut.points[k];
        let [a, b] = case.mesh.edge_vertices(point.edge as usize);
        report.same(
            || format!("edge of cut point {k}"),
            (a, b),
            (int(&w[1]), int(&w[2])),
        );
        report.close(|| format!("t of cut point {k}"), point.t, float(&w[4]));
    }
    // Our vertex -> the oracle's node.
    let node = |v: u32| -> usize {
        let v = v as usize;
        if v < nv { v } else { nv + (v - nv) / 2 }
    };
    for w in &out["t"] {
        let t = int(&w[0]);
        report.same(
            || format!("parent of triangle {t}"),
            tri.parent[t] as usize,
            int(&w[1]),
        );
        report.same(
            || format!("region of triangle {t}"),
            tri.region[t] as usize,
            int(&w[2]),
        );
        let nodes = tri.faces[t].map(node);
        report.same(
            || format!("nodes of triangle {t}"),
            nodes,
            [int(&w[3]), int(&w[4]), int(&w[5])],
        );
        for k in 0..3 {
            report.close(
                || format!("edge {k} of triangle {t}"),
                tri.len2[t][k].sqrt(),
                float(&w[6 + k]),
            );
        }
        report.close(
            || format!("area of triangle {t}"),
            tri.area[t],
            float(&w[9]),
        );
    }
    for (i, w) in out.get("g").map_or(&[][..], |l| &l[..]).iter().enumerate() {
        let segment = e.cut.segments[i];
        let mut ours = segment.points.map(|p| p as usize);
        let mut theirs = [int(&w[0]), int(&w[1])];
        ours.sort_unstable();
        theirs.sort_unstable();
        report.same(|| format!("points of segment {i}"), ours, theirs);
        report.same(
            || format!("face of segment {i}"),
            segment.face as usize,
            int(&w[2]),
        );
        report.close(
            || format!("length of segment {i}"),
            tri.segment_length[i],
            float(&w[3]),
        );
    }
    for w in &out["r"] {
        let r = int(&w[0]);
        report.close(
            || format!("J of region {r}"),
            e.energy.distortion[r],
            float(&w[6]),
        );
        let full = e.energy.area[r] + e.report.excluded_area[r];
        report.close(|| format!("area of region {r}"), full, float(&w[7]));
        report.close(
            || format!("solved area of region {r}"),
            e.energy.area[r],
            float(&w[8]),
        );
    }
    for w in &out["u"] {
        let (r, n) = (int(&w[0]) as u32, int(&w[1]));
        let v = our_vertex(e, r, n);
        report.close(
            || format!("u at node {n} of region {r}"),
            e.log_scale[v],
            float(&w[3]),
        );
        if let Some(p) = &e.adjoint {
            report.close(
                || format!("p at node {n} of region {r}"),
                p[v],
                float(&w[4]),
            );
        }
    }
    for w in &out["b"] {
        let (r, n) = (int(&w[0]) as u32, int(&w[1]));
        let v = our_vertex(e, r, n);
        report.close(
            || format!("v at node {n} of region {r}"),
            e.flux[v],
            float(&w[2]),
        );
        report.close(
            || format!("l at node {n} of region {r}"),
            e.dual_length[v],
            float(&w[3]),
        );
        if let Some(q) = &e.adjoint_flux {
            report.close(
                || format!("q at node {n} of region {r}"),
                q[v],
                float(&w[5]),
            );
        }
    }
    for w in out.get("s").map_or(&[][..], |l| &l[..]) {
        let k = int(&w[0]);
        report.same(
            || format!("sliver flag of cut point {k}"),
            tri.sliver_point[k],
            int(&w[3]) == 1,
        );
    }
    let energy = line("energy");
    for r in 0..2 {
        report.close(
            || format!("energy J{r}"),
            e.energy.distortion[r],
            float(&energy[1 + r]),
        );
        report.close(
            || format!("energy A{r}"),
            e.energy.area[r],
            float(&energy[3 + r]),
        );
    }
    let length = case.problem.length_weight * e.energy.length + e.energy.penalty;
    report.close(|| "energy length term".into(), length, float(&energy[5]));
}

/// The quantities that divide by dual lengths, which can be as small as `delta` times an edge
/// (with the nodal recovery that amplifies rounding errors): the recovered normal derivatives,
/// `G` at the cut vertices, `s`, and the total energy.
fn compare_gradient(case: &Case, e: &Evaluation, out: &Lines, report: &mut Report) {
    let nv = case.mesh.n_vertices();
    for w in &out["b"] {
        let (r, n) = (int(&w[0]) as u32, int(&w[1]));
        let v = our_vertex(e, r, n);
        report.close(
            || format!("d at node {n} of region {r}"),
            e.normal_derivative[v],
            float(&w[4]),
        );
        if let Some(dp) = &e.adjoint_normal_derivative {
            report.close(
                || format!("e at node {n} of region {r}"),
                dp[v],
                float(&w[6]),
            );
        }
        if n < nv {
            continue; // the oracle prints G on the boundary of the mesh too; we do not keep it
        }
        let weights = case.problem.regions[r as usize];
        let oracle = weights.distortion * float(&w[7]) + weights.area;
        report.close(
            || format!("G at node {n} of region {r}"),
            e.gradient[v],
            oracle,
        );
    }
    for w in out.get("s").map_or(&[][..], |l| &l[..]) {
        let k = int(&w[0]);
        report.close(
            || format!("s raw at cut point {k}"),
            e.speed_raw[k],
            float(&w[1]),
        );
        report.close(
            || format!("s smoothed at cut point {k}"),
            e.speed[k],
            float(&w[2]),
        );
    }
    report.close(
        || "energy total".into(),
        e.energy.total,
        float(&out["energy"][0][0]),
    );
}

/// What a replayed step says about the step rule of SPEC §4.
struct StepOutcome {
    /// The fixture asks for the robust statistics.
    robust: bool,
    /// For a fixture with the extremes: whether the robust statistics give another step.
    differs: bool,
    /// Cut points (two regions) or boundary nodes (`N` regions) slowed down by the local limit,
    /// or, where the limit acts at the vertices (SPEC D26), the vertices.
    limited: usize,
    /// The local limit acts at the vertices: the robust statistics with the length step on the
    /// mesh.
    vertex_limit: bool,
}

/// The local limit at the vertices (SPEC §4, D26; the oracle's lines `vlimit` and `hv`): `h_i`,
/// `ω_i`, the factor, whether the vertex is limited, and the count of limited cut-adjacent
/// vertices. `omega` is the roughening rate of the state, which the oracle prints in place of
/// an `ω_i` that is absent or below it. Returns whether the limit is in use.
fn compare_vertex_limit(
    out: &Lines,
    report: &mut Report,
    [edges, rates, factors]: [&[f64]; 3],
    adjacent: &[bool],
    omega: f64,
    n_limited: usize,
) -> bool {
    let in_use = !factors.is_empty();
    report.same(
        || "the local limit acts at the vertices".into(),
        in_use,
        out.contains_key("vlimit"),
    );
    if !in_use || !out.contains_key("vlimit") {
        return in_use;
    }
    report.same(
        || "vertices with a factor".into(),
        factors.len(),
        out["hv"].len(),
    );
    for w in &out["hv"] {
        let i = int(&w[0]);
        report.close(|| format!("h_i of vertex {i}"), edges[i], float(&w[1]));
        // The oracle prints the larger of `ω_i` and `ω`: only the excess over `ω` enters.
        report.close(
            || format!("omega_i of vertex {i}"),
            rates[i].max(omega),
            float(&w[2]),
        );
        report.same(
            || format!("whether vertex {i} has a cut point on its edges"),
            rates[i] > f64::NEG_INFINITY,
            adjacent[i],
        );
        report.close(
            || format!("factor of the local limit at vertex {i}"),
            factors[i],
            float(&w[3]),
        );
        report.same(
            || {
                format!(
                    "whether vertex {i} is limited (factors {:e} and {:e})",
                    factors[i],
                    float(&w[3])
                )
            },
            factors[i] < 1.0,
            float(&w[3]) < 1.0,
        );
    }
    let limited = factors
        .iter()
        .zip(adjacent)
        .filter(|(f, a)| **f < 1.0 && **a)
        .count();
    report.same(
        || "limited cut-adjacent vertices".into(),
        limited,
        int(&out["limited"][0][0]),
    );
    report.same(
        || "limited vertices of the report".into(),
        n_limited,
        limited,
    );
    in_use
}

/// One step of the flow with the stages and the statistics that the fixture names.
fn compare_step(case: &Case, out: &Lines, report: &mut Report) -> StepOutcome {
    let options = paper_flow(case);
    let mut flow = Flow::new(&case.mesh, &case.geometry, &case.phi).unwrap();
    let scale = flow.scale();
    for (w, ours) in out["v"].iter().zip(flow.level_set_unit_area()) {
        report.close(
            || format!("scaled phi of vertex {}", w[0]),
            *ours,
            float(&w[7]),
        );
    }
    // The speed after sliver smoothing, which the statistics of the step rule are taken of.
    let before = flow.evaluate(&case.problem, &case.eval).unwrap();
    let adjacent = before.cut.cut_adjacent_vertices(&case.mesh);
    let speed = before.speed;
    let mut vertex_limit = false;
    let (step, stages) = flow.step_with_stages(&case.problem, &options).unwrap();
    let tau = &out["tau"][0];
    report.close(|| "tau".into(), step.tau, float(&tau[0]));
    let robust = case.statistics == StepStatistics::Robust;
    if case.stages.is_none() {
        report.close(|| "largest speed".into(), step.max_speed, float(&tau[2]));
    } else {
        // Format 2: the speed of the step rule, the roughening rate, and the rates on the cut.
        report.close(
            || "shortest crossed edge".into(),
            step.shortest_crossed_edge,
            float(&tau[1]),
        );
        let on_cut = options.length_step == LengthStep::OnCut;
        // ρ: with the length force where the length step is on the cut.
        let rho: Vec<f64> = if on_cut && !stages.length_speed.is_empty() {
            speed
                .iter()
                .zip(&stages.length_speed)
                .map(|(s, l)| s + l)
                .collect()
        } else {
            speed.clone()
        };
        let rho_max = rho.iter().fold(0.0f64, |m, r| m.max(r.abs()));
        report.close(|| "max |rho|".into(), rho_max, float(&tau[2]));
        let omega = if robust {
            step.roughening_rate
        } else {
            report.close(|| "speed of the step rule".into(), step.step_speed, rho_max);
            step.roughening_rate_max
        };
        report.close(|| "omega".into(), omega, float(&out["omega"][0][0]));
        report.close(
            || "max |lambda0|".into(),
            step.max_length_speed,
            float(&out["lammax"][0][0]),
        );
        if !matches!(case.step_size, StepSize::Fixed(_)) {
            let (motion_speed, edge) = if robust {
                (step.step_speed, 1.0)
            } else {
                (step.step_speed, step.shortest_crossed_edge)
            };
            let ours = step_limits(&case.step_size, omega, motion_speed, edge);
            for (k, name) in ["tau_max", "time limit", "motion limit"].iter().enumerate() {
                close_or_infinite(report, name, ours[k], float(&out["limits"][0][k]));
            }
        }
        if let Some(stats) = out.get("stats") {
            // Format 3: the statistics of the step rule (SPEC §1, §4, D22) and the local limit.
            report.same(
                || "statistics".into(),
                stats[0][0].as_str(),
                if robust { "robust" } else { "extremes" },
            );
            report.close(
                || "largest omega_c".into(),
                step.roughening_rate_max,
                float(&out["omegamax"][0][0]),
            );
            let (h, ell) = (&stages.local_edge, &stages.cut_length);
            let quantile = |values: &[f64]| {
                let in_edges: Vec<f64> = values.iter().zip(h).map(|(r, h)| r.abs() / h).collect();
                weighted_quantile(&in_edges, ell, ROUGHENING_QUANTILE)
            };
            let mhat = &out["mhat"][0];
            let first = quantile(&speed);
            report.close(|| "mhat of s".into(), first, float(&mhat[0]));
            // λ⁰ is there with the length step on the cut, and with the length step on the mesh
            // where the robust rule needs it.
            let with_force: Option<Vec<f64>> = (!stages.length_speed.is_empty()).then(|| {
                speed
                    .iter()
                    .zip(&stages.length_speed)
                    .map(|(s, l)| s + l)
                    .collect()
            });
            assert!(
                with_force.is_some() || !robust,
                "no length force for the step rule"
            );
            let second = with_force.as_deref().map(quantile);
            if let Some(second) = second {
                report.close(|| "mhat of s + lambda0".into(), second, float(&mhat[1]));
            }
            let used = if !robust {
                rho_max / step.shortest_crossed_edge
            } else if on_cut {
                second.unwrap()
            } else {
                first.max(second.unwrap())
            };
            report.close(|| "mhat used".into(), used, float(&mhat[2]));
            if robust {
                report.close(|| "speed of the step rule".into(), step.step_speed, used);
                report.close(
                    || "roughening rate of the state".into(),
                    step.roughening_rate,
                    weighted_quantile(&stages.roughening, ell, ROUGHENING_QUANTILE),
                );
            }
            vertex_limit = compare_vertex_limit(
                out,
                report,
                [
                    &stages.vertex_edge,
                    &stages.vertex_roughening,
                    &stages.vertex_factor,
                ],
                &adjacent,
                omega,
                step.n_limited,
            );
            if vertex_limit {
                // No cut point is limited then, and the limited function is given.
                report.same(
                    || "factors below 1 on the cut".into(),
                    stages.limit_factor.iter().filter(|f| **f < 1.0).count(),
                    0,
                );
                for w in &out["hv"] {
                    let i = int(&w[0]);
                    report.close(
                        || format!("limited phi of vertex {i}"),
                        stages.limited[i],
                        float(&w[4]),
                    );
                }
            } else {
                report.same(
                    || "limited cut points".into(),
                    step.n_limited,
                    int(&out["limited"][0][0]),
                );
                report.same(
                    || "factors below 1".into(),
                    stages.limit_factor.iter().filter(|f| **f < 1.0).count(),
                    int(&out["limited"][0][0]),
                );
            }
            for w in &out["h"] {
                let k = int(&w[0]);
                report.close(|| format!("h_c of cut point {k}"), h[k], float(&w[1]));
                report.close(
                    || format!("|s|/h_c at cut point {k}"),
                    speed[k].abs() / h[k],
                    float(&w[2]),
                );
                if let Some(rho) = &with_force {
                    report.close(
                        || format!("|s + lambda0|/h_c at cut point {k}"),
                        rho[k].abs() / h[k],
                        float(&w[3]),
                    );
                }
                report.close(
                    || format!("factor of the local limit at cut point {k}"),
                    stages.limit_factor[k],
                    float(&w[4]),
                );
                report.same(
                    || {
                        format!(
                            "whether cut point {k} is limited (factors {:e} and {:e})",
                            stages.limit_factor[k],
                            float(&w[4])
                        )
                    },
                    stages.limit_factor[k] < 1.0,
                    float(&w[4]) < 1.0,
                );
                report.close(
                    || format!("rate before the limit at cut point {k}"),
                    stages.unscaled_rate[k],
                    float(&w[5]),
                );
                report.close(
                    || format!("rate after the limit at cut point {k}"),
                    stages.unscaled_rate[k] * stages.limit_factor[k],
                    float(&w[6]),
                );
            }
        }
        let scaled = options.speed_scaling == SpeedScaling::MeanGradient;
        for w in &out["c"] {
            let k = int(&w[0]);
            report.close(
                || format!("length of the cut at cut point {k}"),
                stages.cut_length[k],
                float(&w[2]),
            );
            if !stages.length_force.is_empty() {
                report.close(
                    || format!("Theta at cut point {k}"),
                    stages.length_force[k],
                    float(&w[5]),
                );
            }
            if on_cut {
                report.close(
                    || format!("lambda0 at cut point {k}"),
                    stages.length_speed[k],
                    float(&w[6]),
                );
            }
            report.close(
                || format!("rate at cut point {k}"),
                stages.unscaled_rate[k],
                float(&w[7]),
            );
            if scaled {
                report.close(
                    || format!("mean gradient at cut point {k}"),
                    stages.gradient_scale[k],
                    float(&w[8]),
                );
            }
            report.close(
                || format!("scaled rate at cut point {k}"),
                stages.rate[k],
                float(&w[9]),
            );
            report.close(
                || format!("roughening rate at cut point {k}"),
                stages.roughening[k],
                float(&w[11]),
            );
        }
        for w in &out["m"] {
            let i = int(&w[0]);
            report.close(
                || format!("extension at vertex {i}"),
                stages.extension[i],
                float(&w[2]),
            );
        }
    }
    for w in &out["v"] {
        let i = int(&w[0]);
        report.close(
            || format!("sigma of vertex {i}"),
            stages.extension[i],
            float(&w[8]),
        );
    }
    let new = flow.level_set_unit_area().to_vec();
    let new_caller = flow.level_set();
    for w in &out["n"] {
        let i = int(&w[0]);
        report.close(
            || format!("phi~ of vertex {i}"),
            stages.explicit[i],
            float(&w[1]),
        );
        report.close(
            || format!("phi' of vertex {i}"),
            stages.length[i],
            float(&w[2]),
        );
        report.close(
            || format!("projected phi' of vertex {i}"),
            stages.symmetric[i],
            float(&w[3]),
        );
        report.same(
            || format!("cut-adjacency of vertex {i}"),
            stages.adjacent[i],
            int(&w[5]) == 1,
        );
        report.close(|| format!("new phi of vertex {i}"), new[i], float(&w[7]));
        report.close(
            || format!("new phi (caller's unit) of vertex {i}"),
            new_caller[i] * scale.length_scale,
            float(&w[8]) * scale.length_scale,
        );
    }
    let redist = &out["redist"][0];
    report.same(
        || "unreached vertices".into(),
        step.n_unreached,
        int(&redist[2]),
    );
    report.same(
        || "cut-adjacent vertices".into(),
        stages.adjacent.iter().filter(|&&a| a).count(),
        int(&redist[1]),
    );
    // The energy and the cut of the new state.
    let e = flow.evaluate(&case.problem, &case.eval).unwrap();
    let energy = &out["energy_new"][0];
    report.close(|| "new energy".into(), e.energy.total, float(&energy[0]));
    for r in 0..2 {
        report.close(
            || format!("new J{r}"),
            e.energy.distortion[r],
            float(&energy[1 + r]),
        );
        report.close(
            || format!("new A{r}"),
            e.energy.area[r],
            float(&energy[3 + r]),
        );
    }
    report.same(
        || "new cut points".into(),
        e.cut.points.len(),
        int(&out["cuts_new"][0][0]),
    );
    let other_statistics = FlowOptions {
        step_statistics: StepStatistics::Robust,
        ..options
    };
    let mut other = Flow::new(&case.mesh, &case.geometry, &case.phi).unwrap();
    let again = other.step(&case.problem, &other_statistics).unwrap();
    StepOutcome {
        robust,
        differs: !robust && (again.tau != step.tau || again.n_limited > 0),
        limited: stages.limit_factor.iter().filter(|f| **f < 1.0).count()
            + stages.vertex_factor.iter().filter(|f| **f < 1.0).count(),
        vertex_limit,
    }
}

/// What a set of fixtures covered.
#[derive(Default)]
struct SetSummary {
    fixtures: usize,
    /// Fixtures with a full step.
    stepped: usize,
    nodal: usize,
    smoothed: usize,
    with_triples: usize,
    default_stages: usize,
    /// Numbers compared.
    compared: usize,
    /// Fixtures that ask for the robust statistics, and those of them with a limited point.
    robust: usize,
    limited: usize,
    /// Fixtures in which the local limit acts at the vertices (SPEC D26).
    vertex_limit: usize,
}

fn fixture_names(set: &str, at_least: usize) -> (PathBuf, Vec<String>) {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../lean/oracle/fixtures")
        .join(set);
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .filter_map(|entry| {
            let name = entry.unwrap().file_name().into_string().unwrap();
            name.strip_suffix(".in").map(str::to_string)
        })
        .collect();
    names.sort();
    assert!(
        names.len() >= at_least,
        "only {} fixtures in {}",
        names.len(),
        dir.display()
    );
    (dir, names)
}

#[test]
fn fixtures_of_the_oracle_are_reproduced() {
    let set = replay_two_regions("cuts", 20);
    assert!(set.stepped >= 20 && set.nodal > 0 && set.smoothed > 0 && set.compared > 20_000);
}

/// The fixtures of the default stages and of mixed stage sets (format 2 of the oracle's input):
/// two regions, the stages of SPEC §3.8–§3.10 in any combination, the adaptive step. Format 2
/// has the extremes over the cut in the step rule (SPEC §4 before D22), so the replay asks for
/// [`StepStatistics::Extremes`].
#[test]
fn fixtures_of_the_default_stages_are_reproduced() {
    // `m-featurelap-cotan` smooths the feature size with the cotan Laplacian, which this crate
    // replaced by the Laplacian of the non-negative edge weights (SPEC D17). The two are the
    // same matrix on its mesh, which has no negative edge weight, so it agrees like the rest.
    let set = replay_two_regions("cuts-default", 20);
    assert!(set.stepped >= 20 && set.nodal > 0 && set.smoothed > 0 && set.compared > 20_000);
}

/// The fixtures of the robust step rule (format 3 of the oracle's input, SPEC §1, §4, D22):
/// two regions, the statistics that the fixture names, with the roughening rate as a quantile,
/// both speeds of the step rule, `h_c`, and at every cut point the factor of the local limit
/// and the rate before and after it.
#[test]
fn fixtures_of_the_robust_step_rule_are_reproduced() {
    let set = replay_two_regions("cuts-robust", 9);
    assert!(set.stepped >= 9 && set.robust >= 8 && set.limited >= 3 && set.compared > 20_000);
    assert!(set.vertex_limit >= 2);
}

/// Replays the two-region fixtures of `lean/oracle/fixtures/<set>`.
fn replay_two_regions(set: &str, at_least: usize) -> SetSummary {
    let (dir, names) = fixture_names(set, at_least);
    let mut failures = Vec::new();
    let mut robust_differs: Vec<String> = Vec::new();
    let mut summary = SetSummary {
        fixtures: names.len(),
        ..SetSummary::default()
    };
    let mut worst = [0.0f64; 3];
    let mut compared = [0usize; 3];
    for name in &names {
        let case = parse_case(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let out = lines_by_key(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        let status = out["status"][0][0].as_str();
        assert!(
            status == "ok" || status == "no_cut",
            "{name}: status {status}"
        );
        let mut evaluator = Evaluator::new(&case.mesh, &case.geometry).unwrap();
        let e = evaluator
            .evaluate(&case.phi, &case.problem, &case.eval)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(e.cut.points.is_empty(), status == "no_cut", "{name}");
        // Tolerance 1e-9: what the oracle uses for its own check.
        let mut reports = [Report::new(1e-9), Report::new(1e-9), Report::new(1e-9)];
        compare_evaluation(&case, &mut evaluator, &e, &out, &mut reports[0]);
        compare_gradient(&case, &e, &out, &mut reports[1]);
        if status == "ok" {
            let outcome = compare_step(&case, &out, &mut reports[2]);
            if outcome.differs {
                robust_differs.push(name.clone());
            }
            summary.stepped += 1;
            summary.robust += usize::from(outcome.robust);
            summary.limited += usize::from(outcome.robust && outcome.limited > 0);
            summary.vertex_limit += usize::from(outcome.vertex_limit);
        }
        let problems: Vec<String> = reports
            .iter()
            .flat_map(|r| r.problems.iter().cloned())
            .collect();
        if !problems.is_empty() {
            failures.push(format!("{name}:\n    {}", problems.join("\n    ")));
            continue;
        }
        for (k, r) in reports.iter().enumerate() {
            worst[k] = worst[k].max(r.worst);
            compared[k] += r.compared;
        }
        match case.eval.flux_recovery {
            FluxRecovery::Nodal => summary.nodal += 1,
            FluxRecovery::Smoothed { .. } => summary.smoothed += 1,
        }
    }
    summary.compared = compared.iter().sum();
    println!(
        "{set}: {} fixtures ({} with the nodal flux recovery, {} with the smoothed one; {} with \
         a full step; {} with the robust statistics, {} of them with limited cut points). \
         Numbers compared and largest difference relative to max(1, |value|): evaluation {} / \
         {:.1e}; d, e, G, s {} / {:.1e}; step {} / {:.1e}",
        summary.fixtures,
        summary.nodal,
        summary.smoothed,
        summary.stepped,
        summary.robust,
        summary.limited,
        compared[0],
        worst[0],
        compared[1],
        worst[1],
        compared[2],
        worst[2]
    );
    assert!(
        failures.is_empty(),
        "{} of {} fixtures disagree:\n{}",
        failures.len(),
        names.len(),
        failures.join("\n")
    );
    println!(
        "{set}: of the fixtures with the extremes, the robust statistics of the step rule (SPEC \
         D22) give another step in {}: {}",
        robust_differs.len(),
        robust_differs.join(" ")
    );
    summary
}

/// The fixtures for `N` regions (`lean/oracle/fixtures/regions`, SPEC §8): the evaluation of a
/// network with triple points, and one step with the paper's stages (§8.4) or the default ones
/// (§8.5), stage by stage. These fixtures have the extremes over the cut in the step rule, so
/// the replay asks for [`StepStatistics::Extremes`].
#[test]
fn fixtures_for_several_regions_are_reproduced() {
    let set = replay_region_set("regions", 15);
    assert!(set.with_triples >= 10 && set.default_stages >= 5 && set.compared > 50_000);
}

/// The fixtures of the robust step rule for `N` regions (format 2 of the oracle's `regions`
/// input; SPEC §1, §4, §8.4, §8.5, D22): the quantiles over the cut points and over the
/// boundary nodes, `h` and the weights of both, and at every boundary node the factor of the
/// local limit and the rate before and after it.
#[test]
fn fixtures_of_the_robust_step_rule_for_several_regions_are_reproduced() {
    let set = replay_region_set("regions-robust", 5);
    assert!(set.robust >= 5 && set.limited >= 2 && set.default_stages >= 2);
    assert!(set.with_triples >= 5 && set.compared > 20_000);
    assert!(set.vertex_limit >= 4);
}

/// Replays the fixtures for `N` regions of `lean/oracle/fixtures/<set>`.
fn replay_region_set(set: &str, at_least: usize) -> SetSummary {
    let (dir, names) = fixture_names(set, at_least);
    let mut failures = Vec::new();
    let mut worst = 0.0f64;
    let mut summary = SetSummary {
        fixtures: names.len(),
        ..SetSummary::default()
    };
    let mut robust_differs: Vec<String> = Vec::new();
    for name in &names {
        let input = lines_by_key(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let out = lines_by_key(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        assert_eq!(out["status"][0][0], "ok", "{name}");
        let mut report = Report::new(1e-9);
        let (triples, default, outcome) = replay_regions(&input, &out, &mut report);
        if outcome.differs {
            robust_differs.push(name.clone());
        }
        summary.stepped += 1;
        summary.with_triples += usize::from(triples > 0);
        summary.default_stages += usize::from(default);
        summary.robust += usize::from(outcome.robust);
        summary.limited += usize::from(outcome.robust && outcome.limited > 0);
        summary.vertex_limit += usize::from(outcome.vertex_limit);
        if report.problems.is_empty() {
            worst = worst.max(report.worst);
            summary.compared += report.compared;
        } else {
            failures.push(format!("{name}:\n    {}", report.problems.join("\n    ")));
        }
    }
    println!(
        "{set}: {} fixtures ({} with triple points, {} with the default stages; {} with the \
         robust statistics, {} of them with limited boundary nodes). Numbers compared {}, \
         largest difference relative to max(1, |value|) {worst:.1e}",
        summary.fixtures,
        summary.with_triples,
        summary.default_stages,
        summary.robust,
        summary.limited,
        summary.compared
    );
    assert!(
        failures.is_empty(),
        "{} of {} fixtures disagree:\n{}",
        failures.len(),
        names.len(),
        failures.join("\n")
    );
    println!(
        "{set}: of the fixtures with the extremes, the robust statistics of the step rule (SPEC \
         D22) give another step in {}: {}",
        robust_differs.len(),
        robust_differs.join(" ")
    );
    summary
}

/// One fixture for `N` regions. Returns the number of triple points, whether the stages are
/// the default ones, and what the step says about the step rule.
fn replay_regions(input: &Lines, out: &Lines, report: &mut Report) -> (usize, bool, StepOutcome) {
    let one = |key: &str| -> &Vec<String> {
        let all = input
            .get(key)
            .unwrap_or_else(|| panic!("missing line {key}"));
        assert_eq!(all.len(), 1, "line {key} occurs {} times", all.len());
        &all[0]
    };
    let floats = |words: &[String]| -> Vec<f64> { words.iter().map(|w| float(w)).collect() };
    let nv = int(&one("mesh")[0]);
    let faces: Vec<[u32; 3]> = input["f"]
        .iter()
        .map(|w| [0, 1, 2].map(|k| int(&w[k]) as u32))
        .collect();
    let mesh = TriMesh::new(nv, faces).unwrap();
    let mut lengths = vec![f64::NAN; mesh.n_edges()];
    for w in &input["e"] {
        lengths[mesh.find_edge(int(&w[0]), int(&w[1])).unwrap()] = float(&w[2]);
    }
    let geometry = Geometry::from_edge_lengths(&mesh, lengths).unwrap();
    let n = int(&one("regions")[0]);
    let mut phis = vec![Vec::new(); n];
    for w in &input["phi"] {
        phis[int(&w[0])] = floats(&w[1..]);
    }
    let energy_line = one("energy");
    let energy = match energy_line[0].as_str() {
        "dirichlet" => Energy::Dirichlet,
        "hencky" => Energy::Hencky,
        "scaled_a" => Energy::ScaledDirichlet(ScaleWeight::PerVertex(floats(&energy_line[1..]))),
        // The feature size smoothed with the cotan Laplacian is not in this crate (SPEC D17:
        // the Laplacian of the non-negative edge weights replaced it). Where the fixture asks
        // for it, the weight `a` that the oracle prints is taken as given, and the rest of the
        // case is compared.
        "scaled_curv" if one("featurelap")[0] == "cotan" => {
            let mut a = vec![f64::NAN; nv];
            for w in &out["v"] {
                a[int(&w[0])] = float(&w[5]);
            }
            Energy::ScaledDirichlet(ScaleWeight::PerVertex(a))
        }
        "scaled_curv" => Energy::ScaledDirichlet(ScaleWeight::FeatureSize {
            r_eps: float(&energy_line[1]),
            smoothing: float(&energy_line[2]),
        }),
        other => panic!("unknown energy {other}"),
    };
    let curvature_line = one("curvature");
    let curvature = match curvature_line[0].as_str() {
        "gaussian" => CurvatureSource::Gaussian,
        "uniform" => CurvatureSource::Uniform(float(&curvature_line[1])),
        "per_face" => CurvatureSource::PerFace(floats(&curvature_line[1..])),
        other => panic!("unknown curvature {other}"),
    };
    let (w, beta) = (floats(one("w")), floats(one("beta")));
    let problem = CutProblem {
        energy,
        curvature,
        regions: (0..n)
            .map(|r| RegionWeights {
                distortion: w[r],
                area: beta[r],
            })
            .collect(),
        length_weight: float(&one("alpha")[0]),
        penalty: input.get("eta").map(|l| floats(&l[0])),
    };
    let flux = one("flux");
    let step = one("step");
    let default = one("stages")[0] == "default";
    let base = if default {
        FlowOptions::new(StepSize::Fixed(1.0))
    } else {
        FlowOptions::for_regions(StepSize::Fixed(1.0))
    };
    let options = FlowOptions {
        eval: EvalOptions {
            min_crossing_offset: float(&one("delta")[0]),
            flux_recovery: match flux[0].as_str() {
                "nodal" => FluxRecovery::Nodal,
                _ => FluxRecovery::Smoothed {
                    length: float(&flux[1]),
                },
            },
            smoothing_passes: 0,
            parallelism: Parallelism::Sequential,
            ..EvalOptions::default()
        },
        step_size: match step[0].as_str() {
            "fixed" => StepSize::Fixed(float(&step[1])),
            "capped" => StepSize::Capped {
                tau_max: float(&step[1]),
                max_motion: float(&step[2]),
            },
            _ => StepSize::Adaptive {
                time_scale: float(&step[1]),
                max_motion: float(&step[2]),
                tau_max: float(&step[3]),
            },
        },
        penalty_rule: match one("penalty")[0].as_str() {
            "gradient" => PenaltyRule::Gradient,
            _ => PenaltyRule::PaperDiagonal,
        },
        symmetry: input
            .get("classes")
            .map(|l| l[0].iter().map(|w| int(w) as u32).collect()),
        step_statistics: statistics_of(input),
        ..base
    };
    let robust = options.step_statistics == StepStatistics::Robust;

    let mut flow = ddg_variational_cuts::RegionFlow::new(&mesh, &geometry, &phis).unwrap();
    let scale = flow.scale();
    report.close(
        || "scale S".into(),
        scale.length_scale,
        float(&out["scale"][0][0]),
    );
    let e = flow.evaluate(&problem, &options.eval).unwrap();
    let (cut, tri) = (&e.cut, &e.triangles);
    let n_points = cut.points.len();

    // The network.
    let counts = &out["counts"][0];
    report.same(|| "cut points".into(), n_points, int(&counts[0]));
    report.same(
        || "triple points".into(),
        cut.triple_points.len(),
        int(&counts[1]),
    );
    report.same(|| "segments".into(), cut.segments.len(), int(&counts[2]));
    report.same(
        || "cut-mesh triangles".into(),
        tri.n_faces(),
        int(&counts[3]),
    );
    report.same(|| "moved crossings".into(), cut.n_clamped, int(&counts[4]));
    report.same(
        || "diagonal ties".into(),
        tri.n_diagonal_near_ties,
        int(&counts[5]),
    );
    if !report.problems.is_empty() {
        let outcome = StepOutcome {
            robust,
            differs: false,
            limited: 0,
            vertex_limit: false,
        };
        return (cut.triple_points.len(), default, outcome);
    }
    for w in &out["v"] {
        let i = int(&w[0]);
        report.same(
            || format!("region of vertex {i}"),
            cut.vertex_region[i] as usize,
            int(&w[1]),
        );
    }
    for w in &out["p"] {
        let i = int(&w[0]);
        for k in 0..n {
            report.close(
                || format!("scaled phi_{k} of vertex {i}"),
                flow.level_sets_unit_area()[k][i],
                float(&w[1 + k]),
            );
        }
    }
    for w in &out["x"] {
        let p = int(&w[0]);
        let [a, b] = mesh.edge_vertices(cut.points[p].edge as usize);
        report.same(
            || format!("edge of cut point {p}"),
            (a, b),
            (int(&w[1]), int(&w[2])),
        );
        report.close(
            || format!("t of cut point {p}"),
            cut.points[p].t,
            float(&w[4]),
        );
    }
    for w in out.get("q").map_or(&[][..], |l| &l[..]) {
        let m = int(&w[0]);
        let triple = &cut.triple_points[m];
        report.same(
            || format!("face of triple point {m}"),
            triple.face as usize,
            int(&w[1]),
        );
        for k in 0..3 {
            report.same(
                || format!("corner {k} of triple point {m}"),
                mesh.faces()[triple.face as usize][k] as usize,
                int(&w[2 + k]),
            );
            report.same(
                || format!("cut point {k} of triple point {m}"),
                triple.points[k] as usize,
                int(&w[5 + k]),
            );
            report.close(
                || format!("barycentric coordinate {k} of triple point {m}"),
                triple.barycentric[k],
                float(&w[8 + k]),
            );
        }
    }
    for (i, w) in out["g"].iter().enumerate() {
        let segment = &cut.segments[i];
        let mut ours = segment.points.map(|x| x as usize);
        let mut theirs = [int(&w[0]), int(&w[1])];
        ours.sort_unstable();
        theirs.sort_unstable();
        report.same(|| format!("ends of segment {i}"), ours, theirs);
        report.same(
            || format!("face of segment {i}"),
            segment.face as usize,
            int(&w[2]),
        );
        report.close(
            || format!("length of segment {i}"),
            tri.segment_length[i],
            float(&w[3]),
        );
    }
    for (t, w) in out["t"].iter().enumerate() {
        report.same(
            || format!("parent of triangle {t}"),
            tri.parent[t] as usize,
            int(&w[1]),
        );
        report.same(
            || format!("region of triangle {t}"),
            tri.region[t] as usize,
            int(&w[2]),
        );
        report.same(
            || format!("nodes of triangle {t}"),
            tri.faces[t].map(|v| v as usize),
            [int(&w[3]), int(&w[4]), int(&w[5])],
        );
        report.close(
            || format!("area of triangle {t}"),
            tri.area[t],
            float(&w[9]),
        );
    }

    // The solves, the normal derivatives, the speed.
    for w in &out["u"] {
        let (x, node) = (int(&w[0]), int(&w[1]));
        if int(&w[2]) == 1 {
            report.close(
                || format!("u at node {node} of region {x}"),
                e.log_scale[node],
                float(&w[3]),
            );
            if let Some(p) = &e.adjoint {
                report.close(
                    || format!("p at node {node} of region {x}"),
                    p[node],
                    float(&w[4]),
                );
            }
        }
    }
    for w in &out["b"] {
        let (x, node) = (int(&w[0]), int(&w[1]));
        report.close(
            || format!("v at node {node} of region {x}"),
            e.flux[node],
            float(&w[2]),
        );
        report.close(
            || format!("ell at node {node} of region {x}"),
            e.dual_length[node],
            float(&w[3]),
        );
        report.close(
            || format!("d at node {node} of region {x}"),
            e.normal_derivative[node],
            float(&w[4]),
        );
        if let Some(dp) = &e.adjoint_normal_derivative {
            report.close(
                || format!("e at node {node} of region {x}"),
                dp[node],
                float(&w[6]),
            );
        }
        if node >= nv {
            let weights = problem.regions[x];
            report.close(
                || format!("G at node {node} of region {x}"),
                e.gradient[node],
                weights.distortion * float(&w[7]) + weights.area,
            );
        }
    }
    for w in &out["s"] {
        let p = int(&w[0]);
        report.close(|| format!("s at cut point {p}"), e.speed[p], float(&w[1]));
    }
    let energy = &out["energy"][0];
    report.close(|| "energy".into(), e.energy.total, float(&energy[0]));
    // The oracle's second number is the length term of the energy.
    let length_term = |e: &Evaluation| problem.length_weight * e.energy.length + e.energy.penalty;
    report.close(|| "length term".into(), length_term(&e), float(&energy[1]));
    for w in &out["er"] {
        let x = int(&w[0]);
        report.close(
            || format!("J of region {x}"),
            e.energy.distortion[x],
            float(&w[1]),
        );
        report.close(
            || format!("area of region {x}"),
            e.energy.area[x],
            float(&w[2]),
        );
    }

    // The step.
    let (step, stages) = flow.step_with_stages(&problem, &options).unwrap();
    let tau = &out["tau"][0];
    report.close(|| "tau".into(), step.tau, float(&tau[0]));
    report.close(
        || "shortest crossed edge".into(),
        step.shortest_crossed_edge,
        float(&tau[1]),
    );
    let mut vertex_limit = false;
    // ρ: `s` at the cut points for the paper's stages, `s_v + λ⁰_v` at the boundary nodes for
    // the default ones.
    let with_force: Option<Vec<f64>> = (!stages.length_speed.is_empty()).then(|| {
        stages
            .length_speed
            .iter()
            .enumerate()
            .map(|(i, l)| stages.boundary_speed[nv + i] + l)
            .collect()
    });
    let largest = |values: &[f64]| values.iter().fold(0.0f64, |m, r| m.max(r.abs()));
    let rho_max = match (&with_force, default) {
        (Some(rho), true) => largest(rho),
        _ => largest(&e.speed),
    };
    report.close(|| "max |rho|".into(), rho_max, float(&tau[2]));
    let omega = if robust {
        step.roughening_rate
    } else {
        report.close(|| "speed of the step rule".into(), step.step_speed, rho_max);
        step.roughening_rate_max
    };
    report.close(|| "omega".into(), omega, float(&out["omega"][0][0]));
    if !matches!(options.step_size, StepSize::Fixed(_)) {
        let edge = if robust {
            1.0
        } else {
            step.shortest_crossed_edge
        };
        let ours = step_limits(&options.step_size, omega, step.step_speed, edge);
        for (k, name) in ["tau_max", "time limit", "motion limit"].iter().enumerate() {
            close_or_infinite(report, name, ours[k], float(&out["limits"][0][k]));
        }
    }
    for w in &out["s"] {
        let p = int(&w[0]);
        report.close(
            || format!("omega at cut point {p}"),
            stages.roughening[p],
            float(&w[2]),
        );
    }
    if let Some(stats) = out.get("stats") {
        // Format 2 of `regions`: the statistics of the step rule (SPEC §1, §4, §8.4, D22) and
        // the local limit at every boundary node.
        report.same(
            || "statistics".into(),
            stats[0][0].as_str(),
            if robust { "robust" } else { "extremes" },
        );
        report.close(
            || "largest omega_c".into(),
            step.roughening_rate_max,
            float(&out["omegamax"][0][0]),
        );
        for w in &out["hp"] {
            let p = int(&w[0]);
            report.close(
                || format!("h_c of cut point {p}"),
                stages.local_edge[p],
                float(&w[1]),
            );
            report.close(
                || format!("length of the cut at cut point {p}"),
                stages.cut_length[p],
                float(&w[2]),
            );
            report.close(
                || format!("omega_c of cut point {p}"),
                stages.roughening[p],
                float(&w[3]),
            );
        }
        let mhat = &out["mhat"][0];
        if robust {
            let in_edges: Vec<f64> = e
                .speed
                .iter()
                .zip(&stages.local_edge)
                .map(|(s, h)| s.abs() / h)
                .collect();
            let first = weighted_quantile(&in_edges, &stages.cut_length, ROUGHENING_QUANTILE);
            report.close(|| "mhat of s".into(), first, float(&mhat[0]));
            let rho = with_force
                .as_ref()
                .expect("the length force for the step rule");
            let in_edges: Vec<f64> = rho
                .iter()
                .zip(&stages.node_edge)
                .map(|(r, h)| r.abs() / h)
                .collect();
            let second = weighted_quantile(&in_edges, &stages.node_length, ROUGHENING_QUANTILE);
            report.close(|| "mhat of s + lambda0".into(), second, float(&mhat[1]));
            let used = if default { second } else { first.max(second) };
            report.close(|| "mhat used".into(), used, float(&mhat[2]));
            report.close(|| "speed of the step rule".into(), step.step_speed, used);
            report.close(
                || "roughening rate of the state".into(),
                step.roughening_rate,
                weighted_quantile(&stages.roughening, &stages.cut_length, ROUGHENING_QUANTILE),
            );
        }
        vertex_limit = compare_vertex_limit(
            out,
            report,
            [
                &stages.vertex_edge,
                &stages.vertex_roughening,
                &stages.vertex_factor,
            ],
            &cut.cut_adjacent_vertices(&mesh),
            omega,
            step.n_limited,
        );
        if vertex_limit {
            report.same(
                || "factors below 1 at the boundary nodes".into(),
                stages.limit_factor.iter().filter(|f| **f < 1.0).count(),
                0,
            );
        } else {
            report.same(
                || "limited boundary nodes".into(),
                stages.limit_factor.iter().filter(|f| **f < 1.0).count(),
                int(&out["limited"][0][0]),
            );
            // The report counts the cut points and triple points with a limited copy.
            let n_triples = cut.triple_points.len();
            let mut limited = vec![false; n_points + n_triples];
            for (i, factor) in stages.limit_factor.iter().enumerate() {
                if *factor < 1.0 {
                    let point = if i < 2 * n_points {
                        i / 2
                    } else {
                        n_points + (i - 2 * n_points) / 3
                    };
                    limited[point] = true;
                }
            }
            report.same(
                || "limited cut points and triple points".into(),
                step.n_limited,
                limited.iter().filter(|l| **l).count(),
            );
        }
        for w in &out["hc"] {
            let node = int(&w[1]);
            report.close(
                || format!("h of node {node}"),
                stages.node_edge[node - nv],
                float(&w[2]),
            );
            report.close(
                || format!("factor of the local limit at node {node}"),
                stages.limit_factor[node - nv],
                float(&w[3]),
            );
            report.same(
                || {
                    format!(
                        "whether node {node} is limited (factors {:e} and {:e}, 1 - factor {:e})",
                        stages.limit_factor[node - nv],
                        float(&w[3]),
                        1.0 - stages.limit_factor[node - nv]
                    )
                },
                stages.limit_factor[node - nv] < 1.0,
                float(&w[3]) < 1.0,
            );
            report.close(
                || format!("rate before the limit at node {node}"),
                stages.unlimited_rate[node],
                float(&w[4]),
            );
            report.close(
                || format!("rate after the limit at node {node}"),
                stages.rate[node],
                float(&w[5]),
            );
        }
    }
    for w in &out["sn"] {
        let node = int(&w[1]);
        report.close(
            || format!("s_v at node {node}"),
            stages.boundary_speed[node],
            float(&w[2]),
        );
    }
    for w in out.get("d").map_or(&[][..], |l| &l[..]) {
        let node = int(&w[1]);
        report.close(
            || format!("length of the cut at node {node}"),
            stages.node_length[node - nv],
            float(&w[2]),
        );
        report.close(
            || format!("Theta at node {node}"),
            stages.length_force[node - nv],
            float(&w[6]),
        );
        report.close(
            || format!("lambda0 at node {node}"),
            stages.length_speed[node - nv],
            float(&w[7]),
        );
        report.close(
            || format!("rate at node {node}"),
            stages.unlimited_rate[node],
            float(&w[8]),
        );
        if node < nv + 2 * n_points {
            report.close(
                || format!("mean gradient at node {node}"),
                stages.gradient_scale[(node - nv) / 2],
                float(&w[9]),
            );
        }
        report.close(
            || format!("scaled rate at node {node}"),
            stages.scaled_rate[node],
            float(&w[10]),
        );
    }
    for w in out.get("ev").map_or(&[][..], |l| &l[..]) {
        let i = int(&w[0]);
        report.close(
            || format!("extension at vertex {i}"),
            stages.extension[i],
            float(&w[1]),
        );
    }
    for w in &out["pd"] {
        let (k, i) = (int(&w[0]), int(&w[1]));
        // The rate is recovered from the explicit step, so it carries the rounding of phi/tau.
        let (ours, theirs) = (stages.function_rates[k][i], float(&w[2]));
        let slack = 1e-9 * (1.0 + flow.level_sets_unit_area()[k][i].abs() / step.tau);
        report.compared += 1;
        if (ours - theirs).abs() > slack.max(1e-9 * theirs.abs()) {
            report.note(format!(
                "rate of phi_{k} at vertex {i}: {ours:e} vs oracle {theirs:e}"
            ));
        }
    }
    let new = flow.level_sets_unit_area().to_vec();
    for w in &out["n"] {
        let (k, i) = (int(&w[0]), int(&w[1]));
        report.close(
            || format!("phi~_{k} of vertex {i}"),
            stages.explicit[k][i],
            float(&w[2]),
        );
        report.close(
            || format!("phi'_{k} of vertex {i}"),
            stages.length[k][i],
            float(&w[3]),
        );
        report.close(
            || format!("symmetric phi_{k} of vertex {i}"),
            stages.symmetric[k][i],
            float(&w[4]),
        );
        report.close(
            || format!("projected phi_{k} of vertex {i}"),
            stages.projected[k][i],
            float(&w[5]),
        );
        report.close(
            || format!("new phi_{k} of vertex {i}"),
            new[k][i],
            float(&w[7]),
        );
    }
    for w in &out["a"] {
        let i = int(&w[0]);
        report.same(
            || format!("cut-adjacency of vertex {i}"),
            stages.adjacent[i],
            int(&w[2]) == 1,
        );
    }
    report.same(
        || "unreached vertices".into(),
        step.n_unreached,
        int(&out["redist"][0][1]),
    );
    report.close(
        || "cut motion".into(),
        step.cut_motion,
        float(&out["motion"][0][0]),
    );
    let after = flow.evaluate(&problem, &options.eval).unwrap();
    let energy = &out["energy_new"][0];
    report.close(
        || "new energy".into(),
        after.energy.total,
        float(&energy[0]),
    );
    report.close(
        || "new length term".into(),
        length_term(&after),
        float(&energy[1]),
    );
    let cuts = &out["cuts_new"][0];
    report.same(
        || "new cut points".into(),
        after.cut.points.len(),
        int(&cuts[0]),
    );
    report.same(
        || "new triple points".into(),
        after.cut.triple_points.len(),
        int(&cuts[1]),
    );
    let other_statistics = FlowOptions {
        step_statistics: StepStatistics::Robust,
        ..options
    };
    let mut other = ddg_variational_cuts::RegionFlow::new(&mesh, &geometry, &phis).unwrap();
    let again = other.step(&problem, &other_statistics).unwrap();
    let outcome = StepOutcome {
        robust,
        differs: !robust && (again.tau != step.tau || again.n_limited > 0),
        limited: stages.limit_factor.iter().filter(|f| **f < 1.0).count()
            + stages.vertex_factor.iter().filter(|f| **f < 1.0).count(),
        vertex_limit,
    };
    (cut.triple_points.len(), default, outcome)
}

/// A case that the oracle's own fixtures do not have: a region that is a single vertex whose
/// neighbours alternate between two other regions, so that its boundary is a closed chain in
/// which every segment ends at a triple point. The input is written here, answered by the
/// oracle's executable (`lean/oracle/.lake/build/bin/ddg-oracle`, built with `lake build`),
/// and replayed like a fixture, with the robust statistics and either stage set.
#[test]
#[ignore = "needs the built oracle; run with --ignored"]
fn a_boundary_through_triple_points_only_agrees_with_the_oracle() {
    let oracle = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../lean/oracle/.lake/build/bin/ddg-oracle");
    assert!(oracle.exists(), "{} is not built", oracle.display());
    let sphere = ddg_core::generators::icosphere(2).unwrap();
    let nv = sphere.positions.len();
    let mesh = TriMesh::new(nv, sphere.faces.clone()).unwrap();
    let geometry = Geometry::from_positions(&mesh, &sphere.positions).unwrap();
    let faces = mesh.faces();
    let around = |v: u32| faces.iter().filter(|f| f.contains(&v)).count();
    let centre = (0..nv as u32).find(|&v| around(v) == 6).unwrap();
    let mut next = HashMap::new();
    for f in faces.iter().filter(|f| f.contains(&centre)) {
        let k = f.iter().position(|&v| v == centre).unwrap();
        next.insert(f[(k + 1) % 3], f[(k + 2) % 3]);
    }
    let mut labels = vec![1u32; nv];
    labels[centre as usize] = 0;
    let mut v = *next.keys().min().unwrap();
    for k in 0..6 {
        labels[v as usize] = 1 + (k % 2);
        v = next[&v];
    }
    let phis = ddg_variational_cuts::init::regions_from_vertex_labels(&mesh, &geometry, &labels, 3)
        .unwrap();
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    for stages in ["default", "paper"] {
        let name = format!("triple-chain-{stages}");
        let mut text = format!(
            "case {name}\ntask regions\nmesh {nv} {} {}\n",
            mesh.n_faces(),
            mesh.n_edges()
        );
        for f in faces {
            text += &format!("f {} {} {}\n", f[0], f[1], f[2]);
        }
        for (edge, length) in mesh.edges().iter().zip(geometry.edge_lengths()) {
            text += &format!("e {} {} {}\n", edge[0], edge[1], hex(*length));
        }
        text += "regions 3\n";
        for (k, phi) in phis.iter().enumerate() {
            let values: Vec<String> = phi.iter().map(|x| hex(*x)).collect();
            text += &format!("phi {k} {}\n", values.join(" "));
        }
        text += &format!(
            "energy dirichlet\ncurvature gaussian\nalpha {}\npenalty gradient\nw {one} {one} {one}\n\
             beta {zero} {zero} {zero}\nstep adaptive {} {} {one}\nflux smoothed {}\ndelta {}\n\
             stages {stages}\nfeaturelap positive\nformat 2\nstatistics robust\nend\n",
            hex(1.0),
            hex(0.05),
            hex(0.25),
            hex(0.5),
            hex(1e-3),
            one = hex(1.0),
            zero = hex(0.0),
        );
        let (input_path, out_path) = (
            dir.join(format!("{name}.in")),
            dir.join(format!("{name}.out")),
        );
        std::fs::write(&input_path, &text).unwrap();
        let status = std::process::Command::new(&oracle)
            .arg("regions")
            .arg(&input_path)
            .arg(&out_path)
            .status()
            .unwrap();
        assert!(status.success());
        let input = lines_by_key(&text);
        let out = lines_by_key(&std::fs::read_to_string(&out_path).unwrap());
        assert_eq!(out["status"][0][0], "ok", "{name}: {:?}", out["status"][0]);
        let mut report = Report::new(1e-9);
        let (triples, _, outcome) = replay_regions(&input, &out, &mut report);
        assert_eq!(triples, 6);
        println!(
            "{name}: {} numbers compared, largest difference {:.1e}, {} boundary nodes limited; \
             input and output in {}",
            report.compared,
            report.worst,
            outcome.limited,
            dir.display()
        );
        assert!(
            report.problems.is_empty(),
            "{name}:\n    {}",
            report.problems.join("\n    ")
        );
    }
}
