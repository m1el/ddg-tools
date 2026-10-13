//! Differential test against the Lean oracle (`lean/oracle`, an independent dense implementation
//! written from `SPEC.md` alone): the committed fixtures `lean/oracle/fixtures/cuts/*.in|out`
//! are replayed.
//!
//! **Evaluation** (SPEC §3.1–§3.7, §6): the scaling, the face curvatures, the crossings, every
//! triangle of the cut mesh, the segments, `u`, `p`, `v`, `ℓ`, `q`, the recovered `d`, `e` at
//! every node, `G`, the speed `s` before and after sliver smoothing, `J` and the areas per
//! region, and the energy. Every case states its flux recovery (`flux nodal` or
//! `flux smoothed C`).
//!
//! **One step** (SPEC §3.8–§3.12, §4), with the stages as the oracle implements them: the
//! harmonic extension with the cotan Laplacian, the explicit step, the length step on the mesh
//! with either penalty rule, the symmetry projection, the redistancing. These are the options
//! [`paper_flow`] below; the defaults of this crate differ in the extension and in the length
//! step (SPEC D14, D15) and are tested against exact solutions instead. Compared: `τ`, `σ`,
//! `φ̃`, `φ′`, the projected `φ′`, the cut-adjacent vertices, the new level-set function, the
//! number of unreached vertices, and the energy and number of cut points of the new state.

use std::collections::HashMap;
use std::path::PathBuf;

use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CurvatureSource, CutProblem, Energy, EvalOptions, Evaluation, Evaluator, Extension, Flow,
    FlowOptions, FluxRecovery, LengthStep, Parallelism, PenaltyRule, Redistance, RegionWeights,
    ScaleWeight, SpeedScaling, StepSize, StepStatistics, face_curvatures,
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
        step_statistics: StepStatistics::Extremes,
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

/// One step of the flow with the stages of the oracle.
/// Returns whether the step is another one with the robust statistics of SPEC §4 (D22).
fn compare_step(case: &Case, out: &Lines, report: &mut Report) -> bool {
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
    let (step, stages) = flow.step_with_stages(&case.problem, &options).unwrap();
    let tau = &out["tau"][0];
    report.close(|| "tau".into(), step.tau, float(&tau[0]));
    if case.stages.is_none() {
        report.close(|| "largest speed".into(), step.max_speed, float(&tau[2]));
    } else {
        // Format 2: the speed of the step rule, the roughening rate, and the rates on the cut.
        report.close(
            || "shortest crossed edge".into(),
            step.shortest_crossed_edge,
            float(&tau[1]),
        );
        report.close(|| "max |rho|".into(), step.step_speed, float(&tau[2]));
        report.close(
            || "omega".into(),
            step.roughening_rate_max,
            float(&out["omega"][0][0]),
        );
        report.close(
            || "max |lambda0|".into(),
            step.max_length_speed,
            float(&out["lammax"][0][0]),
        );
        let on_cut = !stages.length_force.is_empty();
        let scaled = options.speed_scaling == SpeedScaling::MeanGradient;
        for w in &out["c"] {
            let k = int(&w[0]);
            if on_cut {
                report.close(
                    || format!("Theta at cut point {k}"),
                    stages.length_force[k],
                    float(&w[5]),
                );
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
    let robust = FlowOptions {
        step_statistics: StepStatistics::Robust,
        ..options
    };
    let mut other = Flow::new(&case.mesh, &case.geometry, &case.phi).unwrap();
    let again = other.step(&case.problem, &robust).unwrap();
    again.tau != step.tau || again.n_limited > 0
}

#[test]
fn fixtures_of_the_oracle_are_reproduced() {
    let stepped = replay_two_regions("cuts", 20, &[]);
    assert!(stepped >= 20);
}

/// The fixtures of the default stages and of mixed stage sets (format 2 of the oracle's input):
/// two regions, the stages of SPEC §3.8–§3.10 in any combination, the adaptive step. The oracle
/// takes the extremes over the cut in the step rule (SPEC §4 before D22), so the replay asks
/// for [`StepStatistics::Extremes`].
#[test]
fn fixtures_of_the_default_stages_are_reproduced() {
    // `m-featurelap-cotan` smooths the feature size with the cotan Laplacian, which this crate
    // replaced by the Laplacian of the non-negative edge weights (SPEC D17). The two are the
    // same matrix on its mesh, which has no negative edge weight, so it agrees like the rest.
    let stepped = replay_two_regions("cuts-default", 20, &[]);
    assert!(stepped >= 20);
}

/// Replays the two-region fixtures of `lean/oracle/fixtures/<set>`; returns how many had a
/// full step. `expected_to_differ` names fixtures whose disagreement is known and explained.
fn replay_two_regions(set: &str, at_least: usize, expected_to_differ: &[&str]) -> usize {
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
    let mut failures = Vec::new();
    let mut robust_differs: Vec<String> = Vec::new();
    let (mut nodal, mut smoothed, mut stepped) = (0, 0, 0);
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
            if compare_step(&case, &out, &mut reports[2]) {
                robust_differs.push(name.clone());
            }
            stepped += 1;
        }
        let problems: Vec<String> = reports
            .iter()
            .flat_map(|r| r.problems.iter().cloned())
            .collect();
        if !problems.is_empty() {
            if expected_to_differ.contains(&name.as_str()) {
                println!(
                    "{name} differs as expected ({} disagreements, first: {})",
                    problems.len(),
                    problems[0]
                );
            } else {
                failures.push(format!("{name}:\n    {}", problems.join("\n    ")));
            }
            continue;
        }
        assert!(
            !expected_to_differ.contains(&name.as_str()),
            "{name} was expected to differ and agrees"
        );
        for (k, r) in reports.iter().enumerate() {
            worst[k] = worst[k].max(r.worst);
            compared[k] += r.compared;
        }
        match case.eval.flux_recovery {
            FluxRecovery::Nodal => nodal += 1,
            FluxRecovery::Smoothed { .. } => smoothed += 1,
        }
    }
    println!(
        "{set}: {} fixtures ({nodal} with the nodal flux recovery, {smoothed} with the smoothed one; \
         {stepped} with a full step). Numbers compared and largest difference relative to \
         max(1, |value|): evaluation {} / {:.1e}; d, e, G, s {} / {:.1e}; step {} / {:.1e}",
        names.len(),
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
        "{set}: with the robust statistics of the step rule (SPEC D22) the step differs in {} \
         fixtures: {}",
        robust_differs.len(),
        robust_differs.join(" ")
    );
    assert!(nodal > 0 && smoothed > 0);
    assert!(compared.iter().sum::<usize>() > 20_000);
    stepped
}

/// The fixtures for `N` regions (`lean/oracle/fixtures/regions`, SPEC §8): the evaluation of a
/// network with triple points, and one step with the paper's stages (§8.4) or the default ones
/// (§8.5), stage by stage. The oracle takes the extremes over the cut in the step rule, so the
/// replay asks for [`StepStatistics::Extremes`].
#[test]
fn fixtures_for_several_regions_are_reproduced() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/regions");
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .filter_map(|entry| {
            let name = entry.unwrap().file_name().into_string().unwrap();
            name.strip_suffix(".in").map(str::to_string)
        })
        .collect();
    names.sort();
    assert!(names.len() >= 15, "only {} fixtures", names.len());
    let mut failures = Vec::new();
    let (mut worst, mut compared, mut with_triples, mut default_stages) = (0.0f64, 0usize, 0, 0);
    let mut robust_differs: Vec<String> = Vec::new();
    for name in &names {
        let input = lines_by_key(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let out = lines_by_key(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        assert_eq!(out["status"][0][0], "ok", "{name}");
        let mut report = Report::new(1e-9);
        let (triples, default, differs) = replay_regions(&input, &out, &mut report);
        if differs {
            robust_differs.push(name.clone());
        }
        with_triples += usize::from(triples > 0);
        default_stages += usize::from(default);
        if report.problems.is_empty() {
            worst = worst.max(report.worst);
            compared += report.compared;
        } else {
            failures.push(format!("{name}:\n    {}", report.problems.join("\n    ")));
        }
    }
    println!(
        "regions: {} fixtures ({with_triples} with triple points, {default_stages} with the \
         default stages). Numbers compared {compared}, largest difference relative to \
         max(1, |value|) {worst:.1e}",
        names.len()
    );
    assert!(
        failures.is_empty(),
        "{} of {} fixtures disagree:\n{}",
        failures.len(),
        names.len(),
        failures.join("\n")
    );
    println!(
        "regions: with the robust statistics of the step rule (SPEC D22) the step differs in {} \
         fixtures: {}",
        robust_differs.len(),
        robust_differs.join(" ")
    );
    assert!(with_triples >= 10 && default_stages >= 5 && compared > 50_000);
}

/// One fixture for `N` regions. Returns the number of triple points, whether the stages are
/// the default ones, and whether the step is another one with the robust statistics of the
/// step rule.
fn replay_regions(input: &Lines, out: &Lines, report: &mut Report) -> (usize, bool, bool) {
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
        step_statistics: StepStatistics::Extremes,
        ..base
    };

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
        return (cut.triple_points.len(), default, false);
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
    report.close(|| "max |rho|".into(), step.step_speed, float(&tau[2]));
    report.close(
        || "omega".into(),
        step.roughening_rate_max,
        float(&out["omega"][0][0]),
    );
    for w in &out["s"] {
        let p = int(&w[0]);
        report.close(
            || format!("omega at cut point {p}"),
            stages.roughening[p],
            float(&w[2]),
        );
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
            stages.rate[node],
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
    let robust = FlowOptions {
        step_statistics: StepStatistics::Robust,
        ..options
    };
    let mut other = ddg_variational_cuts::RegionFlow::new(&mesh, &geometry, &phis).unwrap();
    let again = other.step(&problem, &robust).unwrap();
    (
        cut.triple_points.len(),
        default,
        again.tau != step.tau || again.n_limited > 0,
    )
}
