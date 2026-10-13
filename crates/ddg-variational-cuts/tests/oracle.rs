//! Differential test against the Lean oracle (`lean/oracle`, an independent dense implementation
//! written from `SPEC.md` alone): the committed fixtures `lean/oracle/fixtures/cuts/*.in|out`
//! are replayed, and everything the static evaluation computes is compared with the oracle's
//! output: the scaling, the face curvatures, the crossings, every triangle of the cut mesh, the
//! segments, `u`, `p`, `v`, `ℓ`, `q` at every node, `J` and the areas per region, and the energy.
//!
//! The gradient density `G` and the speed `s` depend on the flux recovery (SPEC §3.5). The
//! fixtures state none, so they are compared for both recoveries and must agree with one of
//! them: the paper's nodal quotient (what the SPEC said when the fixtures were generated) or
//! the smoothed default.
//!
//! Lines of the flow step (`tau`, `n`, `redist`, `energy_new`, `cuts_new`, column `SIGMA`) are
//! not read here.

use std::collections::HashMap;
use std::path::PathBuf;

use ddg_core::{Geometry, TriMesh};
use ddg_variational_cuts::{
    CurvatureSource, CutProblem, Energy, EvalOptions, Evaluation, Evaluator, FluxRecovery,
    Parallelism, RegionWeights, ScaleWeight, face_curvatures,
};

fn float(token: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(token, 16).unwrap_or_else(|_| panic!("bad float {token}")))
}

fn int(token: &str) -> usize {
    token
        .parse()
        .unwrap_or_else(|_| panic!("bad integer {token}"))
}

/// The lines of a block, grouped by their first word.
fn lines_by_key(text: &str) -> HashMap<String, Vec<Vec<String>>> {
    let mut out: HashMap<String, Vec<Vec<String>>> = HashMap::new();
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
    delta: f64,
    rho0: f64,
    passes: usize,
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
    Case {
        mesh,
        geometry,
        phi: floats(one("phi")),
        problem,
        feature,
        delta: float(&one("delta")[0]),
        rho0: float(&one("rho0")[0]),
        passes: int(&one("passes")[0]),
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

    fn close(&mut self, what: impl Fn() -> String, ours: f64, oracle: f64) {
        let error = (ours - oracle).abs() / oracle.abs().max(1.0);
        self.compared += 1;
        if error.is_nan() || error > self.tolerance {
            if self.problems.len() < 8 {
                self.problems
                    .push(format!("{}: {ours:e} vs oracle {oracle:e}", what()));
            } else if self.problems.len() == 8 {
                self.problems.push("...".to_string());
            }
        }
        if error.is_finite() {
            self.worst = self.worst.max(error);
        } else {
            self.worst = f64::INFINITY;
        }
    }

    fn same<T: PartialEq + std::fmt::Debug>(
        &mut self,
        what: impl Fn() -> String,
        ours: T,
        oracle: T,
    ) {
        if ours != oracle && self.problems.len() < 8 {
            self.problems
                .push(format!("{}: {ours:?} vs oracle {oracle:?}", what()));
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

/// Everything that does not depend on the flux recovery.
fn compare_evaluation(
    case: &Case,
    evaluator: &mut Evaluator,
    e: &Evaluation,
    out: &HashMap<String, Vec<Vec<String>>>,
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
                float(&w[4]),
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

/// The quantities that depend on the flux recovery: `G` at the cut vertices, `s`, and the total
/// energy (which does not depend on it, but is checked with the rest of its line).
fn compare_gradient(
    case: &Case,
    e: &Evaluation,
    out: &HashMap<String, Vec<Vec<String>>>,
    report: &mut Report,
) {
    let nv = case.mesh.n_vertices();
    for w in &out["b"] {
        let (r, n) = (int(&w[0]) as u32, int(&w[1]));
        if n < nv {
            continue; // the oracle prints G on the boundary of the mesh too; we do not keep it
        }
        let weights = case.problem.regions[r as usize];
        let oracle = weights.distortion * float(&w[5]) + weights.area;
        report.close(
            || format!("G at node {n} of region {r}"),
            e.gradient[our_vertex(e, r, n)],
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

#[test]
fn fixtures_of_the_oracle_are_reproduced() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../lean/oracle/fixtures/cuts");
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .filter_map(|entry| {
            let name = entry.unwrap().file_name().into_string().unwrap();
            name.strip_suffix(".in").map(str::to_string)
        })
        .collect();
    names.sort();
    assert!(
        names.len() >= 20,
        "only {} fixtures in {}",
        names.len(),
        dir.display()
    );
    let mut failures = Vec::new();
    let (mut nodal, mut smoothed) = (0, 0);
    let (mut worst, mut worst_gradient, mut compared) = (0.0f64, 0.0f64, 0usize);
    for name in &names {
        let case = parse_case(&std::fs::read_to_string(dir.join(format!("{name}.in"))).unwrap());
        let out = lines_by_key(&std::fs::read_to_string(dir.join(format!("{name}.out"))).unwrap());
        let status = out["status"][0][0].as_str();
        assert!(
            status == "ok" || status == "no_cut",
            "{name}: status {status}"
        );
        let mut evaluator = Evaluator::new(&case.mesh, &case.geometry).unwrap();
        let options = |flux_recovery: FluxRecovery| EvalOptions {
            min_crossing_offset: case.delta,
            flux_recovery,
            sliver_aspect_ratio: case.rho0,
            smoothing_passes: case.passes,
            parallelism: Parallelism::Sequential,
        };
        let e = evaluator
            .evaluate(&case.phi, &case.problem, &options(FluxRecovery::Nodal))
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(e.cut.points.is_empty(), status == "no_cut", "{name}");
        // The well-conditioned quantities, to the tolerance the oracle uses for its own check.
        let mut report = Report::new(1e-9);
        compare_evaluation(&case, &mut evaluator, &e, &out, &mut report);
        if !report.problems.is_empty() {
            failures.push(format!("{name}:\n    {}", report.problems.join("\n    ")));
            continue;
        }
        worst = worst.max(report.worst);
        compared += report.compared;
        // G and s divide by dual lengths that can be as small as delta times an edge, which
        // amplifies rounding errors; measured, the two implementations still agree to 3e-12.
        let mut as_nodal = Report::new(1e-9);
        compare_gradient(&case, &e, &out, &mut as_nodal);
        let default = evaluator
            .evaluate(&case.phi, &case.problem, &options(FluxRecovery::default()))
            .unwrap();
        let mut as_smoothed = Report::new(1e-9);
        compare_gradient(&case, &default, &out, &mut as_smoothed);
        if as_nodal.problems.is_empty() {
            nodal += 1;
            worst_gradient = worst_gradient.max(as_nodal.worst);
            compared += as_nodal.compared;
        } else if as_smoothed.problems.is_empty() {
            smoothed += 1;
            worst_gradient = worst_gradient.max(as_smoothed.worst);
            compared += as_smoothed.compared;
        } else {
            failures.push(format!(
                "{name}: G and s agree with neither recovery.\n  nodal:\n    {}\n  smoothed:\n    {}",
                as_nodal.problems.join("\n    "),
                as_smoothed.problems.join("\n    ")
            ));
        }
    }
    println!(
        "{} fixtures, {compared} numbers compared: largest difference {worst:.1e} (relative to \
         max(1, |value|)); G and s agree with the nodal recovery in {nodal} fixtures and with the \
         smoothed one in {smoothed}, largest difference {worst_gradient:.1e}",
        names.len()
    );
    assert!(compared > 20_000, "{compared} numbers compared");
    assert!(
        failures.is_empty(),
        "{} of {} fixtures disagree:\n{}",
        failures.len(),
        names.len(),
        failures.join("\n")
    );
    // The fixtures were generated with one recovery, not a mixture.
    assert!(
        nodal == 0 || smoothed == 0,
        "nodal {nodal}, smoothed {smoothed}"
    );
}
