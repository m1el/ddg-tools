//! Verification layer 3: the Lean oracle (`lean/oracle`), an independent dense implementation
//! of the canonical problem.
//!
//! * `fixtures_are_reproduced` replays the committed golden cases of `lean/oracle/fixtures`.
//! * `differential_fuzz_against_the_oracle_binary` generates random cases, runs them through
//!   the built `ddg-oracle` executable and compares. It is skipped, with a message, when the
//!   executable does not exist (`cd lean/oracle && lake build`, or set `DDG_ORACLE`).
//!
//! The text format is specified in `lean/oracle/README.md`: floats are the 16 hex digits of
//! their IEEE-754 bits, indices are 0-based.
//!
//! Run with `--nocapture` to see the tables.

mod common;

use std::path::{Path, PathBuf};
use std::process::Command;

use common::{FeatureMask, Features, N_KINDS, TestMesh, random_features, random_mesh, summary};
use ddg_core::operators::{
    angle_defects, boundary_geodesic_curvatures, cotan_laplacian, lumped_mass_normalized,
};
use ddg_core::{Geometry, TriMesh};
use ddg_optimal_cones::{Cone, ConeError, ConeProblem, Solution, SolveOptions, Status, solve};

/// Absolute agreement required of `u`, `φ`, `µ` and the objective (relative to
/// `max(1, |value|)`). The oracle's fixtures reproduce to `1e-9` across platforms.
const TOLERANCE: f64 = 1e-8;

fn hex(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn unhex(token: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(token, 16).unwrap_or_else(|_| panic!("bad float {token}")))
}

/// One input case of the oracle's text format.
#[derive(Debug, Clone)]
struct Case {
    name: String,
    positions: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    lambda: f64,
    /// Per vertex: `wE`, `wR`, `lo`, `hi`.
    weights: Vec<[f64; 4]>,
    cones: Vec<Cone>,
}

impl Case {
    fn from_test(name: &str, mesh: &TestMesh, f: &Features) -> Case {
        let nv = mesh.mesh.n_vertices();
        Case {
            name: name.to_string(),
            positions: mesh.positions.clone(),
            faces: mesh.mesh.faces().to_vec(),
            lambda: f.lambda,
            weights: (0..nv)
                .map(|v| {
                    [
                        f.distortion_weights.as_ref().map_or(1.0, |w| w[v]),
                        f.cone_costs.as_ref().map_or(1.0, |w| w[v]),
                        f.lower_bounds.as_ref().map_or(f64::NEG_INFINITY, |w| w[v]),
                        f.upper_bounds.as_ref().map_or(f64::INFINITY, |w| w[v]),
                    ]
                })
                .collect(),
            cones: f.prescribed_cones.clone(),
        }
    }

    fn features(&self) -> Features {
        let column = |k: usize| Some(self.weights.iter().map(|w| w[k]).collect::<Vec<f64>>());
        Features {
            lambda: self.lambda,
            distortion_weights: column(0),
            cone_costs: column(1),
            lower_bounds: column(2),
            upper_bounds: column(3),
            prescribed_cones: self.cones.clone(),
        }
    }

    fn write(&self, out: &mut String) {
        use std::fmt::Write;
        writeln!(out, "case {}", self.name).unwrap();
        writeln!(out, "mesh {} {}", self.positions.len(), self.faces.len()).unwrap();
        for p in &self.positions {
            writeln!(out, "v {} {} {}", hex(p[0]), hex(p[1]), hex(p[2])).unwrap();
        }
        for f in &self.faces {
            writeln!(out, "f {} {} {}", f[0], f[1], f[2]).unwrap();
        }
        writeln!(out, "lambda {}", hex(self.lambda)).unwrap();
        for w in &self.weights {
            writeln!(
                out,
                "w {} {} {} {}",
                hex(w[0]),
                hex(w[1]),
                hex(w[2]),
                hex(w[3])
            )
            .unwrap();
        }
        writeln!(out, "cones {}", self.cones.len()).unwrap();
        for c in &self.cones {
            writeln!(out, "c {} {}", c.vertex, hex(c.angle)).unwrap();
        }
        writeln!(out, "end").unwrap();
    }
}

fn content_lines(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .map(|line| line.split_whitespace().collect::<Vec<&str>>())
        .filter(|tokens| !tokens.is_empty() && !tokens[0].starts_with('#'))
}

/// Parses the input cases of a file. Vertex indices of faces are kept as written, so that an
/// invalid mesh reaches the mesh validation.
fn parse_cases(text: &str) -> Vec<Case> {
    let mut cases = Vec::new();
    let mut current: Option<Case> = None;
    for tokens in content_lines(text) {
        match tokens[0] {
            "case" => {
                current = Some(Case {
                    name: tokens[1].to_string(),
                    positions: Vec::new(),
                    faces: Vec::new(),
                    lambda: f64::NAN,
                    weights: Vec::new(),
                    cones: Vec::new(),
                })
            }
            "end" => cases.push(current.take().expect("end without case")),
            tag => {
                let case = current.as_mut().expect("line outside a case");
                match tag {
                    "v" => {
                        case.positions
                            .push([unhex(tokens[1]), unhex(tokens[2]), unhex(tokens[3])])
                    }
                    "f" => case.faces.push([
                        tokens[1].parse().unwrap(),
                        tokens[2].parse().unwrap(),
                        tokens[3].parse().unwrap(),
                    ]),
                    "lambda" => case.lambda = unhex(tokens[1]),
                    "w" => case.weights.push([
                        unhex(tokens[1]),
                        unhex(tokens[2]),
                        unhex(tokens[3]),
                        unhex(tokens[4]),
                    ]),
                    "c" => case.cones.push(Cone {
                        vertex: tokens[1].parse().unwrap(),
                        angle: unhex(tokens[2]),
                    }),
                    "mesh" | "cones" | "task" => {}
                    other => panic!("unknown tag {other}"),
                }
            }
        }
    }
    cases
}

#[derive(Debug, Clone, PartialEq)]
enum OracleStatus {
    Converged,
    NotConverged(String),
    Infeasible,
    Error(String),
}

/// One block of the oracle's `solve` output.
#[derive(Debug, Clone)]
struct OracleSolve {
    name: String,
    status: OracleStatus,
    /// Total, distortion, cone cost.
    objective: [f64; 3],
    phi_shift: (f64, f64),
    /// Per vertex: `u`, `φ`, `µ`.
    values: Vec<[f64; 3]>,
}

fn parse_solve_output(text: &str) -> Vec<OracleSolve> {
    let mut blocks = Vec::new();
    let mut current: Option<OracleSolve> = None;
    for tokens in content_lines(text) {
        match tokens[0] {
            "case" => {
                current = Some(OracleSolve {
                    name: tokens[1].to_string(),
                    status: OracleStatus::Error("no status line".to_string()),
                    objective: [f64::NAN; 3],
                    phi_shift: (0.0, 0.0),
                    values: Vec::new(),
                })
            }
            "end" => blocks.push(current.take().expect("end without case")),
            tag => {
                let block = current.as_mut().expect("line outside a block");
                match tag {
                    "status" => {
                        let message = tokens[2..].join(" ");
                        block.status = match tokens[1] {
                            "converged" => OracleStatus::Converged,
                            "not_converged" => OracleStatus::NotConverged(message),
                            "infeasible" => OracleStatus::Infeasible,
                            "error" => OracleStatus::Error(message),
                            other => panic!("unknown status {other}"),
                        }
                    }
                    "objective" => {
                        block.objective = [unhex(tokens[1]), unhex(tokens[2]), unhex(tokens[3])]
                    }
                    "phi_shift" => block.phi_shift = (unhex(tokens[1]), unhex(tokens[2])),
                    "s" => {
                        block
                            .values
                            .push([unhex(tokens[2]), unhex(tokens[3]), unhex(tokens[4])])
                    }
                    _ => {}
                }
            }
        }
    }
    blocks
}

/// What a comparison found.
#[derive(Debug, Default, Clone, Copy)]
struct Differences {
    objective: f64,
    measure: f64,
    log_scale: f64,
    adjoint: f64,
    /// Whether `u`, `µ`, `φ` were compared (the minimizer is unique).
    unique: bool,
}

enum Ours {
    Solved(Box<(TestMesh, Features, Solution)>),
    Infeasible,
    Error(String),
}

fn solve_case(case: &Case) -> Ours {
    let mesh = match TriMesh::new(case.positions.len(), case.faces.clone()) {
        Ok(mesh) => mesh,
        Err(e) => return Ours::Error(e.to_string()),
    };
    let geometry = match Geometry::from_positions(&mesh, &case.positions) {
        Ok(geometry) => geometry,
        Err(e) => return Ours::Error(e.to_string()),
    };
    let mesh = TestMesh {
        name: case.name.clone(),
        positions: case.positions.clone(),
        mesh,
        geometry,
    };
    let features = case.features();
    match solve(&features.problem(&mesh), &SolveOptions::default()) {
        Ok(solution) => Ours::Solved(Box::new((mesh, features, solution))),
        Err(ConeError::Infeasible { .. }) => Ours::Infeasible,
        Err(e) => Ours::Error(e.to_string()),
    }
}

/// Compares our result for a case with the oracle's block.
fn compare(case: &Case, oracle: &OracleSolve) -> Result<Differences, String> {
    assert_eq!(case.name, oracle.name);
    let ours = solve_case(case);
    let (mesh, features, solution) = match (&ours, &oracle.status) {
        (Ours::Error(_), OracleStatus::Error(_)) => return Ok(Differences::default()),
        (Ours::Infeasible, OracleStatus::Infeasible) => return Ok(Differences::default()),
        (Ours::Solved(solved), OracleStatus::Converged) => (&solved.0, &solved.1, &solved.2),
        (Ours::Solved(solved), OracleStatus::NotConverged(message)) => {
            // The oracle gave up; nothing to compare, but we must still be certified.
            return if solved.2.report.status == Status::Converged {
                Ok(Differences::default())
            } else {
                Err(format!(
                    "neither converged (oracle: {message}): {}",
                    summary(&solved.2)
                ))
            };
        }
        (ours, status) => {
            let ours = match ours {
                Ours::Solved(s) => summary(&s.2),
                Ours::Infeasible => "infeasible".to_string(),
                Ours::Error(e) => format!("error: {e}"),
            };
            return Err(format!("oracle says {status:?}, we say: {ours}"));
        }
    };
    if solution.report.status != Status::Converged {
        return Err(format!("not converged: {}", summary(solution)));
    }
    let nv = mesh.mesh.n_vertices();
    let relative = |a: f64, b: f64| (a - b).abs() / b.abs().max(1.0);
    let mut d = Differences {
        objective: relative(solution.objective, oracle.objective[0])
            .max(relative(solution.distortion_optimal, oracle.objective[1]))
            .max(relative(solution.cone_cost, oracle.objective[2])),
        unique: features.strictly_weighted(mesh),
        ..Differences::default()
    };
    if d.objective > TOLERANCE {
        return Err(format!(
            "objective (total, distortion, cone cost): ours ({}, {}, {}), oracle {:?}",
            solution.objective, solution.distortion_optimal, solution.cone_cost, oracle.objective
        ));
    }
    if !d.unique {
        // Only the objective and wE·u are determined.
        for v in 0..nv {
            let w = case.weights[v][0];
            d.log_scale = d
                .log_scale
                .max(w * (solution.log_scale_optimal[v] - oracle.values[v][0]).abs());
        }
        if d.log_scale > 1e-5 {
            return Err(format!("wE·u differs by {}", d.log_scale));
        }
        return Ok(d);
    }
    for v in 0..nv {
        d.log_scale = d
            .log_scale
            .max(relative(solution.log_scale_optimal[v], oracle.values[v][0]));
        d.measure = d
            .measure
            .max((solution.measure[v] - oracle.values[v][2]).abs());
    }
    // The adjoint: equal on a mesh with boundary; on a closed mesh equal up to a constant
    // that lies in the interval the oracle reports (and our interval is the same set).
    let shift = if mesh.is_closed() {
        solution.adjoint[0] - oracle.values[0][1]
    } else {
        0.0
    };
    for v in 0..nv {
        d.adjoint = d
            .adjoint
            .max((solution.adjoint[v] - shift - oracle.values[v][1]).abs());
    }
    if d.measure > TOLERANCE || d.log_scale > TOLERANCE || d.adjoint > TOLERANCE {
        return Err(format!(
            "differences: measure {:.3e}, u {:.3e}, adjoint {:.3e} (after the constant {shift:.3e})",
            d.measure, d.log_scale, d.adjoint
        ));
    }
    let (lo, hi) = oracle.phi_shift;
    if shift < lo - TOLERANCE || shift > hi + TOLERANCE {
        return Err(format!(
            "our adjoint is the oracle's plus {shift}, outside its admissible interval [{lo}, {hi}]"
        ));
    }
    let (our_lo, our_hi) = solution.adjoint_shift;
    if (our_lo + shift - lo).abs() > TOLERANCE || (our_hi + shift - hi).abs() > TOLERANCE {
        return Err(format!(
            "admissible adjoint constants: ours [{our_lo}, {our_hi}] + {shift}, oracle [{lo}, {hi}]"
        ));
    }
    Ok(d)
}

fn repository_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// Replays the committed fixtures: same status, objective, `u`, `µ` and `φ` as the oracle.
#[test]
fn fixtures_are_reproduced() {
    let dir = repository_root().join("lean/oracle/fixtures");
    let manifest = match std::fs::read_to_string(dir.join("MANIFEST.tsv")) {
        Ok(text) => text,
        Err(e) => {
            eprintln!(
                "SKIPPED: {} cannot be read ({e}); no fixture was checked",
                dir.join("MANIFEST.tsv").display()
            );
            return;
        }
    };
    println!(
        "{:<34} {:>4} {:<12} {:>10} {:>10} {:>10} {:>10}  stages I/II/safeguard",
        "fixture", "V", "status", "objective", "measure", "u", "adjoint"
    );
    let mut failures = Vec::new();
    let mut checked = 0;
    for line in manifest.lines().filter(|l| !l.trim().is_empty()) {
        let name = line.split('\t').next().unwrap().trim();
        let input = std::fs::read_to_string(dir.join(format!("{name}.in")))
            .unwrap_or_else(|e| panic!("{name}.in: {e}"));
        let output = std::fs::read_to_string(dir.join(format!("{name}.solve.out")))
            .unwrap_or_else(|e| panic!("{name}.solve.out: {e}"));
        let cases = parse_cases(&input);
        let blocks = parse_solve_output(&output);
        assert_eq!(
            (cases.len(), blocks.len()),
            (1, 1),
            "{name}: one case per fixture"
        );
        let (case, block) = (&cases[0], &blocks[0]);
        checked += 1;
        let status = match &block.status {
            OracleStatus::Converged => "converged",
            OracleStatus::NotConverged(_) => "not_converged",
            OracleStatus::Infeasible => "infeasible",
            OracleStatus::Error(_) => "error",
        };
        match compare(case, block) {
            Ok(d) => {
                let stages = match solve_case(case) {
                    Ours::Solved(s) => format!(
                        "{}/{}/{}",
                        s.2.report.regularized_iterations,
                        s.2.report.exact_iterations,
                        s.2.report.safeguard_rounds
                    ),
                    _ => "-".to_string(),
                };
                let cell = |x: f64, shown: bool| {
                    if shown {
                        format!("{x:.2e}")
                    } else {
                        "-".to_string()
                    }
                };
                let solved = block.status == OracleStatus::Converged;
                println!(
                    "{:<34} {:>4} {:<12} {:>10} {:>10} {:>10} {:>10}  {}",
                    name,
                    case.positions.len(),
                    status,
                    cell(d.objective, solved),
                    cell(d.measure, solved && d.unique),
                    cell(d.log_scale, solved),
                    cell(d.adjoint, solved && d.unique),
                    stages
                );
            }
            Err(e) => {
                println!("{name:<34} {:>4} {status:<12} FAILED", case.positions.len());
                failures.push(format!("{name}: {e}"));
            }
        }
    }
    assert!(checked > 0, "the manifest lists no fixture");
    assert!(
        failures.is_empty(),
        "{} of {checked} fixtures disagree:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn oracle_binary() -> Option<PathBuf> {
    let path = match std::env::var_os("DDG_ORACLE") {
        Some(path) => PathBuf::from(path),
        None => repository_root().join("lean/oracle/.lake/build/bin/ddg-oracle"),
    };
    path.is_file().then_some(path)
}

fn run_oracle(binary: &Path, command: &str, input: &str, tag: &str) -> String {
    let dir = std::env::temp_dir().join(format!("ddg-oracle-fuzz-{}-{tag}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (in_path, out_path) = (dir.join("cases.in"), dir.join("cases.out"));
    std::fs::write(&in_path, input).unwrap();
    let status = Command::new(binary)
        .arg(command)
        .arg(&in_path)
        .arg(&out_path)
        .status()
        .unwrap_or_else(|e| panic!("cannot run {}: {e}", binary.display()));
    assert!(status.success(), "ddg-oracle {command} failed: {status}");
    let output = std::fs::read_to_string(&out_path).unwrap();
    std::fs::remove_dir_all(&dir).ok();
    output
}

/// Random meshes and features through the oracle executable and through the solver.
/// `DDG_ORACLE_FUZZ_CASES` sets the number of cases (default 240).
#[test]
fn differential_fuzz_against_the_oracle_binary() {
    let Some(binary) = oracle_binary() else {
        eprintln!(
            "SKIPPED: the ddg-oracle executable was not found (build it with `lake build` in \
             lean/oracle, or set DDG_ORACLE); nothing was compared"
        );
        return;
    };
    let n_cases: u64 = std::env::var("DDG_ORACLE_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(240);
    let mut cases = Vec::new();
    let mut input = String::new();
    for i in 0..n_cases {
        let kind = (i as usize) % N_KINDS;
        let mesh = random_mesh(kind, 7000 + i);
        let mask = if i % 3 == 0 {
            FeatureMask::ALL
        } else {
            FeatureMask::UNIQUE
        };
        let features = random_features(&mesh, 31 * i + 5, mask);
        let case = Case::from_test(&format!("fuzz-{i}-{}", mesh.name), &mesh, &features);
        case.write(&mut input);
        cases.push(case);
    }
    let blocks = parse_solve_output(&run_oracle(&binary, "solve", &input, "solve"));
    assert_eq!(blocks.len(), cases.len(), "one output block per case");
    let mut failures = Vec::new();
    let mut worst = Differences::default();
    let (mut converged, mut infeasible, mut oracle_gave_up, mut unique) = (0, 0, 0, 0);
    for (case, block) in cases.iter().zip(&blocks) {
        match &block.status {
            OracleStatus::Converged => converged += 1,
            OracleStatus::Infeasible => infeasible += 1,
            OracleStatus::NotConverged(_) => oracle_gave_up += 1,
            OracleStatus::Error(e) => failures.push(format!("{}: oracle error: {e}", case.name)),
        }
        match compare(case, block) {
            Ok(d) => {
                unique += usize::from(d.unique);
                worst.objective = worst.objective.max(d.objective);
                worst.measure = worst.measure.max(d.measure);
                worst.adjoint = worst.adjoint.max(d.adjoint);
                if d.unique {
                    worst.log_scale = worst.log_scale.max(d.log_scale);
                }
            }
            Err(e) => failures.push(format!("{}: {e}", case.name)),
        }
    }
    println!(
        "oracle fuzz: {n_cases} cases, {converged} converged ({unique} with a unique minimizer), \
         {infeasible} infeasible, {oracle_gave_up} not converged by the oracle"
    );
    println!(
        "largest differences: objective {:.2e} (relative), measure {:.2e}, u {:.2e} (relative), \
         adjoint {:.2e}",
        worst.objective, worst.measure, worst.log_scale, worst.adjoint
    );
    assert!(
        failures.is_empty(),
        "{} of {n_cases} cases disagree, first ones:\n{}",
        failures.len(),
        failures
            .iter()
            .take(5)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// The oracle's `operators` output against the operators of `ddg-core` that the solver uses,
/// on the same random meshes. The oracle evaluates the literal formulas, `ddg-core` Kahan's
/// arrangements, so agreement is to rounding, not bit for bit.
#[test]
fn operators_agree_with_the_oracle_binary() {
    let Some(binary) = oracle_binary() else {
        eprintln!("SKIPPED: the ddg-oracle executable was not found; no operator was compared");
        return;
    };
    let mut meshes = Vec::new();
    let mut input = String::new();
    for i in 0..60u64 {
        let mesh = random_mesh((i as usize) % N_KINDS, 9000 + i);
        Case::from_test(&format!("ops-{i}"), &mesh, &Features::plain(0.5)).write(&mut input);
        meshes.push(mesh);
    }
    let output = run_oracle(&binary, "operators", &input, "operators");
    let mut index = usize::MAX;
    let mut worst = [0.0f64; 3];
    let mut laplacians = Vec::new();
    for mesh in &meshes {
        laplacians.push(cotan_laplacian(&mesh.mesh, &mesh.geometry).unwrap());
    }
    for tokens in content_lines(&output) {
        match tokens[0] {
            "case" => index = index.wrapping_add(1),
            "status" => assert_eq!(tokens[1], "ok", "{}", tokens.join(" ")),
            "o" => {
                let mesh = &meshes[index];
                let v: usize = tokens[1].parse().unwrap();
                let boundary = tokens[2] == "1";
                assert_eq!(boundary, mesh.mesh.is_boundary_vertex(v));
                let curvature = if boundary {
                    boundary_geodesic_curvatures(&mesh.mesh, &mesh.geometry).unwrap()[v]
                } else {
                    angle_defects(&mesh.mesh, &mesh.geometry).unwrap()[v]
                };
                let mass = lumped_mass_normalized(&mesh.mesh, &mesh.geometry).unwrap()[v];
                worst[0] = worst[0].max((curvature - unhex(tokens[3])).abs());
                worst[1] = worst[1].max((mass / unhex(tokens[4]) - 1.0).abs());
            }
            "l" => {
                let (i, j): (usize, usize) =
                    (tokens[1].parse().unwrap(), tokens[2].parse().unwrap());
                let ours = laplacians[index].get(i, j);
                let theirs = unhex(tokens[3]);
                worst[2] = worst[2].max((ours - theirs).abs() / theirs.abs().max(1.0));
            }
            _ => {}
        }
    }
    assert_eq!(index + 1, meshes.len());
    println!(
        "operators on {} meshes: largest difference of a curvature {:.2e}, of a mass {:.2e} \
         (relative), of a Laplacian entry {:.2e} (relative to max(1, |entry|))",
        meshes.len(),
        worst[0],
        worst[1],
        worst[2]
    );
    assert!(
        worst[0] < 1e-12 && worst[1] < 1e-12 && worst[2] < 1e-11,
        "{worst:?}"
    );
}

/// The problem builder used above must not drop anything: a case written and parsed back
/// gives the same problem.
#[test]
fn text_format_round_trip() {
    let mesh = random_mesh(3, 11);
    let features = random_features(&mesh, 5, FeatureMask::ALL);
    let case = Case::from_test("round-trip", &mesh, &features);
    let mut text = String::new();
    case.write(&mut text);
    let parsed = &parse_cases(&text)[0];
    assert_eq!(parsed.positions, case.positions);
    assert_eq!(parsed.faces, case.faces);
    assert_eq!(parsed.lambda.to_bits(), case.lambda.to_bits());
    assert_eq!(parsed.cones, case.cones);
    for (a, b) in parsed.weights.iter().zip(&case.weights) {
        assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    }
    let problem: ConeProblem<'_> = parsed.features().problem(&mesh);
    assert_eq!(problem.prescribed_cones, features.prescribed_cones);
}
